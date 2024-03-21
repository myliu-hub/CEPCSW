/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "Findlet.h"
#include "Vector2D.h"
#include "CDCTrajectory3D.h"
#include "CDCTrajectory2D.h"
#include "CDCTrajectorySZ.h"

#include "Algorithms.h"
#include "Functional.h"
#include "WeightComperator.h"

#include "CDCCKFState.h"
#include "CDCCKFPath.h"
#include "CDCStateFilterFactory.h"

#include "ChooseableFilter.h"

#include "StringManipulation.h"

namespace Belle2 {
  /// A stack of pre-, helix-extrapolation- , Kalman-extrapolation- and Kalman-update-filters.
  class CDCCKFStateFilter : public TrackFindingCDC::Findlet<const CDCCKFState, CDCCKFState> {
  public:
    /// Add all sub findlets
    CDCCKFStateFilter()
    {
      addProcessingSignalListener(&m_preFilter);
      addProcessingSignalListener(&m_basicFilter);
      addProcessingSignalListener(&m_extrapolationFilter);
      addProcessingSignalListener(&m_finalSelection);
    }



    /// Apply the findlet and do the state selection
    void apply(const CDCCKFPath& path, std::vector<CDCCKFState>& nextStates) override
    {
      const CDCCKFState& lastState = path.back();
      //if((!lastState.isSeed()) && (lastState.getWireHit())){
      //    std::cout << " CDCCKFStateFilter lastState layerID = " << lastState.getWireHit()->getWire().getILayer() << " cellID = " << lastState.getWireHit()->getWire().getIWire() << std::endl;
      //}
      //以lastState为起始点构建3D trajectory
      const TrackFindingCDC::CDCTrajectory3D& trajectory = lastState.getTrajectory();//此时构建的三维径迹以此时的State为圆心

      TrackFindingCDC::Weight weight;

      //B2DEBUG(29, "On layer: " << (lastState.isSeed() ? -1 : lastState.getWireHit()->getWire().getICLayer()));

      int iState =0;
      for (CDCCKFState& nextState : nextStates) {
        //B2DEBUG(29, "Checking layer: " << nextState.getWireHit()->getWire().getICLayer());

        iState++;

        m_preFilter.setFilterName("all");
        //m_preFilter.setFilterName("distance");
        m_preFilter.initialize();
        weight = m_preFilter({&path, &nextState});
        nextState.setWeight(weight);
        if (std::isnan(weight)) {
          std::cout << "CDCCKFStateFilter Fails PreFilter, because weight is nan!! " << std::endl;
          continue;
        }

        // Do a reconstruction based on the helix extrapolation from the last hit
        reconstruct(nextState, trajectory, lastState.getArcLength());
        //m_basicFilter.setFilterName("distance");
        m_basicFilter.setFilterName("rough");
        //m_basicFilter.setFilterName("all");
        m_basicFilter.initialize();
        weight = m_basicFilter({&path, &nextState});
        nextState.setWeight(weight);
        if (std::isnan(weight)) {
          std::cout << "CDCCKFStateFilter Fails Basic Filter, because weight is nan!! " << std::endl;
          continue;
        }

        // Extrapolate and update
        m_extrapolationFilter.setFilterName("extrapolate_and_update");
        m_extrapolationFilter.initialize();
        weight = m_extrapolationFilter({&path, &nextState});
        nextState.setWeight(weight);
        if (std::isnan(weight)) {
          std::cout << "CDCCKFStateFilter Fails Extrapolation, becaus weight is nan!! " << std::endl;
          continue;
        }

        // Do a final hit selection based on the new state
        const TrackFindingCDC::CDCTrajectory3D& thisTrajectory = nextState.getTrajectory();
        reconstruct(nextState, thisTrajectory, nextState.getArcLength());

        m_finalSelection.setFilterName("extrapolate_and_update");
        //m_finalSelection.setFilterName("distance");
        m_finalSelection.initialize();
        weight = m_finalSelection({&path, &nextState});
        nextState.setWeight(weight);
        if (std::isnan(weight)) {
            std::cout << "CDCCKFStateFilter Fails FinalFilter, because weight is nan!! " << std::endl;
            continue;
        }
      }

      //删除nextStates里面weight为NAN的state
      TrackFindingCDC::erase_remove_if(nextStates,
              TrackFindingCDC::Composition<TrackFindingCDC::IsNaN, TrackFindingCDC::GetWeight>());

      //根据weight对state由高到低进行排序
      std::sort(nextStates.begin(), nextStates.end(), TrackFindingCDC::GreaterWeight());

      //从nextStates中找出前m_maximalHitCandidates个state，其他的删除
      TrackFindingCDC::only_best_N(nextStates, m_maximalHitCandidates);
    }

  private:
    /// Parameter: max number of candidates
    //myliu:保留最大的候选击中
    size_t m_maximalHitCandidates = 10;
    /// Pre Filter
    TrackFindingCDC::ChooseableFilter<CDCStateFilterFactory> m_preFilter;
    /// Basic Filter (uses helix extrapolation)
    TrackFindingCDC::ChooseableFilter<CDCStateFilterFactory> m_basicFilter;
    /// Extrapolation Filter  (after Kalman extrapolation)
    TrackFindingCDC::ChooseableFilter<CDCStateFilterFactory> m_extrapolationFilter;
    /// Final Selection Filter (after Kalman update)
    TrackFindingCDC::ChooseableFilter<CDCStateFilterFactory> m_finalSelection;

    /// Helper function to reconstruct the arc length and the hit distance of a state according to the trajectory
    void reconstruct(CDCCKFState& state, const TrackFindingCDC::CDCTrajectory3D& trajectory, const double lastArcLength) const
    {
        // TODO: actually we do not need to do any trajectory creation here. We could save some computing time!
        const TrackFindingCDC::CDCTrajectory2D& trajectory2D = trajectory.getTrajectory2D(); //径迹在XY平面的投影圆
        const TrackFindingCDC::CDCTrajectorySZ& trajectorySZ = trajectory.getTrajectorySZ(); // 径迹沿Z方向的直线

        const TrackFindingCDC::CDCWireHit* wireHit = state.getWireHit();

        TrackFindingCDC::Vector2D recoPos2D;
        TrackFindingCDC::Vector2D refPos2D;//直斜丝有区别
        if (wireHit->isAxial()) {
            //recoPos2D为径迹上的最近点--point on track
            recoPos2D = wireHit->reconstruct2D(trajectory2D);
            refPos2D = wireHit->getRefPos2D();
        } else {
            const TrackFindingCDC::CDCWire& wire = wireHit->getWire();
            const TrackFindingCDC::Vector2D& posOnXYPlane = wireHit->reconstruct2D(trajectory2D);

            const double arcLength = trajectory2D.calcArcLength2D(posOnXYPlane);
            const double z = trajectorySZ.mapSToZ(arcLength);

            const TrackFindingCDC::Vector2D& wirePos2DAtZ = wire.getWirePos2DAtZ(z);

            const TrackFindingCDC::Vector2D& recoPosOnTrajectory = trajectory2D.getClosest(wirePos2DAtZ);
            const double driftLength = wireHit->getRefDriftLength();
            TrackFindingCDC::Vector2D disp2D = recoPosOnTrajectory - wirePos2DAtZ;
            disp2D.normalizeTo(driftLength);
            recoPos2D = wirePos2DAtZ + disp2D;
        }

        const double arcLength = trajectory2D.calcArcLength2D(recoPos2D);
        const double z = trajectorySZ.mapSToZ(arcLength);
        const double distanceToHit = trajectory2D.getDist2D(recoPos2D);
        const double doca = recoPos2D.distance(refPos2D);
        const double doca_test = recoPos2D.distance_VV(refPos2D);

        state.setArcLength(lastArcLength + arcLength);
        state.setHitDistance(distanceToHit); 
        state.setDoca(doca);
        state.setReconstructedZ(z);
    }
  };
}
