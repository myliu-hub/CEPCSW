/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include <tracking/trackFindingCDC/findlets/base/Findlet.h>

#include <tracking/trackFindingCDC/numerics/EForwardBackward.h>
#include <tracking/ckf/general/utilities/SearchDirection.h>

#include <tracking/ckf/cdc/entities/CDCCKFState.h>
#include <tracking/ckf/cdc/entities/CDCCKFPath.h>

#include <tracking/trackFindingCDC/topology/CDCWire.h>
#include <tracking/trackFindingCDC/topology/CDCWireTopology.h>

#include <tracking/trackFindingCDC/utilities/StringManipulation.h>
#include <tracking/trackFindingCDC/numerics/Angle.h>


namespace Belle2 {

  /// Create CKF states, based on the current path. Perform some basic selection at this stage (based on phi, max. jump of layers)
  class CDCCKFStateCreator
    : public TrackFindingCDC::Findlet<CDCCKFState, const CDCCKFState,
      const TrackFindingCDC::CDCWireHit* const > {

    /// Parent class
    using Super = TrackFindingCDC::Findlet<CDCCKFState, const CDCCKFState, const TrackFindingCDC::CDCWireHit* const>;

    /// Store basic wire info for faster access
    struct CDCCKFWireHitCache {
      /// layer index
      int     icLayer;
      /// azimuthal coordinate
      double  phi;
    };


  public:

    /// Expose the parameters of the sub findlets.
    //void exposeParameters( const std::string& prefix) override
    //{
    //  moduleParamList->addParameter(TrackFindingCDC::prefixed(prefix, "maximalLayerJump"),
    //                                m_maximalLayerJump, "Maximal jump over N layers", m_maximalLayerJump);
    //  moduleParamList->addParameter(TrackFindingCDC::prefixed(prefix, "maximalLayerJumpBackwardSeed"),
    //                                m_maximalLayerJump_backwardSeed, "Maximal jump over N layers", m_maximalLayerJump_backwardSeed);
    //  moduleParamList->addParameter(TrackFindingCDC::prefixed(prefix, "maximalDeltaPhi"),
    //                                m_maximalDeltaPhi, "Maximal distance in phi between wires for Z=0 plane", m_maximalDeltaPhi);
    //  moduleParamList->addParameter(TrackFindingCDC::prefixed(prefix, "hitFindingDirection"),
    //                                m_param_writeOutDirectionAsString, "Start from innermost/outermost CDC layers", m_param_writeOutDirectionAsString);
    //}

    /// Clear the wireHit cache
    void beginEvent() override
    {
      Super::beginEvent();
      m_wireHitCache.clear();

      // Determine direction of track building
      m_param_writeOutDirection = fromString(m_param_writeOutDirectionAsString);

      if (m_param_writeOutDirection == TrackFindingCDC::EForwardBackward::c_Forward) {
        doForward = true;
      } else if (m_param_writeOutDirection == TrackFindingCDC::EForwardBackward::c_Backward) {
        doForward = false;
      } else {
        //B2FATAL("CDCCKFStateCreator: No valid direction specified. Please use forward/backward.");
      }
    }

    /// Main method of the findlet. Select + create states (output parameter nextStates) suitable for the input path, based on input wireHits
    void apply(std::vector<CDCCKFState>& nextStates, const CDCCKFPath& path,
               const std::vector<const TrackFindingCDC::CDCWireHit*>& wireHits) override
    {
      // TODO: as we do not need any information on the current state (track state) of the path, we could in principle
      // TODO: precalculate everything in here

      // Create cache over wirehits, if empty:
        std::cout << " wireHits size = " << wireHits.size() << std::endl;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << " m_wireHitCache size = " << m_wireHitCache.size() << std::endl;
      if (m_wireHitCache.empty()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        const size_t nHits = wireHits.size();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        m_wireHitCache.reserve(nHits);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        for (auto  hitPtr : wireHits) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << " hitPtr = " << hitPtr << std::endl;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          // to speed things up, don't consider background/taken hits at all (and not just in the loop below).
          // I can't just remove them from the list, otherwise the relation to the wireHits is broken
          // so set the layer index to a high number.
          if (hitPtr->getAutomatonCell().hasBackgroundFlag() || hitPtr->getAutomatonCell().hasTakenFlag()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            m_wireHitCache.push_back(CDCCKFWireHitCache{99999, 0.});
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          } else {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << " hitPtr->getWire() = " << hitPtr->getWire() << std::endl;
        std::cout << " hitPtr->getWire().getICLayer() = " << hitPtr->getWire().getICLayer() << std::endl;
        std::cout << " hitPtr->getRefPos2D().phi() = " << hitPtr->getRefPos2D().phi() << std::endl;
            m_wireHitCache.push_back(CDCCKFWireHitCache{hitPtr->getWire().getICLayer(), hitPtr->getRefPos2D().phi()});
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

      // Cache last-on-the-path state info too:
      const auto& lastState = path.back();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      double lastPhi = 0;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      double lastICLayer = -1;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      if (lastState.isSeed()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        if (doForward) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          lastICLayer = 0;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        } else {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          const auto& wireTopology = TrackFindingCDC::CDCWireTopology::getInstance();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          const auto& wires = wireTopology.getWires();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          const float maxForwardZ = wires.back().getForwardZ();     // 157.615
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          const float maxBackwardZ = wires.back().getBackwardZ();   // -72.0916
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

          const TrackFindingCDC::Vector3D seedPos(lastState.getSeed()->getPositionSeed());
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          const float seedPosZ = seedPos.z();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

          if (seedPosZ < maxForwardZ && seedPosZ > maxBackwardZ) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            lastICLayer = 56;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          } else {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            // do straight extrapolation of seed momentum to CDC outer walls
            TrackFindingCDC::Vector3D seedMomZOne(lastState.getSeed()->getMomentumSeed());
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            seedMomZOne = seedMomZOne / seedMomZOne.z();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            // const float maxZ = seedPosZ > 0 ? maxForwardZ : maxBackwardZ;
            // const TrackFindingCDC::Vector3D extrapolatedPos = seedPos - seedMom / seedMom.norm() * (seedPosZ - maxZ);

            // find closest iCLayer
            float minDist = 99999;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            for (const auto& wire : wires) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
              const float maxZ = seedPosZ > 0 ? wire.getForwardZ() : wire.getBackwardZ();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
              const TrackFindingCDC::Vector3D extrapolatedPos = seedPos - seedMomZOne * (seedPosZ - maxZ);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

              const auto distance = wire.getDistance(extrapolatedPos);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
              if (distance < minDist) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
                minDist = distance;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
                lastICLayer = wire.getICLayer();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
              }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            //B2DEBUG(29, lastICLayer << " (d=" << minDist << ")");
          }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      } else {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        lastPhi = lastState.getWireHit()->getRefPos2D().phi();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        lastICLayer = lastState.getWireHit()->getWire().getICLayer();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

      // Get sorted vector of wireHits on the path for faster search
      std::vector<const TrackFindingCDC::CDCWireHit*> wireHitsOnPath;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      for (auto const& state : path) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        if (! state.isSeed()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          wireHitsOnPath.push_back(state.getWireHit());
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      std::sort(wireHitsOnPath.begin(), wireHitsOnPath.end());
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

      size_t nHits = wireHits.size();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      for (size_t i = 0; i < nHits; i++) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        // adjust direction of loop (minimal speed gain)
        int idx = doForward ? i : nHits - i - 1;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        const auto iCLayer =  m_wireHitCache[idx].icLayer; // wireHit->getWire().getICLayer();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        if (m_param_writeOutDirection == TrackFindingCDC::EForwardBackward::c_Backward && lastState.isSeed()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          if (std::abs(lastICLayer - iCLayer) > m_maximalLayerJump_backwardSeed) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            continue;
          }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        } else if (std::abs(lastICLayer - iCLayer) > m_maximalLayerJump) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          continue;
        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        if (! lastState.isSeed()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          double deltaPhi = TrackFindingCDC::AngleUtil::normalised(lastPhi - m_wireHitCache[idx].phi);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          if (fabs(deltaPhi)  > m_maximalDeltaPhi)  {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
            continue;
          }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        const TrackFindingCDC::CDCWireHit* wireHit = wireHits[idx];
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        if (std::binary_search(wireHitsOnPath.begin(), wireHitsOnPath.end(), wireHit)) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
          continue;
        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        nextStates.emplace_back(wireHit);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
    }

  private:
    /// Maximum allowed step over layers
    int m_maximalLayerJump = 2;
    /// Maximum allowed step over layers (if outside->in CKF) for first step after seed (e.g. ECLShower)
    int m_maximalLayerJump_backwardSeed = 3;
    /// Maximal distance in phi between the path last hit/seed and the candidate hit
    double m_maximalDeltaPhi =  TMath::Pi() / 8;
    /// Parameter for the direction in which the tracks are built
    std::string m_param_writeOutDirectionAsString = "forward";
    /// Direction parameter converted from the string parameters
    TrackFindingCDC::EForwardBackward m_param_writeOutDirection = TrackFindingCDC::EForwardBackward::c_Unknown;
    /// Direction parameter converted to boolean for convenience
    bool doForward = true;

    /// Cache to store frequently used information
    std::vector<CDCCKFWireHitCache> m_wireHitCache = {};

  };
}
