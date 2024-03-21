/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "Findlet.h"

#include "CDCCKFSeedCreator.h"
#include "StackTreeSearcher.h"
#include "CDCCKFResultFinalizer.h"
#include "CDCCKFResultStorer.h"

#include "CDCCKFPath.h"

#include "CDCWireHit.h"
#include "CDCWire.h"

#include <array>
#include <vector>
#include <iostream>
#include <string>

namespace Belle2 {
  class RecoTrack;


  /// Main findlet of the ToCDCCKF module
  class CKFToCDCFindlet : public TrackFindingCDC::Findlet<const TrackFindingCDC::CDCWireHit> {
    /// Parent class
    using Super = TrackFindingCDC::Findlet<const TrackFindingCDC::CDCWireHit>;

  public:
    /// Constructor, for setting module description and parameters.
    CKFToCDCFindlet();

    /// Default desctructor
    ~CKFToCDCFindlet() override;

    /// Expose the parameters of the sub findlets.
    //void exposeParameters( const std::string& prefix) override;

    double getPathArclength(std::vector<CDCCKFResult> results);
    double getSumDistances(std::vector<CDCCKFResult> results);
    double getChi2(std::vector<CDCCKFResult> results);

    /// Do the track/hit finding/merging.
    void apply(const std::vector<TrackFindingCDC::CDCWireHit>& wireHits) override;

    /// Clear the object pools
    void beginEvent() override;

    //TODO:二维数组初始化
    void getResult(std::vector< std::vector<unsigned short>> & output,
            std::vector<std::vector<double>>& trackParams,std::vector<double>& radius,
            std::vector<std::vector<double>>& circleCenter,std::vector<double>& deltaArcLength,
            std::vector<double>& hitDistance,std::vector<double>& doca,
            std::vector<double>& reconstructedZ,std::vector<TVector3>& recoPos){
        for(int i=0;i<m_results.size();i++){
            //for(const CDCCKFState & state: m_results[i]){

            int nState = 0;
            std::vector<double> arcLength;
            for(const auto& state: m_results[i]){

                if(state.isSeed()) continue;
                const auto* wireHit = state.getWireHit();
                const auto& wire = wireHit->getWire();

                unsigned short layerID = wire.getICLayer();
                unsigned short cellID = wire.getIWire();

                std::vector< unsigned short > ID;
                ID.clear();
                ID.push_back(layerID);
                ID.push_back(cellID);

                output.push_back(ID);
               
                TrackFindingCDC::PerigeeCircle perigeeCircle =
                    state.getTrajectory().getLocalHelix().helix().circleXY();
                TrackFindingCDC::SZLine szLine =
                    state.getTrajectory().getLocalHelix().helix().szLine();

                std::vector<double> trackParam;
                trackParam.push_back(perigeeCircle.impact());
                trackParam.push_back(perigeeCircle.phi0());
                trackParam.push_back(szLine.tanLambda());
                trackParam.push_back(perigeeCircle.curvature());
                trackParam.push_back(szLine.z0());

                trackParams.push_back(trackParam);
                trackParam.clear();

                radius.push_back(perigeeCircle.radius());

                std::vector<double> circleCenterXY;
                circleCenterXY.push_back(perigeeCircle.center().x());
                circleCenterXY.push_back(perigeeCircle.center().y());

                circleCenter.push_back(circleCenterXY);

                circleCenterXY.clear();
                arcLength.push_back(state.getArcLength());
                hitDistance.push_back(state.getHitDistance());
                doca.push_back(state.getDoca());
                reconstructedZ.push_back(state.getReconstructedZ());

                recoPos.push_back(state.getStatePos());
                nState++;
            }
            for(int i =1;i<arcLength.size();i++){
                deltaArcLength.push_back(arcLength[i] - arcLength[i-1]);
            }
        }
    }

    static void addSeedRecoTrack(RecoTrack* recoTrack){
        m_vxdRecoTrackVector.push_back(recoTrack);
    }

    static void clearSeedRecoTrack(){
        m_vxdRecoTrackVector.clear();
    }

    static void setpathSelectFilterName(std::string pathSelectFilterName)
    {
        m_pathSelectFilterName = pathSelectFilterName;
    }

    static void setMaximalLayerJump(int MaximalLayerJump)
    {
        m_MaximalLayerJump = MaximalLayerJump;
    }

    static void addMeasuredStateOnPlane(genfit::MeasuredStateOnPlane state)
    {
        m_state = state;
    }

  private:
    // Findlets
    /// Findlet for retrieving the vxd tracks and writing the result out
    //TrackLoader m_trackHandler;
    /// Seed Creator
    CDCCKFSeedCreator m_seedCreator;
    /// Tree Searcher
    StackTreeSearcher m_treeSearcher;
    /// Result Finalizer
    CDCCKFResultFinalizer m_resultFinalizer;
    /// Result Storer
    CDCCKFResultStorer m_resultStorer;

    // Object pools
    /// Pointers to the CDC Reco tracks as a vector
    static std::vector<RecoTrack*> m_vxdRecoTrackVector;

    static std::string m_pathSelectFilterName;
    static int m_MaximalLayerJump;
    static genfit::MeasuredStateOnPlane m_state;
    /// Current list of paths
    std::vector<CDCCKFPath> m_paths;
    /// Current list of seeds
    std::vector<CDCCKFPath> m_seeds;
    /// Current list of results
    std::vector<CDCCKFResult> m_results;
  };
}
