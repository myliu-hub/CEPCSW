/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/

#include "CKFToCDCFindlet.h"

#include "Algorithms.h"

using namespace Belle2;

CKFToCDCFindlet::~CKFToCDCFindlet() = default;

CKFToCDCFindlet::CKFToCDCFindlet()
{
  //addProcessingSignalListener(&m_trackHandler);
  addProcessingSignalListener(&m_seedCreator);
  addProcessingSignalListener(&m_treeSearcher);
  addProcessingSignalListener(&m_resultFinalizer);
  addProcessingSignalListener(&m_resultStorer);
}

//void CKFToCDCFindlet::exposeParameters( const std::string& prefix)
//{
//  Super::exposeParameters(moduleParamList, prefix);
//
//  //m_trackHandler.exposeParameters(moduleParamList, prefix);
//  m_seedCreator.exposeParameters(moduleParamList, prefix);
//  m_treeSearcher.exposeParameters(moduleParamList, prefix);
//  m_resultFinalizer.exposeParameters(moduleParamList, prefix);
//  m_resultStorer.exposeParameters(moduleParamList, prefix);
//
//  moduleParamList->getParameter<std::string>("statePreFilter").setDefaultValue("all");
//  moduleParamList->getParameter<std::string>("stateBasicFilter").setDefaultValue("rough");
//  moduleParamList->getParameter<std::string>("stateExtrapolationFilter").setDefaultValue("extrapolate_and_update");
//  moduleParamList->getParameter<std::string>("stateFinalFilter").setDefaultValue("distance");
//}

void CKFToCDCFindlet::beginEvent()
{
  Super::beginEvent();

  //m_vxdRecoTrackVector.clear();
  m_paths.clear();
  m_seeds.clear();
  m_results.clear();
}

std::vector<RecoTrack*> CKFToCDCFindlet::m_vxdRecoTrackVector;
std::string CKFToCDCFindlet::m_pathSelectFilterName;
int CKFToCDCFindlet::m_MaximalLayerJump;
genfit::MeasuredStateOnPlane CKFToCDCFindlet::m_state; 

double CKFToCDCFindlet::getPathArclength(std::vector<CDCCKFResult> results)
{
    if(results.size() > 1) return 1e-10;
    int num = results.size()-1;
    const auto& LastState = results[num].back();
    const auto ArcLength = LastState.getArcLength();
    return ArcLength;
}

double CKFToCDCFindlet::getSumDistances(std::vector<CDCCKFResult> results)
{
    if(results.size() > 1) return 1e-10;
    int n = results.size()-1;

    double sum = 0;
    int num =0;
    for (auto const& state : results[n]) {
        double dist = (state.getHitDistance());
        sum += dist * dist;
        num++;
    }
    return sum;
}

double CKFToCDCFindlet::getChi2(std::vector<CDCCKFResult> results)
{
    if(results.size() > 1) return 1e-10;
    int n = results.size()-1;

    return results[n].back().getChi2();
}

void CKFToCDCFindlet::apply(const std::vector<TrackFindingCDC::CDCWireHit>& wireHits)
{
    //m_trackHandler.apply(m_vxdRecoTrackVector);
    //从recoTrack构造CDCCKFPath(m_seeds)为seed
    m_seedCreator.addMeasuredStateOnPlane(m_state);
    m_seedCreator.apply(m_vxdRecoTrackVector, m_seeds);

    const auto& wireHitPtrs = TrackFindingCDC::as_pointers<const TrackFindingCDC::CDCWireHit>(wireHits);

    int i =0;
    //对seeds做循环
    for (const auto& seed : m_seeds) {
        //B2DEBUG(29, "Starting new seed");
        m_paths.clear();
        //Current list of paths,Play a transitional role in the middle
        m_paths.push_back(seed);
        //寻迹算法的核心
        m_treeSearcher.setPathSelectFilterName(m_pathSelectFilterName);
        m_treeSearcher.setMaximalLayerJump(m_MaximalLayerJump);
        m_treeSearcher.apply(m_paths, wireHitPtrs);
        m_resultFinalizer.apply(m_paths, m_results);
        i++;
    }

    //m_resultStorer.initialize();
    //m_resultStorer.apply(m_results);
    
    for(int i=0;i<m_results.size();i++){
        for(int j=0;j<m_results[i].size();j++){ 
        TrackFindingCDC::PerigeeCircle perigeeCircle = m_results[i][j].getTrajectory().getLocalHelix().helix().circleXY();
        TrackFindingCDC::SZLine szLine = m_results[i][j].getTrajectory().getLocalHelix().helix().szLine();
        }
    }


}
