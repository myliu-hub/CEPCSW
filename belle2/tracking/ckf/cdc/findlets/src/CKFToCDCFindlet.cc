/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/

#include <tracking/ckf/cdc/findlets/CKFToCDCFindlet.h>

#include <tracking/trackFindingCDC/utilities/Algorithms.h>


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


void CKFToCDCFindlet::apply(const std::vector<TrackFindingCDC::CDCWireHit>& wireHits)
{
    //m_trackHandler.apply(m_vxdRecoTrackVector);
    std::cout << " m_vxdRecoTrackVector size = " << m_vxdRecoTrackVector.size() << std::endl;
    std::cout << " m_vxdRecoTrackVector[0] = " << m_vxdRecoTrackVector[0] << std::endl;
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    m_seedCreator.apply(m_vxdRecoTrackVector, m_seeds);
    std::cout << " m_vxdRecoTrackVector size = " << m_vxdRecoTrackVector.size() << std::endl;
    std::cout << __FILE__ << " " << __LINE__ << std::endl; 

    const auto& wireHitPtrs = TrackFindingCDC::as_pointers<const TrackFindingCDC::CDCWireHit>(wireHits);
    std::cout << __FILE__ << " " << __LINE__ << std::endl; 

    for (const auto& seed : m_seeds) {
    std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        //B2DEBUG(29, "Starting new seed");
        m_paths.clear();
    std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        m_paths.push_back(seed);
    std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        m_treeSearcher.apply(m_paths, wireHitPtrs);
    std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        m_resultFinalizer.apply(m_paths, m_results);
    std::cout << __FILE__ << " " << __LINE__ << std::endl; 
    }

    std::cout << __FILE__ << " " << __LINE__ << std::endl; 
    m_resultStorer.initialize();
    std::cout << __FILE__ << " " << __LINE__ << std::endl; 
    m_resultStorer.apply(m_results);
    std::cout << __FILE__ << " " << __LINE__ << std::endl; 
}
