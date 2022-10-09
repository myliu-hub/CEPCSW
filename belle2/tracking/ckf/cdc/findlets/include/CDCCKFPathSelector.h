/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include <tracking/trackFindingCDC/findlets/base/Findlet.h>
#include <tracking/trackFindingCDC/utilities/Algorithms.h>
#include <tracking/ckf/cdc/entities/CDCCKFPath.h>

#include <tracking/ckf/cdc/filters/pathPairs/CDCPathPairFilterFactory.h>
#include <tracking/trackFindingCDC/filters/base/ChooseableFilter.icc.h>
#include <tracking/trackFindingCDC/utilities/StringManipulation.h>

namespace Belle2 {
  /// Select the m_maximalCandidatesInFlight paths for further processing
  class CDCCKFPathSelector : public TrackFindingCDC::Findlet<CDCCKFPath> {
  public:
    CDCCKFPathSelector()
    {
      addProcessingSignalListener(&m_filter);
    }

    /// Expose the parameters of the sub findlets.
    //void exposeParameters( const std::string& prefix) override
    //{
    //  moduleParamList->addParameter(TrackFindingCDC::prefixed(prefix, "maximalCandidatesInFlight"),
    //                                m_maximalCandidatesInFlight,
    //                                "Maximal candidates in flight", m_maximalCandidatesInFlight);
    //  m_filter.exposeParameters(moduleParamList, prefix);
    //}

    /// main method of the findlet, out of all paths "newPaths" select the best N=m_maximalCandidatesInFlight
    void apply(std::vector<CDCCKFPath>& newPaths) override
    {
        m_filter.setFilterName("arc_length");
        m_filter.initialize();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        const auto pathComparison = [&](const CDCCKFPath & lhs, const CDCCKFPath & rhs) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        return m_filter({&lhs, &rhs}) > 0;
      };
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      std::sort(newPaths.begin(), newPaths.end(), pathComparison);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

      TrackFindingCDC::only_best_N(newPaths, m_maximalCandidatesInFlight);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
    }

  private:
    /// Maximum number of paths to select
    size_t m_maximalCandidatesInFlight = 3;
    /// Filter to order paths
    TrackFindingCDC::ChooseableFilter<CDCPathPairFilterFactory> m_filter;
  };
}
