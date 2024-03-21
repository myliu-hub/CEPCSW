/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "Findlet.h"
#include "Algorithms.h"
#include "CDCCKFPath.h"

#include "CDCPathPairFilterFactory.h"
#include "ChooseableFilter.icc.h"
#include "StringManipulation.h"

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
        m_filter.setFilterName(m_pathSelectFilterName);

        m_filter.initialize();
        const auto pathComparison = [&](const CDCCKFPath & lhs, const CDCCKFPath & rhs) {
            return m_filter({&lhs, &rhs}) > 0;
      };
      std::sort(newPaths.begin(), newPaths.end(), pathComparison);

      TrackFindingCDC::only_best_N(newPaths, m_maximalCandidatesInFlight);
    }

    void SetPathSelectFilterName(std::string pathSelectFilterName)
    {
        m_pathSelectFilterName = pathSelectFilterName; 
    }
  private:
    std::string m_pathSelectFilterName;
    /// Maximum number of paths to select
    size_t m_maximalCandidatesInFlight = 4;
    /// Filter to order paths
    TrackFindingCDC::ChooseableFilter<CDCPathPairFilterFactory> m_filter;
  };
}
