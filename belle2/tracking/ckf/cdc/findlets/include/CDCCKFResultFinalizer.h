/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include <tracking/trackFindingCDC/findlets/base/Findlet.h>

#include <tracking/trackFindingCDC/filters/base/ChooseableFilter.h>
#include <tracking/ckf/cdc/filters/paths/CDCPathFilterFactory.h>

#include <tracking/trackFindingCDC/numerics/Weight.h>

#include <tracking/ckf/cdc/entities/CDCCKFPath.h>
#include <tracking/ckf/cdc/entities/CDCCKFResult.h>

#include <vector>

namespace Belle2 {

  /// Findlet to finalize CKF Paths in terms of final result.
  class CDCCKFResultFinalizer : public TrackFindingCDC::Findlet<const CDCCKFPath, CDCCKFResult> {
  public:
    CDCCKFResultFinalizer()
    {
      addProcessingSignalListener(&m_filter);
    }

    /// Expose the parameters of the sub findlets.
    //void exposeParameters( const std::string& prefix) override
    //{
    //  m_filter.exposeParameters(moduleParamList, prefix);
    //}

    /// main method of the findlet, for a list of paths create a list of results.
    void apply(const std::vector<CDCCKFPath>& paths, std::vector<CDCCKFResult>& results) override
    {
        std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        if (paths.empty()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl; 
            return;
        }

        std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        const CDCCKFPath* bestElement = nullptr;
        std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        TrackFindingCDC::Weight bestWeight = -NAN;
        std::cout << __FILE__ << " " << __LINE__ << std::endl; 


        for (const CDCCKFPath& path : paths) {
            std::cout << __FILE__ << " " << __LINE__ << std::endl; 
            m_filter.initialize();
            std::cout << __FILE__ << " " << __LINE__ << std::endl; 
            const TrackFindingCDC::Weight weight = m_filter(path);
            std::cout << __FILE__ << " " << __LINE__ << std::endl; 
            if (weight <= bestWeight) {
                std::cout << __FILE__ << " " << __LINE__ << std::endl; 
                continue;
            }
            std::cout << __FILE__ << " " << __LINE__ << std::endl; 
            bestWeight = weight;
            std::cout << __FILE__ << " " << __LINE__ << std::endl; 
            bestElement = &path;
            std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        }

        std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        if (bestElement and not std::isnan(bestWeight)) {
            std::cout << __FILE__ << " " << __LINE__ << std::endl; 
            results.push_back(*bestElement);
            std::cout << __FILE__ << " " << __LINE__ << std::endl; 
        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl; 
    }

      private:
    /// Filter to weigth the best path
    TrackFindingCDC::ChooseableFilter<CDCPathFilterFactory> m_filter;
  };
}
