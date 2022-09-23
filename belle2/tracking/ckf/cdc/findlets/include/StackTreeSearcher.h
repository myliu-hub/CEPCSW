/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include <tracking/trackFindingCDC/findlets/base/Findlet.h>

#include <tracking/ckf/cdc/findlets/CDCCKFStateCreator.h>
#include <tracking/ckf/cdc/findlets/CDCCKFStateFilter.h>
#include <tracking/ckf/cdc/findlets/CDCCKFPathMerger.h>
#include <tracking/ckf/cdc/findlets/CDCCKFPathSelector.h>

#include <tracking/ckf/cdc/entities/CDCCKFState.h>
#include <tracking/ckf/cdc/entities/CDCCKFPath.h>

#include <tracking/trackFindingCDC/utilities/StringManipulation.h>


namespace Belle2 {
  /// CKF tree searcher which traces several best paths.
  class StackTreeSearcher : public
    TrackFindingCDC::Findlet<CDCCKFPath, const TrackFindingCDC::CDCWireHit* const> {
  public:
    StackTreeSearcher()
    {
      addProcessingSignalListener(&m_stateCreator);
      addProcessingSignalListener(&m_stateFilter);
      addProcessingSignalListener(&m_pathMerger);
      addProcessingSignalListener(&m_pathSelector);
    }

    /// Expose the parameters of the sub findlets.
    //void exposeParameters( const std::string& prefix) override
    //{
    //  m_stateCreator.exposeParameters(moduleParamList, prefix);
    //  m_stateFilter.exposeParameters(moduleParamList, TrackFindingCDC::prefixed("state", prefix));
    //  m_pathMerger.exposeParameters(moduleParamList, prefix);
    //  m_pathSelector.exposeParameters(moduleParamList, TrackFindingCDC::prefixed("path", prefix));
    //}

    /// Main method to update the paths. Input: vector of the selected paths and a vector of CDC wirehits to be considered.
    void apply(std::vector<CDCCKFPath>& paths,
               const std::vector<const TrackFindingCDC::CDCWireHit*>& wireHits) override
    {

        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        if (paths.empty()) {
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            return;
        }

        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        std::vector<CDCCKFPath> newPaths;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        std::vector<CDCCKFState> nextStates;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        for (CDCCKFPath& path : paths) {
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            //B2DEBUG(29, "Testing one path from " << path.back());
            m_stateCreator.apply(nextStates, path, wireHits);
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            m_stateFilter.apply(path, nextStates);
            std::cout << __FILE__ << " " << __LINE__ << std::endl;

            // TODO: Attention: if there is no hit anymore, the path will not be added to the final set!
            for (const auto& nextState : nextStates) {
                std::cout << __FILE__ << " " << __LINE__ << std::endl;
                path.push_back(nextState);
                std::cout << __FILE__ << " " << __LINE__ << std::endl;

                //B2DEBUG(29, "will go to " << nextState);
                newPaths.push_back(path);
                std::cout << __FILE__ << " " << __LINE__ << std::endl;
                path.pop_back();
                std::cout << __FILE__ << " " << __LINE__ << std::endl;
            }
            //B2DEBUG(29, "Now having " << newPaths.size() << " in flight");
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            nextStates.clear();
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        //B2DEBUG(29, "Having found " << newPaths.size() << " new paths");
        for (const auto& path : newPaths) {
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            //B2DEBUG(29, path);
        }

        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        m_pathMerger.apply(newPaths);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        //B2DEBUG(29, "Having found " << newPaths.size() << " new paths after merging");
        for (const auto& path : newPaths) {
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            //B2DEBUG(29, path);
        }

        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        m_pathSelector.apply(newPaths);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        //B2DEBUG(29, "Having found " << newPaths.size() << " new paths after selection");
        for (const auto& path : newPaths) {
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            //B2DEBUG(29, path);
        }

        if (newPaths.empty()) {
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            return;
        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        paths.swap(newPaths);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        newPaths.clear();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        apply(paths, wireHits);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
    }

        private:
    /// algorthim to create CDC-CDF states while traversing the path
    CDCCKFStateCreator m_stateCreator;
    /// algorithm to perform state filtering
    CDCCKFStateFilter m_stateFilter;
    /// algorithm to merge similar paths
    CDCCKFPathMerger m_pathMerger;
    /// algorithm to select N best paths, for further processing.
    CDCCKFPathSelector m_pathSelector;
    };
}
