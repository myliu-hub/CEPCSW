/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "Findlet.h"

#include "CDCCKFStateCreator.h"
#include "CDCCKFStateFilter.h"
#include "CDCCKFPathMerger.h"
#include "CDCCKFPathSelector.h"

#include "CDCCKFState.h"
#include "CDCCKFPath.h"

#include "StringManipulation.h"


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

        if (paths.empty()) {
            return;
        }

        std::vector<CDCCKFPath> newPaths;
        std::vector<CDCCKFState> nextStates;

        for (CDCCKFPath& path : paths) {
            //根据jumpLayer排除比较远距离的击中
            m_stateCreator.beginEvent();
            m_stateCreator.setMaximalLayerJump(m_MaximalLayerJum);
            m_stateCreator.apply(nextStates, path, wireHits);
            m_stateFilter.apply(path, nextStates);

            // TODO: Attention: if there is no hit anymore, the path will not be added to the final set!
            for (const auto& nextState : nextStates) {

                path.push_back(nextState);

                newPaths.push_back(path);
                path.pop_back();
            }
            nextStates.clear();
        }

        // TODO:m_pathMerger not working
        m_pathMerger.apply(newPaths);

        m_pathSelector.SetPathSelectFilterName(m_pathSelectFilterName);
        m_pathSelector.apply(newPaths);

        if (newPaths.empty()) {
            return;
        }

        //交换两个vector的内容
        paths.swap(newPaths);
        newPaths.clear();

        apply(paths, wireHits);
    }

    void setPathSelectFilterName(std::string pathSelectFilterName)
    {
        m_pathSelectFilterName = pathSelectFilterName;
    }

    void setMaximalLayerJump(int MaximalLayerJum)
    {
        m_MaximalLayerJum = MaximalLayerJum;
    }
        private:
    std::string m_pathSelectFilterName;
    int m_MaximalLayerJum;
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
