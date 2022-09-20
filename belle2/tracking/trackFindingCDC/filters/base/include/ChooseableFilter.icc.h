/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include <tracking/trackFindingCDC/filters/base/ChooseableFilter.dcl.h>

#include <tracking/trackFindingCDC/filters/base/FilterFactory.icc.h>


#include <tracking/trackFindingCDC/utilities/StringManipulation.h>

#include <memory>
#include <string>
#include <iostream>

namespace Belle2 {
  namespace TrackFindingCDC {

    template <class AFilter>
    Chooseable<AFilter>::Chooseable(std::unique_ptr<FilterFactory<AFilter>> filterFactory)
      : m_param_filterName(filterFactory ? filterFactory->getDefaultFilterName() : "")
      //, m_param_filterParameters()
      , m_filterFactory(std::move(filterFactory))
    {
      //B2ASSERT("Constructing a chooseable filter with no factory", m_filterFactory);
    }

    template <class AFilter>
    Chooseable<AFilter>::Chooseable(std::unique_ptr<FilterFactory<AFilter>> filterFactory,
                                    const std::string& filterName)
      : m_param_filterName(filterName)
      //, m_param_filterParameters()
      , m_filterFactory(std::move(filterFactory))
    {
      //B2ASSERT("Constructing a chooseable filter with no factory", m_filterFactory);
    }

    //template <class AFilter>
    //void Chooseable<AFilter>::exposeParameters( const std::string& prefix)
    //{
    //    m_filterFactory->createFiltersNameDescription();
      //  Super::exposeParameters(moduleParamList, prefix);
      //if (m_param_filterName == "") {
      //  /// Make a force parameter in case no default was given
      //  moduleParamList->addParameter(prefixed(prefix, "filter"),
      //                                m_param_filterName,
      //                                m_filterFactory->createFiltersNameDescription());
      //} else {
      //  /// Make a normal parameter in case default is known
      //  moduleParamList->addParameter(prefixed(prefix, "filter"),
      //                                m_param_filterName,
      //                                m_filterFactory->createFiltersNameDescription(),
      //                                m_param_filterName);
      //}

      //m_param_filterParameters.addParameter(moduleParamList,
      //                                      prefixed(prefix, "filterParameters"),
      //                                      m_filterFactory->createFiltersParametersDescription());
    //}

    template <class AFilter>
    void Chooseable<AFilter>::initialize()
    {
      m_filter = m_filterFactory->create(m_param_filterName);
      std::cout << __FILE__  << " if m_filter, m_param_filterName = " <<  m_param_filterName << std::endl;
      if (not m_filter) {
          std::cout << __FILE__  << " if no m_filter, m_param_filterName = " <<  m_param_filterName << std::endl;
        //B2ERROR("Could not create filter with name " << m_param_filterName);
        return;
      }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

      /// Transfer parameters
      //ModuleParamList filterModuleParamList;
      //const std::string prefix = "";
      //m_filter->exposeParameters(&filterModuleParamList, prefix);
      //m_param_filterParameters.assignTo(&filterModuleParamList);
      this->addProcessingSignalListener(m_filter.get());
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      Super::initialize();
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
    }

    template <class AFilter>
    bool Chooseable<AFilter>::needsTruthInformation()
    {
      return m_filter->needsTruthInformation();
    }

    template <class AFilter>
    Weight Chooseable<AFilter>::operator()(const Object& object)
    {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      return (*m_filter)(object);
    }

    template <class AFilter>
    Weight Chooseable<AFilter>::operator()(const Object& object) const
    {
      return (*m_filter)(object);
    }

    template <class AFilterFactory>
    ChooseableFilter<AFilterFactory>::ChooseableFilter()
      : Super(std::make_unique<AFilterFactory>())
    {
    }

    template <class AFilterFactory>
    ChooseableFilter<AFilterFactory>::ChooseableFilter(const std::string& filterName)
      : Super(std::make_unique<AFilterFactory>(), filterName)
    {
    }
  }
}
