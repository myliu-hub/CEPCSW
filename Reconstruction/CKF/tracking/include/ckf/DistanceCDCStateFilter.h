/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "BaseCDCStateFilter.h"
#include "Weight.h"

namespace Belle2 {
  /// Give a weight based on the distance from the hit to the path
  class DistanceCDCStateFilter : public BaseCDCStateFilter {
  public:
    /// Return the weight based on the distance
    TrackFindingCDC::Weight operator()(const BaseCDCStateFilter::Object& pair) final;
    /// Expose the parameters
    //void exposeParameters( const std::string& prefix) override;
  private:
    /// Cut value for maximal distance
    double m_maximalDistance = 1000.;
    //double m_maximalDistance = 2.;
    //double m_maximalDistance = 0.5;
  };
}
