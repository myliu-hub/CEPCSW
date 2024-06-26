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
#include <string>

namespace Belle2 {

  /// A very rough filter for all CDC states.
  class RoughCDCStateFilter : public BaseCDCStateFilter {
  public:
    /// return 1 if distance < m_maximalHitDistance, NAN otherwise
    TrackFindingCDC::Weight operator()(const BaseCDCStateFilter::Object& pair) final;

    /// Expose the parameters of the sub findlets.
    //void exposeParameters( const std::string& prefix) override;

  private:
    /// maximal distance from track to trajectory (in XY)
    double m_maximalHitDistance = 500;
  };
}
