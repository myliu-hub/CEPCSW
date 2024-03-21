/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#include "DistanceCDCStateFilter.h"

#include "CDCCKFState.h"
#include "StringManipulation.h"

using namespace Belle2;

TrackFindingCDC::Weight DistanceCDCStateFilter::operator()(const BaseCDCStateFilter::Object& pair)
{
    const CDCCKFState& state = *(pair.second);

  double dist = std::abs(state.getHitDistance());

  if (dist > m_maximalDistance) {
    return NAN;
  }

  return 1 / dist;
}
