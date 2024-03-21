/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#include "RoughCDCStateFilter.h"

#include "CDCCKFState.h"
#include "CDCCKFPath.h"

#include "StringManipulation.h"

using namespace Belle2;

TrackFindingCDC::Weight RoughCDCStateFilter::operator()(const BaseCDCStateFilter::Object& pair)
{
  const CDCCKFPath* path = pair.first;
  const CDCCKFState& state = *(pair.second);
  const CDCCKFState& lastState = path->back();

  const double& arcLength = state.getArcLength() - lastState.getArcLength();
  // TODO: magic number
  // TODO: myliu:set ArcLength range
  if (arcLength <= 0 or arcLength > 1000) {
    return NAN;
  }

  const double& hitDistance = state.getHitDistance();
  if (std::abs(hitDistance) > m_maximalHitDistance) {
    return NAN;
  }

  return 1;
}


//void RoughCDCStateFilter::exposeParameters( const std::string& prefix)
//{
//  moduleParamList->addParameter(TrackFindingCDC::prefixed(prefix, "maximalHitDistance"),
//                                m_maximalHitDistance,
//                                "Maximal allowed hit distance",
//                                m_maximalHitDistance);
//}
