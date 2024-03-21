/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "Filter.icc.h"
#include "CDCCKFPath.h"

#include <utility>

namespace Belle2 {
  /// Base filter for CKF CDC paths
  using BaseCDCPathPairFilter = TrackFindingCDC::Filter<std::pair<const CDCCKFPath*, const CDCCKFPath*>>;
}
