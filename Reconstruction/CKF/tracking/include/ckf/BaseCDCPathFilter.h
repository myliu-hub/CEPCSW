/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "Filter.h"
#include "CDCCKFPath.h"

namespace Belle2 {
  /// Base filter for CKF CDC paths
  using BaseCDCPathFilter = TrackFindingCDC::Filter<CDCCKFPath>;
}
