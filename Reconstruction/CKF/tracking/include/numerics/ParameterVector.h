/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "Matrix.h"

namespace Belle2 {
  namespace TrackFindingCDC {

    /// Vector type for n related parameters
    template <int N>
    using ParameterVector = Matrix<double, N, 1>;

  }
}
