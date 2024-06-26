/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#include "SZParameters.h"

#include "UncertainParameters.icc.h"

using namespace Belle2;
using namespace TrackFindingCDC;

template struct TrackFindingCDC::UncertainParametersUtil<SZUtil, ESZParameter>;
