/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#include "CDCStateFilterFactory.h"

#include "FilterFactory.icc.h"
#include "AllFilter.icc.h"
#include "NoneFilter.icc.h"
#include "AndFilter.icc.h"

#include "VariadicUnionVarSet.h"

#include "CDCStateBasicVarSet.h"
#include "RoughCDCStateFilter.h"
#include "ExtrapolateAndUpdateCDCStateFilter.h"
#include "DistanceCDCStateFilter.h"

using namespace Belle2;
using namespace TrackFindingCDC;

namespace {
  /// And filter for cdc states
  using AndCDCStateFilter = AndFilter<BaseCDCStateFilter>;
}


CDCStateFilterFactory::CDCStateFilterFactory(const std::string& defaultFilterName)
  : Super(defaultFilterName)
{
}

CDCStateFilterFactory::~CDCStateFilterFactory() = default;

std::string CDCStateFilterFactory::getIdentifier() const
{
  return "CDCState";
}

std::string CDCStateFilterFactory::getFilterPurpose() const
{
  return "Reject CDC CKF states. ";
}

std::map<std::string, std::string> CDCStateFilterFactory::getValidFilterNamesAndDescriptions() const
{
  return {
    {"none", "no track combination is valid"},
    {"all", "set all track combinations as good"},
    {"mc_truth", "filtering based on the mc truth information"},
    {"mc_truth_eclSeed", "filtering based on the mc truth information"},
    {"rough", "very rough filtering"},
    {"rough_eclSeed", "very rough filtering, seed created from ECL shower"},
    {"extrapolate_and_update", "Extrapolation and update"},
    {"distance", "Give a weight based on the distance"},
    {"recording", "record variables to a TTree"},
    {"recording_eclSeed", "record variables to a TTree"},
    {"rough_and_recording", "very rough filtering, seed created from SVD track"},
    {"rough_and_recording_eclSeed", "very rough filtering, seed created from ECL shower"},
    {"distance_and_recording_eclSeed", "Give a weight based on the distance"},

  };
}

std::unique_ptr<BaseCDCStateFilter>
CDCStateFilterFactory::create(const std::string& filterName) const
{
  if (filterName == "none") {
    return std::make_unique<TrackFindingCDC::NoneFilter<BaseCDCStateFilter>>();
  } else if (filterName == "all") {
    return std::make_unique<TrackFindingCDC::AllFilter<BaseCDCStateFilter>>();
  } else if (filterName == "rough") {
    return std::make_unique<RoughCDCStateFilter>();
  //} else if (filterName == "rough_eclSeed") {
  //  return std::make_unique<RoughCDCfromEclStateFilter>();
  //} else if (filterName == "mc_truth") {
  //  return std::make_unique<MCTruthCDCStateFilter>();
  //} else if (filterName == "mc_truth_eclSeed") {
  //  return std::make_unique<MCTruthEclSeedFilter>();
  } else if (filterName == "extrapolate_and_update") {
    return std::make_unique<ExtrapolateAndUpdateCDCStateFilter>();
  } else if (filterName == "distance") {
    return std::make_unique<DistanceCDCStateFilter>();
  //} else if (filterName == "recording") {
  //  return std::make_unique<RecordingCDCStateFilter>("CDCStateFilter.root");
  //} else if (filterName == "rough_and_recording") {
  //  return std::make_unique<AndCDCStateFilter>(
  //           std::make_unique<RecordingCDCStateFilter>("CDCStateFilter.root"),
  //           std::make_unique<RoughCDCStateFilter>()
  //         );
  //} else if (filterName == "recording_eclSeed") {
  //  return std::make_unique<RecordingCDCfromEclStateFilter>("CDCfromECLStateFilter.root");
  //} else if (filterName == "rough_and_recording_eclSeed") {
  //  return std::make_unique<AndCDCStateFilter>(
  //           std::make_unique<RecordingCDCfromEclStateFilter>("CDCfromECLStateFilter.root"),
  //           std::make_unique<RoughCDCfromEclStateFilter>()
  //         );
  //} else if (filterName == "distance_and_recording_eclSeed") {
  //  return std::make_unique<AndCDCStateFilter>(
  //           std::make_unique<RecordingCDCfromEclStateFilter>("CDCfromECLStateFilter.root"),
  //           std::make_unique<DistanceCDCStateFilter>()
  //         );
  } else {
    return Super::create(filterName);
  }
}
