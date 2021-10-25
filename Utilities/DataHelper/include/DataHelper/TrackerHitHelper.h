#ifndef TrackerHitHelper_H
#define TrackerHitHelper_H

#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitConst.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include <array>
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"


namespace CEPC{
  std::array<float, 6> GetCovMatrix(edm4hep::TrackerHit& hit, bool useSpacePointerBuilderMethod = false);
  float                GetResolutionRPhi(edm4hep::TrackerHit& hit);
  float                GetResolutionZ(edm4hep::TrackerHit& hit);
  std::array<float, 6> ConvertToCovXYZ(float dU, float thetaU, float phiU, float dV, float thetaV, float phiV, bool useSpacePointBuilderMethod = false);

  const edm4hep::ConstSimTrackerHit getAssoClosestSimTrackerHit(
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        const edm4hep::ConstTrackerHit trackerHit,
        const dd4hep::DDSegmentation::GridDriftChamber* segmentation,
        int docaMehtod);
}

#endif
