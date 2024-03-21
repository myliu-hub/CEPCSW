/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#pragma once

#include "Findlet.h"

#include "CDCCKFPath.h"
#include "CDCCKFState.h"
#include "RecoTrack.h"

#include <vector>

namespace genfit {
    class Track;
    class AbsTrackRep;
}

namespace Belle2 {

  /// Create a CKF seed based on RecoTrack (presumably from VXDTF2)
  class CDCCKFSeedCreator : public TrackFindingCDC::Findlet<RecoTrack* const, CDCCKFPath> {
  public:
    /// Main method of the findlet, loop over reco tracks, create seeds for each of them.
    void apply(const std::vector<RecoTrack*>& recoTracks, std::vector<CDCCKFPath>& seeds) override
    {
        for (RecoTrack* recoTrack : recoTracks) {
            //CDCCKFState seedState(recoTrack, recoTrack->getMeasuredStateOnPlaneFromLastHit()); // last fitted hit
            //CDCCKFState seedState(recoTrack, recoTrack->getMeasuredStateOnPlaneFromFirstHit()); //first fitted hit
            TVector3 pos,mom;
            TMatrixDSym cov;
            m_measuredStateOnPlane.getPosMomCov(pos,mom,cov);
            CDCCKFState seedState(recoTrack, m_measuredStateOnPlane); //first fitted hit
            seeds.push_back({seedState});
        }
    }
    void addMeasuredStateOnPlane(genfit::MeasuredStateOnPlane state)
    {
        m_measuredStateOnPlane = state;
    }

  private:
    genfit::MeasuredStateOnPlane m_measuredStateOnPlane;
  };

}
