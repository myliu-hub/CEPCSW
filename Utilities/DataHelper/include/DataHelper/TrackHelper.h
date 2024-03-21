#ifndef TRACKHELPER_H
#define TRACKHELPER_H
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackerHit.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "TMatrixDSym.h"
#include "TVector3.h"
#include "DataHelper/HelixClass.h"

namespace CEPC{
    //get track position and momentum from TrackState
    void getPosMomFromTrackState(const edm4hep::TrackState& trackState,
            double Bz, TVector3& pos,TVector3& mom,double& charge,
            TMatrixDSym& covMatrix_6);

    //Set track state from position, momentum and charge
    void getTrackStateFromPosMom(edm4hep::TrackState& trackState,double Bz,
            TVector3 pos,TVector3 mom,double charge,TMatrixDSym covMatrix_6);

    void getHelixFromPosMom(HelixClass& helix,double& xc,double& yc,
            double& R,double Bz,TVector3 seedPos,TVector3 seedMom,double charge);
    void getAssoMCParticle(
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            edm4hep::TrackerHit trackerHit,
            edm4hep::MCParticle& mcParticle,edm4hep::SimTrackerHit& simTrackerHit);
}

#endif
