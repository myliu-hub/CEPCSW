/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#include "RecoTrack.h"

#include "TrackCand.h"
#include "AbsTrackRep.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitStatus.h"
#include "WireTrackCandHit.h"
#include "RKTrackRep.h"
#include "MplTrackRep.h"

using namespace Belle2;

RecoTrack::RecoTrack(const TVector3& seedPosition, const TVector3& seedMomentum, const short int seedCharge,
                     const std::string& storeArrayNameOfCDCHits,
                     const std::string& storeArrayNameOfRecoHitInformation) :
  m_charge(seedCharge),
  m_storeArrayNameOfCDCHits(storeArrayNameOfCDCHits),
  m_storeArrayNameOfRecoHitInformation(storeArrayNameOfRecoHitInformation)
{
  m_genfitTrack.setStateSeed(seedPosition, seedMomentum);
  // TODO Set the covariance seed (that should be done by the tracking package)
  TMatrixDSym covSeed(6);
  covSeed(0, 0) = 1e-3;
  covSeed(1, 1) = 1e-3;
  covSeed(2, 2) = 4e-3;
  covSeed(3, 3) = 0.01e-3;
  covSeed(4, 4) = 0.01e-3;
  covSeed(5, 5) = 0.04e-3;
  m_genfitTrack.setCovSeed(covSeed);
}

genfit::TrackCand RecoTrack::createGenfitTrackCand() const
{
  genfit::TrackCand createdTrackCand;

  // Set the trajectory parameters
  createdTrackCand.setPosMomSeed(getPositionSeed(), getMomentumSeed(), getChargeSeed());
  createdTrackCand.setCovSeed(getSeedCovariance());
  createdTrackCand.setTimeSeed(getTimeSeed());
  createdTrackCand.sortHits();

  return createdTrackCand;
}

const genfit::TrackPoint* RecoTrack::getCreatedTrackPoint(const RecoHitInformation* recoHitInformation) const
{
  int createdTrackPointID = recoHitInformation->getCreatedTrackPointID();
  if (createdTrackPointID == -1) {
    return nullptr;
  }

  return m_genfitTrack.getPoint(createdTrackPointID);
}


bool RecoTrack::wasFitSuccessful(const genfit::AbsTrackRep* representation) const
{
  checkDirtyFlag();

  if (getRepresentations().empty()) {
    return false;
  }

  if (not hasTrackFitStatus(representation)) {
    return false;
  }

  const genfit::FitStatus* fs = getTrackFitStatus(representation);
  if (not fs) {
    return false;
  }
  if (not fs->isFitConverged()) {
    return false;
  }

  // make sure we only consider fitted if the Kalman method was used
  if (not dynamic_cast<const genfit::KalmanFitStatus*>(fs)) {
    return false;
  }

  // make sure there is at least one hit with a valid mSoP
  const unsigned int trackSize = m_genfitTrack.getNumPoints();
  for (unsigned int i = 0; i < trackSize; i++) {
    try {
      m_genfitTrack.getFittedState(i, representation);
      return true;
    } catch (const genfit::Exception& exception) {
      //B2DEBUG(100, "Can not get mSoP because of: " << exception.what());
    }
  }

  return false;
}

genfit::Track& RecoTrackGenfitAccess::getGenfitTrack(RecoTrack& recoTrack)
{
  return recoTrack.m_genfitTrack;
}

bool RecoTrackGenfitAccess::InsertTrackPoint(RecoTrack& recoTrack,genfit::TrackPoint* trackPoint)
{

    if(nullptr==trackPoint) return false;
    RecoTrackGenfitAccess::getGenfitTrack(recoTrack).insertPoint(trackPoint);
    return true;
}

genfit::AbsTrackRep* RecoTrackGenfitAccess::createOrReturnRKTrackRep(RecoTrack& recoTrack, int PDGcode)
{
  // try to get the trackRep, if it has already been added
  genfit::AbsTrackRep* trackRepresentation = recoTrack.getTrackRepresentationForPDG(std::abs(PDGcode));

  // not available? create one
  trackRepresentation = new genfit::RKTrackRep(PDGcode);
  RecoTrackGenfitAccess::getGenfitTrack(recoTrack).addTrackRep(trackRepresentation);
  return trackRepresentation;
}

bool RecoTrackGenfitAccess::AddTrackRep(RecoTrack& recoTrack, genfit::AbsTrackRep* rep)
{
    if (rep != nullptr) {
        RecoTrackGenfitAccess::getGenfitTrack(recoTrack).addTrackRep(rep);
        return true;
    }
    return false;
}

const genfit::MeasuredStateOnPlane& RecoTrack::getMeasuredStateOnPlaneClosestTo(const TVector3& closestPoint,
    const genfit::AbsTrackRep* representation)
{
  checkDirtyFlag();
  const unsigned int numberOfPoints = m_genfitTrack.getNumPointsWithMeasurement();

  assert(numberOfPoints > 0);

  const genfit::MeasuredStateOnPlane* nearestStateOnPlane = nullptr;
  double minimalDistance2 = 0;
  for (unsigned int hitIndex = 0; hitIndex < numberOfPoints; hitIndex++) {
    try {
      const genfit::MeasuredStateOnPlane& measuredStateOnPlane = m_genfitTrack.getFittedState(hitIndex, representation);

      const double currentDistance2 = (measuredStateOnPlane.getPos() - closestPoint).Mag2();

      if (not nearestStateOnPlane or currentDistance2 < minimalDistance2) {
        nearestStateOnPlane = &measuredStateOnPlane;
        minimalDistance2 = currentDistance2;
      }
    } catch (const genfit::Exception& exception) {
      //B2DEBUG(50, "Can not get mSoP because of: " << exception.what());
      continue;
    }
  }
  return *nearestStateOnPlane;
}


void RecoTrack::deleteFittedInformation()
{
  // Delete all fitted information for all representations
  for (const genfit::AbsTrackRep* rep : getRepresentations()) {
    deleteFittedInformationForRepresentation(rep);
  }
}

void RecoTrack::deleteFittedInformationForRepresentation(const genfit::AbsTrackRep* rep)
{
  m_genfitTrack.deleteFittedState(rep);
}

genfit::AbsTrackRep* RecoTrack::getTrackRepresentationForPDG(int pdgCode)
{

  const std::vector<genfit::AbsTrackRep*>& trackRepresentations = getRepresentations();

  for (genfit::AbsTrackRep* trackRepresentation : trackRepresentations) {
    // Check if the track representation is a RKTrackRep.
    const genfit::RKTrackRep* rkTrackRepresenation = dynamic_cast<const genfit::RKTrackRep*>(trackRepresentation);
    if (rkTrackRepresenation != nullptr) {
      // take the aboslute value of the PDG code as the TrackRep holds the PDG code including the charge (so -13 or 13)
      if (std::abs(rkTrackRepresenation->getPDG()) == pdgCode) {
        return trackRepresentation;
      }
    }
  }

  return nullptr;
}


/// Helper function to get the seed or the measured state on plane from a track
std::tuple<TVector3, TVector3, short> RecoTrack::extractTrackState() const
{
  if (not wasFitSuccessful()) {
    return std::make_tuple(getPositionSeed(), getMomentumSeed(), getChargeSeed());
  } else {
    const auto& measuredStateOnPlane = getMeasuredStateOnPlaneFromFirstHit();
    return std::make_tuple(measuredStateOnPlane.getPos(), measuredStateOnPlane.getMom(), measuredStateOnPlane.getCharge());
  }
}

bool RecoTrack::hasTrackFitStatus(const genfit::AbsTrackRep* representation) const
{
  checkDirtyFlag();

  // there might be the case, where the genfit track has no trackreps, even not the cardinal
  // one because no fit attempt was performed. In this case, the "hasFitStatus" call to genfit
  // will fail with an access violation. To prevent that, check for the number of reps here before
  // actually calling genfit's hasFitStatus(...)
  if (m_genfitTrack.getNumReps() == 0)
    return false;

  return m_genfitTrack.hasFitStatus(representation);
}


const genfit::MeasuredStateOnPlane& RecoTrack::getMeasuredStateOnPlaneFromRecoHit(const RecoHitInformation* recoHitInfo,
    const genfit::AbsTrackRep* representation) const
{
  checkDirtyFlag();

  const auto* hitTrackPoint = getCreatedTrackPoint(recoHitInfo);

  const auto* fittedResult = hitTrackPoint->getFitterInfo(representation);

  return fittedResult->getFittedState();
}

const genfit::MeasuredStateOnPlane& RecoTrack::getMeasuredStateOnPlaneFromFirstHit(const genfit::AbsTrackRep* representation) const
{
  const unsigned int trackSize = m_genfitTrack.getNumPoints();
  for (unsigned int i = 0; i < trackSize; i++) {
    try {
      return m_genfitTrack.getFittedState(i, representation);
    } catch (const genfit::Exception& exception) {
      //B2DEBUG(50, "Can not get mSoP because of: " << exception.what());
    }
  }
}

const genfit::MeasuredStateOnPlane& RecoTrack::getMeasuredStateOnPlaneFromLastHit(const genfit::AbsTrackRep* representation) const
{
  int trackSize = m_genfitTrack.getNumPoints();
  for (int i = -1; i >= -1*trackSize; i--) {
  //for (int i = 0; i <= trackSize; i++) {
    try {
  //representation->Print();
      return m_genfitTrack.getFittedState(i, representation);
    } catch (const genfit::Exception& exception) {
        //B2DEBUG(50, "Can not get mSoP because of: " << exception.what());
    }
  }
}

const genfit::MeasuredStateOnPlane& RecoTrack::getMeasuredStateOnPlaneExtrapolation(int id) const
{
    int repSize = m_genfitTrack.getTrackReps().size();
    genfit::AbsTrackRep* rep = m_genfitTrack.getTrackRep(id);
    
    return genfit::MeasuredStateOnPlane(rep);
}
