/* Copyright 2008-2010, Technische Universitaet Muenchen,
             2014, Ludwig-Maximilians-Universität München
   Authors: Tobias Schlüter

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
/** @addtogroup genfit
 * @{
 */

#include "WireMeasurementDC.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/TrackerHit.h"

#include <iostream>
#include <algorithm>

#include <Exception.h>
#include <RKTrackRep.h>
#include <HMatrixU.h>

#include <cassert>




WireMeasurementDC::WireMeasurementDC()
  : AbsMeasurement(1), maxDistance_(2), leftRight_(0)
{
  memset(wireEndPoint1_, 0, 3*sizeof(double));
  memset(wireEndPoint2_, 0, 3*sizeof(double));
}

WireMeasurementDC::WireMeasurementDC(double driftDistance, double driftDistanceError, const TVector3& endPoint1, const TVector3& endPoint2, int detId, int hitId, TrackPoint* trackPoint)
  : AbsMeasurement(1), maxDistance_(2), leftRight_(0)
{
  TVectorD coords(1);
  coords[0] = driftDistance;
  this->setRawHitCoords(coords);

  TMatrixDSym cov(1);
  cov(0,0) = driftDistanceError*driftDistanceError;
  this->setRawHitCov(cov);

  this->setWireEndPoints(endPoint1, endPoint2);

  this->setDetId(detId);
  this->setHitId(hitId);
  this->setTrackPoint(trackPoint);
}

SharedPlanePtr WireMeasurementDC::constructPlane(const StateOnPlane& state) const {

  // copy state. Neglect covariance.
  StateOnPlane st(state);

  TVector3 wire1(wireEndPoint1_);
  TVector3 wire2(wireEndPoint2_);

  // unit vector along the wire (V)
  TVector3 wireDirection = wire2 - wire1; 
  wireDirection.SetMag(1.);

  // point of closest approach
  const AbsTrackRep* rep = state.getRep();
  rep->extrapolateToLine(st, wire1, wireDirection);
  const TVector3& poca = rep->getPos(st);
  TVector3 dirInPoca = rep->getMom(st);
  dirInPoca.SetMag(1.);
  const TVector3& pocaOnWire = wire1 + wireDirection.Dot(poca - wire1)*wireDirection;

  // check if direction is parallel to wire
  if (fabs(wireDirection.Angle(dirInPoca)) < 0.01){
    Exception exc("WireMeasurementDC::detPlane(): Cannot construct detector plane, direction is parallel to wire", __LINE__,__FILE__);
    throw exc;
  }

  // construct orthogonal vector
  TVector3 U = wireDirection.Cross(dirInPoca);
  // U.SetMag(1.); automatically assured

  return SharedPlanePtr(new DetPlane(pocaOnWire, U, wireDirection));
}


std::vector<MeasurementOnPlane*> WireMeasurementDC::constructMeasurementsOnPlane(const StateOnPlane& state) const
{
  double mR = getRawHitCoords()(0);
  double mL = -mR;

  MeasurementOnPlane* mopL = new MeasurementOnPlane(TVectorD(1, &mL),
      getRawHitCov(),
      state.getPlane(), state.getRep(), constructHMatrix(state.getRep()));

  MeasurementOnPlane* mopR = new MeasurementOnPlane(TVectorD(1, &mR),
      getRawHitCov(),
      state.getPlane(), state.getRep(), constructHMatrix(state.getRep()));

  // set left/right weights
  if (leftRight_ < 0) {
    mopL->setWeight(1);
    mopR->setWeight(0);
  }
  else if (leftRight_ > 0) {
    mopL->setWeight(0);
    mopR->setWeight(1);
  }
  else {
    double val = 0.5 * pow(std::max(0., 1 - mR/maxDistance_), 2.);
    mopL->setWeight(val);
    mopR->setWeight(val);
  }

  std::vector<MeasurementOnPlane*> retVal;
  retVal.push_back(mopL);
  retVal.push_back(mopR);
  return retVal;
}

const AbsHMatrix* WireMeasurementDC::constructHMatrix(const AbsTrackRep* rep) const {
  if (dynamic_cast<const RKTrackRep*>(rep) == nullptr) {
    Exception exc("WireMeasurementDC default implementation can only handle state vectors of type RKTrackRep!", __LINE__,__FILE__);
    throw exc;
  }

  return new HMatrixU();
}

void WireMeasurementDC::setWireEndPoints(const TVector3& endPoint1, const TVector3& endPoint2)
{
  wireEndPoint1_[0] = endPoint1.X();
  wireEndPoint1_[1] = endPoint1.Y();
  wireEndPoint1_[2] = endPoint1.Z();

  wireEndPoint2_[0] = endPoint2.X();
  wireEndPoint2_[1] = endPoint2.Y();
  wireEndPoint2_[2] = endPoint2.Z();
}

void WireMeasurementDC::setLeftRightResolution(int lr){
  if (lr==0) leftRight_ = 0;
  else if (lr<0) leftRight_ = -1;
  else leftRight_ = 1;
}

void WireMeasurementDC::print(){
  std::cout<<"("<<layer_<<","<<cell_
    <<") wire(" <<wireEndPoint1_[0] <<wireEndPoint1_[1] <<wireEndPoint1_[2]
    <<wireEndPoint2_[0] <<wireEndPoint2_[1] <<wireEndPoint2_[2]
    <<")cm ptr: "<< trackerHit_
    << " dt "<<time_
    <<"ns dd "<<getRawHitCoords()[0]
    <<"cm dd err "<< getRawHitCov()(0,0)
    <<"cm lr "<<leftRight_
  <<std::endl;
  //<<" ddSm "<<driftDistanceSmeared
}


WireMeasurementDC::WireMeasurementDC(const GenfitHit* genfitHit,int iHit):
    genfitHit_(genfitHit),
    maxDistance_(genfitHit->getMaxDistance()),
    leftRight_(genfitHit->getLeftRightAmbig()),
    trackerHit_(genfitHit->getTrackerHit()),
    simTrackerHit_(genfitHit->getSimTrackerHit()),
    layer_(genfitHit->getLayer()),
    cell_(genfitHit->getCell())
{
    /// New a WireMeasurement
    try{
        new (this)WireMeasurementDC(
                genfitHit->getDriftDistance(),
                genfitHit->getDriftDistanceErr(),
                genfitHit->getEnd0(),genfitHit->getEnd1(),
                genfitHit->getCellID(),iHit,nullptr);
    }catch(genfit::Exception& e){
        std::cout<<"New WireMeasurementDC exception"<<std::endl;
        e.what();
    }
}

/** @} */
