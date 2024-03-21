/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#include "NonIPCrossingStateFilter.dcl.h"

#include "CDCTrajectory2D.h"
#include "Vector3D.h"
#include "StringManipulation.h"

#include "SearchDirection.h"

#include "SpacePoint.h"
#include "RecoTrack.h"


namespace Belle2 {
  template <class AllStateFilter>
  TrackFindingCDC::Weight NonIPCrossingStateFilter<AllStateFilter>::operator()(const Object& pair)
  {
    if (std::isnan(AllStateFilter::operator()(pair))) {
      return NAN;
    }

    const auto& previousStates = pair.first;
    auto* state = pair.second;

    const RecoTrack* cdcTrack = previousStates.front()->getSeed();
    //B2ASSERT("Path without seed?", cdcTrack);

    const SpacePoint* spacePoint = state->getHit();
    //B2ASSERT("Path without hit?", spacePoint);

    const genfit::MeasuredStateOnPlane& firstMeasurement = [&state, &previousStates]() {
      if (state->mSoPSet()) {
        return state->getMeasuredStateOnPlane();
      } else {
        //B2ASSERT("Previous state was not fitted?", previousStates.back()->mSoPSet());
        return previousStates.back()->getMeasuredStateOnPlane();
      }
    }();

    const TrackFindingCDC::Vector3D& position = static_cast<TrackFindingCDC::Vector3D>(firstMeasurement.getPos());
    const TrackFindingCDC::Vector3D& momentum = static_cast<TrackFindingCDC::Vector3D>(firstMeasurement.getMom());

    const TrackFindingCDC::CDCTrajectory2D trajectory2D(position.xy(), 0, momentum.xy(), cdcTrack->getChargeSeed());

    const TrackFindingCDC::Vector3D& hitPosition = static_cast<TrackFindingCDC::Vector3D>(spacePoint->getPosition());
    const TrackFindingCDC::Vector2D origin(0, 0);

    const double deltaArcLengthHitOrigin = trajectory2D.calcArcLength2DBetween(hitPosition.xy(), origin);
    const double deltaArcLengthTrackHit = trajectory2D.calcArcLength2D(hitPosition.xy());

    if (not arcLengthInRightDirection(deltaArcLengthTrackHit, m_param_direction) or
        not arcLengthInRightDirection(deltaArcLengthHitOrigin, m_param_direction)) {
      return NAN;
    }

    return 1.0;
  }

//  template <class AllStateFilter>
//  void NonIPCrossingStateFilter<AllStateFilter>::exposeParameters( const std::string& prefix)
//  {
//    moduleParamList->addParameter(TrackFindingCDC::prefixed(prefix, "direction"), m_param_directionAsString,
//                                  "The direction where the extrapolation will happen.");
//  }

  template <class AllStateFilter>
  void NonIPCrossingStateFilter<AllStateFilter>::initialize()
  {
    Super::initialize();
    m_param_direction = fromString(m_param_directionAsString);
  }
}
