/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)                           *
 * Author: The Belle II Collaboration                                     *
 *                                                                        *
 * See git log for contributors and copyright holders.                    *
 * This file is licensed under LGPL-3.0, see LICENSE.md.                  *
 **************************************************************************/
#include "CDCWireHit.h"

#include "CDCTrajectory2D.h"

#include "CDCWireSuperLayer.h"
#include "CDCWire.h"
#include "EStereoKind.h"

#include "Circle2D.h"
#include "Vector3D.h"
#include "Vector2D.h"

#include "ERightLeft.h"
#include "Index.h"

#include "TDCCountTranslatorBase.h"
#include "ADCCountTranslatorBase.h"

#include "CDCHit.h"
#include "WireID.h"

#include <ostream>

using namespace Belle2;
using namespace CDC;
using namespace TrackFindingCDC;

TDCCountTranslatorBase& CDCWireHit::getTDCCountTranslator()
{
//  static CDC::RealisticTDCCountTranslator s_tdcCountTranslator;
//  return s_tdcCountTranslator;
}

ADCCountTranslatorBase& CDCWireHit::getADCCountTranslator()
{
//  static CDC::LinearGlobalADCCountTranslator s_adcCountTranslator;
//  return s_adcCountTranslator;
}

CDCWireHit::CDCWireHit(const CDCHit* const ptrHit,
                       const double driftLength,
                       const double driftLengthVariance,
                       const double chargeDeposit,
                       const double driftTime)
  : m_wireID(ptrHit->getID())
  , m_wire(CDCWire::getInstance(*ptrHit))
  , m_automatonCell(1)
  , m_refDriftLength(driftLength)
  , m_refDriftLengthVariance(driftLengthVariance)
  , m_refChargeDeposit(chargeDeposit)
  , m_refDriftTime(driftTime)
  , m_hit(ptrHit)
{
}

//CDCWireHit::~CDCWireHit(){
//
////    assert(m_wire);
////    assert(m_hit);
//
// //   delete m_wire;
//    //delete m_hit;
//
//}

CDCWireHit::CDCWireHit(const CDCHit* const ptrHit,
                       TDCCountTranslatorBase* ptrTDCCountTranslator,
                       ADCCountTranslatorBase* ptrADCCountTranslator)
  : m_wireID(ptrHit->getID())
  , m_wire(ptrHit ? CDCWire::getInstance(*ptrHit) : nullptr)
  , m_automatonCell(1)
  , m_hit(ptrHit)
{
  if (not ptrHit) {
    //B2ERROR("CDCWireHit constructor invoked with nullptr CDCHit");
    return;
  }
  const CDCHit& hit = *ptrHit;

  TDCCountTranslatorBase& tdcCountTranslator =
    ptrTDCCountTranslator ? *ptrTDCCountTranslator : getTDCCountTranslator();
  ADCCountTranslatorBase& adcCountTranslator =
    ptrADCCountTranslator ? *ptrADCCountTranslator : getADCCountTranslator();

  float initialTOFEstimate = 0;

  float refDriftLengthRight = tdcCountTranslator.getDriftLength(hit.getTDCCount(),
                              getWireID(),
                              initialTOFEstimate,
                              false, // bool leftRight
                              getWire().getRefZ());

  float refDriftLengthLeft = tdcCountTranslator.getDriftLength(hit.getTDCCount(),
                             getWireID(),
                             initialTOFEstimate,
                             true, // bool leftRight
                             getWire().getRefZ());

  m_refDriftTime = tdcCountTranslator.getDriftTime(hit.getTDCCount(),
                                                   getWireID(),
                                                   initialTOFEstimate,
                                                   getWire().getRefZ(),
                                                   hit.getADCCount());

  m_refDriftLength = (refDriftLengthLeft + refDriftLengthRight) / 2.0;

  m_refDriftLengthVariance = tdcCountTranslator.getDriftLengthResolution(m_refDriftLength,
                             getWireID(),
                             false, // bool leftRight ?
                             getWire().getRefZ());

  m_refChargeDeposit = adcCountTranslator.getCharge(hit.getADCCount(),
                                                    getWireID(),
                                                    false, // bool leftRight
                                                    getWire().getRefZ(),
                                                    0); // theta
}

CDCWireHit::CDCWireHit(const WireID& wireID,
                       const double driftLength,
                       const double driftLengthVariance,
                       const double chargeDeposit)
  : m_wireID(wireID)
  , m_wire(CDCWire::getInstance(wireID))
  , m_automatonCell(1)
  , m_refDriftLength(driftLength)
  , m_refDriftLengthVariance(driftLengthVariance)
  , m_refChargeDeposit(chargeDeposit)
  , m_hit(nullptr)
{
}

bool CDCWireHit::operator<(const CDCHit& hit)
{
  return this->getWireID().getEWire() < hit.getID();
}

bool TrackFindingCDC::operator<(const CDCWireHit& wireHit, const CDCWireSuperLayer& wireSuperLayer)
{
  return wireHit.getISuperLayer() < wireSuperLayer.getISuperLayer();
}

bool TrackFindingCDC::operator<(const CDCWireSuperLayer& wireSuperLayer, const CDCWireHit& wireHit)
{
  return wireSuperLayer.getISuperLayer() < wireHit.getISuperLayer();
}

bool TrackFindingCDC::operator<(const CDCWireHit& wireHit, const CDCHit& hit)
{
  return wireHit.getWireID().getEWire() < hit.getID();
}

bool TrackFindingCDC::operator<(const CDCHit& hit, const CDCWireHit& wireHit)
{
  return hit.getID() < wireHit.getWireID().getEWire();
}

const CDCWire& CDCWireHit::attachWire() const
{
  m_wire = CDCWire::getInstance(m_wireID);
  assert(m_wire);
  return *m_wire;
}

Vector2D CDCWireHit::reconstruct2D(const CDCTrajectory2D& trajectory2D) const
{
  const Vector2D& refPos2D = getRefPos2D(); // the position of wirehit
  Vector2D recoPos2D = trajectory2D.getClosest(refPos2D);

  const Vector2D& wirePos2D = getWire().getRefPos2D(); // the position of wirehit
  const double driftLength = getRefDriftLength();


  Vector2D disp2D = recoPos2D - wirePos2D;

  // Fix the displacement to lie on the drift circle.
  disp2D.normalizeTo(driftLength);
  return wirePos2D + disp2D;
}

Vector3D CDCWireHit::reconstruct3D(const CDCTrajectory2D& trajectory2D,
                                   const ERightLeft rlInfo,
                                   const double z) const
{
  const EStereoKind stereoType = getStereoKind();

  if (stereoType == EStereoKind::c_StereoV or stereoType == EStereoKind::c_StereoU) {
    const WireLine& wireLine = getWire().getWireLine();
    const double signedDriftLength = isValid(rlInfo) ? rlInfo * getRefDriftLength() : 0.0;
    return trajectory2D.reconstruct3D(wireLine, signedDriftLength, z);

  } else { /*if (stereoType == EStereoKind::c_Axial)*/
    const Vector2D recoPos2D = reconstruct2D(trajectory2D);
    // for axial wire we can not determine the z coordinate by looking at the xy projection only
    // we set it the basic assumption.
    return Vector3D(recoPos2D, z);
  }
}

Circle2D CDCWireHit::conformalTransformed(const Vector2D& relativeTo) const
{
  Circle2D driftCircle(getRefPos2D() - relativeTo, getRefDriftLength());
  driftCircle.conformalTransform();
  return driftCircle;
}

//Index CDCWireHit::getStoreIHit() const
//{
//  return getHit() ? getHit()->getArrayIndex() : c_InvalidIndex;
//}

const Vector2D& CDCWireHit::getRefPos2D() const
{
  return getWire().getRefPos2D();
}

const Vector3D& CDCWireHit::getRefPos3D() const
{
  return getWire().getRefPos3D();
}

double CDCWireHit::getRefCylindricalR() const
{
  return getWire().getRefCylindricalR();
}

std::ostream& TrackFindingCDC::operator<<(std::ostream& output, const CDCWireHit& wirehit)
{
  return output << "CDCWireHit(" << wirehit.getWireID()
         << ", drift length=" << wirehit.getRefDriftLength() << ")";
}
