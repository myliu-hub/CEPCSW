#include "DetSegmentation/GridDriftChamber.h"
#include <map>

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
GridDriftChamber::GridDriftChamber(const std::string& cellEncoding) : Segmentation(cellEncoding) {
  // define type and description
  _type = "GridDriftChamber";
  _description = "Drift chamber segmentation in the global coordinates";

  registerParameter("cell_size", "cell size", m_cellSize, 0., SegmentationParameter::LengthUnit);
  registerParameter("detector_length", "Length of the wire", m_detectorLength, 1., SegmentationParameter::LengthUnit);
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "cellID");
  registerIdentifier("layerID", "layer id", layer_id, "layer");
  registerParameter("DC_inner_rmin", "DC_inner_rmin", m_DC_inner_rmin, 0., SegmentationParameter::LengthUnit);
  registerParameter("DC_inner_rmax", "DC_inner_rmax", m_DC_inner_rmax, 0., SegmentationParameter::LengthUnit);
  registerParameter("DC_outer_rmin", "DC_outer_rmin", m_DC_outer_rmin, 0., SegmentationParameter::LengthUnit);
  registerParameter("DC_outer_rmax", "DC_outer_rmax", m_DC_outer_rmax, 0., SegmentationParameter::LengthUnit);
}

GridDriftChamber::GridDriftChamber(const BitFieldCoder* decoder) : Segmentation(decoder) {

  _type = "GridDriftChamber";
  _description = "Drift chamber segmentation in the global coordinates";

  registerParameter("cell_size", "cell size", m_cellSize, 1., SegmentationParameter::LengthUnit);
  registerParameter("epsilon0", "epsilon", m_epsilon0, 0., SegmentationParameter::AngleUnit, true);
  registerParameter("detector_length", "Length of the wire", m_detectorLength, 1., SegmentationParameter::LengthUnit);
  registerIdentifier("identifier_phi", "Cell ID identifier for phi", m_phiID, "cellID");
  registerIdentifier("layerID", "layer id", layer_id, "layer");
  registerParameter("DC_inner_rmin", "DC_inner_rmin", m_DC_inner_rmin, 0., SegmentationParameter::LengthUnit);
  registerParameter("DC_inner_rmax", "DC_inner_rmax", m_DC_inner_rmax, 0., SegmentationParameter::LengthUnit);
  registerParameter("DC_outer_rmin", "DC_outer_rmin", m_DC_outer_rmin, 0., SegmentationParameter::LengthUnit);
  registerParameter("DC_outer_rmax", "DC_outer_rmax", m_DC_outer_rmax, 0., SegmentationParameter::LengthUnit);
}

Vector3D GridDriftChamber::position(const CellID& /*cID*/) const {
  Vector3D cellPosition = {0, 0, 0};
  return cellPosition;
}

CellID GridDriftChamber::cellID(const Vector3D& /*localPosition*/, const Vector3D& globalPosition,
                                const VolumeID& vID) const {

  CellID cID = vID;

  double posx = globalPosition.X;
  double posy = globalPosition.Y;
  double radius = sqrt(posx*posx+posy*posy);

  int layerid;
  if(radius > m_DC_inner_rmin && radius < m_DC_inner_rmax) layerid = floor(radius-m_DC_inner_rmin);
  else if (radius > m_DC_outer_rmin && radius < m_DC_outer_rmax)
           layerid = floor(radius-m_DC_outer_rmin)+floor(m_DC_inner_rmax-m_DC_inner_rmin);
  else return -1;

  updateParams(layerid);

  double phi_hit = phiFromXY(globalPosition);
  double offsetphi= m_offset;
  int _lphi;

  if(phi_hit >= offsetphi) {
    _lphi = (int) ((phi_hit - offsetphi)/ _currentLayerphi);
  }
  else {
    _lphi = (int) ((phi_hit - offsetphi + 2 * M_PI)/ _currentLayerphi);
  }

  int lphi = _lphi;
  _decoder->set(cID, layer_id, layerid);
  _decoder->set(cID, m_phiID, lphi);

  return cID;
}

double GridDriftChamber::phi(const CellID& cID) const {
  CellID phiValue = _decoder->get(cID, m_phiID);
  return binToPosition(phiValue, _currentLayerphi, m_offset);
}

void GridDriftChamber::cellposition(const CellID& cID, TVector3& Wstart,
                                    TVector3& Wend) const {

  auto layerIndex = _decoder->get(cID, "layer");
  updateParams(layerIndex);

  double phi_start = phi(cID);
  double phi_mid = phi_start + _currentLayerphi/2.;
  double phi_end = phi_mid + returnAlpha();

  Wstart = returnWirePosition(phi_mid, -1);
  Wend = returnWirePosition(phi_end, 1);
}



double GridDriftChamber::distanceTrackWire(const CellID& cID, const TVector3& hit_start,
                                           const TVector3& hit_end) const {

  TVector3 Wstart = {0,0,0};
  TVector3 Wend = {0,0,0};
  cellposition(cID,Wstart,Wend);

  TVector3 a = hit_end - hit_start;
  TVector3 b = Wend - Wstart;
  TVector3 c = Wstart - hit_start;

  double num = std::abs(c.Dot(a.Cross(b)));
  double denum = (a.Cross(b)).Mag();

  double DCA = 0;

   if (denum) {
    DCA = num / denum;
  }

  return DCA;
}


REGISTER_SEGMENTATION(GridDriftChamber)
}
}
