#include "DataHelper/GeomeryHelper.h"
#include "Identifier/CEPCConf.h"
#include "DetInterface/IGeomSvc.h"


//Externals
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Segmentations.h"
#include "DDRec/ISurface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Vector3D.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "GaudiKernel/SmartIF.h"

#include "TMatrixF.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include <bitset>

// must initialize before using Get()!!!!
const GeomeryWire* GeomeryWire::Instance = nullptr;
const GeomeryWire* GeomeryWire::Get(){
    assert(Instance);
    return Instance;
}

GeomeryWire::GeomeryWire(SmartIF<IGeomSvc> geom):
    m_layerID(0),m_wireID(0),m_geomSvc(geom)
{
    m_dd4hepDetector=m_geomSvc->lcdd();
    dd4hep::Readout readout=m_dd4hepDetector->readout("DriftChamberHitsCollection");
    m_gridDriftChamber=dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>(readout.segmentation().segmentation());
    if(nullptr==m_gridDriftChamber){
        std::cout  << "Failed to get the GridDriftChamber" << std::endl;
    }
    Instance = this;
}

GeomeryWire::~GeomeryWire()
{
}

void GeomeryWire::setCellID(int layerID,int wireID)
{
    m_layerID=layerID;
    m_wireID=wireID;
}

void GeomeryWire::getWirePos(int layerID,int wireID,TVector3& Wstart,TVector3& Wend) const
{

    assert(m_gridDriftChamber);
    m_gridDriftChamber->cellposition2(0,layerID,wireID,Wstart,Wend);
}
