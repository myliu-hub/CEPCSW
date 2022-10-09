#ifndef GeomeryHelper_H
#define GeomeryHelper_H


#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitConst.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include <array>
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "TVector3.h"
#include "GaudiKernel/SmartIF.h"

class IGeomSvc;
namespace dd4hep {
    class Detector;
    namespace DDSegmentation{
        class GridDriftChamber;
    }
}
namespace dd4hep {
    class Detector;
    namespace DDSegmentation{
        class GridDriftChamber;
    }
    namespace rec{
        class ISurface;
    }
}

class GeomeryWire{

    public:
        static const GeomeryWire* Instance;
        static const GeomeryWire* Get();
        GeomeryWire(SmartIF<IGeomSvc> geom);
        virtual ~GeomeryWire();

        void setCellID(int layerID,int wireID);
        void getWirePos(int layerID,int wireID,TVector3& Wstart,TVector3& Wend) const;

    private:

    SmartIF<IGeomSvc> m_geomSvc;

    dd4hep::Detector* m_dd4hepDetector;
    dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
    int m_layerID;
    int m_wireID;

};

#endif
