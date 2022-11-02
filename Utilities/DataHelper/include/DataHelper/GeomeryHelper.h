#ifndef GeomeryHelper_H
#define GeomeryHelper_H


#include <array>
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "TVector3.h"

#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/Bootstrap.h"
#include "DetInterface/IGeomSvc.h"

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

        void getWirePos(int layerID,int wireID,TVector3& Wstart,TVector3& Wend) const;
        dd4hep::DDSegmentation::GridDriftChamber * getSeg() const;
    private:

        GeomeryWire();
        virtual ~GeomeryWire();

        GeomeryWire(const GeomeryWire&);
        GeomeryWire& operator=(const GeomeryWire&);

        SmartIF<IGeomSvc> m_geomsvc;
        dd4hep::Detector* m_dd4hepDetector;
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;

};

#endif
