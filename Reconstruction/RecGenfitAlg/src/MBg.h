#ifndef MAKEBACKGROUND_h
#define MAKEBACKGROUND_h

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/NTuple.h"

#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCParticleConst.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"

#include "DetInterface/IGeomSvc.h"

#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"

#include "DD4hep/Detector.h"

#include "TVector3.h"
#include "TRandom3.h"

class GenfitField;
class IGeomSvc;

namespace edm4hep{
    class EventHeaderCollection;
    class MCParticleCollection;
    class SimTrackerHitCollection;
    class TrackCollection;
    class TrackerHitCollection;
    class MCRecoTrackerAssociationCollection;
    class ReconstructedParticle;
    class ReconstructedParticleCollection;
    class ConstTrackerHit;
    class ConstMCParticle;
    class Track;
}

class MakeBackground : public GaudiAlgorithm
{

    public:
        MakeBackground(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private:

        DataHandle<edm4hep::SimTrackerHitCollection> r_SimDCHCol{"DriftChamberHitsCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> r_DigiDCHCol{"DigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection> r_AssociationCol{"DCHitAssociationCollection", Gaudi::DataHandle::Reader, this};

        DataHandle<edm4hep::SimTrackerHitCollection> w_SimNoiseHCol{"NoiseSimHitsCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackerHitCollection> w_NoiseHitCol{"NoiseDCHitsCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>  w_NoiseAssociationCol{"NoiseDCHitAssociationCollection", Gaudi::DataHandle::Writer, this};

        const GenfitField* m_genfitField;


        int _nEvt;
        std::string m_thisName;
        SmartIF<IGeomSvc> m_geoSvc;
        dd4hep::OverlayedField m_dd4hepField;

    protected:

        
        Gaudi::Property<std::string> m_readout_name{ this, "readout",
            "DriftChamberHitsCollection"};
        Gaudi::Property<double> m_fHitPurity{this,"fHitPurity",0.1};
        Gaudi::Property<float> m_pocaTime  { this, "pocaTime", 225};// ns

        dd4hep::DDSegmentation::GridDriftChamber* m_segmentation;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

        TRandom3 fRandom;

};

#endif
