#ifndef DCHDigiSMAlg_h
#define DCHDigiSMAlg_h

#include "k4FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/NTuple.h"

#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitCollection.h"

#include "DetInterface/IGeomSvc.h"

#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"

#include "DD4hep/Detector.h"

#include "TVector3.h"

class GenfitField;
class IGeomSvc;

namespace edm4hep {
    class TrackCollection;
    class TrackerHitCollection;
    class MCRecoTrackerAssociationCollection;

}

class DCHDigiSMAlg : public GaudiAlgorithm
{

    public:
        DCHDigiSMAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize();
        virtual StatusCode execute();
        virtual StatusCode finalize();

    private:
        const GenfitField* m_genfitField;
        DataHandle<edm4hep::TrackCollection>  m_dcTrackCol{"DCTrackCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection>    m_digiDCHitsCol{"DigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_dcHitAssociationCol{"DCHitAssociationCollection",
                Gaudi::DataHandle::Reader, this};

        dd4hep::DDSegmentation::GridDriftChamber* m_segmentation;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
        Gaudi::Property<std::string> m_readout_name{ this, "readout", "DriftChamberHitsCollection"};//readout for getting segmentation

        int _nEvt;
        std::string m_thisName;
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::OverlayedField m_dd4hepField;

        NTuple::Tuple*  m_tuple;
        NTuple::Item<long> m_nEvt;
        NTuple::Item<long>   m_nTracks;
        NTuple::Item<long>   m_nTrackHits;
        NTuple::Item<int> m_nDCDigi;
        NTuple::Array<int> m_chamber;
        NTuple::Array<int> m_layer;
        NTuple::Array<int> m_cellID;
        NTuple::Array<double> m_dcHitTime;
        NTuple::Array<double> m_dcHitDoca;
        NTuple::Array<float> m_trackhitx;
        NTuple::Array<float> m_trackhity;
        NTuple::Array<float> m_trackhitz;
        NTuple::Array<float> m_distance;
        NTuple::Array<float> m_Doca;
        NTuple::Array<float> m_Doca_digi;
        NTuple::Array<float> m_d0;
        NTuple::Array<float> m_phi0;
        NTuple::Array<float> m_omega;
        NTuple::Array<float> m_z0;
        NTuple::Array<float> m_tanLambda;
        NTuple::Array<float> m_xc;
        NTuple::Array<float> m_yc;
        NTuple::Array<float> m_R;
        NTuple::Array<float> m_px;
        NTuple::Array<float> m_py;
        NTuple::Array<float> m_pz;
        NTuple::Array<float> m_pocax;
        NTuple::Array<float> m_pocay;

        Gaudi::Property<double> m_field{this, "Field", 3.0};


};

#endif
