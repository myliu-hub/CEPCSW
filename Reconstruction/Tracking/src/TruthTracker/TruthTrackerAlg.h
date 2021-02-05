#ifndef TruthTrackerAlg_h
#define TruthTrackerAlg_h

#include "GaudiAlg/GaudiAlgorithm.h"
#include "k4FWCore/DataHandle.h"
#include "DD4hep/Fields.h"
#include "GaudiKernel/NTuple.h"

class IGeomSvc;
namespace dd4hep {
    class Detector;
    namespace DDSegmentation{
        class GridDriftChamber;
        class BitFieldCoder;
    }
}
namespace edm4hep {
    class MCParticleCollection;
    class TrackerHitCollection;
    class SimTrackerHitCollection;
    class TrackCollection;
    class MCRecoTrackerAssociationCollection;
    class ReconstructedParticleCollection;
    class MCRecoParticleAssociationCollection;
    class TrackState;
}

class TruthTrackerAlg: public GaudiAlgorithm
{
    public:
        TruthTrackerAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize() override;
        virtual StatusCode execute() override;
        virtual StatusCode finalize() override;

    private:
        void getTrackStateFromMcParticle(const edm4hep::MCParticleCollection*
                mcParticleCol, edm4hep::TrackState& stat);
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::Detector* m_dd4hep;
        dd4hep::OverlayedField m_dd4hepField;
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
        void debugEvent();

        //reader
        DataHandle<edm4hep::MCParticleCollection> m_mcParticleCol{
            "MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_DCDigiCol{
            "DigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_DCHitAssociationCol{ "DCHitAssociationCollection",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackCollection>
            m_siSubsetTrackCol{ "SubsetTracks",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_SITSpacePointCol{
            "SITSpacePoints" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_SETSpacePointCol{
            "SETSpacePoints" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_FTDSpacePointCol{
            "FTDSpacePoints" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_VXDTrackerHits{
            "VXDTrackerHits" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_SETTrackerHits{
            "SETTrackerHits" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_SITTrackerHits{
            "SITTrackerHits" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_FTDTrackerHits{
            "FTDTrackerHits" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_VXDCollection{
            "VXDCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_SETCollection{
            "SETCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_SITCollection{
            "SITCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_FTDCollection{
            "FTDCollection" , Gaudi::DataHandle::Reader, this};
        //writer
        DataHandle<edm4hep::TrackCollection> m_DCTrackCol{
            "DCTrackCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackCollection> m_SDTTrackCol{
            "SDTTrackCollection", Gaudi::DataHandle::Writer, this};

        //readout for getting segmentation
        Gaudi::Property<std::string> m_readout_name{this, "readout",
            "DriftChamberHitsCollection"};
        Gaudi::Property<bool> m_useSET{this,"useSET",true};
        Gaudi::Property<bool> m_useTruthTrack{this,"useTruthTrack",false};
        Gaudi::Property<bool> m_useSiTruthHit{this,"useSiTruthHit",false};
        Gaudi::Property<bool> m_useSiSpacePoint{this,"useSiSpacePoint",true};
        Gaudi::Property<float> m_resPT{this,"resPT",0};//ratio
        Gaudi::Property<float> m_resPz{this,"resPz",0};//ratio
        Gaudi::Property<float> m_resMomPhi{this,"resMomPhi",0};//radian
        Gaudi::Property<float> m_resMomTheta{this,"resMomTheta",0};//radian
        Gaudi::Property<float> m_resVertexX{this,"resVertexX",0.003};//3um
        Gaudi::Property<float> m_resVertexY{this,"resVertexY",0.003};//3um
        Gaudi::Property<float> m_resVertexZ{this,"resVertexZ",0.003};//3um
        Gaudi::Property<int> m_maxDCDigiCut{this,"maxDCDigiCut",1e6};

        NTuple::Tuple*  m_tuple;
        NTuple::Item<int> m_run;
        NTuple::Item<int> m_evt;
        NTuple::Array<double> m_siMom;
        NTuple::Array<double> m_siPos;
        NTuple::Array<double> m_mcMom;
        NTuple::Array<double> m_mcPos;
        NTuple::Item<int> m_nSimTrackerHitVXD;
        NTuple::Item<int> m_nSimTrackerHitSIT;
        NTuple::Item<int> m_nSimTrackerHitSET;
        NTuple::Item<int> m_nSimTrackerHitFTD;
        NTuple::Item<int> m_nTrackerHitVXD;
        NTuple::Item<int> m_nTrackerHitSIT;
        NTuple::Item<int> m_nTrackerHitSET;
        NTuple::Item<int> m_nTrackerHitFTD;
        NTuple::Item<int> m_nTrackerHitErrVXD;
        NTuple::Item<int> m_nTrackerHitErrSIT;
        NTuple::Item<int> m_nTrackerHitErrSET;
        NTuple::Item<int> m_nTrackerHitErrFTD;
        NTuple::Item<int> m_nSpacePointSIT;
        NTuple::Item<int> m_nSpacePointSET;
        NTuple::Item<int> m_nSpacePointFTD;
        NTuple::Item<int> m_nSpacePointErrVXD;
        NTuple::Item<int> m_nSpacePointErrSIT;
        NTuple::Item<int> m_nSpacePointErrSET;
        NTuple::Item<int> m_nSpacePointErrFTD;
};

#endif
