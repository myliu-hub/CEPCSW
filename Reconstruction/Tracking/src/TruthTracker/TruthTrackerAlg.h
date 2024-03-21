#ifndef TruthTrackerAlg_h
#define TruthTrackerAlg_h

#include "GaudiAlg/GaudiAlgorithm.h"
#include "k4FWCore/DataHandle.h"
#include "DD4hep/Fields.h"
#include "GaudiKernel/NTuple.h"

#include "TRandom3.h"

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
    class MCParticle;
    class SimTrackerHitCollection;
    class SimTrackerHit;
    class TrackerHitCollection;
    class TrackerHit;
    class TrackCollection;
    class Track;
    class MutableTrack;
    class TrackState;
    class ReconstructedParticleCollection;
    class MCRecoTrackerAssociationCollection;
    class MCRecoParticleAssociationCollection;
}

class TruthTrackerAlg: public GaudiAlgorithm
{
    public:
        TruthTrackerAlg(const std::string& name, ISvcLocator* svcLoc);

        virtual StatusCode initialize() override;
        virtual StatusCode execute() override;
        virtual StatusCode finalize() override;

    private:
        int eventNo;
        void getTrackStateFromMcParticle(edm4hep::MCParticle mcParticleCol,
                edm4hep::TrackState& stat);
        int addSimHitsToTk(DataHandle<edm4hep::SimTrackerHitCollection>&
                colHandle, edm4hep::TrackerHitCollection*& truthTrackerHitCol,
                edm4hep::MutableTrack& track, const char* msg,int nHitAdded);
        int smearDCTkhit(DataHandle<edm4hep::TrackerHitCollection>&
                colHandle,DataHandle<edm4hep::TrackerHitCollection>& smearCol,
                DataHandle<edm4hep::SimTrackerHitCollection>& SimDCHitCol,
                DataHandle<edm4hep::SimTrackerHitCollection>& SimSmearDCHitCol,
                DataHandle<edm4hep::MCRecoTrackerAssociationCollection>& AssoDCHitCol,
                DataHandle<edm4hep::MCRecoTrackerAssociationCollection>& AssoSmearDCHitCol,
                float resX, float resY, float resZ);
        int addHitsToTk(DataHandle<edm4hep::TrackerHitCollection>&
                colHandle, edm4hep::MutableTrack& track, const char* msg,int nHitAdded);
        int addOneLayerHitsToTk(DataHandle<edm4hep::TrackerHitCollection>&
                colHandle, edm4hep::MutableTrack& track, const char* msg,int nHitAdded);
        int addIdealHitsToTk(DataHandle<edm4hep::TrackerHitCollection>&
                colHandle, edm4hep::TrackerHitCollection*& truthTrackerHitCol,
                edm4hep::MutableTrack& track, const char* msg,int nHitAdded);

        int addHotsToTk(edm4hep::Track& sourceTrack,edm4hep::MutableTrack&
                targetTrack, int hitType,const char* msg,int nHitAdded);
        int nHotsOnTrack(edm4hep::Track& track, int hitType);
        int trackerHitColSize(DataHandle<edm4hep::TrackerHitCollection>& hitCol);
        int simTrackerHitColSize(DataHandle<edm4hep::SimTrackerHitCollection>& hitCol);
        bool getTrackStateFirstHit(
                DataHandle<edm4hep::SimTrackerHitCollection>& dcSimTrackerHitCol,
                float charge,edm4hep::TrackState& trackState);
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::Detector* m_dd4hep;
        dd4hep::OverlayedField m_dd4hepField;
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
        void debugEvent();
        bool debugVertex(edm4hep::Track sourceTrack,
                const edm4hep::MCRecoTrackerAssociationCollection* vxdAsso,
                const edm4hep::MCRecoTrackerAssociationCollection* sitAsso,
                const edm4hep::MCRecoTrackerAssociationCollection* setAsso, 
                const edm4hep::MCRecoTrackerAssociationCollection* ftdAsso);
        //unit length is mm
        void getCircleFromPosMom(float pos[3],float mom[3],
                float Bz,float q,float& helixRadius,float& helixXC,float& helixYC);
        int makeNoiseHit(edm4hep::SimTrackerHitCollection* SimVec,
                edm4hep::TrackerHitCollection* Vec,
                edm4hep::MCRecoTrackerAssociationCollection* AssoVec,
                const edm4hep::TrackerHitCollection* digiDCHitsCol,
                const edm4hep::MCRecoTrackerAssociationCollection* assoHits);
        bool debugNoiseHitsCol(edm4hep::TrackerHitCollection* Vec);
        int getSiMCParticle(edm4hep::Track siTack,
                const edm4hep::TrackerHitCollection* dcTrackerHits,
                const edm4hep::MCRecoTrackerAssociationCollection* vxdAsso,
                const edm4hep::MCRecoTrackerAssociationCollection* sitAsso,
                const edm4hep::MCRecoTrackerAssociationCollection* setAsso,
                const edm4hep::MCRecoTrackerAssociationCollection* ftdAsso,
                const edm4hep::MCRecoTrackerAssociationCollection* dcAsso);
        edm4hep::SimTrackerHit getAssoSimTrackerHit(const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
                edm4hep::TrackerHit trackerHit) const;
        void getAssoMCParticle(const edm4hep::MCRecoTrackerAssociationCollection* assoHits,edm4hep::TrackerHit trackerHit,
                edm4hep::MCParticle& mcParticle) const;

            //reader
            //        DataHandle<edm4hep::TrackerHitCollection> m_NoiseHitCol{
            //            "NoiseDCHitsCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCParticleCollection> m_mcParticleCol{
            "MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_DCSimTrackerHitCol{
            "DriftChamberHitsCollection", Gaudi::DataHandle::Reader, this};
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
        //        DataHandle<edm4hep::TrackerHitCollection> m_SETSpacePointCol{
        //            "SETSpacePoints" , Gaudi::DataHandle::Reader, this};
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
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_VXDHitAssociationCol{ "VXDTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_SETCollection{
            "SETCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_SETHitAssociationCol{ "SETTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_SITCollection{
            "SITCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_SITHitAssociationCol{ "SITTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_FTDCollection{
            "FTDCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_FTDHitAssociationCol{ "FTDTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        //writer
        DataHandle<edm4hep::TrackCollection> m_DCTrackCol{
            "DCTrackCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackCollection> m_SDTTrackCol{
            "SDTTrackCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackerHitCollection> m_truthTrackerHitCol{
            "TruthTrackerHitCollection", Gaudi::DataHandle::Writer, this};

        //readout for getting segmentation
        Gaudi::Property<std::string> m_readout_name{this, "readout",
            "DriftChamberHitsCollection"};
        Gaudi::Property<bool> m_hist{this,"hist",false};
        Gaudi::Property<bool> m_useDC{this,"useDC",true};
        Gaudi::Property<bool> m_useSi{this,"useSi",true};
        Gaudi::Property<bool> m_useTruthTrack{this,"useTruthTrack",false};
        Gaudi::Property<bool> m_useSiTruthHit{this,"useSiTruthHit",false};
        Gaudi::Property<bool> m_useSiSimHit{this,"useSiSimHit",false};
        Gaudi::Property<bool> m_skipSecondaryHit{this,"skipSecondaryHit",true};
        Gaudi::Property<bool> m_useFirstHitForDC{this,"useFirstHitForDC",false};
        Gaudi::Property<bool> m_useSiSpacePoint{this,"useSiSpacePoint",false};
        Gaudi::Property<bool> m_useIdealHit{this,"useIdealHit",false};

        Gaudi::Property<bool> m_useMultiTrack{this,"useMultiTrack",false};
        Gaudi::Property<bool> m_setSiTrackConditions{this,"SetSiTrackConditions",true};

        Gaudi::Property<bool> m_useTrackFinding{this,"DCTrackFinding",false};
        Gaudi::Property<bool> m_useSelectHitsOneLayer{this,"useSelectHitsOneLayer",false};

        Gaudi::Property<float> m_momentumCut{this,"momentumCut",0.1};//momentum cut for the first hit
        Gaudi::Property<float> m_resPT{this,"resPT",0};//ratio
        Gaudi::Property<float> m_resPz{this,"resPz",0};//ratio
        Gaudi::Property<float> m_resX{this,"resX",0.11};//mm
        Gaudi::Property<float> m_resY{this,"resY",0.11};//mm
        Gaudi::Property<float> m_resZ{this,"resZ",0.11};//mm
        Gaudi::Property<float> m_driftVelocity{this,"driftVelocity",40};// um/us
        Gaudi::Property<float> m_resMomPhi{this,"resMomPhi",0};//radian
        Gaudi::Property<float> m_resMomTheta{this,"resMomTheta",0};//radian
        Gaudi::Property<float> m_resVertexX{this,"resVertexX",0.003};//3um
        Gaudi::Property<float> m_resVertexY{this,"resVertexY",0.003};//3um
        Gaudi::Property<float> m_resVertexZ{this,"resVertexZ",0.003};//3um
        Gaudi::Property<int> m_maxDCDigiCut{this,"maxDigiCut",1e6};
        Gaudi::Property<int> m_debugPID{this,"debugPID",0};// {0,1,2,3,4} = {e-,mu-,pi-,K-,p}
        Gaudi::Property<std::vector<float> > m_resVXD{this,"resVXD",{0.003,0.003,0.003}};//mm
        Gaudi::Property<std::vector<float> > m_resSIT{this,"resSIT",{0.003,0.003,0.003}};//mm
        Gaudi::Property<std::vector<float> > m_resSET{this,"resSET",{0.003,0.003,0.003}};//mm
        Gaudi::Property<std::vector<float> > m_resFTDPixel{this,"resFTDPixel",{0.003,0.003,0.003}};//mm
        Gaudi::Property<std::vector<float> > m_resFTDStrip{this,"resFTDStrip",{0.003,0.003,0.003}};//mm
        Gaudi::Property<std::vector<float> > m_cutMomHit{this,"cutMomHit",{0.65,0.04,0.2,0.05,0.14}};//e,mu,pi,K,p

        float m_helixRadius,m_helixXC,m_helixYC;
        float m_helixRadiusFirst,m_helixXCFirst,m_helixYCFirst;

        NTuple::Tuple*  m_tuple;
        NTuple::Item<int> m_nSDTTrack;
        NTuple::Item<int> m_nDCTrackHit;
        NTuple::Item<int> m_nSimTrackerHitVXD;
        NTuple::Item<int> m_nSimTrackerHitSIT;
        NTuple::Item<int> m_nSimTrackerHitSET;
        NTuple::Item<int> m_nSimTrackerHitFTD;
        NTuple::Item<int> m_nSimTrackerHitDC;
        NTuple::Item<int> m_nTrackerHitVXD;
        NTuple::Item<int> m_nTrackerHitSIT;
        NTuple::Item<int> m_nTrackerHitSET;
        NTuple::Item<int> m_nTrackerHitFTD;
        NTuple::Item<int> m_nTrackerHitDC;

        NTuple::Array<int> m_nHitOnSiTkVXD;
        NTuple::Array<int> m_nHitOnSiTkSIT;
        NTuple::Array<int> m_nHitOnSiTkSET;
        NTuple::Array<int> m_nHitOnSiTkFTD;
        NTuple::Array<int> m_nHitOnSdtTkVXD;
        NTuple::Array<int> m_nHitOnSdtTkSIT;
        NTuple::Array<int> m_nHitOnSdtTkSET;
        NTuple::Array<int> m_nHitOnSdtTkFTD;
        NTuple::Array<int> m_nHitOnSdtTkDC;
        NTuple::Array<int> m_nHitOnSdtTk;

        NTuple::Array<float> m_siXc;
        NTuple::Array<float> m_siYc;
        NTuple::Array<float> m_siR;

        TRandom3 fRandom;
};

#endif
