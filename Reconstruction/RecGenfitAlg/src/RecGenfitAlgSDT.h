//////////////////////////////////////////////////////////////////////
///
/// This is an algorithm for track fitting for CEPC track with genfit.
///
/// In this file, including:
///   An algorithm for combined silicon and drift chamber track fitting
///   with genfit for 5 particle hypothesis
///
///   Units are following DD4hepUnits
///
/// Authors:
///   Yao ZHANG(zhangyao@ihep.ac.cn)
///
/////////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_RECGENFITALGSDT_H
#define RECGENFITALG_RECGENFITALGSDT_H

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/NTuple.h"
#include "k4FWCore/DataHandle.h"
#include "DD4hep/Fields.h"
#include <string>
#include <vector>
#include "AbsKalmanFitter.h"

class GenfitFitter;
class GenfitField;
class GenfitTrack;
class IGeomSvc;
class time;
namespace genfit{
    class EventDisplay;
}
namespace dd4hep {
    class Detector;
    namespace DDSegmentation{
        class GridDriftChamber;
        class BitFieldCoder;
    }
}
namespace edm4hep{
    class EventHeaderCollection;
    class MCParticleCollection;
    class MutableReconstructedParticle;
    class SimTrackerHitCollection;
    class TrackCollection;
    class TrackerHitCollection;
    class MCRecoTrackerAssociationCollection;
    class ReconstructedParticle;
    class ReconstructedParticleCollection;
    class MutableReconstructedParticleCollection;
    class TrackerHit;
    class Track;
}

/////////////////////////////////////////////////////////////////////////////

class RecGenfitAlgSDT:public GaudiAlgorithm {
    public:
        RecGenfitAlgSDT (const std::string& name, ISvcLocator* pSvcLocator);
        StatusCode initialize() override; StatusCode execute() override; StatusCode finalize() override;

    private:
        GenfitFitter* m_genfitFitter;//The pointer to a GenfitFitter
        const GenfitField* m_genfitField;//The pointer to a GenfitField

        void debugTrack(int iStrack,int pidType,const GenfitTrack* genfitTrack,
                TVector3 pocaToOrigin_pos,TVector3 pocaToOrigin_mom,
                TMatrixDSym pocaToOrigin_cov);
        void debugEvent(const edm4hep::TrackCollection* sdtTrackCol,
                const edm4hep::TrackCollection* sdtRecTrackCol,
                double eventStartTime,int nFittedSDT);

        void selectHits(const edm4hep::Track&, std::vector<edm4hep::TrackerHit*>& dcDigiSelected);
        bool addHitsToTk(edm4hep::TrackerHit hit,edm4hep::Track& track,const char* msg) const;
        void addSETHitsToTk(edm4hep::Track& track,
                const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
                const edm4hep::TrackerHitCollection* SETHits);

        DataHandle<edm4hep::EventHeaderCollection> m_headerCol{
            "EventHeaderCol", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_DCDigiCol{
            "DigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_TruthDCDigiCol{
            "TruthDigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        //Mc truth
        DataHandle<edm4hep::SimTrackerHitCollection> m_simVXDHitCol{
            "VXDCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_simSETHitCol{
            "SETCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_simSITHitCol{
            "SITCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_simFTDHitCol{
            "FTDCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCParticleCollection> m_mcParticleCol{
            "MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_simDCHitCol{
            "DriftChamberHitsCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_DCHitAssociationCol{"DCHitAssociationCollection",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackCollection> m_dcTrackCol{
            "DCTrackCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_VXDHitAssociationCol{"VXDTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_SITHitAssociationCol{"SITTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_SETHitAssociationCol{"SETTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_FTDHitAssociationCol{"FTDTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};

        DataHandle<edm4hep::TrackerHitCollection> m_SETHitCol{
            "SETTrackerHits",Gaudi::DataHandle::Reader, this};

        //Track from silicon detectors
        DataHandle<edm4hep::TrackCollection> m_SDTTrackCol{"SDTTrackCollection",
            Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackCollection> m_SDTTrackFindCol{"SDTTrackFindCollection",
            Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackCollection> m_SDTRecTrackCol{"SDTRecTrackCollection",
            Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackCollection> m_SDTDebugRecTrackCol{"SDTDebugRecTrackCollection",
            Gaudi::DataHandle::Writer, this};

        //Track Finding
        DataHandle<edm4hep::TrackerHitCollection> m_DCTrackFindingCol{
            "DCTrackFindingHitCollection",Gaudi::DataHandle::Reader, this};

        //Output hits and particles
        DataHandle<edm4hep::ReconstructedParticleCollection> m_SDTRecParticleCol{
            "SDTRecParticleCollection", Gaudi::DataHandle::Writer, this};

        const unsigned int m_nPDG;//5:e,mu,pi,K,proton
        int m_eventNo;
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::OverlayedField m_dd4hepField;
        dd4hep::Detector* m_dd4hepDetector;
        double m_cell_width;
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
        Gaudi::Property<std::string> m_readout_name{this,
            "readout", "DriftChamberHitsCollection"};
        Gaudi::Property<int> m_iterFit{this,"iterFit",3};
        Gaudi::Property<int> m_debug{this,"debug",0};
        Gaudi::Property<int> m_debugGenfit{this,"debugGenfit",0};
        Gaudi::Property<int> m_debugPid{this,"debugPid",-99};
        Gaudi::Property<int> m_eventNoSelection{this,"eventNoSelection",1e9};
        Gaudi::Property<std::vector<float> > m_sigmaHitU{this,
            "sigmaHitU",{0.11, // DC z mm
                0.0028,0.006,0.004,0.004,0.004,0.004, //VXD U
                0.0072, //SIT U 4 layers same resolusion
                0.0072, //SET U
                0.003,0.003,0.0072,0.0072,0.0072,0.0072,0.0072}};//FTD V
        //mm, 0:DC, 1~7:VXD, 8:SIT, 9:SET, FTD:10~16
        Gaudi::Property<std::vector<float> > m_sigmaHitV{this,
            "sigmaHitV",{0.11, // DC z mm
                0.0028,0.006,0.004,0.004,0.004,0.004, //VXD V
                0.086, //SIT V
                0.086, //SET V
                0.003,0.003,0.0072,0.0072,0.0072,0.0072,0.0072}};//FTD V
        Gaudi::Property<float> m_cutSETPos{this,"cutSETPos",0.43}; //mm :2024.01.25之前为50mm
        Gaudi::Property<int> m_measurementTypeSi{this,"measurementTypeSi",0};
        //-1: not use, 0, space point, 1, pixel or planer measurement
        Gaudi::Property<int> m_measurementTypeDC{this,"measurementTypeDC",0};
        //-1: not use, 0, space point, 1, wire measurement
        //Gaudi::Property<bool> m_smearHit{this,"smearHit",true};
        //Gaudi::Property<float> m_nSigmaHit{this,"nSigmaHit",5};
        Gaudi::Property<int> m_sortMethod{this,"sortMethod",0};
        Gaudi::Property<bool> m_truthAmbig{this,"truthAmbig",true};
        Gaudi::Property<bool> m_addSET{this,"addSET",true};
        Gaudi::Property<bool> m_useSalHits{this,"useSalHits",true};
        Gaudi::Property<bool> m_checkGenfitTrack{this,"checkGenfitTrack",false};
        Gaudi::Property<float> m_skipCorner{this,"skipCorner",1200.};
        Gaudi::Property<float> m_skipNear{this,"skipNear",0.};
        Gaudi::Property<bool> m_isUseCovTrack{this,"isUseCovTrack",false};
        //Fitter type default is DAFRef.
        //Candidates are DAF,DAFRef,KalmanFitter and KalmanFitterRefTrack.
        Gaudi::Property<std::string> m_fitterType{this,"fitterType","DAF"};
        Gaudi::Property<bool> m_correctBremsstrahlung{this,
            "correctBremsstrahlung",false};
        Gaudi::Property<bool> m_noMaterialEffects{this,
            "noMaterialEffects",false};
        Gaudi::Property<bool> m_skipWireMaterial{this,
            "skipWireMaterial",true};
        Gaudi::Property<int> m_maxIteration{this,"maxIteration",20};
        Gaudi::Property<int> m_minIteration{this,"minIteration",4};
        Gaudi::Property<int> m_resortHits{this,"resortHits",true};
        Gaudi::Property<int> m_chi2FitCut{this,"chi2FitCut",30};
        Gaudi::Property<int> m_dcFitNum{this,"dcFitNum",5};
        Gaudi::Property<double> m_bStart{this,"bStart",100};
        Gaudi::Property<double> m_bFinal{this,"bFinal",0.01};
        Gaudi::Property<double> m_DCCornerCuts{this,"dcCornerCuts",-999};
        Gaudi::Property<double> m_ndfCut{this,"ndfCut",1e9};
        Gaudi::Property<double> m_chi2Cut{this,"chi2Cut",1e9};
        //-1,chargedGeantino;0,1,2,3,4:e,mu,pi,K,proton
        Gaudi::Property<bool> m_useTruthTrack{this,"useTruthTrack",false};
        Gaudi::Property<std::string> m_genfitHistRootName{this,
            "genfitHistRootName",""};
        Gaudi::Property<bool> m_showDisplay{this,"showDisplay",false};
        Gaudi::Property<bool> m_fitSiliconOnly{this,"fitSiliconOnly",false};
        Gaudi::Property<bool> m_isUseFixedSiHitError{this,"isUseFixedSiHitError",false};
        Gaudi::Property<std::vector<float> > m_hitError{this,"hitError",
            {0.007,0.007,0.03}};
        Gaudi::Property<double> m_extMinDistCut{this,"extMinDistCut",1e-4};
        Gaudi::Property<int> m_multipleMeasurementHandling{this,
            "multipleMeasurementHandling",
            (int) genfit::eMultipleMeasurementHandling::unweightedClosestToPredictionWire};
        Gaudi::Property<double> m_driftVelocity{this,"drift_velocity",40};//um/ns
        Gaudi::Property<bool> m_selectDCHit{this,"selectDCHit",false};
        Gaudi::Property<bool> m_useNoiseDCHit{this,"useNoiseDCHit",false};
        Gaudi::Property<double> m_docaCut{this,"docaCut",3.3};//mm

        Gaudi::Property<bool> m_useTrackFinding{this,"DCTrackFinding",false};
        Gaudi::Property<double> m_extraCut{this,"extraCut",0.055}; //5*sigma cm
        Gaudi::Property<double> m_extraDocaCut{this,"extraDocaCut",1.2}; //cm
        int m_fitSuccess[5];
        int m_nRecTrack;
        bool m_firstTuple;

        genfit::EventDisplay* m_genfitDisplay;

        /// tuples
        NTuple::Tuple*  m_tuple;
        NTuple::Item<int> m_run;
        NTuple::Item<int> m_evt;
        NTuple::Item<int> m_tkId;
        NTuple::Item<int> m_mcIndex;//number of navigated mcParicle
        NTuple::Matrix<float> m_truthPocaMc;//2 dim matched particle and 3 pos.
        NTuple::Array<float> m_seedMomP;//for some track
        NTuple::Array<float> m_seedMomPt;
        NTuple::Array<int> m_seedMomQ;
        NTuple::Matrix<float> m_seedMom;
        NTuple::Matrix<float> m_seedPos;
        NTuple::Matrix<float> m_pocaPosMc;//2 dim matched particle and 3 pos.
        NTuple::Matrix<float> m_pocaMomMc;//2 dim matched particle and 3 mom.
        NTuple::Array<float> m_pocaMomMcP;//2 dim matched particle and p
        NTuple::Array<float> m_pocaMomMcPt;//2 dim matched particle and pt
        NTuple::Matrix<float> m_pocaPosMdc;//pos 0:x,1:y,2:z
        NTuple::Matrix<float> m_pocaMomMdc;//mom. 0:px,1:py,2:pz
        NTuple::Item<int> m_pidIndex;
        NTuple::Matrix<float> m_firstPosKal;//5 hyposis and pos. at first
        NTuple::Array<float> m_firstMomKalP;//5 hyposis and mom. at first
        NTuple::Array<float> m_firstMomKalPt;//5 hyposis and mom. at first

        NTuple::Matrix<float> m_ErrorcovMatrix6;
        NTuple::Array<float> m_posx;
        NTuple::Array<float> m_posy;
        NTuple::Array<float> m_posz;

        NTuple::Array<float> m_momx;
        NTuple::Array<float> m_momy;
        NTuple::Array<float> m_momz;

        NTuple::Array<float> m_fittedXc;
        NTuple::Array<float> m_fittedYc;
        NTuple::Array<float> m_fittedR;

        NTuple::Array<float> m_PosMcX;
        NTuple::Array<float> m_PosMcY;
        NTuple::Array<float> m_PosMcZ;

        NTuple::Array<float> m_MomMcX;
        NTuple::Array<float> m_MomMcY;
        NTuple::Array<float> m_MomMcZ;


        NTuple::Array<float> m_PocaPosX;
        NTuple::Array<float> m_PocaPosY;
        NTuple::Array<float> m_PocaPosZ;

        NTuple::Array<float> m_PocaMomX;
        NTuple::Array<float> m_PocaMomY;
        NTuple::Array<float> m_PocaMomZ;

        NTuple::Matrix<float> m_McErrCov;
        NTuple::Matrix<float> m_PocaErrCov;

        NTuple::Matrix<float> m_ErrorcovMatrix;
        NTuple::Array<float> m_D0;
        NTuple::Array<float> m_phi;
        NTuple::Array<float> m_omega;
        NTuple::Array<float> m_Z0;
        NTuple::Array<float> m_tanLambda;

        NTuple::Array<float> m_trackXc;
        NTuple::Array<float> m_trackYc;
        NTuple::Array<float> m_trackR;

        NTuple::Matrix<float> m_ErrorcovMatrix_Origin;
        NTuple::Array<float> m_D0_Origin;
        NTuple::Array<float> m_phi_Origin;
        NTuple::Array<float> m_omega_Origin;
        NTuple::Array<float> m_Z0_Origin;
        NTuple::Array<float> m_tanLambda_Origin;

        NTuple::Array<float> mcP_D0;
        NTuple::Array<float> mcP_phi;
        NTuple::Array<float> mcP_omega;
        NTuple::Array<float> mcP_Z0;
        NTuple::Array<float> mcP_tanLambda;
        NTuple::Array<float> mcP_Xc;
        NTuple::Array<float> mcP_Yc;
        NTuple::Array<float> mcP_R;

        NTuple::Matrix<float> m_pocaPosKal;//5 hyposis and 3 mom.
        NTuple::Matrix<float> m_pocaMomKal;//5 hyposis and 3 mom.
        NTuple::Matrix<float> m_pocaMomKalP;//5 hyposis and p
        NTuple::Matrix<float> m_pocaMomKalPt;//5 hyposis and pt
        NTuple::Matrix<int> m_chargeKal;
        NTuple::Matrix<float> m_chi2Kal;
        NTuple::Matrix<float> m_nDofKal;
        NTuple::Matrix<int> m_isFitConverged;
        NTuple::Matrix<int> m_isFitConvergedFully;
        NTuple::Matrix<int> m_isFitted;
        NTuple::Array<int> m_PDG;//5 hyposis and mom. at first
        NTuple::Array<int> m_numDCOnTrack;//5 hyposis and mom. at first
        NTuple::Array<int> m_salDCHits;//5 hyposis and mom. at first
        NTuple::Array<int> m_DCHitsCol;//5 hyposis and mom. at first
        NTuple::Matrix<int> m_fittedState;
        NTuple::Item<int> m_nDCDigi;
        NTuple::Item<int> m_nTruthDCDigi;
        NTuple::Item<int> m_nHitMc;
        NTuple::Item<int> m_nSdtTrack;
        NTuple::Item<int> m_nSdtTrackHit;

        NTuple::Item<int> m_nSdtRecTrack;

        NTuple::Matrix<int> m_nHitWithFitInfo;
        NTuple::Item<int> m_nHitKalInput;
        NTuple::Array<int> m_nHitDetType;
        NTuple::Matrix<int> m_nHitFailedKal;
        NTuple::Matrix<int> m_nHitFitted;

        NTuple::Array<float> m_dcDigiChamber;

        NTuple::Array<float> m_docaTrack;
        NTuple::Array<int> m_isNoise;

        NTuple::Array<float> m_dcDigiLayer;
        NTuple::Array<float> m_dcDigiCell;
        NTuple::Array<float> m_dcDigiTime;
        NTuple::Array<float> m_dcDigiDrift;
        NTuple::Array<float> m_dcDigiDocaMC;
        NTuple::Array<float> m_dcDigiPocaOnWireMCX;
        NTuple::Array<float> m_dcDigiPocaOnWireMCY;
        NTuple::Array<float> m_dcDigiWireStartX;
        NTuple::Array<float> m_dcDigiWireStartY;
        NTuple::Array<float> m_dcDigiWireStartZ;
        NTuple::Array<float> m_dcDigiWireEndX;
        NTuple::Array<float> m_dcDigiWireEndY;
        NTuple::Array<float> m_dcDigiWireEndZ;
        NTuple::Array<float> m_dcDigiMcPosX;
        NTuple::Array<float> m_dcDigiMcPosY;
        NTuple::Array<float> m_dcDigiMcPosZ;
        NTuple::Array<float> m_dcDigiMcMomX;
        NTuple::Array<float> m_dcDigiMcMomY;
        NTuple::Array<float> m_dcDigiMcMomZ;
        NTuple::Array<float> m_dcDigiDocaExt;
        NTuple::Array<float> m_dcDigiPocaExtX;
        NTuple::Array<float> m_dcDigiPocaExtY;
        NTuple::Array<float> m_dcDigiPocaExtZ;
        NTuple::Item<float> m_firstMomMc;

        NTuple::Array<int> m_nTrackerHitSDT;
        NTuple::Array<int> m_nTrackerHitDC;
        NTuple::Array<int> m_nTrackerFitHitDC;

};
#endif
