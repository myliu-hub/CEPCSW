//////////////////////////////////////////////////////////////////////
///
/// This is an algorithm for track fitting for CEPC track with genfit.
///
/// In this file, including:
///   An algorithm for drift chamber track fitting with genfit with 5 hypothesis
///
///   Units are following DD4hepUnits
///
/// Authors:
///   Yao ZHANG(zhangyao@ihep.ac.cn)
///
/////////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_RECGENFITALGDC_H
#define RECGENFITALG_RECGENFITALGDC_H

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/NTuple.h"
#include "k4FWCore/DataHandle.h"
#include "DD4hep/Fields.h"
#include <string>
#include "TVector3.h"

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
    //class rec::CellIDPositionConverter;
    namespace DDSegmentation{
        class GridDriftChamber;
        class BitFieldCoder;
    }
}
namespace edm4hep{
    class EventHeaderCollection;
    class MCParticleCollection;
    class SimTrackerHitCollection;
    class TrackCollection;
    class Track;
    class TrackState;
    class TrackerHitCollection;
    class MCRecoTrackerAssociationCollection;
    class ReconstructedParticle;
    class ReconstructedParticleCollection;
}

/////////////////////////////////////////////////////////////////////////////

class RecGenfitAlgDC:public GaudiAlgorithm {
    public:
        RecGenfitAlgDC (const std::string& name, ISvcLocator* pSvcLocator);
        StatusCode initialize() override;
        StatusCode execute() override;
        StatusCode finalize() override;

    private:
        GenfitFitter* m_genfitFitter;//The pointer to a GenfitFitter
        const GenfitField* m_genfitField;//The pointer to a GenfitField

        void debugTrack(int pidType,const GenfitTrack* genfitTrack,
                const edm4hep::Track dcTrack);
        void debugInitTrack(const GenfitTrack* genfitTrack);
        void debugEvent(const edm4hep::TrackCollection* sdtTrackCol,
                const edm4hep::TrackCollection* sdtRecTrackCol,
                double eventStartTime);
        //void debugdEdx();
        void getCircleFromPosMom(double pos[3],double mom[3],
                double Bz,double q,double& r,double& xc,double& yc);
        //get track position and momentum from TrackState
        void getCircleFromTrackState(const edm4hep::TrackState& trackState,
                double& r, double& xc, double& yc,double& charge);
        void lsFit(bool smearTrack,bool smearHit,bool firstHit);

        DataHandle<edm4hep::EventHeaderCollection> _headerCol{
            "EventHeaderCol", Gaudi::DataHandle::Reader, this};
        //Drift chamber rec hit and trac
        DataHandle<edm4hep::TrackerHitCollection> m_digiDCHitsCol{
            "DigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackCollection> m_dcTrackCol{
            "DCTrackCollection", Gaudi::DataHandle::Reader, this};
        //Mc truth
        DataHandle<edm4hep::MCParticleCollection> m_mcParticleCol{
            "MCParticle", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_simDCHitCol{
            "DriftChamberHitsCollection" , Gaudi::DataHandle::Reader, this};

        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_DCHitAssociationCol{"DCHitAssociationCollection",
                Gaudi::DataHandle::Reader, this};
        //output hits and particles
        DataHandle<edm4hep::TrackerHitCollection> m_dcFitRecHitCol{
            "DCFitRecHitsCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::ReconstructedParticleCollection> m_dcRecParticleCol{
            "DCRecParticleCollection", Gaudi::DataHandle::Writer, this};
        DataHandle<edm4hep::TrackCollection> m_DCRecTrackCol{"DCRecTrackCollection",
            Gaudi::DataHandle::Writer, this};

        const unsigned int m_nPDG;//5:e,mu,pi,K,proton
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::OverlayedField m_dd4hepField;
        dd4hep::Detector* m_dd4hep;
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
        Gaudi::Property<std::string> m_readout_name{this,
            "readout", "DriftChamberHitsCollection"};
        Gaudi::Property<int> m_iEventSelect{this,"iEventSelect",-1};
        Gaudi::Property<int> m_debug{this,"debug",false};
        Gaudi::Property<int> m_debugGenfit{this,"debugGenfit",0};
        Gaudi::Property<int> m_debugLocal{this,"debugLocal",0};
        Gaudi::Property<bool> m_useMcParticleSeed{this,"useMcParticleSeed",false};
        Gaudi::Property<bool> m_debugLsFit{this,"debugLsFit",false};
        Gaudi::Property<bool> m_useFirstHitAsSeed{this,"useFirstHitAsSeed",false};
        Gaudi::Property<bool> m_smearTrack{this,"smearTrack",false};
        Gaudi::Property<bool> m_smearHit{this,"smearHit",true};
        Gaudi::Property<std::vector<float> > m_sigmaHitU{this,"sigmaHit",{0.11,0.003,0.003,0.003,0.003}};//mm,0:DC,...TODO
        Gaudi::Property<std::vector<float> > m_sigmaHitV{this,"sigmaHit",{10,0.003,0.003,0.003,0.003}};//mm,0:DC,...TODO
        Gaudi::Property<float> m_sigmaDrift{this,"sigmaDrift",0.11};//mm
        Gaudi::Property<float> m_nSigmaHit{this,"nSigmaHit",5};
        Gaudi::Property<int> m_sortMethod{this,"sortMethod",0};
        Gaudi::Property<bool> m_truthAmbig{this,"truthAmbig",true};
        Gaudi::Property<float> m_skipCorner{this,"skipCorner",999.};
        Gaudi::Property<float> m_skipNear{this,"skipNear",0.};
        Gaudi::Property<bool> m_skipWireMaterial{this,"skipWireMaterial",false};
        Gaudi::Property<double> m_initCovResPos{this,"initCovResPos",1};
        Gaudi::Property<double> m_initCovResMom{this,"initCovResMom",0.1};
        Gaudi::Property<bool> m_isUseCovTrack{this,"isUseCovTrack",false};
        Gaudi::Property<double> m_driftVelocity{this,"drift_velocity",40};
        Gaudi::Property<std::vector<float> > m_hitError{this,"hitError",
            {0.007,0.007,0.03}};
        //Fitter type default is DAFRef.
        //Candidates are DAF,DAFRef,KalmanFitter and KalmanFitterRefTrack.
        Gaudi::Property<std::string> m_fitterType{this,"fitterTyep","DAFRef"};
        Gaudi::Property<bool> m_correctBremsstrahlung{this,
            "correctBremsstrahlung",false};
        Gaudi::Property<bool> m_noMaterialEffects{this,
            "noMaterialEffects",false};
        Gaudi::Property<int> m_maxIteration{this,"maxIteration",20};
        Gaudi::Property<int> m_resortHits{this,"resortHits",true};
        Gaudi::Property<double> m_bStart{this,"bStart",100};
        Gaudi::Property<double> m_bFinal{this,"bFinal",0.01};
        Gaudi::Property<double> m_ndfCut{this,"ndfCut",1e9};
        Gaudi::Property<double> m_chi2Cut{this,"chi2Cut",1e9};
        //-1,chargedGeantino;0,1,2,3,4:e,mu,pi,K,proton
        Gaudi::Property<int> m_debugPid{this,"debugPid",-99};
        Gaudi::Property<int> m_debugMaterial{this,"debugMaterial",0};
        Gaudi::Property<bool> m_useTruthTrack{this,"useTruthTrack",true};
        Gaudi::Property<bool> m_useTruthHit{this,"useTruthHit",true};
        Gaudi::Property<std::string> m_genfitHistRootName{this,
            "genfitHistRootName",""};
        Gaudi::Property<bool> m_showDisplay{this,"showDisplay",false};
        Gaudi::Property<double> m_extMinDistCut{this,"extMinDistCut",1e-4};
        Gaudi::Property<double> m_blowUpFactor{this,"blowUpFactor",500};
        Gaudi::Property<double> m_resetOffDiagonals{this,"resetOffDiagonals",true};
        Gaudi::Property<double> m_blowUpMaxVal{this,"blowUpMaxVal",1e6};


        int m_fitSuccess[5];
        int m_nDCTrack;
        //bool m_useRecLRAmbig;

        genfit::EventDisplay* m_genfitDisplay;
        clock_t m_timer;
        int m_iEvent;

        /// tuples
        NTuple::Tuple*  m_tuple;
        NTuple::Item<int> m_run;
        NTuple::Item<int> m_evt;
        NTuple::Item<int> m_tkId;
        NTuple::Item<int> m_mcIndex;//number of navigated mcParicle

        NTuple::Item<double> m_seedMomP;//for single track
        NTuple::Item<double> m_seedMomPt;
        NTuple::Item<int> m_seedMomQ;
        NTuple::Array<double> m_seedMom;
        NTuple::Array<double> m_seedPos;
        NTuple::Item<double> m_seedCenterX;
        NTuple::Item<double> m_seedCenterY;
        NTuple::Item<double> m_seedR;

        NTuple::Matrix<double> m_truthPocaMc;//2 dim matched particle and 3 pos.
        NTuple::Matrix<double> m_pocaPosMc;//2 dim matched particle and 3 pos.
        NTuple::Matrix<double> m_pocaMomMc;//2 dim matched particle and 3 mom.
        NTuple::Array<double> m_pocaMomMcP;//2 dim matched particle and p
        NTuple::Array<double> m_pocaMomMcPt;//2 dim matched particle and pt
        NTuple::Array<double> m_pocaPosMdc;//pos 0:x,1:y,2:z
        NTuple::Array<double> m_pocaMomMdc;//mom. 0:px,1:py,2:pz
        NTuple::Item<int> m_pidIndex;
        NTuple::Matrix<double> m_firstPosKal;//5 hyposis and pos. at first
        NTuple::Array<double> m_firstMomKalP;//5 hyposis and mom. at first
        NTuple::Array<double> m_firstMomKalPt;//5 hyposis and mom. at first

        NTuple::Array<double> m_ErrorcovMatrix;
        NTuple::Item<double> m_D0;
        NTuple::Item<double> m_phi;
        NTuple::Item<double> m_omega;
        NTuple::Item<double> m_Z0;
        NTuple::Item<double> m_tanLambda;

        NTuple::Item<double> mcP_D0;
        NTuple::Item<double> mcP_phi;
        NTuple::Item<double> mcP_omega;
        NTuple::Item<double> mcP_Z0;
        NTuple::Item<double> mcP_tanLambda;

        NTuple::Item<double> m_fitPosx;
        NTuple::Item<double> m_fitPosy;
        NTuple::Item<double> m_fitPosz;

        NTuple::Item<double> m_fitMomx;
        NTuple::Item<double> m_fitMomy;
        NTuple::Item<double> m_fitMomz;
        NTuple::Item<double> m_fitMom;
        NTuple::Item<double> m_firstMomMc;

        NTuple::Array<double> m_extraPos;
        NTuple::Array<double> m_extraMom;
        NTuple::Array<double> m_Error6;

        NTuple::Item<int> m_nSdtTrack;
        NTuple::Item<int> m_nDcTrack;
        NTuple::Item<int> m_nDCDigi;

        NTuple::Matrix<double> m_pocaPosKal;//5 hyposis and 3 mom.
        NTuple::Matrix<double> m_pocaMomKal;//5 hyposis and 3 mom.
        NTuple::Array<double> m_pocaMomKalPFirst;//5 hyposis and p
        NTuple::Array<double> m_pocaMomKalP;//5 hyposis and p
        NTuple::Array<double> m_pocaMomKalPt;//5 hyposis and pt
        NTuple::Array<int> m_chargeKal;
        NTuple::Array<double> m_chi2Kal;
        NTuple::Array<double> m_nDofKal;
        NTuple::Array<int> m_isFitConverged;
        NTuple::Array<int> m_isFitConvergedFully;
        NTuple::Array<int> m_isFitted;
        NTuple::Item<int> m_nDigi;
        NTuple::Item<int> m_nHitMc;

        NTuple::Item<int> m_nSimDCHit;
        NTuple::Array<int> m_nHitWithFitInfo;
        NTuple::Item<int> m_nHitKalInput;
        NTuple::Array<double> m_dcDigiChamber;
        NTuple::Array<double> m_dcDigiLayer;
        NTuple::Array<double> m_dcDigiCell;
        NTuple::Array<double> m_dcDigiTime;
        NTuple::Array<double> m_dcDigiDoca;
        NTuple::Array<double> m_dcDigiDocaMC;
        NTuple::Array<double> m_dcDigiDocaExt;
        NTuple::Array<double> m_dcDigiDocaIdeal;
        NTuple::Array<double> m_dcDigiPocaExtX;
        NTuple::Array<double> m_dcDigiPocaExtY;
        NTuple::Array<double> m_dcDigiPocaExtZ;
        NTuple::Array<double> m_dcDigiWireStartX;
        NTuple::Array<double> m_dcDigiWireStartY;
        NTuple::Array<double> m_dcDigiWireStartZ;
        NTuple::Array<double> m_dcDigiWireEndX;
        NTuple::Array<double> m_dcDigiWireEndY;
        NTuple::Array<double> m_dcDigiWireEndZ;
        NTuple::Array<double> m_dcDigiMcPosX;
        NTuple::Array<double> m_dcDigiMcPosY;
        NTuple::Array<double> m_dcDigiMcPosZ;
        NTuple::Array<double> m_dcDigiMcMomX;
        NTuple::Array<double> m_dcDigiMcMomY;
        NTuple::Array<double> m_dcDigiMcMomZ;

        NTuple::Array<double> m_dcHitDriftT;
        NTuple::Array<double> m_dcHitDriftDl;
        NTuple::Array<double> m_dcHitDriftDr;
        NTuple::Array<int> m_dcHitLr;
        NTuple::Array<int> m_dcHitLayer;
        NTuple::Array<int> m_dcHitWire;
        NTuple::Array<double> m_dcHitExpDoca;
        NTuple::Array<double> m_dcHitExpMcDoca;
        NTuple::Array<double> m_dcHitErr;
        NTuple::Array<int> m_nHitFailedKal;
        NTuple::Array<int> m_nHitFitted;
        NTuple::Array<double> m_time;
        //truth
        NTuple::Array<int> m_dcHitMcLr;
        NTuple::Array<int> m_dcHitMcTkId;
        NTuple::Array<double> m_dcHitMcDrift;
        NTuple::Array<double> m_dcHitMcX;
        NTuple::Array<double> m_dcHitMcY;
        NTuple::Array<double> m_dcHitMcZ;
        NTuple::Array<double> m_dcHitMcDoca;
        NTuple::Array<double> m_dcHitMcWireX;
        NTuple::Array<double> m_dcHitMcWireY;
        NTuple::Array<double> m_dcHitExpMcPocaX;
        NTuple::Array<double> m_dcHitExpMcPocaY;
        NTuple::Array<double> m_dcHitExpMcPocaZ;
        NTuple::Array<double> m_dcHitExpMcPocaWireX;
        NTuple::Array<double> m_dcHitExpMcPocaWireY;
        NTuple::Array<double> m_dcHitExpMcPocaWireZ;

        //fit
        NTuple::Item<int> m_genfitTrackNumPoint;
        NTuple::Item<int> m_genfitTrackNumPointsWithMeas;
        NTuple::Item<int> m_genfitNHit;
        NTuple::Array<int> m_genfitHitLayer;
        NTuple::Array<int> m_genfitHitCell;
        NTuple::Array<double> m_genfitHitEndX;
        NTuple::Array<double> m_genfitHitEndY;
        NTuple::Array<double> m_genfitHitEndZ;
        NTuple::Array<double> m_genfitHitDrift;
        NTuple::Array<double> m_genfitHitDriftErr;
        NTuple::Array<double> m_genfitTrackPos;
        NTuple::Array<double> m_genfitTrackMom;
        NTuple::Item<double> m_genfitTimeSeed;
        NTuple::Item<double> m_momLsMCP;
        NTuple::Item<double> m_momLsMCPFirst;
        NTuple::Item<double> m_momLsP;
        NTuple::Item<double> m_momLsPt;
        NTuple::Item<double> m_momLsPz;

};
#endif
