#include "RecGenfitAlgSDT.h"
#include "GenfitTrack.h"
#include "GenfitFitter.h"
#include "GenfitField.h"
#include "GenfitUnit.h"

//genfit
#include "EventDisplay.h"
#include "ReferenceStateOnPlane.h"

//cepcsw
#include "DetInterface/IGeomSvc.h"
#include "DataHelper/HelixClass.h"
#include "DataHelper/TrackHelper.h"
#include "DataHelper/TrackerHitHelper.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "UTIL/ILDConf.h"

//externals
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "UTIL/BitField64.h"
#include "DDSegmentation/Segmentation.h"
#include "TRandom.h"
#include "TLorentzVector.h"

//stl
#include <chrono>
#include "time.h"
#include <stdlib.h>
#include <thread>
#include <iostream>

#include <TTimeStamp.h>

#include <ctime>
#include <cstdlib>

#include "map"

DECLARE_COMPONENT( RecGenfitAlgSDT )

    /////////////////////////////////////////////////////////////////////
    RecGenfitAlgSDT::RecGenfitAlgSDT(const std::string& name,
            ISvcLocator* pSvcLocator):
        GaudiAlgorithm(name, pSvcLocator),m_nPDG(5),m_dd4hepDetector(nullptr),
        m_gridDriftChamber(nullptr),m_decoder(nullptr)
{
    declareProperty("EventHeaderCollection", m_headerCol);
    declareProperty("MCParticleCollection", m_mcParticleCol,
            "Handle of the input MCParticle collection");
    declareProperty("DigiDCHitCollection", m_DCDigiCol,
            "Handle of DC digi(TrakerHit) collection");
    declareProperty("TruthDigiDCHitCollection", m_TruthDCDigiCol,
            "Handle of Truth DC digi(TrakerHit) collection");
    declareProperty("DCHitAssociationCollection", m_DCHitAssociationCol,
            "Handle of DCsimTrackerHit and DCTrackerHit association collection");
    declareProperty("SDTTrackCollection", m_SDTTrackCol,
            "Handle of input silicon track collection");
    declareProperty("SDTTrackFindCollection", m_SDTTrackFindCol,
            "Handle of input silicon track finding collection");
    declareProperty("SDTRecTrackCollection",m_SDTRecTrackCol,
            "Handle of input silicon rec. track collection");
    declareProperty("DCTrackCollection", m_dcTrackCol,
            "Handle of DC track collection");
    declareProperty("SDTRecParticleCollection", m_SDTRecParticleCol,
            "Handle of silicon+drift chamber rec. particle collection");

    declareProperty("SimTrackerHitCollection",m_simVXDHitCol,
            "Handle of the VXDsimTrackerHit collection");
    declareProperty("SimTrackerHitCollection",m_simSETHitCol,
            "Handle of the SETsimTrackerHit collection");
    declareProperty("TrackerHitCollection",m_SETHitCol,
            "Handle of the SETTrackerHit collection");
    declareProperty("MCRecoTrackerAssociationCollection",m_SETHitAssociationCol,
            "Handle of the SETTrackerHitAssociation collection");
    declareProperty("SimTrackerHitCollection",m_simSITHitCol,
            "Handle of the SITsimTrackerHit collection");
    declareProperty("SimTrackerHitCollection",m_simFTDHitCol,
            "Handle of the FTDsimTrackerHit collection");
    declareProperty("SimTrackerHitCollection",m_simDCHitCol,
            "Handle of the DCsimTrackerHit collection");
    declareProperty("SimTrackerHitCollection",m_simVXDHitCol,
             "Handle of the VXDsimTrackerHit collection");

    declareProperty("DCTrackFindingHitCollection", m_DCTrackFindingCol,
            "Handle of output TrackFinding TrackerHit collection");

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgSDT::initialize()
{
    MsgStream log(msgSvc(), name());
    info()<<" RecGenfitAlgSDT initialize()"<<endmsg;

    m_eventNo=0;

    ///Get GeomSvc
    m_geomSvc=Gaudi::svcLocator()->service("GeomSvc");
    if (nullptr==m_geomSvc) {
        std::cout<<"Failed to find GeomSvc"<<std::endl;
        return StatusCode::FAILURE;
    }
    ///Get Detector
    m_dd4hepDetector=m_geomSvc->lcdd();

    ///Get Field
    m_dd4hepField=m_geomSvc->lcdd()->field();

    /// New a genfit fitter
    m_genfitFitter=new GenfitFitter(m_fitterType.toString().c_str());
    m_genfitField=new GenfitField(m_dd4hepField);
    m_genfitFitter->setField(m_genfitField);
    m_genfitFitter->setGeoMaterial(m_geomSvc->lcdd(),m_extMinDistCut,
            m_skipWireMaterial);
    m_genfitFitter->setEnergyLossBrems(m_correctBremsstrahlung);
    m_genfitFitter->setNoiseBrems(m_correctBremsstrahlung);
    if(m_debug>10) m_genfitFitter->setDebug(m_debug-10);
    if(m_noMaterialEffects) m_genfitFitter->setNoEffects(true);
    if(-1==m_debugPid) m_genfitFitter->setNoEffects(true);
    if(-1==m_debugPid) m_debugPid=0;//charged geantino with electron pid
    if(m_fitterType=="DAF"||m_fitterType=="DafRef"){
        m_genfitFitter->setMaxIterationsBetas(m_bStart,m_bFinal,m_maxIteration);
    } else {
        m_genfitFitter->setMaxIterations(m_maxIteration);
    }
    //print genfit parameters
    if(m_debug) m_genfitFitter->print();
    if(""!=m_genfitHistRootName) m_genfitFitter->initHist(m_genfitHistRootName);

    //initialize member vairables
    for(int i=0;i<5;i++) m_fitSuccess[i]=0;
    m_nRecTrack=0;
    ///Get Readout
    dd4hep::Readout readout=m_dd4hepDetector->readout(m_readout_name);
    ///Get Segmentation
    m_gridDriftChamber=dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>
        (readout.segmentation().segmentation());
    if(nullptr==m_gridDriftChamber){
        error() << "Failed to get the GridDriftChamber" << endmsg;
        return StatusCode::FAILURE;
    }

    m_cell_width = 0.5*(m_gridDriftChamber->cell_Size()); //dd4hep::cm
    debug() << " m_cell_width = " << m_cell_width << " dd4hep::cm"<< endmsg;
    m_skipCorner = m_cell_width/dd4hep::cm;

    ///Get Decoder
    m_decoder = m_geomSvc->getDecoder(m_readout_name);
    if (nullptr==m_decoder) {
        error() << "Failed to get the decoder" << endmsg;
        return StatusCode::FAILURE;
    }


    ///book tuple
    NTuplePtr nt(ntupleSvc(), "RecGenfitAlgSDT/recGenfitAlgSDT");
    if(nt){
        m_tuple=nt;
    }else{
        m_tuple=ntupleSvc()->book("RecGenfitAlgSDT/recGenfitAlgSDT",
                CLID_ColumnWiseTuple,"RecGenfitAlgSDT");
        if(m_tuple){
            StatusCode sc;
            sc=m_tuple->addItem("run",m_run);
            sc=m_tuple->addItem("evt",m_evt);
            sc=m_tuple->addItem("tkId",m_tkId);
            sc=m_tuple->addItem("nStdTrack",m_nSdtTrack);
            sc=m_tuple->addItem("nStdTrackHit",m_nSdtTrackHit,0,1000);

            sc=m_tuple->addItem("nSdtRecTrack",m_nSdtRecTrack);


            sc=m_tuple->addItem("mcIndex",m_mcIndex,0,100);//max. 100 particles
            sc=m_tuple->addItem("seedMomP",m_mcIndex,m_seedMomP);//for some track debug
            sc=m_tuple->addItem("seedMomPt",m_mcIndex,m_seedMomPt);
            sc=m_tuple->addItem("seedMomQ",m_mcIndex,m_seedMomQ);
            sc=m_tuple->addItem("seedPos",m_mcIndex,m_seedPos,3);
            sc=m_tuple->addItem("seedMom",m_mcIndex,m_seedMom,3);
            sc=m_tuple->addItem("truthPocaMc",m_mcIndex,m_truthPocaMc,3);
            sc=m_tuple->addItem("pocaPosMc",m_mcIndex,m_pocaPosMc,3);
            sc=m_tuple->addItem("pocaMomMc",m_mcIndex,m_pocaMomMc,3);
            sc=m_tuple->addItem("pocaMomMcP",m_mcIndex,m_pocaMomMcP);
            sc=m_tuple->addItem("pocaMomMcPt",m_mcIndex,m_pocaMomMcPt);
            sc=m_tuple->addItem("pocaPosMdc",m_mcIndex,m_pocaPosMdc,3);
            sc=m_tuple->addItem("pocaMomMdc",m_mcIndex,m_pocaMomMdc,3);
            sc=m_tuple->addItem("index",m_pidIndex, 0, 5);

            sc=m_tuple->addItem("ErrorcovMatrix6",m_mcIndex,m_ErrorcovMatrix6,6);
            sc=m_tuple->addItem("McErrCov",m_mcIndex,m_McErrCov,6);
            sc=m_tuple->addItem("posx",m_mcIndex,m_posx);
            sc=m_tuple->addItem("posy",m_mcIndex,m_posy);
            sc=m_tuple->addItem("posz",m_mcIndex,m_posz);

            sc=m_tuple->addItem("momx",m_mcIndex,m_momx);
            sc=m_tuple->addItem("momy",m_mcIndex,m_momy);
            sc=m_tuple->addItem("momz",m_mcIndex,m_momz);

            sc=m_tuple->addItem("fittedXc",m_mcIndex,m_fittedXc);
            sc=m_tuple->addItem("fittedYc",m_mcIndex,m_fittedYc);
            sc=m_tuple->addItem("fittedR",m_mcIndex,m_fittedR);

            sc=m_tuple->addItem("PosMcX",m_mcIndex,m_PosMcX);
            sc=m_tuple->addItem("PosMcY",m_mcIndex,m_PosMcY);
            sc=m_tuple->addItem("PosMcZ",m_mcIndex,m_PosMcZ);

            sc=m_tuple->addItem("MomMcX",m_mcIndex,m_MomMcX);
            sc=m_tuple->addItem("MomMcY",m_mcIndex,m_MomMcY);
            sc=m_tuple->addItem("MomMcZ",m_mcIndex,m_MomMcZ);

            sc=m_tuple->addItem("PocaPosX",m_mcIndex,m_PocaPosX);
            sc=m_tuple->addItem("PocaPosY",m_mcIndex,m_PocaPosY);
            sc=m_tuple->addItem("PocaPosZ",m_mcIndex,m_PocaPosZ);

            sc=m_tuple->addItem("PocaMomX",m_mcIndex,m_PocaMomX);
            sc=m_tuple->addItem("PocaMomY",m_mcIndex,m_PocaMomY);
            sc=m_tuple->addItem("PocaMomZ",m_mcIndex,m_PocaMomZ);

            sc=m_tuple->addItem("PocaErrCov",m_mcIndex,m_PocaErrCov,6);

            sc=m_tuple->addItem("ErrorcovMatrix",m_mcIndex,m_ErrorcovMatrix,15);
            sc=m_tuple->addItem("D0",m_mcIndex,m_D0);
            sc=m_tuple->addItem("phi",m_mcIndex,m_phi);
            sc=m_tuple->addItem("omega",m_mcIndex,m_omega);
            sc=m_tuple->addItem("Z0",m_mcIndex,m_Z0);
            sc=m_tuple->addItem("tanLambda",m_mcIndex,m_tanLambda);

            sc=m_tuple->addItem("trackXc",m_mcIndex,m_trackXc);
            sc=m_tuple->addItem("trackYc",m_mcIndex,m_trackYc);
            sc=m_tuple->addItem("trackR",m_mcIndex,m_trackR);

            sc=m_tuple->addItem("ErrorcovMatrix_Origin",m_mcIndex,m_ErrorcovMatrix_Origin,15);
            sc=m_tuple->addItem("D0_Origin",m_mcIndex,m_D0_Origin);
            sc=m_tuple->addItem("phi_Origin",m_mcIndex,m_phi_Origin);
            sc=m_tuple->addItem("omega_Origin",m_mcIndex,m_omega_Origin);
            sc=m_tuple->addItem("Z0_Origin",m_mcIndex,m_Z0_Origin);
            sc=m_tuple->addItem("tanLambda_Origin",m_mcIndex,m_tanLambda_Origin);

            sc=m_tuple->addItem("mcP_D0",m_mcIndex,mcP_D0);
            sc=m_tuple->addItem("mcP_phi",m_mcIndex,mcP_phi);
            sc=m_tuple->addItem("mcP_omega",m_mcIndex,mcP_omega);
            sc=m_tuple->addItem("mcP_Z0",m_mcIndex,mcP_Z0);
            sc=m_tuple->addItem("mcP_tanLambda",m_mcIndex,mcP_tanLambda);

            sc=m_tuple->addItem("mcP_Xc",m_mcIndex,mcP_Xc);
            sc=m_tuple->addItem("mcP_Yc",m_mcIndex,mcP_Yc);
            sc=m_tuple->addItem("mcP_R",m_mcIndex,mcP_R);

            sc=m_tuple->addItem("pocaPosKal",5,3,m_pocaPosKal);
            sc=m_tuple->addItem("pocaMomKal",5,3,m_pocaMomKal);
            sc=m_tuple->addItem("pocaMomKalP",m_mcIndex,m_pocaMomKalP,5);
            sc=m_tuple->addItem("pocaMomKalPt",m_mcIndex,m_pocaMomKalPt,5);
            sc=m_tuple->addItem("chargeKal",m_mcIndex,m_chargeKal,5);
            sc=m_tuple->addItem("nDofKal",m_mcIndex,m_nDofKal,5);
            sc=m_tuple->addItem("chi2Kal",m_mcIndex,m_chi2Kal,5);
            sc=m_tuple->addItem("isFitted",m_mcIndex,m_isFitted,5);
            sc=m_tuple->addItem("PDG",m_mcIndex,m_PDG);
            sc=m_tuple->addItem("numDCOnTrack",m_mcIndex,m_numDCOnTrack);
            sc=m_tuple->addItem("salDCHits",m_mcIndex,m_salDCHits);
            sc=m_tuple->addItem("DCHitsCol",m_mcIndex,m_DCHitsCol);
            sc=m_tuple->addItem("isFitConverged",m_mcIndex,m_isFitConverged,5);
            sc=m_tuple->addItem("isFitConvergedFully",m_mcIndex,
                    m_isFitConvergedFully,5);
            sc=m_tuple->addItem("fittedState",m_mcIndex,m_fittedState,5);
            sc=m_tuple->addItem("nHitFailedKal",m_mcIndex,m_nHitFailedKal,5);
            sc=m_tuple->addItem("nHitFitted",m_mcIndex,m_nHitFitted,5);
            sc=m_tuple->addItem("nDCDigi",m_nDCDigi,0,50000);
            sc=m_tuple->addItem("nTruthDCDigi",m_nTruthDCDigi,0,50000);
            sc=m_tuple->addItem("nHitKalInput",m_nHitKalInput,0,300000);
            //10 is greater than # of tracking detectors
            sc=m_tuple->addItem("hitDetID",10,m_nHitDetType);
            sc=m_tuple->addItem("nHitWithFitInfo",m_mcIndex,m_nHitWithFitInfo,5);

            sc=m_tuple->addItem("dcDigiChamber",m_nDCDigi,m_dcDigiChamber);
            sc=m_tuple->addItem("dcDigiLayer",m_nDCDigi,m_dcDigiLayer);
            sc=m_tuple->addItem("dcDigiCell",m_nDCDigi,m_dcDigiCell);
            sc=m_tuple->addItem("dcDigiTime",m_nDCDigi,m_dcDigiTime);
            sc=m_tuple->addItem("dcDigiDrift",m_nDCDigi,m_dcDigiDrift);
            sc=m_tuple->addItem("dcDigiDocaMC",m_nDCDigi,m_dcDigiDocaMC);
            sc=m_tuple->addItem("dcDigiPocaOnWireMCX",m_nDCDigi,m_dcDigiPocaOnWireMCX);
            sc=m_tuple->addItem("dcDigiPocaOnWireMCY",m_nDCDigi,m_dcDigiPocaOnWireMCY);
            sc=m_tuple->addItem("dcDigiWireStartX",m_nDCDigi,m_dcDigiWireStartX);
            sc=m_tuple->addItem("dcDigiWireStartY",m_nDCDigi,m_dcDigiWireStartY);
            sc=m_tuple->addItem("dcDigiWireStartZ",m_nDCDigi,m_dcDigiWireStartZ);
            sc=m_tuple->addItem("dcDigiWireEndX",m_nDCDigi,m_dcDigiWireEndX);
            sc=m_tuple->addItem("dcDigiWireEndY",m_nDCDigi,m_dcDigiWireEndY);
            sc=m_tuple->addItem("dcDigiWireEndZ",m_nDCDigi,m_dcDigiWireEndZ);
            sc=m_tuple->addItem("dcDigiMcMomX",m_nDCDigi,m_dcDigiMcMomX);
            sc=m_tuple->addItem("dcDigiMcMomY",m_nDCDigi,m_dcDigiMcMomY);
            sc=m_tuple->addItem("dcDigiMcMomZ",m_nDCDigi,m_dcDigiMcMomZ);
            sc=m_tuple->addItem("dcDigiMcPosX",m_nDCDigi,m_dcDigiMcPosX);
            sc=m_tuple->addItem("dcDigiMcPosY",m_nDCDigi,m_dcDigiMcPosY);
            sc=m_tuple->addItem("dcDigiMcPosZ",m_nDCDigi,m_dcDigiMcPosZ);
            sc=m_tuple->addItem("docaTrack",m_nDCDigi,m_docaTrack);
            sc=m_tuple->addItem("isNoise",m_nDCDigi,m_isNoise);
            sc=m_tuple->addItem("firstMomMc",m_firstMomMc);

            sc=m_tuple->addItem("dcDigiDocaExt",m_nDCDigi,m_dcDigiDocaExt);
            sc=m_tuple->addItem("dcDigiPocaExtX",m_nDCDigi,m_dcDigiPocaExtX);
            sc=m_tuple->addItem("dcDigiPocaExtY",m_nDCDigi,m_dcDigiPocaExtY);
            sc=m_tuple->addItem("dcDigiPocaExtZ",m_nDCDigi,m_dcDigiPocaExtZ);

            sc=m_tuple->addItem("nTrackerHitDC",100,m_nTrackerHitDC);
            sc=m_tuple->addItem("nTrackerFitHitDC",100,m_nTrackerFitHitDC);
            sc=m_tuple->addItem("nTrackerHitSDT",100,m_nTrackerHitSDT);
            debug()<< "Book tuple RecGenfitAlgSDT/recGenfitAlgSDT" << endmsg;
        }else{
            warning()<<"Tuple RecGenfitAlgSDT/recGenfitAlgSDT not booked"<<endmsg;
        }
    }//end of book tuple

    //init genfit event display
    if(m_showDisplay) m_genfitDisplay = genfit::EventDisplay::getInstance();

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgSDT::execute()
{
    info()<<"RecGenfitAlgSDT in execute()"<<endmsg;

    edm4hep::ReconstructedParticleCollection* sdtRecParticleCol=
        m_SDTRecParticleCol.createAndPut();

    edm4hep::TrackCollection* sdtRecTrackCol=
        m_SDTRecTrackCol.createAndPut();

    StatusCode sc=StatusCode::SUCCESS;

    if(m_debug&&(abs(m_eventNoSelection)<1e8)&&m_eventNo!=m_eventNoSelection){
        ++m_eventNo;
        return sc;
    }

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    if(m_tuple) start=std::chrono::high_resolution_clock::now();

    ///retrieve silicon Track and TrackHits
    const edm4hep::TrackCollection* sdtTrackCol=nullptr;
    if(!m_useTrackFinding){
        if(m_SDTTrackCol.exist())sdtTrackCol=m_SDTTrackCol.get();
    } else {
        if(m_SDTTrackFindCol.exist())sdtTrackCol=m_SDTTrackFindCol.get();
    }
    if(nullptr==sdtTrackCol || sdtTrackCol->size()<=0) {
        debug()<<"TrackCollection not found or sdtTrackCol size=0"<<endmsg;
        return StatusCode::SUCCESS;
    }

    auto assoDCHitsCol=m_DCHitAssociationCol.get();
    double eventStartTime=0;

    const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
    dCDigiCol=m_DCDigiCol.get();

    const edm4hep::TrackCollection* dcTrackCol=nullptr;
    if(m_dcTrackCol.exist()) dcTrackCol=m_dcTrackCol.get();
    if(nullptr==dcTrackCol) {
        debug()<<"TrackCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
    mcParticleCol=m_mcParticleCol.get();//FIXME get error when call exist()
    if(nullptr==mcParticleCol){
        debug()<<"MCParticleCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }

    //SET SimTrackerHits
    const edm4hep::TrackerHitCollection* SETHits = m_SETHitCol.get();

    ///----------------------------------------------------
    ///Loop over Track and do fitting for each track
    ///----------------------------------------------------
    m_firstTuple=true;
    debug()<<"SDTTrackCol size="<<sdtTrackCol->size()<<endmsg;
    int iSdtTrack = 0;
    int nDC,nSDT,nFittedSDT,nFittedDC,ngenfitHit;
    std::vector<double> trackL;
    std::vector<double> hitMom;
    std::vector<float> truthMomEdep;
    std::vector<double> driftDis;
    std::vector<double> FittedDoca;
    std::vector<double> truthDoca;
    std::vector<double> Res;
    std::vector<double> truthRes;
    std::vector<int> isRecNoise;
    std::vector<int> isRecFitNoise;
    std::vector<int> layerFit;
    std::vector<int> wireIDFit;
    std::vector<int> layerIDRec;
    std::vector<int> wireIDRec;
    std::vector<double> extrapoLen;

    std::vector<double> docaTrack;
    std::vector<int> isNoise;
    for(auto sdtTrack: *sdtTrackCol)
    {
        if(sdtTrack.trackerHits_size()<1e-9) continue;

        //------------------------------------
        //Add SET hit
        //------------------------------------
        if(m_addSET) addSETHitsToTk(sdtTrack,m_SETHitAssociationCol.get(),SETHits);

        ///Loop over 5 particle hypothesis(0-4): e,mu,pi,K,p
        ///-1 for chargedgeantino
        for(unsigned int pidType=0;pidType<m_nPDG;pidType++)
        {
            if((m_debugPid>=0) && (m_debugPid!=pidType)) continue;
            debug()<<"processing pidType "<<pidType<<endmsg;
            ///-----------------------------------
            ///Create a GenFit track
            ///-----------------------------------
            GenfitTrack* genfitTrack=new GenfitTrack(m_genfitField,
                    m_gridDriftChamber,m_geomSvc);
            genfitTrack->setDebug(m_debug);

            if(!genfitTrack->createGenfitTrackFromEDM4HepTrack(pidType,
                        sdtTrack, eventStartTime,m_isUseCovTrack)){
                debug()<<"createGenfitTrackFromEDM4HepTrack from SDT track failed!"<<endmsg;
                delete genfitTrack;
                continue;
            }
            //}

            ///-----------------------------------
            ///Add hits on track
            ///-----------------------------------
            if(m_debug) std::cout<<" m_measurementTypeSi "<<m_measurementTypeSi<<" "<<m_measurementTypeDC<<" "<<std::endl;
            int nHitAdded=0;
            //add silicon hits
            if(0==m_measurementTypeSi.value()){
                nHitAdded+=genfitTrack->addSpacePointsSi(sdtTrack,
                        m_sigmaHitU,m_sigmaHitV);
            }else if(1==m_measurementTypeSi.value()){
                nHitAdded+=genfitTrack->addSiliconMeasurements(sdtTrack,
                        m_sigmaHitU,m_sigmaHitV);
            }

            //add DC hits
            if(0==m_measurementTypeDC.value()){
                nHitAdded+=genfitTrack->addSpacePointsDC(sdtTrack,
                        assoDCHitsCol,m_sigmaHitU,m_sigmaHitV);
            }else if(1==m_measurementTypeDC.value()){
                if(m_selectDCHit){
                    std::vector<edm4hep::TrackerHit*> selectedHits;
                    selectHits(sdtTrack,selectedHits);
                    nHitAdded+=genfitTrack->addWireMeasurementsFromList(selectedHits,
                            m_sigmaHitU[0],assoDCHitsCol,m_sortMethod,m_truthAmbig,
                            m_skipCorner,m_skipNear);//mm
                    std::vector<edm4hep::TrackerHit*> tmp;
                    selectedHits.swap(tmp);
                }else{
                    nHitAdded+=genfitTrack->addWireMeasurementsOnTrack(sdtTrack,
                            m_sigmaHitU[0],assoDCHitsCol,m_sortMethod,m_truthAmbig,
                            m_skipCorner,m_skipNear);//mm
                }
            }


            // skip events w.o hits
            if(0==nHitAdded){
                debug()<<m_eventNo<<" No hit added to track!"<<endmsg;
                delete genfitTrack;
                continue;
            }
            if(m_debug) genfitTrack->printSeed();

            ///-----------------------------------
            ///call genfit fitting procedure
            ///-----------------------------------
            m_genfitFitter->setDebug(m_debug);
            m_genfitFitter->setDebugGenfit(m_debugGenfit);
            m_genfitFitter->processTrack(genfitTrack,m_resortHits.value());

            //------
            //salvage hits
            //single particle
            //------
            int nAddHits = 0;
            int pdg = genfitTrack->getPDG();
            m_PDG[iSdtTrack] = pdg;
            int numDC;

            if(m_useSalHits){
                nAddHits = genfitTrack->salvageHits(pdg,mcParticleCol,dCDigiCol,
                        assoDCHitsCol,m_sigmaHitU[0],m_truthAmbig,m_skipCorner,m_skipNear,
                        m_extraCut,m_extraDocaCut,docaTrack,isNoise,numDC);
                bool status = genfitTrack->sortedHitOnTrack(pdg);
                m_genfitFitter->processTrackWithRep(genfitTrack,0,m_resortHits);

                int iFiter =0;
                int dcFitNum;
                while(iFiter<m_iterFit){
                    double chi2Fit = genfitTrack->checkGenfitTrack(dcFitNum);
                    if(chi2Fit<m_chi2FitCut && dcFitNum > m_dcFitNum) break;
                    m_genfitFitter->processTrackWithRep(genfitTrack,0,m_resortHits);
                    iFiter++;
                }
            }

            m_numDCOnTrack[iSdtTrack] = numDC;
            m_salDCHits[iSdtTrack] = nAddHits;
            m_DCHitsCol[iSdtTrack] = dCDigiCol->size();

            ///-----------------------------------
            ///Store track
            ///-----------------------------------
            auto dcRecParticle=sdtRecParticleCol->create();
            auto dcRecTrack=sdtRecTrackCol->create();

            TVector3 pocaToOrigin_pos,pocaToOrigin_mom;
            TMatrixDSym pocaToOrigin_cov;
            edm4hep::TrackState pocaToOrigin_trackState;
            if(!genfitTrack->storeTrack(dcRecParticle,dcRecTrack,
                        pocaToOrigin_pos,pocaToOrigin_mom,pocaToOrigin_cov,
                        pidType,m_ndfCut,m_chi2Cut,nDC,nSDT,nFittedDC,nFittedSDT,
                        ngenfitHit,trackL,hitMom,truthMomEdep,assoDCHitsCol,
                        driftDis,FittedDoca,truthDoca,Res,truthRes,isRecNoise,isRecFitNoise,
                        layerFit,wireIDFit,layerIDRec,wireIDRec,extrapoLen)){
                debug()<<"Fitting failed!"<<std::endl;
            }else{
                ++m_fitSuccess[pidType];
            }

            if(m_tuple) debugTrack(iSdtTrack,pidType,genfitTrack,pocaToOrigin_pos,
                    pocaToOrigin_mom,pocaToOrigin_cov);
            if(m_showDisplay) {
                m_genfitDisplay->addEvent(genfitTrack->getTrack());
                m_genfitDisplay->open();

                using namespace std::chrono_literals;
                std::this_thread::sleep_for(1000000000ms);
                system("pause");
            }else{
                delete genfitTrack;
            }
        }//end loop over particle type
        if(m_tuple) {
            m_nTrackerHitDC[iSdtTrack] = nDC;
            m_nTrackerFitHitDC[iSdtTrack] = nFittedDC;
            m_nTrackerHitSDT[iSdtTrack] = nSDT;
        }
        ++iSdtTrack;
    }//end loop over a track
    m_nRecTrack++;
    
    if(m_tuple) {
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        debug() << "Elapsed time: " << elapsed.count() << " s"<<endmsg;
        debugEvent(sdtTrackCol,sdtRecTrackCol,eventStartTime,nFittedSDT);
    }
    
    if(m_tuple) sc=m_tuple->write();
    
    ++m_eventNo;
    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgSDT::finalize()
{
    MsgStream log(msgSvc(), name());
    info()<< " RecGenfitAlgSDT in finalize()" << endmsg;

    m_genfitFitter->writeHist();
    delete m_genfitFitter;
    info()<<"RecGenfitAlgSDT nRecTrack="<<m_nRecTrack<<" success e "
        <<m_fitSuccess[0]<<" mu "<<m_fitSuccess[1]<<" pi "<<m_fitSuccess[2]
        <<" K "<<m_fitSuccess[3]<<" p "<<m_fitSuccess[4]<<std::endl;
    if(m_nRecTrack>0){
        std::cout<<"RecGenfitAlgSDT Success rate = "<<std::endl;
        for (int i=0;i<5;i++){
            std::cout<<Form("%d: %d/%d= %2.2f",i,m_fitSuccess[i],m_nRecTrack,
                    ((float) m_fitSuccess[i])/m_nRecTrack)<<std::endl;
        }
    }
    return StatusCode::SUCCESS;
}

void RecGenfitAlgSDT::debugTrack(int iStrack,int pidType,const GenfitTrack* genfitTrack,
        TVector3 pocaToOrigin_pos,TVector3 pocaToOrigin_mom,
        TMatrixDSym pocaToOrigin_cov)
{
    /// Get fit status
    const genfit::FitStatus* fitState = genfitTrack->getFitStatus();
    int charge= fitState->getCharge();

    if(m_firstTuple){
        m_nHitKalInput=genfitTrack->getNumPoints();
        debug()<<"m_nHitKalInput "<<m_nHitKalInput<<endmsg;
        //FIXME read from config file
        if(m_debug) {
            debug()<<"detType nHot: ";
        }
        int detIDs[5]={1,2,3,4,5};//VXD=1,SIT=2,SET=5;FTD=3,
        for(int i=0;i<5;i++){
            m_nHitDetType[detIDs[i]]=genfitTrack->getNumPointsDet(detIDs[i]);
            if(m_debug){
                debug()<<" "<<detIDs[i]<<"="<<m_nHitDetType[detIDs[i]]<<", ";
            }
        }
        if(m_debug) { debug()<<endmsg; }
        m_firstTuple=false;
    }
    m_chargeKal[iStrack][pidType]= charge;
    m_nHitWithFitInfo[iStrack][pidType]=genfitTrack->getNumPointsWithFittedInfo();
    m_chi2Kal[iStrack][pidType]=fitState->getChi2();
    m_nDofKal[iStrack][pidType]=fitState->getNdf();
    m_isFitted[iStrack][pidType]=(int)fitState->isFitted();
    m_isFitConverged[iStrack][pidType]=(int) fitState->isFitConverged();
    m_isFitConvergedFully[iStrack][pidType]=(int) fitState->isFitConvergedFully();

    ///get fitted state of track
    TMatrixDSym fittedCov;
    TLorentzVector fittedPos;
    TVector3 fittedMom;
    int fittedState=genfitTrack->getFittedState(fittedPos,fittedMom,fittedCov);
    const TLorentzVector seedPos=genfitTrack->getSeedStatePos();
    const TVector3 seedMom=genfitTrack->getSeedStateMom();
    m_fittedState[iStrack][pidType]=fittedState;
    HelixClass helix;//mm and GeV
    HelixClass helix_origin;//mm and GeV


    double pos[3]={(fittedPos.X()/dd4hep::mm),(fittedPos.Y()/dd4hep::mm),
        (fittedPos.Z()/dd4hep::mm)};
    double mom[3]={(fittedMom.X()),(fittedMom.Y()),(fittedMom.Z())};

    m_posx[iStrack] = fittedPos.X();
    m_posy[iStrack] = fittedPos.Y();
    m_posz[iStrack] = fittedPos.Z();

    m_momx[iStrack] = fittedMom.X();
    m_momy[iStrack] = fittedMom.Y();
    m_momz[iStrack] = fittedMom.Z();

    m_PosMcX[iStrack] = seedPos.X();
    m_PosMcY[iStrack] = seedPos.Y();
    m_PosMcZ[iStrack] = seedPos.Z();

    m_MomMcX[iStrack] = seedMom.X();
    m_MomMcY[iStrack] = seedMom.Y();
    m_MomMcZ[iStrack] = seedMom.Z();

    m_PocaPosX[iStrack] = pocaToOrigin_pos.X()*dd4hep::mm;
    m_PocaPosY[iStrack] = pocaToOrigin_pos.Y()*dd4hep::mm;
    m_PocaPosZ[iStrack] = pocaToOrigin_pos.Z()*dd4hep::mm;

    m_PocaMomX[iStrack] = pocaToOrigin_mom.X();
    m_PocaMomY[iStrack] = pocaToOrigin_mom.Y();
    m_PocaMomZ[iStrack] = pocaToOrigin_mom.Z();

    for(int i=0;i<6;i++)
    {
        m_ErrorcovMatrix6[iStrack][i] = fittedCov(i,i);
        m_PocaErrCov[iStrack][i] = pocaToOrigin_cov(i,i);
    }

    double pocaToOrigin_Pos[3] = {pocaToOrigin_pos.X(),pocaToOrigin_pos.Y(),pocaToOrigin_pos.Z()};
    double pocaToOrigin_Mom[3] = {pocaToOrigin_mom.X(),pocaToOrigin_mom.Y(),pocaToOrigin_mom.Z()};
    TLorentzVector pocaToOrigin_POS;
    pocaToOrigin_POS.SetXYZT(pocaToOrigin_pos.X()*dd4hep::mm,pocaToOrigin_pos.Y()*dd4hep::mm,
            pocaToOrigin_pos.Z()*dd4hep::mm,999);
    helix.Initialize_VP(pos,mom,charge,m_genfitField->getBz(fittedPos.Vect())/GenfitUnit::tesla);
    helix_origin.Initialize_VP(pocaToOrigin_Pos,pocaToOrigin_Mom,charge,m_genfitField->getBz(pocaToOrigin_POS.Vect())/GenfitUnit::tesla);
    m_pocaMomKalP[iStrack][pidType]=fittedMom.Mag();
    m_pocaMomKalPt[iStrack][pidType]=fittedMom.Perp();

    TMatrixDSym covMatrix_6=pocaToOrigin_cov;
    for(int i=0;i<5;i++){
        covMatrix_6[0][i]=pocaToOrigin_cov[0][i]/dd4hep::mm;//d0 column
        covMatrix_6[1][i]=pocaToOrigin_cov[1][i]/dd4hep::mm;//omega column
        covMatrix_6[2][i]=pocaToOrigin_cov[2][i]/dd4hep::mm;//z0 column
        covMatrix_6[i][0]=pocaToOrigin_cov[i][0]/dd4hep::mm;//d0 row
        covMatrix_6[i][1]=pocaToOrigin_cov[i][1]/dd4hep::mm;//omega row
        covMatrix_6[i][2]=pocaToOrigin_cov[i][2]/dd4hep::mm;//z0 row
    }
    edm4hep::TrackState trackState_Origin;
    CEPC::getTrackStateFromPosMom(trackState_Origin,m_genfitField->getBz(pocaToOrigin_POS.Vect())/GenfitUnit::tesla,pocaToOrigin_pos,
            pocaToOrigin_mom,charge,covMatrix_6);
    std::array<float,21> errorCov_Origin;
    errorCov_Origin = trackState_Origin.covMatrix;
    for(int j=0; j<15; j++) {
        m_ErrorcovMatrix_Origin[iStrack][j] = errorCov_Origin[j];
    }
    m_D0_Origin[iStrack] = helix_origin.getD0();
    m_phi_Origin[iStrack] = helix_origin.getPhi0();
    m_omega_Origin[iStrack] = helix_origin.getOmega();
    m_Z0_Origin[iStrack] = helix_origin.getZ0();
    m_tanLambda_Origin[iStrack] = helix_origin.getTanLambda();

    m_fittedXc[iStrack] = helix.getXC();
    m_fittedYc[iStrack] = helix.getYC();
    m_fittedR[iStrack] = helix.getRadius();

    m_evt=m_eventNo;
    /// Get fit status
    if((0!=fittedState)||(!m_isFitted[pidType])||(m_nDofKal[iStrack][pidType]>m_ndfCut)){
        debug()<<"evt "<<m_evt<<" fit FAILED !!"
            <<pidType<<" fittedState "<<fittedState<<" isFitted "
            <<m_isFitted[pidType]<<" isConverged "<<m_isFitConverged[pidType]
            <<" isFitConvergedFully "<<m_isFitConvergedFully[pidType]<<endmsg;
    }else{
        debug()<<"==fit result evt "<<m_evt<<" pidType "<<pidType<<" pos("<<
            fittedPos.X()<<" "<<
            fittedPos.Y()<<" "<<
            fittedPos.Z()<<")cm mom("<<
            fittedMom.X()<<" "<<
            fittedMom.Y()<<" "<<
            fittedMom.Z()<<") p_tot "<<
            fittedMom.Mag()<<" pt "<<
            fittedMom.Perp()
            <<" fittedState "<<fittedState<<" isFitted "
            <<m_isFitted[pidType]<<" isConverged "<<m_isFitConverged[pidType]
            <<" isFitConvergedFully "<<m_isFitConvergedFully[pidType]
            <<" ndf "<<m_nDofKal[iStrack][pidType]
            <<" chi2 "<<m_chi2Kal[pidType]<<endmsg;
    }
}

void RecGenfitAlgSDT::debugEvent(const edm4hep::TrackCollection* sdtTrackCol,
        const edm4hep::TrackCollection* sdtRecTrackCol,
        double eventStartTime,int nFittedSDT)
{
    int iSdtTrack=0;
    m_nSdtTrack=sdtTrackCol->size();
    for(auto sdtTrack: *sdtTrackCol){
        m_nSdtTrackHit = sdtTrack.trackerHits_size();
        if(sdtTrack.trackerHits_size()<1e-9) continue;
        for(int ihit=0;ihit<sdtTrack.trackerHits_size();ihit++){
            edm4hep::TrackerHit sdtTrackHit = sdtTrack.getTrackerHits(ihit);
        }

        //if(iSdtTrack>0) break;//TODO debug for some track only
        edm4hep::TrackState trackStat=sdtTrack.getTrackStates(0);//FIXME?
        HelixClass helixClass;
        helixClass.Initialize_Canonical(trackStat.phi,trackStat.D0,
                trackStat.Z0,trackStat.omega,trackStat.tanLambda,
                m_genfitField->getBz({0.,0.,0.})/GenfitUnit::tesla);

        TLorentzVector posInit(helixClass.getReferencePoint()[0],
                helixClass.getReferencePoint()[1],
                helixClass.getReferencePoint()[2],eventStartTime);
        m_seedPos[iSdtTrack][0]=posInit.X();
        m_seedPos[iSdtTrack][1]=posInit.Y();
        m_seedPos[iSdtTrack][2]=posInit.Z();
        TVector3 momInit(helixClass.getMomentum()[0],
                helixClass.getMomentum()[1],helixClass.getMomentum()[2]);
        m_seedMomP[iSdtTrack]=momInit.Mag();
        m_seedMomPt[iSdtTrack]=momInit.Perp();
        m_seedMom[iSdtTrack][0]=momInit.X();
        m_seedMom[iSdtTrack][1]=momInit.Y();
        m_seedMom[iSdtTrack][2]=momInit.Z();
        TVector3 pos,mom;
        TMatrixDSym cov(6);
        double charge;
        CEPC::getPosMomFromTrackState(trackStat,
                m_genfitField->getBz({0.,0.,0.})/GenfitUnit::tesla,pos,mom,charge,cov);
        for(int i =0;i<6;i++)
        {
            m_McErrCov[iSdtTrack][i] = cov(i,i);
        }
        m_seedMomQ[iSdtTrack]=charge;
        iSdtTrack++;
    }

    const edm4hep::MCParticleCollection* mcParticleCol = nullptr;
    const edm4hep::SimTrackerHitCollection* simHitCol=nullptr;

    m_pidIndex=5;

    if((m_mcParticleCol.get())!=nullptr){
        mcParticleCol=m_mcParticleCol.get();
        int iMcParticle=0;
        HelixClass helix_mcP;
        for(auto mcParticle : *mcParticleCol){
            edm4hep::Vector3f mcPocaMom = mcParticle.getMomentum();//GeV
            edm4hep::Vector3d mcPocaPos = mcParticle.getVertex();

            double mcPos[3]={(mcPocaPos.x),(mcPocaPos.y),(mcPocaPos.z)};
            double mcMom[3]={(mcPocaMom.x),(mcPocaMom.y),(mcPocaMom.z)};
            for(int i=0;i<3;i++){debug()<<"mcPos "<<mcPos[i]<<endmsg;}
            for(int i=0;i<3;i++){debug()<<"mcMom "<<mcMom[i]<<endmsg;}
            float mcCharge = mcParticle.getCharge();
            helix_mcP.Initialize_VP(mcPos,mcMom,mcCharge,
                    m_genfitField->getBz(mcPos)/GenfitUnit::tesla);

            mcP_D0[iMcParticle] = helix_mcP.getD0();
            mcP_phi[iMcParticle] = helix_mcP.getPhi0();
            mcP_omega[iMcParticle] = helix_mcP.getOmega();
            mcP_Z0[iMcParticle] = helix_mcP.getZ0();
            mcP_tanLambda[iMcParticle] = helix_mcP.getTanLambda();

            mcP_Xc[iMcParticle] = helix_mcP.getXC();
            mcP_Yc[iMcParticle] = helix_mcP.getYC();
            mcP_R[iMcParticle] = helix_mcP.getRadius();

            debug()<< "debugEvent Bz " << m_genfitField->getBz(mcPos)/GenfitUnit::tesla
                << "Tesla mc d0= " << mcP_D0
                << " phi0= " << mcP_phi
                << " omega= " << mcP_omega
                << " Z0= " << mcP_Z0
                << " tanLambda= " << mcP_tanLambda << endmsg;

            float px=mcPocaMom.x;
            float py=mcPocaMom.y;
            float pz=mcPocaMom.z;
            debug()<<"mc pos("<<mcPos[0]<<","<<mcPos[1]<<","<<mcPos[2]
                <<") pxyz("<<px<<","<<py<<","<<pz<<") p"<<sqrt(px*px+py*py+pz*pz)
                <<endmsg;
            m_pocaMomMcP[iMcParticle]=sqrt(px*px+py*py+pz*pz);
            m_pocaMomMcPt[iMcParticle]=sqrt(px*px+py*py);
            m_pocaMomMc[iMcParticle][0]=px;
            m_pocaMomMc[iMcParticle][1]=py;
            m_pocaMomMc[iMcParticle][2]=pz;
            iMcParticle++;
        }
        m_mcIndex=iMcParticle;
    }

    m_nSdtRecTrack=sdtRecTrackCol->size();
    int isdttrack=0;
    for(auto sdtRecTrack: *sdtRecTrackCol){
        for(int iHit=0;iHit<sdtRecTrack.trackerHits_size();iHit++)
        {
            edm4hep::TrackerHit sdtRecTrackHit = sdtRecTrack.getTrackerHits(iHit);
        }
        for(unsigned int i=0; i<sdtRecTrack.trackStates_size(); i++) {
            edm4hep::TrackState trackStat=sdtRecTrack.getTrackStates(i);
            std::array<float,21> errorCov;
            errorCov = trackStat.covMatrix;
            for(int j=0; j<15; j++) {
                m_ErrorcovMatrix[isdttrack][j] = errorCov[j];
                if(m_debug)debug()<<"errorCov "<<j<<" "<<errorCov[j]<<endmsg;
            }
            m_D0[isdttrack] = trackStat.D0;
            m_phi[isdttrack] = trackStat.phi;
            m_omega[isdttrack] = trackStat.omega;
            m_Z0[isdttrack] = trackStat.Z0;
            m_tanLambda[isdttrack] = trackStat.tanLambda;

            HelixClass helix;
            helix.Initialize_Canonical(trackStat.phi,trackStat.D0,
                    trackStat.Z0,trackStat.omega,trackStat.tanLambda,
                    m_genfitField->getBz({0.,0.,0.})/GenfitUnit::tesla);

            m_trackXc[isdttrack] = helix.getXC();
            m_trackYc[isdttrack] = helix.getYC();
            m_trackR[isdttrack] = helix.getRadius();

        }
        ++isdttrack;
    }

    //truth digi hit
    const edm4hep::TrackerHitCollection* TruthdCDigiCol=m_TruthDCDigiCol.get();
    m_nTruthDCDigi=TruthdCDigiCol->size();

    //debug digi
    const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
    dCDigiCol=m_DCDigiCol.get();
    if(nullptr!=dCDigiCol){ m_nDCDigi=dCDigiCol->size(); }
    int iDCDigi=0;
    for(auto dcDigi: *dCDigiCol){
        m_dcDigiChamber[iDCDigi]=m_decoder->get(dcDigi.getCellID(),"chamber");
        m_dcDigiLayer[iDCDigi]=m_decoder->get(dcDigi.getCellID(),"layer");
        m_dcDigiCell[iDCDigi]=m_decoder->get(dcDigi.getCellID(),"cellID");
        m_dcDigiTime[iDCDigi]=dcDigi.getTime();
        m_dcDigiDrift[iDCDigi]=dcDigi.getTime()*m_driftVelocity.value()/10000.; //cm
        TVector3 endPointStart(0,0,0);
        TVector3 endPointEnd(0,0,0);
        m_gridDriftChamber->cellposition(dcDigi.getCellID(),endPointStart,
                endPointEnd);
        m_dcDigiWireStartX[iDCDigi]=endPointStart.X();
        m_dcDigiWireStartY[iDCDigi]=endPointStart.Y();
        m_dcDigiWireStartZ[iDCDigi]=endPointStart.Z();
        m_dcDigiWireEndX[iDCDigi]=endPointEnd.X();
        m_dcDigiWireEndY[iDCDigi]=endPointEnd.Y();
        m_dcDigiWireEndZ[iDCDigi]=endPointEnd.Z();



        //get information from associated simTrackerHit
        auto dcSimTrackerHit=CEPC::getAssoClosestSimTrackerHit(m_DCHitAssociationCol.get(),dcDigi,m_gridDriftChamber,0);
        m_dcDigiMcMomX[iDCDigi]=dcSimTrackerHit.getMomentum().x*GenfitUnit::GeV;
        m_dcDigiMcMomY[iDCDigi]=dcSimTrackerHit.getMomentum().y*GenfitUnit::GeV;
        m_dcDigiMcMomZ[iDCDigi]=dcSimTrackerHit.getMomentum().z*GenfitUnit::GeV;
        m_dcDigiMcPosX[iDCDigi]=dcSimTrackerHit.getPosition().x*GenfitUnit::mm;
        m_dcDigiMcPosY[iDCDigi]=dcSimTrackerHit.getPosition().y*GenfitUnit::mm;
        m_dcDigiMcPosZ[iDCDigi]=dcSimTrackerHit.getPosition().z*GenfitUnit::mm;
        m_dcDigiDocaMC[iDCDigi]=dcDigi.getTime()*m_driftVelocity.value()/10000.;//cm
        TVector3 pocaOnWire=m_gridDriftChamber->wirePos_vs_z(dcDigi.getCellID(),
                dcSimTrackerHit.getPosition().z*dd4hep::mm);
        m_dcDigiPocaOnWireMCX[iDCDigi]=pocaOnWire.X();
        m_dcDigiPocaOnWireMCY[iDCDigi]=pocaOnWire.Y();
        double firstMom=sqrt(dcSimTrackerHit.getMomentum().x*
                dcSimTrackerHit.getMomentum().x+dcSimTrackerHit.getMomentum().y
                *dcSimTrackerHit.getMomentum().y+dcSimTrackerHit.getMomentum().z
                *dcSimTrackerHit.getMomentum().z);
        if(0==m_decoder->get(dcDigi.getCellID(),"layer")){
            m_firstMomMc=firstMom;
            if(m_debug.value()>0){
                std::cout<<" firstMomMc "<<firstMom<<" ("
                    <<dcSimTrackerHit.getMomentum().x<<","
                    <<dcSimTrackerHit.getMomentum().y
                    <<","<<dcSimTrackerHit.getMomentum().z<<")"<<std::endl;
            }
        }
        if(m_debug) std::cout<<"digi "<<iDCDigi<<" ("
            <<m_decoder->get(dcDigi.getCellID(),"layer")<<","
                <<m_decoder->get(dcDigi.getCellID(),"cellID")<<") truth mom GeV"
                <<firstMom<<" "<<dcSimTrackerHit.getMomentum().x<<" "
                <<dcSimTrackerHit.getMomentum().y<<" "
                <<dcSimTrackerHit.getMomentum().z<<" truth pos "
                <<dcSimTrackerHit.getPosition().x*GenfitUnit::mm<<" "
                <<dcSimTrackerHit.getPosition().y*GenfitUnit::mm<<" "
                <<dcSimTrackerHit.getPosition().z*GenfitUnit::mm<<"cm truth doca "
                <<dcDigi.getTime()*m_driftVelocity.value()/10000.<<" cm"
                <<std::endl;//yzhang debug

        iDCDigi++;
    }
}

void RecGenfitAlgSDT::selectHits(const edm4hep::Track&,
        std::vector<edm4hep::TrackerHit*>& dcDigiSelected)
{
    //for single track only, FIXME
    double eventStartTime=0;
    unsigned int pidType=1;//mu
    GenfitTrack* genfitTrack=new GenfitTrack(m_genfitField,
            m_gridDriftChamber,m_geomSvc);
    genfitTrack->setDebug(m_debug);
    const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
    mcParticleCol=m_mcParticleCol.get();//FIXME get error when call exist()
    if(genfitTrack->createGenfitTrackFromMCParticle(
                pidType,*(mcParticleCol->begin()),
                eventStartTime)){
        const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
        dCDigiCol=m_DCDigiCol.get();
        int iDCDigi=0;
        for(auto dcDigi:*dCDigiCol){
            TVector3 poca,pocaDir,pocaOnWire;

            double docaExt=1e9;
            bool stopAtBoundary=false;
            bool calcJacobianNoise=true;
            edm4hep::MCParticle mcParticle=*(mcParticleCol->begin());//FIXME single track only

            genfitTrack->extrapolateToHit(poca,pocaDir,pocaOnWire,docaExt,
                    mcParticle,dcDigi.getCellID(),pidType,stopAtBoundary,calcJacobianNoise);
            m_dcDigiDocaExt[iDCDigi]=docaExt;
            m_dcDigiPocaExtX[iDCDigi]=poca.X();
            m_dcDigiPocaExtY[iDCDigi]=poca.Y();
            m_dcDigiPocaExtZ[iDCDigi]=poca.Z();
            double docaMC=dcDigi.getTime()*m_driftVelocity.value()/10000.;//cm
            auto dcSimTrackerHit=CEPC::getAssoClosestSimTrackerHit(
                    m_DCHitAssociationCol.get(),dcDigi,m_gridDriftChamber,0);

            debug()<<" CellID "<<dcDigi.getCellID()<<"select hit ("
                <<m_decoder->get(dcDigi.getCellID(),"layer")
                <<","<<m_decoder->get(dcDigi.getCellID(),"cellID")
                <<") by extdoca:"<<docaExt<<"- docaMC:"<<docaMC<<" deltaDoca "
                <<fabs(docaExt-docaMC)<<" cut "<<m_docaCut
                <<" pocaOnWire ("<<pocaOnWire.X()<<","<<pocaOnWire.Y()
                <<","<<pocaOnWire.Z()<<") truthPos("
                <<dcSimTrackerHit.getPosition().x/10.<<","
                <<dcSimTrackerHit.getPosition().y/10.<<","
                <<dcSimTrackerHit.getPosition().z/10.<<") diffZ "
                <<dcSimTrackerHit.getPosition().z/10.-pocaOnWire.Z()<<endmsg;
            if(fabs(docaExt-docaMC)>m_docaCut*GenfitUnit::mm){
                debug()<<"Skip hit delta doca "<<fabs(docaExt-docaMC)<<endmsg;
                continue;
            }
            edm4hep::TrackerHit* thisDigi = new edm4hep::TrackerHit(dcDigi);
            dcDigiSelected.push_back(thisDigi);
            iDCDigi++;
        }//end loop over digi
        if(m_debug>0){
            std::cout<<"selectHits "<<dCDigiCol->size()
                <<" after "<<iDCDigi<<std::endl;
        }
    }//end loop over track
    delete genfitTrack;
}//end of select hit

bool RecGenfitAlgSDT::addHitsToTk(edm4hep::TrackerHit hit,edm4hep::Track& track,
        const char* msg) const
{
    edm4hep::MutableTrack _track = track.clone();
    _track.addToTrackerHits(hit);
    track = _track;
    return true;
}

void RecGenfitAlgSDT::addSETHitsToTk(edm4hep::Track& track,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        const edm4hep::TrackerHitCollection* SETHits)
{
    TVector3 SETPos,SETMom;

    edm4hep::TrackState trackStat=track.getTrackStates(0);//FIXME?
    edm4hep::TrackerHit trackHit  = track.getTrackerHits(0);
    HelixClass helixClass;
    helixClass.Initialize_Canonical(trackStat.phi,trackStat.D0,
            trackStat.Z0,trackStat.omega,trackStat.tanLambda,
            m_genfitField->getBz({0.,0.,0.})/GenfitUnit::tesla);
    double SETRadius = 1815.1; // mm;FIXME
    float refPos[3],Point[6];
    float helixMom[3];
    helixClass.getPointOnCircle(SETRadius,refPos,Point); //mm
    helixClass.getExtrapolatedMomentum(Point,helixMom);//mm
    SETPos.SetXYZ(Point[0],Point[1],Point[2]);
    SETMom.SetXYZ(helixMom[0],helixMom[1],helixMom[2]);
    for(int i=0;i<SETHits->size();i++){
        edm4hep::TrackerHit SETHit = SETHits->at(i);
        TVector3 SETHitPos(SETHit.getPosition()[0],
                SETHit.getPosition()[1],SETHit.getPosition()[2]);
        if(abs(SETPos.X()-SETHitPos.X())<m_cutSETPos &&
                abs(SETPos.Y()-SETHitPos.Y())<m_cutSETPos &&
                abs(SETPos.Z()-SETHitPos.Z())<m_cutSETPos ){
            const edm4hep::SimTrackerHit simTrackerhit =
                CEPC::getAssoSimTrackerHit(assoHits,SETHit);
            if(addHitsToTk(SETHit,track,"SET digi hit")){
                debug()<<"add SET digi hit 1  trackerHit"<<endmsg;
            }

        }
    }
}

