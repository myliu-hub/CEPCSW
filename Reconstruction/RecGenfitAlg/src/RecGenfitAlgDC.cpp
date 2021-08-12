#include "RecGenfitAlgDC.h"
#include "GenfitTrack.h"
#include "GenfitFitter.h"
#include "GenfitField.h"
#include "AbsMeasurement.h"

//genfit
#include "EventDisplay.h"

//cepcsw
#include "DetInterface/IGeomSvc.h"
#include "DataHelper/HelixClass.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "DataHelper/TrackHelper.h"
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
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "UTIL/BitField64.h"
#include "DDSegmentation/Segmentation.h"
#include "TRandom.h"
#include "TLorentzVector.h"

//stl
#include <chrono>
#include "time.h"


DECLARE_COMPONENT( RecGenfitAlgDC )

/////////////////////////////////////////////////////////////////////
RecGenfitAlgDC::RecGenfitAlgDC(const std::string& name, ISvcLocator* pSvcLocator):
        GaudiAlgorithm(name, pSvcLocator),m_nPDG(5),m_dd4hep(nullptr),
        m_gridDriftChamber(nullptr),m_decoder(nullptr)
{
    //declareProperty("EventHeaderCol", _headerCol);
    declareProperty("MCParticle", m_mcParticleCol,
            "Handle of the input MCParticle collection");
    declareProperty("DriftChamberHitsCollection", m_simDCHitCol,
            "Handle of the input SimHit collection");
    declareProperty("DigiDCHitCollection", m_digiDCHitsCol,
            "Handle of digi DCHit collection");
    declareProperty("DCTrackCollection", m_dcTrackCol,
            "Handle of DC track collection");
    declareProperty("DCHitAssociationCollection", m_dcHitAssociationCol,
            "Handle of simTrackerHit and TrackerHit association collection");
    declareProperty("DCRecParticleCollection", m_dcRecParticleCol,
            "Handle of drift chamber reconstructed particle collection");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgDC::initialize()
{
    MsgStream log(msgSvc(), name());
    info()<<"RecGenfitAlgDC initialize()"<<endmsg;

    ///Get GeomSvc
    m_geomSvc=Gaudi::svcLocator()->service("GeomSvc");
    if (nullptr==m_geomSvc) {
        std::cout<<"Failed to find GeomSvc"<<std::endl;
        return StatusCode::FAILURE;
    }
    ///Get Detector
    m_dd4hep = m_geomSvc->lcdd();
    ///Get Field
    m_dd4hepField=m_geomSvc->lcdd()->field();

    /// New a genfit fitter
    m_genfitFitter=new GenfitFitter(m_fitterType.toString().c_str());
    m_genfitField=new GenfitField(m_dd4hepField);
    m_genfitFitter->setField(m_genfitField);
    m_genfitFitter->setGeoMaterial(m_geomSvc->lcdd(),m_extMinDistCut);
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
    m_nDCTrack=0;
    ///Get Readout
    dd4hep::Readout readout=m_dd4hep->readout(m_readout_name);
    ///Get Segmentation
    m_gridDriftChamber=dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>
        (readout.segmentation().segmentation());
    if(nullptr==m_gridDriftChamber){
        error() << "Failed to get the GridDriftChamber" << endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Decoder
    m_decoder = m_geomSvc->getDecoder(m_readout_name);
    if (nullptr==m_decoder) {
        error() << "Failed to get the decoder" << endmsg;
        return StatusCode::FAILURE;
    }


    ///book tuple
    NTuplePtr nt(ntupleSvc(), "RecGenfitAlgDC/recGenfitAlgDC");
    if(nt){
        m_tuple=nt;
    }else{
        m_tuple=ntupleSvc()->book("RecGenfitAlgDC/recGenfitAlgDC",
                CLID_ColumnWiseTuple,"RecGenfitAlgDC");
        if(m_tuple){
            StatusCode sc;
            sc=m_tuple->addItem("run",m_run);
            sc=m_tuple->addItem("evt",m_evt);
            sc=m_tuple->addItem("tkId",m_tkId);
            sc=m_tuple->addItem("nSdtTrack",m_nSdtTrack);
            sc=m_tuple->addItem("nDcTrack",m_nDcTrack);

            sc=m_tuple->addItem("mcIndex",m_mcIndex,0,100);//max. 100 particles
            sc=m_tuple->addItem("seedMomP",m_seedMomP);//for single track debug
            sc=m_tuple->addItem("seedMomPt",m_seedMomPt);
            sc=m_tuple->addItem("seedMomQ",m_seedMomQ);
            sc=m_tuple->addItem("seedMom",3,m_seedMom);
            sc=m_tuple->addItem("seedPos",3,m_seedPos);
            sc=m_tuple->addItem("seedCenterX",m_seedCenterX);
            sc=m_tuple->addItem("seedCenterY",m_seedCenterY);
            sc=m_tuple->addItem("seedR",m_seedR);
            sc=m_tuple->addItem("truthPocaMc",m_mcIndex,m_truthPocaMc,3);
            sc=m_tuple->addItem("pocaPosMc",m_mcIndex,m_pocaPosMc,3);
            sc=m_tuple->addItem("pocaMomMc",m_mcIndex,m_pocaMomMc,3);
            sc=m_tuple->addItem("pocaMomMcP",m_mcIndex,m_pocaMomMcP);
            sc=m_tuple->addItem("pocaMomMcPt",m_mcIndex,m_pocaMomMcPt);
            sc=m_tuple->addItem("pocaPosMdc",3,m_pocaPosMdc);
            sc=m_tuple->addItem("pocaMomMdc",3,m_pocaMomMdc);
            sc=m_tuple->addItem("index",m_pidIndex, 0, 5);
            sc=m_tuple->addItem("firstPosKalP",5,3,m_firstPosKal);
            sc=m_tuple->addItem("firstMomKalP",5,m_firstMomKalP);
            sc=m_tuple->addItem("firstMomKalPt",5,m_firstMomKalPt);

            sc=m_tuple->addItem("ErrorcovMatrix",15,m_ErrorcovMatrix);
            sc=m_tuple->addItem("D0",m_D0);
            sc=m_tuple->addItem("phi",m_phi);
            sc=m_tuple->addItem("omega",m_omega);
            sc=m_tuple->addItem("Z0",m_Z0);
            sc=m_tuple->addItem("tanLambda",m_tanLambda);

            sc=m_tuple->addItem("mcP_D0",mcP_D0);
            sc=m_tuple->addItem("mcP_phi",mcP_phi);
            sc=m_tuple->addItem("mcP_omega",mcP_omega);
            sc=m_tuple->addItem("mcP_Z0",mcP_Z0);
            sc=m_tuple->addItem("mcP_tanLambda",mcP_tanLambda);

            sc=m_tuple->addItem("fitPosx",m_fitPosx);
            sc=m_tuple->addItem("fitPosy",m_fitPosy);
            sc=m_tuple->addItem("fitPosz",m_fitPosz);

            sc=m_tuple->addItem("fitMomx",m_fitMomx);
            sc=m_tuple->addItem("fitMomy",m_fitMomy);
            sc=m_tuple->addItem("fitMomz",m_fitMomz);

            sc=m_tuple->addItem("extraPos",3,m_extraPos);
            sc=m_tuple->addItem("extraMom",3,m_extraMom);

            sc=m_tuple->addItem("Errorcov6",6,m_Error6);

            sc=m_tuple->addItem("nDCDigi",m_nDCDigi,0,50000);

            sc=m_tuple->addItem("pocaPosKal",5,3,m_pocaPosKal);
            sc=m_tuple->addItem("pocaMomKal",5,3,m_pocaMomKal);
            sc=m_tuple->addItem("pocaMomKalPFirst",5,m_pocaMomKalPFirst);
            sc=m_tuple->addItem("pocaMomKalP",5,m_pocaMomKalP);
            sc=m_tuple->addItem("pocaMomKalPt",5,m_pocaMomKalPt);
            sc=m_tuple->addItem("chargeKal",5,m_chargeKal);
            sc=m_tuple->addItem("nDofKal",5,m_nDofKal);
            sc=m_tuple->addItem("chi2Kal",5,m_chi2Kal);
            sc=m_tuple->addItem("isFitted",5,m_isFitted);
            sc=m_tuple->addItem("isFitConverged",5,m_isFitConverged);
            sc=m_tuple->addItem("isFitConvergedFully",5,
                    m_isFitConvergedFully);
            sc=m_tuple->addItem("nHitFailedKal",5,m_nHitFailedKal);
            sc=m_tuple->addItem("nHitFitted",5,m_nHitFitted);
            sc=m_tuple->addItem("nDigi",m_nDigi);
            sc=m_tuple->addItem("nHitMc",m_nHitMc);
            sc=m_tuple->addItem("nHitKalInput",m_nHitKalInput);
            sc=m_tuple->addItem("nHitWithFitInfo",5,m_nHitWithFitInfo);

            sc=m_tuple->addItem("dcDigiChamber",m_nDCDigi,m_dcDigiChamber);
            sc=m_tuple->addItem("dcDigiLayer",m_nDCDigi,m_dcDigiLayer);
            sc=m_tuple->addItem("dcDigiCell",m_nDCDigi,m_dcDigiCell);
            sc=m_tuple->addItem("dcDigiTime",m_nDCDigi,m_dcDigiTime);
            sc=m_tuple->addItem("dcDigiDocaExt",m_nDCDigi,m_dcDigiDocaExt);
            sc=m_tuple->addItem("dcDigiDoca",m_nDCDigi,m_dcDigiDoca);
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

            sc=m_tuple->addItem("nSimDCHit",m_nSimDCHit,0,50000);
            sc=m_tuple->addItem("dcHitDriftT",m_nSimDCHit,m_dcHitDriftT);
            sc=m_tuple->addItem("dcHitDriftDl",m_nSimDCHit,m_dcHitDriftDl);
            sc=m_tuple->addItem("dcHitDriftDr",m_nSimDCHit,m_dcHitDriftDr);
            sc=m_tuple->addItem("dcHitLr",m_nSimDCHit,m_dcHitLr);
            sc=m_tuple->addItem("dcHitLayer",m_nSimDCHit,m_dcHitLayer);
            sc=m_tuple->addItem("dcHitWire",m_nSimDCHit,m_dcHitWire);
            sc=m_tuple->addItem("dcHitExpDoca",m_nSimDCHit,m_dcHitExpDoca);
            sc=m_tuple->addItem("dcHitExpMcDoca",m_nSimDCHit,m_dcHitExpMcDoca);
            sc=m_tuple->addItem("dcHitErr",m_nSimDCHit,m_dcHitErr);
            sc=m_tuple->addItem("time",5,m_time);
            sc=m_tuple->addItem("dcHitMcTkId",m_nSimDCHit,m_dcHitMcTkId);
            sc=m_tuple->addItem("dcHitMcLr",m_nSimDCHit,m_dcHitMcLr);
            sc=m_tuple->addItem("dcHitMcDrift",m_nSimDCHit,m_dcHitMcDrift);
            sc=m_tuple->addItem("dcHitMcX",m_nSimDCHit,m_dcHitMcX);
            sc=m_tuple->addItem("dcHitMcY",m_nSimDCHit,m_dcHitMcY);
            sc=m_tuple->addItem("dcHitMcZ",m_nSimDCHit,m_dcHitMcZ);
            sc=m_tuple->addItem("dcHitMcDoca",m_nSimDCHit,m_dcHitMcDoca);
            sc=m_tuple->addItem("dcHitMcWireX",m_nSimDCHit,m_dcHitMcWireX);
            sc=m_tuple->addItem("dcHitMcWireY",m_nSimDCHit,m_dcHitMcWireY);
            //sc=m_tuple->addItem("dcHitMcLayer",m_nSimDCHit,m_dcHitMcLayer);
            sc=m_tuple->addItem("mcPocaX",m_nSimDCHit,m_dcHitExpMcPocaX);
            sc=m_tuple->addItem("mcPocaY",m_nSimDCHit,m_dcHitExpMcPocaY);
            sc=m_tuple->addItem("mcPocaZ",m_nSimDCHit,m_dcHitExpMcPocaZ);
            sc=m_tuple->addItem("mcPocaWireX",m_nSimDCHit,m_dcHitExpMcPocaWireX);
            sc=m_tuple->addItem("mcPocaWireY",m_nSimDCHit,m_dcHitExpMcPocaWireY);
            sc=m_tuple->addItem("mcPocaWireZ",m_nSimDCHit,m_dcHitExpMcPocaWireZ);
            sc=m_tuple->addItem("genfitTrackNumPoint",m_genfitTrackNumPoint);
            sc=m_tuple->addItem("genfitTrackNumPointsWithMeas",
                    m_genfitTrackNumPointsWithMeas);
            sc=m_tuple->addItem("genfitNHit",m_genfitNHit,0,50000);
            sc=m_tuple->addItem("genfitHitLayer",m_genfitNHit,m_genfitHitLayer);
            sc=m_tuple->addItem("genfitHitCell",m_genfitNHit,m_genfitHitCell);
            sc=m_tuple->addItem("genfitHitDrift",m_genfitNHit,m_genfitHitDrift);
            sc=m_tuple->addItem("genfitHitDriftErr",m_genfitNHit,m_genfitHitDriftErr);
            sc=m_tuple->addItem("genfitTrackPos",3,m_genfitTrackPos);
            sc=m_tuple->addItem("genfitTrackMom",3,m_genfitTrackMom);
            sc=m_tuple->addItem("genfitTimeSeed",m_genfitTimeSeed);
            debug()<< "Book tuple RecGenfitAlgDC/genfit" << endmsg;
        }else{
            error()<< "Cannot book tuple RecGenfitAlgDC/genfit" << endmsg;
        }
    }//end of book tuple

    //init genfit event display
    //if(m_showDisplay) m_genfitDisplay = genfit::EventDisplay::getInstance();

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgDC::execute()
{
    edm4hep::ReconstructedParticleCollection* dcRecParticleCol=
        m_dcRecParticleCol.createAndPut();
    edm4hep::TrackCollection* dcRecParticleTrackCol=
        m_DCRecTrackCol.createAndPut();

    StatusCode sc;
    m_timer=clock();
    info()<<"RecGenfitAlgDC in execute()"<<endmsg;

    /////retrieve EventHeader
    //auto header = _headerCol.get()->at(0);
    //int evtNo = header.getEventNumber();
    //int runNo = header.getRunNumber();

    //std::cout<<"run "<<header.getEventNumber()
    //  <<" "<<header.getRunNumber()<<std::endl;

    ///retrieve Track and TrackHits
    const edm4hep::TrackCollection* dcTrackCol=nullptr;
    if(m_dcTrackCol.exist()) dcTrackCol=m_dcTrackCol.get();
    if(nullptr==dcTrackCol) {
        debug()<<"TrackCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }

    const edm4hep::TrackerHitCollection* didiDCHitsCol=nullptr;
    didiDCHitsCol=m_digiDCHitsCol.get();
    if(nullptr==didiDCHitsCol) {
        debug()<<"DigiDCHitCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    ///retrieve DC Hit Association
    auto assoDCHitsCol=m_dcHitAssociationCol.get();

    double eventStartTime=0;

    ///----------------------------------------------------
    ///Loop over Track and do fitting for each track
    ///----------------------------------------------------
    debug()<<"DCTrackCol size="<<dcTrackCol->size()<<endmsg;
    m_nDcTrack=dcTrackCol->size();
    for(auto dcTrack: *dcTrackCol){
        ///Loop over 5 particle hypothesis(0-4): e,mu,pi,K,p
        ///-1 for chargedgeantino
        for(unsigned int pidType=0;pidType<m_nPDG;pidType++){
            if((m_debugPid>=0) && (m_debugPid!=pidType)) continue;
            ///-----------------------------------
            ///Create a GenFit track
            ///-----------------------------------
            GenfitTrack* genfitTrack=new GenfitTrack(m_genfitField,
                    m_gridDriftChamber,m_geomSvc);
            genfitTrack->setDebug(m_debug);

            if(m_useMcParticleSeed){
                ///Retrieve MC particle(s)
                const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
                mcParticleCol=m_mcParticleCol.get();
                if(nullptr==mcParticleCol){
                    debug()<<"MCParticleCollection not found"<<endmsg;
                    return StatusCode::SUCCESS;
                }
                //for single track only, FIXME
                if(!genfitTrack->createGenfitTrackFromMCParticle(
                            pidType,*(mcParticleCol->begin()),
                            eventStartTime)){
                    debug()<<"createGenfitTrackFromMCParticle failed!"<<endmsg;
                    return StatusCode::SUCCESS;
                }
            }else{
                if(!genfitTrack->createGenfitTrackFromEDM4HepTrack(pidType,dcTrack,
                            eventStartTime,m_isUseCovTrack)){
                    debug()<<"createGenfitTrackFromEDM4HepTrack failed!"<<endmsg;
                    return StatusCode::SUCCESS;
                }
            }


            if(m_useTruthHit){
                bool fitSiliconOnly=false;
                bool isUseFixedSiHitError=true;
                if(0==genfitTrack->addHitsOnEdm4HepTrack(dcTrack,assoDCHitsCol,
                            m_sigmaHit,m_smearHit,fitSiliconOnly
                            , isUseFixedSiHitError)){
                    debug()<<"addSimTrackerHits failed!"<<endmsg;
                    return StatusCode::FAILURE;
                }
            }else{
                if(0==genfitTrack->addWireMeasurementOnTrack(dcTrack,
                            m_sigmaDrift.value(),m_smearHit,assoDCHitsCol)){
                    debug()<<"addWireMeasurementOnTrack failed!"<<endmsg;
                    return StatusCode::FAILURE;
                }
            }
            if(m_tuple) debugInitTrack(pidType,dcTrack,genfitTrack);
            if(m_debug){
                debug()<<"print seed track"<<endmsg;
                genfitTrack->printSeed();
                debug()<<"before processTrack"<<endmsg;
            }

            ///-----------------------------------
            ///call genfit fitting procedure
            ///-----------------------------------
            m_genfitFitter->processTrack(genfitTrack,m_resortHits.value());

            ///-----------------------------------
            ///Store track
            ///-----------------------------------
            auto dcRecParticle=dcRecParticleCol->create();
            auto dcRecTrack=dcRecParticleTrackCol->create();
            genfitTrack->storeTrack(dcRecParticle,dcRecTrack,pidType,m_ndfCut,
                    m_chi2Cut);
            if(m_debug){
                debug()<<"after storeTrack"<<endmsg;
                genfitTrack->printSeed();
            }

            if(m_tuple) debugTrack(pidType,genfitTrack);
            if(m_showDisplay) {
                //m_genfitDisplay->addEvent(genfitTrack->getTrack());
                //m_genfitDisplay->open();
            }else{
                delete genfitTrack;
            }
            ++m_fitSuccess[pidType];
        }//end loop over particle type
    }//end loop over a track

    if(m_tuple) debugEvent(dcTrackCol,dcRecParticleTrackCol,eventStartTime);

    //if(m_genfitDisplay) while(1){
    //    std::cout<<"Press any key to finish..."<<std::endl;
    //    //system ("pause");
    //}

    if(m_tuple) sc=m_tuple->write();

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgDC::finalize()
{
    MsgStream log(msgSvc(), name());
    info()<< "RecGenfitAlgDC in finalize()" << endmsg;

    m_genfitFitter->writeHist();
    delete m_genfitFitter;
    info()<<"RecGenfitAlgDC nRecTrack="<<m_nDCTrack<<" success e "
        <<m_fitSuccess[0]<<" mu "<<m_fitSuccess[1]<<" pi "<<m_fitSuccess[2]
        <<" K "<<m_fitSuccess[3]<<" p "<<m_fitSuccess[4]<<std::endl;
    if(m_nDCTrack>0){
        std::cout<<"RecGenfitAlgDC Success rate = "<<std::endl;
        for (int i=0;i<5;i++){
            std::cout<<Form("%d %2.2f",i,((float) m_fitSuccess[i])/m_nDCTrack)
                <<std::endl;
        }
    }
    return StatusCode::SUCCESS;
}

void RecGenfitAlgDC::debugInitTrack(int pidType,
        const edm4hep::Track& track,
        const GenfitTrack* genfitTrack)
{

    TVectorD seed=genfitTrack->getTrack()->getStateSeed();
    debug()<<"debug seed track"<<endmsg;
    if(m_debug) seed.Print();
    m_genfitTrackPos[0]=seed[0];
    m_genfitTrackPos[1]=seed[1];
    m_genfitTrackPos[2]=seed[2];
    m_genfitTrackMom[0]=seed[3];
    m_genfitTrackMom[1]=seed[4];
    m_genfitTrackMom[2]=seed[5];
    m_genfitTimeSeed=genfitTrack->getTrack()->getTimeSeed();
    debug()<<"NumPoints "<<genfitTrack->getTrack()->getNumPoints()
        <<" NumPointsWithMeasurement "
        <<genfitTrack->getTrack()->getNumPointsWithMeasurement()<<endmsg;
    m_genfitTrackNumPoint=genfitTrack->getTrack()->getNumPoints();
    m_genfitTrackNumPointsWithMeas=genfitTrack->getTrack()
        ->getNumPointsWithMeasurement();
    int iMeas=0;
    for(auto trackPoint:genfitTrack->getTrack()->getPoints()){
        int nMeas=trackPoint->getNumRawMeasurements();
        debug()<<"nMeas "<<nMeas<<endmsg;
        for(int i=0;i<nMeas;i++){
            genfit::AbsMeasurement* rawMeas=trackPoint->getRawMeasurement(i);
            //drift chamber FIXME
            debug()<<"nDim "<<rawMeas->getDim()<<endmsg;
            if(m_debug>0) rawMeas->getRawHitCoords().Print();
            if(m_debug>0) rawMeas->getRawHitCov().Print();
            if(1==trackPoint->getRawMeasurement(i)->getDim()){
                debug()<<" detId "<<rawMeas->getDetId()<<endmsg;
                m_genfitHitLayer[iMeas]=
                    m_decoder->get(rawMeas->getDetId(),"layer");
                m_genfitHitCell[iMeas]=
                    m_decoder->get(rawMeas->getDetId(),"cellID");
                m_genfitHitDrift[iMeas]=rawMeas->getRawHitCoords()[0];
                m_genfitHitDriftErr[iMeas]=rawMeas->getRawHitCov()(0,0);
                if(m_debug>0){
                    std::cout<<"layer "
                        <<m_decoder->get(rawMeas->getDetId(),"layer")<<" cellID "
                        <<m_decoder->get(rawMeas->getDetId(),"cellID")
                        <<" m_genfitHitDriftErr "<<rawMeas->getRawHitCov()(0,0)
                        <<std::endl;
                }
            }
            iMeas++;
        }
    }
    m_genfitNHit=iMeas;
}

void RecGenfitAlgDC::debugTrack(int pidType,const GenfitTrack* genfitTrack)
{
    /// Get fit status
    const genfit::FitStatus* fitState = genfitTrack->getFitStatus();
    int charge= fitState->getCharge();

    m_chargeKal[pidType]= charge;
    m_nHitWithFitInfo[pidType]=genfitTrack->getNumPointsWithFittedInfo();
    m_chi2Kal[pidType]=fitState->getChi2();
    m_nDofKal[pidType]=fitState->getNdf();
    m_isFitted[pidType]=(int)fitState->isFitted();
    m_isFitConverged[pidType]=(int) fitState->isFitConverged();
    m_isFitConvergedFully[pidType]=(int) fitState->isFitConvergedFully();

    ///get fitted state of track
    TMatrixDSym fittedCov;
    TLorentzVector fittedPos;
    TVector3 fittedMom;
    int fittedState=genfitTrack->getFittedState(fittedPos,fittedMom,fittedCov);

    m_fitPosx = fittedPos.X();
    m_fitPosy = fittedPos.Y();
    m_fitPosz = fittedPos.Z();
    m_fitMomx = fittedMom.X();
    m_fitMomy = fittedMom.Y();
    m_fitMomz = fittedMom.Z();

    TVector3 fittedPos_ext2Origin(1e9,1e9,1e9);
    TVector3 fittedMom_ext2Origin(1e9,1e9,1e9);
    TMatrixDSym fittedCov_ext2Origin;
    const TVector3 referencePoint(0.,0.,0.);
    double trackLength = genfitTrack->extrapolateToPoint(fittedPos_ext2Origin,
            fittedMom_ext2Origin,fittedCov_ext2Origin,referencePoint);

    if(m_debug>0){
        debug()<<"Fitted pos and mom "<<endmsg;
        fittedPos.Print();
        fittedMom.Print();
        debug()<<"After extrapolate to origin trackLength "<<trackLength<<endmsg;
        fittedPos_ext2Origin.Print();
        fittedMom_ext2Origin.Print();
    }

    m_extraPos[0]=fittedPos_ext2Origin.X();
    m_extraPos[1]=fittedPos_ext2Origin.Y();
    m_extraPos[2]=fittedPos_ext2Origin.Z();
    m_extraMom[0]=fittedMom_ext2Origin.X();
    m_extraMom[1]=fittedMom_ext2Origin.Y();
    m_extraMom[2]=fittedMom_ext2Origin.Z();

    m_pocaMomKalP[pidType]=fittedMom_ext2Origin.Mag();
    m_pocaMomKal[pidType][0]=fittedMom_ext2Origin.X();
    m_pocaMomKal[pidType][1]=fittedMom_ext2Origin.Y();
    m_pocaMomKal[pidType][2]=fittedMom_ext2Origin.Z();
    m_pocaPosKal[pidType][0]=fittedPos_ext2Origin.X();
    m_pocaPosKal[pidType][1]=fittedPos_ext2Origin.Y();
    m_pocaPosKal[pidType][2]=fittedPos_ext2Origin.Z();

    for(int i=0; i<6;i++) {
        m_Error6[i] = fittedCov[i][i];
    }

    HelixClass helix;//mm and GeV
    float pos[3]={float(fittedPos.X()/dd4hep::mm),float(fittedPos.Y()/dd4hep::mm),
        float(fittedPos.Z()/dd4hep::mm)};
    float mom[3]={float(fittedMom.X()),float(fittedMom.Y()),float(fittedMom.Z())};
    helix.Initialize_VP(pos,mom,charge,m_genfitField->getBzTesla(fittedPos.Vect()));
    m_pocaMomKalPFirst[pidType]=fittedMom.Mag();
    m_pocaMomKalPt[pidType]=fittedMom.Perp();

    if(m_debug>0){
        /// Get fit status
        debug()<<"evt "<<m_evt<<" fit result: get status OK? pidType "
            <<pidType<<" fittedState "<<fittedState<<" isFitted "
            <<m_isFitted[pidType]<<" isConverged "<<m_isFitConverged[pidType]
            <<" isFitConvergedFully "<<m_isFitConvergedFully[pidType]
            <<" ndf "<<m_nDofKal[pidType]
            <<" chi2 "<<m_chi2Kal[pidType]<<endmsg;
        if((0!=fittedState)||(!m_isFitted[pidType])||(m_nDofKal[pidType]>m_ndfCut)){
            debug()<<"fitting failed"<<endmsg;
        }else{
            debug()<<"evt "<<m_evt<<" fit result: Pos("<<
                fittedPos.X()<<" "<<
                fittedPos.Y()<<" "<<
                fittedPos.Z()<<") mom("<<
                fittedMom.X()<<" "<<
                fittedMom.Y()<<" "<<
                fittedMom.Z()<<") p_tot "<<
                fittedMom.Mag()<<" pt "<<
                fittedMom.Perp()<<" ext pos "<<
                fittedPos_ext2Origin.X()<<" "<<
                fittedPos_ext2Origin.Y()<<" "<<
                fittedPos_ext2Origin.Z()<<") ext mom("<<
                fittedMom_ext2Origin.X()<<" "<<
                fittedMom_ext2Origin.Y()<<" "<<
                fittedMom_ext2Origin.Z()<<") ext p_tot "<<
                fittedMom_ext2Origin.Mag()<<" ext pt "<<
                fittedMom_ext2Origin.Perp()
                <<endmsg;
        }
    }

    //debug extrapolation
    ///Get track parameters
    TVectorD seed=genfitTrack->getTrack()->getStateSeed();
    TLorentzVector posInit;
    TVector3 pos_t(seed[0],seed[1],seed[2]);
    TVector3 mom_t(seed[3],seed[4],seed[5]);
    //double charge(0);
    int iDCDigi=0;
    const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
    dCDigiCol=m_digiDCHitsCol.get();
    for(auto dcDigi: *dCDigiCol){
        TVector3 poca,pocaDir,pocaOnWire;
        double doca;
        TVector3 endPointStart(0,0,0);
        TVector3 endPointEnd(0,0,0);
        m_gridDriftChamber->cellposition(dcDigi.getCellID(),endPointStart,
                endPointEnd);//cm
        endPointStart.SetX(endPointStart.x()*dd4hep::cm);
        endPointStart.SetY(endPointStart.y()*dd4hep::cm);
        endPointStart.SetZ(endPointStart.z()*dd4hep::cm);
        endPointEnd.SetX(endPointEnd.x()*dd4hep::cm);
        endPointEnd.SetY(endPointEnd.y()*dd4hep::cm);
        endPointEnd.SetZ(endPointEnd.z()*dd4hep::cm);

        bool stopAtBoundary=false;
        bool calcJacobianNoise=false;
        genfitTrack->extrapolateToHit(poca,pocaDir,pocaOnWire,doca,pos_t,mom_t,
                endPointStart,endPointEnd,m_debug.value(),0,stopAtBoundary,
                calcJacobianNoise);
        if(m_debug){
            std::cout<<" endPointStart ";
            endPointStart.Print();
            std::cout<<" endPointEnd";
            endPointEnd.Print();
            std::cout<<" ext origin pos "; pos_t.Print();
            std::cout<<" mom "; mom_t.Print();
            std::cout<<" poca ";
            poca.Print();
            std::cout<<" pocaDir ";
            pocaDir.Print();
            std::cout<<" pocaOnWire ";
            pocaOnWire.Print();
            std::cout<<" doca "<<doca<<std::endl;
            std::cout<<" pos "<<std::endl;
            pos_t.Print();
        }
        m_dcDigiDocaExt[iDCDigi]=doca;
        iDCDigi++;
    }
}

void RecGenfitAlgDC::debugEvent(const edm4hep::TrackCollection* sdtTrackCol,
        const edm4hep::TrackCollection* sdtRecTrackCol,
        double eventStartTime)
{
    int iSdtTrack=0;
    m_nSdtTrack=sdtTrackCol->size();
    for(auto sdtTrack: *sdtTrackCol){
        if(iSdtTrack>0) break;//TODO debug for single track only
        edm4hep::TrackState trackStat=sdtTrack.getTrackStates(0);//FIXME?
        HelixClass helixClass;//mm and GeV
        helixClass.Initialize_Canonical(trackStat.phi,trackStat.D0,
                trackStat.Z0,trackStat.omega,trackStat.tanLambda,
                m_genfitField->getBzTesla({0.,0.,0.}));

        TLorentzVector posInit(helixClass.getReferencePoint()[0],
                helixClass.getReferencePoint()[1],
                helixClass.getReferencePoint()[2],eventStartTime);
        m_seedPos[0]=posInit.X();
        m_seedPos[1]=posInit.Y();
        m_seedPos[2]=posInit.Z();
        TVector3 momInit(helixClass.getMomentum()[0],
                helixClass.getMomentum()[1],helixClass.getMomentum()[2]);
        m_seedMomP=momInit.Mag();
        m_seedMomPt=momInit.Perp();
        m_seedMom[0]=momInit.X();
        m_seedMom[1]=momInit.Y();
        m_seedMom[2]=momInit.Z();
        iSdtTrack++;
        TVector3 pos,mom;
        TMatrixDSym cov(6);
        double charge;
        CEPC::getPosMomFromTrackState(trackStat,m_genfitField->getBzTesla({0.,0.,0.}),
                pos,mom,charge,cov);
        m_seedMomQ=charge;
        debug()<<" sdtTrack charge "<<charge<<" seed mom "<<momInit.X()<<" "<<
            momInit.Y()<<" "<<momInit.Z()<<endmsg;
        debug()<<" sdtTrack "<<charge<<" seed pos "<<posInit.X()<<" "<<
            posInit.Y()<<" "<<posInit.Z()<<endmsg;
        if(m_debug>0){
            debug()<<"pos"<<endmsg;
            pos.Print();
            debug()<<"mom"<<endmsg;
            mom.Print();
            debug()<<"cov"<<endmsg;
            cov.Print();
        }
    }

    const edm4hep::MCParticleCollection* mcParticleCol = nullptr;
    const edm4hep::SimTrackerHitCollection* simDCHitCol=nullptr;

    m_pidIndex=5;

    mcParticleCol=m_mcParticleCol.get();
    int iMcParticle=0;
    HelixClass helix_mcP;
    for(auto mcParticle : *mcParticleCol){
        edm4hep::Vector3f mcPocaMom = mcParticle.getMomentum();//GeV
        edm4hep::Vector3d mcPocaPos = mcParticle.getVertex();
        double mcPos[3]={(mcPocaPos.x)*dd4hep::mm,(mcPocaPos.y)*dd4hep::mm,(mcPocaPos.z)*dd4hep::mm};
        double mcMom[3]={(mcPocaMom.x),(mcPocaMom.y),(mcPocaMom.z)};
        for(int i=0;i<3;i++){debug()<<"mcPos "<<mcPos[i]<<endmsg;}
        for(int i=0;i<3;i++){debug()<<"mcMom "<<mcMom[i]<<endmsg;}
        float mcCharge = mcParticle.getCharge();
        helix_mcP.Initialize_VP(mcPos,mcMom,mcCharge,
                m_genfitField->getBzTesla(mcPos));
        mcP_D0 = helix_mcP.getD0();
        mcP_phi = helix_mcP.getPhi0();
        mcP_omega = helix_mcP.getOmega();
        mcP_Z0 = helix_mcP.getZ0();
        mcP_tanLambda = helix_mcP.getTanLambda();
        debug()<< " debugEvent Bz " << m_genfitField->getBzTesla(mcPos)
            << " mc d0= " << mcP_D0
            << " phi0= " << mcP_phi
            << " omega= " << mcP_omega
            << " Z0= " << mcP_Z0
            << " tanLambda= " << mcP_tanLambda << endmsg;
        m_seedCenterX=helix_mcP.getXC();
        m_seedCenterY=helix_mcP.getYC();
        m_seedR=helix_mcP.getRadius();

        float px=mcPocaMom.x;
        float py=mcPocaMom.y;
        float pz=mcPocaMom.z;
        debug()<<"mc pxyz   "<<px<<" "<<py<<" "<<pz<<endmsg;
        debug()<<"mc circle cx "<<helix_mcP.getXC()<<" "<<helix_mcP.getYC()<<" "<<helix_mcP.getRadius()<<endmsg;
        m_pocaMomMcP[iMcParticle]=sqrt(px*px+py*py+pz*pz);
        m_pocaMomMcPt[iMcParticle]=sqrt(px*px+py*py);
        m_pocaMomMc[iMcParticle][0]=px;
        m_pocaMomMc[iMcParticle][1]=py;
        m_pocaMomMc[iMcParticle][2]=pz;
        m_pocaPosMc[iMcParticle][0]=mcPocaPos.x;
        m_pocaPosMc[iMcParticle][1]=mcPocaPos.y;
        m_pocaPosMc[iMcParticle][2]=mcPocaPos.z;
        iMcParticle++;
    }
    m_mcIndex=iMcParticle;

    int iHit=0;
    simDCHitCol=m_simDCHitCol.get();
    for(auto simDCHit: *simDCHitCol){
        edm4hep::Vector3d pos=simDCHit.position();
        TVectorD p(3);
        p[0]=pos.x;//no unit conversion here
        p[1]=pos.y;
        p[2]=pos.z;
        m_dcHitMcX[iHit]=pos.x;
        m_dcHitMcY[iHit]=pos.y;
        m_dcHitMcZ[iHit]=pos.z;
        m_dcHitMcDoca[iHit]=simDCHit.getTime()*40.*dd4hep::um*dd4hep::mm;
        TVector3 endPointStart(0,0,0);
        TVector3 endPointEnd(0,0,0);
        m_gridDriftChamber->cellposition(simDCHit.getCellID(),endPointStart,
                endPointEnd);
        m_dcHitMcWireX[iHit]=endPointStart.X();
        m_dcHitMcWireY[iHit]=endPointStart.Y();
        iHit++;
    }
    m_nSimDCHit=simDCHitCol->size();
    const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
    dCDigiCol=m_digiDCHitsCol.get();
    if(nullptr!=dCDigiCol){ m_nDCDigi=dCDigiCol->size(); }
    int iDCDigi=0;
    for(auto dcDigi: *dCDigiCol){
        m_dcDigiChamber[iDCDigi]=m_decoder->get(dcDigi.getCellID(),"chamber");
        m_dcDigiLayer[iDCDigi]=m_decoder->get(dcDigi.getCellID(),"layer");
        m_dcDigiCell[iDCDigi]=m_decoder->get(dcDigi.getCellID(),"cellID");
        m_dcDigiTime[iDCDigi]=dcDigi.getTime();
        m_dcDigiDoca[iDCDigi]=dcDigi.getTime()*m_driftVelocity/10000.; //cm
        TVector3 endPointStart(0,0,0);
        TVector3 endPointEnd(0,0,0);
        m_gridDriftChamber->cellposition(dcDigi.getCellID(),endPointStart,
                endPointEnd);
        m_dcDigiWireStartX[iDCDigi]=endPointStart.X();
        m_dcDigiWireStartY[iDCDigi]=endPointStart.Y();
        m_dcDigiWireStartZ[iDCDigi]=endPointStart.Y();
        m_dcDigiWireEndX[iDCDigi]=endPointEnd.X();
        m_dcDigiWireEndY[iDCDigi]=endPointEnd.Y();
        m_dcDigiWireEndZ[iDCDigi]=endPointEnd.Y();

        //get information from associated simTrackerHit
        const edm4hep::ConstSimTrackerHit dcSimTrackerHit=
            m_dcHitAssociationCol.get()->at(iDCDigi).getSim();
        m_dcDigiMcMomX[iDCDigi]=dcSimTrackerHit.getMomentum().x*dd4hep::mm;
        m_dcDigiMcMomY[iDCDigi]=dcSimTrackerHit.getMomentum().y*dd4hep::mm;
        m_dcDigiMcMomZ[iDCDigi]=dcSimTrackerHit.getMomentum().z*dd4hep::mm;
        m_dcDigiMcPosX[iDCDigi]=dcSimTrackerHit.getPosition().x*dd4hep::mm;
        m_dcDigiMcPosY[iDCDigi]=dcSimTrackerHit.getPosition().y*dd4hep::mm;
        m_dcDigiMcPosZ[iDCDigi]=dcSimTrackerHit.getPosition().z*dd4hep::mm;

        iDCDigi++;
    }

    m_nSdtTrack=sdtRecTrackCol->size();
    for(auto sdtTrack: *sdtRecTrackCol){
        for(unsigned int i=0; i<sdtTrack.trackStates_size(); i++) {
            edm4hep::TrackState trackStat=sdtTrack.getTrackStates(i);
            std::array<float,15> errorCov;
            errorCov = trackStat.covMatrix;
            for(int j=0; j<15; j++) {
                m_ErrorcovMatrix[j] = errorCov[j];
                debug()<<"debugEvent2 errorCov "<<j<<" "<<errorCov[j]<<endmsg;
            }
            m_D0 = trackStat.D0;
            m_phi = trackStat.phi;
            m_omega = trackStat.omega;
            m_Z0 = trackStat.Z0;
            m_tanLambda = trackStat.tanLambda;
        }
    }
}

/*void RecGenfitAlgDC::debugEvent()
  {
  const edm4hep::MCParticleCollection* mcParticleCol = nullptr;
  const edm4hep::SimTrackerHitCollection* simDCHitCol=nullptr;

  m_pidIndex=5;

  mcParticleCol=m_mcParticleCol.get();
  simDCHitCol=m_simDCHitCol.get();
  m_nSimDCHit=simDCHitCol->size();
  int iMcParticle=0;
  int iHit=0;
  for(auto mcParticle : *mcParticleCol){
  for(auto simDCHit: *simDCHitCol){
  edm4hep::Vector3d pos=simDCHit.position();
  TVectorD p(3);
  p[0]=pos.x;//no unit conversion here
  p[1]=pos.y;
  p[2]=pos.z;
  m_dcHitMcX[iHit]=pos.x;
  m_dcHitMcY[iHit]=pos.y;
  m_dcHitMcZ[iHit]=pos.z;
  iHit++;
  }
  edm4hep::Vector3f mcPocaMom = mcParticle.getMomentum();//GeV
  float px=mcPocaMom.x;
  float py=mcPocaMom.y;
  float pz=mcPocaMom.z;
  debug()<<"   "<<px<<" "<<py<<" "<<pz<<endmsg;
  m_pocaMomMcP[iMcParticle]=sqrt(px*px+py*py+pz*pz);
  iMcParticle++;
  }
  m_mcIndex=iHit;

  }*/
