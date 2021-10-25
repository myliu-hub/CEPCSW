#include "RecGenfitAlgDC.h"
#include "GenfitTrack.h"
#include "GenfitFitter.h"
#include "GenfitField.h"

//genfit
#include "EventDisplay.h"
#include "AbsMeasurement.h"
#include "WireMeasurementNew.h"
#include "MaterialEffects.h"

//cepcsw
#include "DetInterface/IGeomSvc.h"
#include "DataHelper/HelixClass.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "DataHelper/TrackHelper.h"
#include "DataHelper/TrackerHitHelper.h"
#include "UTIL/ILDConf.h"
#include "GenfitMaterialInterface.h"

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
#include "TDatabasePDG.h"
#include "TH1D.h"
#include "LSFitting.h"


//stl
#include <chrono>
#include "time.h"
//#undef GENFIT_MY_DEBUG
#define GENFIT_MY_DEBUG 1

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
  declareProperty("DCHitAssociationCollection", m_DCHitAssociationCol,
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
  m_genfitFitter=new GenfitFitter(m_fitterType.toString().c_str(),m_debug.value());
  m_genfitField=new GenfitField(m_dd4hepField);
  m_genfitFitter->setField(m_genfitField);
  m_genfitFitter->setGeoMaterial(m_geomSvc->lcdd(),m_extMinDistCut.value(),
      m_skipWireMaterial.value());
  m_genfitFitter->setEnergyLossBrems(m_correctBremsstrahlung.value());
  m_genfitFitter->setNoiseBrems(m_correctBremsstrahlung.value());
  m_genfitFitter->setDebug(m_debug.value());
  m_genfitFitter->setDebugGenfit(m_debugGenfit.value());
  m_genfitFitter->setDebugLocal((unsigned int) m_debugLocal.value());
  m_genfitFitter->setBlowUpFactor(m_blowUpFactor);
  m_genfitFitter->setResetOffDiagonals(m_resetOffDiagonals);
  m_genfitFitter->setBlowUpMaxVal(m_blowUpMaxVal);
  m_genfitFitter->setMaterialDebugLvl(m_debugMaterial);
  if(m_noMaterialEffects.value()) m_genfitFitter->setNoEffects(true);
  if(-1==m_debugPid.value()) m_genfitFitter->setNoEffects(true);
  if(-1==m_debugPid.value()) m_debugPid=0;//charged geantino with electron pid
  if(m_fitterType=="DAF"||m_fitterType=="DafRef"){
    m_genfitFitter->setMaxIterationsBetas(m_bStart.value(),m_bFinal.value(),m_maxIteration.value());
#ifdef GENFIT_MY_DEBUG
    m_genfitFitter->setDebugLocal(m_debugLocal.value());
#endif
  } else {
    m_genfitFitter->setMaxIterations(m_maxIteration.value());
  }
  //print genfit parameters
  if(m_debug.value()>0) m_genfitFitter->print();
  if(""!=m_genfitHistRootName) m_genfitFitter->initHist(m_genfitHistRootName.value());

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
      sc=m_tuple->addItem("momLsP",m_momLsP);

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
      sc=m_tuple->addItem("fitMom",m_fitMom);
      sc=m_tuple->addItem("firstMomMc",m_firstMomMc);

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
      sc=m_tuple->addItem("dcDigiPocaExtX",m_nDCDigi,m_dcDigiPocaExtX);
      sc=m_tuple->addItem("dcDigiPocaExtY",m_nDCDigi,m_dcDigiPocaExtY);
      sc=m_tuple->addItem("dcDigiPocaExtZ",m_nDCDigi,m_dcDigiPocaExtZ);
      sc=m_tuple->addItem("dcDigiDocaMC",m_nDCDigi,m_dcDigiDocaMC);
      sc=m_tuple->addItem("dcDigiDocaExt",m_nDCDigi,m_dcDigiDocaExt);
      sc=m_tuple->addItem("dcDigiDocaIdeal",m_nDCDigi,m_dcDigiDocaIdeal);
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

      sc=m_tuple->addItem("nSimDCHit",m_nSimDCHit,0,500000);
      sc=m_tuple->addItem("dcHitDriftT",m_nSimDCHit,m_dcHitDriftT);
      sc=m_tuple->addItem("dcHitDriftDl",m_nSimDCHit,m_dcHitDriftDl);
      sc=m_tuple->addItem("dcHitDriftDr",m_nSimDCHit,m_dcHitDriftDr);
      sc=m_tuple->addItem("dcHitLr",m_nSimDCHit,m_dcHitLr);
      sc=m_tuple->addItem("dcHitLayer",m_nSimDCHit,m_dcHitLayer);
      sc=m_tuple->addItem("dcHitWire",m_nSimDCHit,m_dcHitWire);
      sc=m_tuple->addItem("dcHitErr",m_nSimDCHit,m_dcHitErr);
      sc=m_tuple->addItem("time",5,m_time);
      sc=m_tuple->addItem("dcHitMcTkId",m_nSimDCHit,m_dcHitMcTkId);
      sc=m_tuple->addItem("dcHitMcLr",m_nSimDCHit,m_dcHitMcLr);
      sc=m_tuple->addItem("dcHitMcDrift",m_nSimDCHit,m_dcHitMcDrift);
      sc=m_tuple->addItem("dcHitMcX",m_nSimDCHit,m_dcHitMcX);
      sc=m_tuple->addItem("dcHitMcY",m_nSimDCHit,m_dcHitMcY);
      sc=m_tuple->addItem("dcHitMcZ",m_nSimDCHit,m_dcHitMcZ);
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
      sc=m_tuple->addItem("genfitHitEndX",m_genfitNHit,m_genfitHitEndX);
      sc=m_tuple->addItem("genfitHitEndY",m_genfitNHit,m_genfitHitEndY);
      sc=m_tuple->addItem("genfitHitEndZ",m_genfitNHit,m_genfitHitEndZ);
      sc=m_tuple->addItem("genfitHitDrift",m_genfitNHit,m_genfitHitDrift);
      sc=m_tuple->addItem("genfitHitDriftErr",m_genfitNHit,m_genfitHitDriftErr);
      sc=m_tuple->addItem("genfitTrackPos",3,m_genfitTrackPos);
      sc=m_tuple->addItem("genfitTrackMom",3,m_genfitTrackMom);
      sc=m_tuple->addItem("genfitTimeSeed",m_genfitTimeSeed);
      sc=m_tuple->addItem("momLsMCP",m_momLsMCP);
      sc=m_tuple->addItem("momLsMCPFirst",m_momLsMCPFirst);
      sc=m_tuple->addItem("momLsP",m_momLsP);
      sc=m_tuple->addItem("momLsPt",m_momLsPt);
      sc=m_tuple->addItem("momLsPz",m_momLsPz);
      debug()<< "Book tuple RecGenfitAlgDC/genfit" << endmsg;
    }else{
      error()<< "Cannot book tuple RecGenfitAlgDC/genfit" << endmsg;
    }
  }//end of book tuple

  //init genfit event display
  //if(m_showDisplay.value()) m_genfitDisplay = genfit::EventDisplay::getInstance();

  m_iEvent=-1;
  //if(m_debug) {
  //    m_genfitFitter->getMaterialEffects()->drawdEdx(11);
  //    m_genfitFitter->getMaterialEffects()->drawdEdx(13);
  //    m_genfitFitter->getMaterialEffects()->drawdEdx(211);
  //    m_genfitFitter->getMaterialEffects()->drawdEdx(321);
  //    m_genfitFitter->getMaterialEffects()->drawdEdx(2212);
  //}

  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode RecGenfitAlgDC::execute()
{
  m_iEvent++;
  edm4hep::ReconstructedParticleCollection* dcRecParticleCol=
    m_dcRecParticleCol.createAndPut();
  edm4hep::TrackCollection* dcRecParticleTrackCol=
    m_DCRecTrackCol.createAndPut();

  if(m_iEventSelect>0&&m_iEvent!=m_iEventSelect.value()){
    std::cout<<"Skip event "<<m_iEvent<<std::endl;
    return StatusCode::SUCCESS;
  }
  StatusCode sc;
  m_timer=clock();
  std::cout<<"RecGenfitAlgDC in execute() event No."<<m_iEvent<<std::endl;
  //info()<<"RecGenfitAlgDC in execute() event "<<m_iEvent<<endmsg;

  if(m_debugLsFit){
    lsFit(m_smearTrack,m_smearHit,m_useFirstHitAsSeed);
    if(m_tuple) sc=m_tuple->write();
    return StatusCode::SUCCESS;
  }
  /////retrieve EventHeader
  //auto header = _headerCol.get()->at(0);
  //int evtNo = header.getEventNumber();
  //int runNo = header.getRunNumber();

  //std::cout<<"run "<<header.getEventNumber()
  //  <<" "<<header.getRunNumber()<<std::endl;

  ///retrieve Track and TrackHits
  auto dcTrackCol=m_dcTrackCol.get();
  if(nullptr==dcTrackCol) {
    debug()<<"TrackCollection not found"<<endmsg;
    return StatusCode::SUCCESS;
  }

  ///retrieve Digi
  auto didiDCHitsCol=m_digiDCHitsCol.get();
  if(nullptr==didiDCHitsCol) {
    debug()<<"DigiDCHitCollection not found"<<endmsg;
    return StatusCode::SUCCESS;
  }
  ///retrieve DC Hit Association
  auto assoDCHitsCol=m_DCHitAssociationCol.get();

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
      if((m_debugPid.value()>=0) && (m_debugPid.value()!=(int) pidType)) continue;
      ///-----------------------------------
      ///Create a GenFit track
      ///-----------------------------------
      GenfitTrack* genfitTrack=new GenfitTrack(m_genfitField,
          m_gridDriftChamber,m_geomSvc);
      genfitTrack->setDebug(m_debug.value());

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


      if(m_useTruthHit.value()){
        if(0==genfitTrack->addSpacePointsDC(dcTrack,assoDCHitsCol,m_sigmaHitU,m_sigmaHitV)){
          debug()<<"addSimTrackerHits failed!"<<endmsg;
          return StatusCode::FAILURE;
        }
      }else{
        if(0==genfitTrack->addWireMeasurements(dcTrack,m_sigmaDrift,
                    assoDCHitsCol,m_sortMethod,m_truthAmbig,
                    m_skipCorner,m_skipNear)){
          debug()<<"no hits in addWireMeasurements!"<<endmsg;
          return StatusCode::SUCCESS;
        }
      }
      if(m_tuple||m_debug.value()>0) debugInitTrack(genfitTrack);

#ifdef GENFIT_MY_DEBUG
      m_genfitFitter->SetRunEvent(m_iEvent);
#endif
      m_genfitFitter->setDebugGenfit(m_debugGenfit.value());
      m_genfitFitter->setDebugLocal(m_debugLocal.value());
      ///-----------------------------------
      ///call genfit fitting procedure
      ///-----------------------------------
      m_genfitFitter->processTrack(genfitTrack,m_resortHits.value());

      ///-----------------------------------
      ///Store track
      ///-----------------------------------
      auto dcRecParticle=dcRecParticleCol->create();
      auto dcRecTrack=dcRecParticleTrackCol->create();
      genfitTrack->storeTrack(dcRecParticle,dcRecTrack,pidType,m_ndfCut.value(),
          m_chi2Cut.value());
      if(m_debug.value()>0){
        debug()<<"print seed after storeTrack"<<endmsg;
        genfitTrack->printSeed();
      }

      if(m_debug.value()>0 || nullptr!=m_tuple){
        debugTrack(pidType,genfitTrack,dcTrack);
      }
      if(m_showDisplay) {
        //m_genfitDisplay->addEvent(genfitTrack->getTrack());
        //m_genfitDisplay->open();
        //TODO
      }
      delete genfitTrack;
      ++m_fitSuccess[pidType];

    }//end loop over particle type

  }//end loop over a track

  if(m_tuple||m_debug) debugEvent(dcTrackCol,dcRecParticleTrackCol,eventStartTime);

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

#ifdef GENFIT_MY_DEBUG
  m_genfitFitter->writeHist();
#endif
  delete m_genfitFitter;
  delete m_genfitField;
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

void RecGenfitAlgDC::debugInitTrack(const GenfitTrack* genfitTrack)
{

  if(m_debug.value()>0)genfitTrack->printSeed();
  if(m_tuple){
    TVectorD seed=genfitTrack->getTrack()->getStateSeed();
    m_genfitTrackPos[0]=seed[0];
    m_genfitTrackPos[1]=seed[1];
    m_genfitTrackPos[2]=seed[2];
    m_genfitTrackMom[0]=seed[3];
    m_genfitTrackMom[1]=seed[4];
    m_genfitTrackMom[2]=seed[5];
    m_genfitTimeSeed=genfitTrack->getTrack()->getTimeSeed();
    m_genfitTrackNumPoint=genfitTrack->getTrack()->getNumPoints();
    m_genfitTrackNumPointsWithMeas=genfitTrack->getTrack()
      ->getNumPointsWithMeasurement();
    int iMeas=0;
    for(auto trackPoint:genfitTrack->getTrack()->getPoints()){
      int nMeas=trackPoint->getNumRawMeasurements();
      //debug()<<"nMeas "<<nMeas<<endmsg;
      for(int i=0;i<nMeas;i++){
        //genfit::AbsMeasurement* rawMeas=trackPoint->getRawMeasurement(i);
        genfit::WireMeasurementNew* rawMeas=
          dynamic_cast<genfit::WireMeasurementNew*> (trackPoint->getRawMeasurement(i));
        //drift chamber FIXME
        //debug()<<"nDim "<<rawMeas->getDim()<<endmsg;
        //if(m_debug>0) rawMeas->getRawHitCoords().Print();
        //if(m_debug>0) rawMeas->getRawHitCov().Print();
        if(1==trackPoint->getRawMeasurement(i)->getDim()){
          //debug()<<" detId "<<rawMeas->getDetId()<<endmsg;
          m_genfitHitLayer[iMeas]=
            m_decoder->get(rawMeas->getDetId(),"layer");
          m_genfitHitCell[iMeas]=
            m_decoder->get(rawMeas->getDetId(),"cellID");

          TVector3 endPointStart(0,0,0);
          TVector3 endPointEnd(0,0,0);
          m_gridDriftChamber->cellposition(rawMeas->getDetId(),endPointStart,
              endPointEnd);//cm
          m_genfitHitEndX[iMeas]=endPointStart.X();
          m_genfitHitEndY[iMeas]=endPointStart.X();
          m_genfitHitEndZ[iMeas]=endPointStart.X();
          m_genfitHitDrift[iMeas]=rawMeas->getRawHitCoords()[0];
          m_genfitHitDriftErr[iMeas]=rawMeas->getRawHitCov()(0,0);
          //if(m_debug>0){
          //    std::cout<<"layer "
          //        <<m_decoder->get(rawMeas->getDetId(),"layer")<<" cellID "
          //        <<m_decoder->get(rawMeas->getDetId(),"cellID")
          //        <<" m_genfitHitDriftErr "<<rawMeas->getRawHitCov()(0,0)
          //        <<" m_genfitHitEndX "<<rawMeas->getWireEndPoints(0).X()
          //        <<" m_genfitHitEndY "<<rawMeas->getWireEndPoints(0).Y()
          //        <<" m_genfitHitEndZ "<<rawMeas->getWireEndPoints(0).Z()
          //        <<std::endl;
          //}
        }
        iMeas++;
      }
    }
    m_genfitNHit=iMeas;
  }
}

void RecGenfitAlgDC::debugTrack(int pidType,const GenfitTrack* genfitTrack,
    const edm4hep::Track dcTrack)
{

  debug()<<"=========================="<<endmsg;
  debug()<<"debug track after fitting"<<endmsg;
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
  m_fitMom = fittedMom.Mag();

  TVector3 fittedPos_ext2Origin(1e9,1e9,1e9);
  TVector3 fittedMom_ext2Origin(1e9,1e9,1e9);
  TMatrixDSym fittedCov_ext2Origin;
  const TVector3 referencePoint(0.,0.,0.);
  double trackLength = genfitTrack->extrapolateToPoint(fittedPos_ext2Origin,
      fittedMom_ext2Origin,fittedCov_ext2Origin,referencePoint);

  if(m_debug.value()>0){
    debug()<<"Fitted pos and mom "<<fittedPos.X()<<" "<<fittedPos.Y()
      <<" "<<fittedPos.Z()<<") mom"<<fittedMom.X()<<" "<<fittedMom.Y()
      <<" "<<fittedMom.Z() <<endmsg;
    debug()<<"After extrapolate to origin trackLength "<<trackLength<<endmsg;
    debug()<<fittedPos_ext2Origin.X()<<" "<<fittedPos_ext2Origin.Y()
      <<" "<<fittedPos_ext2Origin.Z()<<") mom"<<fittedMom_ext2Origin.X()
      <<" "<<fittedMom_ext2Origin.Y()
      <<" "<<fittedMom_ext2Origin.Z() <<endmsg;
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

  if(m_debug.value()>0){
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

  const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
  mcParticleCol=m_mcParticleCol.get();
  double xc,yc,r;
  double Bz=3.0;
  double pos_t[3];
  double mom_t[3];
  for(auto mcParticle : *mcParticleCol){
    const edm4hep::Vector3d mcParticleVertex=mcParticle.getVertex();//mm
    const edm4hep::Vector3f mcParticleMom=mcParticle.getMomentum();//GeV
    pos_t[0]=mcParticleVertex.x;
    pos_t[1]=mcParticleVertex.y;
    pos_t[2]=mcParticleVertex.z;//mm
    mom_t[0]=mcParticleMom.x;
    mom_t[1]=mcParticleMom.y;
    mom_t[2]=mcParticleMom.z;//GeV
    getCircleFromPosMom(pos_t,mom_t,Bz,mcParticle.getCharge(),r,xc,yc);
  }
  const TLorentzVector seedPos=genfitTrack->getSeedStatePos();
  const TVector3 seedMom=genfitTrack->getSeedStateMom();
  TLorentzVector posInit;
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
    endPointStart.SetX(endPointStart.X()/dd4hep::mm);
    endPointStart.SetY(endPointStart.Y()/dd4hep::mm);
    endPointStart.SetZ(endPointStart.Z()/dd4hep::mm);
    endPointEnd.SetX(endPointEnd.X()/dd4hep::mm);
    endPointEnd.SetY(endPointEnd.Y()/dd4hep::mm);
    endPointEnd.SetZ(endPointEnd.Z()/dd4hep::mm);
    bool stopAtBoundary=false;
    bool calcJacobianNoise=true;
    genfitTrack->extrapolateToHit(poca,pocaDir,pocaOnWire,doca,pos_t,mom_t,
        endPointStart,endPointEnd,0,stopAtBoundary,calcJacobianNoise);
    m_dcDigiDocaExt[iDCDigi]=doca;
    m_dcDigiPocaExtX[iDCDigi]=poca.X();
    m_dcDigiPocaExtY[iDCDigi]=poca.Y();
    m_dcDigiPocaExtZ[iDCDigi]=poca.Z();
    double q;
    TVector3 center(xc,yc,0);
    TVector3 wire(endPointStart.X(),endPointStart.Y(),0);
    double docaIdeal=(center-wire).Mag()-r;

    double r_tk,xc_tk,yc_tk;
    getCircleFromTrackState(dcTrack.getTrackStates(0),r_tk,xc_tk,yc_tk,q);
    TVector3 center_tk(xc_tk,yc_tk,0);
    double docaIdeal_tk=(center_tk-wire).Mag()-r_tk;
    m_dcDigiDocaIdeal[iDCDigi]=docaIdeal;
    if(m_debug.value()>0){
      std::cout.precision(13);
      std::cout<<"truth EXT to hit";
      std::cout<<"("<<m_decoder->get(dcDigi.getCellID(),"layer")
        <<","<<m_decoder->get(dcDigi.getCellID(),"cellID")<<")"<<std::endl;
      std::cout<<" wire "<<wire.X()<<" "<<wire.Y()
        <<" seed("<<seedPos.X()<<" "<<seedPos.Y()
        <<" "<<seedPos.Z()<<" "<<seedMom.X()<<" "<<seedMom.Y() <<" "<<seedMom.Z()
        <<") xc,yc,r MC "<<xc<<" "<<yc<<" "<<r
        <<" xc_tk,yc_tk,r_tk "<<xc_tk<<" "<<yc_tk<<" "<<r_tk
        <<" doca ideal_tk "<<docaIdeal_tk<<" mm "
        <<" doca ideal "<<docaIdeal<<" mm "<<std::endl;//yzhang debug
    }
    iDCDigi++;
  }
}

void RecGenfitAlgDC::debugEvent(const edm4hep::TrackCollection* sdtTrackCol,
    const edm4hep::TrackCollection* sdtRecTrackCol,
    double eventStartTime)
{
  m_evt=m_iEvent;
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
    if(m_debug.value()>0){
      debug()<<"sdtTrack charge "<<charge<<" seed mom "<<momInit.X()<<" "<<
        momInit.Y()<<" "<<momInit.Z()<<endmsg;
      debug()<<"sdtTrack "<<charge<<" seed pos "<<posInit.X()<<" "<<
        posInit.Y()<<" "<<posInit.Z()<<endmsg;
      debug()<<"fitted pos sdtTrack? "<<pos.X()<<" "<<pos.Y()<<" "<<pos.Z()<<endmsg;
      debug()<<"fitted mom sdtTrack? "<<mom.X()<<" "<<mom.Y()<<" "<<mom.Z()<<endmsg;
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
    float mcCharge = mcParticle.getCharge();
    helix_mcP.Initialize_VP(mcPos,mcMom,mcCharge,
        m_genfitField->getBzTesla(mcPos));
    mcP_D0 = helix_mcP.getD0();
    mcP_phi = helix_mcP.getPhi0();
    mcP_omega = helix_mcP.getOmega();
    mcP_Z0 = helix_mcP.getZ0();
    mcP_tanLambda = helix_mcP.getTanLambda();
    m_seedCenterX=helix_mcP.getXC();
    m_seedCenterY=helix_mcP.getYC();
    m_seedR=helix_mcP.getRadius();

    float px=mcPocaMom.x;
    float py=mcPocaMom.y;
    float pz=mcPocaMom.z;
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
    TVector3 endPointStart(0,0,0);
    TVector3 endPointEnd(0,0,0);
    //std::cout<<"RecGenfitAlgDC debugEvent iHit "<<iHit<<" simDCHit "<<simDCHit.getCellID()<<std::endl;
    m_gridDriftChamber->cellposition(simDCHit.getCellID(),endPointStart,
        endPointEnd);
    m_dcHitMcWireX[iHit]=endPointStart.X();
    m_dcHitMcWireY[iHit]=endPointStart.Y();
    iHit++;
  }
  //std::cout<<" yzhang debug RecGenfitAlgDC "<<__LINE__<<" simDCHitCol->size()"<<simDCHitCol->size()<<std::endl;
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
    //m_dcDigiDrift[iDCDigi]=dcDigi.getTime()*m_driftVelocity.value()/10000.; //cm
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
    //edm4hep::ConstSimTrackerHit dcSimTrackerHit;
    auto dcSimTrackerHit=CEPC::getAssoClosestSimTrackerHit(m_DCHitAssociationCol.get(),dcDigi,m_gridDriftChamber,0);
    //const edm4hep::MCRecoTrackerAssociationCollection* assoHits=m_DCHitAssociationCol.get();
    m_dcDigiMcMomX[iDCDigi]=dcSimTrackerHit.getMomentum().x*dd4hep::mm;
    m_dcDigiMcMomY[iDCDigi]=dcSimTrackerHit.getMomentum().y*dd4hep::mm;
    m_dcDigiMcMomZ[iDCDigi]=dcSimTrackerHit.getMomentum().z*dd4hep::mm;
    m_dcDigiMcPosX[iDCDigi]=dcSimTrackerHit.getPosition().x*dd4hep::mm;
    m_dcDigiMcPosY[iDCDigi]=dcSimTrackerHit.getPosition().y*dd4hep::mm;
    m_dcDigiMcPosZ[iDCDigi]=dcSimTrackerHit.getPosition().z*dd4hep::mm;
    m_dcDigiDocaMC[iDCDigi]=dcDigi.getTime()*m_driftVelocity.value()/10000.;//cm
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
      <<m_decoder->get(dcDigi.getCellID(),"layer")<<" "
        <<m_decoder->get(dcDigi.getCellID(),"cellID")<<") truth mom "
        <<firstMom<<" "<<dcSimTrackerHit.getMomentum().x<<" "
        <<dcSimTrackerHit.getMomentum().y<<" "
        <<dcSimTrackerHit.getMomentum().z<<" "
        <<std::endl;//yzhang debug

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
        if(m_debug)debug()<<"debugEvent2 errorCov "<<j<<" "<<errorCov[j]<<endmsg;
      }
      m_D0 = trackStat.D0;
      m_phi = trackStat.phi;
      m_omega = trackStat.omega;
      m_Z0 = trackStat.Z0;
      m_tanLambda = trackStat.tanLambda;
    }
  }

}

//void RecGenfitAlgDC::debugdEdx()
//{
//    //const int c_monopolePDGCode = 99666;
//    int pdg[5]={11,13,211,321,2212};
//    for(int ipdg=0;ipdg<5;ipdg++){
//        int pdg_ = pdg[ipdg];
//        TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(pdg_);
//        double charge_ = part->Charge() / 3.;  // We only ever use the square
//        double mass_ = part->Mass(); // GeV
//
//        double stepSize_ = 1;
//
//        GenfitMaterialInterface* materialInterface_=m_genfitFitter->getGeoMaterial();
//        materialInterface_->initTrack(82, 0, 0, 100, 0, 0);
//        auto currentMaterial = materialInterface_->getMaterialParameters();
//        genfit::MaterialEffects* materialEffects_=m_genfitFitter->getMaterialEffects();
//        std::cout<<" print currentMaterial "<<std::endl;
//        currentMaterial.Print();
//        double matDensity_ = currentMaterial.density;
//        double matZ_ = currentMaterial.Z;
//        double matA_ = currentMaterial.A;
//        double radiationLength_ = currentMaterial.radiationLength;
//        double mEE_ = currentMaterial.mEE;
//
//        double minMom = 0.00001;
//        double maxMom = 10000;
//        int nSteps(10000);
//        double logStepSize = (log10(maxMom) - log10(minMom)) / (nSteps-1);
//
//        TH1D hdEdxBethe("dEdxBethe", "dEdxBethe; log10(mom)", nSteps, log10(minMom), log10(maxMom));
//        TH1D hdEdxBrems("dEdxBrems", "dEdxBrems; log10(mom)", nSteps, log10(minMom), log10(maxMom));
//
//        for (int i=0; i<nSteps; ++i) {
//            double mom = pow(10., log10(minMom) + i*logStepSize);
//            double E = hypot(mom, mass_);
//            //if (pdg_ == c_monopolePDGCode) {
//            //    charge_ = mag_charge_ * mom / E; //effective charge for monopoles
//            //}
//
//            //energyLossBrems_ = false;
//            //energyLossBetheBloch_ = true;
//
//            try {
//                //hdEdxBethe.Fill(log10(mom), materialEffects_->dEdx(E));
//            }
//            catch (...) {
//
//            }
//
//
//            //debugOut<< "E = " << E << "; dEdx = " << dEdx(E) <<"\n";
//
//            //energyLossBrems_ = true;
//            //energyLossBetheBloch_ = false;
//            try {
//                //hdEdxBrems.Fill(log10(mom), materialEffects_->dEdx(E));
//            }
//            catch (...) {
//
//            }
//        }
//
//        //energyLossBrems_ = true;
//        //energyLossBetheBloch_ = true;
//
//        std::string Result;//string which will contain the result
//        std::stringstream convert; // stringstream used for the conversion
//        convert << pdg;//add the value of Number to the characters in the stream
//        Result = convert.str();//set Result to the content of the stream
//        TFile outfile("dEdx_" + TString(Result) + ".root", "recreate");
//
//        outfile.cd();
//        hdEdxBethe.Write();
//        hdEdxBrems.Write();
//        outfile.Close();
//    }
//
//} /* End of namespace genfit */



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
////get track position and momentum from TrackState unit cm
//void RecGenfitAlgDC::getCircleFromPosMom(TVector3 pos,TVector3 mom,double charge,
//        double Bz, double& r, double& xc, double& yc)
//{
//    double referencePoint[3];
//    double momentum[3];
//    referencePoint[0] = (double) pos[0];
//    referencePoint[1] = (double) pos[1];
//    referencePoint[2] = (double) pos[2];
//    momentum[0] = (double) mom[0];
//    momentum[1] = (double) mom[1];
//    momentum[2] = (double) mom[2];
//    double pxy = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
//    double FCT = 2.99792458E-4;
//    r= pxy / (FCT*Bz);
//    double phiMomRefPoint = atan2(mom[1],mom[0]);
//    double const_pi2 = 0.5*M_PI;
//    xc= pos[0] + r*cos(phiMomRefPoint-const_pi2*charge);
//    yc= pos[1] + r*sin(phiMomRefPoint-const_pi2*charge);
//
//}
//get track position and momentum from TrackState
void RecGenfitAlgDC::getCircleFromTrackState(const edm4hep::TrackState& trackState,
    double& r, double& xc, double& yc,double& charge)
{
  double omega=trackState.omega;
  double phi=trackState.phi;
  double referencePointX=trackState.referencePoint[0];
  double referencePointY=trackState.referencePoint[1];
  charge = omega/fabs(omega);
  r = 1./fabs(omega);
  xc = referencePointX+r*cos(phi-0.5*M_PI*charge);
  yc = referencePointY+r*sin(phi-0.5*M_PI*charge);
}
//unit length is mm
void RecGenfitAlgDC::getCircleFromPosMom(double pos[3],double mom[3],
    double Bz,double q,double& helixRadius,double& helixXC,double& helixYC)
{
  double FCT = 2.99792458E-4;//mm
  double pxy = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
  helixRadius = pxy / (FCT*Bz);
  double phiMomRefPoint = atan2(mom[1],mom[0]);
  helixXC= pos[0] + helixRadius*cos(phiMomRefPoint-M_PI*0.5*q);
  helixYC= pos[1] + helixRadius*sin(phiMomRefPoint-M_PI*0.5*q);
}

void RecGenfitAlgDC::lsFit(bool smearTrack,bool smearHit,bool useFirstHit){
  ///get helix parameter from MC
  double xc,yc,r;
  double Bz=3.0;
  double pos_t[3];
  double pos_t_first[3];
  double mom_t[3];
  double mom_t_first[3];
  double charge=-1;

  const edm4hep::SimTrackerHitCollection* simDCHitCol=nullptr;
  simDCHitCol=m_simDCHitCol.get();
  for(auto simDCHit: *simDCHitCol){
      pos_t_first[0]=simDCHit.getPosition().x;
      pos_t_first[1]=simDCHit.getPosition().y;
      pos_t_first[2]=simDCHit.getPosition().z;
      mom_t_first[0]=simDCHit.getMomentum().x;
      mom_t_first[1]=simDCHit.getMomentum().y;
      mom_t_first[2]=simDCHit.getMomentum().z;
      charge=-1;//FIXME
      break;
  }
  const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
  mcParticleCol=m_mcParticleCol.get();
  if(nullptr==mcParticleCol)return;
  for(auto mcParticle : *mcParticleCol){
      const edm4hep::Vector3d mcParticleVertex=mcParticle.getVertex();//mm
      const edm4hep::Vector3f mcParticleMom=mcParticle.getMomentum();//GeV
      pos_t[0]=mcParticleVertex.x;
      pos_t[1]=mcParticleVertex.y;
      pos_t[2]=mcParticleVertex.z;//mm
      mom_t[0]=mcParticleMom.x;
      mom_t[1]=mcParticleMom.y;
      mom_t[2]=mcParticleMom.z;
      charge=mcParticle.getCharge();
  }
  if(useFirstHit){
      getCircleFromPosMom(pos_t_first,mom_t_first,Bz,charge,r,xc,yc);
  }else{
      getCircleFromPosMom(pos_t,mom_t,Bz,charge,r,xc,yc);
  }
  debug()<<" getCircleFromPosMom radius "<<r<<" xc "<<xc<<" yc "<<yc<<std::endl;
  if(m_tuple){
      m_momLsMCP=sqrt(mom_t[0]*mom_t[0]+mom_t[1]*mom_t[1]+mom_t[2]*mom_t[2]);
      m_momLsMCPFirst=sqrt(mom_t_first[0]*mom_t_first[0]+mom_t_first[1]*mom_t_first[1]+mom_t_first[2]*mom_t_first[2]);
  }

  LSFitting* aLSFitting=new LSFitting();
  aLSFitting->clear();
  const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
  dCDigiCol=m_digiDCHitsCol.get();
  debug()<<" getCircleFromPosMom radius "<<r<<" xc "<<xc<<" yc "<<yc<<std::endl;
  for(auto dcDigi: *dCDigiCol){
      TVector3 endPointStart(0,0,0);
      TVector3 endPointEnd(0,0,0);
      m_gridDriftChamber->cellposition(dcDigi.getCellID(),endPointStart,
              endPointEnd);//cm
      endPointStart.SetX(endPointStart.X()/dd4hep::mm);
      endPointStart.SetY(endPointStart.Y()/dd4hep::mm);
      endPointStart.SetZ(endPointStart.Z()/dd4hep::mm);
      aLSFitting->setWire(endPointStart.X(),endPointStart.Y());//mm
      double drift=dcDigi.getTime()*m_driftVelocity.value()/1000.; //mm
      if(smearHit) drift+=gRandom->Gaus(m_sigmaDrift.value());//mm
      aLSFitting->setDrift(drift);
      TVector3 center(xc,yc,0);
  }
  if(m_debug>0)aLSFitting->print();
  double lsXc,lsYc,lsRadius;
  lsXc=xc;
  lsYc=yc;
  lsRadius=r;
  if(smearTrack){
      lsXc+=gRandom->Gaus(0,xc*0.05);//mm
      lsYc+=gRandom->Gaus(0,yc*0.05);
      lsRadius+=gRandom->Gaus(0,r*0.05);
  }
  if(m_debug) std::cout<<" ls before fit "<<lsXc<<" "<<lsYc<<" "<<lsRadius<<std::endl;
  aLSFitting->Fitting(lsXc,lsYc,lsRadius);
  if(m_debug) std::cout<<" ls after fit "<<lsXc<<" "<<lsYc<<" "<<lsRadius<<std::endl;
  HelixClass helix;
  double bZ=1;
  double phi0=0;
  double B=3.;
  double signPz=-1;
  double zBegin=0;
  helix.Initialize_BZ(lsXc,lsYc,lsRadius,bZ,phi0,B,signPz,zBegin);
  const float* mom;
  const float* pos;
  mom=helix.getMomentum();
  pos=helix.getReferencePoint();
  m_momLsP=sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
  m_momLsPt=sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
  m_momLsPz=mom[2];
  if(m_debug) std::cout<<" lsFit mom "<<m_momLsP<<" "<<m_momLsPt<<" "<<m_momLsPz<<std::endl;
  if(m_debug) std::cout<<" lsFit pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
  delete aLSFitting;
}
