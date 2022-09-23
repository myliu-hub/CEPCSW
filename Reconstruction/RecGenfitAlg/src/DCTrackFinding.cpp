#include "DCTrackFinding.h"


#include "CDCCKFPath.h"
#include "CDCWireHit.h"
#include "CKFToCDCFindlet.h"

//genfit
#include "AbsKalmanFitter.h"
#include "KalmanFitterRefTrack.h"
#include "EventDisplay.h"
#include "ReferenceStateOnPlane.h"
#include "AbsMeasurement.h"
#include "AbsTrackRep.h"
#include "FitStatus.h"
#include "KalmanFitterInfo.h"
#include "KalmanFittedStateOnPlane.h"
#include "RKTrackRep.h"
//#include "PlanarMeasurement.h"
#include "PlanarMeasurementSDT.h"
#include "SpacepointMeasurement.h"
#include "StateOnPlane.h"
#include "RectangularFinitePlane.h"
#include "Track.h"
#include "TrackPoint.h"
#include "MeasuredStateOnPlane.h"
#include "WireMeasurementNew.h"
#include "AbsMeasurement.h"
#include "TrackPoint.h"

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

//ROOT
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TRandom.h"
#include "TVector3.h"

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

DECLARE_COMPONENT( DCTrackFinding )

    /////////////////////////////////////////////////////////////////////
    DCTrackFinding::DCTrackFinding(const std::string& name,
            ISvcLocator* pSvcLocator):GaudiAlgorithm(name, pSvcLocator),
    m_dd4hepDetector(nullptr),m_gridDriftChamber(nullptr),m_decoder(nullptr)
{

    // Rec Si Track
    declareProperty("SubsetTracks", m_siSubsetTrackCol,
            "Handle of silicon subset track collection");
    //DC DigiHit SimHit DCAsso
    declareProperty("DigiDCHitCollection", m_DCDigiCol,
            "Handle of DC digi(TrakerHit) collection");
    declareProperty("DCHitAssociationCollection", m_DCHitAssociationCol,
            "Handle of DCsimTrackerHit and DCTrackerHit association collection");
    declareProperty("SimTrackerHitCollection",m_simDCHitCol,
            "Handle of the DCsimTrackerHit collection");
    //SIT SimTrackHit DCHitAssociation
    declareProperty("SITCollection",m_simSITHitCol,
            "Handle of the SITsimTrackerHit collection");
    declareProperty("SITTrackerHitAssociation",m_SITHitAssociationCol,
            "Handle of SITsimTrackerHit and SITTrackerHit association collection");
    //VXD SimTrackHit DCHitAssociation
    declareProperty("VXDCollection",m_simVXDHitCol,
            "Handle of the VXDsimTrackerHit collection");
    //declareProperty("VXDTrackerHitAssociation",m_VXDHitAssociationCol,
    //        "Handle of VXDsimTrackerHit and VXDTrackerHit association collection");
    //FTD SimTrackHit DCHitAssociation
    declareProperty("FTDCollection",m_simFTDHitCol,
            "Handle of the FTDsimTrackerHit collection");
    declareProperty("FTDTrackerHitAssociation",m_FTDHitAssociationCol,
            "Handle of FTDsimTrackerHit and FTDTrackerHit association collection");
    //SET SimTrackHit DCHitAssociation
    declareProperty("SETCollection",m_simSETHitCol,
            "Handle of the SETsimTrackerHit collection");
    declareProperty("SETTrackerHitAssociation",m_SETHitAssociationCol,
            "Handle of SETsimTrackerHit and SETTrackerHit association collection");

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode DCTrackFinding::initialize()
{
    MsgStream log(msgSvc(), name());
    info()<<" DCTrackFinding initialize()"<<endmsg;

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

    ///Get Readout
    dd4hep::Readout readout=m_dd4hepDetector->readout(m_readout_name);
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

    /// New a genfit fitter
    //m_genfitFitter = new genfit::KalmanFitterRefTrack();
    m_genfitFitter=new GenfitFitter(m_fitterType.toString().c_str()); 
    m_genfitField=new GenfitField(m_dd4hepField);
    m_genfitFitter->setField(m_genfitField);
    m_genfitFitter->setGeoMaterial(m_geomSvc->lcdd(),m_extMinDistCut,
            m_skipWireMaterial);
    m_genfitFitter->setEnergyLossBrems(m_correctBremsstrahlung);
    m_genfitFitter->setNoiseBrems(m_correctBremsstrahlung);
    if(m_noMaterialEffects) m_genfitFitter->setNoEffects(true);
    if(-1==m_debugPid) m_genfitFitter->setNoEffects(true);
    if(-1==m_debugPid) m_debugPid=0;
    if(m_fitterType=="DAF"||m_fitterType=="DafRef"){
        m_genfitFitter->setMaxIterationsBetas(m_bStart,m_bFinal,m_maxIteration);
    } else {
        m_genfitFitter->setMaxIterations(m_maxIteration);
    }
    //print genfit parameters
    m_genfitFitter->print();
    if(""!=m_genfitHistRootName) m_genfitFitter->initHist(m_genfitHistRootName);


    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode DCTrackFinding::execute()
{
    info()<<"DCTrackFinding in execute()"<<endmsg;

    const edm4hep::TrackCollection* siTrackCol=nullptr;
    siTrackCol=m_siSubsetTrackCol.get();
    if(nullptr==siTrackCol){
        debug()<<"SDTTrackCollection is empty"<<endmsg;
        return StatusCode::SUCCESS;
    }else{
        debug()<<"SiSubsetTrackCol size "<<siTrackCol->size()<<endmsg;
        std::cout <<"SiSubsetTrackCol size "<<siTrackCol->size()<<std::endl;;
        if(0==siTrackCol->size()){
            return StatusCode::SUCCESS;
        }
    }

    //DC
    auto assoDCHitsCol=m_DCHitAssociationCol.get();
    std::cout << " assoDCHitsCol size = " << assoDCHitsCol->size() << std::endl;
    const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
    dCDigiCol=m_DCDigiCol.get();
    std::cout << " dCDigiCol size = " << dCDigiCol->size() << std::endl;

    //SIT
    auto assoSITHitsCol=m_SITHitAssociationCol.get();
    std::cout << " assoSITHitsCol size = " << assoSITHitsCol->size() << std::endl;
    const edm4hep::SimTrackerHitCollection* simSITHitCol=nullptr;
    simSITHitCol=m_simSITHitCol.get();
    std::cout << " simSITHitCol size = " << simSITHitCol->size() << std::endl;

    for(auto siTk:*siTrackCol){

        edm4hep::TrackState trackStat=siTk.getTrackStates(1);
        TVector3 seedPos,seedMom;
        TMatrixDSym seedCov;
        double charge_double;
        TVectorD seedState(6);
        TMatrixDSym covMInit_6(6);
        GenfitTrack::getTrackFromEDMTrackFinding(siTk,charge_double,seedState,covMInit_6,
                seedPos,seedMom);

        int charge = (int) charge_double;
        std::cout << " seed Pos = " << std::endl;
        seedPos.Print();
        std::cout << " seed MOm = " << std::endl;
        seedMom.Print();
        std::cout << " seed State= " << std::endl;
        seedState.Print();

        Belle2::RecoTrack *newRecoTrack = new Belle2::RecoTrack(seedPos,seedMom,charge);
        Belle2::RecoTrackGenfitAccess *recoTrackGenfitAccess = new Belle2::RecoTrackGenfitAccess();
        genfit::AbsTrackRep* RecoRep = recoTrackGenfitAccess->createOrReturnRKTrackRep(*newRecoTrack,-11);

        std::vector<Belle2::TrackFindingCDC::CDCWireHit> bestElement2;

        edm4hep::ConstMCParticle mcParticle;
        bool flag = false;
        for(unsigned int iHit=0;iHit<siTk.trackerHits_size();iHit++)
        {
            std::cout << "No." << iHit << " is running " << std::endl;
            edm4hep::ConstTrackerHit* trackerHit=
                new edm4hep::ConstTrackerHit(siTk.getTrackerHits(iHit));
            unsigned long long cellID=trackerHit->getCellID();
            ///getDetTypeID
            UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
            encoder.setValue(cellID);
            int detTypeID=encoder[lcio::ILDCellID0::subdet];

            if(!(detTypeID==lcio::ILDDetID::SIT)) continue;
            std::string detectorName = "SIT";
            for(int iSim=0; iSim <(int) assoSITHitsCol->size();iSim++)
            {
                if(assoSITHitsCol->at(iSim).getRec()== *trackerHit)
                {
                    mcParticle = assoSITHitsCol->at(iSim).getSim().getMCParticle();
                    flag = true;
                    break;
                }
            }
            if(flag) break;
        }
        std::cout << " MCParticle = " << mcParticle << std::endl;

        // DC digi hit find MCParticle
        int ndigi =0;
        for(auto dcDigi: *dCDigiCol){
            TVectorD hitpos(3);
            TMatrixDSym hitCov(3);

            unsigned short tdcCount = dcDigi.getTime();
            unsigned short adcCount = 0;
            unsigned short charge;
            unsigned short iSuperLayer=0;
            
            unsigned short iLayer=0;
            unsigned short iWire=0;

            genfit::TrackPoint* trackPoint =nullptr;
            bool insertPoint =false;
            for(int iSim=0; iSim <(int) assoDCHitsCol->size();iSim++)
            {
                if(assoDCHitsCol->at(iSim).getRec()== dcDigi &&
                        (assoDCHitsCol->at(iSim).getSim().getMCParticle() == mcParticle))
                {
                    std::cout << "DC  MCParticle = " << assoDCHitsCol->at(iSim).getSim().getMCParticle() << std::endl;
                    std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    unsigned long long dccellID=dcDigi.getCellID();
                    std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    iLayer = m_decoder->get(dccellID,"layer");
                    std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    iWire = m_decoder->get(dccellID,"cellID");
                    std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    charge = assoDCHitsCol->at(iSim).getSim().getMCParticle().getCharge();
                    std::cout << __FILE__ << " " << __LINE__ << std::endl;

                    hitpos[0]=dcDigi.getPosition()[0];
                    hitpos[1]=dcDigi.getPosition()[1];
                    hitpos[2]=dcDigi.getPosition()[2];

                    hitCov(0,0)=dcDigi.getCovMatrix().at(0);
                    hitCov(1,1)=dcDigi.getCovMatrix().at(2);
                    hitCov(2,2)=dcDigi.getCovMatrix().at(5);

                    std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    trackPoint = new genfit::TrackPoint(
                            new genfit::SpacepointMeasurement(hitpos,hitCov,dccellID,ndigi,nullptr),
                            &(recoTrackGenfitAccess->getGenfitTrack(*newRecoTrack)));
                    std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    insertPoint = recoTrackGenfitAccess->InsertTrackPoint(*newRecoTrack,trackPoint);
                    std::cout << __FILE__ << " " << __LINE__ << std::endl;
                    std::cout << " insertPoint = " << insertPoint << std::endl;
                    ndigi++;
                    break;
                    std::cout << __FILE__ << " " << __LINE__ << std::endl;
                }//loop if
                std::cout << __FILE__ << " " << __LINE__ << std::endl;
            }// loop for Sim Hit over
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            Belle2::CDCHit * cdchit = new Belle2::CDCHit(tdcCount,adcCount,
                    iSuperLayer,iLayer,iWire);
            std::cout << " getIWire() = " << cdchit->getIWire() << " getILayer() = " << cdchit->getILayer() << " getICLayer() = " << cdchit->getICLayer() << " getISuperLayer() = " << cdchit->getISuperLayer() << " getID() = " << cdchit->getID() << std::endl;
            std::cout << " tdcCount = " << tdcCount << " adcCount = " << adcCount << " iSuperLayer = " << iSuperLayer << " iLayer = " << iLayer << " iWire = " << iWire << std::endl;
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            //Belle2::CDCHit cdchit(tdcCount,adcCount,iSuperLayer,iLayer,iWire);
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            Belle2::TrackFindingCDC::CDCWireHit * cdcWireHit =
                new Belle2::TrackFindingCDC::CDCWireHit(cdchit,
                        m_driftVelocity*dcDigi.getTime(),m_sigma,0,dcDigi.getTime());
            std::cout << __FILE__ << " " << __LINE__ << std::endl;

            bestElement2.push_back(*cdcWireHit);
            std::cout << __FILE__ << " " << __LINE__ << std::endl;

        }
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << " ndigi = " << ndigi << std::endl;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        if(ndigi<1e-5) return StatusCode::FAILURE;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        m_genfitFitter->processTrack(&(recoTrackGenfitAccess->getGenfitTrack(*newRecoTrack)));
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        Belle2::CKFToCDCFindlet::addSeedRecoTrack(newRecoTrack);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        const Belle2::CDCCKFPath* bestElement = nullptr;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

        Belle2::CKFToCDCFindlet m_CKFToCDC;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        m_CKFToCDC.apply(bestElement2);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

    }

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode DCTrackFinding::finalize()
{
    MsgStream log(msgSvc(), name());
    info()<< " DCTrackFinding in finalize()" << endmsg;

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
