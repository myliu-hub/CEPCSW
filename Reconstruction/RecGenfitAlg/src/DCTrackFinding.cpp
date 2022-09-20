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
            ISvcLocator* pSvcLocator):GaudiAlgorithm(name, pSvcLocator)
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
    declareProperty("VXDTrackerHitAssociation",m_VXDHitAssociationCol,
            "Handle of VXDsimTrackerHit and VXDTrackerHit association collection");
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

        std::cout << "trackerHits_size = " << siTk.trackerHits_size() << std::endl;
        for(unsigned int iHit=1;iHit<siTk.trackerHits_size();iHit++)
        {

            std::cout << "No." << iHit << " is running " << std::endl;
            edm4hep::ConstTrackerHit* trackerHit=
                new edm4hep::ConstTrackerHit(siTk.getTrackerHits(iHit));
            unsigned long long cellID=trackerHit->getCellID();

            TVectorD hitpos(3);
            hitpos[0]=trackerHit->getPosition()[0];
            hitpos[1]=trackerHit->getPosition()[1];
            hitpos[2]=trackerHit->getPosition()[2];

            TMatrixDSym hitCov(3);
            hitCov(0,0)=trackerHit->getCovMatrix().at(0);
            hitCov(1,1)=trackerHit->getCovMatrix().at(2);
            hitCov(2,2)=trackerHit->getCovMatrix().at(5);

            genfit::TrackPoint* trackPoint = new genfit::TrackPoint(
                    new genfit::SpacepointMeasurement(hitpos,hitCov,cellID,iHit,nullptr),&(recoTrackGenfitAccess->getGenfitTrack(*newRecoTrack)));
            //std::cout << __FILE__ << " genfitTrack rep = "  << std::endl;
            //genfitTrack->getTrackRep(iHit)->Print();
            bool insertPoint = recoTrackGenfitAccess->InsertTrackPoint(*newRecoTrack,trackPoint);
            std::cout << " insertPoint = " << insertPoint << std::endl;

            //delete trackPoint;
        }

        m_genfitFitter->processTrack(&(recoTrackGenfitAccess->getGenfitTrack(*newRecoTrack)));

        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << "newRecoTrack =" << newRecoTrack <<std::endl;
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        Belle2::CKFToCDCFindlet::addSeedRecoTrack(newRecoTrack);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;

    }
    std::cout << __FILE__ << " " << __LINE__ << std::endl;

    const Belle2::CDCCKFPath* bestElement = nullptr;
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    std::vector<Belle2::TrackFindingCDC::CDCWireHit> bestElement2;
    std::cout << __FILE__ << " " << __LINE__ << std::endl;

    Belle2::CKFToCDCFindlet m_CKFToCDC;
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    m_CKFToCDC.apply(bestElement2);
    std::cout << __FILE__ << " " << __LINE__ << std::endl;

    // genfit::RKTrackRep* rep = new genfit::RKTrackRep(pdg);
    // genfit::MeasuredStateOnPlane stateInit(rep);
    // rep->setPosMomCov(stateInit, pos.Vect(), mom, covM);

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
