#include<stdlib.h>
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
#include <fstream>

#include <TTimeStamp.h>

#include <ctime>
#include <cstdlib>

#include <array>

#include "map"
#include "DataHelper/GeomeryHelper.h"

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
    declareProperty("TruthDigiDCHitCollection", m_TruthDCDigiCol,
            "Handle of Truth DC digi(TrakerHit) collection");
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
    declareProperty("VXDTrackerHitAssociation", m_VXDHitAssociationCol,
            "Handle of VXD association collection");
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

    // Output collections
    declareProperty("DCTrackFindingHitCollection",w_DCTrackFindingCol,
            "Handle of TrackFinding DCHit collection");

    //SDTTrackCollection
    declareProperty("SDTTrackCollection", m_SDTTrackCol,
            "Handle of input silicon track collection");
    //SDTTrackFindCollection
    declareProperty("SDTTrackFindCollection", m_SDTTrackFindCol,
            "Handle of input silicon track finding collection");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode DCTrackFinding::initialize()
{
    MsgStream log(msgSvc(), name());
    info()<<" DCTrackFinding initialize()"<<endmsg;

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
    if(m_fitterType=="DAF"||m_fitterType=="DafRef"){
        m_genfitFitter->setMaxIterationsBetas(m_bStart,m_bFinal,m_maxIteration);
    } else {
        m_genfitFitter->setMaxIterations(m_maxIteration);
    }
    //print genfit parameters
    m_genfitFitter->print();
    if(""!=m_genfitHistRootName) m_genfitFitter->initHist(m_genfitHistRootName);

    ///book tuple
    NTuplePtr nt(ntupleSvc(), "DCTrackFinding/dcTrackFinding");
    if(nt){
        m_tuple=nt;
    }else{
        m_tuple=ntupleSvc()->book("DCTrackFinding/dcTrackFinding",
                CLID_ColumnWiseTuple,"DCTrackFinding");
        if(m_tuple){
            StatusCode sc;
            sc=m_tuple->addItem("N_SiTrack",m_n_N_SiTrack,0,100);
            sc=m_tuple->addItem("N_DCTrackFinding",m_n_N_SiTrack,m_n_DCTrackFinding);
            sc=m_tuple->addItem("CPUTime",m_n_N_SiTrack,m_CPUTime);
            sc=m_tuple->addItem("layer",m_n_N_SiTrack,m_layer,100);
            sc=m_tuple->addItem("cellID",m_n_N_SiTrack,m_cellID,100);
            sc=m_tuple->addItem("NoiseFlag",m_n_N_SiTrack,m_NoiseFlag,100);
            sc=m_tuple->addItem("isNoOverlap",m_n_N_SiTrack,m_isNoOverlap);
            sc=m_tuple->addItem("N_digiOverlap",m_n_N_SiTrack,m_digiOverlap);
            sc=m_tuple->addItem("hit_wx",m_n_N_SiTrack,m_hit_wx,100);
            sc=m_tuple->addItem("hit_wy",m_n_N_SiTrack,m_hit_wy,100);
            sc=m_tuple->addItem("hit_wz",m_n_N_SiTrack,m_hit_wz,100);
            sc=m_tuple->addItem("hit_ex",m_n_N_SiTrack,m_hit_ex,100);
            sc=m_tuple->addItem("hit_ey",m_n_N_SiTrack,m_hit_ey,100);
            sc=m_tuple->addItem("hit_ez",m_n_N_SiTrack,m_hit_ez,100);
            sc=m_tuple->addItem("recoPosX",m_n_N_SiTrack,m_recoPosX,100);
            sc=m_tuple->addItem("recoPosY",m_n_N_SiTrack,m_recoPosY,100);
            sc=m_tuple->addItem("recoPosZ",m_n_N_SiTrack,m_recoPosZ,100);
            sc=m_tuple->addItem("cellPosX",m_n_N_SiTrack,m_cellPosX,100);
            sc=m_tuple->addItem("cellPosY",m_n_N_SiTrack,m_cellPosY,100);
            sc=m_tuple->addItem("cellPosZ",m_n_N_SiTrack,m_cellPosZ,100);
            sc=m_tuple->addItem("radius",m_n_N_SiTrack,m_radius,100);
            sc=m_tuple->addItem("arcLength",m_n_N_SiTrack,m_arcLength,100);
            sc=m_tuple->addItem("hitDistance",m_n_N_SiTrack,m_hitDistance,100);
            sc=m_tuple->addItem("doca",m_n_N_SiTrack,m_doca,100);
            sc=m_tuple->addItem("ReconstructedZ",m_n_N_SiTrack,m_ReconstructedZ,100);
            sc=m_tuple->addItem("trackParD0",m_n_N_SiTrack,m_trackParD0,100);
            sc=m_tuple->addItem("trackParPhi0",m_n_N_SiTrack,m_trackParPhi0,100);
            sc=m_tuple->addItem("trackParOmega",m_n_N_SiTrack,m_trackParOmega,100);
            sc=m_tuple->addItem("trackParZ0",m_n_N_SiTrack,m_trackParZ0,100);
            sc=m_tuple->addItem("trackParTanLambda",m_n_N_SiTrack,m_trackParTanLambda,100);
            sc=m_tuple->addItem("CircleCenterX",m_n_N_SiTrack,m_centerX,100);
            sc=m_tuple->addItem("CircleCenterY",m_n_N_SiTrack,m_centerY,100);
            sc=m_tuple->addItem("N_Digi",m_n_Digi,0,100000);
            sc=m_tuple->addItem("N_TruthDigi",m_n_TruthDigi,0,100000);
            sc=m_tuple->addItem("N_match",m_n_N_SiTrack,m_n_match);
            sc=m_tuple->addItem("match_hit_wx",m_n_N_SiTrack,m_match_hit_wx,100);
            sc=m_tuple->addItem("match_hit_wy",m_n_N_SiTrack,m_match_hit_wy,200);
            sc=m_tuple->addItem("match_hit_wz",m_n_N_SiTrack,m_match_hit_wz,200);
            sc=m_tuple->addItem("match_hit_ex",m_n_N_SiTrack,m_match_hit_ex,200);
            sc=m_tuple->addItem("match_hit_ey",m_n_N_SiTrack,m_match_hit_ey,200);
            sc=m_tuple->addItem("match_hit_ez",m_n_N_SiTrack,m_match_hit_ez,200);
        }else{
            warning()<<"Tuple DCTrackFinding/dcTrackFinding not booked"<<endmsg;
        }
    } //end of book tuple 

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode DCTrackFinding::execute()
{
    info()<<"DCTrackFinding in execute()"<<endmsg;

    StatusCode sc=StatusCode::SUCCESS;

    edm4hep::TrackerHitCollection* TrF = w_DCTrackFindingCol.createAndPut();
    
    edm4hep::TrackCollection* sdtTkFCol= m_SDTTrackFindCol.createAndPut();

    const edm4hep::TrackCollection* siTrackCol=nullptr;
    siTrackCol=m_siSubsetTrackCol.get();
    if(nullptr==siTrackCol){
        debug()<<"SDTTrackCollection is empty"<<endmsg;
        return StatusCode::SUCCESS;
    }else{
        debug()<<"SiSubsetTrackCol size "<<siTrackCol->size()<<endmsg;
        if(0==siTrackCol->size()){
            return StatusCode::SUCCESS;
        }
    }

    //DC
    auto assoDCHitsCol=m_DCHitAssociationCol.get();
    const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
    dCDigiCol=m_DCDigiCol.get();

    // Truth DC hit
    const edm4hep::TrackerHitCollection* TruthdCDigiCol=m_TruthDCDigiCol.get();
    int truthNumdc = TruthdCDigiCol->size();
    m_n_TruthDigi = truthNumdc;
   //VXDAssoHits
    const edm4hep::MCRecoTrackerAssociationCollection* VXDAssoVec = m_VXDHitAssociationCol.get();
    //SIT
    auto assoSITHitsCol=m_SITHitAssociationCol.get();
    const edm4hep::SimTrackerHitCollection* simSITHitCol=nullptr;
    simSITHitCol=m_simSITHitCol.get();

    //FTD
    auto assoFTDHitsCol=m_FTDHitAssociationCol.get();
    const edm4hep::SimTrackerHitCollection* simFTDHitCol=nullptr;
    simFTDHitCol=m_simFTDHitCol.get();

    //SDTtrackCollection
    const edm4hep::TrackCollection* sdtTrackCol=nullptr;
    if(m_SDTTrackCol.exist())sdtTrackCol=m_SDTTrackCol.get(); 
    if(nullptr==sdtTrackCol || sdtTrackCol->size()<=0) {
        std::cout<<"TrackCollection not found or sdtTrackCol size=0"<<std::endl;
        return StatusCode::SUCCESS;
    }

    std::vector<std::vector<int>> truthDChit;
    for(auto trudchit: *TruthdCDigiCol){
        unsigned short iLayer=0;
        unsigned short iWire=0;

        unsigned long long dccellID=trudchit.getCellID();
        iLayer = m_decoder->get(dccellID,"layer");
        iWire = m_decoder->get(dccellID,"cellID");

        std::vector<int> mid;
        mid.push_back(iLayer);
        mid.push_back(iWire);

        truthDChit.push_back(mid);

        mid.clear();

    }

    int siNum = 0;
    for(auto sdtTk:*sdtTrackCol){

        edm4hep::MutableTrack sdtTFk=sdtTkFCol->create();

        edm4hep::TrackState trackStat=sdtTk.getTrackStates(0);
        edm4hep::TrackerHit trackHit = sdtTk.getTrackerHits(0);
        edm4hep::MCParticle mcParticle;
        edm4hep::SimTrackerHit simTrackerHit;
        CEPC::getAssoMCParticle(VXDAssoVec,trackHit,mcParticle,simTrackerHit);
        TVector3 seedPos,seedMom;
        TMatrixDSym seedCov;
        double charge_double;
        TVectorD seedState(6);
        TMatrixDSym covMInit_6(6);
        GenfitTrack::getTrackFromEDMTrackFinding(sdtTk,charge_double,seedState,covMInit_6,
                seedPos,seedMom);

        int charge = (int) charge_double;
        int pdg = mcParticle.getPDG();

        Belle2::RecoTrack *newRecoTrack = new Belle2::RecoTrack(seedPos,seedMom,charge);
        Belle2::RecoTrackGenfitAccess *recoTrackGenfitAccess = new Belle2::RecoTrackGenfitAccess();
        genfit::AbsTrackRep* RecoRep = recoTrackGenfitAccess->createOrReturnRKTrackRep(*newRecoTrack,pdg);

        genfit::AbsTrackRep* rep;
        rep = new genfit::RKTrackRep(pdg);
        genfit::MeasuredStateOnPlane state(rep);
        rep->setPosMom(state, seedPos, seedMom);

        const TVector3 linePoint(0,0,0);
        const TVector3 lineDirection(0,0,1);
        const double R = 80.5;//cm

        try { 
            double L = rep->extrapolateToCylinder(state, R, linePoint, lineDirection, false);
        }
        catch (genfit::Exception& e) {
            std::cerr << e.what();
        }

        std::vector<Belle2::TrackFindingCDC::CDCWireHit> bestElement2;

        if(m_tuple) m_n_Digi = dCDigiCol->size();

        // DC digi hit find MCParticle
        int ndigi =0;
        int detID = 7;
        std::map < std::vector<int> , unsigned long long > findCellID;
        edm4hep::TrackerHit digihits[55][1000];
        for(auto dcDigi: *dCDigiCol){

            if(0==dcDigi.getQuality()) m_digiOverlap[siNum]++;
            TVectorD hitpos(3);
            TMatrixDSym hitCov(3);

            hitpos[0]=dcDigi.getPosition()[0];
            hitpos[1]=dcDigi.getPosition()[1];
            hitpos[2]=dcDigi.getPosition()[2];

            hitCov(0,0)=dcDigi.getCovMatrix().at(0);
            hitCov(1,1)=dcDigi.getCovMatrix().at(2);
            hitCov(2,2)=dcDigi.getCovMatrix().at(5);

            unsigned short tdcCount = dcDigi.getTime();
            unsigned short adcCount = 0;
            unsigned short iSuperLayer=0;

            unsigned short iLayer=0;
            unsigned short iWire=0;

            genfit::TrackPoint* trackPoint =nullptr;
            bool insertPoint =false;
            unsigned long long dccellID=dcDigi.getCellID();
            iLayer = m_decoder->get(dccellID,"layer");
            iWire = m_decoder->get(dccellID,"cellID");

            //if(1==dcDigi.getQuality()){
            //    testLayer.push_back(iLayer);
            //}

            digihits[iLayer][iWire] = dcDigi;

            std::vector<int> id;
            id.push_back(iLayer);
            id.push_back(iWire);

            findCellID.insert(std::pair< std::vector<int> , int>(id,dccellID));

            id.clear();

            double driftDistance = 1e-4*m_driftVelocity*dcDigi.getTime(); //cm

            TVector3 Wstart(0,0,0);
            TVector3 Wend  (0,0,0);
            m_gridDriftChamber->cellposition(dcDigi.getCellID(), Wstart, Wend);

            if(1==dcDigi.getQuality()){
                double hitphi = atan(hitpos[1]/hitpos[0]);
                double wirephi = atan(Wstart.Y()/Wstart.X());
                //testPhi.push_back(hitphi);
                //testWirePhi.push_back(wirephi);
            }

            ndigi++;

            Belle2::CDCHit * cdchit = new Belle2::CDCHit(tdcCount,adcCount,
                    iSuperLayer,iLayer,iWire);

            Belle2::TrackFindingCDC::CDCWireHit cdcWireHit(cdchit,
                    //m_driftVelocity*dcDigi.getTime(),m_sigma,0,dcDigi.getTime());
                driftDistance,m_sigma,0,dcDigi.getTime());

            bestElement2.push_back(cdcWireHit);
        }
        if(ndigi<1e-5) return StatusCode::FAILURE;

        Belle2::CKFToCDCFindlet::clearSeedRecoTrack();
        Belle2::CKFToCDCFindlet::addSeedRecoTrack(newRecoTrack);

        Belle2::CKFToCDCFindlet m_CKFToCDC;
        std::vector< std::vector<unsigned short> > output;
        std::vector<double> radius;
        std::vector<double> arcLength;
        std::vector<double> hitDistance;
        std::vector<double> doca;
        std::vector<double> reconstructedZ;
        std::vector<TVector3> recoPos;
        std::vector<std::vector<double>> trackParams;
        std::vector<std::vector<double>> circleCenter;
        m_CKFToCDC.beginEvent();
        m_CKFToCDC.addMeasuredStateOnPlane(state);
        m_CKFToCDC.setpathSelectFilterName(m_pathSelectFilterName);
        m_CKFToCDC.setMaximalLayerJump(m_maximalLayerJump);


        auto start = std::chrono::high_resolution_clock::now(); 
        m_CKFToCDC.apply(bestElement2);

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;

        m_CPUTime[siNum] = elapsed.count();

        m_CKFToCDC.getResult(output,trackParams,
                radius,circleCenter,arcLength,hitDistance,doca,
                reconstructedZ,recoPos);
        m_n_DCTrackFinding[siNum] = output.size();

        for(int i=0;i<radius.size();i++) m_radius[siNum][i] = radius[i];
        for(int i=0;i<arcLength.size();i++) m_arcLength[siNum][i] = arcLength[i];
        for(int i=0;i<hitDistance.size();i++) m_hitDistance[siNum][i] = hitDistance[i];
        for(int i=0;i<doca.size();i++) m_doca[siNum][i] = doca[i];
        for(int i=0;i<reconstructedZ.size();i++) m_ReconstructedZ[siNum][i] = reconstructedZ[i];

        for(int i=0;i<recoPos.size();i++){
            m_recoPosX[siNum][i] = recoPos[i].X();
            m_recoPosY[siNum][i] = recoPos[i].Y();
            m_recoPosZ[siNum][i] = recoPos[i].Z();
        }

        for(int i=0;i<trackParams.size();i++){
            m_trackParD0[siNum][i] = trackParams[i][0];
            m_trackParPhi0[siNum][i] = trackParams[i][1];
            m_trackParOmega[siNum][i] = trackParams[i][2];
            m_trackParZ0[siNum][i] = trackParams[i][3];
            m_trackParTanLambda[siNum][i] = trackParams[i][4];
        }
        for(int i=0;i<circleCenter.size();i++){
            m_centerX[siNum][i] = circleCenter[i][0];
            m_centerY[siNum][i] = circleCenter[i][1];
        }

        int num_match =0;
        //add sdtTk To sdtTFk
        sdtTFk.addToTrackStates(sdtTk.getTrackStates(0));
        for(int i=0;i<sdtTk.trackerHits_size();i++){
            sdtTFk.addToTrackerHits(sdtTk.getTrackerHits(i));
        }

        for(int i = 0;i<output.size();i++)
        {
            auto trkHit = TrF->create();
            edm4hep::TrackerHit digihit = digihits[output[i].at(0)][output[i].at(1)];

            TVector3 Wstart(0,0,0);
            TVector3 Wend  (0,0,0);
            m_gridDriftChamber->cellposition(digihit.getCellID(), Wstart, Wend);

            TVector3 cellPos =  m_gridDriftChamber->wirePos_vs_z2(output[i].at(0),output[i].at(1),
                    reconstructedZ[i]);

            m_cellPosX[siNum][i] = cellPos.X();
            m_cellPosY[siNum][i] = cellPos.Y();
            m_cellPosZ[siNum][i] = cellPos.Z();

            m_layer[siNum][i] = output[i].at(0);
            m_cellID[siNum][i] = output[i].at(1);
            m_NoiseFlag[siNum][i] = digihit.getQuality();

            //match是指跟truth hit匹配上的，不排除有被overlap的可能
            for(int j=0;j<truthDChit.size();j++){
                if((output[i].at(0) == truthDChit[j].at(0)) && 
                        (output[i].at(1) == truthDChit[j].at(1))){

                    m_match_hit_wx[siNum][num_match] = Wstart.X();
                    m_match_hit_wy[siNum][num_match] = Wstart.Y();
                    m_match_hit_wz[siNum][num_match] = Wstart.Z();

                    m_match_hit_ex[siNum][num_match] = Wend.X();
                    m_match_hit_ey[siNum][num_match] = Wend.Y();
                    m_match_hit_ez[siNum][num_match] = Wend.Z();
                    if(1==digihit.getQuality())m_isNoOverlap[siNum]++;

                    num_match++;
                }
            }

            m_hit_wx[siNum][i] = Wstart.X();
            m_hit_wy[siNum][i] = Wstart.Y();
            m_hit_wz[siNum][i] = Wstart.Z();

            m_hit_ex[siNum][i] = Wend.X();
            m_hit_ey[siNum][i] = Wend.Y();
            m_hit_ez[siNum][i] = Wend.Z();

            trkHit.setCellID(digihit.getCellID());
            trkHit.setQuality(digihit.getQuality());
            trkHit.setTime(digihit.getTime());
            trkHit.setEDep(digihit.getEDep());
            trkHit.setPosition(digihit.getPosition());
            trkHit.setCovMatrix(digihit.getCovMatrix());

            sdtTFk.addToTrackerHits(trkHit);
        }

        m_n_match[siNum] = num_match;

        delete newRecoTrack;
        delete recoTrackGenfitAccess;

        for(auto& it : bestElement2){
            delete it.getHit();
        }
        bestElement2.clear();
        siNum++;
    }

    m_n_N_SiTrack = siNum;
    if(m_tuple) sc=m_tuple->write();

    m_eventNo++;

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode DCTrackFinding::finalize()
{
    MsgStream log(msgSvc(), name());
    info()<< " DCTrackFinding in finalize()" << endmsg;
    delete m_genfitFitter;

    std::cout << " End of DCTrackFinding .... " << std::endl;
    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

int DCTrackFinding::getNumPointsWithFittedInfo(genfit::Track genfitTrack,
        int repID) const 
{
    int nHitWithFittedInfo = 0;
    //check number of hit with fitter info
    int nHit = genfitTrack.getNumPointsWithMeasurement();
    //check number of hit with fitter info
    for(int i=0; i<nHit; i++){
        if(nullptr != genfitTrack.getPointWithFitterInfo(i,genfitTrack.getTrackRep(repID))){
            nHitWithFittedInfo++;
        }
    }
    return nHitWithFittedInfo;
}

int DCTrackFinding::getFittedState(genfit::Track genfitTrack,
        TLorentzVector& pos, TVector3& mom,
        TMatrixDSym& cov, int trackPointId, int repID, bool biased) const
{
    //check number of hit with fitter info
    if(getNumPointsWithFittedInfo(genfitTrack,repID)<=2) return 1;

    //get track rep
    genfit::AbsTrackRep* rep = genfitTrack.getTrackRep(repID);
    if(nullptr == rep) return 2;

    //get first or last measured state on plane
    genfit::MeasuredStateOnPlane mop;
    try{
        mop = genfitTrack.getFittedState(trackPointId,rep,biased);
    }catch(genfit::Exception& e){
        std::cout<<" getNumPointsWithFittedInfo="
            <<getNumPointsWithFittedInfo(genfitTrack,repID)
            <<" no TrackPoint with fitted info "<<std::endl;
        return 3;
    }

    //get state
    TVector3 p;
    mop.getPosMomCov(p,mom,cov);
    pos.SetVect(p);
    pos.SetT(9999);//FIXME

    return 0;//success
}

bool DCTrackFinding::getMOP(int hitID,genfit::MeasuredStateOnPlane& mop,
        genfit::AbsTrackRep* trackRep,genfit::Track* track) const
{
    if(nullptr == trackRep) trackRep = getRep(0,track);
    try{
        mop = track->getFittedState(hitID,trackRep);
    }catch(genfit::Exception& e){
        e.what();
        return false;
    }
    return true;
}

genfit::AbsTrackRep* DCTrackFinding::getRep(int id,genfit::Track* track) const
{
    return track->getTrackRep(id);
}

int DCTrackFinding::addHitsToTk(edm4hep::TrackerHitCollection *
        //colHandle, edm4hep::Track& track, const char* msg,int nHitAdded)
col, edm4hep::Track& track, const char* msg) const
{
    //if(nHitAdded>0) return nHitAdded;
    int nHit=0;
    debug()<<"add "<<msg<<" "<<col->size()<<" trackerHit"<<endmsg;
    //sort,FIXME
    edm4hep::MutableTrack sdtTrack = track.clone();
    for(auto hit:*col){
        sdtTrack.addToTrackerHits(hit);
        ++nHit;
    }
    track = sdtTrack;
    return nHit;
}
