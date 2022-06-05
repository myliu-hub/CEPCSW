#include "GenfitTrack.h"
#include "GenfitHit.h"
#include "GenfitField.h"

//CEPCSW
#include "DataHelper/HelixClass.h"
#include "DataHelper/TrackHelper.h"
#include "DataHelper/TrackerHitHelper.h"
#include "DetInterface/IGeomSvc.h"
#include "GenfitUnit.h"
#include "UTIL/ILDConf.h"
#include "WireMeasurementDC.h"

//Externals
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Segmentations.h"
#include "DDRec/ISurface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Vector3D.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackerHitConst.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/Vector3d.h"
#include "GaudiKernel/SmartIF.h"

//genfit
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

//ROOT
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TRandom.h"
#include "TVector3.h"

//cpp
#include <cfloat>

#undef GENFIT_MY_DEBUG
//#define GENFIT_MY_DEBUG 1

const int GenfitTrack::s_PDG[2][5]
={{-11,-13,211,321,2212},{11,13,-211,-321,-2212}};

    bool
sortDCSimHit(edm4hep::ConstSimTrackerHit hit1,edm4hep::ConstSimTrackerHit hit2)
{
    //std::cout<<"hit1"<<hit1<<std::endl;
    //std::cout<<"hit2"<<hit2<<std::endl;
    bool isEarly=hit1.getTime()<hit2.getTime();
    return isEarly;
}
    bool
sortDCDigi(std::pair<double,edm4hep::ConstTrackerHit*> hitPair1,std::pair<double,edm4hep::ConstTrackerHit*> hitPair2)
{
    bool isEarly=hitPair1.first<hitPair2.first;
    return isEarly;
}
//bool sortDCDigiLayer(std::pair<int,edm4hep::ConstTrackerHit> hitPair1,std::pair<int,edm4hep::ConstTrackerHit> hitPair2)
//{
//    bool isEarly=hitPair1.first<hitPair2.first;
//    return isEarly;
//}


GenfitTrack::GenfitTrack(const GenfitField* genfitField,
        const dd4hep::DDSegmentation::GridDriftChamber* seg,
        SmartIF<IGeomSvc> geom, const char* name)
:m_name(name),m_track(nullptr),m_debug(0),m_debugLocal(0),m_geomSvc(geom),
    m_genfitField(genfitField),m_gridDriftChamber(seg),
    m_decoderDC(geom->getDecoder("DriftChamberHitsCollection"))
{
}

GenfitTrack::~GenfitTrack()
{
    clearGenfitHitVec();
    ///Note: track reps and points will be deleted in the destructor of track
    ///implemented in genfit::Track::Clear()
    delete m_track;
}

void GenfitTrack::clearGenfitHitVec(){
    for(auto h:m_genfitHitVec){delete h;}
    m_genfitHitVec.clear();
    std::vector<GenfitHit*>().swap(m_genfitHitVec);
}

/// create a Genfit track from track state
/// Initialize track with seed state and cov
/// unit conversion here!!!
bool GenfitTrack::createGenfitTrack(int pdgType,int charge,
        TVectorD trackParam, TMatrixDSym covMInit_6)
{
    if(m_debug){
        std::cout<<"createGenfitTrack pdgType "<<pdgType<<" charge "<<charge
          <<" trackParam "<<std::endl; trackParam.Print();
        std::cout<<"Cov"<<std::endl; covMInit_6.Print();
    }
    ///new a track and set seed state and cov
    if(nullptr!=m_track) {
        delete m_track;
        clearGenfitHitVec();
    }
    m_track=new genfit::Track();
    m_track->setStateSeed(trackParam);
    m_track->setCovSeed(covMInit_6);

#ifdef GENFIT_MY_DEBUG
    //m_track->setDebugLvlLocal(m_debugLocal);
#endif

    ///new a track representation and add to the track
    addTrackRep(pdgType,charge);

    if(m_debug>0) printSeed();
    return true;
}

///Create a Genfit track with MCParticle, unit conversion here
bool GenfitTrack::createGenfitTrackFromMCParticle(int pidType,
        const edm4hep::MCParticle& mcParticle, double eventStartTime)
{
    ///Get error matrix of seed track
    TMatrixDSym covMInit_6(6);
    getSeedCov(covMInit_6);
    TVectorD seedState(6);
    getTrackFromMCPartile(mcParticle,seedState,covMInit_6);

    return createGenfitTrack(pidType,mcParticle.getCharge(),seedState,covMInit_6);
}//end of createGenfitTrackFromMCParticle

///Create a Genfit track with edm4hep Track
bool GenfitTrack::createGenfitTrackFromEDM4HepTrack(int pidType,
        const edm4hep::Track& track, double eventStartTime, bool isUseCovTrack)
{

    ///Skip track w.o. hit
    if(track.trackerHits_size()<=0) {
        std::cout << " track.trackerHits_size = " << track.trackerHits_size() << std::endl;
        if(m_debug){
          std::cout<<"createGenfitTrackFromEDM4HepTrack skip track w/o hit"<<std::endl;
        }
        return false;
    }

    //pivotToFirstLayer(mcPocaPos,mcPocaMom,firstLayerPos,firstLayerMom);
    ///Get track parameters
    double charge=0;
    TVectorD seedState(6);
    TMatrixDSym covMInit_6(6);

    getTrackFromEDMTrack(track,charge,seedState,covMInit_6);///unit conversion
    if(!isUseCovTrack) getSeedCov(covMInit_6);

    bool status=createGenfitTrack(pidType,charge,seedState,covMInit_6);
    return status;
}

/// Add a 3d SpacepointMeasurement
bool GenfitTrack::addSpacePointMeasurement(const TVector3& pos,
        std::vector<float> sigmaUVec,std::vector<float> sigmaVVec,int cellID,int hitID)
{
    TVectorD pos_smeared(3);
    for(int i=0;i<3;i++) pos_smeared[i]=pos(i)*GenfitUnit::mm;

    /// New a SpacepointMeasurement
    /// space point error matrix and smear, unit cm
    TMatrixDSym hitCov(3);
    hitCov.Zero();

    //get sigmas
    float sigmaU,sigmaV;
    getSigmas(cellID,sigmaUVec,sigmaVVec,sigmaU,sigmaV);
    float sigmaX=sigmaU;//sigmaU*cos(atan2(pos_smeared[1],pos_smeared[0]));
    float sigmaY=sigmaU;//*sin(atan2(pos_smeared[1],pos_smeared[0]));
    float sigmaZ=sigmaV;

    //smear 3d track position
    bool smear=(sigmaUVec[0]>0);
    if(smear){
        pos_smeared[0]+=gRandom->Gaus(0,sigmaX);
        pos_smeared[1]+=gRandom->Gaus(0,sigmaX);
        pos_smeared[2]+=gRandom->Gaus(0,sigmaX);
    }
    hitCov(0,0)=sigmaX*sigmaX;//FIXME
    hitCov(1,1)=sigmaX*sigmaX;
    hitCov(2,2)=sigmaX*sigmaX;

    if(m_debug>=2){
        std::cout<<" addSpacePointMeasurement pos "<<std::endl;
        pos.Print();
        std::cout<<m_name<<" hit "<<hitID
            <<" pos_smeared "<<pos_smeared[0]<<" "<<pos_smeared[1]
            <<" "<<pos_smeared[2]<<" "<<" hitCov (0,0) "<<hitCov(0,0)<<" cm"<<std::endl;
        hitCov.Print();
    }

    /// new space point and TrackPoint, unit in cm
    genfit::TrackPoint* trackPoint = new genfit::TrackPoint(
            new genfit::SpacepointMeasurement(pos_smeared,hitCov,
                cellID ,hitID, nullptr), m_track);
    if(m_debug>=2){
        std::cout<<m_name<<" add TrackPoint \n";
        trackPoint->Print();
    }
    m_track->insertPoint(trackPoint);
    return true;
}//end of addSpacePointMeasurement

/// Return isurface of a silicon hit
const dd4hep::rec::ISurface*
GenfitTrack::getISurface(edm4hep::ConstTrackerHit* hit){
    dd4hep::rec::SurfaceManager surfaceManager(*m_geomSvc->lcdd());

    std::string detectorName;
    unsigned long long cellID=hit->getCellID();
    int detTypeID=getDetTypeID(hit->getCellID());
    if(detTypeID==lcio::ILDDetID::VXD){
        detectorName="VXD";
        if(hit->getPosition().z>0){
            cellID+=0x100000000;
        }else{
            cellID+=0x300000000;
        }
    }else if(detTypeID==lcio::ILDDetID::SIT){
        detectorName="SIT";
    }else if(detTypeID==lcio::ILDDetID::SET){
        detectorName="SET";
    }else if(detTypeID==lcio::ILDDetID::FTD){
        detectorName="FTD";
    }else{
        std::cout << "ERROR:getISurface  iSurface = NULL!" << std::endl;
        return nullptr;
    }
    if(m_debug>2) std::cout<<m_name<<" detectorName  "<<detectorName
        <<" hit position mm "<<hit->getPosition()<<" hit.cellID "<<hit->getCellID()
            <<" cellId "<<cellID<<" detTypeID "<<detTypeID<<std::endl;
    const dd4hep::rec::SurfaceMap* surfaceMap=surfaceManager.map(detectorName);
    auto iter=surfaceMap->find(cellID);
    dd4hep::rec::ISurface* iSurface=nullptr;
    if(iter!=surfaceMap->end()){iSurface=(*iter).second;}

    //std::multimap< unsigned long, dd4hep::rec::ISurface*>::const_iterator it,itend;
    //it=surfaceMap->begin();
    //itend= surfaceMap->end();
    //std::cout<<" print map "<<detectorName<<std::endl;
    //for(; it!=itend; it++){
    //    dd4hep::rec::ISurface* surf = it->second;
    //    dd4hep::rec::Vector3D origin = surf->origin();
    //    std::cout <<"surf id "<< surf->id() << " origin xyz " << origin.x()
    //        << " " << origin.y() << " " << origin.z() << std::endl;
    //}
    return iSurface;
}

/// Add a 1d strip or 2d pixel smeared by sigma
    bool
GenfitTrack::addSiliconMeasurement(edm4hep::ConstTrackerHit* hit,
        float sigmaU,float sigmaV,int cellID,int hitID)
{
    if(m_debug>0) std::cout<<"addSiliconMeasurement "<<*hit<<std::endl;

    ///get surface by cellID
    const dd4hep::rec::ISurface* iSurface = getISurface(hit);
    if(nullptr==iSurface){
        std::cout<<m_name<<" addSiliconMeasurement get surface ERROR!"<<std::endl;
        return false;
    }
    TVector3 o,u,v;
    double lengthU,lengthV;
    getISurfaceOUV(iSurface,o,u,v,lengthU,lengthV);

    ///Get detector plane parameter
    //VXD
    //SIT
    //SET

    ///Get measurement and cov
    TVectorD hitCoords(2);
    TVector3 p;
    TMatrixDSym hitCov(2);
    getMeasurementAndCov(hit,p,hitCov);
    hitCoords(0)=(p-o).Dot(u);
    hitCoords(1)=(p-o).Dot(v);
    if(m_debug) std::cout<<"yzhang debug hitCoords cm "<<hitCoords(0)<<" "<<hitCoords(1)<<std::endl;
    hitCov.Zero();
    hitCov(0,0)=sigmaU*sigmaU;
    hitCov(1,1)=sigmaV*sigmaV;

    ///Create planer finite detector plane, measurement and TrackPoint
    genfit::RectangularFinitePlane* pixelOrStripPlane=
        new genfit::RectangularFinitePlane(
                -lengthU/2.,lengthU/2.,-lengthV/2.,lengthV/2.);
    genfit::SharedPlanePtr plane(new genfit::DetPlane(o,u,v));
    plane->setFinitePlane(pixelOrStripPlane);
    genfit::PlanarMeasurementSDT* planarMeasurement=new genfit::PlanarMeasurementSDT(
            hitCoords,hitCov,cellID,hitID,nullptr);
    planarMeasurement->setTrackerHit(hit);
    //genfit::PlanarMeasurement* planarMeasurement=new genfit::PlanarMeasurement(
    //        hitCoords,hitCov,cellID,hitID,nullptr);
    planarMeasurement->setPlane(plane);
    m_track->insertPoint(new genfit::TrackPoint(planarMeasurement,m_track));

    if(m_debug>2){
        std::cout<<"hitID "<<hitID<<" unit cm"<<std::endl;
        std::cout<<"sigmaU "<<sigmaU<<" sigmaV "<<sigmaV<<std::endl;
        std::cout<<"lengthU "<<lengthU<<" lengthV "<<lengthV<<std::endl;
        std::cout<<"u "<<u.x()<<" "<<u.y()<<" "<<u.z()<<std::endl;
        std::cout<<"v "<<v.x()<<" "<<v.y()<<" "<<v.z()<<std::endl;
        std::cout<<"o "<<o.x()<<" "<<o.y()<<" "<<o.z()<<std::endl;
        std::cout<<"p "<<p.x()<<" "<<p.y()<<" "<<p.z()<<std::endl;
        std::cout<<"hitCoords "<<hitCoords(0)<<" "<<hitCoords(1)<<std::endl;
        std::cout<<"hitCov "<<hitCov(0,0)<<" "<<hitCov(1,1)<<std::endl;
    }

    return true;
}

int GenfitTrack::addSiliconMeasurements(edm4hep::Track& track,
        std::vector<float> sigmaUVec,std::vector<float> sigmaVVec)
{
    ///Get TrackerHit on Track
    int nHitAdd=0;
    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
        edm4hep::ConstTrackerHit* trackerHit=
          new edm4hep::ConstTrackerHit(track.getTrackerHits(iHit));
        unsigned long long cellID=trackerHit->getCellID();
        float sigmaU,sigmaV;
        int sigmaUID=getSigmas(cellID,sigmaUVec,sigmaVVec,sigmaU,sigmaV);
        if(0==sigmaUID) continue;//skip DC track
        addSiliconMeasurement(trackerHit,sigmaU,sigmaV,cellID,nHitAdd++);
        GenfitHit* genfitHit=
            new GenfitHit(trackerHit,nullptr,nullptr,
            nullptr,0.,0.);
        if(nullptr==genfitHit) continue;
        m_genfitHitVec.push_back(genfitHit);
    }
    return nHitAdd;
}

//Add wire measurement on wire, unit conversion here
//int GenfitTrack::addWireMeasurementsFromList(podio::RelationRange<edm4hep::ConstTrackerHit> hits,float sigma,
int GenfitTrack::addWireMeasurementsFromList(std::vector<edm4hep::ConstTrackerHit*>& hits,float sigma,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        int sortMethod, bool truthAmbig,float skipCorner,float skipNear)
{
    if(m_debug>0){ std::cout<<"addWireMeasurementsFromList"<<std::endl; }
    //podio::RelationRange<edm4hep::ConstTrackerHit> hits_t=track.getTrackerHits();
    std::vector<edm4hep::ConstTrackerHit*> sortedTrackerHits;
//    sortedTrackerHits.reserve(100);
    getSortedTrackerHits(hits,assoHits,sortedTrackerHits,sortMethod);

    if(m_debug>0){
      std::cout<<"n sortedTrackerHits "<<sortedTrackerHits.size()<<std::endl;
    }
    int nHitAdd=0;
    for(auto trackerHit : sortedTrackerHits){
        edm4hep::ConstSimTrackerHit simTrackerHitAsso;
        getAssoSimTrackerHit(assoHits,trackerHit,simTrackerHitAsso);
        //FIXME driftVelocity
        GenfitHit* genfitHit=
            makeAGenfitHit(trackerHit,&simTrackerHitAsso,sigma,truthAmbig,
                    skipCorner,skipNear);
        if(nullptr==genfitHit) continue;
        m_genfitHitVec.push_back(genfitHit);

        ///New a TrackPoint,create connection between meas. and trackPoint
        WireMeasurementDC* wireMeasurementDC=
            new WireMeasurementDC(genfitHit,nHitAdd++);
        m_track->insertPoint(new genfit::TrackPoint(wireMeasurementDC,m_track));
        if(m_debug>=2){ std::cout<<nHitAdd-1; wireMeasurementDC->Print(); }
    }//end of loop over sotred hits
    return nHitAdd;
}//end of addWireMeasurementsFromList

//Add wire measurement on wire, unit conversion here
int GenfitTrack::addWireMeasurementsOnTrack(edm4hep::Track& track,float sigma,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        int sortMethod, bool truthAmbig,float skipCorner,float skipNear)
{
    if(m_debug>0){ std::cout<<"addWireMeasurementsOnTrack"<<std::endl; }

    if(m_debug>0){
      std::cout<<"n trackerHits size"<<track.trackerHits_size()<<std::endl;
    }
    int nHitAdd=0;
    if(0==track.trackerHits_size()){
        if(m_debug>0) std::cout<<"No hit on track"<<std::endl;
        return nHitAdd;
    }

    podio::RelationRange<edm4hep::ConstTrackerHit> hits_t=track.getTrackerHits();
    std::vector<edm4hep::ConstTrackerHit*> hits;
    //hits.reserve(1000);
    for(auto &h:hits_t){
        edm4hep::ConstTrackerHit* hit = const_cast<edm4hep::ConstTrackerHit*>(&h);
        hits.push_back(hit); 
    }
    nHitAdd=addWireMeasurementsFromList(hits,sigma,assoHits,sortMethod,truthAmbig,skipCorner,skipNear);

    return nHitAdd;
}//end of addWireMeasurements

/// Get MOP
bool GenfitTrack::getMOP(int hitID,genfit::MeasuredStateOnPlane& mop,
        genfit::AbsTrackRep* trackRep) const
{

    if(nullptr == trackRep) trackRep = getRep();
    try{
        mop = m_track->getFittedState(hitID,trackRep);
    }catch(genfit::Exception& e){
        e.what();
        return false;
    }
    return true;
}

/// New and add a track representation to track
genfit::RKTrackRep* GenfitTrack::addTrackRep(int pdgType,int charge)
{
    int chargeId=0;
    charge>0 ? chargeId=0 : chargeId=1;//s_PDG[0]: positive particle

    genfit::RKTrackRep* rep=new genfit::RKTrackRep(s_PDG[chargeId][pdgType]);
    m_track->addTrackRep(rep);
    return rep;
}

/// Get the position from genfit::Track::getStateSeed
const TLorentzVector GenfitTrack::getSeedStatePos()const
{
    TVectorD seedState(6);
    seedState = m_track->getStateSeed();
    TVector3 p(seedState[0],seedState[1],seedState[2]);
    p = p*dd4hep::cm;
    TLorentzVector pos(p[0],p[1],p[2],9999);//FIXME
    return pos;
}

/// Get the momentum from genfit::Track::getStateSeed
const TVector3 GenfitTrack::getSeedStateMom() const
{
    TVectorD seedState(6); seedState = m_track->getStateSeed();
    TVector3 mom(seedState[3],seedState[4],seedState[5]);
    return mom*dd4hep::GeV;
}

/// Get the seed states of momentum and position
void GenfitTrack::getSeedStateMom(TLorentzVector& pos, TVector3& mom) const
{
    TVectorD seedState(6); seedState = m_track->getStateSeed();
    mom = TVector3(seedState[3],seedState[4],seedState[5])*dd4hep::GeV;
    seedState = m_track->getStateSeed();
    TVector3 p = TVector3(seedState[0],seedState[1],seedState[2])*dd4hep::cm;
    pos.SetXYZT(p[0],p[1],p[2],9999);//FIXME time
}

unsigned int GenfitTrack::getNumPoints() const
{
    return m_track->getNumPoints();
}

GenfitHit* GenfitTrack::GetHit(long unsigned int i) const {
  if(i>=m_genfitHitVec.size())return nullptr;
  return m_genfitHitVec[i];
}

unsigned int GenfitTrack::getNumPointsDet(int detTypeID) const
{
    unsigned int nHit=0;
    const std::vector<genfit::TrackPoint*> tps=m_track->getPoints();
    for(auto tp:tps){
        const genfit::AbsMeasurement* m=nullptr;
        if(tp->hasRawMeasurements()) m=tp->getRawMeasurement();
        if(nullptr!=m && detTypeID==getDetTypeID(m->getDetId())) nHit++;
    }
    return nHit;
}

/// Test the fit result FIXME
bool GenfitTrack::fitSuccess(int repID) const
{

    genfit::FitStatus* fitStatus = m_track->getFitStatus(getRep(repID));

    /// Test fitting converged
    if (!fitStatus->isFitted()||!fitStatus->isFitConverged()
            ||fitStatus->isFitConvergedFully()) {
        if(m_debug>=2)std::cout<<m_name<< "Fitting is failed... isFitted"
            <<fitStatus->isFitted()<<" , isFitConverged "
                <<fitStatus->isFitConverged()<<", isFitConvergedFully "
                <<fitStatus->isFitConvergedFully()<<std::endl;
        return false;
    }

    double chi2 = fitStatus->getChi2();
    double ndf  = fitStatus->getNdf();
    if(m_debug>=2)std::cout<< "Fit result: chi2 "<<chi2 <<" ndf "<<ndf
        << " chi2/ndf = " << chi2/ndf<<std::endl;

    /// Test fitting chi2
    if (chi2<= 0) {
        if(m_debug>=2)std::cout<<m_name<< "Fit chi2<0 (chi2,ndf) = (" <<
            chi2 << "," << ndf  << ")"<<std::endl;
        return false;
    }
    return true;
}

void GenfitTrack::setDebugGenfit(int debug)
{
    if(m_track){
        for(unsigned int i=0;i<m_track->getNumReps();i++){
            m_track->getTrackRep(i)->setDebugLvl(debug);
        }
    }
}
void GenfitTrack::setDebugLocal(int debug){
#ifdef GENFIT_MY_DEBUG
    if(m_track){
        for(unsigned int i=0;i<m_track->getNumReps();i++){
            //m_track->getTrackRep(i)->setDebugLvlLocal(debug);
        }
    }
    m_debugLocal=debug;
#endif
    if(m_debug)std::cout<<"GenfitTrack:setDebugLvlLocal "<<debug<<std::endl;
}

void GenfitTrack::printSeed() const
{
    TLorentzVector pos = getSeedStatePos();
    TVector3 mom = getSeedStateMom();
    std::cout << " Seed pos "  << std::endl;
    pos.Print();
    std::cout << " Seed Mom =  " << std::endl;
    mom.Print();
    print(pos,mom);
    TMatrixDSym covSeed=m_track->getCovSeed();
    //if(m_debug>1)
    std::cout << " covSeed = " << std::endl;
    covSeed.Print();
    //std::cout<<" pdg "<<0<<getRep(0)<<std::endl;//FIXME
}

void GenfitTrack::printFitted(int repID) const
{
    TLorentzVector fittedPos;
    TVector3 fittedMom;
    TMatrixDSym cov;

    if(m_debug>=2)std::cout<<m_name<< "printFitted nHit="
        <<m_track->getNumPoints()<<std::endl;
    for(unsigned int iHit=0; iHit<m_track->getNumPoints(); iHit++){
        if (getPosMomCovMOP((int) iHit, fittedPos, fittedMom, cov, repID)){
            //print(fittedPos,fittedMom,to_string(iHit).c_str());//TODO
        }else{
            if(m_debug>=2)std::cout<<m_name<<"Hit "<<iHit
                <<" have no fitter info"<<std::endl;
        }
    }
}

/// Print track information
void GenfitTrack::print( TLorentzVector pos, TVector3 mom,
        const char* comment) const
{
    TVector3 pglo = pos.Vect();
    TVector3 mglo = mom;

    // TODO
    if(m_debug>=2)std::cout<<m_name<<" "<<comment<<std::endl;

    if(m_debug>1){
        for(unsigned int i=0;i<m_track->getNumReps();i++){
            m_track->getTrackRep(i)->Print();
        }
    }
    //for(unsigned int i=0; i<m_track->getNumPoints(); i++){
    //  m_track->getPoint(i)->print();
    //}
}

/// Get position, momentum, cov on plane of hitID-th hit
bool GenfitTrack::getPosMomCovMOP(int hitID, TLorentzVector& pos,
        TVector3& mom, TMatrixDSym& cov,int repID) const
{
    TVector3 p;
    genfit::MeasuredStateOnPlane mop;
    if(!getMOP(hitID,mop,getRep(repID))) return false;
    mop.getPosMomCov(p,mom,cov);
    pos.SetVect(p*dd4hep::cm);
    pos.SetT(9999);//FIXME
    mom = mom*(dd4hep::GeV);
    //FIXME
    //TrackingUtils::CovConvertUnit(cov, dd4hep::cm, dd4hep::GeV);
    return true;
}

int GenfitTrack::getNumPointsWithFittedInfo(int repID) const
{
    int nHitWithFittedInfo = 0;
    //check number of hit with fitter info
    int nHit = m_track->getNumPointsWithMeasurement();
    //check number of hit with fitter info
    for(int i=0; i<nHit; i++){
        if(nullptr != m_track->getPointWithFitterInfo(i,getRep(repID))){
            nHitWithFittedInfo++;
        }
    }
    return nHitWithFittedInfo;
}

int GenfitTrack::getFittedState(TLorentzVector& pos, TVector3& mom,
        TMatrixDSym& cov, int trackPointId, int repID, bool biased) const
{
    //check number of hit with fitter info
    if(getNumPointsWithFittedInfo(repID)<=2) return 1;

    //get track rep
    genfit::AbsTrackRep* rep = getRep(repID);
    if(nullptr == rep) return 2;

    //get first or last measured state on plane
    genfit::MeasuredStateOnPlane mop;
    try{
        mop = m_track->getFittedState(trackPointId,rep,biased);
    }catch(genfit::Exception& e){
        std::cout<<" getNumPointsWithFittedInfo="
            <<getNumPointsWithFittedInfo(repID)
            <<" no TrackPoint with fitted info "<<std::endl;
        if(m_debug>=2)std::cout<<m_name
            <<"Exception in getFittedState"<<std::endl;
        return 3;
    }

    //get state
    TVector3 p;
    mop.getPosMomCov(p,mom,cov);
    pos.SetVect(p*dd4hep::cm);
    pos.SetT(9999);//FIXME
    mom = mom*(dd4hep::GeV);

    return 0;//success
}

// Get point with fitter info
int GenfitTrack::getDetIDWithFitterInfo(int hitID, int idRaw) const
{
    return m_track->getPointWithFitterInfo(hitID)->
        getRawMeasurement(idRaw)->getDetId();
}

int GenfitTrack::getPDG(int id) const
{
    return m_track->getTrackRep(id)->getPDG();
}

int GenfitTrack::getPDGCharge(int id) const
{
    return m_track->getTrackRep(id)->getPDGCharge();
}

const genfit::FitStatus*
GenfitTrack::getFitStatus(int repID) const
{
    return m_track->getFitStatus(getRep(repID));
}

/// Extrapolate track to the cloest point of approach(POCA) to the wire of hit,
/// return doca, poca on wire and track
/// outputs poca [cm]
/// unit conversion here
double GenfitTrack::extrapolateToHit(TVector3& poca, TVector3& pocaDir,
        TVector3& pocaOnWire, double& doca,edm4hep::MCParticle mcParticle,
        int cellID, int repID, bool stopAtBoundary, bool calcJacobianNoise)const
{
    ///Get MCParticle and convert unit
    TVector3 pos;
    TVector3 mom;
    getPosMomFromMCPartile(mcParticle,pos,mom);

    TVector3 end0(0,0,0);
    TVector3 end1(0,0,0);
    getEndPointsOfWire(cellID,end0,end1);

    int chargeId;
    mcParticle.getCharge() >0 ? chargeId=0 : chargeId=1;//s_PDG[0]: positive particle
    genfit::RKTrackRep* rep = new genfit::RKTrackRep(s_PDG[chargeId][repID]);
    //genfit::MeasuredStateOnPlane state(rep);
    genfit::StateOnPlane state(rep);
    rep->setPosMom(state, pos, mom);

    if(m_debug){
        std::cout<<"GenfitterTrack::extrapolateToHit before pos and mom "
            <<" charge "<<(int)mcParticle.getCharge()<<" repID "<<repID
            <<" pdg "<<s_PDG[chargeId][repID]<<std::endl;
        pos.Print();
        mom.Print();
        std::cout<<" end point"<<std::endl;
        end0.Print();
        end1.Print();
        std::cout<<" state "<<std::endl;
        state.Print();
        std::cout<<" stopAtBoundary "<<stopAtBoundary<<std::endl;
        std::cout<<" calcJacobianNoise "<<calcJacobianNoise<<std::endl;
    }


    double extrapoLen=0.;
    try {
        extrapoLen = rep->extrapolateToLine(state,end0,end0-end1,
                poca,pocaDir,pocaOnWire,stopAtBoundary,calcJacobianNoise);
    }catch(genfit::Exception& e) {
        if(m_debug>=3)std::cout<<
            "Exception in GenfitterTrack::extrapolateToHit"<<e.what()<<std::endl;
        return extrapoLen/GenfitUnit::mm*dd4hep::mm;
    }

    doca = (pocaOnWire-poca).Mag();
    if(m_debug>=2){
        std::cout<< " poca "<<poca.x()<<","<<poca.y()<<" "<<poca.z()
            <<" pocaDir "<<pocaDir.x()<<","<<pocaDir.y()<<" "<<pocaDir.z()
            <<" pocaOnWire "<<pocaOnWire.x()<<","<<pocaOnWire.y()<<" "<<pocaOnWire.z()
            <<" doca "<<doca<<" cm"<<std::endl;
    }

    delete rep;
    return extrapoLen/GenfitUnit::mm*dd4hep::mm;
}//end of extrapolateToHit

/////Add space point measurement from edm4hep::Track to genfit track
//int GenfitTrack::addHitsOnEdm4HepTrack(const edm4hep::Track& track,
//        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
//        std::vector<float> sigma,bool smear,
//        bool measurementTypeSi, bool measurementTypeDC){
//    ///Get TrackerHit on Track
//    int hitID=0;
//    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
//        edm4hep::ConstTrackerHit hit=track.getTrackerHits(iHit);
//        ///Get hit type
//        int detTypeID=getDetTypeID(hit.getCellID());
//        if(m_debug>=2)std::cout<<"addHitsOnEdm4HepTrack "<<iHit<<" hit "<<hit
//            <<" detTypeID "<<detTypeID<<" type "<<hit.getType()<<std::endl;
//
//        bool hitIsSpapcePoint=UTIL::BitSet32(hit.getType())[
//                              UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT];
//        bool hitIsPlanar=UTIL::BitSet32(hit.getType())[
//                         UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL];
//        if(m_debug>2){
//            std::cout<<detTypeID<<" COMPOSITE_SPACEPOINT "<<hitIsSpapcePoint
//                <<std::endl;
//            std::cout<<detTypeID<<" ONE_DIMENSIONAL "<<hitIsPlanar<<std::endl;
//        }
//
//    }
//
//    return 1;
//}

///Add space point measurement of silicon from edm4hep::Track to genfit track
int GenfitTrack::addSpacePointsSi(const edm4hep::Track& track,
        std::vector<float> sigmaU,std::vector<float> sigmaV)
{
    ///Get TrackerHit on Track
    int nHitAdd=0;
    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
        int detTypeID=getDetTypeID(track.getTrackerHits(iHit).getCellID());
        if(m_geomSvc->lcdd()->constant<int>("DetID_DC")==detTypeID) continue;
        edm4hep::ConstTrackerHit hit=track.getTrackerHits(iHit);
        edm4hep::Vector3d pos=hit.getPosition();

        //edm4hep::Vector3d pos=simTrackerHitAsso.getPosition();

        TVector3 p(pos.x,pos.y,pos.z);

        p.Print();
        unsigned long long cellID = hit.getCellID();
        if(addSpacePointMeasurement(p,sigmaU,sigmaV,cellID,nHitAdd)){
            if(m_debug>=2)std::cout<<"add silicon space point"<<std::endl;
            nHitAdd++;
        }else{
            if(m_debug>=2)std::cout<<"addSpacePointMeasurement"
                <<cellID<<" faieled" <<std::endl;
        }
    }
    if(m_debug>=2) std::cout<<m_name<<" Si addHitAdd="<<nHitAdd<<std::endl;
    return nHitAdd;
}//end of addSpacePointsSi

///Add drift chamber space point from edm4hep::track
int GenfitTrack::addSpacePointsDC(const edm4hep::Track& track,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        std::vector<float> sigmaU,std::vector<float> sigmaV)
{
    if(m_debug>=2){
    std::cout<<m_name<<" 5. addSpacePointsDC nTrackerHit="
        <<track.trackerHits_size()<<std::endl;
    }

    ///Get TrackerHit with min. time in Cell
    std::vector<edm4hep::ConstSimTrackerHit> sortedDCTrackHitCol;
    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
        edm4hep::ConstTrackerHit hit=track.getTrackerHits(iHit);
        int detTypeID=getDetTypeID(track.getTrackerHits(iHit).getCellID());
        if(m_geomSvc->lcdd()->constant<int>("DetID_DC")!=detTypeID) continue;

        ///Select the SimTrakerHit with least time
        float minTime=FLT_MAX;
        edm4hep::ConstSimTrackerHit minTimeSimHit;
        for(int iSimHit=0;iSimHit<(int) assoHits->size();iSimHit++){
            if(assoHits->at(iSimHit).getRec()==hit &&
                    assoHits->at(iSimHit).getSim().getTime()<minTime){
                minTimeSimHit=assoHits->at(iSimHit).getSim();
                minTime=assoHits->at(iSimHit).getSim().getTime();
            }
        }
        if(!minTimeSimHit.isProducedBySecondary()){
            sortedDCTrackHitCol.push_back(minTimeSimHit);
        }
    }//end loop over TrackerHit

    ///Sort DC SimTrackerHit by time
    std::sort(sortedDCTrackHitCol.begin(),sortedDCTrackHitCol.end(),sortDCSimHit);

    ///Add DC hits to track
    int nHitAdd=0;
    for(auto dCTrackerHit: sortedDCTrackHitCol){
        edm4hep::Vector3d pos=dCTrackerHit.getPosition();
        TVector3 p(pos.x,pos.y,pos.z);

        if(m_debug) std::cout<<"6. addSpacePointsDC hit "<<dCTrackerHit<<std::endl;
        unsigned long long cellID=dCTrackerHit.getCellID();
        if(addSpacePointMeasurement(p,sigmaU,sigmaV,dCTrackerHit.getCellID()
                    ,nHitAdd)){
            if(m_debug>=2)std::cout<<"add DC space point"<<std::endl;
            nHitAdd++;
        }else{
            if(m_debug>=2)std::cout<<"addSpacePointMeasurement"
                <<cellID<<" faieled" <<std::endl;
        }
    }

    if(m_debug>=2) std::cout<<m_name<<" DC addHitAdd="<<nHitAdd<<std::endl;
    return nHitAdd;
}//end of addSpacePointsDC

//double GenfitTrack::extrapolateToPoint(TVector3& pos, TVector3& mom,
//        TMatrixDSym& cov,
//        const TVector3& point,
//        int repID,// same with pidType
//        bool stopAtBoundary,
//        bool calcJacobianNoise) const
//{
//    return extrapolateToPoint(pos,mom,cov,point,repID,stopAtBoundary,
//            calcJacobianNoise);
//
//}//end of extrapolateToPoint

double GenfitTrack::extrapolateToPoint(TVector3& pos, TVector3& mom,
        TMatrixDSym& cov, const TVector3& point,
        int repID,// same with pidType
        bool stopAtBoundary,
        bool calcJacobianNoise) const
{
    double trackLength(1e9*dd4hep::cm);
    if(!getFitStatus(repID)->isFitted()) return trackLength;
    try{
        // get track rep
        genfit::AbsTrackRep* rep = getRep(repID);
        if(nullptr == rep) {
            if(m_debug>=2)std::cout<<"In extrapolateToPoint rep "
                <<repID<<" not exist!"<<std::endl;
            return trackLength*dd4hep::cm;
        }

        /// extrapolate to point
        //genfit::StateOnPlane state(*(&(track->getTrack()->getFittedState(0,rep))));

        // get track point with fitter info
        genfit::TrackPoint* tp = getTrack()->getPointWithFitterInfo(0,rep);
        if(nullptr == tp) {
            if(m_debug>=2)std::cout<<
                "In extrapolateToPoint TrackPoint is null"<<std::endl;
            return trackLength*dd4hep::cm;//FIXME unit
        }

        // get fitted state on plane of this track point
        genfit::KalmanFittedStateOnPlane* state =
            static_cast<genfit::KalmanFitterInfo*>(
                    tp->getFitterInfo(rep))->getBackwardUpdate();
        genfit::StateOnPlane orignalState(*state);
        if(m_debug>3){
            tp->Print();
            std::cout<<" original state before extrapolate "<<std::endl;
            state->Print();
        }

        if(nullptr == state) {
            if(m_debug>=2)std::cout<<
                "In extrapolateToPoint KalmanFittedStateOnPlane is null"<<std::endl;
            return trackLength*dd4hep::cm;//FIXME unit
        }
        //rep->setDebugLvl(10);
        trackLength = rep->extrapolateToPoint(*state,
                point*(1/dd4hep::cm),stopAtBoundary, calcJacobianNoise);

        rep->getPosMomCov(*state,pos,mom,cov);//FIXME exception exist

        pos = pos*dd4hep::cm;//FIXME unit
        mom = mom*dd4hep::GeV;//FIXME unit
        if(m_debug>3){
            std::cout<<" original state before extrapolate "<<std::endl;
            orignalState.Print();
            std::cout<<" extrapolated state "<<std::endl;
            state->Print();
        }

    } catch(genfit::Exception& e){
        if(m_debug>=3)std::cout
            <<"Exception in GenfitTrack::extrapolateToPoint"
                << e.what()<<std::endl;
        trackLength = 1e9*dd4hep::cm;
    }
    return trackLength*dd4hep::cm;
}//end of extrapolateToPoint

/// Extrapolate the track to the cyliner at fixed raidus
/// position & momentum as starting point
/// position and momentum at global coordinate in dd4hepUnit
/// return trackLength in dd4hepUnit
    double
GenfitTrack::extrapolateToCylinder(TVector3& pos, TVector3& mom,
        double radius, const TVector3 linePoint,
        const TVector3 lineDirection, int hitID, int repID,
        bool stopAtBoundary, bool calcJacobianNoise)
{
    double trackLength(1e9*dd4hep::cm);//FIXME unit
    if(!getFitStatus(repID)->isFitted()) return trackLength;
    try{
        // get track rep
        genfit::AbsTrackRep* rep = getRep(repID);
        if(nullptr == rep) {
            if(m_debug>=2)std::cout<<"In extrapolateToCylinder rep is null"
                <<std::endl;
            return trackLength*dd4hep::cm;//FIXME unit
        }

        // get track point with fitter info
        genfit::TrackPoint* tp =
            getTrack()->getPointWithFitterInfo(hitID,rep);
        if(nullptr == tp) {
            if(m_debug>=2)std::cout<<
                "In extrapolateToCylinder TrackPoint is null"<<std::endl;
            return trackLength*dd4hep::cm;//FIXME unit
        }

        // get fitted state on plane of this track point
        genfit::KalmanFittedStateOnPlane* state =
            static_cast<genfit::KalmanFitterInfo*>(
                    tp->getFitterInfo(rep))->getBackwardUpdate();

        if(nullptr == state){
            if(m_debug>=2)std::cout<<"In extrapolateToCylinder "
                <<"no KalmanFittedStateOnPlane in backwardUpdate"<<std::endl;
            return trackLength*dd4hep::cm;
        }
        rep->setPosMom(*state, pos*(1/dd4hep::cm), mom*(1/dd4hep::GeV));//FIXME unit

        /// extrapolate
        trackLength = rep->extrapolateToCylinder(*state,
                radius/dd4hep::cm, linePoint*(1/dd4hep::cm), lineDirection,
                stopAtBoundary, calcJacobianNoise);
        // get pos&mom at extrapolated point on the cylinder
        rep->getPosMom(*state,pos,mom);//FIXME exception exist
        pos = pos*dd4hep::cm;
        mom = mom*dd4hep::GeV;
    } catch(genfit::Exception& e){
        if(m_debug>=3)std::cout
            <<"Exception in GenfitTrack::extrapolateToCylinder "
                << e.what()<<std::endl;
        trackLength = 1e9*dd4hep::cm;
    }
    return trackLength*dd4hep::cm;
}

bool GenfitTrack::debugDistance(const edm4hep::TrackerHitCollection* dCDigiCol,
        int& nFittedDC,int& nFittedSDT,int& ngenfitHit,
        std::vector<double>& smearDistance,
        std::vector<double>& truthDistance,double driftVelocity)
{

    int DCHit = 0;
    int SDTHit = 0;
    unsigned int nPoints = m_track->getNumPoints();
    for(unsigned int i = 0; i<nPoints; i++){
      genfit::TrackPoint* point = m_track->getPoint(i);
      genfit::AbsMeasurement* absMea = point->getRawMeasurement();
      genfit::PlanarMeasurementSDT* sdtMea =
          dynamic_cast<genfit::PlanarMeasurementSDT*>(absMea);
      if(sdtMea){
        const edm4hep::ConstTrackerHit* TrackerHit_ = sdtMea->getTrackerHit();
        SDTHit++;
      }else{
        WireMeasurementDC* dcMea =
            dynamic_cast<WireMeasurementDC*>(absMea);
        if(dcMea){
          const edm4hep::ConstTrackerHit* TrackerHit_ = dcMea->getTrackerHit();
          smearDistance.push_back(1e-3*driftVelocity*TrackerHit_->getTime());
          DCHit++;
          for(auto dcDigi: *dCDigiCol){
            if(dcDigi.getCellID() == TrackerHit_->getCellID())
            {
                truthDistance.push_back(1e-3*driftVelocity*(dcDigi.getTime()));
            }
          }
        }
      }
    }

    DCHit = nFittedDC;
    SDTHit = nFittedSDT;
    ngenfitHit = nFittedDC+nFittedDC;

    return true;
}

/// Get doca on plane of hitID-th hit
bool GenfitTrack::GetDocaRes(int hitID, double& DriftDis,double& fittedDoca,
        double& res, int repID, bool biased) const
{
    genfit::TrackPoint* trackPoint = m_track->getPointWithFitterInfo(hitID, getRep(repID));
    if(!trackPoint)return false;
    genfit::AbsMeasurement* thismea = trackPoint->getRawMeasurement(0);
    const TVectorD measereddoca = thismea->getRawHitCoords();
    DriftDis = abs(measereddoca[0]);
    genfit::MeasuredStateOnPlane mop;
    if(!getMOP(hitID,mop,getRep(repID))) return false;
    TVectorD state = mop.getState();
    if(state.GetNrows() != 5)return false;
    fittedDoca = abs(state[3]);
    res = 10*(fittedDoca - DriftDis);
    return true;
}

bool GenfitTrack::storeTrack(edm4hep::ReconstructedParticle& recParticle,
        edm4hep::Track& track,TVector3& pocaToOrigin_pos,
        TVector3& pocaToOrigin_mom,TMatrixDSym& pocaToOrigin_cov,
        // edm4hep::TrackState& trackState,
        int pidType, int ndfCut, double chi2Cut,
        int& nFittedDC, int& nFittedSDT, int& ngenfitHit,
        std::vector<double>& trackL, std::vector<double>& hitMom,
        std::vector<float>& truthMomEdep,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        std::vector<double>& driftDis,
        std::vector<double>& FittedDoca,
        std::vector<double>& Res)
{

    int id = 0;
    int fitid = 0;
    int DCHit = 0;
    int dcFit = 0;
    int sdtFit = 0;
    int SDTHit = 0;

    //  getNumRawMeasurements
    unsigned int nPoints = m_track->getNumPoints();
    std::cout << " nPoints =  " << nPoints << std::endl;

    std::vector<double> hitMomMag;

    while(1){
        genfit::TrackPoint* point = m_track->getPointWithFitterInfo(id);
        if(!point)break;
        id++;
        genfit::AbsMeasurement* absMea = point->getRawMeasurement();

        // getPoint FitterInfo
        genfit::AbsFitterInfo* FitterInfo = point->getFitterInfo();
        genfit::KalmanFitterInfo *KalmanfitterInfo =
            dynamic_cast<genfit::KalmanFitterInfo*>(FitterInfo);

        unsigned int nNumMea = KalmanfitterInfo->getNumMeasurements();

        bool flag = false;
        for(int i=0;i<nNumMea;i++)
        {
            genfit::MeasurementOnPlane* MeaOnPlane =
                KalmanfitterInfo->getMeasurementOnPlane(i);
            double weight = MeaOnPlane->getWeight();
            if(weight>0.8) flag = true;
        }
        if(flag) fitid++;

        genfit::MeasurementOnPlane * Mea =
            KalmanfitterInfo->getMeasurementOnPlane();

        genfit::PlanarMeasurementSDT* sdtMea =
            dynamic_cast<genfit::PlanarMeasurementSDT*>(absMea);
        if(sdtMea){
            SDTHit++;
            if(flag) sdtFit++;
            const edm4hep::ConstTrackerHit* TrackerHit_ = sdtMea->getTrackerHit();
            track.addToTrackerHits(*TrackerHit_);
        }else{
            WireMeasurementDC* dcMea =
                dynamic_cast<WireMeasurementDC*>(absMea);
            if(dcMea){
                const edm4hep::ConstTrackerHit* TrackerHit_ = dcMea->getTrackerHit();
                track.addToTrackerHits(*TrackerHit_);


                DCHit++;
                if(flag)
                {

                    truthMomEdep.push_back(TrackerHit_->getEDep());
                    double DriftDis,fittedDoca,res = 0;
                    GetDocaRes(dcFit,DriftDis,fittedDoca,res);
                    driftDis.push_back(DriftDis);
                    FittedDoca.push_back(fittedDoca);
                    Res.push_back(res);

                    //                    TMatrixDSym fittedHitCov(6);//cm, GeV
                    //                    TLorentzVector fittedHitPos;
                    //                    TVector3 fittedHitMom;
                    //                    int fittedState=getFittedState(fittedHitPos,fittedHitMom,fittedHitCov,dcFit);
                    //                    hitMomMag.push_back(fittedHitMom.Mag());
                    //
                    //                    for(int iSim=0;iSim<assoHits->size();iSim++)
                    //                    {
                    //                        //if(*TrackerHit_ == assoHits->at(iSim).getRec())
                    //                        if(TrackerHit_->getCellID() == assoHits->at(iSim).getRec().getCellID())
                    //                        {
                    //                           // std::cout << " if  TrackerHit_ = " << TrackerHit_->getCellID() << std::endl;
                    //                           // std::cout << " if assoHits->at(iSim).getRec() = " << assoHits->at(iSim).getRec().getCellID() << std::endl;
                    //                           // std::cout << " Sim Mom = " << assoHits->at(iSim).getSim().getMomentum() << std::endl;
                    //                           // std::cout << " Sim MomMag = " << sqrt(assoHits->at(iSim).getSim().getMomentum()[0]*assoHits->at(iSim).getSim().getMomentum()[0]+assoHits->at(iSim).getSim().getMomentum()[1]*assoHits->at(iSim).getSim().getMomentum()[1]+assoHits->at(iSim).getSim().getMomentum()[2]*assoHits->at(iSim).getSim().getMomentum()[2]) << std::endl;
                    //
                    //                            break;
                    //                        }
                    //                    }

                    //std::cout << " i = " << dcFit << "fittedMom = " << fittedHitMom.X() << " " << fittedHitMom.Y() << " " << fittedHitMom.Z() << std::endl;
                    //std::cout << " i = " << dcFit << "fittedMomMag = " << fittedHitMom.Mag() << std::endl;
                    dcFit++;
                }
            }
        }
    }

    std::cout << " id = " << id << std::endl;
    std::cout << " fitid = " << fitid << std::endl;

    nFittedDC = dcFit;
    nFittedSDT = SDTHit;

    std::cout<<"nDC: "<<nFittedDC<<", nSDT: "<<nFittedSDT<<std::endl;
    std::cout<<"nFittedDC: "<<dcFit<<", nFittedSDT: "<<sdtFit<<std::endl;
    if(m_debug>0)std::cout<<m_name<<" store track ndfCut "<<ndfCut<<" chi2Cut "
        <<chi2Cut<<std::endl;

    /// Get fit status
    const genfit::FitStatus* fitState = getFitStatus();
    double ndf = fitState->getNdf();
    if(ndf>ndfCut){
        if(m_debug>0){
            std::cout<<m_name<<" cut by ndf="<<ndf<<">"<<ndfCut<<std::endl;
        }
        return false;
    }
    double chi2 = fitState->getChi2();
    if(chi2>chi2Cut){
        if(m_debug>0){
            std::cout<<m_name<<" cut by chi2="<<chi2<<">"<<chi2Cut<<std::endl;
        }
        return false;
    }
    double charge = fitState->getCharge();
    int isFitted = fitState->isFitted();
    int isConverged = fitState->isFitConverged();
    int isConvergedFully = fitState->isFitConvergedFully();

    TMatrixDSym fittedCov(6);//cm, GeV
    TLorentzVector fittedPos;
    TVector3 fittedMom;
    int fittedState=getFittedState(fittedPos,fittedMom,fittedCov);

    //get dx
    TVector3 pos;
    TVector3 mom;
    pos.SetXYZ(fittedPos.X(),fittedPos.Y(),fittedPos.Z());
    mom.SetXYZ(fittedMom.X(),fittedMom.Y(),fittedMom.Z());
    //std::cout << "fitted momx = " << mom.X() << "momy = " << mom.Y() << "momz = " << mom.Z() << " mom.Mag = " << mom.Mag() << std::endl;
    //std::cout << "fitted posx = " << pos.X() << "posy = " << pos.Y() << "posz = " << pos.Z() << std::endl;
    double radius = 80.022;  // cm
    const TVector3 linePoint(0,0,0);
    const TVector3 lineDirection(0,0,1);

    for(int i = 0;i<fitid;i++){
        int repID=1;
        radius+=1.741321;
        double tracklength =
            extrapolateToCylinder(pos,mom,radius,linePoint,lineDirection,repID);
        hitMomMag.push_back(mom.Mag());
        //std::cout << "momx = " << mom.X() << "momy = " << mom.Y() << "momz = " << mom.Z() << " mom.Mag = " << mom.Mag() << std::endl;
        //std::cout << "posx = " << pos.X() << "posy = " << pos.Y() << "posz = " << pos.Z() << std::endl;
        trackL.push_back(tracklength);
    }

    for(int j=0;j<hitMomMag.size()-1;j++)
    {
        hitMom.push_back(hitMomMag[j]-hitMomMag[j+1]);
        //std::cout << " Error Dep = " << hitMomMag[j]-hitMomMag[j+1] << std::endl;
    }

    if(m_debug>0)std::cout<<m_name<<" fit result: get status OK? pidType "
        <<pidType<<" fittedState"<<fittedState<<"==0? "<<(0==fittedState)
            <<" isFitted "<<isFitted
            <<" isConverged "<<isConverged<<" ndf "<<ndf<<std::endl;
    if((0!=fittedState)||(!isFitted)||(!isConvergedFully)||(ndf>ndfCut)){
        if(m_debug>0)std::cout<<m_name<<" fitting FAILED!=========="<<std::endl;
    }else{
        if(m_debug>0){
            const TLorentzVector seedPos=getSeedStatePos();
            const TVector3 seedMom=getSeedStateMom();
            std::cout<<m_name<<"===fit-seed:delta pos("
                <<fittedPos.X()-seedPos.X()<<","
                <<fittedPos.Y()-seedPos.Y()<<","
                <<fittedPos.Z()-seedPos.Z()<<") delta mom("
                <<fittedMom.X()-seedMom.X()<<","
                <<fittedMom.Y()-seedMom.Y()<<","
                <<fittedMom.Z()-seedMom.Z()<<")"<<" delta mag "<<fittedMom.Mag()-seedMom.Mag()<<std::endl;
            std::cout<<m_name<<" seed pos("
                <<seedPos.X()<<","
                <<seedPos.Y()<<","
                <<seedPos.Z()<<") seed mom("
                <<seedMom.X()<<","
                <<seedMom.Y()<<","
                <<seedMom.Z()<<std::endl;
            std::cout<<m_name<<" fit result: Pos("<<
                fittedPos.X()<<" "<<
                fittedPos.Y()<<" "<<
                fittedPos.Z()<<") mom("<<
                fittedMom.X()<<" "<<
                fittedMom.Y()<<" "<<
                fittedMom.Z()<<") p_tot "<<
                fittedMom.Mag()<<" pt "<<
                fittedMom.Perp()<<
                " chi2 "<<chi2<<
                " ndf "<<ndf
                <<std::endl;
            std::cout<<"fittedCov "<<std::endl;
            fittedCov.Print();
        }
    }


    //TVectorD seedState5(5);
    //getSeedState5(seedState5);

    ///track status at POCA to referencePoint: origin
    const TVector3 referencePoint(0,0,0);
    //TVector3 
    pocaToOrigin_pos.SetXYZ(1e9*dd4hep::cm,1e9*dd4hep::cm,1e9*dd4hep::cm);
    //TVector3 
    pocaToOrigin_mom.SetXYZ(1e9*dd4hep::GeV,1e9*dd4hep::GeV,1e9*dd4hep::GeV);
    //TMatrixDSym pocaToOrigin_cov;
    if(extrapolateToPoint(pocaToOrigin_pos,pocaToOrigin_mom,pocaToOrigin_cov,
                referencePoint) > 1e6*dd4hep::cm){
        if(m_debug>0)std::cout<<m_name<<" extrapolate to origin failed"<<std::endl;
        return false;
    }

    //Unit conversion of position and momentum cm->mm
    pocaToOrigin_pos.SetXYZ(pocaToOrigin_pos.X()/dd4hep::mm,
            pocaToOrigin_pos.Y()/dd4hep::mm,
            pocaToOrigin_pos.Z()/dd4hep::mm);
    pocaToOrigin_mom.SetXYZ(pocaToOrigin_mom.X(),
            pocaToOrigin_mom.Y(),
            pocaToOrigin_mom.Z());

    //unit conversion of error matrix
    TMatrixDSym covMatrix_6=fittedCov;
    for(int i=0;i<5;i++){
        covMatrix_6[0][i]=fittedCov[0][i]/dd4hep::mm;//d0 column
        covMatrix_6[2][i]=fittedCov[1][i]/dd4hep::mm;//omega column
        covMatrix_6[3][i]=fittedCov[2][i]/dd4hep::mm;//z0 column
        covMatrix_6[i][0]=fittedCov[i][0]/dd4hep::mm;//d0 row
        covMatrix_6[i][2]=fittedCov[i][1]/dd4hep::mm;//omega row
        covMatrix_6[i][3]=fittedCov[i][2]/dd4hep::mm;//z0 row
        //covMatrix_6[0][i]=fittedCov[0][i]/GenfitUnit::cm*dd4hep::cm;//d0 column
        //covMatrix_6[2][i]=fittedCov[2][i]*GenfitUnit::cm/dd4hep::cm;//omega column
        //covMatrix_6[3][i]=fittedCov[3][i]/GenfitUnit::cm*dd4hep::cm;//z0 column
        //covMatrix_6[i][0]=fittedCov[i][0]/GenfitUnit::cm*dd4hep::cm;//d0 row
        //covMatrix_6[i][2]=fittedCov[i][2]*GenfitUnit::cm/dd4hep::cm;//omega row
        //covMatrix_6[i][3]=fittedCov[i][3]/GenfitUnit::cm*dd4hep::cm;//z0 row
    }

    double Bz=m_genfitField->getBz(referencePoint)/GenfitUnit::tesla;

    if(m_debug>0){
        std::cout<<m_name<<" fit result poca: Pos"<<std::endl;
        pocaToOrigin_pos.Print();
        pocaToOrigin_mom.Print();
        std::cout<<" chi2 "<<chi2<< " ndf "<<ndf <<std::endl;
        std::cout<<"fittedCov "<<std::endl;
        fittedCov.Print();
        std::cout<<"covMatrix_6 "<<std::endl;
        covMatrix_6.Print();
        std::cout<<__LINE__<<" debug Bz tesla"<<Bz<<std::endl;
    }

    edm4hep::TrackState trackState;
    CEPC::getTrackStateFromPosMom(trackState,Bz,pocaToOrigin_pos,
            pocaToOrigin_mom,charge,covMatrix_6);

    trackState.location=0;//edm4hep::AtIP;//FIXME
    if(m_debug>2){
        std::cout<<m_name<<" trackState "<<trackState<<std::endl;
    }
    track.addToTrackStates(trackState);

    //ngenfitHit = m_genfitHitVec.size();
    ngenfitHit = nPoints;
    std::cout << " m_genfitHitVec size = " << m_genfitHitVec.size() << std::endl;
    //    for(long unsigned int i=0; i<m_genfitHitVec.size();i++)
    //    {
    //        GenfitHit * genfitHit = GetHit(i);
    //        edm4hep::ConstTrackerHit* trackHit =
    //            const_cast<edm4hep::ConstTrackerHit*>(genfitHit->getTrackerHit());
    //
    //        track.addToTrackerHits(*trackHit);
    //    }
    //track.setType();
    track.setChi2(chi2);
    track.setNdf(ndf);
    //track.setDEdx();
    //track.setRadiusOfInnermostHit();//FIXME
    //track.addToTrackerHits();

    //new ReconstructedParticle
    //recParticle->setType();
    //dcRecParticle->setEnergy();

    recParticle.setMomentum(edm4hep::Vector3f(pocaToOrigin_mom.X(),
                pocaToOrigin_mom.Y(),pocaToOrigin_mom.Z()));
    recParticle.setReferencePoint(edm4hep::Vector3f(referencePoint.X(),
                referencePoint.Y(),referencePoint.Z()));
    recParticle.setCharge(charge);
    //recParticle->setMass();
    //recParticle.setCovMatrix(trackState->covMatrix);
    //recParticle->setStartVertex();
    recParticle.addToTracks(track);
    if(m_debug>2){
        std::cout<<m_name<<" storeTrack trackState "<<trackState<<std::endl;
        std::cout<<m_name<<" storeTrack track "<<track<<std::endl;
    }

    return true;
}

void GenfitTrack::pivotToFirstLayer(const edm4hep::Vector3d& pos,
        const edm4hep::Vector3f& mom, edm4hep::Vector3d& firstPos,
        edm4hep::Vector3f& firstMom)
{
    //FIXME, TODO
    firstPos=pos;
    firstMom=mom;
}

int GenfitTrack::getDetTypeID(unsigned long long cellID) const
{
    UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
    encoder.setValue(cellID);
    return encoder[lcio::ILDCellID0::subdet];
}

void GenfitTrack::getAssoSimTrackerHit(
        const edm4hep::MCRecoTrackerAssociationCollection*& assoHits,
        edm4hep::ConstTrackerHit* trackerHit,
        edm4hep::ConstSimTrackerHit& simTrackerHit) const
{
    for(auto assoHit: *assoHits){
        if(assoHit.getRec()==*trackerHit)
        {
            simTrackerHit=assoHit.getSim();

        }
    }
}

genfit::AbsTrackRep* GenfitTrack::getRep(int id) const
{
    return m_track->getTrackRep(id);
}

void GenfitTrack::getSeedCov(TMatrixDSym& cov){
    cov.Zero();
    for(int i=0;i<3;i++) {
        double posResolusion=1.;
        cov(i,i)=posResolusion*posResolusion; //seed position
        double momResolusion=5.;
        cov(i+3,i+3)=momResolusion*momResolusion; //seed momentum
    }
}

//unit conversion here
void GenfitTrack::getEndPointsOfWire(int cellID,TVector3& end0,TVector3& end1) const{
    m_gridDriftChamber->cellposition(cellID,end0,end1);//dd4hep unit
    end0*=(1./dd4hep::cm*GenfitUnit::cm);
    end1*=(1./dd4hep::cm*GenfitUnit::cm);
}
//unit conversion here
void GenfitTrack::getTrackFromMCPartile(const edm4hep::MCParticle mcParticle,
        TVectorD& trackParam, TMatrixDSym& cov) const{
    TVector3 pos;
    TVector3 mom;
    getPosMomFromMCPartile(mcParticle,pos,mom);
    trackParam[0]=pos[0];
    trackParam[1]=pos[1];
    trackParam[2]=pos[2];
    trackParam[3]=mom[0];
    trackParam[4]=mom[1];
    trackParam[5]=mom[2];
    //cov TODO
}
//unit conversion here
void GenfitTrack::getPosMomFromMCPartile(const edm4hep::MCParticle mcParticle,
        TVector3& pos,TVector3& mom) const{
    const edm4hep::Vector3d mcParticleVertex=mcParticle.getVertex();//mm
    const edm4hep::Vector3f mcParticleMom=mcParticle.getMomentum();//GeV
    pos[0]=mcParticleVertex.x*GenfitUnit::mm;
    pos[1]=mcParticleVertex.y*GenfitUnit::mm;
    pos[2]=mcParticleVertex.z*GenfitUnit::mm;
    mom[0]=mcParticleMom.x*GenfitUnit::GeV;
    mom[1]=mcParticleMom.y*GenfitUnit::GeV;
    mom[2]=mcParticleMom.z*GenfitUnit::GeV;
    if(m_debug>2){
        std::cout<<"getPosMomFromTrackState pos("<<pos.X()<<" "<<pos.Y()
            <<" "<<pos.Z();
        std::cout<<") "<<" ("<<mom.X()<<" "<<mom.Y()<<" "<<mom.Z()<<")"
            <<" mom "<<mom.Mag()<<" pt "<<mom.Perp()<<std::endl;
    }
}

//unit conversion here
void GenfitTrack::getTrackFromEDMTrack(const edm4hep::Track& edm4HepTrack,
        double& charge, TVectorD& trackParam, TMatrixDSym& cov) const{
    TVector3 pos;
    TVector3 mom;
    double Bz=m_genfitField->getBz(TVector3{0.,0.,0.})/GenfitUnit::tesla;
    double charge_double;
    CEPC::getPosMomFromTrackState(edm4HepTrack.getTrackStates(0),Bz,pos,mom,charge_double,cov);

    //std::cout<<__LINE__<<" Bz "<<Bz<<" charge "<<charge_double<<std::endl;
    //pos.Print();
    //mom.Print();
    charge=(int) charge_double;
    trackParam[0]=pos[0]*GenfitUnit::mm;
    trackParam[1]=pos[1]*GenfitUnit::mm;
    trackParam[2]=pos[2]*GenfitUnit::mm;
    trackParam[3]=mom[0]*GenfitUnit::GeV;
    trackParam[4]=mom[1]*GenfitUnit::GeV;
    trackParam[5]=mom[2]*GenfitUnit::GeV;

    //cov unit conversion
    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            cov(i,j) = cov(i,j)*GenfitUnit::mm;
        }
    }
}

//to genfit unit
void GenfitTrack::getISurfaceOUV(const dd4hep::rec::ISurface* iSurface,TVector3& o,
        TVector3& u,TVector3& v, double& lengthU,double& lengthV){
    o.SetXYZ(iSurface->origin().x()/dd4hep::mm*GenfitUnit::mm,
            iSurface->origin().y()/dd4hep::mm*GenfitUnit::mm,
            iSurface->origin().z()/dd4hep::mm*GenfitUnit::mm);
    u.SetXYZ(iSurface->u().x()/dd4hep::mm*GenfitUnit::mm,
            iSurface->u().y()/dd4hep::mm*GenfitUnit::mm,
            iSurface->u().z()/dd4hep::mm*GenfitUnit::mm);
    v.SetXYZ(iSurface->v().x()/dd4hep::mm*GenfitUnit::mm,
            iSurface->v().y()/dd4hep::mm*GenfitUnit::mm,
            iSurface->v().z()/dd4hep::mm*GenfitUnit::mm);

    lengthU=iSurface->length_along_u();
    lengthV=iSurface->length_along_v();
}

void GenfitTrack::getMeasurementAndCov(edm4hep::ConstTrackerHit* hit,TVector3& pos,TMatrixDSym& cov){

    pos.SetXYZ(hit->getPosition().x*GenfitUnit::mm,
            hit->getPosition().y*GenfitUnit::mm,
            hit->getPosition().z*GenfitUnit::mm);
}


int GenfitTrack::getSigmas(int cellID,std::vector<float> sigmaUVec,
        std::vector<float> sigmaVVec,float& sigmaU,float& sigmaV)const{
    int detTypeID=getDetTypeID(cellID);
    int sigmaUID=0;
    int sigmaVID=0;
    int layer=0;
    dd4hep::DDSegmentation::BitFieldCoder* decoder;
    if(m_geomSvc->lcdd()->constant<int>("DetID_DC")==detTypeID){
        sigmaUID=0;
        sigmaVID=0;
    }else if(detTypeID==lcio::ILDDetID::VXD){
        decoder=m_geomSvc->getDecoder("VXDCollection");//FIXME
        layer=decoder->get(cellID,"layer");//FIXME
        sigmaUID=layer+1;
        sigmaVID=layer+1;
    }else if(detTypeID==lcio::ILDDetID::SIT){
        sigmaUID=7;
        sigmaVID=7;
    }else if(detTypeID==lcio::ILDDetID::SET){
        sigmaUID=8;
        sigmaVID=8;
    }else if(detTypeID==lcio::ILDDetID::FTD){
        decoder=m_geomSvc->getDecoder("FTDCollection");//FIXME
        layer=decoder->get(cellID,"layer");//FIXME
        sigmaUID=layer+9;
        sigmaVID=layer+9;
    }else{
        if(m_debug>=0) std::cout<<m_name<<" Error: no detType!"<<std::endl;
    }
    if(m_debug){
        std::cout<<"sigmaUID "<<sigmaUID<<" sigmaVID "<<sigmaVID<<std::endl;
        //std::cout<<"pos "<<pos_smeared[0]<<" "<<pos_smeared[1]<<" "<<pos_smeared[2]<<std::endl;
        //std::cout<<"angle "<<atan2(pos_smeared[1],pos_smeared[0])<<std::endl;
        std::cout<<"sigmaU "<<sigmaUVec[sigmaUID]*GenfitUnit::mm
            <<" sigmaV "<<sigmaVVec[sigmaVID]*GenfitUnit::mm<<std::endl;
    }
    sigmaU=sigmaUVec[sigmaUID]*GenfitUnit::mm;
    sigmaV=sigmaVVec[sigmaVID]*GenfitUnit::mm;
    return sigmaUID;
}

bool GenfitTrack::isCDCHit(edm4hep::ConstTrackerHit* hit){
    return m_geomSvc->lcdd()->constant<int>("DetID_DC")==
        getDetTypeID(hit->getCellID());
}

void GenfitTrack::getSortedTrackerHits(
        std::vector<edm4hep::ConstTrackerHit*>& trackerHits,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        std::vector<edm4hep::ConstTrackerHit*>& sortedDCTrackerHits,
        int sortMethod){

    std::vector<std::pair<double,edm4hep::ConstTrackerHit*> > sortedDCTrackerHitPair;
    for(auto trackerHit : trackerHits){
        edm4hep::ConstTrackerHit* thisHit = trackerHit;
        if(!isCDCHit(thisHit))continue;//skip non-DC trackerHit

        edm4hep::ConstSimTrackerHit simTrackerHitAsso;
        getAssoSimTrackerHit(assoHits,trackerHit,simTrackerHitAsso);
        double time=simTrackerHitAsso.getTime();
        if(0==sortMethod){
            //by time
            sortedDCTrackerHitPair.push_back(std::make_pair(time,trackerHit));
            if(m_debug>0){ std::cout<<"sorted DC digi by time"<<std::endl;}
        }else if(1==sortMethod){
            //by layer
            sortedDCTrackerHitPair.push_back(std::make_pair(
                        m_decoderDC->get(trackerHit->getCellID(),"layer"),trackerHit));
            if(m_debug>0){ std::cout<<"sorted DC digi by layer"<<std::endl;}
        }else{
            sortedDCTrackerHits.push_back(trackerHit);
        }
    }
    if(0==sortMethod || 1==sortMethod){
        std::sort(sortedDCTrackerHitPair.begin(),sortedDCTrackerHitPair.end(),
                sortDCDigi);
        for(auto trackerHit:sortedDCTrackerHitPair){
            sortedDCTrackerHits.push_back(trackerHit.second);
        }
    }
    if(m_debug>0){
        std::cout<<"trackerHits on track after sort\n";
        for(auto trackerHit:sortedDCTrackerHits){
            std::cout<<trackerHit->getCellID()<<"("<<std::setw(2)
                <<m_decoderDC->get(trackerHit->getCellID(),"layer")
                <<","<<std::setw(3)
                <<m_decoderDC->get(trackerHit->getCellID(),"cellID")<<")\n";
        }
        std::cout<<"\n"; std::cout.unsetf(std::ios_base::floatfield);
    }
    return;
}//end of getSortedTrackerHits

//make a genfit hit and do hits selection
GenfitHit* GenfitTrack::makeAGenfitHit(edm4hep::ConstTrackerHit* trackerHit,
        edm4hep::ConstSimTrackerHit* simTrackerHitAsso,
        double sigma,bool truthAmbig,double skipCorner,double skipNear){

    //TODO truthAmbig
    double driftVelocity=40.;//FIXME
    GenfitHit* genfitHit=new GenfitHit(trackerHit,simTrackerHitAsso,m_decoderDC,
            m_gridDriftChamber,driftVelocity,sigma*GenfitUnit::mm);
    //skip corner hit
    if(fabs(genfitHit->getDriftDistance())>skipCorner){
        if(m_debug) std::cout<<" skip hit dd > skipCorner"<<std::endl;
        delete genfitHit;
        return nullptr;
    }
    if(genfitHit->getDriftDistance()<skipNear){
        if(m_debug) std::cout<<" skip hit dd < skipCorner"<<std::endl;
        delete genfitHit;
        return nullptr;
    }
    return genfitHit;
}
