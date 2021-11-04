#include "GenfitTrack.h"
#include "GenfitField.h"

//CEPCSW
#include "DataHelper/HelixClass.h"
#include "DataHelper/TrackHelper.h"
#include "DataHelper/TrackerHitHelper.h"
#include "DetInterface/IGeomSvc.h"
#include "UTIL/ILDConf.h"

//Externals
#include "GaudiKernel/SmartIF.h"
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

//genfit
#include "AbsMeasurement.h"
#include "AbsTrackRep.h"
#include "FitStatus.h"
#include "KalmanFitterInfo.h"
#include "KalmanFittedStateOnPlane.h"
#include "RKTrackRep.h"
#include "PlanarMeasurement.h"
#include "SpacepointMeasurement.h"
#include "StateOnPlane.h"
#include "RectangularFinitePlane.h"
#include "Track.h"
#include "TrackPoint.h"
#include "MeasuredStateOnPlane.h"
#include "WireMeasurementNew.h"

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
sortDCDigi(std::pair<double,edm4hep::ConstTrackerHit> hitPair1,std::pair<double,edm4hep::ConstTrackerHit> hitPair2)
{
    bool isEarly=hitPair1.first<hitPair2.first;
    return isEarly;
}
    bool
sortDCDigiLayer(std::pair<int,edm4hep::ConstTrackerHit> hitPair1,std::pair<int,edm4hep::ConstTrackerHit> hitPair2)
{
    bool isEarly=hitPair1.first<hitPair2.first;
    return isEarly;
}


GenfitTrack::GenfitTrack(const GenfitField* genfitField,
        const dd4hep::DDSegmentation::GridDriftChamber* seg,
        SmartIF<IGeomSvc> geom, const char* name)
:m_name(name),m_track(nullptr),m_debug(0),m_debugLocal(0),
    m_genfitField(genfitField),m_gridDriftChamber(seg),m_geomSvc(geom)
{

}

GenfitTrack::~GenfitTrack()
{
    ///Note: track reps and points will be deleted in the destructor of track
    ///implemented in genfit::Track::Clear()
    delete m_track;
}

/// create a Genfit track from track state
/// Initialize track with seed state and cov
/// NO unit conversion here
bool GenfitTrack::createGenfitTrack(int pdgType,int charge,
        TLorentzVector posInit, TVector3 momInit, TMatrixDSym covMInit_6)
{
    ///Set track seed state
    TVectorD seedState(6);
    for(int i=0;i<3;i++) {
        seedState(i)=posInit[i]; //seed position
        seedState(i+3)=momInit[i]; //seed momentum
    }

    ///new a track and set seed state and cov
    if(nullptr==m_track) m_track=new genfit::Track();
    m_track->setStateSeed(seedState);
    m_track->setCovSeed(covMInit_6);
#ifdef GENFIT_MY_DEBUG
    //m_track->setDebugLvlLocal(m_debugLocal);
#endif

    ///new a track representation and add to the track
    int chargeId=0;
    charge>0 ? chargeId=0 : chargeId=1;//s_PDG[0]: positive particle
    addTrackRep(s_PDG[chargeId][pdgType]);

    if(m_debug>0){
        std::cout<<m_name<<" CreateGenfitTrack seed pos("
            <<seedState[0]<<" "<<seedState[1]<<" "<<seedState[2]<<")cm ("
            <<seedState[3]<<" "<<seedState[4]<<" "<<seedState[5]<<")GeV charge "
            <<charge<<" pdg "<<s_PDG[chargeId][pdgType]<<std::endl;
        std::cout<<"seedCov "<<std::endl;
        covMInit_6.Print();
    }

    return true;
}

///Create a Genfit track with MCParticle, unit conversion here
bool GenfitTrack::createGenfitTrackFromMCParticle(int pidType,
        const edm4hep::MCParticle& mcParticle, double eventStartTime)
{
    if(m_debug>=2)std::cout<<"createGenfitTrackFromMCParticle "<<std::endl;

    TLorentzVector posInit;
    TVector3 posInit_vector3(mcParticle.getVertex().x,mcParticle.getVertex().y,
            mcParticle.getVertex().z);
    TVector3 momInit(mcParticle.getMomentum().x,mcParticle.getMomentum().y,
            mcParticle.getMomentum().z);

    ///Pivot to first layer to avoid correction of beam pipe
    //edm4hep::Vector3d firstLayerPos(1e9,1e9,1e9);
    //edm4hep::Vector3f firstLayerMom(1e9,1e9,1e9);
    ///Convert MC particle position to first layer mm and GeV
    //pivotToFirstLayer(posInit_vector3,momInit,firstLayerPos,firstLayerMom);

    ///Unit conversion
    posInit.SetX(posInit.X()*dd4hep::mm);
    posInit.SetY(posInit.Y()*dd4hep::mm);
    posInit.SetZ(posInit.Z()*dd4hep::mm);
    momInit.SetX(momInit.X()*dd4hep::GeV);
    momInit.SetY(momInit.Y()*dd4hep::GeV);
    momInit.SetZ(momInit.Z()*dd4hep::GeV);

    ///Get seed position and momentum
    TLorentzVector seedPos(posInit.X(),posInit.Y(),posInit.Z(),eventStartTime);
    TVector3 seedMom(momInit.X(),momInit.Y(),momInit.Z());

    ///Get error matrix of seed track
    TMatrixDSym covMInit_6(6);
    for(int i=0;i<3;i++) {
        double posResolusion=1.;
        covMInit_6(i,i)=posResolusion*posResolusion; //seed position
        double momResolusion=5.;
        covMInit_6(i+3,i+3)=momResolusion*momResolusion; //seed momentum
    }
    if(m_debug>=2){
        std::cout<<"mcPos " << mcParticle.getVertex().x<<" charge "
            <<mcParticle.getCharge()<<" vertex "
            <<mcParticle.getVertex().y<<" "<<mcParticle.getVertex().z<<"\n mcMom "
            <<mcParticle.getMomentum().x<<" "<<mcParticle.getMomentum().y<<" "
            <<mcParticle.getMomentum().z<<"\n ";
    }

    ///Create a genfit track with seed
    bool status=GenfitTrack::createGenfitTrack(pidType,mcParticle.getCharge(),
            seedPos,seedMom,covMInit_6);
    if(!status&&m_debug>=0){std::cout<<m_name
        <<" createGenfitTrackFromMCParticle failed!!!" <<std::endl;}
    return status;
}//end of createGenfitTrackFromMCParticle

///Create a Genfit track with MCParticle, unit conversion here
bool GenfitTrack::createGenfitTrackFromEDM4HepTrack(int pidType,
        const edm4hep::Track& track, double eventStartTime, bool isUseCovTrack)
{
    ///Skip track w.o. hit
    if(track.trackerHits_size()<=0) {
        if(m_debug>=2) std::cout<<m_name<<" skip track n hit=0"<<std::endl;
        return false;
    }
    //TODO
    //pivotToFirstLayer(mcPocaPos,mcPocaMom,firstLayerPos,firstLayerMom);

    ///Get track parameters
    TLorentzVector posInit;
    TVector3 posInit_vector3;
    TVector3 momInit;
    double charge(0);
    TMatrixDSym covMInit_6(6);
    CEPC::getPosMomFromTrackState(track.getTrackStates(0),
            m_genfitField->getBzTesla(TVector3{0.,0.,0.}),
            posInit_vector3,momInit,charge,
            covMInit_6);
    if(m_debug>=2){
        std::cout<<" trackState "<<track.getTrackStates(0)<<std::endl;
        std::cout<<m_name<<" posInit mm" <<std::endl;
        posInit_vector3.Print();
        std::cout<<m_name<<" momInit " <<std::endl;
        momInit.Print();
        std::cout<<m_name<<" covMInit_6 from edm4hep Track" <<std::endl;
        covMInit_6.Print();
        std::cout<<m_name<<" charge "<<charge
            <<" Bz "<<m_genfitField->getBzTesla(TVector3{0.,0.,0.}) <<std::endl;
    }

    ///unit conversion
    posInit.SetX(posInit_vector3.X()*dd4hep::mm);
    posInit.SetY(posInit_vector3.Y()*dd4hep::mm);
    posInit.SetZ(posInit_vector3.Z()*dd4hep::mm);
    posInit.SetT(eventStartTime);
    momInit.SetX(momInit.X()*dd4hep::GeV);
    momInit.SetY(momInit.Y()*dd4hep::GeV);
    momInit.SetZ(momInit.Z()*dd4hep::GeV);
    //unit conversion of error matrix //TODO
    //covMInit_6=

    ///set user defined error matrix
    if(!isUseCovTrack){
        covMInit_6.Zero();
        for(int i = 0; i < 3; ++i) {
            double posResolusion=1.;
            covMInit_6(i,i)=posResolusion*posResolusion; //seed position
            double momResolusion=5.;
            covMInit_6(i+3,i+3)=momResolusion*momResolusion; //seed momentum
        }
    }

    if(m_debug>=2){
        std::cout<<m_name<<" createGenfitTrackFromEDM4HepTrack charge "<<charge
            <<std::endl;
        std::cout<<m_name<<" posInit cm " <<std::endl;
        posInit.Print();
        std::cout<<m_name<<" momInit GeV " <<std::endl;
        momInit.Print();
        std::cout<<m_name<<" covMInit_6 for genfit track" <<std::endl;
        covMInit_6.Print();
        std::cout<<m_name<<" createGenfitTrackFromEDM4HepTrack "
            <<" Bz "<<m_genfitField->getBzTesla({0.,0.,0.})
            <<" n trackerHit "<<track.trackerHits_size()
            <<" TrackState "<<track.getTrackStates(0)
            <<" track "<<track<<std::endl;
    }

    bool status=createGenfitTrack(pidType,charge,posInit,momInit,covMInit_6);
    if(!status && m_debug>=2){
        std::cout<<m_name<<" createGenfitTrackFromEDM4HepTrack failed!!!"
            <<std::endl;
    }
    return status;
}

/// Add a 3d SpacepointMeasurement
bool GenfitTrack::addSpacePointMeasurement(const TVector3& pos,
        std::vector<float> sigmaU,std::vector<float> sigmaV,int cellID,int hitID)
{
    TVectorD pos_smeared(3);
    if(m_debug){
        std::cout<<"addSpacePointMeasurement pos "<<std::endl;
        pos.Print();
    }
    for(int i=0;i<3;i++) pos_smeared[i]=pos(i)*dd4hep::mm;

    /// New a SpacepointMeasurement
    /// space point error matrix and smear, unit cm
    bool smear=(sigmaU[0]>0);
    TMatrixDSym hitCov(3);
    hitCov.Zero();
    //smear 3d track position
    int detTypeID=getDetTypeID(cellID);

    int sigmaUID=0;
    int sigmaVID=0;
    int layer=0;
    dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
    if(m_geomSvc->lcdd()->constant<int>("DetID_DC")==detTypeID){
        sigmaUID=0;
        sigmaVID=0;
    }else if(detTypeID==lcio::ILDDetID::VXD){
        m_decoder=m_geomSvc->getDecoder("VXDCollection");//FIXME
        layer=m_decoder->get(cellID,"layer");//FIXME
        sigmaUID=layer+1;
        sigmaVID=layer+1;
    }else if(detTypeID==lcio::ILDDetID::SIT){
        sigmaUID=7;
        sigmaVID=7;
    }else if(detTypeID==lcio::ILDDetID::SET){
        sigmaUID=8;
        sigmaVID=8;
    }else if(detTypeID==lcio::ILDDetID::FTD){
        m_decoder=m_geomSvc->getDecoder("FTDCollection");//FIXME
        layer=m_decoder->get(cellID,"layer");//FIXME
        sigmaUID=layer+9;
        sigmaVID=layer+9;
    }else{
        if(m_debug>=0) std::cout<<m_name<<" Error: no detType!"<<std::endl;
        return false;
    }

    if(m_debug) std::cout<<" sigmaUID "<<sigmaUID<<" sigmaVID "<<sigmaVID<<std::endl;
    std::cout<<" "<<__LINE__<<" pos "<<pos_smeared[0]<<" "<<pos_smeared[1]<<" "<<pos_smeared[2]<<std::endl;
    std::cout<<" "<<__LINE__<<" angle "<<atan2(pos_smeared[1],pos_smeared[0])<<std::endl;
    float sigmaX=sigmaU[sigmaUID]*dd4hep::mm*cos(atan2(pos_smeared[1],pos_smeared[0]));
    float sigmaY=sigmaU[sigmaUID]*dd4hep::mm*sin(atan2(pos_smeared[1],pos_smeared[0]));
    float sigmaZ=sigmaV[sigmaVID]*dd4hep::mm;
    if(smear){
        pos_smeared[0]+=gRandom->Gaus(0,sigmaX);
        pos_smeared[1]+=gRandom->Gaus(0,sigmaY);
        pos_smeared[2]+=gRandom->Gaus(0,sigmaZ);
    }
    hitCov(0,0)=sigmaX*sigmaX;//FIXME
    hitCov(1,1)=sigmaX*sigmaX;//FIXME
    hitCov(2,2)=sigmaX*sigmaX;//FIXME

    if(m_debug>=2){
        std::cout<<m_name<<" hit "<<hitID<<" detTypeID "<<detTypeID
            <<" pos_smeared "<<pos_smeared[0]<<" "<<pos_smeared[1]
            <<" "<<pos_smeared[2]<<" "<<" hitCov cm"<<std::endl;
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
GenfitTrack::getISurface(edm4hep::ConstTrackerHit hit){
    dd4hep::rec::SurfaceManager surfaceManager(*m_geomSvc->lcdd());

    std::string detectorName;
    unsigned long long cellID=hit.getCellID();
    int detTypeID=getDetTypeID(hit.getCellID());
    if(detTypeID==lcio::ILDDetID::VXD){
        detectorName="VXD";
        if(hit.getPosition().z>0){
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
        <<" hit position "<<hit.getPosition()<<" cellID "<<hit.getCellID()
            <<" cellId+ "<<cellID<<" detTypeID "<<detTypeID<<std::endl;
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
GenfitTrack::addSiliconMeasurement(edm4hep::ConstTrackerHit& hit,
        float sigmaU,float sigmaV,int cellID,int hitID)
{
    if(m_debug>0)std::cout<<"addSiliconMeasurement "<<std::endl;

    ///get surface by cellID
    const dd4hep::rec::ISurface* iSurface = getISurface(hit);
    if(nullptr==iSurface){
        std::cout<<m_name<<" addSiliconMeasurement get surface ERROR!"<<std::endl;
        return false;
    }

    ///Get detector plane parameter
    //VXD
    //SIT
    //SET
    TVector3 o(iSurface->origin().x()*dd4hep::mm,
            iSurface->origin().y()*dd4hep::mm,iSurface->origin().z()*dd4hep::mm);
    TVector3 u(iSurface->u().x()*dd4hep::mm,iSurface->u().y()*dd4hep::mm,
            iSurface->u().z()*dd4hep::mm);
    TVector3 v(iSurface->v().x()*dd4hep::mm,iSurface->v().y()*dd4hep::mm,
            iSurface->v().z()*dd4hep::mm);

    ///Get measurement and cov
    TVectorD hitCoords(2);
    hitCoords(0)=atan2(o.y(),o.x());
    hitCoords(1)=o.z();
    TMatrixDSym hitCov(2);
    hitCov.Zero();
    hitCov(0,0)=sigmaU*sigmaU;
    hitCov(1,1)=sigmaV*sigmaV;

    ///Create planer finite detector plane, measurement and TrackPoint
    genfit::RectangularFinitePlane* pixelOrStripPlane=
        new genfit::RectangularFinitePlane(
                -sigmaU/2.,sigmaU/2.,-sigmaV/2.,sigmaV/2.);
    genfit::SharedPlanePtr plane(new genfit::DetPlane(o,u,v));
    plane->setFinitePlane(pixelOrStripPlane);
    genfit::PlanarMeasurement* planarMeasurement=new genfit::PlanarMeasurement(
            hitCoords,hitCov,cellID,hitID,nullptr);
    planarMeasurement->setPlane(plane);
    m_track->insertPoint(new genfit::TrackPoint(planarMeasurement,m_track));

    if(m_debug>2){
        std::cout<<"hitID "<<hitID<<" unit cm"<<std::endl;
        std::cout<<"u "<<u.x()<<" "<<u.y()<<" "<<u.z()<<std::endl;
        std::cout<<"v "<<v.x()<<" "<<v.y()<<" "<<v.z()<<std::endl;
        std::cout<<"o "<<o.x()<<" "<<o.y()<<" "<<o.z()<<std::endl;
        std::cout<<"hitCoords "<<hitCoords(0)<<" "<<hitCoords(1)<<std::endl;
        std::cout<<"hitCov "<<hitCov(0,0)<<" "<<hitCov(1,1)<<std::endl;
    }

    return true;
}

int GenfitTrack::addSiliconMeasurements(edm4hep::Track& track,
        std::vector<float> sigmaU,std::vector<float> sigmaV)
{
    dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
    ///Get TrackerHit on Track
    int nHitAdd=0;
    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
        unsigned long long cellID=track.getTrackerHits(iHit).getCellID();
        int detTypeID=getDetTypeID(cellID);
        if(m_geomSvc->lcdd()->constant<int>("DetID_DC")==detTypeID) continue;
        edm4hep::ConstTrackerHit hit=track.getTrackerHits(iHit);
        int sigmaUID=0;
        int sigmaVID=0;
        if(detTypeID==lcio::ILDDetID::VXD){
            m_decoder=m_geomSvc->getDecoder("VXDCollection");//FIXME
            int layer=m_decoder->get(cellID,"layer");//FIXME
            sigmaUID=layer+1;
            sigmaVID=layer+1;
        }else if(detTypeID==lcio::ILDDetID::SIT){
            sigmaUID=7;
            sigmaVID=7;
        }else if(detTypeID==lcio::ILDDetID::SET){
            sigmaUID=8;
            sigmaVID=8;
        }else if(detTypeID==lcio::ILDDetID::FTD){
            m_decoder=m_geomSvc->getDecoder("FTDCollection");//FIXME
            int layer=m_decoder->get(cellID,"layer");//FIXME
            sigmaUID=layer+9;
            sigmaVID=layer+9;
        }

        if(m_debug>0){
            std::cout<<sigmaU[sigmaUID]<<" "<<sigmaV[sigmaVID]<<"mm "<<std::endl;
        }
        addSiliconMeasurement(hit,sigmaU[sigmaUID]*dd4hep::mm,
                sigmaV[sigmaVID]*dd4hep::mm,cellID,nHitAdd++);
    }
    return nHitAdd;
}

//Add wire measurement on wire, unit conversion here
int GenfitTrack::addWireMeasurements(edm4hep::Track& track,float sigma,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        int sortMethod, bool truthAmbig,float skipCorner,float skipNear)
{
    int nHitAdd=0;
    if(m_debug>2) std::cout<<"in GenfitTrack::addWireMeasurements"<<std::endl;
    if(0==track.trackerHits_size()){
        if(m_debug>0) std::cout<<"No hit on track"<<std::endl;
        return nHitAdd;
    }

    dd4hep::DDSegmentation::BitFieldCoder* m_decoder
        = m_geomSvc->getDecoder("DriftChamberHitsCollection");

    if(m_debug>0){
        std::cout<<"sort sim hit layer sortMethod "<<sortMethod
            <<" simTrackerHit "<<std::endl;
    }

    std::vector<std::pair<double,edm4hep::ConstTrackerHit> > sortedDCTrackerHit;
    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
        const edm4hep::ConstTrackerHit hit=track.getTrackerHits(iHit);
        if(m_geomSvc->lcdd()->constant<int>("DetID_DC")
                !=getDetTypeID(hit.getCellID())) continue;//skip non-DC hit

        edm4hep::ConstSimTrackerHit simTrackerHitAsso;
        getAssoSimTrackerHit(assoHits,hit,simTrackerHitAsso);
        //simTrackerHitAsso=
        //  CEPC::getAssoClosestSimTrackerHit(assoHits,hit,m_gridDriftChamber,0);
        double time=simTrackerHitAsso.getTime();
        if(0==sortMethod){
            //by time
            sortedDCTrackerHit.push_back(std::make_pair(
                        time,track.getTrackerHits(iHit)));
        }else{
            //by layer
            sortedDCTrackerHit.push_back(std::make_pair(
                        m_decoder->get(hit.getCellID(),"layer"),
                        track.getTrackerHits(iHit)));
        }
        if(m_debug>0){
            std::cout<<"("<<std::setw(2)<<iHit<<","
                <<m_decoder->get(hit.getCellID(),"layer")
                <<","<<std::setw(3)<<m_decoder->get(hit.getCellID(),"cellID")
                <<","<<std::setprecision(5)<<time<<") "<<std::endl;
        }
    }
    if(0==sortedDCTrackerHit.size()){
        if(m_debug>0){ std::cout<<"No DC hit on track"<<std::endl;}
        return nHitAdd;
    }
    if(0==sortMethod){
        std::sort(sortedDCTrackerHit.begin(),sortedDCTrackerHit.end(),sortDCDigi);
        if(m_debug>0){ std::cout<<"addWireMeasurements sorted by time"<<std::endl;}
    }else{
        std::sort(sortedDCTrackerHit.begin(),sortedDCTrackerHit.end(),sortDCDigiLayer);
        if(m_debug>0){ std::cout<<"addWireMeasurements sorted by layer"<<std::endl;}
    }
    if(m_debug>0){
        for(auto hit:sortedDCTrackerHit){
            std::cout<<"("<<std::setw(2)
                <<m_decoder->get(hit.second.getCellID(),"layer")
                <<","<<std::setw(3)
                <<m_decoder->get(hit.second.getCellID(),"cellID")<<")\n";
        }
    }
    if(m_debug>0){ std::cout<<"\n"; std::cout.unsetf(std::ios_base::floatfield);}
    for(auto hitPair:sortedDCTrackerHit){
        edm4hep::ConstTrackerHit hit=hitPair.second;
        double driftVelocity=40.;//FIXME, TODO, um/ns
        double driftDistance=hit.getTime()*driftVelocity*dd4hep::um*dd4hep::cm;
        double driftDistanceSmeared=driftDistance;
        bool smear=(sigma>0);
        if(smear) driftDistanceSmeared+=gRandom->Gaus(0,sigma*dd4hep::mm);
        //cm, skip corner hit
        if(fabs(driftDistance)>skipCorner||fabs(driftDistanceSmeared)>skipCorner){
            if(m_debug) std::cout<<" skip dd "<<driftDistance<<" or "
                <<driftDistanceSmeared<<">"<<skipCorner<<std::endl;
            continue;
        }
        if(driftDistance<skipNear*dd4hep::mm
                ||driftDistanceSmeared<skipNear*dd4hep::mm){//cm
            if(m_debug) std::cout<<" skip dd "<<driftDistance<<" or "
                <<driftDistanceSmeared<<"<"<<skipNear<<" cm"<<std::endl;
            continue;
        }

        TVector3 endPointStart(0,0,0);
        TVector3 endPointEnd(0,0,0);
        m_gridDriftChamber->cellposition(hit.getCellID(),endPointStart,
                endPointEnd);//cm
        double lrAmbig=0;
        endPointStart.SetX(endPointStart.x()*dd4hep::cm);
        endPointStart.SetY(endPointStart.y()*dd4hep::cm);
        endPointStart.SetZ(endPointStart.z()*dd4hep::cm);
        endPointEnd.SetX(endPointEnd.x()*dd4hep::cm);
        endPointEnd.SetY(endPointEnd.y()*dd4hep::cm);
        endPointEnd.SetZ(endPointEnd.z()*dd4hep::cm);

        edm4hep::ConstSimTrackerHit simTrackerHitAsso;
        //simTrackerHitAsso=
        //  CEPC::getAssoClosestSimTrackerHit(assoHits,hit,m_gridDriftChamber,0);
        getAssoSimTrackerHit(assoHits,hit,simTrackerHitAsso);

        int lrAmbigFlag=0;
        if(truthAmbig){
            edm4hep::Vector3d pos=simTrackerHitAsso.getPosition();
            edm4hep::Vector3f mom=simTrackerHitAsso.getMomentum();
            TVector3 pocaOnTrack(pos.x*dd4hep::mm,pos.y*dd4hep::mm,pos.z*dd4hep::mm);
            TVector3 trackDir(mom.x,mom.y,mom.z);
            trackDir=trackDir.Unit();
            TVector3 wireDir=(endPointEnd-endPointStart).Unit();
            TVector3 poca(endPointStart.X(),endPointStart.Y(),pos.z);//axial 
            TVector3 pocaDir=(poca-pocaOnTrack).Unit();
            TVector3 a=pocaDir.Cross(trackDir);
            lrAmbig=(pocaDir.Cross(trackDir))*wireDir;
            lrAmbigFlag=fabs(lrAmbig)/lrAmbig;
        }

        /// New a WireMeasurement
        try{
            genfit::WireMeasurementNew* wireMeas=new genfit::WireMeasurementNew(
                    driftDistanceSmeared,sigma*dd4hep::mm,endPointStart,
                    endPointEnd,hit.getCellID(),nHitAdd++,nullptr);
            wireMeas->setMaxDistance(0.6*1.4);//0.5*sqrt(2) cm FIXME
            wireMeas->setLeftRightResolution(lrAmbigFlag);

#ifdef GENFIT_MY_DEBUG
            //dd4hep::DDSegmentation::BitFieldCoder* m_decoder
            //= m_geomSvc->getDecoder("DriftChamberHitsCollection");

            //wireMeas->setLayerCell(m_decoder->get(cellID,"layer"),m_decoder->get(cellID,"cellID"));
            const edm4hep::Vector3d pos=simTrackerHitAsso.getPosition();
            const edm4hep::Vector3f mom=simTrackerHitAsso.getMomentum();
            wireMeas->SetTruth(
                    TVector3(pos.x*dd4hep::mm,pos.y*dd4hep::mm,pos.z*dd4hep::mm),
                    TVector3(mom.x,mom.y,mom.z),simTrackerHitAsso.getTime()*40/10000.);//cm FIXME
            //std::cout<<" wireMeas set truth pos "<<pos.x*dd4hep::mm<<" "
            //<<pos.y*dd4hep::mm<<" "<<pos.z*dd4hep::mm<<" pos "<<pos<<std::endl;
            //std::cout<<" wireMeas set truth mom "<<mom.x<<" "<<mom.y<<" "
            //<<mom.z<<" mom "<<mom<<std::endl;
#endif
            ///New a TrackPoint,create connection between meas. and trackPoint
            m_track->insertPoint(new genfit::TrackPoint(wireMeas,m_track));
        }catch(genfit::Exception& e){
            if(m_debug>=2)std::cout<<m_name
                <<"Add wire measurement exception"<<std::endl;
            e.what();
        }

        if(m_debug>=2){
            std::cout<<nHitAdd<<"("<<m_decoder->get(hit.getCellID(),"layer")
                <<","<<m_decoder->get(hit.getCellID(),"cellID")
                <<") wire(" <<endPointStart.X()
                <<","<<endPointStart.Y()<<"," <<endPointStart.Z()<<") ("
                <<endPointEnd.X()<<"," <<endPointEnd.Y()<<","
                <<endPointEnd.Z()<<") dt "<<hit.getTime()
                <<" time "<<simTrackerHitAsso.getTime()
                <<" dd "<<driftDistance
                <<" ddSm "<<driftDistanceSmeared
                <<" sigma "<<sigma*dd4hep::mm<<"cm lr "<<lrAmbigFlag<<std::endl;
        }
    }
    return nHitAdd;
}//end of addWireMeasurements

/// Get MOP
bool GenfitTrack::getMOP(int hitID,
        genfit::MeasuredStateOnPlane& mop, genfit::AbsTrackRep* trackRep) const
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
genfit::RKTrackRep* GenfitTrack::addTrackRep(int pdg)
{
    /// create a new track representation
    genfit::RKTrackRep* rep = new genfit::RKTrackRep(pdg);
    //m_reps.push_back(rep);
    //rep->setDebugLvl(m_debug-10);
    m_track->addTrackRep(rep);
    //try{
    //  genfit::MeasuredStateOnPlane stateInit(rep);
    //  rep->setPosMomCov(stateInit, pos, mom, covM);
    //}catch(genfit::Exception e){
    //  if(m_debug>=2)std::cout<<m_name<<" Exception in set track status"
    //  <<std::endl     ;
    //  std::cout<<e.what()<<std::endl;
    //  return false;
    //}
    return rep;
}

/// Get the position from genfit::Track::getStateSeed
const TLorentzVector GenfitTrack::getSeedStatePos()const
{
    TVectorD seedStat(6);
    seedStat = m_track->getStateSeed();
    TVector3 p(seedStat[0],seedStat[1],seedStat[2]);
    p = p*dd4hep::cm;
    TLorentzVector pos(p[0],p[1],p[2],9999);//FIXME
    return pos;
}

/// Get the momentum from genfit::Track::getStateSeed
const TVector3 GenfitTrack::getSeedStateMom() const
{
    TVectorD seedStat(6); seedStat = m_track->getStateSeed();
    TVector3 mom(seedStat[3],seedStat[4],seedStat[5]);
    return mom*dd4hep::GeV;
}

/// Get the seed states of momentum and position
void GenfitTrack::getSeedStateMom(TLorentzVector& pos, TVector3& mom) const
{
    TVectorD seedStat(6); seedStat = m_track->getStateSeed();
    mom = TVector3(seedStat[3],seedStat[4],seedStat[5])*dd4hep::GeV;
    seedStat = m_track->getStateSeed();
    TVector3 p = TVector3(seedStat[0],seedStat[1],seedStat[2])*dd4hep::cm;
    pos.SetXYZT(p[0],p[1],p[2],9999);//FIXME time
}

unsigned int GenfitTrack::getNumPoints() const
{
    return m_track->getNumPoints();
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
    std::cout<<"GenfitTrack::printSeed "<<std::endl;
    TLorentzVector pos = getSeedStatePos();
    TVector3 mom = getSeedStateMom();
    print(pos,mom);
    std::cout<<"GenfitTrack NumPoints "<<m_track->getNumPoints()
        <<" NumPointsWithMeasurement "
        <<m_track->getNumPointsWithMeasurement()<<std::endl;
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
    std::cout<<"pos("<<pos.X()<<" "<<pos.Y()<<" "<<pos.Z();
    std::cout<<") "<<" ("<<mom.X()<<" "<<mom.Y()<<" "<<mom.Z()<<")"
        <<" mom "<<mom.Mag()<<" pt "<<mom.Perp()<<std::endl;

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
        TVector3& mom, TMatrixDSym& cov, int repID) const
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
    int nHit = m_track->getNumPointsWithMeasurement();
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
        std::cout<<e.what()<<std::endl;
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
/// return StateOnPlane of this POCA
/// inputs
///  pos,mom ... position & momentum at starting point, unit= [mm]&[GeV/c]
///              (converted to [cm]&[GeV/c] inside this function)
///  hit ... destination
/// outputs poca [mm] (converted from [cm] in this function) ,global
double GenfitTrack::extrapolateToHit( TVector3& poca, TVector3& pocaDir,
        TVector3& pocaOnWire, double& doca, TVector3 pos, TVector3 mom,
        TVector3 end0,//one end of the hit wire
        TVector3 end1,//the orhter end of the hit wire
        int repID,
        bool stopAtBoundary,
        bool calcJacobianNoise)const
{

    //genfit::MeasuredStateOnPlane state = getMOP(iPoint); // better?
    //genfit::MeasuredStateOnPlane state = getMOP(0);      // better?
    //To do the extrapolation without IHitSelection,above 2 lines are not used.
    pos = pos*dd4hep::cm;//FIXME
    mom = mom*dd4hep::GeV;

    if(m_debug>2){
        std::cout<<"GenfitterTrack::extrapolateToHit before pos and mom "
            <<" pdg "<<getRep(repID)->getPDG()<<std::endl;
        pos.Print();
        mom.Print();
    }

    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(getRep(repID)->getPDG());
    //rep->setDebugLvl(debug-10);
    genfit::MeasuredStateOnPlane state(rep);
    rep->setPosMom(state, pos, mom);


    //genfit::MeasuredStateOnPlane state(m_track->getRep(repID));
    //m_track->getRep(repID)->setPosMom(state, pos, mom);

    //m_track->PrintSeed();
    double extrapoLen(0);
    //std::cout<<" wire1 "<<std::endl;
    //end0.Print();
    //std::cout<<" wire0 "<<std::endl;
    //end1.Print();
    if(m_debug){
        std::cout<<" state "<<std::endl;
        state.Print();
    }


    // forth
    try {
        //genfit::RKTrackRep* rkRep =
        //  dynamic_cast<genfit::RKTrackRep*> (m_track->getRep(repID));
        //extrapoLen = rkRep->extrapolateToLine(state, end0, end1, poca,
        //    pocaDir, pocaOnWire, stopAtBoundary, calcJacobianNoise);
        extrapoLen = rep->extrapolateToLine(state,end0*dd4hep::mm,
                end1*dd4hep::mm,poca,pocaDir,pocaOnWire,stopAtBoundary,
                calcJacobianNoise);
    }catch(genfit::Exception& e) {
        if(m_debug>=3)std::cout<<
            "Exception in GenfitterTrack::extrapolateToHit"<<e.what()<<std::endl;
        return extrapoLen;
    }

    //poca = poca*(dd4hep::cm);
    //pocaOnWire = pocaOnWire*(dd4hep::cm);
    poca.SetZ(0);
    pocaOnWire.SetZ(0);
    doca = (pocaOnWire-poca).Mag();
    //std::cout<<" debug poca yzhang "<<std::endl;
    //poca.Print();
    //pocaOnWire.Print();
    //TODO: debug pocaOnWire
    if(m_debug>=2)std::cout<< " poca "<<poca.x()<<","<<poca.y()
        <<" "<<poca.z()<<" doca "<<doca<<std::endl;
    if(m_debug>=2)std::cout<< " pocaDir "<<pocaDir.x()
        <<","<<pocaDir.y()<<" "<<pocaDir.z()<<std::endl;
    if(m_debug>=2)std::cout<< " pocaOnWire "<<pocaOnWire.x()
        <<","<<pocaOnWire.y()<<" "<<pocaOnWire.z()<<" doca "<<doca<<std::endl;

    delete rep;
    return extrapoLen*(dd4hep::cm);
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

        std::cout<<" "<<__LINE__<<" addSpacePointsSi hit "<<hit<<std::endl;
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
        std::cout<<m_name<<" addSpacePointsDC nTrackerHit="
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

        std::cout<<" "<<__LINE__<<" addSpacePointsDC hit "<<dCTrackerHit<<std::endl;
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


double GenfitTrack::extrapolateToPoint(TVector3& pos, TVector3& mom,
        const TVector3& point,
        int repID,// same with pidType
        bool stopAtBoundary,
        bool calcJacobianNoise) const
{
    TMatrixDSym cov;
    return extrapolateToPoint(pos,mom,cov,point,repID,stopAtBoundary,
            calcJacobianNoise);

}//end of extrapolateToPoint


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
            return trackLength*dd4hep::cm;
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
            return trackLength*dd4hep::cm;
        }
        trackLength = rep->extrapolateToPoint(*state,
                point*(1/dd4hep::cm),stopAtBoundary, calcJacobianNoise);
        rep->getPosMomCov(*state,pos,mom,cov);//FIXME exception exist
        pos = pos*dd4hep::cm;
        mom = mom*dd4hep::GeV;
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
    double trackLength(1e9*dd4hep::cm);
    if(!getFitStatus(repID)->isFitted()) return trackLength;
    try{
        // get track rep
        genfit::AbsTrackRep* rep = getRep(repID);
        if(nullptr == rep) {
            if(m_debug>=2)std::cout<<"In extrapolateToCylinder rep is null"
                <<std::endl;
            return trackLength*dd4hep::cm;
        }

        // get track point with fitter info
        genfit::TrackPoint* tp =
            getTrack()->getPointWithFitterInfo(hitID,rep);
        if(nullptr == tp) {
            if(m_debug>=2)std::cout<<
                "In extrapolateToCylinder TrackPoint is null"<<std::endl;
            return trackLength*dd4hep::cm;
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
        rep->setPosMom(*state, pos*(1/dd4hep::cm), mom*(1/dd4hep::GeV));

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

bool GenfitTrack::storeTrack(edm4hep::ReconstructedParticle& recParticle,
        edm4hep::Track& track,
        int pidType, int ndfCut, double chi2Cut)
{

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

    ///track status at POCA to referencePoint: origin
    const TVector3 referencePoint(0,0,0);
    TVector3 pocaToOrigin_pos(1e9*dd4hep::cm,1e9*dd4hep::cm,1e9*dd4hep::cm);
    TVector3 pocaToOrigin_mom(1e9*dd4hep::GeV,1e9*dd4hep::GeV,1e9*dd4hep::GeV);
    if(extrapolateToPoint(pocaToOrigin_pos,pocaToOrigin_mom,referencePoint)
            > 1e6*dd4hep::cm){
        if(m_debug>0)std::cout<<m_name<<" extrapolate to origin failed"<<std::endl;
        return false;
    }

    //Unit conversion of position and momentum
    pocaToOrigin_pos.SetXYZ(pocaToOrigin_pos.X()*dd4hep::mm,
            pocaToOrigin_pos.Y()*dd4hep::mm,pocaToOrigin_pos.Z()*dd4hep::mm);
    pocaToOrigin_mom.SetXYZ(pocaToOrigin_mom.X()*dd4hep::GeV,
            pocaToOrigin_mom.Y()*dd4hep::GeV,pocaToOrigin_mom.Z()*dd4hep::GeV);
    //unit conversion of error matrix
    TMatrixDSym covMatrix_6=fittedCov;
    for(int i=0;i<5;i++){
        covMatrix_6[0][i]=fittedCov[0][i]*dd4hep::cm;//d0 column
        covMatrix_6[2][i]=fittedCov[2][i]/dd4hep::cm;//omega column
        covMatrix_6[3][i]=fittedCov[3][i]*dd4hep::cm;//z0 column
        covMatrix_6[i][0]=fittedCov[i][0]*dd4hep::cm;//d0 row
        covMatrix_6[i][2]=fittedCov[i][2]/dd4hep::cm;//omega row
        covMatrix_6[i][3]=fittedCov[i][3]*dd4hep::cm;//z0 row
    }

    if(m_debug>0){
        std::cout<<m_name<<" fit result poca: Pos"<<std::endl;
        pocaToOrigin_pos.Print();
        pocaToOrigin_mom.Print();
        std::cout<<" chi2 "<<chi2<< " ndf "<<ndf <<std::endl;
        std::cout<<"fittedCov "<<std::endl;
        fittedCov.Print();
        std::cout<<"covMatrix_6 "<<std::endl;
        covMatrix_6.Print();
    }

    double Bz=m_genfitField->getBzTesla(referencePoint);
    edm4hep::TrackState trackState;
    CEPC::getTrackStateFromPosMom(trackState,Bz,pocaToOrigin_pos,
            pocaToOrigin_mom,charge,covMatrix_6);
    trackState.location=0;//edm4hep::AtIP;//FIXME
    if(m_debug>2){std::cout<<m_name<<" trackState "<<trackState<<std::endl;}
    track.addToTrackStates(trackState);


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
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        edm4hep::ConstTrackerHit trackerHit,
        edm4hep::ConstSimTrackerHit& simTrackerHit) const
{
    for(auto assoHit: *assoHits){
        if(assoHit.getRec()==trackerHit){ simTrackerHit=assoHit.getSim(); }
    }
}

genfit::AbsTrackRep* GenfitTrack::getRep(int id) const
{
    return m_track->getTrackRep(id);
}
