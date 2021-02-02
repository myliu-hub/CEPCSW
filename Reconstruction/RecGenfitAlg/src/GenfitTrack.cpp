#include "GenfitTrack.h"
#include "GenfitField.h"

//CEPCSW
#include "DataHelper/HelixClass.h"
#include "UTIL/ILDConf.h"

//Externals
#include "DD4hep/DD4hepUnits.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackerHitConst.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/Vector3d.h"
#include "DetSegmentation/GridDriftChamber.h"

//genfit
#include "Track.h"
#include "MeasuredStateOnPlane.h"
#include "RKTrackRep.h"
#include "TrackPoint.h"
#include "StateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "KalmanFittedStateOnPlane.h"
#include "AbsTrackRep.h"
#include "FitStatus.h"
#include "SpacepointMeasurement.h"
#include "WireMeasurementNew.h"

//ROOT
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"

//cpp
#include <cfloat>

const int GenfitTrack::s_PDG[2][5]
={{-11,-13,211,321,2212},{11,13,-211,-321,-2212}};

    bool
sortDCHit(edm4hep::ConstSimTrackerHit hit1,edm4hep::ConstSimTrackerHit hit2)
{
    //std::cout<<"hit1"<<hit1<<std::endl;
    //std::cout<<"hit2"<<hit2<<std::endl;
    bool isEarly=hit1.getTime()<hit2.getTime();
    return isEarly;
}

    GenfitTrack::GenfitTrack(const GenfitField* genfitField, const dd4hep::DDSegmentation::GridDriftChamber* seg, const char* name)
:m_name(name),m_track(nullptr),m_reps(),m_debug(0),
m_genfitField(genfitField),m_gridDriftChamber(seg)
{

}

GenfitTrack::~GenfitTrack()
{
    ///Note: track reps and points will be deleted in the destructor of track
    ///implemented in genfit::Track::Clear()
    delete m_track;
}

/// create a Genfit track from track state, without trackRep
/// Initialize track with seed states
/// NO unit conversion here
bool GenfitTrack::createGenfitTrack(int pdgType,int charge,
        TLorentzVector posInit, TVector3 momInit, TMatrixDSym covMInit_6)
{
    TVectorD seedState(6);
    //TMatrixDSym seedCov(6);

    ////yzhang FIXME
    for(int i = 0; i < 3; ++i) {
        //seed position
        seedState(i)=posInit[i];
        //seed momentum
        seedState(i+3)=momInit[i];
    }

    if(nullptr==m_track) m_track=new genfit::Track();
    m_track->setStateSeed(seedState);
    m_track->setCovSeed(covMInit_6);

    /// new a track representation and add to the track
    int chargeId=0;
    charge>0 ? chargeId=0 : chargeId=1;

    if(m_debug>=2)std::cout<<m_name<<" CreateGenfitTrack seed pos("
        <<seedState[0]<<" "<<seedState[1]<<" "<<seedState[2]<<")cm ("
            <<seedState[3]<<" "<<seedState[4]<<" "<<seedState[5]<<")GeV charge "
            <<charge<<" pdg "<<s_PDG[chargeId][pdgType]<<std::endl;
    if(m_debug>=2)std::cout<<"seedCov "<<std::endl;
    if(m_debug>0) covMInit_6.Print();

    addTrackRep(s_PDG[chargeId][pdgType]);

    return true;
}

///Create a Genfit track with MCParticle, unit conversion here
bool GenfitTrack::createGenfitTrackFromMCParticle(int pidType,
        const edm4hep::MCParticle& mcParticle, double eventStartTime)
{
    if(m_debug>=2)std::cout<<"createGenfitTrackFromMCParticle "<<std::endl;
    ///get track parameters from McParticle
    edm4hep::Vector3d mcPocaPos = mcParticle.getVertex();//mm
    edm4hep::Vector3f mcPocaMom = mcParticle.getMomentum();//GeV
    if(m_debug>=2)std::cout<<"seedPos poca "<< mcPocaPos.x
        <<" "<<mcPocaPos.y<<" "<<mcPocaPos.z<<" mm "<<std::endl;
    if(m_debug>=2)std::cout<<"seedMom poca "<< mcPocaMom.x
        <<" "<<mcPocaMom.y<<" "<<mcPocaMom.z<<" GeV "<<std::endl;

    ///Pivot to first layer to avoid correction of beam pipe
    edm4hep::Vector3d firstLayerPos(1e9,1e9,1e9);
    edm4hep::Vector3f firstLayerMom(1e9,1e9,1e9);
    pivotToFirstLayer(mcPocaPos,mcPocaMom,firstLayerPos,firstLayerMom);

    //TODO convert unit
    ///Get seed position and momentum
    TLorentzVector seedPos(firstLayerPos.x,firstLayerPos.y,firstLayerPos.z,
            eventStartTime);
    TVector3 seedMom(firstLayerMom.x,firstLayerMom.y,firstLayerMom.z);
    if(m_debug>=2)std::cout<<"seedPos "<< firstLayerPos.x
        <<" "<<firstLayerPos.y<<" "<<firstLayerPos.z<<std::endl;
    if(m_debug>=2)std::cout<<"seedMom "<< firstLayerMom.x
        <<" "<<firstLayerMom.y<<" "<<firstLayerMom.z<<std::endl;

    ///Get error matrix of seed track
    TMatrixDSym covMInit_6(6);//FIXME, TODO
    ////yzhang FIXME
    for(int i = 0; i < 3; ++i) {
        double posResolusion=1.;
        //seed position
        covMInit_6(i,i)=posResolusion*posResolusion;
        //seed momentum
        double momResolusion=5.;
        covMInit_6(i+3,i+3)=momResolusion*momResolusion;
    }
    if(m_debug>=2){
        std::cout<<"covMInit_6"<<std::endl;
        covMInit_6.Print();
    }

    ///Create a genfit track with seed
    if(!GenfitTrack::createGenfitTrack(pidType,mcParticle.getCharge(),
                seedPos,seedMom,covMInit_6)){
        if(m_debug>=2)std::cout<<"GenfitTrack"
            <<" Error in createGenfitTrackFromMCParticle" <<std::endl;
        return false;
    }
    if(m_debug>=2)std::cout<<"GenfitTrack \
        createGenfitTrackFromMCParticle success." <<std::endl;

    return true;
}//end of createGenfitTrackFromMCParticle

///Create a Genfit track with MCParticle, unit conversion here
bool GenfitTrack::createGenfitTrackFromEDM4HepTrack(int pidType,
        const edm4hep::Track& track, double eventStartTime)
{
    edm4hep::TrackState trackState=track.getTrackStates(0);//FIXME?
    if(m_debug>=2){
        std::cout<<m_name<<" createGenfitTrackFromEDM4HepTrack "
            <<" Bz "<<m_genfitField->getBz({0.,0.,0.})
            *dd4hep::kilogauss/dd4hep::tesla
            <<" edm4hep::TrackState " <<trackState<<" track "<<track<<std::endl;
    }
    //TODO
    //pivotToFirstLayer(mcPocaPos,mcPocaMom,firstLayerPos,firstLayerMom);
    //Get track parameters
    HelixClass helixClass;
    double D0=trackState.D0;
    double phi=trackState.phi;
    double omega=trackState.omega;
    double Z0=trackState.Z0;
    double tanLambda=trackState.tanLambda;
    double Bz=m_genfitField->getBz({0.,0.,0.})*dd4hep::kilogauss/dd4hep::tesla;
    helixClass.Initialize_Canonical(phi,D0,Z0,omega,tanLambda,Bz);
    TLorentzVector posInit(helixClass.getReferencePoint()[0],
            helixClass.getReferencePoint()[1],
            helixClass.getReferencePoint()[2],eventStartTime);
    posInit.SetX(posInit.X()*dd4hep::mm);
    posInit.SetY(posInit.Y()*dd4hep::mm);
    posInit.SetZ(posInit.Z()*dd4hep::mm);
    TVector3 momInit(helixClass.getMomentum()[0],
            helixClass.getMomentum()[1],helixClass.getMomentum()[2]);
    momInit.SetX(momInit.x()*dd4hep::GeV);
    momInit.SetY(momInit.y()*dd4hep::GeV);
    momInit.SetZ(momInit.z()*dd4hep::GeV);

    if(m_debug>=2){
        std::cout<<m_name<<" createGenfitTrackFromEDM4HepTrack posInit "
            <<std::endl;
        posInit.Print();
        std::cout<<m_name<<" createGenfitTrackFromEDM4HepTrack momInit "
            <<std::endl;
        momInit.Print();
    }
    float charge = helixClass.getCharge();
    //TMatrixDSym covMInit_5(5);//FIXME??? UNIT
    /////< lower triangular covariance matrix of the track parameters.
    /////  the order of parameters is  d0, phi, omega, z0, tan(lambda).
    //std::array<float, 15> covMatrix=trackState.covMatrix;
    //covMInit_5(0,0)=covMatrix[0];
    //covMInit_5(1,0)=covMatrix[1];
    //covMInit_5(1,1)=covMatrix[2];
    //covMInit_5(2,0)=covMatrix[3];
    //covMInit_5(2,1)=covMatrix[4];
    //covMInit_5(2,2)=covMatrix[5];
    //covMInit_5(3,0)=covMatrix[6];
    //covMInit_5(3,1)=covMatrix[7];
    //covMInit_5(3,2)=covMatrix[8];
    //covMInit_5(3,3)=covMatrix[9];
    //covMInit_5(4,0)=covMatrix[10];
    //covMInit_5(4,1)=covMatrix[11];
    //covMInit_5(4,2)=covMatrix[12];
    //covMInit_5(4,3)=covMatrix[13];
    //covMInit_5(4,4)=covMatrix[14];

    TMatrixDSym covMInit_6(6);
    /////Error propagation
    ////V(Y)=S * V(X) * ST , mS = S , mVy = V(Y) , helix.covariance() = V(X)
    //TMatrix mS(covMInit_6.GetNrows(),covMInit_5.GetNrows());
    //mS.Zero();
    //double FCT = 2.99792458E-4;
    //mS[0][0]=-sin(phi);
    //mS[0][1]=-1*D0*cos(phi);
    //mS[1][0]=cos(phi);
    //mS[1][1]=-1*D0*sin(phi);
    //mS[2][3]=1;
    //mS[3][1]=FCT*Bz*(1/omega)*sin(phi);
    //mS[3][2]=charge*FCT*Bz*(1/(omega*omega))*cos(phi);
    //mS[4][1]=-1*FCT*Bz*(1/omega)*cos(phi);
    //mS[4][2]=charge*FCT*Bz*(1/(omega*omega))*sin(phi);
    //mS[5][2]=charge*tanLambda*Bz*(1/(omega*omega));
    //mS[5][4]=-FCT*Bz/omega*(1+tanLambda*tanLambda);

    //covMInit_6= covMInit_5.Similarity(mS);
    //if(m_debug>=2){
    //    std::cout<<m_name<<" covMInit_5 " <<std::endl;
    //    covMInit_5.Print();
    //    std::cout<<m_name<<" mS " <<std::endl;
    //    mS.Print();
    //    std::cout<<m_name<<" covMInit_6 " <<std::endl;
    //    covMInit_6.Print();
    //}

    for(int i = 0; i < 3; ++i) {
        double posResolusion=10.;
        //seed position
        covMInit_6(i,i)=posResolusion*posResolusion;
        //seed momentum
        double momResolusion=50.;
        covMInit_6(i+3,i+3)=momResolusion*momResolusion;
    }
    //TODO ini cov with trackState
    if(!createGenfitTrack(pidType,helixClass.getCharge(),posInit,momInit,
                covMInit_6)){
        if(m_debug>=2){
            std::cout<<m_name<<" Failed in \
                createGenfitTrackFromEDM4HepTrack "<<std::endl;
        }
        return false;
    }
    return true;
}

/// Add a 3d SpacepointMeasurement on TrackerHit
bool GenfitTrack::addSpacePointFromTrakerHit(edm4hep::ConstTrackerHit& hit,
        int hitID)
{
    edm4hep::Vector3d pos=hit.getPosition();
    TVectorD p(3);
    p[0]=pos.x*dd4hep::mm;
    p[1]=pos.y*dd4hep::mm;
    p[2]=pos.z*dd4hep::mm;

    if(m_debug>=2)std::cout<<m_name<<" addSpacePointFromTrakerHit"<<hitID
        <<"pos "<<p[0]<<" "<<p[1]<<" "<<p[2]<<" cm"<<std::endl;
    /// New a SpacepointMeasurement
    double cov[6];
    for(int i=0;i<6;i++) {
        cov[i]=hit.getCovMatrix(i);
        if(m_debug>=2)std::cout<<"cov "<<cov[i]<<std::endl;
    }

    TMatrixDSym hitCov_3(3);
    UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
    encoder.setValue(hit.getCellID());
    int detID=encoder[lcio::ILDCellID0::subdet];

    if(m_debug>=2){
        std::cout<<detID<<" COMPOSITE_SPACEPOINT "<<UTIL::BitSet32(hit.getType())
            [UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]<<std::endl;
        std::cout<<detID<<" ONE_DIMENSIONAL "<<UTIL::BitSet32(hit.getType())
            [UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL]<<std::endl;
    }
    if(UTIL::BitSet32(hit.getType())[UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL]){
        if(m_debug>=2){
            std::cout<<m_name<<" detID "<<detID<<" create Planer hit err"<<std::endl;
        }
        // //in SimpleDigi/src/PlanarDigiAlg.cpp
        // //cov[0] = u_direction[0];//theta
        // //cov[1] = u_direction[1];//phi
        // //cov[2] = resU;
        // //cov[3] = v_direction[0];
        // //cov[4] = v_direction[1];
        // //cov[5] = resV;
        // TMatrix mS(3,2);
        // mS[0][0]=-1*sin(cov[1]);//sin(phi_u)
        // mS[0][1]=0;
        // mS[1][0]=cos(cov[1]);//cos(phi_u)
        // mS[1][1]=0;
        // mS[2][0]=0;
        // mS[2][1]=1;
        // TMatrixDSym covUV_2(2);
        // covUV_2[0][0]=cov[2]*cov[2]*dd4hep::mm*dd4hep::mm;//resU^2
        // //covUV_2[1][1]=cov[5]*cov[5]*dd4hep::mm*dd4hep::mm;//resV^2 ,0 FIXME
        // covUV_2[1][1]=1e9*dd4hep::mm*dd4hep::mm;//resU^2//FIXME TODO
        // covUV_2[0][1]=covUV_2[1][0]=0;
        // hitCov_3=covUV_2.Similarity(mS);
        // if(m_debug>=2){
        //     std::cout<<m_name<<" mS "<<std::endl;
        //     mS.Print();
        //     std::cout<<m_name<<" covUV_2 "<<std::endl;
        //     covUV_2.Print();
        // }
        hitCov_3[0][0]=0.003*0.003;
        hitCov_3[1][1]=0.003*0.003;
        hitCov_3[2][2]=0.003*0.003;
    }else if(UTIL::BitSet32(hit.getType())
            [UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]){
        if(m_debug>=2){
            std::cout<<m_name<<" detID "<<detID<<" create point hit err"<<std::endl;
        }
        //space point error matrix, lower triangle
        hitCov_3[0][0]=cov[0];
        hitCov_3[1][0]=cov[1];
        hitCov_3[1][1]=cov[2];
        hitCov_3[2][0]=cov[3];
        hitCov_3[2][1]=cov[4];
        hitCov_3[2][2]=cov[5];
    }else{
        hitCov_3[0][0]=0.003*0.003;
        hitCov_3[1][1]=0.003*0.003;
        hitCov_3[2][2]=0.003*0.003;
    }
    if(m_debug>=2){
        std::cout<<m_name<<" hitCov_3 "<<std::endl;
        hitCov_3.Print();
    }

    genfit::SpacepointMeasurement* sMeas =
        new genfit::SpacepointMeasurement(p,hitCov_3,(int) hit.getCellID(),hitID,
                nullptr);
    genfit::TrackPoint* trackPoint = new genfit::TrackPoint(sMeas,m_track);
    m_track->insertPoint(trackPoint);

    if(m_debug>=2)std::cout<<"end of addSpacePointFromTrakerHit"<<std::endl;
    return true;
}

/// Add a 3d SpacepointMeasurement with MC truth position smeared by sigma
bool GenfitTrack::addSpacePointMeasurement(const TVectorD& pos,
        double sigma, int detID, int hitID, bool smear)
{
    double sigma_t=sigma*dd4hep::mm;
    /// Convert from CEPCSW unit to genfit unit, cm
    TVectorD pos_t(3);
    pos_t(0)=pos(0)*dd4hep::mm;
    pos_t(1)=pos(1)*dd4hep::mm;
    pos_t(2)=pos(2)*dd4hep::mm;

    /// smear hit position with same weight
    TVectorD pos_smeared(3);
    for (int i=0;i<3;i++){
        pos_smeared[i]=pos_t(i);
        if(smear) pos_smeared[i]+=gRandom->Gaus(0,sigma_t/TMath::Sqrt(3.));
    }

    /// New a SpacepointMeasurement
    TMatrixDSym hitCov(3);
    hitCov(0,0)=sigma_t*sigma_t;
    hitCov(1,1)=sigma_t*sigma_t;
    hitCov(2,2)=sigma_t*sigma_t;

    if(m_debug>=2)std::cout<<m_name<<" addSpacePointMeasurement detID "
        <<detID<<" hitId "<<hitID<<" " <<pos_t[0]<<" "<<pos_t[1]<<" "<<pos_t[2]
            <<" cm smeared "<<pos_smeared[0]<<" "<<pos_smeared[1]<<" "
            <<pos_smeared[2]<<" sigma_t "<<sigma_t<<" cm"<<std::endl;

    genfit::SpacepointMeasurement* sMeas =
        new genfit::SpacepointMeasurement(pos_smeared,hitCov,detID,hitID,nullptr);
    genfit::TrackPoint* trackPoint = new genfit::TrackPoint(sMeas,m_track);
    m_track->insertPoint(trackPoint);

    return true;
}


/// Add a WireMeasurement, no Unit conversion here
void GenfitTrack::addWireMeasurement(double driftDistance,
        double driftDistanceError, const TVector3& endPoint1,
        const TVector3& endPoint2, int lrAmbig, int detID, int hitID)
{
    try{
        /// New a WireMeasurement
        genfit::WireMeasurementNew* wireMeas = new genfit::WireMeasurementNew(
                driftDistance, driftDistanceError, endPoint1, endPoint2, detID,
                hitID, nullptr);
        wireMeas->setMaxDistance(0.5);//0.5 cm FIXME
        wireMeas->setLeftRightResolution(lrAmbig);

        if(m_debug>=2)std::cout<<m_name<<" Add wire measurement(cm) "<<hitID
            <<" ep1("<<endPoint1[0]<<" "<<endPoint1[1]<<" "<<endPoint1[2]
                <<") ep2("<<endPoint2[0]<<" "<<endPoint2[1]<<" "<<endPoint2[2]
                <<") drift "<<driftDistance<<" driftErr "<<driftDistanceError
                <<" lr "<<lrAmbig<<" detId "<<detID << " hitId "<< hitID
                <<std::endl;

        ///New a TrackPoint,create connection between meas. and trackPoint
        genfit::TrackPoint* trackPoint=new genfit::TrackPoint(wireMeas,m_track);
        wireMeas->setTrackPoint(trackPoint);

        m_track->insertPoint(trackPoint);

    }catch(genfit::Exception& e){
        if(m_debug>=2)std::cout<<m_name
            <<"Add wire measurement exception"<<std::endl;
        e.what();
    }
}//end of addWireMeasurementOnTrack

//Add wire measurement on wire, unit conversion here
bool GenfitTrack::addWireMeasurementOnTrack(edm4hep::Track& track,double sigma)
{
    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
        edm4hep::ConstTrackerHit hit=track.getTrackerHits(iHit);

        double driftVelocity=40;//FIXME, TODO, um/ns
        double driftDistance=hit.getTime()*driftVelocity*dd4hep::um;
        TVector3 endPointStart(0,0,0);
        TVector3 endPointEnd(0,0,0);
        m_gridDriftChamber->cellposition(hit.getCellID(),endPointStart,
                endPointEnd);
        int lrAmbig=0;
        if(m_debug>=2)std::cout<<m_name<<" time "<<hit.getTime()
            <<" driftVelocity " <<driftVelocity<<std::endl;
        if(m_debug>=2)std::cout<<m_name<<" wire pos " <<endPointStart.X()
            <<" "<<endPointStart.Y()<<" " <<endPointStart.Z()<<" "
                <<endPointEnd.X()<<" " <<endPointEnd.Y()<<" "
                <<endPointEnd.Z()<<" " <<std::endl;
        endPointStart.SetX(endPointStart.x()*dd4hep::cm);
        endPointStart.SetY(endPointStart.y()*dd4hep::cm);
        endPointStart.SetZ(endPointStart.z()*dd4hep::cm);
        endPointEnd.SetX(endPointEnd.x()*dd4hep::cm);
        endPointEnd.SetY(endPointEnd.y()*dd4hep::cm);
        endPointEnd.SetZ(endPointEnd.z()*dd4hep::cm);
        addWireMeasurement(driftDistance,sigma*dd4hep::cm,endPointStart,
                endPointEnd,lrAmbig,hit.getCellID(),iHit);
    }
    return true;
}//end of addWireMeasurementOnTrack of Track

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
    m_reps.push_back(rep);
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

void GenfitTrack::setDebug(int debug)
{
    m_debug = debug;
    for(unsigned int i=0;i<m_reps.size();i++){ m_reps[i]->setDebugLvl(debug); }
}

void GenfitTrack::printSeed() const
{
    TLorentzVector pos = getSeedStatePos();
    TVector3 mom = getSeedStateMom();
    print(pos,mom);
    if(m_debug>=2)std::cout<<m_name<<" NumPoints "<<getNumPoints()<<std::endl;
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
        for(unsigned int i=0;i<m_reps.size();i++){ m_reps[i]->Print(); }
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
        TMatrixDSym& cov, int repID, bool biased) const
{
    //check number of hit with fitter info
    if(getNumPointsWithFittedInfo(repID)<=2) return 1;

    //get track rep
    genfit::AbsTrackRep* rep = getRep(repID);
    if(nullptr == rep) return 2;

    //get first or last measured state on plane
    genfit::MeasuredStateOnPlane mop;
    try{
        mop = m_track->getFittedState(biased);
    }catch(genfit::Exception& e){
        std::cout<<" getNumPointsWithFittedInfo "
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
    return m_reps[id]->getPDG();
}

int GenfitTrack::getPDGCharge(int id) const
{
    return m_reps[id]->getPDGCharge();
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
        int debug,
        int repID,
        bool stopAtBoundary,
        bool calcJacobianNoise)const
{

    //genfit::MeasuredStateOnPlane state = getMOP(iPoint); // better?
    //genfit::MeasuredStateOnPlane state = getMOP(0);      // better?
    //To do the extrapolation without IHitSelection,above 2 lines are not used.
    pos = pos*dd4hep::cm;//FIXME
    mom = mom*dd4hep::GeV;

    //std::cout<<__LINE__<<" extrapolate to Hit pos and mom"<<std::endl;
    //pos.Print();
    //mom.Print();

    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(getRep(repID)->getPDG());
    rep->setDebugLvl(debug);
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
    //std::cout<<" state "<<std::endl;
    //state.Print();


    // forth
    try {
        //genfit::RKTrackRep* rkRep =
        //  dynamic_cast<genfit::RKTrackRep*> (m_track->getRep(repID));
        //extrapoLen = rkRep->extrapolateToLine(state, end0, end1, poca,
        //    pocaDir, pocaOnWire, stopAtBoundary, calcJacobianNoise);
        extrapoLen = rep->extrapolateToLine(state, end0, end1, poca,
                pocaDir, pocaOnWire, stopAtBoundary, calcJacobianNoise);
    }catch(genfit::Exception& e) {
        if(m_debug>=3)std::cout<<
            "Exception in GenfitterTrack::extrapolateToHit"<<e.what()<<std::endl;
        return extrapoLen;
    }

    //poca = poca*(dd4hep::cm);
    //pocaOnWire = pocaOnWire*(dd4hep::cm);
    pocaOnWire = pocaOnWire;
    doca = (pocaOnWire-poca).Mag();
    //TODO: debug pocaOnWire
    if(m_debug>=2)std::cout<< " poca "<<poca.x()<<","<<poca.y()
        <<" "<<poca.z()<<" doca "<<doca<<std::endl;
    if(m_debug>=2)std::cout<< " pocaOnWire "<<pocaOnWire.x()
        <<","<<pocaOnWire.y()<<" "<<pocaOnWire.z()<<" doca "<<doca<<std::endl;

    return extrapoLen*(dd4hep::cm);
}//end of extrapolateToHit


///Add space point measurement from edm4hep::Track to genfit track
int GenfitTrack::addSimTrackerHitsOnTrack(const edm4hep::Track& track,
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        float sigma,bool smear, bool fitSiliconOnly){

    //A TrakerHit collection
    std::vector<edm4hep::ConstSimTrackerHit> sortedDCTrackHitCol;

    if(m_debug>=2)std::cout<<m_name<<" addSimTrackerHitsOnTrack VXD "
        <<lcio::ILDDetID::VXD<<" SIT "
            <<lcio::ILDDetID::SIT<<" SET "
            <<lcio::ILDDetID::SET<<" FTD "
            <<lcio::ILDDetID::FTD<<" "<<std::endl;

    ///Get TrackerHit on Track
    int hitID=0;
    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
        edm4hep::ConstTrackerHit hit=track.getTrackerHits(iHit);
        UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
        encoder.setValue(hit.getCellID());
        int detID=encoder[lcio::ILDCellID0::subdet];
        if(m_debug>=2)std::cout<<m_name<<" "<<iHit<<" hit "<<hit
            <<" detID "<<detID<<std::endl;

        if(detID==lcio::ILDDetID::VXD || detID==lcio::ILDDetID::SIT ||
                detID==lcio::ILDDetID::SET || detID==lcio::ILDDetID::FTD){
            if(addSpacePointFromTrakerHit(hit,hitID)){
                if(m_debug>=2){
                    std::cout<<"silicon addSpacePointFromTrakerHit ";
                    if(detID==lcio::ILDDetID::VXD)std::cout<<"VXD"<<std::endl;
                    if(detID==lcio::ILDDetID::SIT)std::cout<<"SIT"<<std::endl;
                    if(detID==lcio::ILDDetID::FTD)std::cout<<"FTD"<<std::endl;
                    if(detID==lcio::ILDDetID::SET)std::cout<<"SET"<<std::endl;
                }
                hitID++;
            }else{
                if(m_debug>=2)std::cout<<"silicon addSpacePointFromTrakerHit"
                    <<detID<<" "<<hit.getCellID()<<" failed"<<std::endl;
            }
        }else if(7==detID){
            //if(addSpacePointMeasurement(p,sigma,hit.getCellID(),hitID)){
            //    if(m_debug>=2)std::cout<<"add DC space point"<<std::endl;
            //    hitID++;
            //}else{
            //    if(m_debug>=2)std::cout<<"addSpacePointMeasurement"
            //        <<detID<<" faieled" <<std::endl;
            //}
            float minTime=FLT_MAX;
            edm4hep::ConstSimTrackerHit minTimeSimHit;
            //Select the SimTrakerHit with least time
            for(int iSimHit=0;iSimHit<(int) assoHits->size();iSimHit++){
                if(assoHits->at(iSimHit).getRec()==hit &&
                        assoHits->at(iSimHit).getSim().getTime()<minTime){
                    minTimeSimHit=assoHits->at(iSimHit).getSim();
                    minTime=assoHits->at(iSimHit).getSim().getTime();
                }
            }
            //std::cout<<"minTimeSimHit "<<minTimeSimHit<<std::endl;
            sortedDCTrackHitCol.push_back(minTimeSimHit);
        }else{
            if(m_debug>=2)std::cout<<"addSimTrackerHitsOnTrack Skip add this hit!"
                <<std::endl;
        }
    }//end loop over hit on track

    if(m_debug>=2)std::cout<<" addSimTrakerHits trackerHits_size="
        <<track.trackerHits_size()<<std::endl;

    ///Add DC hits to track
    //Sort sim DC hits by time
    //std::sort(sortedDCTrackHitCol.begin(),sortedDCTrackHitCol.end(),sortDCHit);
    if(!fitSiliconOnly){
        for(auto dCTrackerHit: sortedDCTrackHitCol){
            edm4hep::Vector3d pos=dCTrackerHit.getPosition();
            TVectorD p(3);
            p[0]=pos.x;
            p[1]=pos.y;
            p[2]=pos.z;
            unsigned long long detID = dCTrackerHit.getCellID();
            if(addSpacePointMeasurement(p,sigma,detID,hitID,smear)){
                if(m_debug>=2)std::cout<<"add DC space point"<<std::endl;
                hitID++;
            }else{
                if(m_debug>=2)std::cout<<"addSpacePointMeasurement"
                    <<detID<<" faieled" <<std::endl;
            }
        }
    }

    if(m_debug>=2){
        std::cout<<"GenfitTrack addSimTrackerHitsOnTrack="<<hitID<<std::endl;
    }
    return hitID;
}

double GenfitTrack::extrapolateToPoint(TVector3& pos, TVector3& mom,
        const TVector3& point,
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

        if(nullptr == state) {
            if(m_debug>=2)std::cout<<
                "In extrapolateToPoint KalmanFittedStateOnPlane is null"<<std::endl;
            return trackLength*dd4hep::cm;
        }
        trackLength = rep->extrapolateToPoint(*state,
                point*(1/dd4hep::cm),stopAtBoundary, calcJacobianNoise);
        rep->getPosMom(*state,pos,mom);//FIXME exception exist
        pos = pos*dd4hep::cm;
        mom = mom*dd4hep::GeV;
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
    TMatrixDSym fittedCov;
    TLorentzVector fittedPos;
    TVector3 fittedMom;
    int fittedState=getFittedState(fittedPos,fittedMom,fittedCov);
    if(m_debug>0)std::cout<<m_name<<" fit result: get status OK? pidType "
        <<pidType<<" fittedState==0 " <<(0==fittedState)<<" isFitted "<<isFitted
            <<" isConverged "<<isConverged<<" ndf "<<ndf<<std::endl;
    if((0!=fittedState)||(!isFitted)||(!isConvergedFully)||(ndf>ndfCut)){
        if(m_debug>0)std::cout<<m_name<<" fitting failed"<<std::endl;
    }else{
        if(m_debug>0)std::cout<<m_name<<" fit result: Pos("<<
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
    }

    ///track status at POCA to origin
    const TVector3 origin(0,0,0);
    TVector3 pocaToOrigin_pos(1e9*dd4hep::cm,1e9*dd4hep::cm,1e9*dd4hep::cm);
    TVector3 pocaToOrigin_mom(1e9*dd4hep::GeV,1e9*dd4hep::GeV,1e9*dd4hep::GeV);
    if(extrapolateToPoint(pocaToOrigin_pos,pocaToOrigin_mom,origin)
            > 1e6*dd4hep::cm){
        if(m_debug>0)std::cout<<m_name<<" extrapolate to origin failed"<<std::endl;
        return false;
    }

    float pos[3]={float(pocaToOrigin_pos.X()),float(pocaToOrigin_pos.Y()),
        float(pocaToOrigin_pos.Z())};
    float mom[3]={float(pocaToOrigin_mom.X()),float(pocaToOrigin_mom.Y()),float(pocaToOrigin_mom.Z())};

    if(m_debug>0){std::cout<<m_name<<" fit result poca: Pos("<<
        pocaToOrigin_pos.X()<<" "<<
            pocaToOrigin_pos.Y()<<" "<<
            pocaToOrigin_pos.Z()<<") mom("<<
            pocaToOrigin_mom.X()<<" "<<
            pocaToOrigin_mom.Y()<<" "<<
            pocaToOrigin_mom.Z()<<") p_tot "<<
            pocaToOrigin_mom.Mag()<<" pt "<<
            pocaToOrigin_mom.Perp()<<
            " chi2 "<<chi2<< " ndf "<<ndf <<std::endl;
    }

    HelixClass helix;
    helix.Initialize_VP(pos,mom,charge,m_genfitField->getBz(origin));


    //new TrackState
    edm4hep::TrackState* trackState = new edm4hep::TrackState();
    trackState->location=0;//edm4hep::AtIP;
    trackState->D0=helix.getD0();
    trackState->phi=helix.getPhi0();
    trackState->omega=helix.getOmega();
    trackState->Z0=helix.getZ0();
    trackState->tanLambda=helix.getTanLambda();
    trackState->referencePoint=helix.getReferencePoint();

    //    std::array<float,15> covMatrix;
    //    int k=0;
    //    for(int i=0;i<5;i++){
    //        for(int j=0;j<5;j++){
    //            if(i<=j) covMatrix[k]=;//FIXME
    //        }
    //    }
    //    trackState.covMatrix=

    //new Track
    edm4hep::Track* track = new edm4hep::Track();
    //track->setType();
    track->setChi2(fitState->getChi2());
    track->setNdf(fitState->getNdf());
    //track->setDEdx();
    //track->setRadiusOfInnermostHit();//FIXME
    //track->addToTrackerHits();

    //new ReconstructedParticle
    //recParticle->setType();
    //dcRecParticle->setEnergy();

    edm4hep::Vector3f momVec3(helix.getMomentum()[0],
            helix.getMomentum()[1],helix.getMomentum()[2]);
    recParticle.setMomentum(momVec3);
    //recParticle->setReferencePoint(referencePoint);
    recParticle.setCharge(helix.getCharge());
    //    recParticle->setMass();
    //    recParticle->setCovMatrix();
    //    rcRecParticle->setStartVertex();
    //recParticle->addToTracks(track);
    if(m_debug>2){
        std::cout<<m_name<<" storeTrack trackState "<<*trackState<<std::endl;
        std::cout<<m_name<<" storeTrack track "<<*track<<std::endl;
    }

    return true;
}

void GenfitTrack::pivotToFirstLayer(edm4hep::Vector3d& pos,
        edm4hep::Vector3f& mom, edm4hep::Vector3d& firstPos,
        edm4hep::Vector3f& firstMom)
{
    //FIXME, TODO
    firstPos=pos;
    firstMom=mom;
}
