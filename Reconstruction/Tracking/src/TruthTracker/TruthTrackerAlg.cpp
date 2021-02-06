#include "TruthTrackerAlg.h"
#include "DataHelper/HelixClass.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"
#include "DD4hep/DD4hepUnits.h"
#include "UTIL/ILDConf.h"

//external
#include "CLHEP/Random/RandGauss.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"

DECLARE_COMPONENT(TruthTrackerAlg)

TruthTrackerAlg::TruthTrackerAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc),m_dd4hep(nullptr),m_gridDriftChamber(nullptr),
    m_decoder(nullptr)
{
    declareProperty("MCParticle", m_mcParticleCol,
            "Handle of the input MCParticle collection");
    declareProperty("DigiDCHitCollection", m_DCDigiCol,
            "Handle of DC digi(TrackerHit) collection");
    declareProperty("DCHitAssociationCollection", m_DCHitAssociationCol,
            "Handle of association collection");
    declareProperty("DCTrackCollection", m_DCTrackCol,
            "Handle of DC track collection");
    declareProperty("SiSubsetTrackCollection", m_siSubsetTrackCol,
            "Handle of silicon subset track collection");
    declareProperty("SDTTrackCollection", m_SDTTrackCol,
            "Handle of SDT track collection");
    declareProperty("VXDTrackerHits", m_VXDTrackerHits,
            "Handle of input VXD tracker hit collection");
    declareProperty("SITTrackerHits", m_SITTrackerHits,
            "Handle of input SIT tracker hit collection");
    declareProperty("SETTrackerHits", m_SETTrackerHits,
            "Handle of input SET tracker hit collection");
    declareProperty("FTDTrackerHits", m_FTDTrackerHits,
            "Handle of input FTD tracker hit collection");
    declareProperty("SITSpacePoints", m_SITSpacePointCol,
            "Handle of input SIT hit collection");
    declareProperty("SETSpacePoints", m_SETSpacePointCol,
            "Handle of input SET hit collection");
    declareProperty("FTDSpacePoints", m_FTDSpacePointCol,
            "Handle of input FTD hit collection");
    declareProperty("VXDCollection", m_VXDCollection,
            "Handle of input VXD hit collection");
    declareProperty("SITCollection", m_SITCollection,
            "Handle of input SIT hit collection");
    declareProperty("SETCollection", m_SETCollection,
            "Handle of input SET hit collection");
    declareProperty("FTDCollection", m_FTDCollection,
            "Handle of input FTD hit collection");
}


StatusCode TruthTrackerAlg::initialize()
{
    ///Get geometry
    m_geomSvc=service<IGeomSvc>("GeomSvc");
    if (!m_geomSvc) {
        error()<<"Failed to get GeomSvc."<<endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Detector
    m_dd4hep=m_geomSvc->lcdd();
    if (nullptr==m_dd4hep) {
        error()<<"Failed to get dd4hep::Detector."<<endmsg;
        return StatusCode::FAILURE;
    }
    //Get Field
    m_dd4hepField=m_geomSvc->lcdd()->field();
    ///Get Readout
    dd4hep::Readout readout=m_dd4hep->readout(m_readout_name);
    ///Get Segmentation
    m_gridDriftChamber=dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>
        (readout.segmentation().segmentation());
    if(nullptr==m_gridDriftChamber){
        error()<<"Failed to get the GridDriftChamber"<<endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Decoder
    m_decoder = m_geomSvc->getDecoder(m_readout_name);
    if (nullptr==m_decoder) {
        error()<<"Failed to get the decoder"<<endmsg;
        return StatusCode::FAILURE;
    }

    m_tuple=nullptr;
    ///book tuple
    NTuplePtr nt(ntupleSvc(), "TruthTrackerAlg/truthTrackerAlg");
    if(nt){
        m_tuple=nt;
    }else{
        m_tuple=ntupleSvc()->book("TruthTrackerAlg/truthTrackerAlg",
                CLID_ColumnWiseTuple,"TruthTrackerAlg");
        if(m_tuple){
            StatusCode sc;
            sc=m_tuple->addItem("run",m_run);
            sc=m_tuple->addItem("evt",m_evt);
            sc=m_tuple->addItem("siMom",3,m_siMom);
            sc=m_tuple->addItem("siPos",3,m_siPos);
            sc=m_tuple->addItem("mcMom",3,m_mcMom);
            sc=m_tuple->addItem("mcPos",3,m_mcPos);
            sc=m_tuple->addItem("nSimTrackerHitVXD",m_nSimTrackerHitVXD);
            sc=m_tuple->addItem("nSimTrackerHitSIT",m_nSimTrackerHitSIT);
            sc=m_tuple->addItem("nSimTrackerHitSET",m_nSimTrackerHitSET);
            sc=m_tuple->addItem("nSimTrackerHitFTD",m_nSimTrackerHitFTD);
            sc=m_tuple->addItem("nTrackerHitVXD",m_nTrackerHitVXD);
            sc=m_tuple->addItem("nTrackerHitSIT",m_nTrackerHitSIT);
            sc=m_tuple->addItem("nTrackerHitSET",m_nTrackerHitSET);
            sc=m_tuple->addItem("nTrackerHitFTD",m_nTrackerHitFTD);
            sc=m_tuple->addItem("nSpacePointSIT",m_nSpacePointSIT);
            sc=m_tuple->addItem("nSpacePointSET",m_nSpacePointSET);
            sc=m_tuple->addItem("nSpacePointFTD",m_nSpacePointFTD);
            //sc=m_tuple->addItem("nSimTrackerHitErrVXD",3,m_nSimTrackerHitErrVXD);
            //sc=m_tuple->addItem("nSimTrackerHitErrSIT",3,m_nSimTrackerHitErrSIT);
            //sc=m_tuple->addItem("nSimTrackerHitErrSET",3,m_nSimTrackerHitErrSET);
            //sc=m_tuple->addItem("nSimTrackerHitErrFTD",3,m_nSimTrackerHitErrFTD);
            //sc=m_tuple->addItem("nSimSpacePoinErrVXD",3,m_nSimSpacePoinErrVXD);
            //sc=m_tuple->addItem("nSimSpacePoinErrSIT",3,m_nSimSpacePoinErrSIT);
            //sc=m_tuple->addItem("nSimSpacePoinErrSET",3,m_nSimSpacePoinErrSET);
            //sc=m_tuple->addItem("nSimSpacePoinErrFTD",3,m_nSimSpacePoinErrFTD);
        }
    }
    return GaudiAlgorithm::initialize();
}

StatusCode TruthTrackerAlg::execute()
{
    info()<<"In execute()"<<endmsg;

    ///Output Track collection
    edm4hep::TrackCollection* dcTrackCol=m_DCTrackCol.createAndPut();
    edm4hep::TrackCollection* sdtTrackCol=m_SDTTrackCol.createAndPut();
    ///Retrieve MC particle(s)
    const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
    mcParticleCol=m_mcParticleCol.get();
    //if(m_mcParticleCol.exist()){mcParticleCol=m_mcParticleCol.get();}
    if(nullptr==mcParticleCol){
        debug()<<"MCParticleCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }

    ///Retrieve DC digi
    const edm4hep::TrackerHitCollection* digiDCHitsCol=nullptr;
    //if(m_DCDigiCol.exist()){digiDCHitsCol=m_DCDigiCol.get();}//FIXME DEBUG
    digiDCHitsCol=m_DCDigiCol.get();//FIXME DEBUG
    if(nullptr==digiDCHitsCol){
        debug()<<"TrackerHitCollection not found"<<endmsg;
        //return StatusCode::SUCCESS;//FIXME return when no hits in DC + silicon
    }else{
        debug()<<"digiDCHitsCol size "<<digiDCHitsCol->size()<<endmsg;
        if((int) digiDCHitsCol->size()>m_maxDCDigiCut){
            debug()<<"Track cut by m_maxDCDigiCut "<<m_maxDCDigiCut<<endmsg;
            return StatusCode::SUCCESS;
        }
    }

    ////TODO
    //Output MCRecoTrackerAssociationCollection collection
    //const edm4hep::MCRecoTrackerAssociationCollection*
    //    mcRecoTrackerAssociationCol=nullptr;
    //if(nullptr==mcRecoTrackerAssociationCol){
    //    log<<MSG::DEBUG<<"MCRecoTrackerAssociationCollection not found"
    //        <<endmsg;
    //    return StatusCode::SUCCESS;
    //}
    //mcRecoTrackerAssociationCol=m_mcRecoParticleAssociation.get();

    double radiusOfInnermostHit=1e9;
    bool isVXDAdded=false;
    bool isSITAdded=false;
    bool isSETAdded=false;
    bool isFTDAdded=false;

    int nVXDHit=0;
    int nSITHit=0;
    int nSETHit=0;
    int nFTDHit=0;
    int nDCHit=0;

    ///Retrieve silicon Track
    const edm4hep::TrackCollection* siTrackCol=nullptr;
    if(m_siSubsetTrackCol.exist()){
        siTrackCol=m_siSubsetTrackCol.get();
        if(nullptr!=siTrackCol){
            debug()<<"SDTTrackCollection size "<<siTrackCol->size()
                <<endmsg;
            if(0==siTrackCol->size()) return StatusCode::SUCCESS;
        }else{
            debug()<<"SDTTrackCollection is empty"<<endmsg;
            return StatusCode::SUCCESS;
        }
    }
    ///New SDT and DC track
    edm4hep::Track sdtTrack=sdtTrackCol->create();
    edm4hep::Track dcTrack=dcTrackCol->create();
    ///Create track with mcParticle
    edm4hep::TrackState trackStateMc;
    getTrackStateFromMcParticle(mcParticleCol,trackStateMc);
    if(m_useTruthTrack.value()) sdtTrack.addToTrackStates(trackStateMc);
    if((!m_useTruthTrack.value()||(!m_useSiTruthHit)) && nullptr!=siTrackCol){
        for(auto siTrack:*siTrackCol){
            if(!m_useTruthTrack){
                debug()<<"siTrack: "<<siTrack<<endmsg;
                edm4hep::TrackState sdtTrackState;
                edm4hep::TrackState siTrackStat=siTrack.getTrackStates(0);//FIXME?
                sdtTrackState.location=siTrackStat.location;
                sdtTrackState.D0=siTrackStat.D0;
                sdtTrackState.phi=siTrackStat.phi;
                sdtTrackState.omega=siTrackStat.omega;
                sdtTrackState.Z0=siTrackStat.Z0;
                sdtTrackState.tanLambda=siTrackStat.tanLambda;
                sdtTrackState.referencePoint=siTrackStat.referencePoint;
                for(int k=0;k<15;k++){
                    sdtTrackState.covMatrix[k]=siTrackStat.covMatrix[k];
                }
                sdtTrack.addToTrackStates(sdtTrackState);
                sdtTrack.setType(siTrack.getType());
                sdtTrack.setChi2(siTrack.getChi2());
                sdtTrack.setNdf(siTrack.getNdf());
                sdtTrack.setDEdx(siTrack.getDEdx());
                sdtTrack.setDEdxError(siTrack.getDEdxError());
                sdtTrack.setRadiusOfInnermostHit(
                        siTrack.getRadiusOfInnermostHit());
            }
            if(!m_useSiTruthHit){
                for(unsigned int iSiHit=0;iSiHit<siTrack.trackerHits_size();
                        iSiHit++){
                    edm4hep::ConstTrackerHit hit=siTrack.getTrackerHits(iSiHit);
                    UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
                    encoder.setValue(hit.getCellID());
                    int detID=encoder[lcio::ILDCellID0::subdet];
                    if(detID==lcio::ILDDetID::VXD){
                        isVXDAdded=true;
                        debug()<<"add VXD from track "<<endmsg;
                        ++nVXDHit;
                    }
                    if(detID==lcio::ILDDetID::SIT){
                        isSITAdded=true;
                        debug()<<"add SIT from track "<<endmsg;
                        ++nSITHit;
                    }
                    if(detID==lcio::ILDDetID::SET){
                        isSETAdded=true;
                        debug()<<"add SET from track "<<endmsg;
                        ++nSETHit;
                    }
                    if(detID==lcio::ILDDetID::FTD){
                        isFTDAdded=true;
                        debug()<<"add FTD from track "<<endmsg;
                        ++nFTDHit;
                    }
                    debug()<<endmsg<<"siHit "<<iSiHit<<" "<<hit<<endmsg;//got error
                    sdtTrack.addToTrackerHits(hit);

                }//end of loop over hits on siTrack
            }
        }//end of loop over siTrack
    }

    if(m_useSiSpacePoint.value()){
        ///Add vxd silicon trackerHit
        debug()<<"use Si SpacePoint Hit"<<endmsg;
        if(!isVXDAdded){
            const edm4hep::TrackerHitCollection* vxdTrakerHits
                =m_VXDTrackerHits.get();
            debug()<<"add VXD trackerHit "<<vxdTrakerHits->size()<<endmsg;
            for(auto vxdTrakerHit:*vxdTrakerHits){
                sdtTrack.addToTrackerHits(vxdTrakerHit);
                nVXDHit++;
            }
            isVXDAdded=true;
        }
        ///Use silicon spacePoint
        if(!isSITAdded){
            const edm4hep::TrackerHitCollection* sitSpacePointCol
                =m_SITSpacePointCol.get();
            debug()<<"add SIT spacePoint "<<sitSpacePointCol->size()<<endmsg;
            for(auto sitSpacePoint:*sitSpacePointCol){
                sdtTrack.addToTrackerHits(sitSpacePoint);
                nSITHit++;
            }
            isSITAdded=true;
        }
        if(!isSETAdded&&m_useSET){
            const edm4hep::TrackerHitCollection* setSpacePointCol
                =m_SETSpacePointCol.get();
            debug()<<"add SET spacePoint "<<setSpacePointCol->size()<<endmsg;
            for(auto setSpacePoint:*setSpacePointCol){
                sdtTrack.addToTrackerHits(setSpacePoint);
                nSETHit++;
            }
            isSETAdded=true;
        }
        if(!isFTDAdded){
            //TODO
            //const edm4hep::TrackerHitCollection* ftdSpacePointCol
            //    =m_FTDSpacePointCol.get();
            //debug()<<"add FTD spacePoint "<<ftdSpacePointCol->size()<<endmsg;
            //for(auto ftdSpacePoint:*ftdSpacePointCol){
            //    sdtTrack.addToTrackerHits(ftdSpacePoint);
            //    nFTDHit++;
            //}
        }
    }//end of use space point

    ///Add silicon trackerHit
    debug()<<"use Si TrackerHit"<<endmsg;
    if(!isVXDAdded){
        const edm4hep::TrackerHitCollection* vxdTrakerHits
            =m_VXDTrackerHits.get();
        debug()<<"add VXD trackerHit "<<vxdTrakerHits->size()<<endmsg;
        for(auto vxdTrakerHit:*vxdTrakerHits){
            sdtTrack.addToTrackerHits(vxdTrakerHit);
            nVXDHit++;
        }
        isVXDAdded=true;
    }
    if(!isSITAdded){
        const edm4hep::TrackerHitCollection* sitTrakerHits
            =m_SITTrackerHits.get();
        debug()<<"add SIT trackerHit "<<sitTrakerHits->size()<<endmsg;
        for(auto sitTrakerHit:*sitTrakerHits){
            sdtTrack.addToTrackerHits(sitTrakerHit);
            nSITHit++;
        }
        isSITAdded=true;
    }
    if(!isSETAdded){
        const edm4hep::TrackerHitCollection* setTrakerHits
            =m_SETTrackerHits.get();
        debug()<<"add SET trackerHit "<<setTrakerHits->size()<<endmsg;
        for(auto setTrakerHit:*setTrakerHits){
            sdtTrack.addToTrackerHits(setTrakerHit);
            nSETHit++;
        }
        isSETAdded=true;
    }
    if(!isFTDAdded){
        //TODO
        //const edm4hep::TrackerHitCollection* ftdTrakerHits
        //    =m_FTDTrackerHits.get();
        //debug()<<"add FTD trackerHit "<<ftdTrakerHits->size()<<endmsg;
        //for(auto ftdTrakerHit:*ftdTrakerHits){
        //    sdtTrack.addToTrackerHits(ftdTrakerHit);
        //    nFTDHit++;
        //}
        isFTDAdded=true;
    }

    ///Add DC hits to track
    for(auto digiDC:*digiDCHitsCol){
        //if(Sim->MCParti!=current) continue;//TODO
        edm4hep::Vector3d digiPos=digiDC.getPosition();
        double r=sqrt(digiPos.x*digiPos.x+digiPos.y*digiPos.y);
        if(r<radiusOfInnermostHit) radiusOfInnermostHit=r;
        dcTrack.addToTrackerHits(digiDC);
        sdtTrack.addToTrackerHits(digiDC);
        nDCHit++;
    }

    sdtTrack.setNdf(sdtTrack.trackerHits_size()-5);

    ///Create DC Track
    edm4hep::TrackState dcTrackState=trackStateMc;
    dcTrack.addToTrackStates(dcTrackState);
    dcTrack.setNdf(dcTrack.trackerHits_size()-5);
    //track.setType();//TODO
    //track.setChi2(gauss(digiDCHitsCol->size-5(),1));//FIXME
    //track.setDEdx();//TODO
    //track.setRadiusOfInnermostHit(radiusOfInnermostHit);//TODO

    debug()<<"nVXDHit "<<nVXDHit<<" nSITHit "<<nSITHit<<" nSETHit "<<nSETHit
        <<" nFTDHit "<<nFTDHit<<" nDCHit "<<nDCHit<<endmsg;
    debug()<<"sdtTrack nHit "<<sdtTrack.trackerHits_size()<<sdtTrack<<endmsg;
    debug()<<"dcTrack nHit "<<dcTrack.trackerHits_size()<<dcTrack<<endmsg;
    //debug()<<"Output track size="<<trackCol->size()<<endmsg;

    if(m_tuple) debugEvent();

    return StatusCode::SUCCESS;
}

StatusCode TruthTrackerAlg::finalize()
{
    return GaudiAlgorithm::finalize();
}

void TruthTrackerAlg::getTrackStateFromMcParticle(
        const edm4hep::MCParticleCollection* mcParticleCol,
        edm4hep::TrackState& trackState){
    ///Convert MCParticle to DC Track and ReconstructedParticle
    debug()<<"MCParticleCol size="<<mcParticleCol->size()<<endmsg;
    for(auto mcParticle : *mcParticleCol){
        /// skip mcParticleVertex do not have enough associated hits TODO
        ///Vertex
        const edm4hep::Vector3d mcParticleVertex=mcParticle.getVertex();//mm
        edm4hep::Vector3f mcParticleVertexSmeared;//mm
        mcParticleVertexSmeared.x=
            CLHEP::RandGauss::shoot(mcParticleVertex.x,m_resVertexX);
        mcParticleVertexSmeared.y=
            CLHEP::RandGauss::shoot(mcParticleVertex.y,m_resVertexY);
        mcParticleVertexSmeared.z=
            CLHEP::RandGauss::shoot(mcParticleVertex.z,m_resVertexZ);
        ///Momentum
        edm4hep::Vector3f mcParticleMom=mcParticle.getMomentum();//GeV
        double mcParticlePt=sqrt(mcParticleMom.x*mcParticleMom.x+
                mcParticleMom.y*mcParticleMom.y);
        //double mcParticlePtSmeared=
        //    CLHEP::RandGauss::shoot(mcParticlePt,m_resPT);
        double mcParticleMomPhi=atan2(mcParticleMom.y,mcParticleMom.x);
        double mcParticleMomPhiSmeared=
            CLHEP::RandGauss::shoot(mcParticleMomPhi,m_resMomPhi);
        edm4hep::Vector3f mcParticleMomSmeared;
        mcParticleMomSmeared.x=mcParticlePt*cos(mcParticleMomPhiSmeared);
        mcParticleMomSmeared.y=mcParticlePt*sin(mcParticleMomPhiSmeared);
        mcParticleMomSmeared.z=
            CLHEP::RandGauss::shoot(mcParticleMom.z,m_resPz);

        ///Converted to Helix
        double B[3]={1e9,1e9,1e9};
        m_dd4hepField.magneticField({0.,0.,0.},B);
        HelixClass helix;
        //float pos[3]={mcParticleVertexSmeared.x,
        //    mcParticleVertexSmeared.y,mcParticleVertexSmeared.z};
        //float mom[3]={mcParticleMomSmeared.x,mcParticleMomSmeared.y,
        //    mcParticleMomSmeared.z};
        ////FIXME DEBUG
        float pos[3]={(float)mcParticleVertex.x,
            (float)mcParticleVertex.y,(float)mcParticleVertex.z};
        float mom[3]={(float)mcParticleMom.x,(float)mcParticleMom.y,
            (float)mcParticleMom.z};
        helix.Initialize_VP(pos,mom,mcParticle.getCharge(),B[2]/dd4hep::tesla);
        if(m_tuple) {
            for(int ii=0;ii<3;ii++) {
                m_mcMom[ii]=mom[ii];
                m_mcPos[ii]=pos[ii];
            }
        }

        ///new Track
        trackState.D0=helix.getD0();
        trackState.phi=helix.getPhi0();
        trackState.omega=helix.getOmega();
        trackState.Z0=helix.getZ0();
        trackState.tanLambda=helix.getTanLambda();
        trackState.referencePoint=helix.getReferencePoint();
        std::array<float,15> covMatrix;
        for(int i=0;i<15;i++){covMatrix[i]=100.;}//FIXME
        trackState.covMatrix=covMatrix;

        debug()<<"mcParticle "<<mcParticle
            <<" momPhi "<<mcParticleMomPhi
            <<" mcParticleVertex("<<mcParticleVertex<<")mm "
            <<" mcParticleVertexSmeared("<<mcParticleVertexSmeared<<")mm "
            <<" mcParticleMom("<<mcParticleMom<<")GeV "
            <<" mcParticleMomSmeared("<<mcParticleMomSmeared<<")GeV "
            <<" Bxyz "<<B[0]/dd4hep::tesla<<" "<<B[1]/dd4hep::tesla
            <<" "<<B[2]/dd4hep::tesla<<" tesla"<<endmsg;
        debug()<<"D0 "<<trackState.D0<<" phi "<<trackState.phi<<" omega "
            <<trackState.omega<<" Z0 "<<trackState.Z0<<" tanLambda "
            <<trackState.tanLambda<<" referencePoint "
            <<trackState.referencePoint
            <<trackState.covMatrix<<" Bz "<<B[2]/dd4hep::tesla<<endmsg;
        debug()<<"DC trackState:location,D0,phi,omega,Z0,tanLambda"
            <<",referencePoint,cov\n"<<trackState<<endmsg;
    }//end loop over MCParticleCol
}//end of getTrackStateFromMcParticle

void TruthTrackerAlg::debugEvent()
{
    ///Retrieve silicon Track
    const edm4hep::TrackCollection* siTrackCol=nullptr;
    if(m_siSubsetTrackCol.exist()){
        siTrackCol=m_siSubsetTrackCol.get();
        if(nullptr!=siTrackCol){
            debug()<<"SDTTrackCollection size "<<siTrackCol->size()
                <<endmsg;
        }else{
            debug()<<"SDTTrackCollection is empty"<<endmsg;
        }
    }
    if(nullptr!=siTrackCol){
        for(auto siTrack:*siTrackCol){
            debug()<<"siTrack: "<<siTrack<<endmsg;
            edm4hep::TrackState trackStat=siTrack.getTrackStates(0);//FIXME?
            double B[3]={1e9,1e9,1e9};
            m_dd4hepField.magneticField({0.,0.,0.},B);
            HelixClass helix;
            helix.Initialize_Canonical(trackStat.phi, trackStat.D0, trackStat.Z0,
                    trackStat.omega, trackStat.tanLambda, B[2]/dd4hep::tesla);

            if(m_tuple){
                m_siMom[0]=helix.getMomentum()[0];
                m_siMom[1]=helix.getMomentum()[1];
                m_siMom[2]=helix.getMomentum()[2];
                m_siPos[0]=helix.getReferencePoint()[0];
                m_siPos[1]=helix.getReferencePoint()[1];
                m_siPos[2]=helix.getReferencePoint()[2];
            }
        }
    }
    const edm4hep::SimTrackerHitCollection* vxdCollection=m_VXDCollection.get();
    m_nSimTrackerHitVXD=vxdCollection->size();
    const edm4hep::SimTrackerHitCollection* sitCollection=m_SITCollection.get();
    m_nSimTrackerHitSIT=sitCollection->size();
    const edm4hep::SimTrackerHitCollection* setCollection=m_SETCollection.get();
    m_nSimTrackerHitSET=setCollection->size();
    const edm4hep::SimTrackerHitCollection* ftdCollection=m_FTDCollection.get();
    m_nSimTrackerHitFTD=ftdCollection->size();
    const edm4hep::TrackerHitCollection* vxdTrakerHits=m_VXDTrackerHits.get();
    m_nTrackerHitVXD=vxdTrakerHits->size();
    const edm4hep::TrackerHitCollection* sitTrakerHits=m_SITTrackerHits.get();
    m_nTrackerHitSIT=sitTrakerHits->size();
    const edm4hep::TrackerHitCollection* setTrakerHits=m_SETTrackerHits.get();
    m_nTrackerHitSET=setTrakerHits->size();
    const edm4hep::TrackerHitCollection* ftdTrakerHits=m_FTDTrackerHits.get();
    m_nTrackerHitFTD=ftdTrakerHits->size();
    const edm4hep::TrackerHitCollection* sitSpacePointCol=m_SITSpacePointCol.get();
    if(nullptr!=sitSpacePointCol) m_nSpacePointSIT=sitSpacePointCol->size();
    if(m_useSET){
        //if(m_SETSpacePointCol.exist()){
        const edm4hep::TrackerHitCollection* setSpacePointCol=
            m_SETSpacePointCol.get();
        if(nullptr!=setSpacePointCol) m_nSpacePointSET=setSpacePointCol->size();
        //}
    }
    //if(m_FTDSpacePointCol.exist()){
    //    const edm4hep::TrackerHitCollection* ftdSpacePointCol=m_FTDSpacePointCol.get();
    //    if(nullptr!=ftdSpacePointCol) m_nSpacePointSET=ftdSpacePointCol->size();
    //}
    //m_nTrackerHitErrVXD=;
    //m_nTrackerHitErrSIT=;
    //m_nTrackerHitErrSET=;
    //m_nTrackerHitErrFTD=;
    //m_nSpacePointErrVXD=;
    //m_nSpacePointErrSIT=;
    //m_nSpacePointErrSET=;
    //m_nSpacePointErrFTD=;
    if(nullptr!=m_tuple) m_tuple->write();
}
