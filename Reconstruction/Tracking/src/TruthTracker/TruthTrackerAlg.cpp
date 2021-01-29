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
#include "edm4hep/TrackerHitCollection.h"
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
    declareProperty("SITTrackerHits", m_SITTrackerHitCol,
            "Handle of input SIT hit collection");
    declareProperty("SETTrackerHits", m_SETTrackerHitCol,
            "Handle of input SET hit collection");
}

StatusCode TruthTrackerAlg::initialize()
{
    ///Get geometry
    m_geomSvc=service<IGeomSvc>("GeomSvc");
    if (!m_geomSvc) {
        error() << "Failed to get GeomSvc." << endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Detector
    m_dd4hep=m_geomSvc->lcdd();
    if (nullptr==m_dd4hep) {
        error() << "Failed to get dd4hep::Detector." << endmsg;
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
        error() << "Failed to get the GridDriftChamber" << endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Decoder
    m_decoder = m_geomSvc->getDecoder(m_readout_name);
    if (nullptr==m_decoder) {
        error() << "Failed to get the decoder" << endmsg;
        return StatusCode::FAILURE;
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
    }
    if((int) digiDCHitsCol->size()>m_maxDCDigiCut){
        debug()<<"Cut by m_maxDCDigiCut "<<m_maxDCDigiCut<<endmsg;
        return StatusCode::SUCCESS;
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
    bool isAddSITSeperately=true;
    bool isAddSETSeperately=true;
    if(nullptr!=siTrackCol){
        ///New SDT track
        for(auto siTrack:*siTrackCol){
            debug()<<"siTrack: "<<siTrack<<endmsg;
            edm4hep::Track sdtTrack=sdtTrackCol->create();
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
            sdtTrack.setRadiusOfInnermostHit(siTrack.getRadiusOfInnermostHit());
            for(unsigned int iSiHit=0;iSiHit<siTrack.trackerHits_size();
                    iSiHit++){
                edm4hep::ConstTrackerHit hit=siTrack.getTrackerHits(iSiHit);
                UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
                encoder.setValue(hit.getCellID());
                int detID=encoder[lcio::ILDCellID0::subdet];
                if(detID==lcio::ILDDetID::SIT) isAddSITSeperately=false;
                if(detID==lcio::ILDDetID::SET) isAddSETSeperately=false;
                //debug()<<"siHit "<<iSiHit<<" "<<hit<<endmsg;
                std::cout<<"hit "<<hit<<std::endl;
                sdtTrack.addToTrackerHits(hit);
            }
            //TODO For single track only
            int nSITHit=0;
            if(isAddSITSeperately){
                const edm4hep::TrackerHitCollection* sitTrackerHitCol
                    =m_SITTrackerHitCol.get();
                for(auto sitTrackerHit:*sitTrackerHitCol){
                    sdtTrack.addToTrackerHits(sitTrackerHit);
                    nSITHit++;
                }
            }
            int nSETHit=0;
            if(isAddSETSeperately){
                const edm4hep::TrackerHitCollection* setTrackerHitCol
                    =m_SETTrackerHitCol.get();
                for(auto setTrackerHit:*setTrackerHitCol){
                    sdtTrack.addToTrackerHits(setTrackerHit);
                    nSETHit++;
                }
            }
            int nDCHit=0;
            //TODO tracks
            for(auto digiDC:*digiDCHitsCol){
                //if(Sim->MCParti!=current) continue;//TODO
                sdtTrack.addToTrackerHits(digiDC);
                nDCHit++;
            }
            debug()<<"siTrack trackerHits_size="<<siTrack.trackerHits_size()<<
                " nSITHit_sp "<<nSITHit<<" nSETHit_sp "<<nSETHit<<
                " nDCHit "<<nDCHit<<endmsg;
            debug()<<"sdtTrack nHit "<<sdtTrack.trackerHits_size()<<sdtTrack
                <<endmsg;
        }//end of loop over siTrack
    }

    ///Convert MCParticle to DC Track and ReconstructedParticle
    debug()<<"MCParticleCol size="<<mcParticleCol->size()<<endmsg;
    for(auto mcParticle : *mcParticleCol){
        /// skip mcParticleVertex do not have enough associated hits TODO

        debug()<<"MCParticleCol "<<mcParticle<<endmsg;
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

        ///new Track
        edm4hep::Track dcTrack=dcTrackCol->create();
        edm4hep::TrackState trackState;
        trackState.D0=helix.getD0();
        trackState.phi=helix.getPhi0();
        trackState.omega=helix.getOmega();
        trackState.Z0=helix.getZ0();
        trackState.tanLambda=helix.getTanLambda();
        trackState.referencePoint=helix.getReferencePoint();
        std::array<float,15> covMatrix;
        for(int i=0;i<15;i++){covMatrix[i]=999.;}//FIXME
        trackState.covMatrix=covMatrix;
        dcTrack.addToTrackStates(trackState);
        //dcTrack.setType();//TODO
        //dcTrack.setChi2(gauss(digiDCHitsCol->size-5(),1));//FIXME
        dcTrack.setNdf(digiDCHitsCol->size()-5);
        //dcTrack.setDEdx();//TODO

        debug()<<"D0 "<<trackState.D0<<" phi "<<trackState.phi<<" omega "
            <<trackState.omega<<" Z0 "<<trackState.Z0<<" tanLambda "
            <<trackState.tanLambda<<" referencePoint "
            <<trackState.referencePoint
            <<trackState.covMatrix<<" Bz "<<B[2]/dd4hep::tesla<<endmsg;
        //set hits
        double radiusOfInnermostHit=1e9;
        debug()<<"digiDCHitsCol size "<<digiDCHitsCol->size()<<endmsg;
        for(auto digiDC : *digiDCHitsCol){
            //if(Sim->MCParti!=current) continue;//TODO
            edm4hep::Vector3d digiPos=digiDC.getPosition();
            double r=sqrt(digiPos.x*digiPos.x+digiPos.y*digiPos.y);
            if(r<radiusOfInnermostHit) radiusOfInnermostHit=r;
            dcTrack.addToTrackerHits(digiDC);
        }
        dcTrack.setRadiusOfInnermostHit(radiusOfInnermostHit);//TODO
        debug()<<"DC trackState:location,D0,phi,omega,Z0,tanLambda"
            <<",referencePoint,cov\n"<<trackState<<endmsg;
        debug()<<"dcTrack"<<dcTrack<<endmsg;

        debug()<<"mcParticle "<<mcParticle
            <<" momPhi "<<mcParticleMomPhi
            <<" mcParticleVertex("<<mcParticleVertex<<")mm "
            <<" mcParticleVertexSmeared("<<mcParticleVertexSmeared<<")mm "
            <<" mcParticleMom("<<mcParticleMom<<")GeV "
            <<" mcParticleMomSmeared("<<mcParticleMomSmeared<<")GeV "
            <<" Bxyz "<<B[0]/dd4hep::tesla<<" "<<B[1]/dd4hep::tesla
            <<" "<<B[2]/dd4hep::tesla<<" tesla"<<endmsg;
    }//end loop over MCParticleCol

    debug()<<"Output DCTrack size="<<dcTrackCol->size()<<endmsg;
    return StatusCode::SUCCESS;
}

StatusCode TruthTrackerAlg::finalize()
{
    return GaudiAlgorithm::finalize();
}
