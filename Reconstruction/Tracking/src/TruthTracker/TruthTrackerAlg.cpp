#include "TruthTrackerAlg.h"
#include "DataHelper/HelixClass.h"
#include "DataHelper/TrackHelper.h"
#include "DataHelper/TrackerHitHelper.h"
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
    declareProperty("DriftChamberHitsCollection", m_DCSimTrackerHitCol,
            "Handle of DC SimTrackerHit collection");
    declareProperty("DigiDCHitCollection", m_DCDigiCol,
            "Handle of DC digi(TrackerHit) collection");
    declareProperty("DCHitAssociationCollection", m_DCHitAssociationCol,
            "Handle of association collection");
    declareProperty("DCTrackCollection", m_DCTrackCol,
            "Handle of DC track collection");
    declareProperty("SubsetTracks", m_siSubsetTrackCol,
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
    declareProperty("FTDSpacePoints", m_FTDSpacePointCol,
            "Handle of input FTD hit collection");
    declareProperty("VXDCollection", m_VXDCollection,
            "Handle of input VXD hit collection");
    declareProperty("VXDTrackerHitAssociation", m_VXDHitAssociationCol,
            "Handle of VXD association collection");
    declareProperty("SETTrackerHitAssociation", m_SETHitAssociationCol,
            "Handle of SET association collection");
    declareProperty("SITTrackerHitAssociation", m_SITHitAssociationCol,
            "Handle of SIT association collection");
    declareProperty("FTDTrackerHitAssociation", m_FTDHitAssociationCol,
            "Handle of FTD association collection");
    declareProperty("SITCollection", m_SITCollection,
            "Handle of input SIT hit collection");
    declareProperty("SETCollection", m_SETCollection,
            "Handle of input SET hit collection");
    declareProperty("FTDCollection", m_FTDCollection,
            "Handle of input FTD hit collection");
    declareProperty("TruthTrackerHitCollection", m_truthTrackerHitCol,
            "Handle of output truth TrackerHit collection");
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

    eventNo = 0;

    ///book tuple
    if(m_hist){
        NTuplePtr nt(ntupleSvc(), "TruthTrackerAlg/truthTrackerAlg");
        if(nt){
            m_tuple=nt;
        }else{
            m_tuple=ntupleSvc()->book("TruthTrackerAlg/truthTrackerAlg",
                    CLID_ColumnWiseTuple,"TruthTrackerAlg");
            if(m_tuple){
                StatusCode sc;
                sc=m_tuple->addItem("nDCTrackHit",m_nDCTrackHit,0,100000);

                sc=m_tuple->addItem("nSDTTrack",m_nSDTTrack,0,100000);
                sc=m_tuple->addItem("nSimTrackerHitVXD",m_nSimTrackerHitVXD);
                sc=m_tuple->addItem("nSimTrackerHitSIT",m_nSimTrackerHitSIT);
                sc=m_tuple->addItem("nSimTrackerHitSET",m_nSimTrackerHitSET);
                sc=m_tuple->addItem("nSimTrackerHitFTD",m_nSimTrackerHitFTD);
                sc=m_tuple->addItem("nSimTrackerHitDC",m_nSimTrackerHitDC);
                sc=m_tuple->addItem("nTrackerHitVXD",m_nTrackerHitVXD);
                sc=m_tuple->addItem("nTrackerHitSIT",m_nTrackerHitSIT);
                sc=m_tuple->addItem("nTrackerHitSET",m_nTrackerHitSET);
                sc=m_tuple->addItem("nTrackerHitFTD",m_nTrackerHitFTD);
                sc=m_tuple->addItem("nTrackerHitDC",m_nTrackerHitDC);
                sc=m_tuple->addItem("nHitOnSiTkXVD",m_nSDTTrack,m_nHitOnSiTkVXD);
                sc=m_tuple->addItem("nHitOnSiTkSIT",m_nSDTTrack,m_nHitOnSiTkSIT);
                sc=m_tuple->addItem("nHitOnSiTkSET",m_nSDTTrack,m_nHitOnSiTkSET);
                sc=m_tuple->addItem("nHitOnSiTkFTD",m_nSDTTrack,m_nHitOnSiTkFTD);
                sc=m_tuple->addItem("nHitOnSdtTkVXD",m_nSDTTrack,m_nHitOnSdtTkVXD);
                sc=m_tuple->addItem("nHitOnSdtTkSIT",m_nSDTTrack,m_nHitOnSdtTkSIT);
                sc=m_tuple->addItem("nHitOnSdtTkSET",m_nSDTTrack,m_nHitOnSdtTkSET);
                sc=m_tuple->addItem("nHitOnSdtTkFTD",m_nSDTTrack,m_nHitOnSdtTkFTD);
                sc=m_tuple->addItem("nHitOnSdtTkDC",m_nSDTTrack,m_nHitOnSdtTkDC);
                sc=m_tuple->addItem("nHitSdt",m_nSDTTrack,m_nHitOnSdtTk);
                sc=m_tuple->addItem("siXc",m_nSDTTrack,m_siXc);
                sc=m_tuple->addItem("siYc",m_nSDTTrack,m_siYc);
                sc=m_tuple->addItem("siR",m_nSDTTrack,m_siR);
            }
        }
    }
    return GaudiAlgorithm::initialize();
}

StatusCode TruthTrackerAlg::execute()
{
    info()<<"In execute()"<<endmsg;

    ///Output DC Track collection
    edm4hep::TrackCollection* dcTrackCol= m_DCTrackCol.createAndPut();

    ///Output SDT Track collection
    edm4hep::TrackCollection* sdtTkCol= m_SDTTrackCol.createAndPut();

    ///Output Hit collection
    auto truthTrackerHitCol = m_truthTrackerHitCol.createAndPut();

    //VXD MCRecoTrackerAssociation
    const edm4hep::MCRecoTrackerAssociationCollection* VXDAssoVec = m_VXDHitAssociationCol.get();
    //SIT MCRecoTrackerAssociation
    const edm4hep::MCRecoTrackerAssociationCollection* SITAssoVec = m_SITHitAssociationCol.get();
    //SET MCRecoTrackerAssociation
    const edm4hep::MCRecoTrackerAssociationCollection* SETAssoVec = m_SETHitAssociationCol.get();
    //FTD MCRecoTrackerAssociation
    const edm4hep::MCRecoTrackerAssociationCollection* FTDAssoVec = m_FTDHitAssociationCol.get();

    ///Retrieve MC particle(s)
    const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
    mcParticleCol=m_mcParticleCol.get();
    const int mcPNum = mcParticleCol->size();
    float mcP_Mom[mcPNum]={1e-9};

    if(nullptr==mcParticleCol){
        debug()<<"MCParticleCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    for(int i=0;i<mcPNum;i++){
        mcP_Mom[i] = sqrt(mcParticleCol->at(i).getMomentum()[0]*mcParticleCol->at(i).getMomentum()[0]
                +mcParticleCol->at(i).getMomentum()[1]*mcParticleCol->at(i).getMomentum()[1]);
    }

    ///Retrieve DC digi
    const edm4hep::TrackerHitCollection* digiDCHitsCol=nullptr;
    digiDCHitsCol=m_DCDigiCol.get();
    if(digiDCHitsCol->size()<1e-9){
        debug() << " TruthTrackerAlg digiDCHitsCol size =0 !!" << endmsg;
        return StatusCode::SUCCESS;
    }
    const edm4hep::SimTrackerHitCollection* dcSimHitCol
        =m_DCSimTrackerHitCol.get();
    const edm4hep::MCRecoTrackerAssociationCollection* assoHits
        = m_DCHitAssociationCol.get();
    if(m_useDC){
        if(nullptr==dcSimHitCol){
            debug()<<"DC SimTrackerHitCollection not found"<<endmsg;
        }else{
            debug()<<"DriftChamberHitsCollection size "
                <<dcSimHitCol->size()<<endmsg;
            if(m_tuple)m_nSimTrackerHitDC=dcSimHitCol->size();
        }
        std::cout << " digiDCHitsCol size = " << digiDCHitsCol->size() << std::endl;
        if(m_tuple) m_nDCTrackHit = digiDCHitsCol->size();
        if(nullptr==digiDCHitsCol){
            debug()<<"TrackerHitCollection not found"<<endmsg;
        }else{
            debug()<<"digiDCHitsCol size "<<digiDCHitsCol->size()<<endmsg;
            if((int) digiDCHitsCol->size()>m_maxDCDigiCut){
                debug()<<"Track cut by m_maxDCDigiCut "<<m_maxDCDigiCut<<endmsg;
                return StatusCode::SUCCESS;
            }
        }

        //Test DCHit Mom
        if(m_useCutMom){
            const int dchitNum = digiDCHitsCol->size();
            float dchitMom[dchitNum] = {1e-9};
            int nloop =0;
            for(auto dcHit:*digiDCHitsCol){

                edm4hep::SimTrackerHit simDCHit = getAssoSimTrackerHit(assoHits,dcHit);
                dchitMom[nloop] = sqrt(simDCHit.getMomentum()[0]*simDCHit.getMomentum()[0]+
                        simDCHit.getMomentum()[1]*simDCHit.getMomentum()[1]);
                float dchitMomP = sqrt(simDCHit.getMomentum()[0]*simDCHit.getMomentum()[0]+
                        simDCHit.getMomentum()[1]*simDCHit.getMomentum()[1]+
                        simDCHit.getMomentum()[2]*simDCHit.getMomentum()[2]);
                nloop++;
            }

            float cutMom = 0;
            for(int i=0;i<nloop;i++){
                if(dchitMom[i]>1e-1){
                    cutMom = dchitMom[i];
                    break;
                }
            }
            float cutMomLast = 0;
            for(int i=nloop-1;i>-1;i--){
                if(dchitMom[i]>1e-9){
                    cutMomLast = dchitMom[i];
                    break;
                }
            }
            if(((std::abs(mcP_Mom[0]-cutMom))/mcP_Mom[0])>m_cutMomHit[m_debugPID])
                return StatusCode::SUCCESS;
            if((std::abs(cutMom-cutMomLast)/mcP_Mom[0])>m_cutMomHit[m_debugPID])
                return StatusCode::SUCCESS;
        }// end loop of use cutMom
    }//end loop of useDC

    ///Retrieve silicon Track
    const edm4hep::TrackCollection* siTrackCol=nullptr;
    siTrackCol=m_siSubsetTrackCol.get();
    if(nullptr==siTrackCol){
        debug()<<"SDTTrackCollection is empty"<<endmsg;
        if(!m_useTruthTrack.value()) return StatusCode::SUCCESS;
    }else{
        debug()<<"SiSubsetTrackCol size "<<siTrackCol->size()<<endmsg;
        if(!m_useTruthTrack.value()&&0==siTrackCol->size()){
            return StatusCode::SUCCESS;
        }
        for(int iSitrk=0;iSitrk<siTrackCol->size();iSitrk++){
            auto siTk = siTrackCol->at(iSitrk);
            if(nHotsOnTrack(siTk,lcio::ILDDetID::VXD)<1e-9 &&
                    nHotsOnTrack(siTk,lcio::ILDDetID::SIT)<1e-9){
                debug() << " Because siTrack have no VXDHits and SITHits !" << endmsg;
                return StatusCode::SUCCESS;
            }
            if(nHotsOnTrack(siTk,lcio::ILDDetID::FTD)>1e-9){
                debug() << " Because siTrack have FTDHits !" << endmsg;
                return StatusCode::SUCCESS;
            }

        }
    }//end loop of silicon Track

    m_nSDTTrack = siTrackCol->size();
    for(int iSitrk=0;iSitrk<siTrackCol->size();iSitrk++){
        auto siTk = siTrackCol->at(iSitrk);
        ///New SDT track
        edm4hep::MutableTrack sdtTk=sdtTkCol->create();

        int nVXDHit=0;
        int nSITHit=0;
        int nSETHit=0;
        int nFTDHit=0;
        int nDCHitDCTk=0;
        int nDCHitSDTTk=0;

        ///Create track with mcParticle
        edm4hep::TrackerHit trackHit = siTk.getTrackerHits(0);
        edm4hep::MCParticle mcParticle;
        edm4hep::SimTrackerHit simTrackerHit;
        CEPC::getAssoMCParticle(VXDAssoVec,trackHit,mcParticle,simTrackerHit);
        edm4hep::TrackState trackStateMc;
        getTrackStateFromMcParticle(mcParticle,trackStateMc);

        if(m_useTruthTrack.value()||!m_useSi){
            sdtTk.addToTrackStates(trackStateMc);
        }

        if(m_useSi){
            if(!m_useTruthTrack.value()){
                debug()<<"siTk: "<<siTk<<endmsg;
                edm4hep::TrackState siTrackStat=siTk.getTrackStates(0);//FIXME?
                float Bz = 3;//T

                HelixClass helixClass;
                helixClass.Initialize_Canonical(siTrackStat.phi,siTrackStat.D0,
                        siTrackStat.Z0,siTrackStat.omega,siTrackStat.tanLambda,Bz);
                m_siXc[iSitrk] = helixClass.getXC();
                m_siYc[iSitrk] = helixClass.getYC();
                m_siR[iSitrk] = helixClass.getRadius();

                TVector3 siPos,siMom;
                TMatrixDSym siCov(6);
                double charge_test;

                CEPC::getPosMomFromTrackState(siTrackStat,Bz,siPos,siMom,charge_test,siCov);

                if((0==mcParticle.getCharge()) || (charge_test==mcParticle.getCharge())){
                    sdtTk.addToTrackStates(siTrackStat);
                }else{
                    sdtTk.addToTrackStates(trackStateMc);
                }

                sdtTk.setType(siTk.getType());
                sdtTk.setChi2(siTk.getChi2());
                sdtTk.setNdf(siTk.getNdf());
                sdtTk.setDEdx(siTk.getDEdx());
                sdtTk.setDEdxError(siTk.getDEdxError());
                sdtTk.setRadiusOfInnermostHit(
                        siTk.getRadiusOfInnermostHit());
            }
            if(!m_useSiTruthHit){
                debug()<<"use Si hit on track"<<endmsg;
                nVXDHit=addHotsToTk(siTk,sdtTk,lcio::ILDDetID::VXD,"VXD",nVXDHit);
                nSITHit=addHotsToTk(siTk,sdtTk,lcio::ILDDetID::SIT,"SIT",nSITHit);
                //end of loop over hits on siTk
            }else{
                if(m_useSiSimHit){
                    ///Add silicon SimTrackerHit
                    debug()<<"Add silicon SimTrackerHit"<<endmsg;
                    nVXDHit=addSimHitsToTk(m_VXDCollection,truthTrackerHitCol,sdtTk,"VXD",nVXDHit);
                    nSITHit=addSimHitsToTk(m_SITCollection,truthTrackerHitCol,sdtTk,"SIT",nSITHit);
                }else{  
                    ///Add reconstructed hit or digi
                    debug()<<"Add VXD TrackerHit"<<endmsg;
                    nVXDHit=addHitsToTk(m_VXDTrackerHits,sdtTk,"VXD digi",nVXDHit);
                    nSITHit=addHitsToTk(m_SITTrackerHits,sdtTk,"SIT digi",nSITHit);
                    if(m_useSiSpacePoint.value()){
                        ///Add silicon SpacePoint
                        debug()<<"Add silicon SpacePoint"<<endmsg;
                        if(m_useSiSpacePoint){
                            nSITHit=addHitsToTk(m_SITSpacePointCol,sdtTk,"SIT sp",nSITHit);
                        }
                    }//end of use space point
                } // end of Add reconstructed hit or digi
            } // end of truth hits
        }//end of use silicon

        if(m_useDC){
            ///Create DC Track
            edm4hep::MutableTrack dcTrack=dcTrackCol->create();

            //Create TrackState
            edm4hep::TrackState trackStateFirstDCHit;
            float charge=trackStateMc.omega/fabs(trackStateMc.omega);
            if(m_useFirstHitForDC&&getTrackStateFirstHit(m_DCSimTrackerHitCol,
                        charge,trackStateFirstDCHit)){
                dcTrack.addToTrackStates(trackStateFirstDCHit);
                dcTrack.addToTrackStates(trackStateMc);
            }else{
                dcTrack.addToTrackStates(trackStateMc);
                dcTrack.addToTrackStates(trackStateFirstDCHit);
            }
            sdtTk.addToTrackStates(trackStateFirstDCHit);

            ///Add other track properties
            dcTrack.setNdf(dcTrack.trackerHits_size()-5);

            if(!m_useTrackFinding){
                if(m_useIdealHit){
                    nDCHitDCTk=addIdealHitsToTk(m_DCDigiCol,truthTrackerHitCol,dcTrack,
                            "DC digi",nDCHitDCTk);
                }else{
                    nDCHitDCTk=addHitsToTk(m_DCDigiCol,dcTrack,"DC digi(dcTrack)",nDCHitDCTk);
                }
                debug()<<"dcTrack nHit "<<dcTrack.trackerHits_size()<<dcTrack<<endmsg;
            }
        }

        debug()<<"sdtTk nHit "<<sdtTk.trackerHits_size()<<sdtTk<<endmsg;
        debug()<<"nVXDHit "<<nVXDHit<<" nSITHit "<<nSITHit<<" nSETHit "<<nSETHit
            <<" nFTDHit "<<nFTDHit<<" nDCHitSDTTk "<<nDCHitSDTTk<<endmsg;
        std::cout<<"iSitrk = " << iSitrk << " nVXDHit "<<nVXDHit<<" nSITHit "<<nSITHit<<" nSETHit "<<nSETHit
            <<" nFTDHit "<<nFTDHit<<" nDCHitSDTTk "<<nDCHitSDTTk<<std::endl;;

        if(m_tuple){
            m_nHitOnSdtTkVXD[iSitrk]=nVXDHit;
            m_nHitOnSdtTkSIT[iSitrk]=nSITHit;
            m_nHitOnSdtTkDC[iSitrk]=nDCHitSDTTk;
            m_nHitOnSdtTk[iSitrk]=sdtTk.trackerHits_size();
            StatusCode sc=m_tuple->write();
        }
    }//end of loop over siTk

    debugEvent();
    StatusCode sc=m_tuple->write();

    return StatusCode::SUCCESS;
}

StatusCode TruthTrackerAlg::finalize()
{
    return GaudiAlgorithm::finalize();
}

void TruthTrackerAlg::getTrackStateFromMcParticle(
        edm4hep::MCParticle mcParticle,
        edm4hep::TrackState& trackState)
{
    ///Convert MCParticle to DC Track and ReconstructedParticle
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
    const edm4hep::Vector3f mcParticleMom=mcParticle.getMomentum();//GeV
    float mcParticlePt=sqrt(mcParticleMom.x*mcParticleMom.x+
            mcParticleMom.y*mcParticleMom.y);
    //float mcParticlePtSmeared=
    //    CLHEP::RandGauss::shoot(mcParticlePt,m_resPT);
    float mcParticleMomPhi=atan2(mcParticleMom.y,mcParticleMom.x);
    float mcParticleMomPhiSmeared=
        CLHEP::RandGauss::shoot(mcParticleMomPhi,m_resMomPhi);
    edm4hep::Vector3f mcParticleMomSmeared;
    mcParticleMomSmeared.x=mcParticlePt*cos(mcParticleMomPhiSmeared);
    mcParticleMomSmeared.y=mcParticlePt*sin(mcParticleMomPhiSmeared);
    mcParticleMomSmeared.z=CLHEP::RandGauss::shoot(mcParticleMom.z,m_resPz);

    ///Converted to Helix
    double B[3]={1e9,1e9,1e9};
    m_dd4hepField.magneticField({0.,0.,0.},B);
    HelixClass helix;
    ////FIXME DEBUG
    float pos[3]={mcParticleVertex.x,mcParticleVertex.y,mcParticleVertex.z};//mm
    float mom[3]={mcParticleMom.x,mcParticleMom.y,mcParticleMom.z};//mm
    helix.Initialize_VP(pos,mom,mcParticle.getCharge(),B[2]/dd4hep::tesla);

    ///new Track
    trackState.D0=helix.getD0();
    trackState.phi=helix.getPhi0();
    trackState.omega=helix.getOmega();
    trackState.Z0=helix.getZ0();
    trackState.tanLambda=helix.getTanLambda();
    trackState.referencePoint=helix.getReferencePoint();
    decltype(trackState.covMatrix) covMatrix;
    for(int i=0;i<covMatrix.size();i++){covMatrix[i]=999.;}//FIXME
    trackState.covMatrix=covMatrix;

    getCircleFromPosMom(pos,mom,B[2]/dd4hep::tesla,mcParticle.getCharge(),m_helixRadius,m_helixXC,m_helixYC);

    debug()<<"dd4hep::mm "<<dd4hep::mm<<" dd4hep::cm "<<dd4hep::cm<<endmsg;
    debug()<<"mcParticle "<<mcParticle
        <<" helix radius "<<helix.getRadius()<<" "<<helix.getXC()<<" "
        <<helix.getYC()<<" mm "
        <<" myhelix radius "<<m_helixRadius<<" "<<m_helixXC<<" "
        <<m_helixYC<<" mm "
        <<" momMC "<<mom[0]<<" "<<mom[1]<<" "<<mom[2]<<"GeV"
        <<" posMC "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<"mm"
        <<" momPhi "<<mcParticleMomPhi
        <<" mcParticleVertex("<<mcParticleVertex<<")mm "
        <<" mcParticleVertexSmeared("<<mcParticleVertexSmeared<<")mm "
        <<" mcParticleMom("<<mcParticleMom<<")GeV "
        <<" mcParticleMomSmeared("<<mcParticleMomSmeared<<")GeV "
        <<" Bxyz "<<B[0]/dd4hep::tesla<<" "<<B[1]/dd4hep::tesla
        <<" "<<B[2]/dd4hep::tesla<<" tesla"<<endmsg;
}//end of getTrackStateFromMcParticle

bool TruthTrackerAlg::getTrackStateFirstHit(
        DataHandle<edm4hep::SimTrackerHitCollection>& dcSimTrackerHitCol,
        float charge,edm4hep::TrackState& trackState)
{

    const edm4hep::SimTrackerHitCollection* col=nullptr;
    col=dcSimTrackerHitCol.get();
    debug()<<"TruthTrackerAlg::getTrackStateFirstHit"<<endmsg;
    debug()<<"simTrackerHitCol size "<<col->size()<<endmsg;
    float minHitTime=1e9;
    if(nullptr!=col||0==col->size()){
        edm4hep::SimTrackerHit firstHit;
        for(auto dcSimTrackerHit:*col){
            const edm4hep::Vector3f mom=dcSimTrackerHit.getMomentum();
            if(abs(sqrt(mom[0]*mom[0]+mom[1]*mom[1]))<m_momentumCut)continue;//yzhang TEMP skip hits with momentum <0.5GeV/c
            if(dcSimTrackerHit.getTime()<minHitTime) {
                minHitTime=dcSimTrackerHit.getTime();
                firstHit=dcSimTrackerHit;
            }
        }
        const edm4hep::Vector3d pos=firstHit.getPosition();
        const edm4hep::Vector3f mom=firstHit.getMomentum();
        debug()<<"first Hit pos "<<pos<<" mom "<<mom<<" time "<<minHitTime<<endmsg;
        float pos_t[3]={(float)pos[0],(float)pos[1],(float)pos[2]};
        float mom_t[3]={(float)mom[0],(float)mom[1],(float)mom[2]};
        ///Converted to Helix
        double B[3]={1e9,1e9,1e9};
        m_dd4hepField.magneticField({0.,0.,0.},B);
        HelixClass helix;
        helix.Initialize_VP(pos_t,mom_t,charge,B[2]/dd4hep::tesla);
        m_helixRadiusFirst=helix.getRadius();
        m_helixXCFirst=helix.getXC();
        m_helixYCFirst=helix.getYC();

        ///new Track
        trackState.D0=helix.getD0();
        trackState.phi=helix.getPhi0();
        trackState.omega=helix.getOmega();
        trackState.Z0=helix.getZ0();
        trackState.tanLambda=helix.getTanLambda();
        trackState.referencePoint=helix.getReferencePoint();
        std::array<float,21> covMatrix;
        for(int i=0;i<21;i++){covMatrix[i]=100.;}//FIXME
        trackState.covMatrix=covMatrix;
        debug()<<"first hit trackState "<<trackState<<endmsg;
        return true;
    }
    return false;
}//end of getTrackStateFirstHit

void TruthTrackerAlg::debugEvent()
{
    if(m_useSi){
        ///Retrieve silicon Track
        const edm4hep::TrackCollection* siTrackCol=nullptr;
        siTrackCol=m_siSubsetTrackCol.get();
        int nSiTrk = 0;
        if(nullptr!=siTrackCol){
            for(auto siTk:*siTrackCol){
                debug()<<"siTk: "<<siTk<<endmsg;
                edm4hep::TrackState trackStat=siTk.getTrackStates(0);//FIXME?
                double B[3]={1e9,1e9,1e9};
                m_dd4hepField.magneticField({0.,0.,0.},B);
                HelixClass helix;
                helix.Initialize_Canonical(trackStat.phi, trackStat.D0, trackStat.Z0,
                        trackStat.omega, trackStat.tanLambda, B[2]/dd4hep::tesla);

                if(m_tuple){
                    m_nHitOnSiTkVXD[nSiTrk]=nHotsOnTrack(siTk,lcio::ILDDetID::VXD);
                    m_nHitOnSiTkSIT[nSiTrk]=nHotsOnTrack(siTk,lcio::ILDDetID::SIT);
                    m_nHitOnSiTkSET[nSiTrk]=nHotsOnTrack(siTk,lcio::ILDDetID::SET);
                    m_nHitOnSiTkFTD[nSiTrk]=nHotsOnTrack(siTk,lcio::ILDDetID::FTD);
                }
                nSiTrk++;
            }//end of loop over siTk
        }
        if(m_tuple){
            //SimTrackerHits
            m_nSimTrackerHitVXD=simTrackerHitColSize(m_VXDCollection);
            m_nSimTrackerHitSIT=simTrackerHitColSize(m_SITCollection);
            m_nSimTrackerHitSET=simTrackerHitColSize(m_SETCollection);
            m_nSimTrackerHitFTD=simTrackerHitColSize(m_FTDCollection);

            //TrackerHits
            m_nTrackerHitVXD=trackerHitColSize(m_VXDTrackerHits);
            m_nTrackerHitSIT=trackerHitColSize(m_SITTrackerHits);
            m_nTrackerHitSET=trackerHitColSize(m_SETTrackerHits);
            m_nTrackerHitFTD=trackerHitColSize(m_FTDTrackerHits);
            m_nTrackerHitDC=trackerHitColSize(m_DCDigiCol);
        }
    }
}

int TruthTrackerAlg::addIdealHitsToTk(
        DataHandle<edm4hep::TrackerHitCollection>& colHandle,
        edm4hep::TrackerHitCollection*& truthTrackerHitCol,
        edm4hep::MutableTrack& track, const char* msg,int nHitAdded)
{
    if(nHitAdded>0) return nHitAdded;
    int nHit=0;
    const edm4hep::TrackerHitCollection* col=colHandle.get();
    debug()<<"add "<<msg<<" "<<col->size()<<" trackerHit"<<endmsg;
    debug()<<track<<endmsg;
    for(auto hit:*col){
        //get end point of this wire
        TVector3 endPointStart(0,0,0);
        TVector3 endPointEnd(0,0,0);
        m_gridDriftChamber->cellposition(hit.getCellID(),endPointStart,
                endPointEnd);//cm

        //calc. doca of helix to wire
        TVector3 wire(endPointStart.X()/dd4hep::mm,endPointStart.Y()/dd4hep::mm,0);//to mm
        TVector3 center(m_helixXC,m_helixYC,0);//mm
        float docaIdeal=(center-wire).Mag()-m_helixRadius;//mm
        TVector3 centerFirst(m_helixXCFirst,m_helixYCFirst,0);//mm
        float docaIdealFirst=(centerFirst-wire).Mag()-m_helixRadiusFirst;//mm

        //add modified hit
        auto tmpHit = truthTrackerHitCol->create();
        tmpHit=hit.clone();
        //tmpHit=hit;
        tmpHit.setTime(fabs(docaIdeal)*1e3/40.);//40#um/ns, drift time in ns
        track.addToTrackerHits(tmpHit);

        long long int detID=hit.getCellID();
        debug()<<" addIdealHitsToTk "<<m_helixRadius<<" center "<<m_helixXC
            <<" "<<m_helixYC<<" mm wire("<<m_decoder->get(detID,"layer")<<","
            <<m_decoder->get(detID,"cellID")<<") "<<wire.X()<<" "<<wire.Y()
            <<"mm docaIdeal "<<docaIdeal<<" docaIdealFirst "<<docaIdealFirst<<"mm "
            <<"hit.Time orignal "<<hit.getTime()<<" new Time "
            <<fabs(docaIdeal)*1e3/40.<<endmsg;
        ++nHit;
    }
    return nHit;
}

int TruthTrackerAlg::addHitsToTk(DataHandle<edm4hep::TrackerHitCollection>&
colHandle, edm4hep::MutableTrack& track, const char* msg,int nHitAdded)
{
    if(nHitAdded>0) return nHitAdded;
    int nHit=0;
    const edm4hep::TrackerHitCollection* col=colHandle.get();
    debug()<<"add "<<msg<<" "<<col->size()<<" trackerHit"<<endmsg;
    //sort,FIXME
    for(auto hit:*col){
        track.addToTrackerHits(hit);
        ++nHit;
    }
    return nHit;
}

int TruthTrackerAlg::addOneLayerHitsToTk(DataHandle<edm4hep::TrackerHitCollection>&
colHandle, edm4hep::MutableTrack& track, const char* msg,int nHitAdded)
{
    if(nHitAdded>0) return nHitAdded;
    int nHit=0;
    const edm4hep::TrackerHitCollection* col=colHandle.get();
    debug()<<"add "<<msg<<" "<<col->size()<<" trackerHit"<<endmsg;
    //sort,FIXME
    bool isHaveHit[55]={0};
    for(auto hit:*col){
        unsigned long long cellID=hit.getCellID();
        int layer = m_decoder->get(cellID,"layer");
        if(!(isHaveHit[layer])){
            track.addToTrackerHits(hit);
            ++nHit;
        }
        isHaveHit[layer] = 1;
    }
    return nHit;
}

void TruthTrackerAlg::getAssoMCParticle(
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        edm4hep::TrackerHit trackerHit,
        edm4hep::MCParticle& mcParticle) const
{
    edm4hep::SimTrackerHit simTrackerHit;
    for(auto assoHit: *assoHits){
        if(assoHit.getRec()==trackerHit)
        {
            simTrackerHit=assoHit.getSim();
            mcParticle = simTrackerHit.getMCParticle();

        }
    }
}

edm4hep::SimTrackerHit TruthTrackerAlg::getAssoSimTrackerHit(
        const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        edm4hep::TrackerHit trackerHit) const
{
    edm4hep::SimTrackerHit simTrackerHit;
    for(auto assoHit: *assoHits){
        if(assoHit.getRec()==trackerHit)
        {
            simTrackerHit=assoHit.getSim();
        }
    }
    return simTrackerHit;
}

int TruthTrackerAlg::getSiMCParticle(edm4hep::Track siTrack,
        const edm4hep::TrackerHitCollection* dcTrackerHits,
        const edm4hep::MCRecoTrackerAssociationCollection* vxdAsso,
        const edm4hep::MCRecoTrackerAssociationCollection* sitAsso,
        const edm4hep::MCRecoTrackerAssociationCollection* setAsso,
        const edm4hep::MCRecoTrackerAssociationCollection* ftdAsso,
        const edm4hep::MCRecoTrackerAssociationCollection* dcAsso
        )
{
    int nDCHit =0;

    edm4hep::TrackerHit trackerHit = siTrack.getTrackerHits(0);
    unsigned long long cellID = trackerHit.getCellID();
    UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string) ;
    encoder.setValue(cellID);
    int detID=encoder[lcio::ILDCellID0::subdet] ;
    edm4hep::MCParticle mcParticle;
    if(lcio::ILDDetID::VXD==detID){
        std::cout << "the hit from VXD!!" << std::endl;
        getAssoMCParticle(vxdAsso,trackerHit,mcParticle);
    }else if(lcio::ILDDetID::SIT==detID){
        getAssoMCParticle(sitAsso,trackerHit,mcParticle);
        std::cout << "the hit from SIT!!" << std::endl;
    }else if(lcio::ILDDetID::SET==detID){
        getAssoMCParticle(setAsso,trackerHit,mcParticle);
        std::cout << "the hit from SET!!" << std::endl;
    }else if(lcio::ILDDetID::FTD==detID){
        getAssoMCParticle(ftdAsso,trackerHit,mcParticle);
        std::cout << "the hit from FTD!!" << std::endl;
    }else{
        std::cout << "the hit from ERROR!!" << std::endl;
    }

    for(auto dcHit:*dcTrackerHits){

        edm4hep::SimTrackerHit simDCHit = getAssoSimTrackerHit(dcAsso,dcHit);

        if(mcParticle==simDCHit.getMCParticle()) nDCHit++;
    }
    return nDCHit;
}


int TruthTrackerAlg::addSimHitsToTk(
        DataHandle<edm4hep::SimTrackerHitCollection>& colHandle,
        edm4hep::TrackerHitCollection*& truthTrackerHitCol,
        //edm4hep::Track& track, const char* msg,int nHitAdded)
        edm4hep::MutableTrack& track, const char* msg,int nHitAdded)
{
    if(nHitAdded>0) return nHitAdded;
    int nHit=0;
    const edm4hep::SimTrackerHitCollection* col=colHandle.get();
    for(auto simTrackerHit:*col){
        auto trackerHit=truthTrackerHitCol->create();
        if(m_skipSecondaryHit&&simTrackerHit.isProducedBySecondary()) {
            debug()<<"skip secondary simTrackerHit "<<msg<<endmsg;
            continue;
        }
        auto& pos = simTrackerHit.getPosition();
        debug()<<" addSimHitsToTk "<<msg<<" "<<sqrt(pos.x*pos.x+pos.y*pos.y)<<endmsg;
        UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string) ;
        int detID=encoder[lcio::ILDCellID0::subdet] ;
        float resolution[3];
        if(lcio::ILDDetID::VXD==detID){
            for(int i=0;i<3;i++)resolution[i]=m_resVXD[i];
        }else if(lcio::ILDDetID::SIT==detID){
            for(int i=0;i<3;i++)resolution[i]=m_resSIT[i];
        }else if(lcio::ILDDetID::SET==detID){
            for(int i=0;i<3;i++)resolution[i]=m_resSET[i];
        }else if(lcio::ILDDetID::FTD==detID){
            for(int i=0;i<3;i++)resolution[i]=m_resFTDPixel[i];//FIXME
        }else{
            for(int i=0;i<3;i++)resolution[i]=0.003;
        }
        edm4hep::Vector3d posSmeared;//mm
        posSmeared.x=CLHEP::RandGauss::shoot(pos.x,resolution[0]);
        posSmeared.y=CLHEP::RandGauss::shoot(pos.y,resolution[1]);
        posSmeared.z=CLHEP::RandGauss::shoot(pos.z,resolution[2]);
        trackerHit.setPosition(posSmeared) ;
        encoder.setValue(simTrackerHit.getCellID()) ;
        trackerHit.setCellID(encoder.lowWord());//?FIXME
        std::array<float, 6> cov;
        cov[0]=resolution[0]*resolution[0];
        cov[1]=0.;
        cov[2]=resolution[1]*resolution[1];
        cov[3]=0.;
        cov[4]=0.;
        cov[5]=resolution[2]*resolution[2];
        trackerHit.setCovMatrix(cov);
        debug()<<"add simTrackerHit "<<msg<<" trackerHit "<<trackerHit<<endmsg;
        ///Add hit to track
        track.addToTrackerHits(trackerHit);
        trackerHit.setEDep(simTrackerHit.getEDep());
        trackerHit.addToRawHits(simTrackerHit.getObjectID());
        trackerHit.setType(-8);//FIXME?
        ++nHit;
    }
    debug()<<"add simTrackerHit "<<msg<<" "<<nHit<<endmsg;
    return nHit;
}

bool TruthTrackerAlg::debugVertex(edm4hep::Track sourceTrack,
        const edm4hep::MCRecoTrackerAssociationCollection* vxdAsso,
        const edm4hep::MCRecoTrackerAssociationCollection* sitAsso,
        const edm4hep::MCRecoTrackerAssociationCollection* setAsso,
        const edm4hep::MCRecoTrackerAssociationCollection* ftdAsso
        ){
    for(unsigned int iHit=0;iHit<sourceTrack.trackerHits_size();iHit++){
        edm4hep::TrackerHit hit=sourceTrack.getTrackerHits(iHit);
        edm4hep::Vector3d position = hit.getPosition();

        unsigned long long cellID = hit.getCellID();
        UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string) ;
        encoder.setValue(cellID);

        int detID=encoder[lcio::ILDCellID0::subdet] ;

        edm4hep::SimTrackerHit simTrackerHit;
        edm4hep::MCParticle mcParticle;
        if(lcio::ILDDetID::VXD==detID){
            simTrackerHit = getAssoSimTrackerHit(vxdAsso,hit);
            getAssoMCParticle(vxdAsso,hit,mcParticle);
        }else if(lcio::ILDDetID::SIT==detID){
            simTrackerHit = getAssoSimTrackerHit(sitAsso,hit);
            getAssoMCParticle(sitAsso,hit,mcParticle);
        }else if(lcio::ILDDetID::SET==detID){
            simTrackerHit = getAssoSimTrackerHit(setAsso,hit);
            getAssoMCParticle(setAsso,hit,mcParticle);
        }else if(lcio::ILDDetID::FTD==detID){
            simTrackerHit = getAssoSimTrackerHit(ftdAsso,hit);
            getAssoMCParticle(ftdAsso,hit,mcParticle);
        }else{
            std::cout << __LINE__ <<"the hit from ERROR!!" << std::endl;
        }

        edm4hep::Vector3f momentum = simTrackerHit.getMomentum();

        float pos[3] = {position[0],position[1],position[2]};
        float mom[3] = {momentum[0]*1e-1,momentum[1]*1e-1,momentum[2]*1e-1};

        float mcPos[3] = {mcParticle.getVertex().x,mcParticle.getVertex().y,mcParticle.getVertex().z};
        float mcMom[3] = {mcParticle.getMomentum().x,mcParticle.getMomentum().y,mcParticle.getMomentum().z};

        HelixClass helix;
        helix.Initialize_VP(pos,mom,mcParticle.getCharge(),3);

        HelixClass mcHelix;
        mcHelix.Initialize_VP(mcPos,mcMom,mcParticle.getCharge(),3);

        float deltaD0 = helix.getD0()-mcHelix.getD0();
        float deltaZ0 = helix.getZ0()-mcHelix.getZ0();
        if((1000*deltaD0)>30 || (1000*deltaZ0)>30){ //um
            return true;
            break;
        }
    }
    return false;
}

int TruthTrackerAlg::addHotsToTk(edm4hep::Track& sourceTrack,
    edm4hep::MutableTrack& targetTrack, int hitType,const char* msg,int nHitAdded)
{
    int nHit=0;
    for(unsigned int iHit=0;iHit<sourceTrack.trackerHits_size();iHit++){
        edm4hep::TrackerHit hit=sourceTrack.getTrackerHits(iHit);
        UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
        encoder.setValue(hit.getCellID());
        if(encoder[lcio::ILDCellID0::subdet]==hitType){
            targetTrack.addToTrackerHits(hit);
            debug()<<endmsg<<" add siHit "<<msg<<" "<<iHit<<" "<<hit
                <<" pos "<<hit.getPosition().x<<" "<<hit.getPosition().y<<" "
                <<hit.getPosition().z<<" " <<endmsg;
            ++nHit;
        }
    }
    debug()<<endmsg<<" "<<nHit<<" "<<msg<<" hits add on track"<<endmsg;
    return nHit;
}

int TruthTrackerAlg::nHotsOnTrack(edm4hep::Track& track, int hitType)
{
    int nHit=0;
    for(unsigned int iHit=0;iHit<track.trackerHits_size();iHit++){
        edm4hep::TrackerHit hit=track.getTrackerHits(iHit);
        //edm4hep::ConstTrackerHit hit=track.getTrackerHits(iHit);
        UTIL::BitField64 encoder(lcio::ILDCellID0::encoder_string);
        encoder.setValue(hit.getCellID());
        if(encoder[lcio::ILDCellID0::subdet]==hitType){
            ++nHit;
        }
    }
    return nHit;
}

int TruthTrackerAlg::trackerHitColSize(DataHandle<edm4hep::TrackerHitCollection>& col)
{
    const edm4hep::TrackerHitCollection* c=col.get();
    if(nullptr!=c) return c->size();
    return 0;
}

int TruthTrackerAlg::simTrackerHitColSize(DataHandle<edm4hep::SimTrackerHitCollection>& col)
{
    const edm4hep::SimTrackerHitCollection* c=col.get();
    if(nullptr!=c) return c->size();
    return 0;
}
//unit length is mm
void TruthTrackerAlg::getCircleFromPosMom(float pos[3],float mom[3],
        float Bz,float q,float& helixRadius,float& helixXC,float& helixYC)
{
    float FCT = 2.99792458E-4;//mm
    float pxy = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
    helixRadius = pxy / (FCT*Bz);
    float phiMomRefPoint = atan2(mom[1],mom[0]);
    helixXC= pos[0] + helixRadius*cos(phiMomRefPoint-M_PI*0.5*q);
    helixYC= pos[1] + helixRadius*sin(phiMomRefPoint-M_PI*0.5*q);
}
