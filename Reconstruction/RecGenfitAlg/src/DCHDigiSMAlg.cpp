#include "DCHDigiSMAlg.h"

#include "edm4hep/Vector3f.h"

#include "GenfitTrack.h"
#include "GenfitFitter.h"
#include "GenfitField.h"

#include "DetInterface/IGeomSvc.h"
#include "DataHelper/HelixClass.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/Vector3D.h"

#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/MsgStream.h"


#include <array>
#include <math.h>
#include <cmath>
#include <algorithm>

DECLARE_COMPONENT(DCHDigiSMAlg)

DCHDigiSMAlg::DCHDigiSMAlg(const std::string& name, ISvcLocator* svcLoc)
 : GaudiAlgorithm(name, svcLoc),
   _nEvt(0)
{
    declareProperty("DCTrackCollection", m_dcTrackCol,
            "Handle of DC track collection");
    declareProperty("DigiDCHitCollection", m_digiDCHitsCol, "Handle of Digi DCHit collection");
    declareProperty("DCHitAssociationCollection", m_dcHitAssociationCol,
            "Handle of simTrackerHit and TrackerHit association collection");
    m_thisName = name;
}

StatusCode DCHDigiSMAlg::initialize()
{
  info() << "Booking Ntuple" << endmsg;

  m_geomSvc=Gaudi::svcLocator()->service("GeomSvc");
  if (nullptr==m_geomSvc) {
      std::cout<<"Failed to find GeomSvc"<<std::endl;
      return StatusCode::FAILURE;
  }

  m_dd4hepField=m_geomSvc->lcdd()->field();
  m_genfitField=new GenfitField(m_dd4hepField);
  dd4hep::Detector* m_dd4hep = m_geomSvc->lcdd();
  dd4hep::Readout readout = m_dd4hep->readout(m_readout_name);
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>(readout.segmentation().segmentation());
  m_decoder = m_geomSvc->getDecoder(m_readout_name);

  NTuplePtr nt1(ntupleSvc(), "DCHDigiSMAlg/Track");
  if ( nt1 ) {
      m_tuple = nt1;
  }else{
      m_tuple=ntupleSvc()->book("DCHDigiSMAlg/Track",CLID_ColumnWiseTuple,
              "Track");
  }
  if(m_tuple){
      StatusCode sc;
      sc=m_tuple->addItem ("ntrk",      m_nTracks, 0, 500000 );
      sc=m_tuple->addItem ("ntrkhit",      m_nTrackHits, 0, 500000 );
      sc=m_tuple->addItem("nDCDigi",m_nDCDigi,0,50000);
      sc=m_tuple->addItem("dcHitTime",m_nDCDigi,m_dcHitTime);
      sc=m_tuple->addItem("dcHitTime",m_nDCDigi,m_dcHitTime);
      sc=m_tuple->addItem("dcHitDoca",m_nDCDigi,m_dcHitDoca);
      sc=m_tuple->addItem("distance",m_nDCDigi,m_distance);
      sc=m_tuple->addIndexedItem ("trackhitx", m_nTrackHits, m_trackhitx );
      sc=m_tuple->addIndexedItem ("trackhity", m_nTrackHits, m_trackhity );
      sc=m_tuple->addIndexedItem ("trackhitz", m_nTrackHits, m_trackhitz );
      sc=m_tuple->addIndexedItem ("Doca", m_nTrackHits, m_Doca );
      sc=m_tuple->addIndexedItem ("Doca_digi", m_nDCDigi, m_Doca_digi );
      sc=m_tuple->addIndexedItem ("d0",        m_nTracks, m_d0 );
      sc=m_tuple->addIndexedItem ("phi0",      m_nTracks, m_phi0 );
      sc=m_tuple->addIndexedItem ("omega",     m_nTracks, m_omega );
      sc=m_tuple->addIndexedItem ("z0",        m_nTracks, m_z0 );
      sc=m_tuple->addIndexedItem ("tanLambda", m_nTracks, m_tanLambda );
      sc=m_tuple->addIndexedItem ("xc", m_nTracks, m_xc );
      sc=m_tuple->addIndexedItem ("yc",         m_nTracks, m_yc );
      sc=m_tuple->addIndexedItem ("R",         m_nTracks, m_R );
      sc=m_tuple->addIndexedItem ("px",        m_nTracks, m_px );
      sc=m_tuple->addIndexedItem ("py",        m_nTracks, m_py );
      sc=m_tuple->addIndexedItem ("pz",        m_nTracks, m_pz );
      sc=m_tuple->addIndexedItem ("pocax",        m_nTracks, m_pocax );
      sc=m_tuple->addIndexedItem ("pocay",        m_nTracks, m_pocay );
  }else{
      error()<< "Cannot book tuple MyTuples/Track" << endmsg;
  }

  _nEvt = 0;
  
  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode DCHDigiSMAlg::execute()
{
    info() << "Processing " << _nEvt << " events " << endmsg;
    StatusCode sc;

    const edm4hep::TrackCollection* dcTrackCol=nullptr;
    if(m_dcTrackCol.exist()) dcTrackCol=m_dcTrackCol.get();
    if(nullptr==dcTrackCol) {
        debug()<<"TrackCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }

    auto assoHits =m_dcHitAssociationCol.get();

    int iTracks=0;
    float ratio=0;
    TVector3 hitPosition(0,0,0);
    if(m_tuple) m_nTracks = dcTrackCol->size();
    for(auto dcTrack: *dcTrackCol){
       // for(unsigned int i=0;i < dcTrack.trackStates_size(); i++) {
        edm4hep::TrackState trackStat=dcTrack.getTrackStates(0);

        HelixClass helixClass;
        helixClass.Initialize_Canonical(trackStat.phi,trackStat.D0,
                trackStat.Z0,trackStat.omega,trackStat.tanLambda,
                m_genfitField->getBzTesla({0.,0.,0.}));
        m_d0        [iTracks]  = trackStat.D0;
        m_phi0      [iTracks]  = trackStat.phi;
        m_omega     [iTracks]  = trackStat.omega;
        m_z0        [iTracks]  = trackStat.Z0;
        m_tanLambda [iTracks]  = trackStat.tanLambda;
        m_xc [iTracks]= helixClass.getXC();
        m_yc [iTracks]= helixClass.getYC();
        m_R [iTracks]= helixClass.getRadius();
        TVector3 momInit(helixClass.getMomentum()[0],
                helixClass.getMomentum()[1],helixClass.getMomentum()[2]);

        m_px [iTracks]= momInit.X();
        m_py [iTracks]= momInit.Y();
        m_pz [iTracks]= momInit.Z();

        ratio = m_d0[iTracks]/(m_d0[iTracks]+m_R[iTracks]);
        m_pocax [iTracks] = ratio*m_xc[iTracks];
        m_pocay [iTracks] = ratio*m_yc[iTracks];

        const edm4hep::TrackerHitCollection* dCDigiCol=nullptr;
        dCDigiCol=m_digiDCHitsCol.get();
        if(nullptr!=dCDigiCol){ m_nDCDigi=dCDigiCol->size(); }
        int iDCDigi=0;
        for(auto dcDigi: *dCDigiCol){
            m_dcHitTime[iDCDigi]=dcDigi.getTime();
            m_dcHitDoca[iDCDigi]=dcDigi.getTime()*40./1000.; //mm

            unsigned long long wcellid = dcDigi.getCellID();


            TVector3 Wstart_digi(0,0,0);
            TVector3 Wend_digi  (0,0,0);
            m_segmentation->cellposition(wcellid, Wstart_digi, Wend_digi);

            float disx_digi = Wstart_digi.x()/dd4hep::mm-m_xc[iTracks];
            float disy_digi = Wstart_digi.y()/dd4hep::mm-m_yc[iTracks];
            float sqrta_digi = sqrt(disx_digi*disx_digi+disy_digi*disy_digi);
            m_Doca_digi[iDCDigi] = fabs(m_R[iTracks] - sqrta_digi);

            edm4hep::Vector3d pos=dcDigi.position();
            TVector3 p(pos.x,pos.y,pos.z);

//            for(int iSimHit=0;iSimHit<(int) assoHits->size();iSimHit++){
//                if(assoHits->at(iSimHit).getRec()== dcDigi) {
//std::cout << "size= " << iSimHit << std::endl;
            const edm4hep::ConstSimTrackerHit sim = assoHits->at(iDCDigi).getSim();
            edm4hep::Vector3f Mom = sim.getMomentum();
            TVector3 m(Mom.x,Mom.y,Mom.z);

            float Steplength = sim.getPathLength();
//std::cout << " step= " << Steplength
//          << " p[0] " << p[0]
//          << " p[1] " << p[1] 
//          << " p[2] " << p[2]
//          << std::endl;
            TVector3  pos_start = p - 0.5 * Steplength * m.Unit();
            TVector3  pos_end = p + 0.5 * Steplength * m.Unit();

            m_distance[iDCDigi] = (m_segmentation->Distance(wcellid,pos_start,pos_end,hitPosition))/dd4hep::mm;
            //                }
            //            }

            iDCDigi++;
        }


        int iTrackHits = 0;
        if(m_tuple) m_nTrackHits = m_nTracks + dcTrack.trackerHits_size();
        for(unsigned int i=0;i < dcTrack.trackerHits_size(); i++) {
            edm4hep::ConstTrackerHit hit = dcTrack.getTrackerHits(i);

            unsigned long long id = hit.getCellID();
            int chamber = m_decoder->get(id, "chamber");
            int layer   = m_decoder->get(id, "layer"  );
            int cellID  = m_decoder->get(id, "cellID" );

            TVector3 Wstart(0,0,0);
            TVector3 Wend  (0,0,0);
            m_segmentation->cellposition(id, Wstart, Wend);
            m_trackhitx [iTrackHits]= hit.getPosition()[0]/dd4hep::mm;
            m_trackhity [iTrackHits]= hit.getPosition()[1]/dd4hep::mm;
            m_trackhitz [iTrackHits]= hit.getPosition()[2]/dd4hep::mm;
            float disx = Wstart.x()/dd4hep::mm-m_xc[iTracks];
            float disy = Wstart.y()/dd4hep::mm-m_yc[iTracks];
            float sqrta = sqrt(disx*disx+disy*disy);
            m_Doca[iTrackHits] = fabs(m_R[iTracks] - sqrta);

            iTrackHits++;
        }
        iTracks++;
    }
    //    }
    if(m_tuple) sc=m_tuple->write();
    _nEvt++;

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode DCHDigiSMAlg::finalize()
{
    info() << "Processed " << _nEvt << " events " << endmsg;
    return GaudiAlgorithm::finalize();
}


