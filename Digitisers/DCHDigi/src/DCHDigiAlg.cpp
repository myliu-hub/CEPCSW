/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "DCHDigiAlg.h"


#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"

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

#include "time.h"
#include <unordered_map>

#include "DataHelper/HelixClass.h"

DECLARE_COMPONENT( DCHDigiAlg )

DCHDigiAlg::DCHDigiAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
  
  // Input collections
  declareProperty("SimDCHitCollection", r_SimDCHCol, "Handle of the Input SimHit collection");
  
  // Output collections
  declareProperty("DigiDCHitCollection", w_DigiDCHCol, "Handle of Digi DCHit collection");
  declareProperty("TruthDigiDCHitCollection", w_TruthDigiDCHCol, "Handle of Truth Digi DCHit collection");
  
  declareProperty("AssociationCollection", w_AssociationCol, "Handle of Association collection");
   
}

StatusCode DCHDigiAlg::initialize()
{
  m_geosvc = service<IGeomSvc>("GeomSvc");
  if ( !m_geosvc )  throw "DCHDigiAlg :Failed to find GeomSvc ...";
  dd4hep::Detector* m_dd4hep = m_geosvc->lcdd();
  if ( !m_dd4hep )  throw "DCHDigiAlg :Failed to get dd4hep::Detector ...";
  dd4hep::Readout readout = m_dd4hep->readout(m_readout_name);
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>(readout.segmentation().segmentation());
  m_decoder = m_geosvc->getDecoder(m_readout_name);
  if (!m_decoder) {
      error() << "Failed to get the decoder. " << endmsg;
      return StatusCode::FAILURE;
  }

  fRandom.SetSeed(105105);//FIXME: set by users

  if(m_WriteAna){

      NTuplePtr nt( ntupleSvc(), "MyTuples/DCH_digi_evt" );
      if ( nt ) m_tuple = nt;
      else {
          m_tuple = ntupleSvc()->book( "MyTuples/DCH_digi_evt", CLID_ColumnWiseTuple, "DCH_digi_evt" );
          if ( m_tuple ) {
            m_tuple->addItem( "N_evt" , m_evt ).ignore();
            m_tuple->addItem( "N_sim" , m_n_sim , 0, 1000000 ).ignore();
            m_tuple->addItem( "N_digi", m_n_digi, 0, 1000000 ).ignore();
            m_tuple->addItem( "run_time" , m_time ).ignore();
            m_tuple->addItem( "simhit_x", m_n_sim, m_simhit_x).ignore();
            m_tuple->addItem( "simhit_y", m_n_sim, m_simhit_y).ignore();
            m_tuple->addItem( "simhit_z", m_n_sim, m_simhit_z).ignore();
            m_tuple->addItem( "Simdca"      , m_n_sim, m_Simdca       ).ignore();
            m_tuple->addItem( "simhitT", m_n_sim, m_simhitT    ).ignore();
            m_tuple->addItem( "simhitmom", m_n_sim, m_simhitmom    ).ignore();
            m_tuple->addItem( "simPDG", m_n_sim, m_simPDG    ).ignore();
            m_tuple->addItem( "chamber" , m_n_digi, m_chamber  ).ignore();
            m_tuple->addItem( "layer"   , m_n_digi, m_layer    ).ignore();
            m_tuple->addItem( "cell"    , m_n_digi, m_cell     ).ignore();
            m_tuple->addItem( "cell_x"  , m_n_digi, m_cell_x   ).ignore();
            m_tuple->addItem( "cell_y"  , m_n_digi, m_cell_y   ).ignore();
            m_tuple->addItem( "cell1_x"  , m_n_digi, m_cell1_x   ).ignore();
            m_tuple->addItem( "cell1_y"  , m_n_digi, m_cell1_y   ).ignore();
            m_tuple->addItem( "hit_x"    , m_n_digi,m_hit_x     ).ignore();
            m_tuple->addItem( "hit_y"    , m_n_digi,m_hit_y     ).ignore();
            m_tuple->addItem( "hit_z"    , m_n_digi,m_hit_z     ).ignore();
            m_tuple->addItem( "mom_x"    , m_n_digi,m_mom_x     ).ignore();
            m_tuple->addItem( "mom_y"    , m_n_digi,m_mom_y     ).ignore();
            m_tuple->addItem( "dca"      , m_n_digi, m_dca       ).ignore();
            m_tuple->addItem( "poca_x"   , m_n_digi, m_poca_x    ).ignore();
            m_tuple->addItem( "poca_y"   , m_n_digi, m_poca_y    ).ignore();
            m_tuple->addItem( "hit_dE"   , m_n_digi,m_hit_dE    ).ignore();
            m_tuple->addItem( "truthlength", m_n_digi,m_truthlength).ignore();
            m_tuple->addItem( "N_noise", m_n_noise, 0, 1000000 ).ignore();
            m_tuple->addItem( "layer_noise", m_n_noise,m_noise_layer).ignore();
            m_tuple->addItem( "cell_noise", m_n_noise,m_noise_cell).ignore();
            m_tuple->addItem( "NoiseT", m_n_noise,m_noise_time).ignore();
            m_tuple->addItem( "N_noiseCover", m_n_noiseCover, 0, 1000000 ).ignore();
            m_tuple->addItem( "layer_noiseCover", m_n_noiseCover,m_noiseCover_layer).ignore();
            m_tuple->addItem( "cell_noiseCover", m_n_noiseCover,m_noiseCover_cell).ignore();
          } else { // did not manage to book the N tuple....
            info() << "    Cannot book N-tuple:" << long( m_tuple ) << endmsg;
          }
      }
  }
  info()<<"DCHDigiAlg::initialized"<<endmsg;
  return GaudiAlgorithm::initialize();
}


StatusCode DCHDigiAlg::execute()
{
  m_start = clock();

  info() << "Processing " << _nEvt << " events " << endmsg;
  if(m_WriteAna) m_evt = _nEvt;
  edm4hep::TrackerHitCollection* Vec   = w_DigiDCHCol.createAndPut();
  edm4hep::TrackerHitCollection* TruthVec   = w_TruthDigiDCHCol.createAndPut();
  edm4hep::MCRecoTrackerAssociationCollection* AssoVec   = w_AssociationCol.createAndPut();
  const edm4hep::SimTrackerHitCollection* SimHitCol =  r_SimDCHCol.get();
  if (SimHitCol->size() == 0) {
    return StatusCode::SUCCESS;
  }
  debug()<<"input sim hit size="<< SimHitCol->size() <<endmsg;

  auto SimHit0 = SimHitCol->at(0);
  std::map<unsigned long long, std::vector<decltype(SimHit0)>> id_hits_map;


  for( int i = 0; i < SimHitCol->size(); i++ ) 
  {
      auto SimHit = SimHitCol->at(i);
      unsigned long long id = SimHit.getCellID();
      float sim_hit_mom = sqrt( SimHit.getMomentum()[0]*SimHit.getMomentum()[0] + SimHit.getMomentum()[1]*SimHit.getMomentum()[1] + SimHit.getMomentum()[2]*SimHit.getMomentum()[2] );//GeV
      if(sim_hit_mom < m_mom_threshold) continue; 
      if(sim_hit_mom > m_mom_threshold_high) continue; 
      if(SimHit.getEDep() <= m_edep_threshold) continue;

      //Wire efficiency
      float hitProb = fRandom.Uniform(1.);
      if(hitProb > m_wireEff) continue;

      if ( id_hits_map.find(id) != id_hits_map.end()) id_hits_map[id].push_back(SimHit);
      else 
      {
          std::vector< decltype(SimHit) > vhit;
          vhit.push_back(SimHit);
          id_hits_map[id] = vhit ;
      }
  }
  if(m_WriteAna && (nullptr!=m_tuple)){
      m_n_sim = 0;
      m_n_digi = 0 ;
  }

  std::vector<TVector2> HitPosition2D;
  bool ismixNoise[55] = {0};
  std::cout << " id_hits_map  size = " << id_hits_map.size() << std::endl;

  for(auto iter = id_hits_map.begin(); iter != id_hits_map.end(); iter++)
  {
      unsigned long long wcellid = iter->first;
      auto trkHit = Vec->create();
      auto TruthtrkHit = TruthVec->create();
      trkHit.setCellID(wcellid);
      TruthtrkHit.setCellID(wcellid);
      float tot_edep   = 0 ;
      float tot_length = 0 ;
      float pos_x = 0 ;
      float pos_y = 0 ;
      float pos_z = 0 ;
      float momx,momy,momz = 0;
      const int simhit_size = iter->second.size();
      for(unsigned int i=0; i< simhit_size; i++)
      {
          tot_edep += iter->second.at(i).getEDep();//GeV
      }
      int chamber = m_decoder->get(wcellid, "chamber");
      int layer   = m_decoder->get(wcellid, "layer"  );
      int cellID  = m_decoder->get(wcellid, "cellID" );
      TVector3 Wstart(0,0,0);
      TVector3 Wend  (0,0,0);
      m_segmentation->cellposition(wcellid, Wstart, Wend);
      float dd4hep_mm = dd4hep::mm;
      if(m_debug) std::cout<<"DCHDigi wcellid ="<<wcellid<< ",chamber="<<chamber<<",layer="<<layer<<",cellID="<<cellID<<",s_x="<<Wstart.x()<<",s_y="<<Wstart.y()<<",s_z="<<Wstart.z()<<",E_x="<<Wend.x()<<",E_y="<<Wend.y()<<",E_z="<<Wend.z()<<std::endl;

      TVector3  denominator = (Wend-Wstart) ;
      float min_distance = 999 ;
      float sim_distance = 999 ;
      float tmp_distance =0;
      TVector3 hitPosition;
      TVector3 PCA;
      for(unsigned int i=0; i< simhit_size; i++)
      {
          float sim_hit_mom = sqrt( iter->second.at(i).getMomentum()[0]*iter->second.at(i).getMomentum()[0] + iter->second.at(i).getMomentum()[1]*iter->second.at(i).getMomentum()[1] + iter->second.at(i).getMomentum()[2]*iter->second.at(i).getMomentum()[2] );//GeV
          if(sim_hit_mom<1e-4) continue;

          TVector3  pos(iter->second.at(i).getPosition()[0]*dd4hep_mm, iter->second.at(i).getPosition()[1]*dd4hep_mm, iter->second.at(i).getPosition()[2]*dd4hep_mm);

          float simpos[3] = {iter->second.at(i).getPosition()[0],iter->second.at(i).getPosition()[1],iter->second.at(i).getPosition()[2]};
          float simmom[3] = {iter->second.at(i).getMomentum()[0],iter->second.at(i).getMomentum()[1],iter->second.at(i).getMomentum()[2]};

          TVector3 sim_mon(iter->second.at(i).getMomentum()[0],iter->second.at(i).getMomentum()[1],iter->second.at(i).getMomentum()[2]);
          float Steplength = iter->second.at(i).getPathLength();
          TVector3  pos_start = pos - 0.5 * Steplength * sim_mon.Unit();
          TVector3  pos_end = pos + 0.5 * Steplength * sim_mon.Unit();
          tmp_distance = m_segmentation->Distance(wcellid,pos_start,pos_end,hitPosition,PCA);
          tmp_distance = tmp_distance/dd4hep_mm; //mm
          sim_distance = tmp_distance;

          if(tmp_distance < min_distance){
              min_distance = tmp_distance;
              pos_x = hitPosition.x();     //pos.x();
              pos_y = hitPosition.y();     //pos.y();
              pos_z = pos.z();
              momx = iter->second.at(i).getMomentum()[0];
              momy = iter->second.at(i).getMomentum()[1];
          }
          tot_length += iter->second.at(i).getPathLength();//mm
          auto asso = AssoVec->create();
          asso.setRec(trkHit);
          asso.setSim(iter->second.at(i));
          asso.setWeight(iter->second.at(i).getEDep()/tot_edep);

          if(m_WriteAna && (nullptr!=m_tuple)) {
              m_simhit_x[m_n_sim] = pos.x();
              m_simhit_y[m_n_sim] = pos.y();
              m_simhit_z[m_n_sim] = pos.z();
              m_Simdca[m_n_sim] = sim_distance;
              m_simhitT[m_n_sim] = iter->second.at(i).getTime();
              m_simhitmom[m_n_sim] = sim_hit_mom;
              m_simPDG[m_n_sim] = iter->second.at(i).getMCParticle().getPDG();
              m_n_sim ++ ;
          }
      }

      //Smear drift distance
      if(m_smear){
          float dd = min_distance + gRandom->Gaus(0,fabs(m_res_x));
          if(dd>=0) min_distance = dd;
          else min_distance = -dd;
      }

      trkHit.setTime(min_distance*1e3/m_velocity);//m_velocity is um/ns, drift time in ns
      trkHit.setQuality(1);//1:singal   0:noise of overlap singal  -1:noise
      trkHit.setEDep(tot_edep);// GeV
      trkHit.setPosition (edm4hep::Vector3d(pos_x, pos_y, pos_z));//position of closest sim hit
      trkHit.setCovMatrix(std::array<float, 6>{m_res_x, 0, m_res_y, 0, 0, m_res_z});//cov(x,x) , cov(y,x) , cov(y,y) , cov(z,x) , cov(z,y) , cov(z,z) in mm

      //add noise every layer
      if(m_mixNoise && !ismixNoise[layer]) mixNoise(layer,cellID,Vec,&trkHit,ismixNoise);

      TruthtrkHit.setEDep(tot_edep);// GeV
      TruthtrkHit.setPosition (edm4hep::Vector3d(pos_x, pos_y, pos_z));//position of closest sim hit
      TruthtrkHit.setCovMatrix(std::array<float, 6>{m_res_x, 0, m_res_y, 0, 0, m_res_z});//cov(x,x) , cov(y,x) , cov(y,y) , cov(z,x) , cov(z,y) , cov(z,z) in mm

      if(m_WriteAna && (nullptr!=m_tuple)) {
          m_chamber  [m_n_digi] = chamber;
          m_layer    [m_n_digi] = layer  ;
          m_cell     [m_n_digi] = cellID;
          m_cell_x   [m_n_digi] = Wstart.x();
          m_cell_y   [m_n_digi] = Wstart.y();
          m_cell1_x   [m_n_digi] = Wend.x();
          m_cell1_y   [m_n_digi] = Wend.y();
          m_hit_x    [m_n_digi] = pos_x;
          m_hit_y    [m_n_digi] = pos_y;
          m_hit_z    [m_n_digi] = pos_z;
          m_mom_x    [m_n_digi] = momx ;
          m_mom_y    [m_n_digi] = momy ;
          m_dca      [m_n_digi] = min_distance;
          m_poca_x   [m_n_digi] = PCA.x();
          m_poca_y   [m_n_digi] = PCA.y();
          m_hit_dE   [m_n_digi] = trkHit.getEDep();
          m_truthlength[m_n_digi] = tot_length ;
          m_n_digi ++ ;
      }

  }
  debug()<<"output digi DCHhit size="<< Vec->size() <<endmsg;
  _nEvt ++ ;

  if(m_WriteAna && (nullptr!=m_tuple)){
      StatusCode status = m_tuple->write();
      if ( status.isFailure() ) {
          error() << "    Cannot fill N-tuple:" << long( m_tuple ) << endmsg;
          return StatusCode::FAILURE;
      }
  }
  m_end = clock();
  if(m_WriteAna){
      m_time = (m_end - m_start);
  }

  return StatusCode::SUCCESS;
}

StatusCode DCHDigiAlg::finalize()
{
    info() << "Processed " << _nEvt << " events " << endmsg;
    return GaudiAlgorithm::finalize();
}

void DCHDigiAlg::mixNoise(int layerID ,int wireID,
        edm4hep::TrackerHitCollection* Vec,
        edm4hep::MutableTrackerHit* trackerHitLayer,
        bool ismixNoise[55]){

    int maxCellID = m_segmentation->maxWireID(0,layerID);
    int nNoiseHits = (int) maxCellID*m_noiseEff;

    //loop over noise hits in this layer
    for(int ic=0;ic<nNoiseHits;ic++){
        int cellid = (int) maxCellID*fRandom.Uniform(1.);
        float time = m_timeWindow*fRandom.Uniform(1.);
        TVector3 Wstart,Wend;
        m_segmentation->cellposition2(0,layerID,cellid,Wstart,Wend);
        //found noise overlap with signal
        if(!(cellid != wireID))
        {
            if(time<trackerHitLayer->getTime()){
                //1. overlap with signal and time > signal time
                if(!m_mixNoiseAndSkipOverlap){
                    //set time
                    trackerHitLayer->setTime(time);
                    trackerHitLayer->setQuality(0);

                    m_noiseCover_layer[m_n_noiseCover] = layerID;
                    m_noiseCover_cell[m_n_noiseCover] = cellid;

                    m_n_noiseCover++;

                    m_noise_layer[m_n_noise] = layerID;
                    m_noise_cell[m_n_noise] = cellid;
                    m_noise_time[m_n_noise] = time;

                    m_n_noise++;
                } // end fo m_mixNoiseAndSkipOverlap
            }else{
                //2. overlap with signal and time > signal time
                // TODO adc sum with signal
                continue;
            }
        }else{ // not found noise overlap with signal
            //3. create a noise hit
            //create another trackerHit
            if(!m_mixNoiseAndSkipNoOverlap){
                auto trkNoiseHit = Vec->create();

                dd4hep::rec::CellID ncellid;
                m_decoder->set(ncellid,"chamber",0);
                m_decoder->set(ncellid,"layer",layerID);
                m_decoder->set(ncellid,"cellID",cellid);

                trkNoiseHit.setCellID(ncellid);
                trkNoiseHit.setTime(time);
                trkNoiseHit.setQuality(-1);

                m_noise_layer[m_n_noise] = layerID;
                m_noise_cell[m_n_noise] = cellid;
                m_noise_time[m_n_noise] = time;

                m_n_noise++;

            }// end of m_mixNoiseAndSkipNoOverlap
        } // end of if


    }// loop over noise hit
    ismixNoise[layerID] = true;
}//end of mixNoise
