#ifndef DCH_DIGI_ALG_H
#define DCH_DIGI_ALG_H

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeomSvc.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"

#include "TVector3.h"



class DCHDigiAlg : public GaudiAlgorithm
{
 
public:
 
  DCHDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual StatusCode initialize() ;
 
  /** Called for every event - the working horse.
   */
  virtual StatusCode execute() ; 
 
  /** Called after data processing for clean up.
   */
  virtual StatusCode finalize() ;
 
protected:

  SmartIF<IGeomSvc> m_geosvc;
  typedef std::vector<float> FloatVec;
  int _nEvt ;

  NTuple::Tuple* m_tuple = nullptr ;
  NTuple::Item<long>   m_n_sim;
  NTuple::Item<long>   m_n_digi;
  NTuple::Array<int  > m_chamber   ;
  NTuple::Array<int  > m_layer     ;
  NTuple::Array<int  > m_cell      ;
  NTuple::Array<float> m_cell_x    ;
  NTuple::Array<float> m_cell_y    ;
  NTuple::Array<float> m_simhit_x  ;
  NTuple::Array<float> m_simhit_y  ;
  NTuple::Array<float> m_simhit_z  ;
  NTuple::Array<float> m_hit_x     ;
  NTuple::Array<float> m_hit_y     ;
  NTuple::Array<float> m_hit_z     ;
  NTuple::Array<float> m_dca       ;
  NTuple::Array<float> m_hit_dE    ;
  NTuple::Array<float> m_hit_dE_dx ;



  dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
  dd4hep::DDSegmentation::GridDriftChamber* m_segmentation;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  
  Gaudi::Property<std::string> m_readout_name{ this, "readout", "DriftChamberHitsCollection"};//readout for getting segmentation
 
  Gaudi::Property<float> m_res_x     { this, "res_x", 0.11};//mm
  Gaudi::Property<float> m_res_y     { this, "res_y", 0.11};//mm
  Gaudi::Property<float> m_res_z     { this, "res_z", 1   };//mm
  Gaudi::Property<float> m_velocity  { this, "drift_velocity", 40};// um/ns
  Gaudi::Property<float> m_mom_threshold { this, "mom_threshold", 0};// GeV
  Gaudi::Property<bool>  m_WriteAna { this, "WriteAna", false};

//  int cellsum[101] = {0,505, 1017, 1535, 2059, 2589, 3126, 3669, 4218, 4774, 5336, 5904, 6478, 7059, 7646, 8239, 8839, 9445, 10057, 10675, 11300, 11931, 12568,13212, 13862, 14518, 15180, 15849, 16524, 17205, 17893,18587, 19287, 19993, 20706, 21425, 22150, 22881, 23619, 24363, 25113, 25870, 26633, 27402, 28177, 28959,29747, 30541, 31342, 32149, 32962, 33781, 34607, 35439, 36277, 37122, 37973, 38830, 39693, 40563, 41439, 42321, 43210, 44105, 45006, 45913, 46827, 47747, 48673,49606, 50545, 51490, 52441, 53399, 54363, 55333, 56310, 57293, 58282, 59277, 60279, 61287, 62301, 63322, 64349, 65382, 66421, 67467, 68519, 69577, 70641, 71712, 72789, 73872, 74962, 76058, 77160, 78268, 79383, 80504, 81631};

int cellnum[100] = {505,512,518,524,530,537,543,549,556,562,568,574,581,593,600,606,612,618,625,631,637,644,650,656,662,669,675,681,688,694,700,706,713,719,725,731,738,744,750,757,763,769,775,782,788,794,801,807,813,819,826,832,838,845,851,857,863,870,876,882,889,895,901,907,914,920,926,933,939,945,951,958,964,970,977,983,989,995,1002,1008,1014,1021,1027,1033,1039,1046,1052,1058,1064,1071,1077,1083,1090,1096,1102,1108,1115,1121,1127};

  // Input collections
  DataHandle<edm4hep::SimTrackerHitCollection> r_SimDCHCol{"DriftChamberHitsCollection", Gaudi::DataHandle::Reader, this};
  // Output collections
  DataHandle<edm4hep::TrackerHitCollection>    w_DigiDCHCol{"DigiDCHitCollection", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoTrackerAssociationCollection>    w_AssociationCol{"DCHitAssociationCollection", Gaudi::DataHandle::Writer, this};
};

#endif
