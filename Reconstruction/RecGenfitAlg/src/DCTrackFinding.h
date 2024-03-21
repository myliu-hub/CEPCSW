//////////////////////////////////////////////////////////////////////
///
/// This is an algorithm for track fitting for CEPC track with genfit.
///
/// In this file, including:
///   An algorithm for combined silicon and drift chamber track fitting
///   with genfit for 5 particle hypothesis
///
///   Units are following DD4hepUnits
///
/// Authors:
///   Mengyao Liu(myliu@ihep.ac.cn)
///
/////////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_DCTRACKFINDING_H
#define RECGENFITALG_DCTRACKFINDING_H

#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/NTuple.h"
#include "k4FWCore/DataHandle.h"
#include "DD4hep/Fields.h"
#include <string>
#include <vector>
#include "AbsKalmanFitter.h"

#include "DataHelper/GeomeryHelper.h"

#include "GenfitTrack.h"
#include "GenfitFitter.h"
#include "GenfitField.h"
#include "GenfitUnit.h"


//ROOT
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TLorentzVector.h"

class GenfitTrack;
class GenfitFitter;
class GenfitField;
class IGeomSvc;
//class GeomeryWire;

namespace dd4hep {
    class Detector;
    namespace DDSegmentation{
        class GridDriftChamber;
        class BitFieldCoder;
    }
}
namespace edm4hep{
    class EventHeaderCollection;
    class MCParticleCollection;
    class SimTrackerHitCollection;
    class TrackCollection;
    class TrackerHitCollection;
    class MCRecoTrackerAssociationCollection;
    class ReconstructedParticle;
    class ReconstructedParticleCollection;
    class TrackerHit;
    class Track;
}

/////////////////////////////////////////////////////////////////////////////

class DCTrackFinding:public GaudiAlgorithm {
    public:
        DCTrackFinding (const std::string& name,ISvcLocator* pSvcLocator);
        StatusCode initialize() override; StatusCode execute() override; StatusCode finalize() override;

    private:
        //genfit::AbsKalmanFitter* m_genfitFitter;
        GeomeryWire* geomery_wire;
        GenfitFitter* m_genfitFitter;
        const GenfitField* m_genfitField;
        SmartIF<IGeomSvc> m_geomSvc;
        dd4hep::OverlayedField m_dd4hepField;
        dd4hep::Detector* m_dd4hepDetector;
        
        dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

        Gaudi::Property<std::string> m_pathSelectFilterName{this,"pathSelectFilterName","distance"};
        Gaudi::Property<int> m_maximalLayerJump{this,"MaximalLayerJump",1};

        Gaudi::Property<double> m_extMinDistCut{this,"extMinDistCut",1e-4};
        Gaudi::Property<bool> m_skipWireMaterial{this,
            "skipWireMaterial",false};
        Gaudi::Property<bool> m_correctBremsstrahlung{this,
            "correctBremsstrahlung",true};
        Gaudi::Property<std::string> m_fitterType{this,"fitterType","DAF"};
        Gaudi::Property<bool> m_noMaterialEffects{this,
            "noMaterialEffects",false};
        Gaudi::Property<int> m_debugPid{this,"debugPid",-99};
        Gaudi::Property<double> m_bStart{this,"bStart",100};
        Gaudi::Property<double> m_bFinal{this,"bFinal",0.01};
        Gaudi::Property<int> m_maxIteration{this,"maxIteration",20};
        Gaudi::Property<std::string> m_genfitHistRootName{this,
            "genfitHistRootName",""};
        Gaudi::Property<double> m_driftVelocity{this,"driftVelocity",40};//um/ns
        Gaudi::Property<bool> m_Smear{this,"Smear",false};
        Gaudi::Property<float> m_resX{this,"resX",0.11};//mm
        Gaudi::Property<double> m_Noiseff{this,"Noiseff",0.5};//um/ns
        Gaudi::Property<double> m_sigma{this,"sigmaL",0.011};//0.011cm=110um
        Gaudi::Property<std::string> m_readout_name{this,
            "readout", "DriftChamberHitsCollection"};

        // Rec Si Track
        DataHandle<edm4hep::TrackCollection> m_siSubsetTrackCol{
            "SubsetTracks" , Gaudi::DataHandle::Reader, this};
        //DC DigiHit SimTrackHit DCHitAssociation
        DataHandle<edm4hep::TrackerHitCollection> m_DCDigiCol{
            "DigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::TrackerHitCollection> m_TruthDCDigiCol{
            "TruthDigiDCHitCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::SimTrackerHitCollection> m_simDCHitCol{
            "DriftChamberHitsCollection" , Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_DCHitAssociationCol{"DCHitAssociationCollection",
                Gaudi::DataHandle::Reader, this};
        //SIT SimTrackHit DCHitAssociation 
        DataHandle<edm4hep::SimTrackerHitCollection> m_simSITHitCol{
            "SITCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_SITHitAssociationCol{"SITTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        //VXD SimTrackHit DCHitAssociation 
        DataHandle<edm4hep::SimTrackerHitCollection> m_simVXDHitCol{
            "VXDCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_VXDHitAssociationCol{"VXDTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        //FTD SimTrackHit DCHitAssociation 
        DataHandle<edm4hep::SimTrackerHitCollection> m_simFTDHitCol{
            "FTDCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_FTDHitAssociationCol{"FTDTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};
        //SET SimTrackHit DCHitAssociation 
        DataHandle<edm4hep::SimTrackerHitCollection> m_simSETHitCol{
            "SETCollection", Gaudi::DataHandle::Reader, this};
        DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
            m_SETHitAssociationCol{"SETTrackerHitAssociation",
                Gaudi::DataHandle::Reader, this};

        DataHandle<edm4hep::TrackCollection> m_SDTTrackCol{"SDTTrackCollection",
            Gaudi::DataHandle::Reader, this};

        DataHandle<edm4hep::TrackCollection> m_SDTTrackFindCol{"SDTTrackFindCollection",
            Gaudi::DataHandle::Writer, this};

        // Output collections
        DataHandle<edm4hep::TrackerHitCollection>    w_DCTrackFindingCol{
            "DCTrackFindingHitCollection", Gaudi::DataHandle::Writer, this};

        int getNumPointsWithFittedInfo(genfit::Track genfitTrack,int repID) const;
        int getFittedState(genfit::Track genfitTrack, 
                TLorentzVector& pos, TVector3& mom, TMatrixDSym& cov,
                int trackPointId=0, int repID=0, bool biased=true) const;
        genfit::AbsTrackRep* getRep(int id,genfit::Track* track) const;
        bool getMOP(int hitID, genfit::MeasuredStateOnPlane& mop,
                genfit::AbsTrackRep* trackRep,genfit::Track* track) const;

        int addHitsToTk(edm4hep::TrackerHitCollection* col, edm4hep::Track& track, const char* msg) const;
        int m_eventNo;

        /// tuples
        NTuple::Tuple*  m_tuple;
        NTuple::Item<int> m_n_N_SiTrack;
        NTuple::Array<int> m_n_DCTrackFinding;
        NTuple::Array<float> m_CPUTime;
        //NTuple::Array<int> m_n_reslayer;
        NTuple::Item<int> m_n_Digi;
        NTuple::Item<int> m_n_TruthDigi;
        NTuple::Array<int> m_n_match;
        NTuple::Matrix<int> m_layer;
        NTuple::Matrix<int> m_cellID;
        NTuple::Matrix<int> m_NoiseFlag;
        NTuple::Array<int> m_isNoOverlap;
        NTuple::Array<int> m_digiOverlap;
        NTuple::Matrix<double> m_hit_wx;
        NTuple::Matrix<double> m_hit_wy;
        NTuple::Matrix<double> m_hit_wz;
        NTuple::Matrix<double> m_hit_ex;
        NTuple::Matrix<double> m_hit_ey;
        NTuple::Matrix<double> m_hit_ez;
        NTuple::Matrix<double> m_recoPosX;
        NTuple::Matrix<double> m_recoPosY;
        NTuple::Matrix<double> m_recoPosZ;
        NTuple::Matrix<double> m_cellPosX;
        NTuple::Matrix<double> m_cellPosY;
        NTuple::Matrix<double> m_cellPosZ;
        NTuple::Matrix<double> m_match_hit_wx;
        NTuple::Matrix<double> m_match_hit_wy;
        NTuple::Matrix<double> m_match_hit_wz;
        NTuple::Matrix<double> m_match_hit_ex;
        NTuple::Matrix<double> m_match_hit_ey;
        NTuple::Matrix<double> m_match_hit_ez;
        //NTuple::Matrix<int> m_reslayer;
        //NTuple::Matrix<double> m_resPhi;
        NTuple::Matrix<double> m_radius;
        NTuple::Matrix<double> m_arcLength;
        NTuple::Matrix<double> m_hitDistance;
        NTuple::Matrix<double> m_doca;
        NTuple::Matrix<double> m_ReconstructedZ;
        NTuple::Matrix<double> m_trackParD0;
        NTuple::Matrix<double> m_trackParPhi0;
        NTuple::Matrix<double> m_trackParOmega;
        NTuple::Matrix<double> m_trackParZ0;
        NTuple::Matrix<double> m_trackParTanLambda;
        NTuple::Matrix<double> m_centerX;
        NTuple::Matrix<double> m_centerY;

};
#endif
