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
    class ConstTrackerHit;
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

        Gaudi::Property<double> m_extMinDistCut{this,"extMinDistCut",1e-4};
        Gaudi::Property<bool> m_skipWireMaterial{this,
            "skipWireMaterial",true};
        Gaudi::Property<bool> m_correctBremsstrahlung{this,
            "correctBremsstrahlung",false};
        Gaudi::Property<std::string> m_fitterType{this,"fitterType","DAF"};
        Gaudi::Property<bool> m_noMaterialEffects{this,
            "noMaterialEffects",false};
        Gaudi::Property<int> m_debugPid{this,"debugPid",-99};
        Gaudi::Property<double> m_bStart{this,"bStart",100};
        Gaudi::Property<double> m_bFinal{this,"bFinal",0.01};
        Gaudi::Property<int> m_maxIteration{this,"maxIteration",20};
        Gaudi::Property<std::string> m_genfitHistRootName{this,
            "genfitHistRootName",""};
        Gaudi::Property<double> m_driftVelocity{this,"drift_velocity",40};//um/ns
        Gaudi::Property<double> m_sigma{this,"sigmaL",0.00011};//0.00011m=110um
        Gaudi::Property<std::string> m_readout_name{this,
            "readout", "DriftChamberHitsCollection"};

        // Rec Si Track
        DataHandle<edm4hep::TrackCollection> m_siSubsetTrackCol{
            "SubsetTracks" , Gaudi::DataHandle::Reader, this};
        //DC DigiHit SimTrackHit DCHitAssociation
        DataHandle<edm4hep::TrackerHitCollection> m_DCDigiCol{
            "DigiDCHitCollection", Gaudi::DataHandle::Reader, this};
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
        //DataHandle<edm4hep::MCRecoTrackerAssociationCollection>
        //    m_VXDHitAssociationCol{"VXDTrackerHitAssociation",
        //        Gaudi::DataHandle::Reader, this};
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
};
#endif
