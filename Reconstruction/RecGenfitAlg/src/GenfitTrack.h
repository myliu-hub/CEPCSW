//////////////////////////////////////////////////////////////////
///
/// This is an interface of to call genfit
/// A genfit track can be created, fitted and extrapolated
/// Track with only track representation(s)
///
/// In this file, including:
///   a genfit track class
///
///   Units are following DD4hepUnits
///
/// Authors:
///   Zhang Yao (zhangyao@ihep.ac.cn)
///   Y.Fujii (yfujii@ihep.ac.cn)
///   Yohei Nakatsugawa (yohei@ihep.ac.cn)
///
//////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_GENFITTRACK_H
#define RECGENFITALG_GENFITTRACK_H

#include "GenfitFitter.h"

//ROOT
#include "TVectorD.h"
#include "TMatrixDSym.h"

//Gaudi
#include "GaudiKernel/SmartIF.h"

//STL
#include <vector>

class TLorentzVector;
class IGeomSvc;

namespace genfit{
    class Track;
    class FitStatus;
    class AbsTrackRep;
    class RKTrackRep;
    class KalmanFittedStateOnPlane;
}
namespace edm4hep{
    class MCParticle;
    class SimTrackerHitCollection;
    class ReconstructedParticle;
    class MCRecoTrackerAssociationCollection;
    class Track;
    class TrackCollection;
    class ConstTrackerHit;
    class ConstSimTrackerHit;
    class Vector3d;
    class Vector3f;
}
namespace dd4hep {
    namespace DDSegmentation{
        class GridDriftChamber;
    }
    namespace rec{
        class ISurface;
    }
}

class GenfitTrack {
    friend int GenfitFitter::processTrack(
            GenfitTrack* track, bool resort=true);

    friend int GenfitFitter::processTrackWithRep(
            GenfitTrack* track, int repID, bool resort=true);

    public:
    GenfitTrack(const GenfitField* field,
            const dd4hep::DDSegmentation::GridDriftChamber* seg,
            SmartIF<IGeomSvc> geom,
            const char* name="GenfitTrack");
    virtual ~GenfitTrack();

    /// ---------Add a Genfit track-------
    virtual bool createGenfitTrack(int pdgType,int charge,TLorentzVector pos,
            TVector3 mom, TMatrixDSym covM);

    ///Create genfit track from MCParticle
    bool createGenfitTrackFromMCParticle(int pidTyep,const edm4hep::MCParticle&
            mcParticle, double eventStartTime=0.);
    bool createGenfitTrackFromEDM4HepTrack(int pidType,const edm4hep::Track& track,
            double eventStartTime,bool isUseCovTrack);

    /// ---------Add measurements---------
    ///Add one space point measurement, return number of hits on track
    virtual bool addSpacePointMeasurement(const TVectorD&,std::vector<float>
        sigmaU,std::vector<float> sigmaV,int cellID,int hitID);

    ///Add silicon space points from edm4hep::track
    int addSpacePointsSi(const edm4hep::Track& track,
            std::vector<float> sigmaU,std::vector<float> sigmaV);

    ///Add drift chamber space points from edm4hep::track
    int addSpacePointsDC(const edm4hep::Track& track,
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
        std::vector<float> sigmaU,std::vector<float> sigmaV);


    ///Add WireMeasurements of hits on track
    virtual int addWireMeasurements(edm4hep::Track& track,float sigma,
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            int sortMethod,bool truthAmbig,float skipCorner, float skipNear);

    ///Add one silicon hits
    bool addSiliconMeasurement(edm4hep::ConstTrackerHit& hit,
            float sigmaU,float sigmaV,int cellID,int hitID);

    ///Add silicon measurements, return number of hits on track
    int addSiliconMeasurements(edm4hep::Track& track,
            std::vector<float> sigmaU,std::vector<float> sigmaV);

    ///Store track to ReconstructedParticle
    bool storeTrack(edm4hep::ReconstructedParticle& recParticle,
            edm4hep::Track& track, int pidType,
            int ndfCut, double chi2Cut);

    ///A tool to convert track to the first layer of DC
    void pivotToFirstLayer(const edm4hep::Vector3d& pos,
            const edm4hep::Vector3f& mom, edm4hep::Vector3d& firstPos,
            edm4hep::Vector3f& firstMom);

    /// Copy a track to event
    //void CopyATrack()const;

    ///Extrapolate to Hit
    /// Extrapolate the track to the drift chamber hit
    /// Output: poca pos and dir and poca distance to the hit wire
    /// Input: genfit track, pos and mom, two ends of a wire
    ///        pos, and mom are position & momentum at starting point
    double extrapolateToHit( TVector3& poca, TVector3& pocaDir,
            TVector3& pocaOnWire, double& doca, TVector3 pos, TVector3 mom,
            TVector3 end0,//one end of the hit wire
            TVector3 end1,//the orhter end of the hit wire
            int repID,
            bool stopAtBoundary=false,
            bool calcJacobianNoise=true) const;
    /// Extrapolate the track to the point
    /// Output: pos and mom of POCA point to point
    /// Input: genfitTrack,point,repID,stopAtBoundary and calcAverageState
    /// repID same with pidType
    double extrapolateToPoint(TVector3& pos, TVector3& mom,
            const TVector3& point, int repID=0, bool stopAtBoundary = false,
            bool calcJacobianNoise = true) const;

    double extrapolateToPoint(TVector3& pos, TVector3& mom, TMatrixDSym& cov,
            const TVector3& point, int repID=0,
            bool stopAtBoundary = false, bool calcJacobianNoise = true) const;

    /// Extrapolate the track to the cyliner at fixed raidus
    /// Output: pos and mom at fixed radius
    /// Input: genfitTrack, radius of cylinder at center of the origin,
    ///        repID, stopAtBoundary and calcAverageState
    double extrapolateToCylinder(TVector3& pos, TVector3& mom,
            double radius, const TVector3 linePoint,
            const TVector3 lineDirection, int hitID =0, int repID=0,
            bool stopAtBoundary=false, bool calcJacobianNoise=true);


    bool getPosMomCovMOP(int hitID, TLorentzVector& pos, TVector3& mom,
            TMatrixDSym& cov, int repID=0) const;

    /// get the seed position and momentum from track
    const TLorentzVector getSeedStatePos() const;
    const TVector3 getSeedStateMom() const;
    void getSeedStateMom(TLorentzVector& pos, TVector3& mom) const;
    unsigned int getNumPoints() const;
    unsigned int getNumPointsDet(int cellID) const;

    /// get the fitted track status
    const genfit::FitStatus* getFitStatus(int repID=0) const;
    int getFittedState(TLorentzVector& pos, TVector3& mom, TMatrixDSym& cov,
            int trackPointId=0, int repID=0, bool biased=true) const;
    int getNumPointsWithFittedInfo(int repID=0) const;
    bool getFirstPointWithFittedInfo(int repID=0) const;
    bool fitSuccess(int repID) const;

    /// get the wire infomation
    int getDetIDWithFitterInfo(int hitID, int idRaw=0) const;

    int getPDG(int id=0) const;
    int getPDGCharge(int id=0) const;

    /// print genfit track
    void printSeed() const;//print seed
    void printFitted(int repID=0) const;//print seed
    void print(TLorentzVector pos, TVector3 mom, const char* comment="") const;

    /// set and get debug level
    void setDebug(int debug){m_debug=debug;}
    void setDebugGenfit(int debug);
    void setDebugLocal(int debug);
    int getDebug(void) const { return m_debug;}

    /// get name of this object
    const char* getName() const {return m_name;}
    genfit::Track* getTrack() const{return m_track;}

    /// Add a track representation
    genfit::RKTrackRep* addTrackRep(int pdg);

    protected:
    //genfit::Track* getTrack() {return m_track;}

    private:

    int getDetTypeID(unsigned long long cellID) const;
    const char* m_name;

    ///Note! private functions are using genfit unit, cm and MeV

    genfit::AbsTrackRep* getRep(int id=0) const;
    bool getMOP(int hitID, genfit::MeasuredStateOnPlane& mop,
            genfit::AbsTrackRep* trackRep=nullptr) const;
    const dd4hep::rec::ISurface* getISurface(edm4hep::ConstTrackerHit hit);
    void getAssoSimTrackerHit(
            const edm4hep::MCRecoTrackerAssociationCollection* assoHits,
            edm4hep::ConstTrackerHit trackerHit,
            edm4hep::ConstSimTrackerHit& simTrackerHit) const;

    genfit::Track* m_track;/// track
    int m_debug;/// debug level
    int m_debugLocal;/// debug level local

    const GenfitField* m_genfitField;//pointer to genfit field
    const dd4hep::DDSegmentation::GridDriftChamber* m_gridDriftChamber;

    static const int s_PDG[2][5];

    SmartIF<IGeomSvc> m_geomSvc;

};

#endif
