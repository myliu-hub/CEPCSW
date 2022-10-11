#include "MBg.h"

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

DECLARE_COMPONENT(MakeBackground)

MakeBackground::MakeBackground(const std::string& name, ISvcLocator* svcLoc)
 : GaudiAlgorithm(name, svcLoc),
   _nEvt(0)
{

    declareProperty("DriftChamberHitsCollection", r_SimDCHCol,
            "Handle of the Input SimHit collection");
    declareProperty("DigiDCHitCollection", r_DigiDCHCol,
            "Handle of DC digi(TrackerHit) collection");
    declareProperty("DCHitAssociationCollection", r_AssociationCol,
            "Handle of Association collection");

    declareProperty("NoiseSimHitsCollection",w_SimNoiseHCol,
            "Handle of DC NoiseHitsCollection collection");
    declareProperty("NoiseDCHitsCollection",w_NoiseHitCol,
            "Handle of DC NoiseHitsCollection collection");
    declareProperty("NoiseDCHitAssociationCollection",w_NoiseAssociationCol,
            "Handle of NoiseDCAssociation collection");
}

StatusCode MakeBackground::initialize()
{
    m_geoSvc = service<IGeomSvc>("GeomSvc");
    if ( !m_geoSvc )  throw "MakeBackground :Failed to find GeomSvc ...";

    dd4hep::Detector* m_dd4hep = m_geoSvc->lcdd();
    if ( !m_dd4hep )  throw "MakeBackground :Failed to get dd4hep::Detector ..."; 

    dd4hep::Readout readout = m_dd4hep->readout(m_readout_name);
    m_segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>(readout.segmentation().segmentation());

    m_decoder = m_geoSvc->getDecoder(m_readout_name);

    if (!m_decoder) {
        error() << "Failed to get the decoder. " << endmsg;
        return StatusCode::FAILURE;
    }

//    fRandom.SetSeed(105105);//FIXME: set by users

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode MakeBackground::execute()
{
    info() << "Processing " << _nEvt << " events " << endmsg;
    StatusCode sc;

    std::map<double,edm4hep::TrackerHit> resortHit;
    resortHit.clear();

    const edm4hep::SimTrackerHitCollection* SimHitCol =  r_SimDCHCol.get();
    edm4hep::SimTrackerHitCollection* SimVec   = w_SimNoiseHCol.createAndPut();

    const edm4hep::TrackerHitCollection* digiDCHitsCol=nullptr;
    edm4hep::TrackerHitCollection* Vec   = w_NoiseHitCol.createAndPut();

    const edm4hep::MCRecoTrackerAssociationCollection* assoHits
            =r_AssociationCol.get();
    edm4hep::MCRecoTrackerAssociationCollection* AssoVec
             =  w_NoiseAssociationCol.createAndPut();

    digiDCHitsCol=r_DigiDCHCol.get();
    int nNoise = 0;

    //loop all hits
    for(unsigned int i = 0; i<(digiDCHitsCol->size()); i++)
    {
        edm4hep::TrackerHit mcHit = digiDCHitsCol->at(i);

        unsigned long long wcellid = mcHit.getCellID();
        int system  = m_decoder->get(wcellid, "system");
        int chamber = m_decoder->get(wcellid, "chamber");
        int layer   = m_decoder->get(wcellid, "layer"  );
        int cellID  = m_decoder->get(wcellid, "cellID" );

        float pocaTime = mcHit.getTime();

        // store trackHit
        auto trkHit = Vec->create();

        trkHit.setCellID(wcellid);
        trkHit.setTime(pocaTime);
        trkHit.setEDep(mcHit.getEDep());
        //trkHit.setEdx(mcHit.getEdx());
        trkHit.setPosition(mcHit.getPosition());
        trkHit.setCovMatrix(mcHit.getCovMatrix());

        for(int iAsso=0;iAsso<(int) assoHits->size();iAsso++)
        {

            if(assoHits->at(iAsso).getRec().getCellID()==wcellid )
            {
                auto SimtrkHit = SimVec->create();

                //SimtrkHit = assoHits->at(iAsso).getSim();
                SimtrkHit.setCellID(assoHits->at(iAsso).getSim().getCellID());
                SimtrkHit.setEDep(assoHits->at(iAsso).getSim().getEDep());
                SimtrkHit.setTime(assoHits->at(iAsso).getSim().getTime());
                SimtrkHit.setPathLength(assoHits->at(iAsso).getSim().getPathLength());
                SimtrkHit.setQuality(assoHits->at(iAsso).getSim().getQuality());
                SimtrkHit.setPosition(assoHits->at(iAsso).getSim().getPosition());
                SimtrkHit.setMomentum(assoHits->at(iAsso).getSim().getMomentum());
                SimtrkHit.setMCParticle(assoHits->at(iAsso).getSim().getMCParticle());

                auto asso = AssoVec->create();
                asso.setRec(trkHit);
                asso.setSim(SimtrkHit);
                asso.setWeight(SimtrkHit.getEDep()/trkHit.getEDep());
            }
        }

        double prob1 = fRandom.Uniform(1.);
        if(prob1 > m_fHitPurity) continue;

        float noiseTime = m_pocaTime*fRandom.Uniform(1.);
        if(noiseTime<pocaTime) continue;
        nNoise++;

        trkHit.setTime(noiseTime);

        for(int iAsso=0;iAsso<(int) assoHits->size();iAsso++)
        {

            if(assoHits->at(iAsso).getRec().getCellID()==wcellid )
            {

                AssoVec->at(iAsso).setRec(trkHit);

            }

        }

    }

std::cout << " nNoise = " << nNoise << std::endl;

    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

StatusCode MakeBackground::finalize()
{
    info() << "Processed " << _nEvt << " events " << endmsg;
    return GaudiAlgorithm::finalize();
}


