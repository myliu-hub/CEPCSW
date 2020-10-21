#include "GearSvc/IGearSvc.h"
#include "gear/GearMgr.h"

#include "MarlinKalTest.h"

#include "TrackSystemSvc.h"

DECLARE_COMPONENT(TrackSystemSvc)

TrackSystemSvc::TrackSystemSvc(const std::string& name, ISvcLocator* svc)
  : base_class(name, svc){
}

TrackSystemSvc::~TrackSystemSvc(){
}

MarlinTrk::IMarlinTrkSystem* TrackSystemSvc::getTrackSystem(void* address){
  std::map<void*, MarlinTrk::IMarlinTrkSystem*>::iterator it=m_trackSystems.find(address);
  if(it==m_trackSystems.end()){
    gear::GearMgr* mgr=0; 
    auto _gear = service<IGearSvc>("GearSvc");
    if ( !_gear ) {
      info() << "Failed to find GearSvc ..." << endmsg;
    }
    else{
      mgr = _gear->getGearMgr();
    }

    auto _geoSvc = service<IGeoSvc>("GeoSvc");
    if ( !_geoSvc ) {
      info() << "Failed to find GeoSvc ..." << endmsg;
    }
    if(mgr==0&&_geoSvc==0){
      fatal() << "Both GearSvc and GeoSvc invalid!" << endmsg;
      return 0;
    }
    debug() << "GearMgr=" << mgr << " GeoSvc=" << _geoSvc << endmsg;
    MarlinTrk::IMarlinTrkSystem* sys = new MarlinTrk::MarlinKalTest( *mgr, _geoSvc );
    m_trackSystems[address] = sys;
    debug() << "Track system created successfully for " << address << endmsg;
    return sys;
  }
  return it->second;
}

StatusCode TrackSystemSvc::initialize(){
  for(std::map<void*, MarlinTrk::IMarlinTrkSystem*>::iterator it=m_trackSystems.begin();it!=m_trackSystems.end();it++){
    delete it->second;
  }
  m_trackSystems.clear();

  m_trackSystems[0] = getTrackSystem(0);
  
  return StatusCode::SUCCESS;
}

void TrackSystemSvc::removeTrackSystem(void* address){
  std::map<void*, MarlinTrk::IMarlinTrkSystem*>::iterator it=m_trackSystems.find(address);
  if ( it!=m_trackSystems.end() ) {
    delete it->second;
    m_trackSystems.erase(it);
  }
  return;
}

StatusCode TrackSystemSvc::finalize(){
  for(std::map<void*, MarlinTrk::IMarlinTrkSystem*>::iterator it=m_trackSystems.begin();it!=m_trackSystems.end();it++){
    delete it->second;
  }
  
  return StatusCode::SUCCESS;
}
