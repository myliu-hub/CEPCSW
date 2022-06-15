#include "SteppingAction.h"

SteppingAction::SteppingAction(ToolHandleArray<IAnaElemTool>& anatools)
    : m_anaelemtools(anatools) {

}

SteppingAction::~SteppingAction() {

}

void
SteppingAction::UserSteppingAction(const G4Step* aStep) {

//    G4Track *currentTrack = aStep->GetTrack();
//    if(currentTrack->GetDefinition()->GetParticleName() != "e-" &&
//       currentTrack->GetDefinition()->GetParticleName() != "e+" &&
//       currentTrack->GetDefinition()->GetParticleName() != "mu-" &&
//       currentTrack->GetDefinition()->GetParticleName() != "mu+" &&
//       currentTrack->GetDefinition()->GetParticleName() != "pi-" &&
//       currentTrack->GetDefinition()->GetParticleName() != "pi+" &&
//       currentTrack->GetDefinition()->GetParticleName() != "pi0" &&
//       currentTrack->GetDefinition()->GetParticleName() != "K-" &&
//       currentTrack->GetDefinition()->GetParticleName() != "K+" &&
//       currentTrack->GetDefinition()->GetParticleName() != "K0" &&
//       currentTrack->GetDefinition()->GetParticleName() != "anti-p-" &&
//       currentTrack->GetDefinition()->GetParticleName() != "p+"){
//
//            currentTrack->SetTrackStatus(fKillTrackAndSecondaries);
//    }

    for (auto ana: m_anaelemtools) {
        ana->UserSteppingAction(aStep);
    }
}

