//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm10/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "Randomize.hh"

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4VProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction), fVerboseLevel(0),
    fStepNumber(0),
    fParticleType(0),
    fEnergy(0.),
    fSecondaryParticleType(0),
    fSecondaryEnergy(0.),
    fMomentumX(0.),
    fMomentumY(0.),
    fMomentumZ(0.),
    fprocessType(0.)
{
  G4RunManager::GetRunManager()->SetPrintProgress(10000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event*)
{
    fStepNumber = 0;
    fParticleType = 0;
    fEnergy = 0.;
    fSecondaryParticleType = 0;
    fSecondaryEnergy = 0.;
    fMomentumX = 0.;
    fMomentumY = 0.;
    fMomentumZ = 0.;
    fprocessType = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* event)
{  
  // save rndm status
  if (fRunAction->GetRndmFreq() == 2) { 
    CLHEP::HepRandom::saveEngineStatus("endOfEvent.rndm");        
 
    // show rndm status
    G4int eventNb = event->GetEventID();
    G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if (eventNb%printModulo == 0) { 
      G4cout << "\n---> End of Event: " << eventNb << G4endl;
      CLHEP::HepRandom::showEngineStatus();
    }
  }
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleIColumn(0, fStepNumber);
  analysisManager->FillNtupleIColumn(1, fParticleType);
  analysisManager->FillNtupleDColumn(2, fEnergy);
  analysisManager->FillNtupleIColumn(3, fSecondaryParticleType);
  analysisManager->FillNtupleDColumn(4, fSecondaryEnergy);
  analysisManager->FillNtupleDColumn(5, fMomentumX);
  analysisManager->FillNtupleDColumn(6, fMomentumY);
  analysisManager->FillNtupleDColumn(7, fMomentumZ);
  analysisManager->FillNtupleDColumn(8, fprocessType);
  analysisManager->AddNtupleRow();
}
void EventAction::RecordStepData(G4int stepNumber, G4int particleType, G4double energy,
    G4int secondaryParticleType, G4double secondaryEnergy,
    G4double momentumX, G4double momentumY, G4double momentumZ, G4double processType)
{
    fStepNumber = stepNumber;
    fParticleType = particleType;
    fEnergy = energy;
    fSecondaryParticleType = secondaryParticleType;
    fSecondaryEnergy = secondaryEnergy;
    fMomentumX = momentumX;
    fMomentumY = momentumY;
    fMomentumZ = momentumZ;
    fprocessType = processType;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
