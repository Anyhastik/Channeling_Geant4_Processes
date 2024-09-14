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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TransitionRadiationPhysics.hh"
#include "DetectorConstruction.hh"
#include "G4ChanRad.hh"
#include "G4ProcessManager.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4ChanRad.hh"
#include "G4RegularXTRadiator.hh"
#include "G4TransparentRegXTRadiator.hh"
#include "G4GammaXTRadiator.hh"
#include "G4StrawTubeXTRadiator.hh"

#include "G4XTRGammaRadModel.hh"
#include "G4XTRRegularRadModel.hh"
#include "CrystalTarget.hh"

G4ThreadLocal 
G4ChanRad* TransitionRadiationPhysics::fXTRProcess = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TransitionRadiationPhysics::TransitionRadiationPhysics(G4int verb,
    DetectorConstruction* ptr) 
  : G4VPhysicsConstructor("XTR"), 
    fDetector(ptr), 
    fVerbose(verb),
    fXTRModel("transpM")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TransitionRadiationPhysics::~TransitionRadiationPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TransitionRadiationPhysics::ConstructProcess()
{
  if("dummy" == fXTRModel) { return; }
  if(0 < fVerbose) {
    G4cout<< "TransitionRadiationPhysics: XTR model <" << fXTRModel
          << ">" <<G4endl;
  }
  RadiatorDescription* rDescription = fDetector->GetRadiatorDescription();

  if(fXTRModel == "transpM" ) 
  { 
    fXTRProcess = new CrystalTarget(rDescription->fLogicalVolume,
                                                rDescription->fRadMaterial,
                                                rDescription->fRadThickness,
												1565.55896, // менять через мессенджер
                                                "RegularXTRadiator");//ChanRad
  }     
  if(!fXTRProcess) { 
    if(0 < fVerbose) {
      G4cout<< "TransitionRadiationPhysics: XTR model <" << fXTRModel
            << "> is not known - no XTR process defined" <<G4endl;
    }
    return; 
  }

  fXTRProcess->SetVerboseLevel(fVerbose);

  G4Electron* elec = G4Electron::Electron();
  G4ProcessManager* manager = elec->GetProcessManager();
  manager->AddDiscreteProcess(fXTRProcess);

  G4Positron* posi = G4Positron::Positron();
  manager = posi->GetProcessManager();
  manager->AddDiscreteProcess(fXTRProcess);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

