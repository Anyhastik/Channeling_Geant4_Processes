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
/// \file electromagnetic/TestEm10/src/DetectorSimpleALICE.cc
/// \brief Implementation of the DetectorSimpleALICE class
//
//
//
//

#include "DetectorSimpleALICE.hh"
#include "SensitiveDetector.hh"
#include "Materials.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"

#include "G4Region.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include <phys_const.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorSimpleALICE::DetectorSimpleALICE()
  : fRadiatorDescription(0) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorSimpleALICE::~DetectorSimpleALICE()
{
  // delete fRadiatorDescription;
        // the description is deleted in detector construction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorSimpleALICE::Construct()
{
  // Geometry parameters
  //

  G4cout << "DetectorSimpleALICE setup" << G4endl;

  G4double worldSizeZ = 400.*cm;
  G4double worldSizeR = 20.*cm;

  // Radiator and detector parameters
  //Старый увеличенный объем
  //G4double crystalThickness = crystaltickness*0.001*mm *10;  //0.020*mm = 20mkm
  G4double crystalThickness = crystaltickness * 0.001 * mm * 10;  //0.020*mm = 20mkm

  G4double absorberThickness = 0.01 * mm;
 // G4double absorberThickness = 0.1 * mm;

  // Materials
  //

  // Change to create materials using NIST
  /*G4Material* air   = Materials::GetInstance()->GetMaterial("Air");
  G4Material* mylar = Materials::GetInstance()->GetMaterial("Mylar");
  G4Material* xe15CO2 = Materials::GetInstance()->GetMaterial("Xe15CO2");

  G4double foilDensity = mylar->GetDensity();
  G4double gasDensity  = air->GetDensity();  
  G4double totDensity  = foilDensity*foilGasRatio 
                       + gasDensity*(1.0-foilGasRatio);

  G4double fractionFoil =  foilDensity*foilGasRatio/totDensity;
  G4double fractionGas  =  gasDensity*(1.0-foilGasRatio)/totDensity;*/


 // ------------- Materials -------------
   G4double a, z, density, temperature, pressure;

 // vacuum
  G4Material* galactic = new G4Material("Galactic", z = 1., a = 1.01 * g / mole, density = 1.e-25 * g / cm3, kStateGas,
      temperature = 0.1 * kelvin, pressure = 1.e-19 * pascal);

 //Change to create materials using NIST
  G4NistManager* nist = G4NistManager::Instance();

  // World parameters
  G4double env_sizeXY = 0.15 * mm, env_sizeZ = 0.1 * mm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;


  // ------------- Volumes -------------
  // World
  G4double world_sizeXY = 1.2 * env_sizeXY * 10; //1.2 * env_sizeXY*10;
  G4double world_sizeZ = 1.2 * env_sizeZ * 10; //1.2 * env_sizeZ*10
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_WATER");//G4_AIR

  G4Box* solidWorld =
      new G4Box("World",                       //its name
          0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld,          //its solid
          world_mat,           //its material
          "World");            //its name

  G4VPhysicalVolume* physicsWorld =//NAMING
      new G4PVPlacement(0,                     //no rotation
          G4ThreeVector(),       //at (0,0,0)
          logicWorld,            //its logical volume
          "World",               //its name
          0,                     //its mother  volume
          false,                 //no boolean operation
          0,                     //copy number
          checkOverlaps);        //overlaps checking

  //------------------------------------------
  // selection of a material of the crystal 
  G4Material* crystalMaterial =  nist->FindOrBuildMaterial("G4_W");//G4_Si // galactic;//
  if (crystaltype == "ge")
	  crystalMaterial = nist->FindOrBuildMaterial("G4_Ge");
  else if (crystaltype == "w")
	  crystalMaterial = nist->FindOrBuildMaterial("G4_W");
  else if (crystaltype == "c")
	  crystalMaterial = nist->FindOrBuildMaterial("G4_C");
  
  //----------Fake-Vacuum-Box--------------------------------
  G4Box* FakeVacuum_box = new G4Box("FakeVacuum", crystalThickness / 2.,
      crystalThickness / 2., 150.0 * um);
      //crystalThickness / 2., 1500.0 * um);


  G4LogicalVolume* FakeVacuum_log
      = new G4LogicalVolume(FakeVacuum_box, galactic, "FakeVacuum", 0, 0, 0);

  G4VPhysicalVolume* FakeVacuum_phys
      = new G4PVPlacement(0, G4ThreeVector(0, 0, -150.0 * um - 0.02 * um - crystalThickness / 2),
     // = new G4PVPlacement(0, G4ThreeVector(0, 0, -1500.0 * um -0.2 * um - crystalThickness / 2),
          FakeVacuum_log, "FakeVacuum", logicWorld,
          0, false, 0);
  //------------------------------------------


  //----------Fake-Crystal--------------------------------
  G4Box* FakeCrystal_box = new G4Box("FakeCrystal", crystalThickness / 2.,
      crystalThickness / 2., 0.01*um);
      //crystalThickness / 2., 0.1 * um);


  G4LogicalVolume* FakeCrystal_log
      = new G4LogicalVolume(FakeCrystal_box, galactic, "FakeCrystal", 0, 0, 0);

  G4VPhysicalVolume* FakeCrystal_phys
      = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.01 * um - crystalThickness / 2), 
          FakeCrystal_log, "FakeCrystal", logicWorld,
          0, false, 0);
  //------------------------------------------
  
  // // Crystal 
  // Radiator description = crystal
  fRadiatorDescription = new RadiatorDescription;
  fRadiatorDescription->fRadMaterial = crystalMaterial;
  fRadiatorDescription->fRadThickness = crystalThickness;
 
  G4Box* Crystal_box = new G4Box("Crystal", crystalThickness / 2., 
      crystalThickness / 2., crystalThickness / 2.);

  G4LogicalVolume* Crystal_log
   = new G4LogicalVolume(Crystal_box, crystalMaterial, "Crystal", 0, 0, 0);

  G4VPhysicalVolume* Crystal_phys
   = new G4PVPlacement(0, G4ThreeVector(0,0,0), Crystal_log, "Crystal", logicWorld, 
       0, false, 0);

  fRadiatorDescription->fLogicalVolume = Crystal_log;

  // create region for window inside windowR for

  G4Region* radRegion = new G4Region("XTRradiator");
  radRegion->AddRootLogicalVolume(Crystal_log);

  
  //--------------------------------------------
  // Absorber = detector

  /*G4VSolid* solidAbsorber
      = new G4Box("Absorber", crystalThickness / 2., crystalThickness / 2.,
          absorberThickness / 2.);*/
  G4VSolid* solidAbsorber
      = new G4Box("Absorber", 0.09*mm, 0.09*mm,
      //= new G4Box("Absorber", 5.0 * mm, 5.0 * mm,

          absorberThickness / 2.);

  G4LogicalVolume* logicAbsorber
      = new G4LogicalVolume(solidAbsorber, galactic, "Absorber");

  G4double absorberZ = crystalThickness / 2 + absorberThickness / 2. + 10e-3*mm;
  //G4double absorberZ = crystalThickness / 2 + absorberThickness / 2. + 10e-2 * mm;

  new G4PVPlacement(0, 
      G4ThreeVector(0, 0, absorberZ),
      "Absorber", logicAbsorber, physicsWorld, false, 0);

  G4Region* regGasDet = new G4Region("XTRdEdxDetector");
  regGasDet->AddRootLogicalVolume(logicAbsorber);

  // Sensitive Detectors: Absorber

  SensitiveDetector* sd = new SensitiveDetector("AbsorberSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  Crystal_log->SetSensitiveDetector(sd); //logicAbsorber->SetSensitiveDetector(sd);
   //---------------------------------------------------------------------------------------------
 /* G4VSolid* solidWorld 
    = new G4Box("World", worldSizeR, worldSizeR, worldSizeZ/2.);
 
  G4LogicalVolume* logicWorld 
    = new G4LogicalVolume(solidWorld,  worldMaterial,  "World");

  G4VPhysicalVolume* physicsWorld 
    = new G4PVPlacement(0, G4ThreeVector(), "World", logicWorld, 0,  false, 0);*/

  // CR radiator envelope

  //----------------------
  // Print geometry parameters

  /*G4cout << "\n The  WORLD   is made of "
         << worldSizeZ/mm << "mm of " << worldMaterial->GetName();
  G4cout << ", the transverse size (R) of the world is " 
         << worldSizeR/mm << " mm. " << G4endl;
  G4cout << " The ABSORBER is made of "
         << absorberThickness/mm << "mm of " << absorberMaterial->GetName();
  G4cout << ", the transverse size (R) is " 
         << absorberRadius/mm << " mm. " << G4endl;
  G4cout << " Z position of the (middle of the) absorber " 
         << absorberZ/mm << "  mm." << G4endl;

  G4cout << "radZ = " << radZ/mm << " mm" << G4endl;
  G4cout << "startZ = " << startZ/mm<< " mm" << G4endl;

  G4cout << "fRadThick = " << radThick/mm << " mm"<<G4endl;
  G4cout << "fFoilNumber = " << foilNumber << G4endl;
  G4cout << "fRadiatorMat = " << radiatorMat->GetName() << G4endl;
  G4cout << "WorldMaterial = " << worldMaterial->GetName() << G4endl;*/
  G4cout << "Hi! It's DetectorSimpleAlice.cc" << G4endl;
  G4cout << G4endl;

  return physicsWorld;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
