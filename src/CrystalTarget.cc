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
/// \file electromagnetic/TestEm10/src/CrystalTarget.cc
/// \brief Implementation of the CrystalTarget class
// 
//

#include <complex>
#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include "CrystalTarget.hh"
#include "Randomize.hh"
//#include "G4integrator.hh"
#include "G4Gamma.hh"
#include "G4PhysicalConstants.hh"





 #include "phys_const.hh"
 //#include "crystal_functions.hh"
 //using namespace std;

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

CrystalTarget::CrystalTarget(G4LogicalVolume *anEnvelope,
                                         G4Material* RadlMat,
                                         G4double a, G4double EGamma,
                                         const G4String& processName) :
  G4ChanRad(anEnvelope, RadlMat, a, EGamma , processName)
{
  G4cout<<"Regular transparent X-ray TR  radiator EM process is called"<<G4endl;

  // Build energy and angular G4integral spectra of X-ray TR photons from
  // a radiator
  // fExitFlux   = true;
  fAlphaPlate = 10000;
  fAlphaGas   = 1000;

  //  BuildTable();
}

///////////////////////////////////////////////////////////////////////////

CrystalTarget::~CrystalTarget()
{
  ;
}

///////////////////////////////////////////////////////////////////////////
//
//

G4double CrystalTarget::SpectralXTRdEdx(G4int i, G4double energy, vector<ld1> emission, ld1 min_TT)
{
  // G4double result, sum = 0., tmp, cof1, cof2, cofMin, cofPHC,aMa, bMb, sigma,mu;
  // G4int k, kMax, kMin;
	G4double result;
  // sigma = 0.01;
  // mu = 0.05;
  // result = (1/sigma/sqrt(2*pi))*exp(-((energy-mu)*(energy-mu))/(2*sigma*sigma));
  // //my changes
    
	//emission.w.clear();
    //rad_counter_w = 0;
    G4double sum_dwde = 0;
    //G4double coeffi = 4*pi*pow(gama, 2) / min_TT;
    //vector<ld1> radi;
   // radi.clear();
    // std::ofstream total("total_rad_new.txt");
    // for (G4int i = 0; i < 100; i++)
    // {
       // emission.w.push_back(rad_counter_w * hbar);
        //rad_counter_w += 0.01 * coeffi;
		cout<<"i = "<<i<<" energy = "<<energy<<"  emis = "<<emission.at(i)<<endl;
        for (G4int j = i; j < nbeam*(100); j+=100)
        {	
            sum_dwde +=  emission.at(j);//*crystaltickness*pow(10,-4)/nbeam;
        }
        //radi.push_back(sum_dwde);
        //total << emission.w.back() << " " << /*fixed << setprecision(6) <<*/ radi.back()/nbeam << "\n";
        //sum_dwde = 0;
    // }
   
    //fout89.close();
    //total.close();
	result = sum_dwde;
  return result;
}
//
//
////////////////////////////////////////////////////////////////////////////








