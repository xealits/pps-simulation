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
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Description: Test of Continuous Process G4Cerenkov
//              and RestDiscrete Process G4Scintillation
//              -- Generation Cerenkov Photons --
//              -- Generation Scintillation Photons --
//              -- Transport of optical Photons --
// Version:     5.0
// Created:     1996-04-30
// Author:      Juliet Armstrong
// mail:        gum@triumf.ca
//     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "G4RunManager.hh"
#include "G4UImanager.hh"

G4long seed;

#include "ExN06PhysicsList.hh"
#include "ExN06PrimaryGeneratorAction.hh"
#include "ExN06DetectorConstruction.hh"
#include "ExN06RunAction.hh"
#include "ExN06StackingAction.hh"
#include "ExN06SteppingVerbose.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


#include "stdio.h"

// #include "seed.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Seed the random number generator manually
  //
  // G4long myseed = 159929;
  // G4long myseed = 144563;
  // CLHEP::HepRandom::setTheSeed(myseed);
  if (argc > 1)
    {
      sscanf( argv[1], "%lu", &seed ); // so unsighed long is casted to G4long !
      // also the seed is global -- it is included to the filename of the output in include/LeadSD.hh (the definition of photodetector behaviour)
      printf("SEED:\n%lu\n---------------------------------------------\n", seed);
    }
  else return 404;
  CLHEP::HepRandom::setTheSeed(seed);

  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new ExN06SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // Run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // UserInitialization classes - mandatory
  //
  G4VUserPhysicsList* physics = new ExN06PhysicsList;
  runManager-> SetUserInitialization(physics);
  //
  G4VUserPrimaryGeneratorAction* gen_action = new ExN06PrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);
  //
  G4VUserDetectorConstruction* detector = new ExN06DetectorConstruction;
  runManager-> SetUserInitialization(detector);
  
#ifdef G4VIS_USE
  // visualization manager
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // UserAction classes
  //
  G4UserRunAction* run_action = new ExN06RunAction;
  runManager->SetUserAction(run_action);
  //
  G4UserStackingAction* stacking_action = new ExN06StackingAction;
  runManager->SetUserAction(stacking_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer(); 
   
  if (argc<=2)   // Define UI session for interactive mode
    // argc should == 2, one argument -- the seed -- is necessary?
    {
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
    }
  else         // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[2]; // the first argument is seed, the macro file goes next
      UImanager->ApplyCommand(command+fileName);
    }
   
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  delete verbosity;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
