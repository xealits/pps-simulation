
///
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

#include "ExN06DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4UnionSolid.hh"
#include "G4Trd.hh"


#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "LeadSD.hh"
#include "AirSD.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::ExN06DetectorConstruction()
{
  expHall_x = expHall_y = expHall_z = 10.0*cm;
  tank_x    = tank_y    = tank_z    =  5.0*cm;
  bubble_x  = bubble_y  = bubble_z  =  0.5*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::~ExN06DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExN06DetectorConstruction::Construct()
{




  /* -------------------------------------------------------------------------------------------
  // Generate & Add Material Properties Table
  // Optical Properties
  // ------------------------------------------------------------------------------------------- */


  // Quartz refraction from:
  // I. H. Malitson. Interspecimen Comparison of the Refractive Index of Fused Silica
  // at http://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson

  const G4int QuartzRefractionEntries = 101;
  G4double QuartzPhotonEnergy[QuartzRefractionEntries] =
  { 5.90476190476191*eV, 5.06122448979592*eV, 4.42857142857143*eV, 3.93650793650794*eV, 3.54285714285714*eV,
    3.22077922077922*eV, 2.95238095238095*eV, 2.72527472527473*eV, 2.53061224489796*eV, 2.36190476190476*eV,
    2.21428571428571*eV, 2.08403361344538*eV, 1.96825396825397*eV, 1.86466165413534*eV, 1.77142857142857*eV,
    1.68707482993197*eV, 1.61038961038961*eV, 1.54037267080745*eV, 1.47619047619048*eV, 1.41714285714286*eV,
    1.36263736263736*eV, 1.31216931216931*eV, 1.26530612244898*eV, 1.22167487684729*eV, 1.18095238095238*eV,
    1.14285714285714*eV, 1.10714285714286*eV, 1.07359307359307*eV, 1.04201680672269*eV, 1.01224489795918*eV,
    0.984126984126984*eV, 0.957528957528958*eV, 0.932330827067669*eV, 0.908424908424908*eV, 0.885714285714286*eV,
    0.86411149825784*eV, 0.843537414965986*eV, 0.823920265780731*eV, 0.805194805194805*eV, 0.787301587301587*eV,
    0.770186335403727*eV, 0.753799392097264*eV, 0.738095238095238*eV, 0.723032069970845*eV, 0.708571428571429*eV,
    0.694677871148459*eV, 0.681318681318681*eV, 0.668463611859838*eV, 0.656084656084656*eV, 0.644155844155844*eV,
    0.63265306122449*eV, 0.621553884711779*eV, 0.610837438423645*eV, 0.600484261501211*eV, 0.59047619047619*eV,
    0.580796252927401*eV, 0.571428571428571*eV, 0.562358276643991*eV, 0.553571428571428*eV, 0.545054945054945*eV,
    0.536796536796537*eV, 0.528784648187633*eV, 0.521008403361345*eV, 0.513457556935818*eV, 0.506122448979592*eV,
    0.498993963782696*eV, 0.492063492063492*eV, 0.4853228962818*eV, 0.478764478764479*eV, 0.472380952380952*eV,
    0.466165413533835*eV, 0.460111317254174*eV, 0.454212454212454*eV, 0.448462929475588*eV, 0.442857142857143*eV,
    0.437389770723104*eV, 0.43205574912892*eV, 0.426850258175559*eV, 0.421768707482993*eV, 0.416806722689076*eV,
    0.411960132890365*eV, 0.407224958949097*eV, 0.402597402597403*eV, 0.398073836276083*eV, 0.393650793650794*eV,
    0.389324960753532*eV, 0.385093167701863*eV, 0.380952380952381*eV, 0.376899696048632*eV, 0.372932330827068*eV,
    0.369047619047619*eV, 0.365243004418262*eV, 0.361516034985423*eV, 0.357864357864358*eV, 0.354285714285714*eV,
    0.350777934936351*eV, 0.34733893557423*eV, 0.343966712898752*eV, 0.340659340659341*eV, 0.337414965986395*eV,
    0.334231805929919*eV };



  // the ones from the site http://refractiveindex.info/
  // Optical constants of SiO2 (Silicon dioxide, Silica, Quartz)
  // Malitson 1965 - Fused silica; n 0.21-3.71 µm

  G4double QuartzRefractiveIndex[QuartzRefractionEntries] =
  { 1.5383576204905, 1.5102724365895, 1.4941636611188, 1.4839008951423, 1.476891413496,
    1.4718556531995, 1.4680936900401, 1.46519309996, 1.4628966820387, 1.4610366660574,
    1.4594995356592, 1.458206104926, 1.4570996888769, 1.4561387969803, 1.4552924662623,
    1.454537192876, 1.4538548630589, 1.4532313266004, 1.4526553936728, 1.4521181167939,
    1.451612268629, 1.4511319566977, 1.4506723353353, 1.4502293877559, 1.4497997593263,
    1.4493806287126, 1.4489696073537, 1.4485646603469, 1.4481640436751, 1.4477662540207,
    1.4473699883562, 1.4469741111889, 1.4465776278426, 1.4461796625343, 1.4457794402848,
    1.4453762719132, 1.4449695415266, 1.4445586960405, 1.4441432363602, 1.4437227099274,
    1.4432967043935, 1.4428648422301, 1.4424267761167, 1.4419821849822, 1.4415307705927,
    1.4410722546016, 1.4406063759896, 1.440132888836, 1.4396515603723, 1.4391621692763,
    1.438664504173, 1.438158362312, 1.4376435483981, 1.4371198735528, 1.4365871543902,
    1.4360452121918, 1.435493872166, 1.4349329627838, 1.4343623151781, 1.4337817626007,
    1.4331911399286, 1.4325902832137, 1.431979029271, 1.4313572152994, 1.4307246785328,
    1.4300812559156, 1.4294267838019, 1.4287610976737, 1.4280840318763, 1.4273954193693,
    1.4266950914905, 1.4259828777305, 1.4252586055186, 1.4245221000158, 1.4237731839158,
    1.4230116772521, 1.4222373972096, 1.4214501579411, 1.4206497703869, 1.4198360420966,
    1.4190087770532, 1.4181677754984, 1.4173128337577, 1.4164437440668, 1.415560294396,
    1.414662268275, 1.4137494446147, 1.412821597528, 1.411878496148, 1.4109199044426,
    1.4099455810268, 1.4089552789703, 1.4079487456017, 1.4069257223074, 1.4058859443263,
    1.4048291405385, 1.4037550332482, 1.402663337961, 1.401553763154, 1.4004260100389,
    1.3992797723176 };


/*
  G4double QuartzRefractiveIndex[QuartzRefractionEntries] =
  { 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3, 3, 3, 3, 3,
    3 };
*/

/*
  G4double QuartzRefractiveIndexRoof[QuartzRefractionEntries] =
  { 1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8 };
*/



  const G4int refractiveEntries = 101;
  const G4int nEntries = 32;

  G4double PhotonEnergy[refractiveEntries] =
  { 5.90476190476191*eV, 5.06122448979592*eV, 4.42857142857143*eV, 3.93650793650794*eV, 3.54285714285714*eV,
    3.22077922077922*eV, 2.95238095238095*eV, 2.72527472527473*eV, 2.53061224489796*eV, 2.36190476190476*eV,
    2.21428571428571*eV, 2.08403361344538*eV, 1.96825396825397*eV, 1.86466165413534*eV, 1.77142857142857*eV,
    1.68707482993197*eV, 1.61038961038961*eV, 1.54037267080745*eV, 1.47619047619048*eV, 1.41714285714286*eV,
    1.36263736263736*eV, 1.31216931216931*eV, 1.26530612244898*eV, 1.22167487684729*eV, 1.18095238095238*eV,
    1.14285714285714*eV, 1.10714285714286*eV, 1.07359307359307*eV, 1.04201680672269*eV, 1.01224489795918*eV,
    0.984126984126984*eV, 0.957528957528958*eV, 0.932330827067669*eV, 0.908424908424908*eV, 0.885714285714286*eV,
    0.86411149825784*eV, 0.843537414965986*eV, 0.823920265780731*eV, 0.805194805194805*eV, 0.787301587301587*eV,
    0.770186335403727*eV, 0.753799392097264*eV, 0.738095238095238*eV, 0.723032069970845*eV, 0.708571428571429*eV,
    0.694677871148459*eV, 0.681318681318681*eV, 0.668463611859838*eV, 0.656084656084656*eV, 0.644155844155844*eV,
    0.63265306122449*eV, 0.621553884711779*eV, 0.610837438423645*eV, 0.600484261501211*eV, 0.59047619047619*eV,
    0.580796252927401*eV, 0.571428571428571*eV, 0.562358276643991*eV, 0.553571428571428*eV, 0.545054945054945*eV,
    0.536796536796537*eV, 0.528784648187633*eV, 0.521008403361345*eV, 0.513457556935818*eV, 0.506122448979592*eV,
    0.498993963782696*eV, 0.492063492063492*eV, 0.4853228962818*eV, 0.478764478764479*eV, 0.472380952380952*eV,
    0.466165413533835*eV, 0.460111317254174*eV, 0.454212454212454*eV, 0.448462929475588*eV, 0.442857142857143*eV,
    0.437389770723104*eV, 0.43205574912892*eV, 0.426850258175559*eV, 0.421768707482993*eV, 0.416806722689076*eV,
    0.411960132890365*eV, 0.407224958949097*eV, 0.402597402597403*eV, 0.398073836276083*eV, 0.393650793650794*eV,
    0.389324960753532*eV, 0.385093167701863*eV, 0.380952380952381*eV, 0.376899696048632*eV, 0.372932330827068*eV,
    0.369047619047619*eV, 0.365243004418262*eV, 0.361516034985423*eV, 0.357864357864358*eV, 0.354285714285714*eV,
    0.350777934936351*eV, 0.34733893557423*eV, 0.343966712898752*eV, 0.340659340659341*eV, 0.337414965986395*eV,
    0.334231805929919*eV };

/*
  G4double PhotonEnergy[nEntries] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
*/

//
// Water
//
/*
  G4double RefractiveIndex1[nEntries] =
            { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
              1.346,  1.3465, 1.347,  1.3475, 1.348,
              1.3485, 1.3492, 1.35,   1.3505, 1.351,
              1.3518, 1.3522, 1.3530, 1.3535, 1.354,
              1.3545, 1.355,  1.3555, 1.356,  1.3568,
              1.3572, 1.358,  1.3585, 1.359,  1.3595,
              1.36,   1.3608};
*/

  G4double RefractiveIndex1[refractiveEntries] =
  { 1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8, 1.8, 1.8, 1.8, 1.8,
    1.8 };

/*
  G4double RefractiveIndex1[nEntries] =
            { 1.9, 1.9, 1.9, 1.9,  1.9,
              1.9, 1.9, 1.9, 1.9,  1.9,
              1.9, 1.9, 1.9, 1.9,  1.9,
              1.9, 1.9, 1.9, 1.9,  1.9,
              1.9, 1.9, 1.9, 1.9,  1.9,
              1.9, 1.9, 1.9, 1.9,  1.9,
              1.9, 1.9};
*/

  G4double Absorption1[nEntries] =
           {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
           15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
           45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
           52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
           30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
           17.500*m, 14.500*m };

  G4double Absorption2[nEntries] =
           {.3*mm,  .4*mm,  .6*mm,  .9*mm, 1.2*mm, 1.3*mm,
           1.5*mm, 1.7*mm, 1.8*mm, 2.0*mm, 2.6*mm, 3.5*mm,
           4.5*mm, 4.7*mm, 5.2*mm, 5.2*mm, 5.5*mm, 5.2*mm,
           5.2*mm, 4.7*mm, 4.5*mm, 4.1*mm, 3.7*mm, 3.3*mm,
           3.0*mm, 2.8*mm, 2.7*mm, 2.4*mm, 2.2*mm, 1.9*mm,
           1.7*mm, 1.4*mm };


  G4double AbsorptionQuartz[nEntries] =
           {10*m, 10*m, 10*m, 10*m, 10*m, 10*m,
            10*m, 10*m, 10*m, 10*m, 10*m, 10*m,
            10*m, 10*m, 10*m, 10*m, 10*m, 10*m,
            10*m, 10*m, 10*m, 10*m, 10*m, 10*m,
            10*m, 10*m, 10*m, 10*m, 10*m, 10*m,
            10*m, 10*m };
  
            /*
           {1*m, 1*m, 1*m, 1*m, 1*m, 1*m,
            1*m, 1*m, 1*m, 1*m, 1*m, 1*m,
            1*m, 1*m, 1*m, 1*m, 1*m, 1*m,
            1*m, 1*m, 1*m, 1*m, 1*m, 1*m,
            1*m, 1*m, 1*m, 1*m, 1*m, 1*m,
            1*m, 1*m };
            */

/*
  G4double ScintilFast[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4double ScintilSlow[nEntries] =
            { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
              7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
              3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
              4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
              7.00, 6.00, 5.00, 4.00 };
*/

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1, refractiveEntries)
       ->SetSpline(true);
  // myMPT1->AddConstProperty("RINDEX", 2); // DOES NOT WORK

  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries)
       ->SetSpline(true);


  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();

  myMPT2->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1, refractiveEntries)
       ->SetSpline(true);
  // myMPT2->AddConstProperty("RINDEX", 2); // DOES NOT WORK

  myMPT2->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption2,     nEntries)
       ->SetSpline(true);





  const G4int NumberDetOpticEntries = 2;
  G4double DetOpticPhotonEnergies[NumberDetOpticEntries] = { 0.1*eV, 6*eV };

  G4double DetAbsorbtions[NumberDetOpticEntries] = { 0.1*mm, 0.1*mm };

  G4double DetRefractions[NumberDetOpticEntries] = { 5, 5 };







  G4MaterialPropertiesTable* SiO2_MPTdet = new G4MaterialPropertiesTable();

  // the Quart refractive index
  SiO2_MPTdet->AddProperty("RINDEX",       QuartzPhotonEnergy, QuartzRefractiveIndex, QuartzRefractionEntries)
       ->SetSpline(true);


  // SiO2_MPTdet->AddConstProperty("RINDEX", 2); // DOES NOT WORK

  // SiO2_MPTdet->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption2,     nEntries)
  SiO2_MPTdet->AddProperty("ABSLENGTH",    DetOpticPhotonEnergies, DetAbsorbtions, NumberDetOpticEntries)
       ->SetSpline(true);
  
  // myMPT2->AddConstProperty("ABSLENGTH",    0.001*mm);



  G4MaterialPropertiesTable* Sapphire_MPTdet = new G4MaterialPropertiesTable();

  // the Sapphire refractive index
  // 1.8 index
  Sapphire_MPTdet->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1, refractiveEntries)
       ->SetSpline(true);

  Sapphire_MPTdet->AddProperty("ABSLENGTH",    DetOpticPhotonEnergies, DetAbsorbtions, NumberDetOpticEntries)
       ->SetSpline(true);




  // Quartz Properties

  G4MaterialPropertiesTable* quartzMPT = new G4MaterialPropertiesTable();


  // refractive index 1.5
  // quartzMPT->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1, refractiveEntries)
       // ->SetSpline(true);

  quartzMPT->AddProperty("RINDEX",       QuartzPhotonEnergy, QuartzRefractiveIndex, QuartzRefractionEntries)
       ->SetSpline(true);
  // quartzMPT->AddConstProperty("RINDEX", 2); // DOES NOT WORK

  quartzMPT->AddProperty("ABSLENGTH",    QuartzPhotonEnergy, AbsorptionQuartz,     nEntries)
       ->SetSpline(true);




  // Sapphire Properties

  G4MaterialPropertiesTable* sapphireMPT = new G4MaterialPropertiesTable();


  // refractive index 1.8
  sapphireMPT->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1, refractiveEntries)
       ->SetSpline(true);

  // sapphireMPT->AddProperty("RINDEX",       saphirePhotonEnergy, saphireRefractiveIndex, saphireRefractionEntries)
       // ->SetSpline(true);
  // sapphireMPT->AddConstProperty("RINDEX", 2); // DOES NOT WORK

  sapphireMPT->AddProperty("ABSLENGTH",    QuartzPhotonEnergy, AbsorptionQuartz,     nEntries)
       ->SetSpline(true);




/*
  myMPT1->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,     nEntries)
        ->SetSpline(true);
  myMPT1->AddProperty("SLOWCOMPONENT",PhotonEnergy, ScintilSlow,     nEntries)
        ->SetSpline(true);
  
  myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);
*/

/*
  const G4int NUMENTRIES_water = 60;

  G4double ENERGY_water[NUMENTRIES_water] = {
     1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
     1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
     1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
     1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
     1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
     2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
     2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
     2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
     2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
     2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
     3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
     3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
     3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
     4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
     5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
  };

  //assume 100 times larger than the rayleigh scattering for now.
  G4double MIE_water[NUMENTRIES_water] = {
     167024.4*m, 158726.7*m, 150742  *m,
     143062.5*m, 135680.2*m, 128587.4*m,
     121776.3*m, 115239.5*m, 108969.5*m,
     102958.8*m, 97200.35*m, 91686.86*m,
     86411.33*m, 81366.79*m, 76546.42*m,
     71943.46*m, 67551.29*m, 63363.36*m,
     59373.25*m, 55574.61*m, 51961.24*m,
     48527.00*m, 45265.87*m, 42171.94*m,
     39239.39*m, 36462.50*m, 33835.68*m,
     31353.41*m, 29010.30*m, 26801.03*m,
     24720.42*m, 22763.36*m, 20924.88*m,
     19200.07*m, 17584.16*m, 16072.45*m,
     14660.38*m, 13343.46*m, 12117.33*m,
     10977.70*m, 9920.416*m, 8941.407*m,
     8036.711*m, 7202.470*m, 6434.927*m,
     5730.429*m, 5085.425*m, 4496.467*m,
     3960.210*m, 3473.413*m, 3032.937*m,
     2635.746*m, 2278.907*m, 1959.588*m,
     1675.064*m, 1422.710*m, 1200.004*m,
     1004.528*m, 833.9666*m, 686.1063*m
  };

  // gforward, gbackward, forward backward ratio
  G4double MIE_water_const[3]={0.99,0.99,0.8};
*/

/*
  myMPT1->AddProperty("MIEHG",ENERGY_water,MIE_water,NUMENTRIES_water)
        ->SetSpline(true);
  myMPT1->AddConstProperty("MIEHG_FORWARD",MIE_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD",MIE_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",MIE_water_const[2]);
*/


  // Set the Birks Constant for the Water scintillator

  //Water->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  //
  // Air
  //

/*
  G4double RefractiveIndex2[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };
*/

  G4double RefractiveIndex2[nEntries] =
            { 1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1 };

  G4MaterialPropertiesTable* airMPT = new G4MaterialPropertiesTable();
  airMPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
  


  const G4int NUMENTRIESmat = 2;
  // G4double mat_PP[NUMENTRIESmat] = {2.07*eV, 3.28*eV};
  G4double mat_PP[NUMENTRIESmat] = {1.07*eV, 6.28*eV};
 
  const G4int nEntriesV = 34;

  G4double PhotonEnergyV[nEntriesV] = {
      0.0001 * eV, 1.0 * eV, 2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV,
      2.177 * eV, 2.216 * eV, 2.256 * eV, 2.298 * eV,
      2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV,
      2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV,
      2.757 * eV, 2.820 * eV, 2.885 * eV, 2.954 * eV,
      3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV,
      3.353 * eV, 3.446 * eV, 3.545 * eV, 3.649 * eV,
      3.760 * eV, 3.877 * eV, 4.002 * eV, 4.136 * eV
  };

  G4double AbsorptionV[nEntriesV] =
         { 3.448*m, 3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
          15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
          45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
          52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
          30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
          17.500*m, 14.500*m };


  G4MaterialPropertiesTable *VikuitiPT = new G4MaterialPropertiesTable();
  G4double RINDEXvikuiti[NUMENTRIESmat] = {1.3, 1.3};
  VikuitiPT -> AddProperty("RINDEX", mat_PP, RINDEXvikuiti, NUMENTRIESmat);
  VikuitiPT -> AddProperty("ABSLENGTH", PhotonEnergyV, AbsorptionV, nEntriesV);












































  /* -------------------------------------------------------------------------------------------
  // Materials
  // ------------------------------------------------------------------------------------------- */


  G4double a, z, density;
  G4int nelements;


  // Air
  // 
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  Air->SetMaterialPropertiesTable(airMPT);







/*
  // Water
  // 
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* Water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  Water->AddElement(H, 2);
  Water->AddElement(O, 1);

  Water->SetMaterialPropertiesTable(myMPT1);
*/






  G4String name, symbol;
  G4int ncomponents, natoms;

  // Quartz with refraction
  //

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 28.09*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a);

  density = 2.200*g/cm3;
  // density = 1.0*g/cm3;
  G4Material* SiO2 = new G4Material(name="quartz", density, ncomponents=2);
  SiO2->AddElement(elSi, natoms=1);
  SiO2->AddElement(elO , natoms=2);

  G4NistManager* nist = G4NistManager::Instance();
  // G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  SiO2->SetMaterialPropertiesTable(quartzMPT);
  // SiO2->SetMaterialPropertiesTable(myMPT1);
  // SiO2->SetMaterialPropertiesTable(myMPT2);


  G4Material* Sapphire = new G4Material(name="Sapphire", density, ncomponents=2);
  Sapphire->AddElement(elSi, natoms=1); // WRONG
  Sapphire->AddElement(elO , natoms=2); // WRONG
  // only optical density is ok

  Sapphire->SetMaterialPropertiesTable(sapphireMPT);




  // Vikuiti

  G4Element *H = new G4Element ("Hydrogen", "H", z = 1 , a = 1.01 * g / mole);
  G4Element *C = new G4Element ("Carbon"  , "C", z = 6 , a = 12.01 * g / mole);
  G4Element *Ox = new G4Element ("Oxigen", "Ox", z = 8 , a = 16 * g / mole);

  G4double nel;

  G4Material * Vikuiti = new G4Material("Vikuiti", density=1*g/cm3, nel = 3);
  Vikuiti->AddElement (C, 1);
  Vikuiti->AddElement (H, 1);
  Vikuiti->AddElement (Ox, 1);

  Vikuiti -> SetMaterialPropertiesTable(VikuitiPT);









































  /* -------------------------------------------------------------------------------------------
  // Geometry
  // ------------------------------------------------------------------------------------------- */

//
// ------------- Volumes --------------

  G4bool checkOverlaps = true;







// The experimental Hall
//
  //G4Box* expHall_box = new G4Box("World",expHall_x,expHall_y,expHall_z);
  // G4Box* expHall_box = new G4Box("World", 15*cm, 15*cm, 15*cm);
  G4Box* expHall_box = new G4Box("World", 15*cm, 15*cm, 35*cm);


  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box, Air, "World"); //,0,0,0);

/*
  G4LogicalVolume* expHall_log
    // = new G4LogicalVolume(expHall_box, Vikuiti, "World"); //,0,0,0);
    = new G4LogicalVolume(expHall_box, Air, "World"); //,0,0,0);
    //= new G4LogicalVolume(expHall_box,Water,"World",0,0,0);
*/

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,
                        G4ThreeVector(),
                        expHall_log,
                        "World",
                        0,
                        false,
                        0,
                        checkOverlaps);





























/*
  // L-SHAPE

  G4double length = 18*mm;
  G4double hight  = 58.8*mm;
  G4double breadth = 3*mm;

  // G4double width = 0.5*mm;
  G4double det_overlap = 0*mm;

  G4double width = breadth/2;
  G4double side = length/2;
  G4double top = hight/2;

  G4double det_width =  width + 0.5*mm;

  // Tank1, Tank2, their Union

  G4Box* Tank1_box = new G4Box("Tank1b", width, width, side - width);
  G4Box* Tank2_box = new G4Box("Tank2b", width, top, width);

  G4RotationMatrix* RotNul = new G4RotationMatrix;
  G4ThreeVector pos2 = G4ThreeVector(0, top - width, side);

  G4UnionSolid* Tank_union =
    new G4UnionSolid("L_shape", Tank1_box, Tank2_box, RotNul, pos2); 

  G4LogicalVolume* TankUni_log
    // = new G4LogicalVolume(Tank_union, SiO2,"TankL"); //,0,0,0);
    = new G4LogicalVolume(Tank_union, Sapphire, "TankL"); //,0,0,0);


  G4VPhysicalVolume* Tank_phys
    = new G4PVPlacement(0,
                        G4ThreeVector(0, 0, side - width),
                        TankUni_log,
                        "TankL",
                        expHall_log,
                        false,
                        0,
                        checkOverlaps);


  // Detector box

  G4ThreeVector posDetector = G4ThreeVector(0,  hight, length - width);

  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_GLASS_LEAD");
  // detector_mat->SetMaterialPropertiesTable(SiO2_MPTdet);
  detector_mat->SetMaterialPropertiesTable(Sapphire_MPTdet);

  G4Box * DetectorBox =
    new G4Box("SolidShapeD",
      det_width, width + det_overlap, det_width);

  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(DetectorBox,         //its solid
                        detector_mat,          //its material
                        "Detector");           //its name

  G4VPhysicalVolume* Det_phys =
    new G4PVPlacement(0,
                    posDetector,
                    logicDetector,
                    "Detector",
                    expHall_log,
                    false,
                    0,
                    checkOverlaps);
*/














/*
// L-SHAPE, separate bars

  G4double length = 100*mm;
  G4double hight  = 10*mm;
  G4double breadth = 3*mm;

  // G4double width = 0.5*mm;
  G4double det_overlap = 0*mm;

  G4double width = breadth/2;
  G4double side = length/2;
  G4double top = hight/2;

  G4double det_width =  width + 0.5*mm;

  // Tank2
  // TankH, TankV

  G4Box* TankH_box = new G4Box("TankHb", width, width, side - width);
  G4Box* TankV_box = new G4Box("TankVb", width, top, width);

  // G4RotationMatrix* RotNul = new G4RotationMatrix;
  G4ThreeVector pos2 = G4ThreeVector(0, top - width, side);

  // G4UnionSolid* Tank_union =
    // new G4UnionSolid("L_shape", TankH_box, TankV_box, RotNul, pos2); 

  G4LogicalVolume* TankH_log
    = new G4LogicalVolume(TankH_box, SiO2,"TankL"); //,0,0,0);

  G4LogicalVolume* TankV_log
    = new G4LogicalVolume(TankV_box, Sapphire,"TankL"); //,0,0,0);
    // = new G4LogicalVolume(TankV_box, SiO2,"TankL"); //,0,0,0);


  G4VPhysicalVolume* Tank_phys
    = new G4PVPlacement(0,
                        G4ThreeVector(0, 0, side - width),
                        TankH_log,
                        "TankL",
                        expHall_log,
                        false,
                        0,
                        checkOverlaps);


  G4VPhysicalVolume* TankV_phys
    = new G4PVPlacement(0,
                        G4ThreeVector(0, top - width, 2*side - width),
                        TankV_log,
                        "TankL",
                        expHall_log,
                        false,
                        0,
                        checkOverlaps);



  // Detector box

  G4ThreeVector posDetector = G4ThreeVector(0,  hight, length - width);

  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_GLASS_LEAD");
  detector_mat->SetMaterialPropertiesTable(SiO2_MPTdet);

  G4Box * DetectorBox =
    new G4Box("SolidShapeD",
      det_width, width + det_overlap, det_width);

  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(DetectorBox,         //its solid
                        detector_mat,          //its material
                        "Detector");           //its name

  G4VPhysicalVolume* Det_phys =
    new G4PVPlacement(0,
                    posDetector,
                    logicDetector,
                    "Detector",
                    expHall_log,
                    false,
                    0,
                    checkOverlaps);

*/















/*
  // HORIZONTAL BOX

  // Tank

  G4double length = 55.6*mm;
  G4double breadth = 3*mm;
  G4double det_overlap = 0*mm;

  G4double width = breadth/2;
  G4double side = (length)/2;

  G4double det_width = width + 0.5*mm;

  G4Box* Tank1_box = new G4Box("Tank1b", width, width, side);
  G4ThreeVector posTank = G4ThreeVector(0,  0*mm, side);
  G4ThreeVector posDetector = G4ThreeVector(0,  0*mm, (2*side + width));


  G4LogicalVolume* Tank1_log
    // = new G4LogicalVolume(Tank1_box, SiO2, "TankL1"); //,0,0,0);
    = new G4LogicalVolume(Tank1_box, Sapphire, "TankL1"); //,0,0,0);



  G4VPhysicalVolume* Tank_phys
    = new G4PVPlacement(0, posTank, Tank1_log, "Tank1",
                        expHall_log, false, 0, checkOverlaps);


  // Detector box

  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_GLASS_LEAD");
  detector_mat->SetMaterialPropertiesTable(SiO2_MPTdet);

  G4Box * DetectorBox =
    new G4Box("SolidShapeD",
      det_width, det_width, width + det_overlap);

  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(DetectorBox,         //its solid
                        detector_mat,          //its material, it has the same refractive index as the quartz -- to absorb all the photons
                        "Detector");           //its name

  new G4PVPlacement(0,
                    posDetector,
                    logicDetector,
                    "Detector",
                    expHall_log,
                    false,
                    0,
                    checkOverlaps);
*/




















  // VERTICAL BOX

  // Tank

  G4double hight = 50*mm;
  G4double det_overlap = 0*mm;
  G4double breadth = 3*mm;

  G4double width = breadth/2;
  G4double forward_width = width;

  G4double top = hight/2;

  G4double det_width = width + 0.5*mm;
  G4double det_forward_width = forward_width + 0.5*mm;

  G4Box* Tank2_box = new G4Box("Tank2b", width, top, forward_width);

  G4ThreeVector posBack = G4ThreeVector(0,  top - width, width);

  G4ThreeVector posDetector = G4ThreeVector(0,  2*top, width);

  G4LogicalVolume* Tank2_log
    = new G4LogicalVolume(Tank2_box, SiO2,"TankL2"); //,0,0,0);
    // = new G4LogicalVolume(Tank2_box,SiO2,"Tank"); //,0,0,0);
    // = new G4LogicalVolume(Tank2_box,Water,"Tank"); //,0,0,0);

  G4VPhysicalVolume* Tank_phys
    = new G4PVPlacement(0,
                        posBack,
                        Tank2_log,
                        "Tank2",
                        expHall_log,
                        false,
                        0,
                        checkOverlaps);

  // Detector box

  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_GLASS_LEAD");
  detector_mat->SetMaterialPropertiesTable(SiO2_MPTdet);

  G4Box * DetectorBox =
    new G4Box("SolidShapeD",
      det_width, width + det_overlap, det_forward_width);

  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(DetectorBox,         //its solid
                        detector_mat,          //its material
                        "Detector");           //its name


  new G4PVPlacement(0,
                    posDetector,
                    logicDetector,
                    "Detector",
                    expHall_log,
                    false,
                    0,
                    checkOverlaps);


























/*
  // The Cut Piramid (aka Trapezoid)

  // Tank

  G4double hight = 40*mm;
  G4double det_overlap = 0*mm;
  G4double width_bot = 1.5*mm;
  G4double width_top = 1.6*mm;

  G4double top = hight/2;

  G4double det_width = width_top + 0.5*mm;

  // G4Box* Tank2_box = new G4Box("Tank2b", width, top, width);

  G4Trd* Tank2_box = new G4Trd("Tank2b",
             width_bot,
             width_top,
             width_bot,
             width_top,
             top);

  G4ThreeVector posBack = G4ThreeVector(0,  top - width_bot, width_top);

  G4ThreeVector posDetector = G4ThreeVector(0,  2*top - width_bot, width_top);

  G4LogicalVolume* Tank2_log
    = new G4LogicalVolume(Tank2_box, SiO2,"TankL2"); //,0,0,0);
    // = new G4LogicalVolume(Tank2_box,SiO2,"Tank"); //,0,0,0);
    // = new G4LogicalVolume(Tank2_box,Water,"Tank"); //,0,0,0);

  G4RotationMatrix * yRot90deg = new G4RotationMatrix();   // Rotates X and Z axes only

  // yRot90deg->rotateY(M_PI/2.*rad);
  yRot90deg->rotateX(90.*deg);


  G4VPhysicalVolume* Tank_phys
    = new G4PVPlacement(yRot90deg,
                        posBack,
                        Tank2_log,
                        "Tank2",
                        expHall_log,
                        false,
                        0,
                        checkOverlaps);

  // Detector box

  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_GLASS_LEAD");
  detector_mat->SetMaterialPropertiesTable(SiO2_MPTdet);

  G4Box * DetectorBox =
    new G4Box("SolidShapeD",
      det_width, width_top + det_overlap, det_width);

  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(DetectorBox,         //its solid
                        detector_mat,          //its material
                        "Detector");           //its name


  new G4PVPlacement(0,
                    posDetector,
                    logicDetector,
                    "Detector",
                    expHall_log,
                    false,
                    0,
                    checkOverlaps);
*/

















/*
  // The FOCUSED L-SHAPE

  // Tank

  G4double length = 20*mm;
  G4double hight = 40*mm;
  G4double det_overlap = 0*mm;
  // G4double breadth = 0.5*mm;  G4double breadth_top = 0.6*mm;
  // G4double breadth = 1*mm;  G4double breadth_top = 1.1*mm;
  // G4double breadth = 2*mm;  G4double breadth_top = 2.2*mm;
  // G4double breadth = 3*mm;  G4double breadth_top = 3.2*mm;
  G4double breadth = 4*mm;  G4double breadth_top = 4.2*mm;

  G4double width = breadth / 2;
  G4double width_bot = width;
  G4double width_top = breadth_top / 2;

  G4double top = hight/2;

  G4double det_width = width_top + 0.5*mm;

  G4double side = length/2;



  // Tank1, Tank2, their Union

  G4Box* Tank1_box = new G4Box("Tank1b", width, width, side );

  G4ThreeVector posTank = G4ThreeVector(0,  0*mm, side);

  G4LogicalVolume* Tank1_log
    = new G4LogicalVolume(Tank1_box, SiO2, "TankL1"); //,0,0,0);


  G4VPhysicalVolume* Tank1_phys
    = new G4PVPlacement(0, posTank, Tank1_log, "Tank1",
                        expHall_log, false, 0, checkOverlaps);




  // G4Box* Tank2_box = new G4Box("Tank2b", width, top, width);


  G4Trd* Tank2_box = new G4Trd("Tank2b",
             width_bot,
             width_top,
             width_bot,
             width_top,
             top);


  G4ThreeVector posBack = G4ThreeVector(0, top - width_bot, 2*side + width);

  G4ThreeVector posDetector = G4ThreeVector(0, 2*top - width_bot + width_top, 2*side + width);

  G4LogicalVolume* Tank2_log
    = new G4LogicalVolume(Tank2_box, SiO2,"TankL2"); //,0,0,0);
    // = new G4LogicalVolume(Tank2_box,SiO2,"Tank"); //,0,0,0);
    // = new G4LogicalVolume(Tank2_box,Water,"Tank"); //,0,0,0);

  G4RotationMatrix * yRot90deg = new G4RotationMatrix();   // Rotates X and Z axes only

  // yRot90deg->rotateY(M_PI/2.*rad);
  yRot90deg->rotateX(90.*deg);


  G4VPhysicalVolume* Tank_phys
    = new G4PVPlacement(yRot90deg,
                        posBack,
                        Tank2_log,
                        "Tank2",
                        expHall_log,
                        false,
                        0,
                        checkOverlaps);

  // Detector box

  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_GLASS_LEAD");
  detector_mat->SetMaterialPropertiesTable(SiO2_MPTdet);

  G4Box * DetectorBox =
    new G4Box("SolidShapeD",
      det_width, width_top + det_overlap, det_width);

  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(DetectorBox,         //its solid
                        detector_mat,          //its material
                        "Detector");           //its name



  new G4PVPlacement(0,
                    posDetector,
                    logicDetector,
                    "Detector",
                    expHall_log,
                    false,
                    0,
                    checkOverlaps);
*/













  ////////////////////// Sensitiveness

  G4SDManager* SDManager = G4SDManager::GetSDMpointer();

  LeadSD* MySD = new LeadSD("Box_Sensitive_Detector");
  SDManager->AddNewDetector(MySD);

  AirSD* MySD2 = new AirSD("Air_Sensitive_Detector");
  SDManager->AddNewDetector(MySD2);

  logicDetector->SetSensitiveDetector(MySD);
  expHall_log->SetSensitiveDetector(MySD2);

  ////////////////////////////////////
















































  /* -------------------------------------------------------------------------------------------
  // Optical surface
  // ------------------------------------------------------------------------------------------- */



  // surface

  G4OpticalSurface* CrysOpSurface = new G4OpticalSurface("CrystalSurface");

  // G4LogicalBorderSurface* SurfaceCrys[64];

  // for(i=0; i<64; i++){
    // SurfaceCrys[i] = new G4LogicalBorderSurface("CrysSurface",crystal_phys,contact_phys[i],CrysOpSurface);
  // }


  new G4LogicalBorderSurface("CrysSurface1",
                             Tank_phys,
                             // TankV_phys,
                             expHall_phys,
                             CrysOpSurface);


/*
  new G4LogicalBorderSurface("CrysSurfaceV",
                             // Tank_phys,
                             TankV_phys,
                             expHall_phys,
                             CrysOpSurface);
*/


  // new G4LogicalBorderSurface("CrysSurface2",
                             // Tank2_phys,
                             // expHall_phys,
                             // CrysOpSurface);




  // G4double crys_sigma_alpha = 1;
  // G4double crys_sigma_alpha = 0.000115;
  // G4double crys_sigma_alpha = 0.00021;
  G4double crys_sigma_alpha = 0.000115;
  // CrysOpSurface -> SetSigmaAlpha(crys_sigma_alpha);

  // CrysOpSurface -> SetType(dielectric_dielectric);
  CrysOpSurface -> SetType(dielectric_metal);
  CrysOpSurface -> SetModel(unified);
  // CrysOpSurface -> SetFinish(ground);
  CrysOpSurface -> SetFinish(polished);
  // CrysOpSurface -> SetFinish(polishedfrontpainted);
  // CrysOpSurface -> SetFinish(polishedbackpainted);


  // G4double specularlobecrys[NUMcrys] = {0.05, 0.05};
  // G4double specularspikecrys[NUMcrys] = {0.85, 0.85};
  // G4double specularlobecrys[NUMcrys] = {0.01, 0.01};
  // G4double specularspikecrys[NUMcrys] = {0.89, 0.89};

  // G4double backscattercrys[NUMcrys] = {0.05, 0.05};


  /*
  const G4int QuartzRefractionEntries = 101;
  G4double QuartzPhotonEnergy[QuartzRefractionEntries] =
  { 5.90476190476191*eV, 5.06122448979592*eV, 4.42857142857143*eV, 3.93650793650794*eV, 3.54285714285714*eV,
    3.22077922077922*eV, 2.95238095238095*eV, 2.72527472527473*eV, 2.53061224489796*eV, 2.36190476190476*eV,
    2.21428571428571*eV, 2.08403361344538*eV, 1.96825396825397*eV, 1.86466165413534*eV, 1.77142857142857*eV,
    1.68707482993197*eV, 1.61038961038961*eV, 1.54037267080745*eV, 1.47619047619048*eV, 1.41714285714286*eV,
    1.36263736263736*eV, 1.31216931216931*eV, 1.26530612244898*eV, 1.22167487684729*eV, 1.18095238095238*eV,
    1.14285714285714*eV, 1.10714285714286*eV, 1.07359307359307*eV, 1.04201680672269*eV, 1.01224489795918*eV,
    0.984126984126984*eV, 0.957528957528958*eV, 0.932330827067669*eV, 0.908424908424908*eV, 0.885714285714286*eV,
    0.86411149825784*eV, 0.843537414965986*eV, 0.823920265780731*eV, 0.805194805194805*eV, 0.787301587301587*eV,
    0.770186335403727*eV, 0.753799392097264*eV, 0.738095238095238*eV, 0.723032069970845*eV, 0.708571428571429*eV,
    0.694677871148459*eV, 0.681318681318681*eV, 0.668463611859838*eV, 0.656084656084656*eV, 0.644155844155844*eV,
    0.63265306122449*eV, 0.621553884711779*eV, 0.610837438423645*eV, 0.600484261501211*eV, 0.59047619047619*eV,
    0.580796252927401*eV, 0.571428571428571*eV, 0.562358276643991*eV, 0.553571428571428*eV, 0.545054945054945*eV,
    0.536796536796537*eV, 0.528784648187633*eV, 0.521008403361345*eV, 0.513457556935818*eV, 0.506122448979592*eV,
    0.498993963782696*eV, 0.492063492063492*eV, 0.4853228962818*eV, 0.478764478764479*eV, 0.472380952380952*eV,
    0.466165413533835*eV, 0.460111317254174*eV, 0.454212454212454*eV, 0.448462929475588*eV, 0.442857142857143*eV,
    0.437389770723104*eV, 0.43205574912892*eV, 0.426850258175559*eV, 0.421768707482993*eV, 0.416806722689076*eV,
    0.411960132890365*eV, 0.407224958949097*eV, 0.402597402597403*eV, 0.398073836276083*eV, 0.393650793650794*eV,
    0.389324960753532*eV, 0.385093167701863*eV, 0.380952380952381*eV, 0.376899696048632*eV, 0.372932330827068*eV,
    0.369047619047619*eV, 0.365243004418262*eV, 0.361516034985423*eV, 0.357864357864358*eV, 0.354285714285714*eV,
    0.350777934936351*eV, 0.34733893557423*eV, 0.343966712898752*eV, 0.340659340659341*eV, 0.337414965986395*eV,
    0.334231805929919*eV };
    */
  


  G4double specularspikecrys[QuartzRefractionEntries] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0 };





  G4double specularlobecrys[QuartzRefractionEntries] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1 };


/*
  G4double backscattercrys[QuartzRefractionEntries] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    // 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05 };
*/


/*
  G4double reflectivitycrys[QuartzRefractionEntries] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8 };
*/


  G4MaterialPropertiesTable *CrysOpSurfaceProperty = new G4MaterialPropertiesTable();


  // CrysOpSurfaceProperty -> AddProperty("SPECULARLOBECONSTANT",ppcrys,specularlobecrys,NUMcrys);
  // CrysOpSurfaceProperty -> AddProperty("SPECULARSPIKECONSTANT",ppcrys,specularspikecrys,NUMcrys);



  const G4int NUM = 2;
  G4double pp[NUM] = {2.038*eV, 4.144*eV};
  // G4double rindex[NUM] = {1.35, 1.40};
  G4double rindex[NUM] = {1, 1};
  // CrysOpSurfaceProperty -> AddProperty("RINDEX", pp, rindex, NUM);


  // CrysOpSurfaceProperty->AddProperty("RINDEX",       QuartzPhotonEnergy, QuartzRefractiveIndex, QuartzRefractionEntries)
       // ->SetSpline(true);
  // CrysOpSurfaceProperty->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);  

  // CrysOpSurfaceProperty -> AddProperty("BACKSCATTERCONSTANT",ppcrys,backscattercrys,NUMcrys);


  // CrysOpSurfaceProperty -> AddProperty("BACKSCATTERCONSTANT",QuartzPhotonEnergy,backscattercrys,QuartzRefractionEntries);

  // CrysOpSurfaceProperty -> AddProperty("REFLECTIVITY",QuartzPhotonEnergy,reflectivitycrys,QuartzRefractionEntries);





  const G4int NUMcrys = 2;
  G4double ppcrys[NUMcrys] = {1.07*eV, 10.28*eV};


  G4double difpart = 0.25;
  G4double specspikeprob[NUM] = {1 - difpart, 1 - difpart};
  // G4double specspikeprob[NUM] = {1, 1};

  // CrysOpSurfaceProperty -> AddProperty("SPECULARLOBECONSTANT",QuartzPhotonEnergy,specularlobecrys,QuartzRefractionEntries);
  // CrysOpSurfaceProperty -> AddProperty("SPECULARSPIKECONSTANT",QuartzPhotonEnergy,specularspikecrys,QuartzRefractionEntries);
  // CrysOpSurfaceProperty -> AddProperty("SPECULARSPIKECONSTANT", ppcrys, specspikeprob, NUMcrys);
  G4double reflectivitycrys[NUMcrys] = {0.999, 0.999};
  // G4double reflectivitycrys[NUMcrys] = {1, 1};
  CrysOpSurfaceProperty -> AddProperty("REFLECTIVITY",ppcrys,reflectivitycrys,NUMcrys);
  // G4double diffuseprob[NUM] = {1, 1};
  G4double diffuseprob[NUM] = {difpart, difpart};
  // CrysOpSurfaceProperty -> AddProperty("DIFFUSELOBECONSTANT", ppcrys, diffuseprob, NUMcrys);

  CrysOpSurface -> SetMaterialPropertiesTable(CrysOpSurfaceProperty);
















//always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
