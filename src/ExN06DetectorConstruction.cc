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



#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "LeadSD.hh"



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

// ------------- Materials -------------

  G4double a, z, density;
  G4int nelements;

// Air
// 
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

// Water
// 
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* Water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  Water->AddElement(H, 2);
  Water->AddElement(O, 1);


  G4String name, symbol;
  G4int ncomponents, natoms;

  // Quartz with refraction
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


  //
  // ------------ Generate & Add Material Properties Table ------------
  //


  // Quartz refraction from:
  // I. H. Malitson. Interspecimen Comparison of the Refractive Index of Fused Silica
  // at http://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson

  const G4int QuartzRefractionEntries = 10;
  G4double QuartzPhotonEnergy[QuartzRefractionEntries] =
  { 5.90476190476191*eV, 5.06122448979592*eV, 4.42857142857143*eV, 3.93650793650794*eV, 3.54285714285714*eV,
    3.22077922077922*eV, 2.95238095238095*eV, 2.72527472527473*eV, 2.53061224489796*eV, 2.36190476190476*eV};/*,
    2.21428571428571*eV, 2.08403361344538*eV, 1.96825396825397*eV, 1.86466165413534*eV, 1.77142857142857*eV,
    1.68707482993197*eV, 1.61038961038961*eV, 1.54037267080745*eV, 1.47619047619048*eV, 1.41714285714286*eV,
    1.36263736263736*eV, 1.31216931216931*eV, 1.26530612244898*eV, 1.22167487684729*eV, 1.18095238095238*eV,
    1.14285714285714*eV, 1.10714285714286*eV, 1.07359307359307*eV, 1.04201680672269*eV, 1.01224489795918*eV};/*,
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
    0.334231805929919*eV };*/

  G4double QuartzRefractiveIndex[QuartzRefractionEntries] =
  { 1.5383576204905, 1.5102724365895, 1.4941636611188, 1.4839008951423, 1.476891413496,
    1.4718556531995, 1.4680936900401, 1.46519309996, 1.4628966820387, 1.4610366660574};/*,
    1.4594995356592, 1.458206104926, 1.4570996888769, 1.4561387969803, 1.4552924662623,
    1.454537192876, 1.4538548630589, 1.4532313266004, 1.4526553936728, 1.4521181167939,
    1.451612268629, 1.4511319566977, 1.4506723353353, 1.4502293877559, 1.4497997593263,
    1.4493806287126, 1.4489696073537, 1.4485646603469, 1.4481640436751, 1.4477662540207};/*,
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
    1.3992797723176 };*/





  const G4int nEntries = 32;

  G4double PhotonEnergy[nEntries] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
//
// Water
//              
  G4double RefractiveIndex1[nEntries] =
            { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
              1.346,  1.3465, 1.347,  1.3475, 1.348,
              1.3485, 1.3492, 1.35,   1.3505, 1.351,
              1.3518, 1.3522, 1.3530, 1.3535, 1.354,
              1.3545, 1.355,  1.3555, 1.356,  1.3568,
              1.3572, 1.358,  1.3585, 1.359,  1.3595,
              1.36,   1.3608};

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

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries)
       ->SetSpline(true);
  // myMPT1->AddConstProperty("RINDEX", 2); // DOES NOT WORK

  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries)
       ->SetSpline(true);


  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();

  myMPT2->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries)
       ->SetSpline(true);
  // myMPT2->AddConstProperty("RINDEX", 2); // DOES NOT WORK

  myMPT2->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption2,     nEntries)
       ->SetSpline(true);

  
  // myMPT2->AddConstProperty("ABSLENGTH",    0.001*mm);





  // Quartz Properties

  G4MaterialPropertiesTable* quartzMPT = new G4MaterialPropertiesTable();

  quartzMPT->AddProperty("RINDEX",       QuartzPhotonEnergy, QuartzRefractiveIndex, QuartzRefractionEntries)
       ->SetSpline(true);
  // quartzMPT->AddConstProperty("RINDEX", 2); // DOES NOT WORK

  quartzMPT->AddProperty("ABSLENGTH",    QuartzPhotonEnergy, Absorption1,     nEntries)
       ->SetSpline(true);



  // SiO2->SetMaterialPropertiesTable(quartzMPT);
  SiO2->SetMaterialPropertiesTable(myMPT1);
  // SiO2->SetMaterialPropertiesTable(myMPT2);

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

  Water->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the Water scintillator

  //Water->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

//
// Air
//
  G4double RefractiveIndex2[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* airMPT = new G4MaterialPropertiesTable();
  airMPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
  
  Air->SetMaterialPropertiesTable(airMPT);









//
// ------------- Volumes --------------

  G4bool checkOverlaps = true;







// The experimental Hall
//
  //G4Box* expHall_box = new G4Box("World",expHall_x,expHall_y,expHall_z);
  G4Box* expHall_box = new G4Box("World", 10*cm, 10*cm, 10*cm);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box, Air, "World"); //,0,0,0);
    //= new G4LogicalVolume(expHall_box,Water,"World",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0, checkOverlaps);





// The Water Tank
//        
  //G4Box* waterTank_box = new G4Box("Tank",tank_x,tank_y,tank_z);
  G4Box* waterTank_box = new G4Box("Tank", 1*cm, 1*cm, 5*cm);

  G4LogicalVolume* waterTank_log
    = new G4LogicalVolume(waterTank_box, SiO2,"Tank"); //,0,0,0);
    // = new G4LogicalVolume(waterTank_box,SiO2,"Tank"); //,0,0,0);
    // = new G4LogicalVolume(waterTank_box,Water,"Tank"); //,0,0,0);


  G4VPhysicalVolume* waterTank_phys
    = new G4PVPlacement(0,G4ThreeVector(),waterTank_log,"Tank",
                        expHall_log,false,0, checkOverlaps);





  // Detector box

  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_GLASS_LEAD");
  detector_mat->SetMaterialPropertiesTable(myMPT2);


  G4Box * DetectorBox =
    new G4Box("SolidShapeD",
      1*cm, 1*cm, 1*cm);

  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(DetectorBox,         //its solid
                        detector_mat,          //its material
                        "Detector");           //its name



  G4ThreeVector posDetector = G4ThreeVector(0,  0*cm, 6*cm);
  new G4PVPlacement(0,
                    posDetector,
                    logicDetector,
                    "Detector",
                    expHall_log,
                    false,
                    0,
                    checkOverlaps);

  ////////////////////// Sensitiveness

  G4SDManager* SDManager = G4SDManager::GetSDMpointer();

  LeadSD* MySD = new LeadSD("Box_Sensitive_Detector");

  SDManager->AddNewDetector(MySD);
  logicDetector->SetSensitiveDetector(MySD);

  ////////////////////////////////////

/*
// The Air Bubble
//   
  G4Box* bubbleAir_box = new G4Box("Bubble",bubble_x,bubble_y,bubble_z);

  G4LogicalVolume* bubbleAir_log
    = new G4LogicalVolume(bubbleAir_box,Air,"Bubble",0,0,0);

//G4VPhysicalVolume* bubbleAir_phys =
      new G4PVPlacement(0,G4ThreeVector(0,2.5*m,0),bubbleAir_log,"Bubble",
                        waterTank_log,false,0);
*/

// ------------- Surfaces --------------
/*
//
// Water Tank
//
  G4OpticalSurface* OpWaterSurface = new G4OpticalSurface("WaterSurface");
  OpWaterSurface->SetType(dielectric_dielectric);
  OpWaterSurface->SetFinish(ground);
  OpWaterSurface->SetModel(unified);


  new G4LogicalBorderSurface("WaterSurface",
                                 waterTank_phys,expHall_phys,OpWaterSurface);


// Air Bubble
//
  G4OpticalSurface* OpAirSurface = new G4OpticalSurface("AirSurface");
  OpAirSurface->SetType(dielectric_dielectric);
  OpAirSurface->SetFinish(polished);
  OpAirSurface->SetModel(glisur);

  G4LogicalSkinSurface* AirSurface = 
          new G4LogicalSkinSurface("AirSurface", bubbleAir_log, OpAirSurface);

  G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
        (AirSurface->GetSurface(bubbleAir_log)->GetSurfaceProperty());

  if (opticalSurface) opticalSurface->DumpInfo();
*/

//
// Generate & Add Material Properties Table attached to the optical surfaces
//
/*
  const G4int num = 2;
  G4double Ephoton[num] = {2.034*eV, 4.136*eV};

  //OpticalWaterSurface 
  G4double RefractiveIndex[num] = {1.35, 1.40};
  G4double SpecularLobe[num]    = {0.3, 0.3};
  G4double SpecularSpike[num]   = {0.2, 0.2};
  G4double Backscatter[num]     = {0.2, 0.2};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
  
  myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);

  OpWaterSurface->SetMaterialPropertiesTable(myST1);

  //OpticalAirSurface
  G4double Reflectivity[num] = {0.3, 0.5};
  G4double Efficiency[num]   = {0.8, 1.0};

  G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();

  myST2->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);
  myST2->AddProperty("EFFICIENCY",   Ephoton, Efficiency,   num);

  OpAirSurface->SetMaterialPropertiesTable(myST2);
*/

//always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
