/// \file B4cDetectorConstruction.cc
/// \brief Implementation of the B4cDetectorConstruction class

#include "B4cDetectorConstruction.hh"
#include "B4cCalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4ThreadLocal
//G4GlobalMagFieldMessenger* B4cDetectorConstruction::fMagFieldMessenger = 0;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::B4cDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),fScoringVolume(0),fScoringVolume1a(0), 
   fScoringVolume1b(0), fScoringVolume2(0), fScoringVolume3(0), fScoringVolume4(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::~B4cDetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Si");
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  nistManager->FindOrBuildMaterial("G4_Galactic");
  nistManager->FindOrBuildMaterial("G4_B");
  nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE"); 
  nistManager->FindOrBuildMaterial("G4_Gd"); 
  
  
  
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  
  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  
  G4double calorSizeXY  = 1.69*mm;

  
    G4double calorThickness = 300*um;
    G4double L6Thickness = 1*um; 
    G4double L4Thickness = 1*um; 
    G4double L3Thickness = 0.2*um; 
    G4double L2Thickness = 400*nm;
    G4double L1Thickness = 0.5*um;
    G4double bufferThickness = 0.0001*um;
    G4double caseThickness = 2*mm; 
    G4double GdThickness = 2*mm; 
   // auto calorThickness = 2*cm;
    auto worldSizeXY = 1*m;
    auto worldSizeZ  = 1*m;
  
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("G4_AIR");
  //auto defaultMaterial = G4Material::GetMaterial("G4_Galactic");
  
  
  auto bufferMaterial = //G4Material::GetMaterial("G4_Si");
  G4Material::GetMaterial("G4_B");
  //G4Material::GetMaterial("G4_POLYETHYLENE");
  
  
    
    //G4UnitDefinition::BuildUnitsTable();
    G4int ncomponents, natoms;
    G4String name, symbol;
    G4double z, a,  density, abundance;
    G4int n, ncomp;
    
    
    a = 69.723*g/mole;
    G4Element* elGa  = new G4Element(name="Gallium",symbol="Ga" , z= 31., a);

    a = 14.01*g/mole;
    G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);
    
    density = 6.15*g/cm3;
    
    G4Material* GaN = new G4Material(name="GaN", density, ncomponents=2);
    GaN->AddElement(elGa, natoms=1);
    GaN->AddElement(elN, natoms=1);
  
    
  // Defining isotopes as  material
    //Boron-11 
    a = 11.009*g/mole; 
    density = 2.370*g/cm3;
    G4Isotope* isoB11 = new G4Isotope(name = "Boron11", z= 5., n = 11, a); 
    G4Element* elB11 = new G4Element(name = "Boron11", symbol = "B11",  ncomp = 1); 
    elB11->AddIsotope(isoB11, abundance = 100.0*perCent); 
    G4Material* B11 = new G4Material(name="B11", density, ncomponents=1);
    B11->AddElement(elB11, natoms=1);
    
    
    //Boron-10 
    a = 10.01*g/mole; 
    density = 2.13*g/cm3;
    G4Isotope* isoB10 = new G4Isotope(name = "Boron10", z= 5., n = 10, a); 
    G4Element* elB10 = new G4Element(name = "Boron10", symbol = "B10",  ncomp = 1); 
    elB10->AddIsotope(isoB10, abundance = 100.0*perCent); 
    G4Material* B10 = new G4Material(name="B10", density, ncomponents=1);
    B10->AddElement(elB10, natoms=1);
    
    //Phosphorus 31  
    a = 30.97*g/mole; 
    density = 1.82*g/cm3;
    G4Isotope* isoP31 = new G4Isotope(name = "Phosphorus31", z= 15., n = 31, a); 
    G4Element* elP31 = new G4Element(name = "Phosphorus31", symbol = "P31",  ncomp = 1); 
    elP31->AddIsotope(isoP31, abundance = 100.0*perCent); 
    G4Material* P31 = new G4Material(name="P31", density, ncomponents=1);
    P31->AddElement(elP31, natoms=1);
    
    //Aluminium 26 
    a = 26.98*g/mole; 
    density = 2.7*g/cm3;
    G4Isotope* isoAl27 = new G4Isotope(name = "Aluminium27", z= 13., n = 27, a); 
    G4Element* elAl27 = new G4Element(name = "Aluminium27", symbol = "Al27",  ncomp = 1); 
    elAl27->AddIsotope(isoAl27, abundance = 100.0*perCent); 
    G4Material* Al27 = new G4Material(name="Al27", density, ncomponents=1);
    Al27->AddElement(elAl27, natoms=1);
    
    //Silicon 28 
    
    a = 27.97*g/mole; 
    density = 2.3*g/cm3;
    G4Isotope* isoSi28 = new G4Isotope(name = "Silicon28", z= 14., n = 28, a); 
    G4Element* elSi28 = new G4Element(name = "Silicon28", symbol = "Si28",  ncomp = 1); 
    elSi28->AddIsotope(isoSi28, abundance = 100.0*perCent); 
    G4Material* Si28 = new G4Material(name="Si28", density, ncomponents=1);
    Si28->AddElement(elSi28, natoms=1);
    
    
    //Defining Materials for the different layers 
    
    
    
    //Layer 5
    a = 26.98*g/mole; 
    density = 5.5e-7*g/cm3;
    G4Element* elB = new G4Element(name = "Boron", symbol = "B",  ncomp = 2); 
    elB->AddIsotope(isoB10, abundance = 20.0*perCent); 
    elB->AddIsotope(isoB11, abundance = 80.0*perCent);
    G4Material* B = new G4Material(name="Boron", density, ncomponents=1);
    B->AddElement(elB, natoms = 1); 
    
    
    G4double densityL = 2.30000000001*g/cm3; 
    G4Material* SiB = new G4Material(name="SiliconBoron", densityL, ncomponents=2);
    
    SiB->AddMaterial(B,  3.91e-10*perCent);
    SiB->AddMaterial(Si28,  99.9999999996*perCent); 
    
    //Layer 4
    densityL = (2.3)*g/cm3 + (5.5e-7)*g/cm3; 
    G4double perB11 = (5.5e-5)/(2.3+(5.5e-5))*perCent;
    
    G4Material* L4 = new G4Material(name="L4", densityL, ncomponents=2);
    L4->AddMaterial(B11,  perB11);
    L4->AddMaterial(Si28,  100*perCent-perB11); 
    
    //Layer 3
    densityL = (2.3)*g/cm3 + (2.6e-3)*g/cm3; 
    perB11 = (2.6e-1)/(2.3+(2.6e-3))*perCent;
    
    G4Material* L3 = new G4Material(name="L3", densityL, ncomponents=2);
    L3->AddMaterial(P31,  perB11);
    L3->AddMaterial(Si28,  100*perCent-perB11);
    
    //Layer 6
    
    densityL = (2.3)*g/cm3 + (1.8e-3)*g/cm3; 
    perB11 = (1.8e-1)/(2.3+(1.8e-3))*perCent;
    
    G4Material* L6 = new G4Material(name="L6", densityL, ncomponents=2);
    L6->AddMaterial(B11,  perB11);
    L6->AddMaterial(Si28,  100*(perCent-perB11));
    
   
    
    
     

    G4cout << "Material Table: "<<*(G4Material::GetMaterialTable()) << G4endl;
    
    /*
    //Boron-10
    a = 10.01*g/mole;
    G4Element* elB10  = new G4Element(name="Boron10",symbol="B10" , z= 5., a);
    density = 2.13*g/cm3;
    G4Material* B10 = new G4Material(name="B10", density, ncomponents=1);
    B10->AddElement(elB10, natoms=1);
    */
    
  
  auto detectorMaterial = 
  //G4Material::GetMaterial("G4_Si");
  //G4Material::GetMaterial("B10");
  //G4Material::GetMaterial("B11");
  //  G4Material::GetMaterial("P31"); 
  //G4Material::GetMaterial("Al27");
  G4Material::GetMaterial("Si28");
  //G4Material::GetMaterial("L4"); 
  //G4Material::GetMaterial(""); 
  
  auto L4Material = G4Material::GetMaterial("L4");
  auto L3Material = G4Material::GetMaterial("L3"); 
  auto L2Material = G4Material::GetMaterial("G4_SILICON_DIOXIDE"); 
  auto L1Material = G4Material::GetMaterial("Al27"); 
  auto caseMaterial = G4Material::GetMaterial("Al27"); 
  auto GdMaterial = G4Material::GetMaterial("G4_Gd"); 
  auto L6Material = G4Material::GetMaterial("L6"); 
   
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY, worldSizeXY, worldSizeZ); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
    
    
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    //Case  
    // Conical section shape       
  G4double shapecase_rmina =  0.*mm, shapecase_rmaxa = calorSizeXY;
  G4double shapecase_rminb =  0.*mm, shapecase_rmaxb = calorSizeXY;
  G4double shapecase_hz = caseThickness;
  G4double shapecase_phimin = 0.*deg, shapecase_phimax = 360.*deg;
  G4Cons* caseS =  new G4Cons("case", shapecase_rmina, shapecase_rmaxa, shapecase_rminb, shapecase_rmaxb, shapecase_hz/2,
    shapecase_phimin, shapecase_phimax);
  G4double offset = -1*cm; 
  G4ThreeVector casepos = G4ThreeVector(0, 0, offset + (calorThickness/2 + L4Thickness + L3Thickness + L2Thickness+L1Thickness));  
  // */
                         
  auto caseLV
    = new G4LogicalVolume(
                 caseS,     // its solid
                 caseMaterial,  // its material
                 "caseLV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 casepos,  // at (0,0,0)
                 caseLV,          // its logical volume
                 "case",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //Gadolinium Layer 
    // Conical section shape       
  G4double shapeGd_rmina =  0.*mm, shapeGd_rmaxa = calorSizeXY;
  G4double shapeGd_rminb =  0.*mm, shapeGd_rmaxb = calorSizeXY;
  G4double shapeGd_hz = GdThickness;
  G4double shapeGd_phimin = 0.*deg, shapeGd_phimax = 360.*deg;
  G4Cons* GdS =  new G4Cons("Gd", shapeGd_rmina, shapeGd_rmaxa, shapeGd_rminb, shapeGd_rmaxb, shapeGd_hz/2,
    shapeGd_phimin, shapeGd_phimax);
  offset = -5*mm; 
  G4ThreeVector Gdpos = G4ThreeVector(0, 0, offset + (calorThickness/2 + L4Thickness + L3Thickness + L2Thickness+L1Thickness));  
  // */
                         
  auto GdLV
    = new G4LogicalVolume(
                 GdS,     // its solid
                 GdMaterial,  // its material
                 "GdLV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 Gdpos,  // at (0,0,0)
                 GdLV,          // its logical volume
                 "Gd",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
    
    
    
    // Layer 1a 
  
   // Conical section shape       
  G4double shapeL1_rmina =  0.*mm, shapeL1_rmaxa = calorSizeXY;
  G4double shapeL1_rminb =  0.*mm, shapeL1_rmaxb = calorSizeXY;
  G4double shapeL1_hz = L1Thickness;
  G4double shapeL1a_phimin = 0.*deg, shapeL1a_phimax = 180.*deg;
  G4Cons* L1aS =  new G4Cons("L1a", shapeL1_rmina, shapeL1_rmaxa, shapeL1_rminb, shapeL1_rmaxb, shapeL1_hz/2,
    shapeL1a_phimin, shapeL1a_phimax);
  G4ThreeVector L1apos = G4ThreeVector(0, 0, calorThickness/2 + L4Thickness + L3Thickness + L2Thickness+L1Thickness/2);  
  // */
                         
  auto L1aLV
    = new G4LogicalVolume(
                 L1aS,     // its solid
                 L1Material,  // its material
                 "L1aLV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 L1apos,  // at (0,0,0)
                 L1aLV,          // its logical volume
                 "L1a",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  // Layer 1b 
  
   // Conical section shape       
  shapeL1_rmina =  0.*mm, shapeL1_rmaxa = calorSizeXY;
  shapeL1_rminb =  0.*mm, shapeL1_rmaxb = calorSizeXY;
  shapeL1_hz = L1Thickness;
  G4double shapeL1b_phimin = -180.*deg, shapeL1b_phimax = 180.*deg;
  G4Cons* L1bS =  new G4Cons("L1b", shapeL1_rmina, shapeL1_rmaxa, shapeL1_rminb, shapeL1_rmaxb, shapeL1_hz/2,
    shapeL1b_phimin, shapeL1b_phimax);
  G4ThreeVector L1bpos = G4ThreeVector(0, 0, calorThickness/2 + L4Thickness + L3Thickness + L1Thickness/2);  
  // */
                         
  auto L1bLV
    = new G4LogicalVolume(
                 L1bS,     // its solid
                 L1Material,  // its material
                 "L1bLV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 L1bpos,  // at (0,0,0)
                 L1bLV,          // its logical volume
                 "L1b",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
    
    // Layer 2
  
   // Conical section shape       
  G4double shapeL2_rmina =  0.*mm, shapeL2_rmaxa = calorSizeXY;
  G4double shapeL2_rminb =  0.*mm, shapeL2_rmaxb = calorSizeXY;
  G4double shapeL2_hz = L2Thickness;
  G4double shapeL2_phimin = 0.*deg, shapeL2_phimax = 180.*deg;
  G4Cons* L2S =  new G4Cons("L2", shapeL2_rmina, shapeL2_rmaxa, shapeL2_rminb, shapeL2_rmaxb, shapeL2_hz/2,
    shapeL2_phimin, shapeL2_phimax);
  G4ThreeVector L2pos = G4ThreeVector(0, 0, calorThickness/2 + L4Thickness + L3Thickness + L2Thickness/2);  
  // */
                         
  auto L2LV
    = new G4LogicalVolume(
                 L2S,     // its solid
                 L2Material,  // its material
                 "L2LV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 L2pos,  // at (0,0,0)
                 L2LV,          // its logical volume
                 "L2",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  // Layer 3
  
   // Conical section shape       
  G4double shapeL3_rmina =  0.*mm, shapeL3_rmaxa = calorSizeXY;
  G4double shapeL3_rminb =  0.*mm, shapeL3_rmaxb = calorSizeXY;
  G4double shapeL3_hz = L3Thickness;
  G4double shapeL3_phimin = 0.*deg, shapeL3_phimax = 360.*deg;
  G4Cons* L3S =  new G4Cons("L3", shapeL3_rmina, shapeL3_rmaxa, shapeL3_rminb, shapeL3_rmaxb, shapeL3_hz/2,
    shapeL3_phimin, shapeL3_phimax);
  G4ThreeVector L3pos = G4ThreeVector(0, 0, calorThickness/2 + L4Thickness + L3Thickness/2);  
  // */
                         
  auto L3LV
    = new G4LogicalVolume(
                 L3S,     // its solid
                 L3Material,  // its material
                 "L3LV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 L3pos,  // at (0,0,0)
                 L3LV,          // its logical volume
                 "L3",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
                                 
  // L4
  
   // Conical section shape       
  G4double shapeL4_rmina =  0.*mm, shapeL4_rmaxa = calorSizeXY;
  G4double shapeL4_rminb =  0.*mm, shapeL4_rmaxb = calorSizeXY;
  G4double shapeL4_hz = L4Thickness;
  G4double shapeL4_phimin = 0.*deg, shapeL4_phimax = 360.*deg;
  G4Cons* L4S =  new G4Cons("L4", shapeL4_rmina, shapeL4_rmaxa, shapeL4_rminb, shapeL4_rmaxb, shapeL4_hz/2,
    shapeL4_phimin, shapeL4_phimax);
  G4ThreeVector L4pos = G4ThreeVector(0, 0, calorThickness/2 + L4Thickness/2);  
  // */
                         
  auto L4LV
    = new G4LogicalVolume(
                 L4S,     // its solid
                 L4Material,  // its material
                 "L4LV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 L4pos,  // at (0,0,0)
                 L4LV,          // its logical volume
                 "L4",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  
  // Layer 5
  
  // Conical section shape       
  G4double shape1_rmina =  0.*mm, shape1_rmaxa = calorSizeXY;
  G4double shape1_rminb =  0.*mm, shape1_rmaxb = calorSizeXY;
  G4double shape1_hz = calorThickness;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Cons* calorimeterS =  new G4Cons("Calorimeter", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz/2,
    shape1_phimin, shape1_phimax);
  // */
                         
  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 detectorMaterial,  // its material
                 "CalorLV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  
  
  // L6
  
   // Conical section shape       
  G4double shapeL6_rmina =  0.*mm, shapeL6_rmaxa = calorSizeXY;
  G4double shapeL6_rminb =  0.*mm, shapeL6_rmaxb = calorSizeXY;
  G4double shapeL6_hz = L6Thickness;
  G4double shapeL6_phimin = 0.*deg, shapeL6_phimax = 360.*deg;
  G4Cons* L6S =  new G4Cons("L6", shapeL6_rmina, shapeL6_rmaxa, shapeL6_rminb, shapeL6_rmaxb, shapeL6_hz/2,
    shapeL6_phimin, shapeL6_phimax);
  G4ThreeVector L6pos = G4ThreeVector(0, 0, -calorThickness/2 - L6Thickness/2);  
  // */
                         
  auto L6LV
    = new G4LogicalVolume(
                 L6S,     // its solid
                 L6Material,  // its material
                 "L6LV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 L6pos,  // at (0,0,0)
                 L6LV,          // its logical volume
                 "L6",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  /* 
  // Measurement box 

  auto rulerS
    = new G4Box("Ruler",     // its name
                 calorSizeXY, calorSizeXY, calorThickness); // its size

                         
  auto rulerLV
    = new G4LogicalVolume(
                 rulerS,     // its solid
                 defaultMaterial,  // its material
                 "rulerLV");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 rulerLV,          // its logical volume
                 "Ruler",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
   
  */
  
  
    // silicon box
    //
     DetectorPos = bufferThickness/2+calorThickness/2; 
     G4ThreeVector positionTarget = G4ThreeVector(0,0, DetectorPos);
    /*
    auto siliconBOXS
    = new G4Box("siliconBOX",     // its name
                calorSizeXY, calorSizeXY, bufferThickness); // its size
    */
    
    /*
   // Conical section shape       
  G4double shape2_rmina =  0.*mm, shape2_rmaxa = calorSizeXY;
  G4double shape2_rminb =  0.*mm, shape2_rmaxb = calorSizeXY;
  G4double shape2_hz = bufferThickness;
  G4double shape2_phimin = 0.*deg, shape2_phimax = 360.*deg;
  G4Cons* siliconBOXS =  new G4Cons("siliconBOX", shape2_rmina, shape2_rmaxa, shape2_rminb, shape2_rmaxb, shape2_hz/2,
    shape2_phimin, shape2_phimax);
   
    
    auto siliconBOXLV
    = new G4LogicalVolume(
                          siliconBOXS,     // its solid
                          bufferMaterial,  // its material
                          "siliconBOXLV");   // its name
    
    new G4PVPlacement(
                      0,                // no rotation
                      positionTarget,  // at (0,0,0)
                      siliconBOXLV,          // its logical volume
                      "siliconBOX",    // its name
                      worldLV,          // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      fCheckOverlaps);  // checking overlaps
*/
  
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
    
    //auto simpleBoxVisAtttmp= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    //worldLV->SetVisAttributes(simpleBoxVisAtttmp);

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

    auto L1aVisAtt2= new G4VisAttributes(G4Colour(1.0,0.0,1.0));
    L1aLV->SetVisAttributes(L1aVisAtt2);
    auto L1bVisAtt2= new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    L1bLV->SetVisAttributes(L1bVisAtt2); 
    auto GdVisAtt2= new G4VisAttributes(G4Colour(1.0,0.5,1.0));
    GdLV->SetVisAttributes(GdVisAtt2);
    auto L6VisAtt2= new G4VisAttributes(G4Colour(0,126,0));
    L6LV->SetVisAttributes(L6VisAtt2);
    
    
  //
  // Always return the physical World
  //
    //G4cout<<fScoringVolume->GetName()<<G4endl; 
    fScoringVolume4 = L4LV;
    fScoringVolume3 = L3LV; 
    fScoringVolume2 = L2LV; 
    fScoringVolume1a = L1aLV; 
    fScoringVolume1b = L1bLV; 
    fScoringVolume = calorLV; 
    
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //
  auto CalorSD
    = new B4cCalorimeterSD("CalorSD", "CalorHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(CalorSD);
  SetSensitiveDetector("CalorLV",CalorSD);
  SetSensitiveDetector("L4LV", CalorSD); 
  SetSensitiveDetector("L3LV", CalorSD); 
  SetSensitiveDetector("L2LV", CalorSD); 
  SetSensitiveDetector("L1aLV", CalorSD); 
  SetSensitiveDetector("L1bLV", CalorSD); 
  

 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
