#include "B4PrimaryGeneratorAction.hh"
#include "B4cDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr)
{

  G4int nofParticles = 1;

  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particleDefinition);
//  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  
  
  
  
  
  fParticleGun->SetParticleEnergy(14*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0.;
  G4double worldXHalflength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength(); 
    worldXHalflength = worldBox->GetXHalfLength();
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
  
  // Set gun position
  
  G4double zpos =  1.5*cm;
  
  //fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));
   fParticleGun->SetParticlePosition(G4ThreeVector(0., 0. , -zpos));
  // fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 62*um));
  // fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  
  
  //Randomize Particle position 
  /*
  G4double envSizeXY = worldXHalflength*2;
  G4double envSizeZ = worldZHalfLength*2;    
  G4double size = 0.8; 
  G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double z0 = -0.5 * envSizeZ;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  */
  
  // Randomize Particle momentum direction 
  /*
  //Sphere coordinates
    G4double twopi = 2*3.1415;
    G4double cosTheta = -1.0+2.0*G4UniformRand();
    G4double phi = twopi * G4UniformRand();
    G4double sinTheta = sqrt(1-cosTheta*cosTheta);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sinTheta * cos(phi) , sinTheta * sin(phi), cosTheta));
  */
  
  
  
  //Cone coordinates momentum 
    G4double theta = 2*3.141592*G4UniformRand();
    G4double cosTheta = cos(theta);
    G4double sinTheta = sin(theta);
    G4double R = (1.69*mm)*G4UniformRand();
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(R*cosTheta , R*sinTheta , zpos));
  

  
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0 , 0 , zpos));
  
  
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

