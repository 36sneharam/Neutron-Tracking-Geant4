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
/// \file B4cCalorimeterSD.cc
/// \brief Implementation of the B4cCalorimeterSD class

#include "B4cCalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "B4Analysis.hh"
#include "G4UnitsTable.hh"
#include "B4Analysis.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorimeterSD::B4cCalorimeterSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName
                            )
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr)

{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorimeterSD::~B4cCalorimeterSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new B4cCalorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  auto hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  // Create hits
  // fNofCells for cells + one more for total sums 
  
    fHitsCollection->insert(new B4cCalorHit());
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B4cCalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{  
  // energy deposit
  //auto edep = step->GetTotalEnergyDeposit();
  //if(edep>4)  G4cout<< " L'energia e' maggiore di 4 MeV per " << edep << " " << step->GetTrack()->GetDefinition()->GetParticleName() << " is secondary ? " << step->GetTrack()->GetParentID()<< G4endl;
    
    // energy deposit
    auto edep = step->GetTotalEnergyDeposit();
    
    // step length
    G4double stepLength = 0.;
    if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
        stepLength = step->GetStepLength();
    }

 // G4String processName =step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  //G4cout << " Il Processo da cui vieni:  " << processName << G4endl;     
  
 if ( edep==0. ) return false;
 if ( stepLength==0. ) return false;
  
 //if step->GetTrack()->GetParentID() == 0
    
 // B4cCalorHit* hit = new B4cCalorHit();
    auto hit
    = (*fHitsCollection)[fHitsCollection->entries()-1];
    
    if ( ! hit ) {
        G4ExceptionDescription msg;
        msg << "Cannot access hit " ;
        G4Exception("B4cCalorimeterSD::ProcessHits()",
                    "MyCode0004", FatalException, msg);
    }
  // Add values
  // G4cout<< " L'energia e' maggiore di 4 MeV per " << edep << " " << step->GetTrack()->GetDefinition()->GetParticleName() << " is secondary ? " << step->GetTrack()->GetParentID()<< G4endl;
  
    hit->Add(edep, stepLength);
 // hit->SetEdepIons(Edep_ions);
  
  //fHitsCollection->insert( hit );
    
   // if(edep>4)  G4cout<< " L'energia e' maggiore di 4 MeV per " << hit->GetEdep() << " " << step->GetTrack()->GetDefinition()->GetParticleName() << " is secondary ? " << step->GetTrack()->GetParentID()<< G4endl;
    
    //auto analysisManager = G4AnalysisManager::Instance();
    
    // fill histograms
    //analysisManager->FillH1(0, edep);
   // G4double TLength = 0;
    /*
    
    G4ParticleDefinition *particle;
  const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();
    for (size_t lp=0; lp<(*secondary).size(); lp++) {
        particle = (*secondary)[lp]->GetDefinition();
       if (particle->GetParticleType() == "nucleus") {
           
           G4cout << " Ion track is " << stepLength << G4endl;
           
       }
    }*/
    /* //good implementation
    auto analysisManager = G4AnalysisManager::Instance();
    G4int parent_ID = step->GetTrack()->GetParentID();
    G4String particleType = step->GetTrack()->GetParticleDefinition()->GetParticleType();
    if (parent_ID>0 && particleType == "nucleus" ) {
        G4double damntrack = step->GetTrack()->GetTrackLength();
        G4cout << " damn track " << G4BestUnit(damntrack,"Length") << " " << particleType << G4endl;
        analysisManager->FillNtupleDColumn(5, damntrack);
        if (damntrack>0)  analysisManager->FillH1(4, damntrack);
    }*/
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  auto nofHits = fHitsCollection->entries();
  if ( verboseLevel>1 ) { 
    
    G4cout
       << G4endl 
       << "-------->Hits Collection: in this event they are " << nofHits 
       << " hits in the tracker chambers: " << G4endl;
    for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
    
  }
  //  for ( G4int i=0; i<nofHits; i++ ) G4cout << " Il Processo da cui vieni:  " <<  (*fHitsCollection)[i]->GetPN() << G4endl;     

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
