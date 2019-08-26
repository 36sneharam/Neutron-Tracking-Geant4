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
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class

#include "B4RunAction.hh"
#include "B4Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::B4RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //
  
  // Creating histograms (MeV)
  analysisManager->CreateH1("Edep","Edep in Detector",1000, 0., (300*10000)*eV, "eV");
  analysisManager->CreateH1("Amass","Atomic Mass", 100, 0., 100);
  analysisManager->CreateH1("LTrack","trackL in detector", 100, 0., 5*mm);
  analysisManager->CreateH1("Ions","Number of ions per event", 100, 0., 10);
  analysisManager->CreateH1("evtID", "event ID", 1000, 0, 100000000);
  analysisManager->CreateH1("IonTrack","Ion track length", 100, 0., 5*um);
  analysisManager->CreateH1("Anum","Atomic Number", 100, 0., 100);
  analysisManager->CreateH1("EdepStep","Edep in Detector for ions", 1000, 0., 2500*eV);
  analysisManager->CreateH1("Edep10bin", " Edep in Detector for ions", 16, 0., 80*MeV);
  analysisManager->CreateH1("SP","SP for events with only 1 ion",100, 0, 30);
  analysisManager->CreateH1("IonTrack1ion","Ion track length with 1 ion only", 100, 0., 10);
  analysisManager->CreateH1("Interact_Num", "Interaction Number", 10, 10., 20);
  analysisManager->CreateH1("ElasvsInelas", "Elastic (1) vs Inelastic (2)", 2, 1, 3); 
  analysisManager->CreateH1("ElasticEdep", "Elastic Edep in Detector", 1000, 0., (300*10000)*eV,"eV" ); 
  analysisManager->CreateH1("InelasEdep", "Inelastic Edep in Detector", 1000, 0., (300*10000)*eV,"eV"); 
  
  analysisManager->CreateH1("Qdep","Qdep in Detector",1000, 0., (300)*10000);
  analysisManager->CreateH1("ElasticQdep", "Elastic Qdep in Detector", 1000, 0., (300)*10000); 
  analysisManager->CreateH1("InelasQdep", "Inelastic Qdep in Detector", 1000, 0., (300)*10000); 
  
  analysisManager->CreateH1("SecondQdep","Qdep in Detector by Secondary particles",1000, 0., (300)*10000);
  analysisManager->CreateH1("SecondElasticQdep", "Elastic Qdep in Detector by Secondary particles", 1000, 0., (300)*10000); 
  analysisManager->CreateH1("SecondInelasQdep", "Inelastic Qdep in Detector by Secondary particles", 1000, 0., (300)*10000); 

  
  analysisManager->CreateH2("TrackMass","Ion track length vs Amass", 100, 0., 100, 100, 0., 10); // first H2
  analysisManager->CreateH2("TrackAnum","Ion track length vs Anum", 100, 0., 100, 100, 0., 10);
  analysisManager->CreateH2("SPAnum","Stopping power vs Anum", 100, 0., 100, 100, 0., 50);
  analysisManager->CreateH2("EDEPvsZ","Edep vs Anum", 100, 0., 100, 100, 0., 80);
  
  
  
  //analysisManager->CreateH1("electrons", "Number of electrons per event", 100, 0., 10); 
  
  //analysisManager->CreateH1("gamma", "Number of gamma particles per event", 100, 0., 10);
  
  //analysisManager->CreateH1("All electrons", "Total number of electrons produced per event", 100, 0, 400); 
  

    
 //analysisManager->CreateH1("EdepStep","EdepStep", 1000, 0., 0.5*MeV);
    // Creating ntuple
  //
  analysisManager->CreateNtuple("B4", "Edep");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("Amass");
  analysisManager->CreateNtupleDColumn("LTrack");
  analysisManager->CreateNtupleDColumn("Ions");
  analysisManager->CreateNtupleDColumn("evtID");
  analysisManager->CreateNtupleDColumn("IonTrack");
  analysisManager->CreateNtupleDColumn("Anum");
  analysisManager->CreateNtupleDColumn("isEventwithIon");
  analysisManager->CreateNtupleDColumn("electrons"); 
  analysisManager->CreateNtupleDColumn("gamma");
    
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::~B4RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "B4-Calor";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
    
    G4cout << " ECalor : mean = "
       << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
       << " rms = " 
	   << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
    G4cout << " LCalor : mean = "
    << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length")
    << " rms = "
    << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;
    
    
  

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
  
  
  
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
