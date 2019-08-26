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
/// \file B4cEventAction.cc
/// \brief Implementation of the B4cEventAction class

#include "B4cEventAction.hh"
#include "B4cCalorimeterSD.hh"
#include "B4cCalorHit.hh"
#include "B4Analysis.hh"
#include "SteppingAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::B4cEventAction()
 : G4UserEventAction(),
   fCalorHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::~B4cEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorHitsCollection* 
B4cEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<B4cCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("B4cEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::PrintEventStatistics(
                              G4double CalorEdep, G4double CalorTrackLength) const
{
  // print event statistics
  G4cout
     << "   Calor: total energy: "
     << std::setw(7) << G4BestUnit(CalorEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(CalorTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    fTotalEnergyDepositStep = 0.;
    fcount_ions = 0;
    fAtomicMass=0;
    fiSIon=false;
    fisElectron = false; 
    fIonLengthStep=0;
    fAtomicNumber=0;
    fcount_electrons = 0; 
    fcount_elastic = 0; 
    fcount_inelastic = 0; 
    forigin_ID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::EndOfEventAction(const G4Event* event)
{  
  // Get hits collections IDs (only once)
  if ( fCalorHCID == -1 ) {
    fCalorHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("CalorHitsCollection");
    
  }

  // Get hits collections
  auto CalorHC = GetHitsCollection(fCalorHCID, event);

  // Get hit with total values
  auto CalorHit = (*CalorHC)[CalorHC->entries()-1];
    
  
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    PrintEventStatistics(CalorHit->GetEdep(),CalorHit->GetTrackLength());
  }  

  // Fill histograms, ntuple
  //

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
 
  // fill histograms
    
    
/*0*/   if (CalorHit->GetEdep()>0) analysisManager->FillH1(0, CalorHit->GetEdep());
        if (CalorHit->GetEdep()>0) {
            analysisManager->FillH1(15, (CalorHit->GetEdep())/(3.6*0.000001));
            if (forigin_ID != 0) {
                analysisManager->FillH1(18, (CalorHit->GetEdep())/(3.6*0.000001));
            }
        }
/*1*/   // mass atomico A to be filled in stepping action
/*2*/   if (CalorHit->GetTrackLength()>0) analysisManager->FillH1(2, CalorHit->GetTrackLength());
/*3*/   if (fcount_ions>0)  analysisManager->FillH1(3, fcount_ions);
/*4*/   analysisManager->FillH1(4,eventID);
/*5*/   // ion track
/*6*/   // numbero atomico Z to be filled in stepping action
/*7*/   if(fTotalEnergyDepositStep>0) analysisManager->FillH1(7,fTotalEnergyDepositStep);
/*8*/   if(fTotalEnergyDepositStep>0) analysisManager->FillH1(8,fTotalEnergyDepositStep); //10 bin
   
    if (fcount_ions==1 && fIonLengthStep>0 ) {
        G4double SP = (fTotalEnergyDepositStep/(fIonLengthStep/10))/(2.32*1000);
        G4cout << "Edep for the ONLY ion " << fTotalEnergyDepositStep << " Track Length for the ONLY ion " << fIonLengthStep << G4endl;
        G4cout << " Mi stampo anche il ratio " << SP << " MeV/cm**2/mg " <<  G4endl;
       
        /*9*/ analysisManager->FillH1(9,SP); // MeV/cm
        /*10*/ analysisManager->FillH1(10,fIonLengthStep*1000);
        G4cout << " e la massa  " << fAtomicMass <<  G4endl;
        /*0-2*/ analysisManager->FillH2(0, fAtomicMass, fIonLengthStep*1000);
        /*1-2*/ analysisManager->FillH2(1, fAtomicNumber, fIonLengthStep*1000);
        /*2-2*/ analysisManager->FillH2(2, fAtomicNumber, SP);
         /*3-2*/ analysisManager->FillH2(3, fAtomicNumber, fTotalEnergyDepositStep);
         
                 
    };
    
    
    
    
    //analysisManager->FillH1(13, fcount_allelectrons); 
    //analysisManager->FillH1(11, fcount_electrons); 
    //analysisManager->FillH1(12, fcount_gamma);
    
    
    

    //fElasticity = GetElasticity(); 
    G4cout<<"Elastic: "<< fcount_elastic<<G4endl; 
    G4cout<<"Inelastic: "<<fcount_inelastic<<G4endl;
    
    if (fcount_elastic != 0) {
        analysisManager->FillH1(13, CalorHit->GetEdep());
        analysisManager->FillH1(16, (CalorHit->GetEdep())/(3.6*0.000001));
        if (forigin_ID != 0) {
            analysisManager->FillH1(18, (CalorHit->GetEdep())/(3.6*0.000001));
            
        }
    }
    
    if (fcount_inelastic != 0) {
        analysisManager->FillH1(14, CalorHit->GetEdep());
        analysisManager->FillH1(17, (CalorHit->GetEdep())/(3.6*0.000001));
        if (forigin_ID != 0) {
            analysisManager->FillH1(18, (CalorHit->GetEdep())/(3.6*0.000001));
        }
    }
    
    
   
    G4cout<<"Origin ID: "<<forigin_ID<<G4endl; 
    

    //G4cout << " Energy for " << fTotalEnergyDepositStep << G4endl;
  //analysisManager->FillH1(1,fAtomicMass);
  //analysisManager->FillH1(4, fTotalEnergyDepositStep);
    // fill ntuple
  
/*0*/   analysisManager->FillNtupleDColumn(0, CalorHit->GetEdep());
/*1*/   //massa atomica to be filled in stepping action
/*2*/   analysisManager->FillNtupleDColumn(2, CalorHit->GetTrackLength());
/*3*/   analysisManager->FillNtupleDColumn(3, fcount_ions);
/*4*/   analysisManager->FillNtupleDColumn(4, eventID );
        //analysisManager->FillNtupleDColumn(14, fcount_electrons); 
        //analysisManager->FillNtupleDColumn(15, fcount_gamma); 
/*5*/   //ion track
/*6*/   // numbero atomico Z to be filled in stepping action
/*7*/   analysisManager->FillNtupleDColumn(7, fiSIon); // Edep if isIon == 1 mi da l'energia dell'evento che ha lo ion
/*8*/   analysisManager->FillNtupleDColumn(8, fisElectron);

    
    
  // if (fcount_ions >0) G4cout << " Atomic number : " << fAtomicMass <<G4endl;
  //analysisManager->FillNtupleDColumn(1, fAtomicMass);
   
  analysisManager->AddNtupleRow();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
