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
/// \file B4cEventAction.hh
/// \brief Definition of the B4cEventAction class

#ifndef B4cEventAction_h
#define B4cEventAction_h 1

#include "G4UserEventAction.hh"

#include "B4cCalorHit.hh"

#include "globals.hh"

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in the hits collections.

class B4cEventAction : public G4UserEventAction
{
public:
  B4cEventAction();
  virtual ~B4cEventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);
    
  void AddEdepStep(G4double EdepStep)    {fTotalEnergyDepositStep += EdepStep;};
  G4double GetEnergyDepositStep()    {return fTotalEnergyDepositStep;};
  
  void AddElectronNumber(G4int count_electrons) {fcount_electrons = count_electrons;}; 
  G4int GetElectronNumber () {return fcount_electrons;};
  
  void AddGammaNumber(G4int count_gamma) {fcount_gamma = count_gamma;}; 
  G4int GetGammaNumber () {return fcount_gamma;};
    
  void AddIonsNumber(G4int count_ions) {fcount_ions = count_ions;};
  G4int GetIonsNumber() {return fcount_ions;};
    
    void AddAtomicMass(G4int AtomicMass) {fAtomicMass = AtomicMass;};
    G4int GetAtomicMass() {return fAtomicMass;};
    
    void AddAtomicNumber(G4int AtomicNumber) {fAtomicNumber = AtomicNumber;};
    G4int GetAtomicNumber() {return fAtomicNumber;};
    
    void AddIsIon(G4bool isIon) {fiSIon = isIon;};
    G4bool GetIsIon() {return fiSIon;};
    
    void AddIsElectron(G4bool isElectron) {fisElectron = isElectron;}; 
    G4bool GetIsElectron() {return fisElectron;}; 
    
    void AddIonLengthStep(G4double IonLengthStep)    {fIonLengthStep += IonLengthStep;};
    G4double GetIonLengthStep()    {return fIonLengthStep;};
    
    void Addallelectron(G4int count_allelectrons) {fcount_allelectrons = count_allelectrons;}; //Add up all electrons generated (Not just secondary particles) 
    G4int GetAllElectronNumber() {return fcount_allelectrons;};
    
    
    void AddElastic(G4int count_elastic) {fcount_elastic = count_elastic; }; 
    void AddInelastic(G4int count_inelastic) {fcount_inelastic = count_inelastic; }; 
    
    G4int GetOriginID(G4int origin_ID) {return forigin_ID;};  
    
private:
  // methods
  B4cCalorHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;
                                            
  //SteppingAction* GetElasticity(G4int hcID, const G4Event* event) const; 
  
  void PrintEventStatistics(G4double CalorEdep, G4double CaloTrackLength) const;
  
  // data members                   
    G4int  fCalorHCID;
    G4double fTotalEnergyDepositStep;
    G4int fcount_ions;
    G4int fcount_electrons;
    G4int fcount_allelectrons;
    G4int fcount_gamma;
    G4int fcount_alpha; 
    G4int fAtomicMass;
    G4int fAtomicNumber;
    G4bool fiSIon;
    G4bool fisElectron;
    G4bool fisgamma; 
    G4double fIonLengthStep;
    G4int forigin_ID; 
    
    G4int fcount_inelastic; 
    G4int fcount_elastic; 
    
    
    
    

};
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
