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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
//#include "Run.hh"
//#include "HistoManager.hh"

#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcess.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4Ions.hh"
#include "B4Analysis.hh"
#include "G4ParticleTable.hh"
#include "G4StepLimiter.hh"
//#include "Run.hh"
#include "B4cEventAction.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "B4cDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4EventManager.hh"
#include "G4NeutronHPThermalScatteringData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(B4cEventAction* event)
: G4UserSteppingAction(),fEventAction(event),fTryLength(0), fScoringVolume(0), fScoringVolume4(0)
{
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ 
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    
    fTryLength=0;
    
    if (!fScoringVolume) {
        const B4cDetectorConstruction* detectorConstruction
        = static_cast<const B4cDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume = detectorConstruction->GetScoringVolume();
        fScoringVolume4 = detectorConstruction->GetScoringVolume4(); 
        fScoringVolume3 = detectorConstruction->GetScoringVolume3(); 
        fScoringVolume2 = detectorConstruction->GetScoringVolume2(); 
        fScoringVolume1a = detectorConstruction->GetScoringVolume1a(); 
        fScoringVolume1b = detectorConstruction->GetScoringVolume1b(); 
        
        
        
        
    }
    
    //G4cout<<"scoring volume: "<<fScoringVolume->GetName()<<G4endl;
    //G4cout<<"scoring volume: "<<fScoringVolume4->GetName()<<G4endl;
    
    //G4cout<<"Volume: "<<detectorConstruction->GetScoringVolume()->GetName()<<G4endl; 
    // get volume of the current step
    G4LogicalVolume* volume
    = aStep->GetPreStepPoint()->GetTouchableHandle()
    ->GetVolume()->GetLogicalVolume();
    
    
    G4Track *aTrack = aStep->GetTrack() ;
    
    
    
    if (volume->GetName() == "CalorLV") {
        G4cout<<"Particle Name: "<<aStep->GetTrack()->GetParticleDefinition()->GetParticleName()<<G4endl;
    }
    
    G4cout<<"Number of Secondaries: "<<aStep->GetNumberOfSecondariesInCurrentStep()<<G4endl;
    
    
    if (volume != fScoringVolume3 and volume != fScoringVolume4 and volume != fScoringVolume2 and volume != fScoringVolume1a and volume != fScoringVolume1b and volume != fScoringVolume) return;
    
    
    
    
    
    //Run* run = static_cast<Run*>( G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    
    //G4double EdepStep = aStep->GetTotalEnergyDeposit();
    
   // if (EdepStep > 0.) {
     //   fEventAction->AddEdepStep(EdepStep);}
    
   G4int fcount_elastic = 0; 
   G4int fcount_inelastic = 0; 
   
   G4int fcount_ions=0;
   G4int fcount_electrons = 0; 
   G4int fcount_gamma = 0; 
   G4int fcount_alpha = 0; 
   G4int count_Si28 = 0; 
   G4int interaction_num = 0; 
   G4int process_num = 0; 
   G4int forigin_ID = 0;
   //G4int event_id_prev; 
   
   G4int Si28_norm = 0;
   
   
   G4int fcount_allelectrons = 0;
    
   
   fcount_allelectrons = total+fcount_allelectrons;
   //G4int fall_particle = 0;
   G4bool isElectron = true; 
   
   
   G4int A = 0;
   G4int Z = 0;
   auto analysisManager = G4AnalysisManager::Instance();
   bool isIon = false;
    

    G4int parent_ID = aStep->GetTrack()->GetParentID();
    G4int track_ID = aStep->GetTrack()->GetTrackID();
    G4String particleType = aStep->GetTrack()->GetParticleDefinition()->GetParticleType();
    G4String particleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    auto event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    auto process = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName(); 
    auto volume_next =  aStep->GetTrack()->GetNextVolume()->GetName();
    auto volume_current = aStep->GetTrack()->GetVolume()->GetName();
     
    G4cout<<"Process: "<<process<<G4endl; 
    
    G4cout<<"Volume Current: "<<volume_current<<G4endl; 
    G4cout<<"Volume Next: "<<volume_next<<G4endl; 
    if (volume_current != "case" and track_ID > 1) {
        
        G4cout<<" Number of Seconds : "<<aStep->GetNumberOfSecondariesInCurrentStep()<<G4endl; 
        G4cout << "Particle Name: " << G4endl; 
        forigin_ID = 1; 
        fEventAction->GetOriginID(forigin_ID);
    }
    
    
     
    
    
    
    
    if (parent_ID>0 && particleName == "e-" ) {
       
//        G4cout<<track_ID == trackid_prev<<G4endl;
        //if (track_ID != trackid_prev) {
            // fcount_ions++;
            total++;  
            //fcount_allelectrons++;
            
            
            fEventAction->Addallelectron(fcount_allelectrons);
       //}
       
    }

    if (parent_ID>0) {
       // fcount_ions++;
        
       G4double damntrack = aStep->GetTrack()->GetTrackLength();
        fTryLength+=damntrack;
        
       G4double EdepStep = aStep->GetTotalEnergyDeposit();
       // if(aStep->GetTrack()->GetTrackID()) // idea e' presteppoint track id = poststeptrack id allora incrementa
       //G4cout << " damn track " << G4BestUnit(damntrack,"Length") << " " << fTryLength << " particle type " << particleType << G4endl;
       fEventAction->AddEdepStep(EdepStep); // is this the energy of the only ions for real? to be check
       //fEventAction->AddIonsNumber(fcount_ions);
        
        if (damntrack>0)  {
            analysisManager->FillH1(5, damntrack);
            analysisManager->FillNtupleDColumn(5, damntrack);}
            fEventAction->AddIonLengthStep(damntrack);
            
        
    
    }
     
    
    
    
    

        
    
    
   
    
        //if (track->GetTrackID()>1) G4cout << "damn track " << damntrack << G4endl;


        const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();
          //G4cout << "Number of secondary particles: "<< (*secondary).size() <<G4endl;
         if ((*secondary).size() == 0 &&  parent_ID==0 ) {
                interaction_num = 20;
                G4cout << "Interaction Number: "<< interaction_num <<G4endl;
                //G4cout<<"process: "<<process<<G4endl;
                analysisManager->FillH1(11, interaction_num);
        }
        
        
        for (size_t lp=0; lp<(*secondary).size(); lp++) {
            
            G4ParticleDefinition *particle;
            particle = (*secondary)[lp]->GetDefinition();
            G4int Parent_ID = (*secondary)[lp]->GetParentID();
            A = particle->GetAtomicMass();
            Z = particle->GetAtomicNumber();
            
            
            G4String name   = particle->GetParticleName();
            
            if (particle->GetParticleType() == "nucleus") {
                
                isIon= true;
               
                fcount_ions++;
               

                fEventAction->AddIsIon(isIon); // this tell me if there is an event with ions
                fEventAction->AddIonsNumber(fcount_ions); // this one is ok, but do not sum it on the header. Keep it =
                fEventAction->AddAtomicMass(A); // gives me only last ion, what if tehre are two or three so it's wrong to fill it in EventAction - to be deleted
                fEventAction->AddAtomicNumber(Z);
                
                //if (EdepStep > 0.) fEventAction->AddEdepStep(EdepStep);
                //G4cout << " sto maledetto A " << A << G4endl;
                
                analysisManager->FillH1(1, A);
                analysisManager->FillNtupleDColumn(1, A);
                analysisManager->FillH1(6, Z);
                analysisManager->FillNtupleDColumn(6, Z);
            
                
            }
            
            
            if (name == "e-") {
                fcount_electrons++; 
                
                fEventAction->AddIsElectron(isElectron); //this tells me if there is an event with electron
                fEventAction->AddElectronNumber(fcount_electrons); // this one is ok, but do not sum it on the header. Keep it =
            
                
                
            }
            
            
            
            //if (name == "gamma") {
                //fcount_gamma++;
                
                //fEventAction->AddIsElectron(isGamma); //this tells me if there is an event with electron
                //fEventAction->AddGammaNumber(fcount_gamma); // this one is ok, but do not sum it on the header. Keep it =
                
                
                
            //}
            //counting number of gamma particles and alpha particles produced by target neutron 
            G4cout<<" Number of secondaries "<<(*secondary).size()<<G4endl;
            for (size_t i=0; i<(*secondary).size(); i++) {
                        G4int Parent_ID_secondary = (*secondary)[i]->GetParentID();
                        G4ParticleDefinition *particle_secondarycount;
                        particle_secondarycount = (*secondary)[i]->GetDefinition();
                        auto name_secondarycount = particle_secondarycount->GetParticleName();
                        
                        
                        
                        
                        
                        
                
                        if (name_secondarycount == "gamma" && Parent_ID_secondary == 1) {
                            fcount_gamma++;
                    
                    
                            }
                        if (name_secondarycount == "alpha" && Parent_ID_secondary == 1) {
                            fcount_alpha++; 
                        }
                        
            } 
            
            //Classification of inelastic vs elastic 
            
        
            
           
            
           // Si Inelastic vs Elastic 
           
           G4int Zmat = 14; 
           
        
            //G4cout<<"Prev volume: "<< volume_prev<<G4endl;
            //G4cout<<"Current volume: "<<volume_current<<G4endl;
            
            
            if (Parent_ID == 1) { //Checking secondaries 
                
                
                
                
                if (process == "hadElastic") {
                    fcount_elastic++;
                    process_num = 1;
                    analysisManager->FillH1(12, process_num);
                    fEventAction->AddElastic(fcount_elastic); 
                    
                    
                     
                }
                else if (process == "neutronInelastic" && Z>2) {
                    fcount_inelastic++;
                    process_num = 2;
                    G4cout<<"Process Classification: "<<process_num<<G4endl;
                    analysisManager->FillH1(12, process_num);
                    fEventAction->AddInelastic(fcount_inelastic); 
                    
                }
                
                
                
                if (Z == Zmat - 2) {                
                        interaction_num = 12;      
                    
                    //G4cout<<"Process Classification: "<<process_num<<G4endl;
                    
                } 
                
                
                else if (Z == Zmat -1) {
                    interaction_num = 13;
                    
                    //G4cout<<"Process Classification: "<<process_num<<G4endl;
                    //analysisManager->FillH1(12, process_num);
                }
                
                else if (Z == Zmat) { 
                    
                
                    if (fcount_gamma != 0 ) 
                    { 
                        
                        interaction_num = 11;
                    }
                    
                    else {
                        
                        interaction_num = 10;
                        
                    }    
                        //G4cout<<"Process Classification: "<<process_num<<G4endl;
                        //analysisManager->FillH1(12, process_num);
                    
        
                }
                 
                
                    
                else if (A == 1 || A == 0 || A == 4) { 
                    interaction_num = 0;
                }
                else {
                        G4cout <<"Atomic Number "<< A << G4endl;
                        G4cout <<"Atomic Mass " << Z << G4endl; 
                        interaction_num = 16;
                }
            } 
            
            
            
            G4cout << "interaction Number: "<< interaction_num <<G4endl;
            analysisManager->FillH1(11, interaction_num);
            

            
            
        }
        
    
    
     /*
           // General G4Isotope Inelastic vs Elastic 
           
           G4int Zmat = 14; 
    
            if (Parent_ID == 1) { //Checking secondaries 
                
                if (process == "hadElastic") {
                    process_num = 1;
                    G4cout<<"Process Classification: "<<process_num<<G4endl;
                    analysisManager->FillH1(12, process_num);
                }
                else if (process == "neutronInelastic" && Z>2) {
                    process_num = 2;
                    G4cout<<"Process Classification: "<<process_num<<G4endl;
                    analysisManager->FillH1(12, process_num);
                }
                
                
                
                if (Z == Zmat-2) {                
                        interaction_num = 12;       //Aluminum + Alpha
                    
                    //G4cout<<"Process Classification: "<<process_num<<G4endl;
                    
                } 
                
                
                else if (Z == Zmat-1) {
                    interaction_num = 13;  //Silicon
                    
                    //G4cout<<"Process Classification: "<<process_num<<G4endl;
                    //analysisManager->FillH1(12, process_num);
                }
                
                else if (Z == Zmat) { 
                    
                
                    if (fcount_gamma != 0 ) 
                    { 
                        
                        interaction_num = 11;
    
                        
                        //G4cout<<"Process Classification: "<<process_num<<G4endl;
                        //analysisManager->FillH1(12, process_num);
                    } 
                    
                    else { 
                        
                        interaction_num = 10;
                        
                        
                        //G4cout<<"Process Classification: "<<process_num<<G4endl;
                        //analysisManager->FillH1(12, process_num);
                        
                        //G4RunManager::GetRunManager()->AbortEvent();  
                    }
                 
                
                    
                    
                    
                }
                else if (A == 1 || A == 0 || A == 4) { 
                    interaction_num = 0;
                }
                else {
                        G4cout <<"Atomic Number "<< A << G4endl;
                        G4cout <<"Atomic Mass " << Z << G4endl; 
                        interaction_num = 16;
                }
            } 
            
            
            
            G4cout << "interaction Number: "<< interaction_num <<G4endl;
            analysisManager->FillH1(11, interaction_num);
            

            
            
        }
        
    */
     
     /*      
     // SiO2 Inelastic vs Elastic
     
            if (Parent_ID == 1) { //Checking secondaries 
                
                if (process == "hadElastic") {
                    process_num = 1;
                    G4cout<<"Process Classification: "<<process_num<<G4endl;
                    analysisManager->FillH1(12, process_num);
                }
                else if (process == "neutronInelastic" && Z>2) {
                    process_num = 2;
                    G4cout<<"Process Classification: "<<process_num<<G4endl;
                    analysisManager->FillH1(12, process_num);
                }
                
                
                G4int Zmat1 = 8;
                G4int Zmat2 = 14; 
                
                if (Z == Zmat1-2) {                
                        interaction_num = 12;      
                    
                    //G4cout<<"Process Classification: "<<process_num<<G4endl;
                    
                } 
                else if (Z == Zmat1-1) {
                    interaction_num = 13; 
                    
                }
                
                else if (Z == Zmat1) { 
                    
                
                    if (fcount_gamma != 0 ) 
                    { 
                        interaction_num = 11; 
                    }
                    else 
                    {
                        interaction_num = 10; 
                    }
                        
                        //G4cout<<"Process Classification: "<<process_num<<G4endl;
                        //analysisManager->FillH1(12, process_num);
                } 
                
                else if (Z == Zmat2-2) {                
                        interaction_num = 16;      
                    
                    //G4cout<<"Process Classification: "<<process_num<<G4endl;
                    
                } 
                else if (Z == Zmat2-1) {
                    interaction_num = 17; 
                    
                }
                
                else if (Z == Zmat2) { 
                    
                
                    if (fcount_gamma != 0 ) 
                    { 
                        interaction_num = 15; 
                    }
                    else 
                    {
                        interaction_num = 14; 
                    }
                        
                        //G4cout<<"Process Classification: "<<process_num<<G4endl;
                        //analysisManager->FillH1(12, process_num);
                } 
                
            

                else if (A == 1 || A == 0 || A == 4) { 
                    interaction_num = 0;
                }
                else {
                        G4cout <<"Atomic Number "<< A << G4endl;
                        G4cout <<"Atomic Mass " << Z << G4endl; 
                        interaction_num = 18;
                }
            } 
            
            
            
            G4cout << "interaction Number: "<< interaction_num <<G4endl;
            analysisManager->FillH1(11, interaction_num);
            

            
            
        }
        
       */  
          
           
           
           
    //G4cout<< " number of ions from steps " << fcount_ions << G4endl;
    //fEventAction->AddIonsNumber(fcount_ions);
   // analysisManager->FillH1(3, fcount_ions);
    //analysisManager->FillNtupleDColumn(3, fcount_ions);
    
    analysisManager->AddNtupleRow();
     
    /* Determining inelastic scattering */ 
    
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*Run* run
 = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 
 // count processes
 //
 const G4StepPoint* endPoint = aStep->GetPostStepPoint();
 G4VProcess* process   =
 const_cast<G4VProcess*>(endPoint->GetProcessDefinedStep());
 run->CountProcesses(process);
 
 // check that an real interaction occured (eg. not a transportation)
 G4StepStatus stepStatus = endPoint->GetStepStatus();
 G4bool transmit = (stepStatus==fGeomBoundary || stepStatus==fWorldBoundary);
 if (transmit) return;
 
 //real processes : sum track length
 //
 G4double stepLength = aStep->GetStepLength();
 run->SumTrack(stepLength);
 
 //energy-momentum balance initialisation
 //
 const G4StepPoint* prePoint = aStep->GetPreStepPoint();
 G4double Q             = - prePoint->GetKineticEnergy();
 G4ThreeVector Pbalance = - prePoint->GetMomentum();
 
 //initialisation of the nuclear channel identification
 //
 G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
 G4String partName = particle->GetParticleName();
 G4String nuclearChannel = partName;
 G4HadronicProcess* hproc = dynamic_cast<G4HadronicProcess*>(process);
 const G4Isotope* target = NULL;
 if (hproc) target = hproc->GetTargetIsotope();
 G4String targetName = "XXXX";
 if (target) targetName = target->GetName();
 nuclearChannel += " + " + targetName + " --> ";
 if (targetName == "XXXX") run->SetTargetXXX(true);
 
 //scattered primary particle (if any)
 //
 G4AnalysisManager* analysis = G4AnalysisManager::Instance();
 G4int ih = 1;
 if (aStep->GetTrack()->GetTrackStatus() == fAlive) {
 G4double energy = endPoint->GetKineticEnergy();
 analysis->FillH1(ih,energy);
 //
 G4ThreeVector momentum = endPoint->GetMomentum();
 Q        += energy;
 Pbalance += momentum;
 //
 nuclearChannel += partName + " + ";
 }
 
 //secondaries
 //
 const std::vector<const G4Track*>* secondary
 = aStep->GetSecondaryInCurrentStep();
 for (size_t lp=0; lp<(*secondary).size(); lp++) {
 particle = (*secondary)[lp]->GetDefinition();
 G4String name   = particle->GetParticleName();
 G4String type   = particle->GetParticleType();
 G4double energy = (*secondary)[lp]->GetKineticEnergy();
 run->ParticleCount(name,energy);
 //energy spectrum
 ih = 0;
 if (particle == G4Gamma::Gamma())       ih = 2;
 else if (particle == G4Neutron::Neutron())   ih = 3;
 else if (particle == G4Proton::Proton())     ih = 4;
 else if (particle == G4Deuteron::Deuteron()) ih = 5;
 else if (particle == G4Alpha::Alpha())       ih = 6;
 else if (type == "nucleus")                  ih = 7;
 else if (type == "meson")                    ih = 8;
 else if (type == "baryon")                   ih = 9;
 if (ih > 0) analysis->FillH1(ih,energy);
 //atomic mass
 if (type == "nucleus") {
 G4int A = particle->GetAtomicMass();
 analysis->FillH1(12, A);
 }
 //energy-momentum balance
 G4ThreeVector momentum = (*secondary)[lp]->GetMomentum();
 Q        += energy;
 Pbalance += momentum;
 //count e- from internal conversion together with gamma
 if (particle == G4Electron::Electron()) particle = G4Gamma::Gamma();
 //particle flag
 fParticleFlag[particle]++;
 }
 
 //energy-momentum balance
 G4double Pbal = Pbalance.mag();
 run->Balance(Pbal);
 ih = 10;
 analysis->FillH1(ih,Q);
 ih = 11;
 analysis->FillH1(ih,Pbal);
 
 // nuclear channel
 const G4int kMax = 16;
 const G4String conver[] = {"0","","2 ","3 ","4 ","5 ","6 ","7 ","8 ","9 ",
 "10 ","11 ","12 ","13 ","14 ","15 ","16 "};
 std::map<G4ParticleDefinition*,G4int>::iterator ip;
 for (ip = fParticleFlag.begin(); ip != fParticleFlag.end(); ip++) {
 particle = ip->first;
 G4String name = particle->GetParticleName();
 G4int nb = ip->second;
 if (nb > kMax) nb = kMax;
 G4String Nb = conver[nb];
 if (particle == G4Gamma::Gamma()) {
 run->CountGamma(nb);
 Nb = "N ";
 name = "gamma or e-";
 }
 if (ip != fParticleFlag.begin()) nuclearChannel += " + ";
 nuclearChannel += Nb + name;
 }
 
 ///G4cout << "\n nuclear channel: " << nuclearChannel << G4endl;
 run->CountNuclearChannel(nuclearChannel, Q);
 
 fParticleFlag.clear();
 
//  kill event after first interaction

 G4RunManager::GetRunManager()->AbortEvent();  */

/*    G4int count_ions=0;
 G4Track * theTrack = aStep  ->  GetTrack();
 G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
 if (dynamic_cast<G4Ions*>(particleDef)) {
 
 count_ions++;
 G4cout <<"\n number of ions: " << count_ions << G4endl;}
 */
//G4String name = "dontcare"



/*   
           // P31 Inelastic vs Elastic 
           
           G4int Zmat = 15; 
    
            if (Parent_ID == 1) { //Checking secondaries 
                
                if (process == "hadElastic") {
                    process_num = 1;
                    G4cout<<"Process Classification: "<<process_num<<G4endl;
                    analysisManager->FillH1(12, process_num);
                }
                else if (process == "neutronInelastic" && Z>2) {
                    process_num = 2;
                    G4cout<<"Process Classification: "<<process_num<<G4endl;
                    analysisManager->FillH1(12, process_num);
                }
                
                
                
                if (Z == 13) {                
                        interaction_num = 12;       //Aluminum + Alpha
                    
                    //G4cout<<"Process Classification: "<<process_num<<G4endl;
                    
                } 
                
                
                else if (Z == 14) {
                    interaction_num = 13;  //Silicon
                    
                    //G4cout<<"Process Classification: "<<process_num<<G4endl;
                    //analysisManager->FillH1(12, process_num);
                }
                
                else if (Z == Zmat) { 
                    
                
                    if (fcount_gamma != 0 ) 
                    { 
                        
                        interaction_num = 11;
    
                        
                        //G4cout<<"Process Classification: "<<process_num<<G4endl;
                        //analysisManager->FillH1(12, process_num);
                    } 
                    
                    else { 
                        
                        interaction_num = 10;
                        
                        
                        //G4cout<<"Process Classification: "<<process_num<<G4endl;
                        //analysisManager->FillH1(12, process_num);
                        
                        //G4RunManager::GetRunManager()->AbortEvent();  
                    }
                 
                
                    
                    
                    
                }
                else if (A == 1 || A == 0 || A == 4) { 
                    interaction_num = 0;
                }
                else {
                        G4cout <<"Atomic Number "<< A << G4endl;
                        G4cout <<"Atomic Mass " << Z << G4endl; 
                        interaction_num = 16;
                }
            } 
            
            
            
            G4cout << "interaction Number: "<< interaction_num <<G4endl;
            analysisManager->FillH1(11, interaction_num);
            

            
            
        } */
         
/*      
     // B10 and B11 Inelastic vs Elastic 
           
            if (Parent_ID == 1) { //Checking secondaries 
                
                if (process == "hadElastic") {
                    process_num = 1;
                    G4cout<<"Process Classification: "<<process_num<<G4endl;
                    analysisManager->FillH1(12, process_num);
                }
                else if (process == "neutronInelastic" && Z>2) {
                    process_num = 2;
                    G4cout<<"Process Classification: "<<process_num<<G4endl;
                    analysisManager->FillH1(12, process_num);
                }
                
                
                G4int Zmat = 5;
                
                if (Z == 3) {                
                        interaction_num = 12;      
                    
                    //G4cout<<"Process Classification: "<<process_num<<G4endl;
                    
                } 
                
                
                
                
                else if (Z == Zmat) { 
                    
                
                    if (fcount_gamma != 0 ) 
                    { 
                        interaction_num = 11; 
                    }
                    else 
                    {
                        interaction_num = 10; 
                    }
                        
                        //G4cout<<"Process Classification: "<<process_num<<G4endl;
                        //analysisManager->FillH1(12, process_num);
                } 
                
                else if (Z == 4 ) {
                    interaction_num = 13; 
                    
                }
                
                else if ( A == 3 and Z == 2) {
                    if (fcount_alpha != 0 ) {
                        interaction_num = 14; 
                    }
                }
                    
                   
                
                    
                    

                else if (A == 1 || A == 0 || A == 4) { 
                    interaction_num = 0;
                }
                else {
                        G4cout <<"Atomic Number "<< A << G4endl;
                        G4cout <<"Atomic Mass " << Z << G4endl; 
                        interaction_num = 16;
                }
            } 
            
            
            
            G4cout << "interaction Number: "<< interaction_num <<G4endl;
            analysisManager->FillH1(11, interaction_num);
            

            
            
        }
        */
