#include "B4cActionInitialization.hh"
#include "B4PrimaryGeneratorAction.hh"
#include "B4RunAction.hh"
#include "B4cEventAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cActionInitialization::B4cActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cActionInitialization::~B4cActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cActionInitialization::BuildForMaster() const
{
  SetUserAction(new B4RunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void B4cActionInitialization::Build() const
{
  SetUserAction(new B4PrimaryGeneratorAction);
  SetUserAction(new B4RunAction);
  //  SetUserAction(new SteppingAction);
  
  //  SetUserAction(new SteppingAction);
  
  //  B4cEventAction* event_action = new B4cEventAction();
  // SteppingAction* stepping_action = new SteppingAction(event_action);
  //  SetUserAction(stepping_action);
  SetUserAction(new B4cEventAction);
  //  SetUserAction(new SteppingAction);

  }*/
 
void B4cActionInitialization::Build() const
{
  
  SetUserAction(new B4PrimaryGeneratorAction);
    
  B4RunAction* runAction = new B4RunAction;
  SetUserAction(runAction);

  B4cEventAction* eventAction = new B4cEventAction;
  SetUserAction(eventAction);
    
    TrackingAction* trackingAction = new TrackingAction();
    SetUserAction(trackingAction);

  SteppingAction* steppingAction = new SteppingAction(eventAction);
  SetUserAction(steppingAction);
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
