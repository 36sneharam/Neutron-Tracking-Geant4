# Neutron-Tracking-Geant4
Tracking the interaction of neutrons with a Silicon detector 

	
 ## 1. GEOMETRY DEFINITION
   The detector consists of 6 cylindrical (or semi-cylindrical) layers. Each layer is defined by three parameters  
 	- the material of the layer,
	- the radius of the layer, 
	- the thickness of the layer 
 	A cylindrical slab of Al and a slab of Gd are also included and are
    positioned parallel to the detector with the same diameter as the detector. The position of these slabs can be changed by editing the variable **offset.** 

All geometrical and material definitions can be found in the B4cDetectorConstruction Class. To change materials and geometrical parameters, edit **DefineVolumes()** or **DefineMaterials()**. To assign a geometry as a sensitive region, edit **ConstructSDandField()**. 

 	
 ## 2. PHYSICS LIST
 
  The **QGSP_BERT_HP** physics list is set as the default physics list. In order to use another physics list edit the parameter **physName** in *exampleB4c.cc*. All the built-in physics lists for this program can be found in *PhysicsList.cc* and can be called in a similar way to QGSP_BERT_HP. 
 	 
 ## 3. AN EVENT : THE PRIMARY GENERATOR
 
   The primary kinematic consists of a single particle located a certain distance from the detector. The type of the particle, its energy, and its position are set in PrimaryGeneratorAction (neutron, 1 MeV, 1.5cm). By default, the particles are fired in a cone configuration whose diameter equals the diameter of the detector. The direction of the neutron can be changed using the function SetParticleMomentumDirection.  
 	
 ## 4. PHYSICS
 
   An event is killed as soon as the incident particle is no longer in the sensitive region defined by the parameter fSensitiveDetector. Several parameters are determined including total energy deposited in the detector, energy deposited by elastic interactions, energy deposited by inelastic interactions etc. The parameters extracted from the simulation are defined in SteppingAction.cc. 

 ## 5. HISTOGRAMS
         
   The test contains 18 built-in 1D histograms, which are managed by
   G4AnalysisManager and are initialized in RunAction.cc. 
   1.	"Total energy deposited in detector"
   2.  "Atomic Mass of particles created"
   3.	"Track length in detector"
   4.	"Number of ions per event"
   5.	"event ID"
   6.	"Ion track length"
   7.	"Atomic Number of particles created"
   8.	"Energy Deposited in detector for ions"
   9.	"Energy Deposited in detector for ions - larger energy range"
   10. "SP for events with only 1 ion"
   11. "Ion track length with 1 ion only"
   12. "Interaction Number" (Defined in Stepping Action)
   13. "Elasic (1) vs Inelastic (2)"
   14. "Elastic Energy Deposited in Detector"
   15. "Inelastic Energy Deposited in Detector"
   16.  "Total Charge Deposited in Detector"
   17.  "Elastic Charge Deposited in Detector"
   18.  "Inelastic Charge Deposited in Detector"

   In addition, there are 4 2-D histograms :
   1.  "Ion Track length vs Atomic Mass"
   2.  "Ion Track length vs Atomic Number"
   3.  "Stopping power vs Atomic Number"
   4.  "Energy Deposited vs Atomic Number "
	    
      
   The histograms are managed by the HistoManager class and its Messenger. The histograms are filled-in in the EventAction class. The histograms are saved in the file B4-Calor.root. 
 	 				
 ## 6.  VISUALIZATION
 
   The Visualization Manager is set in the main().
   The initialisation of the drawing is done via the commands
   /vis/... in the macro vis.mac. To get visualisation:
   
   > /control/execute vis.mac
 	
   
	
 ## 7. HOW TO START ?
 
   Execute B4c in 'batch' mode from macro files :
 	% B4c   run2.mac
 		
   Execute B4c in 'interactive mode' with visualization :
   
 	% B4c
	Idle> control/execute vis.mac
 	....
 	Idle> type your commands
 	....
 	Idle> exit
