#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessTable.hh"
#include "G4VProcess.hh"
#include "phys_const.hh"
#include <G4ChanRad.hh>
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4ParticleDefinition.hh"
#include "EventAction.hh"

class EventAction;

class MySteppingAction : public G4UserSteppingAction {
private:
    EventAction* fEventAction;
public:
    MySteppingAction(EventAction* eventAction) : G4UserSteppingAction(),
        fEventAction(eventAction) 
    {}
    virtual ~MySteppingAction() {}

    void UserSteppingAction(const G4Step* step) {

        //std::cout << "inside UserSteppingAction" << std::endl;
        G4ThreeVector vec(step->GetPostStepPoint()->GetMomentumDirection());
        //std::cout << "got vec" << std::endl;
        G4double x = vec.getX();
        G4double y = vec.getY();
        G4double z = vec.getZ();
        G4double tetta_ang = acos(z / (sqrt(x * x + y * y + z * z) * 1));
        // std::cout << "got xyz; tetta_ang = " << tetta_ang << std::endl;
        G4double thetacrit;

        G4ChanRad* Class = new G4ChanRad;
        Class->GetThettaCrit(thetacrit);
        std::cout << "got thetacrit = " << thetacrit << std::endl << "got tetta_ang = " << tetta_ang;

        if (step->GetTrack()->GetDynamicParticle()->GetKineticEnergy() < 200) //&& step->GetTrack()->GetTrackID() != 1)//(tetta_ang > thetacrit)
        {
           // todo: add actions
        }
        else if (tetta_ang < thetacrit)
        {
            G4VPhysicalVolume* PreVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
            G4String preVolumeName = PreVolume->GetName();
            //std::cout << "got preVolumeName" << std::endl;
            G4VPhysicalVolume* PostVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
            //std::cout << "got PostVolume" << std::endl;
            if (PostVolume) {
                G4String postVolumeName = PostVolume->GetName();
                if ((preVolumeName != "FakeCrystal" && postVolumeName == "FakeCrystal") /*|| (preVolumeName != "Crystal" && postVolumeName == "Crystal")*/)
                {
                    // новое направление частицы
                    G4ThreeVector newDirection(0, 0, 1);
                    //std::cout << "got newDirection" << std::endl;
                    step->GetPostStepPoint()->SetMomentumDirection(newDirection);
                    //step->GetTrack()->SetMomentumDirection(newDirection); //для конкретного объема
                    //std::cout << "got done" << std::endl;
                    G4VProcess* process_eIoni = nullptr;
                    G4VProcess* process_msc = nullptr;
                    G4VProcess* process_eBrem = nullptr;
                    G4ParticleDefinition* particleDefinition = step->GetTrack()->GetDefinition();
                    G4ProcessManager* processManager = particleDefinition->GetProcessManager();
                    G4ProcessVector* processVector = processManager->GetProcessList();
                    // Получаем размер списка процессов
                    G4int nProcesses = processVector->size();
                    G4VPhysicalVolume* PreVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
                    G4String preVolumeName = PreVolume->GetName();
                    //std::cout << "got preVolumeName" << std::endl;
                    G4VPhysicalVolume* PostVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
                    if (PostVolume) {
                        G4String postVolumeName = PostVolume->GetName();
                        if (true)//preVolumeName == "Crystal" || preVolumeName == "FakeCrystal")
                        {
                            for (G4int i = 0; i < nProcesses; ++i) {
                                G4VProcess* process = (*processVector)[i];
                                if (process->GetProcessName() == "eIoni") {
                                    process_eIoni = process;
                                    break;
                                }
                            }
                            if (process_eIoni) {  
                                processManager->SetProcessActivation(processManager->GetProcessIndex(process_eIoni), false);
                            }
                            for (G4int i = 0; i < nProcesses; ++i) {
                                G4VProcess* process = (*processVector)[i];
                                if (process->GetProcessName() == "msc") {
                                    process_msc = process;
                                    break;
                                }
                            }
                            if (process_msc) {
                                processManager->SetProcessActivation(processManager->GetProcessIndex(process_msc), false);
                            }
                            for (G4int i = 0; i < nProcesses; ++i) {
                                G4VProcess* process = (*processVector)[i];
                                if (process->GetProcessName() == "eBrem") {
                                    process_eBrem = process;
                                    break;
                                }
                            }
                            if (process_eBrem) {  
                                processManager->SetProcessActivation(processManager->GetProcessIndex(process_eBrem), false);
                            }
                        }
                    }
                }
            }
        }
        G4VPhysicalVolume* PreVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
        G4String preVolumeName = PreVolume->GetName();
        //std::cout << "got preVolumeName" << std::endl;
        G4VPhysicalVolume* PostVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
        //std::cout << "got PostVolume" << std::endl;
        if (PostVolume) {
            G4String postVolumeName = PostVolume->GetName();
            //std::cout << "got postVolumeName" << std::endl;

            if (preVolumeName == "World" && postVolumeName != "Crystal") /*|| (preVolumeName != "Crystal" && postVolumeName == "Crystal")*/
            {
                G4VProcess* process_eIoni = nullptr;
                G4VProcess* process_msc = nullptr;
                G4VProcess* process_eBrem = nullptr;
                G4ParticleDefinition* particleDefinition = step->GetTrack()->GetDefinition();
                G4ProcessManager* processManager = particleDefinition->GetProcessManager();
                G4ProcessVector* processVector = processManager->GetProcessList();
                // Получаем размер списка процессов
                G4int nProcesses = processVector->size();
                G4VPhysicalVolume* PreVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
                G4String preVolumeName = PreVolume->GetName();
                //std::cout << "got preVolumeName" << std::endl;
                G4VPhysicalVolume* PostVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
                if (PostVolume) {
                    G4String postVolumeName = PostVolume->GetName();
                    if (true)//preVolumeName == "Crystal" || preVolumeName == "FakeCrystal")
                    {
                        for (G4int i = 0; i < nProcesses; ++i) {
                            G4VProcess* process = (*processVector)[i];
                            if (process->GetProcessName() == "eIoni") {
                                process_eIoni = process;
                                break;
                            }
                        }
                        if (process_eIoni) {
                            //processManager->SetProcessActivation(process_eIoni, false); 
                            processManager->SetProcessActivation(processManager->GetProcessIndex(process_eIoni), true);

                        }
                        for (G4int i = 0; i < nProcesses; ++i) {
                            G4VProcess* process = (*processVector)[i];
                            if (process->GetProcessName() == "msc") {
                                process_msc = process;
                                break;
                            }
                        }
                        if (process_msc) {
                            // processManager->SetProcessActivation(process_msc, false);
                            processManager->SetProcessActivation(processManager->GetProcessIndex(process_msc), true);

                        }
                        for (G4int i = 0; i < nProcesses; ++i) {
                            G4VProcess* process = (*processVector)[i];
                            if (process->GetProcessName() == "eBrem") {
                                process_eBrem = process;
                                break;
                            }
                        }
                        if (process_eBrem) {
                            //processManager->SetProcessActivation(process_eBrem, false);
                            processManager->SetProcessActivation(processManager->GetProcessIndex(process_eBrem), true);
                        }
                    }
                }

            }
        }
        delete Class;
        User2SteppingAction(step);
    }
    /**
    * Выполняется на каждом шаге
    */
    void User2SteppingAction(const G4Step* step)
    {
        G4VPhysicalVolume* PreVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
        G4String preVolumeName = PreVolume->GetName();
        G4VPhysicalVolume* PostVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
        if (PostVolume) {
            G4String postVolumeName = PostVolume->GetName();
            if (preVolumeName == "Crystal")
            {
                G4int stepNumber = step->GetTrack()->GetCurrentStepNumber();
                G4String particleType = step->GetTrack()->GetDefinition()->GetParticleName();
                G4double energy = step->GetPreStepPoint()->GetKineticEnergy() / CLHEP::MeV;

                G4String secondaryParticleType = "none";
                G4double secondaryEnergy = 0.0;

                const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
                if (!secondaries->empty()) {
                    const G4Track* secondary = (*secondaries)[0];
                    secondaryParticleType = secondary->GetDefinition()->GetParticleName();
                    secondaryEnergy = secondary->GetKineticEnergy() / CLHEP::MeV;
                }

                G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentumDirection();
                const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
                // Получение информации о родительской частице
                G4int process_type = 9, particleTypeI = 9, secondaryParticleTypeI = 9;
                if (process != nullptr) {
                    if (process->GetProcessName() == "eBrem") process_type = 1;
                    else if (process->GetProcessName() == "eIoni") process_type = 2;
                    else if (process->GetProcessName() == "msc") process_type = 3;
                    else if (process->GetProcessName() == "RegularXTRadiator") process_type = 4;
                    else if (process->GetProcessName() == "compt") process_type = 5;
                    else if (process->GetProcessName() == "CoulombScat") process_type = 6;
                    else if (process->GetProcessName() == "Transportation") process_type = 7;

                    if (particleType == "e-") particleTypeI = 1;
                    else if (particleType == "gamma") particleTypeI = 2;
                    else if (particleType == "e+") particleTypeI = 3;
                    else if (particleType == "proton") particleTypeI = 4;

                    if (secondaryParticleType == "e-") secondaryParticleTypeI = 1;
                    else if (secondaryParticleType == "gamma") secondaryParticleTypeI = 2;
                    else if (secondaryParticleType == "e+") secondaryParticleTypeI = 3;
                    else if (secondaryParticleType == "proton") secondaryParticleTypeI = 4;
                       
                    fEventAction->RecordStepData(stepNumber, particleTypeI, energy, secondaryParticleTypeI, secondaryEnergy, momentum.x(), momentum.y(), momentum.z(), process_type);
                }
            }
        }
    }
};

#endif

