#ifndef LeadSD_h
#define LeadSD_h 1

#include "G4VSensitiveDetector.hh"
//#include "ThompsonEyeHit.hh"


#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4HCtable.hh"
#include "G4UnitsTable.hh"

#include "G4VProcess.hh"

#include <stdio.h>

#include <iostream>
using namespace std;





class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;


class LeadSD : public G4VSensitiveDetector
{
  int number_of_photons_entered;
  std::vector< G4double > arrivalTimes;
  G4String name;
  public:
    LeadSD(G4String SDname);
    ~LeadSD();

    G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);

    void Initialize(G4HCofThisEvent* HCE);

    void EndOfEvent(G4HCofThisEvent* HCE);

  private:
    // ThompsonEyeHitCollection* HitCollection;
};

#endif







LeadSD::LeadSD(G4String SDname)
: G4VSensitiveDetector(SDname)
{
  G4cout << "Creating SD with name: " << SDname << G4endl;

  name = SDname;
  number_of_photons_entered = 0;
  // collectionName.insert("ThompsonEyeHitCollection");

}

LeadSD::~LeadSD()
{
}

G4bool LeadSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )//ROhist
{

/*
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Track* aTrack = aStep->GetTrack();
  G4ParticleDefinition* aParticle = aTrack->GetDefinition();

  if(preStepPoint->GetStepStatus() == fGeomBoundary && aParticle->GetParticleName() == "opticalphoton")
  {
    G4ThreeVector WorldPos = preStepPoint->GetPosition();
    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4ThreeVector LocalPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(WorldPos);
    G4double energy = aTrack->GetTotalEnergy();
    G4ThreeVector direction = aTrack->GetMomentumDirection();
    G4ThreeVector localDirection = theTouchable->GetHistory()->GetTopTransform().TransformAxis(direction);

    ThompsonEyeHit* newHit = new ThompsonEyeHit();
    newHit->SetWorldPos(WorldPos);
    newHit->SetLocalPos(LocalPos);
    newHit->SetEnergy(energy);
    newHit->SetDirection(direction);
    newHit->SetLocalDirection(localDirection);

    HitCollection->insert(newHit);
  }
*/


  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Track* aTrack = aStep->GetTrack();
  G4ParticleDefinition* aParticle = aTrack->GetDefinition();

  G4double arrivalTime = 0;

  if(preStepPoint->GetStepStatus() == fGeomBoundary && aParticle->GetParticleName() == "opticalphoton")
  {
    number_of_photons_entered += 1;
    arrivalTime = aTrack->GetGlobalTime();
    arrivalTimes.push_back( arrivalTime );
  }


  return true;
}

void LeadSD::Initialize(G4HCofThisEvent* HCE)
{
  // HitCollection = new ThompsonEyeHitCollection (GetName(), "ThompsonEyeHitCollection");

  // static G4int HCID = -1;

  // if (HCID<0)
  // {
    // HCID = GetCollectionID(0);
  // }

  // HCE->AddHitsCollection(HCID, HitCollection);
}

void LeadSD::EndOfEvent(G4HCofThisEvent*)
{
  //HitCollection->PrintAllHits();
  cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
  cout << name << "\n";
  cout << number_of_photons_entered << "\n";


  FILE *fp;
  fp = fopen("photon_arrivals", "w");

  fprintf(fp, "arrival_time\n");

  int N = arrivalTimes.size();
  for (int i = 0; i<N; i++)
  {
    fprintf(fp, "%g\n", arrivalTimes[i]);
    // cout << arrivalTimes[i] << "\n";
  }
  fclose(fp);
  cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
}

