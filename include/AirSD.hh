#ifndef AirSD_h
#define AirSD_h 1

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


class AirSD : public G4VSensitiveDetector
{
  // int number_of_photons_entered;

  std::vector< G4int > hitIds;
  std::vector< G4int > hitTimes;

  std::vector< G4double > hitXs;
  std::vector< G4double > hitYs;
  std::vector< G4double > hitZs;

  G4String name;
  public:
    AirSD(G4String SDname);
    ~AirSD();

    G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);

    void Initialize(G4HCofThisEvent* HCE);

    void EndOfEvent(G4HCofThisEvent* HCE);

  private:
    // ThompsonEyeHitCollection* HitCollection;
};

#endif







AirSD::AirSD(G4String SDname)
: G4VSensitiveDetector(SDname)
{
  G4cout << "Creating SD with name: " << SDname << G4endl;

  // name = SDname;
  // number_of_photons_entered = 0;
  // collectionName.insert("ThompsonEyeHitCollection");

}

AirSD::~AirSD()
{
}

G4bool AirSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )//ROhist
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

/*

  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Track* aTrack = aStep->GetTrack();
  G4ParticleDefinition* aParticle = aTrack->GetDefinition();

  G4int hitId = 0;
  G4int hitTime = 0;
  G4ThreeVector hitPosition; // = new G4ThreeVector();

  if ( preStepPoint->GetStepStatus() == fGeomBoundary && aParticle->GetParticleName() == "opticalphoton" )
  {
    // number_of_photons_entered += 1;

    hitId = aTrack->GetTrackID();
    hitIds.push_back( hitId );
    
    hitPosition = aTrack->GetPosition();
    hitXs.push_back( hitPosition.x() );
    hitYs.push_back( hitPosition.y() );
    hitZs.push_back( hitPosition.z() );

    hitTime = aTrack->GetGlobalTime() / ns; // global time in nanoseconds
    hitTimes.push_back( hitTime );

    // arrivalEnergy = aTrack->GetTotalEnergy() / eV; // energy of the arriving optical photon in eV
    // arrivalEnergies.push_back( arrivalEnergy );
  }
*/

  return true;
}

void AirSD::Initialize(G4HCofThisEvent* HCE)
{
  // HitCollection = new ThompsonEyeHitCollection (GetName(), "ThompsonEyeHitCollection");

  // static G4int HCID = -1;

  // if (HCID<0)
  // {
    // HCID = GetCollectionID(0);
  // }

  // HCE->AddHitsCollection(HCID, HitCollection);
}

void AirSD::EndOfEvent(G4HCofThisEvent*)
{

/*
  //HitCollection->PrintAllHits();
  cout << "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR\n";
  cout << name << "\n";
  // cout << number_of_photons_entered << "\n";

  FILE *fp;
  fp = fopen("photon_reflections", "w");

  // fprintf(fp, "№ photons arrived:\n%d\n\n", number_of_photons_entered);

  fprintf(fp, "track_id,hit_time,x,y,z\n");

  int N = hitIds.size();
  for (int i = 0; i<N; i++)
  {
    fprintf(fp, "%d,%d,%g,%g,%g\n", hitIds[i], hitTimes[i], hitXs[i], hitYs[i], hitZs[i]);
    // cout << hitIds[i] << "\n";
  }
  fclose(fp);
  hitIds.clear();
  hitTimes.clear();
  hitXs.clear();
  hitYs.clear();
  hitZs.clear();
  // number_of_photons_entered = 0;
  cout << "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR\n";
*/

}

