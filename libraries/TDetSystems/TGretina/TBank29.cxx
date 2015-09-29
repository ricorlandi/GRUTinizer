#include "TBank29.h"

#include "TGEBEvent.h"
#include "TGRUTOptions.h"

ClassImp(TBank29)

TBank29::TBank29(){
  channels = new TClonesArray("TMode3");
}

TBank29::~TBank29() {
  channels->Delete();
}

void TBank29::Copy(TObject& obj) const {
  TDetector::Copy(obj);

  TBank29& bank = (TBank29&)obj;
  channels->Copy(*(bank.channels));
  bank.raw_data.clear();
}

void TBank29::InsertHit(const TDetectorHit& hit){
  TMode3* new_hit = (TMode3*)channels->ConstructedAt(Size());
  hit.Copy(*new_hit);
}

int TBank29::BuildHits(){
  for(auto& event : raw_data){
    TGEBEvent* geb = (TGEBEvent*)&event;
    SetTimestamp(geb->GetTimestamp());
    TMode3Hit hit;
    hit.BuildFrom(geb->GetPayloadBuffer());
    InsertHit(hit);
  }
  return Size();
}

void TBank29::Print(Option_t *opt) const { }

void TBank29::Clear(Option_t *opt) {
  TDetector::Clear(opt);
  channels->Clear(opt);//("TBank29Hit");
  raw_data.clear();
}
