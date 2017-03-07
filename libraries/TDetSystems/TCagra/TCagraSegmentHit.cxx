#include "TCagraSegmentHit.h"

#include <algorithm>
#include <iostream>

#include "TANLEvent.h"


void TCagraSegmentHit::Copy(TObject& obj) const{
  TDetectorHit::Copy(obj);

  TCagraSegmentHit& cagra = (TCagraSegmentHit&)obj;
  cagra.fTrace = fTrace;
}

void TCagraSegmentHit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
  fTrace.clear();
}

void TCagraSegmentHit::SetTrace(std::vector<Short_t>& trace) {
  fTrace.clear();
  fTrace.swap(trace);
}

void TCagraSegmentHit::Print(Option_t *opt) const {
  std::cout << "TCagraSegmentHit:\n"
            << "\tCharge: " << Charge() << "\n"
            << std::flush;
}

int TCagraSegmentHit::GetDetnum() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return chan->GetArrayPosition();
  } else {
    std::cout << "Unknown address: " << std::hex << fAddress << std::dec
              << std::endl;
    return -1;
  }
}
char TCagraSegmentHit::GetLeaf() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return *chan->GetArraySubposition();
  } else {
    std::cout << "Unknown address: " << std::hex << fAddress << std::dec
              << std::endl;
    return (char)-1;
  }
}

int TCagraSegmentHit::GetSegnum() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return chan->GetSegment();
  } else {
    return -1;
  }
}

// int TCagraSegmentHit::GetCrate() const {
//   return (fAddress&0x00ff0000)>>16;
// }

int TCagraSegmentHit::GetBoardID() const {
  return (fAddress&0x0000ff00)>>8;
}

int TCagraSegmentHit::GetChannel() const {
  return (fAddress&0x000000ff)>>0;
}

Int_t TCagraSegmentHit::Charge() const {
  if(fCharge > 30000) {
    return fCharge - 32768;
  } else {
    return fCharge;
  }
}

Double_t TCagraSegmentHit::GetCorrectedEnergy(Double_t asym_bl) {

  TChannel* chan = TChannel::GetChannel(fAddress);
  Double_t Energy = 0;
  if(!chan){
    std::cout << std::hex << "Channel 0x" << fAddress << " not defined in calibrations file, no corrections are applied." << std::endl;
  } else {
    auto pzE = chan->PoleZeroCorrection(prerise_energy,postrise_energy,TANLEvent::GetShapingTime());
    pzE = chan->BaselineCorrection(pzE,asym_bl);
    Energy = chan->CalEnergy(pzE, fTimestamp);
  }
  return Energy;
}

Double_t TCagraSegmentHit:: GetTraceBaseline() {
  int length = std::min<int>(6,fTrace.size());
  if (length) {
    Double_t bl = 0;
    for (int i=0; i<length; i++) {
      bl+=fTrace[i];
    }
    return bl/length;
  }
  return 0;
}

Double_t TCagraSegmentHit::GetTraceEnergy(const UShort_t& a,const UShort_t& b,const UShort_t& x,const UShort_t& y) const {
  if (!fTrace.size()) { return 0; }

  if (fTrace.size() < y) {
    static int nprint = 0;
    if (nprint < 10) {
      std::cout << "Warning: Trace length less than requested sampling window: " << fTrace.size() <<std::endl;
    } nprint++;
    return 0;
  }

  double baseline = 0;
  for (int i=a; i<b; i++) { baseline+=fTrace[i]; }
  double integral = 0;
  for (int i=x; i<y; i++) { integral+=fTrace[i]; }

  return (integral/(y-x) - baseline/(b-a));
}
