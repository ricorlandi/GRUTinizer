#ifndef _TCagraSEGMENTHIT_H_
#define _TCagraSEGMENTHIT_H_

#include "TDetectorHit.h"

class TCagraSegmentHit : public TDetectorHit {
public:
  TCagraSegmentHit() { }

  virtual void Copy(TObject&) const;
  virtual void Clear(Option_t *opt = "");
  virtual void Print(Option_t *opt = "") const;

  virtual Int_t Charge() const;

  int GetDetnum() const;
  char GetLeaf() const;
  int GetSegnum() const;

  int GetBoardID() const;
  int GetChannel() const;

  std::vector<Short_t>& GetTrace() { return fTrace; }
  void SetTrace(std::vector<Short_t>& trace);


  void SetDiscTime(const Double_t t) { time = t; }
  Double_t GetDiscTime() { return time; }
  Double_t GetCorrectedEnergy(Double_t asym_bl=0.);
  void SetPreRise(Double_t prerise) { prerise_energy = prerise; }
  void SetPostRise(Double_t postrise) { postrise_energy = postrise; }
  Double_t GetPreRise() { return prerise_energy; }
  Double_t GetPostRise() { return postrise_energy; }
  void SetFlags(UShort_t fl) { flags = fl; }
  const UShort_t& GetFlags() const { return flags; }
  void SetBaseSample(UShort_t base) { base_sample = base; }
  const UShort_t& GetBaseSample() const { return base_sample; }
  Double_t GetTraceEnergy(const UShort_t& a,const UShort_t& b,const UShort_t& x,const UShort_t& y) const; //
  Double_t GetTraceBaseline(); //


private:
  std::vector<Short_t> fTrace;
  Double_t time;
  UShort_t flags;
  Double_t prerise_energy;
  Double_t postrise_energy;
  UShort_t base_sample;
  ClassDef(TCagraSegmentHit,1);
};

#endif /* _TCagraSEGMENTHIT_H_ */
