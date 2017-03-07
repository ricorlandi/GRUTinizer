
#include "TGrandRaidenHit.h"
#include "TGRUTOptions.h"
#include "TSmartBuffer.h"
#include "TMath.h"
#include "GValue.h"

ClassImp(TGrandRaidenHit)


std::vector<double> TGrandRaidenHit::acoefs;
std::vector<double> TGrandRaidenHit::bcoefs;
unsigned int TGrandRaidenHit::xdegree = 3, TGrandRaidenHit::adegree = 1, TGrandRaidenHit::ydegree = 1, TGrandRaidenHit::bdegree = 1;

Double_t TGrandRaidenHit::aoffset = std::sqrt(-1);
Double_t TGrandRaidenHit::boffset = std::sqrt(-1);


TGrandRaidenHit::TGrandRaidenHit() : vector(1.,1.,1.) {
  madc1=0; madc2=0; madc3=0;
  tpos1=0; tpos2=0; tpos3=0;
  excitation_energy=0;
  momentum=0;
}
TGrandRaidenHit::TGrandRaidenHit(const TGrandRaidenHit& gr) {
  labr_hits = gr.labr_hits;
  madc1 = gr.madc1;
  madc2 = gr.madc2;
  madc3 = gr.madc3;
  tpos1 = gr.tpos1;
  tpos2 = gr.tpos2;
  tpos3 = gr.tpos3;
  Timestamp = gr.Timestamp;
  rcnp = gr.rcnp;
  vector = gr.vector;
  excitation_energy=gr.excitation_energy;
  momentum=gr.momentum;
}
TGrandRaidenHit::TGrandRaidenHit(RCNPEvent& rcnpevent)
  : rcnp(rcnpevent), vector(1.,1.,1.) {
  madc1=0; madc2=0; madc3=0;
  tpos1=0; tpos2=0; tpos3=0;
  excitation_energy=0;
  momentum=0;
}
TGrandRaidenHit::~TGrandRaidenHit() {
}

void TGrandRaidenHit::BuildFrom(){
#ifdef RCNP
  static bool once = true;
  if (once) {
    RCNPEvent::HistDefCheckSum();
    once = false;
  }
  Clear();

  Timestamp = rcnp.GetTimestamp();

  auto adc = rcnp.GR_ADC();
  auto tdc = rcnp.GR_TDC();

  auto qtc_le_tdc =  rcnp.QTC_LEADING_TDC();
  auto qtc_le_chan =  rcnp.QTC_LEADING_CH();
  auto qtc_tr_tdc =  rcnp.QTC_TRAILING_TDC();
  auto qtc_tr_chan =  rcnp.QTC_TRAILING_CH();

  if (qtc_le_tdc && qtc_tr_tdc) {

    for (auto i=0u; i<qtc_le_chan->size(); i++) {
      for (auto j=0u; j<qtc_tr_chan->size(); j++) {
        if ((*qtc_le_chan)[i]==(*qtc_tr_chan)[j]) {

          LaBrHit temphit;
          temphit.SetTimestamp(Timestamp);
          temphit.channel = (*qtc_le_chan)[i];
          temphit.width = (*qtc_tr_tdc)[j] - (*qtc_le_tdc)[i];
          temphit.qtc_le =(*qtc_le_tdc)[i];
          temphit.qtc_tr = (*qtc_tr_tdc)[j];
          unsigned int address = (2 << 24) + temphit.channel;
          temphit.SetAddress(address);
          temphit.SetCharge(temphit.width);
          labr_hits.push_back(temphit);
        }
      }
    }

  }

  if (adc) {
    madc1 = TMath::Sqrt((*adc)[0]*(*adc)[1]);
    madc2 = TMath::Sqrt((*adc)[2]*(*adc)[3]);
    madc3 = TMath::Sqrt((*adc)[4]*(*adc)[5]);
  }
  if (tdc) {
    tpos1 = TMath::Sqrt((*tdc)[0]*(*tdc)[1]);
    tpos2 = TMath::Sqrt((*tdc)[2]*(*tdc)[3]);
    tpos3 = TMath::Sqrt((*tdc)[4]*(*tdc)[5]);
  }

#endif
}

void TGrandRaidenHit::SetRaytraceParams(std::vector<double> apar, std::vector<double> bpar, size_t xdeg, size_t adeg, size_t ydeg, size_t bdeg) {
  xdegree=xdeg;
  adegree=adeg;
  ydegree=ydeg;
  bdegree=bdeg;
  acoefs = std::move(apar);
  bcoefs = std::move(bpar);
}

std::pair<double,double> TGrandRaidenHit::Raytrace(bool apply_offsets) {
  return raytrace(rcnp.GR_X(0),rcnp.GR_TH(0),rcnp.GR_Y(0),rcnp.GR_PH(0),apply_offsets);
}
TVector3 TGrandRaidenHit::GetEjectileVector(bool apply_offsets) {

  double thetax=0,thetay=0; // A,B
  std::tie(thetax,thetay) = raytrace(rcnp.GR_X(0),rcnp.GR_TH(0),rcnp.GR_Y(0),rcnp.GR_PH(0),apply_offsets);
  thetax/=1000; // convert to radian from mrad
  thetay/=1000; // convert to radian from mrad
  auto phi = TMath::ATan2(TMath::Sin(thetay),TMath::Sin(thetax));
  auto theta = TMath::ATan(TMath::Sqrt( TMath::Power(TMath::Tan(thetax),2) + TMath::Power(TMath::Tan(thetay),2) ));

  vector.SetTheta(theta);
  vector.SetPhi(phi);
  vector.SetMag(1); // probably unnecessary, but ROOT...

  return vector;
}


TVector3 TGrandRaidenHit::ReconstructInvariant(const TVector3& gamma, bool apply_offsets) {
  auto ejectile = GetEjectileVector(apply_offsets);
  ejectile.SetMag(GetMomentum());
  return ejectile + gamma;
}

double TGrandRaidenHit::GetMomentum() {
  if (momentum) { return momentum; }
  momentum = GValue::Value("GrandRaiden_Slope")*rcnp.GR_X(0) + GValue::Value("GrandRaiden_Offset");
  return momentum;
}

std::pair<double,double> TGrandRaidenHit::raytrace(double x, double a, double y, double b, bool apply_offset) {
  double sum = 0;
  double count = 0;
  double A = 0;
  double B = 0;
  // dispersive angle raytrace
  for (auto i=0u; i<=xdegree; i++) {
    sum += acoefs[i]*pow(x,i);
    count = i;
  }
  for (auto i=1u; i<=adegree; i++) {
    sum += acoefs[i+count]*pow(a,i);
  }
  A=sum;
  sum=count=0;

  // hole arrangement is reversed from sieve slit points
  A *= -1;

  // non-dispersive angle raytrace
  for (auto i=0u; i<= xdegree; i++) {
    double sum2 = 0;
    for (auto j=0u; j<= adegree; j++) {
      double sum3 = 0;
      for (auto k=0u; k<= ydegree; k++) {
        double sum4 = 0;
        for (auto l=0u; l<= bdegree; l++) {
          sum4 += bcoefs[count]*pow(b,l);
          count++;
        }
        sum3 += sum4*pow(y,k);
      }
      sum2 += sum3*pow(a,j);
    }
    sum += sum2*pow(x,i);
  }
  B=sum;


  double a_offset = 0, b_offset = 0;
  std::tie(a_offset,b_offset) = TGrandRaidenHit::GetAngleOffsets();
  return std::pair<double,double>(A-a_offset*1000,B-b_offset*1000);
}

void TGrandRaidenHit::Copy(TObject& obj) const {
  TDetectorHit::Copy(obj);



}


void TGrandRaidenHit::Print(Option_t *opt) const { }

void TGrandRaidenHit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
}













