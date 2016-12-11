
#include "TGrandRaidenHit.h"
#include "TGRUTOptions.h"
#include "TSmartBuffer.h"
#include "TMath.h"
#include "GValue.h"

ClassImp(TGrandRaidenHit)


std::vector<double> TGrandRaidenHit::acoefs;
std::vector<double> TGrandRaidenHit::bcoefs;
unsigned int TGrandRaidenHit::xdegree = 2, TGrandRaidenHit::adegree = 2, TGrandRaidenHit::ydegree = 1;


TGrandRaidenHit::TGrandRaidenHit() {
  madc1=0; madc2=0; tpos1=0; tpos2=0;
  vector = nullptr;
  excitation_energy=0;
}
TGrandRaidenHit::TGrandRaidenHit(const TGrandRaidenHit& gr) {
  labr_hits = gr.labr_hits;
  madc1 = gr.madc1;
  madc2 = gr.madc2;
  tpos1 = gr.tpos1;
  tpos2 = gr.tpos2;
  Timestamp = gr.Timestamp;
  rcnp = gr.rcnp;
}
TGrandRaidenHit::TGrandRaidenHit(RCNPEvent& rcnpevent) :
  rcnp(rcnpevent) {
  madc1=0; madc2=0; tpos1=0; tpos2=0;
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
          temphit.channel = (*qtc_le_chan)[i];
          temphit.width = (*qtc_tr_tdc)[j] - (*qtc_le_tdc)[i];
          temphit.qtc_le =(*qtc_le_tdc)[i];
          temphit.qtc_tr = (*qtc_tr_tdc)[j];

          labr_hits.push_back(temphit);

        }
      }
    }

  }

  if (adc) {
    madc1 = TMath::Sqrt((*adc)[0]*(*adc)[1]);
    madc2 = TMath::Sqrt((*adc)[2]*(*adc)[3]);
  }
  if (tdc) {
    tpos1 = TMath::Sqrt((*tdc)[0]*(*tdc)[1]);
    tpos2 = TMath::Sqrt((*tdc)[2]*(*tdc)[3]);
  }

#endif
}

void TGrandRaidenHit::SetRaytraceParams(std::vector<double> apar, std::vector<double> bpar, size_t xdeg, size_t adeg, size_t ydeg) {
  xdegree=xdeg;
  adegree=adeg;
  ydegree=ydeg;
  acoefs = std::move(apar);
  bcoefs = std::move(bpar);
}

TVector3 TGrandRaidenHit::GetEjectileVector() {
  if (vector) { return *vector; }

  double thetax=0,thetay=0; // A,B
  std::tie(thetax,thetay) = raytrace(rcnp.GR_X(0),rcnp.GR_TH(0),rcnp.GR_Y(0));
  auto phi = TMath::ATan2(TMath::Sin(thetay),TMath::Sin(thetax));
  auto theta = TMath::ASin(TMath::Sqrt( TMath::Power(TMath::Sin(thetax),2) + TMath::Power(TMath::Sin(thetay),2) ));

  vector = new TVector3(1.,1.,1.);
  vector->SetTheta(theta);
  vector->SetPhi(phi);
  vector->SetMag(1); // probably unnecessary, but ROOT...

  return *vector;
}


double TGrandRaidenHit::GetExEnergy() {
  if (excitation_energy) { return excitation_energy; }
  excitation_energy = GValue::Value("GrandRaiden_Slope")*rcnp.GR_X(0) + GValue::Value("GrandRaiden_Offset");
  return excitation_energy;
}

std::pair<double,double> TGrandRaidenHit::raytrace(double x, double a, double y) {
  double sum = 0;
  double count = 0;
  double A = 0;
  double B = 0;
  // dispersive angle raytrace
  for (auto i=0u; i<=xdegree; i++) {
    sum += acoefs[i]*pow(x,i);
    count = i;
  }
  for (auto i=0u; i<=adegree; i++) {
    sum += acoefs[i+count+1]*pow(a,i);
  }
  A=sum;
  sum=count=0;

  // non-dispersive angle raytrace
  for (auto i=0u; i<= xdegree; i++) {
    double sum2 = 0;
    for (auto j=0u; j<= adegree; j++) {
      double sum3 = 0;
      for (auto k=0u; k<= ydegree; k++) {
        sum3 += bcoefs[count]*pow(y,k);
        count++;
      }
      sum2 += sum3*pow(a,j);
    }
    sum += sum2*pow(x,i);
  }
  B=sum;

  return std::pair<double,double>(A,B);
}

void TGrandRaidenHit::Copy(TObject& obj) const {
  TDetectorHit::Copy(obj);



}


void TGrandRaidenHit::Print(Option_t *opt) const { }

void TGrandRaidenHit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
}













