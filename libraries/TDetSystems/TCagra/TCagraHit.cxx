#include "TCagraHit.h"

#include "TCagra.h"

#include <algorithm>
#include <iostream>

#include "TString.h"
#include "TRandom.h"
#include "TANLEvent.h"

#include "GCanvas.h"
#include "GValue.h"

#include "TGRUTOptions.h"
#include "TF1.h"
#include "TFitResult.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/RootFinder.h"


ClassImp(TCagraHit)

TCagraHit::TCagraHit() :  charge(0.0), prerise_energy(0), postrise_energy(0), fit_params({}) {
  is_segment = false;
  baseline_fitted=0;
  time=0;
  prev_time=0;
  flags=0;
  prerise_energy=0;
  postrise_energy=0;
  base_sample=0;
  sampled_baseline=0;
  prev_postrise_begin_sample=0;
  prerise_begin=0;
  prerise_end=0;
  corrected_energy=0;
}

TCagraHit::~TCagraHit() {
}
void TCagraHit::Copy(TObject& obj) const {
  TDetectorHit::Copy(obj);
}
void TCagraHit::Print(Option_t *opt) const {
}
void TCagraHit::PrintChannel() const {
  TChannel::GetChannel(fAddress)->Print();
}

void TCagraHit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
  fTrace.clear();
}
bool TCagraHit::HasCore() const {
  return fCharge != -1;
}

int TCagraHit::GetDetnum() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  int output = -1;
  if(chan && fAddress!=-1){
    output = chan->GetArrayPosition();
  } else if(fSegments.size()) {
    output = fSegments[0].GetDetnum();
  } else {
    // std::cout << "Unknown address: " << std::hex << fAddress << std::dec
    //           << std::endl;
    output = -1;
  }

  if(output == -1 && chan){
    // std::cout << "Chan with det=-1: " << chan->GetName() << std::endl;
    // std::cout << "address: " << fAddress << std::endl;
  }

  return output;
}
char TCagraHit::GetLeaf() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  char output = (char)-1;
  if(chan && fAddress!=-1){
    output = *chan->GetArraySubposition();
  } else if(fSegments.size()) {
    output = fSegments[0].GetLeaf();
  } else {
    // std::cout << "Unknown address: " << std::hex << fAddress << std::dec
    //           << std::endl;
    output = (char)-1;
  }

  if(output == -1 && chan){
    // std::cout << "Chan with det=-1: " << chan->GetName() << std::endl;
    // std::cout << "address: " << fAddress << std::endl;
  }

  return output;
}
int TCagraHit::GetSegnum() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return chan->GetSegment();
  } else {
    return -1;
  }
}

char TCagraHit::GetSystem() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return *chan->GetSystem();
  } else {
    return -1;
  }
}

// int TCagraHit::GetCrate() const {
//   return (fAddress&0x00ff0000)>>16;
// }

int TCagraHit::GetBoardID() const {
  return (fAddress&0x0000ff00)>>8;
}

int TCagraHit::GetChannel() const {
  return (fAddress&0x000000ff)>>0;
}

TCagraHit& TCagraHit::MakeSegmentByAddress(unsigned int address){
  // for(auto& segment : fSegments){
  //   if(segment.Address() == address){
  //     return segment;
  //   }
  // }

  fSegments.emplace_back();
  TCagraHit& output = fSegments.back();
  output.SetAddress(address);
  return output;
}

int TCagraHit::GetMainSegnum() const {
  int output = 0;
  double max_energy = -9e99;
  for(auto& segment : fSegments){
    if(segment.GetEnergy() > max_energy){
      output = segment.GetSegnum();
      max_energy = segment.GetEnergy();
    }
  }
  return output;
}

TVector3 TCagraHit::GetMomentumVector(pos opt, bool apply_array_offset){
  auto gamma = GetPosition(opt,apply_array_offset);
  gamma.SetMag(GetCorrectedEnergy()/1000);
  return gamma;
}

TVector3 TCagraHit::GetPosition(pos opt, bool apply_array_offset) const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  auto clover_type = *chan->GetSystem();
  UShort_t slot = chan->GetArrayPosition(); // slot number
  char core = *chan->GetArraySubposition(); // leaf number
  UShort_t seg = 0;
  // seg = 0 gives the position of crystal center
  // seg = 1 gives the more backward segment
  // seg = 2 gives the more forward segment

  if (chan->GetName()[0] == 'B') { return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1)); }

  if (opt == pos::both || opt == pos::seg_only) {
    if (clover_type == 'Y') {

      TCagraHit const* seg_hit = nullptr;

      size_t num_hits = GetNumSegments();
      switch(num_hits) {
      case 0:
        seg = 0;
        break;
      case 1:
        seg_hit = &fSegments[0];
        break;
      case 2:
        seg_hit =
          (fSegments[0].GetEnergy() > fSegments[1].GetEnergy()) ?
          &fSegments[0] : &fSegments[1];
        break;
      case 3:
      default:
        std::cout << "Slot: " << slot << " Core: " << core << " NumHits: " << num_hits << std::endl;
        std::cout << "Core hit Energy:" << GetEnergy() << " Timestamp: " << Timestamp() << std::endl;
        std::cout << "Segment hits: \n";
        for (auto i=0u; i<fSegments.size(); i++) {
          std::cout << fSegments[i].GetDetnum() << " " << fSegments[i].GetLeaf() <<
            " Timestamp: " << fSegments[i].Timestamp() << " Energy: " << fSegments[i].GetEnergy() << std::endl;
        }
        throw std::runtime_error("More than two segment hits in a single Yale crystal.");
        break;
      }

      if (seg_hit) {
        char leaf = seg_hit->GetLeaf();

        // southern hemisphere
        if (slot <=8 && slot >=5) {
          if (core == 'A' || core == 'D') {
            assert(leaf != 'R'); // comment this after checking
            seg = (leaf == 'L') ? 1 : 2;
          }
          else if (core == 'B' || core == 'C') {
            assert(leaf != 'L');
            seg = (leaf == 'M') ? 1 : 2;
          }
        }
        // northern hemisphere
        else {
          if (core == 'A' || core == 'D') {
            assert(leaf != 'R'); // comment this after checking
            seg = (leaf == 'L') ? 2 : 1;
          }
          else if (core == 'B' || core == 'C') {
            assert(leaf != 'L');
            seg = (leaf == 'M') ? 2 : 1;
          }
        }
      }

    } else if (clover_type == 'I') {
      seg = chan->GetSegment();
    } else {
      std::cout << "Segments from unexpected clover type: " << clover_type << std::endl;
      throw std::runtime_error("Segments from unexpected detector system.");
    }

  } // end if pos is seg or both

  if (opt == pos::seg_only && seg == 0) {
    return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }

  TVector3 array_pos = TCagra::GetSegmentPosition(slot, core, seg);
  if(apply_array_offset){
    array_pos += TVector3(GValue::Value("Cagra_X_offset"),
                          GValue::Value("Cagra_Y_offset"),
                          GValue::Value("Cagra_Z_offset"));
  }
  return array_pos;
}

double TCagraHit::GetDoppler(double beta, pos opt, const TVector3& particle_vec, const TVector3& offset) {

  double gamma = 1/(sqrt(1-pow(beta,2)));
  TVector3 pos = GetPosition(opt) + offset;
  double cos_angle = TMath::Cos(pos.Angle(particle_vec));
  double dc_en = GetCorrectedEnergy()*gamma *(1 - beta*cos_angle);
  return dc_en;
}

Int_t TCagraHit::Charge() const {
  if(fCharge > 30000) {
    return fCharge - 32768;
  } else {
    return fCharge;
  }
}

Float_t TCagraHit::GetCharge() const {
  return fCharge;
}


Double_t TCagraHit::GetCorrectedEnergy(Double_t asym_bl) {
  if (!asym_bl && corrected_energy) { return corrected_energy; }

  TChannel* chan = TChannel::GetChannel(fAddress);
  Double_t Energy = 0;
  if(!chan){
    static int count = 0;
    if (count++ < 10) {
      std::cout << std::hex << "Channel 0x" << fAddress << " not defined in calibrations file, no corrections are applied." << std::dec << std::endl;
    } else if (count == 10) {
      std::cout << "Supressing warning." <<std::endl;
    }
  } else {
    int polarity = TANLEvent::GetSignalPolarity();
    // if (is_segment) {
    //   polarity *= -1;
    // }

    auto pzE = chan->PoleZeroCorrection(prerise_energy,postrise_energy,TANLEvent::GetShapingTime(),polarity);
    pzE = chan->BaselineCorrection(pzE,asym_bl,polarity);
    Energy = chan->CalEnergy(pzE, fTimestamp);
  }

  if (!asym_bl) { corrected_energy = Energy; }
  return Energy;
}

Double_t TCagraHit::GetTraceHeightPZ(Double_t asym_bl, size_t size) {
  if(fTrace.size() < 2*size){
    return std::sqrt(-1);
  }

  TChannel* chan = TChannel::GetChannel(fAddress);
  Double_t Energy = 0;
  if(!chan){
    static int count = 0;
    if (count++ < 10) {
      std::cout << std::hex << "Channel 0x" << fAddress << " not defined in calibrations file, no corrections are applied." << std::dec << std::endl;
    } else if (count == 10) {
      std::cout << "Supressing warning." <<std::endl;
    }
  } else {

    double low = 0;
    double high = 0;
    GetSanitizedTrace();
    for(unsigned int i=0; i<size; i++){
      low += fTrace[i];
      high += fTrace[fTrace.size()-i-1];
    }

    auto pzE = chan->PoleZeroCorrection(low,high,size,TANLEvent::GetSignalPolarity());

    pzE = (asym_bl) ? chan->BaselineCorrection(pzE,asym_bl) : chan->BaselineCorrection(pzE);
    Energy = chan->CalEnergy(pzE, fTimestamp);
  }
  return Energy;
}

Double_t TCagraHit:: GetTraceBaseline() {
  std::vector<Short_t>* trace = GetTrace();
  if (trace) {
    int length = std::min<int>(6,trace->size());
    Double_t bl = 0;
    for (int i=0; i<length; i++) {
      bl+=transform_trace_point((*trace)[i]);
    }
    return bl/length;
  }
  return 0;
}

Double_t ExpBaseline(Double_t* x, Double_t* par) {
  return par[0]*TMath::Exp(-x[0]*par[1]);
}

Double_t Line(Double_t* x, Double_t* par) {
  return par[0]*x[0] + par[1];
}


void TCagraHit::DrawBaselineExponential(int segnum) {

  Int_t time[3];
  Double_t sample[3];
  Double_t baseline = 1;//transform_trace_point(base_sample);

  sample[0] = transform_trace_point(prev_postrise_begin_sample)/baseline;
  sample[1] = transform_trace_point(prerise_begin)/baseline;
  sample[2] = transform_trace_point(prerise_end)/baseline;

  Double_t tprev_postrise_enter = prev_time + GValue::Value("d_window")*100 + GValue::Value("k_window")*100;

  Double_t tdiff = Timestamp() - prev_time;
  time[0] = 0;
  time[2] = Timestamp()
    - GValue::Value("d2_window")*100
    - GValue::Value("d3_window")*100
    - tprev_postrise_enter;

  time[1] = time[2] - TANLEvent::GetShapingTime();


  TH1F hist("hist", "", 2.5*int(tdiff), -tdiff/2, 2*tdiff);
  hist.SetStats(false);


  hist.SetTitle(Form("CAGRA Slot %d at %ld ns - Baseline", GetDetnum(), Timestamp()));
  hist.GetXaxis()->SetTitle("Time (clock ticks)");
  hist.GetYaxis()->SetTitle("ADC units");
  for(size_t i=0; i<3; i++) {
    hist.SetBinContent(time[i]+int(tdiff/2),sample[i]);
    //std::cout << time[i] << " " << sample[i] << std::endl;
  }

  Double_t p[2] = {5.0,5.0e-05};
  TF1* func = new TF1("ExpBaseline",ExpBaseline,-100,100,2);

  // To ignore error messages about Minuit2 not existing
  gROOT->ProcessLine( "gErrorIgnoreLevel = 3002;");

  //Double_t perr[2];
  //while ( perr[0] < 10 && perr[1] < 10) {
    //p[0]
  func->SetParameters(p);
  hist.Fit(func,"Q");

  auto baseline_fit = hist.GetFunction("ExpBaseline");
//  delete func;
  //std::cout << baseline_fit->Eval(0) << " " << sample[0] << std::endl;
  //std::cout << baseline_fit->GetParError(0) << " " << baseline_fit->GetParError(1) << std::endl;
  //perr[0] = baseline_fit->GetParError(0);
  //perr[1] = baseline_fit->GetParError(1);

    //break;
    //}
  gROOT->ProcessLine( "gErrorIgnoreLevel = -1;");



  hist.DrawCopy();


}


Double_t TCagraHit::GetBaselineExpCorr(int segnum) {

  Int_t time[3];
  //if (!baseline_fit) {

    Double_t sample[3];
    Double_t baseline = 1;//transform_trace_point(base_sample);

    sample[0] = transform_trace_point(prev_postrise_begin_sample)/baseline;
    sample[1] = transform_trace_point(prerise_begin)/baseline;
    sample[2] = transform_trace_point(prerise_end)/baseline;

    Double_t tprev_postrise_enter = prev_time + GValue::Value("d_window")*100 + GValue::Value("k_window")*100;
    Double_t tdiff = Timestamp() - prev_time;
    time[0] = 0;
    time[2] = Timestamp()
      - GValue::Value("d2_window")*100
      - GValue::Value("d3_window")*100
      - tprev_postrise_enter;

    time[1] = time[2] - TANLEvent::GetShapingTime();


    TH1F hist("hist", "", 2.5*int(tdiff), -tdiff/2, 2*tdiff);
    hist.SetStats(false);


    hist.SetTitle(Form("CAGRA Slot %d at %ld ns - Baseline", GetDetnum(), Timestamp()));
    hist.GetXaxis()->SetTitle("Time (clock ticks)");
    hist.GetYaxis()->SetTitle("ADC units");
    for(size_t i=0; i<3; i++) {
      hist.SetBinContent(time[i]+int(tdiff/2),sample[i]);
      //std::cout << time[i] << " " << sample[i] << std::endl;
    }

    Double_t p[2] = {5.0,5.0e-05};
    TF1* func = new TF1("ExpBaseline",ExpBaseline,-100,100,2);

    // To ignore error messages about Minuit2 not existing
    gROOT->ProcessLine( "gErrorIgnoreLevel = 3002;");

    //Double_t perr[2];
    //while ( perr[0] < 10 && perr[1] < 10) {
    //p[0]

    func->SetParameters(p);
    hist.Fit(func,"0Q");

    auto baseline_fit = hist.GetFunction("ExpBaseline");
    //delete func;
    //std::cout << baseline_fit->Eval(0) << " " << sample[0] << std::endl;
    //std::cout << baseline_fit->GetParError(0) << " " << baseline_fit->GetParError(1) << std::endl;
    //perr[0] = baseline_fit->GetParError(0);
    //perr[1] = baseline_fit->GetParError(1);

    //break;
    //}
    gROOT->ProcessLine( "gErrorIgnoreLevel = -1;");



  // integrate pre-rise region
  auto prerise_corr = baseline_fit->Integral(time[1],time[2])/(time[2]-time[1]);
  auto energy = (postrise_energy - (prerise_energy-prerise_corr))/TANLEvent::GetShapingTime();
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(!chan){
    static int count = 0;
    if (count++ < 10) {
      std::cout << std::hex << "Channel 0x" << fAddress << " not defined in calibrations file, no corrections are applied." << std::dec << std::endl;
    } else if (count == 10) {
      std::cout << "Supressing warning." <<std::endl;
    }
  } else {
    energy = chan->CalEnergy(energy, fTimestamp);
  }

  // if (baseline_fit->GetParError(0) < 3) {
  //   std::cout << "Bad fit" << std::endl;
  //   std::cout <<baseline_fit->GetParError(0) << " " << baseline_fit->GetParError(1) << std::endl;
  //   return GetEnergy();
  // }


  //delete baseline_fit;
  return energy;

  // std::cout << inverse_transform_trace_point(prerise_corr) << std::endl;
  // std::cout << base_sample << std::endl;
  // std::cout << prerise_energy << " " << postrise_energy << std::endl;

  // Double_t begin = time[1];
  // Double_t end = time[2];
  // Double_t par[2] = {baseline_fit->GetParameter(0), baseline_fit->GetParameter(1)};
  // Double_t integral = -par[0]/par[1]*( ExpBaseline(&begin, par) - ExpBaseline(&end, par) );

  //hist.DrawCopy();


}

///////////////////////////////
double LeastSquaresChi2(double ni, double nui) {return TMath::Power(ni-nui,2)/nui;}
double ModifiedLeastSquaresChi2(double ni, double nui) {return TMath::Power(ni-nui,2)/ni;}
Double_t baseline_exp(Double_t x, const Double_t* par) {
  return par[0]*TMath::Exp(-x*par[1]);
}
std::function<Double_t(const Double_t* params)> wrapper;
Double_t BaselineFit(const Double_t* params) {
  return wrapper(params);
}
///////////////////////////////

Double_t TCagraHit::GetBaselineExpCorrFast(int segnum) {
  const size_t size = 3;
  Int_t time[size];

  Double_t sample[size];
  Double_t baseline = 1;//transform_trace_point(base_sample);

  auto shaping_time = TANLEvent::GetShapingTime();

  sample[0] = transform_trace_point(prev_postrise_begin_sample)/baseline;
  sample[1] = transform_trace_point(prerise_begin)/baseline;
  sample[2] = transform_trace_point(prerise_end)/baseline;

  Double_t tprev_postrise_enter = prev_time + GValue::Value("d_window")*100 + GValue::Value("k_window")*100;
  Double_t tdiff = Timestamp() - prev_time;
  time[0] = 0;
  time[2] = Timestamp()
    - GValue::Value("d2_window")*100
    - GValue::Value("d3_window")*100
    - tprev_postrise_enter;

  time[1] = time[2] - shaping_time;

  std:: cout<< time[0] << " " << time[1] << " " << time[2] << std::endl;

  wrapper = [&time,&sample](const Double_t* params) {
    Double_t sum_chi2 = 0.0;
    for (int i=0; i<size; i++) {
      auto sample_fit = baseline_exp(time[i],params);
      sum_chi2 += LeastSquaresChi2(sample[i],sample_fit);
    }
    return sum_chi2;
  };

  auto minuit = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
  ROOT::Math::Functor kernel(&BaselineFit,2);
  minuit->SetMaxFunctionCalls(1e6);
  minuit->SetMaxIterations(1e5);
  minuit->SetTolerance(0.0001);
  minuit->SetFunction(kernel);
  minuit->SetVariable(0,"N",5.0,0.1);
  minuit->SetVariable(1,"lambda",5.0e-5,1.0e-6);
  minuit->Minimize();
  auto params = minuit->X();
  fit_params[0] = params[0];
  fit_params[1] = params[1];

  if (std::isnan(params[0]) || std::isnan(params[1])) {
    //std::cout << "nan in fit" << std::endl <<std::endl;;
    return 0.0;
  }
  // integral
  Double_t prerise_corr = (baseline_exp(time[1],params) - baseline_exp(time[2],params))/params[1];
  prerise_corr = inverse_transform_trace_point(prerise_corr/shaping_time);

  Double_t tpostrise_enter = Timestamp() + GValue::Value("d_window")*100 + GValue::Value("k_window")*100 - tprev_postrise_enter;
  Double_t tpostrise_exit = tpostrise_enter+shaping_time;
  Double_t postrise_corr = (baseline_exp(tpostrise_enter,params) - baseline_exp(tpostrise_exit,params))/params[1];
  postrise_corr = inverse_transform_trace_point(postrise_corr/shaping_time);

  auto energy = std::abs(postrise_corr-postrise_energy/shaping_time) - std::abs(prerise_energy/shaping_time-prerise_corr);
  //std::cout << "E: " << (postrise_energy-prerise_energy)/shaping_time << " Ecorr: "<< energy << std::endl;



  TChannel* chan = TChannel::GetChannel(fAddress);
  if(!chan){
    static int count = 0;
    if (count++ < 10) {
      std::cout << std::hex << "Channel 0x" << fAddress << " not defined in calibrations file, no corrections are applied." << std::dec << std::endl;
    } else if (count == 10) {
      std::cout << "Supressing warning." <<std::endl;
    }
  } else {
    energy = chan->CalEnergy(energy, fTimestamp);
  }

  return energy;
}

void TCagraHit::DrawTraceBaseline(int segnum) {
  std::vector<Short_t>* trace = GetTrace(segnum);
  if(!trace){
    std::cout << "No segment trace found for segment " << segnum << std::endl;
    return;
  }

  // create base line fit
  GetBaselineExpCorrFast();

  Double_t tprev_postrise_enter = prev_time + GValue::Value("d_window")*100 + GValue::Value("k_window")*100;
  auto time = Timestamp()-tprev_postrise_enter;
  time -= 150; // offset to beginning of trace window
  std::cout << time << std::endl;

  TH1I hist("hist", "", (time+450), 0, 10*(time+450));
  TH1I basehist("baseline", "", (time+450), 0, 10*(time+450));


  hist.SetStats(false);

  if(segnum==0){
    hist.SetTitle(Form("CAGRA Detector %d at %ld ns", GetDetnum(), Timestamp()));
    hist.GetXaxis()->SetTitle("Time (ns)");
    hist.GetYaxis()->SetTitle("ADC units");
  }



  for(size_t i=0; i<(time+450); i++) {

    if (i >= time) {
      auto n=i-time;
      auto adc = (((*trace)[n]) < 0) ? (*trace)[n] + std::pow(2,14) : (*trace)[n];
      hist.SetBinContent(i+1,adc);
    }

    auto fitadc = baseline_exp(i,fit_params);
    fitadc = inverse_transform_trace_point(fitadc);
    basehist.SetBinContent(i+1,fitadc);
  }
  hist.DrawCopy();
  basehist.SetLineWidth(2);
  basehist.SetLineColor(kRed);
  basehist.DrawCopy("same");
}


void TCagraHit::DrawTrace(int segnum, bool draw_baseline) {
  std::vector<Short_t>* trace = GetTrace(segnum);
  if(!trace){
    std::cout << "No segment trace found for segment " << segnum << std::endl;
    return;
  }

  Double_t time = 0.0;
  if (draw_baseline) {
    GetBaselineExpCorrFast();
    Double_t tprev_postrise_enter = prev_time + GValue::Value("d_window")*100 + GValue::Value("k_window")*100;
    time = Timestamp()-tprev_postrise_enter;
    time -= 150; // offset to beginning of trace window
  }
  TH1I hist("hist", "", trace->size(), 0, 10*trace->size());
  TH1I basehist("baseline", "", trace->size(), 0, 10*trace->size());

  hist.SetStats(false);

  if(segnum==0){
    hist.SetTitle(Form("CAGRA Detector %d at %ld ns", GetDetnum(), Timestamp()));
    hist.GetXaxis()->SetTitle("Time (ns)");
    hist.GetYaxis()->SetTitle("ADC units");
  }



  for(size_t i=0; i<trace->size(); i++) {
    auto adc = (((*trace)[i]) < 0) ? (*trace)[i] + std::pow(2,14) : (*trace)[i];
    hist.SetBinContent(i+1,adc);
    if (draw_baseline && fit_params[0]){
      auto fitadc = baseline_exp(time+i,fit_params);
      fitadc = inverse_transform_trace_point(fitadc);
      basehist.SetBinContent(i+1,fitadc);
    }
  }

  hist.DrawCopy();
  if (draw_baseline && fit_params[0]){
    basehist.SetLineColor(kRed);
    basehist.DrawCopy("same");
  }
}

void TCagraHit::DrawTrace(const std::vector<double>& trace, const std::vector<double>& error) {


  Double_t time = 0.0;
  TH1I hist("hist", "", trace.size(), 0, 10*trace.size());

  hist.SetStats(false);

  hist.GetXaxis()->SetTitle("Time (ns)");
  hist.GetYaxis()->SetTitle("ADC units");

  if (error.size()) {
    for(size_t i=0; i<trace.size(); i++) {
      hist.SetBinContent(i+1,trace[i]);
      hist.SetBinError(i+1,error[i]);
    }
  } else {
      for(size_t i=0; i<trace.size(); i++) {
      hist.SetBinContent(i+1,trace[i]);
    }
  }

  hist.DrawCopy();
}

void TCagraHit::DrawTraceSamples(int segnum) {
  std::vector<Short_t> trace(3);
  trace[0] = prev_postrise_begin_sample;
  trace[1] = prerise_begin;
  trace[2] = prerise_end;

  TH1I hist("hist", "", trace.size(), 0, 10*trace.size());
  hist.SetStats(false);

  if(segnum==0){
    hist.SetTitle(Form("CAGRA Detector %d at %ld ns", GetDetnum(), Timestamp()));
    //hist.GetXaxis()->SetTitle("Time (ns)");
    //hist.GetYaxis()->SetTitle("ADC units");
  }

  for(size_t i=0; i<trace.size(); i++) {
    hist.SetBinContent(i+1,trace[i]);
  }
  hist.DrawCopy();
}

void TCagraHit::SetTrace(std::vector<Short_t>& trace) {
  fTrace.clear();
  fTrace.swap(trace);
}

std::vector<Short_t>* TCagraHit::GetTrace(int segnum) {
  if(segnum == 0){
    return &fTrace;
  }
  for(auto& seg : fSegments) {
    if(seg.GetSegnum() == segnum) {
      return seg.GetTrace();
    }
  }
  return NULL;
}

std::vector<Short_t>* TCagraHit::GetSanitizedTrace(int segnum) {
  if(segnum == 0){
    for(size_t i=0; i<fTrace.size(); i++) {
      auto adc = (fTrace[i] < 0) ? fTrace[i] + std::pow(2,14) : fTrace[i];
      fTrace[i] = adc;
    }
    return &fTrace;
  }
  for(auto& seg : fSegments) {
    if(seg.GetSegnum() == segnum) {

      auto segtrace = seg.GetTrace();
      for(size_t i=0; i<segtrace->size(); i++) {
        auto adc = (((*segtrace)[i]) < 0) ? (*segtrace)[i] + std::pow(2,14) : (*segtrace)[i];
        ((*segtrace)[i]) = adc;
      }
      return segtrace;
    }
  }
  return NULL;
}


double TCagraHit::GetTraceHeight(size_t size) const {
  if(fTrace.size() < 2*size){
    return std::sqrt(-1);
  }

  double low = 0;
  double high = 0;
  for(unsigned int i=0; i<size; i++){
    low += fTrace[i];
    high += fTrace[fTrace.size()-i-1];
  }

  return (high-low)/size;
}

double TCagraHit::GetTraceHeightDoppler(double beta,const TVector3& vec) const {
  if(GetNumSegments()<1) {
    return std::sqrt(-1);
  }

  double gamma = 1/(sqrt(1-pow(beta,2)));
  TVector3 pos = GetPosition();
  double cos_angle = TMath::Cos(pos.Angle(vec));
  double dc_en = GetTraceHeight()*gamma *(1 - beta*cos_angle);
  return dc_en;
}

Double_t TCagraHit::GetTraceEnergy(const UShort_t& a,const UShort_t& b,UShort_t x, UShort_t y) const {
  if (fTrace.size() < b) { return 0; }

  if (fTrace.size() < y) {
    static int nprint = 0;
    if (nprint < 10) {
      std::cout << "Warning: Trace length less than requested sampling window: " << fTrace.size() <<std::endl;
    } nprint++;
    return 0;
  }
  if (y==0) { y = fTrace.size(); }
  if ( y > fTrace.size()) {
    std::cout << "Request upper bound on trace integration is larger than the trace size." << std::endl;
    std::cout << "Trace Info: " << a << " " << b << " " << x << " "<<y << " " << fTrace.size() <<std::endl;
    return 0;
  }


  double baseline = 0;
  for (int i=a; i<b; i++) { baseline+=fTrace[i]; }
  double integral = 0;
  for (int i=x; i<y; i++) { integral+=fTrace[i]; }

  return (integral/(y-x) - baseline/(b-a));
}
