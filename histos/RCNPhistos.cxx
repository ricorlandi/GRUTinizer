
#include "TRuntimeObjects.h"

#include <iostream>
#include <chrono>
#include <thread>
#include <sstream>
#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TObject.h"
#include "TFile.h"
#include "TCutG.h"
#include "TPreserveGDirectory.h"
#include "TCagra.h"
#include "TGrandRaiden.h"

#define BAD_NUM -441441

#include "TANLEvent.h"
//#include "GValue.h"

#define PRINT(x) std::cout << #x" = " << x << std::endl
#define STR(x) #x << " = " <<() x

using namespace std;


static string name;
static stringstream stream;
static TCutG* aGraphicalCut = nullptr;

void MakeCAGRAHistograms(TRuntimeObjects&, TCagra&);
void MakeGrandRaidenHistograms(TRuntimeObjects&, TGrandRaiden&);
void MakeCoincidenceHistograms(TRuntimeObjects&, TCagra&, TGrandRaiden&);
void LoadCuts();

// ----------------------------------------------------------------------
// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  LoadCuts();

  TCagra* cagra = obj.GetDetector<TCagra>();
  TGrandRaiden* gr = obj.GetDetector<TGrandRaiden>();

  if (gr) {
    MakeGrandRaidenHistograms(obj,*gr);
  }
  if (cagra) {
    MakeCAGRAHistograms(obj,*cagra);
  }
  if (cagra && gr) {
    MakeCoincidenceHistograms(obj, *cagra, *gr);
  }

}
// ----------------------------------------------------------------------

void MakeCAGRAHistograms(TRuntimeObjects& obj, TCagra& cagra) {

  for (auto& core_hit : cagra) {

    auto boardid = core_hit.GetBoardID();
    auto chan = core_hit.GetChannel();

    // clover hit, do normal energy analysis
    stream.str("");
    stream << "Ge_" << boardid << "_" << chan;
    obj.FillHistogram("CAGRA_Raw", stream.str(),10000,0,10000,core_hit.Charge());
    for (auto& segment : core_hit) {
      stream.str("");
      stream << "Ge_" << boardid << "_" << segment.GetChannel();
      obj.FillHistogram("CAGRA_Raw",stream.str(),10000,0,10000,segment.Charge());
    }

    stream.str("");
    stream << "Ge_" << boardid << "_" << chan;
    obj.FillHistogram("CAGRA_Calibrated",stream.str(),10000,0,10000,core_hit.GetEnergy());
    for (auto& segment : core_hit) {
      stream.str("");
      stream << "Ge_" << boardid << "_" << segment.GetChannel();
      obj.FillHistogram("CAGRA_Calibrated",stream.str(),10000,0,10000,segment.GetEnergy());
    }

    // skip polezero for LaBr3 system
    if ( core_hit.GetSystem() != 'L') {



      Double_t prerise_base = core_hit.GetPreRise()/TANLEvent::GetShapingTime();

      stream.str("");
      stream << "E_BL" << boardid << "_" << chan;
      obj.FillHistogram("CAGRA_Baseline", stream.str(),1000,0,10000,core_hit.Charge(),1000,0,3000,prerise_base);

      stream.str("");
      stream << "E_BL_scale" << boardid << "_" << chan;
      obj.FillHistogram("CAGRA_Baseline", stream.str(),1000,0,10000,core_hit.Charge()-(1.0/-11.21)*prerise_base,1000,0,3000,prerise_base);

      stream.str("");
      stream << "Ge_PZ_" << boardid << "_" << chan;
      obj.FillHistogram("CAGRA_Corrected", stream.str(),10000,0,10000,core_hit.GetCorrectedEnergy(0));
      stream.str("");
      stream << "Ge_PZ_AsymBL_" << boardid << "_" << chan;
      obj.FillHistogram("CAGRA_Corrected", stream.str(),10000,0,10000,core_hit.GetCorrectedEnergy(core_hit.GetBaseSample()));

      stream.str("");
      stream << "E_cor_BL" << boardid << "_" << chan;
      obj.FillHistogram("CAGRA_Baseline", stream.str(),1000,0,10000,core_hit.GetCorrectedEnergy(core_hit.GetBaseSample()),1000,0,3000,prerise_base);
    }

    static ULong_t first_ts = 0;
    if (first_ts <= 1e6){  first_ts = core_hit.Timestamp(); std::cout << first_ts << std::endl; }
    else {
      obj.FillHistogram("NumEvents","cagra_hits_time",1000,0,8000,(core_hit.Timestamp()-first_ts)*10/1.0e9);
    }



    // Manual trace analysis
    if (false && core_hit.GetDetnum() <= 104) {
      //LaBr3 hit, do custom trace analysis
      stream.str(""); stream << "LaBr3_traces_" <<boardid << "_" << chan;
      auto labr_E = core_hit.GetTraceEnergy(0,57);
      obj.FillHistogram("CAGRA_Raw",stream.str(),10000,0,10000,labr_E);
    }




  }


}
void MakeGrandRaidenHistograms(TRuntimeObjects& obj, TGrandRaiden& gr) {

  std::function<void(std::string)> fp_corrections;

  for (auto& hit : gr) {

    auto& rcnp = hit.GR();

    if (rcnp.GR_MYRIAD(0) != BAD_NUM) {
      obj.FillHistogram("Timing","MyriadTimestamp",10000,1e9,5e12,hit.GetTimestamp());
    }

    static ULong_t prev_ts = 0;
    if (prev_ts) {
      obj.FillHistogram("Timing","GR_EventPeriod",5000,100,50000,hit.GetTimestamp()-prev_ts);
    }
    prev_ts = hit.GetTimestamp();

    auto rf = rcnp.GR_RF(0);
    if (rf != BAD_NUM) {
      obj.FillHistogram("GR","GR_RF",1000,0,0,rf);
    }

    obj.FillHistogram("GR","RayID",64,-16,48, rcnp.GR_RAYID(0));
    if (rcnp.GR_RAYID(0) == 0) { // if track reconstruction successfull
      obj.FillHistogram("GR","GR_X",1200,-600,600, rcnp.GR_X(0));
      obj.FillHistogram("GR","GR_X_cal",1000,0,20, rcnp.GR_X(0)*0.0109+7.6324);
      obj.FillHistogram("GR","GR_Y",200,-100,100, rcnp.GR_Y(0));
      obj.FillHistogram("GR","GR_Theta",100,-1,1, rcnp.GR_TH(0)); // need to learn
      obj.FillHistogram("GR","GR_Phi",100,-1,1, rcnp.GR_PH(0)); // from hist.def
      obj.FillHistogram("GR","X_TH",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
    }

    if (rcnp.GR_ADC()) {
      auto& adc = *rcnp.GR_ADC();
      for (int i=0; i<4; i++) {
        stream.str(""); stream << "GR_ADC" << i;
        obj.FillHistogram("GR",stream.str().c_str(), 1000,0,2000, adc[i]);
      }
      obj.FillHistogram("GR","MeanPlastE1", 2000,0,2000, hit.GetMeanPlastE1());
      obj.FillHistogram("GR","MeanPlastE2", 2000,0,2000, hit.GetMeanPlastE2());
      if (rf != BAD_NUM && rcnp.GR_RAYID(0)==0) {
        obj.FillHistogram("GR","dE1[RF]",1000,0,0,rf,2000,0,2000, hit.GetMeanPlastE1());
        obj.FillHistogram("GR","dE2[RF]",1000,0,0,rf,2000,0,2000, hit.GetMeanPlastE2());
        obj.FillHistogram("GR","dE1[dE2]",2000,0,2000, hit.GetMeanPlastE2(),2000,0,2000, hit.GetMeanPlastE1());
        obj.FillHistogram("GR","Y[X]",1200,-600,600,rcnp.GR_X(0),200,-100,100,rcnp.GR_Y(0));
        obj.FillHistogram("GR","RF[X]",1200,-600,600, rcnp.GR_X(0),500,700,1200,rf);
        obj.FillHistogram("GR","RF[TH]",1000,-1,1, rcnp.GR_TH(0),500,700,1200,rf);
        obj.FillHistogram("GR","dE1[X]",1200,-600,600, rcnp.GR_X(0),2000,0,2000, hit.GetMeanPlastE1());
        obj.FillHistogram("GR","dE1[TH]",1000,-1,1, rcnp.GR_TH(0),2000,0,2000, hit.GetMeanPlastE1());
        obj.FillHistogram("GR","dE2[X]",1200,-600,600, rcnp.GR_X(0),2000,0,2000, hit.GetMeanPlastE2());
        obj.FillHistogram("GR","dE2[TH]",1000,-1,1, rcnp.GR_TH(0),2000,0,2000, hit.GetMeanPlastE2());

        // e.g. use of TCutG
        //if  (cut0->IsInside(rcnp.GR_X(0),rcnp.GR_TH(0))) { ; }

      }
    }
    if (rcnp.GR_TDC()) {
      auto& tdc = *rcnp.GR_TDC();
      for (int i=0; i<4; i++) {
        stream.str(""); stream << "GR_TDC" << i;
        obj.FillHistogram("GR",stream.str(), 1000,-40000,40000, tdc[i]);
      }
      obj.FillHistogram("GR","MeanPlastPos1", 1000, 0, 40000, hit.GetMeanPlastPos1());
      obj.FillHistogram("GR","MeanPlastPos2", 1000, 0, 40000, hit.GetMeanPlastPos2());
    }


    for (auto const& labr_hit : hit.GetLaBr()) {

      int channum = labr_hit.channel;
      stream.str(""); stream << "LaBrLeading" << channum;
      obj.FillHistogram("GR",stream.str(), 10000,-40000, 40000, labr_hit.qtc_le);


      stream.str(""); stream << "LaBr" << channum << "_LE[LaBr_width]";
      obj.FillHistogram("GR",stream.str(), 1000, -5000, 15000,labr_hit.width, 1000,-40000, 40000, labr_hit.qtc_le);


      stream.str(""); stream << "LaBrWidth" << channum;
      obj.FillHistogram("GR",stream.str(), 10000, -5000, 15000, labr_hit.width);

      if (rcnp.GR_X(0) != BAD_NUM) {
        obj.FillHistogram("GR","X_LaBr",1200,-600,600,rcnp.GR_X(0),10000,-5000,15000,labr_hit.width);
      }

      /* cut on prompt timing peak */ // Needs to be updated for campaign !!!!!!!!!!!
      if (false && labr_hit.qtc_le>=-5100 && labr_hit.qtc_le <=-4300) {


        obj.FillHistogram("GR_Prompt","RayID",64,-16,48, rcnp.GR_RAYID(0));
        if (rcnp.GR_RAYID(0) == 0) { // if track reconstruction successfull
          obj.FillHistogram("GR_Prompt","GR_X",1200,-600,600, rcnp.GR_X(0));
          obj.FillHistogram("GR_Prompt","GR_Y",200,-100,100, rcnp.GR_Y(0));
          obj.FillHistogram("GR_Prompt","GR_Theta",100,-1,1, rcnp.GR_TH(0));
          obj.FillHistogram("GR_Prompt","GR_Phi",100,-1,1, rcnp.GR_PH(0));
          obj.FillHistogram("GR_Prompt","X_TH",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));

          obj.FillHistogram("GR_Prompt","GR_Theta_Phi",100,-1,1, rcnp.GR_TH(0),100,-1,1, rcnp.GR_PH(0));
          obj.FillHistogram("GR_Prompt","GR_X_Y",1200,-600,600, rcnp.GR_X(0),200,-100,100, rcnp.GR_Y(0));

        }

        stream.str(""); stream << "LaBrWidth_LEGate" << channum;
        obj.FillHistogram("GR_Prompt",stream.str(), 10000, -5000, 15000, labr_hit.width);

        stream.str(""); stream << "X_LaBrWidth" << channum;
        obj.FillHistogram("GR_Prompt",stream.str(),1200,-600,600,rcnp.GR_X(0),10000,-5000,15000,labr_hit.width);

      }



    }


  }

}
void MakeCoincidenceHistograms(TRuntimeObjects& obj, TCagra& cagra, TGrandRaiden& gr) {

  for (auto& hit : gr) {

    auto& rcnp = hit.GR();
    auto grtime = hit.GetTimestamp();

    auto x = rcnp.GR_X(0);
    auto Ex = x*0.0109738+7.65621;

    // coincidence rate
    static ULong_t first_timestamp = grtime;
    if (first_timestamp) {
      auto rate = (grtime-first_timestamp)/1e8;
      //cout << grtime << " " << first_timestamp << endl;
      obj.FillHistogram("COIN","Rate",3000,0,30000, rate);
    }

    // coincidence time difference
    for (auto& core_hit : cagra) {

      auto cagratime = core_hit.Timestamp();
      auto tdiff = cagratime-grtime;
      auto boardid = core_hit.GetBoardID();
      auto channel = core_hit.GetChannel();

      obj.FillHistogram("COIN","Diff",1000,-1000,1000,cagratime-grtime);
      stream.str("");
      stream << "TimeDiff_" << boardid << "_" <<channel;
      obj.FillHistogram("COIN",stream.str().c_str(),1200,-300,800,tdiff);
      if ( (tdiff > 235) && (tdiff < 245) ){

        stream.str("");
        stream << "Ge_"<< boardid << "_" <<channel;
        obj.FillHistogram("COIN_Raw",stream.str(),10000,0,10000,core_hit.Charge());
        stream.str("");
        stream << "Ge_" << boardid << "_" <<channel;
        obj.FillHistogram("COIN_Calibrated",stream.str(),10000,0,10000,core_hit.GetEnergy());
        for (auto& segment : core_hit) {
          stream.str("");
          stream << "Ge_" << boardid << "_" << segment.GetChannel();
          obj.FillHistogram("COIN_Calibrated",stream.str(),10000,0,10000,segment.GetEnergy());
        }

        stream.str("");
        stream << "Ex_Ge_" << boardid << "_" <<channel;
        obj.FillHistogram("COIN_Calibrated",stream.str(),1000,0,20,Ex,2500,0,10,core_hit.GetEnergy()*0.0010552+0.0636885);


        stream << "Ge_PZ_AsymBL_" << boardid << "_" << channel;
        obj.FillHistogram("COIN_Calibrated", stream.str().c_str(),10000,0,10000,core_hit.GetCorrectedEnergy(core_hit.GetBaseSample()));

        stream.str("");
        stream << "Ex_CAGRACorrected_" << boardid << "_" <<channel;
        obj.FillHistogram("COIN_Calibrated", stream.str().c_str(),1000,0,20,Ex,2500,0,10000,core_hit.GetCorrectedEnergy(core_hit.GetBaseSample()));
      }

    }


  }



}


void LoadCuts() {


  if(!aGraphicalCut) {
    /*

      TPreserveGDirectory Preserve;
      TFile fcut("./path/to/cut.root");
      aGraphicalCut = (TCutG*)fcut.Get("_cut0");
      std::cout << "Loaded ___ gate." << std::endl;

    */
  }


  // etc

}
