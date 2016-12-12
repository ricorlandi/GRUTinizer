#include "TRuntimeObjects.h"
#include "e441.h"

// ----------------------------------------------------------------------
// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  LoadCuts();
  LoadRaytraceParams();

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

void PoleZeroHistos(TRuntimeObjects& obj, TCagraHit& core_hit, string local_dirname = "") {
  auto flags = core_hit.GetFlags();
  if (TANLEvent::PileUpFlag(flags) || TANLEvent::PileUpOnlyFlag(flags)) { return; }

  int detector = core_hit.GetDetnum();
  string chan = core_hit.GetLeaf() + std::to_string(core_hit.GetSegnum());

  Double_t prerise = core_hit.GetPreRise()/TANLEvent::GetShapingTime();
  stream.str("");  stream << "Prerise[Q]_" << detector << "_" << chan;
  obj.FillHistogram(local_dirname, stream.str(),2000,0,10000,core_hit.GetCharge(),1250,6000,8500,prerise);
  for (auto& seg_hit : core_hit) {
    string seg_chan = seg_hit.GetLeaf() + std::to_string(seg_hit.GetSegnum());
    stream.str("");  stream << "Prerise[Q]_" << detector << "_" << seg_chan;
    obj.FillHistogram(local_dirname, stream.str(),2000,-10000,0,seg_hit.GetCharge(),1250,6000,8500,prerise);
  }
  // stream.str("");  stream << "Q[Prerise]_" << detector << "_" << chan;
  // obj.FillHistogram(local_dirname, stream.str(),1250,6000,8500,prerise,3000,0,6000,core_hit.GetCharge());
  // for (auto& seg : core_hit) {
  //   string seg_chan = seg_hit.GetLeaf() + std::to_string(seg_hit.GetSegnum());
  //   stream.str("");  stream << "Q[Prerise]Seg_" << detector << "_" << seg_chan;
  //   obj.FillHistogram(local_dirname, stream.str(),1250,6000,8500,prerise,3000,0,6000,seg.GetCharge());
  // }

  stream.str("");  stream << "Prerise[E_pzcor_basesample]_" << detector << "_" << chan;
  obj.FillHistogram(local_dirname, stream.str(),3000,0,6000,core_hit.GetCorrectedEnergy(core_hit.GetBaseSample()),1250,6000,8500,prerise);
  stream.str("");  stream << "Prerise[E_pzcor_constant]_" << detector << "_" << chan;
  obj.FillHistogram(local_dirname, stream.str(),3000,-2000,4000,core_hit.GetCorrectedEnergy(),1250,6000,8500,prerise);



  stream.str(""); stream << "E_pzcor_basesample" << detector << "_" << chan;
  obj.FillHistogram(local_dirname,stream.str(),4000,0,0,core_hit.GetCorrectedEnergy(core_hit.GetBaseSample()));

  stream.str(""); stream << "E_pzcor_constant" << detector << "_" << chan;
  auto pzchan = core_hit.GetCorrectedEnergy();
  obj.FillHistogram(local_dirname,stream.str(),4000,0,0,pzchan);
  for (auto& seg_hit : core_hit) {
    string seg_chan = seg_hit.GetLeaf() + std::to_string(seg_hit.GetSegnum());
    stream.str(""); stream << "E_pzcor_constant" << detector << "_" << seg_chan;
    obj.FillHistogram(local_dirname,stream.str(),4000,0,0,seg_hit.GetCorrectedEnergy());
  }

}

void PileUp (TRuntimeObjects& obj, TCagraHit& core_hit) {
  dirname = "PileUp";

  auto flags = core_hit.GetFlags();
  auto pileup = TANLEvent::PileUpFlag(flags);
  if (pileup) {
    obj.FillHistogram(dirname,"Summary",10,1,10,3);
    PoleZeroHistos(obj,core_hit,dirname);
  } else {
    obj.FillHistogram(dirname,"Summary",10,1,10,7);
    PoleZeroHistos(obj,core_hit,"NoPileUp");
  }
}


void MakeCAGRAHistograms(TRuntimeObjects& obj, TCagra& cagra) {

  for (auto& core_hit : cagra) {

    //PileUp(obj,core_hit);

    int detector = core_hit.GetDetnum();
    char core_leaf = core_hit.GetLeaf();
    string chan = core_leaf + std::to_string(core_hit.GetSegnum());
    size_t crystal_id = TCagra::GetCrystalId(detector,core_leaf);

    // cagra core energy summary
    obj.FillHistogram("CrystalEnergySummary",
                      6000,-2000,22000,core_hit.GetCorrectedEnergy(),
                      49,0,49,crystal_id);
    obj.FillHistogram("EgammaSum",7000,0,21000,core_hit.GetCorrectedEnergy());


    static ULong_t first_ts = 0;
    if (first_ts <= 1e6){  first_ts = core_hit.Timestamp(); std::cout << "Timestamp: " << first_ts << "\n" << std::endl; }
    else {
      // cagra core time summary
      obj.FillHistogram("CrystalTimeSummary",
                        1000,0,4000,(core_hit.Timestamp()-first_ts)*10/1.0e9, // in seconds - 1 bin = 4 seconds
                        49,0,49,crystal_id);

      obj.FillHistogram("NumEvents","cagra_hits_time",1000,0,8000,(core_hit.Timestamp()-first_ts)*10/1.0e9);
    }


    auto position = core_hit.GetPosition(pos::core_only);
    obj.FillHistogram("ArrayHits",
                      180,0,180,position.Theta()*180/TMath::Pi(),
                      360,-180,180,position.Phi()*180/TMath::Pi());


    PoleZeroHistos(obj,core_hit,"PoleZero");


    // central contact signals
    // name = "Det_" + std::to_string(detector) + "_" + chan;
    // obj.FillHistogram("CAGRA_Raw", name,2000,0,0,core_hit.GetCharge());

    // // segment (side channel) signals
    // for (auto& segment : core_hit) {
    //   string seg_chan = segment.GetLeaf() + std::to_string(segment.GetSegnum());
    //   stream.str(""); stream << "Det_" << detector << "_" << seg_chan;
    //   obj.FillHistogram("CAGRA_Raw",stream.str(),2000,0,0,segment.GetCharge());
    // }


    // same but for calibrated energies
    // stream.str("");
    // stream << "Det_" << detector << "_" << chan;
    // obj.FillHistogram("CAGRA_Calibrated",stream.str(),2000,0,10000,core_hit.GetEnergy());
    // for (auto& segment : core_hit) {
    //   string seg_chan = segment.GetLeaf() + std::to_string(segment.GetSegnum());
    //   stream.str(""); stream << "Det_" << detector << "_" << seg_chan;
    //   obj.FillHistogram("CAGRA_Calibrated",stream.str(),2000,0,10000,segment.GetEnergy());
    // }


  } // end loop over cagra hits

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


    // X, A, Y, B, RF, DE1, DE2
    // X[A],X[Y],X[B],X[RF],DE1[X],DE2[X]
    // Y[A],Y[B],Y[RF],DE1[Y],DE2[Y]
    // A[B], A[RF], DE1[A], DE2[A]
    // B[RF], DE1[B], DE2[B]
    // DE1[RF], DE2[RF], DE1[DE2]

    obj.FillHistogram("GR","RayID",64,-16,48, rcnp.GR_RAYID(0));
    if (rcnp.GR_RAYID(0) == 0) { // if track reconstruction successful
      obj.FillHistogram("GR","X",1200,-600,600, rcnp.GR_X(0));
      obj.FillHistogram("GR","X_cal",1000,0,20, rcnp.GR_X(0)*0.01074+6.872);
      obj.FillHistogram("GR","Y",200,-100,100, rcnp.GR_Y(0));
      obj.FillHistogram("GR","A",100,-1,1, rcnp.GR_TH(0)); // need to learn
      obj.FillHistogram("GR","B",100,-1,1, rcnp.GR_PH(0)); // from hist.def
      obj.FillHistogram("GR","A[X]",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
      obj.FillHistogram("GR","Y[A]",600,-0.15,0.15,rcnp.GR_TH(0),500,-50,50,rcnp.GR_Y(0));
      obj.FillHistogram("GR","Y[B]",500,-0.1,0.1,rcnp.GR_PH(0),500,-50,50,rcnp.GR_Y(0));


      dirname = "GR_new";
      obj.FillHistogram(dirname,"X[A]",300,-0.15,0.15,rcnp.GR_TH(0),1200,-600,600,rcnp.GR_X(0));
      obj.FillHistogram(dirname,"X[Y]",200,-100,100,rcnp.GR_Y(0),1200,-600,600,rcnp.GR_X(0));
      obj.FillHistogram(dirname,"X[B]",250,-0.1,0.1,rcnp.GR_PH(0),1200,-600,600,rcnp.GR_X(0));
      //obj.FillHistogram(dirname,"X[RF]",500,0,0,rcnp.GR_RF(0),1200,-600,600,rcnp.GR_X(0));
      obj.FillHistogram(dirname,"RF[A]",1000,-1,1,rcnp.GR_TH(0),500,0,0,rcnp.GR_RF(0));

      auto rf_Acor = rcnp.GR_RF(0)-(-1914.5*rcnp.GR_TH(0));
      obj.FillHistogram(dirname,"RF_Acor[A]",1000,-1,1,rcnp.GR_TH(0),500,0,0,rf_Acor);
      obj.FillHistogram(dirname,"RF[X]",1200,-600,600,rcnp.GR_X(0),500,0,0,rcnp.GR_RF(0));
      obj.FillHistogram(dirname,"RF_Acor[X]",1200,-600,600,rcnp.GR_X(0),500,0,0,rf_Acor);
      auto rf_Acor_Xcor = rf_Acor - (0.17205*rcnp.GR_X(0));
      obj.FillHistogram(dirname,"RF_Acor_Xcor[X]",1200,-600,600,rcnp.GR_X(0),500,0,0,rf_Acor_Xcor);
      obj.FillHistogram(dirname,"DE1[X]",1200,-600,600,rcnp.GR_X(0),2000,0,2000, hit.GetMeanPlastE1());
      obj.FillHistogram(dirname,"DE2[X]",1200,-600,600,rcnp.GR_X(0),2000,0,2000, hit.GetMeanPlastE2());

      obj.FillHistogram(dirname,"Y[A]",300,-0.15,0.15,rcnp.GR_TH(0),200,-100,100,rcnp.GR_Y(0));
      obj.FillHistogram(dirname,"Y[B]",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,rcnp.GR_Y(0));
      obj.FillHistogram(dirname,"Y[X]",1200,-600,600,rcnp.GR_X(0),200,-100,100,rcnp.GR_Y(0));
      obj.FillHistogram(dirname,"Y[X]",1200,-600,600,rcnp.GR_X(0),200,-100,100,rcnp.GR_Y(0));
      obj.FillHistogram(dirname,"Y[X]",1200,-600,600,rcnp.GR_X(0),200,-100,100,rcnp.GR_Y(0));
      obj.FillHistogram(dirname,"Y[X]",1200,-600,600,rcnp.GR_X(0),200,-100,100,rcnp.GR_Y(0));

      obj.FillHistogram(dirname,"A[RF]",500,700,1200,rf,1000,-1,1, rcnp.GR_TH(0));


      auto ycor = rcnp.GR_Y(0)+892.46*rcnp.GR_PH(0);
      obj.FillHistogram(dirname,"Y[B]cor",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,ycor);
      if (ycor<11 && ycor >-8) {
        obj.FillHistogram("GR","A[X]_gateYcor",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
        obj.FillHistogram(dirname,"Y[B]cor",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,ycor);

      }

      obj.FillHistogram(dirname,"DE1[RF]",1000,0,0,rf,2000,0,2000, hit.GetMeanPlastE1());
      obj.FillHistogram(dirname,"DE2[RF]",1000,0,0,rf,2000,0,2000, hit.GetMeanPlastE2());
      obj.FillHistogram(dirname,"DE1[RF_Acor_Xcor]",500,0,0,rf_Acor_Xcor,1000,0,2000, hit.GetMeanPlastE1());
      obj.FillHistogram(dirname,"DE2[RF_Acor_Xcor]",500,0,0,rf_Acor_Xcor,1000,0,2000, hit.GetMeanPlastE2());
      obj.FillHistogram(dirname,"DE1[dE2]",2000,0,2000, hit.GetMeanPlastE2(),2000,0,2000, hit.GetMeanPlastE1());
      obj.FillHistogram(dirname,"dE1[A]",1000,-1,1, rcnp.GR_TH(0),2000,0,2000, hit.GetMeanPlastE1());
      obj.FillHistogram(dirname,"dE2[A]",1000,-1,1, rcnp.GR_TH(0),2000,0,2000, hit.GetMeanPlastE2());


      // raytracing
      // double A=0,B=0;
      // std::tie(A,B) = raytrace(rcnp.GR_X(0),rcnp.GR_TH(0),rcnp.GR_Y(0));
      // obj.FillHistogram(dirname,"B[A]",500,0,0,A,500,0,0,B);

      // for (auto& cut : xcuts) {
      //   if (x < cut+xwidth && x >= cut) {
      //     stream.str(""); stream << "_x[" << cut << "," << cut+xwidth << ")";
      //     obj.FillHistogram(dirname,"B[A]"+stream.str(),500,0,0,A,500,0,0,B);
      //   }
      // }



    } // end rayid == 0 (good reconstruction)

    if (rcnp.GR_ADC()) {
      auto& adc = *rcnp.GR_ADC();
      for (int i=0; i<4; i++) {
        stream.str(""); stream << "GR_ADC" << i;
        obj.FillHistogram("GR",stream.str().c_str(), 1000,0,2000, adc[i]);
      }
    }


    for (auto const& labr_hit : hit.GetLaBr()) {

      int channum = labr_hit.channel;
      stream.str(""); stream << "LaBrLeading" << channum;
      obj.FillHistogram("GR", stream.str(), 10000,-40000, 40000, labr_hit.qtc_le);

      stream.str(""); stream << "LaBr" << channum << "_LE[LaBr_width]";
      obj.FillHistogram("GR", stream.str(), 1000, -5000, 15000, labr_hit.width, 1000,-40000, 40000, labr_hit.qtc_le);


      stream.str(""); stream << "LaBrWidth" << channum;
      obj.FillHistogram("GR", stream.str(), 10000, -5000, 15000, labr_hit.width);

      if (rcnp.GR_X(0) != BAD_NUM) {
        obj.FillHistogram("GR","X_LaBr",1200,-600,600,rcnp.GR_X(0),10000,-5000,15000,labr_hit.width);
      }

      /* !!!!!!!!! cut on prompt timing peak */ // Needs to be updated for e441 !!!!!!!!!!!
      if ((labr_hit.qtc_le>=-4000) && (labr_hit.qtc_le<=-3000)) {

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
    auto rf = rcnp.GR_RF(0);
    auto x = rcnp.GR_X(0);
    auto Ex = x*0.01074+6.872; // calibrated for 12C @ commissioning on 04.10.16

    // coincidence rate
    static ULong_t first_timestamp = grtime;
    if (first_timestamp) {
      auto rate = (grtime-first_timestamp)/1e8;
      //cout << grtime << " " << first_timestamp << endl;
      obj.FillHistogram("Coincident","Rate",3000,0,30000, rate);
    }

    auto ycor = rcnp.GR_Y(0)+892.46*rcnp.GR_PH(0);
    auto a_low = 0.00007835619*rcnp.GR_X(0) - 0.01;
    auto a_mid_low = 0.00007835619*rcnp.GR_X(0) + 0.03;
    auto a_mid_high = 0.00007835619*rcnp.GR_X(0) + 0.07;
    auto a_high = 0.00007835619*rcnp.GR_X(0) + 0.11;
    auto a = rcnp.GR_TH(0);



    // coincidence time difference
    for (auto& core_hit : cagra) {

      int detector = core_hit.GetDetnum();
      char core_leaf = core_hit.GetLeaf();
      string chan = core_leaf + std::to_string(core_hit.GetSegnum());
      size_t crystal_id = TCagra::GetCrystalId(detector,core_leaf);

      auto cagratime = core_hit.Timestamp();
      auto tdiff = cagratime-grtime;

      obj.FillHistogram("Coincident","Diff_CAGRA_GR", 1000,-500,1500,cagratime-grtime);

      stream.str("");
      stream << "TimeDiff_" << detector << "_" << chan;
      obj.FillHistogram("Coincident",stream.str().c_str(),1000,-500,1500,tdiff);


      // stream.str(""); stream << "Egam[tdiff]_" <<detector << "_" << chan;
      // obj.FillHistogram(stream.str(),5000,0,10000,core_hit.GetCorrectedEnergy(),1000,-500,1500,tdiff);


      // doppler reconstruction
      static const double ke_projectile = 600.; // MeV
      static const double m_projectile = li6.GetMass();
      static const double e_projectile = m_projectile + ke_projectile;
      static const double p_projectile = TMath::Sqrt(e_projectile*e_projectile-m_projectile*m_projectile);
      static const double beta = p_projectile/e_projectile;

      dirname = "ParticleGamma";
      auto Ecm = core_hit.GetDoppler(beta,pos::core_only);
      auto Elab = core_hit.GetCorrectedEnergy();

      obj.FillHistogram(dirname,"Ecm[Elab]",
                        2500,0,10000,Elab,
                        2500,0,10000,Ecm);
      obj.FillHistogram(dirname,"EdopplerSummary",
                        6000,-2000,22000,Ecm,
                        49,0,49,crystal_id);
      obj.FillHistogram(dirname,"EdopplerSum",7000,0,21000,Ecm);


      if (ycor>12 && ycor<31) {
        dirname = "PG_random";
        obj.FillHistogram(dirname,"Diff_CAGRA_GR", 1000,-500,1500,cagratime-grtime);
        obj.FillHistogram(dirname,"A[X]_gateYcor",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
        //stream.str(""); stream << "TimeDiff_" << detector << "_" << chan;
        obj.FillHistogram(dirname,stream.str().c_str(),1000,-500,1500,tdiff);
      }

      obj.FillHistogram(dirname,"Y[B]cor",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,ycor);
      //obj.FillHistogram(dirname,"Y[B]cor",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,ycor);

      if (ycor<11 && ycor >-8) {
        obj.FillHistogram(dirname,"Diff_CAGRA_GR", 1000,-500,1500,cagratime-grtime);
        obj.FillHistogram(dirname,"A[X]_gateYcor",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
        //stream.str(""); stream << "TimeDiff_" << detector << "_" << chan;
        obj.FillHistogram(dirname,stream.str().c_str(),1000,-500,1500,tdiff);

        if (a > a_low && a <= a_mid_low) {
          obj.FillHistogram(dirname+"_gates","A[X]_gateYcor_lowgateAX", 300,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
          obj.FillHistogram(dirname+"_gates","Diff_CAGRA_GR_lowgateAX", 1000,-500,1500,cagratime-grtime);
        } else if ( a > a_mid_low && a <= a_mid_high) {
          obj.FillHistogram(dirname+"_gates","A[X]_gateYcor_midgateAX", 300,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
          obj.FillHistogram(dirname+"_gates","Diff_CAGRA_GR_midgateAX", 1000,-500,1500,tdiff);

          // gate on prompt event timing peak
          if (tdiff>=252 && tdiff<=280) {
            int detector = core_hit.GetDetnum();
            char core_leaf = core_hit.GetLeaf();
            string chan = core_leaf + std::to_string(core_hit.GetSegnum());
            size_t crystal_id = TCagra::GetCrystalId(detector,core_leaf);

            // cagra core energy summary
            obj.FillHistogram("CrystalEnergySummaryPrompt",
                              6000,-2000,22000,core_hit.GetCorrectedEnergy(),
                              49,0,49,crystal_id);

            // cagra core energy summary vs. rcnp.GR_X(0)
            obj.FillHistogram("CrystalEnergySummaryPrompt_vs_X",
                            300,-600,600,rcnp.GR_X(0),
                            6000,-2000,22000,core_hit.GetCorrectedEnergy());
            obj.FillHistogram("EgammaSumPrompt",7000,0,21000,core_hit.GetCorrectedEnergy());

          }

          // gate on random events side band of timing peak
          if (tdiff>=212 && tdiff<=240) {
            int detector = core_hit.GetDetnum();
            char core_leaf = core_hit.GetLeaf();
            string chan = core_leaf + std::to_string(core_hit.GetSegnum());
            size_t crystal_id = TCagra::GetCrystalId(detector,core_leaf);

            // cagra core energy summary
            obj.FillHistogram("Summary","rand_CrystalEnergySummaryPrompt",
                              6000,-2000,22000,core_hit.GetCorrectedEnergy(),
                              49,0,49,crystal_id);

            // cagra core energy summary vs. rcnp.GR_X(0)
            obj.FillHistogram("Summary","rand_CrystalEnergySummaryPrompt_vs_X",
                            300,-600,600,rcnp.GR_X(0),
                            6000,-2000,22000,core_hit.GetCorrectedEnergy());

          }
        } else if (a > a_mid_high && a <=a_high) {
          obj.FillHistogram(dirname+"_gates","A[X]_gateYcor_highgateAX",300,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
          obj.FillHistogram(dirname+"_gates","Diff_CAGRA_GR_highgateAX", 1000,-500,1500,cagratime-grtime);
        }
      }
    } // end cagra analysis

    // LaBr3 prompt analysis
    for (auto const& labr_hit : hit.GetLaBr()) {
      int channum = labr_hit.channel;
      auto labr_rf = labr_hit.qtc_le-rf;
      if (a > a_low && a <= a_mid_low) {
        stream.str(""); stream << "LaBrLeading_lowgateAX";
        obj.FillHistogram(dirname+"_gates", stream.str(), 10000,-40000, 40000, labr_hit.qtc_le);
      } else if ( a > a_mid_low && a <= a_mid_high) {
        stream.str(""); stream << "LaBrLeading_RF_midgateAX_" << channum;
        obj.FillHistogram(dirname+"_gates", stream.str(), 1000, -25000, 35000, labr_rf);
        //stream.str(""); stream << "LaBrLeading_RF_midgateAX_" << channum;
        //obj.FillHistogram(dirname+"_gates", stream.str(), 1000, 0, 0, labr_rf);
        if ((labr_rf>=-2300) && (labr_rf<=-1600)) {
          stream.str(""); stream << "X_LaBrWidth_" << channum;
          obj.FillHistogram(dirname+"_gates", stream.str(),300,-600,600,rcnp.GR_X(0),2500,0,20000,labr_hit.width);
          stream.str(""); stream << "X_LaBrWidth_sum";
          obj.FillHistogram(dirname+"_gates", stream.str(),300,-600,600,rcnp.GR_X(0),2500,0,20000,labr_hit.width);
        }
        if ((labr_hit.qtc_le>=-3000) && (labr_hit.qtc_le<=-2300)) {
          stream.str(""); stream << "rand_X_LaBrWidth_" << channum;
          obj.FillHistogram(dirname+"_gates", stream.str(),300,-600,600,rcnp.GR_X(0),2500,0,20000,labr_hit.width);
          stream.str(""); stream << "rand_X_LaBrWidth_sum";
          obj.FillHistogram(dirname+"_gates", stream.str(),300,-600,600,rcnp.GR_X(0),2500,0,20000,labr_hit.width);
        }
      }else if (a > a_mid_high && a <=a_high) {
        stream.str(""); stream << "LaBrLeading_highgateAX";
        obj.FillHistogram(dirname+"_gates", stream.str(), 10000,-40000, 40000, labr_hit.qtc_le);
      }
    } // end labr3 analysis
  }
}

void LoadCuts() {
  // Example of how to load a cut once.
  // statically define the cut in e441.h
  // if(!your_fav_cut) {
  //   TPreserveGDirectory Preserve;
  //   TFile fcut("./cuts/newHe3cut.root");
  //   your_fav_cut = (TCutG*)fcut.Get("_cut0");
  //   std::cout << "Loaded he3 gate." << std::endl;
  // }
}

void LoadRaytraceParams(size_t xdeg, size_t adeg, size_t ydeg) {
  static bool once = false;
  if (once) { return; }
  once = true;
  std::cout << "Loading Raytrace Parameters. " << std::endl;

  // sieve slit transformation coefficients
  // output from sieveslit.py
  TGrandRaidenHit::SetRaytraceParams(
    { // a fit parameters
      25478959.64415743201971054,
        -0.01771996157436579,
        -0.00000375686641008,
        -25478962.08887941017746925,
        339.32639498225580610,
        285.41500265652570079
    },
    { // b fit parameters
      -8.77979136468771948,
        3.57167199982826045,
        144.77038232928563843,
        -102.62252113966886213,
        -1669.00084319835968927,
        1060.28866587772699859,
        0.00281259943988080,
        0.00390660581203093,
        -0.57697427700803183,
        -0.03624075884543333,
        10.73431026156209711,
        -0.86971847190018392,
        0.00002737725798545,
        -0.00000000249027871,
        -0.00140344227013718,
        0.00008863729153273,
        0.01262661852581874,
        -0.00031269175012557
        },
    xdeg, // polynomial degree in fit for x
    adeg,
    ydeg
    );
}
