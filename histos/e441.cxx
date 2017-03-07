#include "TRuntimeObjects.h"
#include "e441.h"

template<typename ...Args>
void hist(bool conditional, TRuntimeObjects& obj, Args&&... args) {
  if (conditional) { obj.FillHistogram(std::forward<Args>(args)...); }
}


static double m_target = 0;

// ----------------------------------------------------------------------
// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  InitTarget();
  LoadCuts();
  LoadRaytraceParams(3,1,1,1);

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

void InitTarget() {
  static bool once = false;
  if (once) { return; }
  once = true;
  std::cout << "Target " << int(target) << " selected." << std::endl;

  switch(target) {
  case Target::n12C:
    m_target = c12.GetMass();
    break;
  case Target::n56Fe:
    m_target = fe56.GetMass();
    break;
  case Target::n24Mg:
    m_target = mg24.GetMass();
    break;
  case Target::n93Nb:
  case Target::n93Nb2:
    m_target = nb93.GetMass();
    break;
  default:
    throw std::runtime_error("Unknown target");
    m_target = 0;
    break;
  }
}

void PoleZeroHistos(TRuntimeObjects& obj, TCagraHit& core_hit, string local_dirname = "") {
  auto flags = core_hit.GetFlags();
  if (TANLEvent::PileUpFlag(flags) || TANLEvent::PileUpOnlyFlag(flags)) { return; }

  int detector = core_hit.GetDetnum();
  string chan = core_hit.GetLeaf() + std::to_string(core_hit.GetSegnum());

  Double_t prerise = core_hit.GetPreRise()/TANLEvent::GetShapingTime();
  stream.str("");  stream << "Prerise[Q]_" << detector << "_" << chan;
  hist(false,obj,local_dirname, stream.str(),2000,0,10000,core_hit.GetCharge(),1250,6000,8500,prerise);
  for (auto& seg_hit : core_hit) {
    string seg_chan = seg_hit.GetLeaf() + std::to_string(seg_hit.GetSegnum());
    stream.str("");  stream << "Prerise[Q]_" << detector << "_" << seg_chan;
    hist(false,obj,local_dirname, stream.str(),2000,0,0,seg_hit.GetCharge(),1250,6000,8500,prerise);
  }
  // stream.str("");  stream << "Q[Prerise]_" << detector << "_" << chan;
  // hist(false,obj,local_dirname, stream.str(),1250,6000,8500,prerise,3000,0,6000,core_hit.GetCharge());
  // for (auto& seg : core_hit) {
  //   string seg_chan = seg_hit.GetLeaf() + std::to_string(seg_hit.GetSegnum());
  //   stream.str("");  stream << "Q[Prerise]Seg_" << detector << "_" << seg_chan;
  //   hist(false,obj,local_dirname, stream.str(),1250,6000,8500,prerise,3000,0,6000,seg.GetCharge());
  // }

  stream.str("");  stream << "Prerise[E_pzcor_basesample]_" << detector << "_" << chan;
  hist(false,obj,local_dirname, stream.str(),3000,0,6000,core_hit.GetCorrectedEnergy(core_hit.GetBaseSample()),1250,6000,8500,prerise);
  stream.str("");  stream << "Prerise[E_pzcor_constant]_" << detector << "_" << chan;
  hist(false,obj,local_dirname, stream.str(),3000,-2000,4000,core_hit.GetCorrectedEnergy(),1250,6000,8500,prerise);


  stream.str(""); stream << "E_pzcor_constant" << detector << "_" << chan;
  auto pzchan = core_hit.GetCorrectedEnergy();
  hist(false,obj,local_dirname,stream.str(),4000,0,12000,pzchan);
  for (auto& seg_hit : core_hit) {
    string seg_chan = seg_hit.GetLeaf() + std::to_string(seg_hit.GetSegnum());
    stream.str(""); stream << "E_pzcor_constant" << detector << "_" << seg_chan;
    hist(false,obj,local_dirname,stream.str(),4000,0,12000,seg_hit.GetCorrectedEnergy());
    stream.str("");  stream << "Prerise[E_pzcor_constant]_" << detector << "_" << seg_chan;
    hist(false,obj,local_dirname, stream.str(),4000,0,12000,seg_hit.GetCorrectedEnergy(),1250,6000,8500,prerise);
  }

}

void PileUp (TRuntimeObjects& obj, TCagraHit& core_hit) {
  dirname = "PileUp";

  auto flags = core_hit.GetFlags();
  auto pileup = TANLEvent::PileUpFlag(flags);
  if (pileup) {
    hist(false,obj,dirname,"Summary",10,1,10,3);
    PoleZeroHistos(obj,core_hit,dirname);
  } else {
    hist(false,obj,dirname,"Summary",10,1,10,7);
    PoleZeroHistos(obj,core_hit,"NoPileUp");
  }
}


void MakeCAGRAHistograms(TRuntimeObjects& obj, TCagra& cagra) {

  for (auto& core_hit : cagra) {

    //PileUp(obj,core_hit);

    int detector = core_hit.GetDetnum();
    char core_leaf = core_hit.GetLeaf();
    string chan = core_leaf + std::to_string(core_hit.GetSegnum());
    auto crystal_id = TCagra::GetCrystalId(detector,core_leaf);

    // cagra core energy summary
    hist(true,obj,"CrystalEnergySummary",
         //6000,-2000,22000,core_hit.GetCorrectedEnergy(),
         18000,-2000,22000,core_hit.GetCorrectedEnergy(),
         49,0,49,crystal_id);
    hist(true,obj,"EgammaSum",7000,0,21000,core_hit.GetCorrectedEnergy());


    static int num_segments = TCagra::NumSegments();
    for (auto& seg_hit : core_hit) {
      auto segment_id = TCagra::GetSegmentId(detector,core_leaf,seg_hit.GetSegnum(),seg_hit.GetSystem());
      hist(true,obj,"SegmentEnergySummary",
           18000,-2000,22000,seg_hit.GetCorrectedEnergy(),
           num_segments+1,0,num_segments+1,segment_id);
    }



    static RateOffset rate_offset(1 /* second */, target);
    auto first_ts = rate_offset.Update(core_hit.Timestamp(),(target == Target::n93Nb || target == Target::n93Nb2) ? crystal_id : 0);
    if (first_ts > 1e6) {

      // cagra core time summary
      hist(true,obj,"CrystalTimeSummary",
                        1000,0,4000,(core_hit.Timestamp()-first_ts)*10/1.0e9, // in seconds - 1 bin = 4 seconds
                        49,0,49,crystal_id);


      // cagra energy vs time to study rate dependence (js)
      obj.FillHistogram("EnergyVsTime",
                        1000,0,4000,(core_hit.Timestamp()-first_ts)*10/1.0e9, // 1 bin = 4 seconds
                        18000,-2000,22000,core_hit.GetCorrectedEnergy());

      // cagra energy vs time to study rate dependence (js)
      obj.FillHistogram("EnergyVsTime_baseline",
                        1000,0,4000,(core_hit.Timestamp()-first_ts)*10/1.0e9, // 1 bin = 4 seconds
                        18000,-2000,22000,core_hit.GetCorrectedEnergy(core_hit.GetBaseSample()));

      obj.FillHistogram("EnergyVsTime_withRateCorrection",
                        1000,0,4000,(core_hit.Timestamp()-first_ts)*10/1.0e9, // 1 bin = 4 seconds
                        18000,-2000,22000,core_hit.GetCorrectedEnergy()+rate_offset());


      hist(false,obj,"NumEvents","cagra_hits_time",1000,0,8000,(core_hit.Timestamp()-first_ts)*10/1.0e9);
      hist(false,obj,"NumEvents","prerise[time]",1000,0,8000,(core_hit.Timestamp()-first_ts)*10/1.0e9,1250,6000,8500,core_hit.GetPreRise()/TANLEvent::GetShapingTime());
      hist(false,obj,"NumEvents","postrise[time]",1000,0,8000,(core_hit.Timestamp()-first_ts)*10/1.0e9,1250,6000,8500,core_hit.GetPostRise()/TANLEvent::GetShapingTime());
      hist(false,obj,"NumEvents","postrise[prerise]",1250,0,0,core_hit.GetPreRise()/TANLEvent::GetShapingTime(),1250,0,0,core_hit.GetPostRise()/TANLEvent::GetShapingTime());

      for (auto& seg_hit : core_hit) {
        auto segment_id = TCagra::GetSegmentId(detector,core_leaf,seg_hit.GetSegnum(),seg_hit.GetSystem());
        hist(true,obj,"SegmentTimeSummary",
             1000,0,4000,(seg_hit.Timestamp()-first_ts)*10/1.0e9, // in seconds - 1 bin = 4 seconds
             num_segments+1,0,num_segments+1,segment_id);
      }

    }


    auto position = core_hit.GetPosition(pos::core_only);
    hist(false,obj,"ArrayHits",
                      180,0,180,position.Theta()*180/TMath::Pi(),
                      360,-180,180,position.Phi()*180/TMath::Pi());


    PoleZeroHistos(obj,core_hit,"PoleZero");


    // central contact signals
    // name = "Det_" + std::to_string(detector) + "_" + chan;
    // hist(false,obj,"CAGRA_Raw", name,2000,0,0,core_hit.GetCharge());

    if (core_hit.GetChannel() == 4 && core_hit.GetSystem()=='Y') {
      name = "BGO_" + std::to_string(detector) + "_" + chan;
      hist(false,obj,"BGO", name,2000,0,0,core_hit.GetCharge());
    }


    // // segment (side channel) signals
    // for (auto& segment : core_hit) {
    //   string seg_chan = segment.GetLeaf() + std::to_string(segment.GetSegnum());
    //   stream.str(""); stream << "Det_" << detector << "_" << seg_chan;
    //   hist(false,obj,"CAGRA_Raw",stream.str(),2000,0,0,segment.GetCharge());
    // }


    // same but for calibrated energies
    // stream.str("");
    // stream << "Det_" << detector << "_" << chan;
    // hist(false,obj,"CAGRA_Calibrated",stream.str(),2000,0,10000,core_hit.GetEnergy());
    // for (auto& segment : core_hit) {
    //   string seg_chan = segment.GetLeaf() + std::to_string(segment.GetSegnum());
    //   stream.str(""); stream << "Det_" << detector << "_" << seg_chan;
    //   hist(false,obj,"CAGRA_Calibrated",stream.str(),2000,0,10000,segment.GetEnergy());
    // }


  } // end loop over cagra hits

}

void MakeGRCorrections(TRuntimeObjects& obj, TGrandRaiden& gr, TCagra* cagra=nullptr, std::string dirname="FP_corections");

void MakeGrandRaidenHistograms(TRuntimeObjects& obj, TGrandRaiden& gr) {

  MakeGRCorrections(obj,gr);

  for (auto& hit : gr) {

    auto& rcnp = hit.GR();



    if (rcnp.GR_MYRIAD(0) != BAD_NUM) {
      hist(false,obj,"Timing","MyriadTimestamp",10000,1e9,5e12,hit.GetTimestamp());
    }

    static ULong_t prev_ts = 0;
    if (prev_ts) {
      hist(false,obj,"Timing","GR_EventPeriod",5000,100,50000,hit.GetTimestamp()-prev_ts);
    }
    prev_ts = hit.GetTimestamp();

    auto rf = rcnp.GR_RF(0);
    if (rf != BAD_NUM) {
      hist(false,obj,"GR","GR_RF",1000,0,0,rf);
    }


    // X, A, Y, B, RF, DE1, DE2
    // X[A],X[Y],X[B],X[RF],DE1[X],DE2[X]
    // Y[A],Y[B],Y[RF],DE1[Y],DE2[Y]
    // A[B], A[RF], DE1[A], DE2[A]
    // B[RF], DE1[B], DE2[B]
    // DE1[RF], DE2[RF], DE1[DE2]

    hist(false,obj,"GR","RayID",64,-16,48, rcnp.GR_RAYID(0));
    if (rcnp.GR_RAYID(0) == 0) { // if track reconstruction successful

      hist(false,obj,"GR","x",1200,-600,600, rcnp.GR_X(0));
      hist(false,obj,"GR","x_cal",1000,0,20, rcnp.GR_X(0)*0.01074+6.872);
      hist(false,obj,"GR","y",200,-100,100, rcnp.GR_Y(0));
      hist(false,obj,"GR","a",100,-1,1, rcnp.GR_TH(0)); // need to learn
      hist(false,obj,"GR","b",100,-1,1, rcnp.GR_PH(0)); // from hist.def
      hist(true,obj,"GR","a[x]",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
      hist(true,obj,"GR","y[a]",600,-0.15,0.15,rcnp.GR_TH(0),500,-50,50,rcnp.GR_Y(0));
      hist(true,obj,"GR","y[b]",500,-0.1,0.1,rcnp.GR_PH(0),500,-50,50,rcnp.GR_Y(0));


      if (rcnp.GR_TH(0)>-0.026 && rcnp.GR_TH(0) < -0.01) {
        hist(true,obj,"sieveslit_x1_a1","y[b]",500,-0.1,0.1,rcnp.GR_PH(0),500,-50,50,rcnp.GR_Y(0));
        hist(true,obj,"sieveslit_x1_a1","y[a]",600,-0.15,0.15,rcnp.GR_TH(0),500,-50,50,rcnp.GR_Y(0));
      }

      if (rcnp.GR_TH(0)>-0.01 && rcnp.GR_TH(0) < 0.006) {
        hist(true,obj,"sieveslit_x1_a2","y[b]",500,-0.1,0.1,rcnp.GR_PH(0),500,-50,50,rcnp.GR_Y(0));
        hist(true,obj,"sieveslit_x1_a2","y[a]",600,-0.15,0.15,rcnp.GR_TH(0),500,-50,50,rcnp.GR_Y(0));
      }

      dirname = "GR_new";
      hist(false,obj,dirname,"x[a]",300,-0.15,0.15,rcnp.GR_TH(0),1200,-600,600,rcnp.GR_X(0));
      hist(false,obj,dirname,"x[y]",200,-100,100,rcnp.GR_Y(0),1200,-600,600,rcnp.GR_X(0));
      hist(false,obj,dirname,"x[b]",250,-0.1,0.1,rcnp.GR_PH(0),1200,-600,600,rcnp.GR_X(0));
      //hist(false,obj,dirname,"X[RF]",500,0,0,rcnp.GR_RF(0),1200,-600,600,rcnp.GR_X(0));
      hist(false,obj,dirname,"RF[a]",1000,-1,1,rcnp.GR_TH(0),500,0,0,rcnp.GR_RF(0));

      auto rf_Acor = rcnp.GR_RF(0)-(-1914.5*rcnp.GR_TH(0));
      hist(false,obj,dirname,"RF_acor[a]",1000,-1,1,rcnp.GR_TH(0),500,0,0,rf_Acor);
      hist(false,obj,dirname,"RF[x]",1200,-600,600,rcnp.GR_X(0),500,0,0,rcnp.GR_RF(0));
      hist(false,obj,dirname,"RF_acor[x]",1200,-600,600,rcnp.GR_X(0),500,0,0,rf_Acor);
      auto rf_Acor_Xcor = rf_Acor - (0.17205*rcnp.GR_X(0));
      hist(false,obj,dirname,"RF_Acor_Xcor[x]",1200,-600,600,rcnp.GR_X(0),500,0,0,rf_Acor_Xcor);
      hist(false,obj,dirname,"DE1[x]",1200,-600,600,rcnp.GR_X(0),2000,0,2000, hit.GetMeanPlastE1());
      hist(false,obj,dirname,"DE2[x]",1200,-600,600,rcnp.GR_X(0),2000,0,2000, hit.GetMeanPlastE2());
      hist(false,obj,dirname,"DE3[x]",1200,-600,600,rcnp.GR_X(0),2000,0,2000, hit.GetMeanPlastE3());

      hist(true,obj,dirname,"y[a]",300,-0.15,0.15,rcnp.GR_TH(0),200,-100,100,rcnp.GR_Y(0));
      hist(true,obj,dirname,"y[b]",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,rcnp.GR_Y(0));
      hist(true,obj,dirname,"y[x]",1200,-600,600,rcnp.GR_X(0),200,-100,100,rcnp.GR_Y(0));
      //hist(false,obj,dirname,"Y[X]",1200,-600,600,rcnp.GR_X(0),200,-100,100,rcnp.GR_Y(0));
      //hist(false,obj,dirname,"Y[X]",1200,-600,600,rcnp.GR_X(0),200,-100,100,rcnp.GR_Y(0));
      //hist(false,obj,dirname,"Y[X]",1200,-600,600,rcnp.GR_X(0),200,-100,100,rcnp.GR_Y(0));

      hist(false,obj,dirname,"a[RF]",500,700,1200,rf,1000,-1,1, rcnp.GR_TH(0));


      auto ycor = rcnp.GR_Y(0)+892.46*rcnp.GR_PH(0);
      hist(false,obj,dirname,"y[b]cor",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,ycor);
      if (ycor<22 && ycor >-13) {
        hist(false,obj,"GR","a[x]_gateYcor",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
        hist(false,obj,dirname,"y[b]cor",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,ycor);
      }

      hist(false,obj,dirname,"DE1[RF]",1000,0,0,rf,2000,0,2000, hit.GetMeanPlastE1());
      hist(false,obj,dirname,"DE2[RF]",1000,0,0,rf,2000,0,2000, hit.GetMeanPlastE2());
      hist(false,obj,dirname,"DE3[RF]",1000,0,0,rf,2000,0,2000, hit.GetMeanPlastE3());
      
      hist(false,obj,dirname,"DE1[RF_Acor_Xcor]",500,0,0,rf_Acor_Xcor,1000,0,2000, hit.GetMeanPlastE1());
      hist(false,obj,dirname,"DE2[RF_Acor_Xcor]",500,0,0,rf_Acor_Xcor,1000,0,2000, hit.GetMeanPlastE2());
      hist(false,obj,dirname,"DE1[dE2]",2000,0,2000, hit.GetMeanPlastE2(),2000,0,2000, hit.GetMeanPlastE1());
      hist(false,obj,dirname,"DE2[DE3]",2000,0,2000, hit.GetMeanPlastE2(),2000,0,2000, hit.GetMeanPlastE3());
      
      hist(false,obj,dirname,"dE1[a]",1000,-1,1, rcnp.GR_TH(0),2000,0,2000, hit.GetMeanPlastE1());
      hist(false,obj,dirname,"dE2[a]",1000,-1,1, rcnp.GR_TH(0),2000,0,2000, hit.GetMeanPlastE2());
      hist(false,obj,dirname,"dE3[a]",1000,-1,1, rcnp.GR_TH(0),2000,0,2000, hit.GetMeanPlastE3());


      {
        auto ejectile = hit.GetEjectileVector(true);
        // missing mass
        ejectile.SetMag(hit.GetMomentum());
        auto p_ejectile = ejectile.Mag();
        auto e_ejectile = TMath::Sqrt(p_ejectile*p_ejectile+m_projectile*m_projectile);
        auto ke_ejectile = e_ejectile - m_projectile;
        double theta_cm=0,Ex=0,J_L=0;

        std::tie(theta_cm,Ex,J_L) = kine_2b(m_projectile,m_target,m_projectile,m_target,ke_projectile,ejectile.Theta(), ke_ejectile);

        hist(true,obj,"MissingMass","ThetaLab[MMEx]",600,5,65,Ex,300,0,0.15,ejectile.Theta());
      }

      // raytracing
      double A=0,B=0;
      std::tie(A,B) = hit.Raytrace(true);
      hist(true,obj,dirname,"B[A]",500,60,100,A,500,-65,65,B);
      hist(true,obj,dirname,"B[b]",500,-0.1,0.1,rcnp.GR_PH(0),500,-65,65,B);





      // missing mass
      auto p_ejectile = hit.GetMomentum();
      auto ejectile = hit.GetEjectileVector(true);
      // create ejectile vector


      hist(true,obj,"MissingMass","ThetaLab[p_ejectile]",600,0,0,p_ejectile,300,0,0.15,ejectile.Theta());
      hist(true,obj,"MissingMass","A[p_ejectile]",600,5,65,p_ejectile,500,0,0,A);
      hist(true,obj,"MissingMass","B[p_ejectile]",600,5,65,p_ejectile,500,0,0,B);
      hist(true,obj,"MissingMass","afp[p_ejectile]",600,5,65,p_ejectile,500,0,0,rcnp.GR_TH(0));
      hist(true,obj,"MissingMass","yfp[p_ejectile]",600,5,65,p_ejectile,500,0,0,rcnp.GR_Y(0));

      auto p_kine = p_ejectile - (-9.407e+01*std::pow(ejectile.Theta(),2) - 5.783e-03*ejectile.Theta() + 3.658e-05);
      hist(true,obj,"MissingMass","ThetaLab[p_kine]",600,0,0,p_kine,300,0,0.15,ejectile.Theta());
      hist(true,obj,"MissingMass","A[p_kine]",600,5,65,p_kine,500,0,0,A);
      hist(true,obj,"MissingMass","B[p_kine]",600,5,65,p_kine,500,0,0,B);
      hist(true,obj,"MissingMass","afp[p_kine]",600,5,65,p_kine,500,0,0,rcnp.GR_TH(0));
      hist(true,obj,"MissingMass","yfp[p_kine]",600,5,65,p_kine,500,0,0,rcnp.GR_Y(0));

      auto p_aber_kine = p_kine - (3.70e-03*std::pow(A,2) - 5.22e-01*A);
      hist(true,obj,"MissingMass","ThetaLab[p_aber_kine]",600,0,0,p_aber_kine,300,0,0.15,ejectile.Theta());
      hist(true,obj,"MissingMass","A[p_aber_kine]",600,5,65,p_aber_kine,500,0,0,A);
      hist(true,obj,"MissingMass","B[p_aber_kine]",600,5,65,p_aber_kine,500,0,0,B);
      hist(true,obj,"MissingMass","afp[p_aber_kine]",600,5,65,p_aber_kine,500,0,0,rcnp.GR_TH(0));
      hist(true,obj,"MissingMass","yfp[p_aber_kine]",600,5,65,p_aber_kine,500,0,0,rcnp.GR_Y(0));

      // applying abberation corrected momentum
      //p_ejectile -= (3.70e-03*std::pow(A,2) - 5.22e-01*A + 2.56e+03);
      p_ejectile -= (3.70e-03*std::pow(A,2) - 5.22e-01*A);
      ejectile.SetMag(p_ejectile);
      auto e_ejectile = TMath::Sqrt(p_ejectile*p_ejectile+m_projectile*m_projectile);
      auto ke_ejectile = e_ejectile - m_projectile;

      hist(false,obj,"MissingMass","ThetaLab",300,0,0,ejectile.Theta());
      hist(false,obj,"MissingMass","PhiLab",300,0,0,ejectile.Phi());
      if (ejectile.Theta() < 0.01) {
        hist(false,obj,"MissingMass","PhiLab_10mrad",300,0,0,ejectile.Phi());
      }

      double theta_cm=0,Ex=0,J_L=0;

      std::tie(theta_cm,Ex,J_L) = kine_2b( m_projectile, m_target, m_projectile, m_target, ke_projectile, ejectile.Theta(), ke_ejectile);

      hist(true,obj,"MissingMass","ThetaCM",300,0,0.15,theta_cm);
      hist(true,obj,"MissingMass","ReconstructedEx",600,5,65,Ex);

      hist(true,obj,"MissingMass","ThetaLab[Ex]",600,5,65,Ex,300,0,0.15,ejectile.Theta());
      hist(true,obj,"MissingMass","A[Ex]",600,5,65,Ex,500,0,0,A);
      hist(true,obj,"MissingMass","B[Ex]",600,5,65,Ex,500,0,0,B);
      hist(true,obj,"MissingMass","afp[Ex]",600,5,65,Ex,500,0,0,rcnp.GR_TH(0));
      hist(true,obj,"MissingMass","yfp[Ex]",600,5,65,Ex,500,0,0,rcnp.GR_Y(0));


      if (ejectile.Theta() < 0.01) {
        hist(false,obj,"MissingMass","ReconstructedEx_0_10_mrad",1024,0,0,Ex);
      }
      hist(false,obj,"MissingMass","KE_ejectile",1024,0,0,ke_ejectile);

      hist(false,obj,"MissingMass","KE_projectile",1024,0,0,ke_projectile);
      hist(false,obj,"MissingMass","ejectile_theta",1000,-.1,.1,ejectile.Theta());


      // LaBr3 prompt analysis
      dirname = "LaBr3_prompt";
      for (auto const& labr_hit : hit.GetLaBr()) {
        int channum = labr_hit.channel;
        if ((labr_hit.qtc_le >= -3000) && (labr_hit.qtc_le < -2000)) { // prompt  // 1000 wide
          stream.str(""); stream << "X_LaBrE_" << channum;
          hist(false,obj,dirname, stream.str(),300,-600,600,rcnp.GR_X(0),500,0,20000,labr_hit.GetEnergy());
          stream.str(""); stream << "X_LaBrE_sum";
          hist(false,obj,dirname, stream.str(),300,-600,600,rcnp.GR_X(0),500,0,20000,labr_hit.GetEnergy());

          // old
          hist(false,obj,dirname,"GR_X",1200,-600,600, rcnp.GR_X(0));
          hist(false,obj,dirname,"GR_Y",200,-100,100, rcnp.GR_Y(0));
          hist(false,obj,dirname,"GR_Theta",100,-1,1, rcnp.GR_TH(0));
          hist(false,obj,dirname,"GR_Phi",100,-1,1, rcnp.GR_PH(0));
          hist(false,obj,dirname,"X_TH",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));

          hist(false,obj,dirname,"GR_Theta_Phi",100,-1,1, rcnp.GR_TH(0),100,-1,1, rcnp.GR_PH(0));
          hist(false,obj,dirname,"GR_X_Y",1200,-600,600, rcnp.GR_X(0),200,-100,100, rcnp.GR_Y(0));


          stream.str(""); stream << "LaBrE_LEGate" << channum;
          hist(false,obj,dirname,stream.str(), 10000, -5000, 15000, labr_hit.GetEnergy());
        }
        else if ((labr_hit.qtc_le >= -2000) && (labr_hit.qtc_le < 28000)) { // random // 30000 wide
          stream.str(""); stream << "rand_X_LaBrE_" << channum;
          hist(false,obj,dirname, stream.str(),300,-600,600,rcnp.GR_X(0),500,0,20000,labr_hit.GetEnergy());
          stream.str(""); stream << "rand_X_LaBrE_sum";
          hist(false,obj,dirname, stream.str(),300,-600,600,rcnp.GR_X(0),500,0,20000,labr_hit.GetEnergy());
        }

        stream.str(""); stream << "LaBrLeading" << channum;
        hist(false,obj,"GR", stream.str(), 10000,-40000, 40000, labr_hit.qtc_le);

        stream.str(""); stream << "LaBr" << channum << "_LE[LaBr_E]";
        hist(false,obj,"GR", stream.str(), 2000, -1000, 19000, labr_hit.GetEnergy(), 1000,-4000, 1000, labr_hit.qtc_le);
      } // end labr3 analysis
      // end rayid == 0 (good reconstruction)
    }

    auto labr_hits = hit.GetLaBr();
    for (int i=0; i<labr_hits.size(); i++) {
      stream.str(""); stream << "LaBrE" << labr_hits[i].channel;
      hist(false,obj,"GR", stream.str(), 12500, -5000, 20000, labr_hits[i].GetEnergy());
      for (int j=0; j<labr_hits.size(); j++) {
        if (j!=i && labr_hits[i].channel == labr_hits[j].channel) {
          stream.str(""); stream << "LaBrE_Gamma_Gamma" << labr_hits[i].channel;
          hist(false,obj,"GR", stream.str(),
               1250, -5000, 20000, labr_hits[i].GetEnergy(),
               1250, -5000, 20000, labr_hits[j].GetEnergy());
        }
      }
    }

    if (rcnp.GR_ADC()) {
      auto& adc = *rcnp.GR_ADC();
      for (int i=0; i<6; i++) {
        stream.str(""); stream << "GR_ADC" << i;
        hist(false,obj,"GR",stream.str().c_str(), 1000,0,2000, adc[i]);
      }
    }







  }


}
void MakeGRCorrections(TRuntimeObjects& obj, TGrandRaiden& gr, TCagra* cagra, std::string dirname) {

  for (auto& hit : gr) {

    auto& rcnp = hit.GR();

    auto grtime = hit.GetTimestamp();
    auto rf = rcnp.GR_RF(0);
    auto x = rcnp.GR_X(0);
    auto a = rcnp.GR_TH(0);
    auto rf_cor = rf - -1121.882*a;

    auto hist_vec = [](
      bool conditional,
      TRuntimeObjects& obj,
      std::string dirname,
      std::string histname,
      int xbins,
      int xlow,
      int xhigh,
      std::vector<double>* xvals) {
      if (xvals) {
        for (auto& val : *xvals) {
          hist(conditional,obj,dirname,histname,xbins,xlow,xhigh,val);
        }
      }
    };

    // PID gate (it's pretty clean so only an RF gate is really needed)
    if ((rf_cor >= 900 && rf_cor < 975 ) || (rf_cor >=1700 && rf_cor < 1775)) {
      hist_vec(true,obj,dirname,"GR_TDCR_X1",500,0,500,rcnp.GR_TDCR_X1());
      hist_vec(true,obj,dirname,"GR_TDCR_U1",500,0,500,rcnp.GR_TDCR_U1());
      hist_vec(true,obj,dirname,"GR_TDCR_X2",500,0,500,rcnp.GR_TDCR_X2());
      hist_vec(true,obj,dirname,"GR_TDCR_U2",500,0,500,rcnp.GR_TDCR_U2());

      // These histograms are required if you want to run the
      // drift time to length conversion script (./util/rcnp/gr_tdc.py)
      hist_vec(true,obj,dirname,"GR_TDC_X1",500,0,500,rcnp.GR_TDC_X1());
      hist_vec(true,obj,dirname,"GR_TDC_U1",500,0,500,rcnp.GR_TDC_U1());
      hist_vec(true,obj,dirname,"GR_TDC_X2",500,0,500,rcnp.GR_TDC_X2());
      hist_vec(true,obj,dirname,"GR_TDC_U2",500,0,500,rcnp.GR_TDC_U2());

      hist_vec(true,obj,dirname,"GR_DRIFT_X1",500,-10,10,rcnp.GR_DRIFT_X1());
      hist_vec(true,obj,dirname,"GR_DRIFT_U1",500,-10,10,rcnp.GR_DRIFT_U1());
      hist_vec(true,obj,dirname,"GR_DRIFT_X2",500,-10,10,rcnp.GR_DRIFT_X2());
      hist_vec(true,obj,dirname,"GR_DRIFT_U2",500,-10,10,rcnp.GR_DRIFT_U2());

      for (auto& tdc : *rcnp.GR_TDC_X1()) {
        obj.FillHistogram(dirname,"GR_TDC_X1",500,0,500,tdc);
      }

    }

    if (rcnp.GR_RAYID(0) == 0) {

      hist(true,obj,dirname,"RF_acor[A]",1000,-1,1,a,500,600,1800,rf_cor);
      hist(true,obj,dirname,"RF_acor[X]",1200,-600,600,x,500,600,1800,rf_cor);

      rf_cor -= 0.119966*x;
      hist(true,obj,dirname,"RFcor",1000,0,0,rf_cor);
      hist(true,obj,dirname,"RF_axcor[X]",1200,-600,600,x,500,600,1800,rf_cor);
      hist(true,obj,dirname,"RF_axcor[A]",1000,-1,1,a,500,600,1800,rf_cor);

      hist(true,obj,dirname,"DE1[RFcor]",900,200,2000,rf_cor,350,0,350, hit.GetMeanPlastE1());
      hist(true,obj,dirname,"DE2[RFcor]",900,200,2000,rf_cor,350,0,350, hit.GetMeanPlastE2());
      hist(true,obj,dirname,"DE3[RFcor]",900,200,2000,rf_cor,2000,0,1300, hit.GetMeanPlastE3());

      hist(true,obj,dirname,"DE1",350,0,350, hit.GetMeanPlastE1());
      hist(true,obj,dirname,"DE2",350,0,350, hit.GetMeanPlastE2());
      hist(true,obj,dirname,"DE3",2000,0,0, hit.GetMeanPlastE3());

      // PID gate (it's pretty clean so only an RF gate is really needed)
      if ((rf_cor >= 900 && rf_cor < 975 ) || (rf_cor >=1700 && rf_cor < 1775)) {
        hist(true,obj,dirname,"DE1_rfgate",350,0,350, hit.GetMeanPlastE1());
        hist(true,obj,dirname,"DE2_rfgate",350,0,350, hit.GetMeanPlastE2());
        hist(true,obj,dirname,"DE3_rfgate",2000,0,0, hit.GetMeanPlastE3());
        hist(true,obj,dirname,"DE1[RFcor]rfgate",900,200,2000,rf_cor,350,0,350, hit.GetMeanPlastE1());
        hist(true,obj,dirname,"DE2[RFcor]rfgate",900,200,2000,rf_cor,350,0,350, hit.GetMeanPlastE2());
        hist(true,obj,dirname,"DE3[RFcor]rfgate",900,200,2000,rf_cor,2000,0,1300, hit.GetMeanPlastE3());



        if (cagra) {

          auto ycor = rcnp.GR_Y(0)+892.46*rcnp.GR_PH(0);

          for (auto& core_hit : *cagra) {

            int detector = core_hit.GetDetnum();
            char core_leaf = core_hit.GetLeaf();
            string chan = core_leaf + std::to_string(core_hit.GetSegnum());
            auto crystal_id = TCagra::GetCrystalId(detector,core_leaf);
            if (crystal_id > 50) { continue; }


            bool bgo_hit = false;
            auto cagratime = core_hit.Timestamp();
            auto tdiff = cagratime-grtime;
            hist(true,obj,dirname,"Diff_CAGRA_GR", 1000,-500,1500,tdiff);

            // doppler reconstruction
            auto Ecm = core_hit.GetDoppler(beta,pos::core_only);
            auto ejectile = hit.GetEjectileVector(true);
            ejectile.SetMag(hit.GetMomentum());
            auto Ecm_particle = core_hit.GetDoppler(beta,pos::core_only,ejectile);
            auto Elab = core_hit.GetCorrectedEnergy();
            auto p_ejectile = ejectile.Mag();
            auto e_ejectile = TMath::Sqrt(p_ejectile*p_ejectile+m_projectile*m_projectile);
            auto ke_ejectile = e_ejectile - m_projectile;
            double theta_cm=0,mmEx=0,J_L=0;
            std::tie(theta_cm,mmEx,J_L) = kine_2b(m_projectile,m_target,m_projectile,m_target,ke_projectile,ejectile.Theta(), ke_ejectile);


            hist(true,obj,dirname,"Ecm_particle[Elab]",
                 2500,0,10000,Elab,
                              2500,0,10000,Ecm_particle);
            hist(true,obj,dirname,"EdopplerSum",7000,0,21000,Ecm);
            hist(true,obj,dirname,"EdopplerParticleSum",7000,0,21000,Ecm_particle);

            if (ycor<22 && ycor >-13) {
              hist(true,obj,dirname,"A[X]_gateYcor",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
              // gate on prompt event timing peak
              if (tdiff>=252 && tdiff<=260) { // 9

                // cagra core energy summary vs. rcnp.GR_X(0)
                hist(true,obj,dirname,"EgamLab_vs_Ex",
                     600,5,65,mmEx,
                     6000,-2000,22000,core_hit.GetCorrectedEnergy());
                hist(true,obj,dirname,"EgamLab",7000,0,21000,core_hit.GetCorrectedEnergy());

                hist(true,obj,dirname,"EgamCM",7000,0,21000,Ecm);
                hist(true,obj,dirname,"EgamCMparticle",7000,0,21000,Ecm_particle);
                // cagra core energy summary vs. rcnp.GR_X(0)
                hist(true,obj,dirname,"EgamCMparticle_vs_Ex",
                     600,5,65,mmEx,
                     6000,-2000,22000,Ecm_particle);

                for (auto& other_hit : *cagra) {
                  if (other_hit.GetDetnum() == detector && other_hit.GetSystem() == 'Y' && other_hit.GetChannel() == 4) {
                    auto bgotdiff = core_hit.Timestamp() - other_hit.Timestamp();
                    hist(true,obj,dirname,"BGOTime_diff",5000,0,0, bgotdiff);
                    if (bgotdiff < 10 && bgotdiff > -10) {
                      bgo_hit = true;
                      break;
                      //hist(true,obj,dirname,"BGOTime_diff_accepted",1000,0,0,bgotdiff);
                    }

                  }
                }

                if (core_hit.GetSystem()=='Y') {
                  if (bgo_hit) {
                    hist(true,obj,dirname,"BGO_counts",10,0,10,3);
                  } else {
                    hist(true,obj,dirname,"BGO_counts",10,0,10,7);
                    hist(true,obj,dirname,"EgamCMparticleBGOveto_vs_Ex",
                         600,5,65,mmEx,
                         6000,-2000,22000,Ecm_particle);
                    hist(true,obj,dirname,"EgamCMparticleBGOveto",7000,0,21000,Ecm_particle);
                  }
                }

                bgo_hit = false;


                hist(true,obj,dirname,"MissingMassThetaCM",300,0,0,theta_cm);
                hist(true,obj,dirname,"MissingMassEx",600,5,65,mmEx);



                // raytracing
                double A=0,B=0;
                std::tie(A,B) = hit.Raytrace(true);
                hist(true,obj,dirname,"B[A]",500,0,0,A,500,0,0,B);
                hist(true,obj,dirname,"Bcor[Acor]",500,0,0,A,500,0,0,B);

                hist(true,obj,dirname,"Momentum",1000,2400,2800,hit.GetMomentum());


                // auto p_gamma = core_hit.GetMomentumVector(pos::core_only);
                // auto p_invariant = hit.ReconstructInvariant(p_gamma,true);
                // auto m_invariant = m_projectile + Li6Ex;
                // auto e_invariant = TMath::Sqrt(p_invariant.Mag()*p_invariant.Mag() + m_invariant*m_invariant);
                // auto ke_invariant = e_invariant - m_invariant;

                auto p_gamma = core_hit.GetMomentumVector(pos::core_only);
                auto p_invariant = hit.ReconstructInvariant(p_gamma,true);
                auto m_invariant = m_projectile + Li6Ex;
                auto e_invariant = TMath::Sqrt(p_invariant.Mag()*p_invariant.Mag() + m_invariant*m_invariant);
                auto ke_invariant = e_invariant - m_invariant;


                auto e_gamma = core_hit.GetCorrectedEnergy()/1000;
                hist(true,obj,"inv_tests","e_gamma",1000,0,0,e_gamma);
                auto p_ejectile = hit.GetEjectileVector(true); p_ejectile.SetMag(hit.GetMomentum());
                hist(true,obj,"inv_tests","p_ejectile",1000,0,0,p_ejectile.Mag());
                auto p_inv = p_ejectile + p_gamma;
                hist(true,obj,"inv_tests","p_inv",1000,0,0,p_inv.Mag());


                auto e_inv = std::sqrt(p_ejectile.Mag2() + std::pow(m_projectile,2)) + e_gamma;
                hist(true,obj,"inv_tests","e_inv",1000,0,0,e_inv);
                auto m_inv = std::sqrt(std::pow(e_inv,2)-p_inv.Mag2());
                hist(true,obj,"inv_tests","m_inv",1000,0,0,m_inv);
                auto ex_inv = m_inv - m_projectile;
                hist(true,obj,dirname,"inv_mass_ex",1000,0,10,ex_inv);




                hist(true,obj,dirname,"KE_invariant",1024,500,600,ke_invariant);
                hist(true,obj,dirname,"p_invariant",1024,2400,2800,p_invariant.Mag());
                auto theta_lab = p_invariant.Theta();
                hist(true,obj,dirname,"InvariantThetaLab",300,0,0,theta_lab);

                double theta_cm=0,Ex=0,J_L=0;
                std::tie(theta_cm,Ex,J_L) = kine_2b(m_projectile,m_target,m_invariant,m_target,ke_projectile,theta_lab, ke_invariant);


                if (Ecm_particle>=3400 && Ecm_particle<=3800) {
                  hist(true,obj,dirname+"_6LiEx_gate","ThetaCM",300,0,0.15,theta_cm);
                  hist(true,obj,dirname+"_6LiEx_gate","ReconstructedEx",600,5,65,Ex);
                  hist(true,obj,dirname+"_6LiEx_gate","ThetaLab[Ex]",600,5,65,Ex,300,0,0.15,theta_lab);
                  hist(true,obj,dirname+"_6LiEx_gate","EgamLab[Ex]",
                                    600,5,65,Ex,
                                    5500,0,22000,core_hit.GetCorrectedEnergy());

                  hist(true,obj,dirname+"_6LiEx_gate","MissingMassEx",600,5,65,mmEx);
                  hist(true,obj,dirname+"_6LiEx_gate","inv_mass_ex",200,0,10,ex_inv);


                  double gammafac = 1/(sqrt(1-pow(beta,2)));
                  auto Elab_gam = Li6Ex/(gammafac*(1-beta*p_gamma.CosTheta()));
                  auto adj_p_gamma = p_gamma;
                  adj_p_gamma.SetMag(Elab_gam);
                  auto adj_p_inv = p_ejectile + adj_p_gamma;
                  hist(true,obj,"inv_tests","adj_p_inv",1000,0,0,adj_p_inv.Mag());
                  auto adj_e_inv = std::sqrt(p_ejectile.Mag2() + std::pow(m_projectile,2)) + Elab_gam;
                  hist(true,obj,"inv_tests","adj_e_inv",1000,0,0,adj_e_inv);
                  auto adj_m_inv = std::sqrt(std::pow(adj_e_inv,2)-adj_p_inv.Mag2());
                  hist(true,obj,"inv_tests","adj_m_inv",1000,0,0,adj_m_inv);
                  auto adj_ex_inv = adj_m_inv - m_projectile;
                  hist(true,obj,"inv_tests","adj_inv_mass_ex",1000,0,10,adj_ex_inv);







                  if (ejectile.Theta() < 0.01) {
                    hist(true,obj,dirname+"_6LiEx_gate","ReconstructedEx_0_10_mrad",600,5,65,Ex);
                  }
                  hist(true,obj,dirname+"_6LiEx_gate","KE_projectile",100,0,0,ke_projectile);
                  hist(true,obj,dirname+"_6LiEx_gate","ejectile_theta",1000,-.1,.1,ejectile.Theta());

                  if (Ex >= 17.0 && Ex < 20.5) {
                    hist(true,obj,dirname+"_6LiEx_gate","B[A]_15.1MeVgate",25,-25,25,A,25,-80,80,B);
                  }

                } else if (Ecm_particle>=3900 && Ecm_particle <=4300) {
                  hist(true,obj,dirname+"_6LiEx_sideband","ThetaCM",300,0,0.15,theta_cm);
                  hist(true,obj,dirname+"_6LiEx_sideband","ReconstructedEx",600,5,65,Ex);
                  hist(true,obj,dirname+"_6LiEx_sideband","EgamLab[Ex]",
                       600,5,65,Ex,
                       5500,0,22000,core_hit.GetCorrectedEnergy());
                  if (ejectile.Theta() < 0.01) {
                    hist(true,obj,dirname+"_6LiEx_sideband","ReconstructedEx_0_10_mrad",600,5,65,Ex);
                  }
                  hist(true,obj,dirname+"_6LiEx_sideband","KE_projectile",100,0,0,ke_projectile);
                  hist(true,obj,dirname+"_6LiEx_sideband","ejectile_theta",1000,-.1,.1,ejectile.Theta());
                }
              } // end prompt timing gate

              // gate on random events side band of timing peak
              if (tdiff>=-155 && tdiff<=240) { // 396
                int detector = core_hit.GetDetnum();
                char core_leaf = core_hit.GetLeaf();
                string chan = core_leaf + std::to_string(core_hit.GetSegnum());
                auto crystal_id = TCagra::GetCrystalId(detector,core_leaf);

                // cagra core energy summary
                hist(true,obj,dirname,"rand_EgamLab",
                     6000,-2000,22000,core_hit.GetCorrectedEnergy(),
                     49,0,49,crystal_id);

                // cagra core energy summary vs. rcnp.GR_X(0)
                hist(true,obj,dirname,"rand_EgamLab_vs_Ex",
                     600,5,65,mmEx,
                     6000,-2000,22000,core_hit.GetCorrectedEnergy());



                auto p_gamma = core_hit.GetMomentumVector(pos::core_only);
                auto p_invariant = hit.ReconstructInvariant(p_gamma,true);
                auto m_invariant = m_projectile + Li6Ex;
                auto e_invariant = TMath::Sqrt(p_invariant.Mag()*p_invariant.Mag() + m_invariant*m_invariant);
                auto ke_invariant = e_invariant - m_invariant;

                double theta_cm=0,Ex=0,J_L=0;
                std::tie(theta_cm,Ex,J_L) = kine_2b(m_projectile,m_target,m_invariant,m_target,ke_projectile,p_invariant.Theta(), ke_invariant);


                if (Ecm_particle>=3400 && Ecm_particle<=3800) {
                  hist(true,obj,dirname+"_6LiEx_gate","rand_EgamLab[Ex]",
                                    600,5,65,Ex,
                                    5500,0,22000,core_hit.GetCorrectedEnergy());
                } else if (Ecm_particle>=3900 && Ecm_particle <=4300) {
                  hist(true,obj,dirname+"_6LiEx_sideband","rand_EgamLab[Ex]",
                       600,5,65,Ex,
                       5500,0,22000,core_hit.GetCorrectedEnergy());
                }







              } // end random timing gate

            } // end ycor cut
          } // end cagra loop
        } // end if cagra
      } // end rf PID gate
    } // ray id
  } // gr hits

}
void MakeCoincidenceHistograms(TRuntimeObjects& obj, TCagra& cagra, TGrandRaiden& gr) {
  MakeGRCorrections(obj,gr,&cagra,"PID_gated");


  for (auto& hit : gr) {

    auto& rcnp = hit.GR();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (rcnp.GR_RAYID(0) == 0) { // HARD CUT TO ONLY TAKE EVENTS WHICH HAD A GOOD GR RECONSTRUCTION !!!!!//
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto grtime = hit.GetTimestamp();
      auto rf = rcnp.GR_RF(0);
      auto x = rcnp.GR_X(0);
      //auto Ex = x*0.01074+6.872; // calibrated for 12C @ commissioning on 04.10.16

      // coincidence rate
      static ULong_t first_timestamp = grtime;
      if (first_timestamp) {
        auto rate = (grtime-first_timestamp)/1e8;
        //cout << grtime << " " << first_timestamp << endl;
        hist(false,obj,"Coincident","Rate",3000,0,30000, rate);
      }

      auto ycor = rcnp.GR_Y(0)+892.46*rcnp.GR_PH(0);
//    auto a = rcnp.GR_TH(0);


      // coincidence time difference
      for (auto& core_hit : cagra) {

        int detector = core_hit.GetDetnum();
        char core_leaf = core_hit.GetLeaf();
        string chan = core_leaf + std::to_string(core_hit.GetSegnum());
        auto crystal_id = TCagra::GetCrystalId(detector,core_leaf);
        bool bgo_hit = false;

        auto cagratime = core_hit.Timestamp();
        auto tdiff = cagratime-grtime;

        hist(false,obj,"Coincident","Diff_CAGRA_GR", 1000,-500,1500,cagratime-grtime);

        stream.str("");
        stream << "TimeDiff_" << detector << "_" << chan;
        hist(false,obj,"Coincident",stream.str().c_str(),1000,-500,1500,tdiff);


        // stream.str(""); stream << "Egam[tdiff]_" <<detector << "_" << chan;
        // hist(false,obj,stream.str(),5000,0,10000,core_hit.GetCorrectedEnergy(),1000,-500,1500,tdiff);


        // doppler reconstruction

        dirname = "ParticleGamma";
        auto Ecm = core_hit.GetDoppler(beta,pos::core_only);
        auto ejectile = hit.GetEjectileVector(true);
        auto Ecm_particle = core_hit.GetDoppler(beta,pos::core_only,ejectile);
        auto Elab = core_hit.GetCorrectedEnergy();

        hist(false,obj,dirname,"Ecm_particle[Elab]",
                          2500,0,10000,Elab,
                          2500,0,10000,Ecm_particle);
        hist(false,obj,dirname,"EdopplerSummary",
                          6000,-2000,22000,Ecm,
                          49,0,49,crystal_id);
        hist(false,obj,dirname,"EdopplerSum",7000,0,21000,Ecm);
        hist(false,obj,dirname,"EdopplerParticleSummary",
                          6000,-2000,22000,Ecm_particle,
                          49,0,49,crystal_id);
        hist(false,obj,dirname,"EdopplerParticleSum",7000,0,21000,Ecm_particle);


        hist(false,obj,dirname,"Y[B]cor",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,ycor);
        //hist(false,obj,dirname,"Y[B]cor",250,-0.1,0.1,rcnp.GR_PH(0),200,-100,100,ycor);

        if (ycor<22 && ycor >-13) {
          hist(false,obj,dirname,"Diff_CAGRA_GR", 1000,-500,1500,cagratime-grtime);
          hist(false,obj,dirname,"A[X]_gateYcor",1200,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
          //stream.str(""); stream << "TimeDiff_" << detector << "_" << chan;
          hist(false,obj,dirname,stream.str().c_str(),1000,-500,1500,tdiff);

          hist(false,obj,dirname+"_gates","A[X]_gateYcor_midgateAX", 300,-600,600,rcnp.GR_X(0),1000,-1,1,rcnp.GR_TH(0));
          hist(false,obj,dirname+"_gates","Diff_CAGRA_GR_midgateAX", 1000,-500,1500,tdiff);

          // gate on prompt event timing peak
          if (tdiff>=252 && tdiff<=260) {
            int detector = core_hit.GetDetnum();
            char core_leaf = core_hit.GetLeaf();
            string chan = core_leaf + std::to_string(core_hit.GetSegnum());
            auto crystal_id = TCagra::GetCrystalId(detector,core_leaf);

            // cagra core energy summary
            hist(false,obj,"CrystalEnergySummaryPrompt",
                              6000,-2000,22000,core_hit.GetCorrectedEnergy(),
                              49,0,49,crystal_id);

            // cagra core energy summary vs. rcnp.GR_X(0)
            hist(false,obj,"CrystalEnergySummaryPrompt_vs_X",
                              300,-600,600,rcnp.GR_X(0),
                              6000,-2000,22000,core_hit.GetCorrectedEnergy());
            hist(false,obj,"EgammaSumPrompt",7000,0,21000,core_hit.GetCorrectedEnergy());

            hist(false,obj,dirname,"EdopplerSummaryPrompt",
                              6000,-2000,22000,Ecm,
                              49,0,49,crystal_id);
            hist(false,obj,dirname,"EdopplerSum",7000,0,21000,Ecm);
            hist(false,obj,dirname,"EdopplerParticleSum",7000,0,21000,Ecm_particle);
            hist(false,obj,dirname,"EdopplerParticleSummaryPrompt",
                              6000,-2000,22000,Ecm_particle,
                              49,0,49,crystal_id);
            hist(false,obj,dirname,"EdopplerParticleSumPrompt",7000,0,21000,Ecm_particle);
            // cagra core energy summary vs. rcnp.GR_X(0)
            hist(false,obj,"EdopplerParticle_vs_X",
                              300,-600,600,rcnp.GR_X(0),
                              6000,-2000,22000,Ecm_particle);

            for (auto& other_hit : cagra) {
              if (other_hit.GetDetnum() == detector && other_hit.GetSystem() == 'Y' && other_hit.GetChannel() == 4) {
                auto bgotdiff = core_hit.Timestamp() - other_hit.Timestamp();
                hist(false,obj,dirname,"BGOTime_diff",5000,0,0, bgotdiff);
                if (bgotdiff < 10 && bgotdiff > -10) {
                  bgo_hit = true;
                  break;
                  //hist(false,obj,dirname,"BGOTime_diff_accepted",1000,0,0,bgotdiff);
                }

              }
            }

            if (core_hit.GetSystem()=='Y') {
              if (bgo_hit) {
                hist(false,obj,dirname,"BGO_counts",10,0,10,3);
              } else {
                hist(false,obj,dirname,"BGO_counts",10,0,10,7);
                hist(false,obj,"EdopplerParticleBGOveto_vs_X",
                                  300,-600,600,rcnp.GR_X(0),
                                  6000,-2000,22000,Ecm_particle);
                hist(false,obj,dirname,"EdopplerParticleSumBGOveto",7000,0,21000,Ecm_particle);
              }
            }




            bgo_hit = false;


            hist(false,obj,"EgammaSumPrompt",7000,0,21000,core_hit.GetCorrectedEnergy());

            if (rcnp.GR_X(0) > -370 && rcnp.GR_X(0) < -270) {
              hist(false,obj,"Doppler_15MeV","EdopplerSum",7000,0,21000,Ecm);
              hist(false,obj,"Doppler_15MeV","EdopplerParticleSum",7000,0,21000,Ecm_particle);
              hist(false,obj,"Doppler_15MeV","EdopplerParticleSummaryPrompt",
                                6000,-2000,22000,Ecm_particle,
                                49,0,49,crystal_id);
            }


            //hist(false,obj,dirname,"Momentum",1000,2500,2700,hit.GetMomentum());
            hist(false,obj,dirname,"Momentum",1000,2400,2800,hit.GetMomentum());
            auto p_gamma = core_hit.GetMomentumVector(pos::core_only);
            auto p_invariant = hit.ReconstructInvariant(p_gamma,true);
            auto m_invariant = m_projectile + Li6Ex;
            auto e_invariant = TMath::Sqrt(p_invariant.Mag()*p_invariant.Mag() + m_invariant*m_invariant);
            auto ke_invariant = e_invariant - m_invariant;
            hist(false,obj,dirname,"KE_invariant",1024,500,600,ke_invariant);
            hist(false,obj,dirname,"p_invariant",1024,2400,2800,p_invariant.Mag());

            double theta_cm=0,Ex=0,J_L=0;
            std::tie(theta_cm,Ex,J_L) = kine_2b(m_projectile,m_target,m_invariant,m_target,ke_projectile,p_invariant.Theta(), ke_invariant);


            if (Ecm_particle>=3400 && Ecm_particle<=3800) {
              hist(false,obj,"6LiEx_gate","ThetaCM",300,0,0.15,theta_cm);
              hist(false,obj,"6LiEx_gate","ReconstructedEx",200,0,100,Ex);
              hist(false,obj,"6LiEx_gate","Egam[Ex]",
                                200,0,100,Ex,
                                5500,0,22000,core_hit.GetCorrectedEnergy());

              if (ejectile.Theta() < 0.01) {
                hist(false,obj,"6LiEx_gate","ReconstructedEx_0_10_mrad",200,0,100,Ex);
              }
              hist(false,obj,"6LiEx_gate","KE_projectile",100,0,0,ke_projectile);
              hist(false,obj,"6LiEx_gate","ejectile_theta",1000,-.1,.1,ejectile.Theta());
            } else if (Ecm_particle>=3900 && Ecm_particle <=4300) {
              hist(false,obj,"6LiEx_sideband","ThetaCM",300,0,0.15,theta_cm);
              hist(false,obj,"6LiEx_sideband","ReconstructedEx",200,0,100,Ex);
              if (ejectile.Theta() < 0.01) {
                hist(false,obj,"6LiEx_sideband","ReconstructedEx_0_10_mrad",200,0,100,Ex);
              }
              hist(false,obj,"6LiEx_sideband","KE_projectile",100,0,0,ke_projectile);
              hist(false,obj,"6LiEx_sideband","ejectile_theta",1000,-.1,.1,ejectile.Theta());
            }




          }

          // gate on random events side band of timing peak
          if (tdiff>=212 && tdiff<=240) {
            int detector = core_hit.GetDetnum();
            char core_leaf = core_hit.GetLeaf();
            string chan = core_leaf + std::to_string(core_hit.GetSegnum());
            auto crystal_id = TCagra::GetCrystalId(detector,core_leaf);

            // cagra core energy summary
            hist(false,obj,"Summary","rand_CrystalEnergySummaryPrompt",
                              6000,-2000,22000,core_hit.GetCorrectedEnergy(),
                              49,0,49,crystal_id);

            // cagra core energy summary vs. rcnp.GR_X(0)
            hist(false,obj,"Summary","rand_CrystalEnergySummaryPrompt_vs_X",
                              300,-600,600,rcnp.GR_X(0),
                              6000,-2000,22000,core_hit.GetCorrectedEnergy());

          }
        }
      } // end cagra analysis
    }
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

void LoadRaytraceParams(size_t xdeg, size_t adeg, size_t ydeg, size_t bdeg) {
  static bool once = false;
  if (once) { return; }
  once = true;
  std::cout << "Loading Raytrace Parameters (" << xdeg << ", " << adeg << ", " << ydeg << ", " << bdeg << ")"<< std::endl;

  // sieve slit transformation coefficients
  // output from sieveslit.py
  TGrandRaidenHit::SetRaytraceParams(
    { // a fit parameters
      -16.46324544230171227,
        -0.02122094932284203,
        -0.00000175692614590,
        -0.00000000019248113,
        355.31589759280973340
    },
    { // b fit parameters

      // 3,1,1,1
      // 8.41837609679958909,
      //   -3607.17949748588671355,
      //   1.90171045803883931,
      //   36.26905443987910616,
      //   -377.31312093195958823,
      //   106966.41170441261783708,
      //   79.10896860670438002,
      //   -1129.69514732780385202,
      //   0.01795698529338584,
      //   -5.32080946070484728,
      //   -0.00532266760289275,
      //   0.18648872716102186,
      //   -0.14599997559969397,
      //   -13.07166978169148663,
      //   0.11687091398322876,
      //   -5.47943181532317070,
      //   -0.00007213194571246,
      //   -0.00593886782913727,
      //   -0.00002109368367002,
      //   -0.00027282209845065,
      //   0.00345501796013203,
      //   -1.11420619753071093,
      //   -0.00062148851409233,
      //   0.00670120998171145,
      //   -0.00000025345727323,
      //   0.00004319814907155,
      //   -0.00000000660292514,
      //   -0.00000085877678875,
      //   0.00000338879604594,
      //   -0.00086992922376217,
      //   -0.00000108370755612,
      //   0.00002594650904564

      // new 3,1,1,1 yayayay
      2.71579588833537944,
        -1746.85262809635059966,
        3.53327355622339967,
        8.15257960703027962,
        -126.19599585126820784,
        30571.41901346024314989,
        15.87646844285069925,
        -193.21437922048644964,
        0.00744103702906201,
        -0.89649047942650717,
        -0.00081059552774759,
        0.08457191245207539,
        -0.14916067272338165,
        -28.03494651190073128,
        0.04082308294975286,
        -4.20611225299385705,
        -0.00002531249611258,
        0.00343833829349489,
        -0.00000932843326045,
        0.00025781765104780,
        -0.00011719780912808,
        -0.27058010317939912,
        0.00000688556974219,
        -0.00434125793925212,
        -0.00000005785646263,
        0.00000036615425737,
        -0.00000001197465549,
        0.00000029289566844,
        0.00000048003272123,
        -0.00024986622710858,
        0.00000003669306537,
        0.00001138850753476

      // 3,1,1,0
      // -1.15083365480030220,
      //   5.35167483518167764,
      //   -49.39598493923793399,
      //   -13.62981275316127672,
      //   0.00846642194952134,
      //   0.00170202243207763,
      //   -0.23158408095867616,
      //   0.04798202827122525,
      //   0.00001434639259220,
      //   -0.00001595369862982,
      //   -0.00151696633766249,
      //   0.00037967592175940,
      //   -0.00000000756710322,
      //   -0.00000002285216941,
      //   -0.00000157023814891,
      //   0.00000057104358734

      // 3,2,1,0
      // 3.06120625217588982,
      //   5.17162599163630610,
      //   -340.41455859951196317,
      //   -0.40399534244789076,
      //   3239.25987934644717825,
      //   -149.09742406499594836,
      //   0.01682807972099215,
      //   -0.00236637879905502,
      //   -0.23857337886064786,
      //   0.32276769594154087,
      //   -3.34930433942084882,
      //   -3.11597380176068750,
      //   -0.00003218483418950,
      //   -0.00002974118837506,
      //   0.00085695301413138,
      //   0.00082206417044957,
      //   -0.02121327508777294,
      //   -0.00227951088543225,
      //   -0.00000010461966499,
      //   -0.00000003539426411,
      //   0.00000030870446224,
      //   0.00000053596095634,
      //   0.00001007076634915,
      //   0.00000489147736253

      //3,2,2,0
      // 2.69530522292509289,
      //   5.18785247049013520,
      //   0.01704192844817459,
      //   -274.55663145212349718,
      //   2.61954079020696495,
      //   -3.37052929740243146,
      //   2236.08165691897875149,
      //   -208.09810111366567753,
      //   51.79953653178802142,
      //   0.01811466414640416,
      //   -0.00235235051422705,
      //   0.00000322586654669,
      //   -0.01012397311721290,
      //   0.34896623475566130,
      //   -0.01154570528127622,
      //   -8.63098403252137025,
      //   -3.64132864335810646,
      //   0.20776318449515729,
      //   -0.00000848892936632,
      //   -0.00003050056281987,
      //   -0.00000066825232710,
      //   -0.00108215986107756,
      //   0.00080272100286088,
      //   0.00007065632349954,
      //   0.00740905538053109,
      //   -0.00118144246444490,
      //   -0.00111538332183897,
      //   -0.00000006973083691,
      //   -0.00000003702694993,
      //   -0.00000000120724884,
      //   -0.00000453670516506,
      //   0.00000037821978831,
      //   0.00000017792107931,
      //   0.00009901865475540,
      //   0.00000986283664030,
      //   -0.00000305991562225

        },
    xdeg, // polynomial degree in fit for x
    adeg,
    ydeg,
    bdeg
    );
}
