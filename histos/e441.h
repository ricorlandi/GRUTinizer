#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <algorithm>
#include <cassert>
#include <cmath>


#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TObject.h"
#include "TFile.h"
#include "TCutG.h"
#include "TVector3.h"
#include "TPreserveGDirectory.h"
#include "TCagra.h"
#include "TGrandRaiden.h"
#include "TANLEvent.h"
#include "GValue.h"
#include "TNucleus.h"
#include <tuple>

#define BAD_NUM -441441
#define PRINT(x) std::cout << #x" = " << x << std::endl
#define STR(x) #x << " = " <<() x

using namespace std;

void MakeCAGRAHistograms(TRuntimeObjects&, TCagra&);
void MakeGrandRaidenHistograms(TRuntimeObjects&, TGrandRaiden&);
void MakeCoincidenceHistograms(TRuntimeObjects&, TCagra&, TGrandRaiden&);
void LoadCuts();
void InitTarget();

enum class Target { n12C=1, n56Fe=2, n24Mg=3, n93Nb=4, n93Nb2=5  };

static const Target target = Target(GValue::Value("Target"));

static string name;
static string dirname="";
static stringstream stream;

static TNucleus li6("6Li");
static TNucleus c12("12C");
static TNucleus mg24("24Mg");
static TNucleus fe56("56Fe");
static TNucleus nb93("93Nb");
static TNucleus sn124("124Sn");

static const double a_offset = std::isnan(GValue::Value("A_Offset")) ? 0 : GValue::Value("A_Offset");
static const double b_offset = std::isnan(GValue::Value("B_Offset")) ? 0 : GValue::Value("B_Offset");

static const double m_projectile = li6.GetMass();
static const double eloss = GValue::Value("EnergyLoss");
static const double ke_projectile = GValue::Value("BeamKineticE") - eloss;
static const double e_projectile = m_projectile + ke_projectile;
static const double p_projectile = TMath::Sqrt(e_projectile*e_projectile-m_projectile*m_projectile);
static const double beta = p_projectile/e_projectile;
static const TVector3 projectile(0.,0.,p_projectile);


//static TNucleus li6ex(3,3,li6.GetMass()+Li6Ex,"6Li*");

///=============Two Body Kinematics===========
//Kinematics based on the Mandelstam invariant variables
double omega(double x, double y, double z){
  return sqrt(x*x + y*y + z*z -2*x*y -2*y*z -2*x*z);
}

std::tuple<double,double,double> kine_2b(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject){
  // m1(projectile) - m2(target) - m3(ejectile) - and m4(recoil)

  double Et1 = K_proj + m1;
  double Et2 = m2;
  double Et3 = K_eject + m3;
  double Et4  = Et1 + Et2 - Et3;
  double m4_ex, Ex, theta_cm;
  double s,t,u; //---Mandelstam variables
  double p1, p3;
  double J_LtoCM; //jacobian Lab to CM

  s = pow(m1,2) + pow(m2,2) +2*m2*Et1;
  u = pow(m2,2) + pow(m3,2) - 2*m2*Et3;

  m4_ex = sqrt(  (cos(thetalab) * omega(s,pow(m1,2),pow(m2,2)) * omega(u,pow(m2,2),pow(m3,2)) - (s - pow(m1,2) - pow(m2,2))*(pow(m2,2) + pow(m3,2) - u) )/(2*pow(m2,2)) + s + u - pow(m2,2)  );
  Ex = m4_ex - m4;

  t =   pow(m2,2) + pow(m4_ex,2) - 2*m2*Et4;

  //for normal kinematics
  theta_cm = acos( ( pow(s,2) +s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2)) + (pow(m1,2) - pow(m2,2))*(pow(m3,2) - pow(m4_ex,2)) )/( omega(s,pow(m1,2),pow(m2,2))*omega(s,pow(m3,2),pow(m4_ex,2))) ) ;

  //for inverse kinematics Note: this angle corresponds to the recoil
  //theta_cm = TMath::Pi() - acos( ( pow(s,2) +s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2)) + (pow(m1,2) - pow(m2,2))*(pow(m3,2) - pow(m4_ex,2)) )/( omega(s,pow(m1,2),pow(m2,2))*omega(s,pow(m3,2),pow(m4_ex,2))) ) ;

  p1= sqrt(pow(Et1,2)-pow(m1,2));
  p3 = sqrt(pow(Et3,2)-pow(m3,2));

  J_LtoCM = abs( ((omega(s,pow(m1,2),pow(m2,2))*omega(s,pow(m3,2),pow(m4,2)))/(4*s*p1*p3))*(1.+Et1/m2 - cos(thetalab)*(Et3*p1)/(m2*p3)) );

  return std::tuple<double,double,double>(theta_cm,Ex,J_LtoCM);
}
///=============Brho to Kinetic energy transform===========
double BrhoToTKE(double  brho, double  mass, double Z) {
  //     convert brho (magnetic rigidity) to tke (total kinetic energy)
  //     input:  brho (Tm)
  //             mass (MeV)
  //             Z = q/e (no-dim)
  //     output: TKE (MeV)
  return mass * (sqrt(1. + TMath::Power(((1.e2 * brho * TMath::C() / 1.e8 * Z) / mass), 2)) - 1.);
}

void LoadRaytraceParams(size_t xdeg=2, size_t adeg=2, size_t ydeg=1, size_t bdeg=1);

double VectorAverage(const std::vector<short int>& vec) {
  double sum = 0;
  for (auto const& el : vec) {
    sum += el;
  }
  return (float)sum/vec.size();
}

double VectorStdDev(const std::vector<short int>& vec) {
  double avg = VectorAverage(vec);
  double sum = 0;
  for (auto const& el : vec) {
    sum += std::pow((el-avg),2);
  }
  return sum/vec.size();
}

std::pair<double,double> VectorStats(const std::vector<short int>& vec) {
  double avg = VectorAverage(vec);
  double sum = 0;
  for (auto const& el : vec) {
    sum += std::pow((el-avg),2);
  }
  return std::pair<double,double>(avg,sqrt(sum/vec.size()));
}


void DrawAverageTrace(TCagraHit& core_hit) {
// Average trace analysis
  //////////////////
  static std::vector<std::vector<short int>> traces;
  auto trace = core_hit.GetSanitizedTrace();
  if(trace->size() > 50) {
    traces.push_back(*trace);
  }

  if (traces.size() >= 100) {
    // do a transpose to convert to vector<short int> for each bin
    std::vector<std::vector<short int>> bins(traces[0].size());
    size_t prev_size = traces[0].size();
    for (auto& atrace : traces) {
      assert(atrace.size() == prev_size);
      prev_size = atrace.size();
      for (auto nbin = 0u; nbin<atrace.size(); nbin++) {
        bins[nbin].push_back(atrace[nbin]);
      }
    }
    // for each bin
    // build avg histo with errors
    std::vector<double> avgs, errors;
    for (auto& bin : bins) {
      double avg = 0, stddev = 0;
      std::tie(avg, stddev) = VectorStats(bin);
      avgs.push_back(avg);
      errors.push_back(stddev);
    }
    TCagraHit::DrawTrace(avgs,errors);
    std::cin.get();

    traces.clear();
  }
}


class RateOffset {
public:
  RateOffset(ULong_t interval, Target target)
    : interval(interval) {
    switch (target) {
    case Target::n12C:
    case Target::n56Fe:
    case Target::n24Mg:
     num_crystals = 1;
     break;
    case Target::n93Nb:
    case Target::n93Nb2:
      num_crystals = 49;
      break;
    default:
      throw std::runtime_error("Unknown target");
    }
    counts.reserve(num_crystals);
    rate.reserve(num_crystals);
    std::fill(counts.begin(), counts.end(), 0);
    std::fill(rate.begin(), rate.end(), 0);
  }


  inline ULong_t Update(const ULong_t& ts, size_t crystal_id=0) {
    if (crystal_id == (size_t)-1) { return first_ts; }

    if (first_ts<=1e6) {
      first_ts = ts;
      std::cout << "Timestamp: " << first_ts << "\n" << std::endl;
    } else {

      assert(crystal_id < num_crystals);

      // calculate rate every interval (js)
      counts[crystal_id]++;
      if ((ts-first_ts)*10/1.0e9 > interval_end) {
        interval_end += interval;

        for(int i=0; i<num_crystals; i++) { rate[i] = counts[i]/interval; }
        for(int i=0; i<num_crystals; i++) { counts[i] = 0; }
        //rate = 8*counts/interval;
        //counts = 0;
      }


      switch (target) {
      case Target::n12C:
        if(rate[0]*8 > 10000) { ecal_offset = -11.511178+0.000120007*(rate[crystal_id])*8; }
        break;
      case Target::n56Fe:
        if(rate[crystal_id] > 10000) { ecal_offset = 0; }
        break;
      case Target::n24Mg:
        if(rate[0]*8 > 10000) { ecal_offset = -16.921206+0.000136664*(rate[crystal_id]*8); }
        break;
      case Target::n93Nb:
        assert(num_crystals==49);
	ecal_offset = ratecorr_intercepts[crystal_id]+ratecorr_slopes[crystal_id]*rate[crystal_id];
        break;
      case Target::n93Nb2:
        assert(num_crystals==49);
	if (crystal_id == 44) {
	  ecal_offset = -2.942630+0.168145*rate[crystal_id];
	} else {
          ecal_offset = ratecorr_intercepts[crystal_id]+ratecorr_slopes[crystal_id]*rate[crystal_id];
	}
        break;
      default:
        throw std::runtime_error("Unknown target");
      }




    }
    return first_ts;
  }
  inline Double_t operator()() { return ecal_offset; }


private:
  ULong_t  first_ts = 0;
  const ULong_t interval = 1; // interval for rate calculation
  size_t num_crystals;
  ULong_t  interval_end = 1;
  std::vector<Double_t> counts;
  std::vector<Double_t> rate;
  Double_t ecal_offset = 0;


  std::vector<double> ratecorr_slopes = {0.080527,0.083746,0.074751,0.077786,
					 0.076732,0.098430,0.093351,0.083438,
					 0.085293,0.101123,0.096509,0.083015,
					 0.082756,0.087294,0.083884,0.075692,
					 0.080683,0.085537,0.080275,0.087235,
					 0,0,0,0,
					 0.131977,0.086933,0.088750,0.105513,
					 0.094898,0.088654,0.087073,0.097078,
					 0.075507,0.077302,0.075756,0.073714,
					 0.075501,0.074556,0.074845,0.071820,
					 0.082464,0.080693,0.082462,0.083240,
					 0.065678,0.081520,0.081546,0.081205,0};

  std::vector<double> ratecorr_intercepts = {17.174516,18.352194,15.124436,13.700136,
					     16.116238,23.141675,23.135089,20.002620,
					     18.941009,22.406043,21.815457,17.134956,
					     20.875168,19.613706,21.087828,18.793253,
					     16.577259,13.321167,15.279376,14.972158,
					     0,0,0,0,
					     22.339411,17.500503,17.637294,22.298484,
					     17.500876,18.257787,15.642874,21.772555,
					     12.990928,11.156014,12.811979,13.767984,
					     11.158173,11.734640,12.296375,11.779592,
					     22.084610,20.376215,20.625643,19.627033,
					     6.867173,21.474127,20.714631,19.395711,0};
};


///=============Mass definition============
#define aum 931.494043    // in MeV
#define e_mass  548.57990945e-6  // in aum
#define n_amu 1.0086654
#define n_ch 0
#define n_mass n_amu*aum //in MeV
#define Li6Ex 3.5628 // MeV
#define H1  1.00782503223 // in aum
#define H1_ch 1
#define H1_mass (H1 -  H1_ch*e_mass )*aum //in MeV
#define H3  3.01604927791 // in aum
#define H3_ch 1
#define H3_mass (H3 -  H3_ch*e_mass )*aum //in MeV
#define He3  3.01602932008 // in aum
#define He3_ch 2
#define He3_mass (He3 -  He3_ch*e_mass )*aum //in MeV
#define C12  12.0000000000 // in aum
#define C12_ch 6
#define C12_mass (C12 -  C12_ch*e_mass )*aum //in MeV
#define C14 14.003255
#define C14_ch 6
#define C14_mass (C14 -  C14_ch*e_mass )*aum //in MeV
#define N14 14.0030869
#define N14_ch 7
#define N14_mass (N14 -  N14_ch*e_mass )*aum //in MeV
#define N12 12.018624
#define N12_ch 7
#define N12_mass (N12 -  N12_ch*e_mass )*aum //in MeV
#define B12  12.014352658 // in aum
#define B12_ch 5
#define B12_mass (B12 -  B12_ch*e_mass )*aum //in MeV
#define O16  15.9949146195// in aum
#define O16_ch  8
#define O16_mass (O16 -  O16_ch*e_mass )*aum //in MeV
#define N16    16.00610192 // in aum
#define N16_ch  7
#define N16_mass (N16 -  N16_ch*e_mass )*aum //in MeV
#define Sc45  44.955908275 // in aum
#define Sc45_ch 21
#define Sc45_mass (Sc45 -  Sc45_ch*e_mass )*aum //in MeV
#define Ca45  44.956186350 // in aum
#define Ca45_ch 20
#define Ca45_mass (Ca45 -  Ca45_ch*e_mass )*aum //in MeV
#define Ti46 45.952627718 // in aum
#define Ti46_ch  22
#define Ti46_mass (Ti46 -  Ti46_ch*e_mass )*aum //in MeV
#define Sc46 45.955168257 // in aum
#define Sc46_ch  21
#define Sc46_mass (Sc46 -  Sc46_ch*e_mass )*aum //in MeV
#define Kr86  85.91061062693 // in aum
#define Kr86_ch 36
#define Kr86_mass (Kr86 -  Kr86_ch*e_mass )*aum //in MeV
#define Br86  85.918805433 // in aum#define O16  15.9949146195// in aum
#define O16_ch  8
#define O16_mass (O16 -  O16_ch*e_mass )*aum //in MeV
#define N16    16.00610192 // in aum
#define N16_ch  7
#define N16_mass (N16 -  N16_ch*e_mass )*aum //in MeV
#define Br86_ch 35
#define Br86_mass (Br86 -  Br86_ch*e_mass )*aum //in MeV
#define Brho2  2.3215  //in Tm
