#ifndef TGRANDRAIDENHIT_H
#define TGRANDRAIDENHIT_H

#include "TDetector.h"
#include "TDetectorHit.h"
#include "GValue.h"
#include <vector>

#ifdef RCNP
#include "RCNPEvent.h"
#else
// ---------------------
class RCNPEvent : public TObject { };
// ---------------------
#endif

struct LaBrHit;

class TGrandRaidenHit : public TDetectorHit {
  public:
    TGrandRaidenHit(); // TODO: move to private
    TGrandRaidenHit(RCNPEvent& rcnpevent);
    TGrandRaidenHit(const TGrandRaidenHit& gr);
    ~TGrandRaidenHit();

    virtual void Copy(TObject& obj) const;
    //virtual void Compare(TObject &obj) const;
    virtual void Print(Option_t *opt = "") const;
    virtual void Clear(Option_t *opt = "");

    void     BuildFrom();

    TVector3 GetEjectileVector(bool apply_offsets);
    TVector3 ReconstructInvariant(const TVector3& gamma, bool apply_offsets);
    double   GetMomentum();

    std::pair<double,double> Raytrace(bool apply_offsets);
    const std::vector<LaBrHit>& GetLaBr() { return labr_hits; }
    const Double_t& GetMeanPlastE1() { return madc1; }
    const Double_t& GetMeanPlastE2() { return madc2; }
    const Double_t& GetMeanPlastE3() { return madc3; }
    const Double_t& GetMeanPlastPos1() { return tpos1; }
    const Double_t& GetMeanPlastPos2() { return tpos2; }
    const Double_t& GetMeanPlastPos3() { return tpos3; }
    const Long_t& GetTimestamp() { return Timestamp; }
    RCNPEvent& GR() { return rcnp; }
    static void SetRaytraceParams(std::vector<double> apar, std::vector<double> bpar,
                                  size_t xdeg=3, size_t adeg=1, size_t ydeg=1, size_t bdeg=1);

    static std::pair<Double_t,Double_t> GetAngleOffsets() {
        if (std::isnan(TGrandRaidenHit::aoffset) || std::isnan(TGrandRaidenHit::boffset)) {
            TGrandRaidenHit::aoffset = std::isnan(GValue::Value("A_Offset")) ? 0 : GValue::Value("A_Offset");
            TGrandRaidenHit::boffset = std::isnan(GValue::Value("B_Offset")) ? 0 : GValue::Value("B_Offset");
        }
        return std::pair<Double_t,Double_t>(TGrandRaidenHit::aoffset,TGrandRaidenHit::boffset);
    }
private:

    // angle reconstruction (energy reconstruction still needs implementation)
    // apply sieveslit transform to raytrace from focal plane to target position
    std::pair<double,double> raytrace(double x, double a, double y, double b, bool apply_offsets = false);

    // Double_t ADC[4];
    // Double_t RF;
    // Double_t QTCLead;

    std::vector<LaBrHit> labr_hits;
    Double_t madc1;
    Double_t madc2;
    Double_t madc3;
    Double_t tpos1;
    Double_t tpos2;
    Double_t tpos3;
    Long_t Timestamp;
    RCNPEvent rcnp;
    Double_t excitation_energy;
    Double_t momentum;

    static unsigned int xdegree;  //!
    static unsigned int adegree; //!
    static unsigned int ydegree; //!
    static unsigned int bdegree; //!
    static std::vector<double> acoefs; //!
    static std::vector<double> bcoefs; //!
    static Double_t aoffset; //!
    static Double_t boffset; //!



    TVector3 vector;


    ClassDef(TGrandRaidenHit,1);
};

struct LaBrHit : public TDetectorHit {
    Int_t channel;
    Double_t qtc_le;
    Double_t qtc_tr;
    Double_t width;
    ClassDef(LaBrHit,1);
};


#endif
