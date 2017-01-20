#ifndef TCagraHIT_H
#define TCagraHIT_H

#include "TDetectorHit.h"

enum pos { seg_only, core_only, both };

class TCagraHit : public TDetectorHit {
  public:
    TCagraHit();
    ~TCagraHit();

    virtual void Copy(TObject& obj) const;
    virtual void Print(Option_t *opt = "") const;
    virtual void Clear(Option_t *opt = "");
    //virtual void Draw(Option_t* opt = "");

    virtual Int_t Charge() const;
    Float_t GetCharge() const;
    void SetCharge(Float_t ch) { fCharge=ch; }

    int GetDetnum() const;
    char GetLeaf() const;
    int GetMainSegnum() const;
    int GetSegnum() const;
    char GetSystem() const;

    bool HasCore() const;

    void AddSegmentHit(const TCagraHit& hit) { fSegments.push_back(hit); }
    unsigned int GetNumSegments() const { return fSegments.size(); }
    TCagraHit& GetSegment(int i) { return fSegments.at(i); }
    void RemoveSegment(int i) { fSegments.erase(fSegments.begin()+i); }
    void ClearSegments() { fSegments.clear(); }
    unsigned long GetSegmentTimestamp() {
        if(fSegments.size()){
            return fSegments[0].Timestamp();
        } else {
            return -1;
        }
    }
    void PrintChannel() const;

    TCagraHit& MakeSegmentByAddress(unsigned int address);

    int GetBoardID() const;
    int GetChannel() const;

    TVector3 GetPosition(pos opt = pos::both, bool apply_array_offset = false) const;
    TVector3 GetMomentumVector(pos opt = pos::both, bool apply_array_offset = false);

    double GetDoppler(double beta, pos opt = pos::both,
                      const TVector3& particle_vec = TVector3(0,0,1),
                      const TVector3& cagra_offset = TVector3(0,0,0));

    std::vector<TCagraHit>::iterator begin() { return fSegments.begin(); }
    std::vector<TCagraHit>::iterator end() { return fSegments.end(); }

    void SetDiscTime(const Double_t t) { time = t; }
    Double_t GetDiscTime() { return time; }
    void SetPrevDiscTime(const Double_t t) { prev_time = t; }
    Double_t GetPrevDiscTime() { return prev_time; }

    Double_t GetCorrectedEnergy(Double_t asym_bl=0.);
    void SetPreRise(Double_t prerise) { prerise_energy = prerise; }
    void SetPostRise(Double_t postrise) { postrise_energy = postrise; }
    Double_t GetPreRise() { return prerise_energy; }
    Double_t GetPostRise() { return postrise_energy; }
    Short_t GetPrevPostRiseBeginSample() { return prev_postrise_begin_sample; }
    Short_t GetPreRiseBeginSample() { return prerise_begin; }
    Short_t GetPreRiseEndSample() { return prerise_end; }
    void SetPrevPostRiseBeginSample(UShort_t val) { prev_postrise_begin_sample = val; }
    void SetPreRiseBeginSample(UShort_t val) { prerise_begin = val; }
    void SetPreRiseEndSample(UShort_t val) { prerise_end = val; }

    void SetFlags(UShort_t fl) { flags = fl; }
    const UShort_t& GetFlags() const { return flags; }
    void SetBaseSample(UShort_t base) { base_sample = base; }
    const UShort_t& GetBaseSample() const { return base_sample; }
    void SetSampledBaseline(UShort_t base) { sampled_baseline = base; }
    const UShort_t& GetSampledBaseline() const { return sampled_baseline; }
    std::vector<Short_t>* GetTrace(int segnum=0);
    std::vector<Short_t>* GetSanitizedTrace(int segnum=0);
    void SetTrace(std::vector<Short_t>& trace);
    void DrawTrace(int segnum, bool draw_baseline=false);
    static void DrawTrace(const std::vector<double>& trace, const std::vector<double>& error);
    void DrawTraceBaseline(int segnum);

    void DrawTraceSamples(int segnum);
    double GetTraceHeight(size_t size=10) const;
    double GetTraceHeightDoppler(double beta,const TVector3& vec = TVector3(0,0,1)) const;
    Double_t GetTraceHeightPZ(Double_t asym_bl=0, size_t size=10);
    Double_t GetTraceEnergy(const UShort_t& a,const UShort_t& b,UShort_t x = 0,UShort_t y=0) const;
    Double_t GetTraceBaseline();
    void SetTimingMarks(std::vector<UShort_t>& marks) { fMarks.swap(marks); }

    Double_t GetBaselineExpCorr(int segnum=0);
    Double_t GetBaselineExpCorrFast(int segnum=0);
    void DrawBaselineExponential(int segnum=0);
    inline void MarkAsSegmentHit() { is_segment = true; }
    bool IsSegment() { return is_segment; }

  private:
    bool is_segment;

    Double_t charge;

    std::vector<Short_t> fTrace;
    std::vector<UShort_t> fMarks;
    std::vector<TCagraHit> fSegments;
    bool baseline_fitted;

    Double_t time;
    Double_t prev_time;

    UShort_t flags;

    Double_t prerise_energy;
    Double_t postrise_energy;
    UShort_t base_sample; // running sum
    UShort_t sampled_baseline; // event by event

    // trace points for base line determination
    Short_t prev_postrise_begin_sample;
    Short_t prerise_begin;
    Short_t prerise_end;
    Double_t corrected_energy;

    Double_t fit_params[2]; //!

    inline Double_t transform_trace_point(Double_t point) {
        //auto adc = (point < 0) ? point + std::pow(2,14) : point;
        //return adc;
        auto adc = (point < 0) ? point + std::pow(2,14) : point;
        adc = 8192 - adc;
        return adc + 1000;
    }

    inline Double_t inverse_transform_trace_point(Double_t point) {
        point -= 1000;
        return 8192 - point;
    }

  ClassDef(TCagraHit,2);
};


#endif
