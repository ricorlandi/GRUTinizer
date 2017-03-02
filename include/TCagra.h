#ifndef TCagra_H
#define TCagra_H

#include <TObject.h>
#include <TMath.h>

#include "TDetector.h"
#include "TCagraHit.h"

#include "TVector3.h"

#include <vector>

class TCagra : public TDetector {

public:
  TCagra();
  ~TCagra();

  virtual void Copy(TObject& obj) const;
  virtual void Print(Option_t *opt = "") const;
  virtual void Clear(Option_t *opt = "");

  virtual void          InsertHit(const TDetectorHit& hit);
  virtual TDetectorHit& GetHit(int i)            { return cagra_hits.at(i); }

  const TCagraHit& GetCagraHit(int i) { return cagra_hits.at(i); }
  void PrintHit(int i){ cagra_hits.at(i).Print(); }

  static TVector3 GetSegmentPosition(int slot, char core, int seg);
  static int GetCrystalId(int slot, char core);
  static int GetSegmentId(int slot, char core, int segnum, char system);
  static int NumSegments() { return segment_ids.size(); }

  // Allows for looping over all hits with for(auto& hit : cagra) { }
  std::vector<TCagraHit>::iterator begin() { return cagra_hits.begin(); }
  std::vector<TCagraHit>::iterator end() { return cagra_hits.end(); }
  size_t size() const { return cagra_hits.size(); }

private:
  virtual int BuildHits(std::vector<TRawEvent>& raw_data);
  static void LoadDetectorPositions();

  std::vector<TCagraHit> cagra_hits;
  static std::map<int,TVector3> detector_positions;
  static std::map<int,size_t> crystal_ids;
  static std::map<int,size_t> segment_ids;
  static bool positions_loaded;

  ClassDef(TCagra,1);
};


unsigned int GetRandomCAGRAChannel();


#endif
