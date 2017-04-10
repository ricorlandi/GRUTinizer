#include "TCagra.h"

#include <iostream>
#include <fstream>

#include "TANLEvent.h"
#include "TChannel.h"

#include "TGEBEvent.h"

#include <random>
#include <chrono>
#include <set>
#include <unordered_map>
#include <cassert>

std::map<int,TVector3> TCagra::detector_positions;
std::map<int,size_t> TCagra::crystal_ids;
std::map<int,size_t> TCagra::segment_ids;
bool TCagra::positions_loaded = false;

TCagra::TCagra(){
  Clear();
  LoadDetectorPositions();

  static bool load_segments = true;
  if (load_segments) {
    load_segments = false; // only perform this loop assigment once

    std::vector<int> detnums(16);
    std::iota(detnums.begin(),detnums.end(),1);
    std::vector<int> yale_seg = { 1, 2, 3 };
    std::vector<int> imp_core = { 1, 2, 3, 4 };
    std::vector<int> imp_seg = { 1, 2, 3, 4 };

    // TODO: Fix.
    // Currently hard code detectors 1-14 as yale
    // 15 and 16 as IMP. Probably not true for _all_ exp.

    int segment_count = 0;

    // yale type
    for (auto i=0u; i<detnums.size()-2; i++) {
      for (auto& seg : yale_seg) {
        int index = detnums[i]*10000 + seg;
        if (segment_ids.count(index) == 0) {
          segment_ids[index] = segment_count++;
        } else {
          throw std::runtime_error("Logic error in building segment ids.");
        }
      }
    }
    //imp type
    for (auto i=detnums.size()-2; i<detnums.size(); i++) {
      for (auto& core : imp_core) {
        for (auto& seg : imp_seg) {
          int index = detnums[i]*10000 + core*100 + seg;
          if (segment_ids.count(index) == 0) {
            segment_ids[index] = segment_count++;
          } else {
            throw std::runtime_error("Logic error in building segment ids.");
          }
        }
      }
    }



  }
}

TCagra::~TCagra() {

}

void TCagra::Copy(TObject& obj) const {
  TDetector::Copy(obj);

  TCagra& detector = (TCagra&)obj;
  detector.cagra_hits = cagra_hits;
}

void TCagra::InsertHit(const TDetectorHit& hit){
  cagra_hits.emplace_back((TCagraHit&)hit);
  fSize++;
}

unsigned int GetRandomCAGRAChannel(const int detnum, const char leaf) {
  // Board ids  101-108 90deg Yale clovers + BGO
  //            109,110 135deg Yale clovers + BGO
  //            111-114 135deg IMP clovers + BGO
  // Channel ids
  //            Yale    1 digitizer per clover
  //                    0-3 core signals
  //                    4 - L, 5 - M, 6 - R
  //                    7 - BGO
  //
  //            IMP     2 digitizers per clover
  //                    0,5 - core signals
  //                    1-4,6-9, 4 signals
  static std::mt19937 mt(std::chrono::system_clock::now().time_since_epoch().count());
  static std::uniform_int_distribution<> board(101,114);
  static std::uniform_real_distribution<float> channel(0,1);


  unsigned int board_id = 0;
  unsigned int chan_id = 0; // channel(mt);
  unsigned int address = 0xeeeeeeee;

  static char parent_type = 0;

  if (detnum == -1 && leaf == -1) {
    // asking for a random detector (clover and leaf)
    do {
      board_id = board(mt);
      if (board_id > 110) {
        chan_id = channel(mt)*9;
      } else {
        chan_id = channel(mt)*7;
      }
      address = ((1<<24) + (board_id << 8) + chan_id);


      auto _system = *TChannel::Get(address)->GetSystem();
      //auto _detnum = TChannel::Get(address)->GetArrayPosition();
      auto _leaf = *TChannel::Get(address)->GetArraySubposition();
      auto _segnum = TChannel::Get(address)->GetSegment();
      // std::cout << "Parent: ";
      // std::cout << _system << " - ";
      // std::cout << _detnum << " ";
      // std::cout << _leaf << " ";
      // std::cout << _segnum << std::endl;
      if (_segnum == 0 && _leaf != 'X') {
        parent_type = _system;
        break;
      }

    } while (true);
    //std::cin.get();



  } else {

    // asking for a segment of specified detector
    do {
      board_id = board(mt);
      if (board_id > 110) {
        chan_id = channel(mt)*9;
      } else {
        chan_id = channel(mt)*7;
      }
      address = ((1<<24) + (board_id << 8) + chan_id);

      //auto _system = *TChannel::Get(address)->GetSystem();
      auto _detnum = TChannel::Get(address)->GetArrayPosition();
      auto _leaf = *TChannel::Get(address)->GetArraySubposition();
      auto _segnum = TChannel::Get(address)->GetSegment();
      // std::cout << "  Segment: ";
      // std::cout << _system << " - ";
      // std::cout << _detnum << " ";
      // std::cout << _leaf << " ";
      // std::cout << _segnum;
      // std:: cout << "  Need to match parent_type: " << parent_type << std::endl;

      if (parent_type == 'Y') {
        if (_detnum == detnum && _segnum != 0) { break; }
      }else {
        if (_detnum == detnum && _leaf == leaf && _segnum != 0) { break; }
      }
    } while (true);

    //std::cout << "Found!" << std::endl;
    //std::cin.get();

  }



  //std::cout << board_id << " " << chan_id << std::endl;
  //std::cout << std::hex << ((1<<24) + (board_id << 8) + chan_id) << std::endl;
  return address;


}

void RemoveExcessHits(std::vector<TCagraHit>& seghits) {
  std::vector<unsigned int> erase_list;
  for (auto i=0u; i<seghits.size();i++) {
    if (seghits[i].GetCharge() < 0) {
      erase_list.push_back(i);
    }
  }
  for (auto& idx : erase_list) {
    seghits.erase(seghits.begin()+idx);
  }
  if (seghits.size() > 2) {
    seghits.erase(std::min_element(seghits.begin(),seghits.end(),[](TCagraHit& one, TCagraHit& two){ return one.GetCharge() < two.GetCharge(); } ));
  }
}

int TCagra::BuildHits(std::vector<TRawEvent>& raw_data){
  // --- uncomment to simulate array data ------------------------ //
  // static std::mt19937 mt(std::chrono::system_clock::now().time_since_epoch().count());
  // static std::uniform_real_distribution<float> random(0,1);
  // bool first_time = true;
  // --- uncomment to simulate array data ------------------------ //

  std::unordered_map<int, std::vector<TCagraHit> > cc_hits;
  std::unordered_map<int, std::vector<TCagraHit> > seg_hits;

  for (auto& event : raw_data) {
    SetTimestamp(event.GetTimestamp());

    auto buf = event.GetPayloadBuffer();
    TANLEvent anl(buf);

    unsigned int address = ( (1<<24) +
                             (anl.GetBoardID() << 8) +
                             anl.GetChannel() );

    // --- uncomment to simulate array data ------------------------ //
    // static int prev_detnum = -1;
    // static char prev_leaf = -1;
    // if (random(mt) < 0.2 || first_time) {
    //   first_time = false;
    //   // pick new clover and leaf
    //   address = GetRandomCAGRAChannel(-1,-1);
    //   prev_detnum = TChannel::Get(address)->GetArrayPosition();
    //   prev_leaf = *TChannel::Get(address)->GetArraySubposition();

    // } else {
    //   // use previous clover and leaf
    //   address = GetRandomCAGRAChannel(prev_detnum,prev_leaf);
    // }
    // --- uncomment to simulate array data ------------------------ //

    TCagraHit* hit = nullptr;
    TChannel* chan = TChannel::GetChannel(address);
    if (chan) {
      int detnum = chan->GetArrayPosition(); // clover number
      //char leaf = *chan->GetArraySubposition(); // leaf number
      int segnum = chan->GetSegment(); // segment number

      // seperate out central contact hits, from segment hits
      if (segnum == 0) {
        cc_hits[detnum].emplace_back();
        hit = &cc_hits[detnum].back();
      } else {
        seg_hits[detnum].emplace_back();
        hit = &seg_hits[detnum].back();
        hit->MarkAsSegmentHit();
      }

      if (*chan->GetSystem() == 'L') {
        // do trace analysis for LaBr3
        hit->SetTrace(anl.GetTrace());
        hit->SetCharge(hit->GetTraceEnergy(0,57));
      } else {
        hit->SetTrace(anl.GetTrace());
        // set clover charge from pre/post rise charges
        // both segments and central contacts will have the same energy range
        // either both (+) or both (-)
        if (segnum) {
          hit->SetCharge(-1*TANLEvent::GetSignalPolarity()*anl.GetEnergy());
        } else {
          hit->SetCharge(TANLEvent::GetSignalPolarity()*anl.GetEnergy());
        }
      }
    } else {
      // no channel map for address exists
      // cc_hits[999].emplace_back();
      // hit = &cc_hits[999].back();
      continue;
    }

    hit->SetAddress(address);
    hit->SetTimestamp(event.GetTimestamp());
    hit->SetDiscTime(anl.GetCFD());
    hit->SetPrevDiscTime(anl.GetPrevDisc());
    hit->SetPreRise(anl.GetPreE());
    hit->SetPostRise(anl.GetPostE());
    hit->SetFlags(anl.GetFlags());
    hit->SetBaseSample(anl.GetBaseSample());
    hit->SetSampledBaseline(anl.GetBaseline());
    hit->SetPrevPostRiseBeginSample(anl.GetPrevPostBegin());
    hit->SetPreRiseEndSample(anl.GetPreEnd());
    hit->SetPreRiseBeginSample(anl.GetPreBegin());
    hit->SetTimingMarks(anl.GetTimingMarks());

  }

  // // if there is a segment hit without a central contact, ignore it
  // if (cc_hits.size() == 0) {
  //   return 0;
  // }

  // now let's do some work with organizing the segments with the central contacts
  for (auto& det : seg_hits) { // there is only one instance of each detector
    if (cc_hits.count(det.first)==0) { continue; } // skip any segments which dont have a core hit in the same detector

    auto& seghits = det.second;
    auto& cchits = cc_hits.at(det.first);
    auto nSegments = seghits.size();
    auto nCores = cchits.size();

    if (seghits[0].GetSystem() == 'Y') { // Yale clover
      switch(nCores) {
      case 1: {
        // if there is only 1 cc, then add in segments according to their respective cc
        for (auto& seg_hit : seghits) {
          if (seg_hit.GetLeaf() == 'M') {
            cchits[0].AddSegmentHit(seg_hit);
          }
          else if ((cchits[0].GetLeaf() == 'A' || cchits[0].GetLeaf() == 'D') &&
                   seg_hit.GetLeaf() == 'L') {
            cchits[0].AddSegmentHit(seg_hit);
          }
          else if ((cchits[0].GetLeaf() == 'B' || cchits[0].GetLeaf() == 'C') &&
                   seg_hit.GetLeaf() == 'R') {
            cchits[0].AddSegmentHit(seg_hit);
          }
        }

        break;
      }
      case 2:
      case 3:
      case 4: {
        if (nSegments == 1) {
          auto& seg_hit = seghits[0];
          auto seg_leaf = seg_hit.GetLeaf();
          //assert(seg_hit.GetLeaf() == 'M');
          for (auto& cc_hit : cchits) {
            auto cc_leaf = cc_hit.GetLeaf();
            if (seg_leaf == 'M') {
              cc_hit.AddSegmentHit(seg_hit);
            }
            else if ((cc_leaf == 'A' || cc_leaf == 'D') &&
                     seg_leaf == 'L') {
              cc_hit.AddSegmentHit(seg_hit);
            }
            else if ((cc_leaf == 'B' || cc_leaf == 'C') &&
                     seg_leaf == 'R') {
              cc_hit.AddSegmentHit(seg_hit);
            }
          }
        } else if (nSegments == 2) {

          // special case
          if (nCores == 2) {
            // if central contacts are in same column (same theta)
            if(((cchits[0].GetLeaf() == 'A' || cchits[0].GetLeaf() == 'D') &&
                (cchits[1].GetLeaf() == 'A' || cchits[1].GetLeaf() == 'D'))
               ||
               ((cchits[0].GetLeaf() == 'B' || cchits[0].GetLeaf() == 'C') &&
                (cchits[1].GetLeaf() == 'B' || cchits[1].GetLeaf() == 'C')))
            {
              // don't do anything with the segments as it's ambiguous.
              // could possibly consider doing energy matching in the future.
              break;
            }
          }

          // first pass only add L and R segments
          for (auto& seg_hit : seghits) {
            auto seg_leaf = seg_hit.GetLeaf();
            for (auto& cc_hit : cchits) {
              auto cc_leaf = cc_hit.GetLeaf();
              if ((cc_leaf == 'A' || cc_leaf == 'D') &&
                  seg_leaf == 'L') {
                cc_hit.AddSegmentHit(seg_hit);
              }
              else if ((cc_leaf == 'B' || cc_leaf == 'C') &&
                       seg_leaf == 'R') {
                cc_hit.AddSegmentHit(seg_hit);
              }
            }
          }

          // second pass add in M segments
          for (auto& seg_hit : seghits) {
            if (seg_hit.GetLeaf() != 'M') { continue; }
            for (auto& cc_hit : cchits) {
              // only add in M segment if it doesn't have a L or R segment
              if (cc_hit.GetNumSegments() == 0) {
                cc_hit.AddSegmentHit(seg_hit);
              }
            }
          }

        } else if (nSegments == 3) {
          // do nothing with them as the segment assignment is ambiguous
        }
        break;
      }

      }

    } else { // IMP clover
      // loop over segments for the given detector
      for (auto& seg_hit : seghits) {
        // loop over core hits for the given detector
        for (auto& cc_hit : cchits) {
          // if it's the same crystal, add segment hits to cc hit
          if (cc_hit.GetLeaf() == seg_hit.GetLeaf()) {
            cc_hit.AddSegmentHit(seg_hit);
          }
        }

      }
    }
  }

  // add resulting cc hits into cagra_hits
  for (auto& det : cc_hits) {
    for (auto const& core : det.second) {
      cagra_hits.push_back(core);
      fSize++;
    }
  }

  return Size();
}

TVector3 TCagra::GetSegmentPosition(int slot, char core, int seg) {
  if(slot < 1 || slot > 16 || seg < 0 || seg > 4) {
    return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }
  LoadDetectorPositions();

  int index = (slot << 16) + (((int)core) << 8) + seg;

  if (detector_positions.count(index)>0) {
    return detector_positions.at(index);
  } else {
    return TVector3(0.,0.,0.);
  }
}

void TCagra::LoadDetectorPositions() {
  if (positions_loaded) { return; }

  std::string filename = std::string(getenv("GRUTSYS")) + "/config/CAGRA_positions.txt";

  //Read the locations from file.
  std::ifstream infile(filename);

  if(!infile){
    std::cout << "Cagra positions file \"" << filename << "\""
              << " does not exist, skipping" << std::endl;
    return;
  }

  std::string line;
  while (std::getline(infile,line)) {
    //Parse the line
    int nclover, nsegment;
    char ncrystal;
    double rho,theta,phi;
    int extracted = sscanf(line.c_str(),"slot.%02d.%c.%01d.vec_sph: %lf %lf %lf",
                           &nclover,&ncrystal,&nsegment,&rho,&theta,&phi);
    if (extracted!=6) {
      continue;
    }

    int index = (nclover << 16) + (((int)ncrystal) << 8);

    // assign a unique id to each clover crystal
    if (ncrystal != 'X') {
      static size_t core_count = 0;
      if (crystal_ids.count(index) == 0) {
        crystal_ids[index] = core_count++;
      }
    }

    index += nsegment;
    //Pack into the vector of transformations.
    TVector3 vec = TVector3(1.,1.,1.);
    vec.SetTheta(theta*TMath::Pi()/180.);
    vec.SetPhi(phi*TMath::Pi()/180.);
    vec.SetMag(rho);
    detector_positions[index] = vec;
  }

  positions_loaded = true;
}
int TCagra::GetCrystalId(int slot, char core) {
  if (core == 'X') { return -1; }
  int index = (slot << 16) + (((int)core) << 8);
  return crystal_ids.at(index);
}
int TCagra::GetSegmentId(int slot, char core, int segnum, char system) {
  if (segnum == 0 || core == 'X') { return -1; } // bgo or central contact

  int index = 0;
  if (system == 'Y') {
    index = slot*10000 + segnum;
  } else if (system == 'I') {
    int coreval = ((int)core) - 0x41 + 1; // normalize to A=1,B=2,etc.
    index = slot*10000 + coreval*100 + segnum;
  }

  return segment_ids.at(index);
}

void TCagra::Print(Option_t *opt) const { }

void TCagra::Clear(Option_t *opt) {
  cagra_hits.clear();
}
