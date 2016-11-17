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


TCagra::TCagra(){
  Clear();
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
      }

      if (*chan->GetSystem() == 'L') {
        // do trace analysis for LaBr3
        hit->SetTrace(anl.GetTrace());
        hit->SetCharge(hit->GetTraceEnergy(0,57));
      } else {
        // set clover charge from pre/post rise charges
        hit->SetTrace(anl.GetTrace());
        hit->SetCharge(TANLEvent::GetSignalPolarity()*anl.GetEnergy());
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
    if (cc_hits.count(det.first)==0) { continue; } // skip ay segments which dont have a core hit in the same detector

    auto& seghits = det.second;
    auto& cchits = cc_hits.at(det.first);
    auto nSegments = seghits.size();
    auto nCores = cchits.size();
    //std::cout << "Detector: " << det.first << " System: " << seghits[0].GetSystem() << " nCores: " << nCores << " nSegments: " << nSegments << std::endl;

    if (seghits[0].GetSystem() == 'Y') { // Yale clover
      // don't need to check since segments must come with a core hit
      //auto nCores = 0u;
      //if (cc_hits.count(detnum)) {
      //}
      switch(nCores) {
      case 1: {
        // if there is only 1 cc, then add all segment hits to it
        //assert(nSegments <= 2);
        for (auto& seg_hit : seghits) {
          cchits[0].AddSegmentHit(seg_hit);
        }
        break;
      }
      case 2:
      case 3:
      case 4: {
        if (nSegments == 1) {
          auto& seg_hit = seghits[0];
          //assert(seg_hit.GetLeaf() == 'M');
          for (auto& cc_hit : cchits) {
            cc_hit.AddSegmentHit(seg_hit);
          }
        } else if (nSegments == 2) {

          // special case
          if (nCores == 2) {
            // if central contacts are in same column (same theta)
            if(((cchits[0].GetLeaf() == 'A' || cchits[0].GetLeaf() == 'C') &&
                (cchits[1].GetLeaf() == 'A' || cchits[1].GetLeaf() == 'C'))
               ||
               ((cchits[0].GetLeaf() == 'B' || cchits[0].GetLeaf() == 'D') &&
                (cchits[1].GetLeaf() == 'B' || cchits[1].GetLeaf() == 'D')))
            {
              // don't do anything with the segments as it's ambiguous.
              // could possibly consider doing energy matching in the future.
              break;
            }
          }

          // first pass only add L and R segments
          for (auto& seg_hit : seghits) {
            for (auto& cc_hit : cchits) {
              if ((cc_hit.GetLeaf() == 'A' || cc_hit.GetLeaf() == 'C') &&
                  seg_hit.GetLeaf() == 'L') {
                cc_hit.AddSegmentHit(seg_hit);
              }
              else if ((cc_hit.GetLeaf() == 'B' || cc_hit.GetLeaf() == 'D') &&
                       seg_hit.GetLeaf() == 'R') {
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

  //std::cin.get();
  return Size();
}

/*
  int TCagra::BuildHits(std::vector<TRawEvent>& raw_data){
  static std::mt19937 mt(std::chrono::system_clock::now().time_since_epoch().count());
  static std::uniform_real_distribution<float> random(0,1);

  std::vector<unsigned int> yale_segments;
  std::vector<unsigned int> yale_addresses;

  unsigned int count = 0;
  for (auto& event : raw_data) {
    SetTimestamp(event.GetTimestamp());

    auto buf = event.GetPayloadBuffer();
    TANLEvent anl(buf);

    unsigned int address = ( (1<<24) +
                             (anl.GetBoardID() << 8) +
                             anl.GetChannel() );



    static int prev_detnum = -1;
    static char prev_leaf = -1;

    if (random(mt) < 0.2 || count == 0) {
      // pick new clover and leaf
      address = GetRandomCAGRAChannel(-1,-1);
      prev_detnum = TChannel::Get(address)->GetArrayPosition();
      prev_leaf = *TChannel::Get(address)->GetArraySubposition();

    } else {
      // use previous clover and leaf
      address = GetRandomCAGRAChannel(prev_detnum,prev_leaf);



    }


    TChannel* chan = TChannel::GetChannel(address);

    auto boardnum = (address&0x0000ff00)>>8;
    auto channum = (address&0x000000ff)>>0;

    static int lines_displayed = 0;
    if(!chan){
      if(lines_displayed < 10) {
        // std::cout << "Unknown (board id, channel): ("
        //           << anl.GetBoardID() << ", " << anl.GetChannel()
        //           << ")" << std::endl;
        std::cout << "Unknown (board id, channel): ("
                  << boardnum << ", " << channum
                  << ")" << std::endl;
      } else if(lines_displayed==1000){
        std::cout << "I'm going to stop telling you that the channel was unknown,"
                  << " you should probably stop the program." << std::endl;
      }
      lines_displayed++;

      continue;
    }


    int detnum = chan->GetArrayPosition(); // clover number
    char leaf = *chan->GetArraySubposition(); // clover number
    int segnum = chan->GetSegment(); // segment number
    char detector_type = *chan->GetSystem();


    // select the correct clover if it exists
    TCagraHit* hit = NULL;
    for(auto& ihit : cagra_hits) {
      if(ihit.GetDetnum() == detnum && ihit.GetLeaf() == leaf){
        hit = &ihit;
        break;
      }
    }

    // otherwise make a new clover hit
    if(hit == NULL){
      cagra_hits.emplace_back();
      hit = &cagra_hits.back();
      fSize++;
    }

    // if this hit is a Yale L, M, R  segment
    if (detector_type == 'Y' && segnum != 0) {
      //defer until yale loop below
      yale_addresses.push_back(address);
      yale_segments.push_back(count);
      continue;
    }
    // otherwise, add the central contact as normal


    if(segnum==0){
      hit->SetAddress(address);
      hit->SetTimestamp(event.GetTimestamp());
      hit->SetDiscTime(anl.GetCFD());
      hit->SetCharge(anl.GetEnergy());
      hit->SetTrace(anl.GetTrace());
      hit->SetPreRise(anl.GetPreE());
      hit->SetPostRise(anl.GetPostE());
      hit->SetFlags(anl.GetFlags());
      hit->SetBaseSample(anl.GetBaseSample());
      std::cout << "Count: " << count << " Address: " << std::hex << hit->Address() << std::dec << std::endl;
    } else {
      TCagraSegmentHit& seg = hit->MakeSegmentByAddress(address);
      det.SetCharge(anl.GetEnergy());
      det.SetTimestamp(event.GetTimestamp());
      std::cout << "Count: " << count << " Seg. Address: " << std :: hex << seg.Address() << std::dec << std::endl;
      // TODO: the following need implementation
      //seg.SetDiscTime(anl.GetCFD());
      //seg->SetPreRise(anl.GetPreE());
      //seg->SetPostRise(anl.GetPostE());
      //seg->SetFlags(anl.GetFlags());
      //seg->SetBaseSample(anl.GetBaseSample());
      seg.SetTrace(anl.GetTrace());
    }
    count++;
  }

  // loop over the data once more and pick out the yale segments
  // and add them in to each possible central contact
  std::set<unsigned int> problem_leaves;

  auto i = 0u;
  for (auto& event_idx : yale_segments) {
    auto& event = raw_data[event_idx];
    auto buf = event.GetPayloadBuffer();
    TANLEvent anl(buf);



    unsigned int seg_address = ( (1<<24) +
                                 (anl.GetBoardID() << 8) +
                                 anl.GetChannel() );


    seg_address = yale_addresses[i++]; // only for simulation !!!

    auto seg_chan = TChannel::GetChannel(seg_address);

    auto nhit = 0u;
    for (auto& hit : cagra_hits) {
      auto chan = TChannel::GetChannel(hit.Address());
      // skip non yale central contacts
      if (*chan->GetSystem() != 'Y') { continue; }

      int detnum = chan->GetArrayPosition();


      // if detnum of current hit matches detnum of current segment hit
      if (detnum == seg_chan->GetArrayPosition()) {
        char leaf = *chan->GetArraySubposition(); //  A/B/C/D
        char seg_leaf = *seg_chan->GetArraySubposition(); //  L/M/R
        if ( ((seg_leaf == 'L' || seg_leaf == 'M') && (leaf == 'A' || leaf == 'C')) ||
             ((seg_leaf == 'R' || seg_leaf == 'M') && (leaf == 'B' || leaf == 'D')) ) {

          // if this hit already has a segment mark it for later review
          if (hit.GetNumSegments() != 0) {
            problem_leaves.insert(nHit);
          }

          TCagraSegmentHit& seg = hit->MakeSegmentByAddress(address);
          seg.SetCharge(anl.GetEnergy());
          seg.SetTimestamp(event.GetTimestamp());
          //std::cout << "Count: " << count << " Seg. Address: " << std :: hex << seg.Address() << std::dec << std::endl;
          // TODO: the following need implementation
          //seg.SetDiscTime(anl.GetCFD());
          //seg->SetPreRise(anl.GetPreE());
          //seg->SetPostRise(anl.GetPostE());
          //seg->SetFlags(anl.GetFlags());
          //seg->SetBaseSample(anl.GetBaseSample());
          seg.SetTrace(anl.GetTrace());
        }
      }
      nhit++;
    }
  }

  for (auto& idx1 : problem_leaves) {
    auto& hit1 = cagra_hits[idx];
    for (auto idx2 = 0u; idx2 < cagra_hits.size(); idx2++) {
      if (idx1 == idx2) { continue; }
      auto& hit2 = cagra_hits[idx];
      if (hit1.GetDetnum() == hit2.GetDetnum()) {
        auto leaf1 = hit1.GetLeaf();
        auto leaf2 = hit2.GetLeaf();


        if( // if both hits occured on one side of the clover
          ((leaf1 == 'A' || leaf1 == 'C') && (leaf2 == 'A' || leaf2 == 'C')) ||
          ((leaf1 == 'B' || leaf1 == 'D') && (leaf2 == 'B' || leaf2 == 'D'))
          )
        {
          if (hit2.GetNumSegments()>1) {
            // could compare charges of hits and segments to determine match
            // if (hit1.GetEnergy() >= hit2.GetEnergy()) {
            //   double max_energy = -9e99;
            //   int output = 0;
            //   for (auto n=0u; n < hit1.GetNumSegments(); n++) {
            //     auto& segment = hit1.GetSegment(i);
            //     if(segment.GetEnergy() > max_energy){
            //       output = n;
            //       max_energy = segment.GetEnergy();
            //     }
            //   }
            // } else {
            // }

            // for now we erase both segments since ambiguity exists
            hit1.ClearSegments();
            hit2.ClearSegments();

          } else {
            std::cout << "Error in leaf parsing logic for yale clover" << std::endl;
            exit(0);
          }
        } else { // hits happened on opposite sides of detector


        }


      }
    }
  }


  //TCagraHit hit;
  //hit.BuildFrom(buf);
  //hit.SetTimestamp(event.GetTimestamp());
  //InsertHit(hit);
  std::cout << "Done." << std::endl;
  std::cin.get();

  return Size();
}
*/

TVector3 TCagra::GetSegmentPosition(int detnum, char subpos, int segnum) {
  if(detnum < 0 || detnum > 15 || segnum < 0 || segnum > 2 ||
     subpos < 0x41 || subpos > 0x44){
    return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }
  LoadDetectorPositions();

  int index = (detnum << 16) + (((int)subpos) << 8) + segnum;

  if (detector_positions.count(index)>0) {
    return detector_positions.at(index);
  } else {
    return TVector3(0.,0.,0.);
  }

  // double segment_height = 1.0;
  // double perp_distance = 1.5;

  // // Middle of the segment
  // double segment_phi = 3.1415926535/4.0;
  // double segment_z = segment_height/2.0;

  // double crystal_phi = segment_phi + (segnum-2)*3.1415926/2.0;
  // double crystal_z = segment_z + ((segnum-1)/4)*segment_height;

  // TVector3 crystal_pos(1,0,0);
  // crystal_pos.SetZ(crystal_z);
  // crystal_pos.SetPhi(crystal_phi);
  // crystal_pos.SetPerp(perp_distance);

  // TVector3 global_pos = CrystalToGlobal(detnum, crystal_pos);


  // return global_pos;
}

void TCagra::LoadDetectorPositions() {
  static bool loaded = false;
  if(loaded){
    return;
  }
  loaded = true;

  //std::string filename = std::string(getenv("GRUTSYS")) + "/../config/SeGA_rotations.txt";
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
    double x,y,z;
    int extracted = sscanf(line.c_str(),"clover.%02d.%c.%01d.vec: %lf %lf %lf",
                           &nclover,&ncrystal,&nsegment,&x,&y,&z);
    if (extracted!=6) {
      continue;
    }

    int index = (nclover << 16) + (((int)ncrystal) << 8) + nsegment;


    //Pack into the vector of transformations.
    TVector3 vec = TVector3(x,y,z);
    detector_positions[index] = vec;
  }
}

void TCagra::Print(Option_t *opt) const { }

void TCagra::Clear(Option_t *opt) {
  cagra_hits.clear();
}




      // case 2: {
      //   if (nSegments == 1) {
      //     auto& seg_hit = seghits[0];
      //     assert(seg_hit.GetLeaf() == 'M');
      //     for (auto& cc_hit : cchits) {
      //       cc_hit.AddSegmentHit(seg_hit);
      //     }
      //   } else if (nSegments == 2) {

      //     // if central contacts are in same column (same theta)
      //     if(((cchits[0].GetLeaf() == 'A' || cchits[0].GetLeaf() == 'C') &&
      //         (cchits[1].GetLeaf() == 'A' || cchits[1].GetLeaf() == 'C'))
      //        ||
      //        ((cchits[0].GetLeaf() == 'B' || cchits[0].GetLeaf() == 'D') &&
      //         (cchits[1].GetLeaf() == 'B' || cchits[1].GetLeaf() == 'D')))
      //     {
      //       // don't do anything with the segments as it's ambiguous.
      //       // could possibly consider doing energy matching in the future.
      //     } else { // otherwise they are adjacent

      //       // first pass only add L and R segments
      //       for (auto& seg_hit : seghits) {
      //         for (auto& cc_hit : cchits) {
      //           if ((cc_hit.GetLeaf() == 'A' || cc_hit.GetLeaf() == 'C') &&
      //               seg_hit.GetLeaf() == 'L') {
      //             cc_hit.AddSegmentHit(seg_hit);
      //           }
      //           else if ((cc_hit.GetLeaf() == 'B' || cc_hit.GetLeaf() == 'D') &&
      //                    seg_hit.GetLeaf() == 'R') {
      //             cc_hit.AddSegmentHit(seg_hit);
      //           }
      //         }
      //       }

      //       // second pass add in M segments
      //       for (auto& seg_hit : seghits) {
      //         if (seg_hit.GetLeaf() != 'M') { continue; }
      //         for (auto& cc_hit : cchits) {
      //           // only add in M segment if it doesn't have a L or R segment
      //           if (cc_hit.GetNumSegments() == 0) {
      //             cc_hit.AddSegmentHit(seg_hit);
      //           }
      //         }
      //       }

      //     }
      //   } else if (nSegments == 3) {
      //     // do nothing with them as the segment assignment is ambiguous
      //   }
      //   break;
      // }
