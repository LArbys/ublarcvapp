#ifndef __UBLARCVAPP_MCTOOLS_MCPIXELLABELS_H__
#define __UBLARCVAPP_MCTOOLS_MCPIXELLABELS_H__

#include <set>
#include <array>
#include <vector>

namespace ublarcvapp {
namespace mctools {

class MCPixelLabels {

public:

  MCPixelLabels()
  : index(-1),
    imgcoord({-1,-1,-1,-1,-1}),
    pos({0,0,0}),
    pos_reco({0,0,0}),
    edep({0.0,0.0,0.0}),
    pixval({0,0,0})
  {
      trackids.clear();
      aids.clear();
      pids.clear();
      origin.clear();
  };

  ~MCPixelLabels() {};

  long index; ///< book-keeping label

  std::array<int,5>   imgcoord; //< wire plane image coordinates (u,v,y,row,tick)
  std::array<float,3> pos;      //< true energy deposition 
  std::array<float,3> pos_reco; //< position as it appears when reco'd from wire plane images

  std::array<double,3> edep;    //< energy deposited as seen per wire plane
  std::array<float,3> pixval;   //< pixel value at wire plane

  std::set<long> trackids;      //< track IDs contributing to edep position
  std::set<long> aids;          //< ancestor track ID
  std::set<int> pids;           //< particle ID
  std::set<int> origin;         //< origin ID: 1=neutrino, 2=cosmic, 0=other

  int tick() { return imgcoord[4]; };
  int row()  { return imgcoord[3]; };
  int numTrackIDs() { return trackids.size(); };


};


}
}

#endif