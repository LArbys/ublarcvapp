#include "MRCNNMatchTypes.h"

// larlite
#include "LArUtil/Geometry.h"

namespace ublarcvapp {
namespace dltagger {

  MaskMatchData::MaskMatchData( int _plane, int _index, const larcv::ClusterMask& mask )
    : plane(_plane),
      index(_index)
  {

    // we loop over the points in the mask

    tick_min = mask.meta.max_y();
    tick_max = mask.meta.min_y();
    detz_min = larutil::Geometry::GetME()->DetLength();
    detz_max = 0;
    
    
    std::cout << "meta: " << mask.meta.dump() << std::endl;
    int xoffset = mask.box.min_x();
    int yoffset = mask.box.min_y();
    for ( size_t ipt=0; ipt<mask.points_v.size(); ipt++ ) {
      int row = yoffset + mask.points_v[ipt].y;
      int col = yoffset + mask.points_v[ipt].x;
      float tick = mask.meta.pos_y(row);
      int wire = mask.meta.pos_x(col);

      double xyzstart[3];
      double xyzend[3];
      larutil::Geometry::GetME()->WireEndPoints( plane, wire, xyzstart, xyzend );

      float zmin = (xyzstart[2]<xyzend[2]) ? xyzstart[2] : xyzend[2];
      float zmax = (xyzstart[2]>xyzend[2]) ? xyzstart[2] : xyzend[2];
      
      if ( tick_min>tick ) tick_min = tick;
      if ( tick_max<tick ) tick_max = tick;
      if ( detz_min>zmin ) detz_min = zmin;
      if ( detz_max<zmax ) detz_max = zmax;
    }

    std::cout << "Mask bounds: "
              << " DTICK=[" << tick_min << "," << tick_max << "] "
              << " DZ=[" << detz_min << "," << detz_max << "]"
              << std::endl;

  }
  
}
}
    
