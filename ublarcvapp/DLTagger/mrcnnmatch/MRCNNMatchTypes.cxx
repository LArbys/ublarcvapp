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
    wire_min = larutil::Geometry::GetME()->Nwires(plane);
    wire_max = 0;
    
    int xoffset = mask.box.min_x();
    int yoffset = mask.box.min_y();
    for ( size_t ipt=0; ipt<mask.points_v.size(); ipt++ ) {
      int row = yoffset + mask.points_v[ipt].y;
      int col = xoffset + mask.points_v[ipt].x;
      float tick = mask.meta.pos_y(row);
      int wire = mask.meta.pos_x(col);

      double xyzstart[3];
      double xyzend[3];
      larutil::Geometry::GetME()->WireEndPoints( plane, wire, xyzstart, xyzend );

      float zmin = (xyzstart[2]<xyzend[2]) ? xyzstart[2] : xyzend[2];
      float zmax = (xyzstart[2]>xyzend[2]) ? xyzstart[2] : xyzend[2];
      
      if ( tick_min>tick ) tick_min = tick;
      if ( tick_max<tick ) tick_max = tick;
      if ( wire_min>wire ) wire_min = wire;
      if ( wire_max<wire ) wire_max = wire;
      if ( detz_min>zmin ) detz_min = zmin;
      if ( detz_max<zmax ) detz_max = zmax;
    }

  }

  std::ostream& operator<<(std::ostream& os,const MaskMatchData& m) {
    return os << "MaskMatchData[plane " << m.plane << ", index " << m.index << "] "
              << " tick=(" << m.tick_min << "," << m.tick_max << ")"
              << " wire=(" << m.wire_min << "," << m.wire_max << ")"
              << " detz=(" << m.detz_min << "," << m.detz_max << ")";
  };
  
  
}
}
    
