#include "MRCNNMatchTypes.h"

// larlite
#include "LArUtil/Geometry.h"

namespace ublarcvapp {
namespace dltagger {

  MaskMatchData::MaskMatchData( int _plane, int _index, const larcv::ClusterMask& mask )
    : plane(_plane),
      index(_index),
      used(false)
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
    return os << "MaskMatchData[plane " << m.plane << ", clustermask-index " << m.index << "] "
              << " tick=(" << m.tick_min << "," << m.tick_max << ")"
              << " wire=(" << m.wire_min << "," << m.wire_max << ")"
              << " detz=(" << m.detz_min << "," << m.detz_max << ")"
              << " used="  << m.used;
  };
  

  // =======================================================================================
  // MaskCombo

  void MaskCombo::addMask( const larcv::ClusterMask& mask, const MaskMatchData& data ) {
    int plane = mask.meta.plane();
    if ( indices.at(plane)!=-1 ) {
      throw std::runtime_error("MaskCombo::addMask: combo already has index for this plane");
    }

    // set index
    indices[plane] = data.index;
    pmasks[plane]  = &mask;
    pdata[plane]   = &data;

    // update union
    if ( union_tick[0]==-1 || union_tick[0]>data.tick_min ) union_tick[0] = data.tick_min;
    if ( union_tick[1]==-1 || union_tick[1]<data.tick_max ) union_tick[1] = data.tick_max;
    if ( union_detz[0]==-1 || union_detz[0]>data.detz_min ) union_detz[0] = data.detz_min;
    if ( union_detz[1]==-1 || union_detz[1]<data.detz_max ) union_detz[1] = data.detz_max;
    
    // update intersection
    // only if overlapping
    bool isoverlapping = iscompatible(data);
    for ( size_t p=0; p<3; p++ ) {
      if ( indices[p]==-1 ) continue;
      if ( (intersection_tick[0]<data.tick_max || intersection_tick[1]>data.tick_min )
           && (intersection_detz[0]<data.detz_max || intersection_detz[1]>data.detz_min ) ) {
        isoverlapping = true;
      }
    }
    
    if ( intersection_tick[0]==-1 || (isoverlapping && intersection_tick[0]<data.tick_min ) ) intersection_tick[0] = data.tick_min;
    if ( intersection_tick[1]==-1 || (isoverlapping && intersection_tick[1]>data.tick_max ) ) intersection_tick[1] = data.tick_max;
    if ( intersection_detz[0]==-1 || (isoverlapping && intersection_detz[0]<data.detz_min ) ) intersection_detz[0] = data.detz_min;
    if ( intersection_detz[1]==-1 || (isoverlapping && intersection_detz[1]>data.detz_max ) ) intersection_detz[1] = data.detz_max;

    IOU = calc_iou(false);
  }


  float MaskCombo::calc_iou( bool wdetz_correction ) const {

    float intersection = (intersection_detz[1]-intersection_detz[0])*(intersection_tick[1]-intersection_tick[0]);
    float unionarea    = (union_tick[1]-union_tick[0])*( pdata[2]->detz_max - pdata[2]->detz_min );
    // will have to deal with no detz from plane2 later, but for now, Y-plane is our anchor (as usual)

    // float union_correction = 0.;
    
    // if ( wdetz_correction ) {
    //   // we account for the wire projection angles
    //   // image dz overlap for one point in the detector. wires look like this: \|/ 
    //   // intersection of Z is 1, but union includes the Z-projection of the U,V wires
    //   // we subtract this off. we use geo info to get the z-projection of the wire

    //   // we find the projections
    //   float zproj[3][4] = { 0 };
    //   for ( size_t p=0; p<indices.size(); p++ ) {
    //     if ( p==0 ) std::cout << std::endl;
    //     if ( indices[p]==-1 ) continue;
    //     if ( p==2 ) continue;

    //     // we find the nearest wire to the min and max z-det range

    //     double testpts[4][3] = { { 0,  117.0, intersection_detz[0]  }, // minz top
    //                              { 0, -117.0, intersection_detz[0]  }, // minz bottom
    //                              { 0,  117.0, intersection_detz[1]  }, // maxz top
    //                              { 0, -117.0, intersection_detz[1]  }  // maxz bot
    //     };    

    //     float near_wireco[4];
    //     int   near_wire[4];
    //     for ( int i=0; i<4; i++ ) {

    //       if ( testpts[i][2]>=1036 )
    //         testpts[i][2] = 1036.0;
    //       if ( testpts[i][2]<0.5 )
    //         testpts[i][2] = 0.5;
          
    //       near_wireco[i] = larutil::Geometry::GetME()->WireCoordinate( testpts[i], p );
    //       near_wire[i] = (int)(near_wireco[i]+0.5);
    //       if ( near_wire[i]<0 ) near_wire[i] = 0;
    //       if ( near_wire[i]>=(int)larutil::Geometry::GetME()->Nwires(p) )
    //         near_wire[i] = (int)larutil::Geometry::GetME()->Nwires(p) - 1;
          
    //       double start[3];
    //       double end[3];        
    //       larutil::Geometry::GetME()->WireEndPoints( p, near_wire[i], start, end );
    //       zproj[p][i] = fabs(end[2]-start[2]);
    //       std::cout << "zproj[plane=" << p << "][pt=" << i << "] = " << zproj[p][i]
    //                 << "  near_wire=" << near_wireco[i]
    //                 << "  start=(" << start[1] << "," << start[2] << ")"
    //                 << "  end=(" << end[1] << "," << end[2] << ")"
    //                 << std::endl;
    //     }
    //   }//end of loop over planes
    //   // min-side correction
    //   float zmin_correction = 0.;
    //   float zmax_correction = 0.;
    //   if ( indices[0]!=-1 ) {
    //     zmin_correction = ( zproj[0][0] > zproj[1][1] ) ? zproj[0][0] : zproj[1][1]; // u-plane uses min-top, v-plane uses min-bottom
    //     zmin_correction *= (union_tick[1]-union_tick[0]);
    //   }
    //   if ( indices[1]!=-1 ) {
    //     zmax_correction = ( zproj[0][1] > zproj[1][0] ) ? zproj[0][1] : zproj[1][0]; // u-plane uses max-bot, v-plane uses max-top
    //     zmax_correction *= (union_tick[1]-union_tick[0]);
    //   }
    //   union_correction = (zmin_correction + zmax_correction);
    //   std::cout << " correction=" << union_correction << " = zmin(" << zmin_correction << ") + zmax(" << zmax_correction << ")" << std::endl;
      
    //   unionarea -= union_correction;
    // }
    
    return intersection/unionarea;
    
  }

  bool MaskCombo::iscompatible( const MaskMatchData& data ) const {
    bool isoverlapping = false;
    for ( size_t p=0; p<3; p++ ) {
      if ( indices[p]==-1 ) continue;
      if ( (intersection_tick[0]<data.tick_max && intersection_tick[1]>data.tick_min )
           && (intersection_detz[0]<data.detz_max && intersection_detz[1]>data.detz_min ) ) {
        isoverlapping = true;
      }
    }
    return isoverlapping;
  }

  std::ostream& operator<<(std::ostream& os,const MaskCombo& m) {
    return os << "MaskCombo("
              << " clustmask-idx=[" << m.indices[0] << "," << m.indices[1] << "," << m.indices[2] << "]"
              << " maskdata-idx=[" << m.maskdata_indices[0] << "," << m.maskdata_indices[1] << "," << m.maskdata_indices[2] << "]"      
              << " intersection{ "
              << " tick=[" << m.intersection_tick[0] << "," << m.intersection_tick[1] << "]"
              << " detz=[" << m.intersection_detz[0] << "," << m.intersection_detz[1] << "]"
              << " }"
              << " union{ "
              << " tick=[" << m.union_tick[0] << "," << m.union_tick[1] << "]"
              << " detz=[" << m.union_detz[0] << "," << m.union_detz[1] << "]"
              << " detz[2]=[" << m.pdata[2]->detz_min << "," << m.pdata[2]->detz_max << "]"
              << " }"
              << " iou: " << m.iou();
  }
  
}
}
    
