#include "Gen3DEndpoints.h"

// larlite
#include "LArUtil/Geometry.h"

// ublarcvapp
#include "ublarcvapp/UBWireTool/UBWireTool.h"

namespace ublarcvapp {
namespace dltagger {

  /**
   * find the 3d end points for cropmaskcombo data
   *
   */
  Gen3DEndpoints::Gen3DEndpoints( const FeaturesMaskCombo& features )
    : larcv::larcv_base("Gen3DEndpoints"),
    pfeatures(&features)
  {
    _gather_plane_endpoints(features);
    _gen_3dendpoints();
  }

  void Gen3DEndpoints::_gather_plane_endpoints( const FeaturesMaskCombo& features ) {
    
    int nplanes = (int)features.pcropdata->crops_v.size();
    
    // operate on planes independently first
    mask_endpoints_vv.clear();
    mask_projdist_vv.clear();
    
    for ( size_t p=0; p<nplanes; p++ ) {

      if ( features.pcropdata->istwoplane() && features.pcropdata->badPlane()==(int)p ) {
        mask_endpoints_vv.push_back( std::vector<EndPoint_t>() );
        mask_projdist_vv.push_back(  std::vector<float>() );
        continue;
      }
      
      // find the extreme points along the 1st pca-axis. they are buried deep into the ContourShapeMeta classes
      const std::vector< ublarcvapp::ContourShapeMeta >& contour_v =
       features.combo_mask_contour.m_plane_atomicmeta_v.at(p);

      const larcv::ImageMeta& meta = features.pcropdata->mask_v.at(p).meta();

      // all of this is in (col,row) pixel coordinates
      std::vector<float> pca1_dir = features.pca1_dir_vv.at(p);
      float pcalen = 0.;
      for (int i=0; i<2; i++ )
        pcalen += pca1_dir[i]*pca1_dir[i];
      pcalen = sqrt(pcalen);
      
      auto const& pca_mean        = features.pca_mean_vv.at(p);

      float pca_proj[2] = {0,0};   // min, max projection
      float endpoints[2][2] = {0}; // min, max (x,y) of maximal points
      
      for ( size_t ictr=0; ictr<contour_v.size(); ictr++ ) {
        const ublarcvapp::ContourShapeMeta& contourdata = contour_v[ictr];
        if ( !contourdata.hasValidPCA() ) continue;
        
        auto const& start = contourdata.getFitSegmentStart();
        auto const& end   = contourdata.getFitSegmentEnd();

        float start_proj = ((start.x-pca_mean[0])*pca1_dir[0] + (start.y-pca_mean[1])*pca1_dir[1])/pcalen;
        float end_proj   = ((end.x-pca_mean[0])*pca1_dir[0] + (end.y-pca_mean[1])*pca1_dir[1])/pcalen;
        
        // min projection
        if ( pca_proj[0] > start_proj ) {
          pca_proj[0] = start_proj;
          endpoints[0][0] = start.x;
          endpoints[0][1] = start.y;
        }
        if ( pca_proj[0] > end_proj ) {
          pca_proj[0] = end_proj;
          endpoints[0][0] = end.x;
          endpoints[0][1] = end.y;
        }

        // max projection
        if ( pca_proj[1] < start_proj ) {
          pca_proj[1] = start_proj;
          endpoints[1][0] = start.x;
          endpoints[1][1] = start.y;
        }
        if ( pca_proj[1] < end_proj ) {
          pca_proj[1] = end_proj;
          endpoints[1][0] = end.x;
          endpoints[1][1] = end.y;
        }
        
      }//end of contour loop

      std::vector<float> minpt(2);
      std::vector<float> maxpt(2);

      minpt[0] = endpoints[0][0];
      minpt[1] = endpoints[0][1];
      maxpt[0] = endpoints[1][0];
      maxpt[1] = endpoints[1][1];

      LARCV_DEBUG() << "plane[" << p << "] pca1_dir: (" << pca1_dir[0]*meta.pixel_width() << ", " << pca1_dir[1]*meta.pixel_height() << ")" << std::endl;
      LARCV_DEBUG() << "plane[" << p << "] min contour pt: (" << meta.pos_x((int)minpt[0]) << "," << meta.pos_y((int)minpt[1]) << ") proj=" << pca_proj[0] << std::endl;
      LARCV_DEBUG() << "plane[" << p << "] max contour pt: (" << meta.pos_x((int)maxpt[0]) << "," << meta.pos_y((int)maxpt[1]) << ") proj=" << pca_proj[1] << std::endl;      

      // these are in local pixel coordinates
      std::vector< std::vector<float> > endpoints_v;
      endpoints_v.push_back( minpt );
      endpoints_v.push_back( maxpt );

      std::vector< float > projdist_v;
      projdist_v.push_back( pca_proj[0] );
      projdist_v.push_back( pca_proj[1] );
      
      mask_endpoints_vv.emplace_back( std::move(endpoints_v) );
      mask_projdist_vv.emplace_back(  std::move(projdist_v) );
    }// end of plane loop
    
  }

  void Gen3DEndpoints::_gen_3dendpoints() {
    // ok, one of the easiest steps to screw up.
    // we have to do things differently for
    // vertically extended tracks and horizontally extended (isochronous) tracks.
    // vertical: we loop through the plane, finding the most extreme point in time.
    // horizontal: we use the  y-plane to define true Z-positions and a Z-length.
    //             we assume its complete.
    //             (lacking examples, we haven't implemented this yet)


    // vertical 3d point making
    int max_row = 0;
    int min_row = 10000;
    int max_pt_plane = 0;
    int min_pt_plane = 0;
    std::vector<float> max_pt;
    std::vector<float> min_pt;

    for ( size_t p=0; p<3; p++ ) {

      // get crop mask
      const larcv::ImageMeta& meta = pfeatures->pcropdata->mask_v.at(p).meta();
      if ( pfeatures->pcropdata->istwoplane() && pfeatures->pcropdata->badPlane()==(int)p ) { 
        // empty image. skip
        continue;
      }
      
      for ( size_t i=0; i<2; i++ ) {
        auto const& endpt = mask_endpoints_vv[p][i];
        float row = endpt[1];
        if ( max_row<row ) {
          max_row = row;
          max_pt_plane = p;
          max_pt = endpt;
        }
        if ( min_row>row ) {
          min_row = row;
          min_pt_plane = p;
          min_pt = endpt;
        }
      }
      
    }//end of plane
    
    LARCV_DEBUG() << "target max point: plane=" << max_pt_plane
                  << " tick=" << pfeatures->pcropdata->mask_v.at(max_pt_plane).meta().pos_y((int)max_row) << std::endl;
    LARCV_DEBUG() << "target min point: plane=" << min_pt_plane 
                  << " tick=" << pfeatures->pcropdata->mask_v.at(min_pt_plane).meta().pos_y((int)min_row) << std::endl;
    
    // now, create point in other planes by extending along pca lines
    std::vector< float > extension_pts[3][2]; // [plane][min or max]
    for ( size_t p=0; p<3; p++ ) {

      // in pixel coordinates
      auto& ext_pt_min = extension_pts[p][0];
      auto& ext_pt_max = extension_pts[p][1];
      ext_pt_min.resize(2,0.0);
      ext_pt_max.resize(2,0.0);
      
      const larcv::ImageMeta& meta = pfeatures->pcropdata->mask_v.at(p).meta();
      if ( pfeatures->pcropdata->istwoplane() && pfeatures->pcropdata->badPlane()==(int)p ) {      
        // empty image. skip.
        continue;
      }
      
      // do extension in (col,row) pixel coordinates
      std::vector<float> pca1_dir = pfeatures->pca1_dir_vv.at(p);
      auto const& pca_mean = pfeatures->pca_mean_vv.at(p);

      // max extension: towards max tick
      if ( max_pt_plane==p ) {
        ext_pt_max = max_pt;
      }
      else {
        float dt = (max_row-pca_mean[1])/pca1_dir[1];
        ext_pt_max[0] = pca_mean[0] + pca1_dir[0]*dt;
        ext_pt_max[1] = max_row;
      }

      if ( min_pt_plane==p ) {
        ext_pt_min = min_pt;
      }
      else {
        float dt = (min_row-pca_mean[1])/pca1_dir[1];
        ext_pt_min[0] = pca_mean[0] + pca1_dir[0]*dt;
        ext_pt_min[1] = min_row;        
      }

      LARCV_DEBUG() << "within image " << meta.dump() << ":" << std::endl;
      LARCV_DEBUG() << " ext max pt: (" << ext_pt_max[0] << "," << ext_pt_max[1] << ")" << std::endl;
      LARCV_DEBUG() << " ext min pt: (" << ext_pt_min[0] << "," << ext_pt_min[1] << ")" << std::endl;      
      
      // convert from (col,row) to (wire,tick)
      ext_pt_min[0] = meta.pos_x( (int)ext_pt_min[0] );
      ext_pt_max[0] = meta.pos_x( (int)ext_pt_max[0] );
      ext_pt_min[1] = meta.pos_y( (int)ext_pt_min[1] );
      ext_pt_max[1] = meta.pos_y( (int)ext_pt_max[1] );

      // bound extensions
      if ( ext_pt_max[0]<0 ) ext_pt_max[0] = 0;
      if ( ext_pt_max[0]>=larutil::Geometry::GetME()->Nwires(p) )
        ext_pt_max[0] = larutil::Geometry::GetME()->Nwires(p)-1;
      
      if ( ext_pt_min[0]<0 ) ext_pt_min[0] = 0;
      if ( ext_pt_min[0]>=larutil::Geometry::GetME()->Nwires(p) )
        ext_pt_min[0] = larutil::Geometry::GetME()->Nwires(p)-1;
      
    }

    // now we get 3d points
    for (int i=0; i<2; i++ ) { // loop over min/max
      
      std::vector< int > wireids(3,-1);
      std::vector<float> intersection_zy(2,0);
      double triarea = 0;
      int crosses = 0;

      int num_good_planes = 0;
      
      for ( size_t p=0; p<3; p++ ) {

        const larcv::ImageMeta& meta = pfeatures->pcropdata->mask_v.at(p).meta();
        if ( pfeatures->pcropdata->istwoplane() && pfeatures->pcropdata->badPlane()==(int)p ) {
          // empty image. skip
          continue;
        }

        //if ( extension_pts[p][i].size()==2 ) {
        wireids[p] = (int)extension_pts[p][i][0];
        num_good_planes++;
        //}
      }

      if ( num_good_planes==3 ) {
        // 3 Plane End point creation
        ublarcvapp::UBWireTool::wireIntersection( wireids, intersection_zy, triarea, crosses );
        std::vector<float> point_tyz(3);
        point_tyz[0] = extension_pts[0][i][1]; // tick coordinate
        point_tyz[1] = intersection_zy[1];     // y-position
        point_tyz[2] = intersection_zy[0];     // z-position
        
        endpt_tyz_v.push_back( point_tyz );
        endpt_wid_v.push_back( wireids );
        endpt_tri_v.push_back( triarea );
        endpt_tpc_v.push_back( crosses );
        
        LARCV_INFO() << "3D Point (3 plane): "
                     << " (" << point_tyz[0] << "," << point_tyz[1] << "," << point_tyz[2] << ")"
                     << " triarea=" << triarea
                     << " crosses=" << crosses
                     << std::endl;
      }
      else if ( num_good_planes==2 ) {
        std::vector<int> goodwid;
        std::vector<int> goodplane;
        int missingplane = -1;
        for ( size_t p=0; p<3; p++ ) {
          if ( wireids[p]!=-1 ) {
            goodwid.push_back( wireids[p] );
            goodplane.push_back(p);
          }
          else
            missingplane = (int)p;
        }
        LARCV_DEBUG() << "we have a missing plane: " << missingplane << std::endl;
        LARCV_DEBUG() << " wids={" << wireids[0] << "," << wireids[1] << "," << wireids[2] << "}" << std::endl;
        LARCV_DEBUG() << " goodplanes={" << goodplane[0] << "," << goodplane[1] << "}" << std::endl;
        LARCV_DEBUG() << " wires on good planes={" << goodwid[0] << "," << goodwid[1] << "}" << std::endl;
        
        // get intersection and infer missing plane wire
        triarea = 0.;
        int otherplane = missingplane;
        int otherwire  = 0;
        //ublarcvapp::UBWireTool::wireIntersection( goodplane[0], goodwid[0], goodplane[1], goodwid[1], intersection_zy, crosses );
        ublarcvapp::UBWireTool::getMissingWireAndPlane( goodplane[0], goodwid[0], goodplane[1], goodwid[1],
                                                        otherplane, otherwire, intersection_zy, crosses );
        std::vector<float> point_tyz(3);
        point_tyz[0] = extension_pts[goodplane[0]][i][1]; // tick coordinate
        point_tyz[1] = intersection_zy[1];     // y-position
        point_tyz[2] = intersection_zy[0];     // z-position

        if ( otherwire!=-1 )
          wireids[missingplane] = otherwire;

        endpt_tyz_v.push_back( point_tyz );
        endpt_wid_v.push_back( wireids );
        endpt_tri_v.push_back( triarea );
        endpt_tpc_v.push_back( crosses );
                
        LARCV_INFO() << "3D Point (2 plane): "
                     << " (" << point_tyz[0] << "," << point_tyz[1] << "," << point_tyz[2] << ")"
                     << " triarea=" << triarea
                     << " crosses=" << crosses
                     << " otherwire=" << otherwire
                     << std::endl;        
      }
      else {
        std::stringstream errmsg;
        errmsg << __FILE__ << ":" << __LINE__ << ": not enough good planes (<2) to form 3D point" << std::endl;
        throw std::runtime_error(errmsg.str());
      }
    }
    
  }

  
}
}
