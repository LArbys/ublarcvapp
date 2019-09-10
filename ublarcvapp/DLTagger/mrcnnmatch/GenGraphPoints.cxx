#include "GenGraphPoints.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

#include "ublarcvapp/UBWireTool/UBWireTool.h"

#include <cilantro/kd_tree.hpp>

#include <set>
#include <algorithm>

namespace ublarcvapp {
namespace dltagger {

  GenGraphPoints::GenGraphPoints( const FeaturesMaskCombo& featuredata, larcv::msg::Level_t msglevel )
    : larcv::larcv_base("GenGraphPoints"),
    pfeatures(&featuredata)
  {
    
    set_verbosity(msglevel);
    
    // upstream data
    auto const& mask_v = pfeatures->pcropdata->mask_v;
    auto const& crop_v = pfeatures->pcropdata->crops_v;
    auto const& miss_v = pfeatures->pcropdata->missing_v;    
    
    // get the bounds
    LARCV_DEBUG() << "[ define bounds ]" << std::endl;    
    std::vector<MaskExtrema_t> bounds_v;
    DefineBoundPoints( *pfeatures, bounds_v);


    LARCV_DEBUG() << "[Gen (wire,tick) plane pairs for min-x end point]" << std::endl;

    // test the bounds from the y-plane, min-x, max-x
    std::vector< std::vector<float> > min_wiretickpos_p0_v; // store 3d points from matching to U-plane
    std::vector< std::vector<float> > min_wiretickpos_p1_v; // store 3d points from matching to V-plane
    // min wire coordinate bounds: U-plane
    makePointsFixedZ( *pfeatures,
                      crop_v[2].meta().pos_x( bounds_v[2].points[0][0] )*0.3,  0,
                      crop_v[2].meta().pos_y( bounds_v[2].points[0][1] ),
                      min_wiretickpos_p0_v );
    // min wire coordinate bounds: V-plane    
    makePointsFixedZ( *pfeatures,
                      crop_v[2].meta().pos_x( bounds_v[2].points[0][0] )*0.3,  1,
                      crop_v[2].meta().pos_y( bounds_v[2].points[0][1] ),
                      min_wiretickpos_p1_v );

    // make 3D points
    LARCV_DEBUG() << "make 3d point for min-x extremum" << std::endl;
    //std::vector< std::vector<float> > points3d_v;
    //std::vector< std::vector<float> > twid_v;
    make3Dpoints( crop_v[2].meta().pos_x( bounds_v[2].points[0][0] ),
                  crop_v[2].meta().pos_y( bounds_v[2].points[0][1] ),
                  min_wiretickpos_p0_v, min_wiretickpos_p1_v, m_points3d_v, m_twid_v,
                  true );


    LARCV_DEBUG() << "[Gen (wire,tick) plane pairs for max-x end point]" << std::endl;    
    std::vector< std::vector<float> > max_wiretickpos_p0_v; // store 3d points from matching to U-plane
    std::vector< std::vector<float> > max_wiretickpos_p1_v; // store 3d points from matching to V-plane    
    
    // max wire coodinate bounds: U-plane
    makePointsFixedZ( *pfeatures,
                      crop_v[2].meta().pos_x( bounds_v[2].points[1][0] )*0.3,  0,
                      crop_v[2].meta().pos_y( bounds_v[2].points[1][1] ),
                      max_wiretickpos_p0_v );
    // max wire coordiante boudns: V-plane
    makePointsFixedZ( *pfeatures,
                      crop_v[2].meta().pos_x( bounds_v[2].points[1][0] )*0.3,  1,
                      crop_v[2].meta().pos_y( bounds_v[2].points[1][1] ),
                      max_wiretickpos_p1_v );

    LARCV_DEBUG() << "make 3d point for max-x extremum" << std::endl;
    std::vector< std::vector<float> > maxx_pts_v;
    std::vector< std::vector<float> > maxx_twid_v;
    make3Dpoints( crop_v[2].meta().pos_x( bounds_v[2].points[1][0] ),
                  crop_v[2].meta().pos_y( bounds_v[2].points[1][1] ),
                  max_wiretickpos_p0_v, max_wiretickpos_p1_v, maxx_pts_v, maxx_twid_v,
                  true );

    // scan across x
    std::vector< std::vector<float> > ywire_points_v;
    for ( int icol=bounds_v[2].points[0][0]+1; icol<(int)bounds_v[2].points[1][0]-1; icol+=2 ) {
      std::vector< std::vector<float> > yscan_points_v;
      scanTickDim( *pfeatures, crop_v[2].meta().pos_x(icol), 2, yscan_points_v );
      if ( yscan_points_v.size()>0 ) {
        ywire_points_v.push_back( yscan_points_v.front() );

        std::vector< std::vector<float> > wtpos_p0_v;
        std::vector< std::vector<float> > wtpos_p1_v;

        makePointsFixedZ( *pfeatures,
                          yscan_points_v.front()[0]*0.3,  0,
                          yscan_points_v.front()[1],
                          wtpos_p0_v );
        
        makePointsFixedZ( *pfeatures,
                          yscan_points_v.front()[0]*0.3, 1,
                          yscan_points_v.front()[1],
                          wtpos_p1_v );

        make3Dpoints( yscan_points_v.front()[0],
                      yscan_points_v.front()[1],
                      wtpos_p0_v, wtpos_p1_v, m_points3d_v, m_twid_v,
                      false );
      }//if yscan
    }// column loop
    LARCV_DEBUG() << "number of 3d points defined: " << m_points3d_v.size() << " twid_v=" << m_twid_v.size() << std::endl;

    // add end points
    for (size_t i=0; i<maxx_pts_v.size(); i++) {
      m_points3d_v.push_back( maxx_pts_v[i] );
    }
      
    std::sort( m_points3d_v.begin(), m_points3d_v.end() );
    if ( logger().debug() ) {
      std::stringstream ptlist;
      ptlist << "{ ";
      for ( auto& pt : m_points3d_v ) {
        ptlist << "(" << pt[0] << "," << pt[1] << "," << pt[2] << ") ";
      }
      ptlist << "}";
      LARCV_DEBUG() << ptlist.str() << std::endl;
    }
    
    // make graph components
    std::vector<Eigen::Vector3f > graph_nodes;
    std::map< std::pair<int,int>, float > distmap;
    std::map< std::pair<int,int>, float > pixgapmap;
    makeGraph( m_points3d_v, graph_nodes, distmap, pixgapmap );
  }

  /**
   * find the min and max points in the wire and tick directions on each plane
   *
   * @param[in]    features     data product containing contour information of match-image crops
   * @param[inout] maskbounds_v for each plane, struct that saves the wire and tick bounds
   *
   */
  void GenGraphPoints::DefineBoundPoints( const FeaturesMaskCombo& features,
                                          std::vector< MaskExtrema_t >& maskbounds_v ) {

    // for each plane, for each contour,
    //   loop to find wire and tick bounds

    maskbounds_v.resize(3);

    auto const& mask_v = features.pcropdata->mask_v;    // crop images where charge above threshold and inside mask-rcnn masks    
    
    for ( size_t p=0; p<3; p++ ) {
      // initialize
      auto& maskbounds = maskbounds_v[p];
      const larcv::ImageMeta* meta = &mask_v[p].meta();
      
      //auto const& ctr_v = features.combo_charge_contour.m_plane_atomicmeta_v[p];
      auto const& ctr_v = features.combo_mask_contour.m_plane_atomicmeta_v[p];      
      for ( auto const& ctr : ctr_v ) {

        bool checkctr = false;
        
        // min-x boundary
        if ( maskbounds.bounds[0] > ctr.getBBox().tl().x
             || maskbounds.bounds[1] < ctr.getBBox().br().x
             || maskbounds.bounds[2] > ctr.getBBox().br().y
             || maskbounds.bounds[3] < ctr.getBBox().tl().y ) {
          checkctr = true;
        }
         
        if ( checkctr ) {

          for ( auto const& pt : ctr ) {
            if ( maskbounds.bounds[0] > pt.x ) {
              maskbounds.bounds[0] = pt.x;
              maskbounds.points[0][0] = pt.x;
              maskbounds.points[0][1] = pt.y;
            }
            if ( maskbounds.bounds[1] < pt.x ) {
              maskbounds.bounds[1] = pt.x;
              maskbounds.points[1][0] = pt.x;
              maskbounds.points[1][1] = pt.y;
            }
            if ( maskbounds.bounds[2] > pt.y ) {
              maskbounds.bounds[2] = pt.y;
              maskbounds.points[2][0] = pt.x;
              maskbounds.points[2][1] = pt.y;
            }
            if ( maskbounds.bounds[3] < pt.y ) {
              maskbounds.bounds[3] = pt.y;
              maskbounds.points[3][0] = pt.x;
              maskbounds.points[3][1] = pt.y;
            }
          }
          
        }
        
      }
      if ( logger().debug() ) {
        LARCV_DEBUG() << "bounds on plane[" << p << "]" << std::endl;
        LARCV_DEBUG() << "min-wire=" << maskbounds.bounds[0]
                      << ": pix=(" << maskbounds.points[0][0] << "," << maskbounds.points[0][1] << ")" 
                      << ": img=(" << meta->pos_x(maskbounds.points[0][0]) << "," << meta->pos_y(maskbounds.points[0][1]) << ")"
                      << std::endl;
        LARCV_DEBUG() << "max-wire=" << maskbounds.bounds[1]
                      << ": pix=(" << maskbounds.points[1][0] << "," << maskbounds.points[1][1] << ")" 
                      << ": img=(" << meta->pos_x(maskbounds.points[1][0]) << "," << meta->pos_y(maskbounds.points[1][1]) << ")"
                      << std::endl;
        LARCV_DEBUG() << "min-tick=" << maskbounds.bounds[2]
                      << ": pix=(" << maskbounds.points[2][0] << "," << maskbounds.points[2][1] << ")" 
                      << ": img=(" << meta->pos_x(maskbounds.points[2][0]) << "," << meta->pos_y(maskbounds.points[2][1]) << ")"
                      << std::endl;
        LARCV_DEBUG() << "max-tick=" << maskbounds.bounds[3]
                      << ": pix=(" << maskbounds.points[3][0] << "," << maskbounds.points[3][1] << ")" 
                      << ": img=(" << meta->pos_x(maskbounds.points[3][0]) << "," << meta->pos_y(maskbounds.points[3][1]) << ")"
                      << std::endl;
      }
    }
    
  }

  /** 
   *
   * Scan along fixed detector-z coordinate and find charge in mask cluster
   *
   * @param[in] detz detector-z coordinate
   */
  void GenGraphPoints::makePointsFixedZ( const FeaturesMaskCombo& features,
                                         const float detz, const int plane, const float tick,
                                         std::vector< std::vector<float> >& wiretickpos_v ) {
    
    // plane=2: scan tick
    
    float ystart = 117.0;
    float ystep  = -fabs(0.3*cos(30.0*3.14159/180.0));

    // get inputs from upstream data products
    auto const& mask_v = features.pcropdata->mask_v;    // crop images where charge above threshold and inside mask-rcnn masks
    auto const& crop_v = features.pcropdata->crops_v;   // crop images where charge above threshold
    auto const& miss_v = features.pcropdata->missing_v; // crop images for plane not matched in 2-plane matches

    // 3 plane matches for now
    if ( mask_v[plane].meta().rows()==0 || mask_v[plane].meta().cols()==0 )
      return;
    const larcv::Image2D& maskimg = mask_v[plane]; // get plane we are scanning


    float minx = maskimg.meta().min_x(); // get bounds of image
    float maxx = maskimg.meta().max_x(); // get bounds of image
    float ypos = ystart; // position variable we will increment
    double coord[3] = { 0, 0, detz }; // 3d pos we will incremement
    float wireco = larutil::Geometry::GetME()->WireCoordinate( coord, plane ); // query initial position
    if ( wireco<maskimg.meta().min_x() || wireco>=maskimg.meta().max_x() ) {
      // if within crop test for charge
      return;
    }
    // theck tick bounds
    if ( tick<maskimg.meta().min_y() || tick>=maskimg.meta().max_y() )
      return;

    // initial pixel coordinate
    int col = maskimg.meta().col(wireco);
    int row = maskimg.meta().row(tick);

    // scan along ypos @ tick for charge.
    // to avoid duplicate points in same neighborhood, we find the maximum of a region with charge
    bool inhit = false; // flag indicating we are within a charge region
    int ninhit = 0;     // size of charge region
    float maxpixinhit = 0; // current maximum pix value in charge region
    float pixval = 0.;  // current pix value in charge region
    int maxrow = 0;     // pixel coordinate of max in charge region
    int maxcol = 0;     // pixel coordinate of max in charge region

    std::vector<int> cols_in_hit;  // list of max cols

    // start step loop, end when at bottom of TPC
    while ( ypos > -117.0 ) {

      // update position
      coord[1] = ypos;

      // get wire on plane at this position
      wireco = larutil::Geometry::GetME()->WireCoordinate( coord, plane );

      if (wireco<minx || wireco>=maxx ) {
        // keep going if outside the crop
        ypos += ystep;
        continue;
      }

      // get pixel coordinate of wire
      col = maskimg.meta().col(wireco);

      // get pixel value
      pixval = maskimg.pixel(row,col);

      if ( pixval>10.0 ) {
        // if pixel above threshold, we start a hit region
        if ( maxpixinhit<pixval ) {
          // update the maximum charge pixel
          maxpixinhit = pixval;
          maxrow = row;
          maxcol = col;
        }
        ninhit++;
        inhit = true;
      }
      else {
        // out of hit
        if ( inhit ) {
          // save col
          cols_in_hit.push_back( maxcol );
        }
        // reset hit vars
        inhit = false;
        maxpixinhit = 0.;
        maxrow = 0;
        maxcol = 0;
        ninhit = 0;
      }//end of not in hit
      ypos += ystep;
    }

    // store the (wire,tick) positions of the hit regions
    std::stringstream ss_wires;
    for ( auto const& hitcol : cols_in_hit ) {
      std::vector<float> hit = { (float)maskimg.meta().pos_x(hitcol), tick };
      ss_wires << " " << hit[0];
      wiretickpos_v.push_back( hit );
    }
    LARCV_DEBUG() << "number of 'hits' in scan of plane=" << plane << ": " << cols_in_hit.size() << " wires={" << ss_wires.str() << "}" << std::endl;
    
  }

  /** 
   *
   * Scan along tick direction for some plane
   *
   * @param[in] detz detector-z coordinate
   */
  void GenGraphPoints::scanTickDim( const FeaturesMaskCombo& features,
                                    const float wireco, const int plane, 
                                    std::vector< std::vector<float> >& wiretickpos_v ) {
    
    // get inputs from upstream data products
    auto const& mask_v = features.pcropdata->mask_v;    // crop images where charge above threshold and inside mask-rcnn masks
    auto const& crop_v = features.pcropdata->crops_v;   // crop images where charge above threshold
    auto const& miss_v = features.pcropdata->missing_v; // crop images for plane not matched in 2-plane matches

    // make sure its not empty
    if ( mask_v[plane].meta().rows()==0 || mask_v[plane].meta().cols()==0 )
      return;
    
    const larcv::Image2D& maskimg = mask_v[plane]; // get plane we are scanning

    float minx = maskimg.meta().min_x(); // get bounds of image
    float maxx = maskimg.meta().max_x(); // get bounds of image
    if ( wireco<minx || wireco>=maxx ) {
      // if within crop test for charge
      return;
    }
    
    // initial pixel coordinate
    int col = maskimg.meta().col(wireco);

    // scan along ypos @ tick for charge.
    // to avoid duplicate points in same neighborhood, we find the maximum of a region with charge
    bool inhit = false; // flag indicating we are within a charge region
    int ninhit = 0;     // size of charge region
    float maxpixinhit = 0; // current maximum pix value in charge region
    float pixval = 0.;  // current pix value in charge region
    int maxrow = 0;     // pixel coordinate of max in charge region
    int maxcol = 0;     // pixel coordinate of max in charge region

    std::set<int> rows_in_hit;  // list of max cols

    // start step loop, end when at bottom of TPC
    for ( int row=0; row<(int)maskimg.meta().rows(); row++ ) {

      // get pixel value
      pixval = maskimg.pixel(row,col);

      if ( pixval>10.0 ) {
        // if pixel above threshold, we start a hit region
        if ( maxpixinhit<pixval ) {
          // update the maximum charge pixel
          maxpixinhit = pixval;
          maxrow = row;
          maxcol = col;
        }
        ninhit++;
        inhit = true;
      }
      else {
        // out of hit
        if ( inhit ) {
          // save col
          rows_in_hit.insert( maxrow );
        }
        // reset hit vars
        inhit = false;
        maxpixinhit = 0.;
        maxrow = 0;
        maxcol = 0;
        ninhit = 0;
      }//end of not in hit
      
    }//end of row

    // store the (wire,tick) positions of the hit regions
    std::stringstream ss_wires;
    for ( auto const& hitrow : rows_in_hit ) {
      std::vector<float> hit = { (float)maskimg.meta().pos_x(col), (float)maskimg.meta().pos_y(hitrow) };
      ss_wires << " " << hit[0];
      wiretickpos_v.push_back( hit );
    }
    LARCV_DEBUG() << "number of 'hits' in tick scan of plane=" << plane << ": " << wiretickpos_v.size() << " wires={" << ss_wires.str() << "}" << std::endl;
    
  }
  
  void GenGraphPoints::make3Dpoints( const float ywire,
                                     const float tick,
                                     const std::vector< std::vector<float> >& wiretick_p0_v,
                                     const std::vector< std::vector<float> >& wiretick_p1_v,
                                     std::vector< std::vector<float> >& points3d_v,
                                     std::vector< std::vector<float> >& twid_v,
                                     bool allow_2plane_matches ) {

    if ( wiretick_p0_v.size()==0 && wiretick_p1_v.size()==0 )
      return; // not handling this case yet


    std::set< std::vector<int> > pairs_formed;
    
    // u-wire
    // -------
    // we make at most 2 points, one from the top and one from the bottom
    std::vector<int> uwire_tests;
    if ( wiretick_p0_v.size()>0 ) {
      uwire_tests.push_back( (int)wiretick_p0_v.front()[0] );
      if ( wiretick_p0_v.size()>1 )
        uwire_tests.push_back( (int)wiretick_p0_v.back()[0] );
    }

    std::vector< std::vector<float> > uwire_intersections;
    std::vector< float > uwire_triarea;
    std::vector< int >   uwire_vindex;
    std::vector< std::vector<int> > uwire_wid_v;
    int uwireindex = 0;

    LARCV_DEBUG() << "forming 3d points from both directions of u-plane scan: num uwires=" << uwire_tests.size() << std::endl;
    for ( auto const& uwire : uwire_tests ) {
      int vwire_tests = 0;    
      if ( wiretick_p1_v.size()==0 && uwireindex==0 && allow_2plane_matches ) {
        int crosses = 0;
        std::vector<float> intersection;
        int missingplane = 1;
        int missingwire = 0;
        //ublarcvapp::UBWireTool::wireIntersection( 2, (int)ywire, 0, (int)uwire, intersection, crosses );
        ublarcvapp::UBWireTool::getMissingWireAndPlane( 2, (int)ywire, 0, (int)uwire, missingplane, missingwire, intersection, crosses );
        if ( crosses==1 ) {
          uwire_intersections.push_back( intersection );
          uwire_triarea.push_back(0.0);
          uwire_vindex.push_back(0);
          uwire_wid_v.push_back( std::vector<int>{ (int)uwire, (int)missingwire, (int)ywire} );
        }
        LARCV_DEBUG() << "  2-plane intersection: (z,y)=(" << intersection[0] << "," << intersection[1] << ") crosses=" << crosses << std::endl;
      }
      else if ( wiretick_p1_v.size()>0 ) {
        int vwire_index = (uwireindex==0) ? 0 : (int)wiretick_p1_v.size()-1 ;
        int vwire_step  = (uwireindex==0) ? 1 : -1;
        int vwire_end   = (uwireindex==0) ? (int)wiretick_p1_v.size() : -1;
        float mintriarea = 1e9;
        for (int iv=vwire_index; iv!=vwire_end; iv += vwire_step ) {
          
          std::vector< int > wid_v = { (int)uwire, (int)wiretick_p1_v[iv][0], (int)ywire };
          std::vector<float> intersection_yz;
          double triarea = -1.0;
          int crosses = 0;
          ublarcvapp::UBWireTool::wireIntersection( wid_v, intersection_yz, triarea, crosses );
          LARCV_DEBUG() << " [vindex=" << iv << "] crosses=" << crosses << "  triarea=" << triarea << std::endl;
          if ( triarea>=0 && triarea < mintriarea ) mintriarea = triarea;           
          if ( crosses==1 && triarea >=0 && triarea<50.0 ) {
            // keep this intersection
            uwire_intersections.push_back( intersection_yz );
            uwire_triarea.push_back( (float)triarea );
            uwire_vindex.push_back( iv );
            uwire_wid_v.push_back( wid_v );
            if ( uwireindex==0 ) {
              std::vector<int> idxpair = { 0, iv };
              pairs_formed.insert( idxpair );
            }
            else {
              std::vector<int> idxpair = { (int)wiretick_p0_v.size()-1, iv };
              pairs_formed.insert( idxpair );
            }
            break;
          }
          vwire_tests++;
        }//end of vwire loop
        LARCV_DEBUG() << " vwire-scan: min-triarea=" << mintriarea << std::endl;      
      }
      LARCV_DEBUG() << "  uwireindex=" << uwireindex << " uwire=" << uwire << " tested-nvwires=" << vwire_tests << std::endl;
      uwireindex++;
    }//end of uwire test loop (2 entries at most)
    
    std::stringstream zy_list_u;
    for ( size_t i=0; i<uwire_intersections.size(); i++ ) {
      zy_list_u << "(" << uwire_intersections[i][0] << "," << uwire_intersections[i][1] << "; tri=" << uwire_triarea[i] << " vidx=" << uwire_vindex[i] << ") ";
      std::vector<float> tyz  = { tick, uwire_intersections[i][1], uwire_intersections[i][0] };
      points3d_v.push_back( tyz );
      std::vector<float> twid = { tick, (float)uwire_wid_v[i][0], (float)uwire_wid_v[i][1], (float)uwire_wid_v[i][2] };
      twid_v.push_back( twid );
    }
    LARCV_DEBUG() << "Build 3 plane 3D intersection with U-bounds: { " << zy_list_u.str() << "}" << std::endl;

    // v-wire
    // -------
    // we make at most 2 points, one from the top and one from the bottom
    std::vector<int> vwire_tests;
    if ( wiretick_p1_v.size()>0 ) {
      vwire_tests.push_back( (int)wiretick_p1_v.front()[0] );
      if ( wiretick_p1_v.size()>1 )
        vwire_tests.push_back( (int)wiretick_p1_v.back()[0] );
    }

    std::vector< std::vector<float> > vwire_intersections;
    std::vector< float > vwire_triarea;
    std::vector< int >   vwire_vindex;
    std::vector< std::vector<int> > vwire_wid_v;
    int vwireindex = 0;
    LARCV_DEBUG() << "forming 3d points from both directions of v-plane scan" << std::endl;
    for ( auto const& vwire : vwire_tests ) {
      int uwire_tests = 0;      
      if ( wiretick_p0_v.size()==0 && vwireindex==0 && allow_2plane_matches ) {
        int crosses = 0;
        std::vector<float> intersection;
        //ublarcvapp::UBWireTool::wireIntersection( 2, (int)ywire, 1, (int)vwire, intersection, crosses );
        int missingplane = 0;
        int missingwire = 0;
        ublarcvapp::UBWireTool::getMissingWireAndPlane( 2, (int)ywire, 1, (int)vwire, missingplane, missingwire, intersection, crosses );        
        if ( crosses==1 ) {
          vwire_intersections.push_back( intersection );
          vwire_triarea.push_back(0.0);
          vwire_vindex.push_back(0);
          vwire_wid_v.push_back( std::vector<int>{(int)missingwire,(int)vwire,(int)ywire} );
        }
        LARCV_DEBUG() << "  2-plane intersection: (z,y)=(" << intersection[0] << "," << intersection[1] << ") crosses=" << crosses << std::endl;
      }
      else if ( wiretick_p0_v.size()>0 ) {
      
        int uwire_index = (vwireindex==0) ? 0 : (int)wiretick_p0_v.size()-1 ;
        int uwire_step  = (vwireindex==0) ? 1 : -1;
        int uwire_end   = (vwireindex==0) ? (int)wiretick_p0_v.size() : -1;
        float mintriarea = 1e9;
        for (int iu=uwire_index; iu!=uwire_end; iu += uwire_step ) {

          std::vector<int> pair_idx(2);
          pair_idx[0] = iu;
          if ( vwireindex==0 )
            pair_idx[1] = 0;
          else
            pair_idx[1] = wiretick_p0_v.size()-1;
          if ( pairs_formed.find(pair_idx)!=pairs_formed.end() )
            continue;
          
          std::vector< int > wid_v = { (int)wiretick_p0_v[iu][0], (int)vwire, (int)ywire };
          std::vector<float> intersection_yz;
          double triarea = -1.0;
          int crosses = 0;
          ublarcvapp::UBWireTool::wireIntersection( wid_v, intersection_yz, triarea, crosses );
          if ( triarea>=0 && triarea < mintriarea ) mintriarea = triarea;
          if ( crosses==1 && triarea >=0 && triarea<50.0 ) {
            // keep this intersection
            vwire_intersections.push_back( intersection_yz );
            vwire_triarea.push_back( (float)triarea );
            vwire_vindex.push_back( iu );
            vwire_wid_v.push_back( wid_v );
            break;
          }
          uwire_tests++;
        }//end of vwire loop
        LARCV_DEBUG() << " uwire-scan for vwire: min-triarea=" << mintriarea << std::endl;       
      }
      LARCV_DEBUG() << "  vwireindex=" << vwireindex << " vwire=" << vwire << " tested-nuwires=" << uwire_tests << std::endl;
      vwireindex++;
    }//end of uwire test loop (2 entries at most)


    std::stringstream zy_list_v;
    for ( size_t i=0; i<vwire_intersections.size(); i++ ) {
      zy_list_v << "(" << vwire_intersections[i][0] << "," << vwire_intersections[i][1] << "; tri=" << vwire_triarea[i] << " uidx=" << vwire_vindex[i] << ") ";
      std::vector<float> tyz = { tick, vwire_intersections[i][1], vwire_intersections[i][0] };
      points3d_v.push_back( tyz );
      std::vector<float> twid = { tick, (float)vwire_wid_v[i][0], (float)vwire_wid_v[i][1],  (float)vwire_wid_v[i][2] };
      twid_v.push_back( twid );      
    }
    LARCV_DEBUG() << "Build 3 plane 3D intersection with V-bounds: { " << zy_list_v.str() << "}" << std::endl;
    
  }
  
  void GenGraphPoints::makeGraph( const std::vector< std::vector<float> >& points_tyz,  // input list of (tick,y,z) points
                                  std::vector< Eigen::Vector3f >& points,               // node list: position of (x,y,z)
                                  std::map< std::pair<int,int>, float >& distmap,       // distance between nodes
                                  std::map< std::pair<int,int>, float >& pixgapdist )   // pixel gap between vertices
  { 
                                  
    // we convert points to 'xyz' and load cilantro point cloud so that we can build kD tree
    // then we build edge information

    // map for graph edge valuesa
    points.clear();
    distmap.clear();
    pixgapdist.clear();
    
    for ( auto const& tyz : points_tyz ) {
      float x = (tyz[0]-3200)*0.5*larutil::LArProperties::GetME()->DriftVelocity();
      Eigen::Vector3f point(x, tyz[1], tyz[2] );
      // make sure its not a duplicate
      if ( points.size()>0 && points.back()==point )
        continue;
      points.emplace_back( std::move(point) );
    }

    LARCV_DEBUG() << "Made point list, size=" << points.size() << ", from input list, size=" << points_tyz.size() << "." <<  std::endl;

    // build kdtree
    cilantro::KDTree3f tree(points);

    // references to image products we'll need
    auto const& mask_v  = pfeatures->pcropdata->mask_v;
    auto const& crop_v  = pfeatures->pcropdata->crops_v;
    auto const& miss_v  = pfeatures->pcropdata->missing_v;
    auto const& badch_v = pfeatures->pcropdata->badch_v;


    // for each point, we define edge using
    //  1) distance. limit distance to nearest neighbor
    //  2) gap in pixels to point
    size_t max_n_neighbors = 5;
    float max_neighbor_dist = 50;
    float nedges_per_vertex = 0;
    
    for ( size_t idx=0; idx<points.size(); idx++ ) {
      auto& vertex = points.at(idx);
      cilantro::NeighborSet<float> nn;
      //= tree.kNNInRadiusSearch(vertex, max_n_neighbors, max_neighbor_dist );
      tree.radiusSearch( vertex, max_neighbor_dist, nn );

      int numneighbors = 0;
      for ( size_t inn=0; inn<nn.size(); inn++ ) {
        auto& neighbor = nn.at(inn);
        if ( neighbor.value==0 )
          continue; // no self-connections
        auto& neighbornode = points[ neighbor.index ];
        numneighbors++;
        float pixgap = get_max_pixelgap( crop_v, badch_v, vertex, neighbornode, 0.3, 3 );
        pixgapdist[ std::pair<int,int>( (int)idx, (int)neighbor.index) ] = pixgap;
        distmap[ std::pair<int,int>( (int)idx, (int)neighbor.index) ]    = neighbor.value;
      }
      nedges_per_vertex += (float)numneighbors;
    }
    nedges_per_vertex /= (float)points.size();
    LARCV_DEBUG() << "Defined graph nodes (size=" << points.size() << ") and edges (size=" << pixgapdist.size() << "). ave neighors per node=" << nedges_per_vertex << std::endl;

  }
  
  /**
   * calculate the pixel gap between two 3D points on a plane
   *
   */
  float GenGraphPoints::get_max_pixelgap( const std::vector<larcv::Image2D>& input_v,
                                          const std::vector<larcv::Image2D>& badch_v,
                                          const Eigen::Vector3f& start,
                                          const Eigen::Vector3f& end,
                                          const float max_steplen,
                                          const int pixel_search_width ) {

    //float max_steplen = 0.3; // cm
    //int pixel_search_width = 3;
    const size_t nplanes = input_v.size();

    // we fill the astar_path container with nodes along this line
    std::map< std::vector<int>, int > pixel_map;
    // std::vector< std::vector<int> >   pixel_list;
    // std::vector< std::vector<float> > pixel_tyz;

    // calcualte normalize distance
    std::vector<float> dir(3,0);
    float pathlen = 0.;
    for ( int i=0; i<3; i++ ) {
      dir[i] = end[i]-start[i];
      pathlen += dir[i]*dir[i];
    }
    pathlen = sqrt(pathlen);
    for ( int i=0; i<3; i++ ) dir[i] /= pathlen;

    
    int nsteps = pathlen/max_steplen;
    if ( fabs(pathlen-max_steplen*nsteps)>0.0001 )
      nsteps++;
    float steplen = pathlen/float(nsteps);

    double pos[3] = {0};
    float start_x = start[0];
    float end_x   = end[0];

    // as we step through, we want to track if the segment is good (on charge) or bad (not on charge).
    // so we have to track the state and the distance. we do this with the following struct and function
    struct SegmentState_t {
      int state; // 0=bad; 1=good;
      float curr_seg_len; // current segment length
      float max_good_seg; // max good segment seen
      float max_bad_seg;  // max bad segment seen
      std::vector<float> seg_v; // list of segs we've seen
      std::vector<int>   state_v; // list of states of the seg's we've seen
      SegmentState_t()
        : state(1), // start in good state
          curr_seg_len(0.0),
          max_good_seg(0.0),
          max_bad_seg(0.0)
      {};
      void update( int currentstate, float steplen ) {
        if ( currentstate==state ) {
          // continuation of state
          curr_seg_len += steplen;
          return;
        }
        else {
          // change of state
          if ( state==1 ) {
            // transition in the middle
            curr_seg_len += 0.5*steplen;
            // was good
            if ( curr_seg_len>max_good_seg )
              max_good_seg = curr_seg_len;
          }
          else {
            // was bad
            if ( curr_seg_len>max_bad_seg )
              max_bad_seg = curr_seg_len;
          }
          // store the seg
          seg_v.push_back( curr_seg_len );
          state_v.push_back( state );
          // set the state and start segment halfway through step
          state = currentstate;
          curr_seg_len = 0.5*steplen;
        };
      };
      float get_max_bad_seg_notfirst() {
        float maxbad = 0.;
        // skip the first, because the end point might not very good at first
        for ( size_t s=1; s<seg_v.size(); s++ ) {
          if ( seg_v[s]>maxbad )
            maxbad = seg_v[s];
        }
        return maxbad;
      };
      float get_max_bad_seg() {
        float maxbad = 0.;
        // skip the first, because the end point might not very good at first
        for ( size_t s=0; s<seg_v.size(); s++ ) {
          if ( seg_v[s]>maxbad )
            maxbad = seg_v[s];
        }
        return maxbad;
      };
    } segstate;
    
    for (int istep=1; istep<nsteps-1; istep++ ) {
      for ( int i=0; i<3; i++ )
        pos[i] = start[i] + steplen*istep*dir[i];

      // back to ticks
      float tick = pos[0]/larutil::LArProperties::GetME()->DriftVelocity()/0.5;
      if ( tick<input_v[0].meta().min_y() || tick>=input_v[0].meta().max_y() ) {
        // not in the image. seems bad.
        segstate.update(0,steplen);
        continue;
      }
      
      // get wires
      std::vector<int> rowcol_v; // 4 coordinates (row,col on each plane...)
      rowcol_v.reserve(4);
      rowcol_v.push_back( input_v[0].meta().row(tick) );
      for ( size_t pl=0; pl<nplanes; pl++ ) {
        float wirecoord = larutil::Geometry::GetME()->WireCoordinate( pos, pl );
        if ( wirecoord<input_v[pl].meta().min_x() || wirecoord>=input_v[pl].meta().max_x() ) {
          break;
        }
        rowcol_v.push_back( (int)input_v[pl].meta().col( wirecoord ) );
      }
      
      if ( rowcol_v.size()!=(1+nplanes) ) {
        // bad point by out of bounds step
        segstate.update(0,steplen);
        continue;
      }

      // we have a good point in the image
      auto it = pixel_map.find( rowcol_v );
      if ( it!=pixel_map.end() ) {
        // already in the set, no need to revaluate
        // all points in the set are good (good idea?)
        segstate.update(it->second,steplen);
        continue;
      }

      // new location in the image
      // get the pixel values
      int nplanes_w_charge = 0;
      int nplanes_w_badch  = 0;
      std::vector<float> step_tyz = { tick, (float)pos[1], (float)pos[2] };

      for ( size_t pl=0; pl<nplanes; pl++ ) {

        float maxpixval = 0.;
        for ( int dc=-pixel_search_width; dc<=pixel_search_width; dc++ ) {
          int c = rowcol_v[pl+1]+dc;
          if ( c<0 || c>=(int)input_v[pl].meta().cols() ) continue; // skip it
          float pixval = input_v[pl].pixel( rowcol_v[0], c );
          if ( pixval>maxpixval )
            maxpixval = pixval;
        }
        
        if ( maxpixval>10.0 )
          nplanes_w_charge++;
        else {
          // check badch
          if ( badch_v[pl].pixel( rowcol_v[0], rowcol_v[pl+1]) > 0.5 )
            nplanes_w_badch++;
        }
      }
      
      if ( nplanes_w_charge==3 ) {
        segstate.update(1,steplen);
        pixel_map[rowcol_v] = 1;
        // pixel_list.push_back( rowcol_v );
        // pixel_tyz.push_back( step_tyz );
      }
      else {
        // not charge complete
        if ( nplanes_w_badch==1 ) {
          segstate.update(1,steplen);
          pixel_map[rowcol_v] = 1;
          // pixel_list.push_back( rowcol_v );
          // pixel_tyz.push_back( step_tyz );          
        }
      }

      // everything else is bad
      pixel_map[rowcol_v] = 0;
      segstate.update(0,steplen);
      // pixel_list.push_back( rowcol_v );
      // pixel_tyz.push_back( step_tyz );
    }//end of step loop

    // make astar node list
    return  segstate.get_max_bad_seg();
  }

  /*
  void GenGraphPoints::shortestpath( const std::vector< Eigen::Vector3f >& points,            // node list: position of (x,y,z)
                                     const std::map< std::pair<int,int>, float >& distmap,    // distance between nodes
                                     const std::map< std::pair<int,int>, float >& pixgapdist) // pixel gap between vertices
  {

    struct VertexData {
      int index;
      int pred;
      int dist;
    };

    struct EdgeData {
      float dist;
    };

    // define boost graph
    typedef boost::adjacency_list<boost::vecS, boost::vecS,
                                  boost::undirectedS,
                                  VertexData,
                                  EdgeData> MyGraphType;
    // root node
    int v0 = 0;
    
  }
  */
}
}
