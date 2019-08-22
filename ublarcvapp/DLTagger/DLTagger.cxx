#include "DLTagger.h"

namespace ublarcvapp {
namespace dltagger {

  void DLTagger::runTagger( const std::vector<larcv::Image2D>& wholeview_v,
                            const larcv::EventChStatus& ev_chstatus,
                            const std::vector< std::vector<larcv::ClusterMask> >& clustermask_vv ) {

    // reset/clear algo
    _mask_match_algo.clear();

    // make matches
    _mask_match_algo.matchMasksAcrossPlanes( clustermask_vv, wholeview_v, ev_chstatus );

    // make thrumu image
    _tagPixels( wholeview_v, m_tagged_v );
    hasRun = true;
  }

  /**
   * transfer images to the given container
   *
   */
  void DLTagger::transferImages( std::vector<larcv::Image2D>& out ) {
    for ( auto& img : m_tagged_v ) {
      out.emplace_back( std::move(img) );
    }
    hasRun = false;
  }

  /**
   * tag pixels using matches and Astar tracks above
   *
   */
  void DLTagger::_tagPixels( const std::vector<larcv::Image2D>& wholeview_v,
                             std::vector<larcv::Image2D>& whole_filled_v ) {

    int nplanes = wholeview_v.size();
    
    int ncombos = _mask_match_algo.m_combo_3plane_v.size();
    int nastar  = _mask_match_algo.m_combo_astar_v.size();

    // container for image crops
    std::vector< std::vector<larcv::Image2D> >  plane_fill_images_v( nplanes );

    // container for full images
    whole_filled_v.clear();
    for ( size_t p=0; p<(size_t)nplanes; p++ ) {
      auto const& img = wholeview_v[p];
      auto const& fullmeta = img.meta();
      larcv::Image2D fillme( img.meta() );
      fillme.paint(0);
      whole_filled_v.emplace_back( std::move(fillme) );
    }
    
    // we use the matches with a good AStar track
    for (int icombo=0; icombo<nastar; icombo++ ) {

      LARCV_DEBUG() << "tagging with combo[" << icombo << "] of naster" << std::endl;
      
      auto& astarcombo = _mask_match_algo.m_combo_astar_v[icombo];

      // for now, skip on tracks that did not complete
      if ( astarcombo.astar_completed==0 )
        continue;
      
      // we want to follow the track and pick up pixels within
      //  charge contour clusters close to the path

      // mask
      // auto const& plane_contours_vv
      //   = _mask_match_algo.m_combo_features_v[icombo].combo_mask_contour.m_plane_atomics_v;
      // auto const& plane_contour_metas_vv
      //   = _mask_match_algo.m_combo_features_v[icombo].combo_mask_contour.m_plane_atomicmeta_v;
      // charge
      auto const& plane_contours_vv
        = _mask_match_algo.m_combo_features_v[icombo].combo_charge_contour.m_plane_atomics_v;
      auto const& plane_contour_metas_vv
        = _mask_match_algo.m_combo_features_v[icombo].combo_charge_contour.m_plane_atomicmeta_v;

      std::vector<reco3d::AStar3DNode>& path = astarcombo.astar_path;

      for ( size_t p=0; p<plane_contours_vv.size(); p++ ) {

        const larcv::ImageMeta* cropmeta = &_mask_match_algo.m_combo_crops_v[icombo].crops_v[p].meta();
        if ( cropmeta->cols()==0 || cropmeta->rows()==0 )
          cropmeta = &_mask_match_algo.m_combo_crops_v[icombo].missing_v[p].meta();

        // the image we fill with contour pixels
        larcv::Image2D fillimg(*cropmeta);
        fillimg.paint(0.0);
        LARCV_DEBUG() << "fill image: " << cropmeta->dump() << std::endl;
        
        auto const& contours_v    = plane_contours_vv[p];
        auto const& contourmeta_v = plane_contour_metas_vv[p];
        std::vector<int> contourkept( contours_v.size(), 0 );
        int ncontours = contours_v.size();
        int naccepted = 0;

        LARCV_DEBUG() << "contours in image: " << ncontours << std::endl;
        
        for ( size_t inode=0; inode<path.size()-1; inode++ ) {
          auto& node     = path[inode];
          auto& nextnode = path[inode+1];

          LARCV_DEBUG() << "node[" << inode << "]" << std::endl;          
          // want step end points in pixel coordinates (same coordinates as contours)
          float endpt[2]    = { (float)node.cols[p],     (float)cropmeta->row(node.tyz[0]) };
          float startpt[2]  = { (float)nextnode.cols[p], (float)cropmeta->row(nextnode.tyz[0]) };
          LARCV_DEBUG() << "node[" << inode << "] start=(" << startpt[0] << "," << startpt[1] << ") "
                        << "end=(" << endpt[0] << "," << endpt[1] << ")" << std::endl;
          
          // get dir and length between nodes
          float stepdir[2] = {0};
          float steplen = 0.;
          for (int i=0; i<2; i++ ) {
            stepdir[i] = endpt[i]-startpt[i];
            steplen += stepdir[i]*stepdir[i];
          }
          steplen = sqrt(steplen);
          for ( int i=0; i<2; i++ ) stepdir[i] /= steplen;
          
          // step in the smallest dimension (as long as its not zero
          float step = 1.0;
          if ( stepdir[0]==0 )      step = 1.0/stepdir[1];
          else if ( stepdir[1]==0 ) step = 1.0/stepdir[0];
          else {
            step = ( stepdir[0]>stepdir[1] ) ? 1.0/stepdir[0] : 1.0/stepdir[1];
          }

          // set number of steps
          int nsteps = (steplen/step);

          // because of the image downsampling used by Astar,
          // we could be (or are usually?) far off the path initially
          // we define a bounding box around the step with some padding
          float minx = ( startpt[0] < endpt[0] ) ? startpt[0] : endpt[0];
          float maxx = ( startpt[0] > endpt[0] ) ? startpt[0] : endpt[0];
          float miny = ( startpt[1] < endpt[1] ) ? startpt[1] : endpt[1];
          float maxy = ( startpt[1] > endpt[1] ) ? startpt[1] : endpt[1];

          // expand bounding box by downsample factor
          maxx += 16;
          maxy += 16*6;

          // test against all the contours, against bounding box first
          std::vector<int> test_indices;
          test_indices.reserve( contours_v.size() ); // reserve maximum size
          for ( size_t ictr=0; ictr<contours_v.size(); ictr++ ) {
            if ( contourkept[ictr]==1 ) continue; // already kept, no need to check further
            
            auto const& contour = contours_v[ictr];
            auto const& contourmeta = contourmeta_v[ictr];
            
            // we test bounding box first for quick test
            if ( maxx < contourmeta.getMinX() ||
                 minx > contourmeta.getMaxX() ||
                 maxy < contourmeta.getMinY() ||
                 miny > contourmeta.getMaxY() )
              continue;
            
            test_indices.push_back(ictr);
          }//end of loop over contours

          LARCV_DEBUG() << "node[" << inode << "] contours to test: " << test_indices.size() << std::endl;
          
          // nothing to test
          if ( test_indices.size()==0 )
            continue;
          
          // now we step along line, testing as we go
          std::vector<int> accepted(test_indices.size(),0);
          int ntestok = 0;
          for (int istep=0; istep<nsteps; istep++ ) {
            float x = startpt[0] + istep*step*stepdir[0];
            float y = startpt[1] + istep*step*stepdir[1];

            // loop over contours to test
            for ( size_t itest=0; itest<test_indices.size(); itest++ ) {
              if ( accepted[itest]==1 ) continue; // no need to check further
              float dist = cv::pointPolygonTest( contours_v[ test_indices[itest] ], cv::Point(x,y), true );
              if ( fabs(dist)<16.0 ) {
                accepted[itest] = 1;
                contourkept[ test_indices[itest] ] = 1;
                ntestok++;
                naccepted++;
              }
              // all accepted, stop
              if ( ntestok==test_indices.size() )
                break;
            }
            
            /// all accepted, stop
            if ( ntestok==test_indices.size() )
              break;
          }//end of loop over step between node points

          LARCV_DEBUG() << "node[" << inode << "] contours accepted " << ntestok << std::endl;          

          // we accepted all the contours in this crop. done!
          if ( naccepted==ncontours )
            break;
          
        }// loop over node points

        // got accepted contours
        // mark up image

        // contour metas have the charge pixels conveniently collected
        auto& fillme = whole_filled_v[p];
        auto& fullmeta = fillme.meta();
        
        for ( int ictr=0; ictr<ncontours; ictr++ ) {
          if ( contourkept[ictr]==0 ) continue;

          auto const& contourmeta = contourmeta_v[ictr];
          auto const& qpixels = contourmeta.getChargePixels();

          for ( auto const& pix2d : qpixels ) {
            int col = pix2d.x;
            int row = pix2d.y;
            if ( col<0 || col>=(int)cropmeta->cols() ) continue;
            if ( row<0 || row>=(int)cropmeta->rows() ) continue;

            float wire = cropmeta->pos_x(col);
            float tick = cropmeta->pos_y(row);

            if ( wire<fullmeta.min_x() || wire>=fullmeta.max_x() ) continue;
            if ( tick<fullmeta.min_y() || tick>=fullmeta.max_y() ) continue;

            int xcol = fullmeta.col(wire);
            int xrow = fullmeta.row(tick);

            fillme.set_pixel( xrow, xcol, 1.0 );
            fillimg.set_pixel( row, col,  1.0 );
          }
          
        }

        plane_fill_images_v[p].emplace_back( std::move(fillimg) );
        
      }//end of loop over planes
      
    }//end of loop over matches
    
  }
  
}
}
    
