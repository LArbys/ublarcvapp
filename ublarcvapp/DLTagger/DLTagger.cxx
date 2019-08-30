#include "DLTagger.h"

#include "larcv/core/DataFormat/ROI.h"

#include "ublarcvapp/ubdllee/FixedCROIFromFlashConfig.h"
#include "ublarcvapp/ubdllee/FixedCROIFromFlashAlgo.h"

#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "Geo2D/Core/Geo2D.h"

namespace ublarcvapp {
namespace dltagger {

  void DLTagger::runTagger( const std::vector<larcv::Image2D>& wholeview_v,
                            const larcv::EventChStatus& ev_chstatus,
                            const std::vector<larlite::opflash>& intime_opflash_v,
                            const std::vector< std::vector<larcv::ClusterMask> >& clustermask_vv ) {

    // clear containers
    reset();
    
    // -------------------------------------------
    // CROI
    // -------------------------------------------
    ublarcvapp::ubdllee::FixedCROIFromFlashConfig croi_config;// rely on defaults for now, need proper config
    ublarcvapp::ubdllee::FixedCROIFromFlashAlgo   croi_algo(croi_config);

    LARCV_DEBUG() << "define CROI and make merged CROI" << std::endl;
    int nflashes = intime_opflash_v.size();
    if ( nflashes>0 ) {
      m_croi_v      = croi_algo.findCROIfromFlash( intime_opflash_v.front() ); ///< 3d-consistent regions
      m_croi_merged = _mergeROI( wholeview_v.size(), m_croi_v );
    }
    else {
      m_croi_v.clear();
    }
    
    // -------------------------------------------
    // Mask-RCNN Cluster Matching Across planes
    // -------------------------------------------
    
    // reset/clear algo
    _mask_match_algo.clear();
    _mask_match_algo.set_verbosity( logger().level() );

    // make matches
    LARCV_DEBUG() << "run MRCNNMatch::matchMasksAcrossPlanes" << std::endl;
    bool use_gap_ch = true;
    _mask_match_algo.matchMasksAcrossPlanes( clustermask_vv, wholeview_v, ev_chstatus, use_gap_ch );

    // make pixel clusters
    _makePixelClusters( wholeview_v, _mask_match_algo, m_pixel_cluster_vv, m_pixel_cluster_meta_v );

    // select clusters as cosmic or not cosmic
    _evaluateClusters( wholeview_v, _mask_match_algo, m_croi_merged, m_pixel_cluster_vv, m_iscosmic_v );
    
    // make thrumu image
    _tagPixels( wholeview_v, m_pixel_cluster_vv, m_iscosmic_v, m_tagged_v, m_notcosmic_v );

    // mark that we have collected made all produces
    hasRun = true;

  }

  /**
   * reset all containers
   *
   */
  void DLTagger::reset() {
    m_pixel_cluster_vv.clear();
    m_pixel_cluster_meta_v.clear();
    m_tagged_v.clear();
    m_notcosmic_v.clear();
    m_croi_v.clear();
    m_iscosmic_v.clear();
    m_select_vars_v.clear();
    _mask_match_algo.clear();
  }

  /**
   * transfer images to the given container
   *
   */
  void DLTagger::transferImages( std::vector<larcv::Image2D>& cosmic_v,
                                 std::vector<larcv::Image2D>& notcosmic_v ) {
    
    if (m_tagged_v.size()==0) {
      LARCV_CRITICAL() << "Asking for tagged images, but have not run tagger yet!" << std::endl;
      return;
    }
    
    for ( auto& img : m_tagged_v ) {
      cosmic_v.emplace_back( std::move(img) );
    }
    for ( auto& img : m_notcosmic_v ) {
      notcosmic_v.emplace_back( std::move(img) );
    }
    hasRun = false;
    m_tagged_v.clear();
    m_notcosmic_v.clear();
  }

  /**
   * transfer pixels to event container
   */
  void DLTagger::transferPixelClusters( larcv::EventPixel2D& ev_cosmic,
                                        larcv::EventPixel2D& ev_notcosmic ) {

    if (m_pixel_cluster_vv.size()==0) {
      LARCV_CRITICAL() << "Asking for pixel clusters, but have not run tagger yet!" << std::endl;
      return;
    }

    std::vector<int> nclusters( m_pixel_cluster_vv.size(), 0 );
    for ( size_t p=0; p<m_pixel_cluster_vv.size(); p++ ) {
      for ( size_t icluster=0; icluster<m_pixel_cluster_vv[p].size(); icluster++ ) {
        
        nclusters[p]++;
        
        auto& cluster = m_pixel_cluster_vv[p][icluster];
        if ( m_iscosmic_v[icluster]==1 ) {
          ev_cosmic.Emplace( (larcv::PlaneID_t)p, std::move(cluster), m_pixel_cluster_meta_v[p] );
        }
        else {
          ev_notcosmic.Emplace( (larcv::PlaneID_t)p, std::move(cluster), m_pixel_cluster_meta_v[p] );
        }
        
      }
    }
    
    for ( auto& ncluster : nclusters ) {
      if ( ncluster!=nclusters.front() ) {
        LARCV_CRITICAL() << "The number of clusters in each plane is not the same!" << std::endl;
        throw larcv::larbys( "The number of clusters in each plane is not the same!" );
      }
    }
    hasRun = false;
    m_pixel_cluster_vv.clear();
    m_pixel_cluster_meta_v.clear();
  }

  /**
   * transfer croi to event containers
   */
  void DLTagger::transferCROI( larcv::EventROI& ev_croi_v, larcv::EventROI& ev_croi_merged ) {
    // these can be empty
    if ( m_croi_v.size()==0 )
      return;
    
    ev_croi_v.Emplace( std::move(m_croi_v) );
    ev_croi_merged.Emplace( std::move(m_croi_merged) );
    m_croi_v.clear();
  }

  /**
   * Collect pixels for each matche combo
   *
   */
  void DLTagger::_makePixelClusters( const std::vector<larcv::Image2D>& wholeview_v,
                                     const MRCNNMatch& matchdata,
                                     std::vector< std::vector<larcv::Pixel2DCluster> >& pixel_cluster_vv,
                                     std::vector< larcv::ImageMeta >& pixel_cluster_meta_v ) {

    LARCV_DEBUG() << "make pixel clusters" << std::endl;
    
    int nplanes = wholeview_v.size();
    
    int ncombos = matchdata.m_combo_3plane_v.size();
    int nastar  = matchdata.m_combo_astar_v.size();

    // container for individual image crops
    std::vector< std::vector<larcv::Image2D> >  plane_fill_images_v( nplanes );

    // container for pixel clusters
    pixel_cluster_vv.clear();
    pixel_cluster_meta_v.clear();
    pixel_cluster_vv.resize(wholeview_v.size());
    // save whole image meta for them
    for ( auto const& whole : wholeview_v ) {
      pixel_cluster_meta_v.push_back( whole.meta() );
    }
    
    // we use the matches with a good AStar track
    for (int icombo=0; icombo<nastar; icombo++ ) {

      LARCV_DEBUG() << "--------------------------------------" << std::endl;
      LARCV_DEBUG() << "tagging with combo[" << icombo << "] from AStar list" << std::endl;
      
      auto& astarcombo = matchdata.m_combo_astar_v[icombo];
      auto& endptcombo = matchdata.m_combo_endpt3d_v[icombo];

      // for now, skip on tracks that did not complete
      if ( astarcombo.astar_completed==0 ) {
        LARCV_DEBUG() << "astar not complete. fill empty pixel cluster" << std::endl;
        // create an empty clusters
        for ( size_t p=0; p<wholeview_v.size(); p++ ) {
          pixel_cluster_vv[p].push_back( larcv::Pixel2DCluster() );
        }
        continue;
      }
      
      // we want to follow the track and pick up pixels within
      //  charge contour clusters close to the path

      // mask image: above threshold + mrcnn mask
      // auto const& plane_contours_vv
      //   = _mask_match_algo.m_combo_features_v[icombo].combo_mask_contour.m_plane_atomics_v;
      // auto const& plane_contour_metas_vv
      //   = _mask_match_algo.m_combo_features_v[icombo].combo_mask_contour.m_plane_atomicmeta_v;
      
      // charge image: above threshold
      auto const& plane_contours_vv
        = _mask_match_algo.m_combo_features_v[icombo].combo_charge_contour.m_plane_atomics_v;
      auto const& plane_contour_metas_vv
        = _mask_match_algo.m_combo_features_v[icombo].combo_charge_contour.m_plane_atomicmeta_v;

      const std::vector<reco3d::AStar3DNode>& path = astarcombo.astar_path;

      for ( size_t p=0; p<plane_contours_vv.size(); p++ ) {

        const larcv::ImageMeta* cropmeta = &matchdata.m_combo_crops_v[icombo].crops_v[p].meta();
        if ( cropmeta->cols()==0 || cropmeta->rows()==0 )
          cropmeta = &_mask_match_algo.m_combo_crops_v[icombo].missing_v[p].meta();

        // the image we fill with contour pixels. around crop
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
          
          // step in the smallest dimension (as long as its not zero)
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
        auto& fullmeta = wholeview_v[p].meta();

        larcv::Pixel2DCluster pixcluster;
        
        for ( int ictr=0; ictr<ncontours; ictr++ ) {
          if ( contourkept[ictr]==0 ) continue;

          auto const& contourmeta = contourmeta_v[ictr];
          auto const& qpixels = contourmeta.getChargePixels();


          int nearpath = 0;
          int farpath  = 0;
          for ( auto const& pix2d : qpixels ) {
            int col = pix2d.x;
            int row = pix2d.y;
            if ( col<0 || col>=(int)cropmeta->cols() ) continue;
            if ( row<0 || row>=(int)cropmeta->rows() ) continue;

            // pixels need to stay close to the line (in pixel coords)
            // if astar, find closest segment and then distance to segment
            float distfrompath = 0.;
            if ( astarcombo.used_astar ) {
              // astar method
              int closest_segment = -1;
              float min_seg_dist = 0;
              for ( int inode=0; inode<(int)astarcombo.astar_path.size()-1; inode++ ) {
                auto const& nodeA = astarcombo.astar_path.at(inode);
                auto const& nodeB = astarcombo.astar_path.at(inode+1);
                float pixdistA = sqrt((col-nodeA.cols[p])*(col-nodeA.cols[p]) + (row-nodeA.row)*(row-nodeA.row));
                float pixdistB = sqrt((col-nodeB.cols[p])*(col-nodeB.cols[p]) + (row-nodeB.row)*(row-nodeB.row));
                float pixdist = pixdistA + pixdistB;
                if ( closest_segment<0 || pixdist<min_seg_dist ) {
                  closest_segment = inode;
                  min_seg_dist = pixdist;
                }
              }

              // define line segment and point in Geo2D, then get dist to line
              auto const& closenodeA = astarcombo.astar_path[closest_segment];
              auto const& closenodeB = astarcombo.astar_path[closest_segment+1];              
              geo2d::LineSegment<float> segment( closenodeA.cols[p], closenodeA.row, closenodeB.cols[p], closenodeB.row );
              geo2d::Vector<float> pt( col, row );
              distfrompath = geo2d::Distance( segment, pt );
            }
            else {
              // straight line method
              geo2d::LineSegment<float> segment( astarcombo.astar_path.front().cols[p], astarcombo.astar_path.front().row,
                                                 astarcombo.astar_path.back().cols[p], astarcombo.astar_path.back().row );
              geo2d::Vector<float> pt( col, row );
              distfrompath = geo2d::Distance( segment, pt );
            }

            if ( distfrompath>_max_pixdist_from_path ) {
              farpath++;
              continue;
            }
            else
              nearpath++;
            
            float wire = cropmeta->pos_x(col);
            float tick = cropmeta->pos_y(row);

            if ( wire<fullmeta.min_x() || wire>=fullmeta.max_x() ) continue;
            if ( tick<fullmeta.min_y() || tick>=fullmeta.max_y() ) continue;

            int xcol = fullmeta.col(wire);
            int xrow = fullmeta.row(tick);
            
            fillimg.set_pixel( row, col,  1.0 );

            larcv::Pixel2D pix(xcol,xrow);
            pix.Intensity( wholeview_v[p].pixel( xrow, xcol ) );
            pixcluster += pix;
          }
          LARCV_DEBUG() << "  kept contour[" << ictr << "] "
                        << " (pixels close to path)=" << nearpath
                        << " (pixels far from path)=" << farpath
                        << std::endl;
        }
        LARCV_DEBUG() << "-- pixels in cluster plane[" << p << "]: " << pixcluster.size() << std::endl;
        plane_fill_images_v[p].emplace_back( std::move(fillimg) );
        pixel_cluster_vv[p].emplace_back( std::move(pixcluster) );
        
      }//end of loop over planes
      
    }//end of loop over matches
    
  }
  
  /**
   * tag whole image pixels using pixelclusters built from mrcnn matches
   *
   * fill images for both cosmic- and non-cosmic- selected clusters
   *
   */
  void DLTagger::_tagPixels( const std::vector<larcv::Image2D>& wholeview_v,
                             const std::vector< std::vector<larcv::Pixel2DCluster> >& pixel_cluster_vv,
                             const std::vector<int>& iscosmic_v,
                             std::vector<larcv::Image2D>& whole_cosmic_v,
                             std::vector<larcv::Image2D>& whole_notcosmic_v ) {
    
    whole_cosmic_v.clear();
    whole_notcosmic_v.clear();    
    for ( auto const& img : wholeview_v ) {
      larcv::Image2D cosmic(img.meta());
      larcv::Image2D noncosmic(img.meta());
      cosmic.paint(0);
      noncosmic.paint(0);
      whole_cosmic_v.emplace_back( std::move(cosmic) );
      whole_notcosmic_v.emplace_back( std::move(noncosmic) );
    }

    // sanity/assumption checks
    if ( wholeview_v.size()!=pixel_cluster_vv.size() )
      LARCV_CRITICAL() << "number of planes in whole image vector and mrcnn-mask-matched cluster vector does not agree" << std::endl;
    
    for ( auto const& pixel_cluster_v : pixel_cluster_vv ) {
      if ( pixel_cluster_v.size()!=iscosmic_v.size() ) {
        std::stringstream errmsg;
        errmsg << "number of clusters in pixel_cluster_v (" << pixel_cluster_v.size() << ") "
               << "does not match the number of iscosmic flags (" << iscosmic_v.size() << ")"
               << std::endl;
        LARCV_CRITICAL() << errmsg.str() << std::endl;
        throw larcv::larbys(errmsg.str());
      }
    }

    // loop over pixel clusters, mark up images
    for ( size_t p=0; p<wholeview_v.size(); p++ ) {

      auto const& pixel_cluster_v = pixel_cluster_vv[p];
      
      auto& cosmic    = whole_cosmic_v[p];
      auto& noncosmic = whole_notcosmic_v[p];

      for ( size_t icluster=0; icluster<pixel_cluster_v.size(); icluster++ ) {

        auto const& cluster = pixel_cluster_v[icluster];

        larcv::Image2D* fillimage = (iscosmic_v[icluster]==1) ? &cosmic : &noncosmic;

        for ( auto const& pix : cluster )
          fillimage->set_pixel( pix.Y(), pix.X(), 10.0 );
      }//end of loop over clusters in plane
    }//end of loop over planes
                            
  }  

  /**
   * we select clusters that are consistent with cosmics.
   *
   * we avoid as best we can clusters that are consistent with neutrinos.
   * cuts are based on
   *  1) Detector coordinates containment
   *  2) location out of time
   *  3) location out of the CROI
   *
   */
  void DLTagger::_evaluateClusters( const std::vector<larcv::Image2D>& wholeview_v,
                                    const MRCNNMatch& matchdata,
                                    const larcv::ROI& croimerged,
                                    const std::vector< std::vector<larcv::Pixel2DCluster> >& pixel_cluster_vv,
                                    std::vector<int>& iscosmic_v ) {

    LARCV_DEBUG() << "evaluate clusters" << std::endl;
    
    size_t nmatches = matchdata.numMatches();
    iscosmic_v.resize(nmatches,0);

    LARCV_DEBUG() << "merged CROI given: " << std::endl;
    for ( auto const& bb : croimerged.BB() ) {
      LARCV_DEBUG() << "  " << bb.dump() << std::endl;
    }
    
    // apply selection
    for ( size_t imatch=0; imatch<nmatches; imatch++ ) {

      // store varialbles for analysis
      // ------------------------------
      CosmicSelectVars_t vars;
        
      // for now, skip on tracks that did not complete
      auto& astarcombo = matchdata.m_combo_astar_v[imatch];      
      if ( astarcombo.astar_completed==0 ) {
        m_select_vars_v.push_back( vars );        
        continue;
      }

      vars.astar_complete = 1;      
      
      auto const& endptdata = matchdata.m_combo_endpt3d_v.at(imatch);
      
      // measure of containment
      _calcDwall( endptdata, vars.dwall_outermost, vars.dwall_innermost,
                  vars.outermost_endpt_tyz, vars.innermost_endpt_tyz );

      // measure within intime
      vars.dtick_outoftime = _calcOutOfTimeTicks( endptdata );

      // measure fraction out of CROI
      std::vector< const larcv::Pixel2DCluster* > pixel_v( wholeview_v.size(), nullptr );
      for ( size_t p=0; p<wholeview_v.size(); p++ ) {
        pixel_v[p] = &pixel_cluster_vv[p][imatch];
        LARCV_DEBUG() << "match[" << imatch << "] plane[" << p << "] number of pixels in cluster: " << pixel_v[p]->size() << std::endl;
      }
      if ( croimerged.BB().size()>0 ) {
        _calcOutOfROIfrac( wholeview_v, croimerged, pixel_v, vars.frac_per_plane, vars.numpixels, vars.total_frac );
      }
      else {
        vars.frac_per_plane.resize( wholeview_v.size(), 0.0 );
        vars.total_frac = 0.;
      }

      // store selection vars
      m_select_vars_v.push_back( vars );

      LARCV_DEBUG() << "---------------------------------------" << std::endl;      
      LARCV_DEBUG() << "MatchCluster [" << imatch << "]" << std::endl;
      LARCV_DEBUG() << "---------------------------------------" << std::endl;
      LARCV_DEBUG() << "  dTick Out of Time: " << vars.dtick_outoftime << std::endl;
      LARCV_DEBUG() << "  Total Frac out of ROI: " << vars.total_frac << std::endl;
      for ( size_t p=0; p<vars.frac_per_plane.size(); p++ )
        LARCV_DEBUG() << "     plane[" << p << "] fraction: " << vars.frac_per_plane[p] << std::endl;
      LARCV_DEBUG() << "  DWALL: innermost=" << vars.dwall_innermost << "  outermost=" << vars.dwall_outermost << std::endl;

      
      // make decision
      // ---------------

      if ( vars.dtick_outoftime<_cut_dtick_outoftime ) {
        //iscosmic_v.push_back(1);
        iscosmic_v[imatch] = 1;
        LARCV_DEBUG() << "  cosmic-tagged using [out-of-time]" << std::endl;
        continue;
      }

      if ( vars.total_frac>_cut_frac_out_of_croi ) {
        //iscosmic_v.push_back(1);
        iscosmic_v[imatch] = 1;
        LARCV_DEBUG() << "  cosmic-tagged using [out-of-croi]" << std::endl;        
        continue;
      }

      // containment can save it
      if ( vars.dwall_innermost<_cut_frac_dwall_innermost || vars.dwall_outermost<_cut_frac_dwall_outermost ) {
        //iscosmic_v.push_back(0);
        iscosmic_v[imatch] = 1;
        LARCV_DEBUG() << "  cosmic-tagged using [dwall]" << std::endl;
      }
      else {
        //iscosmic_v.push_back(1);
        iscosmic_v[imatch] = 0;
        LARCV_DEBUG() << "  not-cosmic using [dwall]" << std::endl;        
      }

    }//end of loop over matches/clusters
    
  }
  
  /**
   * evalute dwall using the clusters' 3d endpoints
   *
   */
  void DLTagger::_calcDwall( const Gen3DEndpoints& endpts,
                             float& outermost_dwall,
                             float& innermost_dwall,
                             std::vector<float>& outermost_endpt_tyz,
                             std::vector<float>& innermost_endpt_tyz ) {
    
    const float tpc_halfheight = larutil::Geometry::GetME()->DetHalfHeight();
    const float tpc_length     = larutil::Geometry::GetME()->DetLength();
    LARCV_DEBUG() << "tpc_halfheight=" << tpc_halfheight << " tpc_length=" << tpc_length << std::endl;
    std::vector<float> dwall_endpt_v;
    dwall_endpt_v.reserve(3);
    for ( auto const& tyz : endpts.endpt_tyz_v ) {
      float dwally[2] = { tpc_halfheight-tyz[1], tyz[1]-(-tpc_halfheight) }; // Y top and bottom
      float dwallz[2] = { tyz[2], tpc_length-tyz[2] }; // Z min and max
      float dwall_endpt = dwally[0];
      dwall_endpt = ( dwall_endpt>dwally[1] ) ? dwally[1] : dwall_endpt;
      dwall_endpt = ( dwall_endpt>dwallz[0] ) ? dwallz[0] : dwall_endpt;
      dwall_endpt = ( dwall_endpt>dwallz[1] ) ? dwallz[1] : dwall_endpt;
      LARCV_DEBUG() << "endpoint det. (y,z)=(" << tyz[1] << "," << tyz[2] << "); dwall=" << dwall_endpt << std::endl;
      
      dwall_endpt_v.push_back( dwall_endpt );
    }
    
    //outermost_dwall = ( dwall_endpt_v[0]<dwall_endpt_v[1] )  dwall_endpt_v[0] : dwall_endpt_v[1];
    //innermost_dwall = ( dwall_endpt_v[0]>dwall_endpt_v[1] ) ? dwall_endpt_v[0] : dwall_endpt_v[1];
    if ( dwall_endpt_v[0]<dwall_endpt_v[1] ) {
      outermost_dwall = dwall_endpt_v[0];
      outermost_endpt_tyz = endpts.endpt_tyz_v[0];
      innermost_dwall = dwall_endpt_v[1];
      innermost_endpt_tyz = endpts.endpt_tyz_v[1];
    }
    else {
      outermost_dwall = dwall_endpt_v[1];
      outermost_endpt_tyz = endpts.endpt_tyz_v[1];
      innermost_dwall = dwall_endpt_v[0];
      innermost_endpt_tyz = endpts.endpt_tyz_v[0];
    }
  }
  
  /**
   * calculate the length of ticks the track extends out of the in-time window ticks
   *
   * should be between [3200-pad, 3200+drifttime+pad]
   *
   */
  float DLTagger::_calcOutOfTimeTicks( const Gen3DEndpoints& endpts ) {
    const float tick_pad = 10.0; 
    const float tick_min = 3200-tick_pad;
    const float tick_max = 3200 + larutil::Geometry::GetME()->DetHalfWidth()*2/larutil::LArProperties::GetME()->DriftVelocity()/0.5 + tick_pad;

    float dtick[2] = { endpts.endpt_tyz_v[0][0]-tick_min, tick_max-endpts.endpt_tyz_v[1][0] };
    if ( dtick[0]<dtick[1] ) return dtick[0];
    else return dtick[1];
    return 0; // never gets here
  }

  /**
   * merge vector of ROIs
   *
   */
  larcv::ROI DLTagger::_mergeROI( const int nplanes, const std::vector<larcv::ROI>& croi_v ) {
    
    larcv::ROI mergedroi;
    std::vector<larcv::ImageMeta> meta_v;
    for ( int p=0; p<nplanes; p++ ) {
      larcv::ImageMeta meta;
      for ( auto const& roi : croi_v ) {
        if (meta.cols()==0 || meta.rows()==0 ) {
          // empty, so just set it
          meta = roi.BB(p);
        }
        else {
          meta = meta.inclusive(roi.BB(p)); // union, meta
        }
      }
      LARCV_INFO() << "CROI-merged plane[" << p << "]: " << meta.dump() << std::endl;
      meta_v.push_back( meta );
    }
    mergedroi.SetBB( meta_v );
    return mergedroi;
  }
  
  /** calculate fraction out of CROI
   * 
   */
  void DLTagger::_calcOutOfROIfrac( const std::vector<larcv::Image2D>& wholeview_v,
                                    const larcv::ROI& croi,
                                    const std::vector<const larcv::Pixel2DCluster*>& ppixel_cluster_v,
                                    std::vector<float>& frac_per_plane,
                                    std::vector<int>& npix,
                                    float& total_frac ) {
    size_t nplanes = ppixel_cluster_v.size();
    frac_per_plane.resize(nplanes,0.0);
    total_frac = 0;
    float total_pix = 0;

    npix.resize(nplanes,0);

    for ( size_t p=0; p<nplanes; p++ ) {
      auto const& meta = wholeview_v[p].meta();
      auto const& bb   = croi.BB(p);
      for ( auto const& pix : *ppixel_cluster_v[p] ) {
        npix[p]++;
        try {
          float tick = meta.pos_y(pix.Y());
          float wire = meta.pos_x(pix.X());
          if ( !bb.contains( wire, tick ) ) {
            frac_per_plane[p] += 1.0;
          }
        }
        catch (std::exception& e) {
          LARCV_WARNING() << "exception testing pixel[" << npix[p] << "]: " << e.what() << std::endl;
          continue;
        }
      }

      if ( npix[p]>0 ) {
        total_frac += frac_per_plane[p];        
        frac_per_plane[p] /= (float)npix[p];
        total_pix += (float)npix[p];
      }

    }//end of plane loop
    if ( total_pix>0 )
      total_frac /= total_pix;
  }

  /**
   * add truth-info to reco products for analysis
   *
   *
   */
  void DLTagger::recoTruthMatching( const std::vector<larcv::Image2D>& mcinstance_v ) {
    
    if (!hasRun) {
      LARCV_CRITICAL() << "Need to run tagger first before running reco-truth matching" << std::endl;
      throw larcv::larbys( "Need to run Tagger first" );
    }

    int numclusters = m_pixel_cluster_vv.front().size();
    for ( int icluster=0; icluster<numclusters; icluster++ ) {
      // loop over the pixels and total up pixels with an instance ID. means its a neutrino.
      auto& vars = m_select_vars_v[icluster];
      vars.total_nufrac = 0.;
      vars.nufrac_per_plane.resize( m_pixel_cluster_vv.size(), 0 );
      for ( size_t p=0; p<m_pixel_cluster_vv.size(); p++ ) {
        auto const& instance_img = mcinstance_v[p];
        auto const& pixel_v = m_pixel_cluster_vv.at(p).at(icluster);
        for ( auto const& pix : pixel_v ) {
          if ( (int)pix.Y()>=0 && (int)pix.Y()<(int)instance_img.meta().rows()
               && (int)pix.X()>=0 && (int)pix.X()<(int)instance_img.meta().cols() ) {
            if ( instance_img.pixel( (int)pix.Y(), (int)pix.X() )>0 ) {
              vars.nufrac_per_plane[p]++;
              vars.total_nufrac += 1.0;
            }
          }
        }//loop over pixels
        if ( vars.numpixels[p]>0 ) {
          vars.nufrac_per_plane[p] /= (float)vars.numpixels[p];
        }
        else {
          vars.nufrac_per_plane[p] = 0.;
        }
      }//loop over planes
      float totpix = vars.numpixels[0] + vars.numpixels[1] + vars.numpixels[2];
      if ( totpix>0 ) {
        vars.total_nufrac /= totpix;
      }
      else {
        vars.total_nufrac = 0.;
      }
    }// end of cluster loop
    
  }
  
}
}
    
