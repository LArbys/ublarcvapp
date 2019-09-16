#include "DLTagger.h"

#include "larcv/core/DataFormat/ROI.h"

#include "ublarcvapp/ubdllee/FixedCROIFromFlashConfig.h"
#include "ublarcvapp/ubdllee/FixedCROIFromFlashAlgo.h"
#include "ublarcvapp/ContourTools/ContourClusterAlgo.h"

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

    // make track objects
    _makeTracks( wholeview_v, _mask_match_algo, m_track_v );

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
    m_track_v.clear();
    _mask_match_algo.clear();
    hasRun = false;
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
   * Collect pixels for each match combo
   *
   */
  void DLTagger::_makePixelClusters( const std::vector<larcv::Image2D>& wholeview_v,
                                     const MRCNNMatch& matchdata,
                                     std::vector< std::vector<larcv::Pixel2DCluster> >& pixel_cluster_vv,
                                     std::vector< larcv::ImageMeta >& pixel_cluster_meta_v )
  {

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

      auto& cropcombo  = matchdata.m_combo_crops_v[icombo];
      auto& astarcombo = matchdata.m_combo_astar_v[icombo];
      auto& endptcombo = matchdata.m_combo_endpt3d_v[icombo];
      auto& graphcombo = matchdata.m_combo_graphpts_v[icombo];

      // for now, skip on tracks that did not complete
      if ( matchdata.m_pass[icombo]==0 ) {
        LARCV_DEBUG() << "  reco considered failure. skip." << std::endl;
        for ( size_t p=0; p<wholeview_v.size(); p++ ) {
          pixel_cluster_vv[p].push_back( larcv::Pixel2DCluster() );
        }
        continue;
      }

      // we use the graph path unless there is none and we have the astar path
      // if we have to use the astar path, we get a list of tick,wid
      std::vector< std::vector<float> > astar_twid_v;
      if ( graphcombo.m_path_twid_v.size()==0 ) {
        LARCV_DEBUG() << " fill astar_twid_v" << std::endl;
        const std::vector<reco3d::AStar3DNode>& astar_path = astarcombo.astar_path;
        for ( auto const& node : astar_path ) {
          std::vector<float> node_twid(4,0);
          // convert to wire in crop coord, then to whole image wid
          node_twid[0] = cropcombo.crops_v[0].meta().pos_y( node.row );
          for (int p=0; p<3; p++ ) {
            node_twid[p+1] = cropcombo.crops_v[p].meta().pos_x( node.cols[p] );
          }
          astar_twid_v.push_back( node_twid );
        }
      }

      std::vector< larcv::Pixel2DCluster > pix_v(wholeview_v.size());
      const std::vector< std::vector<float> >* path_twid_v = nullptr;
      if ( graphcombo.m_path_twid_v.size()>0 )
        path_twid_v = &graphcombo.m_path_twid_v;
      else
        path_twid_v = &astar_twid_v;
      

      for ( int ipt=0; ipt<(int)path_twid_v->size()-1; ipt++ ) {

        auto const& twid1_v = (*path_twid_v)[ipt];
        auto const& twid2_v = (*path_twid_v)[ipt+1];
        if ( twid1_v.size()!=wholeview_v.size()+1 || twid2_v.size()!=wholeview_v.size()+1 ) {
          LARCV_CRITICAL() << "number of wire coordinates does not equal the number of planes!" << std::endl;
          throw larcv::larbys();
        }


        float tick1 = twid1_v[0];
        float tick2 = twid2_v[0];
        if ( tick1<wholeview_v[0].meta().min_y() || tick1>=wholeview_v[0].meta().max_y() )
          continue;
        if ( tick2<wholeview_v[0].meta().min_y() || tick2>=wholeview_v[0].meta().max_y() )
          continue;
        
        for (size_t p=0; p<wholeview_v.size(); p++ ) {
          auto const& meta = wholeview_v[p].meta();

          float wid1 = twid1_v[p+1];
          float wid2 = twid2_v[p+1];
          if ( wid1<meta.min_x() || wid1>=meta.max_x() ) continue;
          if ( wid2<meta.min_x() || wid2>=meta.max_x() ) continue;
          
          // move from point to the next
          float dir[2] = { 0, 0 };
          dir[0] = (float)meta.col(twid2_v[p+1])-(float)meta.col(twid1_v[p+1]);
          dir[1] = (float)meta.row(twid2_v[0])-(float)meta.row(twid1_v[0]); // {wire,tick}
          float dirlen = sqrt( dir[0]*dir[0] + dir[1]*dir[1] );
          if (dirlen==0) continue;
          for (int i=0; i<2; i++ ) dir[i] /= dirlen;
          
          int nsteps = dirlen/0.15 + 1;
          float stepsize = dirlen/float(nsteps);
          
          for (int istep=0; istep<=nsteps; istep++ ) {
            float col = (float)meta.col(twid1_v[p+1]) + istep*stepsize*dir[0];
            float row = (float)meta.row(twid1_v[0]) + istep*stepsize*dir[1];
            
            for (int dr=-2; dr<=2; dr++ ) {
              int r = (int)row+dr;
              if ( r<0 || r>=(int)meta.rows() ) continue;
              for (int dc=-2; dc<=2; dc++ ) {
                int c = (int)col+dc;
                if ( c<0 || c>=(int)meta.cols() ) continue;
                
                float pixval = wholeview_v[p].pixel(r,c);
                if ( pixval>10.0 ) {
                  larcv::Pixel2D pix( c, r );
                  pix.Intensity( pixval );
                  pix_v[p] += pix;
                }
              }//end of col neighorhood loop
            }//end of row neighborhood
          }//end of steps
          
        }//end of plane loop
        
      }//end of loop over path segements

      for (size_t p=0; p<3; p++ ) {
        pixel_cluster_vv[p].emplace_back( std::move(pix_v[p]) );
      }
            
    }//end of loop over matches

    LARCV_INFO() << "number of pixel clusters made: " << pixel_cluster_vv.front().size() << std::endl;
    
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
      if ( matchdata.m_pass[imatch]==0 ) {
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
  
  
  // /**
  //  * get export track info from reco'd match clusters as larlite::track.
  //  *
  //  * 
  //  */
  // void DLTagger::_getTracks( const MRCNNMatch& matchdata,
  //                            larlite::event_track& track, int track_type ) {
    
  //   auto& graphdata_v = matchdata.m_combo_graphpts_v;
  //   auto& astardata_v = matchdata.m_combo_astar_v;
    
  //   for ( size_t i=0; matchdata.m_pass.size(); i++ ) {
  //     if ( matchdata.m_pass[i]==0 )
  //       continue; // failed reco
      
  //     if ( track_type==0 && is_cosmic_v[i]==0 ) continue; // only save is-cosmic
  //     nclusters[p]++;
      
  //     auto& cluster = m_pixel_cluster_vv[p][icluster];
  //     if ( m_iscosmic_v[icluster]==1 ) {
  //       ev_cosmic.Emplace( (larcv::PlaneID_t)p, std::move(cluster), m_pixel_cluster_meta_v[p] );
  //     }
  //     else {
  //       ev_notcosmic.Emplace( (larcv::PlaneID_t)p, std::move(cluster), m_pixel_cluster_meta_v[p] );
  //     }
  //   }
    
  //   for ( auto& ncluster : nclusters ) {
  //     if ( ncluster!=nclusters.front() ) {
  //       LARCV_CRITICAL() << "The number of clusters in each plane is not the same!" << std::endl;
  //       throw larcv::larbys( "The number of clusters in each plane is not the same!" );
  //     }
  //   }
  //   hasRun = false;
  //   m_pixel_cluster_vv.clear();
  //   m_pixel_cluster_meta_v.clear();
  // }
    
  
  /**
   * we reconstruct the core of the track
   * 
   * we attempt to the run the Reco3D tracker. (actually later. let's make a simple track for now to eval vertexer.)
   * We pass in an image where a wide neighorhood path along the existing hits from the graph and/or astar
   * is provided to help the track jump the gaps.
   */
  void DLTagger::_makeTracks( const std::vector<larcv::Image2D>& wholeview_v,
                              const MRCNNMatch& matchdata,
                              larlite::event_track& ev_track_v )
  {

    LARCV_DEBUG() << "make pixel clusters" << std::endl;
    
    int nplanes = wholeview_v.size();    
    int ncombos = matchdata.m_combo_3plane_v.size();
    int ntracks = matchdata.m_combo_graphpts_v.size();
    
    // we use the matches with a good AStar track
    for (int icombo=0; icombo<ntracks; icombo++ ) {

      LARCV_DEBUG() << "--------------------------------------" << std::endl;
      LARCV_DEBUG() << "tagging with combo[" << icombo << "] from Graph/AStar Reco" << std::endl;

      auto& cropcombo  = matchdata.m_combo_crops_v[icombo];
      auto& astarcombo = matchdata.m_combo_astar_v[icombo];
      auto& endptcombo = matchdata.m_combo_endpt3d_v[icombo];
      auto& graphcombo = matchdata.m_combo_graphpts_v[icombo];

      // for now, skip on tracks that did not complete
      if ( matchdata.m_pass[icombo]==0 ) {
        LARCV_DEBUG() << "  reco considered failure. skip." << std::endl;
        ev_track_v.push_back( larlite::track() );
        continue;
      }

      // we use the graph path unless there is none and we have the astar path
      // if we have to use the astar path, we get a list of tick,wid
      std::vector< std::vector<float> > astar_twid_v;
      std::vector< std::vector<float> > astar_xyz_v;
      if ( graphcombo.m_path_twid_v.size()==0 ) {
        LARCV_DEBUG() << " fill astar_twid_v" << std::endl;
        const std::vector<reco3d::AStar3DNode>& astar_path = astarcombo.astar_path;
        for ( auto const& node : astar_path ) {
          std::vector<float> node_twid(4,0);
          // convert to wire in crop coord, then to whole image wid
          node_twid[0] = cropcombo.crops_v[0].meta().pos_y( node.row );
          for (int p=0; p<3; p++ ) {
            node_twid[p+1] = cropcombo.crops_v[p].meta().pos_x( node.cols[p] );
          }
          astar_twid_v.push_back( node_twid );
          astar_xyz_v.push_back( node.tyz );
          astar_xyz_v.back()[0] = (astar_xyz_v.back()[0]-3200)*0.5*larutil::LArProperties::GetME()->DriftVelocity();
        }
      }

      const std::vector< std::vector<float> >* path_twid_v = nullptr;
      const std::vector< std::vector<float> >* path_xyz_v  = nullptr;
      if ( graphcombo.m_path_twid_v.size()>0 ) {
        path_twid_v = &graphcombo.m_path_twid_v;
        path_xyz_v  = &graphcombo.m_path_xyz;
      }
      else {
        path_twid_v = &astar_twid_v;
        path_xyz_v  = &astar_xyz_v;
      }


      // make a basic larlite track
      larlite::track track;
      track.set_track_id( icombo );
      track.reserve( path_xyz_v->size() );
      for ( size_t ipos=0; ipos<path_xyz_v->size(); ipos++ ) {
        auto& xyz = (*path_xyz_v)[ipos];
        TVector3 pos( xyz[0], xyz[1], xyz[2] );
        track.add_vertex( pos );        

        float dirlen = 0.;
        std::vector<float> dir(3,0);
        if ( ipos+1<path_xyz_v->size() ) {
          auto& xyznext = (*path_xyz_v)[ipos+1];
          for ( int i=0; i<3; i++ ) {
            dir[i] = xyznext[i]-xyz[i];
            dirlen += dir[i]*dir[i];
          }
        }
        else if ( ipos>=1 ) {
          auto& xyzprev = (*path_xyz_v)[ipos-1];
          for ( int i=0; i<3; i++ ) {
            dir[i] = xyz[i]-xyzprev[i];
            dirlen += dir[i]*dir[i];
          }
          dirlen = sqrt(dirlen);
          if ( dirlen>0 )
            for ( int i=0; i<3; i++ ) dir[i] /= dirlen;
        }
        TVector3 dirv( dir[0], dir[1], dir[2] );
        track.add_direction( dirv );

      }//loop over points

      LARCV_DEBUG() << " track[" << (int)ev_track_v.size()-1 << "] created from combo[" << icombo  << "]" << std::endl;
      ev_track_v.emplace_back( std::move(track) );
      
    }//end of loop over combo
    
    LARCV_INFO() << "prepared " << ev_track_v.size() << " tracks" << std::endl;
  }


  /**
   * get a copy of the larlite tracks for clusters labeled as "cosmic".
   *
   * @param[inout] ev_track larlite::track event container for tracks to be passed into.
   */
  void DLTagger::getCosmicTracks( larlite::event_track& ev_track ) {
    _getTracks(0,ev_track);
  }

  /**
   * get a copy of the larlite tracks for clusters labeled as "not-cosmic".
   *
   * @param[inout] ev_track larlite::track event container for tracks to be passed into.
   */
  void DLTagger::getNotCosmicTracks( larlite::event_track& ev_track ) {
    _getTracks(1,ev_track);
  }

  /**
   * get a copy of the larlite tracks for clusters labeled as "good reco'd".
   *
   *  this is the union of cosmic and not-cosmic tracks.
   *
   * @param[inout] ev_track larlite::track event container for tracks to be passed into.
   */
  void DLTagger::getAllGoodTracks( larlite::event_track& ev_track ) {
    _getTracks(2,ev_track);
  } 
 
  /**
   * get a copy of the larlite tracks for different labels.
   *
   * @param[in]    track_type Type of track to be returned: 0=cosmic, 1=notcosmic, 2=all
   * @param[inout] ev_track larlite::track event container for tracks to be passed into.
   */  
  void DLTagger::_getTracks( int track_type, larlite::event_track& ev_track ) {    

    if ( m_track_v.size()!=m_iscosmic_v.size() ) {
      LARCV_CRITICAL() << "Number of tracks (" << m_track_v.size() << ") "
                       << " does not match the number of iscosmic/notcosmic tags (" << m_iscosmic_v.size() << ")"
                       << std::endl;
      std::stringstream msg;
      msg << "Number of tracks (" << m_track_v.size() << ") "
          << " does not match the number of iscosmic/notcosmic tags (" << m_iscosmic_v.size() << ")"
          << std::endl;
      throw larcv::larbys(msg.str());
    }

    if ( track_type<0 || track_type>2 ) {
      LARCV_CRITICAL() << "Type of track request is unrecognized" << std::endl;
      throw larcv::larbys();
    }

    int ntracks_passed = 0;
    for ( size_t itrack=0; itrack<m_track_v.size(); itrack++ ) {

      auto const& track = m_track_v[itrack];

      if (track.NumberTrajectoryPoints()==0)
        continue;
      
      if ( track_type==0 && m_iscosmic_v[itrack]==1 ) {
        ev_track.push_back( track );
        ntracks_passed++;
      }
      else if ( track_type==1 && m_iscosmic_v[itrack]==0 ) {
        ev_track.push_back( track );
        ntracks_passed++;
      }
      else if (track_type==2) {
        ev_track.push_back( track );        
        ntracks_passed++;
      }
    }
    std::string str_track_type;
    if ( track_type==0 )
      str_track_type = "cosmic";
    else if (track_type==1)
      str_track_type = "notcosmic";
    else
      str_track_type = "all reco'd";
    
    LARCV_DEBUG() << "Number of (larlite) tracks returned as '" << str_track_type << "': "
                  << ntracks_passed
                  << std::endl;
  }
  
}
}
    
