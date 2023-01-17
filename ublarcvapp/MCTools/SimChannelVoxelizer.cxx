#include "SimChannelVoxelizer.h"

#include "larlite/LArUtil/DetectorProperties.h"

#include "MCPixelPGraph.h"

namespace ublarcvapp {
namespace mctools {

  SimChannelVoxelizer::SimChannelVoxelizer()
  {
    auto const& geo  = *(larlite::larutil::Geometry::GetME());
    for (int icryo=0; icryo<(int)geo.Ncryostats(); icryo++) {
      for (int itpc=0; itpc<(int)geo.NTPCs(icryo); itpc++) {
	makeDefaultTPCVoxels( icryo, itpc );
      }
    }
  }

  /**
   * @brief clear data memebers that store event data
   *
   */
  void SimChannelVoxelizer::clear()
  {
    for ( auto& tpcinfo : _tpcdata_v ) {
      tpcinfo.clear_event_data();
    }
      
  }

  void SimChannelVoxelizer::makeDefaultTPCVoxels( int cryoid, int tpcid )
  {
    auto const& geo  = *(larlite::larutil::Geometry::GetME());
    //auto const& detp = *(larutil::DetectorProperties::GetME());
    auto const& tpcgeo = geo.GetTPC( tpcid, cryoid );

    TPCInfo info;
    info.cryoid = cryoid;
    info.tpcid  = tpcid;
    info.driftdir = geo.TPCDriftDir( tpcid, cryoid )[0];

    // default voxel dimensions
    info._voxel_dim_cm_v.resize(3,0.3);

    // replace first dim as tick width
    int dtick = info._voxel_dim_cm_v[0]/larutil::DetectorProperties::GetME()->GetXTicksCoefficient();
    info._voxel_dim_cm_v[0] = dtick;
    info._voxel_dim_cm_v[0] = 1;    

    // define drift direction (always x?)
    // define readout tick bounds and effect x-axis
    float tickmin = 0;
    float tickmax = getReadoutNumSamples( larutil::LArUtilConfig::Detector() );
    float ymin = tpcgeo.fBounds[0][1];
    float ymax = tpcgeo.fBounds[1][1];    
    float zmin = tpcgeo.fBounds[0][2];
    float zmax = tpcgeo.fBounds[1][2];

    // set origin
    info._origin_cm_v.resize(3,0);
    info._origin_cm_v[0] = tickmin;    
    for (int i=1; i<3; i++)
      info._origin_cm_v[i] = tpcgeo.fBounds[0][i];

    // provide x-position of anode
    if ( info.driftdir>0 )
      info.realx_anode = tpcgeo.fBounds[1][0]; // e- moves in +x direction, so x-max is anode
    else
      info.realx_anode = tpcgeo.fBounds[0][0]; // e- moves in -x direction, so x-min is anode

    info.vx_trigger = float(getReadoutTriggerOffset( larutil::LArUtilConfig::Detector() ))/info._voxel_dim_cm_v[0];
    
    // chop into 3 mm blocks
    long nx = fabs(tickmax-tickmin)/info._voxel_dim_cm_v[0]+1;
    long ny = fabs(ymax-ymin)/0.3+1;
    long nz = fabs(zmax-zmin)/0.3+1;
    
    //long totN = nx*ny*nz;
    // std::cout << "[SimChannelVoxelizer::Constructor] total tensor size: " << totN << std::endl;
    // std::cout << "CRYOID[" << cryoid << "] TPC[" << tpcid << "]" << std::endl;
    // std::cout << "NX: " << nx << std::endl;
    // std::cout << "NY: " << ny << std::endl;
    // std::cout << "NZ: " << nz << std::endl;
    // std::cout << "dtick: " << info._voxel_dim_cm_v[0] << " tpc tick" << std::endl;
    // std::cout << "dy: " << info._voxel_dim_cm_v[1] << " cm" << std::endl;
    // std::cout << "dz: " << info._voxel_dim_cm_v[2] << " cm" << std::endl;    
    // std::cout << "tick-min=" << tickmin << " tick-max=" << tickmax << " driftdir=" << info.driftdir << std::endl;
    // std::cout << "y-min=" << ymin << " y-max=" << ymax << std::endl;
    // std::cout << "z-min=" << zmin << " z-max=" << zmax << std::endl;
    // std::cout << "vx-trigger=" << info.vx_trigger << std::endl;
    
    info._num_voxels_v.resize(3,0);
    info._num_voxels_v[0] = nx;
    info._num_voxels_v[1] = ny;
    info._num_voxels_v[2] = nz;
    
    _tpcdata_v.emplace_back( std::move(info) );
    int index = (int)_tpcdata_v.size()-1;
    std::pair<int,int> ctid = {cryoid,tpcid};
    _id2index_v[ ctid ] = index;
  }

  int SimChannelVoxelizer::getReadoutNumSamples( larlite::geo::DetId_t geoid )
  {
    if ( geoid==larlite::geo::kMicroBooNE )
      return 6048;
    return larutil::DetectorProperties::GetME()->NumberTimeSamples();
  }

  int SimChannelVoxelizer::getReadoutTriggerOffset( larlite::geo::DetId_t geoid )
  {
    if ( geoid==larlite::geo::kMicroBooNE )
      return 3200;
    return larutil::DetectorProperties::GetME()->TriggerOffset();
  }

  int SimChannelVoxelizer::getReadoutFirstTick( larlite::geo::DetId_t geoid )
  {
    if ( geoid==larlite::geo::kMicroBooNE )
      return 2400;
    return 0;
  }
  
  void SimChannelVoxelizer::process( larlite::storage_manager& ioll )
  {
    larlite::event_simch* ev_simch
      = (larlite::event_simch*)ioll.get_data( larlite::data::kSimChannel, "simdrift" );
    larlite::event_mctrack* ev_mctrack
      = (larlite::event_mctrack*)ioll.get_data( larlite::data::kMCTrack, "mcreco" );
    larlite::event_mcshower* ev_mcshower
      = (larlite::event_mcshower*)ioll.get_data( larlite::data::kMCShower, "mcreco" );
    larlite::event_mctruth* ev_mctruth
      = (larlite::event_mctruth*)ioll.get_data( larlite::data::kMCTruth, "generator" );

    process( *ev_simch, *ev_mctrack, *ev_mcshower, *ev_mctruth );
  }
  
  void SimChannelVoxelizer::process( const larlite::event_simch& ev_simch,
				     const larlite::event_mctrack& mctrack,
				     const larlite::event_mcshower& mcshower,
				     const larlite::event_mctruth& mctruth)
  {
    // std::cout << "[SimChannelVoxelizer] start process()" << std::endl;
    // std::cout << " mctrack.size()=" << mctrack.size() << std::endl;
    // std::cout << " mcshower.size()=" << mcshower.size() << std::endl;
    // std::cout << " mctruth.size()=" << mctruth.size() << std::endl;    
    
    _scan_IDEs( ev_simch );

    _pdg_labels( mctrack, mcshower, mctruth );
    
    _make_tensors();
  }

  void SimChannelVoxelizer::_scan_IDEs( const larlite::event_simch& ev_simch )
  {

    size_t nbadide = 0;
    size_t ninvalid = 0;
    size_t notpc = 0;
    size_t nvalid = 0;
    for ( auto const& xsimch : ev_simch ) {

      unsigned short chid = xsimch.Channel();
      auto wid = larlite::larutil::Geometry::GetME()->ChannelToWireID( chid );
      int cryoid = wid.Cryostat;
      int tpcid  = wid.TPC;
	
      for ( auto it = xsimch.TDCIDEMap().begin(); it!=xsimch.TDCIDEMap().end(); it++ ) {
	//unsigned short tick = it->first;
	const std::vector<larlite::ide>& idc_v = it->second;
	//std::cout << "IDE len=" << idc_v.size() << std::endl;
	
	for ( auto const& xide : idc_v ) {
	  //std::cout << "IDE pos=(" << xide.x << "," << xide.y << "," << xide.z << ")" << std::endl;
	  if ( fabs(xide.x)>1e100 || fabs(xide.y)>1e100 || fabs(xide.z)>1e100 ) {
	    //std::cout << " bad IDE" << std::endl;
	    nbadide++;
	    continue;
	  }
	  // determine the voxel to populate
	  TVector3 pos( xide.x, xide.y, xide.z );
	  std::pair<int,int> ctid_pair( cryoid, tpcid );

	  auto it_tpc = _id2index_v.find( ctid_pair );
	  if ( it_tpc==_id2index_v.end() ) {
	    //std::cout << "No corresponding TPC?" << std::endl;
	    notpc++;
	    continue;
	  }
	  //std::cout << "  Found info for TPC: index=" << it_tpc->second << std::endl;
	  // if ( xide.trackID!=xide.originID )

	  auto& tpcinfo = _tpcdata_v.at( it_tpc->second );
	  //std::cout << " tpc[" << tpcinfo.cryoid << "," << tpcinfo.tpcid << "]" << std::endl;

	  std::vector<int> voxel(3,0);
	  bool valid = true;
	  for (int i=1; i<3; i++) {
	    voxel[i] = (pos[i]-tpcinfo._origin_cm_v[i])/tpcinfo._voxel_dim_cm_v[i];
	    if ( voxel[i]<0 || voxel[i]>=tpcinfo._num_voxels_v[i] )
	      valid = true;
	  }

	  int tdc = it->first;
	  float tick = tdc-3000;
	  voxel[0] = tick/tpcinfo._voxel_dim_cm_v[0];
	  if ( voxel[0]<0 || voxel[0]>=tpcinfo._num_voxels_v[0] )
	    valid = true;

	  // std::cout << "IDE TDC=" << it->first << " TID=" << abs(xide.trackID)
	  // 	    << " CH[" << chid << "] cryo=" << cryoid << " tpc=" << tpcid
	  // 	    << " valid=" << valid << " voxel[0]=" << voxel[0]
	  // 	    << std::endl;
	  
	  if ( valid ) {
	    //auto const& dim = tpcinfo._charge_v.shape;
	    //int voxindex = voxel[0]*dim[1]*dim[2] + voxel[1]*dim[2] + voxel[2];
	    VoxelCoord_t vcoord;
	    vcoord.vtick = voxel[0];
	    vcoord.vy = voxel[1];
	    vcoord.vz = voxel[2];
	    vcoord.tdc = it->first;
	    //std::cout << "  fill valid voxel: (" << voxel[0] << "," << voxel[1] << "," << voxel[2] << ")" << std::endl;
	    auto it_vox = tpcinfo._voxcoord_2_index.find( vcoord );
	    if ( it_vox==tpcinfo._voxcoord_2_index.end() ) {
	      tpcinfo._voxcoord_2_index[vcoord] = tpcinfo._voxfeat_v.size();
	      it_vox = tpcinfo._voxcoord_2_index.find( vcoord );
	      VoxelFeat_t newfeat;
	      tpcinfo._voxfeat_v.push_back(newfeat);
	    }
	    //std::cout << "  feature index=" << it_vox->second << std::endl;
	    auto& feat = tpcinfo._voxfeat_v.at( it_vox->second );
	    if ( feat.charge<xide.energy ) {
	      feat.charge = xide.energy;
	      feat.trackid = abs(xide.trackID);
	      feat.realpos[0] = xide.x;
	      feat.realpos[1] = xide.y;
	      feat.realpos[2] = xide.z;
	      feat.realpos[3] = 0.0;
	    }
	    //tpcinfo._charge_v.data[voxindex] = xide.energy;
	    nvalid++;
	  }
	  else {
	    ninvalid++;
	  }
	  
	}//end of loop over ide
      }//end of loop over IDE map
    }//end of loop over simch

    // std::cout << "Number of valid IDE: " <<  nvalid << std::endl;
    // std::cout << "Number of bad IDE: " << nbadide << std::endl;
    // std::cout << "Number of Invalid IDE: " << ninvalid << std::endl;
    // std::cout << "Number not in TPC: " << notpc << std::endl;

    return;
  }

  void SimChannelVoxelizer::_make_tensors()
  {

    // Now populate the TPC info arrays
    for ( auto& tpcinfo : _tpcdata_v ) {
     
      int nvoxels = tpcinfo._voxcoord_2_index.size();
      std::cout << "[SimChannelVoxelizer] tpc has " << nvoxels << " filled voxels" << std::endl;
      
      std::vector<int> coord_dims = { nvoxels, 3 };
      std::vector<int> feat_dims = { nvoxels };
      
      tpcinfo._coordpos_v.ndims = 2;
      tpcinfo._coordpos_v.shape = coord_dims;      
      tpcinfo._coordpos_v.data.resize( nvoxels*3, 0 );
      
      tpcinfo._coordindex_v.ndims = 2;
      tpcinfo._coordindex_v.shape = coord_dims;
      tpcinfo._coordindex_v.data.resize( nvoxels*3, 0 );

      tpcinfo._charge_v.ndims = 1;
      tpcinfo._charge_v.shape = feat_dims;
      tpcinfo._charge_v.data.resize( nvoxels, 0 );

      tpcinfo._pdg_v.ndims = 1;
      tpcinfo._pdg_v.shape = feat_dims;
      tpcinfo._pdg_v.data.resize( nvoxels, 0 );
      
      tpcinfo._trackid_v.ndims = 1;
      tpcinfo._trackid_v.shape = feat_dims;
      tpcinfo._trackid_v.data.resize( nvoxels, 0 );
      
      tpcinfo._ancestorid_v.ndims = 1;
      tpcinfo._ancestorid_v.shape = feat_dims;
      tpcinfo._ancestorid_v.data.resize( nvoxels, 0 );
      
      //std::cout << "allocated space for voxels. fill the arrays" << std::endl;

      size_t index = 0;
      for ( auto it=tpcinfo._voxcoord_2_index.begin(); it!=tpcinfo._voxcoord_2_index.end(); it++ ) {
	tpcinfo._coordindex_v.data[3*index + 0] = it->first.vtick;
	tpcinfo._coordindex_v.data[3*index + 1] = it->first.vy;
	tpcinfo._coordindex_v.data[3*index + 2] = it->first.vz;
	for (int i=0; i<3; i++)
	  tpcinfo._coordpos_v.data[3*index + i] = tpcinfo._voxfeat_v.at( it->second ).realpos[i];
	tpcinfo._charge_v.data[index]  = tpcinfo._voxfeat_v.at( it->second ).charge;
	tpcinfo._trackid_v.data[index] = tpcinfo._voxfeat_v.at( it->second ).trackid;
	tpcinfo._ancestorid_v.data[index] = tpcinfo._voxfeat_v.at( it->second ).ancestorid;
	tpcinfo._pdg_v.data[index] = tpcinfo._voxfeat_v.at( it->second ).pdg;		
	index++;
      }
      
    }//end of loop over info for each TPC
    
  }

  void SimChannelVoxelizer::_pdg_labels( const larlite::event_mctrack& mctrack,
					 const larlite::event_mcshower& mcshower,
					 const larlite::event_mctruth& mctruth)
  {
    
    //auto const pdetp = larutil::DetectorProperties::GetME();

    //std::cout << "Parse MCPixelPGraph" << std::endl;
    MCPixelPGraph mcpg;
    //mcpg.set_verbosity(::larcv::msg::kDEBUG);
    mcpg.buildgraphonly( mcshower, mctrack, mctruth );
    //mcpg.printGraph( nullptr, false );
    size_t no_ancestor = 0;
    for ( auto& tpcinfo : _tpcdata_v ) {

      //std::cout << "CRYO[" << tpcinfo.cryoid << "] TPC[" << tpcinfo.tpcid << "] t0-shift and label voxels" << std::endl;
      
      for (auto it=tpcinfo._voxcoord_2_index.begin(); it!=tpcinfo._voxcoord_2_index.end(); it++ ) {

	auto& vox = tpcinfo._voxfeat_v.at( it->second );

	auto pnode = mcpg.findTrackID( vox.trackid );
	if ( pnode==NULL ) {
	  // no node, no change to voxel data
	  vox.pdg = 11;
	  vox.ancestorid = 0;
	  no_ancestor++;
	  //std::cout << "NO NODE: TID=" << vox.trackid << std::endl;
	}
	else {
	  // get ancestor
	  //auto anode = mcpg.findTrackID( pnode->aid );
	  vox.pdg = pnode->pid;
	  vox.ancestorid = pnode->aid;
	}//end of if TrackID found
      }//end of loop over voxelcoord map iterator

      
    }//end of loop over TPC info

    std::cout << "[SimChannelVoxelizer::_pdg_labels] number w/ no ancestor: " << no_ancestor << std::endl;
  }
}
}
