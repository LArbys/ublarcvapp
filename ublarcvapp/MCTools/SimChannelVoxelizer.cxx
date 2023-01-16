#include "SimChannelVoxelizer.h"

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

  void SimChannelVoxelizer::makeDefaultTPCVoxels( int cryoid, int tpcid )
  {
    auto const& geo  = *(larlite::larutil::Geometry::GetME());
    auto const& detp = *(larutil::DetectorProperties::GetME());
    auto const& tpcgeo = geo.GetTPC( tpcid, cryoid );

    TPCInfo info;
    info.cryoid = cryoid;
    info.tpcid  = tpcid;

    // define drift direction (always x?)
    // define readout tick bounds and effect x-axis
    float xmin = detp.ConvertTicksToX( 0, 0, tpcid, cryoid );
    float xmax = detp.ConvertTicksToX( getReadoutNumSamples( larutil::LArUtilConfig::Detector() ), 0, tpcid, cryoid );
    float ymin = tpcgeo.fBounds[0][1];
    float ymax = tpcgeo.fBounds[1][1];    
    float zmin = tpcgeo.fBounds[0][2];
    float zmax = tpcgeo.fBounds[1][2];

    // set origin
    info._origin_cm_v.resize(3,0);
    info._origin_cm_v[0] = xmin;    
    for (int i=1; i<3; i++)
      info._origin_cm_v[i] = tpcgeo.fBounds[0][i];

    // default voxel dimensions
    info._voxel_dim_cm_v.resize(3,0.3);
    
    // chop into 3 mm blocks
    int nx = (xmax-xmin)/0.3+1;
    int ny = (ymax-ymin)/0.3+1;
    int nz = (zmax-zmin)/0.3+1;

    std::cout << "[SimChannelVoxelizer::Constructor] total tensor size: " << nx*ny*nz << std::endl;
    
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
    std::cout << "[SimChannelVoxelizer] start process()" << std::endl;
    _scan_IDEs( ev_simch );

    _pdg_labels( mctrack, mcshower, mctruth );
    
    _make_tensors();
  }

  void SimChannelVoxelizer::_scan_IDEs( const larlite::event_simch& ev_simch )
  {

    size_t nbadide = 0;
    for ( auto const& xsimch : ev_simch ) {

      for ( auto it = xsimch.TDCIDEMap().begin(); it!=xsimch.TDCIDEMap().end(); it++ ) {
	//unsigned short tick = it->first;
	const std::vector<larlite::ide>& idc_v = it->second;
	//std::cout << "IDE len=" << idc_v.size() << std::endl;
	
	for ( auto const& xide : idc_v ) {
	  //std::cout << "IDE pos=(" << xide.x << "," << xide.y << "," << xide.z << ")" << std::endl;
	  if ( xide.x>1e100 || xide.y>1e100 || xide.z>1e100 ) {
	    //std::cout << " bad IDE" << std::endl;
	    nbadide++;
	    continue;
	  }
	  // determine the voxel to populate
	  TVector3 pos( xide.x, xide.y, xide.z );
	  std::vector<int> ctid = larlite::larutil::Geometry::GetME()->GetContainingCryoAndTPCIDs( pos );
	  //std::cout << "  CTID: " << ctid[0] << " " << ctid[1] << std::endl;
	  std::pair<int,int> ctid_pair( ctid[0], ctid[1] );

	  auto it_tpc = _id2index_v.find( ctid_pair );
	  if ( it_tpc==_id2index_v.end() ) {
	    //std::cout << "No corresponding TPC?" << std::endl;
	    continue;
	  }
	  //std::cout << "  Found info for TPC: index=" << it_tpc->second << std::endl;

	  auto& tpcinfo = _tpcdata_v.at( it_tpc->second );
	  //std::cout << " tpc[" << tpcinfo.cryoid << "," << tpcinfo.tpcid << "]" << std::endl;

	  std::vector<int> voxel(3,0);
	  bool valid = true;
	  for (int i=0; i<3; i++) {
	    voxel[i] = (pos[i]-tpcinfo._origin_cm_v[i])/tpcinfo._voxel_dim_cm_v[i];
	    if ( voxel[i]<0 || voxel[i]>=tpcinfo._num_voxels_v[i] )
	      valid = false;
	  }

	  if ( valid ) {
	    //auto const& dim = tpcinfo._charge_v.shape;
	    //int voxindex = voxel[0]*dim[1]*dim[2] + voxel[1]*dim[2] + voxel[2];
	    VoxelCoord_t vcoord;
	    vcoord.vx = voxel[0];
	    vcoord.vy = voxel[1];
	    vcoord.vz = voxel[2];
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
	      feat.trackid = xide.trackID;
	    }
	    //tpcinfo._charge_v.data[voxindex] = xide.energy;
	  }
	  
	}//end of loop over ide
      }//end of loop over IDE map
    }//end of loop over simch

    std::cout << "Number of bad IDE: " << nbadide << std::endl;


    return;
  }

  void SimChannelVoxelizer::_make_tensors()
  {

    // Now populate the TPC info arrays
    for ( auto& tpcinfo : _tpcdata_v ) {
     
      int nvoxels = tpcinfo._voxfeat_v.size();
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
      
      std::cout << "allocated space for voxels. fill the arrays" << std::endl;

      size_t index = 0;
      for ( auto it=tpcinfo._voxcoord_2_index.begin(); it!=tpcinfo._voxcoord_2_index.end(); it++ ) {
	tpcinfo._coordindex_v.data[3*index + 0] = it->first.vx;
	tpcinfo._coordindex_v.data[3*index + 1] = it->first.vy;
	tpcinfo._coordindex_v.data[3*index + 2] = it->first.vz;
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
    // Will need to apply t0 shift
    // How to do this for IDEs without labels?
    // Have to hope that ancestor ID in IDEs is correctly filled.
    
    MCPixelPGraph mcpg;
    mcpg.buildgraphonly( mcshower, mctrack, mctruth );

    for ( auto& tpcinfo : _tpcdata_v ) {
      for ( auto& vox : tpcinfo._voxfeat_v ) {
	auto pnode = mcpg.findTrackID( vox.trackid );
	if ( pnode==NULL ) {
	  vox.pdg = 11;
	  vox.ancestorid = 0;
	}
	else {
	  vox.pdg = pnode->pid;
	  vox.ancestorid = pnode->aid;
	}
      }
    }
  }
}
}
