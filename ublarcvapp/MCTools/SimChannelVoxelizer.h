#ifndef __UBLARCVAPP_MCTOOLS_SIMCHANNEL_VOXELIZER_H__
#define __UBLARCVAPP_MCTOOLS_SIMCHANNEL_VOXELIZER_H__

#include <vector>
#include <map>

#include "larlite/Base/GeoConstants.h"
#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/LArProperties.h"
#include "larlite/LArUtil/DetectorProperties.h"
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/simch.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mctruth.h"

#include "larcv/core/PyUtil/NumpyArray.h"


namespace ublarcvapp {
namespace mctools {


  /**
   * @brief Read in SimChannel information and produce 3D energy deposition truth
   *
   */
  class SimChannelVoxelizer {
  public:

    struct VoxelCoord_t {
      unsigned int vtick;
      unsigned int vy;
      unsigned int vz;
      int tdc;
      int tick;
      
      bool operator<( const VoxelCoord_t& rhs ) const {
	if ( vtick<rhs.vtick )
	  return true;
	if ( vtick==rhs.vtick && vy<rhs.vy )
	  return true;
	if ( vtick==rhs.vtick && vy==rhs.vy && vz<rhs.vz )
	  return true;
	return false;
      };
    };
    struct VoxelFeat_t {
      int pdg;
      long trackid;
      long ancestorid;
      float charge;
      std::array<float,4> realpos;

      VoxelFeat_t()
	: pdg(0),
	  trackid(0),
	  ancestorid(0),
	  charge(0),
	  realpos({0,0,0,0})
      {};
    };      
    
    struct TPCInfo {

      int cryoid;
      int tpcid;
      float driftdir;
      float realx_anode;
      float vx_trigger;
      
      std::vector< float > _origin_cm_v;
      std::vector< float > _voxel_dim_cm_v;
      std::vector< int >   _num_voxels_v;

      larcv::NumpyArrayInt   _coordindex_v;
      larcv::NumpyArrayFloat _coordpos_v;
      larcv::NumpyArrayFloat _charge_v;
      larcv::NumpyArrayInt   _pdg_v;
      larcv::NumpyArrayInt   _trackid_v;
      larcv::NumpyArrayInt   _ancestorid_v;

      std::map< VoxelCoord_t, unsigned long > _voxcoord_2_index;
      std::vector< VoxelFeat_t > _voxfeat_v;

      void clear_event_data()
      {
	_coordindex_v.data.clear();
	_coordpos_v.data.clear();
	_charge_v.data.clear();
	_charge_v.data.clear();
	_trackid_v.data.clear();
	_ancestorid_v.data.clear();
	_voxcoord_2_index.clear();
	_voxfeat_v.clear();
	
	_coordindex_v.shape.clear();
	_coordpos_v.shape.clear();
	_charge_v.shape.clear();
	_charge_v.shape.clear();
	_trackid_v.shape.clear();
	_ancestorid_v.shape.clear();
      };
      
    };
    
  public:
    SimChannelVoxelizer();
    virtual ~SimChannelVoxelizer() {};

    void clear();
    void process( larlite::storage_manager& ioll );
    void process( const larlite::event_simch& ev_simch,
		  const larlite::event_mctrack& mctrack,
		  const larlite::event_mcshower& mcshower,
		  const larlite::event_mctruth& mctruth);    
    // void setOrigin( float x_cm, float y_cm, float z_cm );
    // void setOrigin( const std::vector<float>& origin_cm ) { _origin_cm_v=origin_cm; };
    // void setVoxelDimensions( float x_cm, float y_cm, float x_cm );
    // void setNumVoxels( int n_x, int n_y, int n_z );
    void makeDefaultTPCVoxels( int cryoid, int tpcid );
    int getReadoutNumSamples( larlite::geo::DetId_t geoid );
    int getReadoutTriggerOffset( larlite::geo::DetId_t geoid );
    int getReadoutFirstTick( larlite::geo::DetId_t geoid );
    
  public:

    std::vector< TPCInfo > _tpcdata_v;
    std::map< std::pair<int,int>, int >   _id2index_v; ///< map from (cryoid,tpcid) to index in _tpcdata_v;

  protected:

    void _scan_IDEs( const larlite::event_simch& ev_simch);
    void _make_tensors();
    void _pdg_labels( const larlite::event_mctrack& mctrack,
		      const larlite::event_mcshower& mcshower,
		      const larlite::event_mctruth& mctruth);
      

    
  };

}
}

#endif
