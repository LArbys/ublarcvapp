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

    // The coordinates of a voxel
    struct VoxelCoord_t {
      unsigned int vtick; ///< voxel index in the tick dim
      unsigned int vy;    ///< voxel index in the y dim
      unsigned int vz;    ///< voxel index in the z dim
      int tdc;            ///< Discritized Global Simulation Time
      int tick;           ///< Tick of the TPC clock
      std::vector<float> xyz; ///< coordinate associated with this voxel
      
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
    
    // Features associated to a voxel
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
      std::vector< float > _voxel_dim_v;
      std::vector< int >   _num_voxels_v;

      larcv::NumpyArrayInt   _coordindex_v;
      larcv::NumpyArrayFloat _coordpos_v;
      larcv::NumpyArrayFloat _charge_v;
      larcv::NumpyArrayInt   _pdg_v;
      larcv::NumpyArrayInt   _trackid_v;
      larcv::NumpyArrayInt   _ancestorid_v;
      std::vector< larcv::NumpyArrayFloat > _simch_img_v;

      std::map< VoxelCoord_t, unsigned long > _voxcoord_2_index;
      std::vector< VoxelFeat_t > _voxfeat_v;
      std::map< VoxelCoord_t, unsigned long > _voxcoord_2_sparsearrayindex;

      void clear_event_data()
      {
	_coordindex_v.data.clear();
	_coordpos_v.data.clear();
	_charge_v.data.clear();
	_charge_v.data.clear();
	_trackid_v.data.clear();
	_ancestorid_v.data.clear();

	_simch_img_v.clear();
	
	_voxcoord_2_index.clear();
	_voxfeat_v.clear();
	_voxcoord_2_sparsearrayindex.clear();
		
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
    SimChannelVoxelizer( const std::vector<float>& voxel_dims );
    virtual ~SimChannelVoxelizer() {};

    void set_simch_treename( std::string treename ) { _simch_tree_name=treename; };

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

    VoxelCoord_t makeVoxelCoord( const float tick, const float y, const float z, const int tdc,
				 TPCInfo& tpcinfo );
    VoxelCoord_t makeVoxelCoord( const float tick, const float y, const float z, const int tdc,
				 const int cryoid, const int tpcid );
    void defineTPCVoxels( const std::vector<float>& voxel_dims );
    void makeTPCVoxels( const int cryoid, const int tpcid, const std::vector<float>& voxel_dims );
    int getReadoutNumSamples( larlite::geo::DetId_t geoid );
    int getReadoutTriggerOffset( larlite::geo::DetId_t geoid );
    int getReadoutFirstTick( larlite::geo::DetId_t geoid );

  public:

    std::vector< float > _global_voxel_dim_v;
    std::vector< TPCInfo > _tpcdata_v;
    std::map< std::pair<int,int>, int >   _id2index_v; ///< map from (cryoid,tpcid) to index in _tpcdata_v;
    TPCInfo& getTPCInfo( const int cryoid, const int tpcid );
    bool getVoxelIndex( const float tick, const float y, const float z,
			const int cryoid, const int tpcid,
			std::vector<long>& indices );
    bool getVoxelIndexWithTPCInfo( const float tick, const float y, const float z,
				   TPCInfo& tpcdata, std::vector<long>& indices );
    
    
  protected:

    std::string _simch_tree_name;
    
    void _scan_IDEs( const larlite::event_simch& ev_simch);
    void _make_tensors();
    void _pdg_labels( const larlite::event_mctrack& mctrack,
		      const larlite::event_mcshower& mcshower,
		      const larlite::event_mctruth& mctruth);
    void _fill_simch_images( const larlite::event_simch& ev_simch );
      

    
  };

}
}

#endif
