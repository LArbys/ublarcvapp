#ifndef __UBLARCVAPP_UBPHOTONLIB_UBPHOTONLIB_H__
#define __UBLARCVAPP_UBPHOTONLIB_UBPHOTONLIB_H__

#include <map>
#include <vector>

namespace ublarcvapp {
namespace ubphotonlib {

  /** 
   * @brief Singleton class providing interface to the UB Photon Library
   *
   *
   */
  class UBPhotonLib {

  public:

    static UBPhotonLib* getPhotonLib();
    static void destroy(); ///< to free up memory

    float getVisibility( const std::vector<float>& pos, int opch );
    float getVisibilityTrilinear( const std::vector<float>& pos, int opch );

    std::vector<long> getVoxelCoords( long voxelid ) const;
    std::vector<long> getVoxelCoords( const std::vector<float>& pos ) const;    
    long getVisLibIndex( long voxelid, int opch ) const;
    long getVoxelIndex( const std::vector<float>& pos, int opch ) const;
    long getVisLibIndex( const std::vector<float>& pos, int opch ) const;
    long getVoxelID( const std::vector<long>& voxcoords ) const;

    struct VisData_t {
      long voxelid;
      long opchid;
      int dim_index[3];
      float vis;
      VisData_t()
	: voxelid(0),
	  opchid(0),
	  dim_index{0,0,0},
	  vis(0.0)
      {};
      bool operator<( VisData_t& rhs ) const {
	if (rhs.voxelid < voxelid )
	  return true;
	if ( rhs.voxelid==voxelid && rhs.opchid < opchid )
	  return true;
	return false;
      };
    };
    
  protected:

    UBPhotonLib();
    virtual ~UBPhotonLib();

    
    
    float _tpc_global_origin_cm[3];
    float _cryo_global_origin_cm[3];
    float _cryo_length_cm[3];    
    int _nvoxels_dim[3];
    float _voxel_len_cm[3];
    float _cryo_origin_tpc_coord_cm[3];
    int _nopchs;
    long _nvoxels;
    long _nvis;

    std::map< std::pair<long,int> , float > _voxelopchindex_to_visibility;
    std::vector< VisData_t > _voxelopchindex_v;
    std::map< long, long > _voxel_to_vindex_v;

    // fastest look up is a 4D matrix    
    // we'll use a 1D vector and execute strides
    // will be a lot of zeros. But it's only double the memory.
    std::vector<float> _visibility_v;
    
  private:

    static UBPhotonLib* _ginstance;

  };

  
}
}

#endif
