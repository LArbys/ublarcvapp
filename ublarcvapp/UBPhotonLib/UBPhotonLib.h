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

    long getVoxelIndex( const std::vector<float>& pos );
    
  protected:

    UBPhotonLib();
    virtual ~UBPhotonLib();


    float _tpc_global_origin_cm[3];
    float _cryo_global_origin_cm[3];
    float _cryo_length_cm[3];    
    int _nvoxels_dim[3];
    float _voxel_len_cm[3];
    float _cryo_origin_tpc_coord_cm[3];    

    std::map< std::pair<long,int> , float > _voxelopchindex_to_visibility;
    
  private:

    static UBPhotonLib* _ginstance;

  };

  
}
}

#endif
