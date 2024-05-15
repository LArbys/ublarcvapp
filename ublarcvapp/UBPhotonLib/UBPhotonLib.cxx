#include "UBPhotonLib.h"

#include <iostream>

#include "TFile.h"
#include "TTree.h"

namespace ublarcvapp {
namespace ubphotonlib {

  UBPhotonLib* UBPhotonLib::_ginstance = nullptr;

  UBPhotonLib::UBPhotonLib()
    : _tpc_global_origin_cm{-1.55, -115.53 + 0.5*(117.47+115.53), 0.1},
      _cryo_global_origin_cm{-63.435,-191.61,-92.375},
      _cryo_length_cm{383.22,383.22,1221.75},
      _nvoxels_dim{75,75,400},
      _voxel_len_cm{0,0,0},
      _cryo_origin_tpc_coord_cm{0,0,0}
  {
    
    std::cout << "[UBPhotonLib] Loaded" << std::endl;
    std::cout << "  NVoxels per dim [" << _nvoxels_dim[0] << "," << _nvoxels_dim[1] << "," << _nvoxels_dim[2] << "]" << std::endl;
    for (int i=0; i<3; i++) {
      _voxel_len_cm[i] = _cryo_length_cm[i]/float(_nvoxels_dim[i]);
      _cryo_origin_tpc_coord_cm[i] = _cryo_global_origin_cm[i] - _tpc_global_origin_cm[i];
    }
    std::cout << "  Voxel Length per dim [" << _voxel_len_cm[0] << "," << _voxel_len_cm[1] << "," << _voxel_len_cm[2] << "]" << std::endl;


    std::string dat_folder = std::getenv("UBLARCVAPP_BASEDIR");
    std::string photon_vis_path = dat_folder + "/ublarcvapp/UBPhotonLib/dat/uboone_photon_library_v6_70kV_EnhancedExtraTPCVis.root";
    std::cout << "  PhotonVis File: " << photon_vis_path << std::endl;
    TFile rootfile( photon_vis_path.c_str(), "open");
    TTree* vistree = (TTree*)rootfile.Get("pmtresponse/PhotonLibraryData");
    size_t nentries = vistree->GetEntries();
    std::cout << "  Entries in Vis File: " << nentries << std::endl;
    Int_t Voxel;
    Int_t OpChannel;
    Float_t Visibility;
    vistree->SetBranchAddress("Voxel",&Voxel);
    vistree->SetBranchAddress("OpChannel",&OpChannel);
    vistree->SetBranchAddress("Visibility",&Visibility);
    
    for (size_t i=0; i<nentries; i++) {
      vistree->GetEntry(i);
      std::pair<long,int> voxelkey( (long)Voxel, (int)OpChannel );
      _voxelopchindex_to_visibility[ voxelkey ] = Visibility;
    }

    std::cout << "  Loaded Visibility Info" << std::endl;
  }

  UBPhotonLib::~UBPhotonLib()
  {
  }
      
  UBPhotonLib* UBPhotonLib::getPhotonLib()
  {
    if ( !_ginstance ) {
      _ginstance = new UBPhotonLib();
    }
    return _ginstance;
  }

  void UBPhotonLib::destroy()
  {
    std::cout<< "[UBPhotonLib::destroy()] free the memory used by the photon library" << std::endl;
    delete _ginstance;
    _ginstance = nullptr;
  }

  std::vector<int> UBPhotonLib::getVoxelDimIndices( const std::vector<float>& pos )
  {
    std::vector<int> iindex(3,0);
    for (int i=0; i<3; i++) {
      iindex[i] = (int)((pos[i] - _cryo_origin_tpc_coord_cm[i])/_voxel_len_cm[i]);
    }
    return iindex;
  }

  long UBPhotonLib::getVoxelIndex( const std::vector<float>& pos )
  {
    std::vector<int> iindex = getVoxelDimIndices( pos );
    bool valid = true;
    for (int i=0; i<3; i++) {
      if ( iindex[i]<0 || iindex[i]>=_nvoxels_dim[i] )
	valid = false;
    }

    if (!valid)
      return -1;
    
    long lindex = iindex[0]*_nvoxels_dim[1]*_nvoxels_dim[2] + iindex[1]*_nvoxels_dim[2] + iindex[2];
    
    return lindex;
  }

  float UBPhotonLib::getVisibility( const std::vector<float>& pos, int opch )
  {
    long lindex = getVoxelIndex( pos );
    if ( lindex<0 )
      return 0.0;

    std::pair<long,int> voxopch_key( lindex, opch );
    auto it=_voxelopchindex_to_visibility.find( voxopch_key );
    if ( it==_voxelopchindex_to_visibility.end() ) {
      // not found
      return 0.0;
    }

    // return visibility
    return it->second;
  }

  float UBPhotonLib::getVisibilityTrilinear( const std::vector<float>& pos, int opch )
  {
    long lindex = getVoxelIndex( pos );
    if ( lindex<0 )
      return 0.0;

    std::pair<long,int> voxopch_key( lindex, opch );
    auto it=_voxelopchindex_to_visibility.find( voxopch_key );
    if ( it==_voxelopchindex_to_visibility.end() ) {
      // not found
      return 0.0;
    }

    // return visibility
    return it->second;
  }
  
}
}
