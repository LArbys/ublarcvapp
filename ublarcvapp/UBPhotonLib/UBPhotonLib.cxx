#include "UBPhotonLib.h"

#include <iostream>
#include <chrono>

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
      _cryo_origin_tpc_coord_cm{0,0,0},
      _nopchs(32),
      _nvoxels(0),
      _nvis(0)
  {
    
    std::cout << "[UBPhotonLib] Loaded" << std::endl;
    std::cout << "  NVoxels per dim [" << _nvoxels_dim[0] << "," << _nvoxels_dim[1] << "," << _nvoxels_dim[2] << "]" << std::endl;
    for (int i=0; i<3; i++) {
      _voxel_len_cm[i] = _cryo_length_cm[i]/float(_nvoxels_dim[i]);
      _cryo_origin_tpc_coord_cm[i] = _cryo_global_origin_cm[i] - _tpc_global_origin_cm[i];
    }
    std::cout << "  Voxel Length per dim [" << _voxel_len_cm[0] << "," << _voxel_len_cm[1] << "," << _voxel_len_cm[2] << "]" << std::endl;
    std::cout << "  Cryostat voxel grid origin in TPC coords: ("
	      << _cryo_origin_tpc_coord_cm[0] << ", "
	      << _cryo_origin_tpc_coord_cm[1] << ", " 
	      << _cryo_origin_tpc_coord_cm[2] << ")"
	      << std::endl;
    std::cout << "  Cryostat Max Coordinates: ("
	      << _cryo_origin_tpc_coord_cm[0]+_cryo_length_cm[0] << ", "
	      << _cryo_origin_tpc_coord_cm[1]+_cryo_length_cm[1] << ", "
	      << _cryo_origin_tpc_coord_cm[2]+_cryo_length_cm[2] << ")"
	      << std::endl;

    

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

    auto t1 = std::chrono::high_resolution_clock::now();
    
    //_voxelopchindex_v.resize(nentries);
    std::cout << "   Allocating memory for " << _nvoxels_dim[0] << "x" << _nvoxels_dim[1] << "x" << _nvoxels_dim[2] << "x" << _nopchs << std::endl;
    _nvoxels = _nvoxels_dim[0]*_nvoxels_dim[1]*_nvoxels_dim[2];
    _nvis = _nvoxels*_nopchs;
    _visibility_v.resize( _nvis, 0.0 );

    std::vector<float> vcoord(3,0);    
    long lindex = 0;
    
    for (size_t i=0; i<nentries; i++) {
      vistree->GetEntry(i);
      //std::pair<long,int> voxelkey( (long)Voxel, (int)OpChannel );
      //_voxelopchindex_to_visibility[ voxelkey ] = Visibility;
      //_voxelopchindex_v[i].voxelid = (long)Voxel;
      //_voxelopchindex_v[i].opchid  = (int)OpChannel;
      //_voxelopchindex_v[i].vis     = Visibility;
      //_voxel_to_vindex_v[ (long)Voxel ] = (long)i;
      // need to slice the Voxel

      // copied from getvoxelcoords to maybe trigger vectorization
      //std::vector<long> voxcoord = getVoxelCoords( Voxel );      
      // vcoord[0] = voxelid % _nvoxels_dim[0];
      // long ivox = voxelid-(long)vcoord[0];    
      // vcoord[1] = ((ivox)/_nvoxels_dim[0]) % _nvoxels_dim[1];
      // vcoord[2] = ((ivox - vcoord[1]*_nvoxels_dim[0])/(_nvoxels_dim[1]*_nvoxels_dim[0])) % _nvoxels_dim[2];
      
      //long lindex = getVisLibIndex( Voxel, OpChannel );
      // opchannel is almost certaintly what UB calls "OpDet"      
      lindex = Voxel*_nopchs + OpChannel; 
      
      if ( i>0 && i%10000000==0 ) {
	std::cout << "   loading entry " << i  << ":"
		  << " vis=" << Visibility
		  << " vox=" << Voxel
		  << " opdet=" << OpChannel
		  <<"  lindex=" << lindex
		  << std::endl;
      }
      _visibility_v[lindex] = Visibility;
    }
    
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    double sec = ms_double.count()*0.001;
    std::cout << "  Loaded Visibility Info in " << sec << " secs"  << std::endl;
    
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

  std::vector<long> UBPhotonLib::getVoxelCoords( const std::vector<float>& pos ) const
  {
    std::vector<long> vcoord(3,0);
    for (int i=0; i<3; i++) {
      vcoord[i] = (long) (pos[i] - _cryo_origin_tpc_coord_cm[i])/_voxel_len_cm[i];
    }
    return vcoord;
  }
    
  std::vector<long> UBPhotonLib::getVoxelCoords( long voxelid ) const
  {
    // taken from larsim/Simulation/PhotonVoxels.cxx
    // std::vector<int> ReturnVector;
    // ReturnVector.resize(3);
    // ReturnVector.at(0) =  ID % fxSteps ;
    // ReturnVector.at(1) =  ((ID - ReturnVector.at(0) ) / fxSteps) % fySteps ;
    // ReturnVector.at(2) =  ((ID - ReturnVector.at(0) - (ReturnVector.at(1) * fxSteps)) / (fySteps * fxSteps)) % fzSteps ;
    // return ReturnVector;

    std::vector<long> vcoord(3,0);

    vcoord[0] = voxelid % _nvoxels_dim[0];
    long ivox = voxelid-(long)vcoord[0];    
    vcoord[1] = ((ivox)/_nvoxels_dim[0]) % _nvoxels_dim[1];
    vcoord[2] = ((ivox - vcoord[1]*_nvoxels_dim[0])/(_nvoxels_dim[1]*_nvoxels_dim[0])) % _nvoxels_dim[2];
    return vcoord;
  }

  long UBPhotonLib::getVisLibIndex( long voxelid, int opch ) const
  {
    return voxelid*_nopchs + opch;
  }

  long UBPhotonLib::getVoxelID( const std::vector<long>& voxcoords ) const
  {
    bool valid = true;
    for (int i=0; i<3; i++) {
      if ( voxcoords[i]<0 || voxcoords[i]>=_nvoxels_dim[i] )
	valid = false;
    }
    if (!valid)
      return -1;
    return voxcoords[0] + voxcoords[1]*_nvoxels_dim[0] + voxcoords[2]*_nvoxels_dim[1]*_nvoxels_dim[0];
  }
  
  long UBPhotonLib::getVisLibIndex( const std::vector<float>& pos, int opch ) const
  {
    std::vector<long> iindex = getVoxelCoords( pos );
    long voxelid = getVoxelID(iindex);
    if (voxelid<0)
      return -1;
    
    long lindex = getVisLibIndex( voxelid, opch );
    
    return lindex;
  }

  float UBPhotonLib::getVisibility( const std::vector<float>& pos, int opch )
  {
    long lindex = getVisLibIndex( pos, opch );
    if ( lindex<0 || lindex>(long)_visibility_v.size())
      return 0.0;

    // return visibility
    return _visibility_v[lindex];
  }

  float UBPhotonLib::getVisibilityTrilinear( const std::vector<float>& pos, int opch )
  {
    long lindex = getVoxelIndex( pos, opch );
    if ( lindex<0 || lindex>(long)_visibility_v.size() )
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
