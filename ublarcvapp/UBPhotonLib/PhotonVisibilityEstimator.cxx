#include "PhotonVisibilityEstimator.h"
#include "UBPhotonLib.h"
#include <iostream>

namespace ublarcvapp {
namespace ubphotonlib {

  PhotonVisibilityEstimator::PhotonVisibilityEstimator() 
    : _num_optical_channels(32)  // MicroBooNE has 32 PMTs
  {
    // Ensure photon library is loaded
    UBPhotonLib::getPhotonLib();
  }
  
  PhotonVisibilityEstimator::~PhotonVisibilityEstimator() 
  {
    // Nothing to clean up
  }
  
  void PhotonVisibilityEstimator::clear() 
  {
    _sources.clear();
  }
  
  void PhotonVisibilityEstimator::addPhotonSource(float x, float y, float z, float num_photons) 
  {
    _sources.emplace_back(x, y, z, num_photons);
  }
  
  void PhotonVisibilityEstimator::addPhotonSource(const PhotonSource& source) 
  {
    _sources.push_back(source);
  }
  
  void PhotonVisibilityEstimator::addPhotonSources(const std::vector<PhotonSource>& sources) 
  {
    _sources.insert(_sources.end(), sources.begin(), sources.end());
  }
  
  std::map<int, float> PhotonVisibilityEstimator::calculateDetectedPhotons(bool use_trilinear) 
  {
    std::map<int, float> photons_per_opdet;
    
    // Initialize map with zeros for all optical channels
    for (int opch = 0; opch < _num_optical_channels; ++opch) {
      photons_per_opdet[opch] = 0.0;
    }
    
    // Get photon library instance
    UBPhotonLib* photon_lib = UBPhotonLib::getPhotonLib();
    
    // For each point source
    for (const auto& source : _sources) {
      // Skip if no photons
      if (source.num_photons <= 0) continue;
      
      // For each optical channel
      for (int opch = 0; opch < _num_optical_channels; ++opch) {
        float visibility = 0.0;
        
        // Get visibility from photon library
        if (use_trilinear) {
          visibility = photon_lib->getVisibilityTrilinear(source.position, opch);
        } else {
          visibility = photon_lib->getVisibility(source.position, opch);
        }
        
        // Calculate number of photons detected by this channel from this source
        // Visibility is the probability that a photon from this voxel reaches this opdet
        float detected_photons = source.num_photons * visibility;
        
        // Add to total for this optical channel
        photons_per_opdet[opch] += detected_photons;
      }
    }
    
    return photons_per_opdet;
  }
  
  float PhotonVisibilityEstimator::getTotalDetectedPhotons(bool use_trilinear) 
  {
    auto photons_per_opdet = calculateDetectedPhotons(use_trilinear);
    
    float total = 0.0;
    for (const auto& entry : photons_per_opdet) {
      total += entry.second;
    }
    
    return total;
  }
  
  float PhotonVisibilityEstimator::getTotalEmittedPhotons() const 
  {
    float total = 0.0;
    for (const auto& source : _sources) {
      total += source.num_photons;
    }
    return total;
  }
  
  float PhotonVisibilityEstimator::getCollectionEfficiency(bool use_trilinear) 
  {
    float emitted = getTotalEmittedPhotons();
    if (emitted <= 0) return 0.0;
    
    float detected = getTotalDetectedPhotons(use_trilinear);
    return detected / emitted;
  }

}
}