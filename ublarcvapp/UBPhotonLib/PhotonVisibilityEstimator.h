#ifndef __UBLARCVAPP_UBPHOTONLIB_PHOTONVISIBILITYESTIMATOR_H__
#define __UBLARCVAPP_UBPHOTONLIB_PHOTONVISIBILITYESTIMATOR_H__

#include <vector>
#include <map>
#include <cstddef>

namespace ublarcvapp {
namespace ubphotonlib {

  /**
   * @brief Class to estimate photon detection from multiple 3D point sources
   * 
   * This class uses the UBPhotonLib to calculate how many photons from
   * multiple 3D point sources will be detected by each optical detector.
   */
  class PhotonVisibilityEstimator {
  public:
    
    /**
     * @brief Structure to hold a 3D point source with photon count
     */
    struct PhotonSource {
      std::vector<float> position; ///< 3D position in detector coordinates (x,y,z) in cm
      float num_photons;           ///< Number of photons emitted at this location
      
      PhotonSource() : position(3,0), num_photons(0) {}
      PhotonSource(float x, float y, float z, float nphotons) 
        : position{x,y,z}, num_photons(nphotons) {}
    };
    
    /**
     * @brief Constructor
     */
    PhotonVisibilityEstimator();
    
    /**
     * @brief Destructor
     */
    ~PhotonVisibilityEstimator();
    
    /**
     * @brief Clear all point sources
     */
    void clear();
    
    /**
     * @brief Add a single point source
     * @param x X position in cm (detector coordinates)
     * @param y Y position in cm (detector coordinates)
     * @param z Z position in cm (detector coordinates)
     * @param num_photons Number of photons emitted at this location
     */
    void addPhotonSource(float x, float y, float z, float num_photons);
    
    /**
     * @brief Add a single point source
     * @param source PhotonSource struct containing position and photon count
     */
    void addPhotonSource(const PhotonSource& source);
    
    /**
     * @brief Add multiple point sources
     * @param sources Vector of PhotonSource structs
     */
    void addPhotonSources(const std::vector<PhotonSource>& sources);
    
    /**
     * @brief Calculate photons detected by each optical detector
     * @param use_trilinear If true, use trilinear interpolation for visibility (default: true)
     * @return Map from optical channel ID to number of detected photons
     */
    std::map<int, float> calculateDetectedPhotons(bool use_trilinear = true);
    
    /**
     * @brief Get total number of photons detected across all optical detectors
     * @param use_trilinear If true, use trilinear interpolation for visibility (default: true)
     * @return Total number of detected photons
     */
    float getTotalDetectedPhotons(bool use_trilinear = true);
    
    /**
     * @brief Get number of point sources
     * @return Number of point sources currently stored
     */
    size_t getNumSources() const { return _sources.size(); }
    
    /**
     * @brief Get total number of emitted photons from all sources
     * @return Total number of emitted photons
     */
    float getTotalEmittedPhotons() const;
    
    /**
     * @brief Get the collection efficiency (detected/emitted)
     * @param use_trilinear If true, use trilinear interpolation for visibility (default: true)
     * @return Collection efficiency as a fraction
     */
    float getCollectionEfficiency(bool use_trilinear = true);
    
    /**
     * @brief Get number of optical channels
     * @return Number of optical channels in the detector
     */
    int getNumOpticalChannels() const { return _num_optical_channels; }
    
  private:
    std::vector<PhotonSource> _sources;  ///< Vector of photon point sources
    int _num_optical_channels;            ///< Number of optical detectors
  };

}
}

#endif
