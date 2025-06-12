#ifndef __UBLARCVAPP_PARTICLETOPIXELUTILS_TRACKTOSPACEPOINTS_H__
#define __UBLARCVAPP_PARTICLETOPIXELUTILS_TRACKTOSPACEPOINTS_H__

#include <vector>
#include <map>
#include <set>

#include "larlite/DataFormat/track.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "TVector3.h"

namespace ublarcvapp {
namespace pixelutils {

  /**
   * @brief Convert larlite tracks to collections of 3D space points with charge
   * 
   * This class takes a larlite::track object and extracts charge information
   * from ADC images at the projected 2D positions along the track.
   * It produces a collection of 3D space points with associated charge values.
   */
  class TrackToSpacePoints {
  public:
    
    /**
     * @brief Structure to hold a 3D space point with charge information
     */
    struct SpacePointCharge {
      TVector3 position;      ///< 3D position in detector coordinates (cm)
      float charge;           ///< Total charge from neighboring pixels
      int wire_u;             ///< U-plane wire number
      int wire_v;             ///< V-plane wire number  
      int wire_y;             ///< Y-plane wire number
      int tick;               ///< Time tick
      std::vector<float> plane_charges; ///< Charge on each plane [U,V,Y]
      
      SpacePointCharge() 
        : position(0,0,0), charge(0), 
          wire_u(-1), wire_v(-1), wire_y(-1), tick(-1),
          plane_charges(3,0) {}
    };
    
    /**
     * @brief Structure to hold pixel coordinate and value
     */
    struct PixelData {
      int plane;              ///< Wire plane index
      int row;                ///< Row (tick) index
      int col;                ///< Column (wire) index  
      float value;            ///< Pixel ADC value
      
      PixelData(int p, int r, int c, float v) 
        : plane(p), row(r), col(c), value(v) {}
        
      // Define comparison operator for set
      bool operator<(const PixelData& other) const {
        if (plane != other.plane) return plane < other.plane;
        if (row != other.row) return row < other.row;
        return col < other.col;
      }
    };
    
    /**
     * @brief Constructor
     */
    TrackToSpacePoints();
    
    /**
     * @brief Destructor
     */
    ~TrackToSpacePoints();
    
    /**
     * @brief Convert a track to space points with charge
     * @param track The larlite::track to convert
     * @param adc_v Vector of ADC images for each wire plane
     * @param threshold ADC threshold for pixel selection
     * @param dcol Column window around projected wire position
     * @param drow Row window around projected tick position
     * @param minstepsize Minimum step size along track (cm)
     * @param maxstepsize Maximum step size along track (cm)
     * @return Vector of SpacePointCharge objects along the track
     */
    std::vector<SpacePointCharge> convertTrack(
      const larlite::track& track,
      const std::vector<larcv::Image2D>& adc_v,
      const float threshold = 10.0,
      const int dcol = 3,
      const int drow = 3,
      const float minstepsize = 0.3,
      const float maxstepsize = 0.5
    );
    
    /**
     * @brief Get all unique pixels touched by the track projection
     * @param track The larlite::track to analyze
     * @param adc_v Vector of ADC images for each wire plane
     * @param threshold ADC threshold for pixel selection
     * @param dcol Column window around projected wire position
     * @param drow Row window around projected tick position
     * @param minstepsize Minimum step size along track (cm)
     * @param maxstepsize Maximum step size along track (cm)
     * @return Set of unique PixelData objects
     */
    std::set<PixelData> getUniquePixels(
      const larlite::track& track,
      const std::vector<larcv::Image2D>& adc_v,
      const float threshold = 10.0,
      const int dcol = 3,
      const int drow = 3,
      const float minstepsize = 0.3,
      const float maxstepsize = 0.5
    );
    
    /**
     * @brief Get charge-weighted average position for a collection of pixels
     * @param pixels Set of pixels with charge values
     * @param adc_v Vector of ADC images for wire/tick to position conversion
     * @return Charge-weighted 3D position
     */
    TVector3 getChargeWeightedPosition(
      const std::set<PixelData>& pixels,
      const std::vector<larcv::Image2D>& adc_v
    );
    
    /**
     * @brief Set whether to use charge weighting for space point positions
     * @param use_weighting If true, use charge-weighted positions
     */
    void setUseChargeWeighting(bool use_weighting) { _use_charge_weighting = use_weighting; }
    
    /**
     * @brief Get total number of pixels processed in last conversion
     */
    int getNumPixelsProcessed() const { return _npix_processed; }
    
    /**
     * @brief Get total charge collected in last conversion
     */
    float getTotalCharge() const { return _total_charge; }
    
  private:
    
    bool _use_charge_weighting;  ///< Use charge-weighted positions
    int _npix_processed;         ///< Number of pixels processed
    float _total_charge;         ///< Total charge collected
    
    /**
     * @brief Convert tick to X position
     * @param tick Time tick
     * @return X position in cm
     */
    float tickToX(int tick) const;
    
    /**
     * @brief Convert wire and tick to Y,Z position
     * @param plane Wire plane index
     * @param wire Wire number
     * @param tick Time tick
     * @param y Output Y position
     * @param z Output Z position
     */
    void wireTickToYZ(int plane, int wire, int tick, float& y, float& z) const;
  };

}
}

#endif