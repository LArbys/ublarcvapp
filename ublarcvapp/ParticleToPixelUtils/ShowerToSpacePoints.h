#ifndef __UBLARCVAPP_PIXELUTILS_SHOWERTOSPACEPOINTS_H__
#define __UBLARCVAPP_PIXELUTILS_SHOWERTOSPACEPOINTS_H__

#include <vector>
#include <set>
#include <map>

#include "larcv/core/DataFormat/Image2D.h"
#include "larlite/DataFormat/larflowcluster.h"
#include "larlite/DataFormat/larflow3dhit.h"

#include "TVector3.h"

namespace ublarcvapp {
namespace pixelutils {

  /**
   * @class ShowerToSpacePoints
   * @brief Convert larlite::larflowcluster (shower representation) to space points with charge
   * 
   * This class takes a larflowcluster which contains a collection of larflow3dhit objects
   * (3D space points), projects them into the 2D wire plane images, and collects the 
   * charge from neighboring pixels.
   */
  class ShowerToSpacePoints {
  public:
    
    /**
     * @struct SpacePointCharge
     * @brief Container for 3D position with associated charge information
     */
    struct SpacePointCharge {
      TVector3 position;      ///< 3D position in detector coordinates (cm)
      float charge;           ///< Total charge from neighboring pixels
      std::vector<float> plane_charges; ///< Charge on each plane [U,V,Y]
      
      SpacePointCharge() : charge(0) { plane_charges.resize(3, 0); }
    };
    
    /**
     * @struct PixelData  
     * @brief Container for pixel location and charge
     */
    struct PixelData {
      int plane;
      int row;
      int col;
      float value;
      
      PixelData() : plane(0), row(0), col(0), value(0) {} // Default constructor for ROOT
      
      PixelData(int p, int r, int c, float v) 
        : plane(p), row(r), col(c), value(v) {}
        
      // For std::set uniqueness
      bool operator<(const PixelData& other) const {
        if (plane != other.plane) return plane < other.plane;
        if (row != other.row) return row < other.row;
        return col < other.col;
      }
    };
    
    ShowerToSpacePoints();
    ~ShowerToSpacePoints();
    
    /**
     * @brief Convert a larflowcluster (shower) to space points with charge
     * 
     * @param cluster The larflowcluster containing larflow3dhit objects
     * @param adc_v Vector of ADC images for each plane
     * @param threshold ADC threshold for pixel selection
     * @param dcol Column window around projected point
     * @param drow Row window around projected point
     * @return Vector of SpacePointCharge objects
     */
    std::vector<SpacePointCharge> convertShower(
      const larlite::larflowcluster& cluster,
      const std::vector<larcv::Image2D>& adc_v,
      const float threshold = 10.0,
      const int dcol = 3,
      const int drow = 3);
    
    /**
     * @brief Get the number of pixels processed in the last conversion
     */
    int getNumPixelsProcessed() const { return _npix_processed; }
    
    /**
     * @brief Get the total charge collected in the last conversion
     */
    float getTotalCharge() const { return _total_charge; }
    
    /**
     * @brief Get the pixel data collected for each space point
     * @return Map from space point index to set of pixels
     */
    const std::map<int, std::set<PixelData>>& getPixelMap() const { return _pixel_map; }
    
    /**
     * @brief Set whether to use charge-weighted position calculation
     */
    void setUseChargeWeighting(bool use) { _use_charge_weighting = use; }
    
  protected:
    
    /**
     * @brief Calculate charge-weighted position from pixel data
     */
    TVector3 getChargeWeightedPosition(const std::set<PixelData>& pixels,
                                      const std::vector<larcv::Image2D>& adc_v,
                                      const TVector3& original_position);
    
    bool _use_charge_weighting;
    int _npix_processed;
    float _total_charge;
    std::map<int, std::set<PixelData>> _pixel_map; ///< Map from hit index to pixels
  };

}
}

#endif