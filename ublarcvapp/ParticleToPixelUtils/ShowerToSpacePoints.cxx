#include "ShowerToSpacePoints.h"
#include "larlite/LArUtil/LArProperties.h"
#include "larlite/LArUtil/Geometry.h"
#include <iostream>

namespace ublarcvapp {
namespace pixelutils {

  ShowerToSpacePoints::ShowerToSpacePoints() 
    : _use_charge_weighting(false),
      _npix_processed(0),
      _total_charge(0.0)
  {
  }
  
  ShowerToSpacePoints::~ShowerToSpacePoints()
  {
  }
  
  std::vector<ShowerToSpacePoints::SpacePointCharge> 
  ShowerToSpacePoints::convertShower(
    const larlite::larflowcluster& cluster,
    const std::vector<larcv::Image2D>& adc_v,
    const float threshold,
    const int dcol,
    const int drow)
  {
    std::vector<SpacePointCharge> spacepoints;
    
    // Reset counters
    _npix_processed = 0;
    _total_charge = 0.0;
    _pixel_map.clear();
    
    // Check we have 3 planes
    if (adc_v.size() != 3) {
      std::cerr << "ShowerToSpacePoints: Expected 3 planes, got " 
                << adc_v.size() << std::endl;
      return spacepoints;
    }
    
    // Get detector parameters
    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float usec_per_tick = 0.5;
    float max_tick = adc_v.front().meta().max_y();
    float min_tick = adc_v.front().meta().min_y();
    
    // Keep track of processed pixels to avoid double counting
    std::set<std::tuple<int,int,int>> processed_pixels; // (plane,row,col)
    
    // Process each larflow3dhit in the cluster
    for (size_t ihit = 0; ihit < cluster.size(); ihit++) {
      
      const larlite::larflow3dhit& hit = cluster[ihit];
      
      // Get 3D position from the larflow3dhit
      // larflow3dhit inherits from std::vector<float> which stores [x,y,z]
      if (hit.size() < 3) {
        std::cerr << "ShowerToSpacePoints: larflow3dhit " << ihit 
                  << " has insufficient position data" << std::endl;
        continue;
      }
      
      TVector3 pos(hit[0], hit[1], hit[2]); // x,y,z stored in vector
      
      // Calculate tick from X position  
      float tick = hit[0]/driftv/usec_per_tick + 3200;
      
      if (tick < min_tick || tick > max_tick)
        continue;
        
      int row = adc_v.front().meta().row(tick);
      
      // Collect charge from each wire plane
      _pixel_map[ihit] = std::set<PixelData>();
      
      for (int p = 0; p < 3; p++) {
        int wire = larutil::Geometry::GetME()->WireCoordinate(pos, p);
        
        // Collect charge from neighboring pixels
        for (int dr = -abs(drow); dr <= abs(drow); dr++) {
          int r = row + dr;
          if (r < 0 || r >= (int)adc_v[p].meta().rows())
            continue;
            
          for (int dc = -abs(dcol); dc <= abs(dcol); dc++) {
            int c = wire + dc;
            if (c < 0 || c >= (int)adc_v[p].meta().cols())
              continue;
              
            float pixval = adc_v[p].pixel(r, c);
            if (pixval > threshold) {
              // Check if we've already processed this pixel
              auto pixel_key = std::make_tuple(p, r, c);
              if (processed_pixels.find(pixel_key) == processed_pixels.end()) {
                processed_pixels.insert(pixel_key);
                _pixel_map[ihit].insert(PixelData(p, r, c, pixval));
                _npix_processed++;
                _total_charge += pixval;
              }
            }
          }
        }
      }
    }
    
    // Now create space points from the collected pixel data
    for (size_t ihit = 0; ihit < cluster.size(); ihit++) {
      
      // Skip hits with no associated pixels
      if (_pixel_map.find(ihit) == _pixel_map.end() || _pixel_map[ihit].empty())
        continue;
        
      const larlite::larflow3dhit& hit = cluster[ihit];
      SpacePointCharge sp;
      
      // Use the original 3D position or charge-weighted position
      // Get position from the std::vector<float> base class
      TVector3 original_pos(hit[0], hit[1], hit[2]);
      
      if (_use_charge_weighting && !_pixel_map[ihit].empty()) {
        sp.position = getChargeWeightedPosition(_pixel_map[ihit], adc_v, original_pos);
      } else {
        sp.position = original_pos;
      }
      
      // Calculate total charge and per-plane charges
      sp.charge = 0;
      sp.plane_charges.resize(3, 0);
      
      for (const auto& pix : _pixel_map[ihit]) {
        sp.charge += pix.value;
        sp.plane_charges[pix.plane] += pix.value;
      }
      
      spacepoints.push_back(sp);
    }
    
    return spacepoints;
  }
  
  TVector3 ShowerToSpacePoints::getChargeWeightedPosition(
    const std::set<PixelData>& pixels,
    const std::vector<larcv::Image2D>& adc_v,
    const TVector3& original_position)
  {
    // For now, just return the original position
    // In the future, could implement charge-weighted averaging
    // across the pixels to refine the position
    return original_position;
  }

}
}