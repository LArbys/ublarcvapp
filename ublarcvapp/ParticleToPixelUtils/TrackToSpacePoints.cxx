#include "TrackToSpacePoints.h"
#include "larlite/LArUtil/LArProperties.h"
#include "larlite/LArUtil/Geometry.h"
#include <iostream>

namespace ublarcvapp {
namespace pixelutils {

  TrackToSpacePoints::TrackToSpacePoints() 
    : _use_charge_weighting(false),
      _npix_processed(0),
      _total_charge(0.0)
  {
  }
  
  TrackToSpacePoints::~TrackToSpacePoints()
  {
  }
  
  std::vector<TrackToSpacePoints::SpacePointCharge> 
  TrackToSpacePoints::convertTrack(
    const larlite::track& track,
    const std::vector<larcv::Image2D>& adc_v,
    const float threshold,
    const int dcol,
    const int drow,
    const float minstepsize,
    const float maxstepsize)
  {
    std::vector<SpacePointCharge> spacepoints;
    
    // Reset counters
    _npix_processed = 0;
    _total_charge = 0.0;
    
    int npts = track.NumberTrajectoryPoints();
    if (npts <= 1) {
      std::cerr << "TrackToSpacePoints: Cannot process track with only " 
                << npts << " point(s)" << std::endl;
      return spacepoints;
    }
    
    // Check we have 3 planes
    if (adc_v.size() != 3) {
      std::cerr << "TrackToSpacePoints: Expected 3 planes, got " 
                << adc_v.size() << std::endl;
      return spacepoints;
    }
    
    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float usec_per_tick = 0.5;
    float max_tick = adc_v.front().meta().max_y();
    float min_tick = adc_v.front().meta().min_y();
    
    // Keep track of processed pixels to avoid double counting
    std::set<std::tuple<int,int,int>> processed_pixels; // (plane,row,col)
    
    // Map from 3D position to pixel data
    std::map<int, std::set<PixelData>> step_pixels; // step index -> pixels
    std::map<int, TVector3> step_positions; // step index -> 3D position
    
    int global_step = 0;
    
    // Loop over track segments
    for (int ipt = 0; ipt < npts-1; ipt++) {
      
      TVector3 start = track.LocationAtPoint(ipt);
      TVector3 end = track.LocationAtPoint(ipt+1);
      TVector3 dir = end - start;
      
      double segsize = dir.Mag();
      
      int nsteps = 1;
      if (segsize > minstepsize) {
        nsteps = segsize/maxstepsize + 1;
      }
      
      float stepsize = segsize/float(nsteps);
      
      // Step along the segment
      for (int istep = 0; istep <= nsteps; istep++) {
        
        // Get 3D position along track
        TVector3 pos = start + istep*(stepsize/segsize)*dir;
        step_positions[global_step] = pos;
        
        // Calculate tick from X position
        int tick = pos[0]/driftv/usec_per_tick + 3200;
        
        if (tick < min_tick || tick > max_tick)
          continue;
          
        int row = adc_v.front().meta().row(tick);
        
        // Project into each wire plane
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
                  step_pixels[global_step].insert(PixelData(p, r, c, pixval));
                  _npix_processed++;
                  _total_charge += pixval;
                }
              }
            }
          }
        }
        
        global_step++;
      }
    }
    
    // Now create space points from the collected pixel data
    for (auto& step_entry : step_pixels) {
      int step_idx = step_entry.first;
      auto& pixels = step_entry.second;
      
      if (pixels.empty())
        continue;
        
      SpacePointCharge sp;
      
      // Use the track position or charge-weighted position
      if (_use_charge_weighting && !pixels.empty()) {
        sp.position = getChargeWeightedPosition(pixels, adc_v);
      } else {
        sp.position = step_positions[step_idx];
      }
      
      // Calculate total charge and per-plane charges
      sp.charge = 0;
      sp.plane_charges.resize(3, 0);
      
      int tick_sum = 0;
      int tick_count = 0;
      std::vector<int> wire_nums(3, -1);
      
      for (const auto& pix : pixels) {
        sp.charge += pix.value;
        sp.plane_charges[pix.plane] += pix.value;
        
        if (pix.plane == 0 && wire_nums[0] == -1) {
          wire_nums[0] = pix.col;
        } else if (pix.plane == 1 && wire_nums[1] == -1) {
          wire_nums[1] = pix.col;
        } else if (pix.plane == 2 && wire_nums[2] == -1) {
          wire_nums[2] = pix.col;
        }
        
        tick_sum += pix.row;
        tick_count++;
      }
      
      // Set wire numbers and average tick
      sp.wire_u = wire_nums[0];
      sp.wire_v = wire_nums[1];
      sp.wire_y = wire_nums[2];
      sp.tick = (tick_count > 0) ? tick_sum / tick_count : -1;
      
      spacepoints.push_back(sp);
    }
    
    return spacepoints;
  }
  
  std::set<TrackToSpacePoints::PixelData> 
  TrackToSpacePoints::getUniquePixels(
    const larlite::track& track,
    const std::vector<larcv::Image2D>& adc_v,
    const float threshold,
    const int dcol,
    const int drow,
    const float minstepsize,
    const float maxstepsize)
  {
    std::set<PixelData> unique_pixels;
    
    int npts = track.NumberTrajectoryPoints();
    if (npts <= 1) {
      return unique_pixels;
    }
    
    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float usec_per_tick = 0.5;
    float max_tick = adc_v.front().meta().max_y();
    float min_tick = adc_v.front().meta().min_y();
    
    // Loop over track segments
    for (int ipt = 0; ipt < npts-1; ipt++) {
      
      TVector3 start = track.LocationAtPoint(ipt);
      TVector3 end = track.LocationAtPoint(ipt+1);
      TVector3 dir = end - start;
      
      double segsize = dir.Mag();
      
      int nsteps = 1;
      if (segsize > minstepsize) {
        nsteps = segsize/maxstepsize + 1;
      }
      
      float stepsize = segsize/float(nsteps);
      
      // Step along the segment
      for (int istep = 0; istep <= nsteps; istep++) {
        
        // Get 3D position along track
        TVector3 pos = start + istep*(stepsize/segsize)*dir;
        
        // Calculate tick from X position
        int tick = pos[0]/driftv/usec_per_tick + 3200;
        
        if (tick < min_tick || tick > max_tick)
          continue;
          
        int row = adc_v.front().meta().row(tick);
        
        // Project into each wire plane
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
                unique_pixels.insert(PixelData(p, r, c, pixval));
              }
            }
          }
        }
      }
    }
    
    return unique_pixels;
  }
  
  TVector3 TrackToSpacePoints::getChargeWeightedPosition(
    const std::set<PixelData>& pixels,
    const std::vector<larcv::Image2D>& adc_v)
  {
    if (pixels.empty()) {
      return TVector3(0, 0, 0);
    }
    
    // Group pixels by plane
    std::map<int, std::vector<const PixelData*>> plane_pixels;
    for (const auto& pix : pixels) {
      plane_pixels[pix.plane].push_back(&pix);
    }
    
    // Calculate charge-weighted wire position for each plane
    std::map<int, float> weighted_wires;
    std::map<int, float> total_charges;
    float weighted_tick = 0;
    float total_tick_charge = 0;
    
    for (const auto& plane_pair : plane_pixels) {
      int plane = plane_pair.first;
      float weighted_wire = 0;
      float plane_charge = 0;
      
      for (const auto* pix : plane_pair.second) {
        weighted_wire += pix->col * pix->value;
        plane_charge += pix->value;
        weighted_tick += pix->row * pix->value;
        total_tick_charge += pix->value;
      }
      
      if (plane_charge > 0) {
        weighted_wires[plane] = weighted_wire / plane_charge;
        total_charges[plane] = plane_charge;
      }
    }
    
    if (total_tick_charge > 0) {
      weighted_tick /= total_tick_charge;
    }
    
    // Convert weighted wire/tick to 3D position
    // This is approximate - ideally would use proper 3D reconstruction
    float x = tickToX(weighted_tick);
    
    // Average Y,Z from different plane projections
    float y_sum = 0, z_sum = 0;
    int n_planes = 0;
    
    for (const auto& wire_pair : weighted_wires) {
      int plane = wire_pair.first;
      float wire = wire_pair.second;
      float y, z;
      wireTickToYZ(plane, wire, weighted_tick, y, z);
      y_sum += y;
      z_sum += z;
      n_planes++;
    }
    
    if (n_planes > 0) {
      y_sum /= n_planes;
      z_sum /= n_planes;
    }
    
    return TVector3(x, y_sum, z_sum);
  }
  
  float TrackToSpacePoints::tickToX(int tick) const
  {
    const float driftv = larutil::LArProperties::GetME()->DriftVelocity();
    const float usec_per_tick = 0.5;
    return (tick - 3200) * driftv * usec_per_tick;
  }
  
  void TrackToSpacePoints::wireTickToYZ(int plane, int wire, int tick, 
                                        float& y, float& z) const
  {
    // Get 3D position from wire/tick
    // This is a simplified version - in practice would use proper geometry
    auto const* geom = larutil::Geometry::GetME();
    
    // Get wire start and end points
    double xyz_start[3], xyz_end[3];
    geom->WireEndPoints(plane, wire, xyz_start, xyz_end);
    
    // For now, use wire midpoint for Y,Z
    // A more sophisticated approach would use the tick info
    y = (xyz_start[1] + xyz_end[1]) / 2.0;
    z = (xyz_start[2] + xyz_end[2]) / 2.0;
  }

}
}