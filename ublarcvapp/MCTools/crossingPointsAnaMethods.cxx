#include "crossingPointsAnaMethods.h"

#include <cstring>
#include <stdexcept>

#include "TTree.h"
#include "TLorentzVector.h"

#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"
#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/LArProperties.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/trigger.h"

// ublarcvapp
#include "ublarcvapp/UBWireTool/UBWireTool.h"
#include "ublarcvapp/ubdllee/dwall.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * get tick from mcstep
   *
   *
   */
  float CrossingPointsAnaMethods::getTick( const larlite::mcstep& step, const float trig_time,
                                           const larutil::SpaceChargeMicroBooNE* psce ) {
    std::vector<float> pos(4,0);
    pos[0] = step.T();
    pos[1] = step.X();
    pos[2] = step.Y();
    pos[3] = step.Z();
    return getTick( pos, trig_time, psce );
  }
  

  /**
   * get tick from mcstep
   *
   *
   */  
  float CrossingPointsAnaMethods::getTick( const std::vector<float>& step, const float trig_time,
                                           const larutil::SpaceChargeMicroBooNE* psce ) {    
    // Function returns the tick time of a MC step point
    // if SCE pointer is null, we do not correct for the space charge
    
    std::vector<double> dpos(3,0);
    if ( psce ) {
      std::vector<double> pos_offset = psce->GetPosOffsets( step[1], step[2], step[3] );
      dpos[0] = step[1] - pos_offset[0] + 0.7;
    }
    else {
      dpos[0] = step[1];
    }
    
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;    
    float tick = ( step[0]*1.0e-3 - (trig_time-4050.0) )/0.5 + dpos[0]/cm_per_tick + 3200.0;
    
    return tick;
  }

  /**
   * @brief get TPC tick from MC step info. do not apply drift time.
   *
   */  
  float CrossingPointsAnaMethods::getTrueTick( const std::vector<float>& step, const float trig_time,
					       const larutil::SpaceChargeMicroBooNE* psce ) {    
    // Function returns the tick time of a MC step point
    // if SCE pointer is null, we do not correct for the space charge
    
    std::vector<double> dpos(3,0);
    if ( psce ) {
      std::vector<double> pos_offset = psce->GetPosOffsets( step[1], step[2], step[3] );
      dpos[0] = step[1] - pos_offset[0] + 0.7;
    }
    else {
      dpos[0] = step[1];
    }
    
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;    
    float tick = ( step[0]*1.0e-3 - (trig_time-4050.0) )/0.5 + 3200.0;
    
    return tick;
  }
  
  /**
   * doesTrackCrossImageBoundary
   *
   *
   */  
  int CrossingPointsAnaMethods::doesTrackCrossImageBoundary( const larlite::mctrack& track, const larcv::ImageMeta& meta,
                                                             const float trig_time, const larutil::SpaceChargeMicroBooNE* psce ) {
    if ( track.size()==0 )
      return -1;
    
    float tick_start = getTick( track.front(), trig_time, psce );
    float tick_end   = getTick( track.back(), trig_time, psce );
    if ( tick_start>meta.min_y() && tick_start<meta.max_y() && tick_end>meta.min_y() && tick_end<meta.max_y() )
      return 2;

    if ( tick_start < meta.min_y() && tick_end > meta.min_y() )
      return 0; // start out -> end in
    else if ( tick_start > meta.min_y() && tick_end < meta.min_y())
      return 1; // start in -> end out
    else if ( tick_start < meta.max_y() && tick_end > meta.max_y())
      return 1; // start in -> end out;
    else if ( tick_start > meta.max_y() && tick_end < meta.max_y() )
      return 0; // start out -> end in

    return -1;
  }

  /**
   * get first visible point on track
   *
   */
  std::vector<int>
  CrossingPointsAnaMethods::getFirstStepPosInsideImage( const larlite::mctrack& track, const larcv::ImageMeta& meta,
                                                        const float trig_time,
                                                        const bool startAtstart,
                                                        const float max_step_size, const float fv_border,
                                                        std::vector<float>& endpt3d,
                                                        const larutil::SpaceChargeMicroBooNE* psce,
                                                        bool verbose ) {
    
    // This function returns the (SCE-corrected) position where a MC track first is inside the image bounds
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;    
    int npts = (int)track.size();
    endpt3d.clear();

    for ( int ipt=1; ipt<npts; ipt++ ) {

      int thispt = ipt;
      int lastpt = thispt-1;
      if ( !startAtstart ) {
	thispt = npts-1-ipt;
	lastpt = thispt+1;
      }
      
      const auto& this_step = track.at( thispt );
      const auto& last_step = track.at( lastpt );
      
      float dir[3] = { float(this_step.X()-last_step.X()), float(this_step.Y()-last_step.Y()), float(this_step.Z()-last_step.Z()) };
      float dirnorm = 0;
      for (int i=0; i<3; i++) {
	dirnorm += dir[i]*dir[i];
      }
      dirnorm = sqrt(dirnorm);
      if ( dirnorm<1.0e-3 )
	continue;
	
      for (int i=0; i<3; i++)
	dir[i] /= dirnorm;
      int nsteps=dirnorm/max_step_size+1;
      if ( nsteps<= 0 )
	nsteps = 1;
      float stepsize = dirnorm/float(nsteps);
      for (int istep=0; istep<nsteps; istep++) {
	std::vector<float> pos(4,0);
	std::vector<float> pos4v(4,0);
	pos4v[0] = last_step.T();	
	pos[0] = pos4v[1] = last_step.X() + stepsize*float(istep)*dir[0];
	pos[1] = pos4v[2] = last_step.Y() + stepsize*float(istep)*dir[1];
	pos[2] = pos4v[3] = last_step.Z() + stepsize*float(istep)*dir[2];

	int boundary_type = -1;        
	float fdwall = ublarcvapp::dwall( pos, boundary_type );
        
	std::vector<double> offset = psce->GetPosOffsets( pos[0], pos[1], pos[2] );
	std::vector<float> pos_sce(3);
	pos_sce[0] = pos[0]-(float)offset[0]+0.7;
	pos_sce[1] = pos[1]+(float)offset[1];
	pos_sce[2] = pos[2]+(float)offset[2];
	float tick = getTick( pos4v, trig_time, psce );
        
        if ( verbose ) {
          std::cout << " [step] dwall=" << fdwall << " tick=" << tick
                    << " pos-nosce=(" << pos[0] << "," << pos[1] << "," << pos[2] << ") "
                    << " sce-offset=(" << offset[0] << "," <<  offset[1] << "," << offset[2] << ") "            
            //<< " imgcoords=(row=" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")"
                    << std::endl;
        }
        
        if (fdwall<0.1 || (offset[0]==0.0 && offset[1]==0.0 && offset[2]==0.0 ))
          continue;        

	pos_sce[0] = (tick-3200.0)*cm_per_tick;
	if ( tick<meta.min_y()+1.0 || tick>meta.max_y()-1.0 )
	  continue;
        
	std::vector<int> imgcoords;
	try {
	  imgcoords = ublarcvapp::UBWireTool::getProjectedImagePixel( pos_sce, meta, 3 );
	}
	catch (...) {
	  //std::cout << std::endl;
	  continue;
	}
        
	// std::cout << " [" << thispt << "/" << npts << ":" << istep << "/" << nsteps << "] tick=" << tick;
	// if ( startAtstart )
	//   std::cout << " trackstart=(" << track.front().X() << "," << track.front().Y() << "," << track.front().Z() << ")";
	// else
	//   std::cout << " trackend=(" << track.back().X() << "," << track.back().Y() << "," << track.back().Z() << ")";	  
	// std::cout << " truepos=(" << pos[0] << "," << pos[1] << "," << pos[2] << ") ";
	// std::cout << " intime pos_sce=(" << pos_sce[0] << "," << pos_sce[1] << "," << pos_sce[2] << ") tick=" << tick << " ";
        endpt3d = pos_sce;
        if ( verbose )
          std::cout << "[return endpoint]" << std::endl;
        
	return imgcoords;
      }
      
    }//end of track loop

    // didn't find the crossing boundary
    if ( verbose )
      std::cout << "[return empty]" << std::endl;
    
    std::vector<int> empty;
    return empty;
  }
  
}
}
