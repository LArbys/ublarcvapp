#include "crossingPointsAnaMethods.h"

#include <cstring>
#include <stdexcept>

#include "TTree.h"
#include "TLorentzVector.h"

#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/trigger.h"

// ublarcvapp
#include "ublarcvapp/UBWireTool/UBWireTool.h"
#include "ublarcvapp/ubdllee/dwall.h"

namespace ublarcvapp {

  int doesTrackCrossImageBoundary( const larlite::mctrack& track, const larcv::ImageMeta& meta,
                                   const float trig_time, const larutil::SpaceChargeMicroBooNE* psce ) {
    
    float tick_start = getTick( track.front(), trig_time, psce );
    float tick_end   = getTick( track.back(), trig_time, psce );
    if ( tick_start>meta.min_y() && tick_start<meta.max_y() && tick_end>meta.min_y() && tick_end<meta.max_y() )
      return -1;

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

  std::vector<int> getImageBoundaryCrossingPoint( const larlite::mctrack& track, std::vector<float>& crossingpt,
                                                  const larcv::ImageMeta& meta,
						  const float boundary_tick_buffer, const float trig_time,
                                                  const larutil::SpaceChargeMicroBooNE* psce ) {
    
    if ( doesTrackCrossImageBoundary( track, meta, trig_time, psce )==-1 ) {
      std::stringstream msg;
      msg << __FILE__ << ":" << __LINE__ << " asking for bundary crossing point for track that does not cross the boundary" << std::endl;
      throw std::runtime_error( msg.str() );
    }
    
    crossingpt.resize(3,0);
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    
    float last_tick = getTick( track.front(), trig_time, psce );
    float fpos_last[3] = { (float)track.front().X(), (float)track.front().Y(), (float)track.front().Z() };
    for (int i=1; i<(int)track.size(); i++) {
      const auto& step_now = track[i];
      const auto& step_last = track[i-1];
      float fpos_now[3] = { (float)step_now.X(), (float)step_now.Y(), (float)step_now.Z() };
      float tick_now = getTick( step_now, trig_time, psce );
      float high_tick = tick_now;
      float low_tick  = last_tick;
      if ( low_tick>high_tick ) {
	high_tick = last_tick;
	low_tick = tick_now;
      }
      float boundary_tick = meta.min_y();
      bool crosses_bounds = false;
      if ( low_tick<meta.min_y() && high_tick>meta.min_y() ) {
	crosses_bounds = true;
	boundary_tick = meta.min_y() + boundary_tick_buffer;
      }
      else if ( low_tick<meta.max_y() && high_tick>meta.max_y() ) {
	crosses_bounds = true;
	boundary_tick = meta.max_y() - boundary_tick_buffer;	
      }

      /// go to next step, if not the boundary crossing step
      if ( crosses_bounds ) {

	//std::cout << "found crossing step: tick now=" << tick_now << " last=" << last_tick;
	
	float dir[3] = {0};
	float dirnorm = 0;
	for (int i=0; i<3; i++) {
	  dir[i] = fpos_now[i]-fpos_last[i];
	  dirnorm += dir[i]*dir[i];
	}
	dirnorm = sqrt(dirnorm);
	for (int i=0; i<3; i++) {
	  dir[i] /= dirnorm;
	}
	float dtick = boundary_tick-last_tick;
	float dcm   = cm_per_tick*dtick;
	for (int i=0; i<3; i++) {
	  crossingpt[i] = fpos_last[i] + dcm*dir[i];
	}
	//std::cout << " dtick=" << dtick << " dcm=" << dcm << std::endl;
	// get the image coordinates. Should be within the image now.
	if ( psce ) {
	  std::vector<double> offset = psce->GetPosOffsets( crossingpt[0], crossingpt[1], crossingpt[2] );
	  crossingpt[0] += -offset[0] + 0.7;
	  crossingpt[1] += offset[1];
	  crossingpt[2] += offset[2];
	}
	std::vector<int> crossing_imgcoords = ublarcvapp::UBWireTool::getProjectedImagePixel( crossingpt, meta, 3 );
	float finaltick = getTick( step_last, trig_time, psce ) + dtick;
	crossing_imgcoords[0] = meta.row( finaltick );
	return crossing_imgcoords;
      }//end of if crosses

      // if not, update last tick info and move on
      last_tick = tick_now;
      std::memcpy( fpos_last, fpos_now, sizeof(float)*3 );
    }
    // should not get here
    std::stringstream msg;
    msg << __FILE__ << ":" << __LINE__ << " routine should not get here." << std::endl;
    throw std::runtime_error( msg.str() );
    std::vector<int> empty(4,0);
    return empty;
  }
  

  float getTick( const larlite::mcstep& step, const float trig_time, const larutil::SpaceChargeMicroBooNE* psce ) {
    std::vector<float> pos(4,0);
    pos[0] = step.T();
    pos[1] = step.X();
    pos[2] = step.Y();
    pos[3] = step.Z();
    return getTick( pos, trig_time, psce );
  }
  

  float getTick( const std::vector<float>& step, const float trig_time, const larutil::SpaceChargeMicroBooNE* psce ) {    
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

  std::vector<float> getFirstStepPosInsideImage( const larlite::mctrack& track, const larcv::ImageMeta& meta,
                                                 const float trig_time,
						 const bool startAtstart, const float max_step_size, const float fv_border,
                                                 const larutil::SpaceChargeMicroBooNE* psce ) {
    // This function returns the (SCE-corrected) position where a MC track first is inside the image bounds
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;    
    int npts = (int)track.size();

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
	
	std::vector<double> offset = psce->GetPosOffsets( pos[0], pos[1], pos[2] );
	std::vector<float> pos_sce(3);
	pos_sce[0] = pos[0]-(float)offset[0]+0.7;
	pos_sce[1] = pos[1]+(float)offset[1];
	pos_sce[2] = pos[2]+(float)offset[2];
	int boundary_type = -1;
	float fdwall = ublarcvapp::dwall( pos_sce, boundary_type ); // use apparent distance...
	if ( fdwall<fv_border )
	  continue;

	float tick = getTick( pos4v, trig_time, psce );
	if ( tick<meta.min_y()+20.0 || tick>meta.max_y()-20.0 )
	  continue;

	// std::cout << " [" << thispt << "/" << npts << ":" << istep << "/" << nsteps << "] tick=" << tick;
	// if ( startAtstart )
	//   std::cout << " trackstart=(" << track.front().X() << "," << track.front().Y() << "," << track.front().Z() << ")";
	// else
	//   std::cout << " trackend=(" << track.back().X() << "," << track.back().Y() << "," << track.back().Z() << ")";	  
	// std::cout << " truepos=(" << pos[0] << "," << pos[1] << "," << pos[2] << ") ";
	// std::cout << " intime pos_sce=(" << pos_sce[0] << "," << pos_sce[1] << "," << pos_sce[2] << ") tick=" << tick << " ";
	pos_sce[0] = (tick-3200.0)*cm_per_tick; // we have to give the apparent-x (relative to the trigger) because we need to test the position in the image	
	std::vector<int> imgcoords;
	try {
	  imgcoords = ublarcvapp::UBWireTool::getProjectedImagePixel( pos_sce, meta, 3 );
	  //std::cout << " imgcoords=(row=" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")" << std::endl;
	}
	catch (...) {
	  //std::cout << std::endl;
	  continue;
	}
	
	return pos;
      }
      
    }//end of track loop

    // didn't find the crossing boundary
    std::vector<float> empty;
    return empty;
  }
  
}
