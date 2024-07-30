#include "MCPos2ImageUtils.h"

#include "larlite/LArUtil/DetectorProperties.h"
#include "larlite/LArUtil/LArProperties.h"
#include "larlite/LArUtil/Geometry.h"

namespace ublarcvapp {
namespace mctools {

  MCPos2ImageUtils* MCPos2ImageUtils::_p_singleton = nullptr;
  
  std::vector< std::vector<float> >
  MCPos2ImageUtils::getRecoSpacepoints( const larlite::mctrack& track,
					bool apply_t0_shift,
					bool apply_sce )
  {
    std::vector< std::vector<float> > reco( track.size() );

    size_t istep = 0;
    for (auto const& step : track ) {
      reco[istep] = truepos_to_recopos( step.X(), step.Y(), step.Z(), step.T(),
					apply_t0_shift, apply_sce );
      istep++;
    }


    return reco;
  }

  std::vector<float>
  MCPos2ImageUtils::truepos_to_recopos( const float x_cm, const float y_cm, const float z_cm,
					const float t_ns,
					bool apply_t0_shift,
					bool apply_sce )
  {
    
    float v_cm_per_us    = ::larutil::LArProperties::GetME()->DriftVelocity();
    float us_per_tick   = (::larutil::DetectorProperties::GetME()->SamplingRate()*1.0e-3); // value return is in MHz, convert to usec
    float trig_g4time_ns = 4050.0;
    float dx_trigger_cm  = 0.0;

    // std::cout << "v_cm_per_us: " << v_cm_per_us << std::endl;
    // std::cout << "us_per_tick: " << us_per_tick << std::endl;
    
    std::vector<float> pos_sce = { x_cm, y_cm, z_cm };
    if ( apply_sce ) {
      pos_sce = get_sce_shifted_pos( x_cm, y_cm, z_cm );
    }
    
    std::vector<float> reco_pos = pos_sce;

    if ( apply_t0_shift ) {
      // apparent depth from the anode due to relative time to trigger
      dx_trigger_cm = (t_ns - trig_g4time_ns)*1.0e-3*v_cm_per_us;
    }
    reco_pos[0] += dx_trigger_cm;    

    // append the tick the charge should appear at
    //float tick = larutil::DetectorProperties::GetME()->ConvertXToTicks( reco_pos[0], 0 ); // this is busted
    float tickx = (reco_pos[0])/v_cm_per_us/us_per_tick + 3200.0;

    //std::cout << "tick=" << tick << " vs " << tickx << std::endl; 
    reco_pos.push_back( tickx );

    return reco_pos;
  }


  /**
   * @brief get the reco position due to electric field non-uniformity
   *
   */
  std::vector<float>
  MCPos2ImageUtils::get_sce_shifted_pos( const float x, const float y, const float z )
  {
    
    std::vector<float> out = { x, y, z };
    
    if ( psce ) {
      std::vector<double> sce_shift = psce->GetPosOffsets( x, y, z );
      out[0] += -sce_shift[0]+0.7;
      out[1] += sce_shift[1];
      out[2] += sce_shift[2];
    }

    return out;
  }

  std::vector<float> MCPos2ImageUtils::to_imagepos( const float x, const float y, const float z, const float t_ns )
  {

    float v_cm_per_us    = ::larutil::LArProperties::GetME()->DriftVelocity();
    float us_per_tick   = (::larutil::DetectorProperties::GetME()->SamplingRate()*1.0e-3); // value return is in MHz, convert to usec
    
    TVector3 vpos;
    vpos[0] = x;
    vpos[1] = y;
    vpos[2] = z;
    //std::cout << "vpos: " << vpos[0] << " " << vpos[1] << " " << vpos[2] << std::endl;

    std::vector<float> imgpos(4,0);
    if ( std::fabs(vpos[1])>116.5 || vpos[2]<0.0 || vpos[2]>1036.0 ) 
      return imgpos;
    
    for (int p=0; p<3; p++) {
      float wire = (float)larutil::Geometry::GetME()->NearestWire( vpos, p );
      imgpos[p] = wire;
    }
    float tickx = vpos[0]/v_cm_per_us/us_per_tick + 3200.0;
    imgpos[3] = tickx;
    
    return imgpos;
  }

  std::vector<float> MCPos2ImageUtils::truepos_to_imagepos( const float x, const float y, const float z,
							    const float t_ns, bool apply_sce )
  {
    std::vector<float> recopos = truepos_to_recopos( x, y, z, t_ns, true, apply_sce );
    TVector3 vpos;
    vpos[0] = recopos[0];
    vpos[1] = recopos[1];
    vpos[2] = recopos[2];
    //std::cout << "vpos: " << vpos[0] << " " << vpos[1] << " " << vpos[2] << std::endl;

    std::vector<float> imgpos(4,0);
    if ( std::fabs(vpos[1])>116.5 || vpos[2]<0.0 || vpos[2]>1036.0 ) 
      return imgpos;

    imgpos[3] = recopos[3]; // should be the tick

    //const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    //float tick = ( recopos[3]-4.050 )/0.5 + recopos[0]/cm_per_tick + 3200.0;
    for (int p=0; p<3; p++) {
      float wire = (float)larutil::Geometry::GetME()->NearestWire( vpos, p );
      imgpos[p] = wire;
    }
	
    return imgpos;
  }

  
  
  
  
}
}
