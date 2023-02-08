#include "ConvertMCInfoForLArCV2.h"

#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"
#include "larlite/LArUtil/DetectorProperties.h"

namespace ublarcvapp {
namespace mctools {

  ConvertMCInfoForLArCV2::~ConvertMCInfoForLArCV2()
  {
    if ( _psce )
      delete _psce;
    _psce = nullptr;
  }
  
  std::vector<larcv::Particle>
  ConvertMCInfoForLArCV2::convert( const std::vector<larlite::mctrack>& mctrack_v,
				   const std::vector<larlite::mcshower>& mcshower_v  )
  {

    LARCV_DEBUG() << "convert " << mctrack_v.size() << " tracks and " << mcshower_v.size() << " showers" << std::endl;
    
    if ( kApplySCE && _psce==nullptr ) {
      LARCV_DEBUG() << "Load Space Charge Module" << std::endl;
      _psce = new larutil::SpaceChargeMicroBooNE();
    }
    
    std::vector<larcv::Particle> particle_v;

    // convert track objects
    for ( auto const& mcpart : mctrack_v ) {

      if ( mcpart.size()<2 )
	continue;
      
      larcv::Particle lcvparticle( larcv::kShapeTrack );

      // PPN minimal info needed
      lcvparticle.id( mcpart.TrackID() );
      lcvparticle.track_id( (unsigned int)mcpart.TrackID() );
      lcvparticle.pdg_code( mcpart.PdgCode() );
      lcvparticle.creation_process( mcpart.Process() );
      lcvparticle.momentum( mcpart.front().Px(), mcpart.front().Py(), mcpart.front().Pz() );
      lcvparticle.energy_init( mcpart.front().E() );
      lcvparticle.parent_track_id( (unsigned int)mcpart.MotherTrackID() );
      lcvparticle.parent_pdg_code( mcpart.MotherPdgCode() );
      lcvparticle.parent_creation_process( mcpart.MotherProcess() );            
      lcvparticle.ancestor_track_id( (unsigned int)mcpart.AncestorTrackID() );      
      lcvparticle.ancestor_pdg_code( mcpart.AncestorPdgCode() );
      lcvparticle.ancestor_creation_process( mcpart.AncestorProcess() );            
      
      float edep = 0.0; // total energy deposited. could use mcstep or mcpart.fdEdx vector if filled
      int nvoxel = 0;
      if ( mcpart.dEdx().size()>0 ) {
	for ( auto const& dedx : mcpart.dEdx() )
	  edep += (float)dedx;
      }
      else {
	// sum the energy loss over mcsteps
	edep = mcpart.front().E()-mcpart.back().E();
      }
      
      lcvparticle.energy_deposit( edep );
      lcvparticle.num_voxels( nvoxel ); // not filled in this context

      std::vector<double> first_step = preparePosition( mcpart.front().X(), mcpart.front().Y(), mcpart.front().Z(), mcpart.front().T() );
      std::vector<double> last_step  = preparePosition( mcpart.back().X(), mcpart.back().Y(), mcpart.back().Z(), mcpart.back().T() );
      lcvparticle.position( first_step[0], first_step[1], first_step[2], first_step[3] );
      lcvparticle.first_step( first_step[0], first_step[1], first_step[2], first_step[3] );
      lcvparticle.last_step(  last_step[0], last_step[1], last_step[2], last_step[3] );

      LARCV_DEBUG() << "Converted track: id=" << mcpart.TrackID() << std::endl;
      LARCV_DEBUG() << lcvparticle.dump() << std::endl;
      
      particle_v.emplace_back( std::move(lcvparticle) );
    }

    // convert shower objects
    for ( auto const& mcpart : mcshower_v ) {

      larcv::Particle lcvparticle( larcv::kShapeShower );

      // PPN minimal info needed
      lcvparticle.pdg_code( mcpart.PdgCode() );
      lcvparticle.track_id( (unsigned int)mcpart.TrackID() );
      lcvparticle.creation_process( mcpart.Process() );
      lcvparticle.momentum( mcpart.DetProfile().Px(), mcpart.DetProfile().Py(), mcpart.DetProfile().Pz() );
      lcvparticle.energy_init( mcpart.DetProfile().E() );      
      lcvparticle.parent_track_id( (unsigned int)mcpart.MotherTrackID() );
      lcvparticle.parent_pdg_code( mcpart.MotherPdgCode() );
      lcvparticle.parent_creation_process( mcpart.MotherProcess() );            
      lcvparticle.ancestor_track_id( (unsigned int)mcpart.AncestorTrackID() );      
      lcvparticle.ancestor_pdg_code( mcpart.AncestorPdgCode() );
      lcvparticle.ancestor_creation_process( mcpart.AncestorProcess() );            
      
      // check the content of the DetProfile momentum
      bool in_tpc = true;
      auto const& detprofmom = mcpart.DetProfile().Momentum().Vect();
      if ( detprofmom[0]==0 && detprofmom[1]==0 && detprofmom[2]==0 )
	in_tpc = false;

      float edep = 0.0;
      int nvoxel = 0;
      if ( in_tpc ) {      
	edep = mcpart.DetProfile().E(); // total energy deposited. could use mcstep or mcpart.fdEdx vector if filled
	// set first step as place that shower first starts to deposit energy in TPC
	std::vector<double> first_step = preparePosition( mcpart.DetProfile().X(), mcpart.DetProfile().Y(), mcpart.DetProfile().Z(), mcpart.DetProfile().T() );
	lcvparticle.position( first_step[0], first_step[1], first_step[2], first_step[3] );	
	lcvparticle.first_step( first_step[0], first_step[1], first_step[2], first_step[3] );
	lcvparticle.last_step( 0, 0, 0, 0 );
      }
      else {
	edep = 0.;
	// set first step as creation point
	lcvparticle.position( mcpart.Start().X(), mcpart.Start().Y(), mcpart.Start().Z(), mcpart.Start().T() );	
	lcvparticle.first_step( mcpart.Start().X(), mcpart.Start().Y(), mcpart.Start().Z(), mcpart.Start().T() );
	lcvparticle.last_step( 0, 0, 0, 0 );	
      }
      lcvparticle.energy_deposit( edep );
      lcvparticle.num_voxels( nvoxel ); // not filled in this context

      LARCV_DEBUG() << "Converted track: id=" << mcpart.TrackID() << std::endl;      
      LARCV_DEBUG() << lcvparticle.dump() << std::endl;
      
      particle_v.emplace_back( std::move(lcvparticle) );
    }

    return particle_v;
  }

  std::vector<double> ConvertMCInfoForLArCV2::preparePosition( double x, double y, double z, double t )
  {
    std::vector<double> step = { x, y, z, t };


    if ( kApplySCE && _psce ) {
      std::vector<double> offsets = _psce->GetPosOffsets( step[0], step[1], step[2] );
      step[0] = step[0]+0.7-offsets[0];
      step[1] = step[1]+offsets[1];
      step[2] = step[2]+offsets[2];
    }
    if ( kApplyT0drift ) {
	// relative time to trigger converted to TPC ticks elapsed
      float step_dtick = ( step[3]*1.0e-3 )/(larutil::DetectorProperties::GetME()->SamplingRate()*1.0e-3);
      // drift time
      float driftx_step = step_dtick*larutil::DetectorProperties::GetME()->GetXTicksCoefficient();
      step[0] += driftx_step;
    }
    return step;
  }

}
}
