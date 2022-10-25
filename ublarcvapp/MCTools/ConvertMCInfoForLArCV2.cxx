#include "ConvertMCInfoForLArCV2.h"

namespace ublarcvapp {
namespace mctools {

  std::vector<larcv::Particle>
  ConvertMCInfoForLArCV2::convert( const std::vector<larlite::mctrack>& mctrack_v,
				   const std::vector<larlite::mcshower>& mcshower_v  )
  {

    std::vector<larcv::Particle> particle_v;

    // convert track objects
    for ( auto const& mcpart : mctrack_v ) {
      
      larcv::Particle lcvparticle( larcv::kShapeTrack );

      // PPN minimal info needed
      lcvparticle.id( mcpart.TrackID() );
      lcvparticle.track_id( mcpart.TrackID() );
      lcvparticle.pdg_code( mcpart.PdgCode() );
      lcvparticle.creation_process( mcpart.Process() );
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
      lcvparticle.first_step( mcpart.front().X(), mcpart.front().Y(), mcpart.front().Z(), mcpart.front().T() );
      lcvparticle.last_step(  mcpart.back().X(), mcpart.back().Y(), mcpart.back().Z(), mcpart.back().T() );
      
      particle_v.emplace_back( std::move(lcvparticle) );
    }

    // convert shower objects
    for ( auto const& mcpart : mcshower_v ) {
      
      larcv::Particle lcvparticle( larcv::kShapeShower );

      // PPN minimal info needed
      lcvparticle.pdg_code( mcpart.PdgCode() );
      lcvparticle.creation_process( mcpart.Process() );

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
	lcvparticle.first_step( mcpart.DetProfile().X(), mcpart.DetProfile().Y(), mcpart.DetProfile().Z(), mcpart.DetProfile().T() );
      }
      else {
	edep = 0.;
	// set first step as creation point
	lcvparticle.first_step( mcpart.Start().X(), mcpart.Start().Y(), mcpart.Start().Z(), mcpart.Start().T() );
      }
      lcvparticle.energy_deposit( edep );
      lcvparticle.num_voxels( nvoxel ); // not filled in this context

      particle_v.emplace_back( std::move(lcvparticle) );
    }

    return particle_v;
  }

}
}
