#include "TrackReverser.h"

namespace ublarcvapp {
namespace reco3d {

  larlite::track
  TrackReverser::reverseTrack( const larlite::track& original,
                               bool flip_direction )
  {
    
    int npts = original.NumberTrajectoryPoints();

    larlite::track reversed;
    reversed.reserve(npts);

    bool has_cov    = ( (int)original.NumberCovariance()==npts )  ? true : false;
    bool has_fitmom = ( (int)original.NumberFitMomentum()==npts ) ? true : false;
    bool has_dqdx   = false;
    if ( (int)original.NumberdQdx((larlite::geo::View_t)0)==4 )
      has_dqdx = true;
    
    int nstored = 0;

    std::vector< std::vector<double> > dqdx_vv;
    if ( has_dqdx ) {
      dqdx_vv.resize( 4 );
      for (size_t v=0; v<4; v++)
        dqdx_vv[v].reserve( npts );
    }
    
    for (int ipt=npts-1; ipt>=0; ipt-- ) {
      // set the position
      TVector3 pos = original.LocationAtPoint(ipt);
      reversed.add_vertex( pos );

      // set the direction
      TVector3 dir;      
      if ( ipt-1>=0 ) {
        // next point available
        TVector3 nextpos  = original.LocationAtPoint(ipt-1);
        TVector3 step_dir = nextpos-pos;
        if ( step_dir.Mag()>0 ) {
          float len = step_dir.Mag();
          for (int i=0; i<3; i++)
            step_dir[i] /= len;
        }
        dir = step_dir;
      }
      else if ( ipt==0 ) {
        dir = reversed.DirectionAtPoint( nstored-1 );
      }
      reversed.add_direction(dir);

      // set the covariance
      if ( has_cov ) {
        reversed.add_covariance( original.CovarianceAtPoint(ipt) );
      }

      // set fit momentum: this seems nonsensical?
      if ( has_fitmom ) {
        reversed.add_momentum( original.MomentumAtPoint(ipt) );
      }

      // set dqdx
      for (size_t v=0; v<4; v++)
        dqdx_vv[v].push_back( original.DQdxAtPoint( ipt, (larlite::geo::View_t)v ) );
      
      nstored++;
    }//end of loop over points

    for (auto& dqdx_v : dqdx_vv )
      reversed.add_dqdx( dqdx_v );

    return reversed;
  }
  
  
}
}
