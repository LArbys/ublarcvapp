#include "TruthTrackSCE.h"

namespace ublarcvapp {
namespace mctools {

  /**
   * @brief default constructor
   *
   *  This instance will make its own SpaceChargeMicroBooNE instance
   *
   */
  TruthTrackSCE::TruthTrackSCE()
    : larcv::larcv_base("TruthTrackSCE"),
    _kown_sce(true),
    _p_sce(nullptr)
  {
    _p_sce = new larutil::SpaceChargeMicroBooNE();
  }

  /**
   * @brief constructor using external spacechargemicroboone instance
   *
   * @param psce Pointer to instance of SpaceChargeMicroBooNE
   */  
  TruthTrackSCE::TruthTrackSCE( larutil::SpaceChargeMicroBooNE* psce )
    : larcv::larcv_base("TruthTrackSCE"),
    _kown_sce(false),
    _p_sce(psce)
  {
  }

  TruthTrackSCE::~TruthTrackSCE()
  {
    if ( _kown_sce && _p_sce )
      delete _p_sce;
    _p_sce = nullptr;
  }

  /**
   * @brief apply space charge effect along a true trajectory
   *
   * @param[in] mct True simulation trajectory
   * @return track with space-charge applied points
   */
  larlite::track TruthTrackSCE::applySCE( const larlite::mctrack& mct )
  {
    larlite::track ll;
    ll.reserve( mct.size() );

    TVector3 last_pt;

    int npts = 0;
    for ( auto const& step : mct ) {
      TVector3 start = { step.Position()[0], step.Position()[1], step.Position()[2] };
      std::vector<double> s_offset = _p_sce->GetPosOffsets(start[0],start[1],start[2]);
      TVector3 start_sce;
      start_sce[0] = start[0] - s_offset[0] + 0.7;
      start_sce[1] = start[1] + s_offset[1];
      start_sce[2] = start[2] + s_offset[2];


      LARCV_DEBUG() << "convert (" << start[0] << "," << start[1] << "," << start[2] << ") "
                    << "to (" << start_sce[0] << "," << start_sce[1] << "," << start_sce[2] << ")" << std::endl;      
      ll.add_vertex( start_sce );
      if ( npts>0 ) {
        TVector3 dir = start_sce - last_pt;
        float norm = dir.Mag();
        if ( norm>0 ) {
          for (int i=0; i<3; i++)
            dir[i] /= norm;
        }
        ll.add_direction(dir);
        if ( npts==1 )
          ll.add_direction(dir); // add extra for first point
      }
      last_pt = start_sce;
      npts++;
    }

    LARCV_DEBUG() << "made track with " << ll.NumberTrajectoryPoints() << " points." << std::endl;
    return ll;
  }

  /**
   * @brief find closest distance to series of line segments
   *
   * @param[in] testpt Test 3D point
   * @param[in] track  series of line segments
   * @param[out] min_r smallest distance to line segments
   * @param[in] min_step  index of line segment that point was closest to
   */
  void TruthTrackSCE::dist2track( const std::vector<float>& testpt,
                                  const larlite::track& track,
                                  float& min_r,
                                  int& min_step )
  {

    min_r = 1.0e9;
    min_step = -1;

    LARCV_DEBUG() << "start test" << std::endl;
    
    for (int i=0; i<(int)track.NumberTrajectoryPoints()-1;i++) {
      std::vector<float> start(3,0);
      std::vector<float> end(3,0);
      float s=0;
      float len1 = 0;
      float len2 = 0;
      float dist1 = 0;
      float dist2 = 0;
      for (int v=0; v<3; v++) {
        start[v] = track.LocationAtPoint(i)[v];
        end[v]   = track.LocationAtPoint(i+1)[v];
        s += (testpt[v]-start[v])*(end[v]-start[v]);
        len1 += (testpt[v]-start[v])*(testpt[v]-start[v]);
        len2 += (end[v]-start[v])*(end[v]-start[v]);
        dist1 += (testpt[v]-start[v])*(testpt[v]-start[v]);
        dist2 += (testpt[v]-end[v])*(testpt[v]-end[v]);
      }
      len1 = sqrt(len1);
      len2 = sqrt(len2);
      dist1 = sqrt(dist1);
      dist2 = sqrt(dist2);      
      if ( len2>0 )
        s /= (len2);

      float r = 1.0e9;
      if (len2>0)
        r = pointLineDistance( start, end, testpt );

      LARCV_DEBUG() << "step[" << i << "] len=" << len2 << " r=" << r << " s=" << s << " vs (" << 0 << "," << len2 << ")" << std::endl;

      // determine if along the track
      if ( len2>0 && s>0 && s<len2 ) {
        if ( min_r>r ) {
          min_r = r;
          min_step = i;
        }
      }
      else if ( min_r>dist1 ) {
        min_r = dist1;
        min_step = i;
      }
      else if ( min_r>dist2 ) {
        min_r = dist2;
        min_step = i;
      }
    }

    LARCV_DEBUG() << "result: min_r=" << min_r << " step=" << min_step << std::endl;
    
  }

  /**
   * @brief distance from point to single line segment
   *
   * @param[in] linept1 one end of line segment
   * @param[in] linept2 other end of line segment
   * @param[in] pt  test point
   * @return distance of point from line segment
   */
  float TruthTrackSCE::pointLineDistance( const std::vector<float>& linept1,
                                          const std::vector<float>& linept2,
                                          const std::vector<float>& pt )
  {
    
    
    std::vector<float> d1(3);
    std::vector<float> d2(3);

    float len1 = 0.;
    float linelen = 0.;
    for (int i=0; i<3; i++ ) {
      d1[i] = pt[i] - linept1[i];
      d2[i] = pt[i] - linept2[i];
      len1 += d1[i]*d1[i];
      linelen += (linept1[i]-linept2[i])*(linept1[i]-linept2[i]);
    }
    len1 = sqrt(len1);
    linelen = sqrt(linelen);

    if ( linelen<1.0e-4 ) {
      // short cluster, use distance to end point
      return len1;
    }

    // cross-product
    std::vector<float> d1xd2(3);
    d1xd2[0] =  d1[1]*d2[2] - d1[2]*d2[1];
    d1xd2[1] = -d1[0]*d2[2] + d1[2]*d2[0];
    d1xd2[2] =  d1[0]*d2[1] - d1[1]*d2[0];
    float len1x2 = 0.;
    for ( int i=0; i<3; i++ ) {
      len1x2 += d1xd2[i]*d1xd2[i];
    }
    len1x2 = sqrt(len1x2);
    float r = len1x2/linelen;
    return r;
  }

  

}
}
