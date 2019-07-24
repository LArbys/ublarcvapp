#include "BoundarySpacePoint.h"
#include <cmath>
#include <vector>
#include <sstream>
#include <string>

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

// larcv
#include "larcv/core/DataFormat/ImageMeta.h"
#include "ublarcvapp/UBWireTool/UBWireTool.h"

namespace ublarcvapp {
namespace tagger {

  BoundarySpacePoint::BoundarySpacePoint( BoundaryEnd_t type, const std::vector<float>& pos,
                                          const std::vector<float>& dir, const larcv::ImageMeta& meta  )
  {
    std::vector<int> imgcoords = ublarcvapp::UBWireTool::getProjectedImagePixel( pos, meta, 3 );

    //if ( imgcoords[0]<0 ) imgcoords[0] = 0;
    //if ( imgcoords[0]>=(int)meta.rows() ) imgcoords[0] = meta.rows()-1;

    try {
      for (int p=0; p<3; p++) {
	//int col = imgcoords[p+1];
	(*this).push_back( BoundaryEndPt(imgcoords[0],imgcoords[p+1],type) );
      }
    }
    catch (std::exception& e) {
      std::stringstream msg;
      msg << __FILE__ << "::" << __LINE__ << "Error making BoundaryEndPts: " << e.what() << std::endl;
      msg << " Input position (" << pos[0] << "," << pos[1] << "," << pos[2] << ")" << std::endl;
      msg << " Imgcoords: (" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")" << std::endl;
      throw std::runtime_error( msg.str() );
    }
    
    m_pos = pos;
    m_dir = dir;
    m_empty = false;
    boundary_type = type;
  }

  BoundarySpacePoint::BoundarySpacePoint( BoundaryEnd_t type,
                                          const std::vector<float>& pos, const larcv::ImageMeta& meta  )
  {
    std::vector<int> imgcoords = ublarcvapp::UBWireTool::getProjectedImagePixel( pos, meta, 3 );
    try {
      for (int p=0; p<3; p++)
	(*this).push_back( BoundaryEndPt(imgcoords[0],imgcoords[p+1],type) );
    }
    catch (std::exception& e) {
      std::stringstream msg;
      msg << __FILE__ << "::" << __LINE__ << "Error making BoundaryEndPts: " << e.what() << std::endl;
      msg << " Input position (" << pos[0] << "," << pos[1] << "," << pos[2] << ")" << std::endl;
      throw std::runtime_error( msg.str() );
    }
    m_pos = pos;
    m_dir.resize(3,0);
    m_empty = false;
    boundary_type = type;
  }

  float BoundarySpacePoint::dwall() const {
    float dwall = 0.0;
    float dy,dz;
    switch( type() ) {
    case ublarcvapp::tagger::kTop:
      dwall = 118.0-m_pos[1];
      break;
    case ublarcvapp::tagger::kBottom:
      dwall = 118.0+m_pos[1];
      break;
    case ublarcvapp::tagger::kUpstream:
      dwall = m_pos[2];
      break;
    case ublarcvapp::tagger::kDownstream:
      dwall = 1037-m_pos[2];
      break;
    case ublarcvapp::tagger::kAnode:
    case ublarcvapp::tagger::kCathode:
    case ublarcvapp::tagger::kImageEnd:
      dy = ( fabs(118.0-m_pos[1])<fabs(118.0+m_pos[1]) ) ? 118.0-m_pos[1] : 118.0+m_pos[1];
      dz = ( fabs(m_pos[2]) < fabs(1037-m_pos[2]) ) ? m_pos[2] : 1037-m_pos[2];
      dwall = ( fabs(dy)<fabs(dz ) ) ? dy : dz;
      break;
    default:
      std::runtime_error("BoundarySpacePoint::dwall[error] cannot calculate dwall for undefined boundary type");
      break;
    }
    return dwall;
  }

  void BoundarySpacePoint::setFlashIndex( int ivec, int idx, const larlite::opflash* popfl ) {
    m_flashidx.ivec = ivec;
    m_flashidx.idx  = idx;
    m_flashidx.popflash = popfl;
  }

  void BoundarySpacePoint::setup( const larcv::ImageMeta& meta ) {
    if ( m_pos.size()==0 ) {
      m_pos.resize(3,0);
      m_dir.resize(3,0);

      // warning: unprotected use of meta

      float x = (meta.pos_y( front().row )-3200.0)*larutil::LArProperties::GetME()->DriftVelocity()*0.5;
      int crosses;
      std::vector<float> intersection;
      double triarea;
      std::vector< int > wids;
      for ( size_t p=0; p<size(); p++) {
	wids.push_back( meta.pos_x( at(p).col ) );
      }
      ublarcvapp::UBWireTool::wireIntersection( wids, intersection, triarea, crosses );
      m_pos[0] = x;
      m_pos[1] = intersection[1];
      m_pos[2] = intersection[0];
    }
  }

  int BoundarySpacePoint::tick( const larcv::ImageMeta& meta ) const {
    float tick = pos()[0]/(::larutil::LArProperties::GetME()->DriftVelocity()*0.5) + 3200.0;
    return (int)tick;
  }

  std::vector<int> BoundarySpacePoint::wires( const larcv::ImageMeta& meta ) const {
    std::vector<int> wires;
    Double_t xyz[3] = { pos()[0], pos()[1], pos()[2] };
    for (size_t p=0; p<3; p++) {
      wires.push_back( larutil::Geometry::GetME()->WireCoordinate( xyz, p ) );
    }
    return wires;
  }

  std::string BoundarySpacePoint::printImageCoords( const larcv::ImageMeta& meta ) const {
    int thetick = tick(meta);
    std::vector<int> thewires = wires(meta);
    std::stringstream ss;
    ss << "(" << thetick << ", " << thewires[0] << ", " << thewires[1] << ", " << thewires[2] << ")";
    return ss.str();
  }
}
}
