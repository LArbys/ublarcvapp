#include "PointImageProjection.h"

#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/DetectorProperties.h"
#include "larlite/LArUtil/LArProperties.h"

namespace ublarcvapp {
namespace ubimagemod {

  /** 
   * @brief project a position inside the tpc into a pixel in the wireplane image and sum pixels around it
   *
   * 
   */
  float PointImageProjection::getPixelSumAroundProjPoint( const std::vector<float>& xyz,
							  const larcv::Image2D& img,
							  int pixel_kernel_radius,
							  float pixval_threshold ) const
  {

    // get wire and time given xyz
    TVector3 worldLoc( xyz[0], xyz[1], xyz[2] );

    //LARCV_NORMAL() << img.meta().dump() << std::endl;
    
    UInt_t wireid = larutil::Geometry::GetME()->NearestWire( worldLoc, img.meta().plane() );
    if ( wireid < img.meta().min_x() ) {
      LARCV_NORMAL() << "Wire " << wireid << " below " << img.meta().min_x() << std::endl;
      return 0.0;
    }
    if ( wireid >= img.meta().max_x() ) {
      LARCV_NORMAL() << "Wire " << wireid << " above max=" << img.meta().max_x() << std::endl;      
      return 0.0;
    }

    //double tick = larutil::DetectorProperties::GetME()->ConvertXToTicks( xyz[0], img.meta().plane() );
    float cm_per_tick = larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    double tick  = xyz[0]/cm_per_tick + 3200.0;
    if (  tick < img.meta().min_y() ) {
      LARCV_NORMAL() << "tick " << tick << " below min=" << img.meta().min_y() << std::endl;      
      return 0.0;
    }

    if ( tick >= img.meta().max_y() ) {
      LARCV_NORMAL() << "tick " << tick << " greater than max=" << img.meta().max_y() << std::endl;
      return 0.0;
    }

    int col = (int)img.meta().col( (float)wireid );
    int row = (int)img.meta().row( (float)tick );

    int dkr = abs(pixel_kernel_radius);
    int dkc = abs(pixel_kernel_radius);

    int startcol = ((col-dkc)<0) ? 0 : col-dkc;
    int endcol   = ((col+dkc)>=(int)img.meta().cols()) ? (int)img.meta().cols()-1 : col+dkc;
    int startrow = ((row-dkr)<0) ? 0 : row-dkr;
    int endrow   = ((row+dkr)>=(int)img.meta().rows()) ? (int)img.meta().rows()-1 : row+dkr;

    float pixsum = 0.;
    for (int pixrow=startrow; pixrow<=endrow; pixrow++) {
      for (int pixcol=startcol; pixcol<=endcol; pixcol++) {
	float pixval = img.pixel( pixrow, pixcol );
	// only add to sum if pixel value above threshold
	if ( pixval > pixval_threshold )
	  pixsum += pixval;
      }//end of dc loop
    }//end of dr loop
    
    return pixsum;
    
  }

}
}
