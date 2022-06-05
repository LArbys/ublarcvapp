#include "ublarcvapp/RecoTools/DetUtils.h"

namespace ublarcvapp {
namespace recotools {

  std::vector< const larcv::Image2D* >
  DetUtils::getTPCImages( const std::vector<larcv::Image2D>& adc_v,
			  const int tpcid, const int cryoid )
  {

    auto const geom = larlite::larutil::Geometry::GetME();
    
    // get the ADC images for the TPC: Make sure they are sorted
    struct PlaneImg_t
    {
      int planeid;
      const larcv::Image2D* pimg;
      PlaneImg_t( int plid, const larcv::Image2D* p )
	: planeid(plid),
	  pimg(p)
      {};
      bool operator<( const PlaneImg_t& rhs ) {
	if ( planeid<rhs.planeid )
	  return true;
	return false;
      };
    };
    
    std::vector< PlaneImg_t > sorter;
    std::vector< const larcv::Image2D* > ptpc_adc_v;
    for (int ii=0; ii<(int)adc_v.size(); ii++) {
      auto const& adc = adc_v[ii];
      int simpleindex = adc.meta().id();
      std::vector<int> ctp = geom->GetCTPfromSimplePlaneIndex( simpleindex );
      if ( ctp[0]==cryoid && ctp[1]==tpcid ) {
	sorter.push_back( PlaneImg_t( ctp[2], &adc ) );
      }
    }
    std::sort( sorter.begin(), sorter.end() );
    for ( auto& p_t : sorter ) {
      ptpc_adc_v.push_back( p_t.pimg );
    }
    
    return ptpc_adc_v;
    
  }
  

}
}
