#include "LArliteManager.h"

namespace ublarcvapp {

  /**
   * constructor. wrapper to larlite::storage_manger class
   *
   * @param[in] mode Possible modes: { kREAD, kWRITE, kBOTH }
   *
   */
  LArliteManager::LArliteManager( larlite::storage_manager::IOMode_t mode, std::string name )
    : larlite::storage_manager(mode), larcv::larcv_base(name)
  {
    m_entry_map.clear();
    m_current_rse   = {-1,-1,-1};
    m_current_entry =  -1;
    m_last_rse   = {-1,-1,-1};
    m_last_entry =  -1;
  }
  
  bool LArliteManager::syncEntry( const larcv::IOManager& iolarcv, bool force_reload ) {
    
    rse_t rse_larcv = { (int)iolarcv.event_id().run(),
                        (int)iolarcv.event_id().subrun(),
                        (int)iolarcv.event_id().event() };

    LARCV_DEBUG() << " sync with larcv "
              << "(" << rse_larcv[0] << "," << rse_larcv[1] << "," << rse_larcv[2] << ")"
              << std::endl;

    if ( rse_larcv[0]==-1 || rse_larcv[1]==-1 || rse_larcv[2]==-1 ) {
      LARCV_CRITICAL() << "Syncing to invalid (run,subrun,event)="
                       << "(" << rse_larcv[0] << "," << rse_larcv[1] << "," << rse_larcv[2] << "). "
                       << " typically need to use 'larcv::IOManager::get_data' before calling this."
                       << std::endl;
      throw std::runtime_error("syncing to invalid (run,subrun,event)");
    }

    if ( !go_to( rse_larcv[0], rse_larcv[1], rse_larcv[2] ) ) {
      LARCV_CRITICAL() << "Could not find entry with rse=(" << rse_larcv[0] << "," << rse_larcv[1] << "," << rse_larcv[2] << ")!" << std::endl;
      return false;
    }
    LARCV_INFO() << "Sync'd w/ rse=(" << rse_larcv[0] << "," << rse_larcv[1] << "," << rse_larcv[2] << ")!" << std::endl;    
    return true;
  }
  
}
