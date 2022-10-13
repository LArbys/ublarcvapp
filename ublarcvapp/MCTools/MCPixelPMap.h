#ifndef __MCPIXELPDICT_H__
#define __MCPIXELPDICT_H__

#include <string>
#include <vector>
#include <map>

#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larlite/DataFormat/storage_manager.h"
#include "MCPixelPGraph.h"

namespace ublarcvapp {
namespace mctools {

  class MCPixelPMap : public larcv::larcv_base {

    public:

      MCPixelPMap() : larcv::larcv_base("MCPixelPMap"), adc_tree("wire") {};
      virtual ~MCPixelPMap() {};

      struct PixLoc_t {
        int plane, row, col;
        PixLoc_t() : plane(-1), row(-1), col(-1) {};
        PixLoc_t(int _plane, int _row, int _col) : plane(_plane), row(_row), col(_col) {};
        bool operator==( const PixLoc_t& rhs ) const {
          if(plane == rhs.plane && row == rhs.row && col == rhs.col) return true;
          return false;
        };
        bool operator<( const PixLoc_t& rhs ) const {
          if(plane < rhs.plane) return true;
          if(plane > rhs.plane) return false;
          if(row < rhs.row) return true;
          if(row > rhs.row) return false;
          if(col < rhs.col) return true;
          return false;
        }
      };

      struct PixCont_t {
        int pdg;
        float pixI;
        PixCont_t() : pdg(0), pixI(0.) {};
        PixCont_t(int _pdg, float _pixI) : pdg(_pdg), pixI(_pixI) {};
        bool operator>(const PixCont_t& rhs) const {
          if(pixI > rhs.pixI) return true;
          return false;
        }
      };

      void buildmap( larcv::IOManager& iolcv, larlite::storage_manager& ioll );
      void buildmap( larcv::IOManager& iolcv, ublarcvapp::mctools::MCPixelPGraph& pixGraph );

      std::vector<PixCont_t> getPixContent(int plane, int row, int col);
      PixCont_t getPixParticle(int plane, int row, int col);
      int getPixParticlePDG(int plane, int row, int col);

      std::map<PixLoc_t, std::vector<PixCont_t> > pixMap;
      std::string adc_tree;
      void set_adc_treename( std::string name ) { adc_tree = name; };

  };

}
}

#endif
