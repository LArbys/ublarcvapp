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

      struct PixPart_t {
        int nodeidx;
        int tid;
        int pdg;
        PixPart_t() : nodeidx(-1), tid(-1), pdg(0) {};
        PixPart_t(int _nodeidx, int _tid, int _pdg) : nodeidx(_nodeidx), tid(_tid), pdg(_pdg) {};
        bool operator>(const PixPart_t& rhs) const {
          if(nodeidx > rhs.nodeidx) return true;
          return false;
        }
      };

      struct PixCont_t {
        float pixI;
        std::vector<PixPart_t> particles;
        PixCont_t() : pixI(0.) {};
        PixCont_t(float _pixI) : pixI(_pixI) {};
        bool operator>(const PixCont_t& rhs) const {
          if(pixI > rhs.pixI) return true;
          return false;
        }
      };

      void buildmap( larcv::IOManager& iolcv, larlite::storage_manager& ioll );
      void buildmap( larcv::IOManager& iolcv, ublarcvapp::mctools::MCPixelPGraph& pixGraph );

      PixCont_t getPixContent(int plane, int row, int col);

      std::map<PixLoc_t, PixCont_t> pixMap;
      std::string adc_tree;
      void set_adc_treename( std::string name ) { adc_tree = name; };

  };

}
}

#endif
