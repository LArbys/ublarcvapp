#ifndef __ASTAR_LATTICE_H__
#define __ASTAR_LATTICE_H__

#include <vector>

#include "larcv/core/DataFormat/ImageMeta.h"

#include "AStar3DTypes.h"

namespace ublarcvapp {
namespace reco3d {

  class Lattice : public std::map< A3DPixPos_t, AStar3DNode* > {
  public:
    Lattice( const std::vector<float>& origin, const std::vector<int>& widths, const std::vector<float>& pixel_lengths, 
    	const std::vector<const larcv::ImageMeta* >& meta_v ) {
      m_origin = origin;
      m_widths = widths;
      m_cm_per_pixel = pixel_lengths;
      m_meta_v = meta_v;
    }
    virtual ~Lattice() { cleanup(); };

    A3DPixPos_t getNodePos( const std::vector<float>& pos );
    std::vector<float> getPos( const int u, const int v, const int w);
    std::vector<float> getPos( const A3DPixPos_t& );
    AStar3DNode* getNode( const std::vector<float>& pos );
    AStar3DNode* getNode( const int u, const int v, const int w);
    AStar3DNode* getNode( const A3DPixPos_t& nodeid );
    void getImageCoordinates( const A3DPixPos_t& nodeid, int& row, std::vector<int>& cols, bool& within_image );
  
    const std::vector<int>& widths() const { return m_widths; }
    const std::vector<float>& origin() const { return m_origin; }
    const std::vector<float>& pixel_len() const { return m_cm_per_pixel; }    

    void setStartNodeWithPadding( AStar3DNode* start, const std::vector<int>& start_padding );

  protected:

    std::vector<float> m_origin;
    std::vector<int> m_widths;
    std::vector<float> m_cm_per_pixel;
    std::vector< const larcv::ImageMeta* > m_meta_v;

    // start padding
    AStar3DNode* m_start_node;
    A3DPixPos_t  m_start_nodeid;
    std::vector<int> m_start_padding;
    bool fCheckStartPadding;

    void cleanup(); // destroy nodes allocated on the heap


  };

  
  
}
}

#endif
