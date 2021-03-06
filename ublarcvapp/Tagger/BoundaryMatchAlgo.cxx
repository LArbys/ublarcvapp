#include "BoundaryMatchAlgo.h"

#include <ctime>
#include <iostream>
#include <assert.h>
#include <cstring>

#include "ublarcvapp/UBWireTool/UBWireTool.h"
#include "BoundaryMatchArrays.h"


namespace ublarcvapp {
namespace tagger {

  BoundaryMatchList::~BoundaryMatchList() {
    delete [] uvec;
    delete [] vvec;
    delete [] yvec;
  }

  void BoundaryMatchList::load( BoundaryMatchArrays::Boundary_t t, BoundaryMatchArrays::MatchMode_t mode ) {
    //const clock_t begin_time = clock();
    
    BoundaryMatchArrays match_arrays(mode);
    
    m_combos.resize( match_arrays.nmatches(t) );

    for (int imatch=0; imatch<match_arrays.nmatches(t); imatch++) {
      int u, v, y;
      match_arrays.getMatch( t, imatch, u, v, y );

      // define and store the combo
      float pos[3] = { 0., 0., 0.};
      match_arrays.getPosition( t, imatch, pos[0], pos[1], pos[2] );

      BoundaryCombo combo( u, v, y, pos[0], pos[1], pos[2] );
      m_combos[imatch] = combo;

      // build index set trees
      auto it_u = u_indices.find((idx)u);
      auto it_v = v_indices.find((idx)v);
      auto it_y = y_indices.find((idx)y);
      if ( it_u==u_indices.end() ) {
	u_indices.insert( std::pair< idx, std::set<idx> >(u, std::set<idx>() ) );
      }
      if ( it_v==v_indices.end() ) {
	v_indices.insert( std::pair< idx, std::set<idx> >(v, std::set<idx>() ) );
      }
      if ( it_y==y_indices.end() ) {
	y_indices.insert( std::pair< idx, std::set<idx> >(y, std::set<idx>() ) );
      }
      u_indices[(idx)u].insert((idx)imatch);
      v_indices[(idx)v].insert((idx)imatch);
      y_indices[(idx)y].insert((idx)imatch);
      if ( u>uid_max ) uid_max = u;
      if ( v>vid_max ) vid_max = v;
      if ( y>yid_max ) yid_max = y;
    }
    uvec = new int[m_combos.size()];
    vvec = new int[m_combos.size()];
    yvec = new int[m_combos.size()];

    //float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    //std::cout << "boundary match list loaded " << elapsed_secs << " secs" << std::endl;
    
  }
  
  std::vector<BoundaryCombo> BoundaryMatchList::findCombos( const std::vector<int>& uwires, const std::vector<int>& vwires, const std::vector<int>& ywires, 
							    const std::vector<larcv::Image2D>& badchs, bool use_badchs ) {
    //const clock_t begin_time = clock();

    memset(uvec,0,sizeof(int)*m_combos.size());
    memset(vvec,0,sizeof(int)*m_combos.size());
    memset(yvec,0,sizeof(int)*m_combos.size());

    // get dsfactor
    int dsfactor = (int)badchs[0].meta().pixel_width();
    
    // filter u wires

    // loop over u-wires
    //std::cout << " u: ";
    for (size_t uwid=0; uwid<uwires.size(); uwid++) {
      if ( uwires[uwid]==0 ) continue;	
      idx uidx = (idx)( uwid );
      //std::cout << " " << uidx;
      // copy indices into passu
      for ( auto &iiu : u_indices[ uidx ] ) {
	//passu.insert( iiu );
	uvec[iiu] = 1;
      }
    }
    //std::cout << std::endl;
    
    // loop over v-wires, filter matched indices
    //std::cout << " v: ";
    for (size_t vwid=0; vwid<vwires.size(); vwid++) {
      if ( vwires[vwid]== 0 ) continue;
      idx vidx = (idx)( vwid );
      //std::cout << " " << vidx;
      // copy indices into passv if in passu
      for ( auto &iiv : v_indices[ vidx ] ) {
	vvec[iiv] = 1;
      }
    }
    //std::cout << std::endl;

    // loop over y-wires, filter matched indices
    //std::cout << " y: ";
    for (size_t ywid=0; ywid<ywires.size(); ywid++) {
      if ( ywires[ywid]== 0 ) continue;
      idx yidx = (idx)(ywid);
      //std::cout << " " << yidx;
      for ( auto &iiy : y_indices[ yidx ] ) {
	yvec[iiy] = 1;
      }
    }
    //std::cout << std::endl;

    //float elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    //std::cout << "boundary match searched in " << elapsed_secs << " secs" << std::endl;    

    std::vector< BoundaryCombo > combo_v;
    combo_v.reserve( m_combos.size() );
    int ncombos=0;
    for (size_t idx=0; idx<m_combos.size(); idx++) {
      int votes = uvec[idx]+vvec[idx]+yvec[idx];
      if ( votes>=3 ) {
	combo_v.push_back( m_combos.at(idx) );
	ncombos++;
      }
      else if ( votes==2 && use_badchs ) {
	int badchp = -1;
	int col = -1;
	if ( uvec[idx]==0 ) { 
	  badchp = 0;
	  col = uvec[idx]/dsfactor;
	}
	else if ( vvec[idx]==0 ) {
	  badchp = 1;
	  col = vvec[idx]/dsfactor;
	}
	else if ( yvec[idx]==0 ) {
	  badchp = 2;
	  col = yvec[idx]/dsfactor;
	}
	if ( badchs[badchp].pixel( 0, col ) > 0 ) {
	  combo_v.push_back( m_combos.at(idx) );
	  ncombos++;
	}
      }
    }
    
    //elapsed_secs = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
    //std::cout << "ncombos=" << ncombos << ". boundary match searched/stored in " << elapsed_secs << " secs" << std::endl;    
    
    return combo_v;
  }

  BoundaryMatchAlgo::BoundaryMatchAlgo( BoundaryMatchArrays::MatchMode_t mode ) {
    fMode = mode;
    matchlist[0] = new BoundaryMatchList( BoundaryMatchArrays::kTop, mode );
    matchlist[1] = new BoundaryMatchList( BoundaryMatchArrays::kBottom, mode );
    matchlist[2] = new BoundaryMatchList( BoundaryMatchArrays::kUpstream, mode );
    matchlist[3] = new BoundaryMatchList( BoundaryMatchArrays::kDownstream, mode );
  }

  BoundaryMatchAlgo::~BoundaryMatchAlgo() {
    for (int i=0; i<4; i++) {
      delete matchlist[i];
      matchlist[i] = NULL;
    }
  }

  void BoundaryMatchAlgo::findCombos( const std::vector<int>& uwires, const std::vector<int>& vwires, const std::vector<int>& ywires,
				      const std::vector<larcv::Image2D>& badchs, bool use_badchs,
				      std::vector<  std::vector<BoundaryCombo>  >& boundary_combos ) {
    if ( boundary_combos.size()!=4 ) {
      std::cout << "wrong size of combo types" << std::endl;
      assert(false);
    }
    
    for (size_t i=0; i<4; i++) {
      std::vector<BoundaryCombo> combos = matchlist[i]->findCombos( uwires, vwires, ywires, badchs, use_badchs );
      std::vector<BoundaryCombo>& out_combos = boundary_combos.at(i);
      for (size_t j=0; j<combos.size(); j++) {
	//boundary_combos.push_back( std::move(combos) );
	out_combos.push_back( combos.at(j) );
      }
    }
  }

}  
}
