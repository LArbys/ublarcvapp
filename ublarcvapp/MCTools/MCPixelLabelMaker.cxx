#include "MCPixelLabelMaker.h"

#include "larlite/LArUtil/TimeService.h"
#include "larlite/LArUtil/Geometry.h"
#include "larlite/LArUtil/LArProperties.h"
#include "larlite/DataFormat/simch.h"

#include "larcv/core/DataFormat/EventImage2D.h"

#include "ublarcvapp/MCTools/MCPos2ImageUtils.h"

#ifdef HAVE_HIGHFIVE
#include <highfive/H5Easy.hpp>
#endif

namespace ublarcvapp {
namespace mctools {

  MCPixelLabelMaker::~MCPixelLabelMaker()
  {
      if ( psce ) {
          delete psce;
          psce = nullptr;
      }

      if ( preverse_sce ) {
          delete preverse_sce;
          preverse_sce = nullptr;
      }
  }

  void MCPixelLabelMaker::process( 
    larlite::storage_manager& ioll, 
    larcv::IOManager& iolcv,
    std::string image2d_tree_name )
  {
    
    // moving real position to apparent position
    if ( psce==nullptr ) {
      psce = new larutil::SpaceChargeMicroBooNE(larutil::SpaceChargeMicroBooNE::kMCC9_Forward);
    }

    ublarcvapp::mctools::MCParticleGraph mcpg;
    mcpg.buildgraph(ioll);

    make_truthlabels_fromsimch(image2d_tree_name,ioll,iolcv,mcpg,psce);

  }

  void MCPixelLabelMaker::make_truthlabels_fromsimch(
      std::string image2d_tree_name,
      larlite::storage_manager& ioll, 
      larcv::IOManager& iolcv,
      ublarcvapp::mctools::MCParticleGraph& mcpg,
      larutil::SpaceChargeMicroBooNE* psce )
  {

    LARCV_INFO() << "start" << std::endl;

    // utility to go from simulated electronics TDC to 
    // ticks (tdcs after readout trigger)
    const larutil::TimeService* timeservice = larutil::TimeService::GetME();

    // geometry
    const larutil::Geometry* geom = larutil::Geometry::GetME();


    // drift velocity
    float driftv = larutil::LArProperties::GetME()->DriftVelocity();

    // get the simch product we need
    larlite::event_simch* ev_simch = 
      (larlite::event_simch*)ioll.get_data(larlite::data::kSimChannel,source);

    // get images
    larcv::EventImage2D* ev_img = 
      (larcv::EventImage2D*)iolcv.get_data(larcv::kProductImage2D,image2d_tree_name);

    auto const& img_v = ev_img->as_vector();
    int nplanes = (int)img_v.size();

    auto const& meta0 = img_v.at(0).meta();

    // Clear the triplet info container
    _pixels_v.clear();

    // loop over simch information, making TripletLabels_t
    size_t nsimch = ev_simch->size();
    size_t ide_w_no_t0 = 0;
    size_t ide_outofimg = 0;
    size_t ide_w_badwire = 0;
    size_t num_ide_used = 0;

    for (size_t isimch=0; isimch<nsimch; isimch++) {
        auto& simch = ev_simch->at(isimch);
        auto chid = simch.Channel();

        larlite::geo::WireID wireid = geom->ChannelToWireID(chid);
        int plane = wireid.Plane;

        //std::cout << "(" << isimch << ") chid=" << chid << " plane=" << plane << std::endl;

        auto& idcmap = simch.TDCIDEMap();
        for ( auto it=idcmap.begin(); it!=idcmap.end(); it++ ) {
            long tdc = it->first;
            int tick = int(timeservice->TPCTDC2Tick(tdc));
            size_t nide = it->second.size();
            for (auto& ide : it->second ) {

                std::vector<double> pos = { ide.x, ide.y, ide.z };
                long tid = ide.trackID;
                long xtid = (tid>=0) ? tid : -tid;
                double edep = ide.energy;

                // get wire coordinates for these positions
                std::vector<int> wire_v(nplanes,0);
                bool bad_wire = false;
                for (int iplane=0; iplane<nplanes; iplane++) {
                    try {
                        UInt_t wireid = geom->NearestWire( pos, iplane );
                        wire_v[iplane] = wireid;
                    }
                    catch (...){
                        bad_wire = true;
                    }
                }
                if ( bad_wire ) {
                    ide_w_badwire++;
                    continue;
                }

                // replace low-energy shower trackid label with mother of the shower
                long mtid = mcpg.getShowerMotherID( tid );
                if (mtid>0) {
                    xtid = mtid;
                }

                auto pnode_t = mcpg.findTrackID( xtid );
                if ( pnode_t==nullptr ) {
                    // don't have an alternative for this right now
                    ide_w_no_t0++;
                    continue;
                }

                long aid = mcpg.getAncestorID( xtid );
                int pid  = pnode_t->pid;
                int origin = pnode_t->origin;

                double t0 = pnode_t->start.at(3); // in ns


                std::vector<double> pos_sce(3,0);
                std::vector<double> imgpos(4,0); // (U,V,Y,tick)
                int row = -1;
                std::array<int,4> imgindex = {-1,-1,-1,-1};

                if ( source=="largeant" ) {

                    bool applied = false;
                    pos_sce = psce->ApplySpaceChargeEffect( pos[0], pos[1], pos[2], applied );

                    // get (u,v,y,tick)
                    std::vector<float> imgpos = 
                        ublarcvapp::mctools::MCPos2ImageUtils::Get()->truepos_to_imagepos( pos[0],
                            pos[1],
                            pos[2],
                            t0,
                            true );
                    pos[0] = (tick-3200)*0.5*driftv;
                    tick = (int)imgpos[3];
                    if ( tick<(int)meta0.min_y() || tick>=(int)meta0.max_y() ) {
                        ide_outofimg++;
                        continue;
                    }

                    row = meta0.row( imgpos[3] );

                    imgindex = { 
                        (int)imgpos[0], 
                        (int)imgpos[1], 
                        (int)imgpos[2], 
                        row };
                }
                else if ( source=="driftWC:simpleSC:Detsim") {
                    if ( preverse_sce == nullptr ) {
                        preverse_sce = new larutil::SpaceChargeMicroBooNE(larutil::SpaceChargeMicroBooNE::kMCC9_Backward);
                    }
                    // the (y,z) is post spacecharge effect and x is drifted to the wireplanes
                    float x_t0_offset = (t0*1.0e-3)*driftv;
                    x_t0_offset = 0.0;
                    pos_sce[0] = float(tick-3200)*0.5*driftv - x_t0_offset;
                    pos_sce[1] = pos[1];
                    pos_sce[2] = pos[2];
                    pos[0] = pos_sce[0];

                    for (int p=0; p<3; p++) {
                        imgpos[p]   = wire_v[p];
                        imgindex[p] = wire_v[p];
                    }
                    imgpos[3] = tick;
                    if ( tick<(int)meta0.min_y() || tick>=(int)meta0.max_y() ) {
                        ide_outofimg++;
                        continue;
                    }
                    row = meta0.row( tick );
                    imgindex[3] = row;
                }
                else {
                    std::stringstream errmsg;
                    errmsg << "unrecognized simch source: " << source << std::endl;
                    throw std::runtime_error(errmsg.str());
                }

                auto it_index = _pixels_v._imgcoord_to_tripindex.find( imgindex );
                if ( it_index==_pixels_v._imgcoord_to_tripindex.end() ) {
                    MCPixelLabels trip;
                    trip.index = (long)_pixels_v._triplets_v.size();
                    for (int i=0; i<3; i++)
                        trip.edep[i] = 0.0;
                    trip.imgcoord[0] = imgindex[0];
                    trip.imgcoord[1] = imgindex[1];
                    trip.imgcoord[2] = imgindex[2];
                    trip.imgcoord[3] = (int)tick;
                    trip.imgcoord[4] = row;
                    trip.pos[0] = pos[0];
                    trip.pos[1] = pos[1];
                    trip.pos[2] = pos[2];
                    trip.pos_reco[0] = pos_sce[0];
                    trip.pos_reco[1] = pos_sce[1];
                    trip.pos_reco[2] = pos_sce[2];

                    // make index map
                    _pixels_v._imgcoord_to_tripindex[imgindex] = trip.index;

                    // we also index stuff
                    for (int iu=-dwire; iu<=dwire; iu++) {
                    for (int iv=-dwire; iv<=dwire; iv++) {
                    for (int iy=-dwire; iy<=dwire; iy++) {
                    for (int ir=-drow; ir<=drow; ir++) {
                        std::array<int,4> modindex = imgindex;
                        modindex[0] += iu;
                        modindex[1] += iv;
                        modindex[2] += iy;
                        modindex[3] += ir;
                        auto it_mod = _pixels_v._imgcoord_to_tripindex.find( modindex );
                        if ( it_mod==_pixels_v._imgcoord_to_tripindex.end()) {
                            _pixels_v._imgcoord_to_tripindex[modindex] = trip.index;
                        }
                        else {
                            // found a previous record
                            int new_dindex = std::abs(iu)+std::abs(iv)+std::abs(iy)+std::abs(ir);

                            if ( new_dindex!=0 ) { 
                              auto& oldtrip = _pixels_v._triplets_v.at( it_mod->second );
                              int old_dindex = 0;
                              for (int ii=0; ii<3; ii++) {
                                  old_dindex += std::abs( imgindex[ii]-oldtrip.imgcoord[ii] );
                              }
                              old_dindex += std::abs( imgindex[3]-oldtrip.imgcoord[4] );
  
                              if ( new_dindex<=old_dindex ) {
                                  // replace index
                                  _pixels_v._imgcoord_to_tripindex[modindex] = trip.index;
                              }
                            }
                        }
                    }
                    }
                    }
                    }

                    _pixels_v._triplets_v.emplace_back( std::move(trip) );     
                    it_index = _pixels_v._imgcoord_to_tripindex.find( imgindex );
                }

                auto& tripinfo = _pixels_v._triplets_v.at(it_index->second);
                tripinfo.edep[plane] += edep;
                tripinfo.trackids.insert(xtid);
                tripinfo.aids.insert(aid);
                tripinfo.pids.insert(pid);
                tripinfo.origin.insert(origin);
                num_ide_used++;

            }

        }
    }

    LARCV_INFO() << "Number of Triplets Created: " << _pixels_v._triplets_v.size() << std::endl;
    LARCV_INFO() << "  IDEs with no track ID match and t0: " << ide_w_no_t0 << std::endl;
    LARCV_INFO() << "  IDEs out-of-image: " << ide_outofimg << std::endl;
    LARCV_INFO() << "  IDEs with no nearby-wire: " << ide_w_badwire << std::endl;
    LARCV_INFO() << "  IDEs used: " << num_ide_used << std::endl; 


  }

  void MCPixelLabelMaker::export_as_hdf( std::string hdf_outfile )
  {

#ifdef HAVE_HIGHFIVE

    LARCV_INFO() << "export to " << hdf_outfile << std::endl;

    HighFive::File file(hdf_outfile, HighFive::File::Overwrite);

    file.createGroup("/mcpixel_labels");

    // export different arrays for export
    int ntriplets = _pixels_v._triplets_v.size();

    std::vector<float> pos_x(ntriplets,0);
    std::vector<float> pos_y(ntriplets,0);
    std::vector<float> pos_z(ntriplets,0);

    std::vector<float> pos_x_reco(ntriplets,0);
    std::vector<float> pos_y_reco(ntriplets,0);
    std::vector<float> pos_z_reco(ntriplets,0);

    std::vector< std::array<double,3> > edep(ntriplets);
    std::vector<long>  trackid(ntriplets,0);
    std::vector<int>   pid(ntriplets,0);
    std::vector<int>   aid(ntriplets,0);
    std::vector<int>   origin(ntriplets,0);
    std::vector<int>   uwire(ntriplets,0);
    std::vector<int>   vwire(ntriplets,0);
    std::vector<int>   ywire(ntriplets,0);
    std::vector<int>   tick(ntriplets,0);
    std::vector<int>   row(ntriplets,0);

    for (auto const& triplet : _pixels_v._triplets_v ) {
        long idx = triplet.index;

        pos_x[idx] = triplet.pos[0];
        pos_y[idx] = triplet.pos[1];
        pos_z[idx] = triplet.pos[2];

        pos_x_reco[idx] = triplet.pos_reco[0];
        pos_y_reco[idx] = triplet.pos_reco[1];
        pos_z_reco[idx] = triplet.pos_reco[2];

        edep[idx] = std::array<double,3>{0,0,0};
        for (int i=0; i<3; i++)
          edep[idx][i] = triplet.edep[i];

        for ( auto& tid : triplet.trackids ) {
            trackid[idx] = tid;
            if (trackid[idx]!=-1)
                break;
        }

        for ( auto& xpid : triplet.pids ) {
            pid[idx]     = xpid;
            if (pid[idx]!=-1)
                break;
        }

        for ( auto& xaid : triplet.aids ) {
            aid[idx]     = xaid;
            if (aid[idx]!=-1)
                break;
        }

        for ( auto& xorigin : triplet.origin ) {
            origin[idx]  = xorigin;
            if (origin[idx]!=-1)
                break;
        }

        uwire[idx]   = triplet.imgcoord[0];
        vwire[idx]   = triplet.imgcoord[1];
        ywire[idx]   = triplet.imgcoord[2];
        tick[idx]    = triplet.imgcoord[4];
        row[idx]     = triplet.imgcoord[3];
    }

    H5Easy::dump( file, "/mcpixel_labels/pos_x", pos_x);
    H5Easy::dump( file, "/mcpixel_labels/pos_y", pos_y);
    H5Easy::dump( file, "/mcpixel_labels/pos_z", pos_z);

    H5Easy::dump( file, "/mcpixel_labels/pos_x_reco", pos_x_reco);
    H5Easy::dump( file, "/mcpixel_labels/pos_y_reco", pos_y_reco);
    H5Easy::dump( file, "/mcpixel_labels/pos_z_reco", pos_z_reco);

    H5Easy::dump( file, "/mcpixel_labels/edep",    edep);
    H5Easy::dump( file, "/mcpixel_labels/trackid", trackid);
    H5Easy::dump( file, "/mcpixel_labels/pid",     pid);
    H5Easy::dump( file, "/mcpixel_labels/aid",     aid);
    H5Easy::dump( file, "/mcpixel_labels/origin",  origin);
    H5Easy::dump( file, "/mcpixel_labels/uwire",   uwire);
    H5Easy::dump( file, "/mcpixel_labels/vwire",   vwire);
    H5Easy::dump( file, "/mcpixel_labels/ywire",   ywire);
    H5Easy::dump( file, "/mcpixel_labels/tick",    tick);
    H5Easy::dump( file, "/mcpixel_labels/row",     row);

    file.flush();

#else
    LARCV_CRITICAL() << "Compiled without HDF5 support." << std::endl;
#endif

  }


}
}