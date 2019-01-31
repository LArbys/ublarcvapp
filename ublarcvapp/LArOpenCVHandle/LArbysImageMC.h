#ifndef __LARBYSIMAGEMC_H__
#define __LARBYSIMAGEMC_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

#include <vector>
#include "TTree.h"
#include "LArUtil/PxUtils.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "Geo2D/Core/HalfLine.h"
#include "Geo2D/Core/Line.h"
#include "larcv/core/DataFormat/EventROI.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class LArbysImageMC ... these comments are used to generate
     doxygen documentation!
  */
  class LArbysImageMC : public ProcessBase {

  public:

    /// Default constructor
    LArbysImageMC(const std::string name="LArbysImageMC");

    /// Default destructor
    ~LArbysImageMC(){}

    void configure(const PSet&);
    void initialize();
    bool process(IOManager& mgr);
    void finalize();
    void Clear();

  protected:

    std::string _producer_roi;
    std::string _producer_image2d;
    bool _neutrino_present;

    TTree* _mc_tree;

    /// Event ID
    int _run;
    int _subrun;
    int _event;
    int _entry;

    std::string _rse_producer;

    /// Primary Particle Info
    int _parent_pdg;//primary particle pdg

    double _energy_deposit;
    double _energy_init;

    double _parent_x;
    double _parent_y;
    double _parent_z;
    double _parent_t;
    double _parent_px;
    double _parent_py;
    double _parent_pz;

    int _current_type;
    int _interaction_type;
    float _length_2d;

    int _nprimary;
    int _ntotal;

    int _nproton;
    int _nlepton;
    int _nshower;
    int _nmeson;
    int _nneutron;

    int _hi_lep_pdg;
    double _hi_lep_e;

    double _dep_sum_lepton;
    double _ke_sum_lepton;

    double _dep_sum_proton;
    double _ke_sum_proton;

    double _dep_sum_meson;
    double _ke_sum_meson;

    double _dep_sum_shower;
    double _ke_sum_shower;

    double _dep_sum_neutron;
    double _ke_sum_neutron;

    geo2d::Vector<float> _start; //2d start point
    geo2d::Vector<float> _dir;   //2d dir

    std::vector<int>   _daughter_pdg_v;
    std::vector<int>   _daughter_trackid_v;
    std::vector<int>   _daughter_parenttrackid_v;

    std::vector<double> _daughter_energyinit_v;
    std::vector<double> _daughter_energydep_v;

    std::vector<std::vector<double> > _daughter_length_vv;
    std::vector<std::vector<double> > _daughter_2dstartx_vv;
    std::vector<std::vector<double> > _daughter_2dstarty_vv;
    std::vector<std::vector<double> > _daughter_2dendx_vv;
    std::vector<std::vector<double> > _daughter_2dendy_vv;
    std::vector<std::vector<double> > _daughter_2dcosangle_vv;

    std::vector<double> _daughterPx_v;
    std::vector<double> _daughterPy_v;
    std::vector<double> _daughterPz_v;
    std::vector<double> _daughterX_v;
    std::vector<double> _daughterY_v;
    std::vector<double> _daughterZ_v;
    std::vector<double> _daughter_length3d_v;

    double _true_proton_end_pt_x;
    double _true_proton_end_pt_y;
    double _true_proton_end_pt_z;

    double _proton_1st_pt_x;
    double _proton_1st_pt_y;
    double _proton_1st_pt_z;

    double _proton_last_pt_x;
    double _proton_last_pt_y;
    double _proton_last_pt_z;

    double _true_lepton_end_pt_x;
    double _true_lepton_end_pt_y;
    double _true_lepton_end_pt_z;

    double _lepton_1st_pt_x;
    double _lepton_1st_pt_y;
    double _lepton_1st_pt_z;

    double _lepton_last_pt_x;
    double _lepton_last_pt_y;
    double _lepton_last_pt_z;

  public:
    /// 2D Vertex Info
    std::vector<double> _vtx_2d_w_v;
    std::vector<double> _vtx_2d_t_v;


  protected:
    /// LARCV Image2D data
    std::vector<larcv::Image2D> _image_v;
    ImageMeta _meta;

    /// Filtering
    float _min_nu_dep_e;
    float _max_nu_dep_e;
    float _min_nu_init_e;
    float _max_nu_init_e;

    int _min_n_proton;
    int _min_n_neutron;
    int _min_n_lepton;
    int _min_n_meson;
    int _min_n_shower;

    float _min_proton_init_e;
    float _min_proton_dep;
    float _max_proton_dep;
    float _min_lepton_init_e;

    bool _check_vis;

    bool _do_not_reco;

    bool _selected;

    bool _select_signal;
    bool _select_background;

    struct aparticle{
      int trackid;
      int parenttrackid;
      float depeng;
      bool primary;

      bool daughterof (const aparticle& particle) const
      { return (this->parenttrackid == particle.trackid); }
    };


    struct entry_info{
      int run;
      int subrun;
      int event;
    };

    entry_info _entry_info;
    bool _is_signal;


  };

  /**
     \class larcv::LArbysImageMCFactory
     \brief A concrete factory class for larcv::LArbysImageMC
  */
  class LArbysImageMCProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    LArbysImageMCProcessFactory() { ProcessFactory::get().add_factory("LArbysImageMC",this); }
    /// dtor
    ~LArbysImageMCProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new LArbysImageMC(instance_name); }

  };

}

#endif
/** @} */ // end of doxygen group
