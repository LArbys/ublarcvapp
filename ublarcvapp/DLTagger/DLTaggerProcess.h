#ifndef __DLTAGGER_PROCESS_H__
#define __DLTAGGER_PROCESS_H__

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

#include "ublarcvapp/LArliteHandler/LArliteManager.h"
#include "DLTagger.h"

namespace ublarcvapp {
namespace dltagger {

  class DLTaggerProcess : public larcv::ProcessBase {
  public:

    DLTaggerProcess(std::string instance_name )
      : larcv::ProcessBase(instance_name),
      _larlite_io(nullptr),
      _ana_tree(nullptr)
    {};
    virtual ~DLTaggerProcess() {};

    // require functions
    void configure(const larcv::PSet& );
    void initialize();
    bool process( larcv::IOManager& mgr );
    void finalize();

    void add_larlite_infile( std::string inputfile ) { _larlite_files.push_back(inputfile); };

  protected:

    DLTagger m_tagger;
    ublarcvapp::LArliteManager* _larlite_io;

    // processor configuration parameters
    std::vector<std::string> _larlite_files; //< list of input larlite files

    // producer and file names
    std::string _input_adc_producer;      //< name of producer for input ADC images
    std::string _input_chstatus_producer; //< name of producer for ChStatus
    std::string _input_opflash_producer;  //< name of producer for in-beam window opflash
    std::string _input_ophit_producer;    //< name of producer for in-beam window ophit
    std::string _input_mask_producer;     //< name of producer containing Mask R-CNN output masks
    std::string _output_tagged_image;     //< name to give event container storing whole-view images tagged with clusters labeled as cosmic
    std::string _output_notcosmic_image;  //< name to give event container storing whole-view images tagged with clusters labeled as not-cosmi
    std::string _output_cosmic_clusters;    //< name to give event container storing pixels from found tracks
    std::string _output_notcosmic_clusters; //< name to give event container storing pixels from found tracksa    
    std::string _output_larlite_file;     //< name of output larlite file
    std::string _output_tracks;           //< name to give event container storing found 3D tracks
    std::string _output_croi;             //< name to give container storing collection of CROIs
    std::string _output_croi_merged;      //< name to give container storing merged CROI

    // true-reco matching
    bool _has_mcinstance_img;             //< if true, expect to have mc instance image for truth-reco matching
    std::string _input_instance_image;    //< name to input producer containing MC Instance Image

    // opflash window
    float _inbeam_win_start_tick;
    float _inbeam_win_end_tick;

    void setupAnaTree();
    void clearAnaVars();
    void fillAnaVars();
    int _num_clusters;
    std::vector< int   > _astar_complete;
    std::vector< float > _dwall_outermost;
    std::vector< float > _dwall_innermost;
    std::vector< float > _dtick_outoftime;
    std::vector< float > _frac_in_croi_plane0;
    std::vector< float > _frac_in_croi_plane1;
    std::vector< float > _frac_in_croi_plane2;
    std::vector< float > _frac_in_croi_total;
    std::vector< float > _nufrac_plane0;
    std::vector< float > _nufrac_plane1;
    std::vector< float > _nufrac_plane2;
    std::vector< float > _nufrac_total;
    std::vector< float > _numpixels_plane0;
    std::vector< float > _numpixels_plane1;
    std::vector< float > _numpixels_plane2;
    std::vector< std::vector<float> > _outermost_endpt_v;
    std::vector< std::vector<float> > _innermost_endpt_v;    
    TTree* _ana_tree;
    
  };

  class DLTaggerProcessFactory : public larcv::ProcessFactoryBase {
  public:
    DLTaggerProcessFactory() { larcv::ProcessFactory::get().add_factory("DLTaggerProcess",this); };
    ~DLTaggerProcessFactory() {};
    larcv::ProcessBase* create(const std::string instance_name) { return new DLTaggerProcess(instance_name); };
  };
  
}
}

#endif
