#ifndef __INFILLSPARSIFYIMAGE_CXX__
#define __INFILLSPARSIFYIMAGE_CXX__

#include <sstream>
#include <random>

#include "InfillSparsifyImage.h"

//larcv
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/EventSparseImage.h"
#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/DataFormat/IOManager.h"

//root
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"

namespace ublarcvapp {

  static InfillSparsifyImageProcessFactory __global_InfillSparsifyImageProcessFactory__;

  InfillSparsifyImage::InfillSparsifyImage(const std::string name)
    : larcv::ProcessBase(name)
  {}

  void InfillSparsifyImage::configure(const larcv::PSet& cfg)
  {

    _verbosity_                   = cfg.get<int>("Verbosity");

    _input_adc_producer           = cfg.get<std::string>("InputADCProducer");
    _input_adcmasked_producer     = cfg.get<std::string>("InputADCMaskedProducer");
    _input_labels_producer        = cfg.get<std::string>("InputLabelsProducer");

    _output_adc_producer          = cfg.get<std::string>("OutputADCProducer");
    _output_adcmasked_producer    = cfg.get<std::string>("OutputADCMaskedProducer");
    // _output_labels_producer       = cfg.get<std::string>("OutputLabelsProducer");

    _output_filename              = cfg.get<std::string>("OutputFilename");

    // output file
    foutIO = new larcv::IOManager( larcv::IOManager::kWRITE );
    foutIO->set_out_file( _output_filename );
    foutIO->initialize();

  }

  void InfillSparsifyImage::initialize()
  {}

  bool InfillSparsifyImage::process(larcv::IOManager& mgr)
  {
    // ---------------------------------------------------------------
    // get data

    // input ADC "truth"
    auto ev_in_adc  = (larcv::EventImage2D*)(mgr.get_data(larcv::kProductImage2D, _input_adc_producer));
    if (!ev_in_adc) {
      LARCV_CRITICAL() << "No Input ADC Image2D found with a name: " << _input_adc_producer << std::endl;
      throw larcv::larbys();
    }
    std::vector< larcv::Image2D > img_adc_v = ev_in_adc->Image2DArray();

    // input ADC masked
    auto ev_in_adc_masked  = (larcv::EventImage2D*)(mgr.get_data(larcv::kProductImage2D, _input_adcmasked_producer));
    if (!ev_in_adc_masked) {
      LARCV_CRITICAL() << "No Input ADCMasked Image2D found with a name: " << _input_adcmasked_producer << std::endl;
      throw larcv::larbys();
    }
    std::vector< larcv::Image2D > img_adcmasked_v = ev_in_adc_masked->Image2DArray();

    // input labels image
    auto ev_in_labels  = (larcv::EventImage2D*)(mgr.get_data(larcv::kProductImage2D, _input_labels_producer));
    if (!ev_in_labels) {
      LARCV_CRITICAL() << "No Input Labels Image2D found with a name: " << _input_labels_producer << std::endl;
      throw larcv::larbys();
    }
    std::vector< larcv::Image2D > img_labels_v = ev_in_labels->Image2DArray();

    // ----------------------------------------------------------------

    // Output containers
    larcv::EventSparseImage* ev_out_adc  = (larcv::EventSparseImage*)foutIO->get_data(larcv::kProductSparseImage, _output_adc_producer);
    larcv::EventSparseImage* ev_out_adcmasked  = (larcv::EventSparseImage*)foutIO->get_data(larcv::kProductSparseImage,_output_adcmasked_producer);

    ev_out_adc->clear();
    ev_out_adcmasked->clear();

    // ----------------------------------------------------------------

    // Copies to work with
    std::vector<larcv::Image2D> adc_image_v = ev_in_adc->Image2DArray();
    std::vector<larcv::Image2D> adcmasked_image_v = ev_in_adc_masked->Image2DArray();
    std::vector<larcv::Image2D> labels_image_v = ev_in_labels->Image2DArray();

    // ----------------------------------------------------------------

    // Run, subrun, event
    int run    = ev_in_adc->run();
    int subrun = ev_in_adc->subrun();
    int event  = ev_in_adc->event();

    std::vector<float> threshold_v (1,10.0);

    larcv::SparseImage adc_sparse_tensor;
    larcv::SparseImage adcmasked_sparse_tensor;

    for(int i =0; i<adc_image_v.size(); i++){
      adc_sparse_tensor = larcv::SparseImage(adc_image_v[i],
                                            labels_image_v[i],
                                            threshold_v);
      adcmasked_sparse_tensor = larcv::SparseImage(adcmasked_image_v[i],
                                            labels_image_v[i],
                                            threshold_v);

      ev_out_adcmasked->Append( adcmasked_sparse_tensor );
      ev_out_adc->Append( adc_sparse_tensor );
    }


  	foutIO->set_id( run, subrun,100*event);
  	foutIO->save_entry();

    return true;
  }



// =============================================================================
  void InfillSparsifyImage::finalize()
  {
    foutIO->finalize();
  }

}
#endif
