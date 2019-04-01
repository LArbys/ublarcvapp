/**
 * \file UBCropLArFlow.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class UBCropLArFlow
 *
 * @author twongjirad
 *
 * We use crops from UBSplitDetector to crop out
 * LArFlow images. This does the work we have to do
 * to change the pixel flow variables from one crop into another.
 *
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __UBCROPLARFLOW_H__
#define __UBCROPLARFLOW_H__

#include "larcv/core/Base/PSet.h"
#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/ROI.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

#include "TH2D.h"

namespace ublarcvapp {

  /**
     \class ProcessBase
     User defined class UBCropLArFlow ... these comments are used to generate
     doxygen documentation!
  */
  
  class UBCropLArFlow : public larcv::ProcessBase {

  public:

    /// Default constructor
    UBCropLArFlow(const std::string name = "UBCropLArFlow");

    /// Default destructor
    ~UBCropLArFlow() {}

    void configure(const larcv::PSet&);

    void initialize();

    bool process(larcv::IOManager& mgr);

    void finalize();

    // ----------------------------------------------------------------------
    // functions

    // make the cropped flow and visibilty images
    void make_cropped_flow_images( const int src_plane,
                                   const larcv::ROI& croppedroi,
                                   const std::vector<larcv::Image2D>& wholeadc,
                                   const larcv::EventChStatus& badch,
                                   const std::vector<larcv::Image2D>& wholeflow,
                                   const std::vector<larcv::Image2D>& wholevisi,
                                   const std::vector<float>& thresholds,
                                   std::vector<larcv::Image2D>& cropped_flow,
                                   std::vector<larcv::Image2D>& cropped_visi,
                                   bool has_visi=true );

    // make the cropped flow images only (no visibility)    
    void make_cropped_flow_images( const int src_plane,
                                   const larcv::ROI& croppedroi,
                                   const std::vector<larcv::Image2D>& wholeadc,
                                   const larcv::EventChStatus& badch,                                   
                                   const std::vector<larcv::Image2D>& wholeflow,
                                   const std::vector<float>& thresholds,
                                   std::vector<larcv::Image2D>& cropped_flow );
    
    std::vector<float> check_cropped_images( const int src_plane,
                                             const std::vector<larcv::Image2D>& croppedadc_v,
                                             const larcv::EventChStatus& badch,
                                             const std::vector<float>& thresholds,
                                             const std::vector<larcv::Image2D>& cropped_flow,
                                             const std::vector<larcv::Image2D>& cropped_visi,
                                             std::vector<TH2D>& hvis,
                                             const bool has_visi,
                                             const bool visualize_flow=false );
    
    static void downsample_crops( const std::vector<larcv::Image2D*>& cropped_adc_v,
				  const std::vector<larcv::Image2D>& cropped_flow_v,
				  const std::vector<larcv::Image2D>& cropped_visi_v,
				  std::vector<larcv::Image2D>& downsampled_adc_v,
				  std::vector<larcv::Image2D>& downsampled_flow_v,
				  std::vector<larcv::Image2D>& downsampled_visi_v );

    /* void maxPool( const int row_downsample_factor, const int col_downsample_factor, */
    /*     	  const larcv::Image2D& src_adc, const larcv::Image2D& target_adc, */
    /*     	  const larcv::Image2D& flow, const larcv::Image2D& visi, */
    /*     	  const std::vector<float>& thresholds, */
    /*     	  larcv::Image2D& ds_src_adc, larcv::Image2D& ds_target_adc, */
    /*     	  larcv::Image2D& ds_flow, larcv::Image2D& ds_visi ); */
    
    
    // ----------------------------------------------------------------------
    // we save data ourselves
    
  public:
    
    const larcv::IOManager& getOutputIOMan() const { return *foutIO; };
    
  protected:
    
    larcv::IOManager* foutIO;

    // ----------------------------------------------------------------------    
    
  private:

    // config parameters
    std::string _input_bbox_producer;    
    std::string _input_adc_producer;
    std::string _input_vis_producer;
    std::string _input_flo_producer;
    std::string _input_cropped_producer;
    std::string _output_adc_producer;
    std::string _output_vis_producer;
    std::string _output_flo_producer;
    std::string _output_meta_producer;    
    std::string _output_filename;
    std::vector<float> _thresholds_v;
    int   _max_images;    
    bool  _check_flow;
    bool  _make_check_image;
    bool  _do_maxpool;
    int   _row_downsample_factor;
    int   _col_downsample_factor;
    bool  _limit_overlap;
    bool  _require_min_goodpixels;
    float _max_overlap_fraction;
    int   _verbosity_;
    bool  _save_output;
    bool  _is_mc;
    static bool  _fusevector;

    // The code uses a lot of the std::algorithm to
    // help parallel many of the operations.
    // These are the algorithms
    // -----------------------------------------------
    
    // Apply Offset to all flow values
    struct FlowOffset {
      FlowOffset( float offset ) { _offset = offset; };
      float _offset;
      float operator()( const float& flo_value ) const {return flo_value+_offset; };
	
    };

    // apply flow value, if out of bounds of target crop image, set visibility to -1
    struct ModVisibility {
      ModVisibility( float target_xmin, float target_xmax, float col ) { _xmin = target_xmin; _xmax = target_xmax; _col = col; };
      float _xmin;
      float _xmax;
      float _col;
      float operator()( const float& flo_value, const float& vis_value ) {
	// too much branching?
	// could break up intro separate pieces
	if ( flo_value<=-4000 ) return -1.0; // no-flow value in that pixel
	float targetwire = _col+flo_value;
	if ( targetwire < _xmin || targetwire >=_xmax ) return -1.0; // goes out of bounds
	return vis_value;
      };
    };

    // if adc value below threshold, mask value
    struct MaskBelowThreshold {
      MaskBelowThreshold( float threshold, float maskvalue ) { _threshold = threshold; _maskvalue = maskvalue; };
      float _threshold;
      float _maskvalue;
      float operator()( const float& adc_value, const float& pixvalue ) {
	if ( adc_value<_threshold )
	  return _maskvalue;
	else
	  return pixvalue;
      };
    };

    // get adc value from target image after following flow from source wire to target wire
    struct FollowFlow {
      FollowFlow( const std::vector<float>& data ) : _data(data) {};
      const std::vector<float>& _data; // target adc data. contains value for cols for a given row
      float operator()( const float& flow, const float& col ) {
        // take in flow for src col, add index of col to get target col
        // if out of bounds return 0
        // if now flow
	int targetcol = flow+col;
	if (flow<=-4000 )
	  return -1.0;
	else if ( targetcol<0 || targetcol>=(int)_data.size() )
	  return 0;
	else
	  return _data[targetcol];
      };
    };

    // get adc value from target image after following flow from source wire to target wire
    struct GoodFlow {
      GoodFlow( const std::vector<float>& src_adc,
                const std::vector<float>& target_adc,
                const std::vector<int>& source_status,
                const std::vector<int>& target_status )
      : _src_adc(src_adc),
        _target_adc(target_adc),
        _source_status(source_status),
        _target_status(target_status)
      {};
      const std::vector<float>& _src_adc;       // target adc data. contains value for cols for a given row
      const std::vector<float>& _target_adc;    // target adc data. contains value for cols for a given row
      const std::vector<int>&   _source_status; // target adc data. contains value for cols for a given row
      const std::vector<int>&   _target_status; // target adc data. contains value for cols for a given row      
      float operator()( const float& flow, const float& col ) {
        // take in flow for src col, add index of col to get target col
        // if out of bounds return 0
        // if now flow
        int sourcecol = (int)col;
	int targetcol = (int)flow+col;
	if (flow<=-4000 )
	  return 0;
	if ( targetcol<0 || targetcol>=(int)_target_adc.size() )
	  return 0;
        if ( flow==0 && (_src_adc[sourcecol]<10.0 || _target_adc[targetcol]<10.0) )
          return 0;
        
        // keep: above thresh or bad channel with flow to above thresh or good channel
        if ( (_source_status[sourcecol]==0 || _src_adc[sourcecol]>=10.0 )
             && (_target_status[targetcol]==0 || _target_adc[targetcol]>=10.0 ) )
          return 1;
        return 0;
      };
    };
    
    // if target adc is below threshold, mask the value
    struct MaskFlowToNothing {
      // masks visibility 
      MaskFlowToNothing( float threshold, float maskvalue ) { _threshold = threshold; _maskvalue = maskvalue; };
      float _threshold;
      float _maskvalue;
      float operator()( const float& target_adc_value, const float& pixvalue ) {
	if ( target_adc_value<=-4000 ) return pixvalue; // no flow
	if ( target_adc_value<_threshold )
	  return _maskvalue;
	else
	  return pixvalue;
      };
    };

    // if target adc is below threshold, mask the value
    struct MaskBadFlowPixels {
      // masks flow values
      MaskBadFlowPixels( float target_threshold,
                         float maskvalue ) {
        _threshold = target_threshold;
        _maskvalue = maskvalue;
      };
      float _threshold;
      float _maskvalue;
      float _xmin;
      float _xmax;
      float operator()( const float& target_adc_value, const float& flowvalue ) {
        // input1: adc values for target image
        // input2: flow value
	if ( target_adc_value<_threshold
             || flowvalue<=-4000)
          return _maskvalue;
        else
	  return flowvalue;
      };
    };
    
    
    static int _check_img_counter;
    static const float _NO_FLOW_VALUE_;

    // color palette for visualization
    static void setBWRPalette();
    static void setRainbowPalette();
    static int* _colors;
  };

  /**
     \class larcv::UBCropLArFlowFactory
     \brief A concrete factory class for larcv::UBCropLArFlow
  */
  class UBCropLArFlowProcessFactory : public larcv::ProcessFactoryBase {
  public:
    /// ctor
    UBCropLArFlowProcessFactory() { larcv::ProcessFactory::get().add_factory("UBCropLArFlow", this); }
    /// dtor
    ~UBCropLArFlowProcessFactory() {}
    /// creation method
    larcv::ProcessBase* create(const std::string instance_name) { return new UBCropLArFlow(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group

