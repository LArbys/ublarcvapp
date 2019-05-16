/**
 * \file InfillSparsifyImage.h
 *
 * \ingroup UBImageMod
 *
 * \brief Class def header for a class InfillSparsifyImage
 *
 * @author kmason
 *
 * We sparsify data crops for sparse infill network
 *
 */

/** \addtogroup UBImageMod

    @{*/
#ifndef __INFILLSPARSIFYIMAGE_H__
#define __INFILLSPARSIFYIMAGE_H__

#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventChStatus.h"
#include "larcv/core/DataFormat/EventSparseImage.h"

namespace ublarcvapp {

  /**
     \class ProcessBase
     User defined class InfillSparsifyImage ... these comments are used to generate
     doxygen documentation!
  */
  class InfillSparsifyImage : public larcv::ProcessBase {

  public:

    /// Default constructor
    InfillSparsifyImage (const std::string name = "InfillSparsifyImage");

    /// Default destructor
    ~InfillSparsifyImage() {}

    void configure(const larcv::PSet&);

    void initialize();

    bool process(larcv::IOManager& mgr);

    void finalize();

    // ----------------------------------------------------------------------
    // functions

    // ----------------------------------------------------------------------
    // we save data ourselves

  public:

    const larcv::IOManager& getOutputIOMan() const { return *foutIO; };

  protected:

    larcv::IOManager* foutIO;

    // ----------------------------------------------------------------------

  private:

    // config parameters
    std::string _input_adc_producer;
    std::string _input_labels_producer;
    std::string _input_adcmasked_producer;

    std::string _output_adc_producer;
    std::string _output_adcmasked_producer;

    std::string _output_filename;

    int _verbosity_;

  };

  /**
     \class larcv::InfillSparsifyImageFactory
     \brief A concrete factory class for larcv::InfillSparsifyImage
  */
  class InfillSparsifyImageProcessFactory : public larcv::ProcessFactoryBase {
  public:
    /// ctor
    InfillSparsifyImageProcessFactory() { larcv::ProcessFactory::get().add_factory("InfillSparsifyImage", this); }
    /// dtor
    ~InfillSparsifyImageProcessFactory() {}
    /// creation method
    larcv::ProcessBase* create(const std::string instance_name) { return new InfillSparsifyImage(instance_name); }
  };

}

#endif
/** @} */ // end of doxygen group
