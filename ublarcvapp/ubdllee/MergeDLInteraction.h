#ifndef __MERGE_DL_INTERACTION_H__
#define __MERGE_DL_INTERACTION_H__

#include "TTree.h"

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"

#include "DLInteraction.h"

namespace ublarcvapp {
namespace ubdllee {

  class MergeDLInteraction : public larcv::ProcessBase {

  public:

    MergeDLInteraction( std::string instance_name )
      : larcv::ProcessBase( instance_name ),
      _treename(""),
      _tree(nullptr),
      _event_interaction_v(nullptr)
      {
      };

    virtual ~MergeDLInteraction() {};


    virtual void configure(const larcv::PSet& pset );
    virtual void initialize();
    virtual bool process(larcv::IOManager& mgr);
    virtual void finalize();

  protected:

    std::string _treename;
    TTree* _tree;
    std::vector< DLInteraction >* _event_interaction_v;
    
  };

  class MergeDLInteractionFactory : public larcv::ProcessFactoryBase {
  public:
    MergeDLInteractionFactory() { larcv::ProcessFactory::get().add_factory("MergeDLInteraction",this); };
    ~MergeDLInteractionFactory() {};
    larcv::ProcessBase* create(const std::string instance_name) { return new MergeDLInteraction(instance_name); };
  };
  
  
}
}

#endif
