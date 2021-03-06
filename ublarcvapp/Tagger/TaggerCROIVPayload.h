#ifndef __TAGGER_CROI_VPAYLOAD_H__
#define __TAGGER_CROI_VPAYLOAD_H__

#include <string>
#include <vector>

namespace ublarcvapp {

  namespace tagger {

    class TaggerCROIVPayload {
      
    protected:
    TaggerCROIVPayload() : m_payloadname("undefined") {};
    TaggerCROIVPayload(std::string payloadname)
      : m_payloadname(payloadname) {};
      virtual ~TaggerCROIVPayload() {};
      
    public:
      std::string getPayLoadName() { return m_payloadname; }
      bool isPayloadType( std::string payloadname ) {
        if ( payloadname==m_payloadname ) 
          return true;
        return false;
      };
      
      virtual void saveSpace() = 0;
      
    protected:
      
      std::string m_payloadname;
      
    };
    
    typedef std::vector<TaggerCROIVPayload> TaggerCROIPayloadList_t;
    
  }
}

#endif
