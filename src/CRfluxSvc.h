
#ifndef CRFLUXSVC_H
#define CRFLUXSVC_H 1

/** 
 * @class CRfluxSvc
 * 
 * @author Heather Kelly 
 *
 * $Header$
 */

#include "GaudiKernel/Service.h"

class CRfluxSvc : public Service
{
public:
    
    CRfluxSvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~CRfluxSvc() {}
    
    StatusCode initialize();
    StatusCode finalize();
    

    /// return the Interface ID
    static const InterfaceID& interfaceID() {
        return ICRfluxSvc::interfaceID(); 
    }
    /// return the service type
    const InterfaceID& type() const;


    /// query by name
    StatusCode queryInterface(const InterfaceID& riid, void** ppvUnknown);


private:

};

#endif // CRFLUXSVC_H
