
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/Service.h"

#include "FluxSvc/IRegisterSource.h"
#include "ICRfluxSvc.h"
#include <iostream>


class CRfluxSvc : public Service, virtual public ICRfluxSvc
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
    IToolSvc *m_toolSvc;
};



DECLARE_SERVICE_FACTORY(CRfluxSvc);


//------------------------------------------------------------------------------
/// Service parameters which can be set at run time must be declared.
/// This should be done in the constructor.

CRfluxSvc::CRfluxSvc(const std::string& name, 
                               ISvcLocator* pSvcLocator) 
                              : Service(name, pSvcLocator), m_toolSvc(0)
{   
    return; 
}

StatusCode CRfluxSvc::initialize()
{
   MsgStream  log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    Service::initialize();
    setProperties();

    if ((sc=service("ToolSvc",m_toolSvc, true)).isFailure() ){
        log << MSG::ERROR << "Couldn't fine ToolSvc" << endreq;
        return sc;
    }

    IRegisterSource *registerTool;
    m_toolSvc->retrieveTool("RegisterCRflux", registerTool);

    return StatusCode::SUCCESS;

}


StatusCode CRfluxSvc::finalize()
{
    return StatusCode::SUCCESS;
}


StatusCode CRfluxSvc::queryInterface(const InterfaceID& riid, void **ppvIF)
{
    if (ICRfluxSvc::interfaceID() == riid) {
        *ppvIF = dynamic_cast<ICRfluxSvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}

// access the type of this service
const InterfaceID&  CRfluxSvc::type () const {
    return ICRfluxSvc::interfaceID();
}


			    
