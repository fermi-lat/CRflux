#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "FluxSvc/IRegisterSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/IFluxSvc.h"
// CR includes
#include "CrExample.h"
#include "CrProton.hh"
#include "CrAlpha.hh"
#include "CrElectron.hh"
#include "CrPositron.hh"
#include "CrGamma.hh"

/** @class RegisterCRflux
 *  @brief Register the various CRflux sources
 *  
 *   @author Toby Burnett
 *   $Header$
 */
class RegisterCRflux : public AlgTool, virtual public IRegisterSource {
public:
    
    RegisterCRflux( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~RegisterCRflux() { }
    
    /// implement to define sources: will be called from FluxSvc
    StatusCode registerMe(IFluxSvc* );
    
};


// Static factory for instantiation of algtool objects
static ToolFactory<RegisterCRflux> s_factory;
const IToolFactory& RegisterCRfluxFactory = s_factory;

// Standard Constructor
RegisterCRflux::RegisterCRflux(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : AlgTool( type, name, parent ) {
    
    // Declare additional interface
    declareInterface<IRegisterSource>(this);
        
}


StatusCode RegisterCRflux::registerMe(IFluxSvc* fsvc) 
{
    MsgStream  log(msgSvc(), name());
    log << MSG::INFO << "Register CRflux Spectrum(s)..." << endreq;

    //declare the factories here:
    static RemoteSpectrumFactory<CrExample> CRfactory(fsvc);
    static RemoteSpectrumFactory<CrProton> CRfactory3(fsvc);
    static RemoteSpectrumFactory<CrAlpha> CRfactory4(fsvc);
    static RemoteSpectrumFactory<CrElectron> CRfactory5(fsvc);
    static RemoteSpectrumFactory<CrPositron> CRfactory6(fsvc);
    static RemoteSpectrumFactory<CrGamma> CRfactory2(fsvc);

    return StatusCode::SUCCESS;
} 

