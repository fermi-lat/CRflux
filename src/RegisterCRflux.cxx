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

#include "CrLocation.h"

#include "CLHEP/Random/Random.h"

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
    static RemoteSpectrumFactory<CrProton> CRfactory2(fsvc);
    static RemoteSpectrumFactory<CrAlpha> CRfactory3(fsvc);
    static RemoteSpectrumFactory<CrElectron> CRfactory4(fsvc);
    static RemoteSpectrumFactory<CrPositron> CRfactory5(fsvc);
    static RemoteSpectrumFactory<CrGamma> CRfactory6(fsvc);

	// CRflux needs to use the same random engine as FluxSvc
	HepRandomEngine* engine = fsvc->getRandomEngine();
	HepRandom::setTheEngine(engine);

   // Get the initial location from FluxSvc and store in the CrLocation singleton
   CrLocation::Instance()->setLat(fsvc->location().first);
   CrLocation::Instance()->setLon(fsvc->location().second);
   CrLocation::Instance()->setFluxSvc(fsvc);

    return StatusCode::SUCCESS;
} 





