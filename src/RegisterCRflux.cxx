#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "FluxSvc/IRegisterSource.h"
#include "FluxSvc/ISpectrumFactory.h"
#include "FluxSvc/IFluxSvc.h"
// CR includes
#include "CrCoordinateTransfer.hh"
#include "CrExample.h"
#include "CrProton.hh"
#include "CrAlpha.hh"
#include "CrElectron.hh"
#include "CrPositron.hh"
#include "CrGamma.hh"
#include "CrHeavyIon.h"
#include "CrHeavyIonVertical.h"
//#include "CrHeavyIonVertZ.h"
//#include "CrHeavyIonZ.h"

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
    StatusCode initialize(); // overload
    
    double m_cutoff;
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
        
    // for testing with a given cutoff
    declareProperty("FixedCutoff", m_cutoff=0);

}


StatusCode RegisterCRflux::registerMe(IFluxSvc* fsvc) 
{
    MsgStream  log(msgSvc(), name());
    log << MSG::INFO << "Register CRflux Spectrum(s)..." << endreq;

    //declare the factories here:
//    static RemoteSpectrumFactory<CrExample> CRfactory(fsvc);
    static RemoteSpectrumFactory<CrProton> CRfactory2(fsvc);
    static RemoteSpectrumFactory<CrAlpha> CRfactory3(fsvc);
    static RemoteSpectrumFactory<CrElectron> CRfactory4(fsvc);
    static RemoteSpectrumFactory<CrPositron> CRfactory5(fsvc);
    static RemoteSpectrumFactory<CrGamma> CRfactory6(fsvc);
    static RemoteSpectrumFactory<CrHeavyIon> CRfactory7(fsvc);
    static RemoteSpectrumFactory<CrHeavyIonVertical> CRfactory8(fsvc);
   // static RemoteSpectrumFactory<CrHeavyIonVertZ> CRfactory9(fsvc);
    //static RemoteSpectrumFactory<CrHeavyIonZ> CRfactory10(fsvc);

	// CRflux needs to use the same random engine as FluxSvc
	CLHEP::HepRandomEngine* engine = fsvc->getRandomEngine();
	CLHEP::HepRandom::setTheEngine(engine);

   // Get the initial location from FluxSvc and store in the CrLocation singleton
   CrLocation::instance()->setFluxSvc(fsvc);
   initialize(); //?

   return StatusCode::SUCCESS;
} 
StatusCode RegisterCRflux::initialize()
{
   //Set the properties
    setProperties();
    if( m_cutoff != 0 ) CrCoordinateTransfer::s_fixedCutoff=0;

    return StatusCode::SUCCESS;
}





