// $Header$



// Include files
#include "FluxSvc/IFluxSvc.h"
#include "FluxSvc/IFlux.h"


// GlastEvent for creating the McEvent stuff
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"

// Gaudi system includes
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include <list>
#include <string>
#include <vector>
#include "GaudiKernel/ParticleProperty.h"


// CR includes
#include "../CrExample.h"
#include "../CrProton.hh"
#include "../CrAlpha.hh"
#include "../CrElectron.hh"
#include "../CrPositron.hh"
#include "../CrGamma.hh"

//#include "FluxAlg.h"
/*! \class CRTestAlg
\brief 

*/

class CRTestAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    CRTestAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    

private:
    IFlux* m_flux;
    IFluxSvc* m_fsvc; /// pointer to the flux Service 
    std::string m_source_name;
    IParticlePropertySvc * m_partSvc;
};


static const AlgFactory<CRTestAlg>  Factory;
const IAlgFactory& CRTestAlgFactory = Factory;

//------------------------------------------------------------------------------
//
CRTestAlg::CRTestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){
    
    declareProperty("source_name", m_source_name="default");
}


//------------------------------------------------------------------------------
/*! */
StatusCode CRTestAlg::initialize() {
    
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initializing..." << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // get the service
    StatusCode sc = service("FluxSvc", m_fsvc);

    return sc;
}


//------------------------------------------------------------------------------
StatusCode CRTestAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );    

  std::vector<const char*> arguments;
  arguments.push_back("CrExample");
  //  arguments.push_back("CrProton");
  //  arguments.push_back("CrAlpha");
  //  arguments.push_back("CrElectron");
  //  arguments.push_back("CrPositron");
  arguments.push_back("CrGamma");
  m_fsvc->rootDisplay(arguments);

    return sc;
}


//------------------------------------------------------------------------------
StatusCode CRTestAlg::finalize() {
    
    return StatusCode::SUCCESS;
}






