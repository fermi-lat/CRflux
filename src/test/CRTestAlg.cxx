// $Header$



// Include files
#include "FluxSvc/IFluxSvc.h"
#include "FluxSvc/IFlux.h"

#include "astro/GPS.h"


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
    DoubleProperty m_latitude;
    DoubleProperty m_longitude;
    DoubleProperty m_time;
    StringArrayProperty m_rootplot;
};


//static const AlgFactory<CRTestAlg>  Factory;
//const IAlgFactory& CRTestAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(CRTestAlg);

//------------------------------------------------------------------------------
//
CRTestAlg::CRTestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){

    declareProperty("source_name", m_source_name="default");
    declareProperty("latitude", m_latitude=20); // not useable yet
    declareProperty("longitude", m_longitude=20);
    declareProperty("rootplot", m_rootplot);
    declareProperty("time", m_time=0);
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
    m_fsvc->GPSinstance()->time(m_time); //try somethin
    m_fsvc->GPSinstance()->notifyObservers();
#if 1
    // make the root plots file here
    std::vector<std::string> sargs;

    for( std::vector<std::string>::const_iterator it = m_rootplot.value().begin(); it!=m_rootplot.value().end(); ++it){
        std::string arg(*it);
        sargs.push_back(arg);
    }

    try {
        m_fsvc->rootDisplay(sargs);
    }catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;; 
    }

#endif

    return sc;
}


//------------------------------------------------------------------------------
StatusCode CRTestAlg::execute() {

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );    

#if 0
    std::vector<std::string> arguments;
    //  arguments.push_back("CrExample");

    // The mix sources handle solid angles better

//    arguments.push_back("CrProton");
    arguments.push_back("CrProtonMix"); // new alternative for CrProton
//    arguments.push_back("CrProtonPrimary");
//    arguments.push_back("CrProtonReentrant");
//    arguments.push_back("CrProtonSplash");

    //  arguments.push_back("CrAlpha");

    //  arguments.push_back("CrElectron");
      arguments.push_back("CrElectronMix"); // alternative
    //  arguments.push_back("CrElectronPrimary");
    //  arguments.push_back("CrElectronReentrant");
    //  arguments.push_back("CrElectronSplash");

    //  arguments.push_back("CrPositron");
      arguments.push_back("CrPositronMix"); // alternative
    //  arguments.push_back("CrPositronPrimary");
    //  arguments.push_back("CrPositronReentrant");
    //  arguments.push_back("CrPositronSplash");

    //  arguments.push_back("CrGamma");
    //  arguments.push_back("CrGammaMix"); // alternative
    //  arguments.push_back("CrGammaPrimary");
    //  arguments.push_back("CrGammaSecondaryDownward");
    //  arguments.push_back("CrGammaSecondaryUpward");
    //  arguments.push_back("-no_integrate");

      arguments.push_back("CrNeutron");
    
    m_fsvc->rootDisplay(arguments);
#endif
    return sc;
}


//------------------------------------------------------------------------------
StatusCode CRTestAlg::finalize() {

    return StatusCode::SUCCESS;
}






