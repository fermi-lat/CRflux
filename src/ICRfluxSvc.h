// $Header$
// 
//!  \author: T. Burnett
//

#ifndef _H_IBFEMfluxSvc
#define _H_IBFEMfluxSvc

// includes
#include "GaudiKernel/IInterface.h"
#include <string>
#include <list>
#include <vector>

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IBFEMfluxSvc(910, 1 , 0); 

// forward declaration
class IParticlePropertySvc;

//! Abstract interface for the flux service, FluxSvc.
class  IBFEMfluxSvc : virtual public IInterface {
public:
 
};

#endif  // _H_IFluxSvc
