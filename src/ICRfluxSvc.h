/** @file ICRfluxSvc.h
*
* $Header$
*/


#ifndef _H_ICRfluxSvc_
#define _H_ICRfluxSvc_

#include "GaudiKernel/IInterface.h"
static const InterfaceID IID_ICRfluxSvc(901, 1 , 1); 


class   ICRfluxSvc : virtual public IInterface {
public:
  

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_ICRfluxSvc; }

};

#endif  // _H_ICRfluxSvc_

