#ifndef CRLOCATION
#define CRLOCATION

/**
 * \class CrLocation
 *
 * \brief This is a singleton class used to temporarily store the latitude
 * and longitude values from FluxSvc for use in CrSpectrum.  It also stores
 * a pointer to FluxSvc.
 * \author Theodore Hierath, University of Washington, August 2003.
*/

#include "FluxSvc/IFluxSvc.h"

class CrLocation {
public:
   static CrLocation* instance();
   void setFluxSvc(IFluxSvc* f);
   IFluxSvc* getFluxSvc(void);

protected:
   CrLocation();
   virtual ~CrLocation();
private:
   static CrLocation* _instance;
   IFluxSvc* fsvc;
};


#endif
