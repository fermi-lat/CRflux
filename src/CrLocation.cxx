#include "CrLocation.h"

// Static pointer to CrLocation
CrLocation* CrLocation::_instance = 0;

// Singleton pattern
CrLocation* CrLocation::instance() {
   if( _instance == 0) {
      _instance = new CrLocation();
   }
   return _instance;
}

CrLocation::CrLocation() {}

CrLocation::~CrLocation() {}


void CrLocation::setFluxSvc(IFluxSvc* f)
{
   fsvc = f;
}

IFluxSvc* CrLocation::getFluxSvc(void)
{
   return fsvc;
}
