#include "CrLocation.h"

// Static pointer to CrLocation
CrLocation* CrLocation::_instance = 0;

// Singleton pattern
CrLocation* CrLocation::Instance() {
   if( _instance == 0) {
      _instance = new CrLocation();
   }
   return _instance;
}

CrLocation::CrLocation() {}

CrLocation::~CrLocation() {}

double CrLocation::getLat() const { return latitude; }

double CrLocation::getLon() const { return longitude; }

void CrLocation::setLat(double lat) { latitude = lat; }

void CrLocation::setLon(double lon) { longitude = lon;}

void CrLocation::setFluxSvc(IFluxSvc* f)
{
   fsvc = f;
}

IFluxSvc* CrLocation::getFluxSvc(void)
{
   return fsvc;
}
