/**
 * CrCoordinateTransfer:
 *  This class calculate the geomagnetic latitude and longitude
 *  from geographic ones.
 */

//$Header$

#ifndef CrCoordinateTransfer_H
#define CrCoordinateTransfer_H

#include "math.h"

class CrCoordinateTransfer
{
public:
  CrCoordinateTransfer();
  ~CrCoordinateTransfer();

  // Calculate geomagnetic latitude and longitude from 
  // geographic ones
  // The input parameters (lat_deg and lon_deg) and return values are
  // given in degree, for consistency with a FluxSvc package.
  double geomagneticLatitude(double lat, double lon) const; // [degree]
  double geomagneticLongitude(double lat, double lon) const; // [degree]

private:
  // latitude and longitude of the geomagnetic north pole, in unit of degree
  double latitude_pole; // [deg]
  double longitude_pole; // [deg]

};

#endif CrCoordinateTransfer_H
