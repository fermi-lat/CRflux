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

  // calculate geomagnetic latitude and longitude from 
  // geographic ones
  double geomagneticLatitude(double lat, double lon) const; // [radian]
  double geomagneticLongitude(double lat, double lon) const; // [radian]

private:
  // latitude and longitude of the geomagnetic north pole
  double latitude_pole;
  double longitude_pole;

};

#endif CrCoordinateTransfer_H
