/****************************************************************************
 * CrCoorinateTransfer.cc:
 ****************************************************************************
 * This class calculates the geomagnetic latitude and longitude from
 * the geographical ones.
 ****************************************************************************
 * 2002-05 Written by T. Mizuno (Hiroshima University)
 ****************************************************************************
 */

//$Header$

#include "CrCoordinateTransfer.hh"

CrCoordinateTransfer::CrCoordinateTransfer()
{
  // latitude and longitude of the geomagnetic north pole in 2000
  // is given below. The values are taken from
  // http://swdcdb.kugi.kyoto-u.ac.jp/trans/index.html
  latitude_pole = 1.388; // 79.55 degree
  longitude_pole = -1.249; // -71.57 degree
}


CrCoordinateTransfer::~CrCoordinateTransfer()
{
  ;
}

// calculate geomagnetic latitude from geographic coordinates
double CrCoordinateTransfer::geomagneticLatitude
(double latitude, double longitude) const
{
  double m_geomagneticLatitude;

  // calculate geomagnetic latitude
  m_geomagneticLatitude = 
    asin( sin(latitude)*sin(latitude_pole) + 
	  cos(latitude)*cos(latitude_pole)*cos(longitude - longitude_pole) );

  // return geomagnetic latitude
  return m_geomagneticLatitude;
}

// calculate geomagnetic longitude from geographic coordinates
double CrCoordinateTransfer::geomagneticLongitude
(double latitude, double longitude) const
{
  double m_geomagneticLatitude, m_geomagneticLongitude;

  // Since geomagnetic latitude and longitude are correlated with each other,
  // we need to calculate both values.
  m_geomagneticLatitude = 
    asin( sin(latitude)*sin(latitude_pole) + 
	  cos(latitude)*cos(latitude_pole)*cos(longitude - longitude_pole) );
  m_geomagneticLongitude = 
    acos ( (cos(latitude)*cos(longitude - longitude_pole) - 
            cos(latitude_pole)*sin(m_geomagneticLatitude))
           /(sin(latitude_pole)*cos(m_geomagneticLatitude)) );
  if ( cos(latitude)*sin(longitude - longitude_pole)/cos(m_geomagneticLatitude) < 0.0 ){
    m_geomagneticLongitude = m_geomagneticLongitude + 3.1415926;
  }

  // return geomagnetic longitude
  return m_geomagneticLongitude;
}
