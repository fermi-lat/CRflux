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

// CLHEP
#include <CLHEP/config/CLHEP.h>

typedef double G4double;

CrCoordinateTransfer::CrCoordinateTransfer()
{
  // latitude and longitude of the geomagnetic north pole in 2000
  // is given below. The values are taken from
  // http://swdcdb.kugi.kyoto-u.ac.jp/trans/index.html
  latitude_pole = 79.55; // 1.388 [rad]
  longitude_pole = -71.57; //  -1.249 [rad]
}


CrCoordinateTransfer::~CrCoordinateTransfer()
{
  ;
}

// calculate geomagnetic latitude from geographic coordinates
double CrCoordinateTransfer::geomagneticLatitude
(double latitude, double longitude) const // parameters are given in degree
{
  double m_geomagneticLatitude;

  // calculate geomagnetic latitude in unit of radian
  m_geomagneticLatitude = 
    asin( sin(latitude*M_PI/180.0)*sin(latitude_pole*M_PI/180.0) + 
	  cos(latitude*M_PI/180.0)*cos(latitude_pole*M_PI/180.0)
	  *cos(longitude*M_PI/180.0 - longitude_pole*M_PI/180.0) );

  // in case that sin(m_geomagneticLatitude) is not from -1.0 to 1.0 due to rounding error
  if ( sin(latitude*M_PI/180.0)*sin(latitude_pole*M_PI/180.0) + 
	  cos(latitude*M_PI/180.0)*cos(latitude_pole*M_PI/180.0)
	  *cos(longitude*M_PI/180.0 - longitude_pole*M_PI/180.0)<-0.999 ){
    m_geomagneticLatitude = -M_PI/2;
  }
  if ( sin(latitude*M_PI/180.0)*sin(latitude_pole*M_PI/180.0) + 
	  cos(latitude*M_PI/180.0)*cos(latitude_pole*M_PI/180.0)
	  *cos(longitude*M_PI/180.0 - longitude_pole*M_PI/180.0)>0.999 ){
    m_geomagneticLatitude = M_PI/2;
  }

  // convert radian to degree
  m_geomagneticLatitude = m_geomagneticLatitude*180.0/M_PI;

  // return geomagnetic latitude
  return m_geomagneticLatitude;
}

// calculate geomagnetic longitude from geographic coordinates
double CrCoordinateTransfer::geomagneticLongitude
(double latitude, double longitude) const // parameters are given in degree
{
  double m_geomagneticLatitude, m_geomagneticLongitude;

  // Since geomagnetic latitude and longitude are correlated with each other,
  // we need to calculate both values.
  m_geomagneticLatitude = 
    asin( sin(latitude*M_PI/180.0)*sin(latitude_pole*M_PI/180.0) + 
	  cos(latitude*M_PI/180.0)*cos(latitude_pole*M_PI/180.0)
	  *cos(longitude*M_PI/180.0 - longitude_pole*M_PI/180.0) );
  // convert from radian to degree
  m_geomagneticLatitude =   m_geomagneticLatitude*180.0/M_PI;

  m_geomagneticLongitude = 
    acos ( (cos(latitude*M_PI/180.0)*cos(longitude*M_PI/180.0 - longitude_pole*M_PI/180.0) - 
            cos(latitude_pole*M_PI/180.0)*sin(m_geomagneticLatitude*M_PI/180.0))
           /(sin(latitude_pole*M_PI/180.0)*cos(m_geomagneticLatitude*M_PI/180.0)) );
  // in case that cos(m_geomagneticLongitude) is not from -1.0 to 1.0 due to rounding error
  if ( (cos(latitude*M_PI/180.0)*cos(longitude*M_PI/180.0 - longitude_pole*M_PI/180.0) - 
            cos(latitude_pole*M_PI/180.0)*sin(m_geomagneticLatitude*M_PI/180.0))
    /(sin(latitude_pole*M_PI/180.0)*cos(m_geomagneticLatitude*M_PI/180.0))<-0.999 ){
    m_geomagneticLongitude = M_PI;
  }
  if ( (cos(latitude*M_PI/180.0)*cos(longitude*M_PI/180.0 - longitude_pole*M_PI/180.0) - 
            cos(latitude_pole*M_PI/180.0)*sin(m_geomagneticLatitude*M_PI/180.0))
    /(sin(latitude_pole*M_PI/180.0)*cos(m_geomagneticLatitude*M_PI/180.0))>0.999 ){
    m_geomagneticLongitude = 0.0;
  }

  if ( cos(latitude*M_PI/180.0)*sin(longitude*M_PI/180.0 - longitude_pole*M_PI/180.0)/cos(m_geomagneticLatitude*M_PI/180.0) < 0.0 ){
    m_geomagneticLongitude = m_geomagneticLongitude + M_PI;
  }

  // convert from radian to degree
  m_geomagneticLongitude = m_geomagneticLongitude*180.0/M_PI;

  // return geomagnetic longitude
  return m_geomagneticLongitude;
}
