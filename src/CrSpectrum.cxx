/****************************************************************************
 * CrSpectrum.cc:
 ****************************************************************************
 * Read comments at the head of in CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note No. (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 ****************************************************************************
 * This program gives a base class for each cosmic-ray component
 * (CrProtonPrimary.cc, CrProtonReentrant.cc, and so on)
 * and calculate geomegnetic cutoff rigidity and 
 * force-field approximation potential from epoch and position 
 * (geographical lattitude and longitude).
 * At the moment, this code gives back fixed values of
 * geomagnetic cutoff rigidity and solar potential for
 * GLAST Balloon Experiment (in July 2001 at Palestine).
 ****************************************************************************
 * 2001-04 Written by M. Ozaki (ISAS)
 * 2001-05 Modified by T. Mizuno (Hiroshima Univ.) 
 * 2001-05 Program checked by T. Kamae and H. Mizushima (Hiroshima)
 * 2001-12 Modified by T. Mizuno to construct a `stand-alone' module
 * 2002-05 Modified by T. Mizuno to calculate geomagnetic cutoff rigidity
 *           and solar potential
 * 2003-12 Modified by T. Mizuno
 *           user can set lower and upper energy to generate gammas.
 ****************************************************************************
 */

//$Header$

#include "CrSpectrum.hh"
#include <iostream>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include "CrLocation.h"
#include "CrCoordinateTransfer.hh"

typedef double G4double;

CrSpectrum::CrSpectrum()
{
  // earth radius in km
  m_earthRadius = 6380;

  //
  // initialize satellite position, time and altitude
  // and calculate the cut off rigidity and solar modulation potential
  //

  // elapsed seconds of 2001-11-01 (solar maximum) from 2001-01-01
  //  m_time = 304*86400; 
  // elapsed seconds of 2001-11-01 from 2001-01-01 + 5.5year 
  m_time = 304*86400+5.5*86400*365; // solar minimum

  // altitude of the satellite is assumed to be 500km,
  // typical value for low earth orbit
  m_altitude = 500.0; 

  // (geographic latitude and longitude) = (31.78deg, -95.73deg)
  // (GLAST Balloon Experiment)
  // m_latitude = 31.78;
  // m_longitude = -95.73;

  // set callback to be notified when the position changes
  m_observer.setAdapter( new ActionAdapter<CrSpectrum>(this,&CrSpectrum::askGPS) );
  CrLocation::instance()->getFluxSvc()->GPSinstance()->notification().attach( &m_observer);
  askGPS(); //initial setup

  // set lower and upper energy to generate gammas
  m_gammaLowEnergy = 1.0e-3; // 1 MeV
  m_gammaHighEnergy = 100.0; // 100 GeV
}

CrSpectrum::~CrSpectrum()
{
  ;
}

void CrSpectrum::setGammaLowEnergy(double ene){ 
  using std::cout;
  using std::endl;

  if (ene>m_gammaHighEnergy){
    cout << "low energy limit to generate gamma should be less than high energy limit of " 
         << m_gammaHighEnergy << " GeV" << endl;
    ene = m_gammaHighEnergy;
  }
  m_gammaLowEnergy = ene;
}

void CrSpectrum::setGammaHighEnergy(double ene){ 
  using std::cout;
  using std::endl;

  if (ene<m_gammaLowEnergy){
    cout << "high energy limit to generate gamma should be more than low energy limit of " 
         << m_gammaLowEnergy << " GeV" << endl;
    ene = m_gammaLowEnergy;
  }
  m_gammaHighEnergy = ene;
}

// set observation time, which is the elapsed seconds from 
// 2000-01-01 00:00:00 UT and can be negative value 
void CrSpectrum::setTime(double time)
{
  setPosition(m_latitude, m_longitude, time);
}

// set satellite position in geographic coordinate
void CrSpectrum::setPosition(double latitude, double longitude)
{
  setPosition(latitude, longitude, m_time);
}


// set satellite position in geographic coordinate
void CrSpectrum::setPosition(double latitude, double longitude, double time)
{
  setPosition(latitude, longitude, time, m_altitude);
}


// set satellite position in geographic coordinate
void CrSpectrum::setPosition
(double latitude, double longitude, double time, double altitude)
{
  using std::cout;
  using std::endl;

  m_latitude  = latitude;
  m_longitude  = longitude;
  m_time = time;
  m_altitude = altitude;
  // compute the geomagnetic coordinates

  CrCoordinateTransfer transfer;
  m_geomagneticLatitude = 
    transfer.geomagneticLatitude(m_latitude, m_longitude);
  m_geomagneticLongitude = 
    transfer.geomagneticLongitude(m_latitude, m_longitude);

  // compute the Cut-Off-Rigidity (GV)
  // ------------------------------
  // For particles incoming vertically to the magnetic field,
  // cut off rigidity (Rc) is calculated as
  // (Rc/GV) = 14.9 * (1+h/R)^-2 * (cos(theta_M))^4,
  // where h gives altitude from earth surface, R means an radius of earth,
  // and theta_M is geomagnetic lattitude.
  // Reference: 
  // "Handbook of space astronomy and astrophysics" 2nd edition, p225
  //   Zombeck, 1990, Cambridge University Press 
  // "High Energy Astrophysics" 2nd edition, p325-330
  //   M. S. Longair, 1992, Cambridge University Press
  // ------------------------------

  m_cutOffRigidity = 
#if 1
    14.9 * pow(1+m_altitude/m_earthRadius,-2)* pow(cos(m_geomagneticLatitude*M_PI/180.0),4);
#else
      8.8; // fixed for testing
#endif
  //  m_cutOffRigidity = 4.46;  // temporarily fixed to Palestine value.

  // magnetic cutoff rigidity is restricted in 0.5 < cor < 14.9[GV]
  if (m_cutOffRigidity < 0.5){
    cout << "In this program, geomagnetic cutoff rigidity is" << endl;
    cout << " restricted in 0.5 < cor < 14.9[GV]" << endl;
    m_cutOffRigidity = 0.5; 
    cout << "COR is set at " << m_cutOffRigidity << " [GV]" << endl;
  }
  if (m_cutOffRigidity > 14.9){ 
    cout << "In this program, geomagnetic cutoff rigidity is" << endl;
    cout << " restricted in 0.5 < cor < 14.9[GV]" << endl;
    m_cutOffRigidity = 14.9; 
    cout << "COR is set at " << m_cutOffRigidity << " [GV]"<< endl;
  }

  // compute the force-field approximation potential (MV)
  // solar activity is assumed to be minimum in 1996-05-01,
  //                   assumed to be maximum in 2001-11-01,
  //               and assumed to be minumim in 2007-05-01.
  // Typical value of solar potential is 540 MV for minimum solar activity 
  // and 1100 MV for maximum solar activity.
  // Cosine curve is used to interpolate.
#if 0 // Original version for Palestine
  // elapsed seconds of 2001-11-01 from 2000-01-01
  double time_0 = (304+365)*86400; 
#else // modified for GLAST mission start at 2001-01-01: this is 2001-11-01 - 2001-01-01
    double time_0 = 304*86400;
#endif
  

  m_solarWindPotential = 820+280*cos(2*M_PI*(m_time-time_0)/(11*365*86400));

  // Solar potential is restricted in 500 < phi < 1100[MV]
  if (m_solarWindPotential < 500.0){ 
    cout << "In this program, solar potential is" << endl;
    cout << " restricted in 500 < phi < 1100[MV]" << endl;
    m_solarWindPotential = 500.0; 
    cout << "Phi is set at " << m_solarWindPotential << " [MV]" << endl;
  }
  if (m_solarWindPotential > 1100.0){ 
    cout << "In this program, solar potential is" << endl;
    cout << " restricted in 500 < phi < 1100[MV]" << endl;
    m_solarWindPotential = 1100.0; 
    cout << "Phi is set at " << m_solarWindPotential << " [MV]" << endl;
  }

}


// set solar modulation potential
void CrSpectrum::setSolarWindPotential(double phi)
{
  using std::cout;
  using std::endl;

  m_solarWindPotential = phi;

  // Solar potential is restricted in 500 < phi < 1100[MV]
  if (m_solarWindPotential < 500.0){ 
    cout << "In this program, solar potential is" << endl;
    cout << " restricted in 500 < phi < 1100[MV]" << endl;
    m_solarWindPotential = 500.0; 
    cout << "Phi is set at " << m_solarWindPotential << " [MV]" << endl;
  }
  if (m_solarWindPotential > 1100.0){ 
    cout << "In this program, solar potential is" << endl;
    cout << " restricted in 500 < phi < 1100[MV]" << endl;
    m_solarWindPotential = 1100.0; 
    cout << "Phi is set at " << m_solarWindPotential << " [MV]" << endl;
  }
}

// set cutoff rigidity
void CrSpectrum::setCutOffRigidity(double cor)
{
  using std::cout;
  using std::endl;

  m_cutOffRigidity = cor;
  // magnetic cutoff rigidity is restricted in 0.5 < cor < 14.9[GV]
  if (m_cutOffRigidity < 0.5){
    cout << "In this program, geomagnetic cutoff rigidity is" << endl;
    cout << " restricted in 0.5 < cor < 14.9[GV]" << endl;
    m_cutOffRigidity = 0.5; 
    cout << "COR is set at " << m_cutOffRigidity << " [GV]" << endl;
  }
  if (m_cutOffRigidity > 14.9){ 
    cout << "In this program, geomagnetic cutoff rigidity is" << endl;
    cout << " restricted in 0.5 < cor < 14.9[GV]" << endl;
    m_cutOffRigidity = 14.9; 
    cout << "COR is set at " << m_cutOffRigidity << " [GV]"<< endl;
  }

  // calculate geomagnetic latitude
  double tmp;
  tmp = pow(m_cutOffRigidity/14.9*pow(1+m_altitude/m_earthRadius,2), 0.25);
  if (tmp>1.0){tmp=1.0;} // when m_altitude>0 and m_cutOffRigidity~14.9
  m_geomagneticLatitude = acos(tmp);
  // convert from radian to degree
  m_geomagneticLatitude = m_geomagneticLatitude*180.0/M_PI;
}


// Gives back observation time [s] from 2000-01-01 00:00:00 UT
double CrSpectrum::time() const
{
  return m_time;
}


// Gives back satellite altitude in km
double CrSpectrum::altitude() const
{
  return m_altitude;
}


// Gives back satellite latitude in [rad]
double CrSpectrum::latitude() const
{
  return m_latitude;
}


// Gives back satellite longitude in [rad]
double CrSpectrum::longitude() const
{
  return m_longitude;
}


// Gives back geomagnetic latitude in [rad]
double CrSpectrum::geomagneticLatitude() const
{
  return m_geomagneticLatitude;
}


// Gives back geomagnetic longitude in [rad]
double CrSpectrum::geomagneticLongitude() const
{
  return m_geomagneticLongitude;
}


// Gives back cutoff rigidity in [GV]
double CrSpectrum::cutOffRigidity() const
{
  return m_cutOffRigidity;
}

// Gives back cutoff rigidity of particular direction in[GV]
double CrSpectrum::cutOffRigidityThisDirection(double theta, double phi) const
{
// The complete formula for the Stoermer cutoff rigidity is given as
//   Rc = (M/2r^2) *
//     ( cos(theta_M)^4/(1+(1-cos(theta_M)^3*sin(theta)*sin(phi))^0.5)^2 )
// where theta and phi are zenith angle and azimuth angle, respectively.
// (theta=0 for vertical downward, phi=0 for particles from north and
//  phi=90deg for particles from east)
// (e.g., p. 100 of "High Energy Cosmic Rays" by Todor Stanev)
//
// Thus, the ratio of cor of particular direction to the vertical cor
// is 4/(1+(1-cos(theta_M)^3*sin(theta)*sin(phi))^0.5)^2

    // convert from degree to radian 
    double magLat = m_geomagneticLatitude/180.0*M_PI;
    // calculate geomagnetic cutoff of particle from (theta, phi)
    double cor = m_cutOffRigidity;
    double tmp = 1+pow(1-pow(cos(magLat),3)*sin(theta)*sin(phi),0.5);
    cor = cor*4.0/pow(tmp,2);
/***
    std::cout << "cor= " << cor << " tmp= " << tmp 
              << " theta= " << theta << " phi= " << phi << std::endl;
***/
    return cor;
}


// Gives back solar modulation potential in [MV]
double CrSpectrum::solarWindPotential() const
{
  return m_solarWindPotential;
}

// call back from GPS when position changes
int CrSpectrum::askGPS()
{
    astro::GPS* gps=CrLocation::instance()->getFluxSvc()->GPSinstance();

    astro::EarthCoordinate pos = gps->earthpos();
    setPosition(pos.latitude(), pos.longitude(), gps->time(), pos.altitude());
    
    return 0; // can't be void in observer pattern
}
