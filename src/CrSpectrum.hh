/**
 * CrSpectrum:
 *  The base class for cosmic-ray spectrum generator class.
 *  All the CR spectrum generator should be derived from this class.
 */

//$Header$

#ifndef CrSpectrum_H
#define CrSpectrum_H

#include <string>
#include <utility>
#include "math.h"

#include "CrCoordinateTransfer.hh"

class HepRandomEngine;

class CrSpectrum
{
public:
  CrSpectrum();
  virtual ~CrSpectrum();

  // time is the elapsed seconds from 2000-01-01 00:00:00 UT 
  // and can be negative value.
  virtual void setTime(double time);

  // set satellite position in geographic coordinate
  // The unit of coordinates is in [radian].
  virtual void setPosition(double latitude, double longitude, double time, double altitude);
  virtual void setPosition(double latitude, double longitude, double time);
  virtual void setPosition(double latitude, double longitude);

  // set geomagnetic cutoff and solar potential
  virtual void setSolarWindPotential(double phi);
  virtual void setCutOffRigidity(double cor);

  // Gives back observation time, satellite altitude,
  // position in geographic/geomagnetic coordinate,
  // cutoff rigidity and solar modulation potential
  virtual double time() const; // [s]
  virtual double altitude() const; // [km]
  virtual double latitude() const; // [deg]
  virtual double longitude() const; // [deg]
  virtual double geomagneticLatitude() const; // [deg]
  virtual double geomagneticLongitude() const; // [deg]
  virtual double cutOffRigidity() const; // [GV]
  virtual double solarWindPotential() const; // [MV]

  // Gives back the energy and direction of the particle
  virtual double energySrc(HepRandomEngine* engine) const=0;
  virtual std::pair<double,double> 
  dir(double energy, HepRandomEngine* engine) const=0;

  // Gives back the flux
  virtual double flux() const=0;
  // Gives back the solid angle from which particle comes
  virtual double solidAngle() const=0;

  // Gives back the particle type and component name
  virtual const char* particleName() const=0;
  virtual std::string title() const=0;

protected:
  // Following member variables defines satellite position 
  // in geographic/geomagnetic coordinate, altitude, time of observation, 
  // geomagnetic cutoff rigidity and solar modulation potential.
  double m_time; // [s]
  double m_altitude; // [km]
  double m_latitude; // [radian]
  double m_longitude; // [radian]
  double m_geomagneticLatitude; // [radian]
  double m_geomagneticLongitude; // [radian]
  double m_cutOffRigidity; // [GV]
  double m_solarWindPotential; // [MV]

  double m_earthRadius; // earth radius in km
};

#endif // CrSpectrum_H
