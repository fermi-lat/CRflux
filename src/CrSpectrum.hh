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
#include "facilities/Observer.h"
#include "astro/EarthCoordinate.h"

class HepRandomEngine;

/** @class CrSpectrum 
 *  @brief base class
 *
 * This program gives a base class for each cosmic-ray component
 * (CrProtonPrimary.cc, CrProtonReentrant.cc, and so on)
 * and calculate geomegnetic cutoff rigidity and 
 * force-field approximation potential from epoch and position 
 * (geographical lattitude and longitude).
 * At the moment, this code gives back fixed values of
 * geomagnetic cutoff rigidity and solar potential for
 * GLAST Balloon Experiment (in July 2001 at Palestine).
*/
class CrSpectrum
{
public:
  CrSpectrum();
  virtual ~CrSpectrum();

  /// time is the elapsed seconds from 2000-01-01 00:00:00 UT 
  /// and can be negative value. [THB: GLAST mission start is 2001-01-01]
  virtual void setTime(double time);

  /// set satellite position in geographic coordinate
  /// The unit of coordinates is in [degree].
  virtual void setPosition(double lat, double lon, double time, double altitude);
  virtual void setPosition(double lat, double lon, double time);
  virtual void setPosition(double lat, double lon);

  /// set geomagnetic cutoff and solar potential
  virtual void setSolarWindPotential(double phi);
  virtual void setCutOffRigidity(double cor);

  /// Gives back observation time, satellite altitude,
  /// position in geographic/geomagnetic coordinate,
  /// cutoff rigidity and solar modulation potential
  /// Please note that the unit of angle is [degree] for consistency with FluxSvc.
  virtual double time() const; // [s]
  virtual double altitude() const; // [km]
  virtual double latitude() const; // [deg]
  virtual double longitude() const; // [deg]
  virtual double geomagneticLatitude() const; // [deg]
  virtual double geomagneticLongitude() const; // [deg]
  virtual double cutOffRigidity() const; // [GV]
  virtual double solarWindPotential() const; // [MV]

  /// Gives back the energy and direction of the particle
  virtual double energySrc(HepRandomEngine* engine) const=0;
  virtual std::pair<double,double> 
  dir(double energy, HepRandomEngine* engine) const=0;

  /// Gives back the flux
  virtual double flux() const=0;
  /// Gives back the solid angle from which particle comes
  virtual double solidAngle() const=0;

  /// Gives back the particle type and component name
  virtual const char* particleName() const=0;
  virtual std::string title() const=0;

  /// set lower and upper energy to generate gammas
  virtual void setGammaLowEnergy(double ene);
  virtual void setGammaHighEnergy(double ene);
  /// Gives back lower and upper energy to generate gammas
  inline virtual double gammaLowEnergy() const { return m_gammaLowEnergy;}
  inline virtual double gammaHighEnergy() const { return m_gammaHighEnergy;}

  /// this one asks the GPS for position
  int askGPS();

protected:
  // Following member variables defines satellite position 
  // in geographic/geomagnetic coordinate, altitude, time of observation, 
  // geomagnetic cutoff rigidity and solar modulation potential.
  double m_time; ///< [s]
  double m_altitude; ///< [km]
  double m_latitude; ///< [deg]
  double m_longitude; ///< [deg]
  double m_geomagneticLatitude; ///< [deg]
  double m_geomagneticLongitude; ///< [deg]
  double m_cutOffRigidity; ///< [GV]
  double m_solarWindPotential; ///< [MV]

  // lower and upper energy to generate gammas
  double m_gammaLowEnergy; ///< [GeV]
  double m_gammaHighEnergy; ///< [GeV]

  double m_earthRadius; ///< earth radius in km
private:
   ObserverAdapter< CrSpectrum > m_observer; ///< obsever tag

   //! will be set by the call back from GPS.
   astro::EarthCoordinate m_pos;
};

#endif // CrSpectrum_H
