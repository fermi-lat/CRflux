/**
 * CrProtonReentrant:
 *   The reentrant cosmic-ray proton spectrum (and incident angle) source.
 */

//$Header$

#ifndef CrProtonReentrant_H
#define CrProtonReentrant_H

#include <utility>
#include <string>

#include "CrSpectrum.hh"

// Forward declaration:
class HepRandomEngine;
class CrProtonReentrant_0002;
class CrProtonReentrant_0203;
class CrProtonReentrant_0304;
class CrProtonReentrant_0405;
class CrProtonReentrant_0506;
class CrProtonReentrant_0607;
class CrProtonReentrant_0708;
class CrProtonReentrant_0809;
class CrProtonReentrant_0910;

class CrProtonReentrant : public CrSpectrum
{
public:
  CrProtonReentrant();
  ~CrProtonReentrant();

  // Gives back particle direction in (cos(theta), phi)
  std::pair<double,double> dir(double energy, HepRandomEngine* engine) const;

  // Gives back particle energy
  double energySrc(HepRandomEngine* engine) const;

  // flux() returns the value integrated over whole energy and direction
  // and devided by 4pi sr: then the unit is [c/s/m^2/sr]
  double flux() const;

  // Gives back solid angle from which particle comes
  double solidAngle() const;

  // Gives back particle name
  const char* particleName() const;

  // Gives back the name of the component
  std::string title() const;

  // cr particle generator sorted in theta_M
private:
  CrProtonReentrant_0002* crProtonReentrant_0002;
  CrProtonReentrant_0203* crProtonReentrant_0203;
  CrProtonReentrant_0304* crProtonReentrant_0304;
  CrProtonReentrant_0405* crProtonReentrant_0405;
  CrProtonReentrant_0506* crProtonReentrant_0506;
  CrProtonReentrant_0607* crProtonReentrant_0607;
  CrProtonReentrant_0708* crProtonReentrant_0708;
  CrProtonReentrant_0809* crProtonReentrant_0809;
  CrProtonReentrant_0910* crProtonReentrant_0910;

};
#endif // CrProtonReentrant_H
