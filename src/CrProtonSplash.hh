/**
 * CrProtonSplash:
 *   The splash cosmic-ray proton spectrum (and incident angle) source.
 */

//$Header$

#ifndef CrProtonSplash_H
#define CrProtonSplash_H

#include <utility>
#include <string>

#include "CrSpectrum.hh"

// Forward declaration:
class HepRandomEngine;
class CrProtonSplash_0002;
class CrProtonSplash_0203;
class CrProtonSplash_0304;
class CrProtonSplash_0405;
class CrProtonSplash_0506;
class CrProtonSplash_0607;
class CrProtonSplash_0708;
class CrProtonSplash_0809;
class CrProtonSplash_0910;

class CrProtonSplash : public CrSpectrum
{
public:
  CrProtonSplash();
  ~CrProtonSplash();

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
  CrProtonSplash_0002* crProtonSplash_0002;
  CrProtonSplash_0203* crProtonSplash_0203;
  CrProtonSplash_0304* crProtonSplash_0304;
  CrProtonSplash_0405* crProtonSplash_0405;
  CrProtonSplash_0506* crProtonSplash_0506;
  CrProtonSplash_0607* crProtonSplash_0607;
  CrProtonSplash_0708* crProtonSplash_0708;
  CrProtonSplash_0809* crProtonSplash_0809;
  CrProtonSplash_0910* crProtonSplash_0910;

};
#endif // CrProtonSplash_H
