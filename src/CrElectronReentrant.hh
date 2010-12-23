/**
 * CrElectronReentrant:
 *   The reentrant cosmic-ray electron spectrum (and incident angle) source.
 */

//$Header$

#ifndef CrElectronReentrant_H
#define CrElectronReentrant_H

#include <utility>
#include <string>

#include "CrSpectrum.hh"

// Forward declaration:
class CLHEP::HepRandomEngine;
class CrElectronReentrant_0001;
class CrElectronReentrant_0102;
class CrElectronReentrant_0203;
class CrElectronReentrant_0304;
class CrElectronReentrant_0405;
class CrElectronReentrant_0506;
class CrElectronReentrant_0611;

class CrElectronReentrant : public CrSpectrum
{
public:
  CrElectronReentrant();
  ~CrElectronReentrant();

  // Gives back particle direction in (cos(theta), phi)
  std::pair<double,double> dir(double energy, CLHEP::HepRandomEngine* engine) const;

  // Gives back particle energy
  double energySrc(CLHEP::HepRandomEngine* engine) const;

  // flux() returns the value averaged over the region from which
  // the particle is coming from and the unit is [c/s/m^2/sr]
  double flux() const;

  // Gives back solid angle from which particle comes
  double solidAngle() const;

  // Gives back particle name
  const char* particleName() const;

  // Gives back the name of the component
  std::string title() const;

  // cr particle generator sorted in theta_M
private:
  CrElectronReentrant_0001* crElectronReentrant_0001;
  CrElectronReentrant_0102* crElectronReentrant_0102;
  CrElectronReentrant_0203* crElectronReentrant_0203;
  CrElectronReentrant_0304* crElectronReentrant_0304;
  CrElectronReentrant_0405* crElectronReentrant_0405;
  CrElectronReentrant_0506* crElectronReentrant_0506;
  CrElectronReentrant_0611* crElectronReentrant_0611;


};
#endif // CrElectronReentrant_H
