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
class CrElectronReentrant_0003;
class CrElectronReentrant_0306;
class CrElectronReentrant_0608;
class CrElectronReentrant_0809;
class CrElectronReentrant_0910;
class CrElectronReentrant_1011;

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
  CrElectronReentrant_0003* crElectronReentrant_0003;
  CrElectronReentrant_0306* crElectronReentrant_0306;
  CrElectronReentrant_0608* crElectronReentrant_0608;
  CrElectronReentrant_0809* crElectronReentrant_0809;
  CrElectronReentrant_0910* crElectronReentrant_0910;
  CrElectronReentrant_1011* crElectronReentrant_1011;


};
#endif // CrElectronReentrant_H
