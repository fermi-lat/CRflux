/**
 * CrElectronSplash:
 *   The splash cosmic-ray electron spectrum (and incident angle) source.
 */

//$Header$

#ifndef CrElectronSplash_H
#define CrElectronSplash_H

#include <utility>
#include <string>
#include "CrSpectrum.hh"

// Forward declaration:
class CLHEP::HepRandomEngine;
class CrElectronSplash_0003;
class CrElectronSplash_0306;
class CrElectronSplash_0608;
class CrElectronSplash_0809;
class CrElectronSplash_0910;
class CrElectronSplash_1011;

class CrElectronSplash : public CrSpectrum
{
public:
  CrElectronSplash();
  ~CrElectronSplash();

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
  CrElectronSplash_0003* crElectronSplash_0003;
  CrElectronSplash_0306* crElectronSplash_0306;
  CrElectronSplash_0608* crElectronSplash_0608;
  CrElectronSplash_0809* crElectronSplash_0809;
  CrElectronSplash_0910* crElectronSplash_0910;
  CrElectronSplash_1011* crElectronSplash_1011;

};
#endif // CrElectronSplash_H
