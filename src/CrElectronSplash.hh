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
class CrElectronSplash_0001;
class CrElectronSplash_0102;
class CrElectronSplash_0203;
class CrElectronSplash_0304;
class CrElectronSplash_0405;
class CrElectronSplash_0506;
class CrElectronSplash_0611;

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
  CrElectronSplash_0001* crElectronSplash_0001;
  CrElectronSplash_0102* crElectronSplash_0102;
  CrElectronSplash_0203* crElectronSplash_0203;
  CrElectronSplash_0304* crElectronSplash_0304;
  CrElectronSplash_0405* crElectronSplash_0405;
  CrElectronSplash_0506* crElectronSplash_0506;
  CrElectronSplash_0611* crElectronSplash_0611;


};
#endif // CrElectronSplash_H
