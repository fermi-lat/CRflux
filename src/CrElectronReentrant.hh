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
class HepRandomEngine;
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

  //
  // "flux" package stuff
  //
  /// r in [0,1)
  float operator()(float r);
  double calculate_rate(double old_rate);
  float flux(float latitude, float longitude) const;
  float flux(std::pair<double,double> coords) const;
  std::string title() const;
  /// fraction function doesn't work in this class... :-(
  float fraction(float energy);
  std::pair<float,float> dir(float energy) const;

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
