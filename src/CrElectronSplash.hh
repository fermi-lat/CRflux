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
class HepRandomEngine;
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
  CrElectronSplash_0003* crElectronSplash_0003;
  CrElectronSplash_0306* crElectronSplash_0306;
  CrElectronSplash_0608* crElectronSplash_0608;
  CrElectronSplash_0809* crElectronSplash_0809;
  CrElectronSplash_0910* crElectronSplash_0910;
  CrElectronSplash_1011* crElectronSplash_1011;

};
#endif // CrElectronSplash_H
