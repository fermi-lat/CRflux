/**
 * CrPositronReentrant:
 *   The reentrant cosmic-ray positron spectrum (and incident angle) source.
 */

//$Header$

#ifndef CrPositronReentrant_H
#define CrPositronReentrant_H

#include <utility>
#include <string>

#include "CrSpectrum.hh"

// Forward declaration:
class HepRandomEngine;
class CrPositronReentrant_0003;
class CrPositronReentrant_0306;
class CrPositronReentrant_0608;
class CrPositronReentrant_0809;
class CrPositronReentrant_0910;
class CrPositronReentrant_1011;

class CrPositronReentrant : public CrSpectrum
{
public:
  CrPositronReentrant();
  ~CrPositronReentrant();

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
  CrPositronReentrant_0003* crPositronReentrant_0003;
  CrPositronReentrant_0306* crPositronReentrant_0306;
  CrPositronReentrant_0608* crPositronReentrant_0608;
  CrPositronReentrant_0809* crPositronReentrant_0809;
  CrPositronReentrant_0910* crPositronReentrant_0910;
  CrPositronReentrant_1011* crPositronReentrant_1011;

};
#endif // CrPositronReentrant_H
