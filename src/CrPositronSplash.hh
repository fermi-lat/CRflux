/**
 * CrPositronSplash:
 *   The splash cosmic-ray positron spectrum (and incident angle) source.
 */

//$Header$

#ifndef CrPositronSplash_H
#define CrPositronSplash_H

#include <utility>
#include <string>
#include "CrSpectrum.hh"

// Forward declaration:
class HepRandomEngine;
class CrPositronSplash_0003;
class CrPositronSplash_0306;
class CrPositronSplash_0608;
class CrPositronSplash_0809;
class CrPositronSplash_0910;
class CrPositronSplash_1011;

class CrPositronSplash : public CrSpectrum
{
public:
  CrPositronSplash();
  ~CrPositronSplash();

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
  CrPositronSplash_0003* crPositronSplash_0003;
  CrPositronSplash_0306* crPositronSplash_0306;
  CrPositronSplash_0608* crPositronSplash_0608;
  CrPositronSplash_0809* crPositronSplash_0809;
  CrPositronSplash_0910* crPositronSplash_0910;
  CrPositronSplash_1011* crPositronSplash_1011;

};
#endif // CrPositronSplash_H
