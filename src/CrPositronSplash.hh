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
class CLHEP::HepRandomEngine;
class CrPositronSplash_0001;
class CrPositronSplash_0102;
class CrPositronSplash_0203;
class CrPositronSplash_0304;
class CrPositronSplash_0405;
class CrPositronSplash_0506;
class CrPositronSplash_0611;

class CrPositronSplash : public CrSpectrum
{
public:
  CrPositronSplash();
  ~CrPositronSplash();

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
  CrPositronSplash_0001* crPositronSplash_0001;
  CrPositronSplash_0102* crPositronSplash_0102;
  CrPositronSplash_0203* crPositronSplash_0203;
  CrPositronSplash_0304* crPositronSplash_0304;
  CrPositronSplash_0405* crPositronSplash_0405;
  CrPositronSplash_0506* crPositronSplash_0506;
  CrPositronSplash_0611* crPositronSplash_0611;


};
#endif // CrPositronSplash_H
