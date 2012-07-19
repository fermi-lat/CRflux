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
class CLHEP::HepRandomEngine;
class CrPositronReentrant_0001;
class CrPositronReentrant_0102;
class CrPositronReentrant_0203;
class CrPositronReentrant_0304;
class CrPositronReentrant_0405;
class CrPositronReentrant_0506;
class CrPositronReentrant_0611;

class CrPositronReentrant : public CrSpectrum
{
public:
  CrPositronReentrant();
  ~CrPositronReentrant();

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
  CrPositronReentrant_0001* crPositronReentrant_0001;
  CrPositronReentrant_0102* crPositronReentrant_0102;
  CrPositronReentrant_0203* crPositronReentrant_0203;
  CrPositronReentrant_0304* crPositronReentrant_0304;
  CrPositronReentrant_0405* crPositronReentrant_0405;
  CrPositronReentrant_0506* crPositronReentrant_0506;
  CrPositronReentrant_0611* crPositronReentrant_0611;


};
#endif // CrPositronReentrant_H
