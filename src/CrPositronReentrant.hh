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
  CrPositronReentrant_0003* crPositronReentrant_0003;
  CrPositronReentrant_0306* crPositronReentrant_0306;
  CrPositronReentrant_0608* crPositronReentrant_0608;
  CrPositronReentrant_0809* crPositronReentrant_0809;
  CrPositronReentrant_0910* crPositronReentrant_0910;
  CrPositronReentrant_1011* crPositronReentrant_1011;

};
#endif // CrPositronReentrant_H
