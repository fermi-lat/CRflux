/**
 * CrNeutronSplash:
 *   Atmospheric neutron spectrum source.
 */

//$Header$

#ifndef CrNeutronSplash_H
#define CrNeutronSplash_H

#include <utility>
#include <string>
#include "CrSpectrum.hh"

// Forward declaration:
class CLHEP::HepRandomEngine;

class CrNeutronSplash : public CrSpectrum
{
public:
  CrNeutronSplash();
  ~CrNeutronSplash();
  
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
};
#endif // CrNeutronSplash_H

