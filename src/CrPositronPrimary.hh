/**
 * CrPositronPrimary:
 *   The primary cosmic-ray positron spectrum (and incident angle) source.
 */

//$Header$

#ifndef CrPositronPrimary_H
#define CrPositronPrimary_H

#include <utility>
#include <string>
#include "CrSpectrum.hh"

// Forward declaration:
class HepRandomEngine;

class CrPositronPrimary : public CrSpectrum
{
public:
  CrPositronPrimary();
  ~CrPositronPrimary();
  
  // Set satellite position, altitude and observation time and
  // calculate energies related to COR.
  // These energies will be used to generate particles.
  void setPosition(double latitude, double longitude);
  void setPosition(double latitude, double longitude, double time);
  void setPosition(double latitude, double longitude, double time, double altitude);

  // Set geomagnetic cutoff rigidity and calculate the energies related.
  // These energies are used to generate the particle. 
  void setCutOffRigidity(double cor);

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
};
#endif // CrPositronPrimary_H

