/**
 * CrNeutron:
 *  The class that calls each cosmic-ray alpha components based on
 *  their flux.
 */

//$Header$

#ifndef CrNeutron_H
#define CrNeutron_H

#include <vector>
#include <utility>
#include <string>

// Activate the next line when in CRflux package
#include "flux/ISpectrum.h"
// Activate the next line when in end-to-end simulation framework
// of Tune's group
//#include "ISpectrum.h"

class CrSpectrum;
class CLHEP::HepRandomEngine;

class CrNeutron : public ISpectrum
{
public:
  // params[0] is bit flag for determining which to include (default 4)
  // 1: primary
  // 2: reentrant
  // 4: splash
  CrNeutron(const std::string& params);
  virtual ~CrNeutron();

  // Gives back component in the ratio of the flux
  CrSpectrum* selectComponent();

  // Gives back energy
  virtual double energy(double time);

  // Gives back paticle direction in cos(theta) and phi[rad]
  virtual std::pair<double,double> dir(double energy);

  // Gives back the total flux (summation of each component's flux)
  virtual double flux (double time) const;  // calculate the flux [c/s/m^2/sr]

  // Gives back the particle kind
  virtual const char* particleName() const{ return "neutron"; }

  // Gives back the component name
  virtual std::string title() const{ return "CrNeutron"; }

  // Gives back solid angle from which particles come
  virtual double solidAngle() const;

  // print out the information of each component
  void dump();

  // Gives back the interval to the next event
  virtual double interval(double time);

private:
  std::vector<CrSpectrum*>  m_subComponents;
  CrSpectrum*               m_component;
  CLHEP::HepRandomEngine* m_engine;
};

#endif // CrNeutron_H

