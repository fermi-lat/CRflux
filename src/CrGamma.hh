/**
 * CrGamma:
 *  The class that calls each gamma-ray components based on
 *  their flux.
 */

//$Header$

#ifndef CrGamma_H
#define CrGamma_H

#include <vector>
#include <utility>
#include <string>

// Activate the next line when in CRflux package
#include "flux/ISpectrum.h"
// Activate the next line when in end-to-end simulation framework
// of Tune's group
//#include "ISpectrum.h"

class CrSpectrum;
class HepRandomEngine;

class CrGamma : public ISpectrum
{
public:
  // params[0] is bit flag for determining which to include (default 7)
  // 1: primary
  // 2: reentrant
  // 4: splash
  CrGamma(const std::string& params);
  virtual ~CrGamma();

  // Gives back component in the ratio of the flux
  CrSpectrum* selectComponent();

  // Gives back energy
  virtual double energy(double time);

  // Gives back paticle direction in cos(theta) and phi[rad]
  virtual std::pair<double,double> dir(double energy);

  // Gives back the total flux (summation of each component's flux)
  virtual double flux(double time) const;  // calculate the flux [c/s/m^2/sr]

  // Gives back the particle kind
  virtual const char* particleName() const{ return "gamma"; }

  // Gives back the component name
  virtual std::string title() const{ return "CrGamma"; }

  // Gives back solid angle from which particles come
  virtual double solidAngle() const;

  // print out the information of each component
  void dump();

  // Gives back the interval to the next event
  virtual double interval(double time);

  //parses the string of input sent to the constructor
  void parseParamList(std::string input, std::vector<float>& output);

private:
  std::vector<CrSpectrum*>  m_subComponents;
  CrSpectrum*               m_component;
  HepRandomEngine* m_engine;
};
#endif // CrGamma_H

