//$Header$

#ifndef CrPositron_H
#define CrPositron_H

/**
 * CrPositron:
 *  The class that calls each cosmic-ray positron components based on
 *  their flux.
 */

#include <vector>
#include <utility>
#include <string>

#include "FluxSvc/ISpectrum.h"
//#include "ISpectrum.h"

class CrSpectrum;
class HepRandomEngine;

class CrPositron : public ISpectrum
{
public:
  // params[0] is bit flag for determining which to include (default 7)
  // 1: primary
  // 2: reentrant
  // 4: splash
  CrPositron(const std::string& params);
  virtual ~CrPositron();

  // Gives back component in the ratio of the flux
  CrSpectrum* selectComponent();

  // Gives back energy
  virtual double energy(double time);

  // Gives back paticle direction in cos(theta) and phi[rad]
  virtual std::pair<double,double> dir(double energy);

  // Gives back the total flux (summation of each component's flux)
  virtual double flux(double time) const;  // calculate the flux [c/s/m^2/sr]

  // Gives back the particle kind
  virtual const char* particleName() const{ return "e+"; }

  // Gives back the component name
  virtual std::string title() const { return "CrPositron"; }

  // Gives back solid angle from which particles come
  virtual double solidAngle() const;

  // Gives back the interval to the next event
  virtual double interval(double time);

  //parses the string of input sent to the constructor
  void parseParamList(std::string input, std::vector<float>& output);

private:
  std::vector<CrSpectrum*>  m_subComponents;
  CrSpectrum*               m_component;
  HepRandomEngine* m_engine;
};
#endif // CrPositron_H

