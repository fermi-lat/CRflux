/**************************************************************************
 * CrAlpha.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program interfaces CrAlphaPrimary.cc to CosmicRayGeneratorAction.cc 
 **************************************************************************
 * This program defines the entry-point class for the cosmic ray alpha 
 * generation and interfaces to the primary cosmic-ray alpha generators.
 **************************************************************************
 * 2002-04 Written by Y. Fukazawa (Hiroshima Univ.)
 * 2003-02 Modified by T. Mizuno to construct a `stand-alone' module
 ****************************************************************************
 */

//$Header$

#include <math.h>
#include <map>
#include <vector>

// CLHEP
#include <CLHEP/Random/JamesRandom.h>

#include "CrAlpha.hh"
#include "CrAlphaPrimary.hh"

#include "CrSpectrum.hh"

// define a factory for anonomous instantiation
#include "FluxSvc/ISpectrumFactory.h"

typedef  double G4double;
namespace{
  const G4double pi    = 3.14159265358979323846264339;
}

// Constructor. Includes each component
CrAlpha::CrAlpha(const std::string& paramstring)
: m_component(0)
{
  std::vector<float> params;
  // including each component (primary alphas)...
  m_subComponents.push_back(new CrAlphaPrimary);

  m_engine = new HepJamesRandom;
}


// Destructor. Delete each component
CrAlpha::~CrAlpha()
{
  std::vector<CrSpectrum*>::iterator  i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    delete *i;
  }
}


// Gives back component in the ratio of the flux
CrSpectrum* CrAlpha::selectComponent()
{
  std::map<CrSpectrum*,double>       integ_flux;
  double                             total_flux = 0;
  std::vector<CrSpectrum*>::iterator i;

  // calculate the total flux 
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
    integ_flux[*i] = total_flux;
  }
  // select component based on the flux
  double  rnum = m_engine->flat() * total_flux;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    if (integ_flux[*i] >= rnum) { break; }
  }

  m_component = *i;

  return *i;
}

// Gives back kinetic energy 
double CrAlpha::energy(double time)
{
  selectComponent();
  return m_component->energySrc(m_engine);
}


// Gives back paticle direction in cos(theta) and phi[rad]
std::pair<double,double> CrAlpha::dir(double energy)
{
  if (!m_component){ selectComponent(); }

  return m_component->dir(energy, m_engine);
}

// Gives back the total flux (summation of each component's flux)
G4double CrAlpha::flux(double time) const
{
  G4double          total_flux = 0;
  std::vector<CrSpectrum*>::const_iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
    //    cout << "flux of this component = " << (*i)->flux() << endl;
  }
  return total_flux;
}

// Gives back solid angle from whick particles come
G4double CrAlpha::solidAngle() const
{
  return 4 * pi;
}

// Gives back the interval to the next event
G4double CrAlpha::interval(double time){
  return -1.0;
}


