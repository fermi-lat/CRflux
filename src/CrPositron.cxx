/**************************************************************************
 * CrPositron.cc
 **************************************************************************
 * Read comments at the head of in CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program interfaces 3 spectra generators, CrPositronPrimary.cc
 * CrPositronReentrant.cc, and CrPositronSplash.cc to 
 * CosmicRayGeneratorAction.cc.
 **************************************************************************
 * This program defines the entry-point class for the cosmic ray positron 
 * generation and interfaces to the primary, reentrant, and splash
 * cosmic-ray positron generators.
 **************************************************************************
 * 2001-07 Written by Y. Fukazawa (Hiroshima Univ.)
 * 2001-10 Modified by T. Mizuno and Y. Fukazawa
 * 2001-12 Modified by T. Mizuno to construct a `stand-alone' module
 ****************************************************************************
 */

//$Header$

#include <math.h>
#include <map>
#include <vector>

// CLHEP
#include <CLHEP/Random/JamesRandom.h>

#include "CrPositron.hh"
#include "CrPositronPrimary.hh"
#include "CrPositronSplash.hh"
#include "CrPositronReentrant.hh"

#include "CrSpectrum.hh"

// define a factory for anonomous instantiation
#include "FluxSvc/ISpectrumFactory.h"

typedef  double G4double;
namespace{
  const G4double pi    = 3.14159265358979323846264339;
}

// Constructor. Includes each component
CrPositron::CrPositron(const std::string& paramstring)
: m_component(0)
{
  std::vector<float> params;
  // including each component (primary/re-entrant/splash positrons)...
  m_subComponents.push_back(new CrPositronPrimary);
  m_subComponents.push_back(new CrPositronReentrant);
  m_subComponents.push_back(new CrPositronSplash);

  m_engine = new HepJamesRandom;
}


// Destructor. Delete each component
CrPositron::~CrPositron()
{
  std::vector<CrSpectrum*>::iterator  i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    delete *i;
  }
}


// Gives back component in the ratio of the flux
CrSpectrum* CrPositron::selectComponent()
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


// Gives back energy
double CrPositron::energy(double time)
{
  selectComponent();
  return m_component->energySrc(m_engine);
}


// Gives back paticle direction in cos(theta) and phi[rad]
std::pair<double,double> CrPositron::dir(double energy)
{
  if (!m_component){ selectComponent(); }

  return m_component->dir(energy, m_engine);
}


// Gives back the total flux (summation of each component's flux)
G4double CrPositron::flux(double time) const
{
  G4double          total_flux = 0;
  std::vector<CrSpectrum*>::const_iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
  }
  return total_flux;
}

// Gives back solid angle from whick particles come
G4double CrPositron::solidAngle() const
{
  return 4 * pi;
}

// Gives back the interval to the next event
G4double CrPositron::interval(double time){
  return -1.0;
}

