/**************************************************************************
 * CrGamma.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note No. (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program interfaces 3 spectra generators, CrGammaPrimary.cc,
 * CrGammaSecondaryDownward.cc, and CrGammaSecondaryUpward.cc, 
 * to CosmicRayGeneratorAction.cc 
 **************************************************************************
 * This program defines the entry-point class for the cosmic ray gamma
 * generation and interfaces to the primary and secondary
 * cosmic-ray gamma generators.
 **************************************************************************
 * 2001-07 Written by Y. Fukazawa (Hiroshima Univ)
 * 2001-12 Modified by T. Mizuno in order to construct a `stand-alone' module
 ****************************************************************************
 */

// $Header$

#include <math.h>
#include <map>
#include <vector>

// CLHEP
#include <CLHEP/Random/JamesRandom.h>

#include "CrGamma.hh"
#include "CrGammaPrimary.hh"
#include "CrGammaSecondaryDownward.hh"
#include "CrGammaSecondaryUpward.hh"

#include "CrSpectrum.hh"

// define a factory for anonomous instantiation
#include "FluxSvc/ISpectrumFactory.h"

typedef  double G4double;
namespace{
  const G4double pi    = 3.14159265358979323846264339;
}


// Constructor. Includes each component
CrGamma::CrGamma(const std::string& paramstring)
: m_component(0)
{
  std::vector<float> params;
  // including each component (primary/secondary gamma-rays)...
  m_subComponents.push_back(new CrGammaPrimary);
  // m_subComponents.push_back(new CrGammaSecondaryDownward);
  m_subComponents.push_back(new CrGammaSecondaryUpward);

  m_engine = new HepJamesRandom;
}


// Destructor. Delete each component
CrGamma::~CrGamma()
{
  std::vector<CrSpectrum*>::iterator  i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    delete *i;
  }
}


// Gives back component in the ratio of the flux
CrSpectrum* CrGamma::selectComponent()
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
double CrGamma::energy(double time)
{
  selectComponent();
  return m_component->energySrc(m_engine);
}


// Gives back paticle direction in cos(theta) and phi[rad]
std::pair<double,double> CrGamma::dir(double energy)
{
  if (!m_component){ selectComponent(); }

  return m_component->dir(energy, m_engine);
}


// Gives back the total flux (summation of each component's flux)
G4double CrGamma::flux(double time) const
{
  G4double          total_flux = 0;
  std::vector<CrSpectrum*>::const_iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
  }
  return total_flux;
}

// Gives back solid angle from whick particles come
G4double CrGamma::solidAngle() const
{
  return 4 * pi;
}

// Gives back the interval to the next event
G4double CrGamma::interval(double time){
  return -1.0;
}
