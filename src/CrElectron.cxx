/**************************************************************************
 * CrElectron.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program interfaces 3 spectra generators, CrElectronPrimary.cc
 * CrElectronReentrant.cc, and CrElectronSplash.cc to 
 * CosmicRayGeneratorAction.cc.
 **************************************************************************
 * This program defines the entry-point class for the cosmic ray electron 
 * generation and interfaces to the primary, reentrant, and splash
 * cosmic-ray electron generators.
 **************************************************************************
 * 2001-04 Written by M. Ozaki (ISAS).
 * 2001-05 Modified by T. Mizuno (Hiroshima Univ.)
 * 2001-05 Comments added and unused codes removed by T. Kamae (SLAC)
 * 2001-05 Program checked by T. Kamae and H. Mizushima (Hiroshima)
 * 2001-11 Modified by T. Mizuno to construct a `stand-alone' module
 ****************************************************************************
 */

//$Header$

#include <math.h>
#include <map>
#include <vector>

// CLHEP
#include <CLHEP/Random/JamesRandom.h>

#include "CrElectron.hh"
#include "CrElectronPrimary.hh"
#include "CrElectronSplash.hh"
#include "CrElectronReentrant.hh"

#include "CrSpectrum.hh"

// define a factory for anonomous instantiation
#include "FluxSvc/ISpectrumFactory.h"

typedef  double G4double;
namespace{
  const G4double pi    = 3.14159265358979323846264339;
}

// Constructor. Includes each component
CrElectron::CrElectron(const std::string& paramstring)
: m_component(0)
{
  std::vector<float> params;
  // including each component (primary/re-entrant/splash electrons)...
  m_subComponents.push_back(new CrElectronPrimary);
  m_subComponents.push_back(new CrElectronReentrant);
  m_subComponents.push_back(new CrElectronSplash);

  m_engine = new HepJamesRandom;
}


// Destructor. Delete each component
CrElectron::~CrElectron()
{
  std::vector<CrSpectrum*>::iterator  i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    delete *i;
  }
}


// Gives back component in the ratio of the flux
CrSpectrum* CrElectron::selectComponent()
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
double CrElectron::energy(double time)
{
  selectComponent();
  return m_component->energySrc(m_engine);
}


// Gives back paticle direction in cos(theta) and phi[rad]
std::pair<double,double> CrElectron::dir(double energy)
{
  if (!m_component){ selectComponent(); }

  return m_component->dir(energy, m_engine);
}


// Gives back the total flux (summation of each component's flux)
G4double CrElectron::flux(double time) const
{
  G4double          total_flux = 0;
  std::vector<CrSpectrum*>::const_iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
  }
  return total_flux;
}

// Gives back solid angle from which particles come
G4double CrElectron::solidAngle() const
{
  return 4 * pi;
}

// Gives back the interval to the next event
G4double CrElectron::interval(double time){
  return -1.0;
}


