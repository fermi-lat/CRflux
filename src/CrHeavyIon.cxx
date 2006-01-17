/**************************************************************************
 * CrHeavyIon.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program interfaces CrExamplePrimary.cc to CosmicRayGeneratorAction.cc 
 **************************************************************************
 * This program defines the entry-point class for the cosmic ray alpha 
 * generation and interfaces to the primary cosmic-ray alpha generators.
 **************************************************************************
 * 2002-04 Written by Y. Fukazawa (Hiroshima Univ.)
 * 2003-02 Modified by T. Mizuno to construct a `stand-alone' module
 *         (User no longer need Geant4 to generate cosmic-ray particles)
 *         Instead of CosmicRayGeneratorAction of BalloonTestV13, 
 *         CrGenerator.cxx provides user with an interface 
 *         to all particle models and is used in end-to-end simulation
 *         framework of Tune's group. User can alternatively use 
 *         this CrHeavyIon.cxx to generate cosmic-ray alphas and this class
 *         is used in CRflux package in SLAC CVS. 
 ****************************************************************************
 */

//$Header$

#include <math.h>
#include <map>
#include <vector>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include "CLHEP/Random/Random.h"
//#include <CLHEP/Random/JamesRandom.h>

#include "CrHeavyIon.h"
#include "CrHeavyIonPrimary.hh"

#include "CrSpectrum.hh"

// Define a factory for anonomous instantiation.
// Comment-out the next line when importing into the CRflux package.
//#include "FluxSvc/ISpectrumFactory.h"

typedef double G4double;

// Constructor. Includes each component
CrHeavyIon::CrHeavyIon(const std::string& paramstring)
: m_component(0)
{
  std::vector<float> params;
  // including each component (primary alphas)...
  m_subComponents.push_back(new CrHeavyIonPrimary);
  
  m_engine = HepRandom::getTheEngine(); //new HepJamesRandom;
}


// Destructor. Delete each component
 CrHeavyIon::~CrHeavyIon()
{
  std::vector<CrSpectrum*>::iterator  i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    delete *i;
  }
}


// Gives back component in the ratio of the flux
CrSpectrum* CrHeavyIon::selectComponent()
{
  std::map<CrSpectrum*,G4double> integ_flux;
  double total_flux = 0;
  std::vector<CrSpectrum*>::iterator i;

  // calculate the total flux 
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
    integ_flux[*i] = total_flux;
  }
  // select component based on the flux
  G4double rnum = m_engine->flat() * total_flux;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    if (integ_flux[*i] >= rnum) { break; }
  }

  m_component = *i;

  return *i;
}

// Gives back kinetic energy 
G4double CrHeavyIon::energy(double time)
{
  selectComponent();
  return m_component->energySrc(m_engine);
}


// Gives back paticle direction in cos(theta) and phi[rad]
std::pair<G4double,G4double> CrHeavyIon::dir(G4double energy)
{
  if (!m_component){ selectComponent(); }

  return m_component->dir(energy, m_engine);
}

// Gives back the total flux (summation of each component's flux)
G4double CrHeavyIon::flux(G4double time) const
{
  G4double total_flux = 0;
  std::vector<CrSpectrum*>::const_iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
    //    std::cout << "flux of this component = " << (*i)->flux() << std::endl;
  }
  return total_flux;
}

// Gives back solid angle from whick particles come
G4double CrHeavyIon::solidAngle() const
{
   if(m_subComponents.size() == 1)
   {
      std::vector<CrSpectrum*>::const_iterator i;
      i = m_subComponents.begin();
      return (*i)->solidAngle();
   }
   else
      return 4 * M_PI;
}

// print out the information of each component
void CrHeavyIon::dump()
{
  std::vector<CrSpectrum*>::const_iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    std::cout << "title: " << (*i)->title() << std::endl;
    std::cout << " flux(c/s/m^2/sr)= " << (*i)->flux() << std::endl;
    std::cout << " geographic latitude/longitude(deg)= " 
         << (*i)->latitude() << " " << (*i)->longitude() << std::endl;
    std::cout << " geomagnetic latitude/longitude(deg)= " 
         << (*i)->geomagneticLatitude() << " " 
         << (*i)->geomagneticLongitude() << std::endl;
    std::cout << " time(s)= " << (*i)->time() 
         << " altitude(km)= " << (*i)->altitude() << std::endl;
    std::cout << " cor(GV)= " << (*i)->cutOffRigidity() 
         << " phi(MV)= " << (*i)->solarWindPotential() << std::endl;
  }
}

// Gives back the interval to the next event.
// If interval(time) is negative, FluxSvc will call
// flux(time) to determine the average flux for that time
// and calculate the arrival time for the next particle using
// the poission distribution.
G4double CrHeavyIon::interval(G4double time){
  return -1.0;
}

const char* CrHeavyIon::particleName() const {
  return m_component->particleName();
}

