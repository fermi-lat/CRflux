/**************************************************************************
 * CrProton.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program interfaces 3 spectra generators, CrProtonPrimary.cc,
 * CrProtonReentrant.cc, and CrProtonSplash.cc to 
 * CosmicRayGeneratorAction.cc 
 **************************************************************************
 * This program defines the entry-point class for the cosmic ray proton 
 * generation and interfaces to the primary, reentrant, and splash
 * cosmic-ray proton generators.
 **************************************************************************
 * 2001-04 Written by M. Ozaki (ISAS) and T. Mizuno (Hiroshima Univ.)
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

#include "CrProton.hh"
#include "CrProtonPrimary.hh"
#include "CrProtonSplash.hh"
#include "CrProtonReentrant.hh"

#include "CrSpectrum.hh"

// define a factory for anonomous instantiation
#include "FluxSvc/ISpectrumFactory.h"

typedef  double G4double;
namespace{
  const G4double pi    = 3.14159265358979323846264339;
}

// Constructor. Includes each component
CrProton::CrProton(const std::string& paramstring)
: m_component(0)
{   
  std::vector<float> params;
  //use parseParamList to parse out the input string
  parseParamList(paramstring,params);
  //the first element in the string is the bit field.(defaults to "all on")
  int flag = params.empty() || params[0]==0 ? 7 : params[0];
  // including each component if it is present in the bit field...
  if(flag& 1) m_subComponents.push_back(new CrProtonPrimary);
  if(flag& 2) m_subComponents.push_back(new CrProtonReentrant);
  if(flag& 4) m_subComponents.push_back(new CrProtonSplash);
  
  m_engine = new HepJamesRandom;
}


// Destructor. Delete each component
CrProton::~CrProton()
{
  std::vector<CrSpectrum*>::iterator  i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    delete *i;
  }
}


// Gives back component in the ratio of the flux
CrSpectrum* CrProton::selectComponent()
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
double CrProton::energy(double time)
{
  selectComponent();
  return m_component->energySrc(m_engine);
}


// Gives back paticle direction in cos(theta) and phi[rad]
std::pair<double,double> CrProton::dir(double energy)
{
  if (!m_component){ selectComponent(); }

  return m_component->dir(energy, m_engine);
}


// Gives back the total flux (summation of each component's flux)
G4double CrProton::flux(double time) const
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
G4double CrProton::solidAngle() const
{
  return 4 * pi;
}

// Gives back the interval to the next event
G4double CrProton::interval(double time){
  return -1.0;
}

void CrProton::parseParamList(std::string input, std::vector<float>& output)
{  
  int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    float f = ::atof( input.c_str() );
    output.push_back(f);
    i=input.find_first_of(",");
    input= input.substr(i+1);
  } 
}
