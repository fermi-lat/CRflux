/****************************************************************************
 * CrNeutron.cc
 ****************************************************************************
 * Atmospheric neutron generator for GLAST LAT.
 ****************************************************************************
 * 2007-09 Written by T. Mizuno (Hiroshima Univ.)
 ****************************************************************************
 */

//$Header$

#include <cmath>
#include <map>
#include <vector>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/Random.h>

#include "CrNeutron.hh"
#include "CrNeutronSplash.hh"

#include "CrSpectrum.hh"


typedef double G4double;

// Constructor. Includes each component
CrNeutron::CrNeutron(const std::string& /* paramstring */)
: m_component(0)
{
  std::vector<float> params;
  // including each component (splash alphas)...
  m_subComponents.push_back(new CrNeutronSplash);
  
   m_engine = HepRandom::getTheEngine(); //new HepJamesRandom;
}


// Destructor. Delete each component
CrNeutron::~CrNeutron()
{
  std::vector<CrSpectrum*>::iterator  i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    delete *i;
  }
}


// Gives back component in the ratio of the flux
CrSpectrum* CrNeutron::selectComponent()
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
G4double CrNeutron::energy(double /* time */)
{
  selectComponent();
  return m_component->energySrc(m_engine);
}


// Gives back paticle direction in cos(theta) and phi[rad]
std::pair<G4double,G4double> CrNeutron::dir(G4double energy)
{
  if (!m_component){ selectComponent(); }

  return m_component->dir(energy, m_engine);
}

// Gives back the total flux (summation of each component's flux)
G4double CrNeutron::flux(G4double /* time */ ) const
{
  G4double total_flux = 0;
  std::vector<CrSpectrum*>::const_iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
    //    cout << "flux of this component = " << (*i)->flux() << endl;
  }
  return total_flux;
}

// Gives back solid angle from whick particles come
G4double CrNeutron::solidAngle() const
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
void CrNeutron::dump()
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
G4double CrNeutron::interval(G4double /* time */){
  return -1.0;
}


