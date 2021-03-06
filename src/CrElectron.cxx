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
 *         (User no longer need Geant4 to generate cosmic-ray particles)
 *         Instead of CosmicRayGeneratorAction of BalloonTestV13, 
 *         CrGenerator.cxx provides user with an interface 
 *         to all particle models and is used in end-to-end simulation
 *         framework of Tune's group. User can alternatively use 
 *         this CrElectron.cxx to generate cosmic-ray electrons and this class
 *         is used in CRflux package in SLAC CVS. 
 ****************************************************************************
 */

//$Header$

#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>

// CLHEP
//#include <CLHEP/config/CLHEP.h>
//#include <CLHEP/Random/JamesRandom.h>
#include "CLHEP/Random/Random.h"


#include "CrElectron.hh"
#include "CrElectronPrimary.hh"
#include "CrElectronSplash.hh"
#include "CrElectronReentrant.hh"

#include "CrSpectrum.hh"

typedef double G4double;

// Constructor. Includes each component
CrElectron::CrElectron(const std::string& paramstring)
: m_component(0)
{
  std::vector<float> params;
  //use parseParamList to parse out the input string
  parseParamList(paramstring,params);
  //the first element in the string is the bit field.(defaults to "all on")
  int flag = params.empty() || params[0]==0 ? 7 : params[0];
  // including each component if it is present in the bit field...
  if(flag& 1) m_subComponents.push_back(new CrElectronPrimary);
  if(flag& 2) m_subComponents.push_back(new CrElectronReentrant);
  if(flag& 4) m_subComponents.push_back(new CrElectronSplash);

  if(params.size()>1 && params[1]>0) {   
     std::vector<CrSpectrum*>::iterator i;
     for (i=m_subComponents.begin() ; i != m_subComponents.end(); i++){
	(*i)->setNormalization(params[1]);
     };
  };	   

  m_engine = CLHEP::HepRandom::getTheEngine(); //new HepJamesRandom;
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
  std::map<CrSpectrum*,G4double> integ_flux;
  G4double                       total_flux = 0;
  std::vector<CrSpectrum*>::iterator i;

  // calculate the total flux
  if(m_subComponents.size() > 1)
  {
     for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
        total_flux += (*i)->solidAngle()*(*i)->flux();
        integ_flux[*i] = total_flux; // Don't divide by 4 pi it's not needed here
     }
  }
  else
  {
     for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
        total_flux += (*i)->flux();
        integ_flux[*i] = total_flux;
     }
  }
  // select component based on the flux
  G4double  rnum = m_engine->flat() * total_flux;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    if (integ_flux[*i] >= rnum) { break; }
  }

  m_component = *i;

  return *i;
}


// Gives back energy
G4double CrElectron::energy(G4double /* time */)
{
  selectComponent();
  return m_component->energySrc(m_engine);
}


// Gives back paticle direction in cos(theta) and phi[rad]
std::pair<G4double,G4double> CrElectron::dir(G4double energy)
{
  if (!m_component){ selectComponent(); }

  return m_component->dir(energy, m_engine);
}


// Gives back the total flux (summation of each component's flux)
G4double CrElectron::flux(G4double /* time */) const
{
  G4double          total_flux = 0;
  std::vector<CrSpectrum*>::const_iterator i;
  // if summing over several sources scale by solid angles
  if(m_subComponents.size() > 1)
  {
     for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
        total_flux += (*i)->solidAngle()*(*i)->flux();
     }
     return total_flux / (4 * M_PI);
  }
  else
  {
     for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
        total_flux += (*i)->flux();
        //    cout << "flux of this component = " << (*i)->flux() << endl;
     }
     return total_flux;
  }
}

// Gives back solid angle from which particles come
G4double CrElectron::solidAngle() const
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
void CrElectron::dump()
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
G4double CrElectron::interval(G4double /* time */){
  return -1.0;
}

void CrElectron::parseParamList(std::string input, std::vector<float>& output)
{  
  int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    float f = ::atof( input.c_str() );
    output.push_back(f);
    i=input.find_first_of(",");
    input= input.substr(i+1);
  } 
}


