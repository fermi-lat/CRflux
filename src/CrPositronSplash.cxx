/**************************************************************************
 * CrPositronSplash.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrPositron, the entry-point class for the cosmic-ray positron
 * generation.
 **************************************************************************
 * This program generates the cosmic-ray secondary positron upward flux 
 * at satellite altitude with proper angular distribution and energy spectrum.
 * The absolute flux and spectrum of upward positrons are assumed 
 * to depend on the geomagnetic cutoff energy.
 * The flux is assumed not to depend on the zenith angle since AMS
 * didn't detect significant difference between downward and upward
 * flux in theta_M<0.6.
 * The energy spectrum above 100 MeV is represented by simple analytic
 * functions such as a power-low referring to AMS data.
 * Below 100 MeV, the spectrum is extrapolated down to 10 MeV
 * assuming that the flux is proportional to E^-1.
 * A method reentrantCRenergy returns an energy and CrPositronReentrant::dir 
 * returns a direction in cos(theta) and phi (azimuth angle). 
 * Please note that we don't have splash nor reentrant component at
 * satellite altitude. The class is named "**Splash" due to the
 * historical reason.
 **************************************************************************
 * Definitions:
 * 1) The z-axis points upward (from Calorimeter to Tracker).  
 * 2) Partile of theta=0 points for downward (i.e., comes from zenith)
 *    and that of theta=pi points for upward (comes from nadir).
 * 3) Partile of phi=0 comes along x-axis (from x>0 to x=0) and
 *    that of phi=pi/2 comes along y-axis (from y>0 to y=0).
 * 4) Energy means kinetic energy unless defined otherwise.
 * 5) Magnetic latitude theta_M is in radian.
 * 6) Particle direction is defined by cos(theta) and phi (in radian).
 **************************************************************************
 * 2001-04 Written by M. Ozaki (ISAS) and T. Mizuno (Hiroshima Univ.) 
 * 2001-05 Modified by T. Mizuno (Hiroshima Univ.)
 * 2001-05 Comments added and unused codes removed by T. Kamae (SLAC)
 * 2001-05 Program checked by T. Kamae and H. Mizushima (Hiroshima)
 * 2001-11 Modified by T. Mizuno
 *           angular distribution is changed to be uniform
 *           energy spectrum is extrapolated down to 10 MeV
 * 2003-02 Modified by T. Mizuno to generate flux at any position in orbit.
 * 2004-04 Modified by T. Mizuno to simplify the model functions.
 **************************************************************************
 */

//$Header$

#include <cmath>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrPositronSplash.hh"
#include "CrPositronSubSplash.hh"


typedef double G4double;


// private function definitions.
namespace {
  const G4double restE = 5.11e-4; // rest energy of positron in [GeV]

  // gives back v/c as a function of kinetic Energy
  inline G4double beta(G4double E /* GeV */)
  {
#if 0	// if E ~ restE
    return sqrt(1 - pow(E/restE+1, -2));
#else	// if E >> restE
    return 1.0;
#endif
  }


  // gives back the rigidity (p/Ze where p is the momentum, e means
  // positron charge magnitude, and Z is the atomic number) in units of [GV],
  // as a function of kinetic Energy [GeV].
  inline G4double rigidity(G4double E /* GeV */)
  {
#if 0	// if E ~ restE
    return sqrt(pow(E + restE, 2) - pow(restE, 2));
#else	// if E >> restE
    return E;
#endif
  }


  // gives back the kinetic energy (GeV) as a function of rigidity (GV)
  inline G4double energy(G4double rigidity /* GV */)
  {
#if 0	// if energy ~ restE
    return sqrt(pow(rigidity, 2) + pow(restE, 2)) - restE;
#else	// if energy >> restE
    return rigidity;
#endif
  }

} // End of noname-namespace: private function definitions.


//
//
//

CrPositronSplash::CrPositronSplash()
{
  crPositronSplash_0003 = new CrPositronSplash_0003();
  crPositronSplash_0306 = new CrPositronSplash_0306();
  crPositronSplash_0608 = new CrPositronSplash_0608();
  crPositronSplash_0809 = new CrPositronSplash_0809();
  crPositronSplash_0910 = new CrPositronSplash_0910();
  crPositronSplash_1011 = new CrPositronSplash_1011();
}


CrPositronSplash::~CrPositronSplash()
{
  delete crPositronSplash_0003;
  delete crPositronSplash_0306;
  delete crPositronSplash_0608;
  delete crPositronSplash_0809;
  delete crPositronSplash_0910;
  delete crPositronSplash_1011;
}


// Gives back particle direction in (cos(theta), phi)
std::pair<G4double,G4double> CrPositronSplash::dir(G4double energy, 
					       HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0).
{
  double theta = acos(engine->flat()); // theta is from 0 to pi/2
  theta = M_PI - theta;
  double phi = engine->flat() * 2 * M_PI;

  return  std::pair<G4double,G4double>(cos(theta), phi);
}


// Gives back particle energy
G4double CrPositronSplash::energySrc(HepRandomEngine* engine) const
{

  G4double r1, r2;
  if (m_geomagneticLatitude*M_PI/180.0<0.15){
    return crPositronSplash_0003->energy(engine);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.15 && m_geomagneticLatitude*M_PI/180.0<0.45){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.15;
    r2 = 0.45-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0003->energy(engine);
    } else {
      return crPositronSplash_0306->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.45 && m_geomagneticLatitude*M_PI/180.0<0.7){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.45;
    r2 = 0.7-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0306->energy(engine);
    } else {
      return crPositronSplash_0608->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.7 && m_geomagneticLatitude*M_PI/180.0<0.85){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.7;
    r2 = 0.85-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0608->energy(engine);
    } else {
      return crPositronSplash_0809->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.85 && m_geomagneticLatitude*M_PI/180.0<0.95){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.85;
    r2 = 0.95-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0809->energy(engine);
    } else {
      return crPositronSplash_0910->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.95 && m_geomagneticLatitude*M_PI/180.0<1.05){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.95;
    r2 = 1.05-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0910->energy(engine);
    } else {
      return crPositronSplash_1011->energy(engine);
    }
  } else{
    return crPositronSplash_1011->energy(engine);
  }

}


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".
G4double CrPositronSplash::flux() const
{
  // energy integrated vertically upward flux, [c/s/m^2/sr]
  G4double upwardFlux; 
  G4double r1, r2;
  if (m_geomagneticLatitude*M_PI/180.0<0.15){
    upwardFlux = crPositronSplash_0003->upwardFlux();
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.15 && m_geomagneticLatitude*M_PI/180.0<0.45){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.15;
    r2 = 0.45-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crPositronSplash_0003->upwardFlux()
		   +r1*crPositronSplash_0306->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.45 && m_geomagneticLatitude*M_PI/180.0<0.7){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.45;
    r2 = 0.7-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crPositronSplash_0306->upwardFlux()
                     +r1*crPositronSplash_0608->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.7 && m_geomagneticLatitude*M_PI/180.0<0.85){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.7;
    r2 = 0.85-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crPositronSplash_0608->upwardFlux()
		   +r1*crPositronSplash_0809->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.85 && m_geomagneticLatitude*M_PI/180.0<0.95){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.85;
    r2 = 0.95-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crPositronSplash_0809->upwardFlux()
                     +r1*crPositronSplash_0910->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.95 && m_geomagneticLatitude*M_PI/180.0<1.05){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.95;
    r2 = 1.05-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crPositronSplash_0910->upwardFlux()
                     +r1*crPositronSplash_1011->upwardFlux() )/(r1+r2);
  } else{
    upwardFlux = crPositronSplash_1011->upwardFlux();
  }

  return upwardFlux; // [c/s/m^2/sr]

}


// Gives back solid angle from which particle comes
G4double CrPositronSplash::solidAngle() const
{
  return 2 * M_PI;
}


// Gives back particle name
const char* CrPositronSplash::particleName() const
{
  return "e+";
}


// Gives back the title of the component
std::string CrPositronSplash::title() const
{
  return  "CrPositronSplash";
}


