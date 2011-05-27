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
 * This program generates the cosmic-ray secondary positron downward flux 
 * at satellite altitude with proper angular distribution and energy spectrum.
 * The absolute flux and spectrum of downward positrons are assumed 
 * to depend on the geomagnetic cutoff energy.
 * The flux is assumed not to depend on the zenith angle since AMS
 * didn't detect significant difference between downward and upward flux.
 * The energy spectrum above 100 MeV is represented by simple analytic
 * functions such as a power-low referring to AMS data.
 * Below 100 MeV, the spectrum is extrapolated down to 10 MeV
 * assuming that the flux is proportional to E^-1.
 * A method splashCRenergy returns an energy and CrPositronSplash::dir 
 * returns a direction in cos(theta) and phi (azimuth angle). 
 * Please note that we don't have splash nor splash component at
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
 * 2001-11 Modified by T. Mizuno.
 *           angular distribution is changed to be uniform
 *           energy spectrum is extrapolated down to 10 MeV
 * 2001-12 Modified by T. Mizuno to construct a `stand-alone' module
 * 2003-02 Modified by T. Mizuno to generate flux at any position in orbit.
 * 2004-04 Modified by T. Mizuno to simplify the model functions.
 * 2005-05 Modified by T. Mizuno to calculate the flux when theta_M<0.
 * 2010-10 Modified by T. Mizuno to adjust the flux to that measured by LAT.
 *   (analysis by Melissa).
 **************************************************************************
 */

//$Header$


#include <cmath>

// CLHEP
//#include <CLHEP/config/CLHEP.h>
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


  // gives back the rigidity (p/Ze where p is the momentum and e is
  // the unit charge and Z the atomic number) in units of [GV],
  // as a function of kinetic energy [GeV].
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

CrPositronSplash::CrPositronSplash():CrSpectrum()
{
  crPositronSplash_0001 = new CrPositronSplash_0001();
  crPositronSplash_0102 = new CrPositronSplash_0102();
  crPositronSplash_0203 = new CrPositronSplash_0203();
  crPositronSplash_0304 = new CrPositronSplash_0304();
  crPositronSplash_0405 = new CrPositronSplash_0405();
  crPositronSplash_0506 = new CrPositronSplash_0506();
  crPositronSplash_0611 = new CrPositronSplash_0611();
}


CrPositronSplash::~CrPositronSplash()
{
  delete crPositronSplash_0001;
  delete crPositronSplash_0102;
  delete crPositronSplash_0203;
  delete crPositronSplash_0304;
  delete crPositronSplash_0405;
  delete crPositronSplash_0506;
  delete crPositronSplash_0611;
}


// Gives back particle direction in (cos(theta), phi)
std::pair<G4double,G4double> CrPositronSplash::dir
(G4double energy, CLHEP::HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0).
{
  double theta = acos(engine->flat()); // theta is from 0 to pi/2
  theta = M_PI - theta;
  G4double phi = engine->flat() * 2 * M_PI;

  return std::pair<G4double,G4double>(cos(theta), phi);
}


// Gives back particle energy
G4double CrPositronSplash::energySrc(CLHEP::HepRandomEngine* engine) const
{
  G4double r1, r2;
  if (fabs(m_geomagneticLatitude)*M_PI/180.0<0.05){
    return crPositronSplash_0001->energy(engine);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.05 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.15){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.05;
    r2 = 0.15-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0001->energy(engine);
    } else {
      return crPositronSplash_0102->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.15 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.25){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.15;
    r2 = 0.25-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0102->energy(engine);
    } else {
      return crPositronSplash_0203->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.25 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.35){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.25;
    r2 = 0.35-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0203->energy(engine);
    } else {
      return crPositronSplash_0304->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.35 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.45){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.35;
    r2 = 0.45-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0304->energy(engine);
    } else {
      return crPositronSplash_0405->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.45 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.55){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.45;
    r2 = 0.55-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0405->energy(engine);
    } else {
      return crPositronSplash_0506->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.55 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.65){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.55;
    r2 = 0.65-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0506->energy(engine);
    } else {
      return crPositronSplash_0611->energy(engine);
    }
  } else{
      return crPositronSplash_0611->energy(engine);
  }

}


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "splash" and "splash".
G4double CrPositronSplash::flux() const
{
  // energy integrated vertically downward flux, [c/s/m^2/sr]
  G4double downwardFlux; 
  G4double r1, r2;
  if (fabs(m_geomagneticLatitude)*M_PI/180.0<0.05){
    downwardFlux = crPositronSplash_0001->downwardFlux();
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.05 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.15){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.05;
    r2 = 0.15-fabs(m_geomagneticLatitude)*M_PI/180.0;
    downwardFlux = ( r2*crPositronSplash_0001->downwardFlux()
             +r1*crPositronSplash_0102->downwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.15 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.25){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.15;
    r2 = 0.25-fabs(m_geomagneticLatitude)*M_PI/180.0;
    downwardFlux = ( r2*crPositronSplash_0102->downwardFlux()
             +r1*crPositronSplash_0203->downwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.25 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.35){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.25;
    r2 = 0.35-fabs(m_geomagneticLatitude)*M_PI/180.0;
    downwardFlux = ( r2*crPositronSplash_0203->downwardFlux()
             +r1*crPositronSplash_0304->downwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.35 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.45){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.35;
    r2 = 0.45-fabs(m_geomagneticLatitude)*M_PI/180.0;
    downwardFlux = ( r2*crPositronSplash_0304->downwardFlux()
             +r1*crPositronSplash_0405->downwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.45 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.55){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.45;
    r2 = 0.55-fabs(m_geomagneticLatitude)*M_PI/180.0;
    downwardFlux = ( r2*crPositronSplash_0405->downwardFlux()
             +r1*crPositronSplash_0506->downwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.55 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.65){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.55;
    r2 = 0.65-fabs(m_geomagneticLatitude)*M_PI/180.0;
    downwardFlux = ( r2*crPositronSplash_0506->downwardFlux()
             +r1*crPositronSplash_0611->downwardFlux() )/(r1+r2);
  } else{
    downwardFlux = crPositronSplash_0611->downwardFlux();
  }

  return m_normalization*downwardFlux; // [c/s/m^2/sr]


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


// Gives back the name of the component
std::string CrPositronSplash::title() const
{
  return  "CrPositronSplash";
}


