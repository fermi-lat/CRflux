/**************************************************************************
 * CrProtonSplash.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrProton, the entry-point class for the cosmic-ray proton 
 * generation.
 **************************************************************************
 * This program generates the cosmic-ray secondary proton upward flux 
 * at satellite altitude with proper angular distribution and energy spectrum.
 * The absolute flux and spectrum of upward protons are assumed 
 * to depend on the geomagnetic cutoff energy.
 * The flux is assumed not to depend on the zenith angle since AMS
 * didn't detect significant difference between downward and upward
 * flux in theta_M<0.6.
 * The energy spectrum above 100 MeV is represented by simple analytic
 * functions such as a power-low referring to AMS data.
 * Below 100 MeV, the spectrum is extrapolated down to 10 MeV
 * assuming that the flux is proportional to E^-1.
 * A method reentrantCRenergy returns an energy and CrProtonReentrant::dir 
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
 * in 2000 Originally written by M. Ozaki (ISAS) 
 *         and K. Hirano (Hiroshima Univ.)
 * 2001-04 Spectrum part modified by T. Mizuno (Hiroshima Univ.)
 * 2001-05 Spectrum part modified by T. Mizuno (Hiroshima Univ.)
 * 2001-05 Comments added and unused codes removed by T. Kamae (SLAC)
 * 2001-05 Program checked by T. Kamae and H. Mizushima (Hiroshima)
 * 2001-11 Modified by T. Mizuno.
 *           angular distribution is changed to be uniform
 *           energy spectrum is extrapolated down to 10 MeV
 * 2001-11 Modified by T. Mizuno to construct a `stand-alone' module
 * 2003-02 Modified by T. Mizuno to generate flux at any position in orbit.
 * 2004-04 Modified by T. Mizuno to simplify the model functions.
 * 2005-05 Modified by T. Mizuno to calculate the flux when theta_M<0.
 **************************************************************************
 */

// $Header$

#include <cmath>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrProtonSplash.hh"
#include "CrProtonSubSplash.hh"


typedef double G4double;


// private function definitions.
namespace {
  // rest energy (rest mass) of proton in units of GeV
  const G4double restE = 0.938;

  // gives back v/c as a function of kinetic Energy
  inline G4double beta(G4double E /* GeV */){
    return sqrt(1 - pow(E/restE+1, -2));
  }

  // gives back the rigidity (p/Ze where p is the momentum, e means
  // electron charge magnitude, and Z is the atomic number) in units of [GV],
  // as a function of kinetic Energy [GeV].
  inline G4double rigidity(G4double E /* GeV */){
    return sqrt(pow(E + restE, 2) - pow(restE, 2));
  }

  // gives back the kinetic energy (GeV) as a function of rigidity (GV)
  inline G4double energy(G4double rigidity /* GV */){
    return sqrt(pow(rigidity, 2) + pow(restE, 2)) - restE;
  }

} // End of noname-namespace: private function definitions.


//
//
//

CrProtonSplash::CrProtonSplash()
{
  crProtonSplash_0002 = new CrProtonSplash_0002();
  crProtonSplash_0203 = new CrProtonSplash_0203();
  crProtonSplash_0304 = new CrProtonSplash_0304();
  crProtonSplash_0405 = new CrProtonSplash_0405();
  crProtonSplash_0506 = new CrProtonSplash_0506();
  crProtonSplash_0607 = new CrProtonSplash_0607();
  crProtonSplash_0708 = new CrProtonSplash_0708();
  crProtonSplash_0809 = new CrProtonSplash_0809();
  crProtonSplash_0910 = new CrProtonSplash_0910();
}


CrProtonSplash::~CrProtonSplash()
{
  delete crProtonSplash_0002;
  delete crProtonSplash_0203;
  delete crProtonSplash_0304;
  delete crProtonSplash_0405;
  delete crProtonSplash_0506;
  delete crProtonSplash_0607;
  delete crProtonSplash_0708;
  delete crProtonSplash_0809;
  delete crProtonSplash_0910;
}


// Gives back particle direction in (cos(theta), phi)
std::pair<double,double> CrProtonSplash::dir(double energy, 
					     HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward has plus sign in cos(theta),
  // and phi = 0 for particle comming along x-axis (from x>0 to x=0)
  // and phi = pi/2 for that comming along y-axis (from y>0 to y=0)
{
  double theta = acos(engine->flat()); // theta is from 0 to pi/2
  theta = M_PI - theta;
  double phi = engine->flat() * 2 * M_PI;

  return std::pair<double,double>(cos(theta), phi);
}


// Gives back particle energy
double CrProtonSplash::energySrc(HepRandomEngine* engine) const
{
  double r1, r2;
  if (fabs(m_geomagneticLatitude)*M_PI/180.0<0.15){
    return crProtonSplash_0002->energy(engine);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.15 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.25){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.15;
    r2 = 0.25-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0002->energy(engine);
    } else {
      return crProtonSplash_0203->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.25 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.35){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.25;
    r2 = 0.35-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0203->energy(engine);
    } else {
      return crProtonSplash_0304->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.35 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.45){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.35;
    r2 = 0.45-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0304->energy(engine);
    } else {
      return crProtonSplash_0405->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.45 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.55){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.45;
    r2 = 0.55-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0405->energy(engine);
    } else {
      return crProtonSplash_0506->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.55 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.65){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.55;
    r2 = 0.65-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0506->energy(engine);
    } else {
      return crProtonSplash_0607->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.65 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.75){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.65;
    r2 = 0.75-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0607->energy(engine);
    } else {
      return crProtonSplash_0708->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.75 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.85){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.75;
    r2 = 0.85-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0708->energy(engine);
    } else {
      return crProtonSplash_0809->energy(engine);
    }
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.85 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.95){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.85;
    r2 = 0.95-fabs(m_geomagneticLatitude)*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0809->energy(engine);
    } else {
      return crProtonSplash_0910->energy(engine);
    }
  } else {
    return crProtonSplash_0910->energy(engine);
  }

}


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".
double CrProtonSplash::flux() const
{
  // energy integrated vertically upward flux, [c/s/m^2/sr]
  double upwardFlux; 
  double r1, r2;
  if (fabs(m_geomagneticLatitude)*M_PI/180.0<0.15){
    upwardFlux = crProtonSplash_0002->upwardFlux();
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.15 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.25){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.15;
    r2 = 0.25-fabs(m_geomagneticLatitude)*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0002->upwardFlux()
		     +r1*crProtonSplash_0203->upwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.25 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.35){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.25;
    r2 = 0.35-fabs(m_geomagneticLatitude)*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0203->upwardFlux()
		     +r1*crProtonSplash_0304->upwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.35 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.45){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.35;
    r2 = 0.45-fabs(m_geomagneticLatitude)*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0304->upwardFlux()
		     +r1*crProtonSplash_0405->upwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.45 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.55){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.45;
    r2 = 0.55-fabs(m_geomagneticLatitude)*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0405->upwardFlux()
		     +r1*crProtonSplash_0506->upwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.55 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.65){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.55;
    r2 = 0.65-fabs(m_geomagneticLatitude)*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0506->upwardFlux()
		     +r1*crProtonSplash_0607->upwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.65 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.75){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.65;
    r2 = 0.75-fabs(m_geomagneticLatitude)*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0607->upwardFlux()
		     +r1*crProtonSplash_0708->upwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.75 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.85){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.75;
    r2 = 0.85-fabs(m_geomagneticLatitude)*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0708->upwardFlux()
		     +r1*crProtonSplash_0809->upwardFlux() )/(r1+r2);
  } else if (fabs(m_geomagneticLatitude)*M_PI/180.0>=0.85 && fabs(m_geomagneticLatitude)*M_PI/180.0<0.95){
    r1 = fabs(m_geomagneticLatitude)*M_PI/180.0-0.85;
    r2 = 0.95-fabs(m_geomagneticLatitude)*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0809->upwardFlux()
		     +r1*crProtonSplash_0910->upwardFlux() )/(r1+r2);
  } else {
    upwardFlux = crProtonSplash_0910->upwardFlux();
  }

  return upwardFlux; // [c/s/m^2/sr]

}


// Gives back solid angle from which particle comes
double CrProtonSplash::solidAngle() const
{
  return 2 * M_PI;
}


// Gives back particle name
const char* CrProtonSplash::particleName() const
{
  return "proton";
}


// Gives back the name of the component
std::string CrProtonSplash::title() const
{
  return  "CrProtonSplash";
}

