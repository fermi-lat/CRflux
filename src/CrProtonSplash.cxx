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
 * This program generates the splash cosmic ray proton flux with
 * proper angular distribution and energy spectrum.
 * The absolute flux and spectrum of "splash" protons are assumed 
 * to depend on the geomagnetic cutoff energy (and is fixed for Palestine,
 * Texas in the current codes). The flux is assumed to 
 * depend on zenith angle as 
 *   1 + 0.6*sin(theta) for theta = pi/2 to pi, and zero for theta < pi/2
 * (theta = pi - zenith angle).
 * The energy spectrum is assumed to be of a power-low with exponential
 * cut-off at a low energy. 
 * A method splashCRenergy returns an energy and CrProtonSplash::dir 
 * returns a direction in cos(theta) and phi (azimuth angle).
 * Energy and angular distribution of splash proton is assumed to be
 * the same as those of reentrant proton.
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
 **************************************************************************
 */

// $Header$

#include <math.h>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrProtonSplash.hh"
#include "CrProtonSubSplash.hh"


typedef  double G4double;


// private function definitions.
namespace {
  const G4double pi = 3.14159265358979323846264339;
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
  /***
   * Here we assume that
   *  theta (= pi - zenith angle) dependence: 
   *    the flux (per steradian) is proportional to 1 + 0.6*sin(theta)
   *  azimuth angle (phi) dependence: isotropic
   *
   * Reference:
   *   Tylka, A. J. 2000-05-12, GLAST team internal report 
   *   "A Review of Cosmic-Ray Albedo Studies: 1949-1970" (1 + 0.6 sin(theta))
   */

  double theta;
  while (1){
    theta = acos(engine->flat()); // theta is from 0 to pi/2
    if (engine->flat()*1.6<1+0.6*sin(theta)){break;}
  }
  theta = pi - theta;

  double phi = engine->flat() * 2 * pi;

  return std::pair<double,double>(cos(theta), phi);
}


// Gives back particle energy
double CrProtonSplash::energySrc(HepRandomEngine* engine) const
{
  double r1, r2;
  if (m_geomagneticLatitude*M_PI/180.0<0.15){
    return crProtonSplash_0002->energy(engine);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.15 && m_geomagneticLatitude*M_PI/180.0<0.25){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.15;
    r2 = 0.25-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0002->energy(engine);
    } else {
      return crProtonSplash_0203->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.25 && m_geomagneticLatitude*M_PI/180.0<0.35){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.25;
    r2 = 0.35-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0203->energy(engine);
    } else {
      return crProtonSplash_0304->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.35 && m_geomagneticLatitude*M_PI/180.0<0.45){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.35;
    r2 = 0.45-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0304->energy(engine);
    } else {
      return crProtonSplash_0405->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.45 && m_geomagneticLatitude*M_PI/180.0<0.55){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.45;
    r2 = 0.55-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0405->energy(engine);
    } else {
      return crProtonSplash_0506->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.55 && m_geomagneticLatitude*M_PI/180.0<0.65){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.55;
    r2 = 0.65-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0506->energy(engine);
    } else {
      return crProtonSplash_0607->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.65 && m_geomagneticLatitude*M_PI/180.0<0.75){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.65;
    r2 = 0.75-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0607->energy(engine);
    } else {
      return crProtonSplash_0708->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.75 && m_geomagneticLatitude*M_PI/180.0<0.85){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.75;
    r2 = 0.85-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0708->energy(engine);
    } else {
      return crProtonSplash_0809->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.85 && m_geomagneticLatitude*M_PI/180.0<0.95){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.85;
    r2 = 0.95-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonSplash_0809->energy(engine);
    } else {
      return crProtonSplash_0910->energy(engine);
    }
  } else {
    return crProtonSplash_0910->energy(engine);
  }

}


// flux() returns the value integrated over whole energy and direction
// and devided by 4pi sr: then the unit is [c/s/m^2/sr].
// This value is used as relative normalization among 
// "primary", "reentrant" and "splash".
double CrProtonSplash::flux() const
{
  // Averaged energy-integrated flux over 4 pi solod angle used 
  // in relative normalizaiotn among "primary", "reentrant" and "splash".

  // energy integrated vertically upward flux, [c/s/m^2/sr]
  double upwardFlux; 
  double r1, r2;
  if (m_geomagneticLatitude*M_PI/180.0<0.15){
    upwardFlux = crProtonSplash_0002->upwardFlux();
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.15 && m_geomagneticLatitude*M_PI/180.0<0.25){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.15;
    r2 = 0.25-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0002->upwardFlux()
		     +r1*crProtonSplash_0203->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.25 && m_geomagneticLatitude*M_PI/180.0<0.35){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.25;
    r2 = 0.35-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0203->upwardFlux()
		     +r1*crProtonSplash_0304->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.35 && m_geomagneticLatitude*M_PI/180.0<0.45){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.35;
    r2 = 0.45-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0304->upwardFlux()
		     +r1*crProtonSplash_0405->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.45 && m_geomagneticLatitude*M_PI/180.0<0.55){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.45;
    r2 = 0.55-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0405->upwardFlux()
		     +r1*crProtonSplash_0506->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.55 && m_geomagneticLatitude*M_PI/180.0<0.65){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.55;
    r2 = 0.65-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0506->upwardFlux()
		     +r1*crProtonSplash_0607->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.65 && m_geomagneticLatitude*M_PI/180.0<0.75){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.65;
    r2 = 0.75-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0607->upwardFlux()
		     +r1*crProtonSplash_0708->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.75 && m_geomagneticLatitude*M_PI/180.0<0.85){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.75;
    r2 = 0.85-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0708->upwardFlux()
		     +r1*crProtonSplash_0809->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.85 && m_geomagneticLatitude*M_PI/180.0<0.95){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.85;
    r2 = 0.95-m_geomagneticLatitude*M_PI/180.0;
    upwardFlux = ( r2*crProtonSplash_0809->upwardFlux()
		     +r1*crProtonSplash_0910->upwardFlux() )/(r1+r2);
  } else {
    upwardFlux = crProtonSplash_0910->upwardFlux();
  }

  // We have assumed that the flux is proportional to 1+0.6 sin(theta)
  // Then, the flux integrated over solid angle is
  // (1+0.15pi)*upwardFlux*2pi
  return 0.5 * (1 + 0.15*pi) * upwardFlux; // [c/s/m^2/sr]

}


// Gives back solid angle from which particle comes
double CrProtonSplash::solidAngle() const
{
  return 2 * pi;
}


// Gives back particle name
const char* CrProtonSplash::particleName() const
{
  return "proton";
}


float CrProtonSplash::operator()(float r)
{
  HepJamesRandom  engine;
  engine.setSeed(r * 900000000);
  // 900000000 comes from HepJamesRandom::setSeed function's comment...
  
  return (float)energySrc(&engine);
}


double CrProtonSplash::calculate_rate(double old_rate)
{
  return  old_rate;
}


float CrProtonSplash::flux(float latitude, float longitude) const
{
  return  flux();
}


float CrProtonSplash::flux(std::pair<double,double> coords) const
{
  return  flux();
}


std::string CrProtonSplash::title() const
{
  return  "CrProtonSplash";
}


float CrProtonSplash::fraction(float energy)
// This function doesn't work in this class... :-(
{
  return  0;
}


std::pair<float,float> CrProtonSplash::dir(float energy) const
{
  HepJamesRandom  engine;

  std::pair<double,double>  d = dir(energy, &engine);
    
  return  std::pair<float,float>(d.first, d.second);
}

