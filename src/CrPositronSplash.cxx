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
 * This program generates the splash cosmic ray positron flux with
 * proper angular distribution and energy spectrum.
 * The absolute flux and spectrum of "splash" positron are assumed 
 * to depend on the geomagnetic cutoff energy (and is fixed for Palestine,
 * Texas in the current codes). The flux is assumed to 
 * depend on zenith angle as 
 *   1 + 0.6*sin(theta) for theta = pi/2 to pi, and zero for theta < pi/2
 * (theta = pi - zenith angle).
 * The energy spectrum above 100 MeV is assumed to be of a power-low 
 * referring to AMS data. Below this energy, it is extrapolated down
 * to 10 MeV assuming that the flux is proportional to E^-1.
 * A method splashCRenergy returns an energy and 
 * CrPositronSplash::dir returns a direction 
 * in cos(theta) and phi (azimuth angle). 
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
 * 2001-07 Written by Y. Fukazawa (Hiroshima Univ.)
 * 2001-10 Modified by T. Mizuno and Y. Fukazawa
 * 2001-11 angular distribution is changed to be uniform (T. Mizuno)
 * 2001-12 Modified by T. Mizuno to cunstruct a `stand-alone' module
 * 2003-02 Modified by T. Mizuno to generate flux at any position in orbit.
 **************************************************************************
 */

// $Header$

#include <math.h>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrPositronSplash.hh"
#include "CrPositronSubSplash.hh"


typedef  double G4double;


// private function definitions.
namespace {
  const G4double pi    = 3.14159265358979323846264339;
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
std::pair<double,double> CrPositronSplash::dir(double energy, 
					       HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0).
{
  /***
   * Here we assume that
   *   theta(pi - zenith angle): 
   *     the flux (per steradian) is proportional to 1 + 0.6*sin(theta)
   *  azimuth angle (phi) : isotropic
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

  return  std::pair<double,double>(cos(theta), phi);
}


// Gives back particle energy
double CrPositronSplash::energySrc(HepRandomEngine* engine) const
{
  double r1, r2;
  if (m_geomagneticLatitude<0.15){
    return crPositronSplash_0003->energy(engine);
  } else if (m_geomagneticLatitude>=0.15 && m_geomagneticLatitude<0.45){
    r1 = m_geomagneticLatitude-0.15;
    r2 = 0.45-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0003->energy(engine);
    } else {
      return crPositronSplash_0306->energy(engine);
    }
  } else if (m_geomagneticLatitude>=0.45 && m_geomagneticLatitude<0.7){
    r1 = m_geomagneticLatitude-0.45;
    r2 = 0.7-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0306->energy(engine);
    } else {
      return crPositronSplash_0608->energy(engine);
    }
  } else if (m_geomagneticLatitude>=0.7 && m_geomagneticLatitude<0.85){
    r1 = m_geomagneticLatitude-0.7;
    r2 = 0.85-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0608->energy(engine);
    } else {
      return crPositronSplash_0809->energy(engine);
    }
  } else if (m_geomagneticLatitude>=0.85 && m_geomagneticLatitude<0.95){
    r1 = m_geomagneticLatitude-0.85;
    r2 = 0.95-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0809->energy(engine);
    } else {
      return crPositronSplash_0910->energy(engine);
    }
  } else if (m_geomagneticLatitude>=0.95 && m_geomagneticLatitude<1.05){
    r1 = m_geomagneticLatitude-0.95;
    r2 = 1.05-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crPositronSplash_0910->energy(engine);
    } else {
      return crPositronSplash_1011->energy(engine);
    }
  } else{
    return crPositronSplash_1011->energy(engine);
  }

}


// flux() returns the value integrated over whole energy and direction
// and devided by 4pi sr: then the unit is [c/s/m^2/sr].
// This value is used as relative normalization among 
// "primary", "reentrant" and "splash".
double CrPositronSplash::flux() const
{
  // Averaged energy-integrated flux over 4 pi solod angle used 
  // in relative normalizaiotn among "primary", "reentrant" and "splash".

  // energy integrated vertically downward flux, [c/s/m^2/sr]
  double upwardFlux; 
  double r1, r2;

  if (m_geomagneticLatitude<0.15){
    upwardFlux = crPositronSplash_0003->upwardFlux();
  } else if (m_geomagneticLatitude>=0.15 && m_geomagneticLatitude<0.45){
    r1 = m_geomagneticLatitude-0.15;
    r2 = 0.45-m_geomagneticLatitude;
    upwardFlux = ( r2*crPositronSplash_0003->upwardFlux()
                   +r1*crPositronSplash_0306->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude>=0.45 && m_geomagneticLatitude<0.7){
    r1 = m_geomagneticLatitude-0.45;
    r2 = 0.7-m_geomagneticLatitude;
    upwardFlux = ( r2*crPositronSplash_0306->upwardFlux()
                     +r1*crPositronSplash_0608->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude>=0.7 && m_geomagneticLatitude<0.85){
    r1 = m_geomagneticLatitude-0.7;
    r2 = 0.85-m_geomagneticLatitude;
    upwardFlux = ( r2*crPositronSplash_0608->upwardFlux()
                   +r1*crPositronSplash_0809->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude>=0.85 && m_geomagneticLatitude<0.95){
    r1 = m_geomagneticLatitude-0.85;
    r2 = 0.95-m_geomagneticLatitude;
    upwardFlux = ( r2*crPositronSplash_0809->upwardFlux()
                     +r1*crPositronSplash_0910->upwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude>=0.95 && m_geomagneticLatitude<1.05){
    r1 = m_geomagneticLatitude-0.95;
    r2 = 1.05-m_geomagneticLatitude;
    upwardFlux = ( r2*crPositronSplash_0910->upwardFlux()
                     +r1*crPositronSplash_1011->upwardFlux() )/(r1+r2);
  } else{
    upwardFlux = crPositronSplash_1011->upwardFlux();
  }

  // We have assumed that the flux is proportional to 1+0.6 sin(theta)
  // Then, the flux integrated over solid angle is
  // (1+0.15pi)*upwardFlux*2pi
  return 0.5 * (1 + 0.15*pi) * upwardFlux; // [c/s/m^2/sr]

}


// Gives back solid angle from which particle comes
double CrPositronSplash::solidAngle() const
{
  return 2 * pi;
}


// Gives back particle name
const char* CrPositronSplash::particleName() const
{
  return "e+";
}


//
// "flux" package stuff
//

float CrPositronSplash::operator()(float r)
{
  HepJamesRandom  engine;
  engine.setSeed(r * 900000000);
  // 900000000 comes from HepJamesRandom::setSeed function's comment...

  return (float)energySrc(&engine);
}


double CrPositronSplash::calculate_rate(double old_rate)
{
  return  old_rate;
}


float CrPositronSplash::flux(float latitude, float longitude) const
{
  return  flux();
}


float CrPositronSplash::flux(std::pair<double,double> coords) const
{
  return  flux();
}


std::string CrPositronSplash::title() const
{
  return  "CrPositronSplash";
}


float CrPositronSplash::fraction(float energy)
// This function doesn't work in this class... :-(
{
  return  0;
}


std::pair<float,float> CrPositronSplash::dir(float energy) const
{
  HepJamesRandom  engine;

  std::pair<double,double>  d = dir(energy, &engine);

  return  std::pair<float,float>(d.first, d.second);
}
