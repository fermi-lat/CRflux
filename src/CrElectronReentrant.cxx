/**************************************************************************
 * CrElectronReentrant.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrElectron, the entry-point class for the cosmic-ray electron 
 * generation.
 **************************************************************************
 * This program generates the reentrant cosmic ray electron flux with
 * proper angular distribution and energy spectrum.
 * The absolute flux and spectrum of "reentrant" electrons are assumed 
 * to depend on the geomagnetic cutoff energy (and is fixed for Palestine,
 * Texas in the current codes). The flux is assumed to 
 * depend on zenith angle as 
 *   1 + 0.6*sin(theta) for theta = 0 to pi/2, and zero for theta> pi/2
 * (theta = pi - zenith angle).
 * The energy spectrum above 100 MeV is assumed to be of a combination 
 * of a power-law function and a power-law with exponential cutoff 
 * in lower energy referring to AMS data.
 * Below 100 MeV, the spectrum is extrapolated down to 10 MeV
 * assuming that the flux is proportional to E^-1.
 * A method reentrantCRenergy returns an energy and 
 * CrElectronReentrant::dir returns a direction 
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
 * 2001-11 Modified by T. Mizuno.
 *           angular distribution is changed to be uniform
 *           energy spectrum is extrapolated down to 10 MeV
 * 2001-12 Modified by T. Mizuno to construct a `stand-alone' module
 * 2003-02 Modified by T. Mizuno to generate flux at any position in orbit.
 **************************************************************************
 */

//$Header$


#include <math.h>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrElectronReentrant.hh"
#include "CrElectronSubReentrant.hh"


typedef  double G4double;

// private function definitions.
namespace {
  const G4double pi    = 3.14159265358979323846264339;
  const G4double restE = 5.11e-4; // rest energy of electron in [GeV]

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

CrElectronReentrant::CrElectronReentrant()
{
  crElectronReentrant_0003 = new CrElectronReentrant_0003();
  crElectronReentrant_0306 = new CrElectronReentrant_0306();
  crElectronReentrant_0608 = new CrElectronReentrant_0608();
  crElectronReentrant_0809 = new CrElectronReentrant_0809();
  crElectronReentrant_0910 = new CrElectronReentrant_0910();
  crElectronReentrant_1011 = new CrElectronReentrant_1011();
}


CrElectronReentrant::~CrElectronReentrant()
{
  delete crElectronReentrant_0003;
  delete crElectronReentrant_0306;
  delete crElectronReentrant_0608;
  delete crElectronReentrant_0809;
  delete crElectronReentrant_0910;
  delete crElectronReentrant_1011;
}


// Gives back particle direction in (cos(theta), phi)
std::pair<double,double> CrElectronReentrant::dir
(double energy, HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0).
{
  /***
   * Here we assume that
   *   theta (pi - zenith angle): 
   *     the flux (per steradian) is proportional to 1 + 0.6*sin(theta)
   *   azimuth angle (phi) : isotropic
   *
   * Reference:
   *   Tylka, A. J. 2000-05-12, GLAST team internal report 
   *   "A Review of Cosmic-Ray Albedo Studies: 1949-1970" (1 + 0.6 sin(theta))
   */

  double theta;
  while (1){
    theta = acos(engine->flat());
    if (engine->flat()*1.6<1+0.6*sin(theta)){break;}
  }

  double phi = engine->flat() * 2 * pi;

  return std::pair<double,double>(cos(theta), phi);
}


// Gives back particle energy
double CrElectronReentrant::energySrc(HepRandomEngine* engine) const
{
  double r1, r2;
  if (m_geomagneticLatitude<0.15){
    return crElectronReentrant_0003->energy(engine);
  } else if (m_geomagneticLatitude>=0.15 && m_geomagneticLatitude<0.45){
    r1 = m_geomagneticLatitude-0.15;
    r2 = 0.45-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crElectronReentrant_0003->energy(engine);
    } else {
      return crElectronReentrant_0306->energy(engine);
    }
  } else if (m_geomagneticLatitude>=0.45 && m_geomagneticLatitude<0.7){
    r1 = m_geomagneticLatitude-0.45;
    r2 = 0.7-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crElectronReentrant_0306->energy(engine);
    } else {
      return crElectronReentrant_0608->energy(engine);
    }
  } else if (m_geomagneticLatitude>=0.7 && m_geomagneticLatitude<0.85){
    r1 = m_geomagneticLatitude-0.7;
    r2 = 0.85-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crElectronReentrant_0608->energy(engine);
    } else {
      return crElectronReentrant_0809->energy(engine);
    }
  } else if (m_geomagneticLatitude>=0.85 && m_geomagneticLatitude<0.95){
    r1 = m_geomagneticLatitude-0.85;
    r2 = 0.95-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crElectronReentrant_0809->energy(engine);
    } else {
      return crElectronReentrant_0910->energy(engine);
    }
  } else if (m_geomagneticLatitude>=0.95 && m_geomagneticLatitude<1.05){
    r1 = m_geomagneticLatitude-0.95;
    r2 = 1.05-m_geomagneticLatitude;
    if (engine->flat()*(r1+r2)<r2){
      return crElectronReentrant_0910->energy(engine);
    } else {
      return crElectronReentrant_1011->energy(engine);
    }
  } else{
      return crElectronReentrant_1011->energy(engine);
  }

}


// flux() returns the value integrated over whole energy and direction
// and devided by 4pi sr: then the unit is [c/s/m^2/sr].
// This value is used as relative normalization among 
// "primary", "reentrant" and "splash".
double CrElectronReentrant::flux() const
{
  // Averaged energy-integrated flux over 4 pi solod angle used 
  // in relative normalizaiotn among "primary", "reentrant" and "splash".

  // energy integrated vertically downward flux, [c/s/m^2/sr]
  double downwardFlux; 
  double r1, r2;
  if (m_geomagneticLatitude<0.15){
    downwardFlux = crElectronReentrant_0003->downwardFlux();
  } else if (m_geomagneticLatitude>=0.15 && m_geomagneticLatitude<0.45){
    r1 = m_geomagneticLatitude-0.15;
    r2 = 0.45-m_geomagneticLatitude;
    downwardFlux = ( r2*crElectronReentrant_0003->downwardFlux()
             +r1*crElectronReentrant_0306->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude>=0.45 && m_geomagneticLatitude<0.7){
    r1 = m_geomagneticLatitude-0.45;
    r2 = 0.7-m_geomagneticLatitude;
    downwardFlux = ( r2*crElectronReentrant_0306->downwardFlux()
             +r1*crElectronReentrant_0608->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude>=0.7 && m_geomagneticLatitude<0.85){
    r1 = m_geomagneticLatitude-0.7;
    r2 = 0.85-m_geomagneticLatitude;
    downwardFlux = ( r2*crElectronReentrant_0608->downwardFlux()
             +r1*crElectronReentrant_0809->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude>=0.85 && m_geomagneticLatitude<0.95){
    r1 = m_geomagneticLatitude-0.85;
    r2 = 0.95-m_geomagneticLatitude;
    downwardFlux = ( r2*crElectronReentrant_0809->downwardFlux()
             +r1*crElectronReentrant_0910->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude>=0.95 && m_geomagneticLatitude<1.05){
    r1 = m_geomagneticLatitude-0.95;
    r2 = 1.05-m_geomagneticLatitude;
    downwardFlux = ( r2*crElectronReentrant_0910->downwardFlux()
             +r1*crElectronReentrant_1011->downwardFlux() )/(r1+r2);
  } else{
    downwardFlux = crElectronReentrant_1011->downwardFlux();
  }

  // We have assumed that the flux is proportional to 1+0.6 sin(theta)
  // Then, the flux integrated over solid angle is
  // (1+0.15pi)*downwardFlux*2pi
  return 0.5 * (1 + 0.15*pi) * downwardFlux; // [c/s/m^2/sr]


}


// Gives back solid angle from which particle comes
double CrElectronReentrant::solidAngle() const
{
  return 2 * pi;
}


// Gives back particle name
const char* CrElectronReentrant::particleName() const
{
  return "e-";
}

//
// "flux" package stuff
//

float CrElectronReentrant::operator()(float r)
{
  HepJamesRandom  engine;
  engine.setSeed(r * 900000000);
  // 900000000 comes from HepJamesRandom::setSeed function's comment...

  return (float)energySrc(&engine);
}


double CrElectronReentrant::calculate_rate(double old_rate)
{
  return  old_rate;
}


float CrElectronReentrant::flux(float latitude, float longitude) const
{
  return  flux();
}


float CrElectronReentrant::flux(std::pair<double,double> coords) const
{
  return  flux();
}


std::string CrElectronReentrant::title() const
{
  return  "CrElectronReentrant";
}


float CrElectronReentrant::fraction(float energy)
  // This function doesn't work in this class... :-(
{
  return  0;
}


std::pair<float,float> CrElectronReentrant::dir(float energy) const
{
  HepJamesRandom  engine;

  std::pair<double,double>  d = dir(energy, &engine);

  return  std::pair<float,float>(d.first, d.second);
}