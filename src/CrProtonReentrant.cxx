/**************************************************************************
 * CrProtonReentrant.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrProton, the entry-point class for the cosmic-ray proton 
 * generation.
 **************************************************************************
 * This program generates the reentrant cosmic ray proton flux with
 * proper angular distribution and energy spectrum.
 * The absolute flux and spectrum of "reentrant" protons are assumed 
 * to depend on the geomagnetic cutoff energy (and is fixed for Palestine, 
 * Texas in the current codes). The flux is assumed to 
 * depend on zenith angle as 
 *   1 + 0.6*sin(theta) for theta = 0 to pi/2, and zero for theta > pi/2
 * (theta = pi - zenith angle).
 * The energy spectrum above 100 MeV is assumed to be of a power-low 
 * with exponential cut-off at a low energy referring to AMS data.
 * Below 100 MeV, the spectrum is extrapolated down to 10 MeV
 * assuming that the flux is proportional to E^-1.
  * A method reentrantCRenergy returns an energy and CrProtonReentrant::dir 
 * returns a direction in cos(theta) and phi (azimuth angle). 
 * Energy and angular distribution of reentrant proton is assumed to be
 * the same as those of splash proton.
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
 * in 2000 Originally written by by M. Ozaki (ISAS) 
 *         and K. Hirano (Hiroshima Univ.)
 * 2001-04 Spectrum part is modified by T. Mizuno (Hiroshima Univ.)
 * 2001-05 Spectrum part is modified by T. Mizuno (Hiroshima Univ.)
 * 2001-05 Comments added and unused codes removed by T. Kamae (SLAC)
 * 2001-05 Program checked by T. Kamae and H. Mizushima (Hiroshima)
 * 2001-11 Modified by T. Mizuno.
 *           angular distribution is changed to be uniform
 *           energy spectrum is extrapolated down to 10 MeV
 * 2001-11 Modified by T. Mizuno to construct a `stand-alone' module
 * 2003-02 Modified by T. Mizuno to generate flux at any position in orbit.
 **************************************************************************
 */

//$Header$

#include <math.h>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrProtonReentrant.hh"
#include "CrProtonSubReentrant.hh"


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
  // as a function of kinetic energy [GeV].
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

CrProtonReentrant::CrProtonReentrant()
{
  crProtonReentrant_0002 = new CrProtonReentrant_0002();
  crProtonReentrant_0203 = new CrProtonReentrant_0203();
  crProtonReentrant_0304 = new CrProtonReentrant_0304();
  crProtonReentrant_0405 = new CrProtonReentrant_0405();
  crProtonReentrant_0506 = new CrProtonReentrant_0506();
  crProtonReentrant_0607 = new CrProtonReentrant_0607();
  crProtonReentrant_0708 = new CrProtonReentrant_0708();
  crProtonReentrant_0809 = new CrProtonReentrant_0809();
  crProtonReentrant_0910 = new CrProtonReentrant_0910();
}


CrProtonReentrant::~CrProtonReentrant()
{
  delete crProtonReentrant_0002;
  delete crProtonReentrant_0203;
  delete crProtonReentrant_0304;
  delete crProtonReentrant_0405;
  delete crProtonReentrant_0506;
  delete crProtonReentrant_0607;
  delete crProtonReentrant_0708;
  delete crProtonReentrant_0809;
  delete crProtonReentrant_0910;
}


// Gives back particle direction in (cos(theta), phi)
std::pair<G4double,G4double> CrProtonReentrant::dir(G4double energy, 
						HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward direction has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0)
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

  G4double theta;
  /// Cos(theta) ranges from 1 to -0.4
  while (1){
    theta = acos(1.4*engine->flat()-0.4);
	double temp = engine->flat()*1.6;
    if (temp<1+0.6*sin(theta) && 1+0.6*sin(theta) > 0){break;}
	if (temp<-(1+0.6*sin(theta)) && -(1+0.6*sin(theta)) > 0){break;}

  }

  G4double phi = engine->flat() * 2 * M_PI;

  return std::pair<G4double,G4double>(cos(theta), phi);
}


// Gives back particle energy
G4double CrProtonReentrant::energySrc(HepRandomEngine* engine) const
{
  //  return crProtonReentrant_0002->energy(engine);
  //  return crProtonReentrant_0203->energy(engine);
  //  return crProtonReentrant_0304->energy(engine);
  //  return crProtonReentrant_0405->energy(engine);
  //  return crProtonReentrant_0506->energy(engine);
  //  return crProtonReentrant_0607->energy(engine);
  //  return crProtonReentrant_0708->energy(engine);
  //  return crProtonReentrant_0809->energy(engine);
  //  return crProtonReentrant_0910->energy(engine);

  G4double r1, r2;
  if (m_geomagneticLatitude*M_PI/180.0<0.15){
    return crProtonReentrant_0002->energy(engine);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.15 && m_geomagneticLatitude*M_PI/180.0<0.25){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.15;
    r2 = 0.25-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonReentrant_0002->energy(engine);
    } else {
      return crProtonReentrant_0203->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.25 && m_geomagneticLatitude*M_PI/180.0<0.35){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.25;
    r2 = 0.35-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonReentrant_0203->energy(engine);
    } else {
      return crProtonReentrant_0304->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.35 && m_geomagneticLatitude*M_PI/180.0<0.45){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.35;
    r2 = 0.45-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonReentrant_0304->energy(engine);
    } else {
      return crProtonReentrant_0405->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.45 && m_geomagneticLatitude*M_PI/180.0<0.55){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.45;
    r2 = 0.55-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonReentrant_0405->energy(engine);
    } else {
      return crProtonReentrant_0506->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.55 && m_geomagneticLatitude*M_PI/180.0<0.65){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.55;
    r2 = 0.65-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonReentrant_0506->energy(engine);
    } else {
      return crProtonReentrant_0607->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.65 && m_geomagneticLatitude*M_PI/180.0<0.75){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.65;
    r2 = 0.75-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonReentrant_0607->energy(engine);
    } else {
      return crProtonReentrant_0708->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.75 && m_geomagneticLatitude*M_PI/180.0<0.85){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.75;
    r2 = 0.85-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonReentrant_0708->energy(engine);
    } else {
      return crProtonReentrant_0809->energy(engine);
    }
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.85 && m_geomagneticLatitude*M_PI/180.0<0.95){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.85;
    r2 = 0.95-m_geomagneticLatitude*M_PI/180.0;
    if (engine->flat()*(r1+r2)<r2){
      return crProtonReentrant_0809->energy(engine);
    } else {
      return crProtonReentrant_0910->energy(engine);
    }
  } else {
    return crProtonReentrant_0910->energy(engine);
  }

}


// flux() returns the value integrated over whole energy and direction
// and devided by 4pi sr: then the unit is [c/s/m^2/sr].
// This value is used as relative normalization among 
// "primary", "reentrant" and "splash".
G4double CrProtonReentrant::flux() const
{
  // Averaged energy-integrated flux over 4 pi solod angle used 
  // in relative normalizaiotn among "primary", "reentrant" and "splash".

  // energy integrated vertically downward flux, [c/s/m^2/sr]
  G4double downwardFlux; 
  G4double r1, r2;
  if (m_geomagneticLatitude*M_PI/180.0<0.15){
    downwardFlux = crProtonReentrant_0002->downwardFlux();
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.15 && m_geomagneticLatitude*M_PI/180.0<0.25){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.15;
    r2 = 0.25-m_geomagneticLatitude*M_PI/180.0;
    downwardFlux = ( r2*crProtonReentrant_0002->downwardFlux()
	     +r1*crProtonReentrant_0203->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.25 && m_geomagneticLatitude*M_PI/180.0<0.35){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.25;
    r2 = 0.35-m_geomagneticLatitude*M_PI/180.0;
    downwardFlux = ( r2*crProtonReentrant_0203->downwardFlux()
	     +r1*crProtonReentrant_0304->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.35 && m_geomagneticLatitude*M_PI/180.0<0.45){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.35;
    r2 = 0.45-m_geomagneticLatitude*M_PI/180.0;
    downwardFlux = ( r2*crProtonReentrant_0304->downwardFlux()
	     +r1*crProtonReentrant_0405->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.45 && m_geomagneticLatitude*M_PI/180.0<0.55){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.45;
    r2 = 0.55-m_geomagneticLatitude*M_PI/180.0;
    downwardFlux = ( r2*crProtonReentrant_0405->downwardFlux()
	     +r1*crProtonReentrant_0506->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.55 && m_geomagneticLatitude*M_PI/180.0<0.65){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.55;
    r2 = 0.65-m_geomagneticLatitude*M_PI/180.0;
    downwardFlux = ( r2*crProtonReentrant_0506->downwardFlux()
	     +r1*crProtonReentrant_0607->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.65 && m_geomagneticLatitude*M_PI/180.0<0.75){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.65;
    r2 = 0.75-m_geomagneticLatitude*M_PI/180.0;
    downwardFlux = ( r2*crProtonReentrant_0607->downwardFlux()
	     +r1*crProtonReentrant_0708->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.75 && m_geomagneticLatitude*M_PI/180.0<0.85){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.75;
    r2 = 0.85-m_geomagneticLatitude*M_PI/180.0;
    downwardFlux = ( r2*crProtonReentrant_0708->downwardFlux()
	     +r1*crProtonReentrant_0809->downwardFlux() )/(r1+r2);
  } else if (m_geomagneticLatitude*M_PI/180.0>=0.85 && m_geomagneticLatitude*M_PI/180.0<0.95){
    r1 = m_geomagneticLatitude*M_PI/180.0-0.85;
    r2 = 0.95-m_geomagneticLatitude*M_PI/180.0;
    downwardFlux = ( r2*crProtonReentrant_0809->downwardFlux()
	     +r1*crProtonReentrant_0910->downwardFlux() )/(r1+r2);
  } else {
    downwardFlux = crProtonReentrant_0910->downwardFlux();
  }

  // We have assumed that the flux is proportional to 1+0.6 sin(theta)
  // Then, the flux integrated over solid angle is
  // (1+0.15pi)*downwardFlux*2pi
  return 0.5 * (1 + 0.15*M_PI) * downwardFlux; // [c/s/m^2/sr]
}


// Gives back solid angle from which particle comes
G4double CrProtonReentrant::solidAngle() const
{
   // * 1.4 since Cos(theta) ranges from 1 to -0.4
  return 2 * M_PI * 1.4;
}


// Gives back particle name
const char* CrProtonReentrant::particleName() const
{
  return "proton";
}


// Gives back the name of the component
std::string CrProtonReentrant::title() const
{
  return  "CrProtonReentrant";
}

