/**************************************************************************
 * CrGammaSecondaryDownward.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrGamma, the entry-point class for the cosmic-ray gamma
 * generation.
 **************************************************************************
 * This program generates the downward component of cosmic ray gamma (2ndary)
 * with proper angular distribution and energy spectrum.
 * Although the absolute flux and spectrum could depend
 * on the geomagnetic cutoff energy, 
 * it is is fixed at Palestine, Texas in the current codes. 
 * Angular dependence of the atmospheric gamma-ray flux 
 * has been measured by some authors and known to depend strongly
 * on zenith angle. Here we refer to Schonfelder et al.
 * (1977, ApJ 217, 306) where zenith-angle dependence of 
 * 1.5-10MeV gamma-ray were measured, and use their data
 * to represent angular depedence of 1 MeV gamma-ray for simplicity.
 * Then, the relative flux is expressed as
 * =======================================
 * relative_flux[c/sr]       theta[rad]
 * ---------------------------------------
 * 1/cos(theta)              0--pi/3
 * 0.3673*exp(1.6182*theta)  pi/3--pi/2
 * 0.02752*exp(3.268*theta)  pi/2--2.007
 * 4318.9*exp(-2.693*theta)  2.007--2.443
 * 6                         2.443-pi 
 * =======================================
 * The energy spectrum of downward 2ndary gamma is expressed with 
 * a power-law function, a power-law with exponential cutoff,
 * and 511 keV emission line.
 * A method downwardCRenergy returns an energy and 
 * CrGammaSecondaryDownward::dir 
 * returns a direction in cos(theta) and phi (azimuth angle). 
 **************************************************************************
 * Definitions:
 * 1) The z-axis points upward (from Calorimeter to Tracker).
 * 2) Partile of theta=0 points for downward (i.e., comes from zenith)
 *    and that of theta=pi points for upward (comes from nadir).
 * 3) Partile of phi=0 comes along x-axis (from x>0 to x=0) and
 *    that of phi=pi/2 comes along y-axis (from y>0 to y=0).
 * 4) Particle direction is defined by cos(theta) and phi (in radian).
 **************************************************************************
 * 2001-07 Written by Y. Fukazawa (Hiroshima Univ). 
 * 2001-09 Modified by T. Mizuno  (Hiroshima Univ). 
 * 2001-09 Modified by Y. Fukazawa  (Hiroshima Univ). 
 * 2001-11 angular distribution is changed to be narrower (T. Mizuno)
 * 2001-12 Modified by T. Mizuno to construct a `stand-alone' module
 * 2002-01 Modified by T. Mizuno
 *         angular distribution is changed to be original (broader) one.
 * 2002-05 Angular distribution is modified by T. Mizuno (Hiroshima Univ). 
 * 2003-08 Energy spectrum/angular distribution of gamma upward is modified by T. Mizuno.
 **************************************************************************
 */

// $Header$


#include <math.h>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrGammaSecondaryDownward.hh"

typedef double G4double;

// Private function definitions.
namespace {
  // lower and higher energy limit of secondary gamma-ray in units of GeV
  const G4double lowE_downward  = 30.0e-6; // 30 keV
  const G4double highE_downward = 100.0;  // 100 GeV
  // energy of spectral break in units of GeV
  const G4double lowE_break  = 1.0e-3; // 1 MeV
  const G4double highE_break  = 1.0; // 1 GeV
 
  // The constant defined below ("ENERGY_INTEGRAL_downward") is the straight 
  // downward (theta=0) flux integrated between lowE_* and highE_*.
  // It must be computed in units of [c/s/m^2/sr] and given here.
  //
  // Integrated flux of continuum over energy is 2.10e3 [c/s/m^2/sr].
  // We also added 511keV line of intensity = 78.33 [c/s/m^/sr],
  // hence the total integral becomes 2.178e3 [c/s/m^2/sr]

  //const G4double ENERGY_INTEGRAL_downward = 2.178e3; // [c/s/m^2/sr]

  // we generate gamma above 1 MeV
  const G4double ENERGY_INTEGRAL_downward = 414.2; // [c/s/m^2/sr]


  //============================================================
  /*
   * We assume that the spectral shape of downward component 
   * does not depend on the direction, 
   * and express it with a power-law, a power-law with exponential cutoff,
   * and 511 keV line emission.
   * The vertically downward spectrum is given as
   * 250.0*(E/MeV)^-1.34 [c/s/m^2/sr/MeV](30keV-1MeV)
   * 250.0*(E/MeV)^-1.70 + 1.14e5*(E/MeV)^-2.5*exp(-(E/120MeV)^-1.5)
   *  [c/s/m^2/sr/MeV] (1-1000MeV)
   * 2.15e4*(E/MeV)^-2.20 [c/s/m^2/sr/MeV](1-100GeV)
   * 511 keV (70.5 [c/s/m^2/sr]).
   * Also note that below 10 MeV and above 1 GeV, downward and 
   * upward components have the same spectral shape.
   */
  
  // normalization of incident spectrum
  const G4double A1_downward = 250.0;
  const G4double A2_downward = 250.0;
  const G4double A3_downward = 1.14e5;
  const G4double A4_downward = 2.15e4;
  // differential spectral index
  const G4double a1_downward = 1.34; 
  const G4double a2_downward = 1.70; 
  const G4double a3_downward = 2.50; 
  const G4double a4_downward = 2.20; 
  const G4double Cutoff = 120; // in unit of MeV
  // 511keV line intensity
  const G4double A_511keV = 78.3;  // [c/s/m^2/sr]
  
  const G4double MeVtoGeV = 1e-3;

  // spectral model of low energy component
  inline G4double downwardCRSpec1(G4double E /* GeV */){
    G4double EMeV = E*1e3;
    return A1_downward * pow(EMeV, -a1_downward);
  }

  // envelope function in lower energy
  inline G4double downwardCRenvelope1
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A1_downward * pow(EMeV, -a1_downward);
  }

  // integral of envelope function in lower energy
  inline G4double downwardCRenvelope1_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A1_downward/(-a1_downward+1) * pow(EMeV, -a1_downward+1);
  }

  // inverse function of integral of envelope function in lower energy 
  // this function returns energy obeying envelope function
  inline G4double downwardCRenvelope1_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a1_downward+1)/ A1_downward * value, 
	       1./(-a1_downward+1)) * MeVtoGeV;
  }

  // spectral model in high energy, power-law component
  inline G4double downwardCRSpec2(G4double E /* GeV */){
    G4double EMeV = E*1e3;
    return A2_downward * pow(EMeV, -a2_downward);
  }

  // envelope function in higher energy
  inline G4double downwardCRenvelope2
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A2_downward * pow(EMeV, -a2_downward);
  }

  // integral of envelope function in higher energy
  inline G4double downwardCRenvelope2_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A2_downward/(-a2_downward+1) * 
      pow(EMeV, -a2_downward+1);
  }

  // inverse function of integral of envelope function in higher energy 
  // this function returns energy obeying envelope function
  inline G4double downwardCRenvelope2_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a2_downward+1)/ A2_downward * value,
	       1./(-a2_downward+1)) * MeVtoGeV;
  }

  // spectral model in high energy, power-law with exponential cutoff
  inline G4double downwardCRSpec3(G4double E /* GeV */){
    G4double EMeV = E*1e3;
    return A3_downward * pow(EMeV, -a3_downward) *
      exp(-pow(EMeV/Cutoff, -a3_downward+1));
  }

  // envelope function in higher energy
  inline G4double downwardCRenvelope3
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A3_downward * pow(EMeV, -a3_downward) *
      exp(-pow(EMeV/Cutoff, -a3_downward+1));
  }

  // integral of envelope function in higher energy
  inline G4double downwardCRenvelope3_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A3_downward * pow(Cutoff, -a3_downward+1) / (a3_downward-1) *
      exp(-pow(EMeV/Cutoff, -a3_downward+1));
  }

  // inverse function of integral of envelope function in higher energy 
  // this function returns energy obeying envelope function
  inline G4double downwardCRenvelope3_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return Cutoff * pow( -log((a3_downward-1) * value
			      / (A3_downward * pow(Cutoff, -a3_downward+1))),
			 1./(-a3_downward+1)) * MeVtoGeV;
  }

  // spectral model of the highest energy component
  inline G4double downwardCRSpec4(G4double E /* GeV */){
    G4double EMeV = E*1e3;
    return A4_downward * pow(EMeV, -a4_downward);
  }

  // envelope function in the highest energy
  inline G4double downwardCRenvelope4
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A4_downward * pow(EMeV, -a4_downward);
  }

  // integral of envelope function in the highest energy
  inline G4double downwardCRenvelope4_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A4_downward/(-a4_downward+1) * pow(EMeV, -a4_downward+1);
  }

  // inverse function of integral of envelope function in the highest energy 
  // this function returns energy obeying envelope function
  inline G4double downwardCRenvelope4_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a4_downward+1)/ A4_downward * value, 
	       1./(-a4_downward+1)) * MeVtoGeV;
  }


  // The random number generator for the secondary (downward) component
  G4double downwardCRenergy
  (HepRandomEngine* engine, G4double cor, G4double solarPotential){
    G4double rand_min_1 = 
      downwardCRenvelope1_integral(lowE_downward, cor, solarPotential);
    G4double rand_max_1 = 
      downwardCRenvelope1_integral(lowE_break, cor, solarPotential);
    G4double rand_min_2 =
      downwardCRenvelope2_integral(lowE_break, cor, solarPotential);
    G4double rand_max_2 =
      downwardCRenvelope2_integral(highE_break, cor, solarPotential);
    G4double rand_min_3 =
      downwardCRenvelope3_integral(lowE_break, cor, solarPotential);
    G4double rand_max_3 =
      downwardCRenvelope3_integral(highE_break, cor, solarPotential);
    G4double rand_min_4 =
      downwardCRenvelope4_integral(highE_break, cor, solarPotential);
    G4double rand_max_4 =
      downwardCRenvelope4_integral(highE_downward, cor, solarPotential);

    G4double envelope1_area = rand_max_1 - rand_min_1;
    G4double envelope2_area = rand_max_2 - rand_min_2;
    G4double envelope3_area = rand_max_3 - rand_min_3;
    G4double envelope4_area = rand_max_4 - rand_min_4;

    /*****
    G4double envelope_area = envelope1_area + envelope2_area 
                               + envelope3_area +envelope4_area;
    *****/

    G4double Ernd,r; 
    G4double E; // E means energy in GeV
    //----------------------------------------
    // following lines are when you generate gammas from 1 MeV to 100 GeV
    G4double envelope_area = envelope2_area + envelope3_area +envelope4_area;
    while(1){
      Ernd = engine->flat();
      if (Ernd <= (envelope2_area)/envelope_area){
	// envelope in higher energy
	r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
	E = downwardCRenvelope2_integral_inv(r, cor, solarPotential);
      } else if (Ernd <= (envelope2_area+envelope3_area)/envelope_area){
	// envelope in highest energy
	r = engine->flat() * (rand_max_3 - rand_min_3) + rand_min_3;
	E = downwardCRenvelope3_integral_inv(r, cor, solarPotential);
      } else {
	r = engine->flat() * (rand_max_4 - rand_min_4) + rand_min_4;
	E = downwardCRenvelope4_integral_inv(r, cor, solarPotential);
      }
      break;
    }
    //----------------------------------------
    return E;
  }
  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrGammaSecondaryDownward::CrGammaSecondaryDownward()
{
  ;
}


CrGammaSecondaryDownward::~CrGammaSecondaryDownward()
{
  ;
}


// Gives back particle direction in (cos(theta), phi)
std::pair<G4double,G4double> CrGammaSecondaryDownward::dir(G4double energy, 
					       HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward has plus sign in cos(theta),
  // and phi = 0 for particle comming along x-axis (from x>0 to x=0)
  // and phi = pi/2 for that comming along y-axis (from y>0 to y=0)
{
  // Based on the measurement by Shonfelder et al.,
  // We expressed zenith-angle dependence of the relative_flux [/sr]
  // in low energy region (@1 MeV) as follows;
  // 1/cos(theta) (0--pi/3 [radian], or 0-60[degree])
  // 0.3673*exp(1.6182*theta) (theta=pi/3--pi/2[rad], or 60--90[degree])
  // 0.02752*exp(3.268*theta) (theta=pi/2--2.007[rad], or 90--115[degree])
  // 4318.9*exp(-2.693*theta) (theta=2.007--2.443[rad], or 115--140[degree])
  // 6 (2.443-pi[rad], or 140--180[degree])
  // Here, theta=0 means particle going vertically downward,
  // and theta=pi is the particle going vertically upward.
  // Integrals over solid angle become as follows;
  // 4.355  (theta=0--pi/3[rad])
  // 9,980  (theta=pi/3--pi/2[rad])
  // 25.157 (theta=pi/2--2.007[rad])
  // 25.414 (theta=2.007--2.443[rad])
  // 8.831  (theta=2.443--pi[rad])

  G4double rand = engine->flat();
  G4double theta;
  if (rand*(4.355+9.980)<=4.355){ // from 0 to pi/3 radian
    while(1){
      theta = acos( cos(M_PI/3)+(engine->flat())*(cos(0.)-cos(M_PI/3)) );
      if ( 2*engine->flat()< (1/cos(theta)) ){break;}
    }
  } else { 
    // pi/3 to pi/2 [rad], where the flux [/sr] depends on theta as
    // a*exp(b*theta)
    while(1){
      // zenith angle distribution: flux[c/s/m^2/sr] is proportional to
      // a*exp(b*theta)
      G4double a=0.3673;
      G4double b=1.6182; 
      G4double max = a/b*exp(b*M_PI/2);
      G4double min = a/b*exp(b*M_PI/3);
      G4double r = engine->flat() * (max-min) + min;
      theta = 1/b*log(b*r/a);
      if (engine->flat()<sin (theta)){break;}
    }
  }
  
  G4double phi   = engine->flat() * 2 * M_PI;

  return std::pair<G4double,G4double>(cos(theta), phi);
}


// Gives back particle energy
G4double CrGammaSecondaryDownward::energySrc(HepRandomEngine* engine) const
{
  return downwardCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".
G4double CrGammaSecondaryDownward::flux() const
{
  // "ENERGY_INTEGRAL_downward" is the energy integrated flux 
  // (between lowE_downward and highE_downward) at theta=0 
  // (vertically downward).

  // Integral over solid angle from theta=0 to pi/3[rad] 
  // becomes 4.355*ENERGY_INTEGRAL_downward,
  // and that from pi/3 to pi/2[rad] becomes
  // 9.980*ENERGY_INTEGRAL_downward (see comments at a "dir" method).
  // Hence the total integrated flux becomes
  // 2.281*2pi*ENERGY_INTEGRAL_downward.

  // Integrated over the upper (sky-side) hemisphere and divided by 2pi.
  return  2.281 * ENERGY_INTEGRAL_downward;  // [c/s/m^2/sr]
}

// Gives back solid angle from which particle comes
G4double CrGammaSecondaryDownward::solidAngle() const
{
  return  2 * M_PI;
}


// Gives back particle name
const char* CrGammaSecondaryDownward::particleName() const
{
  return "gamma";
}


// Gives back the name of the component
std::string CrGammaSecondaryDownward::title() const
{
  return  "CrGammaSecondaryDownward";
}

