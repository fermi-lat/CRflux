/**************************************************************************
 * CrGammaSecondaryUpward.cc
 **************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrGamma, the entry-point class for the cosmic-ray gamma
 * generation.
 **************************************************************************
 * This program generates the upward component of cosmic ray gamma (2ndary)
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
 * The energy spectrum of upward 2ndary gamma is expressed with 
 * three power-law functions and 511 keV emission line.
 * A method upwardCRenergy returns an energy and 
 * CrGammaSecondaryUpward::dir 
 * returns a direction in cos(theta) and phi (azimuth angle). 
 * Please note that both the energy spectrum and angular dependence
 * were poorly know and might suffer uncertainty up to factor of 2.
 **************************************************************************
 * Definitions:
 * 1) The z-axis points upward (from Calorimeter to Tracker).
 * 2) Partile of theta=0 points for upward (i.e., comes from zenith)
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
 * 2003-08 Energy spectrum/angular distribution is modified by T. Mizuno.
 * 2003-12 Modified by T. Mizuno
 *           Cutoff rigidity dependence based on Kur'yan et al. (1979)
 *           is implemented.
 **************************************************************************
 */

// $Header$

#include <math.h>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrGammaSecondaryUpward.hh"


typedef double G4double;


// Private function definitions.
namespace {
  // lower and higher energy limit of secondary gamma-ray in units of GeV
  const G4double lowE_upward  = 30.0e-6; // 30 keV
  const G4double highE_upward = 100.0; // 100 GeV
  // energy of spectral break in units of GeV
  const G4double lowE_break  = 10.e-3; // 10 MeV
  const G4double highE_break  = 1.0; // 1 GeV
 
  // The constant defined below ("ENERGY_INTEGRAL_upward") is the straight 
  // upward (theta=pi) flux integrated between lowE_* and highE_*.
  // It must be computed in units of [c/s/m^2/sr] and given here.
  // The value is at Palestine, Texas and the cutoff rigidity dependence
  // based on Kur'yan et al. is taken into account in flux() method.
  // Therefore, if you want to modify the spectral shape, 
  // you also need to modify the value of ENERGY_INTEGRAL_upward.
  // If you want to modify the cutoff rigidity dependence, all you have to do
  // is to modify the return value of flux() method.

  //
  // Integrated flux of continuum is 1.457e4 [c/s/m^2/sr].
  // We also added 511keV line of intensity = 470[c/s/m^/sr],
  // hence the total integral becomes 1.504e4 [c/s/m^2/sr]

  //  const G4double ENERGY_INTEGRAL_upward = 1.504e4; //[c/m^2/s/sr]

  // we generate gamma above 1 MeV
  const G4double ENERGY_INTEGRAL_upward = 3197; //[c/m^2/s/sr]
  const G4double lowE = 1.e-3; // 1 MeV

  //============================================================
  /*
   * We assume that the spectral shape of upward component does not 
   * depend on the direction, 
   * and express it with three power-law functions
   * and 511 keV line emission.
   * The vertically upward flux is expressed as
   * 1500*(E/MeV)^-1.34 [c/s/m^2/sr/MeV](30keV-10MeV)
   * 4854*(E/MeV)^-1.85 [c/s/m^2/sr/MeV] (10-1000MeV)
   * 5.29e4 * (E/MeV)^-2.20 [c/s/m^2/sr/MeV] (1-100GeV)
   * 511 keV (470 [c/s/m^2/sr]).
   */
  
  // normalization of incident spectrum
  const G4double A1_upward = 1500;
  const G4double A2_upward = 4854;
  const G4double A3_upward = 5.29e4;
  // differential spectral index
  const G4double a1_upward = 1.34; 
  const G4double a2_upward = 1.85; 
  const G4double a3_upward = 2.20; 
  // 511keV line intensity
  const G4double A_511keV = 470;  // c/s/m2/str
  

  const double MeVtoGeV = 1e-3;

  // spectral model of low energy component
  inline G4double upwardCRSpec1(G4double E /* GeV */){
    double EMeV = E*1e3;
    return A1_upward * pow(EMeV, -a1_upward);
  }

  // envelope function in lower energy
  inline G4double upwardCRenvelope1
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A1_upward * pow(EMeV, -a1_upward);
  }

  // integral of envelope function in lower energy
  inline G4double upwardCRenvelope1_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A1_upward/(-a1_upward+1) * pow(EMeV, -a1_upward+1);
  }

  // inverse function of integral of envelope function in lower energy 
  // this function returns energy obeying envelope function
  inline G4double upwardCRenvelope1_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a1_upward+1)/ A1_upward * value, 
	       1./(-a1_upward+1)) * MeVtoGeV;
  }

  // spectral model in high energy, power-law component
  inline G4double upwardCRSpec2(G4double E /* GeV */){
    double EMeV = E*1e3;
    return A2_upward * pow(EMeV, -a2_upward);
  }

  // envelope function in higher energy
  inline G4double upwardCRenvelope2
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A2_upward * pow(EMeV, -a2_upward);
  }

  // integral of envelope function in higher energy
  inline G4double upwardCRenvelope2_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A2_upward/(-a2_upward+1) * 
      pow(EMeV, -a2_upward+1);
  }

  // inverse function of integral of envelope function in higher energy 
  // this function returns energy obeying envelope function
  inline G4double upwardCRenvelope2_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a2_upward+1)/ A2_upward * value,
	       1./(-a2_upward+1)) * MeVtoGeV;
  }

  // spectral model in the highest energy, power-law component
  inline G4double upwardCRSpec3(G4double E /* GeV */){
    double EMeV = E*1e3;
    return A3_upward * pow(EMeV, -a3_upward);
  }

  // envelope function in the highest energy
  inline G4double upwardCRenvelope3
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A3_upward * pow(EMeV, -a3_upward);
  }

  // integral of envelope function in the highest energy
  inline G4double upwardCRenvelope3_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A3_upward/(-a3_upward+1) * 
      pow(EMeV, -a3_upward+1);
  }

  // inverse function of integral of envelope function in the highest energy 
  // this function returns energy obeying envelope function
  inline G4double upwardCRenvelope3_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a3_upward+1)/ A3_upward * value,
	       1./(-a3_upward+1)) * MeVtoGeV;
  }


  // The random number generator for the secondary (upward) component
  G4double upwardCRenergy
  (HepRandomEngine* engine, G4double cor, G4double solarPotential){
    G4double rand_min_1 = 
      upwardCRenvelope1_integral(lowE_upward, cor, solarPotential);
    G4double rand_max_1 = 
      upwardCRenvelope1_integral(lowE_break, cor, solarPotential);
    G4double rand_min_2 =
      upwardCRenvelope2_integral(lowE_break, cor, solarPotential);
    G4double rand_max_2 =
      upwardCRenvelope2_integral(highE_break, cor, solarPotential);
    G4double rand_min_3 =
      upwardCRenvelope3_integral(highE_break, cor, solarPotential);
    G4double rand_max_3 =
      upwardCRenvelope3_integral(highE_upward, cor, solarPotential);

    G4double envelope1_area = rand_max_1 - rand_min_1;
    G4double envelope2_area = rand_max_2 - rand_min_2;
    G4double envelope3_area = rand_max_3 - rand_min_3;

    G4double envelope_area = envelope1_area + envelope2_area + envelope3_area;

    double Ernd,r; // E means energy in GeV
    G4double E;
    //----------------------------------------
    // following lines are when you generate gammas from 1 MeV to 100 GeV
    while(1){
      Ernd = engine->flat();
      if (Ernd <= (envelope1_area)/envelope_area){
	// envelop in lower energy
	r = engine->flat() * (rand_max_1 - rand_min_1) + rand_min_1;
	E = upwardCRenvelope1_integral_inv(r, cor, solarPotential);
	if (E>lowE){break;}
      } else if (Ernd <= (envelope1_area+envelope2_area)/envelope_area){
	// envelop in higher energy
	r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
	E = upwardCRenvelope2_integral_inv(r, cor, solarPotential);
	break;
      } else {
	// envelope in highest energy
	r = engine->flat() * (rand_max_3 - rand_min_3) + rand_min_3;
	E = upwardCRenvelope3_integral_inv(r, cor, solarPotential);
	break;
      }
    }
    //----------------------------------------
    return E;
  }
  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrGammaSecondaryUpward::CrGammaSecondaryUpward()
{
  ;
}


CrGammaSecondaryUpward::~CrGammaSecondaryUpward()
{
  ;
}


// Gives back particle direction in (cos(theta), phi)
std::pair<G4double,G4double> CrGammaSecondaryUpward::dir(G4double energy, 
					       HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The upward has plus sign in cos(theta),
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
  if (rand*(27.157+25.414+8.83)<=8.83){ // from 140 to 180 deg.
    theta = acos(-1+(engine->flat())*(cos(2.443)-cos(M_PI)));
  }
  // pi/2 to 2.443 [rad], where the flux [/sr] depends on theta as
  // a*exp(b*theta)
  else if (rand*(27.157+25.414+8.83)<=(25.414+8.82)){ 
    G4double a=4318.9;
    G4double b=-2.693; 
    while(1){
      G4double max = a/b*exp(b*2.443);
      G4double min = a/b*exp(b*2.007);
      G4double r = engine->flat() * (max-min) + min;
      theta = 1/b*log(b*r/a);
      if (engine->flat()<sin (theta)){break;}
    }
  } else {
    G4double a=0.02752;
    G4double b=3.268; 
    while(1){
      G4double max = a/b*exp(b*2.007);
      G4double min = a/b*exp(b*M_PI/2);
      G4double r = engine->flat() * (max-min) + min;
      theta = 1/b*log(b*r/a);
      if (engine->flat()<sin (theta)){break;}
    }
  }

  G4double phi   = engine->flat() * 2 * M_PI;

  return std::pair<G4double,G4double>(cos(theta), phi);
}


// Gives back particle energy
G4double CrGammaSecondaryUpward::energySrc(HepRandomEngine* engine) const
{
  return upwardCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".
G4double CrGammaSecondaryUpward::flux() const
{
  // "ENERGY_INTEGRAL_upward" is the energy integrated flux 
  // (between lowE_upward and upward) at theta=pi 
  // (vertically upward).

  // Integral over solid angle from theta=pi/2 to 2.007[rad] 
  // becomes 4.526*ENERGY_INTEGRAL_upward,
  // that from 2.007 to 2.443 becomes
  // 4.235*ENERGY_INTEGRAL_upward,
  // and that from 2.443 to pi becomes
  // 1.471*ENERGY_INTEGRAL_upward (see comments at a "dir" method).
  // Hence the total integrated flux becomes
  // 1.628*2pi*ENERGY_INTEGRAL_upward.

  // Cutoff rigidity dependence:
  // The atmospheric gamma flux is expected to depend on the cutoff rigidity.
  // The higher the cutoff rigidity is, the lower the primary cosmic-ray flux
  // is, and the lower the atmospheric gamma-ray flux is.
  // The dependence is discussed by several authors, 
  // e.g., Kur'yan et al. (1979) and Dean et al. (1989).
  // Kur'yan et al. (1979) measurement the upward gamma-ray flux
  // above 80 MeV and found the rigidity dependence as R^-1.13.
  // This relation holds for the measurement above 30 MeV by SAS-2 and
  // the balloon observation with the similar apparatus 
  // (Thompson and Simpson 1981).
  // We adopt this relation since most of trigger comes from gamma
  // above 30 MeV.

  G4double Rc_palestine = 4.46; // Cutoff rigidity at Palestine, Texas [GV]
  // Integrated over the lower (earth-side) hemisphere and divided by 2pi.
  return  1.628 * ENERGY_INTEGRAL_upward * 
    pow(m_cutOffRigidity/Rc_palestine, -1.13);  // [c/s/m^2/sr]
}

// Gives back solid angle from which particle comes
G4double CrGammaSecondaryUpward::solidAngle() const
{
  return  2 * M_PI;
}


// Gives back particle name
const char* CrGammaSecondaryUpward::particleName() const
{
  return "gamma";
}


// Gives back the name of the component
std::string CrGammaSecondaryUpward::title() const
{
  return  "CrGammaSecondaryUpward";
}



