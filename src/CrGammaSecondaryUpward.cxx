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
 * The absolute flux and spectrum are assumed 
 * to depend on the geomagnetic cutoff energy (and is fixed for Palestine, 
 * Texas in the current codes). 
 * We assume that the spectral shape below 10 MeV and above 1 GeV 
 * does not depend on the zenith angle, 
 * and express the angular dependence of the 
 * relative_flux (<10 MeV) as follows;
 * =======================================
 * relative_flux[c/sr]       theta[rad]
 * ---------------------------------------
 * 1                         0--0.698
 * 0.2067*exp(2.263*theta)   0.698--2.007
 * 2590.67*exp(-2.441*theta) 2.007--2.443
 * 20/3                      2.443-pi 
 * =======================================
 * The energy spectrum of upward 2ndary gamma is expressed with 
 * two power-law functions and 511 keV emission line.
 * A method upwardCRenergy returns an energy and 
 * CrGammaSecondaryUpward::dir 
 * returns a direction in cos(theta) and phi (azimuth angle). 
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
 **************************************************************************
 */

// $Header$

#include <math.h>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrGammaSecondaryUpward.hh"


typedef  double G4double;


// Private function definitions.
namespace {
  const G4double pi = 3.14159265358979323846264339;
  // lower and higher energy limit of secondary gamma-ray in units of GeV
  const G4double lowE_upward  = 30.0e-6; // 30 keV
  const G4double highE_upward = 100.0; // 100 GeV
  // energy of spectral break in units of GeV
  const G4double lowE_break  = 1.0e-3; // 1 MeV
  const G4double highE_break  = 1.0; // 1 GeV
 
  // The constant defined below ("ENERGY_INTEGRAL_upward") is the straight 
  // upward (theta=pi) flux integrated between lowE_* and highE_*.
  // It must be computed in units of [c/s/m^2/sr] and given here.
  // ** Change the value  when the COR or force-field potential
  // are changed.  Current value is for COR = 4.46 [GV] and Phi = 
  // 1100 [MV], respectively, and lowE_upward and highE_upward are 
  // set as 
  // lowE_upward = 30 keV and highE_upward = 100 GeV.
  // **
  // Integrated flux of continuum is 1.36e4 [c/s/m^2/sr].
  // We also added 511keV line of intensity = 470[c/s/m^/sr],
  // hence the total integral becomes 1.41e4 [c/s/m^2/sr]

  //  const G4double ENERGY_INTEGRAL_upward = 1.41e4; //[c/m^2/s/sr]

  // we generate gamma above 1 MeV
  const G4double ENERGY_INTEGRAL_upward = 2378; //[c/m^2/s/sr]

  //============================================================
  /*
   * We assume that the spectral shape of upward component does not 
   * depend on the direction, 
   * and express it with two power-law functions
   * and 511 keV line emission.
   * The vertically upward flux is expressed as
   * 1670*(E/MeV)^-1.34 [c/s/m^2/sr/MeV](30keV-1MeV)
   * 1670*(E/MeV)^-1.70 [c/s/m^2/sr/MeV] (1-1000MeV)
   * 5.29e4 * (E/MeV)^-2.20 [c/s/m^2/sr/MeV] (1-100GeV)
   * 511 keV (470 [c/s/m^2/sr]).
   * Also note that below 10 MeV, downward and upward components
   * have the same spectral shape.
   */
  
  // normalization of incident spectrum
  const G4double A1_upward = 1670;
  const G4double A2_upward = 1670;
  const G4double A3_upward = 5.29e4;
  // differential spectral index
  const G4double a1_upward = 1.34; 
  const G4double a2_upward = 1.70; 
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
    /*****
    G4double envelope_area = envelope1_area + envelope2_area + envelope3_area;
    *****/
    G4double envelope_area = envelope2_area + envelope3_area;

    double Ernd,r; // E means energy in GeV
    G4double E;
    while(1){
      Ernd = engine->flat();
      if (Ernd <= (envelope2_area)/envelope_area){
	// envelop in higher energy
	r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
	E = upwardCRenvelope2_integral_inv(r, cor, solarPotential);
      } else {
	// envelope in highest energy
	r = engine->flat() * (rand_max_3 - rand_min_3) + rand_min_3;
	E = upwardCRenvelope3_integral_inv(r, cor, solarPotential);
      }
      break;
    }
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
std::pair<double,double> CrGammaSecondaryUpward::dir(double energy, 
					       HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The upward has plus sign in cos(theta),
  // and phi = 0 for particle comming along x-axis (from x>0 to x=0)
  // and phi = pi/2 for that comming along y-axis (from y>0 to y=0)
{
  // Based on the measurement by Shonfelder et al.,
  // We expressed zenith-angle dependence of the relative_flux [/sr]
  // in low energy region (<10 MeV) as follows;
  // 1 (0--0.698 [radian], or 0-40[degree])
  // 0.2067*exp(2.263*theta) (theta=0.698--2.007[rad], or 40--115[degree])
  // 2590.67*exp(-2.441*theta) (theta=2.007--2.443[rad], or 115--140[degree])
  // 20/3 (2.443-pi[rad], or 140--180[degree])
  // Here, theta=0 means particle going vertically downward,
  // and theta=pi is the particle going vertically upward.
  // Integrals over solid angle become as follows;
  // 1.47 (theta=0--0.698[rad])
  // 16.08  (theta=0.698--pi/2[rad])
  // 32.453 (theta=pi/2--2.007[rad])
  // 26.373 (theta=2.007--2.443[rad])
  // 9.811  (theta=2.443--pi[rad])

  double rand = engine->flat();
  double theta;
  if (rand*(32.453+26.373+9.811)<=9.811){ // from 140 to 180 deg.
    theta = acos(-1+(engine->flat())*(cos(2.443)-cos(pi)));
  }
  // pi/2 to 2.443 [rad], where the flux [/sr] depends on theta as
  // a*exp(b*theta)
  else if (rand*(32.453+26.373+9.811)<=(26.373+9.811)){ 
    double a=2590.67;
    double b=-2.441; 
    while(1){
      double max = a/b*exp(b*2.443);
      double min = a/b*exp(b*2.007);
      double r = engine->flat() * (max-min) + min;
      theta = 1/b*log(b*r/a);
      if (engine->flat()<sin (theta)){break;}
    }
  } else {
    double a=0.2067;
    double b=2.263; 
    while(1){
      double max = a/b*exp(b*2.007);
      double min = a/b*exp(b*pi/2);
      double r = engine->flat() * (max-min) + min;
      theta = 1/b*log(b*r/a);
      if (engine->flat()<sin (theta)){break;}
    }
  }

  //  double theta = acos(engine->flat());
  double phi   = engine->flat() * 2 * pi;

  return std::pair<double,double>(cos(theta), phi);
}


// Gives back particle energy
double CrGammaSecondaryUpward::energySrc(HepRandomEngine* engine) const
{
  return upwardCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


// flux() returns the value integrated over whole energy and direction
// and devided by 4pi sr: then the unit is [c/s/m^2/sr].
// This value is used as relative normalization among 
// "primary", "secondary downward" and "secondary upward".
double CrGammaSecondaryUpward::flux() const
{
  // Averaged energy-integrated flux over 4 pi solod angle used 
  // in relative normalizaiotn among "primary", "downward(2ndary)" 
  // and "upward(2ndary)".
  // "ENERGY_INTEGRAL_upward" is the energy integrated flux 
  // (between lowE_upward and upward) at theta=pi 
  // (vertically upward).

  // Integral over solid angle from theta=pi/2 to 2.007[rad] 
  // becomes 4.868*ENERGY_INTEGRAL_upward,
  // that from 2.007 to 2.443 becomes
  // 3.955*ENERGY_INTEGRAL_upward,
  // and that from 2.443 to pi becomes
  // 1.471*ENERGY_INTEGRAL_upward (see comments at a "dir" method).
  // Hence the total integrated flux becomes
  // 1.64*2pi*ENERGY_INTEGRAL_upward.

  // Integrated over the lower (earth-side) hemisphere and divided by 4pi.
  return  1.64 * 0.5 * ENERGY_INTEGRAL_upward;  // [c/s/m^2/sr]
}

// Gives back solid angle from which particle comes
double CrGammaSecondaryUpward::solidAngle() const
{
  return  2 * pi;
}


// Gives back particle name
const char* CrGammaSecondaryUpward::particleName() const
{
  return "gamma";
}


//
// "flux" package stuff
//

float CrGammaSecondaryUpward::operator()(float r)
{
  HepJamesRandom  engine;
  engine.setSeed(r * 900000000);
  // 900000000 comes from HepJamesRandom::setSeed function's comment...

  return (float)energySrc(&engine);
}


double CrGammaSecondaryUpward::calculate_rate(double old_rate)
{
  return  old_rate;
}


float CrGammaSecondaryUpward::flux(float latitude, float longitude) const
{
  return  flux();
}


float CrGammaSecondaryUpward::flux(std::pair<double,double> coords) const
{
  return  flux();
}


std::string CrGammaSecondaryUpward::title() const
{
  return  "CrGammaSecondaryUpward";
}


float CrGammaSecondaryUpward::fraction(float energy)
// This function doesn't work in this class... :-(
{
  return  0;
}


std::pair<float,float> CrGammaSecondaryUpward::dir(float energy) const
{
  HepJamesRandom  engine;

  std::pair<double,double>  d = dir(energy, &engine);
  
  return  std::pair<float,float>(d.first, d.second);
}



