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
 * This program generates the downward component of CR gamma (2ndary) with
 * proper angular distribution and energy spectrum.
 * The absolute flux and spectrum are assumed 
 * to depend on the geomagnetic cutoff energy (and is fixed for Palestine, 
 * Texas in the current codes). 
 * We assume that the spectral shape below 10 MeV and above 1 GeV 
 * does not depend on the
 * zenith angle, and express the angular dependence of the 
 * relative_flux (<10 MeV) as follows;
 * =======================================
 * relative_flux[c/sr]       theta[rad]
 * ---------------------------------------
 * 1                         0--0.698
 * 0.2067*exp(2.263*theta)   0.698--2.007
 * 2590.67*exp(-2.441*theta) 2.007--2.443
 * 20/3                      2.443-pi 
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
 **************************************************************************
 */

// $Header$


#include <math.h>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrGammaSecondaryDownward.hh"


typedef  double G4double;


// Private function definitions.
namespace {
  const G4double pi = 3.14159265358979323846264339;
  // lower and higher energy limit of secondary gamma-ray in units of GeV
  const G4double lowE_downward  = 30.0e-6; // 30 keV
  const G4double highE_downward = 100.0;  // 100 GeV
  // energy of spectral break in units of GeV
  const G4double lowE_break  = 1.0e-3; // 1 MeV
  const G4double highE_break  = 1.0; // 1 GeV
 
  // The constant defined below ("ENERGY_INTEGRAL_downward") is the straight 
  // downward (theta=0) flux integrated between lowE_* and highE_*.
  // It must be computed in units of [c/s/m^2/sr] and given here.
  // ** Change the value  when the COR or force-field potential
  // are changed.  Current value is for COR = 4.46 [GV] and Phi = 
  // 1100 [MV], respectively, and lowE_downward and highE_downward are 
  // set as 
  // lowE_downward = 30 keV and highE_downward = 1 GeV.
  // **
  // Integrated flux of continuum over energy is 2.10e3 [c/s/m^2/sr].
  // We also added 511keV line of intensity = 70.5 [c/s/m^/sr],
  // hence the total integral becomes 2.17e3 [c/s/m^2/sr]

  //const G4double ENERGY_INTEGRAL_downward = 2.17e3; // [c/s/m^2/sr]

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
  const G4double A_511keV = 70.5;  // [c/s/m^2/sr]
  
  const double MeVtoGeV = 1e-3;

  // spectral model of low energy component
  inline G4double downwardCRSpec1(G4double E /* GeV */){
    double EMeV = E*1e3;
    return A1_downward * pow(EMeV, -a1_downward);
  }

  // envelope function in lower energy
  inline G4double downwardCRenvelope1
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A1_downward * pow(EMeV, -a1_downward);
  }

  // integral of envelope function in lower energy
  inline G4double downwardCRenvelope1_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
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
    double EMeV = E*1e3;
    return A2_downward * pow(EMeV, -a2_downward);
  }

  // envelope function in higher energy
  inline G4double downwardCRenvelope2
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A2_downward * pow(EMeV, -a2_downward);
  }

  // integral of envelope function in higher energy
  inline G4double downwardCRenvelope2_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
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
    double EMeV = E*1e3;
    return A3_downward * pow(EMeV, -a3_downward) *
      exp(-pow(EMeV/Cutoff, -a3_downward+1));
  }

  // envelope function in higher energy
  inline G4double downwardCRenvelope3
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A3_downward * pow(EMeV, -a3_downward) *
      exp(-pow(EMeV/Cutoff, -a3_downward+1));
  }

  // integral of envelope function in higher energy
  inline G4double downwardCRenvelope3_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
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
    double EMeV = E*1e3;
    return A4_downward * pow(EMeV, -a4_downward);
  }

  // envelope function in the highest energy
  inline G4double downwardCRenvelope4
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
    return A4_downward * pow(EMeV, -a4_downward);
  }

  // integral of envelope function in the highest energy
  inline G4double downwardCRenvelope4_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    double EMeV = E*1e3;
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
    G4double envelope_area = envelope2_area + envelope3_area +envelope4_area;

    double Ernd,r; 
    G4double E; // E means energy in GeV
      // continuum
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
std::pair<double,double> CrGammaSecondaryDownward::dir(double energy, 
					       HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward has plus sign in cos(theta),
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
  if (rand*(1.47+16.08)<=1.47){ // from 0 to 0.698 radian
    theta = acos(1-(engine->flat())*(cos(0)-cos(0.698)));
  } else { 
    // 0.698 to pi/2 [rad], where the flux [/sr] depends on theta as
    // a*exp(b*theta)
    while(1){
      // zenith angle distribution: flux[c/s/m^2/sr] is proportional to
      // a*exp(b*theta)
      double a=0.2067;
      double b=2.263; 
      double max = a/b*exp(b*pi/2);
      double min = a/b*exp(b*0.698);
      double r = engine->flat() * (max-min) + min;
      theta = 1/b*log(b*r/a);
      if (engine->flat()<sin (theta)){break;}
    }
  }
  
  double phi   = engine->flat() * 2 * pi;

  return std::pair<double,double>(cos(theta), phi);
}


// Gives back particle energy
double CrGammaSecondaryDownward::energySrc(HepRandomEngine* engine) const
{
  return downwardCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


// flux() returns the value integrated over whole energy and direction
// and devided by 4pi sr: then the unit is [c/s/m^2/sr].
// This value is used as relative normalization among 
// "primary", "secondary downward" and "secondary upward".
double CrGammaSecondaryDownward::flux() const
{
  // Averaged energy-integrated flux over 4 pi solod angle used 
  // in relative normalizaiotn among "primary", "downward(2ndary)" 
  // and "upward(2ndary)".
  // "ENERGY_INTEGRAL_downward" is the energy integrated flux 
  // (between lowE_downward and highE_downward) at theta=0 
  // (vertically downward).

  // Integral over solid angle from theta=0 to 0.698[rad] 
  // becomes 1.47*ENERGY_INTEGRAL_downward,
  // and that from 0.698 to pi/2[rad] becomes
  // 16.08*ENERGY_INTEGRAL_downward (see comments at a "dir" method).
  // Hence the total integrated flux becomes
  // 2.79*2pi*ENERGY_INTEGRAL_downward.

  // Integrated over the upper (sky-side) hemisphere and divided by 4pi.
  return  2.79 * 0.5 * ENERGY_INTEGRAL_downward;  // [c/s/m^2/sr]
}

// Gives back solid angle from which particle comes
double CrGammaSecondaryDownward::solidAngle() const
{
  return  2 * pi;
}


// Gives back particle name
const char* CrGammaSecondaryDownward::particleName() const
{
  return "gamma";
}


//
// "flux" package stuff
//

float CrGammaSecondaryDownward::operator()(float r)
{
  HepJamesRandom  engine;
  engine.setSeed(r * 900000000);
  // 900000000 comes from HepJamesRandom::setSeed function's comment...

  return (float)energySrc(&engine);
}


double CrGammaSecondaryDownward::calculate_rate(double old_rate)
{
  return  old_rate;
}


float CrGammaSecondaryDownward::flux(float latitude, float longitude) const
{
  return  flux();
}


float CrGammaSecondaryDownward::flux(std::pair<double,double> coords) const
{
  return  flux();
}


std::string CrGammaSecondaryDownward::title() const
{
  return  "CrGammaSecondaryDownward";
}


float CrGammaSecondaryDownward::fraction(float energy)
// This function doesn't work in this class... :-(
{
  return  0;
}


std::pair<float,float> CrGammaSecondaryDownward::dir(float energy) const
{
  HepJamesRandom  engine;

  std::pair<double,double>  d = dir(energy, &engine);
  
  return  std::pair<float,float>(d.first, d.second);
}

