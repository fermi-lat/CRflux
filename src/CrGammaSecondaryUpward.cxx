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
 * Angular dependence of the atmospheric gamma-ray flux 
 * has been measured by some authors and known to depend strongly
 * on zenith angle. Here we refer to Thompson et al.
 * (1981, JGR 86, 1265) where zenith-angle dependence of were
 * calculated and compared with data by SAS-2.
 * The relative flux is expressed as
 * =======================================
 * relative_flux[c/sr]       cosTheta[rad]
 * ---------------------------------------
 * 1                         -1.0 -- -0.8
 * 39.33*exp(4.59*cosTheta)  -0.8 -- -0.45
 * 343.5*exp(9.40*cosTheta)  -0.45 -- -0.4
 * 8*exp(-40(cosTheta+0.4))  -0.4 -- 0.0
 * =======================================
 * Please note that the relative flux is 1 in cosTheta=-1.0 to -0.8,
 * 5.0 at cosTheta = -0.45, 8.0 at cosTheta = -0.4,
 * and drops off sharpely.
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
 * 2003-12 Modified by T. Mizuno
 *         user can set lower and Upper energy to generate gammas.
 * 2004-05 Modified by T. Mizuno
 *           Spectrum and angular distribution are modified.
 * 2004-10 Modified by T. Mizuno
 *           Angular dependence is modeled based on SAS-2 paper
 *           (Thompsont et al. 1981), to represent the flux
 *           in satellite altitude.
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
  const G4double lowE_break  = 20.e-3; // 10 MeV
  const G4double highE_break  = 1.0; // 1 GeV
 
  // The constant defined below ("ENERGY_INTEGRAL_upward") is the straight 
  // upward (theta=pi) flux integrated between 
  // m_gammaLowEnergy and m_gammaHighEnergy.
  // It will be computed in units of [c/s/m^2/sr] in flux() method
  // from parameters given below
  // (i.e., A*_upward, a*_upward, and A_511keV).
  // The value is at Palestine, Texas and the cutoff rigidity dependence
  // based on Kur'yan et al. is taken into account in flux() method.
  G4double ENERGY_INTEGRAL_upward;

  //============================================================
  /*
   * We assume that the spectral shape of upward component does not 
   * depend on the direction, 
   * and express it with three power-law functions
   * and 511 keV line emission.
   * The vertically upward flux is expressed as
   * 1010*(E/MeV)^-1.34 [c/s/m^2/sr/MeV](30keV-20MeV)
   * 7290*(E/MeV)^-2.0 [c/s/m^2/sr/MeV] (20-1000MeV)
   * 2.9e4 * (E/MeV)^-2.20 [c/s/m^2/sr/MeV] (1-100GeV)
   * 511 keV (470 [c/s/m^2/sr]).
   */
  
  // normalization of incident spectrum
  const G4double A1_upward = 1010;
  const G4double A2_upward = 7290;
  const G4double A3_upward = 2.9e4;
  // differential spectral index
  const G4double a1_upward = 1.34; 
  const G4double a2_upward = 2.00; 
  const G4double a3_upward = 2.20; 
  // 511keV line intensity based on measurement by Imhof et al.
  const G4double A_511keV = 470;  // [c/s/m^2/sr]
  

  const double MeVtoGeV = 1e-3;

  // function that returns the higher value
  inline G4double max(G4double x, G4double y){
    if (x>y){
      return x;
    } else {
      return y;
    }
  }

  // function that returns the lower value
  inline G4double min(G4double x, G4double y){
    if (x<y){
      return x;
    } else {
      return y;
    }
  }

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
  // We refer to the model for SAS-2 data in Thompson et al. 1981
  // (solid line in Fig. 3 of Thompson et al. 1981, JGR 86, 1265)
  // and modled the zenith-angle dependence of the atmospheric gamma
  // with the following functions.
  // 1 (cosTheta = -1 to -0.8)
  // 39.33*exp(4.59*cosTheta) (cosTheta = -0.8 to -0.45)
  // 343.5*exp(9.40*cosTheta) (cosTheta = -0.45 to -0.4)
  // 8*exp(-40(cosTheta+0.4)) (cosTheta = -0.4 to 0)
  // Here, costheta=-1 means particle going vertically upward.
  // Integrals over solid angle become as follows;
  // 0.2*2pi   (cosTheta = -1 to -0.8)
  // 0.868*2pi (cosTheta = -0.8 to -0.45)
  // 0.319*2pi (cosTheta = -0.45 to -0.4)
  // 0.20*2pi  (cosTheta = -0.4 to 0)

  G4double rand = engine->flat();
  G4double theta, cosTheta;
  if (rand*(0.2+0.868+0.319+0.2)<=0.2){ // cosTheta = -1 to -0.8
    theta = acos(-1+(engine->flat())*(-0.8+1));
  }
  // cosTheta = -0.8 to 0, where the flux [/sr] depends on theta as
  // a*exp(b*cosTheta)
  else if (rand*(0.2+0.868+0.319+0.2)<=(0.2+0.868)){ 
    G4double a=39.33;
    G4double b=4.59; 
    G4double min = a/b*exp(b*(-0.8));
    G4double max = a/b*exp(b*(-0.45));
    G4double r = engine->flat() * (max-min) + min;
    cosTheta = 1/b*log(b*r/a);
    theta = acos(cosTheta);
  } else if (rand*(0.2+0.868+0.319+0.2)<=(0.2+0.868+0.319)){
    G4double a=343.5;
    G4double b=9.40; 
    G4double min = a/b*exp(b*(-0.45));
    G4double max = a/b*exp(b*(-0.4));
    G4double r = engine->flat() * (max-min) + min;
    cosTheta = 1/b*log(b*r/a);
    theta = acos(cosTheta);
  } else {
    G4double a=8;
    G4double b=-40; 
    G4double min = a/b*exp(b*(-0.4));
    G4double max = a/b*exp(b*(0.0));
    G4double r = engine->flat() * (max-min) + min;
    cosTheta = 1/b*log(b*r/a);
    theta = acos(cosTheta);
  }

  G4double phi   = engine->flat() * 2 * M_PI;

  return std::pair<G4double,G4double>(cos(theta), phi);
}


// Gives back particle energy
G4double CrGammaSecondaryUpward::energySrc(HepRandomEngine* engine) const
{

  G4double rand_min_1 = 
    upwardCRenvelope1_integral(min(lowE_break, max(m_gammaLowEnergy, lowE_upward)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_1 = 
    upwardCRenvelope1_integral(max(lowE_upward, min(m_gammaHighEnergy, lowE_break)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_2 =
    upwardCRenvelope2_integral(min(highE_break, max(m_gammaLowEnergy,lowE_break)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_2 =
    upwardCRenvelope2_integral(max(lowE_break, min(m_gammaHighEnergy, highE_break)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_3 =
    upwardCRenvelope3_integral(min(highE_upward, max(m_gammaLowEnergy, highE_break)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_3 =
    upwardCRenvelope3_integral(max(highE_break, min(m_gammaHighEnergy, highE_upward)), 
			       m_cutOffRigidity, m_solarWindPotential);
  
  G4double envelope1_area = rand_max_1 - rand_min_1;
  G4double envelope2_area = rand_max_2 - rand_min_2;
  G4double envelope3_area = rand_max_3 - rand_min_3;
  G4double envelope_511keV;
  if (m_gammaHighEnergy>=511.0e-6 && m_gammaLowEnergy<=511.0e-6){
    envelope_511keV = A_511keV;
  } else {
    envelope_511keV = 0.0;
  }
  
  G4double envelope_area = envelope1_area + envelope2_area 
    + envelope3_area + envelope_511keV;
  
  G4double Ernd,r;
  G4double E; // E means energy in GeV

  while(1){
    Ernd = engine->flat();
    if (Ernd <= (envelope1_area)/envelope_area){
      // envelop in lower energy
      r = engine->flat() * (rand_max_1 - rand_min_1) + rand_min_1;
      E = upwardCRenvelope1_integral_inv(r, m_cutOffRigidity, m_solarWindPotential);
    } else if (Ernd <= (envelope1_area+envelope2_area)/envelope_area){
      // envelop in higher energy
      r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
      E = upwardCRenvelope2_integral_inv(r, m_cutOffRigidity, m_solarWindPotential);
    } else if (Ernd <= (envelope1_area+envelope2_area+envelope3_area)/envelope_area){
      // envelope in highest energy
      r = engine->flat() * (rand_max_3 - rand_min_3) + rand_min_3;
      E = upwardCRenvelope3_integral_inv(r, m_cutOffRigidity, m_solarWindPotential);
    } else {
      E = 511.0e-6;
    }
    break;
  }
  
  return E;
}



// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".
G4double CrGammaSecondaryUpward::flux() const
{
  G4double rand_min_1 = 
    upwardCRenvelope1_integral(min(lowE_break, max(m_gammaLowEnergy, lowE_upward)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_1 = 
    upwardCRenvelope1_integral(max(lowE_upward, min(m_gammaHighEnergy, lowE_break)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_2 =
    upwardCRenvelope2_integral(min(highE_break, max(m_gammaLowEnergy,lowE_break)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_2 =
    upwardCRenvelope2_integral(max(lowE_break, min(m_gammaHighEnergy, highE_break)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_3 =
    upwardCRenvelope3_integral(min(highE_upward, max(m_gammaLowEnergy, highE_break)), 
			       m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_3 =
    upwardCRenvelope3_integral(max(highE_break, min(m_gammaHighEnergy, highE_upward)),
                               m_cutOffRigidity, m_solarWindPotential);
  
  G4double envelope1_area = rand_max_1 - rand_min_1;
  G4double envelope2_area = rand_max_2 - rand_min_2;
  G4double envelope3_area = rand_max_3 - rand_min_3;
  G4double envelope_511keV;
  if (m_gammaHighEnergy>=511.0e-6 && m_gammaLowEnergy<=511.0e-6){
    envelope_511keV = A_511keV;
  } else {
    envelope_511keV = 0.0;
  }
  
  G4double envelope_area = envelope1_area + envelope2_area 
    + envelope3_area + envelope_511keV;
  ENERGY_INTEGRAL_upward = envelope_area;

  /***
  cout << "m_gammaLowEnergy: " << m_gammaLowEnergy << endl;
  cout << "m_gammaHighEnergy: " << m_gammaHighEnergy << endl;
  cout << "envelope1_area: " << envelope1_area << endl;
  cout << "envelope2_area: " << envelope2_area << endl;
  cout << "envelope3_area: " << envelope3_area << endl;
  cout << "envelope_511keV: " << envelope_511keV << endl;
  cout << "envelope_area: " << ENERGY_INTEGRAL_upward << endl;
  ***/

  // "ENERGY_INTEGRAL_upward" is the energy integrated flux 
  // (between gammaLowEnergy and gammaHighEnergy) at theta=pi 
  // (vertically upward).

  // Integral over solid angle from theta=pi/2 to pi is
  // 1.587*2pi*ENERGY_INTEGRAL_upward.
  // (please look at the formula given in "dir" method)

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
  return  1.587 * ENERGY_INTEGRAL_upward * 
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



