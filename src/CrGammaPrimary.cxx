/****************************************************************************
 * CrGammaPrimary.cc:
 ****************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrGamma, the entry-point class for the cosmic-ray gamma
 * generation.
 ****************************************************************************
 * This program provides the primary component of cosmic ray gamma.
 * Its angular distribution is assumed to be uniform for downward 
 * (theta = pi - zenith angle is 0 to pi/2) and zero for theta > pi/2.
 * The following itemizes other important features.
 * 1) One intrinsic spectrum is assumed common to all locations 
 *    on earth: three power-law functions.  
 * 2) Note that secondary component is handled in other programs.
 * 3) A method CrGammaPrimary::energySrc returns an energy and 
 *    CrGammaPrimary::dir returns a direction in 
 *    cos(theta) and phi (azimuth angle). 
 ****************************************************************************
 * Definitions:
 * 1) The z-axis points upward (from Calorimeter to Tracker).  
 * 2) Partile of theta=0 points for downward (i.e., comes from zenith)
 *    and that of theta=pi points for upward (comes from nadir).
 * 3) Partile of phi=0 comes along x-axis (from x>0 to x=0) and
 *    that of phi=pi/2 comes along y-axis (from y>0 to y=0).
 * 4) Particle direction is defined by cos(theta) and phi (in radian).
 ****************************************************************************
 * 2001-07 Written by Y. Fukazawa (Hiroshima Univ)
 * 2001-09 Modified by T. Mizuno (Hiroshima Univ)
 * 2001-09 Modified by Y. Fukazawa (Hiroshima Univ)
 * 2001-12 Modified by T. Mizuno to construct a `stand-alone' module
 ****************************************************************************
 */

// $Header$

#include <math.h>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrGammaPrimary.hh"

typedef double G4double;


// private function definitions.
namespace {
  // lower and higher energy limit of primary gamma-ray in units of GeV
  const G4double lowE_primary  = 30.0e-6; // 30 keV
  const G4double highE_primary = 100.0; // 100 GeV
  // middle-lower and middle-higher energy of spectral break in units of GeV
  const G4double lowE_break  = 50.0e-6; // 50 keV
  const G4double highE_break = 1.0e-3; // 1 MeV
 
  // The constant defined below ("ENERGY_INTEGRAL_primary") is the straight 
  // downward (theta=0) flux integrated between 
  // m_gammaLowEnergy and m_gammaHighEnergy.
  // It will be computed in units of [c/s/m^2/sr] in flux() method
  // from parameters given below 
  // (i.e., A*_primary and a*_primary).
  G4double ENERGY_INTEGRAL_primary;


  //============================================================
  /**
   * We represent the spactral shape of CR gamma-ray (primary)
   * with three power-law functions as
   * 570.8*(E/MeV)^-1.86 (30-50keV)    [c/s/m^2/sr/MeV]
   *  40.0*(E/MeV)^-2.75 (50keV-1MeV)  [c/s/m^2/sr/MeV]
   *  40.0*(E/MeV)^-2.15 (1MeV-100GeV) [c/s/m^2/sr/MeV]
   * For more detail, see reports by T. Mizuno.
   */

  // normalization of incident spectrum
  const G4double A1_primary = 570.8;
  const G4double A2_primary = 40.0;
  const G4double A3_primary = 40.0;
  // differential spectral index 
  const G4double a1_primary = 1.86; 
  const G4double a2_primary = 2.75; 
  const G4double a3_primary = 2.15; 

  const G4double MeVtoGeV = 1e-3;

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

  // The model function of CR primary gamma
  inline G4double primaryCRspec
  (G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */){
    G4double Aflux,EMeV;
    EMeV = E*1e3;
    if ( EMeV < lowE_break ) {
      Aflux = A1_primary * pow(EMeV, -a1_primary);
    } else if ( EMeV < highE_break ) {
      Aflux = A2_primary * pow(EMeV, -a2_primary);
    } else {
      Aflux = A3_primary * pow(EMeV, -a3_primary);
    }
    return Aflux;
  }


  // envelope function in lower energy
  inline G4double primaryCRenvelope1
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A1_primary * pow(EMeV, -a1_primary);
  }

  // integral of envelope function in lower energy
  inline G4double primaryCRenvelope1_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3 ;
    return A1_primary/(-a1_primary+1) * pow(EMeV, -a1_primary+1);
  }

  // inverse function of integral of envelope function in lower energy 
  // this function returns energy obeying envelope function
  inline G4double primaryCRenvelope1_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a1_primary+1)/ A1_primary * value, 
	       1./(-a1_primary+1)) * MeVtoGeV;
   }
 
  // envelope function in middle energy
  inline G4double primaryCRenvelope2
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A2_primary * pow(EMeV, -a2_primary);
  }

  // integral of envelope function in middle energy
  inline G4double primaryCRenvelope2_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A2_primary/(-a2_primary+1) * pow(EMeV, -a2_primary+1);
  }

  // inverse function of integral of envelope function in middle energy 
  // this function returns energy obeying envelope function
  inline G4double primaryCRenvelope2_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a2_primary+1)/ A2_primary * value
	       , 1./(-a2_primary+1)) * MeVtoGeV;
  }

  // envelope function in higher energy
  inline G4double primaryCRenvelope3
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A3_primary * pow(EMeV, -a3_primary);
  }

  // integral of envelope function in higher energy
  inline G4double primaryCRenvelope3_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A3_primary/(-a3_primary+1) * pow(EMeV, -a3_primary+1);
  }

  // inverse function of integral of envelope function in higher energy 
  // this function returns energy obeying envelope function
  inline G4double primaryCRenvelope3_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a3_primary+1)/ A3_primary * value
	       , 1./(-a3_primary+1)) * MeVtoGeV;
  }

  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrGammaPrimary::CrGammaPrimary()
{
  ;
}


CrGammaPrimary::~CrGammaPrimary()
{
  ;
}


// Gives back particle direction in (cos(theta), phi)
std::pair<G4double,G4double> CrGammaPrimary::dir(G4double energy, 
					     HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward direction has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0).
{
  // We assume isotropic distribution from the upper hemisphere.
  // After integration over the azimuth angle (phi), 
  // the theta distribution should be sin(theta) for a constant theta width.

  /// Cos(theta) ranges from 1 to -0.4
  double theta = acos(1.4*engine->flat()-0.4);
  double phi   = engine->flat() * 2 * M_PI;

  return std::pair<G4double,G4double>(cos(theta), phi);
}


// Gives back particle energy
G4double CrGammaPrimary::energySrc(HepRandomEngine* engine) const
{

  G4double rand_min_1 = 
    primaryCRenvelope1_integral(min(lowE_break, max(m_gammaLowEnergy, lowE_primary)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_1 = 
    primaryCRenvelope1_integral(max(lowE_primary, min(m_gammaHighEnergy, lowE_break)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_2 =
    primaryCRenvelope2_integral(min(highE_break, max(m_gammaLowEnergy, lowE_break)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_2 =
    primaryCRenvelope2_integral(max(lowE_break, min(m_gammaHighEnergy, highE_break)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_3 =
    primaryCRenvelope3_integral(min(highE_primary, max(m_gammaLowEnergy, highE_break)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_3 =
    primaryCRenvelope3_integral(max(highE_break, min(m_gammaHighEnergy, highE_primary)), 
				m_cutOffRigidity, m_solarWindPotential);
  
  G4double envelope1_area = rand_max_1 - rand_min_1;
  G4double envelope2_area = rand_max_2 - rand_min_2;
  G4double envelope3_area = rand_max_3 - rand_min_3;
  G4double envelope_area = envelope1_area + envelope2_area + envelope3_area;
  
  G4double Ernd,r;
  G4double E; // E means energy in GeV
  
  while(1){
    Ernd = engine->flat();
    if (Ernd <= envelope1_area/envelope_area){
      // use the envelope function in the lower energy range
      r = engine->flat() * (rand_max_1 - rand_min_1) + rand_min_1;
      E = primaryCRenvelope1_integral_inv(r, m_cutOffRigidity, m_solarWindPotential);
    } else if (Ernd <= (envelope1_area + envelope2_area)/envelope_area){
      // use envelope function in middle energy range
      r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
      E = primaryCRenvelope2_integral_inv(r, m_cutOffRigidity, m_solarWindPotential);
    } else {
      // use envelope function in the higher energy range
      r = engine->flat() * (rand_max_3 - rand_min_3) + rand_min_3;
      E = primaryCRenvelope3_integral_inv(r, m_cutOffRigidity, m_solarWindPotential);
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
G4double CrGammaPrimary::flux() const
{

  G4double rand_min_1 = 
    primaryCRenvelope1_integral(min(lowE_break, max(m_gammaLowEnergy, lowE_primary)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_1 = 
    primaryCRenvelope1_integral(max(lowE_primary, min(m_gammaHighEnergy, lowE_break)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_2 =
    primaryCRenvelope2_integral(min(highE_break, max(m_gammaLowEnergy, lowE_break)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_2 =
    primaryCRenvelope2_integral(max(lowE_break, min(m_gammaHighEnergy, highE_break)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_3 =
    primaryCRenvelope3_integral(min(highE_primary, max(m_gammaLowEnergy, highE_break)), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_3 =
    primaryCRenvelope3_integral(max(highE_break, min(m_gammaHighEnergy, highE_primary)), 
				m_cutOffRigidity, m_solarWindPotential);
  
  G4double envelope1_area = rand_max_1 - rand_min_1;
  G4double envelope2_area = rand_max_2 - rand_min_2;
  G4double envelope3_area = rand_max_3 - rand_min_3;
  G4double envelope_area = envelope1_area + envelope2_area + envelope3_area;
  ENERGY_INTEGRAL_primary = envelope_area;

  /***
  cout << "m_gammaLowEnergy: " << m_gammaLowEnergy << endl;
  cout << "m_gammaHighEnergy: " << m_gammaHighEnergy << endl;
  cout << "envelope1_area: " << envelope1_area << endl;
  cout << "envelope2_area: " << envelope2_area << endl;
  cout << "envelope3_area: " << envelope3_area << endl;
  cout << "envelope_area: " << ENERGY_INTEGRAL_primary << endl;
  ***/

  // We assume that the flux is uniform above the earth horizon.
  // Then the average flux is equal to the vertically downward one.

  return  ENERGY_INTEGRAL_primary;  // [c/s/m^2/sr]
}

// Gives back solid angle from which particle comes
G4double CrGammaPrimary::solidAngle() const
{
  // * 1.4 since Cos(theta) ranges from 1 to -0.4
  return  2 * M_PI * 1.4;
}


// Gives back particle name
const char* CrGammaPrimary::particleName() const
{
  return "gamma";
}


std::string CrGammaPrimary::title() const
{
  return  "CrGammaPrimary";
}


