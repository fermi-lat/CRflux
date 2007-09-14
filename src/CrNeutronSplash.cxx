/****************************************************************************
 * CrNeutronSplash.cc
 ****************************************************************************
 * Atmospheric neutron generator for GLAST LAT.
 ****************************************************************************
 * 2007-09 Written by T. Mizuno (Hiroshima Univ.)
 ****************************************************************************
 */

// $Header$

#include <math.h>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrNeutronSplash.hh"

typedef double G4double;


// private function definitions.
namespace {
  // lower and higher energy limit in units of GeV
  const G4double lowE_neutron  = 0.01; // 10 MeV
  const G4double highE_neutron = 1000.0; // 1 TeV
  // energy of spectral break in units of GeV
  const G4double lowE_break  = 0.1; // 100 MeV
 
  // The constant defined below ("ENERGY_INTEGRAL_neutron") is the straight 
  // downward (theta=0) flux integrated between 
  // lowE_neutron and highE_neutron.
  // It will be computed in units of [c/s/m^2/sr] in flux() method
  // from parameters given below 
  // (i.e., A*_neutron and a*_neutron).
  G4double ENERGY_INTEGRAL_neutron;


  //============================================================
  /**
   * We represent the spactral shape of CR neutrons
   * with two power-law functions as
   *  1.e3/(2pi)*(E/1MeV)^-1.05 (1MeV-100MeV)   [c/s/m^2/sr/MeV]
   *    = 159.1*(E/MeV)^-1.05
   *  8/(2pi)*(E/100MeV)^-1.9 (100MeV-1TeV)     [c/s/m^2/sr/MeV]
   *    = 8033*(E/MeV)^-1.9
   */

  // normalization of incident spectrum
  const G4double A0_neutron = 159.1;
  const G4double A1_neutron = 8033;
  // differential spectral index 
  const G4double a0_neutron = 1.05; 
  const G4double a1_neutron = 1.9; 

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
    if ( E < lowE_break ) {
      Aflux = A0_neutron * pow(EMeV, -a0_neutron);
    } else {
      Aflux = A1_neutron * pow(EMeV, -a1_neutron);
    }
    return Aflux;
  }


  // envelope function in the lowest energy
  inline G4double primaryCRenvelope0
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1.e3;
    return A0_neutron * pow(EMeV, -a0_neutron);
  }

  // integral of envelope function in the lowest
  inline G4double primaryCRenvelope0_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A0_neutron/(-a0_neutron+1) * pow(EMeV, -a0_neutron+1);
  }

  // inverse function of integral of envelope function in the lowest energy 
  // this function returns energy obeying envelope function
  inline G4double primaryCRenvelope0_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a0_neutron+1)/ A0_neutron * value
	       , 1./(-a0_neutron+1)) * MeVtoGeV;
   }
 
  // envelope function in lower energy
  inline G4double primaryCRenvelope1
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A1_neutron * pow(EMeV, -a1_neutron);
  }

  // integral of envelope function in lower energy
  inline G4double primaryCRenvelope1_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double EMeV = E*1e3;
    return A1_neutron/(-a1_neutron+1) * pow(EMeV, -a1_neutron+1);
  }

  // inverse function of integral of envelope function in lower energy 
  // this function returns energy obeying envelope function
  inline G4double primaryCRenvelope1_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a1_neutron+1)/ A1_neutron * value
	       , 1./(-a1_neutron+1)) * MeVtoGeV;
   }
 

  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrNeutronSplash::CrNeutronSplash()
{
  ;
}


CrNeutronSplash::~CrNeutronSplash()
{
  ;
}


// Gives back particle direction in (cos(theta), phi)
std::pair<G4double,G4double> CrNeutronSplash::dir(G4double energy, 
					     CLHEP::HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward direction has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0).
{
  // We assume isotropic distribution from the upper hemisphere.
  // After integration over the azimuth angle (phi), 
  // the theta distribution should be sin(theta) for a constant theta width.


  /// Cos(theta) ranges from -1 to -0.4
  G4double theta = acos(0.6*engine->flat()-1); // for low-earth orbit
//  G4double theta = acos(2*engine->flat()-1);
  G4double phi   = engine->flat() * 2 * M_PI;

  return std::pair<G4double,G4double>(cos(theta), phi);
}


// Gives back particle energy
G4double CrNeutronSplash::energySrc(CLHEP::HepRandomEngine* engine) const
{

  G4double rand_min_0 = 
    primaryCRenvelope0_integral(min(lowE_break, lowE_neutron), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_0 = 
    primaryCRenvelope0_integral(max(lowE_neutron, lowE_break), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_1 = 
    primaryCRenvelope1_integral(min(highE_neutron, lowE_break), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_1 = 
    primaryCRenvelope1_integral(max(lowE_break, highE_neutron), 
				m_cutOffRigidity, m_solarWindPotential);
  
  G4double envelope0_area = rand_max_0 - rand_min_0;
  G4double envelope1_area = rand_max_1 - rand_min_1;
  G4double envelope_area = envelope0_area + envelope1_area;
  
  G4double Ernd,r;
  G4double E; // E means energy in GeV
  
  while(1){
    Ernd = engine->flat();
    if (Ernd <= envelope0_area/envelope_area){
      // use the envelope function in the lowest energy range
      r = engine->flat() * (rand_max_0 - rand_min_0) + rand_min_0;
      E = primaryCRenvelope0_integral_inv(r, m_cutOffRigidity, m_solarWindPotential);
    } else {
      // use the envelope function in the lower energy range
      r = engine->flat() * (rand_max_1 - rand_min_1) + rand_min_1;
      E = primaryCRenvelope1_integral_inv(r, m_cutOffRigidity, m_solarWindPotential);
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
G4double CrNeutronSplash::flux() const
{

  G4double rand_min_0 = 
    primaryCRenvelope0_integral(min(lowE_break, lowE_neutron), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_0 = 
    primaryCRenvelope0_integral(max(lowE_neutron, lowE_break), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_min_1 = 
    primaryCRenvelope1_integral(min(highE_neutron, lowE_break), 
				m_cutOffRigidity, m_solarWindPotential);
  G4double rand_max_1 = 
    primaryCRenvelope1_integral(max(lowE_break, highE_neutron), 
				m_cutOffRigidity, m_solarWindPotential);
  
  G4double envelope0_area = rand_max_0 - rand_min_0;
  G4double envelope1_area = rand_max_1 - rand_min_1;
  G4double envelope_area = envelope0_area + envelope1_area;
  
  ENERGY_INTEGRAL_neutron = envelope_area;


  /***
  cout << "m_gammaLowEnergy: " << m_gammaLowEnergy << endl;
  cout << "m_gammaHighEnergy: " << m_gammaHighEnergy << endl;
  cout << "envelope2_area: " << envelope2_area << endl;
  cout << "envelope3_area: " << envelope3_area << endl;
  cout << "envelope_area: " << ENERGY_INTEGRAL_neutron << endl;
  ***/

  // We assume that the flux is uniform in region of cos(theta)>0.4,
  // then the average flux is equal to the vertically downward one.
  // We also assume that the flux scales as exp(-0.152Rc),
  // where Rc is the geomagnetic cutoff rigidity.

  return  exp(-0.152*(m_cutOffRigidity-5))*ENERGY_INTEGRAL_neutron;  // [c/s/m^2/sr]

}

// Gives back solid angle from which particle comes
G4double CrNeutronSplash::solidAngle() const
{
  // Cos(theta) ranges from 1 to 0.4
  return  1.2 * M_PI;
}


// Gives back particle name
const char* CrNeutronSplash::particleName() const
{
  return "neutron";
}


std::string CrNeutronSplash::title() const
{
  return  "CrNeutronSplash";
}


