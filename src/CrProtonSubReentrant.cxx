/**************************************************************************
 * CrProtonSubReentrant.cc
 **************************************************************************
 * This file defines the subclasses of CrProtonReentrant. 
 * They describe the cosmic-ray proton
 * reentrant spectrum sorted on geomagnetic latitude.
 **************************************************************************
 * 2003-02 written by T. Mizuno.
 **************************************************************************
 */

//$Header$

#include <cmath>

// CLHEP
//#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrProtonSubReentrant.hh"

typedef double G4double;

// private function definitions.
namespace {
  // rest energy (rest mass) of proton in units of GeV
  const G4double restE = 0.938;
  // lower and higher energy limit of secondary proton in units of GeV
  const G4double lowE_reent  = 0.01;
  const G4double highE_reent = 20.0;

  //------------------------------------------------------------
  // cutoff power-law function: A*E^-a*exp(-E/cut)

  // The spectral model: cutoff power-law
  inline G4double cutOffPowSpec
  (G4double norm, G4double index, G4double cutOff, G4double E /* GeV */){
    return norm * pow(E, -index) * exp(-E/cutOff);
  }
  
  // The envelope function of cutoff power-law
  // Random numbers are generated to this envelope function for
  // whose integral the inverse function can be found.  The final one 
  // is obtained by throwing away some random picks of energy.
  inline G4double envelopeCutOffPowSpec
  (G4double norm, G4double index, G4double E /* GeV */){
    return norm * pow(E, -index);
  }
  
  // integral of the envelope function of cutoff power-law
  inline G4double envelopeCutOffPowSpec_integral
  (G4double norm, G4double index, G4double E /* GeV */){
    if (index==1){
      return norm*log(E);
    } else {
      return norm * pow(E, -index+1) / (-index+1);
    }
  }
  
  // inverse function of the integral of the envelope
  inline G4double envelopeCutOffPowSpec_integral_inv
  (G4double norm, G4double index, G4double value){
    if (index==1){
      return exp(value/norm);
    } else {
      return pow( (-index+1)*value/norm, -1./(index-1));
    }
  }
  //------------------------------------------------------------

  //------------------------------------------------------------
  // cutoff power-law function: A*E^-a*exp(-(E/cut)^(-a+1)))

  // The spectral model:
  inline G4double cutOffPowSpec2
  (G4double norm, G4double index, G4double cutOff, G4double E /* GeV */){
    return norm * pow(E, -index) * exp(-pow(E/cutOff, -index+1));
  }
  
  // integral of the envelope function of cutoff power-law
  inline G4double cutOffPowSpec2_integral
  (G4double norm, G4double index, G4double cutOff, G4double E /* GeV */){
    return norm * pow(cutOff, -index+1)/(index-1) *
      exp(-pow(E/cutOff, -index+1));
  }
  
  // inverse function of the integral
  inline G4double cutOffPowSpec2_integral_inv
  (G4double norm, G4double index, G4double cutOff, G4double value){
    return cutOff * pow(-log( (index-1)*value/(norm*pow(cutOff, -index+1)) ), 
			1./(-index+1));
  }
  //------------------------------------------------------------

  //------------------------------------------------------------
  // power-law
  
  // The spectral model: power law
  inline G4double powSpec
  (G4double norm, G4double index, G4double E /* GeV */){
    return norm * pow(E, -index);
  }
  
  // integral of the power-law
  inline G4double powSpec_integral
  (G4double norm, G4double index, G4double E /* GeV */){
    if (index==1){
      return norm * log(E);
    } else {
      return norm * pow(E, -index+1) / (-index+1);
    }
  }
  
  // inverse function of the integral of the power-law
  inline G4double powSpec_integral_inv
  (G4double norm, G4double index, G4double value){
    if (index==1){
      return exp(value/norm);
    } else {
      return pow( (-index+1)*value/norm, -1./(index-1));
    }
  }
  //------------------------------------------------------------
 
}
// end of namespace

//------------------------------------------------------------
// The random number generator for the downward component
// in 0<theta_M<0.2.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrProtonReentrant_0002::CrProtonReentrant_0002(){
  /*
   * Below 100 MeV
   *   j(E) = 0.136*(E/100MeV)^0.4 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 0.123*(E/GeV)^-0.155*exp(-(E/0.51GeV)^0.845) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_reent = 0.136*pow(10.0, 0.4);
  a_reent = -0.4;
  // Normalization and spectral index for E>breakE
  B_reent = 0.123;
  b_reent = 0.155;
  cutOff = 0.509;
  // angular distribution constant
  ang_reent = -0.5;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrProtonReentrant_0002::~CrProtonReentrant_0002()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
G4double CrProtonReentrant_0002::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, breakE);
  G4double rand_min_B = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, breakE);
  G4double rand_max_B = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, highE_reent);
  
  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double spec_area = specA_area + specB_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // spectrum below breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_reent, a_reent, r);
      break;
    } else {
      // spectrum above breakE
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = cutOffPowSpec2_integral_inv(B_reent, b_reent, cutOff, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrProtonReentrant_0002::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, breakE);
  G4double rand_min_2 = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, breakE);
  G4double rand_max_2 = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  // 1+2./3.*ang_reent is to take the angular distribution into account
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2))
    *(1+2./3.*ang_reent);
}

// gives back theta 
G4double CrProtonReentrant_0002::theta(CLHEP::HepRandomEngine* engine){
  G4double theta;
  while(1){
    theta = acos(engine->flat());
    if ((engine->flat())*(1+fabs(ang_reent))
        <=1+ang_reent*sin(theta)*sin(theta)){break;}
  }
  return theta;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the downward component
// in 0.2<theta_M<0.3.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrProtonReentrant_0203::CrProtonReentrant_0203(){
  /*
   * Below 100 MeV
   *   j(E) = 0.1*(E/100MeV)^0.4 [c/s/m^2/sr]
   * 100-600 MeV
   *   j(E) = 0.1*(E/100MeV)^-0.87 [c/s/m^2/sr]
   * Above 600 MeV
   *   j(E) = 0.1*pow(6, -0.87)*(E/600MeV)^-2.53 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<lowE_break
  A_reent = 0.1*pow(10.0, 0.4);
  a_reent = -0.4;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_reent = 0.1*pow(10.0, -0.87);
  b_reent = 0.87;
  // Normalization and spectral index for E>highE_break
  C_reent = 0.1*pow(6.0, -0.87)*pow(1.0/0.6, -2.53);
  c_reent = 2.53;
  // angular distribution constant
  ang_reent = 0.0;
  // The spectrum breaks at 100MeV and 600 MeV
  lowE_break = 0.1;
  highE_break = 0.6;
}
    
CrProtonReentrant_0203::~CrProtonReentrant_0203()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
G4double CrProtonReentrant_0203::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double spec_area = specA_area + specB_area + specC_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<highE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrProtonReentrant_0203::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  // 1+2./3.*ang_reent is to take the angular distribution into account

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3))
      *(1+2./3.*ang_reent);

}

// gives back theta 
G4double CrProtonReentrant_0203::theta(CLHEP::HepRandomEngine* engine){
  G4double theta;
  while(1){
    theta = acos(engine->flat());
    if ((engine->flat())*(1+fabs(ang_reent))
        <=1+ang_reent*sin(theta)*sin(theta)){break;}
  }
  return theta;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the downward component
// in 0.3<theta_M<0.4.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrProtonReentrant_0304::CrProtonReentrant_0304(){
  /*
   * Below 100 MeV
   *   j(E) = 0.1*(E/100MeV)^0.4 [c/s/m^2/sr]
   * 100-600 MeV
   *   j(E) = 0.1*(E/100MeV)^-1.09 [c/s/m^2/sr]
   * Above 600 MeV
   *   j(E) = 0.1*pow(6, -1.09)*(E/600MeV)^-2.40 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<lowE_break
  A_reent = 0.1*pow(10.0, 0.4);
  a_reent = -0.4;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_reent = 0.1*pow(10.0, -1.09);
  b_reent = 1.09;
  // Normalization and spectral index for E>highE_break
  C_reent = 0.1*pow(6.0, -1.09)*pow(1.0/0.6, -2.40);
  c_reent = 2.40;
  // angular distribution constant
  ang_reent = 1.0;
  // The spectrum breaks at 100MeV and 600 MeV
  lowE_break = 0.1;
  highE_break = 0.6;
}
    
CrProtonReentrant_0304::~CrProtonReentrant_0304()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
G4double CrProtonReentrant_0304::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double spec_area = specA_area + specB_area + specC_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<highE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrProtonReentrant_0304::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  // 1+2./3.*ang_reent is to take the angular distribution into account
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3))
      *(1+2./3.*ang_reent);
}

// gives back theta 
G4double CrProtonReentrant_0304::theta(CLHEP::HepRandomEngine* engine){
  G4double theta;
  while(1){
    theta = acos(engine->flat());
    if ((engine->flat())*(1+fabs(ang_reent))
        <=1+ang_reent*sin(theta)*sin(theta)){break;}
  }
  return theta;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the downward component
// in 0.4<theta_M<0.5.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrProtonReentrant_0405::CrProtonReentrant_0405(){
  /*
   * Below 100 MeV
   *   j(E) = 0.1*(E/100MeV)^0.4 [c/s/m^2/sr]
   * 100-600 MeV
   *   j(E) = 0.1*(E/100MeV)^-1.19 [c/s/m^2/sr]
   * Above 600 MeV
   *   j(E) = 0.1*pow(6, -1.19)*(E/600MeV)^-2.54 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<lowE_break
  A_reent = 0.1*pow(10.0, 0.4);
  a_reent = -0.4;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_reent = 0.1*pow(10.0, -1.19);
  b_reent = 1.19;
  // Normalization and spectral index for E>highE_break
  C_reent = 0.1*pow(6.0, -1.19)*pow(1.0/0.6, -2.54);
  c_reent = 2.54;
  // angular distribution constant
  ang_reent = 2.0;
  // The spectrum breaks at 100MeV and 600 MeV
  lowE_break = 0.1;
  highE_break = 0.6;
}
    
CrProtonReentrant_0405::~CrProtonReentrant_0405()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
G4double CrProtonReentrant_0405::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double spec_area = specA_area + specB_area + specC_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<highE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrProtonReentrant_0405::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  // 1+2./3.*ang_reent is to take the angular distribution into account
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3))
      *(1+2./3.*ang_reent);
}

// gives back theta 
G4double CrProtonReentrant_0405::theta(CLHEP::HepRandomEngine* engine){
  G4double theta;
  while(1){
    theta = acos(engine->flat());
    if ((engine->flat())*(1+fabs(ang_reent))
        <=1+ang_reent*sin(theta)*sin(theta)){break;}
  }
  return theta;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the downward component
// in 0.5<theta_M<0.6.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrProtonReentrant_0506::CrProtonReentrant_0506(){
  /*
   * Below 100 MeV
   *   j(E) = 0.1*(E/100MeV)^0.0 [c/s/m^2/sr]
   * 100-400 MeV
   *   j(E) = 0.1*(E/100MeV)^-1.18 [c/s/m^2/sr]
   * Above 400 MeV
   *   j(E) = 0.1*pow(4, -1.18)*(E/400MeV)^-2.31 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<lowE_break
  A_reent = 0.1*pow(10.0, 0.0);
  a_reent = 0.0;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_reent = 0.1*pow(10.0, -1.18);
  b_reent = 1.18;
  // Normalization and spectral index for E>highE_break
  C_reent = 0.1*pow(4.0, -1.18)*pow(1.0/0.4, -2.31);
  c_reent = 2.31;
  // angular distribution constant
  ang_reent = 4.0;
  // The spectrum breaks at 100MeV and 400 MeV
  lowE_break = 0.1;
  highE_break = 0.4;
}
    
CrProtonReentrant_0506::~CrProtonReentrant_0506()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
G4double CrProtonReentrant_0506::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double spec_area = specA_area + specB_area + specC_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<highE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrProtonReentrant_0506::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  // 1+2./3.*ang_reent is to take the angular distribution into account
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3))
      *(1+2./3.*ang_reent);
}

// gives back theta 
G4double CrProtonReentrant_0506::theta(CLHEP::HepRandomEngine* engine){
  G4double theta;
  while(1){
    theta = acos(engine->flat());
    if ((engine->flat())*(1+fabs(ang_reent))
        <=1+ang_reent*sin(theta)*sin(theta)){break;}
  }
  return theta;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the downward component
// in 0.6<theta_M<0.7.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrProtonReentrant_0607::CrProtonReentrant_0607(){
  /*
   * Below 100 MeV
   *   j(E) = 0.13*(E/100MeV)^0.0 [c/s/m^2/sr]
   * 100-300 MeV
   *   j(E) = 0.13*(E/100MeV)^-1.1 [c/s/m^2/sr]
   * Above 300 MeV
   *   j(E) = 0.13*pow(3, -1.1)*(E/300MeV)^-2.25 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_reent = 0.13*pow(10.0, 0.0);
  a_reent = 0.0;
  // Normalization and spectral index for breakE<E
  B_reent = 0.13*pow(10.0, -1.1);
  b_reent = 1.1;
  // Normalization and spectral index for E>highE_break
  C_reent = 0.13*pow(3.0, -1.1)*pow(1.0/0.3, -2.25);
  c_reent = 2.25;
  // angular distribution constant
  ang_reent = 4.0;
  // The spectrum breaks at 100MeV and 300 MeV
  lowE_break = 0.1;
  highE_break = 0.3;
}
    
CrProtonReentrant_0607::~CrProtonReentrant_0607()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
G4double CrProtonReentrant_0607::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double spec_area = specA_area + specB_area + specC_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<highE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrProtonReentrant_0607::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  // 1+2./3.*ang_reent is to take the angular distribution into account
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3))
      *(1+2./3.*ang_reent);
}

// gives back theta 
G4double CrProtonReentrant_0607::theta(CLHEP::HepRandomEngine* engine){
  G4double theta;
  while(1){
    theta = acos(engine->flat());
    if ((engine->flat())*(1+fabs(ang_reent))
        <=1+ang_reent*sin(theta)*sin(theta))       {break;}
  }
  return theta;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the downward component
// in 0.7<theta_M<0.8.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrProtonReentrant_0708::CrProtonReentrant_0708(){
  /*
   * Below 100 MeV
   *   j(E) = 0.2*(E/100MeV)^0.0 [c/s/m^2/sr]
   * 100-400 MeV
   *   j(E) = 0.2*(E/100MeV)^-1.5 [c/s/m^2/sr]
   * Above 400 MeV
   *   j(E) = 0.2*pow(4, -1.5)*(E/400MeV)^-1.85 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_reent = 0.2*pow(10.0, 0.0);
  a_reent = 0.0;
  // Normalization and spectral index for breakE<E
  B_reent = 0.2*pow(10.0, -1.5);
  b_reent = 1.5;
  // Normalization and spectral index for E>highE_break
  C_reent = 0.2*pow(4.0, -1.5)*pow(1.0/0.4, -1.85);
  c_reent = 1.85;
  // angular distribution constant
  ang_reent = 4.0;
  // The spectrum breaks at 100MeV and 400 MeV
  lowE_break = 0.1;
  highE_break = 0.4;
}
    
CrProtonReentrant_0708::~CrProtonReentrant_0708()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
G4double CrProtonReentrant_0708::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double spec_area = specA_area + specB_area + specC_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<highE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrProtonReentrant_0708::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, highE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  // 1+2./3.*ang_reent is to take the angular distribution into account
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3))
      *(1+2./3.*ang_reent);
}

// gives back theta 
G4double CrProtonReentrant_0708::theta(CLHEP::HepRandomEngine* engine){
  G4double theta;
  while(1){
    theta = acos(engine->flat());
    if ((engine->flat())*(1+fabs(ang_reent))
        <=1+ang_reent*sin(theta)*sin(theta)){break;}
  }
  return theta;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the downward component
// in 0.8<theta_M<0.9.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrProtonReentrant_0809::CrProtonReentrant_0809(){
  /*
   * Below 100 MeV
   *   j(E) = 0.23*(E/100MeV)^0.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 0.017*(E/GeV)^-1.83*exp(-(E/0.177GeV)^-0.83) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_reent = 0.23*pow(10.0, 0.0);
  a_reent = 0.0;
  // Normalization and spectral index for E>breakE
  B_reent = 0.017;
  b_reent = 1.83;
  cutOff = 0.177;
  // angular distribution constant
  ang_reent = 4.0;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrProtonReentrant_0809::~CrProtonReentrant_0809()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
G4double CrProtonReentrant_0809::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, breakE);
  G4double rand_min_B = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, breakE);
  G4double rand_max_B = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, highE_reent);
  
  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double spec_area = specA_area + specB_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;

  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // spectrum below breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_reent, a_reent, r);
      break;
    } else {
      // spectrum above breakE
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = cutOffPowSpec2_integral_inv(B_reent, b_reent, cutOff, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrProtonReentrant_0809::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, breakE);
  G4double rand_min_2 = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, breakE);
  G4double rand_max_2 = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  // 1+2./3.*ang_reent is to take the angular distribution into account
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2))
      *(1+2./3.*ang_reent);
}

// gives back theta 
G4double CrProtonReentrant_0809::theta(CLHEP::HepRandomEngine* engine){
  G4double theta;
  while(1){
    theta = acos(engine->flat());
    if ((engine->flat())*(1+fabs(ang_reent))
        <=1+ang_reent*sin(theta)*sin(theta)){break;}
  }
  return theta;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the downward component
// in 0.9<theta_M<1.0.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrProtonReentrant_0910::CrProtonReentrant_0910(){
  /*
   * Below 100 MeV
   *   j(E) = 0.44*(E/100MeV)^0.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 0.037*(E/GeV)^-1.98*exp(-(E/0.21GeV)^-0.98) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_reent = 0.44*pow(10.0, 0.0);
  a_reent = 0.0;
  // Normalization and spectral index for E>breakE
  B_reent = 0.037;
  b_reent = 1.98;
  cutOff = 0.21;
  // angular distribution constant
  ang_reent = 4.0;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrProtonReentrant_0910::~CrProtonReentrant_0910()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
G4double CrProtonReentrant_0910::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, breakE);
  G4double rand_min_B = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, breakE);
  G4double rand_max_B = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, highE_reent);
  
  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double spec_area = specA_area + specB_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;

  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // spectrum below breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_reent, a_reent, r);
      break;
    } else {
      // spectrum above breakE
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = cutOffPowSpec2_integral_inv(B_reent, b_reent, cutOff, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrProtonReentrant_0910::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, breakE);
  G4double rand_min_2 = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, breakE);
  G4double rand_max_2 = 
    cutOffPowSpec2_integral(B_reent, b_reent, cutOff, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  // 1+2./3.*ang_reent is to take the angular distribution into account
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2))
      *(1+2./3.*ang_reent);
}

// gives back theta 
G4double CrProtonReentrant_0910::theta(CLHEP::HepRandomEngine* engine){
  G4double theta;
  while(1){
    theta = acos(engine->flat());
    if ((engine->flat())*(1+fabs(ang_reent))
        <=1+ang_reent*sin(theta)*sin(theta)){break;}
  }
  return theta;
}
//------------------------------------------------------------


