/**************************************************************************
 * CrPositronSubSplash.cc
 **************************************************************************
 * This file defines the subclasses of CrPositronSplash. 
 * They describe the cosmic-ray positron
 * splash spectrum sorted on geomagnetic latitude.
 **************************************************************************
 * 2003-02 written by T. Mizuno.
 **************************************************************************
 */

//$Header$

#include <cmath>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrPositronSubSplash.hh"

typedef double G4double;

// private function definitions.
namespace {
  // rest energy (rest mass) of positron in units of GeV
  const G4double restE = 5.11e-4;
  // lower and higher energy limit of secondary positron in units of GeV
  const G4double lowE_splash  = 0.01;
  const G4double highE_splash = 20.0;

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
  
  // integral of the cutoff power-law
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
// The random number generator for the upward component
// in 0<theta_M<0.3.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray upward positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrPositronSplash_0003::CrPositronSplash_0003(){
  /*
   * Below 100 MeV
   *   j(E) = 1.0*(E/100MeV)^-1.0 [c/s/m^2/sr]
   * 100 MeV -3 GeV
   *   j(E) = 1.0*(E/100MeV)^-2.2 [c/s/m^2/sr]
   * Above 3 GeV
   *   j(E) = 1.0*pow(30, -2.2)*(E/300MeV)^-4.0 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization and spectral index for E<lowE_break
  A_splash = 1.0*pow(10.0, -1.0);
  a_splash = 1.0;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_splash = 1.0*pow(10.0, -2.2);
  b_splash = 2.2;
  // Normalization and spectral index for E>highE_break
  C_splash = 1.0*pow(30.0, -2.2)*pow(1.0/3.0, -4.0);
  c_splash = 4.0;
  // The spectrum breaks at 100MeV and 3 GeV
  lowE_break = 0.1;
  highE_break = 3.0;
}
    
CrPositronSplash_0003::~CrPositronSplash_0003()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
G4double CrPositronSplash_0003::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_A = 
    powSpec_integral(A_splash, a_splash, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_splash, b_splash, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_splash, b_splash, highE_break);
  G4double rand_min_C = 
    powSpec_integral(C_splash, c_splash, highE_break);
  G4double rand_max_C = 
    powSpec_integral(C_splash, c_splash, highE_splash);

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
    E = powSpec_integral_inv(A_splash, a_splash, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<highE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_splash, b_splash, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_splash, c_splash, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrPositronSplash_0003::upwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_1 = 
    powSpec_integral(A_splash, a_splash, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_splash, b_splash, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_splash, b_splash, highE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_splash, c_splash, highE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_splash, c_splash, highE_splash);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3));
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the upward component
// in 0.3<theta_M<0.6.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray upward positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrPositronSplash_0306::CrPositronSplash_0306(){
  /*
   * Below 100 MeV
   *   j(E) = 0.5*(E/100MeV)^-1.0 [c/s/m^2/sr]
   * Above 100 GeV
   *   j(E) = 0.5*(E/100MeV)^-2.7 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization and spectral index for E<breakE
  A_splash = 0.5*pow(10.0, -1.0);
  a_splash = 1.0;
  // Normalization and spectral index for lowE_break<E
  B_splash = 0.5*pow(10.0, -2.7);
  b_splash = 2.7;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrPositronSplash_0306::~CrPositronSplash_0306()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
G4double CrPositronSplash_0306::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  G4double rand_min_B = 
    powSpec_integral(B_splash, b_splash, breakE);
  G4double rand_max_B = 
    powSpec_integral(B_splash, b_splash, highE_splash);
  
  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double spec_area = specA_area + specB_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_splash, a_splash, r);
      break;
    } else if(rnd <= (specA_area+specB_area) / spec_area){
      // breakE<E, powe-law component
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = powSpec_integral_inv(B_splash, b_splash, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrPositronSplash_0306::upwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_1 = 
    powSpec_integral(A_splash, a_splash, breakE);
  G4double rand_min_2 = 
    powSpec_integral(B_splash, b_splash, breakE);
  G4double rand_max_2 = 
    powSpec_integral(B_splash, b_splash, highE_splash);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2));
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the upward component
// in 0.6<theta_M<0.8.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray upward positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrPositronSplash_0608::CrPositronSplash_0608(){
  /*
   * Below 100 MeV
   *   j(E) = 0.3*(E/100MeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 0.3*(E/100MeV)^-3.3
   *          + 2.0*10^-4*(E/GeV)^1.5*exp(-(E/2.3GeV)^2.5) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.0.
   */
  
  // Normalization and spectral index for E<breakE
  A_splash = 0.3*pow(10.0, -1.0);
  a_splash = 1.0;
  // Normalization and spectral index for breakE<E
  B1_splash = 0.3*pow(10.0, -3.3);
  b1_splash = 3.3;
  B2_splash = 2.0e-4;
  b2_splash = -1.5; // positive slope
  cutOff = 2.3;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrPositronSplash_0608::~CrPositronSplash_0608()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
G4double CrPositronSplash_0608::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  G4double rand_min_B1 = 
    powSpec_integral(B1_splash, b1_splash, breakE);
  G4double rand_max_B1 = 
    powSpec_integral(B1_splash, b1_splash, highE_splash);
  G4double rand_min_B2 = 
    cutOffPowSpec2_integral(B2_splash, b2_splash, cutOff, breakE);
  G4double rand_max_B2 = 
    cutOffPowSpec2_integral(B2_splash, b2_splash, cutOff, highE_splash);
  
  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB1_area = rand_max_B1 - rand_min_B1;
  G4double specB2_area = rand_max_B2 - rand_min_B2;
  G4double spec_area = specA_area + specB1_area + specB2_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_splash, a_splash, r);
      break;
    } else if(rnd <= (specA_area+specB1_area) / spec_area){
      // breakE<E, powe-law component
      r = engine->flat() * (rand_max_B1 - rand_min_B1) + rand_min_B1;
      E = powSpec_integral_inv(B1_splash, b1_splash, r);
      break;
    } else {
      // breakE<E, cut off powe-law component
      r = engine->flat() * (rand_max_B2 - rand_min_B2) + rand_min_B2;
      E = cutOffPowSpec2_integral_inv(B2_splash, b2_splash, cutOff, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrPositronSplash_0608::upwardFlux(){
  G4double rand_min_1 =
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_1 =
    powSpec_integral(A_splash, a_splash, breakE);
  G4double rand_min_2 =
    powSpec_integral(B1_splash, b1_splash, breakE);
  G4double rand_max_2 =
    powSpec_integral(B1_splash, b1_splash, highE_splash);
  G4double rand_min_3 =
    cutOffPowSpec2_integral(B2_splash, b2_splash, cutOff, breakE);
  G4double rand_max_3 =
    cutOffPowSpec2_integral(B2_splash, b2_splash, cutOff, highE_splash);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3));
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the upward component
// in 0.8<theta_M<0.9.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray upward positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrPositronSplash_0809::CrPositronSplash_0809(){
  /*
   * Below 100 MeV
   *   j(E) = 0.3*(E/100MeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 0.3*(E/100MeV)^-3.5
   *          + 1.6*10^-3*(E/GeV)^2.0*exp(-(E/2.0GeV)^3.0) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.0.
   */
  
  // Normalization and spectral index for E<breakE
  A_splash = 0.3*pow(10.0, -1.0);
  a_splash = 1.0;
  // Normalization and spectral index for breakE<E
  B1_splash = 0.3*pow(10.0, -3.5);
  b1_splash = 3.5;
  B2_splash = 1.6e-3;
  b2_splash = -2.0; // positive slope
  cutOff = 1.6;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrPositronSplash_0809::~CrPositronSplash_0809()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
G4double CrPositronSplash_0809::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  G4double rand_min_B1 = 
    powSpec_integral(B1_splash, b1_splash, breakE);
  G4double rand_max_B1 = 
    powSpec_integral(B1_splash, b1_splash, highE_splash);
  G4double rand_min_B2 = 
    cutOffPowSpec2_integral(B2_splash, b2_splash, cutOff, breakE);
  G4double rand_max_B2 = 
    cutOffPowSpec2_integral(B2_splash, b2_splash, cutOff, highE_splash);
  
  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB1_area = rand_max_B1 - rand_min_B1;
  G4double specB2_area = rand_max_B2 - rand_min_B2;
  G4double spec_area = specA_area + specB1_area + specB2_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_splash, a_splash, r);
      break;
    } else if(rnd <= (specA_area+specB1_area) / spec_area){
      // breakE<E, powe-law component
      r = engine->flat() * (rand_max_B1 - rand_min_B1) + rand_min_B1;
      E = powSpec_integral_inv(B1_splash, b1_splash, r);
      break;
    } else {
      // breakE<E, cut off powe-law component
      r = engine->flat() * (rand_max_B2 - rand_min_B2) + rand_min_B2;
      E = cutOffPowSpec2_integral_inv(B2_splash, b2_splash, cutOff, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrPositronSplash_0809::upwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_1 = 
    powSpec_integral(A_splash, a_splash, breakE);
  G4double rand_min_2 = 
    powSpec_integral(B1_splash, b1_splash, breakE);
  G4double rand_max_2 = 
    powSpec_integral(B1_splash, b1_splash, highE_splash);
  G4double rand_min_3 = 
    cutOffPowSpec2_integral(B2_splash, b2_splash, cutOff, breakE);
  G4double rand_max_3 = 
    cutOffPowSpec2_integral(B2_splash, b2_splash, cutOff, highE_splash);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.
  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3));
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the upward component
// in 0.9<theta_M<1.0.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray upward positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrPositronSplash_0910::CrPositronSplash_0910(){
  /*
   * Below 100 MeV
   *   j(E) = 0.3*(E/100MeV)^-1.0 [c/s/m^2/sr]
   * Above 100 GeV
   *   j(E) = 0.3*(E/100MeV)^-2.5 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization and spectral index for E<breakE
  A_splash = 0.3*pow(10.0, -1.0);
  a_splash = 1.0;
  // Normalization and spectral index for lowE_break<E
  B_splash = 0.3*pow(10.0, -2.5);
  b_splash = 2.5;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrPositronSplash_0910::~CrPositronSplash_0910()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
G4double CrPositronSplash_0910::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  G4double rand_min_B = 
    powSpec_integral(B_splash, b_splash, breakE);
  G4double rand_max_B = 
    powSpec_integral(B_splash, b_splash, highE_splash);
  
  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double spec_area = specA_area + specB_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_splash, a_splash, r);
      break;
    } else if(rnd <= (specA_area+specB_area) / spec_area){
      // breakE<E, powe-law component
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = powSpec_integral_inv(B_splash, b_splash, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrPositronSplash_0910::upwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_1 = 
    powSpec_integral(A_splash, a_splash, breakE);
  G4double rand_min_2 = 
    powSpec_integral(B_splash, b_splash, breakE);
  G4double rand_max_2 = 
    powSpec_integral(B_splash, b_splash, highE_splash);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2));
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the upward component
// in 1.0<theta_M<1.1.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray upward positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrPositronSplash_1011::CrPositronSplash_1011(){
  /*
   * Below 100 MeV
   *   j(E) = 6.6*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 9.11*10^-4*(E/GeV)^-2.86  [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization and spectral index for E<breakE
  A_splash = 6.60e-2;
  a_splash = 1.0;
  // Normalization and spectral index for breakE<E
  B_splash = 9.11e-4;
  b_splash = 2.86;
  // The spectrum breaks at 100 MeV
  breakE = 0.1;
}
    
CrPositronSplash_1011::~CrPositronSplash_1011()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
G4double CrPositronSplash_1011::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  G4double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  G4double rand_min_B = 
    powSpec_integral(B_splash, b_splash, breakE);
  G4double rand_max_B = 
    powSpec_integral(B_splash, b_splash, highE_splash);
  
  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double spec_area = specA_area + specB_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_splash, a_splash, r);
      break;
    } else{
      // breakE<E
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = powSpec_integral_inv(B_splash, b_splash, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrPositronSplash_1011::upwardFlux(){
  return 187.45;
}
//------------------------------------------------------------

