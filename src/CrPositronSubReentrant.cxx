/**************************************************************************
 * CrPositronSubReentrant.cc
 **************************************************************************
 * This file defines the subclasses of CrPositronReentrant. 
 * They describe the cosmic-ray positron
 * reentrant spectrum sorted on geomagnetic latitude.
 **************************************************************************
 * 2003-02 written by T. Mizuno.
 **************************************************************************
 */

//$Header$

#include <math.h>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrPositronSubReentrant.hh"

typedef  double G4double;

// private function definitions.
namespace {
  const G4double pi = 3.14159265358979323846264339;
  // rest energy (rest mass) of positron in units of GeV
  const G4double restE = 5.11e-4;
  // lower and higher energy limit of secondary positron in units of GeV
  const G4double lowE_reent  = 0.01;
  const G4double highE_reent = 20.0;

  //------------------------------------------------------------
  // cutoff power-law function: A*E^-a*exp(-E/cut)

  // The spectral model: cutoff power-law
  inline double cutOffPowSpec
  (double norm, double index, double cutOff, double E /* GeV */){
    return norm * pow(E, -index) * exp(-E/cutOff);
  }
  
  // The envelope function of cutoff power-law
  // Random numbers are generated to this envelope function for
  // whose integral the inverse function can be found.  The final one 
  // is obtained by throwing away some random picks of energy.
  inline double envelopeCutOffPowSpec
  (double norm, double index, double E /* GeV */){
    return norm * pow(E, -index);
  }
  
  // integral of the envelope function of cutoff power-law
  inline double envelopeCutOffPowSpec_integral
  (double norm, double index, double E /* GeV */){
    if (index==1){
      return norm*log(E);
    } else {
      return norm * pow(E, -index+1) / (-index+1);
    }
  }
  
  // inverse function of the integral of the envelope
  inline double envelopeCutOffPowSpec_integral_inv
  (double norm, double index, double value){
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
  inline double cutOffPowSpec2
  (double norm, double index, double cutOff, double E /* GeV */){
    return norm * pow(E, -index) * exp(-pow(E/cutOff, -index+1));
  }
  
  // integral of the cutoff power-law
  inline double cutOffPowSpec2_integral
  (double norm, double index, double cutOff, double E /* GeV */){
    return norm * pow(cutOff, -index+1)/(index-1) *
      exp(-pow(E/cutOff, -index+1));
  }
  
  // inverse function of the integral
  inline double cutOffPowSpec2_integral_inv
  (double norm, double index, double cutOff, double value){
    return cutOff * pow(-log( (index-1)*value/(norm*pow(cutOff, -index+1)) ), 
			1./(-index+1));
  }
  //------------------------------------------------------------

  //------------------------------------------------------------
  // power-law
  
  // The spectral model: power law
  inline double powSpec
  (double norm, double index, double E /* GeV */){
    return norm * pow(E, -index);
  }
  
  // integral of the power-law
  inline double powSpec_integral
  (double norm, double index, double E /* GeV */){
    if (index==1){
      return norm * log(E);
    } else {
      return norm * pow(E, -index+1) / (-index+1);
    }
  }
  
  // inverse function of the integral of the power-law
  inline double powSpec_integral_inv
  (double norm, double index, double value){
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
// The random number generator for the reentrant component
// in 0<theta_M<0.3.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray reentrant positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrPositronReentrant_0003::CrPositronReentrant_0003(){
  /*
   * Below 100 MeV
   *   j(E) = 9.8*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 2 GeV
   *   j(E) = 1.02*10^-2*(E/GeV)^-2.0*exp(-E/2.5GeV) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization and spectral index for E<breakE
  A_reent = 9.8e-2;
  a_reent = 1.0;
  // Normalization, spectral index and cutoff for breakE<E
  B_reent = 1.02e-2;
  b_reent = 2.0;
  cutOff = 2.5;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrPositronReentrant_0003::~CrPositronReentrant_0003()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
double CrPositronReentrant_0003::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  double rand_max_A = 
    powSpec_integral(A_reent, a_reent, breakE);
  double rand_min_B = 
    envelopeCutOffPowSpec_integral(B_reent, b_reent, breakE);
  double rand_max_B = 
    envelopeCutOffPowSpec_integral(B_reent, b_reent, highE_reent);
  
  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double spec_area = specA_area + specB_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_reent, a_reent, r);
      break;
    } else {
      // breakE<E
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = envelopeCutOffPowSpec_integral_inv(B_reent, b_reent, r);
      if (engine->flat() * envelopeCutOffPowSpec(B_reent, b_reent, E)
	  < cutOffPowSpec(B_reent, b_reent, cutOff, E)){
        break;
      }
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrPositronReentrant_0003::downwardFlux(){
  return 311.0;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the reentrant component
// in 0.3<theta_M<0.6.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray reentrant positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrPositronReentrant_0306::CrPositronReentrant_0306(){
  /*
   * Below 100 MeV
   *   j(E) = 5.03*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 1.18*10^-3*(E/GeV)^-2.63 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization and spectral index for E<breakE
  A_reent = 5.03e-2;
  a_reent = 1.0;
  // Normalization and spectral index for breakE<E
  B_reent = 1.18e-3;
  b_reent = 2.63;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrPositronReentrant_0306::~CrPositronReentrant_0306()
{
  ;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
double CrPositronReentrant_0306::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  double rand_max_A = 
    powSpec_integral(A_reent, a_reent, breakE);
  double rand_min_B = 
    powSpec_integral(B_reent, b_reent, breakE);
  double rand_max_B = 
    powSpec_integral(B_reent, b_reent, highE_reent);
  
  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double spec_area = specA_area + specB_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_reent, a_reent, r);
      break;
    } else {
      // breakE<E
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = powSpec_integral_inv(B_reent, b_reent, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrPositronReentrant_0306::downwardFlux(){
  return 146.69;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the reentrant component
// in 0.6<theta_M<0.8.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray reentrant positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrPositronReentrant_0608::CrPositronReentrant_0608(){
  /*
   * Below 100 MeV
   *   j(E) = 2.61*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 2.12*10^-4*(E/GeV)^-3.09
   *          + 3.0*10^-4*(E/GeV)*exp(-(E/2GeV)^2) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization and spectral index for E<breakE
  A_reent = 2.61e-2;
  a_reent = 1.0;
  // Normalization, spectral index and cutoff for breakE<E
  B1_reent = 2.12e-4;
  b1_reent = 3.09;
  B2_reent = 3.0e-4;
  b2_reent = -1.0; // positive slope
  cutOff = 2.0;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrPositronReentrant_0608::~CrPositronReentrant_0608()
{
  ;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
double CrPositronReentrant_0608::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  double rand_max_A = 
    powSpec_integral(A_reent, a_reent, breakE);
  double rand_min_B1 = 
    powSpec_integral(B1_reent, b1_reent, breakE);
  double rand_max_B1 = 
    powSpec_integral(B1_reent, b1_reent, highE_reent);
  double rand_min_B2 = 
    cutOffPowSpec2_integral(B2_reent, b2_reent, cutOff, breakE);
  double rand_max_B2 = 
    cutOffPowSpec2_integral(B2_reent, b2_reent, cutOff, highE_reent);
  
  double specA_area = rand_max_A - rand_min_A;
  double specB1_area = rand_max_B1 - rand_min_B1;
  double specB2_area = rand_max_B2 - rand_min_B2;
  double spec_area = specA_area + specB1_area + specB2_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_reent, a_reent, r);
      break;
    } else if(rnd <= (specA_area+specB1_area) / spec_area){
      // breakE<E, powe-law component
      r = engine->flat() * (rand_max_B1 - rand_min_B1) + rand_min_B1;
      E = powSpec_integral_inv(B1_reent, b1_reent, r);
      break;
    } else {
      // breakE<E, cut off powe-law component
      r = engine->flat() * (rand_max_B2 - rand_min_B2) + rand_min_B2;
      E = cutOffPowSpec2_integral_inv(B2_reent, b2_reent, cutOff, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrPositronReentrant_0608::downwardFlux(){
  return 73.15;
}
//------------------------------------------------------------


//------------------------------------------------------------
// The random number generator for the reentrant component
// in 0.8<theta_M<0.9.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray reentrant positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrPositronReentrant_0809::CrPositronReentrant_0809(){
  /*
   * Below 100 MeV
   *   j(E) = 2.59*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 1.33*10^-4*(E/GeV)^-3.29
   *          + 1.6*10^-3*(E/GeV)^2.0*exp(-(E/1.6GeV)^3) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization and spectral index for E<breakE
  A_reent = 2.59e-2;
  a_reent = 1.0;
  // Normalization and spectral index for E>breakE
  B1_reent = 1.33e-4;
  b1_reent = 3.29;
  B2_reent = 1.6e-3;
  b2_reent = -2.0; // positive slope
  cutOff = 1.6;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrPositronReentrant_0809::~CrPositronReentrant_0809()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
double CrPositronReentrant_0809::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  double rand_max_A = 
    powSpec_integral(A_reent, a_reent, breakE);
  double rand_min_B1 = 
    powSpec_integral(B1_reent, b1_reent, breakE);
  double rand_max_B1 = 
    powSpec_integral(B1_reent, b1_reent, highE_reent);
  double rand_min_B2 = 
    cutOffPowSpec2_integral(B2_reent, b2_reent, cutOff, breakE);
  double rand_max_B2 = 
    cutOffPowSpec2_integral(B2_reent, b2_reent, cutOff, highE_reent);
  
  double specA_area = rand_max_A - rand_min_A;
  double specB1_area = rand_max_B1 - rand_min_B1;
  double specB2_area = rand_max_B2 - rand_min_B2;
  double spec_area = specA_area + specB1_area + specB2_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_reent, a_reent, r);
      break;
    } else if(rnd <= (specA_area+specB1_area) / spec_area){
      // breakE<E, powe-law component
      r = engine->flat() * (rand_max_B1 - rand_min_B1) + rand_min_B1;
      E = powSpec_integral_inv(B1_reent, b1_reent, r);
      break;
    } else {
      // breakE<E, cut off powe-law component
      r = engine->flat() * (rand_max_B2 - rand_min_B2) + rand_min_B2;
      E = cutOffPowSpec2_integral_inv(B2_reent, b2_reent, cutOff, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrPositronReentrant_0809::downwardFlux(){
  return 73.13;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the reentrant component
// in 0.9<theta_M<1.0.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray reentrant positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrPositronReentrant_0910::CrPositronReentrant_0910(){
  /*
   * Below 100 MeV
   *   j(E) = 4.09*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * 100-300 MeV
   *   j(E) = 3.82*10^-4*(E/GeV)^-3.03 [c/s/m^2/sr]
   * Above 300 MeV
   *   j(E) = 5.55*10^-3*(E/GeV)^(-1.0)*exp(-E/1.3GeV) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 484, 10
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization and spectral index for E<lowE_break
  A_reent = 4.09e-2;
  a_reent = 1.0;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_reent = 3.82e-4;
  b_reent = 3.03;
  // Normalization, spectral index and cut off for highE_break<E
  C_reent = 5.55e-3;
  c_reent = 1.0;
  cutOff = 1.3;
  // The spectrum breaks at 100 MeV and 300 MeV
  lowE_break = 0.1;
  highE_break = 0.3;
}
    
CrPositronReentrant_0910::~CrPositronReentrant_0910()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
double CrPositronReentrant_0910::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  double rand_max_B = 
    powSpec_integral(B_reent, b_reent, highE_break);
  double rand_min_C = 
    envelopeCutOffPowSpec_integral(C_reent, c_reent, highE_break);
  double rand_max_C = 
    envelopeCutOffPowSpec_integral(C_reent, c_reent, highE_reent);
  
  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double specC_area = rand_max_C - rand_min_C;
  double spec_area = specA_area + specB_area + specC_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<lowE_break
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_reent, a_reent, r);
      break;
    } else if(rnd <= (specA_area+specB_area) / spec_area){
      // lowE_break<E<highE_break
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = powSpec_integral_inv(B_reent, b_reent, r);
      break;
    } else {
      // highE_breakE<E
      r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
      E = envelopeCutOffPowSpec_integral_inv(C_reent, c_reent, r);
      if (engine->flat() * envelopeCutOffPowSpec(C_reent, c_reent, E)
	    < cutOffPowSpec(C_reent, c_reent, cutOff, E)){
	break;
      }
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrPositronReentrant_0910::downwardFlux(){
  return 120.73;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the reentrant component
// in 1.0<theta_M<1.1.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray reentrant positron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrPositronReentrant_1011::CrPositronReentrant_1011(){
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
  A_reent = 6.60e-2;
  a_reent = 1.0;
  // Normalization and spectral index for breakE<E
  B_reent = 9.11e-4;
  b_reent = 2.86;
  // The spectrum breaks at 100 MeV
  breakE = 0.1;
}
    
CrPositronReentrant_1011::~CrPositronReentrant_1011()
{
;
}

// returns energy obeying re-entrant cosmic-ray positron spectrum
double CrPositronReentrant_1011::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  double rand_max_A = 
    powSpec_integral(A_reent, a_reent, breakE);
  double rand_min_B = 
    powSpec_integral(B_reent, b_reent, breakE);
  double rand_max_B = 
    powSpec_integral(B_reent, b_reent, highE_reent);
  
  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double spec_area = specA_area + specB_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // E<breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_reent, a_reent, r);
      break;
    } else{
      // breakE<E
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = powSpec_integral_inv(B_reent, b_reent, r);
      break;
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrPositronReentrant_1011::downwardFlux(){
  return 187.45;
}
//------------------------------------------------------------
