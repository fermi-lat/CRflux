/**************************************************************************
 * CrProtonSubSplash.cc
 **************************************************************************
 * This file defines the subclasses of CrProtonSplash. 
 * They describe the cosmic-ray proton
 * splash spectrum sorted on geomagnetic latitude.
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

#include "CrProtonSubSplash.hh"

typedef  double G4double;

// private function definitions.
namespace {
  const G4double pi = 3.14159265358979323846264339;
  // rest energy (rest mass) of proton in units of GeV
  const G4double restE = 0.938;
  // lower and higher energy limit of secondary proton in units of GeV
  const G4double lowE_splash  = 0.01;
  const G4double highE_splash = 10.0;

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
  
  // integral of the envelope function of cutoff power-law
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
// The random number generator for the splash component
// in 0<theta_M<0.2.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray splash proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrProtonSplash_0002::CrProtonSplash_0002(){
  /*
   * Below 100 MeV
   *   j(E) = 1.44*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 7.55*10^-2*(E/GeV)^-0.336*exp(-E/0.787GeV) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_splash = 1.44e-2;
  a_splash = 1.0;
  // Normalization and spectral index for E>breakE
  B_splash = 7.55e-2;
  b_splash = 0.336;
  cutOff = 0.787;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrProtonSplash_0002::~CrProtonSplash_0002()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
double CrProtonSplash_0002::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  double rand_min_B = 
    envelopeCutOffPowSpec_integral(B_splash, b_splash, breakE);
  double rand_max_B = 
    envelopeCutOffPowSpec_integral(B_splash, b_splash, highE_splash);
  
  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double spec_area = specA_area + specB_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  while(1){
    rnd = engine->flat();
    if (rnd <= specA_area / spec_area){
      // spectrum below breakE
      r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
      E = powSpec_integral_inv(A_splash, a_splash, r);
      break;
    } else {
      // spectrum above breakE
      r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
      E = envelopeCutOffPowSpec_integral_inv(B_splash, b_splash, r);
      if (engine->flat() * envelopeCutOffPowSpec(B_splash, b_splash, E) 
          < cutOffPowSpec(B_splash, b_splash, cutOff, E)){
        break;
      }
    }
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrProtonSplash_0002::upwardFlux(){
  return 97.21;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the splash component
// in 0.2<theta_M<0.3.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray splash proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrProtonSplash_0203::CrProtonSplash_0203(){
  /*
   * Below 100 MeV
   *   j(E) = 9.30*10^-3*(E/GeV)^-1.0 [c/s/m^2/sr]
   * 100-600 MeV
   *   j(E) = 1.29*10^-2*(E/GeV)^-0.858 [c/s/m^2/sr]
   * Above 600 MeV
   *   j(E) = 5.78*10^-3*(E/GeV)^-2.43 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<lowE_break
  A_splash = 9.30e-3;
  a_splash = 1.0;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_splash = 1.29e-2;
  b_splash = 0.858;
  // Normalization and spectral index for E>highE_break
  C_splash = 5.78e-3;
  c_splash = 2.43;
  // The spectrum breaks at 100MeV and 600 MeV
  lowE_break = 0.1;
  highE_break = 0.6;
}
    
CrProtonSplash_0203::~CrProtonSplash_0203()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
double CrProtonSplash_0203::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  double rand_max_A = 
    powSpec_integral(A_splash, a_splash, lowE_break);
  double rand_min_B = 
    powSpec_integral(B_splash, b_splash, lowE_break);
  double rand_max_B = 
    powSpec_integral(B_splash, b_splash, highE_break);
  double rand_min_C = 
    powSpec_integral(C_splash, c_splash, highE_break);
  double rand_max_C = 
    powSpec_integral(C_splash, c_splash, highE_splash);

  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double specC_area = rand_max_C - rand_min_C;
  double spec_area = specA_area + specB_area + specC_area;

  double r, E; // E means energy in GeV
  double rnd;
  
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
double CrProtonSplash_0203::upwardFlux(){
  return 48.63;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the splash component
// in 0.3<theta_M<0.4.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray splash proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrProtonSplash_0304::CrProtonSplash_0304(){
  /*
   * Below 600 MeV
   *   j(E) = 8.88*10^-3*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 600 MeV
   *   j(E) = 4.34*10^-3*(E/GeV)^-2.4 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_splash = 8.88e-3;
  a_splash = 1.0;
  // Normalization and spectral index for breakE<E
  B_splash = 4.34e-3;
  b_splash = 2.4;
  // The spectrum breaks at 600MeV
  breakE = 0.6;
}
    
CrProtonSplash_0304::~CrProtonSplash_0304()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
double CrProtonSplash_0304::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  double rand_min_B = 
    powSpec_integral(B_splash, b_splash, breakE);
  double rand_max_B = 
    powSpec_integral(B_splash, b_splash, highE_splash);

  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double spec_area = specA_area + specB_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below breakE
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_splash, a_splash, r);
  } else {
    // spectrum above breakE
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_splash, b_splash, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrProtonSplash_0304::upwardFlux(){
  return 42.56;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the splash component
// in 0.4<theta_M<0.5.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray splash proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrProtonSplash_0405::CrProtonSplash_0405(){
  /*
   * Below 100 MeV
   *   j(E) = 9.85*10^-3*(E/GeV)^-1.0 [c/s/m^2/sr]
   * 100-600 MeV
   *   j(E) = 5.54*10^-3*(E/GeV)^-1.25 [c/s/m^2/sr]
   * Above 600 MeV
   *   j(E) = 3.08*10^-3*(E/GeV)^-2.4 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<lowE_break
  A_splash = 9.85e-3;
  a_splash = 1.0;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_splash = 5.54e-3;
  b_splash = 1.25;
  // Normalization and spectral index for E>highE_break
  C_splash = 3.08e-3;
  c_splash = 2.4;
  // The spectrum breaks at 100MeV and 600 MeV
  lowE_break = 0.1;
  highE_break = 0.6;
}
    
CrProtonSplash_0405::~CrProtonSplash_0405()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
double CrProtonSplash_0405::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  double rand_max_A = 
    powSpec_integral(A_splash, a_splash, lowE_break);
  double rand_min_B = 
    powSpec_integral(B_splash, b_splash, lowE_break);
  double rand_max_B = 
    powSpec_integral(B_splash, b_splash, highE_break);
  double rand_min_C = 
    powSpec_integral(C_splash, c_splash, highE_break);
  double rand_max_C = 
    powSpec_integral(C_splash, c_splash, highE_splash);

  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double specC_area = rand_max_C - rand_min_C;
  double spec_area = specA_area + specB_area + specC_area;

  double r, E; // E means energy in GeV
  double rnd;
  
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
double CrProtonSplash_0405::upwardFlux(){
  return 41.31;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the splash component
// in 0.5<theta_M<0.6.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray splash proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrProtonSplash_0506::CrProtonSplash_0506(){
  /*
   * Below 100 MeV
   *   j(E) = 9.57*10^-3*(E/GeV)^-1.0 [c/s/m^2/sr]
   * 100-400 MeV
   *   j(E) = 6.18*10^-3*(E/GeV)^-1.19 [c/s/m^2/sr]
   * Above 400 MeV
   *   j(E) = 2.11*10^-3*(E/GeV)^-2.37 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<lowE_break
  A_splash = 9.57e-3;
  a_splash = 1.0;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_splash = 6.18e-3;
  b_splash = 1.19;
  // Normalization and spectral index for E>highE_break
  C_splash = 2.11e-3;
  c_splash = 2.37;
  // The spectrum breaks at 100MeV and 600 MeV
  lowE_break = 0.1;
  highE_break = 0.6;
}
    
CrProtonSplash_0506::~CrProtonSplash_0506()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
double CrProtonSplash_0506::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  double rand_max_A = 
    powSpec_integral(A_splash, a_splash, lowE_break);
  double rand_min_B = 
    powSpec_integral(B_splash, b_splash, lowE_break);
  double rand_max_B = 
    powSpec_integral(B_splash, b_splash, highE_break);
  double rand_min_C = 
    powSpec_integral(C_splash, c_splash, highE_break);
  double rand_max_C = 
    powSpec_integral(C_splash, c_splash, highE_splash);

  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double specC_area = rand_max_C - rand_min_C;
  double spec_area = specA_area + specB_area + specC_area;

  double r, E; // E means energy in GeV
  double rnd;
  
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
double CrProtonSplash_0506::upwardFlux(){
  return 39.02;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the splash component
// in 0.6<theta_M<0.7.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray splash proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrProtonSplash_0607::CrProtonSplash_0607(){
  /*
   * Below 300 MeV
   *   j(E) = 1.15*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 300 MeV
   *   j(E) = 2.49*10^-3*(E/GeV)^-2.27 [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_splash = 1.15e-2;
  a_splash = 1.0;
  // Normalization and spectral index for breakE<E
  B_splash = 2.49e-3;
  b_splash = 2.27;
  // The spectrum breaks at 600MeV
  breakE = 0.3;
}
    
CrProtonSplash_0607::~CrProtonSplash_0607()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
double CrProtonSplash_0607::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  double rand_min_B = 
    powSpec_integral(B_splash, b_splash, breakE);
  double rand_max_B = 
    powSpec_integral(B_splash, b_splash, highE_splash);

  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double spec_area = specA_area + specB_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below breakE
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_splash, a_splash, r);
  } else {
    // spectrum above breakE
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_splash, b_splash, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrProtonSplash_0607::upwardFlux(){
  return 48.05;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the splash component
// in 0.7<theta_M<0.8.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray splash proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrProtonSplash_0708::CrProtonSplash_0708(){
  /*
   * Below 100 MeV
   *   j(E) = 1.82*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 4.99*10^-3*(E/GeV)^-1.92*exp(-(E/0.081GeV)^-0.92) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_splash = 1.80e-2;
  a_splash = 1.0;
  // Normalization and spectral index for breakE<E
  B_splash = 4.99e-3;
  b_splash = 1.92;
  cutOff = 0.081;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrProtonSplash_0708::~CrProtonSplash_0708()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
double CrProtonSplash_0708::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  double rand_min_B = 
    cutOffPowSpec2_integral(B_splash, b_splash, cutOff, breakE);
  double rand_max_B = 
    cutOffPowSpec2_integral(B_splash, b_splash, cutOff, highE_splash);

  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double spec_area = specA_area + specB_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below breakE
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_splash, a_splash, r);
  } else {
    // spectrum above breakE
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = cutOffPowSpec2_integral_inv(B_splash, b_splash, cutOff, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrProtonSplash_0708::upwardFlux(){
  return 71.52;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the splash component
// in 0.8<theta_M<0.9.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray splash proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrProtonSplash_0809::CrProtonSplash_0809(){
  /*
   * Below 100 MeV
   *   j(E) = 2.29*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 1.70*10^-2*(E/GeV)^-1.83*exp(-(E/0.178GeV)^-0.83) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_splash = 2.29e-2;
  a_splash = 1.0;
  // Normalization and spectral index for breakE<E
  B_splash = 1.70e-2;
  b_splash = 1.83;
  cutOff = 0.178;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrProtonSplash_0809::~CrProtonSplash_0809()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
double CrProtonSplash_0809::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  double rand_min_B = 
    cutOffPowSpec2_integral(B_splash, b_splash, cutOff, breakE);
  double rand_max_B = 
    cutOffPowSpec2_integral(B_splash, b_splash, cutOff, highE_splash);

  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double spec_area = specA_area + specB_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below breakE
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_splash, a_splash, r);
  } else {
    // spectrum above breakE
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = cutOffPowSpec2_integral_inv(B_splash, b_splash, cutOff, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrProtonSplash_0809::upwardFlux(){
  return 118.46;
}
//------------------------------------------------------------

//------------------------------------------------------------
// The random number generator for the splash component
// in 0.9<theta_M<1.0.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray splash proton
// ("energy" method) and the other returns 
// the energy integrated downward flux ("upwardFlux" method) 
CrProtonSplash_0910::CrProtonSplash_0910(){
  /*
   * Below 100 MeV
   *   j(E) = 4.36*10^-2*(E/GeV)^-1.0 [c/s/m^2/sr]
   * Above 100 MeV
   *   j(E) = 3.65*10^-2*(E/GeV)^-1.98*exp(-(E/0.211GeV)^-0.98) [c/s/m^2/sr]
   * reference:
   *   AMS data, Alcaratz et al. 2000, Phys. Let. B 490, 27
   * Above 100 MeV, we modeled AMS data with analytic function.
   * Below 100 MeV, we do not have enouth information and just
   * extrapolated the spectrum down to 10 MeV with E^-1.
   */
  
  // Normalization, spectral index, and cutoff for E<breakE
  A_splash = 4.36e-2;
  a_splash = 1.0;
  // Normalization and spectral index for breakE<E
  B_splash = 3.65e-2;
  b_splash = 1.98;
  cutOff = 0.211;
  // The spectrum breaks at 100MeV
  breakE = 0.1;
}
    
CrProtonSplash_0910::~CrProtonSplash_0910()
{
;
}

// returns energy obeying re-entrant cosmic-ray proton spectrum
double CrProtonSplash_0910::energy(HepRandomEngine* engine){

  double rand_min_A = 
    powSpec_integral(A_splash, a_splash, lowE_splash);
  double rand_max_A = 
    powSpec_integral(A_splash, a_splash, breakE);
  double rand_min_B = 
    cutOffPowSpec2_integral(B_splash, b_splash, cutOff, breakE);
  double rand_max_B = 
    cutOffPowSpec2_integral(B_splash, b_splash, cutOff, highE_splash);

  double specA_area = rand_max_A - rand_min_A;
  double specB_area = rand_max_B - rand_min_B;
  double spec_area = specA_area + specB_area;

  double r, E; // E means energy in GeV
  double rnd;
  
  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below breakE
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_splash, a_splash, r);
  } else {
    // spectrum above breakE
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = cutOffPowSpec2_integral_inv(B_splash, b_splash, cutOff, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
double CrProtonSplash_0910::upwardFlux(){
  return 246.24;
}
//------------------------------------------------------------


