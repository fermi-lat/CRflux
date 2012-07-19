/**************************************************************************
 * CrElectronSubReentrant.cc
 **************************************************************************
 * This file defines the subclasses of CrElectronReentrant. 
 * They describe the cosmic-ray electron
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

#include "CrElectronSubReentrant.hh"

typedef double G4double;

// private function definitions.
namespace {
  // rest energy (rest mass) of electron in units of GeV
  const G4double restE = 5.11e-4;
  // lower and higher energy limit of secondary electron in units of GeV
  const G4double lowE_reent  = 0.01;
  const G4double highE_reent = 10.0;
  
  const G4double PosToEle_0001 = 4.8; // e+/e- of secondary in 0.0<theta_M<0.1
  const G4double PosToEle_0102 = 4.2; // e+/e- of secondary in 0.1<theta_M<0.2
  const G4double PosToEle_0203 = 3.8; // e+/e- of secondary in 0.2<theta_M<0.3
  const G4double PosToEle_0304 = 2.6; // e+/e- of secondary in 0.3<theta_M<0.4
  const G4double PosToEle_0405 = 1.8; // e+/e- of secondary in 0.4<theta_M<0.5
  const G4double PosToEle_0506 = 1.0; // e+/e- of secondary in 0.5<theta_M<0.6
  const G4double PosToEle_0611 = 1.0; // e+/e- of secondary in 0.6<theta_M<1.1

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
// The random number generator for the downward component
// in 0<theta_M<0.1.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward electron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrElectronReentrant_0001::CrElectronReentrant_0001(){
  /*
   * Below 100 MeV
   *   j(E) = 0.45*(E/100MeV)^-2.0 [c/s/m^2/sr/MeV]
   * 100 MeV - 400 MeV
   *   j(E) = 0.45*(E/100MeV)^-1.5 [c/s/m^2/sr/MeV]
   * 400 MeV - 3 GeV
   *   j(E) = 0.056*(E/400MeV)^-2.5 [c/s/m^2/sr/MeV]
   * Above 3 GeV
   *   j(E) = 3.65e-4*pow(E/3000MeV)^-3.6 [c/s/m^2/sr/MeV]
   * reference:
   *   LAT measurement of 2ndary e- + e+
   */
  
  // Normalization and spectral index E<lowE_break
  A_reent = 0.45*pow(1000./100., -2.0)/(1+PosToEle_0001);
  a_reent = 2.0;
  // Normalization and spectral index for lowE_break<E<midE_break
  B_reent = 0.45*pow(1000./100., -1.5)/(1+PosToEle_0001);
  b_reent = 1.5;
  // Normalization and spectral index for midE_break<E<highE_break
  C_reent = 0.056*pow(1000./400., -2.5)/(1+PosToEle_0001);
  c_reent = 2.5;
  // Normalization and spectral index for E>highE_break
  D_reent = 3.65e-4*pow(1000./3000.0, -3.6)/(1+PosToEle_0001);
  d_reent = 3.6;
  // The spectrum breaks at 100 MeV, 400 MeV and 3 GeV
  lowE_break = 0.1;
  midE_break = 0.4;
  highE_break = 3.0;
}
    
CrElectronReentrant_0001::~CrElectronReentrant_0001()
{
;
}

// returns energy obeying re-entrant cosmic-ray electron spectrum
G4double CrElectronReentrant_0001::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_D = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_D = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double specD_area = rand_max_D - rand_min_D;
  G4double spec_area = specA_area + specB_area + specC_area + specD_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;

  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<midE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else if (rnd <= (specA_area+specB_area+specC_area) / spec_area){
    // spectrum in midE_break<E<highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_D - rand_min_D) + rand_min_D;
    E = powSpec_integral_inv(D_reent, d_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrElectronReentrant_0001::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_4 = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_4 = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3)+(rand_max_4-rand_min_4));
}
//------------------------------------------------------------


//------------------------------------------------------------
// The random number generator for the downward component
// in 0.1<theta_M<0.2.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward electron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrElectronReentrant_0102::CrElectronReentrant_0102(){
  /*
   * Below 100 MeV
   *   j(E) = 0.45*(E/100MeV)^-2.0 [c/s/m^2/sr/MeV]
   * 100 MeV - 400 MeV
   *   j(E) = 0.45*(E/100MeV)^-1.5 [c/s/m^2/sr/MeV]
   * 400 MeV - 1 GeV
   *   j(E) = 0.056*(E/400MeV)^-2.5 [c/s/m^2/sr/MeV]
   * Above 1 GeV
   *   j(E) = 0.0056*pow(E/1000MeV)^-2.9 [c/s/m^2/sr/MeV]
   * reference:
   *   LAT measurement of 2ndary e- + e+
   */
  
  // Normalization and spectral index E<lowE_break
  A_reent = 0.45*pow(1000./100., -2.0)/(1+PosToEle_0102);
  a_reent = 2.0;
  // Normalization and spectral index for lowE_break<E<midE_break
  B_reent = 0.45*pow(1000./100., -1.5)/(1+PosToEle_0102);
  b_reent = 1.5;
  // Normalization and spectral index for midE_break<E<highE_break
  C_reent = 0.056*pow(1000./400., -2.5)/(1+PosToEle_0102);
  c_reent = 2.5;
  // Normalization and spectral index for E>highE_break
  D_reent = 0.0056*pow(1000./1000.0, -2.9)/(1+PosToEle_0102);
  d_reent = 2.9;
  // The spectrum breaks at 100 MeV, 400 MeV and 3 GeV
  lowE_break = 0.1;
  midE_break = 0.4;
  highE_break = 1.0;
}
    
CrElectronReentrant_0102::~CrElectronReentrant_0102()
{
;
}

// returns energy obeying re-entrant cosmic-ray electron spectrum
G4double CrElectronReentrant_0102::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_D = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_D = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double specD_area = rand_max_D - rand_min_D;
  G4double spec_area = specA_area + specB_area + specC_area + specD_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;

  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<midE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else if (rnd <= (specA_area+specB_area+specC_area) / spec_area){
    // spectrum in midE_break<E<highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_D - rand_min_D) + rand_min_D;
    E = powSpec_integral_inv(D_reent, d_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrElectronReentrant_0102::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_4 = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_4 = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3)+(rand_max_4-rand_min_4));
}
//------------------------------------------------------------


//------------------------------------------------------------
// The random number generator for the downward component
// in 0.2<theta_M<0.3.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward electron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrElectronReentrant_0203::CrElectronReentrant_0203(){
  /*
   * Below 100 MeV
   *   j(E) = 0.45*(E/100MeV)^-2.0 [c/s/m^2/sr/MeV]
   * 100 MeV - 300 MeV
   *   j(E) = 0.45*(E/100MeV)^-1.5 [c/s/m^2/sr/MeV]
   * 300 MeV - 400 MeV
   *   j(E) = 0.086*(E/300MeV)^-1.8 [c/s/m^2/sr/MeV]
   * Above 400 MeV
   *   j(E) = 0.051*pow(E/400MeV)^-2.8 [c/s/m^2/sr/MeV]
   * reference:
   *   LAT measurement of 2ndary e- + e+
   */
  
  // Normalization and spectral index E<lowE_break
  A_reent = 0.45*pow(1000./100., -2.0)/(1+PosToEle_0203);
  a_reent = 2.0;
  // Normalization and spectral index for lowE_break<E<midE_break
  B_reent = 0.45*pow(1000./100., -1.5)/(1+PosToEle_0203);
  b_reent = 1.5;
  // Normalization and spectral index for midE_break<E<highE_break
  C_reent = 0.086*pow(1000./300., -1.8)/(1+PosToEle_0203);
  c_reent = 1.8;
  // Normalization and spectral index for E>highE_break
  D_reent = 0.051*pow(1000./400.0, -2.8)/(1+PosToEle_0203);
  d_reent = 2.8;
  // The spectrum breaks at 100 MeV, 300 MeV and 400 MeV
  lowE_break = 0.1;
  midE_break = 0.3;
  highE_break = 0.4;
}
    
CrElectronReentrant_0203::~CrElectronReentrant_0203()
{
;
}

// returns energy obeying re-entrant cosmic-ray electron spectrum
G4double CrElectronReentrant_0203::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_D = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_D = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double specD_area = rand_max_D - rand_min_D;
  G4double spec_area = specA_area + specB_area + specC_area + specD_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;

  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<midE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else if (rnd <= (specA_area+specB_area+specC_area) / spec_area){
    // spectrum in midE_break<E<highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_D - rand_min_D) + rand_min_D;
    E = powSpec_integral_inv(D_reent, d_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrElectronReentrant_0203::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_4 = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_4 = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3)+(rand_max_4-rand_min_4));
}
//------------------------------------------------------------


//------------------------------------------------------------
// The random number generator for the downward component
// in 0.3<theta_M<0.4.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward electron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrElectronReentrant_0304::CrElectronReentrant_0304(){
  /*
   * Below 100 MeV
   *   j(E) = 0.45*(E/100MeV)^-2.0 [c/s/m^2/sr/MeV]
   * 100 MeV - 300 MeV
   *   j(E) = 0.45*(E/100MeV)^-1.6 [c/s/m^2/sr/MeV]
   * 300 MeV - 600 MeV
   *   j(E) = 0.078*(E/300MeV)^-2.5 [c/s/m^2/sr/MeV]
   * Above 600 MeV
   *   j(E) = 0.0137*pow(E/600MeV)^-2.8 [c/s/m^2/sr/MeV]
   * reference:
   *   LAT measurement of 2ndary e- + e+
   */
  
  // Normalization and spectral index E<lowE_break
  A_reent = 0.45*pow(1000./100., -2.0)/(1+PosToEle_0304);
  a_reent = 2.0;
  // Normalization and spectral index for lowE_break<E<midE_break
  B_reent = 0.45*pow(1000./100., -1.6)/(1+PosToEle_0304);
  b_reent = 1.6;
  // Normalization and spectral index for midE_break<E<highE_break
  C_reent = 0.078*pow(1000./300., -2.5)/(1+PosToEle_0304);
  c_reent = 2.5;
  // Normalization and spectral index for E>highE_break
  D_reent = 0.0137*pow(1000./600.0, -2.8)/(1+PosToEle_0304);
  d_reent = 2.8;
  // The spectrum breaks at 100 MeV, 300 MeV and 600 MeV
  lowE_break = 0.1;
  midE_break = 0.3;
  highE_break = 0.6;
}
    
CrElectronReentrant_0304::~CrElectronReentrant_0304()
{
;
}

// returns energy obeying re-entrant cosmic-ray electron spectrum
G4double CrElectronReentrant_0304::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_D = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_D = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double specD_area = rand_max_D - rand_min_D;
  G4double spec_area = specA_area + specB_area + specC_area + specD_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;

  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<midE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else if (rnd <= (specA_area+specB_area+specC_area) / spec_area){
    // spectrum in midE_break<E<highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_D - rand_min_D) + rand_min_D;
    E = powSpec_integral_inv(D_reent, d_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrElectronReentrant_0304::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_4 = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_4 = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3)+(rand_max_4-rand_min_4));
}
//------------------------------------------------------------


//------------------------------------------------------------
// The random number generator for the downward component
// in 0.4<theta_M<0.5.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward electron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrElectronReentrant_0405::CrElectronReentrant_0405(){
  /*
   * Below 100 MeV
   *   j(E) = 0.5*(E/100MeV)^-2.0 [c/s/m^2/sr/MeV]
   * 100 MeV - 300 MeV
   *   j(E) = 0.5*(E/100MeV)^-1.7 [c/s/m^2/sr/MeV]
   * Above 300 MeV
   *   j(E) = 0.077*pow(E/300MeV)^-2.8 [c/s/m^2/sr/MeV]
   * reference:
   *   LAT measurement of 2ndary e- + e+
   */
  
  // Normalization and spectral index E<lowE_break
  A_reent = 0.5*pow(1000./100., -2.0)/(1+PosToEle_0405);
  a_reent = 2.0;
  // Normalization and spectral index for lowE_break<E<highE_break
  B_reent = 0.5*pow(1000./100., -1.7)/(1+PosToEle_0405);
  b_reent = 1.7;
  // Normalization and spectral index for E>highE_break
  C_reent = 0.077*pow(1000./300., -2.8)/(1+PosToEle_0405);
  c_reent = 2.8;
  // The spectrum breaks at 100 MeV, 400 MeV and 3 GeV
  lowE_break = 0.1;
  highE_break = 0.3;
}
    
CrElectronReentrant_0405::~CrElectronReentrant_0405()
{
;
}

// returns energy obeying re-entrant cosmic-ray electron spectrum
G4double CrElectronReentrant_0405::energy(CLHEP::HepRandomEngine* engine){

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
G4double CrElectronReentrant_0405::downwardFlux(){
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

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3));
}
//------------------------------------------------------------


//------------------------------------------------------------
// The random number generator for the downward component
// in 0.5<theta_M<0.6.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward electron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrElectronReentrant_0506::CrElectronReentrant_0506(){
  /*
   * Below 100 MeV
   *   j(E) = 0.6*(E/100MeV)^-2.0 [c/s/m^2/sr/MeV]
   * 100 MeV - 300 MeV
   *   j(E) = 0.6*(E/100MeV)^-1.9 [c/s/m^2/sr/MeV]
   * 300 MeV - 1.5 GeV
   *   j(E) = 0.074*(E/300MeV)^-3.0 [c/s/m^2/sr/MeV]
   * Above 1.5 GeV
   *   j(E) = 0.00059*pow(E/1500MeV)^-2.3 [c/s/m^2/sr/MeV]
   * reference:
   *   LAT measurement of 2ndary e- + e+
   */
  
  // Normalization and spectral index E<lowE_break
  A_reent = 0.6*pow(1000./100., -2.0)/(1+PosToEle_0506);
  a_reent = 2.0;
  // Normalization and spectral index for lowE_break<E<midE_break
  B_reent = 0.6*pow(1000./100., -1.9)/(1+PosToEle_0506);
  b_reent = 1.9;
  // Normalization and spectral index for midE_break<E<highE_break
  C_reent = 0.074*pow(1000./300., -3.0)/(1+PosToEle_0506);
  c_reent = 3.0;
  // Normalization and spectral index for E>highE_break
  D_reent = 0.00059*pow(1000./1500.0, -2.3)/(1+PosToEle_0506);
  d_reent = 2.3;
  // The spectrum breaks at 100 MeV, 400 MeV and 3 GeV
  lowE_break = 0.1;
  midE_break = 0.3;
  highE_break = 1.5;
}
    
CrElectronReentrant_0506::~CrElectronReentrant_0506()
{
;
}

// returns energy obeying re-entrant cosmic-ray electron spectrum
G4double CrElectronReentrant_0506::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_D = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_D = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double specD_area = rand_max_D - rand_min_D;
  G4double spec_area = specA_area + specB_area + specC_area + specD_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;

  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<midE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else if (rnd <= (specA_area+specB_area+specC_area) / spec_area){
    // spectrum in midE_break<E<highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_D - rand_min_D) + rand_min_D;
    E = powSpec_integral_inv(D_reent, d_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrElectronReentrant_0506::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_4 = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_4 = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3)+(rand_max_4-rand_min_4));
}
//------------------------------------------------------------


//------------------------------------------------------------
// The random number generator for the downward component
// in 0.6<theta_M<1.1.
// This class has two methods. 
// One returns the kinetic energy of cosmic-ray downward electron
// ("energy" method) and the other returns 
// the energy integrated downward flux ("downwardFlux" method) 
CrElectronReentrant_0611::CrElectronReentrant_0611(){
  /*
   * Below 100 MeV
   *   j(E) = 0.65*(E/100MeV)^-2.0 [c/s/m^2/sr/MeV]
   * 100 MeV - 300 MeV
   *   j(E) = 0.65*(E/100MeV)^-1.9 [c/s/m^2/sr/MeV]
   * 300 MeV - 1.2 GeV
   *   j(E) = 0.08*(E/300MeV)^-3.2 [c/s/m^2/sr/MeV]
   * Above 1.2 GeV
   *   j(E) = 9e-4*pow(E/1200MeV)^-1.8 [c/s/m^2/sr/MeV]
   * reference:
   *   LAT measurement of 2ndary e- + e+
   */
  
  // Normalization and spectral index E<lowE_break
  A_reent = 0.65*pow(1000./100., -2.0)/(1+PosToEle_0611);
  a_reent = 2.0;
  // Normalization and spectral index for lowE_break<E<midE_break
  B_reent = 0.65*pow(1000./100., -1.9)/(1+PosToEle_0611);
  b_reent = 1.9;
  // Normalization and spectral index for midE_break<E<highE_break
  C_reent = 0.08*pow(1000./300., -3.2)/(1+PosToEle_0611);
  c_reent = 3.2;
  // Normalization and spectral index for E>highE_break
  D_reent = 9e-4*pow(1000./1200.0, -1.8)/(1+PosToEle_0611);
  d_reent = 1.8;
  // The spectrum breaks at 100 MeV, 400 MeV and 3 GeV
  lowE_break = 0.1;
  midE_break = 0.3;
  highE_break = 1.2;
}
    
CrElectronReentrant_0611::~CrElectronReentrant_0611()
{
;
}

// returns energy obeying re-entrant cosmic-ray electron spectrum
G4double CrElectronReentrant_0611::energy(CLHEP::HepRandomEngine* engine){

  G4double rand_min_A = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_A = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_B = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_B = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_C = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_C = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_D = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_D = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  G4double specA_area = rand_max_A - rand_min_A;
  G4double specB_area = rand_max_B - rand_min_B;
  G4double specC_area = rand_max_C - rand_min_C;
  G4double specD_area = rand_max_D - rand_min_D;
  G4double spec_area = specA_area + specB_area + specC_area + specD_area;

  G4double r, E; // E means energy in GeV
  G4double rnd;

  rnd = engine->flat();
  if (rnd <= specA_area / spec_area){
    // spectrum below lowE_break
    r = engine->flat() * (rand_max_A - rand_min_A) + rand_min_A;
    E = powSpec_integral_inv(A_reent, a_reent, r);
  } else if (rnd <= (specA_area+specB_area) / spec_area){
    // spectrum in lowE_break<E<midE_break
    r = engine->flat() * (rand_max_B - rand_min_B) + rand_min_B;
    E = powSpec_integral_inv(B_reent, b_reent, r);
  } else if (rnd <= (specA_area+specB_area+specC_area) / spec_area){
    // spectrum in midE_break<E<highE_break
    r = engine->flat() * (rand_max_C - rand_min_C) + rand_min_C;
    E = powSpec_integral_inv(C_reent, c_reent, r);
  } else {
    // spectrum above highE_break
    r = engine->flat() * (rand_max_D - rand_min_D) + rand_min_D;
    E = powSpec_integral_inv(D_reent, d_reent, r);
  }
  return E;
}

// returns energy integrated downward flux in c/s/m^2/sr
G4double CrElectronReentrant_0611::downwardFlux(){
  G4double rand_min_1 = 
    powSpec_integral(A_reent, a_reent, lowE_reent);
  G4double rand_max_1 = 
    powSpec_integral(A_reent, a_reent, lowE_break);
  G4double rand_min_2 = 
    powSpec_integral(B_reent, b_reent, lowE_break);
  G4double rand_max_2 = 
    powSpec_integral(B_reent, b_reent, midE_break);
  G4double rand_min_3 = 
    powSpec_integral(C_reent, c_reent, midE_break);
  G4double rand_max_3 = 
    powSpec_integral(C_reent, c_reent, highE_break);
  G4double rand_min_4 = 
    powSpec_integral(D_reent, d_reent, highE_break);
  G4double rand_max_4 = 
    powSpec_integral(D_reent, d_reent, highE_reent);

  // Original model function is given in "/MeV" and the energy in "GeV".
  // This is why 1000.* is required below.

  return 1000.*((rand_max_1-rand_min_1)+(rand_max_2-rand_min_2)+(rand_max_3-rand_min_3)+(rand_max_4-rand_min_4));
}
//------------------------------------------------------------



