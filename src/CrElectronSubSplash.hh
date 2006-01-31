//
// CrElectronSubSplash.hh
//

//$Header$

#ifndef CrElectronSubSplash_H
#define CrElectronSubSplash_H

#include <utility>
#include <string>

// Forward declaration:
class CLHEP::HepRandomEngine;

class CrElectronSplash_0003
{
public:
  CrElectronSplash_0003();
  ~CrElectronSplash_0003();

private:
  // Normalization and spectral index for E<lowE_break
  double A_splash, a_splash;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_splash, b_splash;
  // Normalization and spectral index for highE_break<E
  double C_splash, c_splash;
  // The energies where spectrum breaks
  double lowE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double upwardFlux();
};

class CrElectronSplash_0306
{
public:
  CrElectronSplash_0306();
  ~CrElectronSplash_0306();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization, spectral index and cut off for breakE<E
  double B_splash, b_splash;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double upwardFlux();
};

class CrElectronSplash_0608
{
public:
  CrElectronSplash_0608();
  ~CrElectronSplash_0608();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization, spectral index and cut off for breakE<E
  double B1_splash, b1_splash;
  double B2_splash, b2_splash;
  double cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double upwardFlux();
};

class CrElectronSplash_0809
{
public:
  CrElectronSplash_0809();
  ~CrElectronSplash_0809();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization, spectral index and cut off for breakE<E
  double B1_splash, b1_splash;
  double B2_splash, b2_splash;
  double cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double upwardFlux();
};

class CrElectronSplash_0910
{
public:
  CrElectronSplash_0910();
  ~CrElectronSplash_0910();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization, spectral index and cut off for breakE<E
  double B_splash, b_splash;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double upwardFlux();
};

class CrElectronSplash_1011
{
public:
  CrElectronSplash_1011();
  ~CrElectronSplash_1011();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization and spectral index for breakE<E
  double B_splash, b_splash;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double upwardFlux();
};


#endif // CrElectronSubSplash_H
