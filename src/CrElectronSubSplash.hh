//
// CrElectronSubSplash.hh
//

//$Header$

#ifndef CrElectronSubSplash_H
#define CrElectronSubSplash_H

#include <utility>
#include <string>

#include "CrSpectrum.hh"

// Forward declaration:
class CLHEP::HepRandomEngine;

class CrElectronSplash_0001
{
public:
  CrElectronSplash_0001();
  ~CrElectronSplash_0001();

private:
  // Normalization and spectral index for E<lowE_break
  double A_splash, a_splash;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_splash, b_splash;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_splash, c_splash;
  // Normalization and spectral index for highE_break<E
  double D_splash, d_splash;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronSplash_0102
{
public:
  CrElectronSplash_0102();
  ~CrElectronSplash_0102();

private:
  // Normalization and spectral index for E<lowE_break
  double A_splash, a_splash;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_splash, b_splash;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_splash, c_splash;
  // Normalization and spectral index for highE_break<E
  double D_splash, d_splash;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronSplash_0203
{
public:
  CrElectronSplash_0203();
  ~CrElectronSplash_0203();

private:
  // Normalization and spectral index for E<lowE_break
  double A_splash, a_splash;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_splash, b_splash;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_splash, c_splash;
  // Normalization and spectral index for highE_break<E
  double D_splash, d_splash;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronSplash_0304
{
public:
  CrElectronSplash_0304();
  ~CrElectronSplash_0304();

private:
  // Normalization and spectral index for E<lowE_break
  double A_splash, a_splash;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_splash, b_splash;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_splash, c_splash;
  // Normalization and spectral index for highE_break<E
  double D_splash, d_splash;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronSplash_0405
{
public:
  CrElectronSplash_0405();
  ~CrElectronSplash_0405();

private:
  // Normalization and spectral index for E<lowE_break
  double A_splash, a_splash;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_splash, b_splash;
  // Normalization and spectral index for highE_break<E
  double C_splash, c_splash;
  double lowE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronSplash_0506
{
public:
  CrElectronSplash_0506();
  ~CrElectronSplash_0506();

private:
  // Normalization and spectral index for E<lowE_break
  double A_splash, a_splash;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_splash, b_splash;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_splash, c_splash;
  // Normalization and spectral index for highE_break<E
  double D_splash, d_splash;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronSplash_0611
{
public:
  CrElectronSplash_0611();
  ~CrElectronSplash_0611();

private:
  // Normalization and spectral index for E<lowE_break
  double A_splash, a_splash;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_splash, b_splash;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_splash, c_splash;
  // Normalization and spectral index for highE_break<E
  double D_splash, d_splash;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};


#endif // CrElectronSubSplash_H
