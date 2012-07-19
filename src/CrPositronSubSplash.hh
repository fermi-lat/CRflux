//
// CrPositronSubSplash.hh
//

//$Header$

#ifndef CrPositronSubSplash_H
#define CrPositronSubSplash_H

#include <utility>
#include <string>

#include "CrSpectrum.hh"

// Forward declaration:
class CLHEP::HepRandomEngine;

class CrPositronSplash_0001
{
public:
  CrPositronSplash_0001();
  ~CrPositronSplash_0001();

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

class CrPositronSplash_0102
{
public:
  CrPositronSplash_0102();
  ~CrPositronSplash_0102();

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

class CrPositronSplash_0203
{
public:
  CrPositronSplash_0203();
  ~CrPositronSplash_0203();

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

class CrPositronSplash_0304
{
public:
  CrPositronSplash_0304();
  ~CrPositronSplash_0304();

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

class CrPositronSplash_0405
{
public:
  CrPositronSplash_0405();
  ~CrPositronSplash_0405();

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

class CrPositronSplash_0506
{
public:
  CrPositronSplash_0506();
  ~CrPositronSplash_0506();

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

class CrPositronSplash_0611
{
public:
  CrPositronSplash_0611();
  ~CrPositronSplash_0611();

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


#endif // CrPositronSubSplash_H
