//
// CrPositronSubSplash.hh
//

//$Header$

#ifndef CrPositronSubSplash_H
#define CrPositronSubSplash_H

#include <utility>
#include <string>

// Forward declaration:
class HepRandomEngine;

class CrPositronSplash_0003
{
public:
  CrPositronSplash_0003();
  ~CrPositronSplash_0003();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization, spectral index and cutoff for breakE<E
  double B_splash, b_splash;
  double cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrPositronSplash_0306
{
public:
  CrPositronSplash_0306();
  ~CrPositronSplash_0306();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization and spectral index for breakE<E
  double B_splash, b_splash;
  double cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrPositronSplash_0608
{
public:
  CrPositronSplash_0608();
  ~CrPositronSplash_0608();

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
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrPositronSplash_0809
{
public:
  CrPositronSplash_0809();
  ~CrPositronSplash_0809();

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
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrPositronSplash_0910
{
public:
  CrPositronSplash_0910();
  ~CrPositronSplash_0910();

private:
  // Normalization and spectral index for E<lowE_break
  double A_splash, a_splash;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_splash, b_splash;
  // Normalization, spectral index and cut off for highE_break<E
  double C_splash, c_splash, cutOff;
  // The energies where spectrum breaks
  double lowE_break, highE_break;

public:
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrPositronSplash_1011
{
public:
  CrPositronSplash_1011();
  ~CrPositronSplash_1011();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization and spectral index for breakE<E
  double B_splash, b_splash;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};


#endif // CrPositronSubSplash_H
