//
// CrProtonSubSplash.hh
//

//$Header$

#ifndef CrProtonSubSplash_H
#define CrProtonSubSplash_H

#include <utility>
#include <string>

// Forward declaration:
class HepRandomEngine;

class CrProtonSplash_0002
{
public:
  CrProtonSplash_0002();
  ~CrProtonSplash_0002();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization, spectral index and cutoff for E>breakE
  double B_splash, b_splash, cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrProtonSplash_0203
{
public:
  CrProtonSplash_0203();
  ~CrProtonSplash_0203();

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
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrProtonSplash_0304
{
public:
  CrProtonSplash_0304();
  ~CrProtonSplash_0304();

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

class CrProtonSplash_0405
{
public:
  CrProtonSplash_0405();
  ~CrProtonSplash_0405();

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
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrProtonSplash_0506
{
public:
  CrProtonSplash_0506();
  ~CrProtonSplash_0506();

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
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrProtonSplash_0607
{
public:
  CrProtonSplash_0607();
  ~CrProtonSplash_0607();

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

class CrProtonSplash_0708
{
public:
  CrProtonSplash_0708();
  ~CrProtonSplash_0708();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization, spectral index and cut off for breakE<E
  double B_splash, b_splash, cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrProtonSplash_0809
{
public:
  CrProtonSplash_0809();
  ~CrProtonSplash_0809();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization, spectral index and cut off for breakE<E
  double B_splash, b_splash, cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};

class CrProtonSplash_0910
{
public:
  CrProtonSplash_0910();
  ~CrProtonSplash_0910();

private:
  // Normalization and spectral index for E<breakE
  double A_splash, a_splash;
  // Normalization, spectral index and cut off for breakE<E
  double B_splash, b_splash, cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(HepRandomEngine* engine);
  double upwardFlux();
};


#endif // CrProtonSubSplash_H
