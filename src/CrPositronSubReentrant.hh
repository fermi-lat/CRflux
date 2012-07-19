//
// CrPositronSubReentrant.hh
//

//$Header$

#ifndef CrPositronSubReentrant_H
#define CrPositronSubReentrant_H

#include <utility>
#include <string>

#include "CrSpectrum.hh"

// Forward declaration:
class CLHEP::HepRandomEngine;

class CrPositronReentrant_0001
{
public:
  CrPositronReentrant_0001();
  ~CrPositronReentrant_0001();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_reent, b_reent;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_reent, c_reent;
  // Normalization and spectral index for highE_break<E
  double D_reent, d_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0102
{
public:
  CrPositronReentrant_0102();
  ~CrPositronReentrant_0102();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_reent, b_reent;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_reent, c_reent;
  // Normalization and spectral index for highE_break<E
  double D_reent, d_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0203
{
public:
  CrPositronReentrant_0203();
  ~CrPositronReentrant_0203();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_reent, b_reent;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_reent, c_reent;
  // Normalization and spectral index for highE_break<E
  double D_reent, d_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0304
{
public:
  CrPositronReentrant_0304();
  ~CrPositronReentrant_0304();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_reent, b_reent;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_reent, c_reent;
  // Normalization and spectral index for highE_break<E
  double D_reent, d_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0405
{
public:
  CrPositronReentrant_0405();
  ~CrPositronReentrant_0405();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_reent, b_reent;
  // Normalization and spectral index for highE_break<E
  double C_reent, c_reent;
  double lowE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0506
{
public:
  CrPositronReentrant_0506();
  ~CrPositronReentrant_0506();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_reent, b_reent;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_reent, c_reent;
  // Normalization and spectral index for highE_break<E
  double D_reent, d_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0611
{
public:
  CrPositronReentrant_0611();
  ~CrPositronReentrant_0611();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<midE_break
  double B_reent, b_reent;
  // Normalization and spectral index for midE_break<E<highE_break
  double C_reent, c_reent;
  // Normalization and spectral index for highE_break<E
  double D_reent, d_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double midE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};


#endif // CrPositronSubReentrant_H
