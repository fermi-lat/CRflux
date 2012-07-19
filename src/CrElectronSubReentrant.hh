//
// CrElectronSubReentrant.hh
//

//$Header$

#ifndef CrElectronSubReentrant_H
#define CrElectronSubReentrant_H

#include <utility>
#include <string>

#include "CrSpectrum.hh"

// Forward declaration:
class CLHEP::HepRandomEngine;

class CrElectronReentrant_0001
{
public:
  CrElectronReentrant_0001();
  ~CrElectronReentrant_0001();

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

class CrElectronReentrant_0102
{
public:
  CrElectronReentrant_0102();
  ~CrElectronReentrant_0102();

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

class CrElectronReentrant_0203
{
public:
  CrElectronReentrant_0203();
  ~CrElectronReentrant_0203();

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

class CrElectronReentrant_0304
{
public:
  CrElectronReentrant_0304();
  ~CrElectronReentrant_0304();

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

class CrElectronReentrant_0405
{
public:
  CrElectronReentrant_0405();
  ~CrElectronReentrant_0405();

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

class CrElectronReentrant_0506
{
public:
  CrElectronReentrant_0506();
  ~CrElectronReentrant_0506();

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

class CrElectronReentrant_0611
{
public:
  CrElectronReentrant_0611();
  ~CrElectronReentrant_0611();

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


#endif // CrElectronSubReentrant_H
