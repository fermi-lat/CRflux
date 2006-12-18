//
// CrProtonSubReentrant.hh
//

//$Header$

#ifndef CrProtonSubReentrant_H
#define CrProtonSubReentrant_H

#include <utility>
#include <string>

#include "CrSpectrum.hh"

// Forward declaration:
class CLHEP::HepRandomEngine;

class CrProtonReentrant_0002
{
public:
  CrProtonReentrant_0002();
  ~CrProtonReentrant_0002();

private:
  // Normalization and spectral index for E<breakE
  double A_reent, a_reent;
  // Normalization, spectral index and cutoff for E>breakE
  double B_reent, b_reent, cutOff;
  // angular distribution constant
  double ang_reent;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double theta(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrProtonReentrant_0203
{
public:
  CrProtonReentrant_0203();
  ~CrProtonReentrant_0203();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_reent, b_reent;
  // Normalization and spectral index for highE_break<E
  double C_reent, c_reent;
  // angular distribution constant
  double ang_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double theta(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrProtonReentrant_0304
{
public:
  CrProtonReentrant_0304();
  ~CrProtonReentrant_0304();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_reent, b_reent;
  // Normalization and spectral index for highE_break<E
  double C_reent, c_reent;
  // angular distribution constant
  double ang_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double theta(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrProtonReentrant_0405
{
public:
  CrProtonReentrant_0405();
  ~CrProtonReentrant_0405();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_reent, b_reent;
  // Normalization and spectral index for highE_break<E
  double C_reent, c_reent;
  // angular distribution constant
  double ang_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double theta(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrProtonReentrant_0506
{
public:
  CrProtonReentrant_0506();
  ~CrProtonReentrant_0506();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_reent, b_reent;
  // Normalization and spectral index for highE_break<E
  double C_reent, c_reent;
  // angular distribution constant
  double ang_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double theta(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrProtonReentrant_0607
{
public:
  CrProtonReentrant_0607();
  ~CrProtonReentrant_0607();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_reent, b_reent;
  // Normalization and spectral index for highE_break<E
  double C_reent, c_reent;
  // angular distribution constant
  double ang_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double theta(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrProtonReentrant_0708
{
public:
  CrProtonReentrant_0708();
  ~CrProtonReentrant_0708();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_reent, b_reent;
  // Normalization and spectral index for highE_break<E
  double C_reent, c_reent;
  // angular distribution constant
  double ang_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double theta(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrProtonReentrant_0809
{
public:
  CrProtonReentrant_0809();
  ~CrProtonReentrant_0809();

private:
  // Normalization and spectral index for E<breakE
  double A_reent, a_reent;
  // Normalization, spectral index and cutoff for E>breakE
  double B_reent, b_reent, cutOff;
  // angular distribution constant
  double ang_reent;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double theta(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrProtonReentrant_0910
{
public:
  CrProtonReentrant_0910();
  ~CrProtonReentrant_0910();

private:
  // Normalization and spectral index for E<breakE
  double A_reent, a_reent;
  // Normalization, spectral index and cutoff for E>breakE
  double B_reent, b_reent, cutOff;
  // angular distribution constant
  double ang_reent;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double theta(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};


#endif // CrProtonSubReentrant_H
