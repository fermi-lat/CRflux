//
// CrPositronSubReentrant.hh
//

//$Header$

#ifndef CrPositronSubReentrant_H
#define CrPositronSubReentrant_H

#include <utility>
#include <string>

// Forward declaration:
class CLHEP::HepRandomEngine;

class CrPositronReentrant_0003
{
public:
  CrPositronReentrant_0003();
  ~CrPositronReentrant_0003();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<middleE_break
  double B_reent, b_reent;
  // Normalization and spectral index for middleE_break<E<highE_break
  double C_reent, c_reent;
  // Normalization and spectral index for highE_break<E
  double D_reent, d_reent;
  // The energies where spectrum breaks
  double lowE_break;
  double middleE_break;
  double highE_break;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0306
{
public:
  CrPositronReentrant_0306();
  ~CrPositronReentrant_0306();

private:
  // Normalization and spectral index for E<breakE
  double A_reent, a_reent;
  // Normalization, spectral index and cut off for breakE<E
  double B_reent, b_reent;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0608
{
public:
  CrPositronReentrant_0608();
  ~CrPositronReentrant_0608();

private:
  // Normalization and spectral index for E<breakE
  double A_reent, a_reent;
  // Normalization, spectral index and cut off for breakE<E
  double B1_reent, b1_reent;
  double B2_reent, b2_reent;
  double cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0809
{
public:
  CrPositronReentrant_0809();
  ~CrPositronReentrant_0809();

private:
  // Normalization and spectral index for E<breakE
  double A_reent, a_reent;
  // Normalization, spectral index and cut off for breakE<E
  double B1_reent, b1_reent;
  double B2_reent, b2_reent;
  double cutOff;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_0910
{
public:
  CrPositronReentrant_0910();
  ~CrPositronReentrant_0910();

private:
  // Normalization and spectral index for E<breakE
  double A_reent, a_reent;
  // Normalization, spectral index and cut off for breakE<E
  double B_reent, b_reent;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};

class CrPositronReentrant_1011
{
public:
  CrPositronReentrant_1011();
  ~CrPositronReentrant_1011();

private:
  // Normalization and spectral index for E<breakE
  double A_reent, a_reent;
  // Normalization and spectral index for breakE<E
  double B_reent, b_reent;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(CLHEP::HepRandomEngine* engine);
  double downwardFlux();
};


#endif // CrPositronSubReentrant_H
