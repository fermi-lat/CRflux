//
// CrElectronSubReentrant.hh
//

//$Header$

#ifndef CrElectronSubReentrant_H
#define CrElectronSubReentrant_H

#include <utility>
#include <string>

// Forward declaration:
class HepRandomEngine;

class CrElectronReentrant_0003
{
public:
  CrElectronReentrant_0003();
  ~CrElectronReentrant_0003();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_reent, b_reent;
  // Normalization, spectral index and cutoff for E>highE_break
  double C1_reent, c1_reent;
  double C2_reent, c2_reent;
  double cutOff;
  // The energies where spectrum breaks
  double lowE_break, highE_break;

public:
  double energy(HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronReentrant_0306
{
public:
  CrElectronReentrant_0306();
  ~CrElectronReentrant_0306();

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
  double energy(HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronReentrant_0608
{
public:
  CrElectronReentrant_0608();
  ~CrElectronReentrant_0608();

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
  double energy(HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronReentrant_0809
{
public:
  CrElectronReentrant_0809();
  ~CrElectronReentrant_0809();

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
  double energy(HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronReentrant_0910
{
public:
  CrElectronReentrant_0910();
  ~CrElectronReentrant_0910();

private:
  // Normalization and spectral index for E<lowE_break
  double A_reent, a_reent;
  // Normalization and spectral index for lowE_break<E<highE_break
  double B_reent, b_reent;
  // Normalization, spectral index and cut off for highE_break<E
  double C_reent, c_reent, cutOff;
  // The energies where spectrum breaks
  double lowE_break, highE_break;

public:
  double energy(HepRandomEngine* engine);
  double downwardFlux();
};

class CrElectronReentrant_1011
{
public:
  CrElectronReentrant_1011();
  ~CrElectronReentrant_1011();

private:
  // Normalization and spectral index for E<breakE
  double A_reent, a_reent;
  // Normalization and spectral index for breakE<E
  double B_reent, b_reent;
  // The energy where spectrum breaks
  double breakE;

public:
  double energy(HepRandomEngine* engine);
  double downwardFlux();
};


#endif // CrElectronSubReentrant_H
