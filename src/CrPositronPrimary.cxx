/****************************************************************************
 * CrPositronPrimary.cc:
 ****************************************************************************
 * Read comments at the head of CosmicRay PrimaryGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayPrimaryGeneratorAction.cc 
 * via CrPositron, the entry-point class for the cosmic-ray positron 
 * generation.
 ****************************************************************************
 * This program provides the primary component of cosmic ray positrons.
 * Its angular distribution is assumed to be uniform for downward 
 * (theta = pi - zenith angle is 0 to pi/2) and zero for theta > pi/2.
 * The following itemizes other important features.
 * 1) One intrinsic spectrum is assumed common to all locations 
 *    on earth: a power-low with the geomagnetic cut-off at lower energy.  
 *    It is modified by the solar modulation, the geomagnetic cutoff, 
 *    and the east-west asymmetry effect (not yet implemented).
 * 2) Note that there are two other components of cosmic ray positrons. 
 *    "reentrant" and "splash".  They are handled in other programs.
 * 3) A method CrPositronPrimary::energySrc returns an energy and 
 *    CrPositronPrimary::dir returns a direction 
 *    in cos(theta) and phi (azimuth angle). 
 ****************************************************************************
 * Definitions:
 * 1) The z-axis points upward (from Calorimeter to Tracker).
 * 2) Partile of theta=0 points for downward (i.e., comes from zenith)
 *    and that of theta=pi points for upward (comes from nadir).
 * 3) Partile of phi=0 comes along x-axis (from x>0 to x=0) and
 *    that of phi=pi/2 comes along y-axis (from y>0 to y=0).
 * 4) Energy means kinetic energy unless defined otherwise. 
 * 5) Magnetic latitude theta_M is in radian. 
 * 6) Particle direction is defined by cos(theta) and phi (in radian).
 ****************************************************************************
 * 2001-07 Written by Y. Fukazawa (Hiroshima Univ).
 * 2001-10 Modified by T. Mizuno and Y. Fukazawa
 * 2001-11 Modified by T. Mizuno to construct a `stand-alone' module
 ****************************************************************************
 */

// $Header$


#include <cmath>

// CLHEP
//#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrPositronPrimary.hh"


typedef  double G4double;


// private function definitions.
namespace {
  // The rest energy (rest mass) of positron in [GeV]
  const G4double restE = 5.11e-4; // rest energy of positron in [GeV]
  // The lower and higher (kinetic) energy limits of primary positrons 
  // generated in this program in units of GeV.  These values are set 
  // in constructor and when the satellite position is set
  G4double lowE_primary;
  G4double highE_primary;
  // The CutE_primary is the kinetic Energy corresponding to 
  // cutoff-rigidity.   The value is set in constructor.
  G4double cutE_primary;

  // Gives back v/c as a function of kinetic Energy
  inline G4double beta(G4double /* E */  /* in GeV */)
  {
#if 0	// if E ~ restE
    return sqrt(1 - pow(E/restE+1, -2));
#else	// if E >> restE
    return 1.0;
#endif
  }


  // Gives back the rigidity (p/e where p is the momentum and e means
  // positron charge magnitude) in units of [GV]
  // as a function of the kinetic energy [GeV].
  inline G4double rigidity(G4double E /* GeV */)
  {
#if 0	// if E ~ restE
    return sqrt(pow(E + restE, 2) - pow(restE, 2));
#else	// if E >> restE
    return E;
#endif
  }


  // Gives back the kinetic energy [GeV] for a rigidity [GV].
  inline G4double energy(G4double rigidity /* GV */)
  {
#if 0	// if energy ~ restE
    return sqrt(pow(rigidity, 2) + pow(restE, 2)) - restE;
#else	// if energy >> restE
    return rigidity;
#endif
  }



  //=====================================================================
  /**
   * Generate a random distribution of primary cosmic ray positrons:
   * j(E) = mod_spec(E, phi) * geomag_cut(E, CutOff)
   *   mod_spec(E, phi) = org_spec(E + phi * 1e-3) * 
   *    ((E+restE)**2 - restE**2) / ((E+restE+phi*1e-3)**2 - restE**2)
   *   org_spec(E) = A * rigidity(E)**-a
   *     A = 0.0564
   *     a = 3.33
   *   rigidity(E) = sqrt((E+restE)**2 - restE**2)
   *   beta(E) = sqrt(1 - (E/restE+1)**-2)
   *   geomag_cut(E, CutOff) = 1/(1 + (rigidity(E)/CutOff)**-6.0)
   *     CutOff = 4.46 for Theta_M = 0.735 and altitude = 35km 
   *                                           (balloon experiment)
   *     phi = 540, 1100 [MV] for Solar minimum, maximum
   *   E: [GeV]
   *   j: [c/s/m^2/sr/MeV]
   *
   * References:
   *  org_spec: Komori, Y. et al. 1999, Proc. of Dai-Kikyu (large balloon) 
   *                                    Sympo. Heisei-11yr, 33-36  (Fig. 2)
   *            R.L.Golden et al. 1996, ApJL 457, 103
   *            Webber 1983, Composition and Origin of Cosmic Rays
   *              (ed. M. M. Shapiro)
   *            Longair, M. S. 1992, High Energy Astrophysics, vol 1
   *              (2nd edition, Cambridge University Press)
   *  mod_spec: Gleeson, L. J. and Axford, W. I. 
   *            1968, ApJ, 154, 1011-1026 (Eq. 11)
   *  geomag_cut formula: an eyeball fitting function to represent 
   *                      the AMS proton data
   *  CutOff: calculated as (Rc/GV) = 14.9 * (1+h/R)^-2 * (cos(theta_M))^4,
   *          where h gives altitude from earth surface, 
   *                R means an radius of earth,
   *            and theta_M is geomagnetic lattitude. See references below.
   *          "Handbook of space astronomy and astrophysics" 2nd edition, p225 
   *            (Zombeck, 1990, Cambridge University Press)
   *          "High Energy Astrophysics" 2nd edition, p325-330
   *            (M. S. Longair, 1992, Cambridge University Press)
   *
   */

  // Normalization factor of the original spec.
  const G4double A_primary = 0.05; // 0.07*0.078
  // Differential spectral index
  const G4double a_primary = 3.3;


  // Gives back the geomagnetic cutoff factor to the intrinsic 
  // primary cosmic ray spectrum for a kinetic energy E (GeV) 
  // and a cut-off rigidity value (GV).
  inline G4double geomag_cut(G4double E /* GeV */, G4double cor /* GV */)
  {
    return 1./(1 + pow(rigidity(E)/cor, -6.0));
  }


  // The unmodulated primary positron spectrum outside the Solar system is 
  // returned.
  inline G4double org_spec(G4double E /* GeV */)
  {
    return A_primary * pow(rigidity(E), -a_primary);
  }


  // The modulated positron flux for a "phi" value is calculated 
  // according to the force-field approximation of the solar modulation.
  inline G4double mod_spec(G4double E /* GeV */, G4double phi /* MV */)
  {
    return org_spec(E + phi * 1e-3) * (pow(E+restE, 2) - pow(restE, 2))
      / (pow(E+restE+phi*1e-3, 2) - pow(restE, 2));
  }


  // The final spectrum for the primary positron.
  inline G4double primaryCRspec
  (G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */)
  {
    return mod_spec(E, phi) * geomag_cut(E, cor);
  }


  // To generate the spectrum (primaryCRspec) efficiently (ie. minimum 
  // call of random numbers), the spectrum is devided into 2 parts:
  // between lowE_primary and cutE_primary and between cutE_primary
  // and highE_primary. For both portion the spectrum has a complicated
  // formula and the inverse function of its integral doesn't come
  // easily. Instead, we use a simpler function that envelopes the 
  // true spectrum and generate a random number that obeys this simple
  // function. We then removed the excess by comparing with the true spectrum.
  // A linear function is used to envelope the lower part of the 
  // spectrum, and power-law function is used to envelope the
  // higher part of the spectrum.

  // Envelope function in the lower energy range
  // (E<Ec, where Ec corresponds to cutoff rigidity) is given.
  // Primary spectrum of cosmic-ray positron is enveloped by a linear function
  // between lowE_primary and cutE_primary
  inline G4double primaryCRenvelope1
  (G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */){
    G4double coeff = 
      ( primaryCRspec(cutE_primary,cor,phi)
	- primaryCRspec(lowE_primary,cor,phi) ) 
      / (cutE_primary-lowE_primary);
    return coeff * (E-lowE_primary) + primaryCRspec(lowE_primary,cor,phi);
  }


  // Integral of the envelope function in the lower energy range
  inline G4double primaryCRenvelope1_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double coeff =
      ( primaryCRspec(cutE_primary,cor,phi) 
	- primaryCRspec(lowE_primary,cor,phi) ) 
      / (cutE_primary-lowE_primary);
    return 0.5 * coeff * pow(E-lowE_primary,2) + 
      primaryCRspec(lowE_primary,cor,phi) * (E-lowE_primary);
  }


  // Envelope function in higher energy
  // (E>Ec, where Ec corresponds to cutoff rigidity) is given.
  // Primary spectrum of cosmic-ray positron is enveloped by power-law function
  // between cutE_primary and highE_primary
  inline G4double primaryCRenvelope2
  (G4double E /* GeV */, G4double /* cor */  /* in GV */, 
   G4double /* phi */  /* MV */)
  {
    return A_primary * pow(E, -a_primary);
  }


  // Integrated envelope function in the higher energy range.
  inline G4double primaryCRenvelope2_integral
  (G4double E  /* in GeV */, 
   G4double /* cor */  /* GV */, G4double /* phi */ /* MV */)
  {
    return A_primary * pow(E, -a_primary+1) / (-a_primary+1);
  }


  // Inverse function of the integrated envelope function
  // in the higher energy range.
  inline G4double primaryCRenvelope2_integral_rev
  (G4double value, G4double /* cor */ /* GV */, G4double /* phi */ /* MV */)
  {
    return pow(value * (-a_primary+1) / A_primary, 1./(-a_primary+1));
  }


  // The rundam number generator for the primary component.
  G4double primaryCRenergy(CLHEP::HepRandomEngine* engine, 
			   G4double cor, G4double solarPotential)
  {
    G4double rand_min_1 = 
      primaryCRenvelope1_integral(lowE_primary, cor, solarPotential);
    G4double rand_max_1 = 
      primaryCRenvelope1_integral(cutE_primary, cor, solarPotential);
    G4double rand_min_2 = 
      primaryCRenvelope2_integral(cutE_primary, cor, solarPotential);
    G4double rand_max_2 = 
      primaryCRenvelope2_integral(highE_primary, cor, solarPotential);

    G4double envelope1_area = rand_max_1 - rand_min_1;
    G4double envelope2_area = rand_max_2 - rand_min_2;

    double r, E; // E means energy in GeV

    while (1){
      if (engine->flat() <= 
	  envelope1_area / (envelope1_area + envelope2_area)){
        // The envelope function for the lower energy part:
        // We enveloped the spectrum by a linear function between
        // lowE_primary and cutE_primary, and assume the flux to be
        // zero below lowE_primary.
	G4double E1, E2;
        E1 = engine->flat() * (cutE_primary-lowE_primary) + lowE_primary;
        E2 = engine->flat() * (cutE_primary-lowE_primary) + lowE_primary;
        if (E1>E2){E=E1;} else {E=E2;}
        if (engine->flat() <= 
	    primaryCRspec(E, cor, solarPotential) 
	    / primaryCRenvelope1(E, cor, solarPotential))
          break;
      } else {
        // Envelope in the higher energy range.
        r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
        E = primaryCRenvelope2_integral_rev(r, cor, solarPotential);
        if (engine->flat() <= 
	    primaryCRspec(E, cor, solarPotential) 
	    / primaryCRenvelope2(E, cor, solarPotential))
          break;
      }
    }
    return E;
  }

  // This array stores vertically downward flux in unit of [c/s/m^s/sr]
  // as a function of COR and phi (integral_array[COR][phi]).
  // The flux is integrated between lowE_primary and highE_primary.
  // COR = 0.5, 1, 2, ..., 15 [GV]
  // phi = 500, 600, ..., 1100 [MV]
  G4double integral_array[16][7] = {
    {9.404, 6.871, 5.194, 4.036, 3.209, 2.6, 2.142}, // COR = 0.5GV
    {5.212, 4.125, 3.328, 2.728, 2.269, 1.91, 1.625}, // COR = 1 GV
    {2.126, 1.82, 1.572, 1.368, 1.2, 1.058, 0.939}, // COR = 2 GV
    {1.109, 0.987, 0.883, 0.794, 0.716, 0.649, 0.59}, // COR = 3 GV
    {0.668, 0.608, 0.556, 0.509, 0.468, 0.432, 0.399}, // COR = 4 GV
    {0.441, 0.407, 0.378, 0.351, 0.327, 0.305, 0.285}, // COR = 5 GV
    {0.31, 0.29, 0.271, 0.255, 0.239, 0.225, 0.212}, // COR = 6 GV
    {0.228, 0.215, 0.203, 0.192, 0.182, 0.172, 0.164}, // COR = 7 GV
    {0.174, 0.165, 0.157, 0.149, 0.142, 0.136, 0.129}, // COR = 8 GV
    {0.137, 0.131, 0.125, 0.119, 0.114, 0.109, 0.105}, // COR = 9 GV
    {0.11, 0.105, 0.101, 0.097, 0.093, 0.089, 0.086}, // COR = 10 GV
    {0.09, 0.087, 0.083, 0.08, 0.077, 0.075, 0.072}, // COR = 11 GV
    {0.075, 0.072, 0.07, 0.067, 0.065, 0.063, 0.061}, // COR = 12 GV
    {0.063, 0.061, 0.059, 0.057, 0.055, 0.054, 0.052}, // COR = 13 GV
    {0.054, 0.052, 0.051, 0.049, 0.048, 0.046, 0.045}, // COR = 14 GV
    {0.047, 0.045, 0.044, 0.043, 0.042, 0.04, 0.039} // COR = 15 GV
  };

  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrPositronPrimary::CrPositronPrimary():CrSpectrum()
{
  // Set lower and higher energy limits of primary positron.
  // At lowE_primary, flux of primary positron can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 1000.0;
  // "cutE_primary" is the kinetic energy corresponding to 
  // cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);
}


CrPositronPrimary::~CrPositronPrimary()
{
  ;
}

// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrPositronPrimary::
setPosition(double latitude, double longitude){
  CrSpectrum::setPosition(latitude, longitude);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary positron can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 1000.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrPositronPrimary::
setPosition(double latitude, double longitude, double time){
  CrSpectrum::setPosition(latitude, longitude, time);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary positron can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 1000.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrPositronPrimary::
setPosition(double latitude, double longitude, double time, double altitude){
  CrSpectrum::setPosition(latitude, longitude, time, altitude);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary positron can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 1000.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set geomagnetic cutoff rigidity and calculate the energies related.
// These energies are used to generate the particle. 
void CrPositronPrimary::setCutOffRigidity(double cor){
  CrSpectrum::setCutOffRigidity(cor);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary proton can be
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 1000.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Gives back particle direction in (cos(theta), phi)
std::pair<double,double> CrPositronPrimary::dir(double energy, 
						CLHEP::HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward direction has plus sign in cos(theta),
  // and phi = 0 for the particle comming from east
  // (smallest flux for positively charged particle)
  // and phi=pi/2 for that comming from north
{

  // CrSpectrum class takes care of direction generation
  double rig = rigidity(energy*0.001);
  double coeff = -6.0;
  double polarity = 1.0; // positively charged particle
  return CrSpectrum::EW_dir(rig, coeff, polarity, engine);

}



// Gives back particle energy
double CrPositronPrimary::energySrc(CLHEP::HepRandomEngine* engine) const
{
  return primaryCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".
double CrPositronPrimary::flux() const
{
  // Straight downward (theta=0) flux integrated over energy,
  // given by integral_array[16][7]
  G4double energy_integral;

  // cutoff rigidity and solar potential to calculate the energy integral.
  G4double cor, phi; 

  // cor is restricted in 0.5 < cor < 14.9[GV] and 
  // phi is restricted in 500 < phi < 1100[MV]
  cor = m_cutOffRigidity;
  if (m_cutOffRigidity < 0.5){ cor = 0.5; }
  if (m_cutOffRigidity > 14.9){ cor = 14.9; }
  phi = m_solarWindPotential;
  if (m_solarWindPotential < 500.0){ phi = 500.0; }
  if (m_solarWindPotential > 1100.0){ phi = 1100.0; }

  if (cor >=1.0){
    phi = phi/100.0 - 5; // 500 MV corresponds to 0, 600 MV corresponds to 1, etc..
    G4double tmp1 = 
      integral_array[int(cor)][int(phi)] +
      (cor - int(cor)) * (integral_array[int(cor)+1][int(phi)]-integral_array[int(cor)][int(phi)]);
    G4double tmp2 = 
      integral_array[int(cor)][int(phi)+1] +
      (cor - int(cor)) * (integral_array[int(cor)+1][int(phi)+1]-integral_array[int(cor)][int(phi)+1]);
    energy_integral = tmp1 + (tmp2-tmp1)*(phi-int(phi));
  }

  if (cor < 1.0){
    phi = phi/100.0 - 5; // 500 MV corresponds to 0, 600 MV corresponds to 1, etc..
    G4double tmp1 = 
      integral_array[0][int(phi)] +
      (2*cor - 1) *(integral_array[1][int(phi)]-integral_array[0][int(phi)]);
    G4double tmp2 = 
      integral_array[0][int(phi)+1] +
      (2*cor - 1) *(integral_array[1][int(phi)+1]-integral_array[0][int(phi)+1]);
    energy_integral = tmp1 + (tmp2-tmp1)*(phi-int(phi));
  }

  // We assume that the flux is uniform above the earth horizon.
  // Then the average flux is equal to the vertically downward one.
  return  energy_integral;  // [c/s/m^2/sr]
}


// Gives back solid angle from which particle comes
double CrPositronPrimary::solidAngle() const
{
  // * 1.4 since Cos(theta) ranges from 1 to -0.4
  return  2 * M_PI * 1.4;
}


// Gives back particle name
const char* CrPositronPrimary::particleName() const
{
  return "e+";
}


// Gives back the name of the component
std::string CrPositronPrimary::title() const
{
  return  "CrPositronPrimary";
}

