/****************************************************************************
 * CrProtonPrimary.cc:
 ****************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrProton, the entry-point class for the cosmic-ray proton 
 * generation.
 ****************************************************************************
 * This program provides the primary component of cosmic ray protons.
 * Its angular distribution is assumed to be uniform for downward 
 * (theta = pi - zenith angle is 0 to pi/2) and zero for theta > pi/2.
 * The following itemizes other important features.
 * 1) One intrinsic spectrum is assumed common to all locations 
 *    on earth: a power-low with the geomagnetic cut-off at lower energy.  
 *    It is modified by the solar modulation, the geomagnetic cutoff, 
 *    and the east-west asymmetry effect (not yet implemented).
 * 2) Note that there are two other components of cosmic ray protons. 
 *    "reentrant" and "splash".  They are handled in other programs.
 * 3) A method CrProtonPrimary::energySrc returns an energy and 
 *    CrProtonPrimary::dir returns a direction in 
 *    cos(theta) and phi (azimuth angle). 
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
 * in 2000 Originally written by M. Ozaki (ISAS) 
 *         and K. Hirano (Hiroshima Univ.)
 * 2001-04 Spectrum part modified by T. Mizuno (Hiroshima Univ.)
 * 2001-05 Spectrum part modified by T. Mizuno (Hiroshima Univ.)
 * 2001-05 Comments added and unused codes removed by T. Kamae (SLAC)
 * 2001-05 Program checked by T. Kamae and H. Mizushima (Hiroshima)
 * 2001-11 Modified by T. Mizuno to construct a `stand-alone' module
 * 2003-02 Parameters are updated by T. Mizuno
 ****************************************************************************
 */

//$Header$


#include <cmath>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrProtonPrimary.hh"

typedef double G4double;

// private function definitions.
namespace {
  // rest energy (rest mass) of proton in units of GeV
  const G4double restE = 0.938;

  // The lower and higher (kinetic) energy limits of primary protons
  // generated in this program in units of GeV.  These values are set 
  // in constructor and when the satellite position is set
  G4double lowE_primary;
  G4double highE_primary;
  // The CutE_primary is the kinetic Energy corresponding to 
  // cutoff-rigidity. The value is set also in constructor and 
  // when the satellite position is set.
  G4double cutE_primary;


  // gives back v/c as a function of kinetic Energy
  inline G4double beta(G4double E /* GeV */){
    return sqrt(1 - pow(E/restE+1, -2));
  }

  // gives back the rigidity (p/Ze where p is the momentum, e means 
  // electron charge magnitude, and Z is the atomic number) in units of [GV],
  // as a function of kinetic Energy [GeV].
  inline G4double rigidity(G4double E /* GeV */){
    return sqrt(pow(E + restE, 2) - pow(restE, 2));
  }

  // gives back the kinetic energy [GeV] as a function of rigidity [GV]
  inline G4double energy(G4double rigidity /* GV */){
    return sqrt(pow(rigidity, 2) + pow(restE, 2)) - restE;
  }

  //============================================================
  /**
   *  Generate a random distribution of primary cosmic ray protons
   *  j(E) = mod_spec(E, phi) * geomag_cut(E, CutOff)
   *    mod_spec(E, phi) = org_spec(E+phi*1e-3) * 
   *     ((E+restE)**2 - restE**2)/((E+restE+phi*1e-3)**2-restE**2)
   *    org_spec(E) = A * rigidity(E)**-a
   *      A = 23.9 and a = 2.83
   *    rigidity(E) = sqrt((E+restE)**2 - restE**2)
   *    beta(E) = sqrt(1 - (E/restE+1)**-2)
   *    geomag_cut(E, CutOff) = 1/(1 + (rigidity(E)/CutOff)**-12.0)
   *      CutOff = 4.46 for Theta_M = 0.735 and altitude = 35km 
   *                                         (balloon experiment)
   *      phi = 540 and 1100 [MV] for Solar minimum and maximum, respectively
   *    E: [GeV]
   *    j: [c/s/m^2/sr/MeV]
   *
   *  References:
   *  org_spec: AMS data (Alcaraz et al. 2000, Phys. Letter B 472, 215)
   *  geomag_cut formula: 
   *    an eyeball fitting function to represent the AMS proton data.
   *  CutOff: calculated as (Rc/GV) = 14.9 * (1+h/R)^-2 * (cos(theta_M))^4,
   *          where h is the altitude from earth surface, 
   *                R is the mean radius of earth,
   *            and theta_M is geomagnetic lattitude. 
   *          References:
   *          "Handbook of space astronomy and astrophysics" 2nd edition, p225 
   *            (Zombeck, 1990, Cambridge University Press)
   *          "High Energy Astrophysics" 2nd edition, p325-330
   *            (M. S. Longair, 1992, Cambridge University Press)
   *  Solar modulation model (mod_spec) is from:
   *     Gleeson, L. J. and Axford, W. I. 1968, ApJ, 154, 1011-1026 (Eq. 11)
   */

  const G4double A_primary = 23.9; // normalization of incident spectrum
  const G4double a_primary = 2.83; // differential spectral index

  // Gives back the geomagnetic cutoff factor to the intrinsic 
  // primary cosmic ray spectrum for a kinetic energy E(GeV) 
  // and a geomagnetic lattitude theta_M(rad)
  inline G4double geomag_cut(G4double E, G4double cor /* GV */){
    return 1./(1 + pow(rigidity(E)/cor, -12.0));
  }

  // The unmodulated primary proton spectrum outside the Solar system is 
  // returned.
  inline G4double org_spec(G4double E /* GeV */)
  {
    return A_primary * pow(rigidity(E), -a_primary);
  }

  // The modulated proton flux for a "phi" value is returned. 
  // Force-field approximation of the Solar modulation is used.
  inline G4double mod_spec(G4double E /* GeV */, G4double phi /* MV */){
    return org_spec(E + phi*1e-3) * (pow(E+restE, 2) - pow(restE, 2))
      / (pow(E+restE+phi*1e-3,2) - pow(restE,2));
  }

  // The final spectrum for the primary proton.
  inline G4double primaryCRspec
  (G4double E /* GeV */, G4double cor /* GV*/, G4double phi /* MV */){
    return mod_spec(E, phi) * geomag_cut(E, cor);
  }

  // To speed up generation of the spectrum below the geomagnetic cut 
  // off, random numbers are generated to "an envelope function" for
  // whose integral the inverse function can be found.  The final one 
  // is obtained by throwing away some random picks of energy.

  // Envelope function in the lower energy range 
  // (E<Ec, where Ec corresponds to cutoff rigidity).
  // Primary spectrum of cosmic-ray proton is enveloped by a linear function
  // between lowE_primary and cutE_primary.
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

  // The envelope function in the higher energy range
  // (E>Ec, where Ec corresponds to cutoff rigidity).
  // Primary spectrum of cosmic-ray proton is enveloped by power-law function
  // between cutE_primary and highE_primary
  inline G4double primaryCRenvelope2
  (G4double E /* GeV */, G4double /* cor */ /* MV */, 
   G4double /* phi */ /* MV */){
    return A_primary * pow(E, -a_primary);
  }

  // The integral of the envelope function in the higher energy range
  inline G4double primaryCRenvelope2_integral
  (G4double E /* GeV */, G4double /* cor */ /* MV */, 
   G4double /* phi */ /* MV */){
    return A_primary/(-a_primary+1) * pow(E, -a_primary+1);
  }

  // The inverse function of the integral of the envelope function
  // in the higher energy range.
  // This function returns energy obeying envelope function.
  inline G4double primaryCRenvelope2_integral_inv
  (G4double value, G4double /* cor */ /* MV */, G4double /* phi */ /* MV */){
    return pow((-a_primary+1)/ A_primary * value , 1./(-a_primary+1));
  }

  // The random number generator for the primary component
  G4double primaryCRenergy(CLHEP::HepRandomEngine* engine, 
			   G4double cor, G4double solarPotential){
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
    while(1){
      if (engine->flat() <= envelope1_area/(envelope1_area + envelope2_area)){
        // Use the envelop function in the lower energy range
	// (E<Ec where Ec corresponds to the cutoff rigidity).
	// We enveloped proton spectrum by linear function between
	// lowE_primary and cutE_primary, and assume that the flux
	// at lowE_primary is 0.
	G4double E1, E2;
	E1 = engine->flat() * (cutE_primary-lowE_primary) + lowE_primary;
	E2 = engine->flat() * (cutE_primary-lowE_primary) + lowE_primary;
	if (E1>E2){E=E1;} else {E=E2;}
        if (engine->flat() <= 
	    primaryCRspec(E, cor, solarPotential) 
	    / primaryCRenvelope1(E, cor, solarPotential))
          break;
      }
      else{
        // Use the envelop function in the higher energy range
        // (E>Ec where Ec corresponds to the cutoff rigidity).
        r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
        E = primaryCRenvelope2_integral_inv(r, cor, solarPotential);
        if (engine->flat() <= primaryCRspec(E, cor, solarPotential) 
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
    {3819, 3086, 2552, 2150, 1838, 1592, 1393}, // COR = 0.5GV
    {3080, 2578, 2191, 1885, 1640, 1440, 1275}, // COR = 1 GV
    {1746, 1549, 1383, 1242, 1120, 1016, 924.9}, // COR = 2 GV
    {1072, 981.1, 901.0, 830.2, 767.1, 710.7, 660.1}, // COR = 3 GV
    {720.4, 671.8, 627.8, 587.9, 551.5, 518.2, 487.7}, // COR = 4 GV
    {517.8, 488.9, 462.1, 437.5, 414.6, 393.5, 373.8}, // COR = 5 GV
    {390.9, 372.2, 354.8, 338.5, 323.2, 308.9, 295.5}, // COR = 6 GV
    {306.2, 293.4, 281.3, 269.9, 259.3, 249.2, 239.6}, // COR = 7 GV
    {246.7, 237.5, 228.9, 220.6, 212.8, 205.4, 198.3}, // COR = 8 GV
    {203.3, 196.5, 190.1, 183.9, 178.0, 172.4, 167.0}, // COR = 9 GV
    {170.6, 165.5, 160.5, 155.8, 151.3, 146.9, 142.7}, // COR = 10 GV
    {145.4, 141.4, 137.5, 133.8, 130.2, 126.7, 123.4}, // COR = 11 GV
    {125.5, 122.3, 119.2, 116.2, 113.7, 110.6, 107.9}, // COR = 12 GV
    {109.5, 106.9, 104.4, 102.0, 99.7, 97.4, 95.2}, // COR = 13 GV
    {96.5, 94.3, 92.3, 90.3, 88.3, 86.5, 84.6}, // COR = 14 GV
    {85.7, 83.9, 82.2, 80.5, 78.9, 77.3, 75.8} // COR = 15 GV
  };
  


  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrProtonPrimary::CrProtonPrimary():CrSpectrum()
{
  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary proton can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 1000.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);
}


CrProtonPrimary::~CrProtonPrimary()
{
  ;
}


// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrProtonPrimary::setPosition(double latitude, double longitude){
  CrSpectrum::setPosition(latitude, longitude);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary proton can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 1000.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrProtonPrimary::setPosition(double latitude, double longitude, double time){
  CrSpectrum::setPosition(latitude, longitude, time);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary proton can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 1000.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrProtonPrimary::
setPosition(double latitude, double longitude, double time, double altitude){
  CrSpectrum::setPosition(latitude, longitude, time, altitude);
  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary proton can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 1000.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set geomagnetic cutoff rigidity and calculate the energies related.
// These energies are used to generate the particle. 
void CrProtonPrimary::setCutOffRigidity(double cor){
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
std::pair<double,double> CrProtonPrimary::dir(double energy, 
					      CLHEP::HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward direction has plus sign in cos(theta),
  // and phi = 0 for the particle comming from east
  // (smallest flux for positively charged particle)
  // and phi=pi/2 for that comming from north
{

  // CrSpectrum class takes care of direction generation
  double rig = rigidity(energy*0.001);
  double coeff = -12.0;
  double polarity = 1.0; // positively charged particle
  return CrSpectrum::EW_dir(rig, coeff, polarity, engine);
}


// Gives back particle energy
double CrProtonPrimary::energySrc(CLHEP::HepRandomEngine* engine) const
{
  return primaryCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".
double CrProtonPrimary::flux() const
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
double CrProtonPrimary::solidAngle() const
{
   // * 1.4 since Cos(theta) ranges from 1 to -0.4
  return  2 * M_PI * 1.4;
}


// Gives back particle name
const char* CrProtonPrimary::particleName() const
{
  return "proton";
}


// Gives back the name of the component
std::string CrProtonPrimary::title() const
{
  return  "CrProtonPrimary";
}


