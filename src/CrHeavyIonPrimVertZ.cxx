/****************************************************************************
 * CrHeavyIonPrimVertZ.cc:  
 ********************************
 * COPIED from CrHeavyIonPrimary.cc, only defining theta=0: C.LAVALLEY
 ****************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrExample, the entry-point class for the cosmic-ray ion generation.
 ****************************************************************************
 * This program provides the primary component of cosmic ray ions.
 * Its angular distribution is assumed to be uniform for downward 
 * (theta = pi - zenith angle is 0 to pi/2) and zero for theta > pi/2.
 * The following itemizes other important features.
 * 1) The ion atomic number is selected according to the relative abundances 
 *    reported by J.J Engelmann et al., A&A 233,96 (1990). 
 *    The mass number is taken as that of the most abundant isotope.
 * 2) One intrinsic spectrum is assumed common to all locations 
 *    on earth: a power-low with the geomagnetic cut-off at lower energy.  
 *    It is modified by the solar modulation, the geomagnetic cutoff, 
 *    and the east-west asymmetry effect (not yet implemented).
 * 3) A method CrHeavyIonPrimary::energySrc returns an energy and 
 *    CrHeavyIonPrimary::dir returns a direction in 
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
 * 2004-11 Adapted by B. Lott from CrAlphaPrimary.cxx written by  Y. Fukazawa
 *         and T. Mizuno
 ****************************************************************************
 */

//$Header: 

#include <math.h>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/Random.h>

#include "CrHeavyIonPrimVertZ.hh"

typedef double G4double;

// private function definitions.
namespace { 

  CLHEP::HepRandomEngine* mm_engine;
// Set A and Z of ion
  double A_ion,z_ion;
  double get_a_ion();
  //double get_z_ion();

  //const G4double z_ion, A_ion; 
  // rest energy of  ion  in units of GeV
  const G4double restE = 0.931*A_ion;

  // The lower and higher (kinetic) energy limits of primary ions
  // generated in this program in units of GeV.  These values are set 
  // in constructor and when the satellite position is set
  G4double lowE_primary;
  G4double highE_primary;
  // The CutE_primary is the kinetic Energy corresponding to 
  // cutoff-rigidity. The value is set also in constructor and 
  // when the satellite position is set.
  G4double cutE_primary;
  // atomic number of  ion 
 /** inline G4double get_z_ion(){   
 
   //CL: to fix a Z: 
   int iz=11;
   std::cout << "IN CrHEavyIonPrimary: selected z = " << iz+3 << std::endl;
   return (G4double) iz+3;
  }*/
  
  
  inline G4double get_a_ion(){
     G4double mass[24]={7.,9.,11.,12.,14.,16.,19.,20.,23.,24.,27.,28.,31., 32.,     35.,40.,39., 40.,45.,48.,51.,52.,55.,56.};
     int iz=(int) z_ion;
     return mass[iz-3];
  }
  // gives back v/c as a function of kinetic Energy
  inline G4double beta(G4double E /* GeV */){
    return sqrt(1 - pow(E/restE+1, -2));
  }

  // gives back the rigidity (p/Ze where p is the momentum, e means 
  // electron charge magnitude, and Z is the atomic number) in units of [GV],
  // as a function of kinetic Energy [GeV].
  inline G4double rigidity(G4double E /* GeV */){
    return sqrt(pow(E + restE, 2) - pow(restE, 2))/z_ion;
  }

  // gives back the kinetic energy [GeV] as a function of rigidity [GV]
  inline G4double energy(G4double rigidity /* GV */){
    return sqrt(pow(rigidity*z_ion, 2) + pow(restE, 2)) - restE;
  }

  //============================================================
  /**
   *  Generate a random distribution of primary cosmic ray ions
   *  j(E) = mod_spec(E, phi) * geomag_cut(E, CutOff)
   *    mod_spec(E, phi) = org_spec(E+2*phi*1e-3) * 
   *     ((E+restE)**2 - restE**2)/((E+restE+2*phi*1e-3)**2-restE**2)
   *    org_spec(E) = A * rigidity(E)**-a
   *      A = 1.50 and a = 2.77 for alphas
   *      A = 0.204 and a =2.77 for ions 
   *    rigidity(E) = sqrt((E+restE)**2 - restE**2) / 2.0
   *    geomag_cut(E, CutOff) = 1/(1 + (rigidity(E)/CutOff)**-12.0)
   *      CutOff = 4.46 for Theta_M = 0.735 and altitude = 35km 
   *                                         (balloon experiment)
   *      phi = 540 and 1100 [MV] for Solar minimum and maximum, respectively
   *    E: [GeV]
   *    j: [c/s/m^2/sr/MeV]
   *
   *  References:
   *  org_spec: Set et al. 1991, ApJ 378, 763 (LEAP)
   *            Sanuki et al. 2000, ApJ 545, 1135 (BESS)
   *  geomag_cut formula: 
   *    an eyeball fitting function to represent the AMS proton data
   *  CutOff: calculated as (Rc/GV) = 14.9 * (1+h/R)^-2 * (cos(theta_M))^4,
   *          where h is the altitude from earth surface, 
   *                R is the mean radius of earth,
   *                and theta_M is geomagnetic lattitude. 
   *          References:
   *          "Handbook of space astronomy and astrophysics" 2nd edition, p225 
   *            (Zombeck, 1990, Cambridge University Press)
   *          "High Energy Astrophysics" 2nd edition, p325-330
   *            (M. S. Longair, 1992, Cambridge University Press)
   *  Solar modulation model (mod_spec) is from:
   *     Gleeson, L. J. and Axford, W. I. 1968, ApJ, 154, 1011-1026 (Eq. 11)
   */

  const G4double A_primary = 0.204; // normalization of incident spectrum,scaled from 1.5 for alphas
  const G4double a_primary = 2.77; // differential spectral index

  // Gives back the geomagnetic cutoff factor to the intrinsic 
  // primary cosmic ray spectrum for a kinetic energy E(GeV) 
  // and a geomagnetic lattitude theta_M(rad)
  inline G4double geomag_cut(G4double E, G4double cor /* GV */){
    return 1./(1 + pow(rigidity(E)/cor, -12.0));
  }

  // The unmodulated primary ion spectrum outside the Solar system is 
  // returned.
  inline G4double org_spec(G4double E /* GeV */)
  {
    return A_primary * pow(rigidity(E), -a_primary);
  }

  // The modulated ion flux for a "phi" value is returned. 
  // Force-field approximation of the Solar modulation is used.
  // The value of phi(potential) is doubled due to the charge of 2.
  inline G4double mod_spec(G4double E /* GeV */, G4double phi /* MV */){
    return org_spec(E + z_ion*phi*1e-3) * (pow(E+restE, 2) - pow(restE, 2))
      / (pow(E+restE+z_ion*phi*1e-3,2) - pow(restE,2));
  }

  // The final spectrum for the primary ion.
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
  // Primary spectrum of cosmic-ray ion is enveloped by a linear function
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
  // Primary spectrum of cosmic-ray ion is enveloped by power-law function
  // between cutE_primary and highE_primary
  inline G4double primaryCRenvelope2
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    return A_primary * pow(E/z_ion, -a_primary);
  }

  // The integral of the envelope function in the higher energy range
  inline G4double primaryCRenvelope2_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    return A_primary*z_ion/(-a_primary+1) * pow(E/z_ion, -a_primary+1);
  }

  // The inverse function of the integral of the envelope function
  // in the higher energy range.
  // This function returns energy obeying envelope function.
  inline G4double primaryCRenvelope2_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return z_ion*pow((-a_primary+1)/ (A_primary*z_ion) * value , 
		       1./(-a_primary+1));
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

    G4double r, E; // E means energy in GeV
    while(1){
      if (engine->flat() <= envelope1_area/(envelope1_area + envelope2_area)){

       // Use the envelop function in the lower energy range
        // (E<Ec where Ec corresponds to the cutoff rigidity).
        // We enveloped ion spectrum by linear function between
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
    {323.1, 270.4, 230.9, 200.5, 175.8, 155.9, 139.9}, // COR = 0.5GV
    {295.5, 251.5, 217.4, 190.3, 168.3, 150.1, 134.9}, // COR = 1 GV
    {205.3, 182.6, 163.7, 147.6, 131.8, 121.9, 111.6}, // COR = 2 GV
    {138.3, 126.7, 116.6, 107.6, 99.6, 92.4, 86.1}, // COR = 3 GV
    {97.4, 91.0, 85.1, 79.7, 75.0, 70.5, 66.4}, // COR = 4 GV
    {72.0, 68.0, 64.4, 61.0, 57.9, 55.0, 52.7}, // COR = 5 GV
    {55.4, 52.8, 50.3, 48.1, 45.9, 43.9, 42.1}, // COR = 6 GV
    {44.0, 42.2, 40.5, 38.8, 37.3, 35.9, 34.5}, // COR = 7 GV
    {35.8, 34.5, 33.3, 32.1, 31.0, 29.0, 28.9}, // COR = 8 GV
    {29.7, 28.8, 27.8, 27.0, 26.1, 25.3, 24.5}, // COR = 9 GV
    {25.1, 24.4, 23.8, 23.0, 22.3, 21.7, 21.1}, // COR = 10 GV
    {21.6, 21.0, 20.4, 19.9, 19.3, 18.8, 18.3}, // COR = 11 GV
    {18.7, 18.2, 17.8, 17.3, 16.9, 16.5, 16.1}, // COR = 12 GV
    {16.4, 16.0, 15.6, 15.3, 14.9, 14.6, 14.3}, // COR = 13 GV
    {14.5, 14.2, 13.9, 13.6, 13.3, 13.0, 12.7}, // COR = 14 GV
    {12.9, 12.6, 12.4, 12.1, 11.9, 11.6, 11.4} // COR = 15 GV
  };
  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrHeavyIonPrimVertZ::CrHeavyIonPrimVertZ(int z)
{
  // Set lower and higher energy limit of the primary ion (GeV).
  // At lowE_primary, flux of primary ion can be 
  // assumed to be 0, due to geomagnetic cutoff 
  mm_engine = HepRandom::getTheEngine(); //new HepJamesRandom;  
  z_ion = z;
}


CrHeavyIonPrimVertZ::~CrHeavyIonPrimVertZ()
{
  ;
}


// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrHeavyIonPrimVertZ::setPosition(G4double latitude, G4double longitude){
  CrSpectrum::setPosition(latitude, longitude);

  // Set lower and higher energy limit of the primary ion (GeV).
  // At lowE_primary, flux of primary ion can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 50.*A_ion;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrHeavyIonPrimVertZ::setPosition
(G4double latitude, G4double longitude, G4double time){
  CrSpectrum::setPosition(latitude, longitude, time);
  // Set lower and higher energy limit of the primary ion (GeV).
  // At lowE_primary, flux of primary ion can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 50.*A_ion ;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);
  
}

// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrHeavyIonPrimVertZ::
setPosition(G4double latitude, G4double longitude, 
	    G4double time, G4double altitude){
  CrSpectrum::setPosition(latitude, longitude, time, altitude);
  // Set lower and higher energy limit of the primary ion (GeV).
  // At lowE_primary, flux of primary ion can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 50.*A_ion;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set geomagnetic cutoff rigidity and calculate the energies related.
// These energies are used to generate the particle. 
void CrHeavyIonPrimVertZ::setCutOffRigidity(G4double cor){
  CrSpectrum::setCutOffRigidity(cor);
  // Set lower and higher energy limit of the primary ion (GeV).
  // At lowE_primary, flux of primary ion can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 50.0*A_ion;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Gives back particle direction in (cos(theta), phi)
std::pair<G4double,G4double> CrHeavyIonPrimVertZ::dir(G4double energy, 
                                              CLHEP::HepRandomEngine* engine) const
 // return: cos(theta) and phi [rad]
  // The downward direction has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0).
{
  // We assume isotropic distribution from the upper hemisphere.
  // After integration over the azimuth angle (phi), 
  // the theta distribution should be sin(theta) for a constant theta width.
  /// Cos(theta) ranges from 1 to -0.4
  //G4double theta = acos(1.4*engine->flat()-0.4);
  //CL: to define Vertical CrHeavyIon, we set theta=0;
  G4double theta = 0;
  G4double phi   = engine->flat() * 2 * M_PI;

  return std::pair<G4double,G4double>(cos(theta), phi);
}


// Gives back particle energy
G4double CrHeavyIonPrimVertZ::energySrc(CLHEP::HepRandomEngine* engine) const
{ 
  return primaryCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".
G4double CrHeavyIonPrimVertZ::flux() const
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
    phi = phi/100.0 - 5; // 500 MV corresponds to 0, 600 MV corresponds to 1, et..
    G4double tmp1 = 
      integral_array[int(cor)][int(phi)] +
      (cor - int(cor)) * (integral_array[int(cor)+1][int(phi)]-integral_array[int(cor)][int(phi)]);
    G4double tmp2 = 
      integral_array[int(cor)][int(phi)+1] +
      (cor - int(cor)) * (integral_array[int(cor)+1][int(phi)+1]-integral_array[int(cor)][int(phi)+1]);
    energy_integral = tmp1 + (tmp2-tmp1)*(phi-int(phi));
  }

  if (cor < 1.0){
    // 500 MV corresponds to 0, 600 MV corresponds to 1, etc..
    phi = phi/100.0 - 5; 
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
  return energy_integral/7.34;  // [c/s/m^2/sr]  7.34 is the flux ratio between alphas and heavy ions
}

// Gives back solid angle from which particle comes
G4double CrHeavyIonPrimVertZ::solidAngle() const
{
  // * 1.4 since Cos(theta) ranges from 1 to -0.4 
  //z_ion=get_z_ion();
  //C.L: to fix z:
  //z_ion = z;
  std::cout << "CrHeavyIonPrimVertZ::solidAngle(), z_ion=" << z_ion << std::endl;  
  A_ion= get_a_ion();
  // std::cout << z_ion << " " << A_ion << std::endl;  
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 50.*A_ion; // corresponds to 100 GeV/n! not any more
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity); 
  return  2 * M_PI * 1.4;
}

// Gives back particle name
const  char* CrHeavyIonPrimVertZ::particleName() const 
{
  char* nameIon[24]={"Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe"};
  
  int zz=(int) z_ion;
  //  std::cout << " in hi name"<< zz << " " << nameIon[zz-3] <<std::endl;
 
  return nameIon[zz-3];
}


// Gives back the name of the component
std::string CrHeavyIonPrimVertZ::title() const
{
  return "CrHeavyIonPrimary";
}


