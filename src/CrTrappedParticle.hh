/**
 * CrProtonPrimary:
 *   The primary cosmic-ray proton spectrum (and incident angle) source.
 */

//$Header$

#ifndef CrTrappedParticle_H
#define CrTrappedParticle_H

#include <utility>
#include <string>
#include <map>

#include "CrSpectrum.hh"

typedef double G4double;


// Forward declaration:
class CLHEP::HepRandomEngine;

class CrTrappedParticle : public CrSpectrum
{

protected:  
  G4double m_spectrumLatitude, m_spectrumLongitude, m_spectrumAltitude;
  std::map<G4double,G4double> m_intSpectrum;
  G4double m_integralFlux;
  G4double m_maxNonzeroFluxEnergy;
  G4double m_thresholdEnergy,m_eStep,m_eMax,m_modelMinEnergy,m_modelMaxEnergy;
  std::string m_model;
  
  enum {invalid,proton,electron} m_particleType;
  
  
private:  
  std::string m_serverAddress;
  std::string m_xmlDirectory;
  int m_socketHandle;


public:
  CrTrappedParticle(const std::string& paramstring,const std::string& ptype);
  ~CrTrappedParticle();


  // Gives back particle direction in (cos(theta), phi)
  std::pair<double,double> dir(double energy, CLHEP::HepRandomEngine* engine) const;

  // Gives back particle energy
  double energySrc(CLHEP::HepRandomEngine* engine) const;

  // flux() returns the value averaged over the region from which
  // the particle is coming from and the unit is [c/s/m^2/sr]
  double flux() const;

  // Gives back solid angle from which particle comes
  double solidAngle() const;

  // Gives back particle name
  const char* particleName() const;

  // Gives back the name of the component
  std::string title() const;
  
protected:
   bool checkModelCompatibility(const std::string& model,const std::string& particle);
   bool coordinatesChanged() const;
   bool requestNewSpectrum(const G4double minE,G4double maxE,G4double stepE);
   bool psb97UpdateSpectrum(const G4double minE,G4double maxE,G4double stepE);
   void connectToServer();
   void disconnectFromServer();

};

// convenience classes .....

class CrTrappedProton : public CrTrappedParticle {
public:
  CrTrappedProton(const std::string& paramstring): CrTrappedParticle(paramstring,"p"){};
  ~CrTrappedProton(){};      
  
  const char* particleName() const {  return "proton"; };
  std::string title() const {  return  "CrTrappedProton"; };


};

class CrTrappedPositron : public CrTrappedParticle {
public:
  CrTrappedPositron(const std::string& paramstring): CrTrappedParticle(paramstring,"e+"){};
  ~CrTrappedPositron(){};   
  
  const char* particleName() const {  return "positron"; };
  std::string title() const {  return  "CrTrappedPositron"; };

     
};

class CrTrappedElectron : public CrTrappedParticle {
public:
  CrTrappedElectron(const std::string& paramstring): CrTrappedParticle(paramstring,"e-"){};
  ~CrTrappedElectron(){};      

  const char* particleName() const {  return "electron"; };
  std::string title() const {  return  "CrTrappedElectron"; };
};



#endif // CrTrappedParticle_H

