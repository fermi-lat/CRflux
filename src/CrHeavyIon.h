//$Header$


#ifndef CrHeavyIon_H
#define CrHeavyIon_H


//!  The class that calls each cosmic-ray Example components based on
//!  their flux.

#include <vector>
#include <utility>
#include <string>
#include "flux/ISpectrum.h"

class CrSpectrum;
class CLHEP::HepRandomEngine;

class CrHeavyIon : public ISpectrum
{
public:
    // params[0] is bit flag for determining which to include (default 7)
    // 1: primary
    // 2: reentrant
    // 4: splash
    
    CrHeavyIon(const std::string& params);
    virtual ~CrHeavyIon();

   // Gives back component in the ratio of the flux
  CrSpectrum* selectComponent();
   
  //virtual double energy();
    std::pair<double,double> dir(double energy);
    
    //! calculate the flux, particles/m^2/sr.
    virtual double    flux (double time ) const;
    
    virtual const char * particleName()const; // BL{ return "p";}
    virtual std::string title()const{return "CrHeavyIon";}
    virtual double solidAngle( )const;
    
    virtual double energy(double time);
    virtual double interval (double time);
    void dump();
private:
 std::vector<CrSpectrum*>  m_subComponents;
  CrSpectrum*               m_component;
  CLHEP::HepRandomEngine* m_engine;    
};
#endif // CrExample_H

