//$Header$


#ifndef CrExample_H
#define CrExample_H


//!  The class that calls each cosmic-ray Example components based on
//!  their flux.

#include <vector>
#include <utility>
#include <string>
#include "FluxSvc/ISpectrum.h"

class CrExample : public /*CrSpectrum*/ISpectrum
{
public:
    // params[0] is bit flag for determining which to include (default 7)
    // 1: primary
    // 2: reentrant
    // 4: splash
    
    CrExample(const std::string& params);
    ~CrExample();
    
    virtual double energy();
    std::pair<double,double> dir(double energy);
    
    //! calculate the flux, particles/m^2/sr.
    virtual double    flux (double time ) const;
    
    virtual const char * particleName()const{ return "p";}
    virtual std::string title()const{return "CrExample";}
    virtual double solidAngle( )const;
    
    virtual double energy(double time);
    
    
    virtual double interval (double time);

private:
    
};
#endif // CrExample_H

