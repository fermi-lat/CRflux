//
// CrExample.cc
//
//   An example flux.
//

#include <math.h>
#include <map>
#include <vector>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>

#include "CrExample.h"

// define a factory for anonomous instantiation
#include "FluxSvc/ISpectrumFactory.h"

CrExample::CrExample(const std::string& paramstring)
{
   std::vector<float> params;
}


CrExample::~CrExample()
{}

double CrExample::energy()
{
  return 1.0;
}


std::pair<double,double> CrExample::dir(double energy)
// return: cos(zenith_angle) and azimuth [rad]
{
    return std::make_pair<double,double>(0.0,1.0);
}

double CrExample::flux (double time ) const
{
  double          total_flux = 0.;
  return total_flux;
}

double CrExample::solidAngle( )const
{
    return 4*M_PI;
}


double CrExample::energy(double time)
{
    return energy();
}

double CrExample::interval (double time)
{      
        return 1.0;
}
