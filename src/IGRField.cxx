#include <cmath>
#include "igrf_sub/igrf_sub.h"
#include "IGRField.h"




IGRField IGRField::m_model;

IGRField::IGRField() {
    IGRFf2c::initize_(); 
    setYear(2005.); 
    compute(0.,0.,0.);
}; 

void IGRField::setYear(const float year){
   if(fabs(year-m_year)<0.001) return;
   m_year=year;
   if (m_year>=2010.) m_year=2010.;
   if (m_year<=1990.) m_year=1990.;
   IGRFf2c::feldcof_(&m_year,&m_dipoleMoment);
};      
      
int IGRField::compute(const float latitude,const float longitude,const float altitude,const float year){

   float earth_radius= 6371.2 ;

   float b0,rr0;
   IGRFf2c::integer icode;
   IGRFf2c::logical value;
   float stps=0.05;
   float bdel=0.001; 

   if(year>0) m_model.setYear(year);
   m_latitude=latitude;
   m_longitude=longitude;
   m_altitude=altitude;
   
   IGRFf2c::feldg_(&m_latitude,&m_longitude,&m_altitude,&m_bNorth,&m_bEast,&m_bDown,&m_bAbs);
   IGRFf2c::shellg_(&m_latitude,&m_longitude,&m_altitude,&m_dipoleMoment,&m_L,&icode,&b0) ;
   
   IGRFf2c::findb0_(&stps, &bdel, &value, &m_bEquator, &rr0);
   if(value==0) m_bEquator=m_dipoleMoment/(m_L*m_L*m_L);
   
   m_B=m_bAbs/m_bEquator;   
   
// conversion from dipole moment in IGRF G*R_earth^3 to rigidity in GV:
//                                                R_earth   km->cm  statvolt -> gigavolt
   float rigidity_const= 0.25 * m_dipoleMoment * earth_radius *  1e5 *  300. / 1e9 ;
   
// the invariant latitude and the rigidity cutoff
// as described in Smart and Shea, Adv. Space Res. 36 (2005) p.2012

   m_invariantLat=m_L>=1 ? acos(sqrt(1./m_L)): 0.;
   m_rigidityCutoff=rigidity_const/(m_L*m_L);
   
// and alternatively the r, lambda coordinates 
// from http://www.spenvis.oma.be/spenvis/help/background/magfield/rlambda.html 
// gives identical rigidity cutoffs if rigidity cutoff is calculated as 
// rigidity_const * cos(m_lambda)^4 / m_R^2

// needs an approximation to the solution rl(b) of the equation b^2*rl^6+3*rl-4==0 in the range b=1..10
   double rl= pow((double)m_B,-0.215108)*(1.-0.020551*((double)m_B-1.)+0.0008148*((double)m_B-1.));
   m_R= rl*m_L;
   m_lambda= (rl<=1) ? acos(sqrt(rl)) : 0.;

   return icode;

};



