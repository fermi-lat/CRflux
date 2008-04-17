#include <cmath>

#include "igrf_sub.h"

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define TRUE_ (1)
#define FALSE_ (0)


namespace IGRFf2c {

 typedef double doublereal;
 doublereal signc_(real *, real *);
 int stoer_(real *, real *, real *);
 int extrashc_(real *, real *, integer *, real *, integer *, real *, integer *, real *);
 int getshc_(real *, integer *, real *, real *, integer *);
 int intershc_(real *, real *, integer *, real *, 
    real *, integer *, real *, integer *, real *), extrashc_(real *, 
    real *, integer *, real *, integer *, real *, integer *, real *);


/* igrf_sub.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

    static real dgrf90[263] = { 10.,6371.2,1990.,
	    1.,0.,-29775.,0.,1.,1., -1848.,5406.,2.,0.,-2131.,0.,2.,1.,3059.,-2279.,2.,2.,
	    1686.,-373.,3.,0.,1314.,0.,3.,1.,-2239.,-284.,3.,2.,1248.,293.,3.,3.,802.,
	    -352.,4.,0.,939.,0.,4.,1.,780.,247.,4.,2.,325.,-240.,4.,3.,-423.,84.,4.,
	    4.,141.,-299.,5.,0.,-214.,0.,5.,1.,353.,46.,5.,2.,245.,154.,5.,3.,
	    -109.,-153.,5.,4.,-165.,-69.,5.,5.,-36.,97.,6.,0.,61.,0.,6.,1.,65.,-16.,
	    6.,2.,59.,82.,6.,3.,-178.,69.,6.,4.,3.,-52.,6.,5.,18.,1.,6.,6.,
	    -96.,24.,7.,0.,77.,0.,7.,1.,-64.,-80.,7.,2.,2.,-26.,7.,3.,26.,0.,
	    7.,4.,-1.,21.,7.,5.,5.,17.,7.,6.,9.,-23.,7.,7.,0.,-4.,8.,0.,
	    23.,0.,8.,1.,5.,10.,8.,2.,-1.,-19.,8.,3.,-10.,6.,8.,4.,-12.,-22.,
	    8.,5.,3.,12.,8.,6.,4.,12.,8.,7.,2.,-16.,8.,8.,-6.,-10.,9.,0.,
	    4.,0.,9.,1.,9.,-20.,9.,2.,1.,15.,9.,3.,-12.,11.,9.,4.,9.,-7.,
	    9.,5.,-4.,-7.,9.,6.,-2.,9.,9.,7.,7.,8.,9.,8.,1.,-7.,9.,9.,
	    -6.,2.,10.,0.,-3.,0.,10.,1.,-4.,2.,10.,2.,2.,1.,10.,3.,-5.,3.,
	    10.,4.,-2.,6.,10.,5.,4.,-4.,10.,6.,3.,0.,10.,7.,1.,-2.,10.,8.,
	    3.,3.,10.,9.,3.,-1.,10.,10.,0.,-6. };
	    
    static real dgrf95[263] = { 10.,6371.2,1995.,
	    1.,0.,-29692.,0.,1.,1.,-1784.,5306.,2.,0.,-2200.,0.,2.,1.,3070.,-2366.,2.,2.,
	    1681.,-413.,3.,0.,1335.,0.,3.,1.,-2267.,-262.,3.,2.,1249.,302.,3.,3.,759.,
	    -427.,4.,0.,940.,0.,4.,1.,780.,262.,4.,2.,290.,-236.,4.,3.,-418.,97.,4.,
	    4.,122.,-306.,5.,0.,-214.,0.,5.,1.,352.,46.,5.,2.,235.,165.,5.,3.,
	    -118.,-143.,5.,4.,-166.,-55.,5.,5.,-17.,107.,6.,0.,68.,0.,6.,1.,67.,-17.,
	    6.,2.,68.,72.,6.,3.,-170.,67.,6.,4.,-1.,-58.,6.,5.,19.,1.,6.,6.,
	    -93.,36.,7.,0.,77.,0.,7.,1.,-72.,-69.,7.,2.,1.,-25.,7.,3.,28.,4.,
	    7.,4.,5.,24.,7.,5.,4.,17.,7.,6.,8.,-24.,7.,7.,-2.,-6.,8.,0.,
	    25.,0.,8.,1.,6.,11.,8.,2.,-6.,-21.,8.,3.,-9.,8.,8.,4.,-14.,-23.,
	    8.,5.,9.,15.,8.,6.,6.,11.,8.,7.,-5.,-16.,8.,8.,-7.,-4.,9.,0.,
	    4.,0.,9.,1.,9.,-20.,9.,2.,3.,15.,9.,3.,-10.,12.,9.,4.,8.,-6.,
	    9.,5.,-8.,-8.,9.,6.,-1.,8.,9.,7.,10.,5.,9.,8.,-2.,-8.,9.,9.,
	    -8.,3.,10.,0.,-3.,0.,10.,1.,-6.,1.,10.,2.,2.,0.,10.,3.,-4.,4.,
	    10.,4.,-1.,5.,10.,5.,4.,-5.,10.,6.,2.,-1.,10.,7.,2.,-2.,10.,8.,
	    5.,1.,10.,9.,1.,-2.,10.,10.,0.,-7. };
	    
    static real dgrf00[263] = { 10.,6371.2,2e3,1.,
	    0.,-29619.4,0.,1.,1.,-1728.2,5186.1,2.,0.,-2267.7,0.,2.,1.,3068.4,-2481.6,2.,
	    2.,1670.9,-458.,3.,0.,1339.6,0.,3.,1.,-2288.,-227.6,3.,2.,1252.1,293.4,3.,3.,
	    714.5,-491.1,4.,0.,932.3,0.,4.,1.,786.8,272.6,4.,2.,250.,-231.9,4.,3.,
	    -403.,119.8,4.,4.,111.3,-303.8,5.,0.,-218.8,0.,5.,1.,351.4,43.8,5.,2.,222.3,
	    171.9,5.,3.,-130.4,-133.1,5.,4.,-168.6,-39.3,5.,5.,-12.9,106.3,6.,0.,72.3,0.,
	    6.,1.,68.2,-17.4,6.,2.,74.2,63.7,6.,3.,-160.9,65.1,6.,4.,-5.9,-61.2,6.,
	    5.,16.9,.7,6.,6.,-90.4,43.8,7.,0.,79.,0.,7.,1.,-74.,-64.6,7.,2.,0.,
	    -24.2,7.,3.,33.3,6.2,7.,4.,9.1,24.,7.,5.,6.9,14.8,7.,6.,7.3,-25.4,7.,
	    7.,-1.2,-5.8,8.,0.,24.4,0.,8.,1.,6.6,11.9,8.,2.,-9.2,-21.5,8.,3.,-7.9,
	    8.5,8.,4.,-16.6,-21.5,8.,5.,9.1,15.5,8.,6.,7.,8.9,8.,7.,-7.9,-14.9,8.,
	    8.,-7.,-2.1,9.,0.,5.,0.,9.,1.,9.4,-19.7,9.,2.,3.,13.4,9.,3.,-8.4,
	    12.5,9.,4.,6.3,-6.2,9.,5.,-8.9,-8.4,9.,6.,-1.5,8.4,9.,7.,9.3,3.8,9.,
	    8.,-4.3,-8.2,9.,9.,-8.2,4.8,10.,0.,-2.6,0.,10.,1.,-6.,1.7,10.,2.,1.7,
	    0.,10.,3.,-3.1,4.,10.,4.,-.5,4.9,10.,5.,3.7,-5.9,10.,6.,1.,-1.2,10.,
	    7.,2.,-2.9,10.,8.,4.2,.2,10.,9.,.3,-2.2,10.,10.,-1.1,-7.4 };
	    
    static real igrf05[263] = { 10.,6371.2,2005.,
	    1.,0.,-29556.8,0.,1.,1.,-1671.8,5080.,2.,0.,-2340.5,0.,2.,1.,3047.,-2594.9,2.,
	    2.,1656.9,-516.7,3.,0.,1335.7,0.,3.,1.,-2305.3,-200.4,3.,2.,1246.8,269.3,3.,3.,
	    674.4,-524.5,4.,0.,919.8,0.,4.,1.,798.2,281.4,4.,2.,211.5,-225.8,4.,3.,-379.5,
	    145.7,4.,4.,100.2,-304.7,5.,0.,-227.6,0.,5.,1.,354.4,42.7,5.,2.,208.8,179.8,
	    5.,3.,-136.6,-123.,5.,4.,-168.3,-19.5,5.,5.,-14.1,103.6,6.,0.,72.9,0.,6.,
	    1.,69.6,-20.2,6.,2.,76.6,54.7,6.,3.,-151.1,63.7,6.,4.,-15.,-63.4,6.,5.,
	    14.7,0.,6.,6.,-86.4,50.3,7.,0.,79.8,0.,7.,1.,-74.4,-61.4,7.,2.,-1.4,
	    -22.5,7.,3.,38.6,6.9,7.,4.,12.3,25.4,7.,5.,9.4,10.9,7.,6.,5.5,-26.4,7.,
	    7.,2.,-4.8,8.,0.,24.8,0.,8.,1.,7.7,11.2,8.,2.,-11.4,-21.,8.,3.,-6.8,
	    9.7,8.,4.,-18.,-19.8,8.,5.,10.,16.1,8.,6.,9.4,7.7,8.,7.,-11.4,-12.8,8.,
	    8.,-5.,-.1,9.,0.,5.6,0.,9.,1.,9.8,-20.1,9.,2.,3.6,12.9,9.,3.,-7.,
	    12.7,9.,4.,5.,-6.7,9.,5.,-10.8,-8.1,9.,6.,-1.3,8.1,9.,7.,8.7,2.9,9.,
	    8.,-6.7,-7.9,9.,9.,-9.2,5.9,10.,0.,-2.2,0.,10.,1.,-6.3,2.4,10.,2.,1.6,
	    .2,10.,3.,-2.5,4.4,10.,4.,-.1,4.7,10.,5.,3.,-6.5,10.,6.,.3,-1.,10.,
	    7.,2.1,-3.4,10.,8.,3.9,-.9,10.,9.,-.1,-2.3,10.,10.,-2.2,-8. };
	    
    static real igrf05s[179] = { 8.,6371.2,2010.,
	    1.,0.,8.8,0.,1.,1.,10.8,-21.3,2.,0.,-15.,0.,2.,1.,-6.9,-23.3,2.,2.,-1.,
	    -14.,3.,0.,-.3,0.,3.,1.,-3.1,5.4,3.,2.,-.9,-6.5,3.,3.,-6.8,-2.,4.,
	    0.,-2.5,0.,4.,1.,2.8,2.,4.,2.,-7.1,1.8,4.,3.,5.9,5.6,4.,4.,-3.2,
	    0.,5.,0.,-2.6,0.,5.,1.,.4,.1,5.,2.,-3.,1.8,5.,3.,-1.2,2.,5.,
	    4.,.2,4.5,5.,5.,-.6,-1.,6.,0.,-.8,0.,6.,1.,.2,-.4,6.,2.,-.2,
	    -1.9,6.,3.,2.1,-.4,6.,4.,-2.1,-.4,6.,5.,-.4,-.2,6.,6.,1.3,.9,7.,
	    0.,-.4,0.,7.,1.,0.,.8,7.,2.,-.2,.4,7.,3.,1.1,.1,7.,4.,.6,
	    .2,7.,5.,.4,-.9,7.,6.,-.5,-.3,7.,7.,.9,.3,8.,0.,-.2,0.,8.,
	    1.,.2,-.2,8.,2.,-.2,.2,8.,3.,.2,.2,8.,4.,-.2,.4,8.,5.,.2,.2,
	    8.,6.,.5,-.3,8.,7.,-.7,.5,8.,8.,.5,.4 };

/* Common Block Declarations */

union {
    struct {
	real umr, era, aquad, bquad;
    } _1;
    struct {
	real umr, erad, aquad, bquad;
    } _2;
} gener_;

#define gener_1 (gener_._1)
#define gener_2 (gener_._2)

struct {
    real sp[3];
} fidb0_;

#define fidb0_1 fidb0_

union {
    struct {
	real x[3], h__[144];
    } _1;
    struct {
	real xi[3], h__[144];
    } _2;
} _BLNK__;

#define _BLNK__1 (_BLNK__._1)
#define _BLNK__2 (_BLNK__._2)

union {
    struct {
	char name__[12];
	integer nmax;
	real time, g[144];
    } _1;
    struct {
	char fil1[12];
	integer nmax;
	real time, gh1[144];
    } _2;
} model_;

#define model_1 (model_._1)
#define model_2 (model_._2)

/* IGRF_SUB.FOR */
/* ********************************************************************* */
/*  SUBROUTINES igrf_sub plus SHELLIG.FOR (see below)                 * */
/* ********************************************************************* */
/* ********************************************************************* */

/* 11/01/91 SHELLG: lowest starting point for B0 search is 2 */
/*  1/27/92 Adopted to IGRF-91 coeffcients model */
/*  2/05/92 Reduce variable names: INTER(P)SHC,EXTRA(P)SHC,INITI(ALI)ZE */
/*  8/08/95 Updated to IGRF-45-95; new coeff. DGRF90, IGRF95, IGRF95S */
/*  5/31/00 Updated to IGRF-45-00; new coeff.: IGRF00, IGRF00s */
/* -Version-mm/dd/yy-Description (Person reporting the correction) */
/* 2000.01 05/07/01 initial version */
/* 2000.02 07/11/01 replace feldi(xi,h) by feldi (P. Wilkinson) */
/* 2000.02 07/11/01 variables EGNR, AGNR,OGNR not used (P. Wilkinson) */
/* 2000.01 10/28/02 replace TAB/6 blanks, enforce 72/line (D. Simpson) */
/* 2000.02 11/08/02 change unit for coefficients to 14 */
/* 2000.03 06/05/03 correct DIPL computation (V. Truhlik) */
/* 2005.00 04/25/05 CALL FELDI and DO 1111 I=1,7 (Alexey Petrov) */
/* 2005.01 11/10/05 added igrf_dip and geodip (MLAT) */
/* 2005.02 11/10/05 updated to IGRF-10 version */
/* 2006.00 12/21/06 GH2(120) -> GH2(144) */

/* ********************************************************************* */
/* Subroutine */ int igrf_sub__(real *xlat, real *xlong, real *year, real *
	height, real *xl, integer *icode, real *dip, real *dec)
{

    /* Local variables */
    static real bab1;
    static integer ibbb;
    static real babs, dimo, lati, alog2;
    int feldg_(real *, real *, real *, real *, real *,
	     real *, real *);
    static real beast, longi, bdown;
    int shellg_(real *, real *, real *, real *, real *
	    , integer *, real *);
    static real bnorth;
    static integer istart;
    int feldcof_(real *, real *), initize_();

/* ---------------------------------------------------------------- */
/*   INPUT: */
/* 	xlat	geodatic latitude in degrees */
/* 	xlong	geodatic longitude in degrees */
/* 	year	decimal year (year+month/12.0-0.5 or year+day-of-year/365 */
/* 		or 366 if leap year) */
/* 	height	height in km */
/*   OUTPUT: */
/* 	xl	L value */
/* 	icode	=1  L is correct; =2  L is not correct; */
/* 		=3  an approximation is used */
/* 	dip	geomagnetic inclination in degrees */
/* 	dec	geomagnetic declination in degress */
/* ---------------------------------------------------------------- */

    initize_();
    ibbb = 0;
    alog2 = log((float)2.);
    istart = 1;
    lati = *xlat;
    longi = *xlong;

/* ----------------CALCULATE PROFILES----------------------------------- */

    feldcof_(year, &dimo);
    feldg_(&lati, &longi, height, &bnorth, &beast, &bdown, &babs);
    shellg_(&lati, &longi, height, &dimo, xl, icode, &bab1);
    *dip = asin(bdown / babs) / gener_1.umr;
    *dec = asin(beast / sqrt(beast * beast + bnorth * bnorth)) / gener_1.umr;
    return 0;
} /* igrf_sub__ */



/* SHELLIG.FOR */

/* ********************************************************************* */
/*  SUBROUTINES FINDB0, SHELLG, STOER, FELDG, FELDCOF, GETSHC,        * */
/*       INTERSHC, EXTRASHC, INITIZE                                  * */
/* ********************************************************************* */
/* ********************************************************************* */


/* Subroutine */ int findb0_(real *stps, real *bdel, logical *value, real *
	bequ, real *rr0)
{

    /* Local variables */
    static real b;
    static integer i__, j, n;
    static real p[32]	/* was [8][4] */, r1, r2, r3, zz, bq1, bq2, bq3, bold,
	     bmin, rold, step;
    static integer irun;
    doublereal signc_(real *, real *);
    static real step12;
    int stoer_(real *, real *, real *);
    static real bdelta;

/* -------------------------------------------------------------------- */
/* FINDS SMALLEST MAGNETIC FIELD STRENGTH ON FIELD LINE */

/* INPUT:   STPS   STEP SIZE FOR FIELD LINE TRACING */
/*       COMMON/FIDB0/ */
/*          SP     DIPOLE ORIENTED COORDINATES FORM SHELLG; P(1,*), */
/*                 P(2,*), P(3,*) CLOSEST TO MAGNETIC EQUATOR */
/*          BDEL   REQUIRED ACCURACY  = [ B(LAST) - BEQU ] / BEQU */
/*                 B(LAST)  IS FIELD STRENGTH BEFORE BEQU */

/* OUTPUT:  VALUE  =.FALSE., IF BEQU IS NOT MINIMAL VALUE ON FIELD LINE */
/*          BEQU   MAGNETIC FIELD STRENGTH AT MAGNETIC EQUATOR */
/*          RR0    EQUATORIAL RADIUS NORMALIZED TO EARTH RADIUS */
/*          BDEL   FINAL ACHIEVED ACCURACY */
/* -------------------------------------------------------------------- */
/* f2py intent(in) stps */
/* f2py intent(in) bdel */
/* f2py intent(out) value */
/* f2py intent(out) bequ */
/* f2py intent(out) rr0 */

    step = *stps;
    irun = 0;
L7777:
    ++irun;
    if (irun > 5) {
	*value = FALSE_;
	goto L8888;
    }
/* *********************FIRST THREE POINTS */
    p[8] = fidb0_1.sp[0];
    p[9] = fidb0_1.sp[1];
    p[10] = fidb0_1.sp[2];
    step = -signc_(&step, &p[10]);
    
    stoer_(&p[8], &bq2, &r2);
    p[16] = p[8] + step * (float).5 * p[11];
    p[17] = p[9] + step * (float).5 * p[12];
    p[18] = p[10] + step * (float).5;
    stoer_(&p[16], &bq3, &r3);
    p[0] = p[8] - step * (p[11] * (float)2. - p[19]);
    p[1] = p[9] - step * (p[12] * (float)2. - p[20]);
    p[2] = p[10] - step;
    stoer_(p, &bq1, &r1);
    p[16] = p[8] + step * (p[19] * (float)20. - p[11] * (float)3. + p[3]) / (
	    float)18.;
    p[17] = p[9] + step * (p[20] * (float)20. - p[12] * (float)3. + p[4]) / (
	    float)18.;
    p[18] = p[10] + step;
    stoer_(&p[16], &bq3, &r3);
/* ******************INVERT SENSE IF REQUIRED */
    if (bq3 <= bq1) {
	goto L2;
    }
    step = -step;
    r3 = r1;
    bq3 = bq1;
    for (i__ = 1; i__ <= 5; ++i__) {
	zz = p[i__ - 1];
	p[i__ - 1] = p[i__ + 15];
/* L1: */
	p[i__ + 15] = zz;
    }
/* ******************INITIALIZATION */
L2:
    step12 = step / (float)12.;
    *value = TRUE_;
    bmin = (float)1e4;
    bold = (float)1e4;
/* ******************CORRECTOR (FIELD LINE TRACING) */
    n = 0;
L5555:
    p[16] = p[8] + step12 * (p[19] * (float)5. + p[11] * (float)8. - p[3]);
    ++n;
    p[17] = p[9] + step12 * (p[20] * (float)5. + p[12] * (float)8. - p[4]);
/* ******************PREDICTOR (FIELD LINE TRACING) */
    p[24] = p[16] + step12 * (p[19] * (float)23. - p[11] * (float)16. + p[3] *
	     (float)5.);
    p[25] = p[17] + step12 * (p[20] * (float)23. - p[12] * (float)16. + p[4] *
	     (float)5.);
    p[26] = p[18] + step;
    stoer_(&p[24], &bq3, &r3);
    for (j = 1; j <= 3; ++j) {
/*        DO 1111 I=1,8 */
	for (i__ = 1; i__ <= 7; ++i__) {
/* L1111: */
	    p[i__ + (j << 3) - 9] = p[i__ + (j + 1 << 3) - 9];
	}
    }
    b = sqrt(bq3);
    if (b < bmin) {
	bmin = b;
    }
    if (b <= bold) {
	bold = b;
	rold = (float)1. / r3;
	fidb0_1.sp[0] = p[24];
	fidb0_1.sp[1] = p[25];
	fidb0_1.sp[2] = p[26];
	goto L5555;
    }
    if (bold != bmin) {
	*value = FALSE_;
    }
    bdelta = (b - bold) / bold;
    if (bdelta > *bdel) {
	step /= (float)10.;
	goto L7777;
    }
L8888:
    *rr0 = rold;
    *bequ = bold;
    *bdel = bdelta;
    return 0;
} /* findb0_ */



/* Subroutine */ int shellg_0_(int n__, real *glat, real *glon, real *alt, 
	real *dimo, real *fl, integer *icode, real *b0, real *v)
{
    /* Initialized data */

    static real rmin = (float).05;
    static real rmax = (float)1.01;
    static real step = (float).2;
    static real steq = (float).03;
    static real u[9]	/* was [3][3] */ = { (float).3511737,(float)-.9148385,
	    (float)-.1993679,(float).9335804,(float).358368,(float)0.,(float)
	    .0714471,(float)-.186126,(float).9799247 };

    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static real d__;
    static integer i__, n;
    static real p[800]	/* was [8][100] */, r__, t, z__, c0, c1, c2, c3, d0, 
	    d1, d2, e0, e1, e2, r1, r2, r3, ff, gg, fi, ct, rq, st, zq, xx, 
	    zz, bq1, bq2, bq3, r3h, hli, stp, arg1, arg2, bequ, rlat;
    static integer iequ;
    static real term, rlon, step2, radik;
    doublereal signc_(real *, real *);
    static real step12, oterm;
    int stoer_(real *, real *, real *);
    static real dimob0, oradik;

/* -------------------------------------------------------------------- */
/* CALCULATES L-VALUE FOR SPECIFIED GEODAETIC COORDINATES, ALTITUDE */
/* AND GEMAGNETIC FIELD MODEL. */
/* REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTER, INTERNAL NOTE */
/*      NO. 67, 1970. */
/*      G. KLUGE, COMPUTER PHYSICS COMMUNICATIONS 3, 31-35, 1972 */
/* -------------------------------------------------------------------- */
/* CHANGES (D. BILITZA, NOV 87): */
/*   - USING CORRECT DIPOL MOMENT I.E.,DIFFERENT COMMON/MODEL/ */
/*   - USING IGRF EARTH MAGNETIC FIELD MODELS FROM 1945 TO 1990 */
/* -------------------------------------------------------------------- */
/*  INPUT:  ENTRY POINT SHELLG */
/*               GLAT  GEODETIC LATITUDE IN DEGREES (NORTH) */
/*               GLON  GEODETIC LONGITUDE IN DEGREES (EAST) */
/*               ALT   ALTITUDE IN KM ABOVE SEA LEVEL */

/*          ENTRY POINT SHELLC */
/*               V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM) */
/*                       X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE */
/*                       Y-AXIS POINTING TO EQUATOR AT 90 LONG. */
/*                       Z-AXIS POINTING TO NORTH POLE */

/*          DIMO       DIPOL MOMENT IN GAUSS (NORMALIZED TO EARTH RADIUS) */

/*          COMMON */
/*               X(3)    NOT USED */
/*               H(144)  FIELD MODEL COEFFICIENTS ADJUSTED FOR SHELLG */
/* ----------------------------------------------------------------------- */
/*  OUTPUT: FL           L-VALUE */
/*          ICODE        =1 NORMAL COMPLETION */
/*                       =2 UNPHYSICAL CONJUGATE POINT (FL MEANINGLESS) */
/*                       =3 SHELL PARAMETER GREATER THAN LIMIT UP TO */
/*                          WHICH ACCURATE CALCULATION IS REQUIRED; */
/*                          APPROXIMATION IS USED. */
/*          B0           MAGNETIC FIELD STRENGTH IN GAUSS */
/* ----------------------------------------------------------------------- */
/* f2py intent(in) glat */
/* f2py intent(in) glon */
/* f2py intent(in) alt */
/* f2py intent(out) fl */
/* f2py intent(out) icode */
/* f2py intent(out) b0 */

/* -- RMIN, RMAX ARE BOUNDARIES FOR IDENTIFICATION OF ICODE=2 AND 3 */
/* -- STEP IS STEP SIZE FOR FIELD LINE TRACING */
/* -- STEQ IS STEP SIZE FOR INTEGRATION */

    /* Parameter adjustments */
    if (v) {
	--v;
	}

    /* Function Body */
    switch(n__) {
	case 1: goto L_shellc;
	}

    bequ = (float)1e10;
/* *****ENTRY POINT  SHELLG  TO BE USED WITH GEODETIC CO-ORDINATES */
    rlat = *glat * gener_1.umr;
    ct = sin(rlat);
    st = cos(rlat);
    d__ = sqrt(gener_1.aquad - (gener_1.aquad - gener_1.bquad) * ct * ct);
    _BLNK__1.x[0] = (*alt + gener_1.aquad / d__) * st / gener_1.era;
    _BLNK__1.x[2] = (*alt + gener_1.bquad / d__) * ct / gener_1.era;
    rlon = *glon * gener_1.umr;
    _BLNK__1.x[1] = _BLNK__1.x[0] * sin(rlon);
    _BLNK__1.x[0] *= cos(rlon);
    goto L9;

L_shellc:
/* *****ENTRY POINT  SHELLC  TO BE USED WITH CARTESIAN CO-ORDINATES */
    _BLNK__1.x[0] = v[1];
    _BLNK__1.x[1] = v[2];
    _BLNK__1.x[2] = v[3];
/* *****CONVERT TO DIPOL-ORIENTED CO-ORDINATES */
L9:
    rq = (float)1. / (_BLNK__1.x[0] * _BLNK__1.x[0] + _BLNK__1.x[1] * 
	    _BLNK__1.x[1] + _BLNK__1.x[2] * _BLNK__1.x[2]);
    r3h = sqrt(rq * sqrt(rq));
    p[8] = (_BLNK__1.x[0] * u[0] + _BLNK__1.x[1] * u[1] + _BLNK__1.x[2] * u[2]
	    ) * r3h;
    p[9] = (_BLNK__1.x[0] * u[3] + _BLNK__1.x[1] * u[4]) * r3h;
    p[10] = (_BLNK__1.x[0] * u[6] + _BLNK__1.x[1] * u[7] + _BLNK__1.x[2] * u[
	    8]) * rq;
/* *****FIRST THREE POINTS OF FIELD LINE */
    step = -signc_(&step, &p[10]);
    stoer_(&p[8], &bq2, &r2);
    *b0 = sqrt(bq2);
    p[16] = p[8] + step * (float).5 * p[11];
    p[17] = p[9] + step * (float).5 * p[12];
    p[18] = p[10] + step * (float).5;
    stoer_(&p[16], &bq3, &r3);
    p[0] = p[8] - step * (p[11] * (float)2. - p[19]);
    p[1] = p[9] - step * (p[12] * (float)2. - p[20]);
    p[2] = p[10] - step;
    stoer_(p, &bq1, &r1);
    p[16] = p[8] + step * (p[19] * (float)20. - p[11] * (float)3. + p[3]) / (
	    float)18.;
    p[17] = p[9] + step * (p[20] * (float)20. - p[12] * (float)3. + p[4]) / (
	    float)18.;
    p[18] = p[10] + step;
    stoer_(&p[16], &bq3, &r3);
/* *****INVERT SENSE IF REQUIRED */
    if (bq3 <= bq1) {
	goto L2;
    }
    step = -step;
    r3 = r1;
    bq3 = bq1;
    for (i__ = 1; i__ <= 7; ++i__) {
	zz = p[i__ - 1];
	p[i__ - 1] = p[i__ + 15];
/* L1: */
	p[i__ + 15] = zz;
    }
/* *****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH */
L2:
    if (bq1 < bequ) {
	bequ = bq1;
	iequ = 1;
    }
    if (bq2 < bequ) {
	bequ = bq2;
	iequ = 2;
    }
    if (bq3 < bequ) {
	bequ = bq3;
	iequ = 3;
    }
/* *****INITIALIZATION OF INTEGRATION LOOPS */
    step12 = step / (float)12.;
    step2 = step + step;
    steq = signc_(&steq, &step);
    fi = (float)0.;
    *icode = 1;
    oradik = (float)0.;
    oterm = (float)0.;
    stp = r2 * steq;
    z__ = p[10] + stp;
    stp /= (float).75;
    p[7] = step2 * (p[0] * p[3] + p[1] * p[4]);
    p[15] = step2 * (p[8] * p[11] + p[9] * p[12]);
/* *****MAIN LOOP (FIELD LINE TRACING) */
    for (n = 3; n <= 3333; ++n) {
/* *****CORRECTOR (FIELD LINE TRACING) */
	p[(n << 3) - 8] = p[(n - 1 << 3) - 8] + step12 * (p[(n << 3) - 5] * (
		float)5. + p[(n - 1 << 3) - 5] * (float)8. - p[(n - 2 << 3) - 
		5]);
	p[(n << 3) - 7] = p[(n - 1 << 3) - 7] + step12 * (p[(n << 3) - 4] * (
		float)5. + p[(n - 1 << 3) - 4] * (float)8. - p[(n - 2 << 3) - 
		4]);
/* *****PREPARE EXPANSION COEFFICIENTS FOR INTERPOLATION */
/* *****OF SLOWLY VARYING QUANTITIES */
	p[(n << 3) - 1] = step2 * (p[(n << 3) - 8] * p[(n << 3) - 5] + p[(n <<
		 3) - 7] * p[(n << 3) - 4]);
/* Computing 2nd power */
	r__1 = p[(n - 1 << 3) - 8];
/* Computing 2nd power */
	r__2 = p[(n - 1 << 3) - 7];
	c0 = r__1 * r__1 + r__2 * r__2;
	c1 = p[(n - 1 << 3) - 1];
	c2 = (p[(n << 3) - 1] - p[(n - 2 << 3) - 1]) * (float).25;
	c3 = (p[(n << 3) - 1] + p[(n - 2 << 3) - 1] - c1 - c1) / (float)6.;
	d0 = p[(n - 1 << 3) - 3];
	d1 = (p[(n << 3) - 3] - p[(n - 2 << 3) - 3]) * (float).5;
	d2 = (p[(n << 3) - 3] + p[(n - 2 << 3) - 3] - d0 - d0) * (float).5;
	e0 = p[(n - 1 << 3) - 2];
	e1 = (p[(n << 3) - 2] - p[(n - 2 << 3) - 2]) * (float).5;
	e2 = (p[(n << 3) - 2] + p[(n - 2 << 3) - 2] - e0 - e0) * (float).5;
/* *****INNER LOOP (FOR QUADRATURE) */
L4:
	t = (z__ - p[(n - 1 << 3) - 6]) / step;
	if (t > (float)1.) {
	    goto L5;
	}
	hli = (((c3 * t + c2) * t + c1) * t + c0) * (float).5;
	zq = z__ * z__;
	r__ = hli + sqrt(hli * hli + zq);
	if (r__ <= rmin) {
	    goto L30;
	}
	rq = r__ * r__;
	ff = sqrt(zq * (float)3. / rq + (float)1.);
	radik = *b0 - ((d2 * t + d1) * t + d0) * r__ * rq * ff;
	if (r__ - rmax <= (float)0.) {
	    goto L44;
	} else {
	    goto L45;
	}
L45:
	*icode = 2;
/* Computing 2nd power */
	r__1 = r__ - rmax;
	radik -= r__1 * r__1 * (float)12.;
L44:
	if (radik + radik <= oradik) {
	    goto L10;
	}
	term = sqrt(radik) * ff * ((e2 * t + e1) * t + e0) / (rq + zq);
	fi += stp * (oterm + term);
	oradik = radik;
	oterm = term;
	stp = r__ * steq;
	z__ += stp;
	goto L4;
/* *****PREDICTOR (FIELD LINE TRACING) */
L5:
	p[(n + 1 << 3) - 8] = p[(n << 3) - 8] + step12 * (p[(n << 3) - 5] * (
		float)23. - p[(n - 1 << 3) - 5] * (float)16. + p[(n - 2 << 3) 
		- 5] * (float)5.);
	p[(n + 1 << 3) - 7] = p[(n << 3) - 7] + step12 * (p[(n << 3) - 4] * (
		float)23. - p[(n - 1 << 3) - 4] * (float)16. + p[(n - 2 << 3) 
		- 4] * (float)5.);
	p[(n + 1 << 3) - 6] = p[(n << 3) - 6] + step;
	stoer_(&p[(n + 1 << 3) - 8], &bq3, &r3);
/* *****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH */
	if (bq3 < bequ) {
	    iequ = n + 1;
	    bequ = bq3;
	}
/* L3: */
    }
L10:
    if (iequ < 2) {
	iequ = 2;
    }
    fidb0_1.sp[0] = p[(iequ - 1 << 3) - 8];
    fidb0_1.sp[1] = p[(iequ - 1 << 3) - 7];
    fidb0_1.sp[2] = p[(iequ - 1 << 3) - 6];
    if (oradik < (float)1e-15) {
	goto L11;
    }
    fi += stp / (float).75 * oterm * oradik / (oradik - radik);

/* -- The minimal allowable value of FI was changed from 1E-15 to 1E-12, */
/* -- because 1E-38 is the minimal allowable arg. for ALOG in our envir. */
/* -- D. Bilitza, Nov 87. */

L11:
    fi = dabs(fi) * (float).5 / sqrt(*b0) + (float)1e-12;

/* *****COMPUTE L FROM B AND I.  SAME AS CARMEL IN INVAR. */

/* -- Correct dipole moment is used here. D. Bilitza, Nov 87. */

    dimob0 = *dimo / *b0;
    arg1 = log(fi);
    arg2 = log(dimob0);
/* 	arg = FI*FI*FI/DIMOB0 */
/* 	if(abs(arg).gt.88.0) arg=88.0 */
    xx = arg1 * 3 - arg2;
    if (xx > (float)23.) {
	goto L776;
    }
    if (xx > (float)11.7) {
	goto L775;
    }
    if (xx > (float)3.) {
	goto L774;
    }
    if (xx > (float)-3.) {
	goto L773;
    }
    if (xx > (float)-22.) {
	goto L772;
    }
/* L771: */
    gg = xx * (float).333338 + (float).30062102;
    goto L777;
L772:
    gg = ((((((((xx * (float)-8.1537735e-14 + (float)8.3232531e-13) * xx + (
	    float)1.0066362e-9) * xx + (float)8.1048663e-8) * xx + (float)
	    3.2916354e-6) * xx + (float)8.2711096e-5) * xx + (float)
	    .0013714667) * xx + (float).015017245) * xx + (float).43432642) * 
	    xx + (float).62337691;
    goto L777;
L773:
    gg = ((((((((xx * (float)2.6047023e-10 + (float)2.3028767e-9) * xx - (
	    float)2.1997983e-8) * xx - (float)5.3977642e-7) * xx - (float)
	    3.3408822e-6) * xx + (float)3.8379917e-5) * xx + (float)
	    .0011784234) * xx + (float).014492441) * xx + (float).43352788) * 
	    xx + (float).6228644;
    goto L777;
L774:
    gg = ((((((((xx * (float)6.3271665e-10 - (float)3.958306e-8) * xx + (
	    float)9.9766148e-7) * xx - (float)1.2531932e-5) * xx + (float)
	    7.9451313e-5) * xx - (float)3.2077032e-4) * xx + (float)
	    .0021680398) * xx + (float).012817956) * xx + (float).43510529) * 
	    xx + (float).6222355;
    goto L777;
L775:
    gg = (((((xx * (float)2.8212095e-8 - (float)3.8049276e-6) * xx + (float)
	    2.170224e-4) * xx - (float).0067310339) * xx + (float).12038224) *
	     xx - (float).18461796) * xx + (float)2.0007187;
    goto L777;
L776:
    gg = xx - (float)3.0460681;
L777:
    *fl = exp(log((exp(gg) + (float)1.) * dimob0) / (float)3.);
    return 0;
/* *****APPROXIMATION FOR HIGH VALUES OF L. */
L30:
    *icode = 3;
    t = -p[(n - 1 << 3) - 6] / step;
    *fl = (float)1. / ((r__1 = ((c3 * t + c2) * t + c1) * t + c0, dabs(r__1)) 
	    + (float)1e-15);
    return 0;
} /* shellg_ */

/* Subroutine */ int shellg_(real *glat, real *glon, real *alt, real *dimo, 
	real *fl, integer *icode, real *b0)
{
    return shellg_0_(0, glat, glon, alt, dimo, fl, icode, b0, (real *)0);
    }

/* Subroutine */ int shellc_(real *v, real *fl, real *b0)
{
    return shellg_0_(1, (real *)0, (real *)0, (real *)0, (real *)0, fl, (
	    integer *)0, b0, v);
    }



/* Subroutine */ int stoer_(real *p, real *bq, real *r__)
{
    /* Initialized data */

    static real u[9]	/* was [3][3] */ = { (float).3511737,(float)-.9148385,
	    (float)-.1993679,(float).9335804,(float).358368,(float)0.,(float)
	    .0714471,(float)-.186126,(float).9799247 };

    /* System generated locals */
    real r__1;

    /* Local variables */
    static real q, dr, dx, dy, dz, rq, xm, ym, zm, wr, fli, dsq, dxm, dym, 
	    dzm;
    int feldi_();

/* ******************************************************************* */
/* * SUBROUTINE USED FOR FIELD LINE TRACING IN SHELLG                * */
/* * CALLS ENTRY POINT FELDI IN GEOMAGNETIC FIELD SUBROUTINE FELDG   * */
/* ******************************************************************* */
/* *****XM,YM,ZM  ARE GEOMAGNETIC CARTESIAN INVERSE CO-ORDINATES */
    /* Parameter adjustments */
    --p;

    /* Function Body */
    zm = p[3];
    fli = p[1] * p[1] + p[2] * p[2] + (float)1e-15;
/* Computing 2nd power */
    r__1 = zm + zm;
    *r__ = (fli + sqrt(fli * fli + r__1 * r__1)) * (float).5;
    rq = *r__ * *r__;
    wr = sqrt(*r__);
    xm = p[1] * wr;
    ym = p[2] * wr;
/* *****TRANSFORM TO GEOGRAPHIC CO-ORDINATE SYSTEM */
    _BLNK__2.xi[0] = xm * u[0] + ym * u[3] + zm * u[6];
    _BLNK__2.xi[1] = xm * u[1] + ym * u[4] + zm * u[7];
    _BLNK__2.xi[2] = xm * u[2] + zm * u[8];
/* *****COMPUTE DERIVATIVES */
/*      CALL FELDI(XI,H) */
    feldi_();
    q = _BLNK__2.h__[0] / rq;
    dx = _BLNK__2.h__[2] + _BLNK__2.h__[2] + q * _BLNK__2.xi[0];
    dy = _BLNK__2.h__[3] + _BLNK__2.h__[3] + q * _BLNK__2.xi[1];
    dz = _BLNK__2.h__[1] + _BLNK__2.h__[1] + q * _BLNK__2.xi[2];
/* *****TRANSFORM BACK TO GEOMAGNETIC CO-ORDINATE SYSTEM */
    dxm = u[0] * dx + u[1] * dy + u[2] * dz;
    dym = u[3] * dx + u[4] * dy;
    dzm = u[6] * dx + u[7] * dy + u[8] * dz;
    dr = (xm * dxm + ym * dym + zm * dzm) / *r__;
/* *****FORM SLOWLY VARYING EXPRESSIONS */
    p[4] = (wr * dxm - p[1] * (float).5 * dr) / (*r__ * dzm);
    p[5] = (wr * dym - p[2] * (float).5 * dr) / (*r__ * dzm);
    dsq = rq * (dxm * dxm + dym * dym + dzm * dzm);
    *bq = dsq * rq * rq;
    p[6] = sqrt(dsq / (rq + zm * (float)3. * zm));
    p[7] = p[6] * (rq + zm * zm) / (rq * dzm);
    return 0;
} /* stoer_ */



/* Subroutine */ int feldg_0_(int n__, real *glat, real *glon, real *alt, 
	real *bnorth, real *beast, real *bdown, real *babs, real *v, real *b)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real d__, f;
    static integer i__, k, m;
    static real s, t, x, y, z__;
    static integer ih;
    static real cp;
    static integer il;
    static real ct;
    static integer is;
    static real sp, rq, st, rho, xxx, yyy, zzz, brho;
    static integer imax;
    static real rlat;
    static integer last;
    static real rlon, bxxx, byyy, bzzz;
    static integer ihmax;

/* ------------------------------------------------------------------- */
/* CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL */
/* REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61, */
/*      1970. */
/* -------------------------------------------------------------------- */
/* CHANGES (D. BILITZA, NOV 87): */
/*   - FIELD COEFFICIENTS IN BINARY DATA FILES INSTEAD OF BLOCK DATA */
/*   - CALCULATES DIPOL MOMENT */
/* -------------------------------------------------------------------- */
/*  INPUT:  ENTRY POINT FELDG */
/*               GLAT  GEODETIC LATITUDE IN DEGREES (NORTH) */
/*               GLON  GEODETIC LONGITUDE IN DEGREES (EAST) */
/*               ALT   ALTITUDE IN KM ABOVE SEA LEVEL */

/*          ENTRY POINT FELDC */
/*               V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM) */
/*                       X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE */
/*                       Y-AXIS POINTING TO EQUATOR AT 90 LONG. */
/*                       Z-AXIS POINTING TO NORTH POLE */

/*          COMMON BLANK AND ENTRY POINT FELDI ARE NEEDED WHEN USED */
/*            IN CONNECTION WITH L-CALCULATION PROGRAM SHELLG. */

/*          COMMON /MODEL/ AND /GENER/ */
/*               UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT> */
/*               ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN */
/*                       COORDINATES (6371.2 KM) */
/*               AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS FOR */
/*                       EARTH ELLIPSOID AS RECOMMENDED BY INTERNATIONAL */
/*                       ASTRONOMICAL UNION (6378.160, 6356.775 KM). */
/*               NMAX    MAXIMUM ORDER OF SPHERICAL HARMONICS */
/*               TIME    YEAR (DECIMAL: 1973.5) FOR WHICH MAGNETIC */
/*                       FIELD IS TO BE CALCULATED */
/*               G(M)    NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF) */
/*                       M=NMAX*(NMAX+2) */
/* ------------------------------------------------------------------------ */
/*  OUTPUT: BABS   MAGNETIC FIELD STRENGTH IN GAUSS */
/*          BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT */
/*                 TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS */
/*                 POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST */
/*                 AND DOWNWARD. */
/* ----------------------------------------------------------------------- */
/* f2py intent(in) glat */
/* f2py intent(in) glon */
/* f2py intent(in) alt */
/* f2py intent(out) babs */
/* f2py intent(out) bnorth */
/* f2py intent(out) beast */
/* f2py intent(out) bdown */

/* -- IS RECORDS ENTRY POINT */

/* *****ENTRY POINT  FELDG  TO BE USED WITH GEODETIC CO-ORDINATES */
    /* Parameter adjustments */
    if (v) {
	--v;
	}
    if (b) {
	--b;
	}

    /* Function Body */
    switch(n__) {
	case 1: goto L_feldc;
	case 2: goto L_feldi;
	}

    is = 1;
    rlat = *glat * gener_1.umr;
    ct = sin(rlat);
    st = cos(rlat);
    d__ = sqrt(gener_1.aquad - (gener_1.aquad - gener_1.bquad) * ct * ct);
    rlon = *glon * gener_1.umr;
    cp = cos(rlon);
    sp = sin(rlon);
    zzz = (*alt + gener_1.bquad / d__) * ct / gener_1.era;
    rho = (*alt + gener_1.aquad / d__) * st / gener_1.era;
    xxx = rho * cp;
    yyy = rho * sp;
    goto L10;

L_feldc:
/* *****ENTRY POINT  FELDC  TO BE USED WITH CARTESIAN CO-ORDINATES */
    is = 2;
    xxx = v[1];
    yyy = v[2];
    zzz = v[3];
L10:
    rq = (float)1. / (xxx * xxx + yyy * yyy + zzz * zzz);
    _BLNK__2.xi[0] = xxx * rq;
    _BLNK__2.xi[1] = yyy * rq;
    _BLNK__2.xi[2] = zzz * rq;
    goto L20;

L_feldi:
/* *****ENTRY POINT  FELDI  USED FOR L COMPUTATION */
    is = 3;
L20:
    ihmax = model_1.nmax * model_1.nmax + 1;
    last = ihmax + model_1.nmax + model_1.nmax;
    imax = model_1.nmax + model_1.nmax - 1;
    i__1 = last;
    for (i__ = ihmax; i__ <= i__1; ++i__) {
/* L8: */
	_BLNK__2.h__[i__ - 1] = model_1.g[i__ - 1];
    }
    for (k = 1; k <= 3; k += 2) {
	i__ = imax;
	ih = ihmax;
L1:
	il = ih - i__;
	f = (float)2. / (real) (i__ - k + 2);
	x = _BLNK__2.xi[0] * f;
	y = _BLNK__2.xi[1] * f;
	z__ = _BLNK__2.xi[2] * (f + f);
	i__ += -2;
	if ((i__1 = i__ - 1) < 0) {
	    goto L5;
	} else if (i__1 == 0) {
	    goto L4;
	} else {
	    goto L2;
	}
L2:
	i__1 = i__;
	for (m = 3; m <= i__1; m += 2) {
	    _BLNK__2.h__[il + m] = model_1.g[il + m] + z__ * _BLNK__2.h__[ih 
		    + m] + x * (_BLNK__2.h__[ih + m + 2] - _BLNK__2.h__[ih + 
		    m - 2]) - y * (_BLNK__2.h__[ih + m + 1] + _BLNK__2.h__[ih 
		    + m - 3]);
/* L3: */
	    _BLNK__2.h__[il + m - 1] = model_1.g[il + m - 1] + z__ * 
		    _BLNK__2.h__[ih + m - 1] + x * (_BLNK__2.h__[ih + m + 1] 
		    - _BLNK__2.h__[ih + m - 3]) + y * (_BLNK__2.h__[ih + m + 
		    2] + _BLNK__2.h__[ih + m - 2]);
	}
L4:
	_BLNK__2.h__[il + 1] = model_1.g[il + 1] + z__ * _BLNK__2.h__[ih + 1] 
		+ x * _BLNK__2.h__[ih + 3] - y * (_BLNK__2.h__[ih + 2] + 
		_BLNK__2.h__[ih - 1]);
	_BLNK__2.h__[il] = model_1.g[il] + z__ * _BLNK__2.h__[ih] + y * 
		_BLNK__2.h__[ih + 3] + x * (_BLNK__2.h__[ih + 2] - 
		_BLNK__2.h__[ih - 1]);
L5:
	_BLNK__2.h__[il - 1] = model_1.g[il - 1] + z__ * _BLNK__2.h__[ih - 1] 
		+ (x * _BLNK__2.h__[ih] + y * _BLNK__2.h__[ih + 1]) * (float)
		2.;
	ih = il;
	if (i__ >= k) {
	    goto L1;
	}
/* L6: */
    }
    if (is == 3) {
	return 0;
    }
    s = _BLNK__2.h__[0] * (float).5 + (_BLNK__2.h__[1] * _BLNK__2.xi[2] + 
	    _BLNK__2.h__[2] * _BLNK__2.xi[0] + _BLNK__2.h__[3] * _BLNK__2.xi[
	    1]) * (float)2.;
    t = (rq + rq) * sqrt(rq);
    bxxx = t * (_BLNK__2.h__[2] - s * xxx);
    byyy = t * (_BLNK__2.h__[3] - s * yyy);
    bzzz = t * (_BLNK__2.h__[1] - s * zzz);
    if (is == 2) {
	goto L7;
    }
    *babs = sqrt(bxxx * bxxx + byyy * byyy + bzzz * bzzz);
    *beast = byyy * cp - bxxx * sp;
    brho = byyy * sp + bxxx * cp;
    *bnorth = bzzz * st - brho * ct;
    *bdown = -bzzz * ct - brho * st;
    return 0;
L7:
    b[1] = bxxx;
    b[2] = byyy;
    b[3] = bzzz;
    return 0;
} /* feldg_ */

/* Subroutine */ int feldg_(real *glat, real *glon, real *alt, real *bnorth, 
	real *beast, real *bdown, real *babs)
{
    return feldg_0_(0, glat, glon, alt, bnorth, beast, bdown, babs, (real *)0,
	     (real *)0);
    }

/* Subroutine */ int feldc_(real *v, real *b)
{
    return feldg_0_(1, (real *)0, (real *)0, (real *)0, (real *)0, (real *)0, 
	    (real *)0, (real *)0, v, b);
    }

/* Subroutine */ int feldi_()
{
    return feldg_0_(2, (real *)0, (real *)0, (real *)0, (real *)0, (real *)0, 
	    (real *)0, (real *)0, (real *)0, (real *)0);
    }



/* Subroutine */ int feldcof_(real *year, real *dimo)
{
    /* Initialized data */

    
    static real dtemod[5] = { 1990.,1995.,2e3,
	    2005.,2010. };

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    int intershc_(real *, real *, integer *, real *, 
	    real *, integer *, real *, integer *, real *), extrashc_(real *, 
	    real *, integer *, real *, integer *, real *, integer *, real *);
    static doublereal f;
    static integer i__, j, l, m, n;
    static doublereal x, f0;
    static integer is, iu;
    static real gh2[144], gha[144];
    static integer ier;
    static real dte1, dte2;
    static integer iyea, nmax1, nmax2;
    static real sqrt2;
    static integer numye;
    int getshc_(real *, integer *, real *, real *, 
	    integer *);

/* ------------------------------------------------------------------------ */
/*  DETERMINES COEFFICIENTS AND DIPOL MOMENT FROM IGRF MODELS */

/*       INPUT:  YEAR    DECIMAL YEAR FOR WHICH GEOMAGNETIC FIELD IS TO */
/*                       BE CALCULATED */
/*       OUTPUT: DIMO    GEOMAGNETIC DIPOL MOMENT IN GAUSS (NORMALIZED */
/*                       TO EARTH'S RADIUS) AT THE TIME (YEAR) */
/*  D. BILITZA, NSSDC, GSFC, CODE 633, GREENBELT, MD 20771, */
/*       (301)286-9536   NOV 1987. */
/*  ### updated to IGRF-2000 version -dkb- 5/31/2000 */
/*  ### updated to IGRF-2005 version -dkb- 3/24/2005 */
/* ----------------------------------------------------------------------- */
/* f2py intent(in) year */
/* f2py intent(out) dimo */
/* ### FILMOD, DTEMOD arrays +1 */

/* ### numye = numye + 1 ; is number of years represented by IGRF */

    numye = 4;

/*  IS=0 FOR SCHMIDT NORMALIZATION   IS=1 GAUSS NORMALIZATION */
/*  IU  IS INPUT UNIT NUMBER FOR IGRF COEFFICIENT SETS */

    iu = 10;
    is = 0;
/* -- DETERMINE IGRF-YEARS FOR INPUT-YEAR */
    model_2.time = *year;
    iyea = (integer) (*year / (float)5.) * 5;
    l = (iyea - 1990) / 5 + 1;
    if (l < 1) {
	l = 1;
    }
    if (l > numye) {
	l = numye;
    }
    dte1 = dtemod[l - 1];
    dte2 = dtemod[l];

/* -- GET IGRF COEFFICIENTS FOR THE BOUNDARY YEARS */
    if (l == 1) {
	getshc_(dgrf90, &nmax1, &gener_2.erad, model_2.gh1, &ier);
	getshc_(dgrf95, &nmax2, &gener_2.erad, gh2, &ier);
    }
    if (l == 2) {
	getshc_(dgrf95, &nmax1, &gener_2.erad, model_2.gh1, &ier);
	getshc_(dgrf00, &nmax2, &gener_2.erad, gh2, &ier);
    }
    if (l == 3) {
	getshc_(dgrf00, &nmax1, &gener_2.erad, model_2.gh1, &ier);
	getshc_(igrf05, &nmax2, &gener_2.erad, gh2, &ier);
    }
    if (l == 4) {
	getshc_(igrf05, &nmax1, &gener_2.erad, model_2.gh1, &ier);
	getshc_(igrf05s, &nmax2, &gener_2.erad, gh2, &ier);
    }
    if (ier != 0) {
	goto L9998;
    }
/* -- DETERMINE IGRF COEFFICIENTS FOR YEAR */
    if (l <= numye - 1) {
	intershc_(year, &dte1, &nmax1, model_2.gh1, &dte2, &nmax2, gh2, &
		model_2.nmax, gha);
    } else {
	extrashc_(year, &dte1, &nmax1, model_2.gh1, &nmax2, gh2, &
		model_2.nmax, gha);
    }
/* -- DETERMINE MAGNETIC DIPOL MOMENT AND COEFFIECIENTS G */
    f0 = 0.;
    for (j = 1; j <= 3; ++j) {
	f = gha[j - 1] * 1e-5;
	f0 += f * f;
/* L1234: */
    }
    *dimo = sqrt(f0);
    model_2.gh1[0] = (float)0.;
    i__ = 2;
    f0 = 1e-5;
    if (is == 0) {
	f0 = -f0;
    }
    sqrt2 = sqrt((float)2.);
    i__1 = model_2.nmax;
    for (n = 1; n <= i__1; ++n) {
	x = (doublereal) n;
	f0 = f0 * x * x / (x * 4. - 2.);
	if (is == 0) {
	    f0 = f0 * (x * 2. - 1.) / x;
	}
	f = f0 * .5;
	if (is == 0) {
	    f *= sqrt2;
	}
	model_2.gh1[i__ - 1] = gha[i__ - 2] * f0;
	++i__;
	i__2 = n;
	for (m = 1; m <= i__2; ++m) {
	    f = f * (x + m) / (x - m + 1.);
	    if (is == 0) {
		f *= sqrt((x - m + 1.) / (x + m));
	    }
	    model_2.gh1[i__ - 1] = gha[i__ - 2] * f;
	    model_2.gh1[i__] = gha[i__ - 1] * f;
	    i__ += 2;
/* L9: */
	}
    }
    return 0;
L9998:
    *dimo = (float)-9.999e102;
    *year = (float)-1.;
    return 0;
} /* feldcof_ */



/* Subroutine */ int getshc_(real *multip, integer *nmax, real *erad, real *
	gh, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static real g, h__;
    static integer i__, m, n, mm;
    static real cx;
    static integer nn;

/* =============================================================== */

/*       Version 1.01 */

/*       Reads spherical harmonic coefficients from the specified */
/*       file into an array. */

/*       Input: */
/*           IU    - Logical unit number */
/*           FSPEC - File specification */

/*       Output: */
/*           NMAX  - Maximum degree and order of model */
/*           ERAD  - Earth's radius associated with the spherical */
/*                   harmonic coefficients, in the same units as */
/*                   elevation */
/*           GH    - Schmidt quasi-normal internal spherical */
/*                   harmonic coefficients */
/*           IER   - Error number: =  0, no error */
/*                                 = -2, records out of order */
/*                                 = FORTRAN run-time error number */

/*       A. Zunde */
/*       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225 */

/* =============================================================== */
/*      INTEGER CX */
    /* Parameter adjustments */
    --gh;
    --multip;

    /* Function Body */
    *nmax = (integer) multip[1];
    *erad = multip[2];
/* --------------------------------------------------------------- */
/*       Read the coefficient file, arranged as follows: */

/*                                       N     M     G     H */
/*                                       ---------------------- */
/*                                   /   1     0    GH(1)  - */
/*                                  /    1     1    GH(2) GH(3) */
/*                                 /     2     0    GH(4)  - */
/*                                /      2     1    GH(5) GH(6) */
/*           NMAX*(NMAX+3)/2     /       2     2    GH(7) GH(8) */
/*              records          \       3     0    GH(9)  - */
/*                                \      .     .     .     . */
/*                                 \     .     .     .     . */
/*           NMAX*(NMAX+2)          \    .     .     .     . */
/*           elements in GH          \  NMAX  NMAX   .     . */

/*       N and M are, respectively, the degree and order of the */
/*       coefficient. */
/* --------------------------------------------------------------- */
    i__ = 0;
    cx = (float)4.;
    i__1 = *nmax;
    for (nn = 1; nn <= i__1; ++nn) {
	i__2 = nn;
	for (mm = 0; mm <= i__2; ++mm) {
	    n = (integer) multip[(integer) cx];
	    m = (integer) multip[(integer) (cx + 1)];
	    g = multip[(integer) (cx + 2)];
	    h__ = multip[(integer) (cx + 3)];
	    cx += 4;
	    if (nn != n || mm != m) {
		*ier = -2;
		goto L999;
	    }
	    ++i__;
	    gh[i__] = g;
	    if (m != 0) {
		++i__;
		gh[i__] = h__;
	    }
/* L2233: */
	}
/* L2211: */
    }
L999:
    return 0;
} /* getshc_ */



/* Subroutine */ int intershc_(real *date, real *dte1, integer *nmax1, real *
	gh1, real *dte2, integer *nmax2, real *gh2, integer *nmax, real *gh)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, l;
    static real factor;

/* =============================================================== */

/*       Version 1.01 */

/*       Interpolates linearly, in time, between two spherical */
/*       harmonic models. */

/*       Input: */
/*           DATE  - Date of resulting model (in decimal year) */
/*           DTE1  - Date of earlier model */
/*           NMAX1 - Maximum degree and order of earlier model */
/*           GH1   - Schmidt quasi-normal internal spherical */
/*                   harmonic coefficients of earlier model */
/*           DTE2  - Date of later model */
/*           NMAX2 - Maximum degree and order of later model */
/*           GH2   - Schmidt quasi-normal internal spherical */
/*                   harmonic coefficients of later model */

/*       Output: */
/*           GH    - Coefficients of resulting model */
/*           NMAX  - Maximum degree and order of resulting model */

/*       A. Zunde */
/*       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225 */

/* =============================================================== */
/* --------------------------------------------------------------- */
/*       The coefficients (GH) of the resulting model, at date */
/*       DATE, are computed by linearly interpolating between the */
/*       coefficients of the earlier model (GH1), at date DTE1, */
/*       and those of the later model (GH2), at date DTE2. If one */
/*       model is smaller than the other, the interpolation is */
/*       performed with the missing coefficients assumed to be 0. */
/* --------------------------------------------------------------- */
    /* Parameter adjustments */
    --gh;
    --gh2;
    --gh1;

    /* Function Body */
    factor = (*date - *dte1) / (*dte2 - *dte1);
    if (*nmax1 == *nmax2) {
	k = *nmax1 * (*nmax1 + 2);
	*nmax = *nmax1;
    } else if (*nmax1 > *nmax2) {
	k = *nmax2 * (*nmax2 + 2);
	l = *nmax1 * (*nmax1 + 2);
	i__1 = l;
	for (i__ = k + 1; i__ <= i__1; ++i__) {
/* L1122: */
	    gh[i__] = gh1[i__] + factor * (-gh1[i__]);
	}
	*nmax = *nmax1;
    } else {
	k = *nmax1 * (*nmax1 + 2);
	l = *nmax2 * (*nmax2 + 2);
	i__1 = l;
	for (i__ = k + 1; i__ <= i__1; ++i__) {
/* L1133: */
	    gh[i__] = factor * gh2[i__];
	}
	*nmax = *nmax2;
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1144: */
	gh[i__] = gh1[i__] + factor * (gh2[i__] - gh1[i__]);
    }
    return 0;
} /* intershc_ */



/* Subroutine */ int extrashc_(real *date, real *dte1, integer *nmax1, real *
	gh1, integer *nmax2, real *gh2, integer *nmax, real *gh)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, l;
    static real factor;

/* =============================================================== */

/*       Version 1.01 */

/*       Extrapolates linearly a spherical harmonic model with a */
/*       rate-of-change model. */

/*       Input: */
/*           DATE  - Date of resulting model (in decimal year) */
/*           DTE1  - Date of base model */
/*           NMAX1 - Maximum degree and order of base model */
/*           GH1   - Schmidt quasi-normal internal spherical */
/*                   harmonic coefficients of base model */
/*           NMAX2 - Maximum degree and order of rate-of-change */
/*                   model */
/*           GH2   - Schmidt quasi-normal internal spherical */
/*                   harmonic coefficients of rate-of-change model */

/*       Output: */
/*           GH    - Coefficients of resulting model */
/*           NMAX  - Maximum degree and order of resulting model */

/*       A. Zunde */
/*       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225 */

/* =============================================================== */
/* --------------------------------------------------------------- */
/*       The coefficients (GH) of the resulting model, at date */
/*       DATE, are computed by linearly extrapolating the coef- */
/*       ficients of the base model (GH1), at date DTE1, using */
/*       those of the rate-of-change model (GH2), at date DTE2. If */
/*       one model is smaller than the other, the extrapolation is */
/*       performed with the missing coefficients assumed to be 0. */
/* --------------------------------------------------------------- */
    /* Parameter adjustments */
    --gh;
    --gh2;
    --gh1;

    /* Function Body */
    factor = *date - *dte1;
    if (*nmax1 == *nmax2) {
	k = *nmax1 * (*nmax1 + 2);
	*nmax = *nmax1;
    } else if (*nmax1 > *nmax2) {
	k = *nmax2 * (*nmax2 + 2);
	l = *nmax1 * (*nmax1 + 2);
	i__1 = l;
	for (i__ = k + 1; i__ <= i__1; ++i__) {
/* L1155: */
	    gh[i__] = gh1[i__];
	}
	*nmax = *nmax1;
    } else {
	k = *nmax1 * (*nmax1 + 2);
	l = *nmax2 * (*nmax2 + 2);
	i__1 = l;
	for (i__ = k + 1; i__ <= i__1; ++i__) {
/* L1166: */
	    gh[i__] = factor * gh2[i__];
	}
	*nmax = *nmax2;
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1177: */
	gh[i__] = gh1[i__] + factor * gh2[i__];
    }
    return 0;
} /* extrashc_ */



/* Subroutine */ int initize_()
{
    /* Local variables */
    static real erequ, erpol;

/* ---------------------------------------------------------------- */
/* Initializes the parameters in COMMON/GENER/ */

/*       UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT> */
/*       ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN */
/*                       COORDINATES (6371.2 KM) */
/*       EREQU   MAJOR HALF AXIS FOR EARTH ELLIPSOID (6378.160 KM) */
/*       ERPOL   MINOR HALF AXIS FOR EARTH ELLIPSOID (6356.775 KM) */
/*       AQUAD   SQUARE OF MAJOR HALF AXIS FOR EARTH ELLIPSOID */
/*       BQUAD   SQUARE OF MINOR HALF AXIS FOR EARTH ELLIPSOID */

/* ERA, EREQU and ERPOL as recommended by the INTERNATIONAL */
/* ASTRONOMICAL UNION . */
/* ----------------------------------------------------------------- */
    gener_1.era = (float)6371.2;
    erequ = (float)6378.16;
    erpol = (float)6356.775;
    gener_1.aquad = erequ * erequ;
    gener_1.bquad = erpol * erpol;
    gener_1.umr = atan((float)1.) * (float)4. / (float)180.;
    return 0;
} /* initize_ */

doublereal signc_(real *v1, real *v2)
{
    /* System generated locals */
    doublereal ret_val;

    ret_val = *v1;
    if (*v2 < (float)0.) {
	ret_val = -fabs(*v1);
    } else {
    	ret_val = fabs(*v1);
    }
    return ret_val;
} /* signc_ */

}

