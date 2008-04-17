namespace IGRFf2c {

  typedef float real;
  typedef long int integer;
  typedef long int logical;

  int igrf_sub__(real *xlat, real *xlong, real *year, real *height, real *xl, integer *icode, real *dip, real *dec);
  int findb0_(real *stps, real *bdel, logical *value, real *bequ, real *rr0);
  int shellg_(real *glat, real *glon, real *alt, real *dimo, real *fl, integer *icode, real *b0);
  int shellc_(real *v, real *fl, real *b0);
  int feldg_(real *glat, real *glon, real *alt, real *bnorth, real *beast, real *bdown, real *babs);
  int feldc_(real *v, real *b);
  int feldi_();
  int feldcof_(real *year, real *dimo);
  int initize_();
}
