/*****************************************************************************************

 Spline function definitions, BG/P optimized version

*****************************************************************************************/

#ifndef _SPLINES_H
#define _SPLINES_H

// #include "vector-bgp.h"


/*****************************************************************************************
vspline_n1

1st order (linear) splines.

w[0] = 1/2 - x
w[1] = 1/2 + x
*****************************************************************************************/


 inline void vspline_s1( vector const dx, vector w[] ) {
  
  vector const c1_2 = __cmplx(0.5,   0.5);
  
  w[0] = __fpsub( c1_2, dx );
  w[1] = __fpadd( c1_2, dx );
}


/*****************************************************************************************
vspline_n2

2nd order (quadratic) splines.


w[0] = (1 - 2*dx)^2/8.
w[1] = 0.75 - dx^2
w[2] = (1 + 2*dx)^2/8.
1
*****************************************************************************************/


 inline void vspline_s2( vector const dx, vector w[] ) {
  
  vector const oneHalf      = __cmplx(0.5,   0.5);
  vector const oneEight     = __cmplx(0.125, 0.125);
  vector const threeQuarter = __cmplx(0.75,  0.75);
  
  w[0] = __fpmadd( oneEight, dx, __fpmsub( oneHalf, dx, oneHalf ));
  w[1] = __fpnmsub( threeQuarter, dx, dx );
  w[2] = __fpmadd( oneEight, dx, __fpmadd( oneHalf, dx, oneHalf ));
}

/*****************************************************************************************
vspline_n3

3nd order (cubic) splines.

w[0] =  1/6 (1/2 - dx)^3
w[1] = 1/6 (4 - 6 (1/2 + dx)^2 + 3 (1/2 + dx)^3)
w[2] = 1/6 (4 - 6 (1/2 - dx)^2 + 3 (1/2 - dx)^3)
w[3] = 1/6 (1/2 + dx)^3
*****************************************************************************************/


 inline void vspline_s3( vector const dx, vector w[] ) {
  
  vector const c1_2 = __cmplx( 0.5    , 0.5     );
  vector const c2_3 = __cmplx( 2.0/3.0, 2.0/3.0 );
  double const c1_6 = 1.0/6.0;
  
  vector t0, t1, t2, t3;
  
  t0 = __fpsub( c1_2, dx );
  t1 = __fpadd( c1_2, dx );
  
  t2 = __fpmul( t0, t0 );
  t3 = __fpmul( t1, t1 );
  
  t0 = __fpmul( t0, t2 );
  t1 = __fpmul( t1, t3 );
  
  w[0] = __fxpmul( t0, c1_6 );
  w[1] = __fpmadd( __fpsub( c2_3, t3 ), c1_2, t1 );
  w[2] = __fpmadd( __fpsub( c2_3, t2 ), c1_2, t0 );
  w[3] = __fxpmul( t1, c1_6 );
}

/*****************************************************************************************
vspline_n4

4th order (quartic) splines.

w[0] = (1/384) (1 - 2*x)^4
w[1] = (1/96) (19 - 44*x + 24*x^2 + 16*x^3 - 16*x^4)
w[2] = (115/192) - (5*x^2)/8 + x^4/4
w[3] = (1/96) (19 + 44*x + 24*x^2 - 16*x^3 - 16*x^4)
w[4] = (1/384) (1 + 2*x)^4
*****************************************************************************************/


 inline void vspline_s4( vector const dx, vector w[] ) {
  
   vector const c1_8      = __cmplx( 0.125       , 0.125        );
   vector const c3        = __cmplx( 3.          , 3.           );
   vector const c2        = __cmplx( 2.          , 2.           );
   vector const c4        = __cmplx( 4.          , 4.           );
   
   vector const c19_96    = __cmplx( 19.0/96.0   , 19.0/96.0    );
   vector const c1_4      = __cmplx( 0.25        , 0.25         );
   vector const c1_6      = __cmplx( 1.0/6.0     , 1.0/6.0      );   
    
   vector const c11_24    = __cmplx( 11.0/24.0   , 11.0/24.0    );
   double const c1_48     = 1.0/48.0;
   vector const c115_192  = __cmplx( 115.0/192.0 , 115.0/192.0  );
   vector const c5_8      = __cmplx( 0.625       , 0.625        );
   
   vector dx2 = __fpmul(dx, dx);
   vector dx3 = __fpmul(dx, dx2);
   vector dx4 = __fpmul(dx2, dx2);
   
   vector t0 = __fpmadd(__fpmadd(c1_8, c3, dx2), c2, dx4);
   vector t1 = __fpmadd(dx, c4, dx3);
   
   vector t2 = __fpnmsub(__fpmadd(c19_96, c1_4, dx2), c1_6, dx4);
   vector t3 = __fpnmsub(__fpmul(c11_24, dx), c1_6, dx3);
   
   w[0] = __fxpmul(__fpsub(t0, t1), c1_48);
   w[1] = __fpsub(t2, t3);
   w[2] = __fpnmsub(__fpmadd( c115_192, c1_4, dx4), c5_8, dx2);
   w[3] = __fpadd(t2, t3);
   w[4] = __fxpmul(__fpadd(t0, t1), c1_48);  
}


#endif /* _SPLINES_H */
