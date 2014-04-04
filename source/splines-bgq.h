/*****************************************************************************************

 Spline function definitions, BG/Q optimized version

*****************************************************************************************/

#ifndef _SPLINES_H
#define _SPLINES_H

// #include "vector-bgq.h"


/*****************************************************************************************
vspline_n1

1st order (linear) splines.

w[0] = 1/2 - x
w[1] = 1/2 + x
*****************************************************************************************/

inline void vspline_s1( vector const dx, vector w[] ) {
  
  vector const c1_2 = vec_splats(0.5);
  
  w[0] = vec_sub( c1_2, dx );
  w[1] = vec_add( c1_2, dx );
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
  
  vector const oneHalf      = vec_splats(0.5);
  vector const oneEight     = vec_splats(0.125);
  vector const threeQuarter = vec_splats(0.75);
  
  w[0] = vec_madd( vec_msub( oneHalf, dx, oneHalf ), dx, oneEight);
  w[1] = vec_nmsub( dx, dx, threeQuarter );
  w[2] = vec_madd( vec_madd( oneHalf, dx, oneHalf ), dx, oneEight);
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
  
  vector const c1_2 = vec_splats( 0.5     );
  vector const c2_3 = vec_splats( 2.0/3.0 );
  vector const c1_6 = vec_splats( 1.0/6.0 );
  
  vector t0, t1, t2, t3;
  
  t0 = vec_sub( c1_2, dx );
  t1 = vec_add( c1_2, dx );
  
  t2 = vec_mul( t0, t0 );
  t3 = vec_mul( t1, t1 );
  
  t0 = vec_mul( t0, t2 );
  t1 = vec_mul( t1, t3 );
  
  w[0] = vec_mul( t0, c1_6 );
  w[1] = vec_madd( t1, c1_2, vec_sub( c2_3, t3 ) );
  w[2] = vec_madd( t0, c1_2, vec_sub( c2_3, t2 ) );
  w[3] = vec_mul( t1, c1_6 );
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
  
   vector const c1_8      = vec_splats( 0.125       );
   vector const c3        = vec_splats( 3.          );
   vector const c2        = vec_splats( 2.          );
   vector const c4        = vec_splats( 4.          );
   
   vector const c19_96    = vec_splats( 19.0/96.0   );
   vector const c1_4      = vec_splats( 0.25        );
   vector const c1_6      = vec_splats( 1.0/6.0     );   
    
   vector const c11_24    = vec_splats( 11.0/24.0   );
   vector const c1_48     = vec_splats( 1.0/48.0 );
   vector const c115_192  = vec_splats( 115.0/192.0 );
   vector const c5_8      = vec_splats( 0.625       );
   
   vector dx2 = vec_mul(dx, dx);
   vector dx3 = vec_mul(dx, dx2);
   vector dx4 = vec_mul(dx2, dx2);
   
   vector t0 = vec_madd(dx4, c2, vec_madd(dx2, c3, c1_8));
   vector t1 = vec_madd(dx3, c4, dx);
   
   vector t2 = vec_nmsub(dx4, c1_6, vec_madd(dx2, c1_4, c19_96));
   vector t3 = vec_nmsub(dx3, c1_6, vec_mul(c11_24, dx));
   
   w[0] = vec_mul(vec_sub(t0, t1), c1_48);
   w[1] = vec_sub(t2, t3);
   w[2] = vec_nmsub(dx2, c5_8, vec_madd( dx4, c1_4, c115_192));
   w[3] = vec_add(t2, t3);
   w[4] = vec_mul(vec_add(t0, t1), c1_48);  
}


#endif /* _SPLINES_H */
