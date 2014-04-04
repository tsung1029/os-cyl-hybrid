/*

 Spline function definitions for x86 AVX extensions

*/

#ifndef _SPLINES_H
#define _SPLINES_H


#include "vector-avx.h"

/********************************* Linear Interpolation *********************************/ 

static inline void vspline_s1( __m256 const dx, __m256 w[] ) {
  
  w[0] = _mm256_sub_ps( _mm256_set1_ps( 0.5f ), dx );
  w[1] = _mm256_add_ps( _mm256_set1_ps( 0.5f ), dx );
}

static inline void vsplined_s1( __m256d const dx, __m256d w[] ) {
  
  w[0] = _mm256_sub_pd( _mm256_set1_pd( 0.5 ), dx );
  w[1] = _mm256_add_pd( _mm256_set1_pd( 0.5 ), dx );
}

/******************************** Quadratic Interpolation *******************************/ 


static inline void vspline_s2( __m256 const dx, __m256 w[] ) {
  
  __m256 tmp;
  __m256 const oneHalf = _mm256_set1_ps( 0.5f );
  
  tmp  = _mm256_sub_ps( oneHalf, dx );
  w[0] = _mm256_mul_ps( oneHalf, _mm256_mul_ps( tmp, tmp ) );

  w[1]  = _mm256_sub_ps( _mm256_set1_ps( 0.75f ), _mm256_mul_ps( dx, dx ) );
  
  tmp  = _mm256_add_ps( oneHalf, dx );
  w[2] = _mm256_mul_ps( oneHalf, _mm256_mul_ps( tmp, tmp ) );  
}

static inline void vsplined_s2( __m256d const dx, __m256d w[] ) {
  
  __m256d tmp;
  __m256d const oneHalf = _mm256_set1_pd( 0.5 );
  
  tmp  = _mm256_sub_pd( oneHalf, dx );
  w[0] = _mm256_mul_pd( oneHalf, _mm256_mul_pd( tmp, tmp ) );

  w[1]  = _mm256_sub_pd( _mm256_set1_pd( 0.75 ), _mm256_mul_pd( dx, dx ) );
  
  tmp  = _mm256_add_pd( oneHalf, dx );
  w[2] = _mm256_mul_pd( oneHalf, _mm256_mul_pd( tmp, tmp ) );  
}

/********************************** Cubic Interpolation *********************************/ 

static inline void vspline_s3( __m256 const dx, __m256 w[] ) {
  
  register __m256 t0, t1, t2, t3;
  register __m256 const c1_2 = _mm256_set1_ps( 1.0f/2.0f );
  register __m256 const c1_6 = _mm256_set1_ps( 1.0f/6.0f );
  register __m256 const c2_3 = _mm256_set1_ps( 2.0f/3.0f );

  t0 = _mm256_sub_ps( c1_2, dx );
  t1 = _mm256_add_ps( c1_2, dx );
  
  t2 = _mm256_mul_ps( t0, t0 );  // t2 = (1/2 - dx)^2
  t3 = _mm256_mul_ps( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm256_mul_ps( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm256_mul_ps( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm256_mul_ps( c1_6, t0 );
  w[1] = _mm256_add_ps( _mm256_sub_ps( c2_3, t3 ), _mm256_mul_ps( c1_2, t1 ) );
  w[2] = _mm256_add_ps( _mm256_sub_ps( c2_3, t2 ), _mm256_mul_ps( c1_2, t0 ) );
  w[3] = _mm256_mul_ps( c1_6, t1 ); 

}

static inline void vsplined_s3( __m256d const dx, __m256d w[] ) {
  
  register __m256d t0, t1, t2, t3;
  register __m256d const c1_2 = _mm256_set1_pd( 1.0/2.0 );
  register __m256d const c1_6 = _mm256_set1_pd( 1.0/6.0 );
  register __m256d const c2_3 = _mm256_set1_pd( 2.0/3.0 );

  t0 = _mm256_sub_pd( c1_2, dx );
  t1 = _mm256_add_pd( c1_2, dx );
  
  t2 = _mm256_mul_pd( t0, t0 );  // t2 = (1/2 - dx)^2
  t3 = _mm256_mul_pd( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm256_mul_pd( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm256_mul_pd( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm256_mul_pd( c1_6, t0 );
  w[1] = _mm256_add_pd( _mm256_sub_pd( c2_3, t3 ), _mm256_mul_pd( c1_2, t1 ) );
  w[2] = _mm256_add_pd( _mm256_sub_pd( c2_3, t2 ), _mm256_mul_pd( c1_2, t0 ) );
  w[3] = _mm256_mul_pd( c1_6, t1 ); 

}



/********************************* Quartic Interpolation ********************************/ 

static inline void vspline_s4( __m256 const dx, __m256 w[] ) {
  
  __m256 const c1_48 = _mm256_set1_ps( 1.0f/48.0f );
  __m256 const c1_6 = _mm256_set1_ps( 1.0f/6.0f );
  
  
  __m256 dx2 = _mm256_mul_ps( dx, dx );
  __m256 dx3 = _mm256_mul_ps( dx, dx2 );
  __m256 dx4 = _mm256_mul_ps( dx2, dx2 );
  
  // t0 = 1/8 + 3 dx^2 + 2 dx^4
  __m256 t0 = _mm256_add_ps( _mm256_add_ps( _mm256_set1_ps( 1.0f/8.0f ), 
                                      _mm256_mul_ps( _mm256_set1_ps(3.0f), dx2 ) ), 
                          _mm256_add_ps( dx4, dx4 ) );
  // t1 = dx + 4 dx^3
  __m256 t1 = _mm256_add_ps( dx, _mm256_mul_ps( _mm256_set1_ps(4.0f), dx3 ) );
  
  // t2 = 19/16 + (3/2) dx^2 - dx^4
  __m256 t2 = _mm256_add_ps( _mm256_set1_ps( 19.0f/16.0f ), 
                          _mm256_sub_ps( _mm256_mul_ps( _mm256_set1_ps(3.0f/2.0f), dx2 ), dx4 ) );
  
  // t3 = dx^3 - (11/4) dx
  __m256 t3 = _mm256_sub_ps( dx3, _mm256_mul_ps( _mm256_set1_ps(11.0f/4.0f), dx ) );
  
  
  w[0] = _mm256_mul_ps( c1_48, _mm256_sub_ps( t0, t1 ) );
  w[1] = _mm256_mul_ps(  c1_6, _mm256_add_ps( t2, t3 ) );
  w[2] = _mm256_add_ps( _mm256_set1_ps( 115.0f / 192.0f ), 
                     _mm256_sub_ps( _mm256_mul_ps( _mm256_set1_ps( 1.0f/4.0f ), dx4 ),
                                 _mm256_mul_ps( _mm256_set1_ps( 5.0f/8.0f ), dx2 ) ) );
  w[3] = _mm256_mul_ps(  c1_6, _mm256_sub_ps( t2, t3 ) );
  w[4] = _mm256_mul_ps( c1_48, _mm256_add_ps( t0, t1 ) );

}

static inline void vsplined_s4( __m256d const dx, __m256d w[] ) {
  
  __m256d const c1_48 = _mm256_set1_pd( 1.0/48.0 );
  __m256d const c1_6 = _mm256_set1_pd( 1.0/6.0 );
  
  
  __m256d dx2 = _mm256_mul_pd( dx, dx );
  __m256d dx3 = _mm256_mul_pd( dx, dx2 );
  __m256d dx4 = _mm256_mul_pd( dx2, dx2 );
  
  // t0 = 1/8 + 3 dx^2 + 2 dx^4
  __m256d t0 = _mm256_add_pd( _mm256_add_pd( _mm256_set1_pd( 1.0/8.0 ), 
                                      _mm256_mul_pd( _mm256_set1_pd(3.0), dx2 ) ), 
                          _mm256_add_pd( dx4, dx4 ) );
  // t1 = dx + 4 dx^3
  __m256d t1 = _mm256_add_pd( dx, _mm256_mul_pd( _mm256_set1_pd(4.0), dx3 ) );
  
  // t2 = 19/16 + (3/2) dx^2 - dx^4
  __m256d t2 = _mm256_add_pd( _mm256_set1_pd( 19.0/16.0 ), 
                          _mm256_sub_pd( _mm256_mul_pd( _mm256_set1_pd(3.0/2.0), dx2 ), dx4 ) );
  
  // t3 = dx^3 - (11/4) dx
  __m256d t3 = _mm256_sub_pd( dx3, _mm256_mul_pd( _mm256_set1_pd(11.0/4.0), dx ) );
  
  
  w[0] = _mm256_mul_pd( c1_48, _mm256_sub_pd( t0, t1 ) );
  w[1] = _mm256_mul_pd(  c1_6, _mm256_add_pd( t2, t3 ) );
  w[2] = _mm256_add_pd( _mm256_set1_pd( 115.0 / 192.0 ), 
                     _mm256_sub_pd( _mm256_mul_pd( _mm256_set1_pd( 1.0/4.0 ), dx4 ),
                                 _mm256_mul_pd( _mm256_set1_pd( 5.0/8.0 ), dx2 ) ) );
  w[3] = _mm256_mul_pd(  c1_6, _mm256_sub_pd( t2, t3 ) );
  w[4] = _mm256_mul_pd( c1_48, _mm256_add_pd( t0, t1 ) );

}

#endif /* _SPLINES_H */
