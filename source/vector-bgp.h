/****************************************************************************************/
/* SIMD utilites.                                                                       */
/* Definitions for BlueGene/P                                                           */
/****************************************************************************************/

#ifndef _VECTOR_BGP_H
#define _VECTOR_BGP_H

#include <stdio.h>
#include <complex.h>

typedef double _Complex vector;

typedef union 
{
   vector   f;
   double f64[2];
   int    s32[4];
} bg_vec;

/*
typedef struct
{
   int a;
   int b;
} int2;
*/


typedef union{
  int i[2];
  struct { int a; int b; };
} int2;


#define DECLARE_ALIGNED_16(v)       v __attribute__ ((aligned (16)))

#include <stdint.h>
#define IS_ALIGNED_16(addr) (((uintptr_t)addr & 0xF) == 0)


/* 
Number of Newton-Raphson refinements for the full precision operations. The maximum
is 2 (default) that for a = M_PI gives exactly the same results as the normal scalar op. */

#ifndef __MATH_PRECISION__
#define __MATH_PRECISION__ 2
#endif

/* Full precision reciprocal, reciprocal square root                                    */

/*   __fpr - Parallel Reciprocal  */

inline vector __fpr( const vector a )
{
   
   // Get the reciprocal estimate
   vector estimate = __fpre( a );

#if __MATH_PRECISION__ > 0   
   vector one     = __cmplx(1.0,1.0);
   estimate = __fpmadd( estimate, estimate, __fpnmsub( one, a, estimate ) );
#endif

#if __MATH_PRECISION__ > 1   
   estimate = __fpmadd( estimate, estimate, __fpnmsub( one, a, estimate ) );
#endif

   return estimate;
}

/* __fprsqrt - Parallel Reciprocal square root */

inline vector __fprsqrt( const vector a )
{
   // Get the reciprocal square root estimate
   vector estimate = __fprsqrte( a );

#if __MATH_PRECISION__ > 0   
   vector one     = __cmplx(1.0,1.0);
  
   // One round of Newton-Raphson refinement
   vector estimateSquared = __fpmul( estimate, estimate );
   vector halfEstimate    = __fxpmul( estimate, 0.5 );
   estimate = __fpmadd( estimate, halfEstimate,  __fpnmsub( one, estimateSquared, a ));
#endif

#if __MATH_PRECISION__ > 1   
   estimateSquared = __fpmul( estimate, estimate );
   halfEstimate    = __fxpmul( estimate, 0.5 );
   estimate = __fpmadd( estimate, halfEstimate,  __fpnmsub( one, estimateSquared, a ));
#endif
   
   return estimate;
}


/* __fpsqrt - Parallel square root */
inline vector __fpsqrt( const vector a )
{  
  return __fpmul( a, __fprsqrt(a) );
  
}


/* Division is implemented as a / b = a * (1/b)                                        */


/* __fpdive - Parallel division estimate  */

inline vector __fpdive( vector a, vector b )
{
   return __fpmul( a, __fpre(b) );
}

inline vector __fxpdive( double a, vector b )
{
   return __fxpmul( __fpre(b), a );
}

/* __fpdiv - Parallel division */

inline vector __fpdiv( vector a, vector b )
{
   return __fpmul( a, __fpr(b) );
}

inline vector __fxpdiv( double a, vector b )
{
   return __fxpmul( __fpr(b), a );
}


/***********************************************************************
__LFPD2v3

Loads 2 x 3 elements vectors stored sequentially in memory into 3 vector
variables:

b = [x1,y1,z1,x2,y2,z2] -> __LFPD2v3(v1,v2,v3,b)

-> v1 = [x1,x2]
-> v2 = [y1,y2]
-> v3 = [z1,z2]

************************************************************************/

#define __LFPD2v3(v1, v2, v3, b) {               \
                                                 \
  register vector a0,a1,a2,s;                    \
                                                 \
  s = __cmplx( -1.0, 1.0 );                      \
                                                 \
  a0 = __lfpd( &b[0] );                          \
  a1 = __lfpd( &b[2] );                          \
  a2 = __lfpd( &b[4] );                          \
                                                 \
  (v1) = __fpsel( s, a0, a1 );                   \
  (v2) = __fxmr( __fpsel( __fpneg(s), a0, a2 )); \
  (v3) = __fpsel( s, a1, a2 );                   \
}

#define __LFPS2v3(v1, v2, v3, b) {               \
                                                 \
  register vector a0,a1,a2,s;                    \
                                                 \
  s = __cmplx( -1.0, 1.0 );                      \
                                                 \
  a0 = __lfps( &b[0] );                          \
  a1 = __lfps( &b[2] );                          \
  a2 = __lfps( &b[4] );                          \
                                                 \
  (v1) = __fpsel( s, a0, a1 );                   \
  (v2) = __fxmr( __fpsel( __fpneg(s), a0, a2 )); \
  (v3) = __fpsel( s, a1, a2 );                   \
}

/***********************************************************************
__STFPD2v3

Stores 3 vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2], v2 = [y1,y2], v3 = [z1,z2]

-> __STFPD2v3(b, v1, v2, v3)

-> b = [x1,y1,z1,x2,y2,z2]

************************************************************************/

#define __STFPD2v3(b, v1, v2, v3) {              \
                                                 \
  register vector t,s;                           \
                                                 \
  s = __cmplx( -1.0, 1.0 );                      \
                                                 \
  t = __fxmr( v2 );                              \
  __stfpd(&b[0], __fpsel( s, v1, t  ));          \
  __stfpd(&b[2], __fpsel( s, v3, v1 ));          \
  __stfpd(&b[4], __fpsel( s, t,  v3 ));          \
}

#define __STFPS2v3(b, v1, v2, v3) {              \
                                                 \
  register vector t,s;                           \
                                                 \
  s = __cmplx( -1.0, 1.0 );                      \
                                                 \
  t = __fxmr( v2 );                              \
  __stfps(&b[0], __fpsel( s, v1, t  ));          \
  __stfps(&b[2], __fpsel( s, v3, v1 ));          \
  __stfps(&b[4], __fpsel( s, t,  v3 ));          \
}

/***********************************************************************
__LFPD2v2

Loads 2 x 2 elements vectors stored sequentially in memory into 2 vector
variables:

b = [x1,y1,x2,y2] -> __LFPD2v2(v1,v2,b)

-> v1 = [x1,x2]
-> v2 = [y1,y2]

************************************************************************/

#define __LFPD2v2(v1, v2, b) {                   \
                                                 \
  register vector a0,a1,s;                       \
                                                 \
  s = __cmplx( -1.0, 1.0 );                      \
                                                 \
  a0 = __lfpd( &b[0] );                          \
  a1 = __lfpd( &b[2] );                          \
                                                 \
  (v1) = __fpsel( s, a0, __fxmr(a1) );           \
  (v2) = __fpsel( s, __fxmr(a0), a1 );           \
}

#define __LFPS2v2(v1, v2, b) {                   \
                                                 \
  register vector a0,a1,s;                       \
                                                 \
  s = __cmplx( -1.0, 1.0 );                      \
                                                 \
  a0 = __lfps( &b[0] );                          \
  a1 = __lfps( &b[2] );                          \
                                                 \
  (v1) = __fpsel( s, a0, __fxmr(a1) );           \
  (v2) = __fpsel( s, __fxmr(a0), a1 );           \
}


/***********************************************************************
__STFPD2v2

Stores 2 vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2], v2 = [y1,y2]

-> __STFPD2v2(b, v1, v2)

-> b = [x1,y1,x2,y2]

*********************************************************************/

#define __STFPD2v2(b, v1, v2) {                  \
                                                 \
  register vector s;                             \
                                                 \
  s = __cmplx( -1.0, 1.0 );                      \
                                                 \
  __stfpd(&b[0], __fpsel( s, v1, __fxmr(v2) ));  \
  __stfpd(&b[2], __fpsel( s, __fxmr(v1), v2 ));  \
}

#define __STFPS2v2(b, v1, v2) {                  \
                                                 \
  register vector s;                             \
                                                 \
  s = __cmplx( -1.0, 1.0 );                      \
                                                 \
  __stfps(&b[0], __fpsel( s, v1, __fxmr(v2) ));  \
  __stfps(&b[2], __fpsel( s, __fxmr(v1), v2 ));  \
}

/***********************************************************************
printfv

Print vector components
***********************************************************************/

static void printfv(char *tag, vector v)
{
  printf("%s = %f, %f \n", tag, __creal(v), __cimag(v) );
}

/***********************************************************************
test vector components

***********************************************************************/
#define TEST_VECTOR( a ) {							\
if ( isinf(__creal(a)) || isinf(__cimag(a)) ||		\
     isnan(__creal(a)) || isnan(__cimag(a)) ) {		\
     printf( "%s:%d %s = %g, %g\n", __FILE__, __LINE__, #a, __creal(a), __cimag(a) );	\
     exit(-1);										\
   }												\
}

#endif
