#ifndef _OS_SPEC_PUSH_H_BGP
#define _OS_SPEC_PUSH_H_BGP

/*****************************************************************************************

Precision (single/double) of PIC algorithm. Defaults to double precision

*****************************************************************************************/

#ifdef PRECISION_SINGLE

  /* Single precision */
  
  typedef float t_real;
  
  #define __LF(b)              __lfps(b)
  #define __STF(b, a)           __stfps(b,a)
  
  #define __LF2v2(v1,v2,b)     __LFPS2v2(v1,v2,b)
  #define __STF2v2(b,v1,v2)    __STFPS2v2(b,v1,v2)
  
  #define __LF2v3(v1,v2,v3,b)  __LFPS2v3(v1,v2,v3,b)
  #define __STF2v3(b,v1,v2,v3) __STFPS2v3(b,v1,v2,v3)
  
  #define LOADFLD2( flda, fldb, shift ) __cmplxf( flda[shift], fldb[shift] )

#else

  /* Double precision */
  
  typedef double t_real;
  
  #define __LF(b)              __lfpd(b)
  #define __STF(b, a)          __stfpd(b,a)
  
  #define __LF2v2(v1,v2,b)     __LFPD2v2(v1,v2,b)
  #define __STF2v2(b,v1,v2)    __STFPD2v2(b,v1,v2)
  
  #define __LF2v3(v1,v2,v3,b)  __LFPD2v3(v1,v2,v3,b)
  #define __STF2v3(b,v1,v2,v3) __STFPD2v3(b,v1,v2,v3)
  
  #define LOADFLD2( flda, fldb, shift ) __cmplx( flda[shift], fldb[shift] )
  
#endif



/* Size of particle buffer */

#define p_cache_size_2D 32
#define p_cache_size_3D 512

#include "vector-bgp.h"

/***********************************************************************
t_vp2D

Virtual particle type for current deposition in 2D. This is organized in
memory as:

 | x0, x1 | y0, y1 | q, vz | (ix, iy), p0 |

The pusher will store one vp per simulation particle; the splitter will
then add new vp as necessary / correct the previous one.

This has the size of 4 vectors / 8 doubles.

************************************************************************/

typedef struct virtualPart2D {
	double x0,  x1;
	double y0,  y1;
	double  q,  vz;
	int ix, iy; double p0;
} t_vp2D;

/***********************************************************************
STORE2VP2D

Stores 2 virtual particles for 2D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
************************************************************************/

#define STORE2VP2D( buf, vx0, vx1, vy0, vy1, vq, vvz, ix, iy ) { 			\
  register vector s = __cmplx(-1.0,1.0);                                    \
  bg_vec t;											                        \
                                                                            \
  /* 1st vp */  											  			    \
  __stfpd( &buf[ 0 ]  , __fpsel( s, vx0, __fxmr( vx1 ) ));                  \
  __stfpd( &buf[ 2 ]  , __fpsel( s, vy0, __fxmr( vy1 ) ));                  \
  __stfpd( &buf[ 4 ]  , __fpsel( s, vq,  __fxmr( vvz ) ));                  \
                                                                            \
  /* store ix.a, iy.a in buf + 6 */                                         \ 
  t.s32[0] = ix.a;                                                          \
  t.s32[1] = iy.a;                                                          \
  __stfpd( &buf[ 6 ]  , t.f );                                              \
                                                                            \
  /* 2nd vp	*/                                                              \
  __stfpd( &buf[  8 ] , __fpsel( s, __fxmr(vx0), vx1 ));                    \
  __stfpd( &buf[ 10 ], __fpsel( s, __fxmr(vy0), vy1 ));                     \
  __stfpd( &buf[ 12 ], __fpsel( s, __fxmr(vq),  vvz ));                     \
                                                                            \
  /* store ix.b, iy.b in buf + 14 */                                        \ 
  t.s32[0] = ix.b;                                                          \
  t.s32[1] = iy.b;                                                          \
  __stfpd( &buf[ 14 ]  , t.f );                                             \
                                                                            \
}

/***********************************************************************
LOAD2VP2D

Loads 2 virtual particles for 2D current deposition. 

************************************************************************/

#define LOAD2VP2D( buf, vx0, vx1, vy0, vy1, vq, vvz, ix, iy ) { 			\
	register vector s = __cmplx(-1.0,1.0);                                  \
	register vector t0, t1, t2, t4;                                         \
	bg_vec t3;								                                \
																			\
	t0    = __lfpd( &buf[ 0 ] );											\
	t1    = __lfpd( &buf[ 2 ] );											\
	t2    = __lfpd( &buf[ 4 ] );											\
	t3.f  = __lfpd( &buf[ 6 ] );											\
    	                                                                    \
	t4  = __lfpd( &buf[ 8 ] );                                              \
	vx0 = __fpsel( s, t0, __fxmr( t4 ) );                                   \
	vx1 = __fpsel( s, __fxmr( t0 ), t4 );                                   \
	                                                                        \
	t0  = __lfpd( &buf[ 10 ] );                                             \
	vy0 = __fpsel( s, t1, __fxmr( t0 ) );                                   \
	vy1 = __fpsel( s, __fxmr( t1 ), t0 );                                   \
	                                                                        \
	t1  = __lfpd( &buf[ 12 ] );                                             \
	vq  = __fpsel( s, t2, __fxmr( t1 ) );                                   \
	vvz = __fpsel( s, __fxmr( t2 ), t1 );                                   \
	                                                                        \
	ix.a = t3.s32[0];                                                       \
	iy.a = t3.s32[1];                                                       \
	                                                                        \
	t3.f  = __lfpd( &buf[ 14 ] );                                           \
                                                                            \
	ix.b = t3.s32[0];                                                       \
	iy.b = t3.s32[1];                                                       \
                                                                            \
 }


/***********************************************************************
t_vpbuf2D

Buffer to hold virtual particles for current deposition in 2D. 

************************************************************************/

typedef struct { 
  // 3 splits maximum
  DECLARE_ALIGNED_16( t_vp2D buf[ 3 * p_cache_size_2D ] );
  double *p;
  unsigned int np;
} t_vpbuf2D;

/***********************************************************************
t_vp3D

Virtual particle type for current deposition. This corresponds to 5
vectors in the form:

 | x0, x1 | y0, y1 | z0, z1 | q, dp0 | (ix, iy) , (iz, ip0) |

dp0 and ip0 are used for padding to improve memory alignment.
Structure size is 5 vectors / 10 doubles.

************************************************************************/

typedef struct virtualPart3D {
	double      x0,       x1;
	double      y0,       y1;
	double      z0,       z1;
	double       q,      dp0;
	int     ix, iy,  iz, ip0;
} t_vp3D;

/***********************************************************************
STORE2VP3D

Stores 2 virtual particles for 3D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
************************************************************************/

#define STORE2VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz ) { 	\
  register vector s = __cmplx(-1.0,1.0);                                    \
  bg_vec t;											                        \
  																			\
  __stfpd( &buf[  0 ] , __fpsel( s, vx0, __fxmr( vx1 ) ));                  \
  __stfpd( &buf[  2 ] , __fpsel( s, vy0, __fxmr( vy1 ) ));                  \
  __stfpd( &buf[  4 ] , __fpsel( s, vz0, __fxmr( vz1 ) ));                  \
  __stfpd( &buf[  6 ] , vq );                                               \
                                                                            \
  t.s32[0] = ix.a;                                                          \
  t.s32[1] = iy.a;                                                          \
  t.s32[2] = iz.a;                                                          \
  __stfpd( &buf[  8 ]  , t.f );                                             \
                                                                            \
  __stfpd( &buf[ 10 ], __fpsel( s, __fxmr( vx0 ), vx1 ));                   \
  __stfpd( &buf[ 12 ], __fpsel( s, __fxmr( vy0 ), vy1 ));                   \
  __stfpd( &buf[ 14 ], __fpsel( s, __fxmr( vz0 ), vz1 ));                   \
  __stfpd( &buf[ 16 ], __fxmr( vq ) );                                      \
                                                                            \
  t.s32[0] = ix.b;                                                          \
  t.s32[1] = iy.b;                                                          \
  t.s32[2] = iz.b;                                                          \
  __stfpd( &buf[ 18 ] , t.f );                                              \
}

/***********************************************************************
LOAD4VP3D

Loads 2 virtual particles for 3D current deposition. 

************************************************************************/

#define LOAD2VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz ) {	\
  register vector s = __cmplx(-1.0,1.0);                                    \
  register vector t0, t1, t2, t3, t5;                                       \
  bg_vec t4;								                                \
                                                                            \
  t0   = __lfpd( &buf[ 0 ]  );                                              \
  t1   = __lfpd( &buf[ 2 ]  );                                              \
  t2   = __lfpd( &buf[ 4 ]  );                                              \
  t3   = __lfpd( &buf[ 6 ]  );                                              \
                                                                            \
  t4.f = __lfpd( &buf[ 8 ]  );                                              \
  ix.a = t4.s32[0];                                                         \
  iy.a = t4.s32[1];                                                         \
  iz.a = t4.s32[2];                                                         \
                                                                            \
  t5  = __lfpd( &buf[ 10 ]  );                                              \
  vx0 = __fpsel( s, t0, __fxmr( t5 ) );                                     \
  vx1 = __fpsel( s, __fxmr( t0 ), t5 );                                     \
                                                                            \
  t0  = __lfpd( &buf[ 12 ]  );                                              \
  vy0 = __fpsel( s, t1, __fxmr( t0 ) );                                     \
  vy1 = __fpsel( s, __fxmr( t1 ), t0 );                                     \
                                                                            \
  t1  = __lfpd( &buf[ 14 ]  );                                              \
  vz0 = __fpsel( s, t2, __fxmr( t1 ) );                                     \
  vz1 = __fpsel( s, __fxmr( t2 ), t1 );                                     \
                                                                            \
  t2  = __lfpd( &buf[ 16 ]  );                                              \
  vq  = __fpsel( s, t3, __fxmr( t2 ) );                                     \
                                                                            \
  t4.f = __lfpd( &buf[ 18 ]  );                                               \
  ix.b = t4.s32[0];                                                         \
  iy.b = t4.s32[1];                                                         \
  iz.b = t4.s32[2];                                                         \
}

// Virtual particle buffer
typedef struct { 
  // 4 splits maximum
  DECLARE_ALIGNED_16( t_vp3D buf[ 4 * p_cache_size_3D ] );
  double *p;
  unsigned int np;
} t_vpbuf3D;


/*****************************************************************************************
VRGAMMA

Get 1 / Lorentz gamma. 
*****************************************************************************************/

/* This version uses FMA operations */

/*
#define VRGAMMA( vrg, vu1, vu2, vu3 ) {				\
  const vector one = __cmplx(1.0,1.0); 				\
  (vrg) = __fpmadd( one,   (vu1), (vu1) );			\
  (vrg) = __fpmadd( (vrg), (vu2), (vu2) );			\
  (vrg) = __fpmadd( (vrg), (vu3), (vu3) );			\
  (vrg) = __fprsqrt( (vrg) );						\
}
*/

#define VRGAMMA( vrg, vu1, vu2, vu3 ) {						\
   const vector one = __cmplx(1.0,1.0);						\
   (vrg) = __fprsqrt( __fpmadd( 							\
                      __fpmadd( 							\
                      __fpmadd( one,   (vu1), (vu1) ), 		\
                                       (vu2), (vu2) ), 		\
                                       (vu3), (vu3) ) ); 	\
}


/***************************************************************************************************
vget_hp_near

Gets the nearest half points of the 4 particles loaded as a vector for
interpolating scattered grids. 
This routine is for positions defined as distance to the nearest grid
points.
***************************************************************************************************/

/*
inline void vget_hp_near( vector dx, int2 ix, vector* dxh, int2* ixh, int delta )
{
  vector s;
  vector oneHalf = __cmplx( 0.5, 0.5 );
  
  // s = (dx < 0 ? -0.5 : 0.5) 
  s   = __fpsel( dx, __fpneg( oneHalf ), oneHalf );
  
  // dxh = dx - s
  *dxh = __fpsub( dx, s );

  // ixh = (dx < 0 ? -1 : 0) 
  s   = __fpsub( s, oneHalf );
  s   = __fpctiw( s );
  __stfpiw( &ixh, s );
  
  // ixh = ix - ixh & delta
  ixh.a = ix.a - (ixh.a & delta);
  ixh.b = ix.b - (ixh.b & delta);
  
}
*/

// The macro version is slightly faster

#define vget_hp_near(dx,ix,dxh,ixh,delta ) 					\
{ vector s;													\
															\
  vector oneHalf = __cmplx( 0.5, 0.5 );						\
  															\
  /* s = (dx < 0 ? -0.5 : 0.5) */							\
  s   = __fpsel( dx, __fpneg( oneHalf ), oneHalf );			\
  dxh = __fpsub( dx, s );									\
  															\
  s   = __fpsub( s, oneHalf );								\
  s   = __fpctiw( s );										\
  __stfpiw( &(ixh.a), s );									\
  															\
  ixh.a = ix.a - (ixh.a & delta);							\
  ixh.b = ix.b - (ixh.b & delta); 							\
}

inline vector vntrim( const vector vx );


#endif
