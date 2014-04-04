#ifndef _OS_SPEC_PUSH_H_SSE
#define _OS_SPEC_PUSH_H_SSE

/* Size of particle buffer */
/* for values between 32 and 8192 this has little impact on performance */
#define p_cache_size_2D 84
#define p_cache_size_3D 1024

#include "vector-sse.h"

/***********************************************************************
t_vp2Dd

Virtual particle type for current deposition in 2D. This corresponds to 4
vectors ( 8 doubles ) in the form:

 | x0, x1 | y0, y1| q, vz | ix, iy, p0, p1 |

The pusher will store one vp per simulation particle; the splitter will
then add new vp as necessary / correct the previous one.

************************************************************************/

typedef struct virtualPart2D {
	double x0, x1;
	double y0, y1;
	double q, vz;
	int    ix, iy, p0, p1;
} t_vp2D;

/***********************************************************************
STORE2VP2D

Stores 2 virtual particles for 2D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
************************************************************************/

#define STORE2VP2D( buf, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { 			    \
  __m128i t0 = _mm_unpacklo_epi32( vix, viy );	                                \
  _mm_store_pd( &buf[0], _mm_unpacklo_pd( vx0, vx1 ) );						    \
  _mm_store_pd( &buf[2], _mm_unpacklo_pd( vy0, vy1 ) );						    \
  _mm_store_pd( &buf[4], _mm_unpacklo_pd( vq, vvz ) );						    \
  _mm_store_pd( &buf[6], _mm_castsi128_pd( t0 ) );						        \
  _mm_store_pd( &buf[8], _mm_unpackhi_pd( vx0, vx1 ) );						    \
  _mm_store_pd( &buf[10], _mm_unpackhi_pd( vy0, vy1 ) );					    \
  _mm_store_pd( &buf[12], _mm_unpackhi_pd( vq, vvz ) );						    \
  _mm_store_pd( &buf[14], _mm_castsi128_pd( _mm_unpackhi_epi64( t0, t0 ) ) );	\
}

/***********************************************************************
LOAD2VP2D

Loads 2 virtual particles for 2D current deposition. 

************************************************************************/

#define LOAD2VP2D( buf, vx0, vx1, vy0, vy1, vq, vvz, ix, iy ) { 			\
    __m128d t0, t1, t2, t4, t5, t6;                                         \
    __m128i t3, t7;                                                         \
    t0 = _mm_load_pd( &buf[0] );                                            \
    t1 = _mm_load_pd( &buf[2] );                                            \
    t2 = _mm_load_pd( &buf[4] );                                            \
    t3 = _mm_castpd_si128( _mm_load_pd( &buf[6]  ));                        \
    t4 = _mm_load_pd( &buf[8] );                                            \
    t5 = _mm_load_pd( &buf[10] );                                           \
    t6 = _mm_load_pd( &buf[12] );                                           \
    t7 = _mm_castpd_si128( _mm_load_pd( &buf[14] ));                        \
                                                                            \
    vx0 = _mm_unpacklo_pd( t0, t4 );                                        \
    vx1 = _mm_unpackhi_pd( t0, t4 );                                        \
    vy0 = _mm_unpacklo_pd( t1, t5 );                                        \
    vy1 = _mm_unpackhi_pd( t1, t5 );                                        \
    vq  = _mm_unpacklo_pd( t2, t6 );                                        \
    vvz = _mm_unpackhi_pd( t2, t6 );                                        \
                                                                            \
    t7  = _mm_unpacklo_epi32( t3, t7 );                                     \
    _mm_store_si128((__m128i*)ix, t7 );                                     \
    _mm_store_si128((__m128i*)iy, _mm_unpackhi_epi64( t7, t7 ));            \
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

#define STORE2VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { 	\
  __m128i t0 =  _mm_unpacklo_epi32( vix, viy );                                 \
  _mm_store_pd( &buf[0], _mm_unpacklo_pd( vx0, vx1 ) );						    \
  _mm_store_pd( &buf[2], _mm_unpacklo_pd( vy0, vy1 ) );						    \
  _mm_store_pd( &buf[4], _mm_unpacklo_pd( vz0, vz1 ) );						    \
  _mm_store_pd( &buf[6], _mm_unpacklo_pd( vq, vq ) );					        \
  _mm_store_pd( &buf[8], _mm_castsi128_pd( _mm_unpacklo_epi64( t0, viz ) ) );	\
  _mm_store_pd( &buf[10], _mm_unpackhi_pd( vx0, vx1 ) );						\
  _mm_store_pd( &buf[12], _mm_unpackhi_pd( vy0, vy1 ) );						\
  _mm_store_pd( &buf[14], _mm_unpackhi_pd( vz0, vz1 ) );						\
  _mm_store_pd( &buf[16], _mm_unpackhi_pd( vq, vq ) );					        \
  _mm_store_pd( &buf[18], _mm_castsi128_pd(                                     \
                     _mm_unpackhi_epi64( t0, _mm_unpacklo_epi32(viz,viz) ) ) ); \
}


/***********************************************************************
LOAD2VP3D

Loads 2 virtual particles for 3D current deposition. 

************************************************************************/
#define LOAD2VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz ) {    \
    __m128d t0, t1, t2, t3, t5, t6, t7, t8;                                 \
    __m128i t4, t9, t10;                                                    \
    t0 = _mm_load_pd( &buf[0] );                                            \
    t1 = _mm_load_pd( &buf[2] );                                            \
    t2 = _mm_load_pd( &buf[4] );                                            \
    t3 = _mm_load_pd( &buf[6] );                                            \
    t4 = _mm_castpd_si128( _mm_load_pd( &buf[8] ));                         \
    t5 = _mm_load_pd( &buf[10] );                                           \
    t6 = _mm_load_pd( &buf[12] );                                           \
    t7 = _mm_load_pd( &buf[14] );                                           \
    t8 = _mm_load_pd( &buf[16] );                                           \
    t9 = _mm_castpd_si128( _mm_load_pd( &buf[18] ));                        \
                                                                            \
    vx0 = _mm_unpacklo_pd( t0, t5 );                                        \
    vx1 = _mm_unpackhi_pd( t0, t5 );                                        \
    vy0 = _mm_unpacklo_pd( t1, t6 );                                        \
    vy1 = _mm_unpackhi_pd( t1, t6 );                                        \
    vz0 = _mm_unpacklo_pd( t2, t7 );                                        \
    vz1 = _mm_unpackhi_pd( t2, t7 );                                        \
    vq  = _mm_unpacklo_pd( t3, t8 );                                        \
                                                                            \
    t10  = _mm_unpacklo_epi32( t4, t9 );                                    \
    _mm_store_si128((__m128i*)ix, t10 );                                    \
    _mm_store_si128((__m128i*)iy, _mm_unpackhi_epi64( t10, t10 ));          \
    _mm_store_si128((__m128i*)iz, _mm_unpackhi_epi32( t4, t9 ));            \
}


// Virtual particle buffer
typedef struct { 
  // 4 splits maximum
  DECLARE_ALIGNED_16( t_vp3D buf[ 4 * p_cache_size_3D ] );
  double *p;
  unsigned int np;
} t_vpbuf3D;


#endif
