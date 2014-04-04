#ifndef _OS_SPEC_PUSH_H
#define _OS_SPEC_PUSH_H

/* Size of particle buffer */
/* for values between 32 and 8192 this has little impact on performance */
#define p_cache_size_2D 84
#define p_cache_size_3D 1024

#include "vector-avx.h"

/***********************************************************************
t_vp2Dd

Virtual particle type for current deposition in 2D. This corresponds to 2
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
STORE4VP2D

Stores 4 virtual particles for 2D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
************************************************************************/

#define STORE4VP2D( b, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
  __m256d t0, t1, t2, t4, t5, t6; \
  __m128i t3, t7; \
  \
  t0 = _mm256_unpacklo_pd( vx0, vx1 ); \
  t1 = _mm256_unpacklo_pd( vy0, vy1 ); \
  t2 = _mm256_unpacklo_pd( vq, vvz ); \
  t3 = _mm_unpacklo_epi32( vix, viy ); \
  t4 = _mm256_unpackhi_pd( vx0, vx1 ); \
  t5 = _mm256_unpackhi_pd( vy0, vy1 ); \
  t6 = _mm256_unpackhi_pd( vq, vvz ); \
  t7 = _mm_unpackhi_epi32( vix, viy ); \
  \
  _mm_store_pd( &b[ 0], _mm256_castpd256_pd128( t0 ) ); \
  _mm_store_pd( &b[ 2], _mm256_castpd256_pd128( t1 ) ); \
  _mm_store_pd( &b[ 4], _mm256_castpd256_pd128( t2 ) ); \
  _mm_store_si128( (__m128i *) &b[ 6], t3 ); \
  \
  _mm_store_pd( &b[ 8], _mm256_castpd256_pd128( t4 ) ); \
  _mm_store_pd( &b[10], _mm256_castpd256_pd128( t5 ) ); \
  _mm_store_pd( &b[12], _mm256_castpd256_pd128( t6 ) ); \
  _mm_store_si128( (__m128i *) &b[14], _mm_unpackhi_epi64( t3, t3 ) ); \
  \
  _mm_store_pd( &b[16], _mm256_extractf128_pd( t0, 1 ) ); \
  _mm_store_pd( &b[18], _mm256_extractf128_pd( t1, 1 ) ); \
  _mm_store_pd( &b[20], _mm256_extractf128_pd( t2, 1 ) ); \
  _mm_store_si128( (__m128i *) &b[22], t7 ); \
  \
  _mm_store_pd( &b[24], _mm256_extractf128_pd( t4, 1 ) ); \
  _mm_store_pd( &b[26], _mm256_extractf128_pd( t5, 1 ) ); \
  _mm_store_pd( &b[28], _mm256_extractf128_pd( t6, 1 ) ); \
  _mm_store_si128( (__m128i *) &b[30], _mm_unpackhi_epi64( t7, t7 ) ); \
}


/***********************************************************************
LOAD4VP2D

Loads 4 virtual particles for 2D current deposition. 

************************************************************************/

#define LOAD4VP2D( buf, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
  __m256d t0, t1, t2, t3, t4, t5, t6, t7; \
  t0  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 0]) ); \
  t1  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 2]) ); \
  t2  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 4]) ); \
  t3  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 6]) ); \
  t4  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 8]) ); \
  t5  = _mm256_castpd128_pd256( _mm_load_pd(&buf[10]) ); \
  t6  = _mm256_castpd128_pd256( _mm_load_pd(&buf[12]) ); \
  t7  = _mm256_castpd128_pd256( _mm_load_pd(&buf[14]) ); \
  \
  t0  = _mm256_insertf128_pd( t0, _mm_load_pd(&buf[16]), 1); \
  t1  = _mm256_insertf128_pd( t1, _mm_load_pd(&buf[18]), 1); \
  t2  = _mm256_insertf128_pd( t2, _mm_load_pd(&buf[20]), 1); \
  t3  = _mm256_insertf128_pd( t3, _mm_load_pd(&buf[22]), 1); \
  t4  = _mm256_insertf128_pd( t4, _mm_load_pd(&buf[24]), 1); \
  t5  = _mm256_insertf128_pd( t5, _mm_load_pd(&buf[26]), 1); \
  t6  = _mm256_insertf128_pd( t6, _mm_load_pd(&buf[28]), 1); \
  t7  = _mm256_insertf128_pd( t7, _mm_load_pd(&buf[30]), 1); \
  \
  vx0 = _mm256_unpacklo_pd( t0, t4 ); \
  vx1 = _mm256_unpackhi_pd( t0, t4 ); \
  \
  vy0 = _mm256_unpacklo_pd( t1, t5 ); \
  vy1 = _mm256_unpackhi_pd( t1, t5 ); \
  \
  vq  = _mm256_unpacklo_pd( t2, t6 ); \
  vvz = _mm256_unpackhi_pd( t2, t6 ); \
  \
  __m256 a = _mm256_castpd_ps( _mm256_unpacklo_pd( t3, t7 ) ); \
  __m128 b = _mm256_castps256_ps128( a ); \
  __m128 c = _mm256_extractf128_ps( a, 1 ); \
  \
  vix = _mm_castps_si128( _mm_shuffle_ps( b, c, _MM_SHUFFLE( 2, 0, 2, 0 ) ) ); \
  viy = _mm_castps_si128( _mm_shuffle_ps( b, c, _MM_SHUFFLE( 3, 1, 3, 1 ) ) ); \
}


/***********************************************************************
t_vpbuf2D

Buffer to hold virtual particles for current deposition in 2D. 

************************************************************************/

typedef struct { 
  // 3 splits maximum
  DECLARE_ALIGNED_32( t_vp2D buf[ 3 * p_cache_size_2D ] );
  double *p;
  unsigned int np;
} t_vpbuf2D;


/***********************************************************************
t_vp3D

Virtual particle type for current deposition. 

 | x0, x1 | y0, y1 | z0, z1 | q, dp0 | (ix, iy) , (iz, ip0) |

dp0 and ip0 are used for padding to improve memory alignment.
Structure size is 10 doubles.

************************************************************************/

typedef struct virtualPart3D {
	double      x0,       x1;	double      y0,       y1;
	double      z0,       z1;	double       q,      dp0;
	int     ix, iy, iz, ip0;
} t_vp3D;

/***********************************************************************
STORE4VP3D

Stores 4 virtual particles for 3D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
************************************************************************/

#define STORE4VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
  __m256d t0, t1, t2, t3, t5, t6, t7, t8; \
  __m128i t4a, t4b, t9a, t9b; \
  \
  t0 = _mm256_unpacklo_pd( vx0, vx1 ); \
  t1 = _mm256_unpacklo_pd( vy0, vy1 ); \
  t2 = _mm256_unpacklo_pd( vz0, vz1 ); \
  t3 = _mm256_unpacklo_pd( vq, vq ); \
  t4a = _mm_unpacklo_epi32( vix, viy ); \
  t4b = _mm_unpacklo_epi32( viz, viz ); \
  t5 = _mm256_unpackhi_pd( vx0, vx1 ); \
  t6 = _mm256_unpackhi_pd( vy0, vy1 ); \
  t7 = _mm256_unpackhi_pd( vz0, vz1 ); \
  t8 = _mm256_unpackhi_pd( vq, vq ); \
  t9a = _mm_unpackhi_epi32( vix, viy ); \
  t9b = _mm_unpackhi_epi32( viz, viz ); \
  \
  _mm_store_pd( &buf[ 0], _mm256_castpd256_pd128( t0 ) ); \
  _mm_store_pd( &buf[ 2], _mm256_castpd256_pd128( t1 ) ); \
  _mm_store_pd( &buf[ 4], _mm256_castpd256_pd128( t2 ) ); \
  _mm_store_pd( &buf[ 6], _mm256_castpd256_pd128( t3 ) ); \
  _mm_store_pd( &buf[ 8], _mm_castsi128_pd( _mm_unpacklo_epi64( t4a, t4b ) ) ); \
  \
  _mm_store_pd( &buf[10], _mm256_castpd256_pd128( t5 ) ); \
  _mm_store_pd( &buf[12], _mm256_castpd256_pd128( t6 ) ); \
  _mm_store_pd( &buf[14], _mm256_castpd256_pd128( t7 ) ); \
  _mm_store_pd( &buf[16], _mm256_castpd256_pd128( t8 ) ); \
  _mm_store_pd( &buf[18], _mm_castsi128_pd( _mm_unpackhi_epi64( t4a, t4b ) ) ); \
  \
  _mm_store_pd( &buf[20], _mm256_extractf128_pd( t0, 1 ) ); \
  _mm_store_pd( &buf[22], _mm256_extractf128_pd( t1, 1 ) ); \
  _mm_store_pd( &buf[24], _mm256_extractf128_pd( t2, 1 ) ); \
  _mm_store_pd( &buf[26], _mm256_extractf128_pd( t3, 1 ) ); \
  _mm_store_pd( &buf[28], _mm_castsi128_pd( _mm_unpacklo_epi64( t9a, t9b ) ) ); \
  \
  _mm_store_pd( &buf[30], _mm256_extractf128_pd( t5, 1 ) ); \
  _mm_store_pd( &buf[32], _mm256_extractf128_pd( t6, 1 ) ); \
  _mm_store_pd( &buf[34], _mm256_extractf128_pd( t7, 1 ) ); \
  _mm_store_pd( &buf[36], _mm256_extractf128_pd( t8, 1 ) ); \
  _mm_store_pd( &buf[38], _mm_castsi128_pd( _mm_unpackhi_epi64( t9a, t9b ) ) ); \
}


/***********************************************************************
LOAD4VP3D

Loads 4 virtual particles for 3D current deposition. 

************************************************************************/
#define LOAD4VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz ) { \
  __m256d t0, t1, t2, t3, t4, t5, t6, t7, t8, t9; \
  t0  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 0]) ); \
  t1  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 2]) ); \
  t2  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 4]) ); \
  t3  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 6]) ); \
  t4  = _mm256_castpd128_pd256( _mm_load_pd(&buf[ 8]) ); \
  \
  t5  = _mm256_castpd128_pd256( _mm_load_pd(&buf[10]) ); \
  t6  = _mm256_castpd128_pd256( _mm_load_pd(&buf[12]) ); \
  t7  = _mm256_castpd128_pd256( _mm_load_pd(&buf[14]) ); \
  t8  = _mm256_castpd128_pd256( _mm_load_pd(&buf[16]) ); \
  t9  = _mm256_castpd128_pd256( _mm_load_pd(&buf[18]) ); \
  \
  t0  = _mm256_insertf128_pd( t0, _mm_load_pd(&buf[20]), 1); \
  t1  = _mm256_insertf128_pd( t1, _mm_load_pd(&buf[22]), 1); \
  t2  = _mm256_insertf128_pd( t2, _mm_load_pd(&buf[24]), 1); \
  t3  = _mm256_insertf128_pd( t3, _mm_load_pd(&buf[26]), 1); \
  t4  = _mm256_insertf128_pd( t4, _mm_load_pd(&buf[28]), 1); \
  \
  t5  = _mm256_insertf128_pd( t5, _mm_load_pd(&buf[30]), 1); \
  t6  = _mm256_insertf128_pd( t6, _mm_load_pd(&buf[32]), 1); \
  t7  = _mm256_insertf128_pd( t7, _mm_load_pd(&buf[34]), 1); \
  t8  = _mm256_insertf128_pd( t8, _mm_load_pd(&buf[36]), 1); \
  t9  = _mm256_insertf128_pd( t9, _mm_load_pd(&buf[38]), 1); \
  \
  vx0 = _mm256_unpacklo_pd( t0, t5 ); \
  vx1 = _mm256_unpackhi_pd( t0, t5 ); \
  vy0 = _mm256_unpacklo_pd( t1, t6 ); \
  vy1 = _mm256_unpackhi_pd( t1, t6 ); \
  vz0 = _mm256_unpacklo_pd( t2, t7 ); \
  vz1 = _mm256_unpackhi_pd( t2, t7 ); \
  vq  = _mm256_unpacklo_pd( t3, t8 ); \
  \
  __m128 a, b; \
  t0 = _mm256_unpacklo_pd( t4, t9 ); \
  a = _mm_castpd_ps( _mm256_castpd256_pd128( t0 ) ); \
  b = _mm_castpd_ps( _mm256_extractf128_pd( t0, 1 ) ); \
  vix = _mm_castps_si128( _mm_shuffle_ps( a, b, _MM_SHUFFLE( 2, 0, 2, 0 ) ) ); \
  viy = _mm_castps_si128( _mm_shuffle_ps( a, b, _MM_SHUFFLE( 3, 1, 3, 1 ) ) ); \
  \
  t1 = _mm256_unpackhi_pd( t4, t9 ); \
  a = _mm_castpd_ps( _mm256_castpd256_pd128( t1 ) ); \
  b = _mm_castpd_ps( _mm256_extractf128_pd( t1, 1 ) ); \
  viz = _mm_castps_si128( _mm_shuffle_ps( a, b, _MM_SHUFFLE( 2, 0, 2, 0 ) ) ); \
}


// Virtual particle buffer
typedef struct { 
  // 4 splits maximum
  DECLARE_ALIGNED_32( t_vp3D buf[ 4 * p_cache_size_3D ] );
  double *p;
  unsigned int np;
} t_vpbuf3D;


#endif
