#ifndef _OS_SPEC_PUSH_H_AVX
#define _OS_SPEC_PUSH_H_AVX

/* Size of particle buffer */
/* for values between 32 and 8192 this has little impact on performance */
#define p_cache_size_2D 64
#define p_cache_size_3D 1024

#include "vector-avx.h"

/*****************************************************************************************
t_vp2D

Virtual particle type for current deposition in 2D. This corresponds to 2 128bit
vectors in the form:

 | x0, x1, y0, y1| q, vz, ix, iy |

The pusher will store one vp per simulation particle; the splitter will
then add new vp as necessary / correct the previous one.

*****************************************************************************************/

typedef struct virtualPart2D {
	float x0, x1, y0, y1;
	float q, vz;
	int ix, iy;
} t_vp2D;

/*****************************************************************************************
STOREVP2D

Stores 1 virtual particle for 2D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
*****************************************************************************************/
#define STOREVP2D( buf, x0, x1, y0, y1, q, vz, ix, iy ) { 			        \
  buf -> x0 = x0;                                                           \
  buf -> x1 = x1;                                                           \
  buf -> y0 = y0;                                                           \
  buf -> y1 = y1;                                                           \
  buf -> q  = q;                                                            \
  buf -> vz = vz;                                                           \
  buf -> ix = ix;                                                           \
  buf -> iy = iy;                                                           \
}

/*****************************************************************************************
STORE8VP2D

Stores 8 virtual particles for 2D current deposition. If the particle did not cross any
boundary, the splitter will not do any additional work. 

Please note that particles are stored in an interleaved order namely:
 0,2,4,6,1,3,5,7 
This means that the "cross" vector must also be shuffled.
*****************************************************************************************/

#if 0

/*
This version minimizes register use.

Please note that particles are stored in an interleaved order namely:
 0,2,4,6,1,3,5,7 
 
This means that the "cross" vector must also be shuffled.
*/

#define STORE8VP2D( f, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { 			             \
   register __m256 a, b, t0, t1, t2, t3;                                                 \
																		                 \
   t0 = _mm256_unpacklo_ps( vx0, vx1 );											         \
   t1 = _mm256_unpacklo_ps( vy0, vy1 );											         \
   t2 = _mm256_unpacklo_ps( vq, vvz );											         \
   t3 = _mm256_unpacklo_ps( _mm256_castsi256_ps(vix), _mm256_castsi256_ps(viy) );		 \
                                                                                         \
   a = _mm256_shuffle_ps( t0, t1, _MM_SHUFFLE( 1, 0, 1, 0 ) );                           \
   b = _mm256_shuffle_ps( t2, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) );                           \
                                                                                         \
   _mm_store_ps( &f[ 0], _mm256_castps256_ps128( a ) );                                  \
   _mm_store_ps( &f[ 4], _mm256_castps256_ps128( b ) );                                  \
                                                                                         \
   _mm_store_ps( &f[ 8], _mm256_extractf128_ps( a, 1 ) );                                \
   _mm_store_ps( &f[12], _mm256_extractf128_ps( b, 1 ) );                                \
                                                                                         \
   a = _mm256_shuffle_ps( t0, t1, _MM_SHUFFLE( 3, 2, 3, 2 ) );                           \
   b = _mm256_shuffle_ps( t2, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) );                           \
                                                                                         \
   _mm_store_ps( &f[16], _mm256_castps256_ps128( a ) );                                  \
   _mm_store_ps( &f[20], _mm256_castps256_ps128( b ) );                                  \
                                                                                         \
   _mm_store_ps( &f[24], _mm256_extractf128_ps( a, 1 ) );                                \
   _mm_store_ps( &f[28], _mm256_extractf128_ps( b, 1 ) );                                \
                                                                                         \
   t0 = _mm256_unpackhi_ps( vx0, vx1 );											         \
   t1 = _mm256_unpackhi_ps( vy0, vy1 );											         \
   t2 = _mm256_unpackhi_ps( vq, vvz );											         \
   t3 = _mm256_unpackhi_ps( _mm256_castsi256_ps(vix), _mm256_castsi256_ps(viy) );		 \
                                                                                         \
   a = _mm256_shuffle_ps( t0, t1, _MM_SHUFFLE( 1, 0, 1, 0 ) );                           \
   b = _mm256_shuffle_ps( t2, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) );                           \
                                                                                         \
   _mm_store_ps( &f[32], _mm256_castps256_ps128( a ) );                                  \
   _mm_store_ps( &f[36], _mm256_castps256_ps128( b ) );                                  \
                                                                                         \
   _mm_store_ps( &f[40], _mm256_extractf128_ps( a, 1 ) );                                \
   _mm_store_ps( &f[44], _mm256_extractf128_ps( b, 1 ) );                                \
                                                                                         \
   a = _mm256_shuffle_ps( t0, t1, _MM_SHUFFLE( 3, 2, 3, 2 ) );                           \
   b = _mm256_shuffle_ps( t2, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) );                           \
                                                                                         \
   _mm_store_ps( &f[48], _mm256_castps256_ps128( a ) );                                  \
   _mm_store_ps( &f[52], _mm256_castps256_ps128( b ) );                                  \
                                                                                         \
   _mm_store_ps( &f[56], _mm256_extractf128_ps( a, 1 ) );                                \
   _mm_store_ps( &f[60], _mm256_extractf128_ps( b, 1 ) );                                \
}

#endif

#define STORE8VP2D( f, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
   __m256 a, b, c, d; \
   __m256 t0, t1, t2, t3, t4, t5, t6, t7; \
   \
   a  = _mm256_unpacklo_ps( vx0, vx1 ); \
   b  = _mm256_unpacklo_ps( vy0, vy1 ); \
   c  = _mm256_unpacklo_ps( vq, vvz ); \
   d  = _mm256_unpacklo_ps( _mm256_castsi256_ps(vix), _mm256_castsi256_ps(viy) ); \
   \
   t0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t1 = _mm256_shuffle_ps( c, d, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t2 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE( 3, 2, 3, 2 ) ); \
   t3 = _mm256_shuffle_ps( c, d, _MM_SHUFFLE( 3, 2, 3, 2 ) ); \
   \
   a  = _mm256_unpackhi_ps( vx0, vx1 ); \
   b  = _mm256_unpackhi_ps( vy0, vy1 ); \
   c  = _mm256_unpackhi_ps( vq, vvz ); \
   d  = _mm256_unpackhi_ps( _mm256_castsi256_ps(vix), _mm256_castsi256_ps(viy) ); \
   \
   t4 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t5 = _mm256_shuffle_ps( c, d, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t6 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE( 3, 2, 3, 2 ) ); \
   t7 = _mm256_shuffle_ps( c, d, _MM_SHUFFLE( 3, 2, 3, 2 ) ); \
   \
   _mm_store_ps( &f[ 0], _mm256_castps256_ps128( t0 ) ); \
   _mm_store_ps( &f[ 4], _mm256_castps256_ps128( t1 ) ); \
   _mm_store_ps( &f[ 8], _mm256_castps256_ps128( t2 ) ); \
   _mm_store_ps( &f[12], _mm256_castps256_ps128( t3 ) ); \
   \
   _mm_store_ps( &f[16], _mm256_castps256_ps128( t4 ) ); \
   _mm_store_ps( &f[20], _mm256_castps256_ps128( t5 ) ); \
   _mm_store_ps( &f[24], _mm256_castps256_ps128( t6 ) ); \
   _mm_store_ps( &f[28], _mm256_castps256_ps128( t7 ) ); \
   \
   _mm_store_ps( &f[32], _mm256_extractf128_ps( t0, 1 ) ); \
   _mm_store_ps( &f[36], _mm256_extractf128_ps( t1, 1 ) ); \
   _mm_store_ps( &f[40], _mm256_extractf128_ps( t2, 1 ) ); \
   _mm_store_ps( &f[44], _mm256_extractf128_ps( t3, 1 ) ); \
   \
   _mm_store_ps( &f[48], _mm256_extractf128_ps( t4, 1 ) ); \
   _mm_store_ps( &f[52], _mm256_extractf128_ps( t5, 1 ) ); \
   _mm_store_ps( &f[56], _mm256_extractf128_ps( t6, 1 ) ); \
   _mm_store_ps( &f[60], _mm256_extractf128_ps( t7, 1 ) ); \
   \
}


/*****************************************************************************************
LOAD8VP2D

Loads 8 virtual particles for 2D current deposition. 

*****************************************************************************************/

#if 0

/*
This version minimizes register use.

Please note that particles are stored in an interleaved order namely:
 0,2,4,6,1,3,5,7 
 
This means that the "cross" vector must also be shuffled.
*/


#define LOAD8VP2D( f, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { 			             \
   register __m256 a, b, t1, t3, t5, t7;                                                 \
		                                                                                 \
   vy0 = _mm256_castps128_ps256(  _mm_load_ps(&f[ 0]));                                  \
   t1  = _mm256_castps128_ps256(  _mm_load_ps(&f[ 4]));                                  \
   vy0 = _mm256_insertf128_ps(vy0,_mm_load_ps(&f[ 8]),1);                                \
   t1  = _mm256_insertf128_ps(t1, _mm_load_ps(&f[12]),1);                                \
                                                                                         \
   vy1 = _mm256_castps128_ps256(  _mm_load_ps(&f[16]));                                  \
   t3  = _mm256_castps128_ps256(  _mm_load_ps(&f[20]));                                  \
   vy1 = _mm256_insertf128_ps(vy1,_mm_load_ps(&f[24]),1);                                \
   t3  = _mm256_insertf128_ps(t3, _mm_load_ps(&f[28]),1);                                \
                                                                                         \
   vq = _mm256_castps128_ps256(   _mm_load_ps(&f[32]));                                  \
   t5 = _mm256_castps128_ps256(   _mm_load_ps(&f[36]));                                  \
   vq = _mm256_insertf128_ps(vq,  _mm_load_ps(&f[40]),1);                                \
   t5 = _mm256_insertf128_ps(t5,  _mm_load_ps(&f[44]),1);                                \
                                                                                         \
   vvz = _mm256_castps128_ps256(  _mm_load_ps(&f[48]));                                  \
   t7  = _mm256_castps128_ps256(  _mm_load_ps(&f[52]));                                  \
   vvz = _mm256_insertf128_ps(vvz,_mm_load_ps(&f[56]),1);                                \
   t7  = _mm256_insertf128_ps(t7, _mm_load_ps(&f[60]),1);                                \
                                                                                         \
   a = _mm256_shuffle_ps( vy0, vy1, _MM_SHUFFLE(1,0,1,0) );                              \
   b = _mm256_shuffle_ps(  vq, vvz, _MM_SHUFFLE(1,0,1,0) );                              \
                                                                                         \
   vx0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(2,0,2,0) );                                \
   vx1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,1,3,1) );                                \
		                                                                                 \
   a = _mm256_shuffle_ps( vy0, vy1, _MM_SHUFFLE(3,2,3,2) );                              \
   b = _mm256_shuffle_ps(  vq, vvz, _MM_SHUFFLE(3,2,3,2) );                              \
                                                                                         \
   vy0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(2,0,2,0) );                                \
   vy1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,1,3,1) );                                \
                                                                                         \
   a = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE(1,0,1,0) );                                \
   b = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE(1,0,1,0) );                                \
                                                                                         \
   vq  =  _mm256_shuffle_ps( a, b, _MM_SHUFFLE(2,0,2,0) );                               \
   vvz =  _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,1,3,1) );                               \
		                                                                                 \
   a = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE(3,2,3,2) );                                \
   b = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE(3,2,3,2) );                                \
                                                                                         \
   vix = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(2,0,2,0) ));          \
   viy = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,1,3,1) ));          \
}

#endif

#define LOAD8VP2D( f, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
   __m256 t0, t1, t2, t3, t4, t5, t6, t7; \
   __m256 a, b; \
   \
   t0  = _mm256_castps128_ps256(  _mm_load_ps(&f[ 0])); \
   t4  = _mm256_castps128_ps256(  _mm_load_ps(&f[ 4])); \
   t1  = _mm256_castps128_ps256(  _mm_load_ps(&f[ 8])); \
   t5  = _mm256_castps128_ps256(  _mm_load_ps(&f[12])); \
   t2  = _mm256_castps128_ps256(  _mm_load_ps(&f[16])); \
   t6  = _mm256_castps128_ps256(  _mm_load_ps(&f[20])); \
   t3  = _mm256_castps128_ps256(  _mm_load_ps(&f[24])); \
   t7  = _mm256_castps128_ps256(  _mm_load_ps(&f[28])); \
   \
   t0  = _mm256_insertf128_ps( t0, _mm_load_ps(&f[32]), 1); \
   t4  = _mm256_insertf128_ps( t4, _mm_load_ps(&f[36]), 1); \
   t1  = _mm256_insertf128_ps( t1, _mm_load_ps(&f[40]), 1); \
   t5  = _mm256_insertf128_ps( t5, _mm_load_ps(&f[44]), 1); \
   t2  = _mm256_insertf128_ps( t2, _mm_load_ps(&f[48]), 1); \
   t6  = _mm256_insertf128_ps( t6, _mm_load_ps(&f[52]), 1); \
   t3  = _mm256_insertf128_ps( t3, _mm_load_ps(&f[56]), 1); \
   t7  = _mm256_insertf128_ps( t7, _mm_load_ps(&f[60]), 1); \
   \
   a = _mm256_unpacklo_ps( t0, t1 ); \
   b = _mm256_unpacklo_ps( t2, t3 ); \
   \
   vx0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(1,0,1,0) ); \
   vx1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,2,3,2) ); \
   \
   a = _mm256_unpackhi_ps( t0, t1 ); \
   b = _mm256_unpackhi_ps( t2, t3 ); \
   \
   vy0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(1,0,1,0) ); \
   vy1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,2,3,2) ); \
   \
   a = _mm256_unpacklo_ps( t4, t5 ); \
   b = _mm256_unpacklo_ps( t6, t7 ); \
   \
   vq = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(1,0,1,0) ); \
   vvz = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,2,3,2) ); \
   \
   a = _mm256_unpackhi_ps( t4, t5 ); \
   b = _mm256_unpackhi_ps( t6, t7 ); \
   \
   vix = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(1,0,1,0) ) ); \
   viy = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,2,3,2) ) ); \
 }

/***********************************************************************
t_vpbuf2D

Buffer to hold virtual particles for current deposition in 2D. 

************************************************************************/

typedef struct { 
  // 3 splits maximum
  DECLARE_ALIGNED_32( t_vp2D buf[ 3 * p_cache_size_2D ] );
  float *p;
  unsigned int np;
} t_vpbuf2D;


/***********************************************************************
t_vp3D

Virtual particle type for current deposition. This corresponds to 3
vectors in the form:

 | x0, x1, y0, y1| z0, z1, q, ix | iy, iz, p0, p1 |

p0 and p1 are used for padding to improve memory alignment

************************************************************************/

typedef struct virtualPart3D {
	float x0,     x1,    y0,     y1;
	float z0,     z1,     q; int ix;
	int   iy,     iz,    p0,     p1;
} t_vp3D;

/***********************************************************************
STORE8VP3D

Stores 8 virtual particles for 3D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
************************************************************************/

#if 0

/*
This version minimizes register use.

Please note that particles are stored in an interleaved order namely:
 0,2,4,6,1,3,5,7 
 
This means that the "cross" vector must also be shuffled.
*/


#define STORE8VP3D( f, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { 	          \
   register __m256 a, b, c, t0, t1, t2, t3, t4;                                       \
																					  \
   t0 = _mm256_unpacklo_ps( vx0, vx1 );											      \
   t1 = _mm256_unpacklo_ps( vy0, vy1 );											      \
   t2 = _mm256_unpacklo_ps( vz0, vz1 );											      \
   t3 = _mm256_unpacklo_ps( vq, _mm256_castsi256_ps(vix) );                           \
   t4 = _mm256_unpacklo_ps( _mm256_castsi256_ps(viy), _mm256_castsi256_ps(viz) );     \
																					  \
   a = _mm256_shuffle_ps( t0, t1, _MM_SHUFFLE( 1, 0, 1, 0 ) );                        \
   b = _mm256_shuffle_ps( t2, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) );                        \
   c = _mm256_shuffle_ps( t4, t4, _MM_SHUFFLE( 1, 0, 1, 0 ) );                        \
																					  \
   _mm_store_ps( &f[ 0], _mm256_castps256_ps128( a ) );                               \
   _mm_store_ps( &f[ 4], _mm256_castps256_ps128( b ) );                               \
   _mm_store_ps( &f[ 8], _mm256_castps256_ps128( c ) );                               \
																					  \
   _mm_store_ps( &f[12], _mm256_extractf128_ps( a, 1 ) );                             \
   _mm_store_ps( &f[16], _mm256_extractf128_ps( b, 1 ) );                             \
   _mm_store_ps( &f[20], _mm256_extractf128_ps( c, 1 ) );                             \
																					  \
   a = _mm256_shuffle_ps( t0, t1, _MM_SHUFFLE( 3, 2, 3, 2 ) );                        \
   b = _mm256_shuffle_ps( t2, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) );                        \
   c = _mm256_shuffle_ps( t4, t4, _MM_SHUFFLE( 3, 2, 3, 2 ) );                        \
																					  \
   _mm_store_ps( &f[24], _mm256_castps256_ps128( a ) );                               \
   _mm_store_ps( &f[28], _mm256_castps256_ps128( b ) );                               \
   _mm_store_ps( &f[32], _mm256_castps256_ps128( c ) );                               \
																					  \
   _mm_store_ps( &f[36], _mm256_extractf128_ps( a, 1 ) );                             \
   _mm_store_ps( &f[40], _mm256_extractf128_ps( b, 1 ) );                             \
   _mm_store_ps( &f[44], _mm256_extractf128_ps( c, 1 ) );                             \
																					  \
   t0 = _mm256_unpackhi_ps( vx0, vx1 );											      \
   t1 = _mm256_unpackhi_ps( vy0, vy1 );											      \
   t2 = _mm256_unpackhi_ps( vz0, vz1 );											      \
   t3 = _mm256_unpackhi_ps( vq, _mm256_castsi256_ps(vix) );                           \
   t4 = _mm256_unpackhi_ps( _mm256_castsi256_ps(viy), _mm256_castsi256_ps(viz) );     \
																					  \
   a = _mm256_shuffle_ps( t0, t1, _MM_SHUFFLE( 1, 0, 1, 0 ) );                        \
   b = _mm256_shuffle_ps( t2, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) );                        \
   c = _mm256_shuffle_ps( t4, t4, _MM_SHUFFLE( 1, 0, 1, 0 ) );                        \
																					  \
   _mm_store_ps( &f[48], _mm256_castps256_ps128( a ) );                               \
   _mm_store_ps( &f[52], _mm256_castps256_ps128( b ) );                               \
   _mm_store_ps( &f[56], _mm256_castps256_ps128( c ) );                               \
																					  \
   _mm_store_ps( &f[60], _mm256_extractf128_ps( a, 1 ) );                             \
   _mm_store_ps( &f[64], _mm256_extractf128_ps( b, 1 ) );                             \
   _mm_store_ps( &f[68], _mm256_extractf128_ps( c, 1 ) );                             \
																					  \
   a = _mm256_shuffle_ps( t0, t1, _MM_SHUFFLE( 3, 2, 3, 2 ) );                        \
   b = _mm256_shuffle_ps( t2, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) );                        \
   c = _mm256_shuffle_ps( t4, t4, _MM_SHUFFLE( 3, 2, 3, 2 ) );                        \
																					  \
   _mm_store_ps( &f[72], _mm256_castps256_ps128( a ) );                               \
   _mm_store_ps( &f[76], _mm256_castps256_ps128( b ) );                               \
   _mm_store_ps( &f[80], _mm256_castps256_ps128( c ) );                               \
																					  \
   _mm_store_ps( &f[84], _mm256_extractf128_ps( a, 1 ) );                             \
   _mm_store_ps( &f[88], _mm256_extractf128_ps( b, 1 ) );                             \
   _mm_store_ps( &f[92], _mm256_extractf128_ps( c, 1 ) );                             \
}

#endif

#define STORE8VP3D( f, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
   __m256 a, b, c, d, e; \
   __m256 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11; \
   \
   a  = _mm256_unpacklo_ps( vx0, vx1 ); \
   b  = _mm256_unpacklo_ps( vy0, vy1 ); \
   c  = _mm256_unpacklo_ps( vz0, vz1 ); \
   d  = _mm256_unpacklo_ps( vq, _mm256_castsi256_ps(vix) ); \
   e =  _mm256_unpacklo_ps( _mm256_castsi256_ps(viy), _mm256_castsi256_ps(viz) ); \
   \
   t0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t1 = _mm256_shuffle_ps( c, d, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t2 = _mm256_shuffle_ps( e, e, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t3 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE( 3, 2, 3, 2 ) ); \
   t4 = _mm256_shuffle_ps( c, d, _MM_SHUFFLE( 3, 2, 3, 2 ) ); \
   t5 = _mm256_shuffle_ps( e, e, _MM_SHUFFLE( 3, 2, 3, 2 ) );\
   \
   a  = _mm256_unpackhi_ps( vx0, vx1 ); \
   b  = _mm256_unpackhi_ps( vy0, vy1 ); \
   c  = _mm256_unpackhi_ps( vz0, vz1 ); \
   d  = _mm256_unpackhi_ps( vq, _mm256_castsi256_ps(vix) ); \
   e =  _mm256_unpackhi_ps( _mm256_castsi256_ps(viy), _mm256_castsi256_ps(viz) ); \
   \
   t6 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t7 = _mm256_shuffle_ps( c, d, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t8 = _mm256_shuffle_ps( e, e, _MM_SHUFFLE( 1, 0, 1, 0 ) ); \
   t9 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE( 3, 2, 3, 2 ) ); \
   t10= _mm256_shuffle_ps( c, d, _MM_SHUFFLE( 3, 2, 3, 2 ) ); \
   t11= _mm256_shuffle_ps( e, e, _MM_SHUFFLE( 3, 2, 3, 2 ) ); \
    \
   _mm_store_ps( &f[ 0], _mm256_castps256_ps128( t0 ) ); \
   _mm_store_ps( &f[ 4], _mm256_castps256_ps128( t1 ) ); \
   _mm_store_ps( &f[ 8], _mm256_castps256_ps128( t2 ) ); \
   \
   _mm_store_ps( &f[12], _mm256_castps256_ps128( t3 ) ); \
   _mm_store_ps( &f[16], _mm256_castps256_ps128( t4 ) ); \
   _mm_store_ps( &f[20], _mm256_castps256_ps128( t5 ) ); \
   \
   _mm_store_ps( &f[24], _mm256_castps256_ps128( t6 ) ); \
   _mm_store_ps( &f[28], _mm256_castps256_ps128( t7 ) ); \
   _mm_store_ps( &f[32], _mm256_castps256_ps128( t8 ) ); \
   \
   _mm_store_ps( &f[36], _mm256_castps256_ps128( t9 ) ); \
   _mm_store_ps( &f[40], _mm256_castps256_ps128( t10 ) ); \
   _mm_store_ps( &f[44], _mm256_castps256_ps128( t11 ) ); \
   \
   _mm_store_ps( &f[48], _mm256_extractf128_ps( t0, 1 ) ); \
   _mm_store_ps( &f[52], _mm256_extractf128_ps( t1, 1 ) ); \
   _mm_store_ps( &f[56], _mm256_extractf128_ps( t2, 1 ) ); \
   \
   _mm_store_ps( &f[60], _mm256_extractf128_ps( t3, 1 ) ); \
   _mm_store_ps( &f[64], _mm256_extractf128_ps( t4, 1 ) ); \
   _mm_store_ps( &f[68], _mm256_extractf128_ps( t5, 1 ) ); \
   \
   _mm_store_ps( &f[72], _mm256_extractf128_ps( t6, 1 ) ); \
   _mm_store_ps( &f[76], _mm256_extractf128_ps( t7, 1 ) ); \
   _mm_store_ps( &f[80], _mm256_extractf128_ps( t8, 1 ) ); \
   \
   _mm_store_ps( &f[84], _mm256_extractf128_ps( t9, 1 ) ); \
   _mm_store_ps( &f[88], _mm256_extractf128_ps( t10, 1 ) ); \
   _mm_store_ps( &f[92], _mm256_extractf128_ps( t11, 1 ) ); \
 }

/***********************************************************************
LOAD8VP3D

Loads 8 virtual particles for 3D current deposition. 

************************************************************************/

#if 0

/*
This version minimizes register use.

Please note that particles are stored in an interleaved order namely:
 0,2,4,6,1,3,5,7 
 
This means that the "cross" vector must also be shuffled.
*/


#define LOAD8VP3D( f, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { 		      \
   register __m256 a, b, c, t2, t4, t5, t7, t8, t10, t11;                             \
																					  \
   vy0 = _mm256_castps128_ps256(  _mm_load_ps(&f[ 0]));                               \
   vq = _mm256_castps128_ps256(  _mm_load_ps(&f[ 4]));                                \
   t2 = _mm256_castps128_ps256(  _mm_load_ps(&f[ 8]));                                \
   vy0 = _mm256_insertf128_ps(vy0, _mm_load_ps(&f[12]),1);                            \
   vq = _mm256_insertf128_ps(vq, _mm_load_ps(&f[16]),1);                              \
   t2 = _mm256_insertf128_ps(t2, _mm_load_ps(&f[20]),1);                              \
																					  \
   vy1 = _mm256_castps128_ps256(  _mm_load_ps(&f[24]));                               \
   t4 = _mm256_castps128_ps256(  _mm_load_ps(&f[28]));                                \
   t5 = _mm256_castps128_ps256(  _mm_load_ps(&f[32]));                                \
   vy1 = _mm256_insertf128_ps(vy1, _mm_load_ps(&f[36]),1);                            \
   t4 = _mm256_insertf128_ps(t4, _mm_load_ps(&f[40]),1);                              \
   t5 = _mm256_insertf128_ps(t5, _mm_load_ps(&f[44]),1);                              \
																					  \
   vz0 = _mm256_castps128_ps256(  _mm_load_ps(&f[48]));                               \
   t7 = _mm256_castps128_ps256(  _mm_load_ps(&f[52]));                                \
   t8 = _mm256_castps128_ps256(  _mm_load_ps(&f[56]));                                \
   vz0 = _mm256_insertf128_ps(vz0, _mm_load_ps(&f[60]),1);                            \
   t7 = _mm256_insertf128_ps(t7, _mm_load_ps(&f[64]),1);                              \
   t8 = _mm256_insertf128_ps(t8, _mm_load_ps(&f[68]),1);                              \
																					  \
   vz1  = _mm256_castps128_ps256(  _mm_load_ps(&f[72]));                              \
   t10 = _mm256_castps128_ps256(  _mm_load_ps(&f[76]));                               \
   t11 = _mm256_castps128_ps256(  _mm_load_ps(&f[80]));                               \
   vz1  = _mm256_insertf128_ps( vz1, _mm_load_ps(&f[84]),1);                          \
   t10 = _mm256_insertf128_ps(t10, _mm_load_ps(&f[88]),1);                            \
   t11 = _mm256_insertf128_ps(t11, _mm_load_ps(&f[92]),1);                            \
																					  \
   a = _mm256_shuffle_ps( vy0, vy1, _MM_SHUFFLE(1,0,1,0) );                           \
   b = _mm256_shuffle_ps( vz0, vz1, _MM_SHUFFLE(1,0,1,0) );                           \
																					  \
   vx0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(2,0,2,0) );                             \
   vx1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,1,3,1) );                             \
																					  \
   a = _mm256_shuffle_ps( vy0, vy1, _MM_SHUFFLE(3,2,3,2) );                           \
   b = _mm256_shuffle_ps( vz0, vz1, _MM_SHUFFLE(3,2,3,2) );                           \
																					  \
   vy0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(2,0,2,0) );                             \
   vy1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,1,3,1) );                             \
																					  \
   a = _mm256_shuffle_ps( vq, t4, _MM_SHUFFLE(1,0,1,0) );                             \
   b = _mm256_shuffle_ps( t7, t10, _MM_SHUFFLE(1,0,1,0) );                            \
																					  \
   vz0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(2,0,2,0) );                             \
   vz1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,1,3,1) );                             \
																					  \
   a = _mm256_shuffle_ps( vq, t4, _MM_SHUFFLE(3,2,3,2) );                             \
   b = _mm256_shuffle_ps( t7, t10, _MM_SHUFFLE(3,2,3,2) );                            \
																					  \
   vq  = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(2,0,2,0) );                             \
   vix = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,1,3,1) ) );      \
																					  \
   a = _mm256_shuffle_ps( t2, t5, _MM_SHUFFLE(1,0,1,0) );                             \
   b = _mm256_shuffle_ps( t8, t11, _MM_SHUFFLE(1,0,1,0) );                            \
																					  \
   viy = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(2,0,2,0) ) );      \
   viz = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,1,3,1) ) );      \
}

#endif

#define LOAD8VP3D( f, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
   __m256 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11; \
   __m256 a, b; \
   \
   t0  = _mm256_castps128_ps256(  _mm_load_ps(&f[ 0])); \
   t4  = _mm256_castps128_ps256(  _mm_load_ps(&f[ 4])); \
   t8  = _mm256_castps128_ps256(  _mm_load_ps(&f[ 8])); \
   t1  = _mm256_castps128_ps256(  _mm_load_ps(&f[12])); \
   t5  = _mm256_castps128_ps256(  _mm_load_ps(&f[16])); \
   t9  = _mm256_castps128_ps256(  _mm_load_ps(&f[20])); \
   t2  = _mm256_castps128_ps256(  _mm_load_ps(&f[24])); \
   t6  = _mm256_castps128_ps256(  _mm_load_ps(&f[28])); \
   t10 = _mm256_castps128_ps256(  _mm_load_ps(&f[32])); \
   t3  = _mm256_castps128_ps256(  _mm_load_ps(&f[36])); \
   t7  = _mm256_castps128_ps256(  _mm_load_ps(&f[40])); \
   t11 = _mm256_castps128_ps256(  _mm_load_ps(&f[44])); \
   \
   t0  = _mm256_insertf128_ps( t0  , _mm_load_ps(&f[48]), 1); \
   t4  = _mm256_insertf128_ps( t4  , _mm_load_ps(&f[52]), 1); \
   t8  = _mm256_insertf128_ps( t8  , _mm_load_ps(&f[56]), 1); \
   t1  = _mm256_insertf128_ps( t1  , _mm_load_ps(&f[60]), 1); \
   t5  = _mm256_insertf128_ps( t5  , _mm_load_ps(&f[64]), 1); \
   t9  = _mm256_insertf128_ps( t9  , _mm_load_ps(&f[68]), 1); \
   t2  = _mm256_insertf128_ps( t2  , _mm_load_ps(&f[72]), 1); \
   t6  = _mm256_insertf128_ps( t6  , _mm_load_ps(&f[76]), 1); \
   t10 = _mm256_insertf128_ps( t10 , _mm_load_ps(&f[80]), 1); \
   t3  = _mm256_insertf128_ps( t3  , _mm_load_ps(&f[84]), 1); \
   t7  = _mm256_insertf128_ps( t7  , _mm_load_ps(&f[88]), 1); \
   t11 = _mm256_insertf128_ps( t11 , _mm_load_ps(&f[92]), 1); \
   \
   a = _mm256_unpacklo_ps( t0, t1 ); \
   b = _mm256_unpacklo_ps( t2, t3 ); \
   \
   vx0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(1,0,1,0) ); \
   vx1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,2,3,2) ); \
   \
   a = _mm256_unpackhi_ps( t0, t1 ); \
   b = _mm256_unpackhi_ps( t2, t3 ); \
   \
   vy0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(1,0,1,0) ); \
   vy1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,2,3,2) ); \
   \
   a = _mm256_unpacklo_ps( t4, t5 ); \
   b = _mm256_unpacklo_ps( t6, t7 ); \
   \
   vz0 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(1,0,1,0) ); \
   vz1 = _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,2,3,2) ); \
   \
   a = _mm256_unpackhi_ps( t4, t5 ); \
   b = _mm256_unpackhi_ps( t6, t7 ); \
   \
   vq  =                      _mm256_shuffle_ps( a, b, _MM_SHUFFLE(1,0,1,0) ) ; \
   vix = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,2,3,2) ) ); \
   \
   a = _mm256_unpacklo_ps( t8, t9 ); \
   b = _mm256_unpacklo_ps( t10, t11 ); \
   \
   viy = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(1,0,1,0) ) ); \
   viz = _mm256_castps_si256( _mm256_shuffle_ps( a, b, _MM_SHUFFLE(3,2,3,2) ) ); \
 }


/*****************************************************************************************
_MM_STORE_INTERLEAVE

Stores a __m256 vector to memory interleaving lower half / upper half elements for
compatibility with the (LOAD/STORE)VP(2/3)D routines. Values are stored in the following
order:

v = [v0, v1, v2, v3, v4, v5, v6, v7 ]

-> _MM_STORE_INTERLEAVE_PS( b, v )
-> b = [ v0, v4, v1, v5, v2, v6, v3, v7 ]

*****************************************************************************************/

/*

// Currently disabled

#define _MM_STORE_INTERLEAVE( b, x ) {   \
  __m128 t0 = _mm256_castps256_ps128( x ); \
  __m128 t1 = _mm256_extractf128_ps( x, 1 ); \
  _mm_store_ps( (float*) &b[0], _mm_unpacklo_ps( t0, t1 ) ); \
  _mm_store_ps( (float*) &b[4], _mm_unpackhi_ps( t0, t1 ) ); \
}

*/

// Virtual particle buffer
typedef struct { 
  // 4 splits maximum
  DECLARE_ALIGNED_32( t_vp3D buf[ 4 * p_cache_size_3D ] );
  float *p;
  unsigned int np;
} t_vpbuf3D;

#endif
