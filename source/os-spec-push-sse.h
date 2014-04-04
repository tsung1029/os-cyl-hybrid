#ifndef _OS_SPEC_PUSH_H_SSE
#define _OS_SPEC_PUSH_H_SSE

#ifdef __SSE__

/* Size of particle buffer */
/* for values between 32 and 8192 this has little impact on performance */
#define p_cache_size_2D 64
#define p_cache_size_3D 1024

#include "vector-sse.h"

/***********************************************************************
t_vp2D

Virtual particle type for current deposition in 2D. This corresponds to 2
vectors in the form:

 | x0, x1, y0, y1| q, vz, ix, iy |

The pusher will store one vp per simulation particle; the splitter will
then add new vp as necessary / correct the previous one.

************************************************************************/

typedef struct virtualPart2D {
	float x0, x1, y0, y1;
	float q, vz;
	int ix, iy;
} t_vp2D;

/***********************************************************************
STOREVP2D

Stores 1 virtual particle for 2D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
************************************************************************/
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


/***********************************************************************
STORE4VP2D

Stores 4 virtual particles for 2D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
************************************************************************/

#define STORE4VP2D( buf, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { 			\
  register __m128 t0, t1, t2, t3;											\
  																			\
  t0 = _mm_unpacklo_ps( vx0, vx1 );											\
  t1 = _mm_unpacklo_ps( vy0, vy1 );											\
  t2 = _mm_unpacklo_ps( vq, vvz );											\
  t3 = _mm_unpacklo_ps( _mm_castsi128_ps(vix), _mm_castsi128_ps(viy) );		\
  _mm_store_ps( buf    , _mm_movelh_ps( t0, t1 ) );							\
  _mm_store_ps( buf + 4, _mm_movelh_ps( t2, t3 ) );							\
  _mm_store_ps( buf + 8, _mm_movehl_ps( t1, t0 ) );							\
  _mm_store_ps( buf +12, _mm_movehl_ps( t3, t2 ) );							\
																			\
  t0 = _mm_unpackhi_ps( vx0, vx1 );											\
  t1 = _mm_unpackhi_ps( vy0, vy1 );											\
  t2 = _mm_unpackhi_ps( vq, vvz );											\
  t3 = _mm_unpackhi_ps( _mm_castsi128_ps(vix), _mm_castsi128_ps(viy) );		\
  _mm_store_ps( buf +16, _mm_movelh_ps( t0, t1 ) );							\
  _mm_store_ps( buf +20, _mm_movelh_ps( t2, t3 ) );							\
  _mm_store_ps( buf +24, _mm_movehl_ps( t1, _mm_movehl_ps(t0,t0) ) );		\
  _mm_store_ps( buf +28, _mm_movehl_ps( t3, _mm_movehl_ps(t2,t2) ) );		\
}

/***********************************************************************
LOAD4VP2D

Loads 4 virtual particles for 2D current deposition. 

************************************************************************/

#define LOAD4VP2D( buf, vx0, vx1, vy0, vy1, vq, vvz, ix, iy ) { 			\
	register __m128 t0, t1, t2, t3, vix, viy;								\
																			\
	vx0  = _mm_load_ps( buf + 0  );											\
	vq   = _mm_load_ps( buf + 4  );											\
	vx1  = _mm_load_ps( buf + 8  );											\
	vvz  = _mm_load_ps( buf + 12 );											\
	vy0  = _mm_load_ps( buf + 16 );											\
	vix  = _mm_load_ps( buf + 20 );											\
	vy1  = _mm_load_ps( buf + 24 );											\
	viy  = _mm_load_ps( buf + 28 );											\
																			\
	t0 = _mm_unpacklo_ps(vx0, vx1);											\
	t2 = _mm_unpacklo_ps(vy0, vy1);											\
	t1 = _mm_unpackhi_ps(vx0, vx1);											\
	t3 = _mm_unpackhi_ps(vy0, vy1);											\
																			\
	vx0  = _mm_movelh_ps(t0, t2);											\
	vx1  = _mm_movehl_ps(t2, t0);											\
	vy0  = _mm_movelh_ps(t1, t3);											\
	vy1  = _mm_movehl_ps(t3, t1);											\
																			\
	t0 = _mm_unpacklo_ps( vq, vvz);											\
	t2 = _mm_unpacklo_ps(vix, viy);											\
	t1 = _mm_unpackhi_ps( vq, vvz);											\
	t3 = _mm_unpackhi_ps(vix, viy);											\
																			\
	vq   = _mm_movelh_ps(t0, t2);											\
	vvz  = _mm_movehl_ps(t2, t0);											\
	vix  = _mm_movelh_ps(t1, t3);											\
	viy  = _mm_movehl_ps(t3, t1);											\
	_mm_store_si128((__m128i*)ix, _mm_castps_si128(vix));					\
	_mm_store_si128((__m128i*)iy, _mm_castps_si128(viy));					\
 }


/***********************************************************************
t_vpbuf2D

Buffer to hold virtual particles for current deposition in 2D. 

************************************************************************/

typedef struct { 
  // 3 splits maximum
  DECLARE_ALIGNED_16( t_vp2D buf[ 3 * p_cache_size_2D ] );
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
STORE4VP3D

Stores 4 virtual particles for 3D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
************************************************************************/

#define STORE4VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { 	\
  register __m128 t0,t1,t2,t3,t4;												\
  																				\
  t0  = _mm_unpacklo_ps(vx0, vx1);												\
  t1  = _mm_unpacklo_ps(vy0, vy1);												\
  t2  = _mm_unpacklo_ps(vz0, vz1);												\
  t3  = _mm_unpacklo_ps( vq, _mm_castsi128_ps(vix));							\
  t4  = _mm_unpacklo_ps(_mm_castsi128_ps(viy), _mm_castsi128_ps(viz));			\
																				\
  _mm_store_ps(buf     , _mm_movelh_ps( t0, t1 ) );								\
  _mm_store_ps(buf +  4, _mm_movelh_ps( t2, t3 ) );								\
  _mm_store_ps(buf +  8, t4 );													\
																				\
  _mm_store_ps(buf + 12, _mm_movehl_ps( t1, t0 ) );								\
  _mm_store_ps(buf + 16, _mm_movehl_ps( t3, t2 ) );								\
  _mm_store_ps(buf + 20, _mm_movehl_ps( t4, t4 ) );								\
																				\
  t0  = _mm_unpackhi_ps(vx0, vx1);												\
  t1  = _mm_unpackhi_ps(vy0, vy1);												\
  t2  = _mm_unpackhi_ps(vz0, vz1);												\
  t3  = _mm_unpackhi_ps( vq, _mm_castsi128_ps(vix));							\
  t4  = _mm_unpackhi_ps(_mm_castsi128_ps(viy), _mm_castsi128_ps(viz));			\
																				\
  _mm_store_ps(buf + 24, _mm_movelh_ps( t0, t1 ) );								\
  _mm_store_ps(buf + 28, _mm_movelh_ps( t2, t3 ) );								\
  _mm_store_ps(buf + 32, t4 );													\
																				\
  _mm_store_ps(buf + 36, _mm_movehl_ps( t1, t0 ) );								\
  _mm_store_ps(buf + 40, _mm_movehl_ps( t3, t2 ) );								\
  _mm_store_ps(buf + 44, _mm_movehl_ps( t4, t4 ) );								\
 }

/***********************************************************************
LOAD4VP3D

Loads 4 virtual particles for 3D current deposition. 

************************************************************************/

#define LOAD4VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz ) { 		\
  register __m128 t0, t1, t2, t3, vix, viy, viz, c10, c11;						\
  vx0  = _mm_load_ps( buf      );												\
  vz0  = _mm_load_ps( buf + 4  );												\
  vix  = _mm_load_ps( buf + 8  );												\
  vx1  = _mm_load_ps( buf + 12 );												\
  vz1  = _mm_load_ps( buf + 16 );												\
  viy  = _mm_load_ps( buf + 20 );												\
  vy0  = _mm_load_ps( buf + 24 );												\
  vq   = _mm_load_ps( buf + 28 );												\
  viz  = _mm_load_ps( buf + 32 );												\
  vy1  = _mm_load_ps( buf + 36 );												\
  c10  = _mm_load_ps( buf + 40 );												\
  c11  = _mm_load_ps( buf + 44 );												\
																				\
  t0 = _mm_unpacklo_ps( vx0, vx1 );												\
  t1 = _mm_unpacklo_ps( vy0, vy1 );												\
  t2 = _mm_unpackhi_ps( vx0, vx1 );												\
  t3 = _mm_unpackhi_ps( vy0, vy1 );												\
																				\
  vx0 = _mm_movelh_ps( t0, t1 );												\
  vx1 = _mm_movehl_ps( t1, t0 );												\
  vy0 = _mm_movelh_ps( t2, t3 );												\
  vy1 = _mm_movehl_ps( t3, t2 );												\
																				\
  t0 = _mm_unpacklo_ps( vz0, vz1 );												\
  t1 = _mm_unpacklo_ps( vq, c10 );												\
  t2 = _mm_unpackhi_ps( vz0, vz1 );												\
  t3 = _mm_unpackhi_ps( vq, c10 );												\
  																				\
  vz0 = _mm_movelh_ps( t0, t1 );												\
  vz1 = _mm_movehl_ps( t1, t0 );												\
  vq  = _mm_movelh_ps( t2, t3 );												\
  t0 = _mm_unpacklo_ps( vix, viy );												\
  t1 = _mm_unpacklo_ps( viz, c11 );												\
  																				\
  vix = _mm_movehl_ps( t3, t2 );												\
  viy = _mm_movelh_ps( t0, t1 );												\
  viz = _mm_movehl_ps( t1, t0 );												\
																				\
  _mm_store_si128((__m128i*)ix, _mm_castps_si128(vix));							\
  _mm_store_si128((__m128i*)iy, _mm_castps_si128(viy));							\
  _mm_store_si128((__m128i*)iz, _mm_castps_si128(viz));							\
}

// Virtual particle buffer
typedef struct { 
  // 4 splits maximum
  DECLARE_ALIGNED_16( t_vp3D buf[ 4 * p_cache_size_3D ] );
  float *p;
  unsigned int np;
} t_vpbuf3D;

#else

#error Target system does not support SSE code.

#endif 


#endif
