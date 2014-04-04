#ifndef _OS_SPEC_PUSH_H_BGQ
#define _OS_SPEC_PUSH_H_BGQ



/* Size of particle buffer */

#define p_cache_size_2D 32
#define p_cache_size_3D 512

#include "vector-bgq.h"

/*****************************************************************************************
LOADFLD4

Loads 4 field values corresponding to 4 particle positions into a vector
variable. 

*****************************************************************************************/


#define LOADFLD4( fp, shift ) \
  (vector) { (fp[0])[shift], (fp[1])[shift], (fp[2])[shift], (fp[3])[shift] };


/*****************************************************************************************
t_vp2D

Virtual particle type for current deposition in 2D. This is organized in
memory as:

 | x0, x1 , y0, y1 | q, vz , (p0, ix), (p1, iy) |

The pusher will store one vp per simulation particle; the splitter will
then add new vp as necessary / correct the previous one.

This has the size of 2 vectors / 8 doubles.

*****************************************************************************************/

typedef struct virtualPart2D {
	double         x0,         x1,         y0,         y1;
	double          q,         vz; int p0, ix,     p1, iy;
} t_vp2D;

/*****************************************************************************************
STORE4VP2D

Stores 4 virtual particles for 2D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
*****************************************************************************************/

#define STORE4VP2D( buf, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) {	\
   																	\
   register vector _vix, _viy, t0, t1, t2, t3;						\
   																	\
   vector const unpacklo = vec_gpci( 00415 );						\
   vector const movelh   = vec_gpci( 00145 );						\
   vector const movehl   = vec_gpci( 06723 );						\
   vector const unpackhi = vec_gpci( 02637 );						\
   																	\
   _vix = vec_ldiaa( 0, vix );										\
   _viy = vec_ldiaa( 0, viy );										\
   																	\
   t0 = vec_perm( vx0, vx1, unpacklo );								\
   t1 = vec_perm( vy0, vy1, unpacklo );								\
   t2 = vec_perm( vq,  vvz, unpacklo );								\
   t3 = vec_perm( _vix, _viy, unpacklo );							\
   																	\
   vec_sta( vec_perm( t0, t1, movelh ), 0x00, buf );				\
   vec_sta( vec_perm( t2, t3, movelh ), 0x20, buf );				\
   vec_sta( vec_perm( t1, t0, movehl ), 0x40, buf );				\
   vec_sta( vec_perm( t3, t2, movehl ), 0x60, buf );				\
																	\
   t0 = vec_perm( vx0, vx1, unpackhi );								\
   t1 = vec_perm( vy0, vy1, unpackhi );								\
   t2 = vec_perm( vq,  vvz, unpackhi );								\
   t3 = vec_perm( _vix, _viy, unpackhi );							\
   																	\
   vec_sta( vec_perm( t0, t1, movelh ), 0x80, buf );				\
   vec_sta( vec_perm( t2, t3, movelh ), 0xA0, buf );				\
   vec_sta( vec_perm( t1, t0, movehl ), 0xC0, buf );				\
   vec_sta( vec_perm( t3, t2, movehl ), 0xE0, buf );				\
}

/*****************************************************************************************
LOAD4VP2D

Loads 4 virtual particles for 2D current deposition. 

*****************************************************************************************/

#define LOAD4VP2D( buf, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
   																\
   register vector _vix, _viy, t0, t1, t2, t3;					\
																\
   vector const unpacklo = vec_gpci( 00415 );					\
   vector const movelh   = vec_gpci( 00145 );					\
   vector const movehl   = vec_gpci( 06723 );					\
   vector const unpackhi = vec_gpci( 02637 );					\
																\
   vx0  = vec_lda( 0x00, buf );									\
   vq   = vec_lda( 0x20, buf );									\
   vx1  = vec_lda( 0x40, buf );									\
   vvz  = vec_lda( 0x60, buf );									\
   vy0  = vec_lda( 0x80, buf );									\
   _vix = vec_lda( 0xA0, buf );									\
   vy1  = vec_lda( 0xC0, buf );									\
   _viy = vec_lda( 0xE0, buf );									\
																\
   t0 = vec_perm(vx0, vx1, unpacklo);							\
   t2 = vec_perm(vy0, vy1, unpacklo);							\
   t1 = vec_perm(vx0, vx1, unpackhi);							\
   t3 = vec_perm(vy0, vy1, unpackhi);							\
																\
   vx0  = vec_perm(t0, t2, movelh);								\
   vx1  = vec_perm(t2, t0, movehl);								\
   vy0  = vec_perm(t1, t3, movelh);								\
   vy1  = vec_perm(t3, t1, movehl);								\
																\
   t0 = vec_perm( vq, vvz, unpacklo);							\
   t2 = vec_perm(_vix, _viy, unpacklo);							\
   t1 = vec_perm( vq, vvz, unpackhi);							\
   t3 = vec_perm(_vix, _viy, unpackhi);							\
																\
   vq   = vec_perm(t0, t2, movelh);								\
   vvz  = vec_perm(t2, t0, movehl);								\
   _vix  = vec_perm(t1, t3, movelh);							\
   _viy  = vec_perm(t3, t1, movehl);							\
   																\
   vec_sta( _vix, 0, vix );										\
   vec_sta( _viy, 0, viy );										\
 }


/*****************************************************************************************
t_vpbuf2D

Buffer to hold virtual particles for current deposition in 2D. 

*****************************************************************************************/

typedef struct { 
  // 3 splits maximum
  DECLARE_ALIGNED_32( t_vp2D buf[ 3 * p_cache_size_2D ] );
  double *p;
  unsigned int np;
} t_vpbuf2D;

/*****************************************************************************************
t_vp3D

Virtual particle type for current deposition. This corresponds to 3
vectors in the form:

 | x0, x1 , y0, y1 | z0, z1 , q, (ix,ip0) | (iy, ip1) , (iz, ip2), dp0, dp1 |

dp0 and ip0 are used for padding to improve memory alignment.
Structure size is 3 vectors / 12 doubles.

*****************************************************************************************/

typedef struct virtualPart3D {
	double      x0,       x1,        y0,        y1;
	double      z0,       z1,         q; int ip0, ix;
	int         ip1, iy,  ip2, iz; double dp0, dp1;  
} t_vp3D;

/*****************************************************************************************
STORE4VP3D

Stores 4 virtual particles for 3D current deposition. If the particle
did not cross any boundary, the splitter will not do any additional 
work.
*****************************************************************************************/

#define STORE4VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ){ \
																			\
   vector const unpacklo = vec_gpci( 00415 );								\
   vector const movelh   = vec_gpci( 00145 );								\
   vector const movehl   = vec_gpci( 06723 );								\
   vector const unpackhi = vec_gpci( 02637 );								\
																			\
   register vector _vix = vec_ldiaa( 0, vix );								\
   register vector _viy = vec_ldiaa( 0, viy );								\
   register vector _viz = vec_ldiaa( 0, viz );								\
   																			\
   register vector t0, t1, t2, t3, t4;										\
																			\
   t0  = vec_perm(vx0, vx1, unpacklo);										\
   t1  = vec_perm(vy0, vy1, unpacklo);										\
   t2  = vec_perm(vz0, vz1, unpacklo);										\
   t3  = vec_perm( vq, _vix, unpacklo);										\
   t4  = vec_perm(_viy, _viz, unpacklo);									\
																			\
   vec_sta( vec_perm( t0, t1, movelh ), 0x000, buf);						\
   vec_sta( vec_perm( t2, t3, movelh ), 0x020, buf);						\
   vec_sta(                         t4, 0x040, buf);						\
																			\
   vec_sta( vec_perm( t1, t0, movehl ), 0x060, buf);						\
   vec_sta( vec_perm( t3, t2, movehl ), 0x080, buf);						\
   vec_sta( vec_perm( t4, t4, movehl ), 0x0A0, buf);						\
																			\
   t0  = vec_perm(vx0, vx1, unpackhi);										\
   t1  = vec_perm(vy0, vy1, unpackhi);										\
   t2  = vec_perm(vz0, vz1, unpackhi);										\
   t3  = vec_perm( vq, _vix, unpackhi);										\
   t4  = vec_perm(_viy, _viz, unpackhi);									\
																			\
   vec_sta( vec_perm( t0, t1, movelh ), 0x0C0, buf );						\
   vec_sta( vec_perm( t2, t3, movelh ), 0x0E0, buf );						\
   vec_sta(                         t4, 0x100, buf );						\
																			\
   vec_sta( vec_perm( t1, t0, movehl ), 0x120, buf );						\
   vec_sta( vec_perm( t3, t2, movehl ), 0x140, buf );						\
   vec_sta( vec_perm( t4, t4, movehl ), 0x160, buf );						\
}

/*****************************************************************************************
LOAD4VP3D

Loads 4 virtual particles for 3D current deposition. 

*****************************************************************************************/

#define LOAD4VP3D( buf, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) {	\
																			\
   vector const unpacklo = vec_gpci( 00415 );								\
   vector const movelh   = vec_gpci( 00145 );								\
   vector const movehl   = vec_gpci( 06723 );								\
   vector const unpackhi = vec_gpci( 02637 );								\
																			\
   register vector t0, t1, t2, t3, _vix, _viy, _viz, c10, c11;				\
																			\
   vx0  = vec_lda( 0x000, buf );											\
   vz0  = vec_lda( 0x020, buf );											\
   _vix  = vec_lda( 0x040, buf );											\
   vx1  = vec_lda( 0x060, buf );											\
   vz1  = vec_lda( 0x080, buf );											\
   _viy  = vec_lda( 0x0A0, buf );											\
   vy0  = vec_lda( 0x0C0, buf );											\
   vq   = vec_lda( 0x0E0, buf );											\
   _viz  = vec_lda( 0x100, buf );											\
   vy1  = vec_lda( 0x120, buf );											\
   c10  = vec_lda( 0x140, buf );											\
   c11  = vec_lda( 0x160, buf );											\
																			\
   t0 = vec_perm( vx0, vx1 , unpacklo);										\
   t1 = vec_perm( vy0, vy1 , unpacklo);										\
   t2 = vec_perm( vx0, vx1 , unpackhi);										\
   t3 = vec_perm( vy0, vy1 , unpackhi);										\
																			\
   vx0 = vec_perm( t0, t1 , movelh);										\
   vx1 = vec_perm( t1, t0 , movehl);										\
   vy0 = vec_perm( t2, t3 , movelh);										\
   vy1 = vec_perm( t3, t2 , movehl);										\
																			\
   t0 = vec_perm( vz0, vz1 , unpacklo);										\
   t1 = vec_perm( vq, c10 , unpacklo);										\
   t2 = vec_perm( vz0, vz1 , unpackhi);										\
   t3 = vec_perm( vq, c10 , unpackhi);										\
																			\
   vz0 = vec_perm( t0, t1 , movelh);										\
   vz1 = vec_perm( t1, t0 , movehl);										\
   vq  = vec_perm( t2, t3 , movelh);										\
   t0 = vec_perm( _vix, _viy , unpacklo);									\
   t1 = vec_perm( _viz, c11 , unpacklo);									\
																			\
   _vix = vec_perm( t3, t2 , movehl);										\
   _viy = vec_perm( t0, t1 , movelh);										\
   _viz = vec_perm( t1, t0 , movehl);										\
																			\
   vec_sta( _vix, 0, vix );													\
   vec_sta( _viy, 0, viy );													\
   vec_sta( _viz, 0, viz );													\
																			\
}

// Virtual particle buffer
typedef struct { 
  // 4 splits maximum
  DECLARE_ALIGNED_32( t_vp3D buf[ 4 * p_cache_size_3D ] );
  double *p;
  unsigned int np;
} t_vpbuf3D;


/*****************************************************************************************
VRGAMMA

Get 1 / Lorentz gamma. 
*****************************************************************************************/

#define VRGAMMA( vrg, vu1, vu2, vu3 ) {                     \
   const vector one = vec_splats( 1.0 );                    \
   (vrg) = vec_rsqrt( vec_madd( (vu3), (vu3),               \
                      vec_madd( (vu2), (vu2),               \
                      vec_madd( (vu1), (vu1), one ) ) ) );  \
}


/*****************************************************************************************
vget_hp_near

Gets the nearest half points of the 4 particles loaded as a vector for
interpolating scattered grids. 
This routine is for positions defined as distance to the nearest grid
points.
*****************************************************************************************/

#define vget_hp_near(dx,ix,dxh,ixh,delta ) 					\
{ vector s;													\
															\
  vector const oneHalf = vec_splats( 0.5 );					\
  															\
  /* s = (dx < 0 ? -0.5 : 0.5) */							\
  s   = vec_sel( vec_splats( -0.5 ), oneHalf, dx );			\
  dxh = vec_sub( dx, s );									\
  															\
  s   = vec_sub( s, oneHalf );								\
  s   = vec_ctiw( s );										\
  vec_sta( s, 0, ixh );										\
  															\
  for( int i = 0; i < 4; i++ ) 								\
    ixh[i] = ix[i] - (ixh[i] & delta);						\
															\
}

inline vector vntrim( const vector vx );


#endif
