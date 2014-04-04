/*****************************************************************************************

Current deposition, SSE optimized version

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-sse.h"

#include <stdio.h>

#include "fortran.h"

#include "vector-sse.h"
#include "splines-sse.h"


/***********************************************************************
vwl_s1

Compute the longitudinal current weight for 1st order current deposition
***********************************************************************/

inline void vwl_s1( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] ) {

  vwl[0] = _mm_mul_ps( vqn, _mm_sub_ps( vx1, vx0 ) ); 
}


/***********************************************************************
vwl_s2

Compute the longitudinal current weight for 2nd order current deposition
***********************************************************************/

inline void vwl_s2( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] ) {

  __m128 vtx = _mm_add_ps( vx1, vx0 );
  __m128 vdx = _mm_sub_ps( vx1, vx0 );
  
  __m128 t0  = _mm_mul_ps( vdx, vtx );
 
  vwl[0] = _mm_mul_ps( vqn, _mm_sub_ps( vdx, t0 ) );
  vwl[1] = _mm_mul_ps( vqn, _mm_add_ps( vdx, t0 ) );  

}



/***********************************************************************
vwl_s3

Compute the longitudinal current weight for 3rd order current deposition
***********************************************************************/

inline void vwl_s3( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] ) {
   
   __m128 const c2 = _mm_set_ps1( 2.0f );
   
   __m128 vdx  = _mm_sub_ps( vx1, vx0 );
   __m128 vdx3 = _mm_mul_ps( _mm_set_ps1( 4.0f/3.0f ),
                             _mm_sub_ps( _mm_mul_ps( vx1, _mm_mul_ps( vx1, vx1 ) ), 
                                         _mm_mul_ps( vx0, _mm_mul_ps( vx0, vx0 ) ) ) );

   __m128 vtx  = _mm_mul_ps( _mm_mul_ps( c2, _mm_add_ps( vx1, vx0 ) ), vdx );

   __m128 t0 = _mm_add_ps( vdx, vdx3 );
   
   vwl[0] = _mm_mul_ps( vqn, _mm_sub_ps( t0, vtx ) );
   
   vwl[1] = _mm_mul_ps( vqn, _mm_mul_ps( c2, 
                                         _mm_sub_ps( _mm_mul_ps( _mm_set_ps1( 3.0f ), vdx ), 
                                                     vdx3 ) ) ) ;
   vwl[2] = _mm_mul_ps( vqn, _mm_add_ps( t0, vtx ) );

}

/***********************************************************************
vwl_s4

Compute the longitudinal current weight for 4th order current deposition
***********************************************************************/
inline void vwl_s4( __m128 const vqn, __m128 const vx0, __m128 const vx1, __m128 vwl[] ) {

  __m128 vx02 = _mm_mul_ps( vx0, vx0 );
  __m128 vx12 = _mm_mul_ps( vx1, vx1 );
  
  __m128 vdx = _mm_sub_ps( vx1, vx0 );
  __m128 vtx = _mm_mul_ps( _mm_set_ps1( 2.0f ), _mm_add_ps( vx1, vx0 ) );
  
  __m128 vdx2 = _mm_sub_ps( vx12, vx02 );
  
  
  __m128 t0 = _mm_add_ps( _mm_set_ps1( 1.0f ), 
                          _mm_mul_ps( _mm_set_ps1(2.0f), _mm_add_ps( vx12, vx02 ) ) );

  __m128 t1 = _mm_add_ps( _mm_mul_ps( _mm_set_ps1(23.0f), vdx ),
                          _mm_mul_ps( _mm_set_ps1(4.0f), 
                                      _mm_sub_ps( _mm_mul_ps( vx02, vx0 ),
                                                  _mm_mul_ps( vx12, vx1 ) ) ) ); 

  __m128 t2 = _mm_add_ps( _mm_mul_ps( _mm_set_ps1(15.0f), vdx2 ),
                          _mm_mul_ps( _mm_set_ps1(6.0f), 
                                      _mm_sub_ps( _mm_mul_ps( vx02, vx02 ),
                                                  _mm_mul_ps( vx12, vx12 ) ) ) ); 
  
  vwl[0] = _mm_mul_ps( vqn,
                       _mm_mul_ps( _mm_sub_ps( vdx, vdx2 ), _mm_sub_ps( t0, vtx ) ) );
  vwl[1] = _mm_mul_ps( vqn, _mm_sub_ps( t1, t2 ) ); 
  vwl[2] = _mm_mul_ps( vqn, _mm_add_ps( t1, t2 ) );
  vwl[3] = _mm_mul_ps( vqn,
                       _mm_mul_ps( _mm_add_ps( vdx, vdx2 ), _mm_add_ps( t0, vtx ) ) );
      
}



/****************************************************************************************

  Generate specific current deposition functions for 2D and 3D, 1st to 4th order

****************************************************************************************/

#define __TEMPLATE__

// These macros will append _s1, _s2, etc to function names.
#define ONAME(f, o) OJOIN(f, o)
#define OJOIN(f, o) f ## _s ## o

/********************************** Linear interpolation ********************************/

// Interpolation order
#define ORDER 1
// Number of interpolation points ( ORDER + 1 )
#define NP 2
// Current grid offset (1 + ORDER/2)
#define OFFSET 1

#include __FILE__

/******************************** Quadratic interpolation *******************************/

#define ORDER 2
#define NP 3
#define OFFSET 2

#include __FILE__

/********************************** Cubic interpolation *********************************/

#define ORDER 3
#define NP 4
#define OFFSET 2

#include __FILE__

/********************************* Quartic interpolation ********************************/

#define ORDER 4
#define NP 5
#define OFFSET 3

#include __FILE__

#else

/****************************************************************************************

  Template function definitions for 2D and 3D current deposition

****************************************************************************************/

/*

Note: In 2D there was a division by 2 removed from the calculation of the perpendicular
      weights, so jnorm must be dx / dt / (2*NORM). In 3D it was a division by 3, so 
      jnorm must be dx / dt / (3*NORM). To avoid repeating this every time the 
      DEP_CURRENT functions are called, this is done in the ADVANCE_DEPOSIT functions.

*/

/***************   Generate Function names based on interpolation order  ****************/

#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

#define SPLINE  ONAME( vspline, ORDER )
#define WL      ONAME( vwl, ORDER )


inline void DEP_CURRENT_2D
(float * const current, int const * const size, int const * const offset, 
 float * const norm, t_vpbuf2D * const part)
{
        
  typedef struct { float j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
    
  fvec j3[NP][NP];
  
  fvec vwp1[NP], vwp2[NP]; 
  fvec vwl1[ORDER], vwl2[ORDER];
  
  __m128 vx0, vx1, vy0, vy1, vq, vvz;
  
  DECLARE_ALIGNED_16( int ix[4] );
  DECLARE_ALIGNED_16( int iy[4] );

		  
  __m128 vqnx, vqny, vqvz;
  __m128 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 

  int np;
  int i, k, k1, k2;

  __m128 const oneThird = _mm_set_ps1( 1.0f/3.0f );

  int const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                                 ( offset[1] - OFFSET ) * 3 * Dy );

  __m128 const vnorm1 = _mm_set_ps1(norm[0]);
  __m128 const vnorm2 = _mm_set_ps1(norm[1]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 4 add dummy particles to the end

  if ( np % 4 != 0 ) {
    for( i = 0; i < 4 - np%4; i ++ ) {
       t_vp2D * vpt = ((t_vp2D *) part -> p) + i;
       
       vpt -> x0 = vpt -> x1 = 0.;
       vpt -> y0 = vpt -> y1 = 0.;
       vpt -> q  = vpt -> vz = 0.;
       vpt -> ix = vpt -> iy  = 1;
    }
  }

  float* vp = (float *) part -> buf;

  for( i = 0; i < np; i+=4, vp+=32 ) {

    // load 4 particles
    LOAD4VP2D( vp, vx0, vx1, vy0, vy1, vq, vvz, ix, iy );
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = _mm_mul_ps( vq, vnorm1);
    vqny = _mm_mul_ps( vq, vnorm2);
    vqvz = _mm_mul_ps( vq, vvz);
    vqvz = _mm_mul_ps( vqvz, oneThird );

 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m128 *) vwl1 );
    WL( vqny, vy0, vy1, (__m128 *) vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k].v4 = _mm_add_ps(vs0y[k], vs1y[k]);
      vwp2[k].v4 = _mm_add_ps(vs0x[k], vs1x[k]);
    }
     
    // get j3 current
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
		__m128 const oneHalf  = _mm_set_ps1( 0.5f );
		__m128 tmp1, tmp2;

		tmp1 = _mm_add_ps( _mm_mul_ps( vs0x[k1], vs0y[k2] ), 
		                   _mm_mul_ps( vs1x[k1], vs1y[k2] ) );
	   
		tmp2 = _mm_add_ps( _mm_mul_ps( vs0x[k1], vs1y[k2] ), 
		       _mm_mul_ps( vs1x[k1], vs0y[k2] ) );
		
		j3[k1][k2].v4 = _mm_mul_ps( vqvz, _mm_add_ps( tmp1, _mm_mul_ps( tmp2, oneHalf ) ) );
      }
    }

    // New version loop by particle on the outside loop
    for ( k = 0; k < 4; k ++ ) {

      pjpart = pj + ix[k] + iy[k] * Dy;

	  // accumulate j1
	  for( k2 = 0; k2 < NP; k2++ ) {
		p0 = pjpart + k2*Dy;
		for ( k1 = 0; k1 < ORDER; k1++ ) {
		  p0[k1].j1 += vwl1[k1].v[k] * vwp1[k2].v[k];
		}
	  }
  
	  // accumulate j2 - making k2 the outside loop gives marginal perf. gain
	  for( k2 = 0; k2 < ORDER; k2++ ) {
		p0 = pjpart + k2*Dy;
		for ( k1 = 0; k1 < NP; k1++ ) {
		   p0[k1].j2 += vwl2[k2].v[k] * vwp2[k1].v[k];
		}
	  }
  
	  // accumulate j3
	  for( k2 = 0; k2 < NP; k2++ ) {
		p0 = pjpart + k2*Dy;
		for ( k1 = 0; k1 < NP; k1++ ) {
		   p0[k1].j3 += (j3[k1][k2]).v[k];
		}
	  }
	  
    }
	  
  }
}


inline void DEP_CURRENT_3D
(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_vpbuf3D * const part)
{
        
  typedef struct { float j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
  float *vp;
    
  fvec vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  fvec vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];
  
  __m128 vx0, vx1, vy0, vy1,  vz0, vz1, vq; 
  DECLARE_ALIGNED_16( int ix[4] );
  DECLARE_ALIGNED_16( int iy[4] );
  DECLARE_ALIGNED_16( int iz[4] );
  
  __m128 vqnx, vqny, vqnz;
  __m128 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 
  
  
  int np;
  int i, k, k1, k2, k3;
  
  int const Dy = size[0];    
  int const Dz = Dy * size[1]; 
  
  t_current* const pj = (t_current *) ( current + ( offset[0] - OFFSET ) * 3 + 
                                                  ( offset[1] - OFFSET ) * 3 * Dy + 
                                                  ( offset[2] - OFFSET ) * 3 * Dz );

  __m128 const vnorm1 = _mm_set_ps1(norm[0]);
  __m128 const vnorm2 = _mm_set_ps1(norm[1]);
  __m128 const vnorm3 = _mm_set_ps1(norm[2]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np % 4 != 0 ) {
    for( i = 0; i < 4 - np%4; i ++ ) {
       t_vp3D * vpt = (t_vp3D *) part -> p + i;
       
       vpt -> x0 = vpt -> x1 = 0.;
       vpt -> y0 = vpt -> y1 = 0.;
       vpt -> z0 = vpt -> z1 = 0.;
       vpt -> q = 0.;
       vpt -> ix = vpt -> iy = vpt -> iz = 1;
    }
  }

  vp = (float *) part -> buf;

  for( i = 0; i < np; i+=4, vp+=48 ) {
      
	// load 4 particles
    LOAD4VP3D( vp, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz )

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = _mm_mul_ps( vq, vnorm1);
    vqny = _mm_mul_ps( vq, vnorm2);
    vqnz = _mm_mul_ps( vq, vnorm3);
    
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m128 *) vwl1 );
    WL( vqny, vy0, vy1, (__m128 *) vwl2 );
    WL( vqnz, vz0, vz1, (__m128 *) vwl3 );

    
    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
      for( k1 = 0; k1 < NP; k1++ ) {
         
         __m128 tmp1, tmp2;        
         const __m128 oneHalf = _mm_set_ps1(0.5f);

         // wp1[k1][k2] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] + 
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )
         
         tmp1 =  _mm_add_ps( _mm_mul_ps( vs0y[k1], vs0z[k2]), _mm_mul_ps( vs1y[k1], vs1z[k2]));
         tmp2 =  _mm_add_ps( _mm_mul_ps( vs0y[k1], vs1z[k2]), _mm_mul_ps( vs1y[k1], vs0z[k2]));
         vwp1[k2][k1].v4 = _mm_add_ps( tmp1, _mm_mul_ps( tmp2, oneHalf));

         tmp1 =  _mm_add_ps( _mm_mul_ps( vs0x[k1], vs0z[k2]), _mm_mul_ps( vs1x[k1], vs1z[k2]));
         tmp2 =  _mm_add_ps( _mm_mul_ps( vs0x[k1], vs1z[k2]), _mm_mul_ps( vs1x[k1], vs0z[k2]));
         vwp2[k2][k1].v4 = _mm_add_ps( tmp1, _mm_mul_ps( tmp2, oneHalf));

         tmp1 =  _mm_add_ps( _mm_mul_ps( vs0x[k1], vs0y[k2]), _mm_mul_ps( vs1x[k1], vs1y[k2]));
         tmp2 =  _mm_add_ps( _mm_mul_ps( vs0x[k1], vs1y[k2]), _mm_mul_ps( vs1x[k1], vs0y[k2]));
         vwp3[k2][k1].v4 = _mm_add_ps( tmp1, _mm_mul_ps( tmp2, oneHalf));

      }
    }

    // looping by particle on the outside loop yields the best performance
    for ( k = 0; k < 4; k ++ ) {

       pjpart = pj + ix[k] + iy[k] * Dy + iz[k] * Dz;
          
	   // accumulate j1
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < ORDER; k1++ ) {
			  p0[k1].j1 += vwl1[k1].v[k] * vwp1[k3][k2].v[k];
		   }
		 }
	   }
   
	   // accumulate j2
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < ORDER; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < NP; k1++ ) {
			 p0[k1].j2 += vwl2[k2].v[k] * vwp2[k3][k1].v[k];
		   }
		 }
	   }
   
	   
	   // accumulate j3
	   for( k3 = 0; k3 < ORDER; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < NP; k1++ ) {
			 p0[k1].j3  += vwl3[k3].v[k] * vwp3[k2][k1].v[k];
		   }
		 }
	   }

    }

  }
}

#undef DEP_CURRENT_2D
#undef DEP_CURRENT_3D
#undef SPLINE
#undef WL

#undef OFFSET
#undef ORDER
#undef NP

#endif
