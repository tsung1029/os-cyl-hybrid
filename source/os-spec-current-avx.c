/*****************************************************************************************

Current deposition, AVX optimized version

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-avx.h"

#include <stdio.h>

#include "fortran.h"

#include "vector-avx.h"
#include "splines-avx.h"


/***********************************************************************
vwl_s1

Compute the longitudinal current weight for 1st order current deposition
***********************************************************************/

inline void vwl_s1( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] ) {

  vwl[0] = _mm256_mul_ps( vqn, _mm256_sub_ps( vx1, vx0 ) ); 
}


/***********************************************************************
vwl_s2

Compute the longitudinal current weight for 2nd order current deposition
***********************************************************************/

inline void vwl_s2( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] ) {

  __m256 vtx = _mm256_add_ps( vx1, vx0 );
  __m256 vdx = _mm256_sub_ps( vx1, vx0 );
  
  __m256 t0  = _mm256_mul_ps( vdx, vtx );
 
  vwl[0] = _mm256_mul_ps( vqn, _mm256_sub_ps( vdx, t0 ) );
  vwl[1] = _mm256_mul_ps( vqn, _mm256_add_ps( vdx, t0 ) );   

}



/***********************************************************************
vwl_s3

Compute the longitudinal current weight for 3rd order current deposition
***********************************************************************/

inline void vwl_s3( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] ) {
   
   __m256 const c2 = _mm256_set1_ps( 2.0f );
   
   __m256 vdx  = _mm256_sub_ps( vx1, vx0 );
   __m256 vdx3 = _mm256_mul_ps( _mm256_set1_ps( 4.0f/3.0f ),
                             _mm256_sub_ps( _mm256_mul_ps( vx1, _mm256_mul_ps( vx1, vx1 ) ), 
                                         _mm256_mul_ps( vx0, _mm256_mul_ps( vx0, vx0 ) ) ) );

   __m256 vtx  = _mm256_mul_ps( _mm256_mul_ps( c2, _mm256_add_ps( vx1, vx0 ) ), vdx );

   __m256 t0 = _mm256_add_ps( vdx, vdx3 );
   
   vwl[0] = _mm256_mul_ps( vqn, _mm256_sub_ps( t0, vtx ) );
   
   vwl[1] = _mm256_mul_ps( vqn, _mm256_mul_ps( c2, 
                                         _mm256_sub_ps( _mm256_mul_ps( _mm256_set1_ps( 3.0f ), vdx ), 
                                                     vdx3 ) ) ) ;
   vwl[2] = _mm256_mul_ps( vqn, _mm256_add_ps( t0, vtx ) );

}

/***********************************************************************
vwl_s4

Compute the longitudinal current weight for 4th order current deposition
***********************************************************************/
inline void vwl_s4( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] ) {

  __m256 vx02 = _mm256_mul_ps( vx0, vx0 );
  __m256 vx12 = _mm256_mul_ps( vx1, vx1 );
  
  __m256 vdx = _mm256_sub_ps( vx1, vx0 );
  __m256 vtx = _mm256_mul_ps( _mm256_set1_ps( 2.0f ), _mm256_add_ps( vx1, vx0 ) );
  
  __m256 vdx2 = _mm256_sub_ps( vx12, vx02 );
  
  
  __m256 t0 = _mm256_add_ps( _mm256_set1_ps( 1.0f ), 
                          _mm256_mul_ps( _mm256_set1_ps(2.0f), _mm256_add_ps( vx12, vx02 ) ) );

  __m256 t1 = _mm256_add_ps( _mm256_mul_ps( _mm256_set1_ps(23.0f), vdx ),
                          _mm256_mul_ps( _mm256_set1_ps(4.0f), 
                                      _mm256_sub_ps( _mm256_mul_ps( vx02, vx0 ),
                                                  _mm256_mul_ps( vx12, vx1 ) ) ) ); 

  __m256 t2 = _mm256_add_ps( _mm256_mul_ps( _mm256_set1_ps(15.0f), vdx2 ),
                          _mm256_mul_ps( _mm256_set1_ps(6.0f), 
                                      _mm256_sub_ps( _mm256_mul_ps( vx02, vx02 ),
                                                  _mm256_mul_ps( vx12, vx12 ) ) ) ); 
  
  vwl[0] = _mm256_mul_ps( vqn,
                       _mm256_mul_ps( _mm256_sub_ps( vdx, vdx2 ), _mm256_sub_ps( t0, vtx ) ) );
  vwl[1] = _mm256_mul_ps( vqn, _mm256_sub_ps( t1, t2 ) ); 
  vwl[2] = _mm256_mul_ps( vqn, _mm256_add_ps( t1, t2 ) );
  vwl[3] = _mm256_mul_ps( vqn,
                       _mm256_mul_ps( _mm256_add_ps( vdx, vdx2 ), _mm256_add_ps( t0, vtx ) ) );
      
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
  
  __m256 vx0, vx1, vy0, vy1, vq, vvz;
  
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
		  
  __m256 vqnx, vqny, vqvz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 
 
  int np;
  int i, k, k1, k2;

  __m256 const oneThird = _mm256_set1_ps( 1.0f/3.0f );

  int const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                                 ( offset[1] - OFFSET ) * 3 * Dy );

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH != 0 ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       t_vp2D * vpt = ((t_vp2D *) part -> p) + i;
       
       vpt -> x0 = vpt -> x1 = 0.;
       vpt -> y0 = vpt -> y1 = 0.;
       vpt -> q  = vpt -> vz = 0.;
       vpt -> ix = vpt -> iy  = 1;
    }
  }

  float* vp = (float *) part -> buf;
  
  // Each vp is 8 floats
  for( i = 0; i < np; i+=VEC_WIDTH, vp+= 8 * VEC_WIDTH ) {
    
    __m256i vix, viy;
    
    // load 8 particles
    LOAD8VP2D( vp, vx0, vx1, vy0, vy1, vq, vvz, vix, viy );
    
    // Store cell index
    viadd_store( idx, vix, vimul_s( viy, Dy ) ) ;
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = _mm256_mul_ps( vq, vnorm1);
    vqny = _mm256_mul_ps( vq, vnorm2);
    vqvz = _mm256_mul_ps( vq, vvz);
    vqvz = _mm256_mul_ps( vqvz, oneThird );
 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 );
    WL( vqny, vy0, vy1, (__m256 *) vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k].v8 = _mm256_add_ps(vs0y[k], vs1y[k]);
      vwp2[k].v8 = _mm256_add_ps(vs0x[k], vs1x[k]);
    }
     
    // get j3 current
    for( k1 = 0; k1 < NP; k1++ ) {
      for ( k2 = 0; k2 < NP; k2++ ) {
  		
  		__m256 const oneHalf  = _mm256_set1_ps( 0.5f );
  		__m256  tmp1, tmp2;
  		
		tmp1 = _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs0y[k2] ), 
		                      _mm256_mul_ps( vs1x[k1], vs1y[k2] ) );
	    
		tmp2 = _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs1y[k2] ), 
		                      _mm256_mul_ps( vs1x[k1], vs0y[k2] ) );

		j3[k1][k2].v8 = _mm256_mul_ps( vqvz, _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf ) ) );
      }
    }

    // Loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

      // This is not vectorially because it's a 64bit op
      pjpart = pj + idx[k];

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
      
  __m256 vx0, vx1, vy0, vy1,  vz0, vz1, vq; 
  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
  
  __m256 vqnx, vqny, vqnz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 
  
  fvec vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  fvec vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];
  
  unsigned int np;
  unsigned int i, k, k1, k2, k3;
     
  int const Dy = size[0];    
  int const Dz = Dy * size[1]; 
    
  t_current* const pj = (t_current *) ( current + ( offset[0] - OFFSET  ) * 3 + 
                                                  ( offset[1] - OFFSET  ) * 3 * Dy + 
                                                  ( offset[2] - OFFSET  ) * 3 * Dz );

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);
  __m256 const vnorm3 = _mm256_set1_ps(norm[2]);
  
  np = part -> np;
    
  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       t_vp3D* vpt = (t_vp3D *) part -> p + i;
       
       vpt -> x0 = vpt -> x1 = 0.;
       vpt -> y0 = vpt -> y1 = 0.;
       vpt -> z0 = vpt -> z1 = 0.;
       vpt -> q = 0.;
       vpt -> ix = vpt -> iy = vpt -> iz = 1;
    }
  }

  float* vp = (float *) part -> buf;

  for( i = 0; i < np; i+=VEC_WIDTH, vp+= 12 * VEC_WIDTH ) {

    __m256i vix, viy, viz;
      
	// load 8 particles
    LOAD8VP3D( vp, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz )
    
    // Do idx[k] = ix[k] + iy[k]*Dy + iz[k]*Dz vectorially 
    {
      __m128i vDy, vDz, a, b;
      vDy = _mm_set1_epi32( Dy );
      vDz = _mm_set1_epi32( Dz );
      
      a = _mm_mullo_epi32( _mm256_castsi256_si128(viz), vDz );
      b = _mm_mullo_epi32( _mm256_extractf128_si256(viz,1), vDz );

      a = _mm_add_epi32( a, _mm_mullo_epi32( _mm256_castsi256_si128(viy), vDy ));
      b = _mm_add_epi32( b, _mm_mullo_epi32( _mm256_extractf128_si256(viy,1), vDy ));
      
      a = _mm_add_epi32( a, _mm256_castsi256_si128(vix) );
      b = _mm_add_epi32( b, _mm256_extractf128_si256(vix,1) );
      
      _mm_store_si128((__m128i *)  &idx[0], a );
      _mm_store_si128((__m128i *)  &idx[4], b );
    }
    
    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = _mm256_mul_ps( vq, vnorm1);
    vqny = _mm256_mul_ps( vq, vnorm2);
    vqnz = _mm256_mul_ps( vq, vnorm3);
    
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 );
    WL( vqny, vy0, vy1, (__m256 *) vwl2 );
    WL( vqnz, vz0, vz1, (__m256 *) vwl3 );
    
    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
      for( k1 = 0; k1 < NP; k1++ ) {
         
         // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] + 
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )
         
         const __m256 oneHalf = _mm256_set1_ps(0.5f);
         __m256 tmp1, tmp2;        
        
         tmp1 =  _mm256_add_ps( _mm256_mul_ps( vs0y[k1], vs0z[k2]), _mm256_mul_ps( vs1y[k1], vs1z[k2]));
         tmp2 =  _mm256_add_ps( _mm256_mul_ps( vs0y[k1], vs1z[k2]), _mm256_mul_ps( vs1y[k1], vs0z[k2]));
         vwp1[k2][k1].v8 = _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf));

         tmp1 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs0z[k2]), _mm256_mul_ps( vs1x[k1], vs1z[k2]));
         tmp2 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs1z[k2]), _mm256_mul_ps( vs1x[k1], vs0z[k2]));
         vwp2[k2][k1].v8 = _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf));

         tmp1 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs0y[k2]), _mm256_mul_ps( vs1x[k1], vs1y[k2]));
         tmp2 =  _mm256_add_ps( _mm256_mul_ps( vs0x[k1], vs1y[k2]), _mm256_mul_ps( vs1x[k1], vs0y[k2]));
         vwp3[k2][k1].v8 = _mm256_add_ps( tmp1, _mm256_mul_ps( tmp2, oneHalf));

      }
    }

    // looping by particle on the outside loop yields the best performance
     
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

       pjpart = pj + idx[k];
          
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
