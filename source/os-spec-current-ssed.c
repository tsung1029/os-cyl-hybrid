/*****************************************************************************************

Current deposition, SSE optimized version

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-ssed.h"

#include <stdio.h>

#include "fortran.h"

#include "vector-sse.h"
#include "splines-sse.h"


/***********************************************************************
vwl_s1

Compute the longitudinal current weight for 1st order current deposition
***********************************************************************/

inline void vwl_s1( __m128d const vqn, __m128d const vx0, __m128d const vx1, __m128d vwl[] ) {

  vwl[0] = _mm_mul_pd( vqn, _mm_sub_pd( vx1, vx0 ) ); 
}


/***********************************************************************
vwl_s2

Compute the longitudinal current weight for 2nd order current deposition
***********************************************************************/

inline void vwl_s2( __m128d const vqn, __m128d const vx0, __m128d const vx1, __m128d vwl[] ) {

  __m128d vtx = _mm_add_pd( vx1, vx0 );
  __m128d vdx = _mm_sub_pd( vx1, vx0 );
  
  __m128d t0  = _mm_mul_pd( vdx, vtx );
 
  vwl[0] = _mm_mul_pd( vqn, _mm_sub_pd( vdx, t0 ) );
  vwl[1] = _mm_mul_pd( vqn, _mm_add_pd( vdx, t0 ) );  

}



/***********************************************************************
vwl_s3

Compute the longitudinal current weight for 3rd order current deposition
***********************************************************************/

inline void vwl_s3( __m128d const vqn, __m128d const vx0, __m128d const vx1, __m128d vwl[] ) {
   
   __m128d const c2 = _mm_set1_pd( 2.0 );
   
   __m128d vdx  = _mm_sub_pd( vx1, vx0 );
   __m128d vdx3 = _mm_mul_pd( _mm_set1_pd( 4.0/3.0 ),
                             _mm_sub_pd( _mm_mul_pd( vx1, _mm_mul_pd( vx1, vx1 ) ), 
                                         _mm_mul_pd( vx0, _mm_mul_pd( vx0, vx0 ) ) ) );

   __m128d vtx  = _mm_mul_pd( _mm_mul_pd( c2, _mm_add_pd( vx1, vx0 ) ), vdx );

   __m128d t0 = _mm_add_pd( vdx, vdx3 );
   
   vwl[0] = _mm_mul_pd( vqn, _mm_sub_pd( t0, vtx ) );
   
   vwl[1] = _mm_mul_pd( vqn, _mm_mul_pd( c2, 
                                         _mm_sub_pd( _mm_mul_pd( _mm_set1_pd( 3.0 ), vdx ), 
                                                     vdx3 ) ) ) ;
   vwl[2] = _mm_mul_pd( vqn, _mm_add_pd( t0, vtx ) );

}

/***********************************************************************
vwl_s4

Compute the longitudinal current weight for 4th order current deposition
***********************************************************************/
inline void vwl_s4( __m128d const vqn, __m128d const vx0, __m128d const vx1, __m128d vwl[] ) {

  __m128d vx02 = _mm_mul_pd( vx0, vx0 );
  __m128d vx12 = _mm_mul_pd( vx1, vx1 );
  
  __m128d vdx = _mm_sub_pd( vx1, vx0 );
  __m128d vtx = _mm_mul_pd( _mm_set1_pd( 2.0 ), _mm_add_pd( vx1, vx0 ) );
  
  __m128d vdx2 = _mm_sub_pd( vx12, vx02 );
  
  
  __m128d t0 = _mm_add_pd( _mm_set1_pd( 1.0 ), 
                          _mm_mul_pd( _mm_set1_pd(2.0), _mm_add_pd( vx12, vx02 ) ) );

  __m128d t1 = _mm_add_pd( _mm_mul_pd( _mm_set1_pd(23.0), vdx ),
                          _mm_mul_pd( _mm_set1_pd(4.0), 
                                      _mm_sub_pd( _mm_mul_pd( vx02, vx0 ),
                                                  _mm_mul_pd( vx12, vx1 ) ) ) ); 

  __m128d t2 = _mm_add_pd( _mm_mul_pd( _mm_set1_pd(15.0), vdx2 ),
                          _mm_mul_pd( _mm_set1_pd(6.0), 
                                      _mm_sub_pd( _mm_mul_pd( vx02, vx02 ),
                                                  _mm_mul_pd( vx12, vx12 ) ) ) ); 
  
  vwl[0] = _mm_mul_pd( vqn,
                       _mm_mul_pd( _mm_sub_pd( vdx, vdx2 ), _mm_sub_pd( t0, vtx ) ) );
  vwl[1] = _mm_mul_pd( vqn, _mm_sub_pd( t1, t2 ) ); 
  vwl[2] = _mm_mul_pd( vqn, _mm_add_pd( t1, t2 ) );
  vwl[3] = _mm_mul_pd( vqn,
                       _mm_mul_pd( _mm_add_pd( vdx, vdx2 ), _mm_add_pd( t0, vtx ) ) );
      
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

#define SPLINE  ONAME( vsplined, ORDER )
#define WL      ONAME( vwl, ORDER )


inline void DEP_CURRENT_2D
(double * const current, int const * const size, int const * const offset, 
 double * const norm, t_vpbuf2D * const part)
{
        
  typedef struct Current { double j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
  
  dvec j3[NP][NP];
  

  dvec vwp1[NP], vwp2[NP]; 
  dvec vwl1[ORDER], vwl2[ORDER];
  
  __m128d vx0, vx1, vy0, vy1, vq, vvz;
  
  DECLARE_ALIGNED_16( int ix[4] );
  DECLARE_ALIGNED_16( int iy[4] );

		  
  __m128d vqnx, vqny, vqvz;
  __m128d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 

  int np;
  int i, k, k1, k2;

  __m128d const oneThird = _mm_set1_pd( 1.0/3.0 );

  const int Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                                 ( offset[1] - OFFSET ) * 3 * Dy );


  __m128d const vnorm1 = _mm_set1_pd(norm[0]);
  __m128d const vnorm2 = _mm_set1_pd(norm[1]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 2 add dummy particles to the end

  if ( np % 2 != 0 ) {
    for( i = 0; i < 2 - np%2; i ++ ) {
       t_vp2D * vpt = ((t_vp2D *) part -> p) + i;
       
       vpt -> x0 = vpt -> x1 = 0.;
       vpt -> y0 = vpt -> y1 = 0.;
       vpt -> q  = vpt -> vz = 0.;
       vpt -> ix = vpt -> iy  = 1;
    }
  }

  double* vp = (double *) part -> buf;

  for( i = 0; i < np; i+=2, vp+=16 ) {

    // load 2 particles
    LOAD2VP2D( vp, vx0, vx1, vy0, vy1, vq, vvz, ix, iy );
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = _mm_mul_pd( vq, vnorm1);
    vqny = _mm_mul_pd( vq, vnorm2);
    vqvz = _mm_mul_pd( vq, vvz);
    vqvz = _mm_mul_pd( vqvz, oneThird );

 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m128d *) vwl1 );
    WL( vqny, vy0, vy1, (__m128d *) vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k].v2 = _mm_add_pd(vs0y[k], vs1y[k]);
      vwp2[k].v2 = _mm_add_pd(vs0x[k], vs1x[k]);
    }
     
    // get j3 current
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
        __m128d const oneHalf  = _mm_set1_pd( 0.5 );
        __m128d tmp1, tmp2;
  		
		tmp1 = _mm_add_pd( _mm_mul_pd( vs0x[k1], vs0y[k2] ), _mm_mul_pd( vs1x[k1], vs1y[k2] ) );
		tmp2 = _mm_add_pd( _mm_mul_pd( vs0x[k1], vs1y[k2] ), _mm_mul_pd( vs1x[k1], vs0y[k2] ) );

		j3[k1][k2].v2 = _mm_mul_pd( vqvz, _mm_add_pd( tmp1, _mm_mul_pd( tmp2, oneHalf ) ) );
      }
    }

    // New version loop by particle on the outside loop
    for ( k = 0; k < 2; k ++ ) {

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
(double * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_vpbuf3D * const part)
{
        
  typedef struct Current { double j1, j2, j3; } t_current;
  t_current *pjpart, *p0; 
  
  dvec vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  dvec vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];
  
  __m128d vx0, vx1, vy0, vy1,  vz0, vz1, vq; 
  DECLARE_ALIGNED_16( int ix[4] );
  DECLARE_ALIGNED_16( int iy[4] );
  DECLARE_ALIGNED_16( int iz[4] );
  
  __m128d vqnx, vqny, vqnz;
  __m128d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 
  
  
  int np;
  int i, k, k1, k2, k3;
  
  int const Dy = size[0];    
  int const Dz = Dy * size[1]; 
  
  t_current* const pj = (t_current *) ( current + ( offset[0] - OFFSET ) * 3 + 
                                                  ( offset[1] - OFFSET ) * 3 * Dy + 
                                                  ( offset[2] - OFFSET ) * 3 * Dz );

  __m128d const vnorm1 = _mm_set1_pd(norm[0]);
  __m128d const vnorm2 = _mm_set1_pd(norm[1]);
  __m128d const vnorm3 = _mm_set1_pd(norm[2]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 2 add dummy particles to the end
  if ( np % 2 != 0 ) {
    for( i = 0; i < 2 - np%2; i ++ ) {
       t_vp3D * vpt = (t_vp3D *) part -> p + i;
       
       vpt -> x0 = vpt -> x1 = 0.;
       vpt -> y0 = vpt -> y1 = 0.;
       vpt -> z0 = vpt -> z1 = 0.;
       vpt -> q = 0.;
       vpt -> ix = vpt -> iy = vpt -> iz = 1;
    }
  }

  double* vp = (double *) part -> buf;

  for( i = 0; i < np; i += 2, vp += 2 * 10 ) {
      
	// load 2 particles
    LOAD2VP3D( vp, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz );

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = _mm_mul_pd( vq, vnorm1);
    vqny = _mm_mul_pd( vq, vnorm2);
    vqnz = _mm_mul_pd( vq, vnorm3);
    
    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m128d *) vwl1 );
    WL( vqny, vy0, vy1, (__m128d *) vwl2 );
    WL( vqnz, vz0, vz1, (__m128d *) vwl3 );

    
    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
      for( k1 = 0; k1 < NP; k1++ ) {
		 __m128d tmp1, tmp2;        
		 const __m128d oneHalf = _mm_set1_pd(0.5);
         
         // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] + 
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )
         
         tmp1 =  _mm_add_pd( _mm_mul_pd( vs0y[k1], vs0z[k2]), _mm_mul_pd( vs1y[k1], vs1z[k2]));
         tmp2 =  _mm_add_pd( _mm_mul_pd( vs0y[k1], vs1z[k2]), _mm_mul_pd( vs1y[k1], vs0z[k2]));
         vwp1[k2][k1].v2 = _mm_add_pd( tmp1, _mm_mul_pd( tmp2, oneHalf));

         tmp1 =  _mm_add_pd( _mm_mul_pd( vs0x[k1], vs0z[k2]), _mm_mul_pd( vs1x[k1], vs1z[k2]));
         tmp2 =  _mm_add_pd( _mm_mul_pd( vs0x[k1], vs1z[k2]), _mm_mul_pd( vs1x[k1], vs0z[k2]));
         vwp2[k2][k1].v2 = _mm_add_pd( tmp1, _mm_mul_pd( tmp2, oneHalf));

         tmp1 =  _mm_add_pd( _mm_mul_pd( vs0x[k1], vs0y[k2]), _mm_mul_pd( vs1x[k1], vs1y[k2]));
         tmp2 =  _mm_add_pd( _mm_mul_pd( vs0x[k1], vs1y[k2]), _mm_mul_pd( vs1x[k1], vs0y[k2]));
         vwp3[k2][k1].v2 = _mm_add_pd( tmp1, _mm_mul_pd( tmp2, oneHalf));

      }
    }
    
    // loop by particle on the outside loop
    for ( k = 0; k < 2; k ++ ) {
       
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
