/*****************************************************************************************

Charge conserving current deposition, BG/Q optimized version

*****************************************************************************************/


/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-bgq.h"

#include "vector-bgq.h"
#include "splines-bgq.h"
#include "os-spec-push-bgq.h"

/***********************************************************************
vwl_s1

Compute the longitudinal current weight for 1st order current deposition
***********************************************************************/

inline void vwl_s1( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {
   
  vwl[0] = vec_mul(vqn, vec_sub(vx1, vx0));
    
}


/***********************************************************************
vwl_s2

Compute the longitudinal current weight for 2nd order current deposition.
This expects that the normalization factor has been multiplied by 1/2 i.e.
vqn = q dx / dt / 2
***********************************************************************/

inline void vwl_s2( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {

  vector vtx = vec_add( vx1, vx0 );
  vector vdx = vec_sub( vx1, vx0 );

  vwl[0] = vec_mul( vqn, vec_nmsub( vtx, vdx, vdx ));
  vwl[1] = vec_mul( vqn, vec_madd(  vtx, vdx, vdx ));
    
}


/***********************************************************************
vwl_s3

Compute the longitudinal current weight for 3rd order current deposition.
This expects that the normalization factor has been multiplied by 1/8 i.e.
vqn = q dx / dt / 8
***********************************************************************/

inline void vwl_s3( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {
  
  vector const c2        = vec_splats( 2.      );
  vector const c4_3      = vec_splats( 4.0/3.0 );
  vector const c6        = vec_splats( 6.      );
  vector const c8_3      = vec_splats( 8.0/3.0 );
  
  vector vdx, vtx, vdx3;
   
  vdx = vec_sub(vx1, vx0);
  vtx = vec_mul(c2, vec_add(vx1, vx0));
  vdx3 = vec_sub(vec_mul(vx1, vec_mul(vx1, vx1)), vec_mul(vx0, vec_mul(vx0, vx0)));
  
  vwl[0] = vec_mul(vqn, vec_madd(vdx3, c4_3, vec_nmsub(vtx, vdx, vdx)));
  vwl[1] = vec_mul(vqn, vec_nmsub(vdx3, c8_3, vec_mul(c6, vdx)));
  vwl[2] = vec_mul(vqn, vec_madd(vdx3, c4_3, vec_madd(vtx, vdx, vdx)));
    
}

/***********************************************************************
vwl_s4

Compute the longitudinal current weight for 4th order current deposition.
This expects that the normalization factor has been multiplied by 1/48 i.e.
vqn = q dx / dt / 48
***********************************************************************/

inline void vwl_s4( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {
  
  vector const c1        = vec_splats( 1.                 );
  vector const c2        = vec_splats( 2.                 );
  vector const c23       = vec_splats( 23.                );
  vector const c4        = vec_splats( 4.                 );
  vector const c6        = vec_splats( 6.                 );
  vector const c15       = vec_splats( 15.                );
  
   
  vector vx12 = vec_mul(vx1, vx1);
  vector vx02 = vec_mul(vx0, vx0);
  
  vector vdx = vec_sub(vx1, vx0);
  vector vtx = vec_add(vx1, vx0);
  
  vector vdx2 = vec_sub(vx12, vx02);
  
  vector t0 = vec_madd(vec_add(vx12, vx02), c2, c1);
  vector t1 = vec_madd(vec_nmsub(vx1, vx12, vec_mul(vx02, vx0)), c4, vec_mul(c23, vdx));
  vector t2 = vec_madd(vec_nmsub(vx12, vx12, vec_mul(vx02, vx02)), c6, vec_mul(c15, vdx2));
  
  vwl[0] = vec_mul( vqn, vec_mul(vec_sub(vdx, vdx2), vec_nmsub(vtx, c2, t0))); 
  vwl[1] = vec_mul( vqn, vec_sub( t1, t2 )) ;
  vwl[2] = vec_mul( vqn, vec_add( t1, t2 )) ;
  vwl[3] = vec_mul( vqn, vec_mul(vec_add(vdx, vdx2), vec_madd(vtx, c2, t0))); 
    
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

void DEP_CURRENT_2D 
(t_real * const current, int const * const size, int const * const offset, 
                       double * const norm, t_vpbuf2D * const part)
{
 
  typedef struct { t_real j1, j2, j3; } t_current;
  t_current *pj, *pjpart, *p0; 
  
  double *vp;
  
  
  vector vnorm1, vnorm2;
  vector vwp1[NP], vwp2[NP]; 
  vector vwl1[ORDER], vwl2[ORDER];

  vector j3[NP][NP];
  
  vector vx0, vx1, vy0, vy1, vq, vvz;
  
  DECLARE_ALIGNED_32( int ix[4] );
  DECLARE_ALIGNED_32( int iy[4] );
		  
  vector vqnx, vqny, vqvz;
  vector vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 

  
  int np;
  int i, k, k1, k2;
  int Dx,Dy;
  
  vector const oneThird = vec_splats( 1.0/3.0 );
  vector const oneHalf  = vec_splats( 0.5 );
  
  
  pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                ( offset[1] - OFFSET ) * 3 * size[0] );

  Dy = size[0];

  // Setting this to 1/4 removes 1 multiplication from wl2
  vnorm1 = vec_splats(norm[0]);
  vnorm2 = vec_splats(norm[1]);
  
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

  vp = (double *) (&part -> buf[0]);
  
  // Each virtual particle uses 8 doubles, so 4 vp = 32 doubles
  for( i = 0; i < np; i += 4, vp += 32 ) {
    
    // load 4 particles
    LOAD4VP2D( vp, vx0, vx1, vy0, vy1, vq, vvz, ix, iy );
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = vec_mul( vq, vnorm1 );
    vqny = vec_mul( vq, vnorm2 );
    vqvz = vec_mul( vq, vvz);
    vqvz = vec_mul( vqvz, oneThird );
 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 );
    WL( vqny, vy0, vy1, vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k] = vec_add(vs0y[k], vs1y[k]);
      vwp2[k] = vec_add(vs0x[k], vs1x[k]);
    }
     
    // get j3 current
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
        vector s00, s10, tmp1, tmp2;        
		
		s00  = vec_mul( vs0x[k1], vs0y[k2] );
		tmp1 = vec_madd( vs1y[k2], vs1x[k1], s00 );  // tmp1 = s0x*s0y + s1x*s1y
	   
		s10  = vec_mul( vs0x[k1], vs1y[k2] );
		tmp2 = vec_madd( vs0y[k2], vs1x[k1], s10 );  // tmp2 = s0x*s1y + s1x*s0y
		
		tmp1 = vec_madd( tmp2, oneHalf, tmp1 );       // tmp1 = tmp1 + 0.5*tmp2
		
		j3[k1][k2] = vec_mul( vqvz, tmp1 );          // j3 = vqvz * tmp1
      }
    }

    // New version loop by particle on the outside loop
    
    for ( k = 0; k < 4; k ++ ) {
    
	   pjpart = pj + ix[k] + iy[k] * Dy;
   
	   // accumulate j1
	   for( k2 = 0; k2 < NP; k2++ ) {
		 for ( k1 = 0; k1 < ORDER; k1++ ) {
		   pjpart[k2*Dy + k1].j1 += (vwl1[k1])[k] * (vwp1[k2])[k];
		 }
	   }
   
	   // accumulate j2 - making k2 the outside loop gives marginal perf. gain
	   for( k2 = 0; k2 < ORDER; k2++ ) {
		 for ( k1 = 0; k1 < NP; k1++ ) {
			pjpart[k2*Dy + k1].j2 += (vwl2[k2])[k] * (vwp2[k1])[k];
		 }
	   }
   
	   // accumulate j3
	   for( k2 = 0; k2 < NP; k2++ ) {
		 for ( k1 = 0; k1 < NP; k1++ ) {
			pjpart[k2*Dy + k1].j3 += (j3[k1][k2])[k];
		 }
	   }
	}

	  
  }
  
}

void DEP_CURRENT_3D
(t_real * const current, int const * const size, int const * const offset, 
                              double * const norm, t_vpbuf3D * const part)
{
        
  typedef struct { t_real j1, j2, j3; } t_current;
  t_current *pj, *pjpart, *p0; 
  double *vp;
  
  vector vnorm1, vnorm2, vnorm3;
  
  vector vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  vector vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];
  
  vector vx0, vx1, vy0, vy1,  vz0, vz1, vq; 
  DECLARE_ALIGNED_32( int ix[4] );
  DECLARE_ALIGNED_32( int iy[4] );
  DECLARE_ALIGNED_32( int iz[4] );
  
  vector vqnx, vqny, vqnz;
  vector vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 

  const vector oneHalf  = vec_splats( 0.5);
  
  
  int np;
  int i, k, k1, k2, k3;
  int Dx,Dy,Dz;
  
  Dy = size[0];    
  Dz = Dy * size[1]; 
  

  // The particle cell indexes in t_vpbuf3D are indexed to 1 	(subtract 1)
  // The deposition starts at current[*][ix-1][ij-1][ik-1] 	    (subtract 1)
  // 												Total :  subtract 2
  

  pj = (t_current *) ( current + ( offset[0] - OFFSET ) * 3 + 
                                 ( offset[1] - OFFSET ) * 3 * Dy + 
                                 ( offset[2] - OFFSET ) * 3 * Dz );


  // norm[i] = dx[i] / dt / 6
  vnorm1 = vec_splats( norm[0]);
  vnorm2 = vec_splats( norm[1]);
  vnorm3 = vec_splats( norm[2]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np % 4 != 0 ) {
    for( i = 0; i < 4 - np%4; i ++ ) {
       t_vp3D * vp = (t_vp3D *) part -> p + i;
       
       vp -> x0 = vp -> x1 = 0.;
       vp -> y0 = vp -> y1 = 0.;
       vp -> z0 = vp -> z1 = 0.;
       vp -> q = 0.;
       vp -> ix = vp -> iy = vp -> iz = 1;
    }
  }

  vp = (double *) part -> buf;
  
  // Each virtual particle uses 12 doubles, so 4 vp = 48 doubles
  for( i = 0; i < np; i += 4, vp += 48 ) {
      
	// load 2 particles
    LOAD4VP3D( vp, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz );

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = vec_mul( vq, vnorm1 );
    vqny = vec_mul( vq, vnorm2 );
    vqnz = vec_mul( vq, vnorm3 );
    
    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 );
    WL( vqny, vy0, vy1, vwl2 );
    WL( vqnz, vz0, vz1, vwl3 );

    
    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
	  for( k1 = 0; k1 < NP; k1++ ) {
 
         vector tmp1, tmp2;        
        
         // wp1[k1][k2] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] + 
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )
         
         tmp1 = vec_mul( vs0y[k1], vs0z[k2]);
         tmp2 = vec_mul( vs0y[k1], vs1z[k2] );
         tmp1 = vec_madd( vs1z[k2], vs1y[k1], tmp1 );
         tmp2 = vec_madd( vs0z[k2], vs1y[k1], tmp2 );

         vwp1[k2][k1] = vec_madd( tmp2, oneHalf, tmp1 );

         // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] + 
         //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

         tmp1 = vec_mul( vs0x[k1], vs0z[k2] );
         tmp2 = vec_mul( vs0x[k1], vs1z[k2] );
         tmp1 = vec_madd( vs1z[k2], vs1x[k1], tmp1 );
         tmp2 = vec_madd( vs0z[k2], vs1x[k1], tmp2 );
         vwp2[k2][k1] = vec_madd( tmp2, oneHalf, tmp1 );

         // wp2[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] + 
         //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

         tmp1 = vec_mul( vs0x[k1], vs0y[k2] );
         tmp2 = vec_mul( vs0x[k1], vs1y[k2] );
         tmp1 = vec_madd( vs1y[k2], vs1x[k1], tmp1 );
         tmp2 = vec_madd(vs0y[k2],  vs1x[k1], tmp2 );
         vwp3[k2][k1] = vec_madd( tmp2, oneHalf, tmp1 );

      }
    }
    
    // loop by particle on the outside loop
   for ( k = 0; k < 4; k ++ ) {
    
	   pjpart = pj + ix[k] + iy[k] * Dy + iz[k] * Dz;
		  
	   // accumulate j1
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < ORDER; k1++ ) {
			  p0[k1].j1 += (vwl1[k1])[k] * (vwp1[k3][k2])[k];
		   }
		 }
	   }
   
	   // accumulate j2
	   for( k3 = 0; k3 < NP; k3++ ) {
		 for( k2 = 0; k2 < ORDER; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < NP; k1++ ) {
			 p0[k1].j2 += (vwl2[k2])[k] * (vwp2[k3][k1])[k];
		   }
		 }
	   }
   
	   
	   // accumulate j3
	   for( k3 = 0; k3 < ORDER; k3++ ) {
		 for( k2 = 0; k2 < NP; k2++ ) {
		   p0 = pjpart + k2*Dy + k3*Dz;
		   for ( k1 = 0; k1 < NP; k1++ ) {
			 p0[k1].j3  += (vwl3[k3])[k] * (vwp3[k2][k1])[k];
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







