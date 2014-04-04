/*****************************************************************************************

Charge conserving current deposition, BG/P optimized version

*****************************************************************************************/


/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-bgp.h"

#include "vector-bgp.h"
#include "splines-bgp.h"
#include "os-spec-push-bgp.h"



/***********************************************************************
vwl_s1

Compute the longitudinal current weight for 1st order current deposition
***********************************************************************/

inline void vwl_s1( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {
   
  vwl[0] = __fpmul(vqn, __fpsub(vx1, vx0));
    
}


/***********************************************************************
vwl_s2

Compute the longitudinal current weight for 2nd order current deposition.
This expects that the normalization factor has been multiplied by 1/2 i.e.
vqn = q dx / dt / 2
***********************************************************************/

inline void vwl_s2( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {

  vector vtx = __fpadd( vx1, vx0 );
  vector vdx = __fpsub( vx1, vx0 );

  vwl[0] = __fpmul( vqn, __fpnmsub( vdx, vdx, vtx ));
  vwl[1] = __fpmul( vqn, __fpmadd(  vdx, vdx, vtx ));
    
}


/***********************************************************************
vwl_s3

Compute the longitudinal current weight for 3rd order current deposition.
This expects that the normalization factor has been multiplied by 1/8 i.e.
vqn = q dx / dt / 8
***********************************************************************/

inline void vwl_s3( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {
  
  vector const c2        = __cmplx( 2.      , 2.      );
  vector const c4_3      = __cmplx( 4.0/3.0 , 4.0/3.0 );
  vector const c6        = __cmplx( 6.      , 6.      );
  vector const c8_3      = __cmplx( 8.0/3.0 , 8.0/3.0 );
  
  vector vdx, vtx, vdx3;
   
  vdx = __fpsub(vx1, vx0);
  vtx = __fpmul(c2, __fpadd(vx1, vx0));
  vdx3 = __fpsub(__fpmul(vx1, __fpmul(vx1, vx1)), __fpmul(vx0, __fpmul(vx0, vx0)));
  
  vwl[0] = __fpmul(vqn, __fpmadd(__fpnmsub(vdx, vdx, vtx), c4_3, vdx3));
  vwl[1] = __fpmul(vqn, __fpnmsub(__fpmul(c6, vdx), c8_3, vdx3));
  vwl[2] = __fpmul(vqn, __fpmadd(__fpmadd(vdx, vdx, vtx), c4_3, vdx3));
    
}

/***********************************************************************
vwl_s4

Compute the longitudinal current weight for 4th order current deposition.
This expects that the normalization factor has been multiplied by 1/48 i.e.
vqn = q dx / dt / 48
***********************************************************************/

inline void vwl_s4( vector const vqn, vector const vx0, vector const vx1, vector vwl[] ) {
  
  vector const c1        = __cmplx( 1.                 , 1.                 );
  vector const c2        = __cmplx( 2.                 , 2.                 );
  vector const c23       = __cmplx( 23.                , 23.                );
  vector const c4        = __cmplx( 4.                 , 4.                 );
  vector const c6        = __cmplx( 6.                 , 6.                 );
  vector const c15       = __cmplx( 15.                , 15.                );
  
   
  vector vx12 = __fpmul(vx1, vx1);
  vector vx02 = __fpmul(vx0, vx0);
  
  vector vdx = __fpsub(vx1, vx0);
  vector vtx = __fpadd(vx1, vx0);
  
  vector vdx2 = __fpsub(vx12, vx02);
  
  vector t0 = __fpmadd(c1, c2, __fpadd(vx12, vx02));
  vector t1 = __fpmadd(__fpmul(c23, vdx), c4, __fpnmsub(__fpmul(vx02, vx0), vx12, vx1));
  vector t2 = __fpmadd(__fpmul(c15, vdx2) , c6, __fpnmsub(__fpmul(vx02, vx02), vx12, vx12));
  
  vwl[0] = __fpmul( vqn, __fpmul(__fpsub(vdx, vdx2), __fpnmsub(t0, c2, vtx))); 
  vwl[1] = __fpmul( vqn, __fpsub( t1, t2 )) ;
  vwl[2] = __fpmul( vqn, __fpadd( t1, t2 )) ;
  vwl[3] = __fpmul( vqn, __fpmul(__fpadd(vdx, vdx2), __fpmadd(t0, c2, vtx))); 
    
}


/****************************************************************************************

  Generate specific current deposition functions for 2D and 3D, 1st to 4th order

****************************************************************************************/

#define __TEMPLATE__

// These macros will append _s1, _s2, etc to function names.
#define ONAME(f, o) OJOIN(f, o)
#define OJOIN(f, o) f ## _s ## o

// This macro handles data type conversion if required
#ifdef PRECISION_SINGLE
#define VPRI(a) __crealf(a)
#define VSEC(a) __cimagf(a)
#else
#define VPRI(a) __creal(a)
#define VSEC(a) __cimag(a)
#endif


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

void extern inline DEP_CURRENT_2D
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
  
  int2 ix, iy;

		  
  vector vqnx, vqny, vqvz;
  vector vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP]; 

  
  int np;
  int i, k, k1, k2;
  int Dx,Dy;
  
  vector const oneThird = __cmplx( 1.0/3.0, 1.0/3.0 );
  vector const oneHalf  = __cmplx( 0.5, 0.5 );
  
  
  pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 + 
                                ( offset[1] - OFFSET ) * 3 * size[0] );

  Dy = size[0];

  // Setting this to 1/4 removes 1 multiplication from wl2
  vnorm1 = __cmplx(norm[0],norm[0]);
  vnorm2 = __cmplx(norm[1],norm[1]);
  
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

  vp = (double *) (&part -> buf[0]);
  
  // Each virtual particle uses 8 doubles, so 2 vp = 16 doubles
  
  for( i = 0; i < np; i+=2, vp+=16 ) {
    
    // load 2 particles
    LOAD2VP2D( vp, vx0, vx1, vy0, vy1, vq, vvz, ix, iy );
    
    // Get splines    
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    
    vqnx = __fpmul( vq, vnorm1 );
    vqny = __fpmul( vq, vnorm2 );
    vqvz = __fpmul( vq, vvz);
    vqvz = __fpmul( vqvz, oneThird );
 
    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 );
    WL( vqny, vy0, vy1, vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k] = __fpadd(vs0y[k], vs1y[k]);
      vwp2[k] = __fpadd(vs0x[k], vs1x[k]);
    }
     
    // get j3 current
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
        vector s00, s10, tmp1, tmp2;        
		
		s00  = __fpmul( vs0x[k1], vs0y[k2] );
		tmp1 = __fpmadd( s00, vs1x[k1], vs1y[k2] );  // tmp1 = s0x*s0y + s1x*s1y
	   
		s10  = __fpmul( vs0x[k1], vs1y[k2] );
		tmp2 = __fpmadd( s10, vs1x[k1], vs0y[k2] );  // tmp2 = s0x*s1y + s1x*s0y
		
		tmp1 = __fpmadd( tmp1, oneHalf, tmp2 );       // tmp1 = tmp1 + 0.5*tmp2
		
		j3[k1][k2] = __fpmul( vqvz, tmp1 );          // j3 = vqvz * tmp1
      }
    }

    // New version loop by particle on the outside loop
    
    // Particle (A)

	pjpart = pj + ix.a + iy.a * Dy;

	// accumulate j1
	for( k2 = 0; k2 < NP; k2++ ) {
	  for ( k1 = 0; k1 < ORDER; k1++ ) {
		pjpart[k2*Dy + k1].j1 += (VPRI(vwl1[k1]) * VPRI(vwp1[k2]));
	  }
	}

	// accumulate j2 - making k2 the outside loop gives marginal perf. gain
	for( k2 = 0; k2 < ORDER; k2++ ) {
	  for ( k1 = 0; k1 < NP; k1++ ) {
		 pjpart[k2*Dy + k1].j2 += (VPRI(vwl2[k2]) * VPRI(vwp2[k1]));
	  }
	}

	// accumulate j3
	for( k2 = 0; k2 < NP; k2++ ) {
	  for ( k1 = 0; k1 < NP; k1++ ) {
		 pjpart[k2*Dy + k1].j3 += (VPRI(j3[k1][k2]));
	  }
	}

    // Particle (B)

	pjpart = pj + ix.b + iy.b * Dy;

	// accumulate j1
	for( k2 = 0; k2 < NP; k2++ ) {
	  for ( k1 = 0; k1 < ORDER; k1++ ) {
		pjpart[k2*Dy + k1].j1 += (VSEC(vwl1[k1]) * VSEC(vwp1[k2]));
	  }
	}

	// accumulate j2 - making k2 the outside loop gives marginal perf. gain
	for( k2 = 0; k2 < ORDER; k2++ ) {
	  for ( k1 = 0; k1 < NP; k1++ ) {
		 pjpart[k2*Dy + k1].j2 += (VSEC(vwl2[k2]) * VSEC(vwp2[k1]));
	  }
	}

	// accumulate j3
	for( k2 = 0; k2 < NP; k2++ ) {
	  for ( k1 = 0; k1 < NP; k1++ ) {
		 pjpart[k2*Dy + k1].j3 += (VSEC(j3[k1][k2]));
	  }
	}
	  
  }
  
}

void extern inline  DEP_CURRENT_3D
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
  int2 ix, iy, iz;
  
  vector vqnx, vqny, vqnz;
  vector vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP]; 

  const vector oneHalf  = __cmplx( 0.5, 0.5 );
  
  
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
  vnorm1 = __cmplx( norm[0], norm[0]);
  vnorm2 = __cmplx( norm[1], norm[1]);
  vnorm3 = __cmplx( norm[2], norm[2]);
  
  np = part -> np;
  
  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np % 2 != 0 ) {
    for( i = 0; i < 2 - np%2; i ++ ) {
       t_vp3D * vp = (t_vp3D *) part -> p + i;
       
       vp -> x0 = vp -> x1 = 0.;
       vp -> y0 = vp -> y1 = 0.;
       vp -> z0 = vp -> z1 = 0.;
       vp -> q = 0.;
       vp -> ix = vp -> iy = vp -> iz = 1;
    }
  }

  vp = (double *) part -> buf;
  
  // Each virtual particle uses 10 doubles, so 2 vp = 20 doubles
  
  for( i = 0; i < np; i+=2, vp+=20 ) {
      
	// load 2 particles
    LOAD2VP3D( vp, vx0, vx1, vy0, vy1, vz0, vz1, vq, ix, iy, iz );

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );
    
    vqnx = __fpmul( vq, vnorm1 );
    vqny = __fpmul( vq, vnorm2 );
    vqnz = __fpmul( vq, vnorm3 );
    
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
         
         tmp1 = __fpmul( vs0y[k1], vs0z[k2]);
         tmp2 = __fpmul( vs0y[k1], vs1z[k2] );
         tmp1 = __fpmadd( tmp1, vs1y[k1], vs1z[k2] );
         tmp2 = __fpmadd( tmp2, vs1y[k1], vs0z[k2] );

         vwp1[k2][k1] = __fpmadd( tmp1, oneHalf, tmp2 );

         // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] + 
         //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

         tmp1 = __fpmul( vs0x[k1], vs0z[k2] );
         tmp2 = __fpmul( vs0x[k1], vs1z[k2] );
         tmp1 = __fpmadd( tmp1, vs1x[k1], vs1z[k2] );
         tmp2 = __fpmadd( tmp2, vs1x[k1], vs0z[k2] );
         vwp2[k2][k1] = __fpmadd( tmp1, oneHalf, tmp2 );

         // wp2[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] + 
         //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

         tmp1 = __fpmul( vs0x[k1], vs0y[k2] );
         tmp2 = __fpmul( vs0x[k1], vs1y[k2] );
         tmp1 = __fpmadd( tmp1, vs1x[k1], vs1y[k2] );
         tmp2 = __fpmadd(tmp2,  vs1x[k1], vs0y[k2] );
         vwp3[k2][k1] = __fpmadd( tmp1, oneHalf, tmp2 );

      }
    }
    
    // loop by particle on the outside loop
    
    // Particle (A)

	pjpart = pj + ix.a + iy.a * Dy + iz.a * Dz;
	   
	// accumulate j1
	for( k3 = 0; k3 < NP; k3++ ) {
	  for( k2 = 0; k2 < NP; k2++ ) {
		p0 = pjpart + k2*Dy + k3*Dz;
		for ( k1 = 0; k1 < ORDER; k1++ ) {
		   p0[k1].j1 += (VPRI(vwl1[k1]) * VPRI(vwp1[k3][k2]));
		}
	  }
	}

	// accumulate j2
	for( k3 = 0; k3 < NP; k3++ ) {
	  for( k2 = 0; k2 < ORDER; k2++ ) {
		p0 = pjpart + k2*Dy + k3*Dz;
		for ( k1 = 0; k1 < NP; k1++ ) {
		  p0[k1].j2 += (VPRI(vwl2[k2]) * VPRI(vwp2[k3][k1]));
		}
	  }
	}

	
	// accumulate j3
	for( k3 = 0; k3 < ORDER; k3++ ) {
	  for( k2 = 0; k2 < NP; k2++ ) {
		p0 = pjpart + k2*Dy + k3*Dz;
		for ( k1 = 0; k1 < NP; k1++ ) {
		  p0[k1].j3  += (VPRI(vwl3[k3]) * VPRI(vwp3[k2][k1]));
		}
	  }
	}

    // Particle (B)

	pjpart = pj + ix.b + iy.b * Dy + iz.b * Dz;
	   
	// accumulate j1
	for( k3 = 0; k3 < NP; k3++ ) {
	  for( k2 = 0; k2 < NP; k2++ ) {
		p0 = pjpart + k2*Dy + k3*Dz;
		for ( k1 = 0; k1 < ORDER; k1++ ) {
		   p0[k1].j1 += (VSEC(vwl1[k1]) * VSEC(vwp1[k3][k2]));
		}
	  }
	}

	// accumulate j2
	for( k3 = 0; k3 < NP; k3++ ) {
	  for( k2 = 0; k2 < ORDER; k2++ ) {
		p0 = pjpart + k2*Dy + k3*Dz;
		for ( k1 = 0; k1 < NP; k1++ ) {
		  p0[k1].j2 += (VSEC(vwl2[k2]) * VSEC(vwp2[k3][k1]));
		}
	  }
	}

	
	// accumulate j3
	for( k3 = 0; k3 < ORDER; k3++ ) {
	  for( k2 = 0; k2 < NP; k2++ ) {
		p0 = pjpart + k2*Dy + k3*Dz;
		for ( k1 = 0; k1 < NP; k1++ ) {
		  p0[k1].j3  += (VSEC(vwl3[k3]) * VSEC(vwp3[k2][k1]));
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







