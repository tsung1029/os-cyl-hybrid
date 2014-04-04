/*****************************************************************************************

Relativistic particle pusher, AVX optimized version (double precision)

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__


#include "os-spec-push-avxd.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "fortran.h"

#include "vector-avx.h"
#include "splines-avx.h"
#include "os-spec-current-avxd.h"


/*****************************************************************************************
Single pass particle advance deposit

This version of the code does the full push in a single pass: field interpolation, du/dt and dx/dt. 
This allows for savings on memory accesses, since particle positions need only be read once and it 
avoids storing temporary field and momenta values. 

Each group of 4 particles is immediatly split to avoid storing old position values in a special
buffer optimized for the current deposition routines.

The routine then deposits the current for all virtual (split) particles.
*****************************************************************************************/


/*****************************************************************************************
vmntrim

Returns the required shift (trim) for particles that must remain in [ -0.5, 0.5 [
to the nearest cell point:
          x >= 0.5 : +1
  -0.5 <= x <  0.5 :  0
  -0.5 <  x        : -1

*****************************************************************************************/

static inline __m256d vmntrim( const __m256d vx )
{
   register __m256d va, vb;
   
   va = _mm256_cmp_pd( vx, _mm256_set1_pd( -0.5 ), _CMP_LT_OS );
   va = _mm256_and_pd( va, _mm256_set1_pd( +1.0 ) );
   
   vb = _mm256_cmp_pd( vx, _mm256_set1_pd( +0.5 ), _CMP_GE_OS );
   vb = _mm256_and_pd( vb, _mm256_set1_pd( +1.0 ) );
   
   return  _mm256_sub_pd( vb, va );
}


/*****************************************************************************************
LOADFLD4

Loads 4 field values corresponding to 4 particle positions into a vector
variable.
*****************************************************************************************/


#define LOADFLD4( fp, shift ) \
  _mm256_set_pd( (fp[3])[shift], (fp[2])[shift], (fp[1])[shift], (fp[0])[shift] )

/*****************************************************************************************
vget_hp_near

Gets the nearest half points of the 2 particles loaded as a vector for interpolating 
scattered grids. This routine is for positions defined as distance to the nearest grid
points.
*****************************************************************************************/


#define vget_hp_near(dx, ix, dxh, ixh, delta ) 	{ \
   __m256d cmp  = _mm256_cmp_pd( (dx), _mm256_setzero_pd(), _CMP_LT_OS ); \
   (dxh) =  _mm256_add_pd( dx, _mm256_blendv_pd( _mm256_set1_pd( -0.5 ), \
											     _mm256_set1_pd(  0.5 ), cmp ) ); \
   cmp = _mm256_castps_pd( _mm256_permute_ps( _mm256_castpd_ps( cmp ), 0x08 ) ); \
   __m128i cmpi = _mm_castpd_si128( \
			  _mm_unpacklo_pd( _mm256_castpd256_pd128( cmp  ), \
							   _mm256_extractf128_pd( cmp , 1 ) ) ); \
   (ixh)  = _mm_sub_epi32( (ix), _mm_and_si128( cmpi,_mm_set1_epi32( delta )) ); \
}


/*****************************************************************************************
vdudt_boris

Use a Boris push to advance particle velocities using interpolated
fields.
*****************************************************************************************/

#define VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, u, vu1, vu2, vu3 )                   \
{                                                                                             \
   __m256d vut1, vut2, vut3;\
   __m256d tmp1, tmp2, tmp3;\
   __m256d vgamma_tem, votsq;\
   \
   /* Perform first half of electric field acceleration.*/\
   ve1 = _mm256_mul_pd( ve1, vtem );\
   ve2 = _mm256_mul_pd( ve2, vtem );\
   ve3 = _mm256_mul_pd( ve3, vtem );\
   \
   _MM256_LOAD4v3_PD( vu1, vu2, vu3, u )\
   \
   vut1 = _mm256_add_pd( vu1, ve1 );\
   vut2 = _mm256_add_pd( vu2, ve2 );\
   vut3 = _mm256_add_pd( vu3, ve3 );\
   \
   /* Perform first half of the rotation and store in u. */\
   tmp1 = _mm256_mul_pd( vut1, vut1 );\
   tmp2 = _mm256_mul_pd( vut2, vut2 );\
   tmp3 = _mm256_mul_pd( vut3, vut3 );\
   \
   vgamma_tem = _mm256_add_pd( tmp1, tmp2 );\
   vgamma_tem = _mm256_add_pd( vgamma_tem, tmp3 );\
   vgamma_tem = _mm256_add_pd( vgamma_tem, _mm256_set1_pd( 1.0 ) );\
   \
   vgamma_tem = _mm256_sqrt_pd( vgamma_tem );\
   vgamma_tem = _mm256_div_pd( vtem, vgamma_tem );\
   \
   vb1 = _mm256_mul_pd( vb1, vgamma_tem );\
   vb2 = _mm256_mul_pd( vb2, vgamma_tem );\
   vb3 = _mm256_mul_pd( vb3, vgamma_tem );\
   \
   vu1 = _mm256_add_pd( vut1, _mm256_sub_pd( _mm256_mul_pd( vut2, vb3 ), _mm256_mul_pd( vut3, vb2 ) ) );\
   vu2 = _mm256_add_pd( vut2, _mm256_sub_pd( _mm256_mul_pd( vut3, vb1 ), _mm256_mul_pd( vut1, vb3 ) ) );\
   vu3 = _mm256_add_pd( vut3, _mm256_sub_pd( _mm256_mul_pd( vut1, vb2 ), _mm256_mul_pd( vut2, vb1 ) ) );\
   \
   /* Perform second half of the rotation. */\
   tmp1 = _mm256_mul_pd( vb1, vb1 );\
   tmp2 = _mm256_mul_pd( vb2, vb2 );\
   tmp3 = _mm256_mul_pd( vb3, vb3 );\
   \
   tmp1 = _mm256_add_pd( tmp1, tmp2 );\
   tmp1 = _mm256_add_pd( tmp1, tmp3 );\
   tmp1 = _mm256_add_pd( tmp1, _mm256_set1_pd(1.0));\
   \
   votsq = _mm256_div_pd( _mm256_set1_pd(2.0), tmp1 );\
   \
   vb1 = _mm256_mul_pd( vb1, votsq );\
   vb2 = _mm256_mul_pd( vb2, votsq );\
   vb3 = _mm256_mul_pd( vb3, votsq );\
   \
   vut1 = _mm256_add_pd( vut1, _mm256_sub_pd( _mm256_mul_pd( vu2, vb3 ), _mm256_mul_pd( vu3, vb2 ) ) );\
   vut2 = _mm256_add_pd( vut2, _mm256_sub_pd( _mm256_mul_pd( vu3, vb1 ), _mm256_mul_pd( vu1, vb3 ) ) );\
   vut3 = _mm256_add_pd( vut3, _mm256_sub_pd( _mm256_mul_pd( vu1, vb2 ), _mm256_mul_pd( vu2, vb1 ) ) );\
   \
   /* Perform second half of electric field acceleration.*/\
   vu1 = _mm256_add_pd( vut1, ve1 );\
   vu2 = _mm256_add_pd( vut2, ve2 );\
   vu3 = _mm256_add_pd( vut3, ve3 );\
\
}

/*****************************************************************************************
VRGAMMA

Get 1 / Lorentz gamma. Using a macro is beneficial because this is 
called from a function that is already inlined.
*****************************************************************************************/

#define VRGAMMA( vrg, vu1, vu2, vu3 ) {  \
  (vrg) = _mm256_add_pd( _mm256_add_pd(_mm256_mul_pd((vu1),(vu1)),   \
                                       _mm256_mul_pd((vu2),(vu2))),  \
                                       _mm256_mul_pd((vu3),(vu3)) ); \
  (vrg) = _mm256_add_pd( (vrg), _mm256_set1_pd(1.0) ); \
  (vrg) = _mm256_sqrt_pd( (vrg) ); \
  (vrg) = _mm256_div_pd( _mm256_set1_pd(1.0), (vrg) ); \
}


/*****************************************************************************************

The current splitters expect that the particle paths stored in the vpbuf array are relative to the
cell where the motion started. This means that for particles crossing a cell boundary the values
will be outside the correct range. This is done to improve numerical accuracy since the total 
distance travelled is used as the denominator in a division. For small distances occuring next to
cell crossings, numerical roundoff could lead to divisions by 0 otherwise.

*****************************************************************************************/

/*****************************************************************************************
vsplit2D

Splits 2 particle trajectories for current deposition and stores segment on a special virtual
particle buffer which is optimized for the vector current deposition routines
*****************************************************************************************/

void vsplit2D( t_vpbuf2D* const vpbuf, const int cross[], const int dix[], const int diy[] )
{

  int k, np;
  t_vp2D *vp;

  np = VEC_WIDTH;
  vp = (t_vp2D *) vpbuf -> p;
  
  for( k = 0 ; k < VEC_WIDTH;  k++ ) {
     double delta, xint, yint, vzint, xint2, yint2, vzint2;
    
     switch ( cross[k] ) {
       case (0) :  // no split

         // no action required
		 break;
         
       case (1) :  //  x cross only

		 xint  = 0.5 * dix[k];
		 delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);
		 yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta;

         // add extra particle (2nd segment)
         vp[np].x0 = -xint; 	vp[np].x1 = vp[k].x1 - dix[k];
         vp[np].y0 = yint; 		vp[np].y1 = vp[k].y1;
         vp[np].ix = vp[k].ix + dix[k]; 
         vp[np].iy = vp[k].iy;
         vp[np].q  = vp[k].q;  
         vp[np].vz = vp[k].vz * (1-delta);
         np++;
         
         // correct existing particle (1st segment)
         vp[k].x1 = xint;
         vp[k].y1 = yint;
         vp[k].vz *= delta;

		 break;
          
       case (2) :  // y cross only

		 yint  = 0.5 * ( diy[k] );
		 delta = ( yint - vp[k].y0 ) / (vp[k].y1-vp[k].y0);
		 xint  = vp[k].x0 + (vp[k].x1 - vp[k].x0) * delta;

         // add extra particle (2nd segment)
         vp[np].x0 =  xint; 	vp[np].x1 = vp[k].x1;
         vp[np].y0 = -yint; 	vp[np].y1 = vp[k].y1  - diy[k];
         vp[np].ix = vp[k].ix; 
         vp[np].iy = vp[k].iy + diy[k];
         vp[np].q  = vp[k].q;  
         vp[np].vz = vp[k].vz * (1-delta);
         np++;

         // correct existing particle (1st segment)         	 
		 vp[k].x1 = xint;
		 vp[k].y1 = yint;
		 vp[k].vz *= delta;
		 
		 break;
       
       case (3) : // x-y cross
		 // split in x direction first
		 xint  = 0.5 * dix[k];
		 delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);
		 yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta;

         if ( ( yint >= -0.5 ) && (yint < 0.5) ) {

           // no y cross on 1st vp
           vzint = vp[k].vz * delta;

           // y split 2nd vp
           vzint2 = vp[k].vz * (1-delta);
           yint2 = 0.5 * diy[k];
		   delta = ( yint2 - yint ) / (vp[k].y1 - yint);
		   xint2 = -xint + ( vp[k].x1 - xint ) * delta;

           // add extra particle (1st segment y-split)
		   vp[np].x0 =  -xint; 	vp[np].x1 = xint2;
		   vp[np].y0 = yint; 	vp[np].y1 = yint2;
		   vp[np].ix = vp[k].ix + dix[k]; 
		   vp[np].iy = vp[k].iy;
		   vp[np].q  = vp[k].q;  
		   vp[np].vz = vzint2 * delta;
		   np++;


           // add extra particle (2nd segment y-split)
		   vp[np].x0 =  xint2; 	vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 = -yint2; 	vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].ix = vp[k].ix + dix[k]; 
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].q  = vp[k].q;  
		   vp[np].vz = vzint2 * (1-delta);
		   np++;

		   // correct existing particle (1st segment)  
		   vp[k].x1 = xint;
		   vp[k].y1 = yint;
		   vp[k].vz = vzint;

         } else {
         
           vzint = vp[k].vz*delta;
		   // y split 1st vp
		   yint2 = 0.5 * diy[k];
		   delta = ( yint2 - vp[k].y0 ) / ( yint - vp[k].y0);
		   xint2 = vp[k].x0 + (xint - vp[k].x0) * delta;

           // add extra particle (2nd segment y-split)
		   vp[np].x0 =  xint2; vp[np].x1 = xint;
		   vp[np].y0 = -yint2; vp[np].y1 = yint - diy[k];
		   vp[np].ix = vp[k].ix; 
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].q  = vp[k].q;  
		   vp[np].vz = vzint * (1-delta);
		   np++;
 
		   // no y cross on last particle
		   vp[np].x0 = -xint;          vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 = yint - diy[k];  vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].q  = vp[k].q;
		   vp[np].vz = vp[k].vz - vzint;
		   np++;

		   // correct existing particle (1st segment y-split)  
		   vp[k].x1 = xint2;
		   vp[k].y1 = yint2;
		   vp[k].vz = vzint * delta;

         }
		 break;
     }
     
  }

  vpbuf -> p = (double *) (vp + np);
  vpbuf -> np += np;
}


/*****************************************************************************************
vsplit3D

Splits 2 particle trajectories for current deposition and stores segment on a special virtual
particle buffer which is optimized for the vector current deposition routines
*****************************************************************************************/

void vsplit3D( t_vpbuf3D* const vpbuf, const int cross[], 
               const int dix[], const int diy[], const int diz[] )
{

  int k, np;
  t_vp3D *vp;

  np = VEC_WIDTH;
  vp = (t_vp3D *) vpbuf -> p;

  for( k = 0; k < VEC_WIDTH; k++ ) {
     double xint, yint, zint;
     double xint2, yint2, zint2;
     double delta;
        
     switch ( cross[k] ) {
       case (0) :  // no split

         // no action required
		 break;
         
       case (1) :  //  x cross only
         
		 xint  = 0.5 * dix[k];
		 delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);
		 yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta;
		 zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta;

         // add extra particle (2nd segment)
         vp[np].x0 = -xint; 		vp[np].x1 = vp[k].x1 - dix[k];
         vp[np].y0 =  yint; 		vp[np].y1 = vp[k].y1;
         vp[np].z0 =  zint; 		vp[np].z1 = vp[k].z1;
         vp[np].ix = vp[k].ix + dix[k];
         vp[np].iy = vp[k].iy;
         vp[np].iz = vp[k].iz;
         vp[np].q  = vp[k].q;
         np ++;

         // correct existing particle (1st segment)
         vp[k].x1 = xint; 
         vp[k].y1 = yint; 
         vp[k].z1 = zint;

		 break;
          
       case (2) :  // y cross only

		 yint  = 0.5 * diy[k];
		 delta = ( yint - vp[k].y0 ) / (vp[k].y1 - vp[k].y0);
		 xint  = vp[k].x0 + (vp[k].x1 - vp[k].x0) * delta;
		 zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta;

         // add extra particle (2nd segment)
         vp[np].x0 =  xint; 		vp[np].x1 = vp[k].x1;
         vp[np].y0 = -yint; 		vp[np].y1 = vp[k].y1 - diy[k];
         vp[np].z0 =  zint; 		vp[np].z1 = vp[k].z1;
         vp[np].ix = vp[k].ix;
         vp[np].iy = vp[k].iy + diy[k];
         vp[np].iz = vp[k].iz;
         vp[np].q  = vp[k].q;
         np ++;

         // correct existing particle (1st segment)
         vp[k].x1 = xint; 
         vp[k].y1 = yint; 
         vp[k].z1 = zint;

		 break;
       
       case (3) : // x-y cross
       
		 // split in x direction first
		 xint  = 0.5 * dix[k];
		 delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);

		 // note that this is different from case(1) because of change of cell in y direction
		 yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta;
		 zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta;

         if ( ( yint >= -0.5 ) && (yint < 0.5) ) {
                      
           // y split 2nd vp
           yint2 = 0.5 * diy[k];
		   delta = ( yint2 - yint ) / (vp[k].y1 - yint);
		   xint2 = -xint + ( vp[k].x1 - xint ) * delta;
		   zint2 =  zint + ( vp[k].z1 - zint ) * delta;
		   
		   vp[np].x0 = -xint; 		vp[np].x1 = xint2;
		   vp[np].y0 =  yint; 		vp[np].y1 = yint2;
		   vp[np].z0 =  zint; 		vp[np].z1 = zint2;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

		   vp[np].x0 =  xint2; 		vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 = -yint2; 		vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].z0 =  zint2; 		vp[np].z1 = vp[k].z1;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

           // no y cross on 1st vp
           vp[k].x1 = xint; vp[k].y1 = yint; vp[k].z1 = zint;
           
         } else {
         
		   // y split 1st vp
		   yint2 = 0.5 * (diy[k]);
		   delta = ( yint2 - vp[k].y0 ) / ( yint - vp[k].y0);
		   xint2 = vp[k].x0 + (xint - vp[k].x0) * delta;
		   zint2 = vp[k].z0 + (zint - vp[k].z0) * delta;

		   vp[np].x0 =  xint2; 		vp[np].x1 = xint;
		   vp[np].y0 = -yint2; 		vp[np].y1 = yint - diy[k];
		   vp[np].z0 =  zint2; 		vp[np].z1 = zint;
		   vp[np].ix = vp[k].ix;
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;
					  
		   // no y cross on 2nd vp
		   vp[np].x0 =       -xint;	vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 = yint-diy[k]; vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].z0 =        zint;	vp[np].z1 = vp[k].z1;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

           // Correct 1st vp
           vp[k].x1 = xint2; vp[k].y1 = yint2; vp[k].z1 = zint2;

         }
  
		 break;

       case (4) :  // z cross only

		 zint  = 0.5 * diz[k];
		 delta = ( zint - vp[k].z0 ) / (vp[k].z1-vp[k].z0);
		 xint  = vp[k].x0 + (vp[k].x1 - vp[k].x0) * delta;
		 yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta;

         // add extra particle (2nd segment)
         vp[np].x0 =  xint; 		vp[np].x1 = vp[k].x1;
         vp[np].y0 =  yint; 		vp[np].y1 = vp[k].y1;
         vp[np].z0 = -zint; 		vp[np].z1 = vp[k].z1 - diz[k];
         vp[np].ix = vp[k].ix;
         vp[np].iy = vp[k].iy;
         vp[np].iz = vp[k].iz + diz[k];
         vp[np].q  = vp[k].q;
         np ++;

         // correct existing particle (1st segment)
         vp[k].x1 = xint; vp[k].y1 = yint; vp[k].z1 = zint;
	 
		 break;

       case (5) :  // x-z cross
       
		 // split in x direction first
		 xint  = 0.5*( dix[k] );
		 delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);

		 // note that this is different from case(1) because of change of cell in z direction
		 yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta; // no y cross
		 zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta; 

         if ( ( zint >= -0.5 ) && (zint < 0.5) ) {
                      
           // z split 2nd vp
           zint2 = 0.5 * (diz[k]);
		   delta = ( zint2 - zint ) / (vp[k].z1 - zint);
		   xint2 = -xint + ( vp[k].x1 - xint ) * delta;
		   yint2 =  yint + ( vp[k].y1 - yint ) * delta;
	
		   vp[np].x0 = -xint; 		vp[np].x1 = xint2;
		   vp[np].y0 =  yint; 		vp[np].y1 = yint2;
		   vp[np].z0 =  zint; 		vp[np].z1 = zint2;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

		   vp[np].x0 =  xint2; 		vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 =  yint2; 		vp[np].y1 = vp[k].y1;
		   vp[np].z0 = -zint2; 		vp[np].z1 = vp[k].z1 - diz[k];
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz + diz[k];
		   vp[np].q  = vp[k].q;
		   np ++;

           // no z cross on 1st vp
		   vp[k].x1 = xint;
		   vp[k].y1 = yint;
		   vp[k].z1 = zint;

         } else {
         
		   // z split 1st vp
		   zint2 = 0.5 * (diz[k]);
		   delta = ( zint2 - vp[k].z0 ) / ( zint - vp[k].z0);
		   xint2 = vp[k].x0 + (xint - vp[k].x0) * delta;
		   yint2 = vp[k].y0 + (yint - vp[k].y0) * delta;

		   vp[np].x0 =  xint2; 		vp[np].x1 = xint;
		   vp[np].y0 =  yint2; 		vp[np].y1 = yint;
		   vp[np].z0 = -zint2; 		vp[np].z1 = zint - diz[k];
		   vp[np].ix = vp[k].ix;
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz + diz[k];
		   vp[np].q  = vp[k].q;
		   np ++;
           					  
		   // no y cross on 2nd vp
		   vp[np].x0 = -xint; 		vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 =  yint; 		vp[np].y1 = vp[k].y1;
		   vp[np].z0 = zint-diz[k];	vp[np].z1 = vp[k].z1 - diz[k];
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz + diz[k];
		   vp[np].q  = vp[k].q;
		   np ++;

           // correct existing particle
		   vp[k].x1 = xint2;
		   vp[k].y1 = yint2;
		   vp[k].z1 = zint2;

         }
		 break;       

       case (6) :  // y-z cross
       
		 // split in x direction first
		 yint  = 0.5*( diy[k] );
		 delta = ( yint - vp[k].y0 ) / (vp[k].y1 - vp[k].y0);

		 xint  = vp[k].x0 + (vp[k].x1 - vp[k].x0) * delta; // no x cross
		 zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta; 

         if ( ( zint >= -0.5 ) && (zint < 0.5) ) {
                      
           // z split 2nd vp
           zint2 = 0.5 * diz[k];
		   delta = ( zint2 - zint ) / (vp[k].z1 - zint);
		   xint2 =   xint + ( vp[k].x1 - xint ) * delta;
		   yint2 =  -yint + ( vp[k].y1 - yint ) * delta;
		   
		   vp[np].x0 =  xint; 		vp[np].x1 = xint2;
		   vp[np].y0 = -yint; 		vp[np].y1 = yint2;
		   vp[np].z0 =  zint; 		vp[np].z1 = zint2;
		   vp[np].ix = vp[k].ix;
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

		   vp[np].x0 =  xint2; 		vp[np].x1 = vp[k].x1;
		   vp[np].y0 =  yint2; 		vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].z0 = -zint2;	 	vp[np].z1 = vp[k].z1 - diz[k];
		   vp[np].ix = vp[k].ix;
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz + diz[k];
		   vp[np].q  = vp[k].q;
		   np ++;

           // no z cross on 1st vp
		   vp[k].x1 = xint; 
		   vp[k].y1 = yint;
		   vp[k].z1 = zint;
           
         } else {
         
		   // z split 1st vp
		   zint2 = 0.5 * diz[k];
		   delta = ( zint2 - vp[k].z0 ) / ( zint - vp[k].z0);
		   xint2 = vp[k].x0 + (xint - vp[k].x0) * delta;
		   yint2 = vp[k].y0 + (yint - vp[k].y0) * delta;

		   vp[np].x0 =  xint2; 		vp[np].x1 = xint;
		   vp[np].y0 =  yint2; 		vp[np].y1 = yint;
		   vp[np].z0 = -zint2; 		vp[np].z1 = zint - diz[k];
		   vp[np].ix = vp[k].ix;
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz + diz[k];
		   vp[np].q  = vp[k].q;
		   np ++;
					  
		   // no y cross on 2nd vp
		   vp[np].x0 =  xint; 		vp[np].x1 = vp[k].x1;
		   vp[np].y0 = -yint; 		vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].z0 = zint-diz[k];	vp[np].z1 = vp[k].z1 - diz[k];
		   vp[np].ix = vp[k].ix;
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz + diz[k];
		   vp[np].q  = vp[k].q;
		   np ++;

           // correct existing particle
		   vp[k].x1 = xint2;
		   vp[k].y1 = yint2;
		   vp[k].z1 = zint2;

         }
		 break;       

       case (7) : // x-y-z cross
       {
        
          // do an x,y cross first 
		  xint  = 0.5 * dix[k];
		  delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);
		  yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta;
		  zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta; 

         if ( ( yint >= -0.5 ) && (yint < 0.5) ) {
                      
           // y split 2nd vp
           yint2 = 0.5 * (diy[k]);
		   delta = ( yint2 - yint ) / (vp[k].y1 - yint);
		   xint2 = -xint + ( vp[k].x1 - xint ) * delta;
		   zint2 =  zint + ( vp[k].z1 - zint ) * delta;
		   
		   vp[np].x0 = -xint; 		vp[np].x1 = xint2;
		   vp[np].y0 =  yint; 		vp[np].y1 = yint2;
		   vp[np].z0 =  zint; 		vp[np].z1 = zint2;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

		   vp[np].x0 =  xint2; 		vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 = -yint2; 		vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].z0 =  zint2; 		vp[np].z1 = vp[k].z1;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

           // no y cross on 1st vp
           vp[k].x1 = xint; 
           vp[k].y1 = yint;
           vp[k].z1 = zint;
           
         } else {
         
		   // y split 1st vp
		   yint2 = 0.5 * (diy[k]);
		   delta = ( yint2 - vp[k].y0 ) / ( yint - vp[k].y0);
		   xint2 = vp[k].x0 + (xint - vp[k].x0) * delta;
		   zint2 = vp[k].z0 + (zint - vp[k].z0) * delta;

		   vp[np].x0 =  xint2; 		vp[np].x1 = xint;
		   vp[np].y0 = -yint2; 		vp[np].y1 = yint - diy[k];
		   vp[np].z0 =  zint2; 		vp[np].z1 = zint;
		   vp[np].ix = vp[k].ix;
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;
					  
		   // no y cross on 2nd vp
		   vp[np].x0 =       -xint;	vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 = yint-diy[k]; vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].z0 =        zint;	vp[np].z1 = vp[k].z1;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

           // Correct 1st vp
           vp[k].x1 = xint2; 
           vp[k].y1 = yint2; 
           vp[k].z1 = zint2;

         }
         
 		  // one of the 3 vp requires an additional z split
          /* int vpidx[3] = {k, np-2, np-1}; */
          int vpidx[3]; vpidx[0] = k; vpidx[1] = np-2; vpidx[2] = np-1;
          int i, j, idx;

		  zint = 0.5 * (diz[k]);
		  for ( i = 0; i < 3; i ++ ) {
			idx = vpidx[i];
			if ( (vp[idx].z1 < -0.5) || (vp[idx].z1 >= +0.5) ) {
			   
			   // z1 and z0 are already indexed to the same cell
			   delta = ( zint - vp[idx].z0) / ( vp[idx].z1-vp[idx].z0 );
			   xint  = vp[idx].x0 + (vp[idx].x1 - vp[idx].x0) * delta;
			   yint  = vp[idx].y0 + (vp[idx].y1 - vp[idx].y0) * delta;
	  
			   // store new vp
			   vp[np].x0 =  xint; 		vp[np].x1 = vp[idx].x1;
			   vp[np].y0 =  yint; 		vp[np].y1 = vp[idx].y1;
			   vp[np].z0 = -zint; 		vp[np].z1 = vp[idx].z1 - diz[k];
			   vp[np].ix = vp[idx].ix;
			   vp[np].iy = vp[idx].iy;
			   vp[np].iz = vp[idx].iz + diz[k];
			   vp[np].q  = vp[idx].q;
			   np ++;
			   		   
		       // correct old vp
			   vp[idx].x1 = xint;
			   vp[idx].y1 = yint;
			   vp[idx].z1 = zint; 
			   
			   for ( j = i+1; j < 3; j++ ) {
			     idx = vpidx[j];
			     vp[idx].z0 -= diz[k];
			     vp[idx].z1 -= diz[k];
			     vp[idx].iz += diz[k];
			   }
			   
			   // No further z splits
			   break;			   
			}
		  }
		          
       }

     }
     
  }

  vpbuf -> p = (double *) (vp + np);
  vpbuf -> np += np;

}


/****************************************************************************************

  Generate specific particle push functions for 2D and 3D, 1st to 4th order

****************************************************************************************/

#define __TEMPLATE__

// These macros will append _s1, _s2, etc to function names.
// These cannot be used to define the fortran function names because of recursive macro
// expansion issues.
#define ONAME(f, o) OJOIN(f, o)
#define OJOIN(f, o) f ## _s ## o


/********************************** Linear interpolation ********************************/

#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s1 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s1 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s1 )

// Interpolation order
#define ORDER 1
// Number of interpolation points ( ORDER + 1 )
#define NP 2
// Grid offset ( 1 + ORDER/2 )
#define OFFSET 1
// Current normalization ( see the wl* definitions in os-spec-current )
#define JNORM 1

#include __FILE__

/******************************** Quadratic interpolation *******************************/

#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s2 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s2 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s2 )

#define ORDER 2
#define NP 3
#define OFFSET 2
#define JNORM 2

#include __FILE__

/********************************** Cubic interpolation *********************************/

#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s3 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s3 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s3 )

#define ORDER 3
#define NP 4
#define OFFSET 2
#define JNORM 8

#include __FILE__

/********************************* Quartic interpolation ********************************/

#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s4 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s4 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s4 )

#define ORDER 4
#define NP 5
#define OFFSET 3
#define JNORM 48

#include __FILE__

/****************************************************************************************/

#else

/****************************************************************************************

  Template function definitions for 2D and 3D particle push

****************************************************************************************/

/***************   Generate Function names based on interpolation order  ****************/

#define SPLINE  ONAME( vsplined, ORDER )

#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

/********************************** 2D advance deposit **********************************/
 
extern void ADVANCE_DEPOSIT_2D
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  double *efield, double *bfield, int *emf_size, int *emf_offset,  
  double *j, int *j_size, int *j_offset, 
  double *dx, double *dt )
{
  int k, i, np, np_total;
  
  DECLARE_ALIGNED_16( double jnorm[2] );
  
  DECLARE_ALIGNED_16( int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_16( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_16( int diy[VEC_WIDTH] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_vpbuf2D vpbuf;

  // Simulation constants
  __m256d const vtem     = _mm256_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m256d const vdt_dx1 = _mm256_set1_pd( (*dt) / dx[0] );
  __m256d const vdt_dx2 = _mm256_set1_pd( (*dt) / dx[1] );
  
  // get pointers to position 0,0 of each field component
  int const deltaX = 3; // 3 field components
  int const deltaY = emf_size[0] * deltaX;

  double* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX + 
                              (emf_offset[1] - OFFSET) * deltaY;
  double* const e2 = e1 + 1;
  double* const e3 = e1 + 2;

  double* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX + 
                              (emf_offset[1] - OFFSET) * deltaY;
  double* const b2 = b1  + 1;
  double* const b3 = b1  + 2;
  
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / (2*JNORM);
  jnorm[1] = dx[1] / (*dt) / (2*JNORM);
  
  
  // jump to 1st particle
  x  += 2*(*i0-1);
  ix += 2*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 2 add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	
	for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH;  i ++ ) {
	  ix[ 2 * (np_total + i)     ] = 1;
	  ix[ 2 * (np_total + i) + 1 ] = 1;
	  x[ 2 * (np_total + i)     ] = 0.;
	  x[ 2 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i)     ] = 0.;
	  u[ 3 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i) + 2 ] = 0.;
	  q[ (np_total + i ) ] = 0.;
	}
    
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       printf("(*error*) Number of particles to push is still not a multiple of %d!\n", VEC_WIDTH);
       exit(1);
    }
    
  }
  
  if ( p_cache_size_2D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_2D is not a multiple of %d!\n", VEC_WIDTH);
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

	// Initialize virtual particle buffer
	vpbuf.np = 0;
	vpbuf.p = (double *) &(vpbuf.buf);
		
	// Push all particles in group
	for( i = 0; i < np; i+=VEC_WIDTH ) {
	  __m256d vx0, vx1, vy0, vy1;
	  __m128i vix0, vix1, viy0, viy1;
  
	  register __m256d ve1, ve2, ve3;
	  register __m256d vb1, vb2, vb3;
  
	  __m256d vu1, vu2, vu3;
	  
	  __m256d vv1, vv2, vv3;
	  
	  __m256d vq;
  
	  // Load particle positions 
	  _MM256_LOAD4v2_PD( vx0, vy0, x );
	  _MM_LOAD4v2_EPI32( vix0, viy0, ix );
  
	  // Interpolate fields
	  {
		 __m256d  vx0h, vy0h;
		 __m128i vidxi, vidxj, vidxih, vidxjh;
		 __m256d vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	 
	     // idxi = 3 * ix0
	     vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
	     // idxj = deltaY * iy0
	     vidxj = _mm_mullo_epi32( viy0, _mm_set1_epi32(deltaY) );
	     	     
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy);
	 
		 vget_hp_near( vx0, vidxi, vx0h, vidxih, 3 ); 
		 vget_hp_near( vy0, vidxj, vy0h, vidxjh, deltaY ); 
	 
		 SPLINE( vx0h, vwxh );
		 SPLINE( vy0h, vwyh );

		 // Interpolate E field
         {
			int k1, k2;
			
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx2[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx3[VEC_WIDTH] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32(  vidxi,  vidxj ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
              fld1[k1] = e1 + idx1[k1];
              fld2[k1] = e2 + idx2[k1];
              fld3[k1] = e3 + idx3[k1];
			}
			
			ve1 = _mm256_setzero_pd();
			ve2 = _mm256_setzero_pd();
			ve3 = _mm256_setzero_pd();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m256d f1line, f2line, f3line;
			  f1line = _mm256_setzero_pd();
			  f2line = _mm256_setzero_pd();
			  f3line = _mm256_setzero_pd();

			  __m256d f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm256_add_pd( f1line, _mm256_mul_pd( f1point[k1], vwxh[k1] ) );
				 f2line = _mm256_add_pd( f2line, _mm256_mul_pd( f2point[k1],  vwx[k1] ) );
				 f3line = _mm256_add_pd( f3line, _mm256_mul_pd( f3point[k1],  vwx[k1] ) );
			  }

			  ve1 = _mm256_add_pd( ve1, _mm256_mul_pd( f1line,  vwy[k2] ) );
			  ve2 = _mm256_add_pd( ve2, _mm256_mul_pd( f2line, vwyh[k2] ) );
			  ve3 = _mm256_add_pd( ve3, _mm256_mul_pd( f3line,  vwy[k2] ) );
			}
			        
         }		 

		 // Interpolate B field
         {
			int k1, k2;
			
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx2[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx3[VEC_WIDTH] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32( vidxih, vidxjh ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
              fld1[k1] = b1 + idx1[k1];
              fld2[k1] = b2 + idx2[k1];
              fld3[k1] = b3 + idx3[k1];
			}
			
			vb1 = _mm256_setzero_pd();
			vb2 = _mm256_setzero_pd();
			vb3 = _mm256_setzero_pd();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m256d f1line, f2line, f3line;
			  f1line = _mm256_setzero_pd();
			  f2line = _mm256_setzero_pd();
			  f3line = _mm256_setzero_pd();

			  __m256d f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm256_add_pd( f1line, _mm256_mul_pd( f1point[k1],  vwx[k1] ) );
				 f2line = _mm256_add_pd( f2line, _mm256_mul_pd( f2point[k1], vwxh[k1] ) );
				 f3line = _mm256_add_pd( f3line, _mm256_mul_pd( f3point[k1], vwxh[k1] ) );
			  }

			  vb1 = _mm256_add_pd( vb1, _mm256_mul_pd( f1line, vwyh[k2] ) );
			  vb2 = _mm256_add_pd( vb2, _mm256_mul_pd( f2line,  vwy[k2] ) );
			  vb3 = _mm256_add_pd( vb3, _mm256_mul_pd( f3line, vwyh[k2] ) );
			}
			        
         }		 
	  }

	   
	  // ---------- advance momenta
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, u, vu1, vu2, vu3 )
	  
	  // Store Momenta
	  _MM256_STORE4v3_PD(u, vu1, vu2, vu3) 
	     
	  // ---------- advance position and get velocities
	  {
		register __m256d vrg;
		register __m256d vtr1, vtr2;
		register __m128i vitr1, vitr2;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );

        vq = _mm256_load_pd( q );
	  
		// get velocities
		vv1 = _mm256_mul_pd( vu1, vrg );
		vv2 = _mm256_mul_pd( vu2, vrg );
		vv3 = _mm256_mul_pd( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		vx1 = _mm256_add_pd( vx0, _mm256_mul_pd( vv1, vdt_dx1 ) );
		vy1 = _mm256_add_pd( vy0, _mm256_mul_pd( vv2, vdt_dx2 ) );

        // Store virtual particles with positions still indexed to the
        // original cell
        STORE4VP2D( vpbuf.p, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );
	  
	    // Trim positions and store results
		vtr1 = vmntrim(vx1);
		vtr2 = vmntrim(vy1);

		vx1  = _mm256_sub_pd( vx1, vtr1 );
		vy1  = _mm256_sub_pd( vy1, vtr2 );
		_MM256_STORE4v2_PD( x, vx1, vy1 );
		
		vitr1 = _mm256_cvttpd_epi32( vtr1 ) ;
		vitr2 = _mm256_cvttpd_epi32( vtr2 ) ;
		_mm_store_si128( (__m128i *) dix, vitr1 );
		_mm_store_si128( (__m128i *) diy, vitr2 );
		
		vix1 = _mm_add_epi32( vix0, vitr1 ); 
		viy1 = _mm_add_epi32( viy0, vitr2 ); 
		
		_MM_STORE4v2_EPI32( ix, vix1, viy1 );

        // calculate crossings
        _mm_store_si128( (__m128i *) cross, 
                    _mm_add_epi32(_mm_andnot_si128(_mm_cmpeq_epi32(vix0, vix1),_mm_set1_epi32(1)), 
                                  _mm_andnot_si128(_mm_cmpeq_epi32(viy0, viy1),_mm_set1_epi32(2))));
	  }
	  
	  
	  // ---------- split trajectories for current deposition
	  
	  // The new split uses positions indexed to the original cell
      vsplit2D( &vpbuf, cross, dix, diy );

	  
	  // ---------- advance pointers
	  x  += 2 * VEC_WIDTH;
	  u  += 3 * VEC_WIDTH;  
	  ix += 2 * VEC_WIDTH;
	  q  += VEC_WIDTH;
  
	}
  
	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );
		
  }
    
}

extern void ADVANCE_DEPOSIT_2D_CYL
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  int *ilb2 , 
  double *efield, double *bfield, int *emf_size, int *emf_offset,  
  double *j, int *j_size, int *j_offset, 
  double *dx, double *dt )
{
  int k, i, np, np_total;
    
  DECLARE_ALIGNED_16( double jnorm[2] );
  
  DECLARE_ALIGNED_16( int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_16( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_16( int diy[VEC_WIDTH] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_vpbuf2D vpbuf;

  // Simulation Constants
  __m256d const vtem     = _mm256_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m256d const vdt_dx1 = _mm256_set1_pd( (*dt) / dx[0] );
  
  // get pointers to position 0,0 of each field component
  int const deltaX = 3; // 3 field components
  int const deltaY = emf_size[0] * deltaX;

  double* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX + 
                              (emf_offset[1] - OFFSET) * deltaY;
  double* const e2 = e1 + 1;
  double* const e3 = e1 + 2;

  double* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX + 
                              (emf_offset[1] - OFFSET) * deltaY;
  double* const b2 = b1  + 1;
  double* const b3 = b1  + 2;
  
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / (2*JNORM);
  jnorm[1] = dx[1] / (*dt) / (2*JNORM);
  
  __m256d const gshift2 = _mm256_set1_pd( (*ilb2 - 2 ) - 0.5 );
  __m256d const dr      = _mm256_set1_pd( dx[1] );
  __m256d const rdr     = _mm256_set1_pd( 1.0/dx[1] );
  __m256d const vdt     = _mm256_set1_pd( *dt );
  
  // jump to 1st particle
  x  += 2*(*i0-1);
  ix += 2*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 2 add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	
	for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH;  i ++ ) {
	  ix[ 2 * (np_total + i)     ] = 1;
	  ix[ 2 * (np_total + i) + 1 ] = 1;
	  x[ 2 * (np_total + i)     ] = 0.;
	  x[ 2 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i)     ] = 0.;
	  u[ 3 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i) + 2 ] = 0.;
	  q[ (np_total + i ) ] = 0.;
	}
    
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       printf("(*error*) Number of particles to push is still not a multiple of 2!\n");
       exit(1);
    }
    
  }
  
  if ( p_cache_size_2D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_2D is not a multiple of 2!\n");
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

	// Initialize virtual particle buffer
	vpbuf.np = 0;
	vpbuf.p = (double *) &(vpbuf.buf);
		
	// Push all particles in group
	for( i = 0; i < np; i+=VEC_WIDTH ) {
	  __m256d vx0, vx1, vy0, vy1;
	  __m128i vix0, vix1, viy0, viy1;
  
	  register __m256d ve1, ve2, ve3;
	  register __m256d vb1, vb2, vb3;
  
	  __m256d vu1, vu2, vu3;
	  
	  __m256d vv1, vv2, vv3;
	  
	  __m256d vq;
  
	  // Load particle positions 
	  _MM256_LOAD4v2_PD( vx0, vy0, x );
	  _MM_LOAD4v2_EPI32( vix0, viy0, ix );
  
	  // Interpolate fields
	  {
		 __m256d  vx0h, vy0h;
		 __m128i vidxi, vidxj, vidxih, vidxjh;
		 __m256d vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	 
	     // idxi = 3 * ix0
	     vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
	     // idxj = deltaY * iy0
	     vidxj = _mm_mullo_epi32( viy0, _mm_set1_epi32(deltaY) );
	     	     
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy);
	 
		 vget_hp_near( vx0, vidxi, vx0h, vidxih, 3 ); 
		 vget_hp_near( vy0, vidxj, vy0h, vidxjh, deltaY ); 
	 
		 SPLINE( vx0h, vwxh );
		 SPLINE( vy0h, vwyh );

		 // Interpolate E field
         {
			int k1, k2;
			
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx2[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx3[VEC_WIDTH] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32(  vidxi,  vidxj ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
              fld1[k1] = e1 + idx1[k1];
              fld2[k1] = e2 + idx2[k1];
              fld3[k1] = e3 + idx3[k1];
			}
			
			ve1 = _mm256_setzero_pd();
			ve2 = _mm256_setzero_pd();
			ve3 = _mm256_setzero_pd();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m256d f1line, f2line, f3line;
			  f1line = _mm256_setzero_pd();
			  f2line = _mm256_setzero_pd();
			  f3line = _mm256_setzero_pd();

			  __m256d f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm256_add_pd( f1line, _mm256_mul_pd( f1point[k1], vwxh[k1] ) );
				 f2line = _mm256_add_pd( f2line, _mm256_mul_pd( f2point[k1],  vwx[k1] ) );
				 f3line = _mm256_add_pd( f3line, _mm256_mul_pd( f3point[k1],  vwx[k1] ) );
			  }

			  ve1 = _mm256_add_pd( ve1, _mm256_mul_pd( f1line,  vwy[k2] ) );
			  ve2 = _mm256_add_pd( ve2, _mm256_mul_pd( f2line, vwyh[k2] ) );
			  ve3 = _mm256_add_pd( ve3, _mm256_mul_pd( f3line,  vwy[k2] ) );
			}
			        
         }		 

		 // Interpolate B field
         {
			int k1, k2;
			
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx2[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx3[VEC_WIDTH] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32( vidxih, vidxjh ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
              fld1[k1] = b1 + idx1[k1];
              fld2[k1] = b2 + idx2[k1];
              fld3[k1] = b3 + idx3[k1];
			}
			
			vb1 = _mm256_setzero_pd();
			vb2 = _mm256_setzero_pd();
			vb3 = _mm256_setzero_pd();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m256d f1line, f2line, f3line;
			  f1line = _mm256_setzero_pd();
			  f2line = _mm256_setzero_pd();
			  f3line = _mm256_setzero_pd();

			  __m256d f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm256_add_pd( f1line, _mm256_mul_pd( f1point[k1],  vwx[k1] ) );
				 f2line = _mm256_add_pd( f2line, _mm256_mul_pd( f2point[k1], vwxh[k1] ) );
				 f3line = _mm256_add_pd( f3line, _mm256_mul_pd( f3point[k1], vwxh[k1] ) );
			  }

			  vb1 = _mm256_add_pd( vb1, _mm256_mul_pd( f1line, vwyh[k2] ) );
			  vb2 = _mm256_add_pd( vb2, _mm256_mul_pd( f2line,  vwy[k2] ) );
			  vb3 = _mm256_add_pd( vb3, _mm256_mul_pd( f3line, vwyh[k2] ) );
			}
			        
         }		 
	  }

	   
	  // ---------- advance momenta
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, u, vu1, vu2, vu3 );
   
	  // ---------- advance position in cylindrical geometry
	  {
		register __m256d vrg;
		register __m256d vtr1, vtr2;
		register __m128i vitr1, vitr2;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );

        vq = _mm256_load_pd( q );
	  
		// get velocities
		vv1 = _mm256_mul_pd( vu1, vrg );
		vv2 = _mm256_mul_pd( vu2, vrg );
		vv3 = _mm256_mul_pd( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		
		// x1 normal push
		vx1 = _mm256_add_pd( vx0, _mm256_mul_pd( vv1, vdt_dx1 ) );
        
        // x2 - 3D like push
        // This further changes u2 and u3
        {  
		   
		   __m256d gix2, tmp;
		   __m256d r_old, r_new;
		   __m256d x2_new, x3_new; 
		   
		   // gshift2 = shift_ix2 - 0.5d0
   
		   // gix2   = viy0 + gshift2;          // global cell
		   gix2   = _mm256_add_pd( _mm256_cvtepi32_pd( viy0 ), gshift2 );
		   r_old  = _mm256_mul_pd( _mm256_add_pd( vy0, gix2 ), dr );

		   x2_new = _mm256_add_pd( r_old, _mm256_mul_pd( vv2, vdt ) );
		   x3_new =                    _mm256_mul_pd( vv3, vdt ) ;

		   r_new  = _mm256_sqrt_pd( _mm256_add_pd( _mm256_mul_pd( x2_new, x2_new ),
											 _mm256_mul_pd( x3_new, x3_new ) ) ); 
		   
		   // This is a protection against roundoff for cold plasmas
		   // vy1 = ( r_old == r_new ) ? vy0 : r_new * rdr - gix2);
		   vy1 = _mm256_blendv_pd( _mm256_sub_pd( _mm256_mul_pd( r_new, rdr ), gix2 ), vy0, 
		                           _mm256_cmp_pd( r_old, r_new, _CMP_EQ_OQ ) );
		   		   
		   // Correct p_r and p_\theta to conserve angular momentum
		   tmp  = _mm256_div_pd( _mm256_set1_pd(1.0), r_new );
		   vu2  = _mm256_mul_pd( _mm256_add_pd( _mm256_mul_pd( vu2, x2_new ), 
										  _mm256_mul_pd( vu3, x3_new ) ), 
										  tmp );
		   vu3  = _mm256_mul_pd( vu3 , _mm256_mul_pd( r_old, tmp ) );
   
		}        
        
		// Store Momenta
		_MM256_STORE4v3_PD(u, vu1, vu2, vu3) 
        
        // Store virtual particles with positions still indexed to the
        // original cell
        STORE4VP2D( vpbuf.p, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );
	  
	    // Trim positions and store results
		vtr1 = vmntrim(vx1);
		vtr2 = vmntrim(vy1);

		vx1  = _mm256_sub_pd( vx1, vtr1 );
		vy1  = _mm256_sub_pd( vy1, vtr2 );
		_MM256_STORE4v2_PD( x, vx1, vy1 );
		
		vitr1 = _mm256_cvttpd_epi32( vtr1 ) ;
		vitr2 = _mm256_cvttpd_epi32( vtr2 ) ;
		_mm_store_si128( (__m128i *) dix, vitr1 );
		_mm_store_si128( (__m128i *) diy, vitr2 );
		
		vix1 = _mm_add_epi32( vix0, vitr1 ); 
		viy1 = _mm_add_epi32( viy0, vitr2 ); 
		
		_MM_STORE4v2_EPI32( ix, vix1, viy1 );

        // calculate crossings
        _mm_store_si128( (__m128i *) cross, 
                    _mm_add_epi32(_mm_andnot_si128(_mm_cmpeq_epi32(vix0, vix1),_mm_set1_epi32(1)), 
                                  _mm_andnot_si128(_mm_cmpeq_epi32(viy0, viy1),_mm_set1_epi32(2))));
	  }
	  
	  
	  // ---------- split trajectories for current deposition
	  
	  // The new split uses positions indexed to the original cell
      vsplit2D( &vpbuf, cross, dix, diy );

	  
	  // ---------- advance pointers
	  x  += 2 * VEC_WIDTH;  
	  ix += 2 * VEC_WIDTH;
	  u  += 3 * VEC_WIDTH; 
	  q  += VEC_WIDTH;
  
	}
  
	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );
		
  }
    
}

/********************************** 3D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_3D
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  double *efield, double *bfield, int *emf_size, int *emf_offset,  
  double *j, int *j_size, int *j_offset, 
  double *dx, double *dt )
{
  int k, i, np, np_total;
    
  DECLARE_ALIGNED_16( double jnorm[3] );

  DECLARE_ALIGNED_16( int cross[ VEC_WIDTH ] );
  DECLARE_ALIGNED_16( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_16( int diy[VEC_WIDTH] );
  DECLARE_ALIGNED_16( int diz[VEC_WIDTH] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_vpbuf3D vpbuf;

  // Simulation constants
  __m256d const vtem    = _mm256_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m256d const vdt_dx1 = _mm256_set1_pd( (*dt) / dx[0] );
  __m256d const vdt_dx2 = _mm256_set1_pd( (*dt) / dx[1] );
  __m256d const vdt_dx3 = _mm256_set1_pd( (*dt) / dx[2] );
  
  // get pointers to position -1,-1,-1 of each field component
  int const deltaX = 3; // 3 field components
  int const deltaY = emf_size[0] * deltaX;
  int const deltaZ = emf_size[1] * deltaY;

  double* const e1 = efield + (emf_offset[0] - OFFSET)*deltaX +
                              (emf_offset[1] - OFFSET)*deltaY + 
                              (emf_offset[2] - OFFSET)*deltaZ;
  double* const e2 = e1 + 1;
  double* const e3 = e1 + 2;

  double* const b1 = bfield + (emf_offset[0] - OFFSET)*deltaX + 
                              (emf_offset[1] - OFFSET)*deltaY + 
                              (emf_offset[2] - OFFSET)*deltaZ;
  double* const b2 = b1  + 1;
  double* const b3 = b1  + 2;
    
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / (3*JNORM);
  jnorm[1] = dx[1] / (*dt) / (3*JNORM);
  jnorm[2] = dx[2] / (*dt) / (3*JNORM);
  
  // jump to 1st particle
  x  += 3*(*i0-1);
  ix += 3*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;
  
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	
	for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH; i ++ ) {
	  ix[ 3 * (np_total + i)     ] = 1;
	  ix[ 3 * (np_total + i) + 1 ] = 1;
	  ix[ 3 * (np_total + i) + 2 ] = 1;
	  x[ 3 * (np_total + i)     ] = 0.;
	  x[ 3 * (np_total + i) + 1 ] = 0.;
	  x[ 3 * (np_total + i) + 2 ] = 0.;
	  u[ 3 * (np_total + i)     ] = 0.;
	  u[ 3 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i) + 2 ] = 0.;
	  q[ (np_total + i ) ] = 0.;
	}
    
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       printf("(*error*) Number of particles to push is still not a multiple of 2!\n");
       exit(1);
    }
    
  }
  
  if ( p_cache_size_3D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_3D is not a multiple of 2!\n");
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_3D ) {
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_3D ) ? np_total - k : p_cache_size_3D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
	vpbuf.p = (double *) &(vpbuf.buf);
		
	// Push all particles in group
	for( i = 0; i < np; i += VEC_WIDTH ) {
	  __m256d vx0, vx1, vy0, vy1, vz0, vz1;
	  __m128i vix0, vix1, viy0, viy1, viz0, viz1;
	  __m256d vq;
  
	  register __m256d ve1, ve2, ve3;
	  register __m256d vb1, vb2, vb3;
  
	  __m256d vu1, vu2, vu3;
	    
	  // Load particle positions 
	  _MM256_LOAD4v3_PD( vx0, vy0, vz0, x );
	  _MM_LOAD4v3_EPI32( vix0, viy0, viz0, ix  );
  
	  // Interpolate fields
	  {
		 __m256d  vx0h, vy0h, vz0h;
		 __m128i vidxi, vidxj, vidxk, vidxih, vidxjh, vidxkh;
		 __m256d vwx[NP], vwy[NP], vwz[NP], vwxh[NP], vwyh[NP], vwzh[NP];

	     // idxi = 3 * ix0
	     vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
	     // idxj = deltaY * iy0
	     vidxj = _mm_mullo_epi32( viy0, _mm_set1_epi32(deltaY) );
	     // idxk = deltaZ * iz0
	     vidxk = _mm_mullo_epi32( viz0, _mm_set1_epi32(deltaZ) );
	 
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy );
		 SPLINE( vz0, vwz );
	 
		 vget_hp_near( vx0, vidxi, vx0h, vidxih, 3 ); 
		 vget_hp_near( vy0, vidxj, vy0h, vidxjh, deltaY ); 
		 vget_hp_near( vz0, vidxk, vz0h, vidxkh, deltaZ ); 
	 
		 SPLINE( vx0h, vwxh );
		 SPLINE( vy0h, vwyh );
		 SPLINE( vz0h, vwzh );
	 
         // Interpolate E field
         {
			int k1, k2, k3;
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx2[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx3[VEC_WIDTH] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, 
			                  _mm_add_epi32( _mm_add_epi32( vidxih, vidxj ), vidxk ) );
			_mm_store_si128( (__m128i *) idx2, 
			                  _mm_add_epi32( _mm_add_epi32( vidxi, vidxjh ), vidxk ) );
			_mm_store_si128( (__m128i *) idx3, 
			                  _mm_add_epi32( _mm_add_epi32( vidxi, vidxj ), vidxkh ) );

            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
              fld1[k1] = e1 + idx1[k1];
              fld2[k1] = e2 + idx2[k1];
              fld3[k1] = e3 + idx3[k1];
			}
			
			ve1 = _mm256_setzero_pd();
			ve2 = _mm256_setzero_pd();
			ve3 = _mm256_setzero_pd();
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register __m256d f1plane, f2plane, f3plane;
			   f1plane = _mm256_setzero_pd();
			   f2plane = _mm256_setzero_pd();
			   f3plane = _mm256_setzero_pd();
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register __m256d f1line, f2line, f3line;
				 f1line = _mm256_setzero_pd();
				 f2line = _mm256_setzero_pd();
				 f3line = _mm256_setzero_pd();

				 __m256d f1point[NP], f2point[NP], f3point[NP];
				 for ( k1 = 0; k1 < NP; k1++ ) {
					int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point[k1] = LOADFLD4( fld1, shift ); 
					f2point[k1] = LOADFLD4( fld2, shift ); 
					f3point[k1] = LOADFLD4( fld3, shift ); 
				 }
				 
				 for ( k1 = 0; k1 < NP; k1 ++ ) {
					f1line = _mm256_add_pd( f1line, _mm256_mul_pd( f1point[k1], vwxh[k1] ) );
					f2line = _mm256_add_pd( f2line, _mm256_mul_pd( f2point[k1],  vwx[k1] ) );
					f3line = _mm256_add_pd( f3line, _mm256_mul_pd( f3point[k1],  vwx[k1] ) );
				 }
				 
			     f1plane = _mm256_add_pd( f1plane, _mm256_mul_pd( f1line,  vwy[k2] ) );
			     f2plane = _mm256_add_pd( f2plane, _mm256_mul_pd( f2line, vwyh[k2] ) );
			     f3plane = _mm256_add_pd( f3plane, _mm256_mul_pd( f3line,  vwy[k2] ) );
			   }
			   
			   ve1 = _mm256_add_pd( ve1, _mm256_mul_pd( f1plane,  vwz[k3] ) );
			   ve2 = _mm256_add_pd( ve2, _mm256_mul_pd( f2plane,  vwz[k3] ) );
			   ve3 = _mm256_add_pd( ve3, _mm256_mul_pd( f3plane, vwzh[k3] ) );
			}
         
         }		 

         // Interpolate B field
	     {
			int k1, k2, k3;
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx2[VEC_WIDTH] );
			DECLARE_ALIGNED_16( int idx3[VEC_WIDTH] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, 
			                  _mm_add_epi32( _mm_add_epi32( vidxi, vidxjh ), vidxkh ) );
			_mm_store_si128( (__m128i *) idx2, 
			                  _mm_add_epi32( _mm_add_epi32( vidxih, vidxj ), vidxkh ) );
			_mm_store_si128( (__m128i *) idx3, 
			                  _mm_add_epi32( _mm_add_epi32( vidxih, vidxjh ), vidxk ) );

            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
              fld1[k1] = b1 + idx1[k1];
              fld2[k1] = b2 + idx2[k1];
              fld3[k1] = b3 + idx3[k1];
			}
			
			vb1 = _mm256_setzero_pd();
			vb2 = _mm256_setzero_pd();
			vb3 = _mm256_setzero_pd();
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register __m256d f1plane, f2plane, f3plane;
			   f1plane = _mm256_setzero_pd();
			   f2plane = _mm256_setzero_pd();
			   f3plane = _mm256_setzero_pd();
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register __m256d f1line, f2line, f3line;
				 f1line = _mm256_setzero_pd();
				 f2line = _mm256_setzero_pd();
				 f3line = _mm256_setzero_pd();

				 __m256d f1point[NP], f2point[NP], f3point[NP];
				 for ( k1 = 0; k1 < NP; k1++ ) {
					int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point[k1] = LOADFLD4( fld1, shift ); 
					f2point[k1] = LOADFLD4( fld2, shift ); 
					f3point[k1] = LOADFLD4( fld3, shift ); 
				 }
				 
				 for ( k1 = 0; k1 < NP; k1 ++ ) {
					f1line = _mm256_add_pd( f1line, _mm256_mul_pd( f1point[k1],  vwx[k1] ) );
					f2line = _mm256_add_pd( f2line, _mm256_mul_pd( f2point[k1], vwxh[k1] ) );
					f3line = _mm256_add_pd( f3line, _mm256_mul_pd( f3point[k1], vwxh[k1] ) );
				 }
				 
			     f1plane = _mm256_add_pd( f1plane, _mm256_mul_pd( f1line, vwyh[k2] ) );
			     f2plane = _mm256_add_pd( f2plane, _mm256_mul_pd( f2line,  vwy[k2] ) );
			     f3plane = _mm256_add_pd( f3plane, _mm256_mul_pd( f3line, vwyh[k2] ) );
			   }
			   
			   vb1 = _mm256_add_pd( vb1, _mm256_mul_pd( f1plane, vwzh[k3] ) );
			   vb2 = _mm256_add_pd( vb2, _mm256_mul_pd( f2plane, vwzh[k3] ) );
			   vb3 = _mm256_add_pd( vb3, _mm256_mul_pd( f3plane,  vwz[k3] ) );
			}
         
         }		 
		 
	  }	  
	  
	  // ---------- advance momenta
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, u, vu1, vu2, vu3 )
	  
	  // Store Momenta
	  _MM256_STORE4v3_PD(u, vu1, vu2, vu3) 

   
	  // ---------- advance position and get velocities
	  {
		register __m256d vrg;
		register __m256d vtrx, vtry, vtrz;
		register __m128i vitrx, vitry, vitrz;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );
	  
		// get velocities
		vu1 = _mm256_mul_pd( vu1, vrg );
		vu2 = _mm256_mul_pd( vu2, vrg );
		vu3 = _mm256_mul_pd( vu3, vrg );
			
		vx1 = _mm256_add_pd( vx0, _mm256_mul_pd( vu1, vdt_dx1 ) );
		vy1 = _mm256_add_pd( vy0, _mm256_mul_pd( vu2, vdt_dx2 ) );
		vz1 = _mm256_add_pd( vz0, _mm256_mul_pd( vu3, vdt_dx3 ) );
	  
		vtrx = vmntrim(vx1);
		vtry = vmntrim(vy1);
		vtrz = vmntrim(vz1);
		
		vq = _mm256_load_pd( q );
		
		STORE4VP3D( vpbuf.p, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix0, viy0, viz0 );

		vx1   = _mm256_sub_pd( vx1, vtrx );
		vy1   = _mm256_sub_pd( vy1, vtry );
		vz1   = _mm256_sub_pd( vz1, vtrz );
        _MM256_STORE4v3_PD( x, vx1, vy1, vz1 );

		vitrx = _mm256_cvttpd_epi32( vtrx );
		vitry = _mm256_cvttpd_epi32( vtry );
		vitrz = _mm256_cvttpd_epi32( vtrz );
		_mm_store_si128( (__m128i *) dix, vitrx );
		_mm_store_si128( (__m128i *) diy, vitry );
		_mm_store_si128( (__m128i *) diz, vitrz );
		
		vix1 = _mm_add_epi32( vix0, vitrx ); 
		viy1 = _mm_add_epi32( viy0, vitry ); 
		viz1 = _mm_add_epi32( viz0, vitrz ); 
		_MM_STORE4v3_EPI32( ix, vix1, viy1, viz1 );

		// use vitrx to save registers
		vitrx =                       _mm_andnot_si128(_mm_cmpeq_epi32(vix0, vix1),_mm_set1_epi32(1));
		vitrx = _mm_add_epi32( vitrx, _mm_andnot_si128(_mm_cmpeq_epi32(viy0, viy1),_mm_set1_epi32(2)));
		vitrx = _mm_add_epi32( vitrx, _mm_andnot_si128(_mm_cmpeq_epi32(viz0, viz1),_mm_set1_epi32(4)));
		_mm_store_si128( (__m128i *) cross, vitrx );                  

	  }
	  
	  
	  // ---------- split trajectories for current deposition

      vsplit3D( &vpbuf, cross, dix, diy, diz );

	  
	  // ---------- advance pointers
	  x  += 3*VEC_WIDTH;  
	  u  += 3*VEC_WIDTH; 
	  ix += 3*VEC_WIDTH;
	  q  += VEC_WIDTH;
  
	}
  
	// Deposit current from all virtual particles
	DEP_CURRENT_3D( j, j_size, j_offset, jnorm, &vpbuf );
		
  }
    
}


/****************************************************************************************
   Clear template definitions
****************************************************************************************/

#undef ORDER
#undef ADVANCE_DEPOSIT_2D
#undef ADVANCE_DEPOSIT_2D_CYL
#undef ADVANCE_DEPOSIT_3D
#undef OFFSET
#undef SPLINE
#undef NP
#undef JNORM
#undef DEP_CURRENT_2D
#undef DEP_CURRENT_3D

#endif
