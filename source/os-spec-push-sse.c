/*****************************************************************************************

Relativistic particle pusher, SSE optimized version (single precision)

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__


#include "os-spec-push-sse.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "fortran.h"

#include "vector-sse.h"
#include "splines-sse.h"
#include "os-spec-current-sse.h"


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

static inline __m128 vmntrim( const __m128 vx )
{
  register __m128 va, vb;

  va = _mm_cmplt_ps( vx, _mm_set_ps1( -0.5f ) );
  va = _mm_and_ps(   va, _mm_set_ps1( -1.0f ) );
  
  vb = _mm_cmpge_ps( vx, _mm_set_ps1( +0.5f ) );
  vb = _mm_and_ps(   vb, _mm_set_ps1( +1.0f ) );
  
  return _mm_add_ps( va, vb );
}


/*
// This version uses one less constant and saves a register
static inline __m128 vmntrim( const __m128 vx )
{
  __m128 const c1 = _mm_set_ps1( +1.0f );

  return _mm_sub_ps( _mm_and_ps( _mm_cmpge_ps( vx, _mm_set_ps1( +0.5f ) ), c1 ), 
                     _mm_and_ps( _mm_cmplt_ps( vx, _mm_set_ps1( -0.5f ) ), c1 ) );
}
*/

/*****************************************************************************************
LOADFLD4

Loads 4 field values corresponding to 4 particle positions into a vector
variable. Using _mm_set_ps yields the most efficient code in all scenarios

*****************************************************************************************/


#define LOADFLD4( fp, shift ) \
  _mm_set_ps( (fp[3])[shift], (fp[2])[shift], (fp[1])[shift], (fp[0])[shift] )

/*
// The above function should expand to something like this, but optimizes better

#define LOADFLD4( fp, shift ) \
  _mm_movelh_ps( _mm_unpacklo_ps( _mm_load_ss( &(fp[0])[shift] ), 
                                  _mm_load_ss( &(fp[1])[shift] ) ),
                 _mm_unpacklo_ps( _mm_load_ss( &(fp[2])[shift] ), 
                                  _mm_load_ss( &(fp[3])[shift] ) ) );
*/

/*****************************************************************************************
vget_hp_near

Gets the nearest half points of the 4 particles loaded as a vector for
interpolating scattered grids. 
This routine is for positions defined as distance to the nearest grid
points.
*****************************************************************************************/


#ifdef __SSE4_1__

// SSE 4.1 version - Use blendv function
#define vget_hp_near(dx, ix, dxh, ixh, delta ) 						    \
{ register __m128 cmp;													\
  cmp  = _mm_cmplt_ps( dx, _mm_setzero_ps() );	    					\
  ixh  = _mm_sub_epi32( ix, _mm_and_si128( _mm_castps_si128(cmp), 		\
                                           _mm_set1_epi32( delta )) );	\
  dxh =  _mm_add_ps( dx, _mm_blendv_ps( _mm_set1_ps( -0.5f ),           \
                                        _mm_set1_ps(  0.5f ), cmp ) );	\
}

#else

#define vget_hp_near(dx,ix,dxh,ixh,delta ) 							    \
{ register __m128 cmp;													\
  cmp  = _mm_cmplt_ps( dx, _mm_setzero_ps() );	    					\
  ixh  = _mm_sub_epi32( ix, _mm_and_si128( _mm_castps_si128(cmp), 		\
                                           _mm_set1_epi32( delta )) );	\
  dxh = _mm_sub_ps( dx, _mm_set_ps1( 0.5 ) );							\
  dxh = _mm_add_ps( dxh, _mm_and_ps( cmp, _mm_set_ps1( 1.0 )) );		\
}

#endif

/*****************************************************************************************
vdudt_boris

Use a Boris push to advance particle velocities using interpolated
fields.
*****************************************************************************************/

#if 0
inline void vdudt_boris( __m128 vtem, __m128 ve1, __m128 ve2, __m128 ve3,
                         __m128 vb1, __m128 vb2, __m128 vb3,
                         float *u, 
                         __m128 *vu1, __m128 *vu2, __m128 *vu3 )
{
   register __m128 vut1, vut2, vut3;
   register __m128 tmp1, tmp2, tmp3;
   register __m128 vgamma_tem, votsq;
   
   // Perform first half of electric field acceleration.
   ve1 = _mm_mul_ps( ve1, vtem );
   ve2 = _mm_mul_ps( ve2, vtem );
   ve3 = _mm_mul_ps( ve3, vtem );
   
   _MM_LOAD4v3_PS( *vu1, *vu2, *vu3, u )
   
   vut1 = _mm_add_ps( *vu1, ve1 );
   vut2 = _mm_add_ps( *vu2, ve2 );
   vut3 = _mm_add_ps( *vu3, ve3 );

   // Perform first half of the rotation and store in u.
   tmp1 = _mm_mul_ps( vut1, vut1 );
   tmp2 = _mm_mul_ps( vut2, vut2 );
   tmp3 = _mm_mul_ps( vut3, vut3 );
	
   vgamma_tem = _mm_add_ps( tmp1, tmp2 );
   vgamma_tem = _mm_add_ps( vgamma_tem, tmp3 );
   vgamma_tem = _mm_add_ps( vgamma_tem, _mm_set1_ps( 1.0 ) );
   
   vgamma_tem = _mm_sqrt_ps( vgamma_tem );
   vgamma_tem = _mm_div_ps( vtem, vgamma_tem );
	
   vb1 = _mm_mul_ps( vb1, vgamma_tem );
   vb2 = _mm_mul_ps( vb2, vgamma_tem );
   vb3 = _mm_mul_ps( vb3, vgamma_tem );
   
   *vu1 = _mm_add_ps( vut1, _mm_sub_ps( _mm_mul_ps( vut2, vb3 ), _mm_mul_ps( vut3, vb2 ) ) );
   *vu2 = _mm_add_ps( vut2, _mm_sub_ps( _mm_mul_ps( vut3, vb1 ), _mm_mul_ps( vut1, vb3 ) ) );
   *vu3 = _mm_add_ps( vut3, _mm_sub_ps( _mm_mul_ps( vut1, vb2 ), _mm_mul_ps( vut2, vb1 ) ) );

   // Perform second half of the rotation.
   tmp1 = _mm_mul_ps( vb1, vb1 );
   tmp2 = _mm_mul_ps( vb2, vb2 );
   tmp3 = _mm_mul_ps( vb3, vb3 );
   
   tmp1 = _mm_add_ps( tmp1, tmp2 );
   tmp1 = _mm_add_ps( tmp1, tmp3 );
   tmp1 = _mm_add_ps( tmp1, _mm_set1_ps(1.0f));
   
   votsq = _mm_div_ps( _mm_set1_ps(2.0f), tmp1 ); 

   vb1 = _mm_mul_ps( vb1, votsq );
   vb2 = _mm_mul_ps( vb2, votsq );
   vb3 = _mm_mul_ps( vb3, votsq );
   
   vut1 = _mm_add_ps( vut1, _mm_sub_ps( _mm_mul_ps( *vu2, vb3 ), _mm_mul_ps( *vu3, vb2 ) ) );
   vut2 = _mm_add_ps( vut2, _mm_sub_ps( _mm_mul_ps( *vu3, vb1 ), _mm_mul_ps( *vu1, vb3 ) ) );
   vut3 = _mm_add_ps( vut3, _mm_sub_ps( _mm_mul_ps( *vu1, vb2 ), _mm_mul_ps( *vu2, vb1 ) ) );
   
   // Perform second half of electric field acceleration.
   *vu1 = _mm_add_ps( vut1, ve1 );
   *vu2 = _mm_add_ps( vut2, ve2 );
   *vu3 = _mm_add_ps( vut3, ve3 );
   
   // Store results
   _MM_STORE4v3_PS(u, *vu1, *vu2, *vu3) 

}
#endif

#define VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, u, vu1, vu2, vu3 )                   \
{                                                                                             \
   __m128 vut1, vut2, vut3;\
   __m128 tmp1, tmp2, tmp3;\
   __m128 vgamma_tem, votsq;\
   \
   /* Perform first half of electric field acceleration.*/\
   ve1 = _mm_mul_ps( ve1, vtem );\
   ve2 = _mm_mul_ps( ve2, vtem );\
   ve3 = _mm_mul_ps( ve3, vtem );\
   \
   _MM_LOAD4v3_PS( vu1, vu2, vu3, u )\
   \
   vut1 = _mm_add_ps( vu1, ve1 );\
   vut2 = _mm_add_ps( vu2, ve2 );\
   vut3 = _mm_add_ps( vu3, ve3 );\
   \
   /* Perform first half of the rotation and store in u. */\
   tmp1 = _mm_mul_ps( vut1, vut1 );\
   tmp2 = _mm_mul_ps( vut2, vut2 );\
   tmp3 = _mm_mul_ps( vut3, vut3 );\
   \
   vgamma_tem = _mm_add_ps( tmp1, tmp2 );\
   vgamma_tem = _mm_add_ps( vgamma_tem, tmp3 );\
   vgamma_tem = _mm_add_ps( vgamma_tem, _mm_set1_ps( 1.0f ) );\
   \
   vgamma_tem = _mm_sqrt_ps( vgamma_tem );\
   vgamma_tem = _mm_div_ps( vtem, vgamma_tem );\
   \
   vb1 = _mm_mul_ps( vb1, vgamma_tem );\
   vb2 = _mm_mul_ps( vb2, vgamma_tem );\
   vb3 = _mm_mul_ps( vb3, vgamma_tem );\
   \
   vu1 = _mm_add_ps( vut1, _mm_sub_ps( _mm_mul_ps( vut2, vb3 ), _mm_mul_ps( vut3, vb2 ) ) );\
   vu2 = _mm_add_ps( vut2, _mm_sub_ps( _mm_mul_ps( vut3, vb1 ), _mm_mul_ps( vut1, vb3 ) ) );\
   vu3 = _mm_add_ps( vut3, _mm_sub_ps( _mm_mul_ps( vut1, vb2 ), _mm_mul_ps( vut2, vb1 ) ) );\
   \
   /* Perform second half of the rotation. */\
   tmp1 = _mm_mul_ps( vb1, vb1 );\
   tmp2 = _mm_mul_ps( vb2, vb2 );\
   tmp3 = _mm_mul_ps( vb3, vb3 );\
   \
   tmp1 = _mm_add_ps( tmp1, tmp2 );\
   tmp1 = _mm_add_ps( tmp1, tmp3 );\
   tmp1 = _mm_add_ps( tmp1, _mm_set1_ps(1.0f));\
   \
   votsq = _mm_div_ps( _mm_set1_ps(2.0f), tmp1 );\
   \
   vb1 = _mm_mul_ps( vb1, votsq );\
   vb2 = _mm_mul_ps( vb2, votsq );\
   vb3 = _mm_mul_ps( vb3, votsq );\
   \
   vut1 = _mm_add_ps( vut1, _mm_sub_ps( _mm_mul_ps( vu2, vb3 ), _mm_mul_ps( vu3, vb2 ) ) );\
   vut2 = _mm_add_ps( vut2, _mm_sub_ps( _mm_mul_ps( vu3, vb1 ), _mm_mul_ps( vu1, vb3 ) ) );\
   vut3 = _mm_add_ps( vut3, _mm_sub_ps( _mm_mul_ps( vu1, vb2 ), _mm_mul_ps( vu2, vb1 ) ) );\
   \
   /* Perform second half of electric field acceleration.*/\
   vu1 = _mm_add_ps( vut1, ve1 );\
   vu2 = _mm_add_ps( vut2, ve2 );\
   vu3 = _mm_add_ps( vut3, ve3 );\
\
}

/*****************************************************************************************
VRGAMMA

Get 1 / Lorentz gamma. Using a macro is beneficial because this is 
called from a function that is already inlined.

Note: using rsqrt_ps did not yield any measurable performance improvement
*****************************************************************************************/

#define VRGAMMA( vrg, vu1, vu2, vu3 ) {  \
  														    \
  (vrg) = _mm_add_ps( _mm_add_ps(_mm_mul_ps((vu1),(vu1)),   \
                                 _mm_mul_ps((vu2),(vu2))),  \
                                 _mm_mul_ps((vu3),(vu3)) ); \
  (vrg) = _mm_add_ps( (vrg), _mm_set1_ps(1.0f) );		    \
  (vrg) = _mm_sqrt_ps( (vrg) );							    \
  (vrg) = _mm_div_ps( _mm_set1_ps(1.0f), (vrg) );		    \
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

Splits 4 particle trajectories for current deposition and stores segment on a special virtual
particle buffer which is optimized for the vector current deposition routines
*****************************************************************************************/

void vsplit2D( t_vpbuf2D* const vpbuf, const int cross[], const int dix[], const int diy[] )
{
 
  int k, np;
  t_vp2D *vp;

  np = 4;
  vp = (t_vp2D *) vpbuf -> p;
  
  for( k = 0 ; k < 4; k++ ) {
     float delta, xint, yint, vzint, xint2, yint2, vzint2;
          
     switch ( cross[k] ) {
       case (0) :  // no split

         // no action required
		 break;
         
       case (1) :  //  x cross only

		 xint  = 0.5f * dix[k];
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

		 yint  = 0.5f * ( diy[k] );
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
		 xint  = 0.5f * dix[k];
		 delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);
		 yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta;

         if ( ( yint >= -0.5f ) && (yint < 0.5f) ) {

           // no y cross on 1st vp
           vzint = vp[k].vz * delta;

           // y split 2nd vp
           vzint2 = vp[k].vz * (1-delta);
           yint2 = 0.5f * diy[k];
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
		   yint2 = 0.5f * diy[k];
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

  vpbuf -> p = (float *) (vp + np);
  vpbuf -> np += np;
}


/*****************************************************************************************
vsplit3D

Splits 4 particle trajectories for current deposition and stores segment on a special virtual
particle buffer which is optimized for the vector current deposition routines
*****************************************************************************************/

void vsplit3D( t_vpbuf3D* const vpbuf, const int cross[], 
               const int dix[], const int diy[], const int diz[] )
{

  int k, np;
  t_vp3D *vp;

  np = 4;
  vp = (t_vp3D *) vpbuf -> p;

  for( k = 0; k < 4; k++ ) {
     float xint, yint, zint;
     float xint2, yint2, zint2;
     float delta;
        
     switch ( cross[k] ) {
       case (0) :  // no split

         // no action required
		 break;
         
       case (1) :  //  x cross only
         
		 xint  = 0.5f * dix[k];
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

		 yint  = 0.5f * diy[k];
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
		 xint  = 0.5f * dix[k];
		 delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);

		 // note that this is different from case(1) because of change of cell in y direction
		 yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta;
		 zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta;

         if ( ( yint >= -0.5f ) && (yint < 0.5f) ) {
                      
           // y split 2nd vp
           yint2 = 0.5f * diy[k];
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
		   yint2 = 0.5f * (diy[k]);
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

		 zint  = 0.5f * diz[k];
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
		 xint  = 0.5f*( dix[k] );
		 delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);

		 // note that this is different from case(1) because of change of cell in z direction
		 yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta; // no y cross
		 zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta; 

         if ( ( zint >= -0.5f ) && (zint < 0.5f) ) {
                      
           // z split 2nd vp
           zint2 = 0.5f * (diz[k]);
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
		   zint2 = 0.5f * (diz[k]);
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
		 yint  = 0.5f*( diy[k] );
		 delta = ( yint - vp[k].y0 ) / (vp[k].y1 - vp[k].y0);

		 xint  = vp[k].x0 + (vp[k].x1 - vp[k].x0) * delta; // no x cross
		 zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta; 

         if ( ( zint >= -0.5f ) && (zint < 0.5f) ) {
                      
           // z split 2nd vp
           zint2 = 0.5f * diz[k];
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
		   zint2 = 0.5f * diz[k];
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
		  xint  = 0.5f * dix[k];
		  delta = ( xint - vp[k].x0 ) / (vp[k].x1 - vp[k].x0);
		  yint  = vp[k].y0 + (vp[k].y1 - vp[k].y0) * delta;
		  zint  = vp[k].z0 + (vp[k].z1 - vp[k].z0) * delta; 

         if ( ( yint >= -0.5f ) && (yint < 0.5f) ) {
                      
           // y split 2nd vp
           yint2 = 0.5f * (diy[k]);
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
		   yint2 = 0.5f * (diy[k]);
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
          
          /*int vpidx[3] = {k, np-2, np-1};*/
          int vpidx[3]; vpidx[0] = k; vpidx[1] = np-2; vpidx[2] = np-1; 
          
          int i, j, idx;

		  zint = 0.5f * (diz[k]);
		  for ( i = 0; i < 3; i ++ ) {
			idx = vpidx[i];
			if ( (vp[idx].z1 < -0.5f) || (vp[idx].z1 >= +0.5f) ) {
			   
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

  vpbuf -> p = (float *) (vp + np);
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

#define ADVANCE_DEPOSIT_2D FNAME( vadvance_deposit_2d_s2 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s2 )
#define ADVANCE_DEPOSIT_3D FNAME( vadvance_deposit_3d_s2 )

#define ORDER 2
#define NP 3
#define OFFSET 2
#define JNORM 2

#include __FILE__

/********************************** Cubic interpolation *********************************/

#define ADVANCE_DEPOSIT_2D FNAME( vadvance_deposit_2d_s3 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s3 )
#define ADVANCE_DEPOSIT_3D FNAME( vadvance_deposit_3d_s3 )

#define ORDER 3
#define NP 4
#define OFFSET 2
#define JNORM 8

#include __FILE__

/********************************* Quartic interpolation ********************************/

#define ADVANCE_DEPOSIT_2D FNAME( vadvance_deposit_2d_s4 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s4 )
#define ADVANCE_DEPOSIT_3D FNAME( vadvance_deposit_3d_s4 )

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

#define SPLINE  ONAME( vspline, ORDER )

#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

/********************************** 2D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_2D
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  float *efield, float *bfield, int *emf_size, int *emf_offset,  
  float *j, int *j_size, int *j_offset, 
  double *dx, double *dt )
{
  int k, i, np, np_total;
  
  DECLARE_ALIGNED_16( float jnorm[2] );
  
  DECLARE_ALIGNED_16( int cross[4] );
  DECLARE_ALIGNED_16( int dix[4] );
  DECLARE_ALIGNED_16( int diy[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_vpbuf2D vpbuf;

  // Simulation constants
  __m128 const vtem     = _mm_set_ps1( 0.5f * (*dt) / (*rqm) );
  __m128 const vdt_dx1 = _mm_set_ps1( (*dt) / dx[0] );
  __m128 const vdt_dx2 = _mm_set_ps1( (*dt) / dx[1] );

  // get pointers to position 0,0 of each field component
  int const deltaX = 3; // 3 field components
  int const deltaY = emf_size[0] * deltaX;

  float* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX + 
                             (emf_offset[1] - OFFSET) * deltaY;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX + 
                             (emf_offset[1] - OFFSET) * deltaY;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;
  
  
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

  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np_total % 4 ) {
	
	for( i = 0; i < 4 - np_total%4; i ++ ) {
	  ix[ 2 * (np_total + i)     ] = 1;
	  ix[ 2 * (np_total + i) + 1 ] = 1;
	  x[ 2 * (np_total + i)     ] = 0.;
	  x[ 2 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i)     ] = 0.;
	  u[ 3 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i) + 2 ] = 0.;
	  q[ (np_total + i ) ] = 0.;
	}
    
    np_total += ( 4 - np_total%4 );
    if ( np_total % 4 ) {
       printf("(*error*) Number of particles to push is still not a multiple of 4!\n");
       exit(1);
    }
    
  }
  
  if ( p_cache_size_2D % 4 ) {
       printf("(*error*) p_cache_size_2D is not a multiple of 4!\n");
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
	vpbuf.p = (float *) &(vpbuf.buf);
		
	// Push all particles in group
	for( i = 0; i < np; i+=4 ) {
	  __m128 vx0, vx1, vy0, vy1;
	  __m128i vix0, vix1, viy0, viy1;
  
	  register __m128 ve1, ve2, ve3;
	  register __m128 vb1, vb2, vb3;
  
	  __m128 vu1, vu2, vu3;
  
	  // Load particle positions 
	  _MM_LOAD4v2_PS( vx0, vy0, x );
	  _MM_LOAD4v2_EPI32( vix0, viy0, ix  );
  
	  // Interpolate fields
	  {
		 __m128  vx0h, vy0h;
		 __m128i vidxi, vidxj, vidxih, vidxjh;
		 __m128 vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	 
	     // idxi = 3 * ix0
	     vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
	     // idxj = deltaY * iy0
	     vidxj = vsmul( viy0, deltaY );
	     	     
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy);
	 
		 vget_hp_near( vx0, vidxi, vx0h, vidxih, 3 ); 
		 vget_hp_near( vy0, vidxj, vy0h, vidxjh, deltaY ); 
	 
		 SPLINE( vx0h, vwxh );
		 SPLINE( vy0h, vwyh );

		 // Interpolate E field
         {
			int k1, k2;
			
			float *fld1[4], *fld2[4], *fld3[4];
			DECLARE_ALIGNED_16( int idx1[4] );
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32(  vidxi,  vidxj ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 4; k1++) {
              fld1[k1] = e1 + idx1[k1];
              fld2[k1] = e2 + idx2[k1];
              fld3[k1] = e3 + idx3[k1];
			}
			
			ve1 = _mm_setzero_ps();
			ve2 = _mm_setzero_ps();
			ve3 = _mm_setzero_ps();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m128 f1line, f2line, f3line;
			  f1line = _mm_setzero_ps();
			  f2line = _mm_setzero_ps();
			  f3line = _mm_setzero_ps();

			  __m128 f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm_add_ps( f1line, _mm_mul_ps( f1point[k1], vwxh[k1] ) );
				 f2line = _mm_add_ps( f2line, _mm_mul_ps( f2point[k1],  vwx[k1] ) );
				 f3line = _mm_add_ps( f3line, _mm_mul_ps( f3point[k1],  vwx[k1] ) );
			  }

			  ve1 = _mm_add_ps( ve1, _mm_mul_ps( f1line,  vwy[k2] ) );
			  ve2 = _mm_add_ps( ve2, _mm_mul_ps( f2line, vwyh[k2] ) );
			  ve3 = _mm_add_ps( ve3, _mm_mul_ps( f3line,  vwy[k2] ) );
			}
			        
         }		 

		 // Interpolate B field
         {
			int k1, k2;
			
			float *fld1[4], *fld2[4], *fld3[4];
			DECLARE_ALIGNED_16( int idx1[4] );
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32( vidxih, vidxjh ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 4; k1++) {
              fld1[k1] = b1 + idx1[k1];
              fld2[k1] = b2 + idx2[k1];
              fld3[k1] = b3 + idx3[k1];
			}
			
			vb1 = _mm_setzero_ps();
			vb2 = _mm_setzero_ps();
			vb3 = _mm_setzero_ps();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m128 f1line, f2line, f3line;
			  f1line = _mm_setzero_ps();
			  f2line = _mm_setzero_ps();
			  f3line = _mm_setzero_ps();

			  __m128 f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm_add_ps( f1line, _mm_mul_ps( f1point[k1],  vwx[k1] ) );
				 f2line = _mm_add_ps( f2line, _mm_mul_ps( f2point[k1], vwxh[k1] ) );
				 f3line = _mm_add_ps( f3line, _mm_mul_ps( f3point[k1], vwxh[k1] ) );
			  }

			  vb1 = _mm_add_ps( vb1, _mm_mul_ps( f1line, vwyh[k2] ) );
			  vb2 = _mm_add_ps( vb2, _mm_mul_ps( f2line,  vwy[k2] ) );
			  vb3 = _mm_add_ps( vb3, _mm_mul_ps( f3line, vwyh[k2] ) );
			}
			        
         }		 
	  }
	   
	  // ---------- advance momenta
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, u, vu1, vu2, vu3 );

      // Store Results
      _MM_STORE4v3_PS(u, vu1, vu2, vu3); 

	  // advance u
	  u  += 12; 
   
	  // ---------- advance position and get velocities
	  {
		register __m128 vrg;
		register __m128 vtr1, vtr2;
		register __m128i vitr1, vitr2;
	    register __m128 vv1, vv2, vv3;
	    register __m128 vq;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );

        vq = _mm_load_ps( q );
	  
		// get velocities
		vv1 = _mm_mul_ps( vu1, vrg );
		vv2 = _mm_mul_ps( vu2, vrg );
		vv3 = _mm_mul_ps( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		vx1 = _mm_add_ps( vx0, _mm_mul_ps( vv1, vdt_dx1 ) );
		vy1 = _mm_add_ps( vy0, _mm_mul_ps( vv2, vdt_dx2 ) );

        // Store virtual particles with positions still indexed to the
        // original cell
        STORE4VP2D( vpbuf.p, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );
	  
	    // Trim positions and store results
		vtr1 = vmntrim(vx1);
		vtr2 = vmntrim(vy1);

		vx1  = _mm_sub_ps( vx1, vtr1 );
		vy1  = _mm_sub_ps( vy1, vtr2 );
		_MM_STORE4v2_PS( x, vx1, vy1 );
		
		vitr1 = _mm_cvttps_epi32( vtr1 ) ;
		vitr2 = _mm_cvttps_epi32( vtr2 ) ;
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
	  x  += 8;  
	  ix += 8;
	  q  += 4;
  
	}
  
	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );
		
  }
    
}

extern void ADVANCE_DEPOSIT_2D_CYL
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  int *ilb2 , 
  float *efield, float *bfield, int *emf_size, int *emf_offset,  
  float *j, int *j_size, int *j_offset, 
  double *dx, double *dt )
{
  int k, i, np, np_total;
  
  DECLARE_ALIGNED_16( float jnorm[2] );
  
  DECLARE_ALIGNED_16( int cross[4] );
  DECLARE_ALIGNED_16( int dix[4] );
  DECLARE_ALIGNED_16( int diy[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_vpbuf2D vpbuf;
  
  // Simulation Constants
  __m128 const vtem     = _mm_set_ps1( 0.5 * (*dt) / (*rqm) );
  __m128 const vdt_dx1 = _mm_set_ps1( (*dt) / dx[0] );
  
  // get pointers to position 0,0 of each field component
  int const deltaX = 3; // 3 field components
  int const deltaY = emf_size[0] * deltaX;

  float* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX + 
                             (emf_offset[1] - OFFSET) * deltaY;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX + 
                             (emf_offset[1] - OFFSET) * deltaY;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;
   
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / (2*JNORM);
  jnorm[1] = dx[1] / (*dt) / (2*JNORM);
  
  __m128d const gshift2 = _mm_set1_pd( (*ilb2 - 2 ) - 0.5 );
  __m128d const dr      = _mm_set1_pd( dx[1] );
  __m128d const rdr     = _mm_set1_pd( 1.0/dx[1] );
  __m128d const vdt     = _mm_set1_pd( *dt );
  
  // jump to 1st particle
  x  += 2*(*i0-1);
  ix += 2*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np_total % 4 ) {
	
	for( i = 0; i < 4 - np_total%4; i ++ ) {
	  ix[ 2 * (np_total + i)     ] = 1;
	  ix[ 2 * (np_total + i) + 1 ] = 1;
	  x[ 2 * (np_total + i)     ] = 0.;
	  x[ 2 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i)     ] = 0.;
	  u[ 3 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i) + 2 ] = 0.;
	  q[ (np_total + i ) ] = 0.;
	}
    
    np_total += ( 4 - np_total%4 );
    if ( np_total % 4 ) {
       printf("(*error*) Number of particles to push is still not a multiple of 4!\n");
       exit(1);
    }
    
  }
  
  if ( p_cache_size_2D % 4 ) {
       printf("(*error*) p_cache_size_2D is not a multiple of 4!\n");
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
	vpbuf.p = (float *) &(vpbuf.buf);
		
	// Push all particles in group
	for( i = 0; i < np; i+=4 ) {
	  __m128 vx0, vx1, vy0, vy1;
	  __m128i vix0, vix1, viy0, viy1;
  
	  register __m128 ve1, ve2, ve3;
	  register __m128 vb1, vb2, vb3;
  
	  __m128 vu1, vu2, vu3;
  
	  // Load particle positions 
	  _MM_LOAD4v2_PS( vx0, vy0, x );
	  _MM_LOAD4v2_EPI32( vix0, viy0, ix  );
  
	  // Interpolate fields
	  {
		 __m128  vx0h, vy0h;
		 __m128i vidxi, vidxj, vidxih, vidxjh;
		 __m128 vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	 
	     // idxi = 3 * ix0
	     vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
	     // idxj = deltaY * iy0
	     vidxj = vsmul( viy0, deltaY );
	     	     
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy);
	 
		 vget_hp_near( vx0, vidxi, vx0h, vidxih, 3 ); 
		 vget_hp_near( vy0, vidxj, vy0h, vidxjh, deltaY ); 
	 
		 SPLINE( vx0h, vwxh );
		 SPLINE( vy0h, vwyh );

		 // Interpolate E field
         {
			int k1, k2;
			
			float *fld1[4], *fld2[4], *fld3[4];
			DECLARE_ALIGNED_16( int idx1[4] );
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32(  vidxi,  vidxj ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 4; k1++) {
              fld1[k1] = e1 + idx1[k1];
              fld2[k1] = e2 + idx2[k1];
              fld3[k1] = e3 + idx3[k1];
			}
			
			ve1 = _mm_setzero_ps();
			ve2 = _mm_setzero_ps();
			ve3 = _mm_setzero_ps();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m128 f1line, f2line, f3line;
			  f1line = _mm_setzero_ps();
			  f2line = _mm_setzero_ps();
			  f3line = _mm_setzero_ps();

			  __m128 f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm_add_ps( f1line, _mm_mul_ps( f1point[k1], vwxh[k1] ) );
				 f2line = _mm_add_ps( f2line, _mm_mul_ps( f2point[k1],  vwx[k1] ) );
				 f3line = _mm_add_ps( f3line, _mm_mul_ps( f3point[k1],  vwx[k1] ) );
			  }

			  ve1 = _mm_add_ps( ve1, _mm_mul_ps( f1line,  vwy[k2] ) );
			  ve2 = _mm_add_ps( ve2, _mm_mul_ps( f2line, vwyh[k2] ) );
			  ve3 = _mm_add_ps( ve3, _mm_mul_ps( f3line,  vwy[k2] ) );
			}
			        
         }		 

		 // Interpolate B field
         {
			int k1, k2;
			
			float *fld1[4], *fld2[4], *fld3[4];
			DECLARE_ALIGNED_16( int idx1[4] );
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32( vidxih, vidxjh ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 4; k1++) {
              fld1[k1] = b1 + idx1[k1];
              fld2[k1] = b2 + idx2[k1];
              fld3[k1] = b3 + idx3[k1];
			}
			
			vb1 = _mm_setzero_ps();
			vb2 = _mm_setzero_ps();
			vb3 = _mm_setzero_ps();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m128 f1line, f2line, f3line;
			  f1line = _mm_setzero_ps();
			  f2line = _mm_setzero_ps();
			  f3line = _mm_setzero_ps();

			  __m128 f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm_add_ps( f1line, _mm_mul_ps( f1point[k1],  vwx[k1] ) );
				 f2line = _mm_add_ps( f2line, _mm_mul_ps( f2point[k1], vwxh[k1] ) );
				 f3line = _mm_add_ps( f3line, _mm_mul_ps( f3point[k1], vwxh[k1] ) );
			  }

			  vb1 = _mm_add_ps( vb1, _mm_mul_ps( f1line, vwyh[k2] ) );
			  vb2 = _mm_add_ps( vb2, _mm_mul_ps( f2line,  vwy[k2] ) );
			  vb3 = _mm_add_ps( vb3, _mm_mul_ps( f3line, vwyh[k2] ) );
			}
			        
         }		 
	  }
	   
	  // ---------- advance momenta
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, u, vu1, vu2, vu3 );
      
   
	  // ---------- advance position in cylindrical geometry
	  {
		__m128 vrg;
		__m128 vtr1, vtr2;
		__m128i vitr1, vitr2;
	    __m128 vv1, vv2, vv3;
	    __m128 vq;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );

        vq = _mm_load_ps( q );
	  
		// get velocities
		vv1 = _mm_mul_ps( vu1, vrg );
		vv2 = _mm_mul_ps( vu2, vrg );
		vv3 = _mm_mul_ps( vu3, vrg ); 
			
		// Push positions
		vx1 = _mm_add_ps( vx0, _mm_mul_ps( vv1, vdt_dx1 ) );

        {  // This section needs to be performed in double precision
		   
		   __m128d gix2a, gix2b;
		   __m128d vy0a, vy0b, vy1a, vy1b;
		   __m128d r_olda, r_oldb, r_newa, r_newb;
		   __m128d x2_newa, x2_newb, x3_newa, x3_newb; 
		   
		   __m128d vv2a, vv2b, vv3a, vv3b;
		   __m128d vu2a, vu2b, vu3a, vu3b;
		   
		   __m128d tmpa, tmpb;
		   
		   // gshift2 = shift_ix2 - 0.5d0
   
		   // gix2   = viy0 + gshift2;          // global cell
		   gix2a = _mm_add_pd( _mm_cvtepi32_pd( viy0 ), gshift2 );
		   gix2b = _mm_add_pd( _mm_cvtepi32_pd( _mm_unpackhi_epi64( viy0, viy0 ) ), gshift2 );
   
		   // r_old = ( vy0 + gix2 ) * dr;      // global radial position
		   _MM_CVTPS2_PD( vy0a, vy0b, vy0 );
		   r_olda = _mm_mul_pd( _mm_add_pd( vy0a, gix2a ), dr );
		   r_oldb = _mm_mul_pd( _mm_add_pd( vy0b, gix2b ), dr );
   
		   // x2_new = r_old + vv2 * dt;
		   // x3_new =         vv3 * dt;

		   _MM_CVTPS2_PD( vv2a, vv2b, vv2 );
		   _MM_CVTPS2_PD( vv3a, vv3b, vv3 );

		   x2_newa = _mm_add_pd( r_olda, _mm_mul_pd( vv2a, vdt ) );
		   x3_newa =                     _mm_mul_pd( vv3a, vdt ) ;

		   x2_newb = _mm_add_pd( r_oldb, _mm_mul_pd( vv2b, vdt ) );
		   x3_newb =                     _mm_mul_pd( vv3b, vdt );
   
		   // r_new = sqrt( x2_new*x2_new + x3_new*x3_new );
		   r_newa = _mm_sqrt_pd( _mm_add_pd( _mm_mul_pd( x2_newa, x2_newa ),
											 _mm_mul_pd( x3_newa, x3_newa ) ) ); 
		   
		   r_newb = _mm_sqrt_pd( _mm_add_pd( _mm_mul_pd( x2_newb, x2_newb ),
											 _mm_mul_pd( x3_newb, x3_newb ) ) ); 
		   
		   // This is a protection against roundoff for cold plasmas
		   // vy1 = ( r_old == r_new ) ? vy0 : r_new * rdr - gix2);
		   vy1a = _mm_blendv_pd( _mm_sub_pd( _mm_mul_pd( r_newa, rdr ), gix2a ), vy0a, 
		                      _mm_cmpeq_pd( r_olda, r_newa ) );
		   vy1b = _mm_blendv_pd( _mm_sub_pd( _mm_mul_pd( r_newb, rdr ), gix2b ), vy0b, 
		                      _mm_cmpeq_pd( r_oldb, r_newb ) );
		   
		   // Convert vy1 to single precision
		   vy1 = _MM_CVTPD2_PS( vy1a, vy1b );
		   
		   // Correct p_r and p_\theta to conserve angular momentum
		   // tmp = 1.0 / r_new;
		   // vu2 = ( vu2 * x2_new + vu3 * x3_new ) * tmp;
		   // vu3 = vu3 * r_old * tmp;
		   
		   // Convert vu2, vu3 to double precision
		   _MM_CVTPS2_PD( vu2a, vu2b, vu2 );
		   _MM_CVTPS2_PD( vu3a, vu3b, vu3 );
		   
		   tmpa  = _mm_div_pd( _mm_set1_pd(1.0), r_newa );
		   vu2a  = _mm_mul_pd( _mm_add_pd( _mm_mul_pd( vu2a, x2_newa ), 
										   _mm_mul_pd( vu3a, x3_newa ) ), 
										   tmpa );
		   vu3a = _mm_mul_pd( vu3a, _mm_mul_pd( r_olda, tmpa ) );
   
		   tmpb  = _mm_div_pd( _mm_set1_pd(1.0), r_newb );
		   vu2b = _mm_mul_pd( _mm_add_pd( _mm_mul_pd( vu2b, x2_newb ), 
										  _mm_mul_pd( vu3b, x3_newb ) ), 
										  tmpb );
		   vu3b = _mm_mul_pd( vu3b, _mm_mul_pd( r_oldb, tmpb ) );
		   
		   vu2 = _MM_CVTPD2_PS( vu2a, vu2b );
		   vu3 = _MM_CVTPD2_PS( vu3a, vu3b );
		}
		
		// Store Momenta
		_MM_STORE4v3_PS(u, vu1, vu2, vu3); 

        // Store virtual particles with positions still indexed to the
        // original cell
        STORE4VP2D( vpbuf.p, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );
	  
	    // Trim positions and store results
		vtr1 = vmntrim(vx1);
		vtr2 = vmntrim(vy1);

		vx1  = _mm_sub_ps( vx1, vtr1 );
		vy1  = _mm_sub_ps( vy1, vtr2 );
		_MM_STORE4v2_PS( x, vx1, vy1 );
		
		vitr1 = _mm_cvttps_epi32( vtr1 ) ;
		vitr2 = _mm_cvttps_epi32( vtr2 ) ;
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
	  x  += 8;  
	  u  += 12; 
	  ix += 8;
	  q  += 4;
  
	}
  
	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );
		
  }
    
}


/********************************** 3D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_3D
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  float *efield, float *bfield, int *emf_size, int *emf_offset,  
  float *j, int *j_size, int *j_offset, 
  double *dx, double *dt )
{
  int k, i, np, np_total;
    
  DECLARE_ALIGNED_16( float jnorm[3] );

  DECLARE_ALIGNED_16( int cross[4] );
  DECLARE_ALIGNED_16( int dix[4] );
  DECLARE_ALIGNED_16( int diy[4] );
  DECLARE_ALIGNED_16( int diz[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_vpbuf3D vpbuf;

  // Simulation constants
  __m128 const vtem    = _mm_set1_ps( 0.5 * (*dt) / (*rqm) );
  __m128 const vdt_dx1 = _mm_set1_ps( (*dt) / dx[0] );
  __m128 const vdt_dx2 = _mm_set1_ps( (*dt) / dx[1] );
  __m128 const vdt_dx3 = _mm_set1_ps( (*dt) / dx[2] );

  // get pointers to position -1,-1,-1 of each field component
  int const deltaX = 3; // 3 field components
  int const deltaY = emf_size[0] * deltaX;
  int const deltaZ = emf_size[1] * deltaY;
  
  float* const e1 = efield + (emf_offset[0] - OFFSET)*deltaX + 
                             (emf_offset[1] - OFFSET)*deltaY + 
                             (emf_offset[2] - OFFSET)*deltaZ;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0] - OFFSET)*deltaX + 
                             (emf_offset[1] - OFFSET)*deltaY + 
                             (emf_offset[2] - OFFSET)*deltaZ;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;
  
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

  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np_total % 4 ) {
	
	// printf("Adding dummy particles to the end of particle buffer\n");
	// printf("Total particles : %i, remainder (np_total/4): %i \n",
	//		  np_total, np_total % 4);

	for( i = 0; i < 4 - np_total%4; i ++ ) {
	  // printf(" idx = %d \n ", (np_total + i) ); 
	  
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
    
    np_total += ( 4 - np_total%4 );
    if ( np_total % 4 ) {
       printf("(*error*) Number of particles to push is still not a multiple of 4!\n");
       exit(1);
    }
    
    // printf("new np_total = %i", np_total ); 
  }
  
  if ( p_cache_size_3D % 4 ) {
       printf("(*error*) p_cache_size_3D is not a multiple of 4!\n");
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_3D ) {
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_3D ) ? np_total - k : p_cache_size_3D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
	vpbuf.p = (float *) &(vpbuf.buf);
		
	// Push all particles in group
	for( i = 0; i < np; i+=4 ) {
	  __m128 vx0, vx1, vy0, vy1, vz0, vz1;
	  __m128i vix0, vix1, viy0, viy1, viz0, viz1;
	  __m128 vq;
  
	  register __m128 ve1, ve2, ve3;
	  register __m128 vb1, vb2, vb3;
  
	  __m128 vu1, vu2, vu3;
	    
	  // Load particle positions 
	  _MM_LOAD4v3_PS( vx0, vy0, vz0, x );
	  _MM_LOAD4v3_EPI32( vix0, viy0, viz0, ix  );
  
	  // Interpolate fields
	  {
		 __m128  vx0h, vy0h, vz0h;
		 __m128i vidxi, vidxj, vidxk, vidxih, vidxjh, vidxkh;
		 __m128 vwx[NP], vwy[NP], vwz[NP], vwxh[NP], vwyh[NP], vwzh[NP];

	     // idxi = 3 * ix0
	     vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
	     // idxj = deltaY * iy0
	     vidxj = vsmul( viy0, deltaY );
	     // idxk = deltaZ * iz0
	     vidxk = vsmul( viz0, deltaZ );
	 
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
			float *fld1[4], *fld2[4], *fld3[4];
			DECLARE_ALIGNED_16( int idx1[4] );
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, 
			                  _mm_add_epi32( _mm_add_epi32( vidxih, vidxj ), vidxk ) );
			_mm_store_si128( (__m128i *) idx2, 
			                  _mm_add_epi32( _mm_add_epi32( vidxi, vidxjh ), vidxk ) );
			_mm_store_si128( (__m128i *) idx3, 
			                  _mm_add_epi32( _mm_add_epi32( vidxi, vidxj ), vidxkh ) );

            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 4; k1++) {
              fld1[k1] = e1 + idx1[k1];
              fld2[k1] = e2 + idx2[k1];
              fld3[k1] = e3 + idx3[k1];
			}
			
			ve1 = _mm_setzero_ps();
			ve2 = _mm_setzero_ps();
			ve3 = _mm_setzero_ps();
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register __m128 f1plane, f2plane, f3plane;
			   f1plane = _mm_setzero_ps();
			   f2plane = _mm_setzero_ps();
			   f3plane = _mm_setzero_ps();
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register __m128 f1line, f2line, f3line;
				 f1line = _mm_setzero_ps();
				 f2line = _mm_setzero_ps();
				 f3line = _mm_setzero_ps();

				 __m128 f1point[NP], f2point[NP], f3point[NP];
				 for ( k1 = 0; k1 < NP; k1++ ) {
					int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point[k1] = LOADFLD4( fld1, shift ); 
					f2point[k1] = LOADFLD4( fld2, shift ); 
					f3point[k1] = LOADFLD4( fld3, shift ); 
				 }
				 
				 for ( k1 = 0; k1 < NP; k1 ++ ) {
					f1line = _mm_add_ps( f1line, _mm_mul_ps( f1point[k1], vwxh[k1] ) );
					f2line = _mm_add_ps( f2line, _mm_mul_ps( f2point[k1],  vwx[k1] ) );
					f3line = _mm_add_ps( f3line, _mm_mul_ps( f3point[k1],  vwx[k1] ) );
				 }
				 
			     f1plane = _mm_add_ps( f1plane, _mm_mul_ps( f1line,  vwy[k2] ) );
			     f2plane = _mm_add_ps( f2plane, _mm_mul_ps( f2line, vwyh[k2] ) );
			     f3plane = _mm_add_ps( f3plane, _mm_mul_ps( f3line,  vwy[k2] ) );
			   }
			   
			   ve1 = _mm_add_ps( ve1, _mm_mul_ps( f1plane,  vwz[k3] ) );
			   ve2 = _mm_add_ps( ve2, _mm_mul_ps( f2plane,  vwz[k3] ) );
			   ve3 = _mm_add_ps( ve3, _mm_mul_ps( f3plane, vwzh[k3] ) );
			}
         
         }		 

         // Interpolate B field
	     {
			int k1, k2, k3;
			float *fld1[4], *fld2[4], *fld3[4];
			DECLARE_ALIGNED_16( int idx1[4] );
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, 
			                  _mm_add_epi32( _mm_add_epi32( vidxi, vidxjh ), vidxkh ) );
			_mm_store_si128( (__m128i *) idx2, 
			                  _mm_add_epi32( _mm_add_epi32( vidxih, vidxj ), vidxkh ) );
			_mm_store_si128( (__m128i *) idx3, 
			                  _mm_add_epi32( _mm_add_epi32( vidxih, vidxjh ), vidxk ) );

            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 4; k1++) {
              fld1[k1] = b1 + idx1[k1];
              fld2[k1] = b2 + idx2[k1];
              fld3[k1] = b3 + idx3[k1];
			}
			
			vb1 = _mm_setzero_ps();
			vb2 = _mm_setzero_ps();
			vb3 = _mm_setzero_ps();
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register __m128 f1plane, f2plane, f3plane;
			   f1plane = _mm_setzero_ps();
			   f2plane = _mm_setzero_ps();
			   f3plane = _mm_setzero_ps();
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register __m128 f1line, f2line, f3line;
				 f1line = _mm_setzero_ps();
				 f2line = _mm_setzero_ps();
				 f3line = _mm_setzero_ps();

				 __m128 f1point[NP], f2point[NP], f3point[NP];
				 for ( k1 = 0; k1 < NP; k1++ ) {
					int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point[k1] = LOADFLD4( fld1, shift ); 
					f2point[k1] = LOADFLD4( fld2, shift ); 
					f3point[k1] = LOADFLD4( fld3, shift ); 
				 }
				 
				 for ( k1 = 0; k1 < NP; k1 ++ ) {
					f1line = _mm_add_ps( f1line, _mm_mul_ps( f1point[k1],  vwx[k1] ) );
					f2line = _mm_add_ps( f2line, _mm_mul_ps( f2point[k1], vwxh[k1] ) );
					f3line = _mm_add_ps( f3line, _mm_mul_ps( f3point[k1], vwxh[k1] ) );
				 }
				 
			     f1plane = _mm_add_ps( f1plane, _mm_mul_ps( f1line, vwyh[k2] ) );
			     f2plane = _mm_add_ps( f2plane, _mm_mul_ps( f2line,  vwy[k2] ) );
			     f3plane = _mm_add_ps( f3plane, _mm_mul_ps( f3line, vwyh[k2] ) );
			   }
			   
			   vb1 = _mm_add_ps( vb1, _mm_mul_ps( f1plane, vwzh[k3] ) );
			   vb2 = _mm_add_ps( vb2, _mm_mul_ps( f2plane, vwzh[k3] ) );
			   vb3 = _mm_add_ps( vb3, _mm_mul_ps( f3plane,  vwz[k3] ) );
			}
         
         }		 
		 
	  }	  
	  
	  // ---------- advance momenta
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, u, vu1, vu2, vu3 );

      // Store Results
      _MM_STORE4v3_PS(u, vu1, vu2, vu3); 

	  // advance u
	  u  += 12; 
   
	  // ---------- advance position and get velocities
	  {
		register __m128 vrg;
		register __m128 vtrx, vtry, vtrz;
		register __m128i vitrx, vitry, vitrz;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );
	  
		// get velocities
		vu1 = _mm_mul_ps( vu1, vrg );
		vu2 = _mm_mul_ps( vu2, vrg );
		vu3 = _mm_mul_ps( vu3, vrg );
			
		vx1 = _mm_add_ps( vx0, _mm_mul_ps( vu1, vdt_dx1 ) );
		vy1 = _mm_add_ps( vy0, _mm_mul_ps( vu2, vdt_dx2 ) );
		vz1 = _mm_add_ps( vz0, _mm_mul_ps( vu3, vdt_dx3 ) );
	  
		vtrx = vmntrim(vx1);
		vtry = vmntrim(vy1);
		vtrz = vmntrim(vz1);
		
		vq = _mm_load_ps( q );
		
		STORE4VP3D( vpbuf.p, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix0, viy0, viz0 );

		vx1   = _mm_sub_ps( vx1, vtrx );
		vy1   = _mm_sub_ps( vy1, vtry );
		vz1   = _mm_sub_ps( vz1, vtrz );
        _MM_STORE4v3_PS( x, vx1, vy1, vz1 );

		vitrx = _mm_cvttps_epi32( vtrx );
		vitry = _mm_cvttps_epi32( vtry );
		vitrz = _mm_cvttps_epi32( vtrz );
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
	  x  += 12;  
	  ix += 12;
	  q  += 4;
  
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
