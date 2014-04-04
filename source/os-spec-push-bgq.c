/*****************************************************************************************

Relativistic particle pusher, BG/Q optimized version

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-push-bgq.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "fortran.h"

#include "vector-bgq.h"
#include "splines-bgq.h"

#include "os-spec-current-bgq.h"


/*****************************************************************************************
vntrim

Returns the required cell shift (trim) for particles that must remain in [ -0.5, 0.5 [
to the nearest cell point.

dx = ( vx < -0.5 ) ? -1 : 0 + ( vx >= +0.5 ) ? +1 : 0

x'  = x  - dx
ix' = ix + dx

*****************************************************************************************/

inline vector vntrim( const vector vx )
{
  const vector oneHalf = vec_splats( 0.5 );
  const vector one     = vec_splats( 1.0 );
  const vector zero    = vec_splats( 0.0 );
  
  vector va = vec_sel( one, zero, vec_add( oneHalf, vx ) );
  vector vb = vec_sel( zero, one, vec_sub( vx, oneHalf ) );
  
  return vec_sub( vb, va );
}


/*****************************************************************************************
vdudt_boris

Use a Boris push to advance particle velocities using interpolated
fields.
*****************************************************************************************/

#define VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, u, vu1, vu2, vu3 ) \
{																			\
   vector const one = vec_splats( 1.0 );									\
   register vector vut1, vut2, vut3;										\
   register vector vgamma_tem, votsq;										\
   																			\
   ve1 = vec_mul( ve1, vtem );												\
   ve2 = vec_mul( ve2, vtem );												\
   ve3 = vec_mul( ve3, vtem );												\
      																		\
   VEC_LD4v3( vu1, vu2, vu3, u );											\
   																			\
   vut1 = vec_add( vu1, ve1 );												\
   vut2 = vec_add( vu2, ve2 );												\
   vut3 = vec_add( vu3, ve3 );												\
   																			\
   vgamma_tem = vec_madd( vut1, vut1, one );								\
   vgamma_tem = vec_madd( vut2, vut2, vgamma_tem );							\
   vgamma_tem = vec_madd( vut3, vut3, vgamma_tem );							\
   vgamma_tem = vec_rsqrt( vgamma_tem );									\
   vgamma_tem = vec_mul( vtem, vgamma_tem );								\
   																			\
   vb1 = vec_mul( vb1, vgamma_tem );										\
   vb2 = vec_mul( vb2, vgamma_tem );										\
   vb3 = vec_mul( vb3, vgamma_tem );										\
   																			\
   vu1 = vec_madd( vb3, vut2, vut1 );										\
   vu2 = vec_madd( vb1, vut3, vut2 );										\
   vu3 = vec_madd( vb2, vut1, vut3 );										\
																			\
   vu1 = vec_nmsub( vb2, vut3, vu1 );										\
   vu2 = vec_nmsub( vb3, vut1, vu2 );										\
   vu3 = vec_nmsub( vb1, vut2, vu3 );										\
																			\
   votsq = vec_madd( vb1,   vb1, one );										\
   votsq = vec_madd( vb2, vb2, votsq );										\
   votsq = vec_madd( vb3, vb3, votsq );										\
   																			\
   votsq = vec_r(votsq) ;													\
   votsq = vec_add( votsq, votsq );											\
																			\
   vb1 = vec_mul( vb1, votsq );												\
   vb2 = vec_mul( vb2, votsq );												\
   vb3 = vec_mul( vb3, votsq );												\
   																			\
   vut1 = vec_madd( vb3, vu2, vut1 );										\
   vut2 = vec_madd( vb1, vu3, vut2 );										\
   vut3 = vec_madd( vb2, vu1, vut3 );										\
																			\
   vut1 = vec_nmsub( vb2, vu3, vut1 );										\
   vut2 = vec_nmsub( vb3, vu1, vut2 );										\
   vut3 = vec_nmsub( vb1, vu2, vut3 );										\
   																			\
   vu1 = vec_add( vut1, ve1 );												\
   vu2 = vec_add( vut2, ve2 );												\
   vu3 = vec_add( vut3, ve3 );												\
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

inline void vsplit2D( t_vpbuf2D* const vpbuf, const int cross[], const int dix[], const int diy[] )
{

  int k, np;
  t_vp2D *vp;

  vp = (t_vp2D *) (vpbuf -> p);

  np = 4;
  for( k = 0 ; k < 4; k++ ) {

     double delta, xint, yint, vzint, xint4, yint4, vzint4;
    
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
           vzint4 = vp[k].vz * (1-delta);
           yint4 = 0.5 * diy[k];
		   delta = ( yint4 - yint ) / (vp[k].y1 - yint);
		   xint4 = -xint + ( vp[k].x1 - xint ) * delta;

           // add extra particle (1st segment y-split)
		   vp[np].x0 =  -xint; 	vp[np].x1 = xint4;
		   vp[np].y0 = yint; 	vp[np].y1 = yint4;
		   vp[np].ix = vp[k].ix + dix[k]; 
		   vp[np].iy = vp[k].iy;
		   vp[np].q  = vp[k].q;  
		   vp[np].vz = vzint4 * delta;
		   np++;


           // add extra particle (2nd segment y-split)
		   vp[np].x0 =  xint4; 	vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 = -yint4; 	vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].ix = vp[k].ix + dix[k]; 
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].q  = vp[k].q;  
		   vp[np].vz = vzint4 * (1-delta);
		   np++;

		   // correct existing particle (1st segment)  
		   vp[k].x1 = xint;
		   vp[k].y1 = yint;
		   vp[k].vz = vzint;

         } else {
         
           vzint = vp[k].vz*delta;
		   // y split 1st vp
		   yint4 = 0.5 * diy[k];
		   delta = ( yint4 - vp[k].y0 ) / ( yint - vp[k].y0);
		   xint4 = vp[k].x0 + (xint - vp[k].x0) * delta;

           // add extra particle (2nd segment y-split)
		   vp[np].x0 =  xint4; vp[np].x1 = xint;
		   vp[np].y0 = -yint4; vp[np].y1 = yint - diy[k];
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
		   vp[k].x1 = xint4;
		   vp[k].y1 = yint4;
		   vp[k].vz = vzint * delta;

         }
		 break;
     }
     
  }
  
  // Increase buffer pointer
  vp += np;
  vpbuf -> p = (double *) vp;
  vpbuf -> np += np;
}

/*****************************************************************************************
vsplit3D

Splits 4 particle trajectories for current deposition and stores segment on a special virtual
particle buffer which is optimized for the vector current deposition routines
*****************************************************************************************/

inline void vsplit3D( t_vpbuf3D* const vpbuf, const int cross[], 
               const int dix[], const int diy[], const int diz[] )
{

  int k, np;
  t_vp3D *vp;

  np = 4;
  vp = (t_vp3D *) (vpbuf -> p);

  for( k = 0; k < 4; k++ ) {
     double xint, yint, zint;
     double xint4, yint4, zint4;
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
           yint4 = 0.5 * diy[k];
		   delta = ( yint4 - yint ) / (vp[k].y1 - yint);
		   xint4 = -xint + ( vp[k].x1 - xint ) * delta;
		   zint4 =  zint + ( vp[k].z1 - zint ) * delta;
		   
		   vp[np].x0 = -xint; 		vp[np].x1 = xint4;
		   vp[np].y0 =  yint; 		vp[np].y1 = yint4;
		   vp[np].z0 =  zint; 		vp[np].z1 = zint4;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

		   vp[np].x0 =  xint4; 		vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 = -yint4; 		vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].z0 =  zint4; 		vp[np].z1 = vp[k].z1;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

           // no y cross on 1st vp
           vp[k].x1 = xint; vp[k].y1 = yint; vp[k].z1 = zint;
           
         } else {
         
		   // y split 1st vp
		   yint4 = 0.5 * (diy[k]);
		   delta = ( yint4 - vp[k].y0 ) / ( yint - vp[k].y0);
		   xint4 = vp[k].x0 + (xint - vp[k].x0) * delta;
		   zint4 = vp[k].z0 + (zint - vp[k].z0) * delta;

		   vp[np].x0 =  xint4; 		vp[np].x1 = xint;
		   vp[np].y0 = -yint4; 		vp[np].y1 = yint - diy[k];
		   vp[np].z0 =  zint4; 		vp[np].z1 = zint;
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
           vp[k].x1 = xint4; vp[k].y1 = yint4; vp[k].z1 = zint4;

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
           zint4 = 0.5 * (diz[k]);
		   delta = ( zint4 - zint ) / (vp[k].z1 - zint);
		   xint4 = -xint + ( vp[k].x1 - xint ) * delta;
		   yint4 =  yint + ( vp[k].y1 - yint ) * delta;
	
		   vp[np].x0 = -xint; 		vp[np].x1 = xint4;
		   vp[np].y0 =  yint; 		vp[np].y1 = yint4;
		   vp[np].z0 =  zint; 		vp[np].z1 = zint4;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

		   vp[np].x0 =  xint4; 		vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 =  yint4; 		vp[np].y1 = vp[k].y1;
		   vp[np].z0 = -zint4; 		vp[np].z1 = vp[k].z1 - diz[k];
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
		   zint4 = 0.5 * (diz[k]);
		   delta = ( zint4 - vp[k].z0 ) / ( zint - vp[k].z0);
		   xint4 = vp[k].x0 + (xint - vp[k].x0) * delta;
		   yint4 = vp[k].y0 + (yint - vp[k].y0) * delta;

		   vp[np].x0 =  xint4; 		vp[np].x1 = xint;
		   vp[np].y0 =  yint4; 		vp[np].y1 = yint;
		   vp[np].z0 = -zint4; 		vp[np].z1 = zint - diz[k];
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
		   vp[k].x1 = xint4;
		   vp[k].y1 = yint4;
		   vp[k].z1 = zint4;

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
           zint4 = 0.5 * diz[k];
		   delta = ( zint4 - zint ) / (vp[k].z1 - zint);
		   xint4 =   xint + ( vp[k].x1 - xint ) * delta;
		   yint4 =  -yint + ( vp[k].y1 - yint ) * delta;
		   
		   vp[np].x0 =  xint; 		vp[np].x1 = xint4;
		   vp[np].y0 = -yint; 		vp[np].y1 = yint4;
		   vp[np].z0 =  zint; 		vp[np].z1 = zint4;
		   vp[np].ix = vp[k].ix;
		   vp[np].iy = vp[k].iy + diy[k];
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

		   vp[np].x0 =  xint4; 		vp[np].x1 = vp[k].x1;
		   vp[np].y0 =  yint4; 		vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].z0 = -zint4;	 	vp[np].z1 = vp[k].z1 - diz[k];
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
		   zint4 = 0.5 * diz[k];
		   delta = ( zint4 - vp[k].z0 ) / ( zint - vp[k].z0);
		   xint4 = vp[k].x0 + (xint - vp[k].x0) * delta;
		   yint4 = vp[k].y0 + (yint - vp[k].y0) * delta;

		   vp[np].x0 =  xint4; 		vp[np].x1 = xint;
		   vp[np].y0 =  yint4; 		vp[np].y1 = yint;
		   vp[np].z0 = -zint4; 		vp[np].z1 = zint - diz[k];
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
		   vp[k].x1 = xint4;
		   vp[k].y1 = yint4;
		   vp[k].z1 = zint4;

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
           yint4 = 0.5 * (diy[k]);
		   delta = ( yint4 - yint ) / (vp[k].y1 - yint);
		   xint4 = -xint + ( vp[k].x1 - xint ) * delta;
		   zint4 =  zint + ( vp[k].z1 - zint ) * delta;
		   
		   vp[np].x0 = -xint; 		vp[np].x1 = xint4;
		   vp[np].y0 =  yint; 		vp[np].y1 = yint4;
		   vp[np].z0 =  zint; 		vp[np].z1 = zint4;
		   vp[np].ix = vp[k].ix + dix[k];
		   vp[np].iy = vp[k].iy;
		   vp[np].iz = vp[k].iz;
		   vp[np].q  = vp[k].q;
		   np ++;

		   vp[np].x0 =  xint4; 		vp[np].x1 = vp[k].x1 - dix[k];
		   vp[np].y0 = -yint4; 		vp[np].y1 = vp[k].y1 - diy[k];
		   vp[np].z0 =  zint4; 		vp[np].z1 = vp[k].z1;
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
		   yint4 = 0.5 * (diy[k]);
		   delta = ( yint4 - vp[k].y0 ) / ( yint - vp[k].y0);
		   xint4 = vp[k].x0 + (xint - vp[k].x0) * delta;
		   zint4 = vp[k].z0 + (zint - vp[k].z0) * delta;

		   vp[np].x0 =  xint4; 		vp[np].x1 = xint;
		   vp[np].y0 = -yint4; 		vp[np].y1 = yint - diy[k];
		   vp[np].z0 =  zint4; 		vp[np].z1 = zint;
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
           vp[k].x1 = xint4; 
           vp[k].y1 = yint4; 
           vp[k].z1 = zint4;

         }
         
 		  // one of the 3 vp requires an additional z split
          int vpidx[3] = {k, np-2, np-1};
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

  vp += np;
  vpbuf -> p = (double *) vp;
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

#define ADVANCE_DEPOSIT_2D FNAME( vadvance_deposit_2d_s1 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s1 )
#define ADVANCE_DEPOSIT_3D FNAME( vadvance_deposit_3d_s1 )

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
( int* const ix, t_real* const x, t_real* const u, t_real* const q, 
  int const* const i0, int const * const i1, double const * const rqm,
  t_real * const e, t_real * const b, 
  int const * const emf_size, int const * const emf_offset,  
  t_real * const j, int const * const j_size, int const * const j_offset, 
  double const * const dx, double const * const dt )
{
  
  t_real *e1, *e2, *e3;
  t_real *b1, *b2, *b3;
  
  int *pix;
  t_real *px, *pu, *pq;
  
  int deltaX, deltaY, k, i, np;
  int np_total;
  
  vector vdt_dx, vdt_dy, vtem;
  
  DECLARE_ALIGNED_32( double jnorm[2] );
  
  DECLARE_ALIGNED_32( int cross[4] );
  DECLARE_ALIGNED_32( int dix[4] );
  DECLARE_ALIGNED_32( int diy[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_vpbuf2D vpbuf;

   // get pointers to position 0,0 of each field component
  deltaX = 3; // 3 field components 
  deltaY = emf_size[0] * deltaX;

  e1 = e   + (emf_offset[0]-OFFSET)*deltaX + (emf_offset[1]-OFFSET)*deltaY;
  e2 = e1 + 1;
  e3 = e1 + 2;

  b1 = b   + (emf_offset[0]-OFFSET)*deltaX + (emf_offset[1]-OFFSET)*deltaY;
  b2 = b1  + 1;
  b3 = b1  + 2;
  
  vtem   = vec_splats( 0.5 * (*dt) / (*rqm) );
  vdt_dx = vec_splats( (*dt) / dx[0] );
  vdt_dy = vec_splats( (*dt) / dx[1] );

  
  // Normalization for currents
  // The factor 2 comes from a simplification on the calculation of parallel current weights
  jnorm[0] = dx[0] / (*dt) / (2*JNORM);
  jnorm[1] = dx[1] / (*dt) / (2*JNORM);
  
  
  // jump to 1st particle
  px  = x  + 2*(*i0-1);
  pix = ix + 2*(*i0-1);
  pu  = u  + 3*(*i0-1);
  pq  = q  + (*i0-1);
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  if ( p_cache_size_2D % 4 ) {
       printf("(*error*) p_cache_size_2D is not a multiple of 4!\n");
       exit(1);
  }
  
  // If number of particles in buffer is not multiple of 2 add dummy particles to the end
  if ( np_total % 4 ) {
	for( i = 0; i < 4 - np_total%4; i ++ ) {
	  pix[ 2 * (np_total + i)     ] = 1;
	  pix[ 2 * (np_total + i) + 1 ] = 1;
	  px[ 2 * (np_total + i)     ] = 0.;
	  px[ 2 * (np_total + i) + 1 ] = 0.;
	  pu[ 3 * (np_total + i)     ] = 0.;
	  pu[ 3 * (np_total + i) + 1 ] = 0.;
	  pu[ 3 * (np_total + i) + 2 ] = 0.;
	  pq[ (np_total + i ) ] = 0.;
	}
	
	np_total += (4 - np_total%4);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
	vpbuf.p = (double *) (&vpbuf.buf[0]);
		
	// Push all particles in group
	for( i = 0; i < np; i+=4 ) {
	  vector vx0, vy0;
	  DECLARE_ALIGNED_32( int vix[4] );
	  DECLARE_ALIGNED_32( int viy[4] );
  
	  vector ve1, ve2, ve3;
	  vector vb1, vb2, vb3;
  
	  vector vu1, vu2, vu3;	  
	  vector vq;
	  
	  int l;
  
	  // Load particle positions 
	  VEC_LD4v2( vx0, vy0, px );
	  
	  for( l = 0; l < 4; l++ ) {
	    vix[l] = pix[2*l    ];
	    viy[l] = pix[2*l + 1];
	  }
 
	  // Interpolate fields
	  {

		 vector  vx0h, vy0h;
		 int idxi[4], idxj[4], idxih[4], idxjh[4];
		 register vector vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	     
	     int l;
	     
	     for( l=0; l<4; l++) {
	       idxi[l] =      3 * vix[l];
	       idxj[l] = deltaY * viy[l];
	     }

		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy );

		 vget_hp_near( vx0, idxi, vx0h, idxih, 3 ); 
		 vget_hp_near( vy0, idxj, vy0h, idxjh, deltaY ); 

		 SPLINE( vx0h, vwxh );
		 SPLINE( vy0h, vwyh );

		 // Interpolate E field
         {
			int k1, k2;
			t_real *fld1[4], *fld2[4], *fld3[4];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < 4; k1++) {
			   fld1[k1] = e1 + idxih[k1] + idxj[k1];
			   fld2[k1] = e2 + idxi[k1] + idxjh[k1];
			   fld3[k1] = e3 + idxi[k1] + idxj[k1];
			}
						 
			ve1 = zero;
			ve2 = zero;
			ve3 = zero;
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register vector f1line, f2line, f3line;
			  f1line = zero;
			  f2line = zero;
			  f3line = zero;

			  register vector f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = vec_madd( vwxh[k1], f1point[k1], f1line );
				 f2line = vec_madd( vwx[k1],  f2point[k1], f2line );
				 f3line = vec_madd( vwx[k1],  f3point[k1], f3line );
			  }

			  ve1 = vec_madd( vwy[k2],  f1line, ve1 );
			  ve2 = vec_madd( vwyh[k2], f2line, ve2 );
			  ve3 = vec_madd( vwy[k2],  f3line, ve3 );
			}
			        
         }		 

		 
		 // Interpolate B field
         {
			int k1, k2;
			t_real *fld1[4], *fld2[4], *fld3[4];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < 4; k1++) {
			   fld1[k1] = b1 + idxi[k1]  + idxjh[k1];
			   fld2[k1] = b2 + idxih[k1] + idxj[k1];
			   fld3[k1] = b3 + idxih[k1] + idxjh[k1];
			}
			
			vb1 = zero;
			vb2 = zero;
			vb3 = zero;
		  
			for ( k2 = 0; k2 < 3; k2++ ) {
			  
			  register vector f1line, f2line, f3line;
			  f1line = zero;
			  f2line = zero;
			  f3line = zero;

			  register vector f1point[NP], f2point[NP], f3point[NP];
			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
              
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = vec_madd( vwx[k1],  f1point[k1], f1line );
				 f2line = vec_madd( vwxh[k1], f2point[k1], f2line );
				 f3line = vec_madd( vwxh[k1], f3point[k1], f3line );
			  }

			  vb1 = vec_madd( vwyh[k2], f1line, vb1 );
			  vb2 = vec_madd( vwy[k2],  f2line, vb2 );
			  vb3 = vec_madd( vwyh[k2], f3line, vb3 );
			}
         
         }

	  }

	  // ---------- advance momenta
	  // results are stored in vu1, vu2, vu3 (temp) and pu (permanent)
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, pu, vu1, vu2, vu3 );
      VEC_ST4v3(pu, vu1, vu2, vu3);

	  // advance u 4 particles * 3 components
	  pu  += 12; 
   
	  // ---------- advance position and get velocities
	  {
		register vector vrg;
		register vector vx1, vy1, vtrx, vtry;
		int l;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );
		
        vq = vec_lda( 0, pq );
	
		// get velocities
		vu1 = vec_mul( vu1, vrg );
		vu2 = vec_mul( vu2, vrg );
		vu3 = vec_mul( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		vx1 = vec_madd( vdt_dx, vu1, vx0 );
		vy1 = vec_madd( vdt_dy, vu2, vy0 );

        // Store virtual particles with positions still indexed to the
        // original cell
        STORE4VP2D( vpbuf.p, vx0, vx1, vy0, vy1, vq, vu3, vix, viy );
	  
	    // Trim positions and store results
		vtrx = vntrim(vx1);
		vtry = vntrim(vy1);

		vx1  = vec_sub( vx1, vtrx );
		vy1  = vec_sub( vy1, vtry );
		VEC_ST4v2( px, vx1, vy1 );
			
		vtrx = vec_ctiw( vtrx ) ;
		vtry = vec_ctiw( vtry ) ;
		vec_sta( vtrx, 0, dix );
		vec_sta( vtry, 0, diy );
		
        for( l=0; l<4; l++) {
          // store new particle cell indexes
          pix[2*l     ] = vix[l] + dix[l];
          pix[2*l + 1 ] = viy[l] + diy[l];
          
          // calculate crossings - this can be done using vtr1 / vtr2
          cross[l] = (( dix[l] ) ? 1 : 0 ) + (( diy[l] ) ? 2 : 0 );
        }

	  }
	  
	  
	  // ---------- split trajectories for current deposition

	  
	  // The new split uses positions indexed to the original cell
      vsplit2D( &vpbuf, cross, dix, diy );
	  
	  // advance x, ix 4 particles * 2 components
	  px  += 8;  
	  pix += 8;
	  // advance q 4 particles * 1 component
	  pq  += 4;
  
	}

	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );
  }
    
}

extern void ADVANCE_DEPOSIT_2D_CYL
( int* const ix, t_real* const x, t_real* const u, t_real* const q, 
  int const* const i0, int const * const i1, double const * const rqm,
  int const *ilb2 , 
  t_real * const e, t_real * const b, 
  int const * const emf_size, int const * const emf_offset,  
  t_real * const j, int const * const j_size, int const * const j_offset, 
  double const * const dx, double const * const dt )
{
  
  t_real *e1, *e2, *e3;
  t_real *b1, *b2, *b3;
  
  int *pix;
  t_real *px, *pu, *pq;
  
  int deltaX, deltaY, k, i, np;
  int np_total;
  
  vector vdt_dx, vdt_dy, vtem;
  
  DECLARE_ALIGNED_32( double jnorm[2] );
  
  DECLARE_ALIGNED_32( int cross[4] );
  DECLARE_ALIGNED_32( int dix[4] );
  DECLARE_ALIGNED_32( int diy[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_vpbuf2D vpbuf;

   // get pointers to position 0,0 of each field component
  deltaX = 3; // 3 field components 
  deltaY = emf_size[0] * deltaX;

  e1 = e   + (emf_offset[0]-OFFSET)*deltaX + (emf_offset[1]-OFFSET)*deltaY;
  e2 = e1 + 1;
  e3 = e1 + 2;

  b1 = b   + (emf_offset[0]-OFFSET)*deltaX + (emf_offset[1]-OFFSET)*deltaY;
  b2 = b1  + 1;
  b3 = b1  + 2;
  
  vtem   = vec_splats( 0.5 * (*dt) / (*rqm) );
  vdt_dx = vec_splats( (*dt) / dx[0] );
  vdt_dy = vec_splats( (*dt) / dx[1] );

  
  // Normalization for currents
  // The factor 2 comes from a simplification on the calculation of parallel current weights
  jnorm[0] = dx[0] / (*dt) / (2*JNORM);
  jnorm[1] = dx[1] / (*dt) / (2*JNORM);

  vector gshift2 = vec_splats( (*ilb2 - 2 ) - 0.5 );
  vector dr      = vec_splats( dx[1] );
  vector rdr     = vec_splats( 1.0/dx[1] );
  vector vdt     = vec_splats( *dt );
  
  
  // jump to 1st particle
  px  = x  + 2*(*i0-1);
  pix = ix + 2*(*i0-1);
  pu  = u  + 3*(*i0-1);
  pq  = q  + (*i0-1);
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  if ( p_cache_size_2D % 4 ) {
       printf("(*error*) p_cache_size_2D is not a multiple of 4!\n");
       exit(1);
  }
  
  // If number of particles in buffer is not multiple of 2 add dummy particles to the end
  if ( np_total % 4 ) {
	for( i = 0; i < 4 - np_total%4; i ++ ) {
	  pix[ 2 * (np_total + i)     ] = 1;
	  pix[ 2 * (np_total + i) + 1 ] = 1;
	  px[ 2 * (np_total + i)     ] = 0.;
	  px[ 2 * (np_total + i) + 1 ] = 0.;
	  pu[ 3 * (np_total + i)     ] = 0.;
	  pu[ 3 * (np_total + i) + 1 ] = 0.;
	  pu[ 3 * (np_total + i) + 2 ] = 0.;
	  pq[ (np_total + i ) ] = 0.;
	}
	
	np_total += (4 - np_total%4);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
	vpbuf.p = (double *) (&vpbuf.buf[0]);
		
	// Push all particles in group
	for( i = 0; i < np; i+=4 ) {
	  vector vx0, vy0;
	  DECLARE_ALIGNED_32( int vix[4] );
	  DECLARE_ALIGNED_32( int viy[4] );
  
	  vector ve1, ve2, ve3;
	  vector vb1, vb2, vb3;
  
	  vector vu1, vu2, vu3;	  
	  vector vq;

	  int l;
  
	  // Load particle positions 
	  VEC_LD4v2( vx0, vy0, px );
	  
	  for( l = 0; l < 4; l++ ) {
	    vix[l] = pix[2*l    ];
	    viy[l] = pix[2*l + 1];
	  }
 
	  // Interpolate fields
	  {

		 vector  vx0h, vy0h;
		 int idxi[4], idxj[4], idxih[4], idxjh[4];
		 register vector vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	 
	     int l;
	     
	     for( l=0; l<4; l++) {
	       idxi[l] =      3 * vix[l];
	       idxj[l] = deltaY * viy[l];
	     }

		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy );

		 vget_hp_near( vx0, idxi, vx0h, idxih, 3 ); 
		 vget_hp_near( vy0, idxj, vy0h, idxjh, deltaY ); 

		 SPLINE( vx0h, vwxh );
		 SPLINE( vy0h, vwyh );

		 // Interpolate E field
         {
			int k1, k2;
			t_real *fld1[4], *fld2[4], *fld3[4];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < 4; k1++) {
			   fld1[k1] = e1 + idxih[k1] + idxj[k1];
			   fld2[k1] = e2 + idxi[k1] + idxjh[k1];
			   fld3[k1] = e3 + idxi[k1] + idxj[k1];
			}
						 
			ve1 = zero;
			ve2 = zero;
			ve3 = zero;
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register vector f1line, f2line, f3line;
			  f1line = zero;
			  f2line = zero;
			  f3line = zero;

			  register vector f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = vec_madd( vwxh[k1], f1point[k1], f1line );
				 f2line = vec_madd( vwx[k1],  f2point[k1], f2line );
				 f3line = vec_madd( vwx[k1],  f3point[k1], f3line );
			  }

			  ve1 = vec_madd( vwy[k2],  f1line, ve1 );
			  ve2 = vec_madd( vwyh[k2], f2line, ve2 );
			  ve3 = vec_madd( vwy[k2],  f3line, ve3 );
			}
			        
         }		 

		 
		 // Interpolate B field
        {
			int k1, k2;
			t_real *fld1[4], *fld2[4], *fld3[4];
			
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1)
			for( k1 = 0; k1 < 4; k1++) {
			   fld1[k1] = b1 + idxi[k1]  + idxjh[k1];
			   fld2[k1] = b2 + idxih[k1] + idxj[k1];
			   fld3[k1] = b3 + idxih[k1] + idxjh[k1];
			}
			
			vb1 = zero;
			vb2 = zero;
			vb3 = zero;
		  
			for ( k2 = 0; k2 < 3; k2++ ) {
			  
			  register vector f1line, f2line, f3line;
			  f1line = zero;
			  f2line = zero;
			  f3line = zero;

			  register vector f1point[NP], f2point[NP], f3point[NP];
			  for ( k1 = 0; k1 < NP; k1++ ) {
				 int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD4( fld1, shift ); 
				 f2point[k1] = LOADFLD4( fld2, shift ); 
				 f3point[k1] = LOADFLD4( fld3, shift ); 
              }
              
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = vec_madd( vwx[k1],  f1point[k1], f1line );
				 f2line = vec_madd( vwxh[k1], f2point[k1], f2line );
				 f3line = vec_madd( vwxh[k1], f3point[k1], f3line );
			  }

			  vb1 = vec_madd( vwyh[k2], f1line, vb1 );
			  vb2 = vec_madd( vwy[k2],  f2line, vb2 );
			  vb3 = vec_madd( vwyh[k2], f3line, vb3 );
			}
         
         }

	  }

	  // ---------- advance momenta
	  // results are stored in vu1, vu2, vu3 (temp) and pu (permanent)
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, pu, vu1, vu2, vu3 );

   
	  // ---------- advance position and get velocities
	  {
		register vector vrg;
		register vector vx1, vy1, vtrx, vtry;
		register vector vv1, vv2, vv3;

		int l;
  
		VRGAMMA( vrg, vu1, vu2, vu3 );
		
        vq = vec_lda( 0, pq );
	
		// get velocities
		vv1 = vec_mul( vu1, vrg );
		vv2 = vec_mul( vu2, vrg );
		vv3 = vec_mul( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		vx1 = vec_madd( vdt_dx, vv1, vx0 );

        // x2 - 3D like push
        // This further changes u2 and u3
        {  
		   
		   vector gix2;
		   vector r_old, r_new;
		   vector x2_new, x3_new; 
		   vector tmp;
		   
		   
		   // gshift2 = shift_ix2 - 0.5d0
   
		   // gix2   = viy0 + gshift2;          // global cell
	       gix2   = vec_add( vec_cfid( vec_ldiaa( 0, viy ) ), gshift2 );
		   r_old  = vec_mul( vec_add( vy0, gix2 ), dr );

		   x2_new = vec_madd( vdt, vv2, r_old );
		   x3_new = vec_mul( vv3, vdt ) ;

		   r_new  = vec_swsqrt_nochk( vec_madd( x3_new, x3_new, vec_mul( x2_new, x2_new ) ) ); 
		   
		   // This is a protection against roundoff for cold plasmas
		   // vy1 = ( r_old == r_new ) ? vy0 : r_new * rdr - gix2);
		   vy1 = vec_sel(  vec_sub( vec_mul( r_new, rdr ), gix2 ), vy0, 
		                   vec_nabs( vec_sub( r_old, r_new ) ) );
		   		   		   
		   // Correct p_r and p_\theta to conserve angular momentum
		   tmp  = vec_r( r_new );
		   vu2  = vec_mul( vec_add( vec_mul( vu2, x2_new ), vec_mul( vu3, x3_new ) ), tmp );
		   vu3  = vec_mul( vu3 , vec_mul( r_old, tmp ) );
   
		}        

        // Store Momenta
        VEC_ST4v3(pu, vu1, vu2, vu3);

        // Store virtual particles with positions still indexed to the
        // original cell
        STORE4VP2D( vpbuf.p, vx0, vx1, vy0, vy1, vq, vv3, vix, viy );
	  
	    // Trim positions and store results
		vtrx = vntrim(vx1);
		vtry = vntrim(vy1);

		vx1  = vec_sub( vx1, vtrx );
		vy1  = vec_sub( vy1, vtry );
		VEC_ST4v2( px, vx1, vy1 );
			
		vtrx = vec_ctiw( vtrx ) ;
		vtry = vec_ctiw( vtry ) ;
		vec_sta( vtrx, 0, dix );
		vec_sta( vtry, 0, diy );
		
        for( l=0; l<4; l++) {
          // store new particle cell indexes
          pix[2*l     ] = vix[l] + dix[l];
          pix[2*l + 1 ] = viy[l] + diy[l];
          
          // calculate crossings - this can be done using vtr1 / vtr2
          cross[l] = (( dix[l] ) ? 1 : 0 ) + (( diy[l] ) ? 2 : 0 );
        }

	  }
	  
	  
	  // ---------- split trajectories for current deposition

	  
	  // The new split uses positions indexed to the original cell
      vsplit2D( &vpbuf, cross, dix, diy );
	  
	  // advance x, ix 4 particles * 2 components
	  px  += 8;  
	  pix += 8;
	  // advance u 4 particles * 3 components
	  pu  += 12; 
	  // advance q 4 particles * 1 component
	  pq  += 4;
  
	}

	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );
  }
    
}


/********************************** 3D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_3D
( int* const ix, t_real * const x, t_real * const u, t_real * const q, 
  int const* const i0, int const* const i1, double const * const rqm,
  t_real  * const e, t_real  * const b, 
  int const * const emf_size, int const * const emf_offset,  
  t_real  * const j, int const * const j_size, int const * const j_offset, 
  double const * const dx, double const * const dt )
{
  t_real  *e1, *e2, *e3;
  t_real  *b1, *b2, *b3;
  
  int *pix;
  t_real  *px, *pu, *pq;
  
  int deltaX, deltaY, deltaZ, k, i, np;
  int np_total;
    
  vector vdt_dx, vdt_dy, vdt_dz, vtem;
    
  DECLARE_ALIGNED_32( double jnorm[3] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_vpbuf3D vpbuf;

  // get pointers to position -1,-1,-1 of each field component
  deltaX = 3; // 3 field components
  deltaY = emf_size[0] * deltaX;
  deltaZ = emf_size[1] * deltaY;


  e1 = e  + (emf_offset[0]-OFFSET)*deltaX +  
            (emf_offset[1]-OFFSET)*deltaY + 
            (emf_offset[2]-OFFSET)*deltaZ;
  e2 = e1 + 1;
  e3 = e1 + 2;

  b1 = b   + (emf_offset[0]-OFFSET)*deltaX + 
             (emf_offset[1]-OFFSET)*deltaY + 
             (emf_offset[2]-OFFSET)*deltaZ;
  b2 = b1  + 1;
  b3 = b1  + 2;
  
  //
  
  vtem   = vec_splats( 0.5 * (*dt) / (*rqm) );
  vdt_dx = vec_splats( (*dt) / dx[0] );
  vdt_dy = vec_splats( (*dt) / dx[1] );
  vdt_dz = vec_splats( (*dt) / dx[2] );

  
  // Normalization for currents
  // The factor 3 comes from a simplification on the calculation of parallel current weights
  jnorm[0] = dx[0] / (*dt) / ( 3 * JNORM );
  jnorm[1] = dx[1] / (*dt) / ( 3 * JNORM );
  jnorm[2] = dx[2] / (*dt) / ( 3 * JNORM );
  
  // jump to 1st particle
  px  = x  + 3*(*i0-1);
  pix = ix + 3*(*i0-1);
  pu  = u  + 3*(*i0-1);
  pq  = q  +   (*i0-1);
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 2 add dummy particles to the end
  if ( np_total % 4 ) {
	
	for( i = 0; i < 4 - np_total%4; i ++ ) {
	  
	  pix[ 3 * (np_total + i)     ] = 1;
	  pix[ 3 * (np_total + i) + 1 ] = 1;
	  pix[ 3 * (np_total + i) + 2 ] = 1;
	   
	   px[ 3 * (np_total + i)     ] = 0.;
	   px[ 3 * (np_total + i) + 1 ] = 0.;
	   px[ 3 * (np_total + i) + 2 ] = 0.;
	   
	   pu[ 3 * (np_total + i)     ] = 0.;
	   pu[ 3 * (np_total + i) + 1 ] = 0.;
	   pu[ 3 * (np_total + i) + 2 ] = 0.;
	   
	   pq[     (np_total + i)     ] = 0.;
	}
    
    np_total += ( 4 - np_total%4 );
    if ( np_total % 4 ) {
       printf("(*error*) Number of particles to push is still not a multiple of 2!\n");
       exit(1);
    }
    
  }
  
  if ( p_cache_size_3D % 4 ) {
       printf("(*error*) p_cache_size_3D is not a multiple of 2!\n");
       exit(1);
  }
  
  for ( k = 0; k < np_total; k += p_cache_size_3D ) {
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_3D ) ? np_total - k : p_cache_size_3D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
	vpbuf.p = (double *) (&vpbuf.buf[0]);
		
	// Push all particles in group
	for( i = 0; i < np; i+=4 ) {
	  vector vx0, vy0, vz0;
	  DECLARE_ALIGNED_32( int vix[4] );
	  DECLARE_ALIGNED_32( int viy[4] );
	  DECLARE_ALIGNED_32( int viz[4] );
  
	  vector ve1, ve2, ve3;
	  vector vb1, vb2, vb3;
  
	  vector vu1, vu2, vu3;

	  // Load particle positions 
	  VEC_LD4v3( vx0, vy0, vz0, px );

	  for( int l = 0; l < 4; l++ ) {
	    vix[l] = pix[3*l    ];
	    viy[l] = pix[3*l + 1];
	    viz[l] = pix[3*l + 2];
	  }

	  // Interpolate fields
	  {
		 vector  vx0h, vy0h, vz0h;
		 int idxi[4], idxj[4], idxk[4], idxih[4], idxjh[4], idxkh[4];
		 vector vwx[NP], vwy[NP], vwz[NP], vwxh[NP], vwyh[NP], vwzh[NP];

	     int l;
	     
	     for( l=0; l<4; l++) {
	       idxi[l] =      3 * vix[l];
	       idxj[l] = deltaY * viy[l];
	       idxk[l] = deltaZ * viz[l];
	     }
	     	 
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy );
		 SPLINE( vz0, vwz );
	 
		 vget_hp_near( vx0, idxi, vx0h, idxih, 3 ); 
		 vget_hp_near( vy0, idxj, vy0h, idxjh, deltaY ); 
		 vget_hp_near( vz0, idxk, vz0h, idxkh, deltaZ ); 
	 
		 SPLINE( vx0h, vwxh );
		 SPLINE( vy0h, vwyh );
		 SPLINE( vz0h, vwzh );
      
         // Interpolate E field
         {
			int k1, k2, k3;
			t_real *fld1[4], *fld2[4], *fld3[4];
			vector const zero = vec_splats( 0.0 );
						  
			// get pointers to fields in particle cells (i-1, j-1, k-1)
			for( k1 = 0; k1 < 4; k1++) {
			   fld1[k1] = e1 + idxih[k1] +  idxj[k1] +  idxk[k1]; 
			   fld2[k1] = e2 +  idxi[k1] + idxjh[k1] +  idxk[k1];
			   fld3[k1] = e3 +  idxi[k1] +  idxj[k1] + idxkh[k1];
            }
            			
			ve1 = zero;
			ve2 = zero;
			ve3 = zero;
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register vector f1plane, f2plane, f3plane;
			   f1plane = zero;
			   f2plane = zero;
			   f3plane = zero;
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register vector f1line, f2line, f3line;
				 f1line = zero;
				 f2line = zero;
				 f3line = zero;

				 register vector f1point[NP], f2point[NP], f3point[NP];
				 for ( k1 = 0; k1 < NP; k1++ ) {
					int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point[k1] = LOADFLD4( fld1, shift ); 
					f2point[k1] = LOADFLD4( fld2, shift ); 
					f3point[k1] = LOADFLD4( fld3, shift ); 
				 }
				 
				 for ( k1 = 0; k1 < NP; k1 ++ ) {
					f1line = vec_madd( vwxh[k1], f1point[k1], f1line );
				    f2line = vec_madd(  vwx[k1], f2point[k1], f2line );
				    f3line = vec_madd(  vwx[k1], f3point[k1], f3line );
				 }
				 
				 f1plane = vec_madd(  vwy[k2], f1line, f1plane );
				 f2plane = vec_madd( vwyh[k2], f2line, f2plane );
				 f3plane = vec_madd(  vwy[k2], f3line, f3plane );
			   }
			   
			   ve1 = vec_madd(  vwz[k3], f1plane, ve1 );
			   ve2 = vec_madd(  vwz[k3], f2plane, ve2 );
			   ve3 = vec_madd( vwzh[k3], f3plane, ve3 );
			}
         
         }		 

         // Interpolate B field
         {
			int k1, k2, k3;
			t_real *fld1[4], *fld2[4], *fld3[4];
			vector const zero = vec_splats( 0.0 );
			  
			// get pointers to fields in particle cells (i-1, j-1, k-1)
			for( k1 = 0; k1 < 4; k1++) {
			   fld1[k1] = b1 +  idxi[k1] + idxjh[k1] + idxkh[k1]; 
			   fld2[k1] = b2 + idxih[k1] +  idxj[k1] + idxkh[k1];
			   fld3[k1] = b3 + idxih[k1] + idxjh[k1] +  idxk[k1];
			}
			
			vb1 = zero;
			vb2 = zero;
			vb3 = zero;
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register vector f1plane, f2plane, f3plane;
			   f1plane = zero;
			   f2plane = zero;
			   f3plane = zero;
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register vector f1line, f2line, f3line;
				 f1line = zero;
				 f2line = zero;
				 f3line = zero;

				 register vector f1point[NP], f2point[NP], f3point[NP];
				 for ( k1 = 0; k1 < NP; k1++ ) {
					int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point[k1] = LOADFLD4( fld1, shift ); 
					f2point[k1] = LOADFLD4( fld2, shift ); 
					f3point[k1] = LOADFLD4( fld3, shift ); 
				 }
				 
				 for ( k1 = 0; k1 < NP; k1 ++ ) {
					f1line = vec_madd(  vwx[k1], f1point[k1], f1line );
					f2line = vec_madd( vwxh[k1], f2point[k1], f2line );
					f3line = vec_madd( vwxh[k1], f3point[k1], f3line );
				 }
				 
				 f1plane = vec_madd( vwyh[k2], f1line, f1plane  );
				 f2plane = vec_madd(  vwy[k2], f2line, f2plane  );
				 f3plane = vec_madd( vwyh[k2], f3line, f3plane  );
			   }
			   
			   vb1 = vec_madd( vwzh[k3], f1plane, vb1 );
			   vb2 = vec_madd( vwzh[k3], f2plane, vb2 );
			   vb3 = vec_madd(  vwz[k3], f3plane, vb3 );
			}
         
         }		 
      
	  }	  

	  // ---------- advance momenta
	  // results are stored in vu1, vu2, vu3 (temp) and u (permanent)
	  
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, pu, vu1, vu2, vu3 );
      VEC_ST4v3(pu, vu1, vu2, vu3);

	  // advance u 4 particles * 3 components
	  pu  += 12; 
   
	  // ---------- advance position and get velocities
	  {
		register vector vrg, vq;
		register vector vtrx, vtry, vtrz;
		register vector vx1, vy1, vz1;
		
	    DECLARE_ALIGNED_32( int cross[4] );
		DECLARE_ALIGNED_32( int dix[4] );
		DECLARE_ALIGNED_32( int diy[4] );
		DECLARE_ALIGNED_32( int diz[4] );

  
		VRGAMMA( vrg, vu1, vu2, vu3 );
	    
	    vq = vec_lda( 0, pq );
	    
		// get velocities
		vu1 = vec_mul( vu1, vrg );
		vu2 = vec_mul( vu2, vrg );
		vu3 = vec_mul( vu3, vrg );
			
		vx1 = vec_madd( vdt_dx, vu1, vx0 );
		vy1 = vec_madd( vdt_dy, vu2, vy0 );
		vz1 = vec_madd( vdt_dz, vu3, vz0 );
	  
		vtrx = vntrim(vx1);
		vtry = vntrim(vy1);
		vtrz = vntrim(vz1);
				
		STORE4VP3D( vpbuf.p, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz );

		vx1   = vec_sub( vx1, vtrx );
		vy1   = vec_sub( vy1, vtry );
		vz1   = vec_sub( vz1, vtrz );
        VEC_ST4v3( px, vx1, vy1, vz1 );

		vtrx = vec_ctiw( vtrx );
		vtry = vec_ctiw( vtry );
		vtrz = vec_ctiw( vtrz );
		vec_sta( vtrx, 0, dix );
		vec_sta( vtry, 0, diy );
		vec_sta( vtrz, 0, diz );
				
        for( int l = 0; l < 4; l++) {
          // store new particle cell indexes
		  pix[3*l    ] = vix[l] + dix[l];
		  pix[3*l + 1] = viy[l] + diy[l];
		  pix[3*l + 2] = viz[l] + diz[l];
        
		  // calculate crossings
		  cross[l] = (( dix[l] )?1:0) + (( diy[l] )?2:0) + (( diz[l] )?4:0);
        }
        
		// ---------- split trajectories for current deposition
  
		vsplit3D( &vpbuf, cross, dix, diy, diz );

	  }

	  // advance x, ix 4 particles * 3 components
	  px  += 12;  
	  pix += 12;
	  // advance q 4 particles * 1 component
	  pq  += 4;
  
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

