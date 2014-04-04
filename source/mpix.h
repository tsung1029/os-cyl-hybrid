/*

Dummy mpix.h file to test BlueGene/Q code

*/

#ifndef __include_mpix_h__
#define __include_mpix_h__


#include <stdio.h>

static const int torus_ndims  = 5;
static const int torus_part[] = { 2, 2, 2, 2, 2, 16 }; 
static const int torus_size = 512;

#define MPIX_TORUS_MAX_DIMS 5 /* This is the maximum physical size of the torus */

 typedef struct
 {
/* These fields will be used on all platforms. */
   unsigned prank;    /**< Physical rank of the node (irrespective of mapping) */
   unsigned psize;    /**< Size of the partition (irrespective of mapping) */
   unsigned ppn;      /**< Processes per node */
   unsigned coreID;   /**< Process id; values monotonically increase from 0..63 */

   unsigned clockMHz; /**< Frequency in MegaHertz */
   unsigned memSize;  /**< Size of the core memory in MB */

/* These fields are only set on torus platforms (i.e. Blue Gene) */
   unsigned torus_dimension;              /**< Actual dimension for the torus */
   unsigned Size[MPIX_TORUS_MAX_DIMS];    /**< Max coordinates on the torus */
   unsigned Coords[MPIX_TORUS_MAX_DIMS];  /**< This node's coordinates */
   unsigned isTorus[MPIX_TORUS_MAX_DIMS]; /**< Do we have wraparound links? */

/* These fields are only set on systems using Blue Gene IO psets. */
   unsigned rankInPset;
   unsigned sizeOfPset;
   unsigned idOfPset;
 } MPIX_Hardware_t;
 
static int MPIX_Hardware(MPIX_Hardware_t *hw)
 {
    return 0;
 }
 
static int MPIX_Torus_ndims(int *numdim)
{
  *numdim = torus_ndims;
}

static int MPIX_Rank2torus(int rank, int *coords) {
    
   if ( rank < 0 || rank >= torus_size ) {
	  fprintf( stderr, "(*error*) Invalid rank : %d, should be in [0,%d[", 
						rank, torus_size );
	  return -1;
   }
   
   for( int i = 0; i < torus_ndims+1; i++ ) {
	 coords[i] = rank % torus_part[i];
	 rank      = rank / torus_part[i];
   }
   
   return 0;
}
 
 
static int MPIX_Torus2rank(int *coords, int *rank) {
   *rank = 0;
   int stride = 1;
   
   for( int i = 0; i < torus_ndims+1; i++ ) {
      if ( coords[i] < 0 || coords[i] >= torus_part[i] ) {
		 fprintf( stderr, "(*error*) Invalid coordinate[%d] : %d, should be in [0,%d[", 
						   i, coords[i], torus_part[i] );
		  return -1;
      }
      
      rank += coords[i] * stride;
      stride *= torus_part[i];
   }
   return 0;
}


#endif
