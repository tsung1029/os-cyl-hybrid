/*****************************************************************************************

Current deposition, SIMD optimized version

*****************************************************************************************/

#ifdef SIMD



/**************************************** SSE *******************************************/

#ifdef SIMD_SSE

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-current-sse.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-current-ssed.c"

#endif 

#endif

/**************************************** SSE *******************************************/

#ifdef SIMD_AVX

#define HAS_SIMD_CODE 1

#if defined( PRECISION_SINGLE )

#include "os-spec-current-avx.c"

#elif defined( PRECISION_DOUBLE )

#include "os-spec-current-avxd.c"

#endif 

#endif

/************************************ BlueGene/P ****************************************/

#ifdef SIMD_BGP
#define HAS_SIMD_CODE 1

#include "os-spec-current-bgp.c"

#endif 

/************************************ BlueGene/Q ****************************************/

#ifdef SIMD_BGQ
#define HAS_SIMD_CODE 1

#include "os-spec-current-bgq.c"

#endif 

/*********************************** Sanity Check ***************************************/

#ifndef HAS_SIMD_CODE

#error When using SIMD code define SIMD = SSE | AVX | BGP | BGQ in the config file

#endif

#endif


