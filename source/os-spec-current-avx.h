#ifndef _OS_SPEC_CURRENT_H
#define _OS_SPEC_CURRENT_H

#include "vector-avx.h"
#include "os-spec-push-avx.h"

void vwl_s1( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] );
void vwl_s2( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] );
void vwl_s3( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] );
void vwl_s4( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] );

/* Current deposition routines for x86 SSE */

void vdepcurrent_2d_s1(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_vpbuf2D * const part);
void vdepcurrent_2d_s2(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_vpbuf2D * const part);
void vdepcurrent_2d_s3(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_vpbuf2D * const part);
void vdepcurrent_2d_s4(float * const current, int const * const size, int const * const offset, 
                       float * const norm, t_vpbuf2D * const part);

void vdepcurrent_3d_s1(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_vpbuf3D * const part);
void vdepcurrent_3d_s2(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_vpbuf3D * const part);
void vdepcurrent_3d_s3(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_vpbuf3D * const part);
void vdepcurrent_3d_s4(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_vpbuf3D * const part);


#endif
