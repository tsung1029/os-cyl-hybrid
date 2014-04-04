#ifndef _OS_SPEC_CURRENT_H
#define _OS_SPEC_CURRENT_H

#include "vector-bgp.h"
#include "os-spec-push-bgp.h"

inline void vwl_s1( vector const vqn, vector const vx0, vector const vx1, vector vwl[] );
inline void vwl_s2( vector const vqn, vector const vx0, vector const vx1, vector vwl[] );
inline void vwl_s3( vector const vqn, vector const vx0, vector const vx1, vector vwl[] );
inline void vwl_s4( vector const vqn, vector const vx0, vector const vx1, vector vwl[] );


/* Current deposition routines for BG/P */

inline void vdepcurrent_2d_s1(t_real * const current, int const * const size, int const * const offset, 
                       double * const norm, t_vpbuf2D * const part);
inline void vdepcurrent_2d_s2(t_real * const current, int const * const size, int const * const offset, 
                       double * const norm, t_vpbuf2D * const part);
inline void vdepcurrent_2d_s3(t_real * const current, int const * const size, int const * const offset, 
                       double * const norm, t_vpbuf2D * const part);
inline void vdepcurrent_2d_s4(t_real * const current, int const * const size, int const * const offset, 
                       double * const norm, t_vpbuf2D * const part);

inline void vdepcurrent_3d_s1(t_real * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_vpbuf3D * const part);
inline void vdepcurrent_3d_s2(t_real * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_vpbuf3D * const part);
inline void vdepcurrent_3d_s3(t_real * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_vpbuf3D * const part);
inline void vdepcurrent_3d_s4(t_real * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_vpbuf3D * const part);

#endif
