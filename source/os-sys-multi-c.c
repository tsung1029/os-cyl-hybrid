/*

 $URL: svn+ssh://exppmaster/svn_repositories/osiris/trunk/source/os-sys-macosx-c.c $
 $Id: os-sys-macosx-c.c 55 2005-11-30 17:31:32Z zamb $

 Interfacing routines for POSIX system calls
 (should work under any POSIX compliant system)

*/

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <sys/resource.h>

/* Macros for fortran functions */
#include "fortran.h"

/* Force stdout and stderr to be line buffered */

void FNAME(linebuf_stdio_f)
(void)
{
  setvbuf(stdout, (char *)NULL, _IOLBF, 0);
  setvbuf(stderr, (char *)NULL, _IOLBF, 0);
}

/* Remove a file */

void FNAME(remove_f)
(const char *path, int *ierr)
{
 if ((*ierr = remove(path))) *ierr = errno;
}


/* Create a directory */

void FNAME(mkdir_f)
(const char *path, int *ierr)
{
 char uppath[256], *p;
    
 if (mkdir(path,S_IRWXU | (S_IRGRP | S_IXGRP ) | (S_IROTH | S_IXOTH) ))
   switch (errno) {
     case ENOENT : /* A component of the path does not exist */
        
        /* get upper path */
        
        strcpy(uppath, path);
        p = uppath + strlen(uppath);
        while(*p!='/') p--;
        *p=0;
        
        /* recursively build the path */
        
        FNAME(mkdir_f) (uppath, ierr);
        if (!*ierr) FNAME(mkdir_f) (path, ierr);
     
     case EEXIST : /* if directory already exists ignore the error */
        *ierr = 0;
        break;
     default: *ierr = errno;
   }
 else *ierr = 0;
 
}

/* change working a directory */

void FNAME(chdir_f)
(const char *path, int *ierr)
{
  if ((*ierr = chdir(path))) *ierr=errno;
}

/* get current working a directory */

void FNAME(getcwd_f)
(char *path, int *len, int *ierr)
{
  char *buf = (char*) malloc( (*len+1) );
  
  if (getcwd(buf, *len)) {
     strcpy(path, buf);
     *ierr = 0;
   } else {
    strcpy(path, "\0");
    *ierr=errno;
   }
  
  free(buf);
}

/* Create a symbolic link */

void FNAME(symlink_f)
(const char *name1, const char *name2, int *ierr)
{
  if ((*ierr = symlink(name1,name2))) *ierr=errno;
}

/* Get the error string */

void FNAME(strerror_f)
(const int *errnum, char *err)
{
 strcpy(err, strerror(*errnum));
}

/* Get the environment variable value */

void FNAME(getenv_f)
(const char *name, char *value, int *ierr)
{
 char *lvalue;
  
 lvalue = (char *) getenv(name);
  
 if (lvalue == NULL) *ierr = -1;
 else 
 {
   *ierr = 0;
   strcpy(value,lvalue);
 }
 
}

/* Get the hostname */

void FNAME(gethostname_f)
(char *hostname, int *ierr)
{
 char lhostname[256];
  
 if (gethostname(lhostname, 255)) {
   strcpy(hostname,"\0");
   *ierr = errno; 
 } else {
   strcpy(hostname,lhostname);
   *ierr = 0;
 }
}

/* Set the umask */
void FNAME(umask_f)
(int *numask)
{
  umask( *numask );
}

/* check if floating point is infinity */

void FNAME(isinf_f)
( double *x, int *res )
{
  *res = isinf( *x);
}

/* check if floating point is nan */

void FNAME(isnan_double_f)
( double *x, int *res )
{
  *res = isnan( *x);
}

void FNAME(isnan_single_f)
( float *x, int *res )
{
  *res = isnan( *x);
}


/* memcpy wrappers */
void FNAME(memcpy_int)
( void *s1, const void *s2, const int *n )
{
  size_t csize = (*n) * sizeof(int);
  memcpy( s1, s2, csize );
}

void FNAME(memcpy_r4)
( void *s1, const void *s2, const int *n )
{
  
  size_t csize = (*n) * sizeof(float);
  memcpy( s1, s2, csize );
}

void FNAME(memcpy_r8)
( void *s1, const void *s2, const int *n )
{
  size_t csize = (*n) * sizeof(double);
  memcpy( s1, s2, csize );
}

/* setnan wrappers */
void FNAME(setnan_r4)
( float s[], const int *n )
{
  int i;
  float v;
  
  v = nanf("");
  
  for( i = 0; i < *n; i++ ) s[i] = v;
}

void FNAME(setnan_r8)
( double s[], const int *n )
{
  int i;
  double v;
  
  v = nan("");
  
  for( i = 0; i < *n; i++ ) s[i] = v;
}

/* zero_buffer wrappers */
void FNAME(zero_buffer_r4)
( void *b, const int *n )
{
  size_t len = (*n) * sizeof(float);
  memset( b, 0, len );
}

void FNAME(zero_buffer_r8)
( void *b, const int *n )
{
  size_t len = (*n) * sizeof(double);
  memset( b, 0, len );
}


/****************************************************************************** 
   
   System dependent high resolution timers.
   
   Note: There is a significant overhead on converting the time from ticks
         to seconds, so for high precision measurements a tick based approach
         is preferrable.
         
         Resolution is defined in terms of the minimum time difference between
         calls to the system clock, without accounting the time for converting
         the ticks to seconds.
   
*******************************************************************************/

#ifdef __AIX_TIMER__
#define HAVE_TIMER

/* 
  Use high resolution timers from the Standard C Library (libc.a) on IBM AIX 
  systems running on PowerPC architecture. The resolution is system dependent
  and is on the order of tens of nanoseconds (~ 70 ns on a JS21 PPC5+ blade).    
*/

#include <sys/time.h>

int tb_flag;

void FNAME(timer_init)
(void)
{
  /* store real time flag for future use */
  
  timebasestruct_t tb;
  read_real_time(&tb, TIMEBASE_SZ);
  tb_flag = tb.flag;
}

uint64_t FNAME(timer_ticks)
(void)
{
  timebasestruct_t tb;
  read_real_time(&tb, TIMEBASE_SZ);
 
  return ((uint64_t)tb.tb_high)<<32 | (uint64_t)tb.tb_low;

}

double FNAME(timer_interval_seconds)
(uint64_t *start, uint64_t *end)
{
  int sec, nsec;
  timebasestruct_t tb_start, tb_end;

  /* Convert start/end values to real time */
  
  tb_start.tb_high = (unsigned int)((*start) >> 32);
  tb_start.tb_low  = (unsigned int)((*start) & 0xFFFFFFFF);
  tb_start.flag = tb_flag;
  time_base_to_time(&tb_start, TIMEBASE_SZ);

  tb_end.tb_high = (unsigned int)((*end) >> 32);
  tb_end.tb_low  = (unsigned int)((*end) & 0xFFFFFFFF);
  tb_end.flag = tb_flag;
  time_base_to_time(&tb_end, TIMEBASE_SZ);
    
  sec  = tb_end.tb_high - tb_start.tb_high;
  nsec = tb_end.tb_low  - tb_start.tb_low;
  if ( nsec < 0 ) {
    sec --;
    nsec += 1000000000;
  } 
  return (double) sec + 1.0e-9 * ((double) nsec);
}

double FNAME(timer_cpu_seconds)
( void )
{
    timebasestruct_t tb;
    double wtime;
    
    read_real_time(&tb, TIMEBASE_SZ);
    time_base_to_time(&tb, TIMEBASE_SZ);

    wtime  =  (double)tb.tb_low * 1.0e-9;
    wtime +=  (double)tb.tb_high;
    return wtime;
}

double FNAME(timer_resolution)
( void )
{
  timebasestruct_t tb1,tb2;
  int res;
   
  read_real_time(&tb1, TIMEBASE_SZ);
  do {
       read_real_time(&tb2, TIMEBASE_SZ);
  } while (tb1.tb_low == tb2.tb_low);
  
  time_base_to_time(&tb1, TIMEBASE_SZ);
  time_base_to_time(&tb2, TIMEBASE_SZ);
  
  res = tb2.tb_low - tb1.tb_low;
  if ( res < 0 ) res += 1000000000;
    
  return (double)(res) * 1e-9;

}

#endif /* __AIX_TIMER__ */


#ifdef __MACH_TIMER__
#define HAVE_TIMER

/* 
  Use high resolution timers from the Mach kernel. The resolution is system
  dependent and is on the order of tens of nanoseconds.
  
  - Resolution in a G4 800Mhz: 60.22 ns
   
  There is a (unlikely) possibility of overflow that is not checked. 
  (this could be fixed by recording mach_absolute_time() when timer_init 
  is called and then checking for overflow on timer_cpu_seconds)
  
*/

#include <mach/mach_time.h>

static mach_timebase_info_data_t    sTimebaseInfo;
static double tick_conv;

void FNAME(timer_init)
(void)
{
  (void) mach_timebase_info(&sTimebaseInfo);
  tick_conv = (double) sTimebaseInfo.numer / (double) sTimebaseInfo.denom * 1e-9;
}

uint64_t FNAME(timer_ticks)
(void)
{
  return mach_absolute_time();
}

double FNAME(timer_interval_seconds)
(uint64_t *start, uint64_t *end)
{
  /* we should check for overflow */
  
  return (*end - *start) * tick_conv;
}

double FNAME(timer_cpu_seconds)
( void )
{
  return (double)mach_absolute_time() * tick_conv;
  
}

double FNAME(timer_resolution)
( void )
{
  uint64_t        start, end;
  start = mach_absolute_time();
  end = mach_absolute_time();
  
  return (end-start)*tick_conv;
}

#endif /* __MACH_TIMER__ */

#ifdef __HRTIME_TIMER__
#define HAVE_TIMER

/* 
  Use high resolution timers based on the gethrtime function,
  available on most linux distros.  
*/

#include <sys/time.h>

void FNAME(timer_init)
(void)
{
 /* No initialization required */
}

hrtime_t FNAME(timer_ticks)
(void)
{
  return gethrtime();
}

double FNAME(timer_interval_seconds)
(hrtime_t *start, hrtime_t *end)
{
  /* we should check for overflow */
  
  return (*end - *start) * 1.0e-9;
}

double FNAME(timer_cpu_seconds)
( void )
{
  return (double)gethrtime() * 1.0e-9;
  
}

double FNAME(timer_resolution)
( void )
{
  hrtime_t        start, end;
  start = gethrtime();
  end = gethrtime();
  
  return (end-start)*1.0e-9;
}

#endif /* __HRTIME_TIMER__ */

/***************************************************************************************************
  BlueGene/P specific code (includes timers)
***************************************************************************************************/

#ifdef __bgp__         

#include "mpi.h"
#include "mpix.h"

/* returns BlueGene hardware topology has an MPI topology ordered has [T,Z,Y,X] */

void FNAME(bgp_mpix_cart_comm_create_f)
(MPI_Fint *comm, int *ierr)
{
  MPI_Comm cart_comm_bg;

  *ierr = MPIX_Cart_comm_create (& cart_comm_bg);
  if (*ierr != MPI_SUCCESS ) cart_comm_bg = MPI_COMM_NULL;

  *comm = MPI_Comm_c2f( cart_comm_bg );
}


/* returns BlueGene memory per core in MBytes */

#include <sys/resource.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#include <spi/kernel_interface.h>
static _BGP_Personality_t mybgp;

void FNAME( bgp_core_memory_f )
( long *bg_coreMB )
{
  unsigned procMB, coreMB;
  Kernel_GetPersonality(&mybgp, sizeof(_BGP_Personality_t));
  procMB = BGP_Personality_DDRSizeMB(&mybgp);
  coreMB = procMB/Kernel_ProcessCount();
  *bg_coreMB = coreMB;
}

#include <unistd.h>

/* returns BlueGene process memory in kbytes */

void FNAME(bgp_used_memory_f)
( long *used_kb )
{
  struct rusage usage;
  if ( getrusage(RUSAGE_SELF, &usage) ) *used_kb = 0;
  else *used_kb = usage.ru_maxrss;

}


/* BG specific timers */
#define HAVE_TIMER

static double tick_conv;

void FNAME(timer_init)
(void)
{
   unsigned clockMHz;
   Kernel_GetPersonality(&mybgp, sizeof(_BGP_Personality_t));
   clockMHz = BGP_Personality_clockMHz(&mybgp);
   tick_conv = 1.0e-6/(double)clockMHz;
}

uint64_t FNAME(timer_ticks)
( void )
{
   return _bgp_GetTimeBase();
}

double FNAME(timer_interval_seconds)
(uint64_t  *start, uint64_t *end)
{
  /* we should check for overflow */
  
  return (*end - *start) * tick_conv;
}

double FNAME(timer_resolution)
( void )
{
  uint64_t        start, end;
  start = _bgp_GetTimeBase();
  end   = _bgp_GetTimeBase();
  
  return (end-start)*tick_conv;
}

double FNAME(timer_cpu_seconds)
( void )
{
    return (double)_bgp_GetTimeBase() * tick_conv;
}


#endif


/*****************************************************************************************
  BlueGene/Q specific code
*****************************************************************************************/

#ifdef __bgq__

/* MPI extensions for BlueGene/Q */
#include "mpi.h"
#include "mpix.h"

/*
 HWPART_MAX_DIMS defines the maximum number of dimensions. The algorithm include the
 T (core+thread ID) coordinate and it also allows for the additional splitting of the one
 dimension (A through T) into 2 dimensions, so this should be set to TORUS_DIMS + 2
*/
#define HWPART_MAX_DIMS MPIX_TORUS_MAX_DIMS + 2

/*
 SIMPART_MAX_DIMS defines the maximum number of dimensions for the simulation partition.
*/
#define SIMPART_MAX_DIMS 3

/* Structure describing the BlueGene/Q network partition */
typedef struct {

  // Number of nodes on partition
  int nnodes;

  // Number of dimensions of network torus
  int hwpart_torusdims;
  // Size of each hardware partition dimension
  int hwpart[ HWPART_MAX_DIMS ];
  // Do we have wraparound links in this dimension
  int hwpart_istorus[ HWPART_MAX_DIMS ];
  // Location on global torus
  int hwpart_coords[ HWPART_MAX_DIMS ];
  // Global MPI rank
  int mpi_rank;

  // Coordinate to split (if any)
  int hwpart_splitdim; 
  // Number of dimensions of mapping partition (includes splitting)
  int mappart_ndims;
  // Size of each hardware partition dimension (includes splitting)
  int mappart[ HWPART_MAX_DIMS ];
  
  // Number of dimensions of simulation partition
  int simpart_ndims;
  // Size of each simulation partition
  int simpart[ SIMPART_MAX_DIMS ];
  
  // Number of hardware dimensions for each simulation dimension
  int nmap[ SIMPART_MAX_DIMS ];
  // Mapping between simulation and hardware partitions
  int map[ SIMPART_MAX_DIMS ][ HWPART_MAX_DIMS ];
   
} bgq_partition_t;

static bgq_partition_t bgq_part;


/*****************************************************************************************
* Find mapping between simulation and hardware topology
*****************************************************************************************/
int bgq_findmap(  int dim,     int ndims,   int sim_part[], 
                  int npart,   int part[],  int part_map[],
                  int nmap[],  int map[][ HWPART_MAX_DIMS ] )
{
  
  for( int i = 0; i < 1 << npart; i ++ ) {
      int a = 1;
      for( int k = 0; k < npart; k++ ) if ( i & 1 << k ) a*= part[k];
      
      if ( a == sim_part[dim] ) {
         
         // Store solution
         nmap[dim] = 0;
         for( int k = 0; k < npart; k++ ) if ( i & 1 << k ) {
            map[dim][nmap[dim]] = part_map[k];
            nmap[dim]++;
         }
         
         if ( dim == ndims-1 ) {
            // If last dimension then we're done!
            return 1;
         } else {
            // Otherwise find mapping for remaining dimensions
             int new_part_map[ HWPART_MAX_DIMS ];
             int new_part[ HWPART_MAX_DIMS ];
             int new_npart;
            
            // Build a hardware partition with the remaining hw dimensions
            new_npart = 0;
            for( int k = 0; k < npart; k++ ) if ( !(i & 1 << k) ) {
               new_part_map[ new_npart ] = part_map[k];
               new_part[ new_npart ] = part[k];
               new_npart++;
            }
            
            // Find a mapping for the remaining dimensions; if successful return 1
            // otherwise keep searching
            if ( bgq_findmap( dim+1, ndims, sim_part,
                              new_npart, new_part, new_part_map, 
                              nmap, map ) > 0 ) return 1;
         }
      }
   }
   
   // All combinations were tested and no mapping was found
   return 0;
}
             

/*****************************************************************************************
* Get simulation coordinates from hardware coordinates
*****************************************************************************************/
void bgq_hw2sim( int hwcoords[], int simcoords[] )
{
   int mappart[ HWPART_MAX_DIMS ];
   int mapcoords[ HWPART_MAX_DIMS ];
   int weight[ HWPART_MAX_DIMS ];
   
   // Account for splitting of a hw dimension
   int _hwcoords[ HWPART_MAX_DIMS ];

   // If no split just copy the torus coordinates 
   if ( bgq_part.hwpart_splitdim < 0 ) {
      for( int i=0; i < bgq_part.hwpart_torusdims+1; i++ ) _hwcoords[i] = hwcoords[i];  
   } else {
      // Otherwise split the corresponding torus coordinate
      int dim = bgq_part.hwpart_splitdim;
 	  for( int i=0; i < dim ; i++ ) _hwcoords[i] = hwcoords[i];
	  _hwcoords[dim]   = hwcoords[dim] / bgq_part.mappart[ dim+1 ];
	  _hwcoords[dim+1] = hwcoords[dim] % bgq_part.mappart[ dim+1 ];
	  for( int i=dim+2; i < bgq_part.mappart_ndims; i++ ) _hwcoords[i] = hwcoords[i-1];
   }

   for( int dim = 0; dim < bgq_part.simpart_ndims; dim++ ) {
      
      // Get hw partition/coordinates that map to the current dimension
      int nmap = bgq_part.nmap[dim];
      int size = 1;
      for ( int i = 0; i < nmap; i++ ) {
         mappart[i] = bgq_part.mappart[ bgq_part.map[dim][i] ];
         mapcoords[i] = _hwcoords[ bgq_part.map[dim][i] ];
         size *= mappart[i];
      }
      
      // Get weight of each dimension
      weight[0] = size / mappart[0];
      for ( int i = 1; i < nmap; i++ ) weight[i] = weight[i-1] / mappart[i];
      
      // Get dim rank from map coordinates
      int rank = 0;
      int ff = 1;
      for ( int i = 0; i < nmap; i++ ) {
        rank += ((ff) ? mapcoords[i] : mappart[i] - mapcoords[i] - 1 ) * weight[i];
        ff = ! ( (rank / weight[i]) % 2 );
      }
      
      // Store result
      simcoords[dim] = rank;
   }
  
}

/*****************************************************************************************
* Get hardware coordinates from simulation coordinates
*****************************************************************************************/
void bgq_sim2hw( int simcoords[], int hwcoords[] )
{
   int mappart[ HWPART_MAX_DIMS ];
   int mapcoords[ HWPART_MAX_DIMS ];
   int weight[ HWPART_MAX_DIMS ];
   int _hwcoords[ HWPART_MAX_DIMS ];
         
   for( int dim = 0; dim < bgq_part.simpart_ndims; dim++ ) {
      
      // Get hw partition/coordinates that map to the current dimension
      int nmap = bgq_part.nmap[dim];
      int size = 1;
      for ( int i = 0; i < nmap; i++ ) {
         mappart[i] = bgq_part.mappart[ bgq_part.map[dim][i] ];
         size *= mappart[i];
      }
      
      // Get weight of each dimension
      weight[0] = size / mappart[0];
      for ( int i = 1; i < nmap; i++ ) weight[i] = weight[i-1] / mappart[i];
      
      // Get map coordinates from dim rank 
      int rank = simcoords[dim];
      int ff = 1;
      for ( int i = 0; i < nmap; i++ ) {
        _hwcoords[ bgq_part.map[dim][i] ] = (ff) ? rank/weight[i] : 
                                                    mappart[i] - rank/weight[i] - 1;
        ff = ! ( (rank / weight[i]) % 2 );
        rank %= weight[i];
      }
   }

   // If no split just copy the torus coordinates 
   if ( bgq_part.hwpart_splitdim < 0 ) {
      for( int i=0; i < bgq_part.hwpart_torusdims+1; i++ ) hwcoords[i] = _hwcoords[i];  
   } else {
      // Otherwise merge split the corresponding coordinate
      int dim = bgq_part.hwpart_splitdim;
 	  for( int i=0; i < dim ; i++ ) hwcoords[i] = _hwcoords[i];
	  hwcoords[dim]   = _hwcoords[dim] * bgq_part.mappart[ dim+1 ] + _hwcoords[dim+1];
	  for( int i=dim+1; i < bgq_part.hwpart_torusdims+1; i++ ) hwcoords[i] = _hwcoords[i+1];
   }
   
}

/*****************************************************************************************
* Get simulation coordinates from global MPI rank
*****************************************************************************************/
void FNAME(bgq_rank2sim)
( int *rank, int *simcoords, int *ierr ) {
   int hwcoords[MPIX_TORUS_MAX_DIMS + 1 ];
   
   // Validate mpi rank
   if ( *rank < 0 || *rank >= bgq_part.nnodes ) {
	  fprintf( stderr, "(*error*) in bgq_rank2sim() invalid rank = %d,"
					   " should be in the range [%d,%d]\n",      
					   *rank, 0, bgq_part.nnodes-1);
	  *ierr = -1;
	  return;
   }
   
   // Convert mpi rank to torus coordinates
   *ierr = MPIX_Rank2torus( *rank, hwcoords );
   
   // If successful convert hardware coordinates to simulation coordinates
   if ( ! *ierr ) bgq_hw2sim( hwcoords, simcoords );
}

/*****************************************************************************************
* Get global MPI rank from simulation coordinates
*****************************************************************************************/
void FNAME(bgq_sim2rank)
( int* simcoords, int *rank, int *ierr ) {
   int hwcoords[MPIX_TORUS_MAX_DIMS + 1 ];
   
   // Validate simcoords
   for( int i=0; i< bgq_part.simpart_ndims; i++) {
     if ( simcoords[i] < 0 || simcoords[i] >=  bgq_part.simpart[i] ) {
        fprintf( stderr, "(*error*) in bgq_sim2rank() invalid simcoords[%d] = %d,"
                         " should be in the range [%d,%d]\n",      
                         i, simcoords[i], 0, bgq_part.simpart[i]-1);
        *ierr = -1;
        return;
     }
   }
   
   // Convert simulation coordinate to hardware coordinates
   bgq_sim2hw( simcoords, hwcoords );
   
   // Debug
   {
      int _simcoords[3];
      bgq_hw2sim(hwcoords, _simcoords);

       for(int i=0; i<bgq_part.simpart_ndims; i++) {
        if ( simcoords[i] != _simcoords[i] ) { 
		  fprintf( stderr, "(*error*) in bgq_sim2rank() Problem in bgq_hw2sim / sim2hw "
		                   "simulation coordinates -> dim %d, in: %d, out:%d\n",
		                   i, simcoords[i], _simcoords[i] );
		  *ierr = -2;
		  return;
        }
      }
     
   }
   
   // Convert hardware coordinates to MPI rank
   *ierr = MPIX_Torus2rank( hwcoords, rank );
   
   //Debug
   {
      int _simcoords[3];
      
      *ierr = MPIX_Rank2torus( *rank, hwcoords );
      if ( *ierr ) {
        fprintf( stderr, "(*error*) in bgq_sim2rank() Unable to convert back to hwcoords");
        return;
      } 
      
      bgq_hw2sim(hwcoords, _simcoords);
      
       for(int i=0; i<bgq_part.simpart_ndims; i++) {
        if ( simcoords[i] != _simcoords[i] ) { 
		  fprintf( stderr, "(*error*) in bgq_sim2rank() Unable to retrieve original "
		                   "simulation coordinates -> dim %d, in: %d, out:%d\n",
		                   i, simcoords[i], _simcoords[i] );
		  *ierr = -2;
		  return;
        }
      }
   
   }
   
}


/*****************************************************************************************
* Initialize bgq_partition structure.
*   - This will be done using MPIX_Inithw
*****************************************************************************************/

void FNAME(bgq_initpart)
( int* simpart_ndims, int simpart[], int* ierr )
{
  
  MPIX_Hardware_t bgq_hw;
  int torusdims, mpi_rank;
  int hw_nnodes, sim_nnodes;
  
  // Get hardware / network information from system
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPIX_Hardware( &bgq_hw );
  MPIX_Torus_ndims( &torusdims ); // this will always be 5
  
  bgq_part.nnodes = bgq_hw.psize;
  bgq_part.hwpart_torusdims = torusdims;
  for( int i=0; i < torusdims; i++ ) {
    bgq_part.hwpart[i] = bgq_hw.Size[i];
    bgq_part.hwpart_coords[i] = bgq_hw.Coords[i];
    bgq_part.hwpart_istorus[i] = bgq_hw.isTorus[i];
  }
  bgq_part.hwpart[torusdims] = bgq_hw.ppn;
  bgq_part.hwpart_coords[torusdims] = bgq_hw.coreID;
  bgq_part.hwpart_istorus[torusdims] = 1; // The T coordinate is always periodic
  bgq_part.mpi_rank = mpi_rank;

  // Simulation partition
  bgq_part.simpart_ndims = *simpart_ndims;
  for( int i=0; i < *simpart_ndims; i++) bgq_part.simpart[i] = simpart[i];
   
  // Check the total number of nodes
  sim_nnodes=1;
  for( int i=0; i < *simpart_ndims; i++ ) 
       sim_nnodes *= bgq_part. simpart[i]; 
   
  if ( bgq_part.nnodes != sim_nnodes ) {
     if ( mpi_rank == 0 ) 
        fprintf( stderr, "(*error*) Simulation and hardware partition sizes don't match");
     *ierr = -1;
     return;
  }

  // Generate mapping between simulation and hardware partitions
  int part_map[ HWPART_MAX_DIMS ];
  for( int i=0; i < HWPART_MAX_DIMS; i++ ) part_map[i] =i;

  // Try without splitting first
  *ierr = 0;
  
  bgq_part.hwpart_splitdim = -1;
  bgq_part.mappart_ndims = bgq_part.hwpart_torusdims + 1;
  for( int i=0; i<bgq_part.mappart_ndims; i++) 
     bgq_part.mappart[i] = bgq_part.hwpart[i];

  int done = 0;

  done = bgq_findmap( 0, bgq_part.simpart_ndims, bgq_part.simpart, 
                         bgq_part.mappart_ndims, bgq_part.mappart, part_map, 
                         bgq_part.nmap, bgq_part.map );
  
  // if unable to find mapping without splitting, try to split coordinates starting from 
  // the last ( T )
  
  for( int dim = bgq_part.hwpart_torusdims + 1; dim >= 0 && !done; dim-- ) {
	  // Only try to split dimensions larger than 2
	  if ( bgq_part.hwpart[dim] > 2 ) {
		 bgq_part.mappart_ndims = bgq_part. hwpart_torusdims + 2;
		 bgq_part.hwpart_splitdim = dim;
		 
		 int n = bgq_part. hwpart[dim];
		 for( int d0 = 2; ( n / d0 ) * d0 == n && !done; d0*=2 ) {
			// Split dimention
			for( int i=0; i < dim-1 ; i++ ) 
			   bgq_part. mappart[i] = bgq_part. hwpart[i];
			bgq_part. mappart[dim]   = d0;
			bgq_part. mappart[dim+1] = bgq_part. hwpart[dim] / d0;
			for( int i=dim+2; i < bgq_part. mappart_ndims; i++ )
			  bgq_part. mappart[i] = bgq_part. hwpart[i-1];
			
			// Try finding map with the coordinate split
			done = bgq_findmap( 0, bgq_part. simpart_ndims, bgq_part. simpart, 
							       bgq_part. mappart_ndims, bgq_part. mappart, part_map, 
							       bgq_part. nmap, bgq_part. map );
		 }
	  
	  
	  }
  }

  if ( done ) {
     if ( mpi_rank == 0 ) {
        // Print partition information     
        char dimLabel[ HWPART_MAX_DIMS ];
        for(int i=0; i < bgq_part.hwpart_torusdims; i++) dimLabel[i] = 'A'+i;
        dimLabel[ bgq_part.hwpart_torusdims ] = 'T';
        
        printf( "\n------------------------- BlueGene/Q Topology -------------------------\n" );
        printf( "Simulation = { %d", bgq_part.simpart[0] );
        for(int i=1; i<bgq_part.simpart_ndims; i++) printf(", %d", bgq_part.simpart[i] );
        printf(" }\n");
         
        printf( "Hardware   = { %d", bgq_part.hwpart[0] );
        for(int i=1; i<bgq_part.hwpart_torusdims+1; i++) printf(", %d", bgq_part.hwpart[i] );
        if ( bgq_part.hwpart_splitdim < 0 ) {
           printf(" }\n");
        } else {
           printf(" }, %c is split into {%d, %d}\n", dimLabel[bgq_part.hwpart_splitdim],
                       bgq_part.mappart[bgq_part.hwpart_splitdim ],
                       bgq_part.mappart[bgq_part.hwpart_splitdim + 1] );
        }
        
        printf( "Coord. mapping");
        for(int i=0; i < bgq_part.simpart_ndims; i++ ) {
           printf( (i)?", ":" " );
           printf( "%c -> {", 'x'+i );
           if ( bgq_part.hwpart_splitdim < 0 ) {
              for( int j=0; j < bgq_part.nmap[i]; j++ ) 
                 printf(" %c", dimLabel[ bgq_part.map[i][j] ] );
           } else {
              for( int j=0; j < bgq_part.nmap[i]; j++ ) {
                 int idx = bgq_part.map[i][j];
                 if (idx > bgq_part.hwpart_splitdim ) idx--;
                 printf(" %c", dimLabel[ idx ] );
                 if ( bgq_part.map[i][j] == bgq_part.hwpart_splitdim ) printf("0");
                 if ( bgq_part.map[i][j] == bgq_part.hwpart_splitdim + 1) printf("1");
              }
           }
           printf(" }");
        }
        printf( "\n" );
        printf( "-----------------------------------------------------------------------\n" );
        
        // Issue a warning if the mapping required a split of a coordinate other than T
        if ( bgq_part.hwpart_splitdim >= 0 && 
             bgq_part.hwpart_splitdim !=  bgq_part.hwpart_torusdims ) {
		   fprintf( stderr, "(*warning*) Mapping between simulation and hardware partitions required splitting\n");
		   fprintf( stderr, "(*warning*) the %c torus coordinate, network topology will be suboptimal.\n", dimLabel[bgq_part.hwpart_splitdim]);
		   fprintf( stderr, "(*warning*) Please consider choosing a different simulation partition, check the\n");
		   fprintf( stderr, "(*warning*) documentation for details");
        }
     }   
     
     *ierr = 0;
  } else {
	 // If we reach this then we were unable to find a mapping between the simulation and
	 // the hardware partition and the code should default to using MPI_Cart_Create  
	 // Note: I don't think this will ever happen...

     bgq_part. mappart_ndims = -1;
     if ( mpi_rank == 0 ) 
        fprintf( stderr, "(*warning*) Unable to find mapping between simulation and "
                         "hardware partitions\n");
     *ierr = 1;
  }
  
}




#include <unistd.h>

/* returns BlueGene/Q process memory in kbytes */

#include <spi/include/kernel/memory.h>

void FNAME(bgq_used_memory_f)
( long *used_kb )
{
  uint64_t allocated_memory = 0;
  
  /* Size in bytes of the heap size*/
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &allocated_memory);  
  *used_kb = allocated_memory / 1024;
}


/* BG/Q specific timers */
#define HAVE_TIMER

#include <hwi/include/bqc/testint_inlines.h>
static double tick_conv;

void FNAME(timer_init)
(void)
{
   /* The core frequency should always be 1600 MHz */
   
   /* This causes a signal 11 */
   /*tick_conv = 1.0e-6 / (double) TI_CoreFrequency();*/

   tick_conv = 1.0e-6 / 1600.0;
}


uint64_t FNAME(timer_ticks)
( void )
{
   return GetTimeBase();
}

double FNAME(timer_interval_seconds)
(uint64_t  *start, uint64_t *end)
{
  /* we should check for overflow */
  
  return (*end - *start) * tick_conv;
}

double FNAME(timer_resolution)
( void )
{
  uint64_t        start, end;
  start = GetTimeBase();
  end   = GetTimeBase();
  
  return (end-start)*tick_conv;
}

double FNAME(timer_cpu_seconds)
( void )
{
    return (double) GetTimeBase() * tick_conv;
}

#endif

/*****************************************************************************************
  MPI Timers
*****************************************************************************************/
#ifdef __MPI_TIMER__
#define HAVE_TIMER

#include <stdint.h>

void FNAME(timer_init)
(void)
{
  /* No initialization required for these timers */
}

uint64_t FNAME(timer_ticks)
( void )
{
  /* Return number of microsseconds */
  return 1.0e6*MPI_Wtime();
}

double FNAME(timer_interval_seconds)
(uint64_t  *start, uint64_t *end)
{
  return (*end - *start) * 1.0e-6;
}

double FNAME(timer_resolution)
( void )
{
  return MPI_Wtick();
}

double FNAME(timer_cpu_seconds)
( void )
{
    return MPI_Wtime();
}

#endif /* __MPI_TIMER__ */


/*****************************************************************************************
  Default timers (use POSIX timers)
*****************************************************************************************/

#ifndef HAVE_TIMER

/*
  Use standard POSIX timers based on gettimeofday. (This is actually how
  LAM and OpenMPI implement MPI_Wtime). Minimum resolution is 1 microssecond.
  
  The actual resolution is measured by finding the minimum difference >0 
  between succsessive gettimeofday calls. 
*/

#include <sys/time.h>
#include <stdint.h>

void FNAME(timer_init)
(void)
{
  /*  No initialization is required for this type of timers */
}

uint64_t FNAME(timer_ticks)
(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  
  return ((uint64_t)tv.tv_sec)*1000000 + (uint64_t)tv.tv_usec;
}

double FNAME(timer_interval_seconds)
(uint64_t *start, uint64_t *end)
{
  return (*end - *start) * 1.0e-6;
}

double FNAME(timer_cpu_seconds)
( void )
{
    struct timeval tv;
    double wtime;
    gettimeofday(&tv, NULL);
    wtime = tv.tv_sec;
    wtime += (double)tv.tv_usec * 1.0e-6;
    return wtime;
}

double FNAME(timer_resolution)
( void )
{
  struct timeval tv1, tv2;
   
  gettimeofday(&tv1, NULL);
  do {
       gettimeofday(&tv2, NULL);
  } while (tv1.tv_usec == tv2.tv_usec);
  
  return (double)(tv2.tv_usec - tv1.tv_usec) * 1.0e-6;
}


#endif 
