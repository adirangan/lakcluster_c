// lakcluster biclustering algorithm: (64-bit) ver18 (dated 060815 -- 071417). ;
// Accounts for sparsity, and allows for row/col controls as well as multiple (categorical/mds) row covariates. ;
// Note that low-rank-update has not yet been implemented for QR_strategy "ZtSWn". ;
// Written by Aaditya Rangan, with thanks to Chris Chang (see https://www.cog-genomics.org/plink2/). ;
// This program is free software: ;
// you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
// See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with this program.  
// If not, see <http://www.gnu.org/licenses/>.

#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <omp.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

#include <emmintrin.h>

#define PI 3.141592653589793
#define CACHELINE 64 // assumed number of bytes per cache line, for alignment
#define CACHELINE_DBL (CACHELINE / sizeof(double))
#define bget_off(bmr,nr) (((bmr[(nr)/BIT8] >> (7-((nr)%BIT8))) & 1) ^ 1)
#define bget__on(bmr,nr) (((bmr[(nr)/BIT8] >> (7-((nr)%BIT8))) & 1) ^ 0)
#define bget____(bmr,nr) ((int)(((bmr[(nr)/BIT8] >> (7-((nr)%BIT8))) & 1) << 1)-(int)1)
#define b_copy_u(bmr,umr,nr) (bmr[(nr)/BIT8] |= (umr[(nr)]>0) ? (1 << (7 - ((nr)%BIT8))) : 0)
#define bset_off(bmr,nr) (bmr[(nr)/BIT8] &= ~(1 << (7 - ((nr)%BIT8))))
#define bset__on(bmr,nr) (bmr[(nr)/BIT8] |=  (1 << (7 - ((nr)%BIT8))))
#define MAX_THREADS 144
#define MAX_THREADS_P1 145
#define RSNSIZE 64
#define FNAMESIZE 2048
/* #define MALLOC_DEFAULT_MB 2176 */
/* #define MALLOC_DEFAULT_MB 4352 */
/* #define MALLOC_DEFAULT_MB 8704 */
/* #define MALLOC_DEFAULT_MB 17408 */
#define MALLOC_DEFAULT_MB 31744
/* #define MALLOC_DEFAULT_MB 34816 */
#define MULTIPLEX_DIST 960
#define MULTIPLEX_2DIST (MULTIPLEX_DIST * 2)
#define BIT8 8
#define BITJ 16

/* 5 --> 0101 */
#define FIVEMASK 0x5555555555555555LU
/* cachealign(val) returns the multiple of 64 geq to val */
#define CACHEALIGN(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - 1)))
/* cachealign_dbl(val) returns the multiple of 8 geq to val */
#define CACHEALIGN_DBL(val) ((val + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1)))

#define RET_SUCCESS 0
#define RET_NOMEM 1
#define RET_READ_FAIL 2
#define RET_INVALID_FORMAT 3
#define RET_CALC_NOT_YET_SUPPORTED 4
#define RET_INVALID_CMDLINE 5

const char ver_str[] =  "lakcluster_ver18: "
  " linear-affine k-dimensional biclustering algorithm "
  " [64-bit] "
  " (July 2017)"
  " Aaditya Rangan, cims nyu edu\n"
  " "
  " Note: low-rank-update for QR_strategy ""ZtSWn"" not yet implemented ";

/* ---------------------------------------------------------------- */

#define GLOBAL_omp_unused -1 /* unused */
#define GLOBAL_omp_off 0 /* parallelize manually */
#define GLOBAL_omp__on 1 /* use omp in addition */
#define TYPE_CAT 0 /* used in dra_normalize */
#define TYPE_CNT 1 /* used in dra_normalize */
#define TYPE_p0 0
#define TYPE_pm 1
#define TYPE_00 2
#define TYPE_pz 3
#define TYPE_0z 4
#define TYPE_p_ 5
#define TYPE_0_ 6
#define TYPE_pZ 7
#define TYPE_0Z 8
const char * TYPE_name[] = { "TYPE_p0","TYPE_pm","TYPE_00","TYPE_pz","TYPE_0z","TYPE_p_","TYPE_0_","TYPE_pZ","TYPE_0Z" };
int addressable_type_p0=TYPE_p0;
int addressable_type_pm=TYPE_pm;
int addressable_type_00=TYPE_00;
int addressable_type_pz=TYPE_pz;
int addressable_type_0z=TYPE_0z;
int addressable_type_p_=TYPE_p_;
int addressable_type_0_=TYPE_0_;
int addressable_type_pZ=TYPE_pZ;
int addressable_type_0Z=TYPE_0Z;
#define SPACING_b 0
#define SPACING_j 1
#define SPACING_a 2
const char * SPACING_name[] = { "SPACING_b","SPACING_j","SPACING_a" };
int addressable_spacing_b=SPACING_b;
int addressable_spacing_j=SPACING_j;
int addressable_spacing_a=SPACING_a;
#define AJDK_TOT 15
#define AJDK_0_0 0
#define AJDK_1_0 1
#define AJDK_2_0 2
#define AJDK_3_0 3
#define AJDK_4_0 4
#define AJDK_0_1 5
#define AJDK_1_1 6
#define AJDK_2_1 7
#define AJDK_3_1 8
#define AJDK_4_1 9
#define AJDK_0_2 10
#define AJDK_1_2 11
#define AJDK_2_2 12
#define AJDK_3_2 13
#define AJDK_4_2 14
#define GLOBAL_gamma__log 0
#define GLOBAL_gamma_sqrt 1
#define GLOBAL_gamma__lin 2
#define GLOBAL_gamma_rows 3

/* global variables used for timing */
#define GLOBAL_NTICKS 8
clock_t GLOBAL_t_start[GLOBAL_NTICKS],GLOBAL_t_final[GLOBAL_NTICKS];
struct timeval GLOBAL_d_start[GLOBAL_NTICKS],GLOBAL_d_final[GLOBAL_NTICKS];
long GLOBAL_l_msec[GLOBAL_NTICKS],GLOBAL_l_ssec[GLOBAL_NTICKS],GLOBAL_l_usec[GLOBAL_NTICKS]; 
double GLOBAL_elct[GLOBAL_NTICKS],GLOBAL_elrt[GLOBAL_NTICKS];

/* These are the kappa-values used for continuous-covariate correction; note that the former are the square of the latter. */
double GLOBAL_mds_loop_scale_factor_[6] = {1.00 , 0.34 , 0.19 , 0.13 , 0.10 , 0.09};
double GLOBAL_mds_half_scale_factor_[6] = {1.0000 , 0.5747 , 0.4338 , 0.3613 , 0.3154 , 0.3000 } ; 

char GLOBAL_CWD[FNAMESIZE];
char GLOBAL_DIR_BASE[FNAMESIZE]="\0";
char GLOBAL_DIR_NAME[FNAMESIZE]="\0";
char GLOBAL_DIR_XPRE[FNAMESIZE]="\0";
int GLOBAL_verbose=0; // set to 1 to see sysconf(_SC_NPROCESSORS_ONLN) ;
int GLOBAL_thread_count=1; // Set >1 to use pthreads to parallelize across nbins ;
int GLOBAL_omp_type=GLOBAL_omp_off; // manual parallelization using pthreads ;
int GLOBAL_1_icache_linesize=64;
int GLOBAL_1_icache_size=1024;
int GLOBAL_1_icache_assoc=16;
int GLOBAL_1_dcache_linesize=64;
int GLOBAL_1_dcache_size=1024;
int GLOBAL_1_dcache_assoc=16;
int GLOBAL_2_cache_linesize=64;
int GLOBAL_2_cache_size=1024;
int GLOBAL_2_cache_assoc=16;
int GLOBAL_3_cache_linesize=64;
int GLOBAL_3_cache_size=1024;
int GLOBAL_3_cache_assoc=16;
int GLOBAL_4_cache_linesize=64;
int GLOBAL_4_cache_size=1024;
int GLOBAL_4_cache_assoc=16;
double GLOBAL_tolerance=0.000000000001;
unsigned int GLOBAL_recursion_limit=1024*32;
int GLOBAL_LBITS=8;
char GLOBAL_TEST_TYPE[FNAMESIZE]="\0";
char GLOBAL_TEST_TYP2[FNAMESIZE]="\0";
int GLOBAL_TEST_sparse=1;
int GLOBAL_TEST_A_n_rows=12;
int GLOBAL_TEST_A_n_cols=13;
int GLOBAL_TEST_Z_n_rows=15;
int GLOBAL_TEST_Y_n_cols=11;
int GLOBAL_TEST_T_n_cols=3;
int GLOBAL_TEST_niter=2;
double GLOBAL_TEST_mrand=0.0; /* 0-fraction of masks (100-percent - 1-fraction) */
char GLOBAL_skip[FNAMESIZE]="\0";
char GLOBAL_QR_strategy[FNAMESIZE]="\0";
char GLOBAL_QC_strategy[FNAMESIZE]="\0";
long long int GLOBAL_D_MLT = 4294967296; /* #define D_MLT 65536 */ /* #define D_MLT 32768 */ /* #define D_MLT 256 */
int GLOBAL_B_MLT = 32;
int GLOBAL_gamma_type = GLOBAL_gamma__log;
double GLOBAL_gamma = 0.125;
int addressable_1=1;
int addressable_0=0;
int addressable_int_length=128;
int addressable_int[128] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127};
char GLOBAL_A_n_name[FNAMESIZE]="\0";
char GLOBAL_A_t_name[FNAMESIZE]="\0";
char **GLOBAL_A_n_name_=NULL;
char **GLOBAL_A_t_name_=NULL;
int  *GLOBAL_A_n_rows_=NULL;
int  GLOBAL_A_n_cols=0;
char GLOBAL_A_n_rind[FNAMESIZE]="\0";
char **GLOBAL_A_n_rind_=NULL;
char GLOBAL_A_n_cind[FNAMESIZE]="\0";
char GLOBAL_Z_n_name[FNAMESIZE]="\0";
char GLOBAL_Z_t_name[FNAMESIZE]="\0";
char **GLOBAL_Z_n_name_=NULL;
char **GLOBAL_Z_t_name_=NULL;
int  *GLOBAL_Z_n_rows_=NULL;
char GLOBAL_Z_n_rind[FNAMESIZE]="\0";
char **GLOBAL_Z_n_rind_=NULL;
char GLOBAL_Y_n_name[FNAMESIZE]="\0";
char GLOBAL_Y_t_name[FNAMESIZE]="\0";
char **GLOBAL_Y_n_name_=NULL;
char **GLOBAL_Y_t_name_=NULL;
int  GLOBAL_Y_n_cols=0;
char GLOBAL_Y_n_cind[FNAMESIZE]="\0";
char GLOBAL_W_n_name[FNAMESIZE]="\0";
char GLOBAL_W_t_name[FNAMESIZE]="\0";
char **GLOBAL_W_n_name_=NULL;
char **GLOBAL_W_t_name_=NULL;
char GLOBAL_T_n_name[FNAMESIZE]="\0";
char GLOBAL_T_t_name[FNAMESIZE]="\0";
char **GLOBAL_T_n_name_=NULL;
char **GLOBAL_T_t_name_=NULL;
int  GLOBAL_T_n_cols=0;
char GLOBAL_T_n_cind[FNAMESIZE]="\0";
char GLOBAL_S_n_name[FNAMESIZE]="\0";
char GLOBAL_S_t_name[FNAMESIZE]="\0";
char **GLOBAL_S_n_name_=NULL;
char **GLOBAL_S_t_name_=NULL;
char GLOBAL_out_name[FNAMESIZE]="\0";
char GLOBAL_scorebox_out_xdrop[FNAMESIZE]="\0";
int GLOBAL_scorebox_row_max=0;
int GLOBAL_scorebox_row_num=0;
int GLOBAL_scorebox_col_max=0;
int GLOBAL_scorebox_col_num=0;

#define rup(A,B) ((A) + !!((A)%(B))*((B) - ((A)%(B))))
#define maximum(A,B) ((A) > (B) ? (A) : (B))
#define minimum(A,B) ((A) < (B) ? (A) : (B))
#define periodize(A,B,C) ((A) < (B) ? (A) + (C) - (B) : ((A) >= (C) ? (A) - (C) + (B) : (A)))
#define crop(A,B,C) ((A) < (B) ? (B) : ((A) > (C) ? (C) : (A)))
#define rand01 ((double)rand()/(double)RAND_MAX)
#define bsize(A) ((rup((A) + ((BITJ - ((A) % BITJ)) % BITJ),POPLENGTH))/BIT8)
#define psize(A) ((rup((A) + ((BITJ - ((A) % BITJ)) % BITJ),POPLENGTH))/POPLENGTH)

#define RGB3 3 
#define RGBA 4 
#define UCHAR_MAX 255 
#define POPLENGTH (MULTIPLEX_2DIST / 128 * sizeof(__m128i) * BIT8)

unsigned char popcount_uchar[256] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};

struct S_handle
{
  char out_xdrop_name[FNAMESIZE]; /* name of out_xdrop file */
  int out_xdrop_nlines; /* number of lines in out_xdrop */
  int out_xdrop_nrows; /* number of row-indices stored in out_xdrop */
  int out_xdrop_ncols; /* number of col-indices stored in out_xdrop */
  int *out_xdrop; /* actual array for out_xdrop */
  int *mr_srt; /* row-indices from out_xdrop */
  int *mc_srt; /* col-indices from out_xdrop */
  int *mr_lnb; /* local nb from row-indices */
  int *mr_lmr; /* local mr from row-indices */
  int row_max; /* maximum row index for score calculation */
  int row_num; /* number of row indices for score calculation */
  int col_max; /* maximum col index for score calculation */
  int col_num; /* number of col indices for score calculation */
  int *rval; /* array of r-values (rows remaining) */
  int *rdrop; /* array of rdrop-values */
  int *cval; /* array of c-values (cols remaining) */
  int *cdrop; /* array of cdrop-values */
  int *A_rpop_j_total; /* output array (size row_num-x-col_num) */
  int *A_cpop_j; /* output array (size row_num-x-col_num) */
  int *Irem; /* output array (size row_num-x-col_num) */
  double *QR_avg; /* output array (size row_num-x-col_num) */
  double *QC_avg; /* output array (size row_num-x-col_num) */
  int nr; /* temporary row index */
  int nc; /* temporary row index */
};

struct L_handle
{
  int length; /* total length of lf -- often of the form row_stride*col_stride*lyr_stride */
  double *lf; /* holds data in memory -- often of the form lf[nr + nc*row_stride + nl*row_stride*col_stride] */
  int spacing_row; /* spacing of row-index */
  int row_stride; /* associated stride */
  int spacing_col; /* spacing of column-index */
  int col_stride; /* associated stride */
  int spacing_lyr; /* spacing of layer-index */
  int lyr_stride; /* associated stride -- often unused */
};

struct M_handle
{
  unsigned char *mark;
  int header_length; /* header length in bytes expected for array data file */
  char A_name[FNAMESIZE]; /* name of array data file */ 
  FILE *A_fp; /* stores array data file pointer */ unsigned char *Ara; /* stores array data in memory */
  int bitj; /* set to 16, this is the chunk-size used when storing rows and columns */
  int nrows; int ncols; /* dimensions of A */
  int nrows_extend; int ncols_extend; /* additional rows and columns required to reach a multiple of bitj */
  long long int rpop_b; long long int rpop_j; /* popcounts for mr_b and mr_j, respectively */
  long long int cpop_b; long long int cpop_j; /* popcounts for mc_b and mc_j, respectively */
  unsigned char *wX; /* holds data in memory, spaced by POPLENGTH */
  unsigned char *mr_b; unsigned char *mr_j; /* bitmasks for the rows; original and current */
  unsigned char *mc_b; unsigned char *mc_j; /* bitmasks for the cols; original and current */
  double *rsum; /* holds the row-sums of A */
  unsigned int *m_a_; /* m_a_[nr_j] = nr_a */ /* position nr_a of mr_b,mr_j corresponds to nr_j^th positive entry of mr_j*/
  unsigned int *m_b_; /* m_b_[nr_j] = nr_b */ /* and nr_b^th positive entry of mr_b */  
  unsigned int *n_a_; /* n_a_[nc_j] = nc_a */ /* position nc_a of mc_b,mc_j corresponds to nc_j^th positive entry of mc_j*/
  unsigned int *n_b_; /* n_b_[nc_j] = nc_b */ /* and nc_b^th positive entry of mc_b */  
  unsigned long mr_length; unsigned long mc_length; /* bsize(nrows) and bsize(ncols), respectively */
  unsigned long cc_length; /* (ncols + ncols_extend)/BIT8 */
  unsigned long wt_length; /* length of tail required to cap lengthof M_handle to a multiple of POPLENGTH */
  unsigned char *wt; /* tail (not used) */
  unsigned long long int length_total; /* total length of M_handle */
  int ncols_per_z; /* for use with xcalc_dra_to_M_x_bitbuffer */
  double max_d; double min_d; double mlt_d; /* for use with xcalc_dra_to_M_x_bitbuffer */
};

#include "./dir_h/bcc.h"
#include "./dir_h/dcc.h"

/* thread management */
int GLOBAL_NBINS=0;
int GLOBAL_Ireq=0;
int *GLOBAL_lnb = NULL;
int GLOBAL_nf=0; int GLOBAL_nf_cur=0; int GLOBAL_nf_ind=0; int GLOBAL_nf_opn=0;
int GLOBAL_tint[MAX_THREADS];
void *GLOBAL_tvp[MAX_THREADS][128];
pthread_t GLOBAL_threads[MAX_THREADS];
unsigned long long int GLOBAL_ops_f_[MAX_THREADS];
unsigned long long int GLOBAL_ops_f_sum=0;
unsigned long long int GLOBAL_ops_b_[MAX_THREADS];
unsigned long long int GLOBAL_ops_b_sum=0;

inline void ping(){ printf(" %% ping\n");}
inline void pong(){ printf(" %% pong\n");}

/* ---------------------------------------------------------------- */

/* manually managed memory stack */
unsigned char* wkspace=NULL;
unsigned char* wkspace_base=NULL;
int GLOBAL_wkspace_point=1;
struct wkspace_point
{
  struct wkspace_point * parent;
  struct wkspace_point * child;
  long long int check;
  long long int *point;
};
struct wkspace_point *wkspace_point_0=NULL,*wkspace_point_t=NULL;
long malloc_size_mb = MALLOC_DEFAULT_MB;
long long int wkspace_left=0;
long long int wkspace_used=0;

/* lakcluster functions */
#include "./dir_c/mda_io.c"
#include "./dir_c/fillzero.c"
#include "./dir_c/wkspace.c" 
#include "./dir_c/GLOBAL_pthread.c" 
#include "./dir_c/popcount.c"
#include "./dir_c/raprintf.c"
#include "./dir_c/binary_read.c"
#include "./dir_c/updateglobals.c"
#include "./dir_c/L_init.c"
#include "./dir_c/M_init.c"
#include "./dir_c/rastats.c"
#include "./dir_c/Quicksort.c"
#include "./dir_c/PNMfile.c"
#include "./dir_c/An_ajdk_v.c"
#include "./dir_c/An_v.c"
#include "./dir_c/AnZt_vv.c"
#include "./dir_c/AnAt_vv.c"
#include "./dir_c/AtTYn_vv.c"
#include "./dir_c/AtTAn_vv.c"
#include "./dir_c/AttYn_vv.c"
#include "./dir_c/AttAn_vv.c"
#include "./dir_c/xcalc.c"
#include "./dir_c/An_ZtSn_ww.c"
#include "./dir_c/An_ZtSWn_Yt_vv.c"
#include "./dir_c/AnZt_S_WnYt_vv.c"
#include "./dir_c/An_ZtSWn_Yt_ww.c"
#include "./dir_c/At_T_YnWt_S_Zn_vv.c"
#include "./dir_c/AtTYn____WtSZn_vv.c"
#include "./dir_c/AttYn____WtsZn_vv.c"
#include "./dir_c/At_T_YnWt_S_Zn_ww.c"
#include "./dir_c/xdrop.c"
#include "./dir_c/bcc.c"
#include "./dir_c/bcc_lf.c"
#include "./dir_c/bcc_flattenloop.c"
#include "./dir_c/bcc_lrup.c"
#include "./dir_c/bcc_sumscores.c"
#include "./dir_c/lakcluster_driver.c"

/* dexcluster functions */
#include "./dir_c/D_AtTn_ZtSn.c"
#include "./dir_c/TAnZtS.c"
#include "./dir_c/dcc.c"
#include "./dir_c/dcc_lf.c"
#include "./dir_c/dcc_sumscores.c"
#include "./dir_c/dexcluster_driver.c"

/* scorebox functions */
#include "./dir_c/bcc_scorebox.c"
#include "./dir_c/S_init.c"
#include "./dir_c/lakcluster_scorebox.c"

/* ---------------------------------------------------------------- */

void set_globals()
{
  int verbose=0;
  int np=0,nb=0;
  struct stat st_tmp = {0};
  if (verbose){ printf(" %% setting GLOBAL_DIR_NAME\n");}
  if (strcmp(GLOBAL_DIR_NAME,"\0")){ /* do nothing */ }
  else /* if GLOBAL_DIR_NAME not defined */{
    if (strcmp(GLOBAL_out_name,"\0")){
      if (!strcmp(GLOBAL_DIR_BASE,"\0")){ sprintf(GLOBAL_DIR_NAME,"%s/dir_%s",GLOBAL_CWD,GLOBAL_out_name);}
      else /* if GLOBAL_DIR_BASE */{ sprintf(GLOBAL_DIR_NAME,"%s/dir_%s",GLOBAL_DIR_BASE,GLOBAL_out_name);}
      /* if (strcmp(GLOBAL_out_name,"\0"))){ } */}
    /* if not defined */}
  if (verbose>-1){ printf(" %% setting GLOBAL_DIR_NAME: %s (%d)\n",GLOBAL_DIR_NAME,strcmp(GLOBAL_DIR_NAME,"\0"));}
  if (strcmp(GLOBAL_DIR_NAME,"\0") && stat(GLOBAL_DIR_NAME, &st_tmp) == -1) {  printf(" %% mkdir %s;\n",GLOBAL_DIR_NAME); mkdir(GLOBAL_DIR_NAME, 0755);}
  if (verbose){ printf(" %% initializing GLOBAL_lnb\n");}
  GLOBAL_lnb = (int *) wkspace_all0c(GLOBAL_NBINS*sizeof(int)); for (nb=0;nb<GLOBAL_NBINS;nb++){ GLOBAL_lnb[nb] = nb;}
  if (verbose){ printf(" %% initializing GLOBAL_tint\n");}
  for (GLOBAL_nf=0;GLOBAL_nf<MAX_THREADS;GLOBAL_nf++){ GLOBAL_tint[GLOBAL_nf]=GLOBAL_nf;}
}

int main(int argc, char** argv) {
  unsigned char* wkspace_ua;
  int cur_arg = 1,ii=0;
  char* argptr;
  char* bubble;
  int retval=0;
  int cmdline_param=0;
  if (omp_get_max_threads()>MAX_THREADS){ printf(" %% Warning! global variable MAX_THREADS set to %d when it should be %d\n",MAX_THREADS,omp_get_max_threads()); exit(0);}
  if (!getcwd(GLOBAL_CWD,FNAMESIZE)){ printf(" %% Warning! GLOBAL_CWD not read in main\n");} else{ if (GLOBAL_verbose>0){ printf(" %% GLOBAL_CWD: %s\n",GLOBAL_CWD);};} 
  if (GLOBAL_verbose>0){ printf(ver_str);}
  GLOBAL_thread_count = sysconf(_SC_NPROCESSORS_ONLN); if (GLOBAL_verbose>-1){ printf("sysconf(_SC_NPROCESSORS_ONLN) = %d\n",GLOBAL_thread_count);}
  GLOBAL_1_icache_linesize = sysconf(_SC_LEVEL1_ICACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_ICACHE_LINESIZE) = %d\n",GLOBAL_1_icache_linesize);}
  GLOBAL_1_icache_size = sysconf(_SC_LEVEL1_ICACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_ICACHE_SIZE) = %d\n",GLOBAL_1_icache_size);}
  GLOBAL_1_icache_assoc = sysconf(_SC_LEVEL1_ICACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_ICACHE_ASSOC) = %d\n",GLOBAL_1_icache_assoc);}
  GLOBAL_1_dcache_linesize = sysconf(_SC_LEVEL1_DCACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_DCACHE_LINESIZE) = %d\n",GLOBAL_1_dcache_linesize);}
  GLOBAL_1_dcache_size = sysconf(_SC_LEVEL1_DCACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_DCACHE_SIZE) = %d\n",GLOBAL_1_dcache_size);}
  GLOBAL_1_dcache_assoc = sysconf(_SC_LEVEL1_DCACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL1_DCACHE_ASSOC) = %d\n",GLOBAL_1_dcache_assoc);}
  GLOBAL_2_cache_linesize = sysconf(_SC_LEVEL2_CACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL2_CACHE_LINESIZE) = %d\n",GLOBAL_2_cache_linesize);}
  GLOBAL_2_cache_size = sysconf(_SC_LEVEL2_CACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL2_CACHE_SIZE) = %d\n",GLOBAL_2_cache_size);}
  GLOBAL_2_cache_assoc = sysconf(_SC_LEVEL2_CACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL2_CACHE_ASSOC) = %d\n",GLOBAL_2_cache_assoc);}
  GLOBAL_3_cache_linesize = sysconf(_SC_LEVEL3_CACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL3_CACHE_LINESIZE) = %d\n",GLOBAL_3_cache_linesize);}
  GLOBAL_3_cache_size = sysconf(_SC_LEVEL3_CACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL3_CACHE_SIZE) = %d\n",GLOBAL_3_cache_size);}
  GLOBAL_3_cache_assoc = sysconf(_SC_LEVEL3_CACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL3_CACHE_ASSOC) = %d\n",GLOBAL_3_cache_assoc);}
  GLOBAL_4_cache_linesize = sysconf(_SC_LEVEL4_CACHE_LINESIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL4_CACHE_LINESIZE) = %d\n",GLOBAL_4_cache_linesize);}
  GLOBAL_4_cache_size = sysconf(_SC_LEVEL4_CACHE_SIZE); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL4_CACHE_SIZE) = %d\n",GLOBAL_4_cache_size);}
  GLOBAL_4_cache_assoc = sysconf(_SC_LEVEL4_CACHE_ASSOC); if (GLOBAL_verbose>2){ printf("sysconf(_SC_LEVEL4_CACHE_ASSOC) = %d\n",GLOBAL_4_cache_assoc);}
  if (GLOBAL_thread_count == -1) { GLOBAL_thread_count = 1; } else if (GLOBAL_thread_count > MAX_THREADS) { GLOBAL_thread_count = MAX_THREADS; }
  while (cur_arg < argc) {
    argptr = argv[cur_arg];
    if (0){ /* do nothing */}
    else if (!strcmp(argptr, "--nthreads")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --nthreads parameter."); exit(RET_INVALID_CMDLINE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid --something parameter.\n"); exit(RET_INVALID_CMDLINE);}
      GLOBAL_thread_count = ii; printf(" %% nthreads %d\n",GLOBAL_thread_count);
      cur_arg += 2; /* else if pminor_cols */}
    else if (!strcmp(argptr, "--cmdline_param")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --cmdline_param parameter."); exit(RET_INVALID_CMDLINE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid parameter.\n"); exit(RET_INVALID_CMDLINE);}
      cmdline_param = ii; printf(" %% cmdline_param %d\n",cmdline_param);
      cur_arg += 2; /* else if pminor_cols */}
    else { printf("Error: Invalid argument (%s).", argv[cur_arg]); exit(RET_INVALID_CMDLINE);}
    /* while (cur_arg < argc) { } */}
  bubble = (char*)malloc(67108864 * sizeof(char));
  if (!bubble) { printf("Error: not enough memory!\n");}
  wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
  if ((malloc_size_mb > MALLOC_DEFAULT_MB) && !wkspace_ua) { printf("%ld MB memory allocation failed.  Using default allocation behavior.\n", malloc_size_mb); malloc_size_mb = MALLOC_DEFAULT_MB;}
  while (!wkspace_ua) { 
    if (GLOBAL_verbose>2){ if (!wkspace_ua){ printf("Could not allocate %ld MB; trying %ld MB instead.\n",malloc_size_mb,malloc_size_mb - 64);}}
    if (malloc_size_mb > 128) { malloc_size_mb -= 64;} else { malloc_size_mb = 64;}
    wkspace_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
    if (wkspace_ua) { printf("Allocated %ld MB successfully.\n", malloc_size_mb);}
    /* while (!wkspace_ua) { } */}
  // force 64-byte align on OS X to make cache line sensitivity work
  wkspace = (unsigned char*)CACHEALIGN((unsigned long)wkspace_ua);
  wkspace_base = wkspace;
  wkspace_left = malloc_size_mb * 1048576 - (unsigned long long int)(wkspace - wkspace_ua);
  free(bubble);
  if (GLOBAL_verbose>0){ wkspace_printf();}
  wkspace_point_0 = wkspace_make_point();
  wkspace_point_0->parent = NULL; wkspace_point_0->child = NULL;
  wkspace_point_0->check = 0; *(wkspace_point_0->point) = wkspace_point_0->check;
  wkspace_point_t = wkspace_point_0;
  GLOBAL_verbose=0;
  read_input(); if (GLOBAL_verbose>1){ printf("sysconf(_SC_NPROCESSORS_ONLN) = %d\n",(int)sysconf(_SC_NPROCESSORS_ONLN));}
  set_globals();
  if (strcmp(GLOBAL_TEST_TYPE,"wrap_M_setup_test")==0){ wrap_M_setup_test_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"An_ajdk_v")==0){ wrap_An_ajdk_v_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"An_v")==0){ wrap_An_v_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"AnZt_vv")==0){ wrap_AnZt_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"AnAt_vv")==0){ wrap_AnAt_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"AtTYn_vv")==0){ wrap_AtTYn_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"AtTAn_vv")==0){ wrap_AtTAn_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"An_ZtSWn_Yt_vv")==0){ wrap_An_ZtSWn_Yt_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"AnZt_S_WnYt_vv")==0){ wrap_AnZt_S_WnYt_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"An_ZtSWn_ww")==0){ wrap_An_ZtSWn_ww_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"An_ZtSWn_Yt_ww")==0){ wrap_An_ZtSWn_Yt_ww_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"At_T_YnWt_S_Zn_vv")==0){ wrap_At_T_YnWt_S_Zn_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"AtTYn____WtSZn_vv")==0){ wrap_AtTYn____WtSZn_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"AttYn____WtsZn_vv")==0){ wrap_AttYn____WtsZn_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"At_T_YnWt_ww")==0){ wrap_At_T_YnWt_ww_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"At_T_YnWt_S_Zn_ww")==0){ wrap_At_T_YnWt_S_Zn_ww_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"An_ZtSn_ww")==0){ wrap_An_ZtSn_ww_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"bcc_init")==0){ bcc_init_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"bcc_lf_AnZtSWnYt")==0){ bcc_lf_AnZtSWnYt_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"bcc_lf_AtTYnWtSZn")==0){ bcc_lf_AtTYnWtSZn_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"bcc_lrup")==0){ bcc_lrup_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"bcc_flattenloop")==0){ bcc_flattenloop_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"bcc_lrup_flattenloop")==0){ bcc_lrup_flattenloop_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"bcc_sumscores")==0){ bcc_sumscores_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"bcc_lrup_sumscores")==0){ bcc_lrup_sumscores_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"lakcluster_driver")==0){ lakcluster_driver();}
  if (strcmp(GLOBAL_TEST_TYPE,"lakcluster_scorebox")==0){ lakcluster_scorebox();}
  if (strcmp(GLOBAL_TEST_TYPE,"D_AtTn_ZtSn_vv")==0){ wrap_D_AtTn_ZtSn_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcc_init")==0){ dcc_init_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcc_lf_D_AtTn_ZtSn")==0){ dcc_lf_D_AtTn_ZtSn_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcc_lf_TAnZtS")==0){ dcc_lf_TAnZtS_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcc_sumscores")==0){ dcc_sumscores_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcc_time_sumscores")==0){ dcc_time_sumscores_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dexcluster_driver")==0){ dexcluster_driver();}
  free(wkspace_ua); exit(RET_SUCCESS);
  return 0;
}
