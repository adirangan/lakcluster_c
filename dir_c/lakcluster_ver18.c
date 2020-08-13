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

#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */
#ifdef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

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
const char * SPACING_name[] = { "SPACING_b","SPACING_j","SPACING_a" };
int addressable_spacing_b=SPACING_b;
int addressable_spacing_j=SPACING_j;
int addressable_spacing_a=SPACING_a;

/* global variables used for timing */
clock_t GLOBAL_t_start[GLOBAL_NTICKS],GLOBAL_t_final[GLOBAL_NTICKS];
struct timeval GLOBAL_d_start[GLOBAL_NTICKS],GLOBAL_d_final[GLOBAL_NTICKS];
long GLOBAL_l_msec[GLOBAL_NTICKS],GLOBAL_l_ssec[GLOBAL_NTICKS],GLOBAL_l_usec[GLOBAL_NTICKS]; 
double GLOBAL_elct[GLOBAL_NTICKS],GLOBAL_elrt[GLOBAL_NTICKS];

/* These are the kappa-values used for continuous-covariate correction; note that the former are the square of the latter. */
double GLOBAL_kappa_squared_loop_scale_factor_[6] = {1.00 , 0.34 , 0.19 , 0.13 , 0.10 , 0.09};
double GLOBAL_kappa_squared_half_scale_factor_[6] = {1.0000 , 0.5747 , 0.4338 , 0.3613 , 0.3154 , 0.3000 } ; 
double GLOBAL_kappa_squared = 0.0;

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
char GLOBAL_A_p_name[FNAMESIZE]="\0";
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
char GLOBAL_Y_p_name[FNAMESIZE]="\0";
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
char GLOBAL_J_n_rind[FNAMESIZE]="\0";
char **GLOBAL_J_n_rind_=NULL;
char GLOBAL_K_n_cind[FNAMESIZE]="\0";
char GLOBAL_out_name[FNAMESIZE]="\0";
char GLOBAL_scorebox_out_xdrop[FNAMESIZE]="\0";
int GLOBAL_scorebox_row_max=0;
int GLOBAL_scorebox_row_num=0;
int GLOBAL_scorebox_col_max=0;
int GLOBAL_scorebox_col_num=0;
char GLOBAL_pca_out_xdrop[FNAMESIZE]="\0";
char GLOBAL_pca_V_[FNAMESIZE]="\0";
char GLOBAL_pca_infix[FNAMESIZE]="\0";
int GLOBAL_pca_iteration_num=0;
int GLOBAL_pca_iteration_max=0;
int GLOBAL_pca_iteration_min=0;
int GLOBAL_pca_rank=0;
double GLOBAL_pca_tolerance=0.000001;
int GLOBAL_scramble_num=0;
char **GLOBAL_scramble_out_xdrop_=NULL;
unsigned long int *GLOBAL_scramble_rseed_=NULL;
/* bed_to_b16 */
char GLOBAL_fname_bed_0in[FNAMESIZE]="\0";
char GLOBAL_fname_b16_out[FNAMESIZE]="\0";
char GLOBAL_fname_bim_0in[FNAMESIZE]="\0";
char GLOBAL_fname_bim_out[FNAMESIZE]="\0";
char GLOBAL_fname_fam_0in[FNAMESIZE]="\0";
char GLOBAL_fname_fam_out[FNAMESIZE]="\0";
char GLOBAL_fname_flip_flag[FNAMESIZE]="\0";
double GLOBAL_snp_mss_threshold=0.01;
double GLOBAL_snp_maf_threshold=0.25;
double GLOBAL_snp_I_opt_threshold=0.01;
double GLOBAL_pat_mss_threshold=0.04;
int GLOBAL_n_fam_char_max=64;
int GLOBAL_n_bim_char_max=64;
/* bed_merge */
int GLOBAL_n_file=0;
char **GLOBAL_fname_b16_0in_=NULL;
char **GLOBAL_fname_bim_0in_=NULL;
char **GLOBAL_fname_fam_0in_=NULL;
/* wrap_transpose */
char GLOBAL_fname_b16_0in[FNAMESIZE]="\0";
int GLOBAL_n_bytes_per_read=0;

/* popcount table */
unsigned char popcount_uchar[256] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
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

/* ---------------------------------------------------------------- */

/* manually managed memory stack */
unsigned char* wkspace=NULL;
unsigned char* wkspace_base=NULL;
int GLOBAL_wkspace_point=1;
struct wkspace_point *wkspace_point_0=NULL,*wkspace_point_t=NULL;
unsigned long long int GLOBAL_memory_gb=GLOBAL_MEMORY_GB_DEFAULT;
unsigned long long int GLOBAL_memory_mb=GLOBAL_MEMORY_GB_DEFAULT*(unsigned long long int)1024;
unsigned long long int GLOBAL_memory_kb=GLOBAL_MEMORY_GB_DEFAULT*(unsigned long long int)1048576;
long long int wkspace_left=0;
long long int wkspace_used=0;

/* ---------------------------------------------------------------- */

/* RAND functions */
unsigned long int POW2RPOWPLUSRADD=35L;
unsigned long int POW22RPOWMINUSONE=2147483647LL;
int RCYCLENUM=7;

/*---------------------------------------------------------------- */

#ifdef _MONOLITH
#include "lakcluster_function.c"
#endif /* _MONOLITH */

inline void ping(){ printf(" %% ping\n");}
inline void pong(){ printf(" %% pong\n");}

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
    else if (!strcmp(argptr, "--GLOBAL_memory_gb")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --GLOBAL_memory_gb parameter."); exit(RET_INVALID_CMDLINE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid --GLOBAL_memory_gb parameter.\n"); exit(RET_INVALID_CMDLINE);}
      GLOBAL_memory_gb = ii; printf(" %% GLOBAL_memory_gb %d\n",GLOBAL_memory_gb);
      cur_arg += 2; /* else if GLOBAL_memory_gb */}
    else if (!strcmp(argptr, "--nthreads")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --nthreads parameter."); exit(RET_INVALID_CMDLINE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid --nthreads parameter.\n"); exit(RET_INVALID_CMDLINE);}
      GLOBAL_thread_count = ii; printf(" %% nthreads %d\n",GLOBAL_thread_count);
      cur_arg += 2; /* else if nthreads */}
    else if (!strcmp(argptr, "--cmdline_param")) {
      if (cur_arg == argc - 1) { printf("Error: Missing --cmdline_param parameter."); exit(RET_INVALID_CMDLINE);}
      ii = atoi(argv[cur_arg + 1]);
      if (ii < 1) { printf("Error: Invalid parameter.\n"); exit(RET_INVALID_CMDLINE);}
      cmdline_param = ii; printf(" %% cmdline_param %d\n",cmdline_param);
      cur_arg += 2; /* else if other */}
    else { printf("Error: Invalid argument (%s).", argv[cur_arg]); exit(RET_INVALID_CMDLINE);}
    /* while (cur_arg < argc) { } */}
  bubble = (char*)malloc((unsigned long long int)67108864 * sizeof(char));
  if (!bubble) { printf("Error: not enough memory!\n");}
  GLOBAL_memory_mb = GLOBAL_memory_gb * (unsigned long long int)1024; GLOBAL_memory_kb = GLOBAL_memory_mb * (unsigned long long int)1024;
  wkspace_ua = (unsigned char*)malloc(GLOBAL_memory_kb * (unsigned long long int)1024 * sizeof(char));
  while (!wkspace_ua) { 
    if (GLOBAL_verbose>-2){ if (!wkspace_ua){ printf("Could not allocate %ld GB; trying %ld GB instead...\n",GLOBAL_memory_gb,GLOBAL_memory_gb - 1);}}
    GLOBAL_memory_gb -= 1; GLOBAL_memory_mb = GLOBAL_memory_gb * (unsigned long long int)1024; GLOBAL_memory_kb = GLOBAL_memory_mb * (unsigned long long int)1024;
    wkspace_ua = (unsigned char*)malloc(GLOBAL_memory_kb * (unsigned long long int)1024 * sizeof(char));
    if (wkspace_ua) { printf("Allocated %ld KB = %ld MB = %ld GB successfully.\n", GLOBAL_memory_kb,GLOBAL_memory_mb,GLOBAL_memory_gb);}
    /* while (!wkspace_ua) { } */}
  // force 64-byte align on OS X to make cache line sensitivity work
  wkspace = (unsigned char*)CACHEALIGN((unsigned long)wkspace_ua);
  wkspace_base = wkspace;
  wkspace_left = GLOBAL_memory_kb * (unsigned long long int)1024 - (unsigned long long int)(wkspace - wkspace_ua);
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
  if (strcmp(GLOBAL_TEST_TYPE,"M_Ax_to_L2")==0){ wrap_M_Ax_to_L2_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"An_ajdk_v")==0){ wrap_An_ajdk_v_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"An_v")==0){ wrap_An_v_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"AnZt_vv")==0){ wrap_AnZt_vv_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"AnZt_mm")==0){ wrap_AnZt_mm_test();}
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
  if (strcmp(GLOBAL_TEST_TYPE,"dexcluster_scorebox")==0){ dexcluster_scorebox();}
  if (strcmp(GLOBAL_TEST_TYPE,"At_T_Xn_ww")==0){ wrap_At_T_Xn_ww_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"An_Xn_ww")==0){ wrap_An_Xn_ww_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"pca_init")==0){ pca_init_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"pca_iter")==0){ pca_iter_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"pca_driver")==0){ pca_driver();}
  if (strcmp(GLOBAL_TEST_TYPE,"pca_proj_test")==0){ pca_proj_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"pca_proj_driver")==0){ pca_proj_driver();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcg_init")==0){ dcg_init_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcg_lf_D_AtTn_ZtSn")==0){ dcg_lf_D_AtTn_ZtSn_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcg_lf_TAnZtS")==0){ dcg_lf_TAnZtS_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcg_sumscores")==0){ dcg_sumscores_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcg_time_sumscores")==0){ dcg_time_sumscores_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"dcgxpander_driver")==0){ dcgxpander_driver();}
  if (strcmp(GLOBAL_TEST_TYPE,"bed_to_b16_test")==0){ bed_to_b16_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"bed_to_b16")==0){ bed_to_b16();}  
  if (strcmp(GLOBAL_TEST_TYPE,"MDA_io_test")==0){ MDA_io_test();}  
  if (strcmp(GLOBAL_TEST_TYPE,"b16_merge_test")==0){ b16_merge_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"b16_merge")==0){ b16_merge();}
  if (strcmp(GLOBAL_TEST_TYPE,"wrap_transpose_test")==0){ wrap_transpose_test();}
  if (strcmp(GLOBAL_TEST_TYPE,"wrap_transpose_driver")==0){ wrap_transpose_driver();}
  free(wkspace_ua); exit(RET_SUCCESS);
  return 0;
}
