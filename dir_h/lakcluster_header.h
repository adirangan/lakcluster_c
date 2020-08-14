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
#include <cblas.h>

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
#define FNAMESIZE (4096+256)
#define GLOBAL_MEMORY_GB_DEFAULT 31
#define MULTIPLEX_DIST 960
#define MULTIPLEX_2DIST (MULTIPLEX_DIST * 2)
#define BIT8 8
#define BITJ 16
#define SNP_AND_TAG 3
#define SNP_XOR_TAG 2
#define SNP_NOR_TAG 0
#define SNP_MSS_TAG 1
#define SNP_NOT_TAG 5

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
extern const char *TYPE_name[9];
extern int addressable_type_p0;
extern int addressable_type_pm;
extern int addressable_type_00;
extern int addressable_type_pz;
extern int addressable_type_0z;
extern int addressable_type_p_;
extern int addressable_type_0_;
extern int addressable_type_pZ;
extern int addressable_type_0Z;
#define SPACING_b 0
#define SPACING_j 1
#define SPACING_a 2
extern const char *SPACING_name[3];
extern int addressable_spacing_b;
extern int addressable_spacing_j;
extern int addressable_spacing_a;
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
extern clock_t GLOBAL_t_start[GLOBAL_NTICKS];
extern clock_t GLOBAL_t_final[GLOBAL_NTICKS];
extern struct timeval GLOBAL_d_start[GLOBAL_NTICKS];
extern struct timeval GLOBAL_d_final[GLOBAL_NTICKS];
extern long GLOBAL_l_msec[GLOBAL_NTICKS],GLOBAL_l_ssec[GLOBAL_NTICKS],GLOBAL_l_usec[GLOBAL_NTICKS]; 
extern double GLOBAL_elct[GLOBAL_NTICKS],GLOBAL_elrt[GLOBAL_NTICKS];

/* These are the kappa-values used for continuous-covariate correction; note that the former are the square of the latter. */
extern double GLOBAL_kappa_squared_loop_scale_factor_[6];
extern double GLOBAL_kappa_squared_half_scale_factor_[6];
extern double GLOBAL_kappa_squared;

/* global variables used for most routines */
extern char GLOBAL_CWD[FNAMESIZE];
extern char GLOBAL_DIR_BASE[FNAMESIZE];
extern char GLOBAL_DIR_NAME[FNAMESIZE];
extern char GLOBAL_DIR_XPRE[FNAMESIZE];
extern int GLOBAL_verbose;
extern int GLOBAL_thread_count;
extern int GLOBAL_omp_type;
extern int GLOBAL_1_icache_linesize;
extern int GLOBAL_1_icache_size;
extern int GLOBAL_1_icache_assoc;
extern int GLOBAL_1_dcache_linesize;
extern int GLOBAL_1_dcache_size;
extern int GLOBAL_1_dcache_assoc;
extern int GLOBAL_2_cache_linesize;
extern int GLOBAL_2_cache_size;
extern int GLOBAL_2_cache_assoc;
extern int GLOBAL_3_cache_linesize;
extern int GLOBAL_3_cache_size;
extern int GLOBAL_3_cache_assoc;
extern int GLOBAL_4_cache_linesize;
extern int GLOBAL_4_cache_size;
extern int GLOBAL_4_cache_assoc;
extern double GLOBAL_tolerance;
extern unsigned int GLOBAL_recursion_limit;
extern int GLOBAL_LBITS;
extern char GLOBAL_TEST_TYPE[FNAMESIZE];
extern char GLOBAL_TEST_TYP2[FNAMESIZE];
extern int GLOBAL_TEST_sparse;
extern int GLOBAL_TEST_A_n_rows;
extern int GLOBAL_TEST_A_n_cols;
extern int GLOBAL_TEST_Z_n_rows;
extern int GLOBAL_TEST_Y_n_cols;
extern int GLOBAL_TEST_T_n_cols;
extern int GLOBAL_TEST_niter;
extern double GLOBAL_TEST_mrand;
extern char GLOBAL_skip[FNAMESIZE];
extern char GLOBAL_QR_strategy[FNAMESIZE];
extern char GLOBAL_QC_strategy[FNAMESIZE];
extern long long int GLOBAL_D_MLT;
extern int GLOBAL_B_MLT;
extern int GLOBAL_gamma_type;
extern double GLOBAL_gamma;
extern int addressable_1;
extern int addressable_0;
extern int addressable_int_length;
extern int addressable_int[128];
extern char GLOBAL_A_n_name[FNAMESIZE];
extern char GLOBAL_A_t_name[FNAMESIZE];
extern char **GLOBAL_A_n_name_;
extern char **GLOBAL_A_t_name_;
extern int  *GLOBAL_A_n_rows_;
extern int  GLOBAL_A_n_cols;
extern char GLOBAL_A_n_rind[FNAMESIZE];
extern char **GLOBAL_A_n_rind_;
extern char GLOBAL_A_n_cind[FNAMESIZE];
extern char GLOBAL_A_p_name[FNAMESIZE];
extern char GLOBAL_Z_n_name[FNAMESIZE];
extern char GLOBAL_Z_t_name[FNAMESIZE];
extern char **GLOBAL_Z_n_name_;
extern char **GLOBAL_Z_t_name_;
extern int  *GLOBAL_Z_n_rows_;
extern char GLOBAL_Z_n_rind[FNAMESIZE];
extern char **GLOBAL_Z_n_rind_;
extern char GLOBAL_Y_n_name[FNAMESIZE];
extern char GLOBAL_Y_t_name[FNAMESIZE];
extern char **GLOBAL_Y_n_name_;
extern char **GLOBAL_Y_t_name_;
extern int  GLOBAL_Y_n_cols;
extern char GLOBAL_Y_n_cind[FNAMESIZE];
extern char GLOBAL_Y_p_name[FNAMESIZE];
extern char GLOBAL_W_n_name[FNAMESIZE];
extern char GLOBAL_W_t_name[FNAMESIZE];
extern char **GLOBAL_W_n_name_;
extern char **GLOBAL_W_t_name_;
extern char GLOBAL_T_n_name[FNAMESIZE];
extern char GLOBAL_T_t_name[FNAMESIZE];
extern char **GLOBAL_T_n_name_;
extern char **GLOBAL_T_t_name_;
extern int  GLOBAL_T_n_cols;
extern char GLOBAL_T_n_cind[FNAMESIZE];
extern char GLOBAL_S_n_name[FNAMESIZE];
extern char GLOBAL_S_t_name[FNAMESIZE];
extern char **GLOBAL_S_n_name_;
extern char **GLOBAL_S_t_name_;
extern char GLOBAL_out_name[FNAMESIZE];
extern int GLOBAL_scramble_num;
extern char **GLOBAL_scramble_out_xdrop_;
extern unsigned long int *GLOBAL_scramble_rseed_;
extern char GLOBAL_scorebox_out_xdrop[FNAMESIZE];
extern int GLOBAL_scorebox_row_max;
extern int GLOBAL_scorebox_row_num;
extern int GLOBAL_scorebox_col_max;
extern int GLOBAL_scorebox_col_num;
extern char GLOBAL_pca_out_xdrop[FNAMESIZE];
extern char GLOBAL_pca_V_[FNAMESIZE];
extern char GLOBAL_pca_infix[FNAMESIZE];
extern int GLOBAL_pca_iteration_num;
extern int GLOBAL_pca_iteration_max;
extern int GLOBAL_pca_iteration_min;
extern int GLOBAL_pca_rank;
extern double GLOBAL_pca_tolerance;
extern char GLOBAL_J_n_rind[FNAMESIZE];
extern char **GLOBAL_J_n_rind_;
extern char GLOBAL_K_n_cind[FNAMESIZE];
/* bed_to_b16 */
extern char GLOBAL_fname_bed_0in[FNAMESIZE];
extern char GLOBAL_fname_b16_out[FNAMESIZE];
extern char GLOBAL_fname_bim_0in[FNAMESIZE];
extern char GLOBAL_fname_bim_out[FNAMESIZE];
extern char GLOBAL_fname_fam_0in[FNAMESIZE];
extern char GLOBAL_fname_fam_out[FNAMESIZE];
extern char GLOBAL_fname_flip_flag[FNAMESIZE];
extern double GLOBAL_snp_mss_threshold;
extern double GLOBAL_snp_maf_threshold;
extern double GLOBAL_snp_I_opt_threshold;
extern double GLOBAL_pat_mss_threshold;
extern int GLOBAL_n_fam_char_max;
extern int GLOBAL_n_bim_char_max;
/* bed_merge */
extern int GLOBAL_n_file;
extern char **GLOBAL_fname_b16_0in_;
extern char **GLOBAL_fname_bim_0in_;
extern char **GLOBAL_fname_fam_0in_;
/* wrap_transpose */
extern char GLOBAL_fname_b16_0in[FNAMESIZE];
extern int GLOBAL_n_bytes_per_read;

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

extern unsigned char popcount_uchar[256];

#include "S_handle.h"
#include "L_handle.h"
#include "M_handle.h"
#include "P_handle.h"
#include "R_handle.h"
#include "bcc.h"
#include "dcc.h"
#include "dcg.h"

/* thread management */
extern int GLOBAL_NBINS;
extern int GLOBAL_Ireq;
extern int *GLOBAL_lnb;
extern int GLOBAL_nf; extern int GLOBAL_nf_cur; extern int GLOBAL_nf_ind; extern int GLOBAL_nf_opn;
extern int GLOBAL_tint[MAX_THREADS];
extern void *GLOBAL_tvp[MAX_THREADS][128];
extern pthread_t GLOBAL_threads[MAX_THREADS];
extern unsigned long long int GLOBAL_ops_f_[MAX_THREADS];
extern unsigned long long int GLOBAL_ops_f_sum;
extern unsigned long long int GLOBAL_ops_b_[MAX_THREADS];
extern unsigned long long int GLOBAL_ops_b_sum;

/* ---------------------------------------------------------------- */

/* manually managed memory stack */
extern unsigned char* wkspace;
extern unsigned char* wkspace_base;
extern int GLOBAL_wkspace_point;
struct wkspace_point
{
  struct wkspace_point * parent;
  struct wkspace_point * child;
  long long int check;
  long long int *point;
};
extern struct wkspace_point *wkspace_point_0,*wkspace_point_t;
extern unsigned long long int GLOBAL_memory_kb;
extern unsigned long long int GLOBAL_memory_mb;
extern unsigned long long int GLOBAL_memory_gb;
extern long long int wkspace_left;
extern long long int wkspace_used;

/* ---------------------------------------------------------------- */

/* RAND functions */
extern unsigned long int POW2RPOWPLUSRADD;
extern unsigned long int POW22RPOWMINUSONE;
extern int RCYCLENUM;
unsigned long int lrand();
double randn();
unsigned long int RGET(unsigned long int *);
double R01GET(unsigned long int *);
double RNGET(unsigned long int *);
double RISIGET(unsigned long int *,double);

/* ---------------------------------------------------------------- */


#ifndef _MONOLITH
void fill_uchar_zero(unsigned char* iarr, size_t size);
void fill_uchar_ones(unsigned char* iarr, size_t size);
void fill_long_zero(long* larr, size_t size);
#endif /* _MONOLITH */
void ping();
void pong();
void calc_A_ajdk(double *A_p,int A_pcols,double **A_ajdk_p);
void wrap_M_setup_test_excerpt_0(char *error_vs_speed,double mrnd,int nrows,int rpop_b,int ncols,int cpop_b,unsigned char *bn,unsigned char *bt);
void wrap_M_setup_test_excerpt_1(char *error_vs_speed,double mrnd,int nrows,int ncols,unsigned char *mr_b,unsigned char *mc_b,unsigned char **bn_ra_p,unsigned char **bt_ra_p);
void wrap_M_setup_test_excerpt_2(int nrows,int ncols,unsigned char *bn,unsigned char *bt,unsigned char *mr_b,unsigned char *mr_j,unsigned char *mc_b,unsigned char *mc_j,struct M_handle **Mn_p,struct M_handle **Mt_p);
void wrap_M_setup_test(char *error_vs_speed,double mrnd,int nbins,int *A_n_rows,int A_n_cols,int *Z_n_rows,int Y_n_cols,int T_n_cols,struct M_handle ***M_An_p,struct M_handle ***M_At_p,struct M_handle ***M_Zn_p,struct M_handle ***M_Zt_p,struct M_handle ***M_Yn_p,struct M_handle ***M_Yt_p,struct M_handle ***M_Wn_p,struct M_handle ***M_Wt_p,struct M_handle ***M_Tn_p,struct M_handle ***M_Tt_p,struct M_handle ***M_Sn_p,struct M_handle ***M_St_p,double **A_p_p,double **A_ajdk_p,struct L_handle ***lf_An_ajdk_p,struct L_handle ***lf_Zn_ajdk_p,double **Y_p_p,double **Y_ajdk_p,struct L_handle ***lf_Yn_ajdk_p,struct L_handle ***lf_Wn_ajdk_p);
void wrap_M_setup_test_test();
void *get_An_ajdk_v(void *vp);
int wrap_An_ajdk_v(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,struct M_handle *M_An,double *A_ajdk,struct L_handle **output_An_ajdk_p);
void *get_An_ajdk_u(void *vp);
int wrap_An_ajdk_u(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,struct M_handle *M_An,double *A_ajdk,struct L_handle **output_An_ajdk_p);
void wrap_An_ajdk_v_test();
int wrap_AnAt_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,int output_spacing_c,struct M_handle *M_An,struct M_handle *M_Zn,double *A_ajdk,struct L_handle *lf_An_ajdk,struct L_handle *lf_Zn_ajdk,struct L_handle **output_AnAt_p);
void wrap_AnAt_vv_test();
int wrap_An_v__run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,struct M_handle *M_An,struct L_handle **output_An_v_p);
void *get_An_u__run(void *vp);
int wrap_An_u__run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,struct M_handle *M_An,struct L_handle **output_An_u_p);
void wrap_An_v_test();
void *get_AnZt_mm__run(void *vp);
int wrap_AnZt_mm__run(int *tidx,void **vpra,pthread_t *thread_in,struct L_handle *lf_An,struct L_handle *lf_Zn,struct L_handle **output_AnZt_p);
void wrap_AnZt_mm_test();
int wrap_An_ZtSn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,struct M_handle *M_An,struct L_handle *lf_ZtSn,struct M_handle *M_ZtSn,struct M_handle *M_St,struct L_handle **output_An_ZtSn_ww_p);
void *get_An_ZtSn_uu(void *vp);
int wrap_An_ZtSn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,struct M_handle *M_An,struct L_handle *lf_ZtSn,struct M_handle *M_St,struct L_handle **output_An_ZtSn_uu_p);
void An_a0d1_ZtSn_excerpt(struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle *lf_ZtSn,struct L_handle *lf_a0d1_ZtSn);
void An_a2d2_ZtSn_excerpt(struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle *lf_ZtSn,struct L_handle *lf_a2d2_ZtSn);
void wrap_An_ZtSn_ww_test();
int wrap_An_ZtSWn_Yt_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_St,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle *lf_ZtSWn,struct L_handle **output_p);
void *get_AnZtSWnYt_uu(void *vp);
int wrap_AnZtSWnYt_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_Zt,struct M_handle *M_St,struct M_handle *M_Wt,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle **output_p);
void wrap_An_ZtSWn_Yt_vv_test();
int wrap_AnZt_S_WnYt_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_St,struct M_handle *M_Zn,struct L_handle *lf_AnZt,struct L_handle *lf_YnWt,struct L_handle **output_p);
void wrap_AnZt_S_WnYt_vv_test();
int wrap_An_ZtSWn_Yt_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,struct M_handle *M_An,struct M_handle *M_St,struct M_handle **M_ZtSWn_,struct M_handle *M_Wt,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle **output_An_ZtSWn_Yt_ww_p);
void wrap_An_ZtSWn_Yt_ww_test();
int wrap_AnZt_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,int output_spacing_c,struct M_handle *M_An,struct M_handle *M_Zn,double *A_ajdk,struct L_handle *lf_An_ajdk,struct L_handle *lf_Zn_ajdk,struct L_handle **output_AnZt_p);
void *get_AnZt_uu__run(void *vp);
int wrap_AnZt_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,int output_spacing_c,struct M_handle *M_An,struct M_handle *M_Zn,double *A_ajdk,struct L_handle **output_AnZt_p);
void wrap_AnZt_vv_test();
int wrap_AttAn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_y,struct M_handle *M_At,struct M_handle *M_Tt,int ns_j,struct M_handle *M_Yt,double *A_ajdk,double *Y_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_YtTn,struct L_handle **output_AttAn_p);
int wrap_AtTAn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_y,int output_spacing_t,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yt,double *A_ajdk,double *Y_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_YtTn,struct L_handle **output_AtTAn_p);
void wrap_AtTAn_vv_test();
int wrap_AttYn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_y,struct M_handle *M_At,struct M_handle *M_Tt,int ns_j,struct M_handle *M_Yt,double *A_ajdk,double *Y_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_YtTn,struct L_handle **output_AttYn_p);
int wrap_AtTYn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_y,int output_spacing_t,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yt,double *A_ajdk,double *Y_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_YtTn,struct L_handle **output_AtTYn_p);
void *get_AtTYn_uu(void *vp);
int wrap_AtTYn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_y,int output_spacing_t,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yt,double *A_ajdk,double *Y_ajdk,struct L_handle **output_AtTYn_p);
void wrap_AtTYn_vv_test();
int wrap_AttYn____WtsZn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,int ns_j,struct M_handle *M_Yt,struct M_handle *M_Wt,struct M_handle *M_St,struct M_handle *M_Zt,double *Y_ajdk,struct L_handle *lf_AttYn,struct L_handle *lf_ZtsWn,struct L_handle **output_p);
void wrap_AttYn____WtsZn_vv_test();
int wrap_At_T_YnWt_S_Zn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,struct L_handle *lf_YnWt,struct L_handle **output_p);
void *get_AtTYnWtSZn_uu(void *vp);
int wrap_AtTYnWtSZn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yn,struct M_handle *M_Wn,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,double *Y_ajdk,struct L_handle **output_p);
void wrap_At_T_YnWt_S_Zn_vv_test();
int wrap_AtTYn____WtSZn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yt,struct M_handle *M_Wt,struct M_handle *M_St,struct M_handle *M_Zt,double *Y_ajdk,struct L_handle *lf_AtTYn,struct L_handle *lf_ZtSWn,struct L_handle **output_p);
void wrap_AtTYn____WtSZn_vv_test();
void *get_At_T_YnWt_S_Zn_ww_precompute(void *vp);
int wrap_At_T_YnWt_S_Zn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_t,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yn,struct M_handle *M_Wn,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,struct M_handle *M_YnWt,struct L_handle **output_At_T_YnWt_S_Zn_ww_p);
void wrap_At_T_YnWt_S_Zn_ww_test();
void bcc_ajdk_copy(struct bcc_ajdk *D,struct bcc_ajdk *D_in);
void bcc_ajdk_load(struct bcc_ajdk *D);
void bcc_ajdk_init(double mrnd,struct bcc_ajdk *D);
void bcc_single_copy_M_An(struct bcc_single *E,struct bcc_single *E_in);
void bcc_single_load_M_An(struct bcc_single *E);
void bcc_single_init_M_An(char *error_vs_speed,double mrnd,struct bcc_single *E);
void bcc_single_copy_lf(struct bcc_single *E,struct bcc_single *E_in);
void bcc_single_init_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,struct M_handle **M_p);
void bcc_single_init_lf(struct bcc_single *E);
void bcc_double_copy_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,int *trm_flag_p,struct M_handle **M_p,struct L_handle **L_in_p,struct M_handle **M_in_p);
void bcc_double_copy_lf(struct bcc_double *F,struct bcc_double *F_in);
void bcc_double_init_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,int *trm_flag_p,struct M_handle **M_p);
void bcc_double_init_lf(struct bcc_double *F);
void bcc_copy_A_p(struct bcc_ajdk *D,struct bcc_ajdk *D_in);
void bcc_X_nrows_total(struct bcc_ajdk *D);
void bcc_load_A_p(struct bcc_ajdk *D);
void bcc_init_A_p(double mrnd,struct bcc_ajdk *D);
void bcc_copy_QX(struct bcc_ajdk *D,struct bcc_ajdk *D_in);
void bcc_init_QX(struct bcc_ajdk *D);
void bcc_M_mxset(struct bcc_ajdk *D);
void bcc_copy(struct bcc_ajdk *D,struct bcc_ajdk *D_in);
void bcc_load(struct bcc_ajdk **D_p,struct bcc_single ***E_p,struct bcc_double ***F_p,char *QR_strategy,char *QC_strategy);
void bcc_init(char *error_vs_speed,double mrnd,struct bcc_ajdk **D_p,struct bcc_single ***E_p,struct bcc_double ***F_p,char *QR_strategy,char *QC_strategy);
void bcc_init_test();
void *get_singlestudy_uu(void *vp);
void wrap_singlestudy_uu(int *tidx,void **vpra,pthread_t *thread_in,struct bcc_single *E);
void bcc_singlestudy_uu(struct bcc_ajdk *D);
void *get_singlestudy_vv(void *vp);
void wrap_singlestudy_vv(int *tidx,void **vpra,pthread_t *thread_in,struct bcc_single *E);
void At_Yn_a1d1_excerpt(struct M_handle *M_Yn,double *Y_ajdk,struct L_handle *lf_Yn_ajdk,struct L_handle *lf_Yn_a1d1);
void bcc_singlestudy_ww(struct bcc_ajdk *D);
void *get_doublestudy_uu(void *vp);
void wrap_doublestudy_uu(int *tidx,void **vpra,pthread_t *thread_in,struct bcc_double *F);
void bcc_doublestudy_uu(struct bcc_ajdk *D);
void bcc_doublestudy_ww(struct bcc_ajdk *D);
void *get_QR_AnZtSWnYt_uu(void *vp);
int wrap_QR_AnZtSWnYt_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,int AZ_flag,int AY_flag,struct M_handle *M_An,struct M_handle *M_Zt,struct M_handle *M_St,struct M_handle *M_Wt,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle **output_p);
void bcc_QR_AnZtSWnYt_uu(struct bcc_ajdk *D);
void *get_QC_AtTYnWtSZn_uu(void *vp);
int wrap_QC_AtTYnWtSZn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,int AZ_flag,int AY_flag,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yn,struct M_handle *M_Wn,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,double *Y_ajdk,struct L_handle **output_p);
void bcc_QC_AtTYnWtSZn_uu(struct bcc_ajdk *D);
void *get_flattenloop(void *vp);
void wrap_flattenloop(int *tidx,void **vpra,pthread_t *thread_in,struct bcc_double *F);
void bcc_flattenloop(struct bcc_ajdk *D);
void bcc_QR_AnZtSWnYt_error(int verbose,struct bcc_ajdk *D);
void bcc_QC_AtTYnWtSZn_error(int verbose,struct bcc_ajdk *D);
void bcc_flattenloop_test();
void bcc_An_ajdk(struct bcc_ajdk *D);
void bcc_lf_ZtSn(struct bcc_ajdk *D);
void bcc_lf_ZtSWn(struct bcc_ajdk *D);
void bcc_M_ZtSWn_excerpt(int trm_flag,struct M_handle *M_St,struct M_handle *M_in1,struct M_handle *M_in2,struct L_handle *L_in,struct M_handle **M_out_);
void bcc_M_ZtSWn(struct bcc_ajdk *D);
void bcc_lf_YnWt(struct bcc_ajdk *D);
void bcc_M_YnWt_excerpt(int trm_flag,struct M_handle *M_in1,struct M_handle *M_in2,struct L_handle *L_in,struct M_handle **M_out_p);
void bcc_M_YnWt(struct bcc_ajdk *D);
void bcc_lf_AnZt_S_WnYt(struct bcc_ajdk *D);
void bcc_lf_An_ZtSWn_Yt_excerpt(int trm_flag,int n_spacing_A,struct M_handle *M_An,struct M_handle *M_St,struct M_handle **M_ZtSWn_,struct M_handle *M_Zt,struct M_handle *M_Wt,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle **L_out_p);
void bcc_lf_An_ZtSWn_Yt(struct bcc_ajdk *D);
void bcc_lf_AnZtSWnYt(struct bcc_ajdk *D);
void bcc_lf_AnZtSWnYt_error(int verbose,struct bcc_ajdk *D);
void bcc_lf_AnZtSWnYt_test();
void bcc_lf_AtTYn____WtSZn(struct bcc_ajdk *D);
void bcc_lf_At_T_YnWt_S_Zn_excerpt(int trm_flag,int n_spacing_A,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yn,struct M_handle *M_Wn,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,struct M_handle *M_YnWt,struct L_handle **L_out_p);
void bcc_lf_At_T_YnWt_S_Zn(struct bcc_ajdk *D);
void bcc_lf_AtTYnWtSZn(struct bcc_ajdk *D);
void bcc_lf_AtTYnWtSZn_error(int verbose,struct bcc_ajdk *D);
void bcc_lf_AtTYnWtSZn_test();
void bcc_lrup_QC_YnWt_stage_a0(struct bcc_ajdk *D);
void bcc_lrup_QC_YnWt_stage_a1(struct bcc_ajdk *D);
void bcc_lrup_QC_YnWt_stage_a2(struct bcc_ajdk *D);
void bcc_lrup_QC_YnWt_stage_a3(struct bcc_ajdk *D);
void bcc_lrup_QC_YnWt_stage_b0(struct bcc_ajdk *D);
void bcc_lrup_QC_YnWt_stage_b1(struct bcc_ajdk *D);
void bcc_lrup_QC_YnWt_stage_b2(struct bcc_ajdk *D);
void bcc_lrup_QC_YnWt_stage_b3(struct bcc_ajdk *D);
void bcc_lrup_QC_YnWt_stage_c(struct bcc_ajdk *D);
void bcc_lrup_QR_YnWt_stage_0(struct bcc_ajdk *D);
void bcc_lrup_QR_YnWt_stage_1(struct bcc_ajdk *D);
void bcc_lrup_QR_YnWt_stage_2(struct bcc_ajdk *D);
void bcc_lrup_QR_YnWt_stage_3(struct bcc_ajdk *D);
void bcc_lrup_QC_ZtSWn_stage_0(struct bcc_ajdk *D);
void bcc_lrup_QC_ZtSWn_stage_1(struct bcc_ajdk *D);
void bcc_lrup_QC_ZtSWn_stage_2(struct bcc_ajdk *D);
void bcc_test_mxcut(struct bcc_ajdk *D);
void bcc_lrup_mxcut(struct bcc_ajdk *D);
void bcc_lrup_mxset(struct bcc_ajdk *D);
void bcc_lrup_mxdup(struct bcc_ajdk *D);
void bcc_lrup_test();
void bcc_lrup_flattenloop_test();
void bcc_scorebox_mxA(struct bcc_ajdk *D,int rdrop,int cdrop);
void bcc_scorebox_svalue(struct bcc_ajdk *D);
void bcc_scorebox_srt(struct bcc_ajdk *D,int nrows,int *mr_index_local_nb,int *mr_index_local_mr,int ncols,int *mc_srt);
void bcc_out_xdrop_lkp(struct bcc_ajdk *D,int nrows,int *mr_index_sort,int **mr_index_local_nb_,int **mr_index_local_mr_);
void xcc_out_xdrop_srt(int nlines,int *out_xdrop,int *nrows_,int **mr_index_sort_,int *ncols_,int **mc_srt_);
void xcc_out_xdrop_get(char *fname,int *nlines_,int **out_xdrop_);
void xcc_scorebox_xdrop_array(int nxols,int x_max,int x_num,int **xval_,int **xdrop_);
void xcc_out_trace_get(char *fname,int *nlines_,int **out_trace_ni_,int **out_trace_nr_,int **out_trace_nc_,double **out_trace_QR_,double **out_trace_QC_,int **out_trace_nb_);
void bcc_sumscores_mxB(struct bcc_ajdk *D);
void bcc_sumscores_mxC(struct bcc_ajdk *D);
void bcc_sumscores_dmp(struct bcc_ajdk *D);
void bcc_sumscores_mxA(struct bcc_ajdk *D);
void xcc_get_Ireq(char *strategy,int nbins,int Irem,int *Ireq_p,int *Irow_p,int *Icol_p);
void bcc_sumscores_xij(struct bcc_ajdk *D);
void bcc_sumscores_cmb(struct bcc_ajdk *D);
void bcc_sumscores_srt(struct bcc_ajdk *D);
void bcc_sumscores_ifT(struct bcc_ajdk *D);
void bcc_sumscores_nrm(struct bcc_ajdk *D);
void bcc_sumscores_test();
void bcc_lrup_sumscores_test();
void MDA_io_test();
void MDA_write_i4(int n_d,int *d_,int *i4_,char *fname);
void MDA_read_i4(int *n_d_p,int **d_p,int **i4_p,char *fname);
void MDA_write_r8(int n_d,int *d_,double *r8_,char *fname);
void MDA_read_r8(int *n_d_p,int **d_p,double **r8_p,char *fname);
void MDA_write_ulli(int n_d,int *d_,unsigned long long int *ulli_,char *fname);
void MDA_read_ulli(int *n_d_p,int **d_p,unsigned long long int **ulli_p,char *fname);
void binary_read(char *filename,int *bitj_p,int *nrows_p,int *ncols_p,unsigned char **A_p);
void binary_read_getsize(char *filename,int *bitj_p,int *nrows_p,int *ncols_p);
void binary_write(char *filename,int bitj,int nrows,int ncols,unsigned char *A);
void uchar_write_as_binary(char *filename,int bitj,int nrows,int ncols,unsigned char *u_);
void wrap_rand(char *filename,int bitj,int nrows,int ncols,int type_set);
void wrap_transpose(char *filename,char *outname,int nbytes_per_read);
void wrap_transpose_test();
void *get_D_AtTn_ZtSn_vv(void *vp);
int wrap_D_AtTn_ZtSn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_ZtSn,struct L_handle **output_p);
void *get_D_AtTn_AtTn_vv(void *vp);
int wrap_D_AtTn_AtTn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_ZtSn,struct L_handle **output_p);
void *get_D_AtTn_ZtSn_uu(void *vp);
int wrap_D_AtTn_ZtSn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle **output_p);
void *get_D_AtTn_AtTn_uu(void *vp);
int wrap_D_AtTn_AtTn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle **output_p);
void wrap_D_AtTn_ZtSn_vv_test();
void dcc_ajdk_copy(struct dcc_ajdk *D,struct dcc_ajdk *D_in);
void dcc_ajdk_load(struct dcc_ajdk *D);
void dcc_ajdk_init(double mrnd,struct dcc_ajdk *D);
void dcc_single_copy_M_An(struct dcc_single *E,struct dcc_single *E_in);
void dcc_single_load_M_An(struct dcc_single *E,struct dcc_single *E_0);
void dcc_single_load_M_An_bkp(struct dcc_single *E);
void dcc_single_init_M_An(char *error_vs_speed,double mrnd,struct dcc_single *E);
void dcc_single_copy_lf(struct dcc_single *E,struct dcc_single *E_in);
void dcc_single_init_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,struct M_handle **M_p);
void dcc_single_init_lf(struct dcc_single *E);
void dcc_double_copy_lf(struct dcc_double *F,struct dcc_double *F_in);
void dcc_double_init_lf(struct dcc_double *F);
void dcc_copy_A_p(struct dcc_ajdk *D,struct dcc_ajdk *D_in);
void dcc_X_nrows_total(struct dcc_ajdk *D);
void dcc_load_A_p(struct dcc_ajdk *D);
void dcc_init_A_p(double mrnd,struct dcc_ajdk *D);
void dcc_copy_QX(struct dcc_ajdk *D,struct dcc_ajdk *D_in);
void dcc_init_QX(struct dcc_ajdk *D);
void dcc_M_mxset(struct dcc_ajdk *D);
void dcc_copy(struct dcc_ajdk *D,struct dcc_ajdk *D_in);
void dcc_load(struct dcc_ajdk **D_p,struct dcc_single ***E_p,struct dcc_double ***F_p);
void dcc_init(char *error_vs_speed,double mrnd,struct dcc_ajdk **D_p,struct dcc_single ***E_p,struct dcc_double ***F_p);
void dcc_init_test();
void dcc_An_ajdk(struct dcc_ajdk *D);
void dcc_lf_ZtSn(struct dcc_ajdk *D);
void dcc_lf_D_AtTn_ZtSn_vv(struct dcc_ajdk *D);
void dcc_lf_D_AtTn_ZtSn_uu(struct dcc_ajdk *D);
void dcc_lf_D_AtTn_ZtSn_error(int verbose,struct dcc_ajdk *D);
void dcc_lf_D_AtTn_ZtSn_test();
void dcc_lf_TAnZtS_ww(struct dcc_ajdk *D);
void dcc_lf_TAnZtS_uu(struct dcc_ajdk *D);
void dcc_lf_TAnZtS_error(int verbose,struct dcc_ajdk *D);
void dcc_lf_TAnZtS_test();
void *get_dcc_halfloop(void *vp);
void wrap_dcc_halfloop(int *tidx,void **vpra,pthread_t *thread_in,struct dcc_double *F);
void dcc_wrap_dcc_halfloop(struct dcc_ajdk *D);
void dcc_lrup_mxdup(struct dcc_ajdk *D);
void dcc_scorebox_mxA(struct dcc_ajdk *D,int rdrop,int cdrop);
void dcc_scorebox_svalue(struct dcc_ajdk *D);
void dcc_scorebox_srt(struct dcc_ajdk *D,int nrows,int *mr_index_local_nb,int *mr_index_local_mr,int ncols,int *mc_srt);
void dcc_out_xdrop_lkp(struct dcc_ajdk *D,int nrows,int *mr_index_sort,int **mr_index_local_nb_,int **mr_index_local_mr_);
void dcc_sumscores_mxB(struct dcc_ajdk *D);
void dcc_sumscores_mxC(struct dcc_ajdk *D);
void dcc_sumscores_dmp(struct dcc_ajdk *D);
void dcc_sumscores_mxA(struct dcc_ajdk *D);
void dcc_sumscores_xij(struct dcc_ajdk *D);
void dcc_sumscores_cmb(struct dcc_ajdk *D);
void dcc_sumscores_srt(struct dcc_ajdk *D);
void dcc_sumscores_ifT(struct dcc_ajdk *D);
void dcc_sumscores_nrm(struct dcc_ajdk *D);
void dcc_sumscores_test();
void dcc_time_sumscores_test();
void dexcluster_driver();
void dexcluster_scorebox_excerpt_0(int verbose,struct dcc_ajdk *D,struct S_handle *S);
void dexcluster_scorebox_excerpt_1(int verbose,struct dcc_ajdk *D,struct S_handle *S);
void dexcluster_scorebox_excerpt_2(int verbose,struct dcc_ajdk *D,struct S_handle *S);
void dexcluster_scorebox_rc();
void dexcluster_scorebox_cr();
void dexcluster_scorebox_xx();
void dexcluster_scorebox();
void GLOBAL_ops_reset_one(int nt);
void GLOBAL_ops_count_one(int nt,unsigned long long int f,unsigned long long int b);
void GLOBAL_ops_addup_all();
void GLOBAL_ops_addup_one(int nt);
void GLOBAL_ops_printf_all(int verbose,char *prefix);
void GLOBAL_tic(int nx);
void GLOBAL_toc(int nx,int verbose,char *prefix);
void GLOBAL_ops_toc(int nt,int nx,int verbose,char *prefix);
void GLOBAL_pthread_tic();
void GLOBAL_pthread_toc();
void GLOBAL_pthread_tuc();
void timing_dump(int nt);
void lakcluster_driver();
void lakcluster_scorebox_excerpt_0(int verbose,struct bcc_ajdk *D,struct S_handle *S);
void lakcluster_scorebox_excerpt_1(int verbose,struct bcc_ajdk *D,struct S_handle *S);
void lakcluster_scorebox_excerpt_2(int verbose,struct bcc_ajdk *D,struct S_handle *S);
void lakcluster_scorebox_rc();
void lakcluster_scorebox_cr();
void lakcluster_scorebox_xx();
void lakcluster_scorebox();
struct L_handle *L_handle_make(unsigned long long int length);
void L_handle_copy(struct L_handle *L,struct L_handle *L_in);
void L_zero(struct L_handle *L);
void L2_transpose(struct L_handle *L1,struct L_handle *L2);
void L2_duplicate(struct L_handle *L1,struct L_handle *L2);
void L3_duplicate(struct L_handle *L1,struct L_handle *L2);
double *L3_get(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a);
double *L2_get(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a);
double *L1_get(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a);
void L3_set(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a,double val);
void L2_set(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,double val);
void L1_set(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,double val);
void L3_plusequals(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a,double val);
void L2_plusequals(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,double val);
void L1_plusequals(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,double val);
double *L3_lf_get(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a);
double *L2_lf_get(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a);
double *L1_lf_get(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a);
void L3_lf_set(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a,double val);
void L2_lf_set(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,double val);
void L1_lf_set(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,double val);
void L3_lf_plusequals(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a,double val);
void L2_lf_plusequals(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,double val);
void L1_lf_plusequals(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,double val);
void L2_clean(struct L_handle *L,unsigned long long int nrows,unsigned char *mr_b,unsigned char *mr_j,unsigned long long int ncols,unsigned char *mc_b,unsigned char *mc_j);
void LL2_plustimesequals(struct L_handle *LA,struct L_handle *LB,double dtmp,unsigned long long int nrows,unsigned char *mr_b,unsigned char *mr_j,unsigned long long int ncols,unsigned char *mc_b,unsigned char *mc_j);
void lfprintf(struct L_handle *L,char *prefix);
void *get_M_An_to_L2_run(void *vp);
int wrap_M_An_to_L2_run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,int output_spacing_c,struct M_handle *M_An,double *A_ajdk,struct L_handle **output_L_p);
void *get_M_At_to_L2_run(void *vp);
int wrap_M_At_to_L2_run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,int output_spacing_c,struct M_handle *M_At,double *A_ajdk,struct L_handle **output_L_p);
void wrap_M_Ax_to_L2_test();
void mda_write_d2_r8(char *fname,int nrows,int ncols,double *ra);
void mda_write_d3_r8(char *fname,int nrows,int ncols,int nlyrs,double *ra);
void mda_write_d2_i4(char *fname,int nrows,int ncols,int *ra);
void mda_read_d1_r8(char *fname,int *nrows_p,double *dra);
void mda_read_d2_r8(char *fname,int *nrows_p,int *ncols_p,double *dra);
void mda_read_d3_r8(char *fname,int *nrows_p,int *ncols_p,int *nlyrs_p,double *dra);
unsigned long long int M_wkspace_copy(struct M_handle *M,struct M_handle *M_in);
unsigned long long int M_wkspace_alloc(struct M_handle *M,int estim_flag);
void M_handle_printf(struct M_handle *M,int verbose,char *prefix);
struct M_handle * M_handle_make(int bitj,int nrows,int ncols,char *A_filename,unsigned char *Ara);
void M_handle_copy(struct M_handle *M,struct M_handle *M_in);
void M_fp_free(struct M_handle *M);
void M_load(struct M_handle *M);
void M_mxget_excerpt(int verbose,int nrows,unsigned char *mr_b,unsigned char *mr_j,long long int *rpop_b_p,long long int *rpop_j_p,unsigned int *m_a_,unsigned int *m_b_);
void M_mxget(struct M_handle *M);
void M_mxget_bkp(struct M_handle *M);
void M_mxset(struct M_handle *M_An,unsigned char *mr_j,unsigned char *mc_j);
struct M_handle * M_handle_v_make(int bitj,int nrows,int ncols,char *A_filename,unsigned char *Ara,unsigned char *mr_b,unsigned char *mc_b);
struct M_handle * M_handle_w_make(int bitj,int b_mlt,int nrows_bZ,int ncols_bZ_tmp);
void hsv2rgb(double h,double s,double v,double *r,double *g,double *b);
void colorscale(double val,double valmin,double valmax,double *rcolor,double *gcolor,double *bcolor);
int WritePNMfile_color(double *ra,int rows,int cols,double min,double max,char *filename);
#ifndef _MONOLITH
unsigned int popcount_uchar_array(unsigned char *wp,unsigned long wl);
long long int popcount(__m128i** mem1p, __m128i** maskp, __m128i** maskp_end);
long long int popcount_and(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end);
long long int popcount_xor(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end);
long long int popcount_notxorxor(__m128i** mem1p, __m128i** memSp,__m128i** mem2p, __m128i** maskp, __m128i** maskp_end);
double popcount_lf(__m128i** mem1p, __m128i** maskp, __m128i** maskp_end, double **dinp);
double popcount_xor_lf(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end, double **dinp);
double popcount_and_lf(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end, double **dinp);
double popcount_pm0(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end);
long long int popcount_pmpm0(__m128i** mem0p,__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end);
double popcount_pm0_lf(__m128i** mem1p, __m128i** mem2p, __m128i** maskp, __m128i** maskp_end,double **dinp);
#endif /* _MONOLITH */
int dPartition(double *dra,int stride,int l,int r);
void dQuickSort(unsigned int nn,double *dra,int stride,int l,int r);
int dPartition_xij(double *dra,int stride,int *ira,int *ira_b,int *ira_j,int *lnb,int *mr_a,int *mr_b,int l,int r);
unsigned int dQuickSort_xij(unsigned int nn,double *dra,int stride,int *ira,int *ira_b,int *ira_j,int *lnb,int *mr_a,int *mr_b,int l,int r);
int ulliPartition_index(unsigned long long int *ulli_,int stride,int *i_,int l,int r) ;
unsigned int ulliQuickSort_index(unsigned int nn,unsigned long long int *ulli_,int stride,int *i_,int l,int r);
void ulliQuickSort_index_driver(int n_ulli,unsigned long long int *ulli_,int stride,unsigned long long int *ulli_workspace_,int *i_);
int dPartition_index(double *d_,int stride,int *i_,int l,int r) ;
unsigned int dQuickSort_index(unsigned int nn,double *d_,int stride,int *i_,int l,int r);
void dQuickSort_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *i_);
void dQuickSort_index_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *i_orig_from_sort_,int *i_sort_from_orig_,int *i_workspace_);
int lPartition_index(int *l_,int stride,int *i_,int l,int r);
unsigned int lQuickSort_index(unsigned int nn,int *l_,int stride,int *i_,int l,int r);
void lQuickSort_index_driver(int n_l,int *l_,int stride,int *l_workspace_,int *i_);
int charpPartition_index(char **charp_,int stride,int *i_,int l,int r);
unsigned int charpQuickSort_index(unsigned int nn,char **charp_,int stride,int *i_,int l,int r);
void charpQuickSort_index_driver(int n_charp,char **charp_,int stride,char **charp_workspace_,int *i_);
void charpQuickSort_index_test();
void charpQuickSort_index_index_driver(int n_d,char **charp_,int stride,char **charp_workspace_,int *i_orig_from_sort_,int *i_sort_from_orig_,int *i_workspace_);
void raprintf(void *vra,char *type,int rows,int cols,char *prefix);
void ra_fprintf(char *fname,void *vra,char *type,int rows,int cols,char *prefix);
void getBinW(unsigned char *w, char *str, int k);
void bprintf(unsigned char *w,int bitj,int nrows,int ncols,char *prefix);
void ra_stats(void *vra,char *type,unsigned long long int length,void *max_p,void *min_p,double *mean_p,double *stdev_p);
void dra_plustimesequals_s___m_m(double *dra1,unsigned long long int nrows,unsigned long long int ncols,double *dra2,double multiplier,unsigned char *bmr_b,unsigned char *bmc_b);
void dra_plusdivequals_s___m_m(double *dra1,unsigned long long int nrows,unsigned long long int ncols,double *dra2,double denominator,unsigned char *bmr_b,unsigned char *bmc_b);
void dra_mds_nrm_s___m_m(double *dra1,unsigned long long int nrows,unsigned long long int ncols,double denominator,unsigned char *bmr_b,unsigned char *bmc_b);
void dra_mds_pow_s___m_m(double *dra1,unsigned long long int nrows,unsigned long long int ncols,double denominator,unsigned char *bmr_b,unsigned char *bmc_b);
void dra_sumx(unsigned long long int nx,double *dra,unsigned long long int stride,unsigned long long int length);
void dra_plustimesequals(double *dra1,unsigned long long int length,double *dra2,double multiplier);
void dra_plusequals(double *dra1,unsigned long long int length,double *dra2);
void dra_subtequals(double *dra1,unsigned long long int length,double *dra2);
void dra_plus(double *dra1,unsigned long long int length,double adder);
void dra_times(double *dra,unsigned long long int length,double multiplier);
void dra_dup(double *dra1,unsigned long long int length,double *dra2);
void ura_dup(unsigned char *ura1,unsigned long long int length,unsigned char *ura2);
double dra_diff(double *dra,double *drb,unsigned long long int length,unsigned long long int stride);
struct S_handle *S_handle_make(char *fname,int row_max,int row_num,int col_max,int col_num);
void S_init_bcc(struct bcc_ajdk *D,struct S_handle *S);
void S_init_dcc(struct dcc_ajdk *D,struct S_handle *S);
void S_handle_printf(int verbose,struct S_handle *S,char *prefix);
void S_handle_dmp(struct S_handle *S);
void *get_dcc_TAnZtS_ww_stage_0(void *vp);
void wrap_dcc_TAnZtS_ww_stage_0(int *tidx,void **vpra,pthread_t *thread_in,struct dcc_single *E);
void dcc_TAnZtS_ww_stage_0(struct dcc_ajdk *D);
void dcc_TAnZtS_ww_stage_1(struct dcc_ajdk *D);
void *get_dcc_TAnZtS_ww_stage_2(void *vp);
void wrap_dcc_TAnZtS_ww_stage_2(int *tidx,void **vpra,pthread_t *thread_in,struct dcc_double *F);
void dcc_TAnZtS_ww_stage_2(struct dcc_ajdk *D);
void *get_TAnZtS_uu(void *vp);
int wrap_TAnZtS_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_Tt,struct M_handle *M_Zn,struct M_handle *M_St,double *A_ajdk,struct L_handle **output_p);
void *get_TAnAtT_uu(void *vp);
int wrap_TAnAtT_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_Tt,struct M_handle *M_Zn,struct M_handle *M_St,double *A_ajdk,struct L_handle **output_p);
void updateglobals(char *vname);
void read_input();
unsigned char *wkspace_alloc_nocheck(unsigned long long int size);
struct wkspace_point * wkspace_make_point();
struct wkspace_point * wkspace_set_point(struct wkspace_point *w);
void wkspace_printf_point(struct wkspace_point *w);
long long int wkspace_check_point(struct wkspace_point *w);
unsigned char* wkspace_alloc(unsigned long long int size);
unsigned char* wkspace_all0c(unsigned long long int size);
void wkspace_reset(void* new_base);
void wkspace_printf();
void xcalc_setmaxmin(struct L_handle *lf_YnWt,double *YnWt,int rpop_j,int rpop_b,int nrows,int cpop_j,int cpop_b,int ncols,unsigned char *mr_j,unsigned char *mr_b,unsigned char *mc_j,unsigned char *mc_b,int trm_flag,double *max_data_p,double *min_data_p);
void *get_xcalc(void *vp);
void wrap_xcalc(int *tidx,void **vpra,pthread_t *thread_in,unsigned char *mr_j,unsigned char *mr_b,unsigned char *mc_j,unsigned char *mc_b,struct L_handle *lf_YnWt,double *YnWt,int *Yn_nrows_p,int *Wn_nrows_p,struct M_handle **M_X_p,int *b_mlt_p,int trm_flag);
void *get_At_T_YnWt_ww(void *vp);
int wrap_At_T_YnWt_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_t,int spacing_w,struct M_handle *M_At,struct M_handle *M_Tt,struct L_handle *lf_YnWt,struct M_handle *M_YnWt,struct M_handle *M_Wn,double *A_ajdk,int trm_flag,struct L_handle **output_At_T_YnWt_uu_p,struct L_handle **output_At_T_YnWt_ww_p);
void wrap_At_T_YnWt_ww_test();
void *get_An_ZtSWn_ww(void *vp);
int wrap_An_ZtSWn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,int spacing_w,struct M_handle *M_An,struct M_handle *M_St,struct L_handle *lf_ZtSWn,struct M_handle **M_ZtSWn_,struct M_handle *M_Wt,double *A_ajdk,int trm_flag,struct L_handle **output_An_ZtSWn_uu_p,struct L_handle **output_An_ZtSWn_ww_p);
void wrap_An_ZtSWn_ww_test();
void get_xdrop(double lrij,double lcij,int *rdrop_p,int *cdrop_p);
int get_xdrop_total(double lrij,double lcij);
void get_xdrop_array(double lrij,double lcij,int *length_p,int **rdrop_p_,int **cdrop_p_,int **rkeep_p_,int **ckeep_p_);
void get_xdrop_array_sub(double lrij,double lcij,int iteration_num,int iteration_max,int iteration_min,int *length_p,int **rdrop_p_,int **cdrop_p_,int **rkeep_p_,int **ckeep_p_);
void *get_At_T_Xn_ww(void *vp);
int wrap_At_T_Xn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_t,int spacing_w,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Xn,struct M_handle *M_Xt,double *A_ajdk,int trm_flag,struct L_handle **output_At_T_Xn_ww_p);
void *get_At_T_Xn_uu(void *vp);
int wrap_At_T_Xn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_t,int spacing_w,struct M_handle *M_At,struct M_handle *M_Tt,struct L_handle *lf_Xn,struct M_handle *M_Xt,double *A_ajdk,int trm_flag,struct L_handle **output_At_T_Xn_uu_p);
void wrap_At_T_Xn_ww_test();
void *get_An_Xn_ww(void *vp);
int wrap_An_Xn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,int spacing_w,struct M_handle *M_An,/* struct M_handle *M_St, */struct M_handle *M_Xn,struct M_handle *M_Xt,double *A_ajdk,int trm_flag,struct L_handle **output_An_Xn_ww_p);
void *get_An_Xn_uu(void *vp);
int wrap_An_Xn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,int spacing_w,struct M_handle *M_An,/* struct M_handle *M_St, */struct L_handle *lf_Xn,struct M_handle *M_Xt,double *A_ajdk,int trm_flag,struct L_handle **output_An_Xn_uu_p);
void wrap_An_Xn_ww_test();
struct P_handle *P_handle_make(char *pca_infix,char *out_xdrop_fname,char *V_fname,struct dcc_ajdk *D,int iteration_num,int iteration_max,int iteration_min,int rank,double tolerance,int b_mlt);
void P_init_dcc(struct dcc_ajdk *D,struct P_handle *P);
void P_handle_printf(int verbose,struct P_handle *P,char *prefix);
void P_handle_dmp(struct P_handle *P);
void pca_driver();
struct R_handle *R_handle_make(char *out_xdrop_fname,unsigned long int rseed);
void R_handle_printf(int verbose,struct R_handle *R,char *prefix);
void dcc_scramble(struct dcc_ajdk *D,struct R_handle *R);
void bcc_scramble(struct bcc_ajdk *D,struct R_handle *R);
void wrap_dcc_scramble(struct dcc_ajdk *D);
void wrap_bcc_scramble(struct bcc_ajdk *D);
/* %%%%%%%%%%%%%%%% */
void dcg_ajdk_copy(struct dcg_ajdk *D,struct dcg_ajdk *D_in);
void dcg_ajdk_load(struct dcg_ajdk *D);
void dcg_ajdk_init(double mrnd,struct dcg_ajdk *D);
void dcg_single_copy_M_An(struct dcg_single *E,struct dcg_single *E_in);
void dcg_single_load_M_An(struct dcg_single *E);
void dcg_single_init_M_An(char *error_vs_speed,double mrnd,struct dcg_single *E);
void dcg_single_copy_lf(struct dcg_single *E,struct dcg_single *E_in);
void dcg_single_init_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,struct M_handle **M_p);
void dcg_single_init_lf(struct dcg_single *E);
void dcg_double_copy_lf(struct dcg_double *F,struct dcg_double *F_in);
void dcg_double_init_lf(struct dcg_double *F);
void dcg_copy_A_p(struct dcg_ajdk *D,struct dcg_ajdk *D_in);
void dcg_X_nrows_total(struct dcg_ajdk *D);
void dcg_load_A_p(struct dcg_ajdk *D);
void dcg_init_A_p(double mrnd,struct dcg_ajdk *D);
void dcg_copy_QX(struct dcg_ajdk *D,struct dcg_ajdk *D_in);
void dcg_init_QX(struct dcg_ajdk *D);
void dcg_M_mxset(struct dcg_ajdk *D);
void dcg_copy(struct dcg_ajdk *D,struct dcg_ajdk *D_in);
void dcg_load(struct dcg_ajdk **D_p,struct dcg_single ***E_p,struct dcg_double ***F_p);
void dcg_init(char *error_vs_speed,double mrnd,struct dcg_ajdk **D_p,struct dcg_single ***E_p,struct dcg_double ***F_p);
void dcg_init_test();
void dcg_An_ajdk(struct dcg_ajdk *D);
void dcg_lf_ZtSn(struct dcg_ajdk *D);
void dcg_lf_D_AtTn_ZtSn_vv(struct dcg_ajdk *D);
void dcg_lf_D_AtTn_ZtSn_uu(struct dcg_ajdk *D);
void dcg_lf_D_AtTn_ZtSn_error(int verbose,struct dcg_ajdk *D);
void dcg_lf_D_AtTn_ZtSn_test();
void dcg_lf_TAnZtS_ww(struct dcg_ajdk *D);
void dcg_lf_TAnZtS_uu(struct dcg_ajdk *D);
void dcg_lf_TAnZtS_error(int verbose,struct dcg_ajdk *D);
void dcg_lf_TAnZtS_test();
void *get_dcg_halfloop(void *vp);
void wrap_dcg_halfloop(int *tidx,void **vpra,pthread_t *thread_in,struct dcg_double *F);
void dcg_wrap_dcg_halfloop(struct dcg_ajdk *D);
void dcg_lrup_mxdup(struct dcg_ajdk *D);
void dcg_scorebox_mxA(struct dcg_ajdk *D,int rdrop,int cdrop);
void dcg_scorebox_svalue(struct dcg_ajdk *D);
void dcg_scorebox_srt(struct dcg_ajdk *D,int nrows,int *mr_index_local_nb,int *mr_index_local_mr,int ncols,int *mc_srt);
void dcg_out_xdrop_lkp(struct dcg_ajdk *D,int nrows,int *mr_index_sort,int **mr_index_local_nb_,int **mr_index_local_mr_);
void dcg_sumscores_mxB(struct dcg_ajdk *D);
void dcg_sumscores_mxC(struct dcg_ajdk *D);
void dcg_sumscores_dmp(struct dcg_ajdk *D);
void dcg_sumscores_mxA(struct dcg_ajdk *D);
void dcg_sumscores_xij(struct dcg_ajdk *D);
void dcg_sumscores_cmb(struct dcg_ajdk *D);
void dcg_sumscores_srt(struct dcg_ajdk *D);
void dcg_sumscores_ifT(struct dcg_ajdk *D);
void dcg_sumscores_nrm(struct dcg_ajdk *D);
void dcg_sumscores_test();
void dcg_time_sumscores_test();
void dcgxpander_driver();
void dcg_TAnZtS_ww_stage_0(struct dcg_ajdk *D);
void dcg_TAnZtS_ww_stage_1(struct dcg_ajdk *D);
void wrap_dcg_TAnZtS_ww_stage_2(int *tidx,void **vpra,pthread_t *thread_in,struct dcg_double *F);
void dcg_TAnZtS_ww_stage_2(struct dcg_ajdk *D);
void dcg_scramble(struct dcg_ajdk *D,struct R_handle *R);
void wrap_dcg_scramble(struct dcg_ajdk *D);
/* %%%%%%%%%%%%%%%% */
unsigned long long int fsize_0(char *);
unsigned long int wc_0(char *);
void bed_to_b16_test();
void bed_to_b16();
