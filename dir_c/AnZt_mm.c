#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */


void *get_AnZt_mm__run(void *vp)
{
  /* This function serves as a module which calculates the matrix-matrix product AnZt, where both An and Zn are passed in as L_handles. ;
     The calculation is performed using matrix-operations. ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct L_handle *lf_At = (struct L_handle *)(vpra[ip++]);
  struct L_handle *lf_Zt = (struct L_handle *)(vpra[ip++]);
  struct L_handle *output_AnZt = (struct L_handle *)(vpra[ip++]);
  int aper = maximum(1,minimum(lf_At->col_stride,maximum(1,/*262144*/ GLOBAL_2_cache_size / maximum(1,lf_At->row_stride*sizeof(double)) * 3/4 )));
  int zper = maximum(1,minimum(lf_Zt->col_stride,maximum(1,/*262144*/ GLOBAL_2_cache_size / maximum(1,lf_Zt->row_stride*sizeof(double)) * 3/4 )));
  int ma_bgn=0,ma_end=0,ma_tot=0,ma_block=0,ma_block_max=0;
  int mz_bgn=0,mz_end=0,mz_tot=0,mz_block=0,mz_block_max=0;
  int mx=0;
  ma_block_max = rup(lf_At->col_stride,aper)/aper;
  mz_block_max = rup(lf_Zt->col_stride,zper)/zper;
  if (verbose>1){ printf(" %% [entering get_AnZt_mm__run] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% lf_At (%d,%d) aper = (l2 %d)/%d*3/4 = %d, ma_block_max %d\n",lf_At->row_stride,lf_At->col_stride,(int)GLOBAL_2_cache_size,(int)lf_At->col_stride,aper,ma_block_max);}
  if (verbose>1){ printf(" %% lf_Zt (%d,%d) zper = (l2 %d)/%d*3/4 = %d, mz_block_max %d\n",lf_Zt->row_stride,lf_Zt->col_stride,(int)GLOBAL_2_cache_size,(int)lf_Zt->col_stride,zper,mz_block_max);}
  if (verbose>1){ printf(" %% calculating output_AnZt\n");}
  output_AnZt->spacing_row = lf_At->spacing_col; output_AnZt->row_stride = lf_At->col_stride;
  output_AnZt->spacing_col = lf_Zt->spacing_col; output_AnZt->col_stride = lf_Zt->col_stride;
  fill_uchar_zero((unsigned char *)(output_AnZt->lf),output_AnZt->row_stride*output_AnZt->col_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"AnZt_mm")){ goto skip_AnZt_mm;}
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zt->row_stride*lf_Zt->col_stride>0){
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,lf_At->col_stride,lf_Zt->col_stride,lf_At->row_stride,1.0,lf_At->lf,lf_At->row_stride,lf_Zt->lf,lf_Zt->row_stride,0.0,output_AnZt->lf,output_AnZt->row_stride);
    /* if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zt->row_stride*lf_Zt->col_stride>0){ } */}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)lf_At->col_stride*(unsigned long long int)lf_Zt->col_stride*(unsigned long long int)lf_At->row_stride,0);
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp__on){
#pragma omp parallel private(mx,ma_bgn,ma_end,ma_tot,ma_block)
    { /* begin omp parallel */
      mx=0;
      ma_bgn=0,ma_end=0,ma_tot=0,ma_block=0;
#pragma omp for schedule(dynamic)
      for (mx=0;mx<ma_block_max;mx++){
	ma_block = mx;
	ma_bgn = minimum(lf_At->col_stride,(ma_block+0)*aper);
	ma_end = minimum(lf_At->col_stride,(ma_block+1)*aper);
	ma_tot = ma_end-ma_bgn;
	if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zt->row_stride*lf_Zt->col_stride>0){
	  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,ma_tot,lf_Zt->col_stride,lf_At->row_stride,1.0,&(lf_At->lf[ma_bgn*lf_At->row_stride]),lf_At->row_stride,lf_Zt->lf,lf_Zt->row_stride,0.0,&(output_AnZt->lf[ma_bgn]),output_AnZt->row_stride);
	  /* if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zt->row_stride*lf_Zt->col_stride>0){ } */}
	/* for (mx=0;mx<ma_block_max;mx++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)lf_At->col_stride*(unsigned long long int)lf_Zt->col_stride*(unsigned long long int)lf_At->row_stride,0);
    /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp_unused){
#pragma omp parallel private(mx,ma_bgn,ma_end,ma_tot,ma_block,mz_bgn,mz_end,mz_tot,mz_block)
    { /* begin omp parallel */
      mx=0;
      ma_bgn=0,ma_end=0,ma_tot=0,ma_block=0;
      mz_bgn=0,mz_end=0,mz_tot=0,mz_block=0;
#pragma omp for schedule(dynamic)
      for (mx=0;mx<ma_block_max*mz_block_max;mx++){
	ma_block = mx%ma_block_max;
	mz_block = (mx-ma_block)/ma_block_max;
	ma_bgn = minimum(lf_At->col_stride,(ma_block+0)*aper);
	ma_end = minimum(lf_At->col_stride,(ma_block+1)*aper);
	ma_tot = ma_end-ma_bgn;
	mz_bgn = minimum(lf_Zt->col_stride,(mz_block+0)*zper);
	mz_end = minimum(lf_Zt->col_stride,(mz_block+1)*zper);
	mz_tot = mz_end-mz_bgn;
	if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zt->row_stride*lf_Zt->col_stride>0){
	  cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,ma_tot,mz_tot,lf_At->row_stride,1.0,&(lf_At->lf[ma_bgn*lf_At->row_stride]),lf_At->row_stride,&(lf_Zt->lf[mz_bgn*lf_Zt->row_stride]),lf_Zt->row_stride,0.0,&(output_AnZt->lf[ma_bgn + mz_bgn*lf_At->col_stride]),output_AnZt->row_stride);
	  /* if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zt->row_stride*lf_Zt->col_stride>0){ } */}
	/* for (mx=0;mx<ma_block_max*mz_block_max;mx++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)lf_At->col_stride*(unsigned long long int)lf_Zt->col_stride*(unsigned long long int)lf_At->row_stride,0);
    /* else if (GLOBAL_omp_type==GLOBAL_omp_unused){ } */}
 skip_AnZt_mm:
  if (verbose>2){ raprintf(output_AnZt->lf,"double",output_AnZt->row_stride,output_AnZt->col_stride," %% output_AnZt: ");}
  if (verbose>1){ printf(" %% [finished get_AnZt_mm__run] tidx %d \n",tidx);}
  return NULL;
}

int wrap_AnZt_mm__run(int *tidx,void **vpra,pthread_t *thread_in,struct L_handle *lf_At,struct L_handle *lf_Zt,struct L_handle **output_AnZt_p)
{
  /* This function uses the L_handles lf_At and lf_Zt to run a series of parallel calls to get_AnZt_mm__run ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 3)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_r=0,length_c=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_AnZt_mm__run] tidx %d\n",*tidx);}
  if (verbose){ raprintf(lf_At->lf,"double",lf_At->row_stride,lf_At->col_stride," %% lf_At: ");}
  if (verbose){ raprintf(lf_Zt->lf,"double",lf_Zt->row_stride,lf_Zt->col_stride," %% lf_Zt: ");}
  length_r = lf_At->col_stride; length_c = lf_Zt->col_stride;
  length = length_r*length_c; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_r,length_c,length);}
  if (*output_AnZt_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_AnZt_p = L_handle_make(length);}
  length = length_r*length_c; if ((*output_AnZt_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_AnZt_mm__run\n",(*output_AnZt_p)->length,length);} memset((*output_AnZt_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = lf_At; vpra[ip++] = lf_Zt; vpra[ip++] = *output_AnZt_p;
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AnZt_mm__run,vpra)){ printf("Warning! cannot create thread %d in wrap_AnZt_mm__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AnZt_mm__run(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_AnZt_p)->lf,"double",length_r,length_c," %% (*output_AnZt_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_AnZt_mm__run] tidx %d\n",*tidx);}
  return length;
}  

void wrap_AnZt_mm_test()
{
  /* test for errors with input file: AnZt_mm_error.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= AnZt_mm;
     GLOBAL_TEST_TYP2= error;
     GLOBAL_NBINS= 4;
     GLOBAL_TEST_mrand= 0.5;
     GLOBAL_TEST_A_n_rows= 19;
     GLOBAL_TEST_A_n_cols= 2002;
     GLOBAL_TEST_Z_n_rows= 21;
     GLOBAL_TEST_Y_n_cols= 1950;
     GLOBAL_TEST_T_n_cols= 1;
     GLOBAL_TEST_niter= 1;
     END= 0;
  */
  /* test for speed with input file: AnZt_mm_speed.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= AnZt_mm;
     GLOBAL_TEST_TYP2= speed;
     GLOBAL_NBINS= 4;
     GLOBAL_TEST_mrand= 0.05;
     GLOBAL_TEST_A_n_rows= 1920;
     GLOBAL_TEST_A_n_cols= 19200;
     GLOBAL_TEST_Z_n_rows= 1920;
     GLOBAL_TEST_Y_n_cols= 19200;
     GLOBAL_TEST_T_n_cols= 1;
     GLOBAL_TEST_niter= 1;
     END= 0;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_At_ajdk=NULL,**lf_Zt_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_At=NULL,**lf_Zt=NULL,**lf_Yn=NULL,**lf_Wn=NULL;
  struct L_handle **lf_AnZt_mm=NULL,**lf_YnWt_mm=NULL; int *length_AnZt_mm=NULL,*length_YnWt_mm=NULL;
  struct L_handle **lf_AnZt_vv=NULL,**lf_YnWt_vv=NULL; int *length_AnZt_vv=NULL,*length_YnWt_vv=NULL;
  struct L_handle **lf_AnZt_uu=NULL,**lf_YnWt_uu=NULL; int *length_AnZt_uu=NULL,*length_YnWt_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_AnZt_vv_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_At_ajdk,&lf_Zt_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  lf_At = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_Zt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_Yn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_Wn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AnZt_mm = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnZt_mm = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_YnWt_mm = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_YnWt_mm = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_At[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_An[nb]->ncols);
    lf_Zt[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->nrows*(unsigned long long int)M_An[nb]->ncols);
    lf_Yn[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Yn[nb]->ncols);
    lf_Wn[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->nrows*(unsigned long long int)M_Yn[nb]->ncols);
    lf_AnZt_mm[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    lf_YnWt_mm[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_AnZt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnZt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_YnWt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_YnWt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_AnZt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    lf_YnWt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
  	if (verbose){ printf(" %% %s; B: %s; A: %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
	  wrap_M_At_to_L2_run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_B,M_At[nb],A_ajdk,&(lf_At[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  wrap_M_At_to_L2_run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_B,M_Zt[nb],A_ajdk,&(lf_Zt[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  wrap_M_At_to_L2_run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_B,M_Yt[nb],Y_ajdk,&(lf_Yn[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  wrap_M_At_to_L2_run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_B,M_Wt[nb],Y_ajdk,&(lf_Wn[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% M_Ax_to_L2: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (verbose && error_check){
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% M_An[%d] (%d,%d,%d)-x-(%d,%d,%d) ; M_At[%d] (%d,%d,%d)-x-(%d,%d,%d) ; lf_At[%d] (%d,%d)\n",nb,M_An[nb]->rpop_j,M_An[nb]->rpop_b,M_An[nb]->nrows,M_An[nb]->cpop_j,M_An[nb]->cpop_b,M_An[nb]->ncols,nb,M_At[nb]->rpop_j,M_At[nb]->rpop_b,M_At[nb]->nrows,M_At[nb]->cpop_j,M_At[nb]->cpop_b,M_At[nb]->ncols,nb,lf_At[nb]->row_stride,lf_At[nb]->col_stride);
	    printf(" %% M_Zn[%d] (%d,%d,%d)-x-(%d,%d,%d) ; M_Zt[%d] (%d,%d,%d)-x-(%d,%d,%d) ; lf_Zt[%d] (%d,%d)\n",nb,M_Zn[nb]->rpop_j,M_Zn[nb]->rpop_b,M_Zn[nb]->nrows,M_Zn[nb]->cpop_j,M_Zn[nb]->cpop_b,M_Zn[nb]->ncols,nb,M_Zt[nb]->rpop_j,M_Zt[nb]->rpop_b,M_Zt[nb]->nrows,M_Zt[nb]->cpop_j,M_Zt[nb]->cpop_b,M_Zt[nb]->ncols,nb,lf_Zt[nb]->row_stride,lf_Zt[nb]->col_stride);
	    printf(" %% M_Yn[%d] (%d,%d,%d)-x-(%d,%d,%d) ; M_Yt[%d] (%d,%d,%d)-x-(%d,%d,%d) ; lf_Yn[%d] (%d,%d)\n",nb,M_Yn[nb]->rpop_j,M_Yn[nb]->rpop_b,M_Yn[nb]->nrows,M_Yn[nb]->cpop_j,M_Yn[nb]->cpop_b,M_Yn[nb]->ncols,nb,M_Yt[nb]->rpop_j,M_Yt[nb]->rpop_b,M_Yt[nb]->nrows,M_Yt[nb]->cpop_j,M_Yt[nb]->cpop_b,M_Yt[nb]->ncols,nb,lf_Yn[nb]->row_stride,lf_Yn[nb]->col_stride);
	    printf(" %% M_Wn[%d] (%d,%d,%d)-x-(%d,%d,%d) ; M_Wt[%d] (%d,%d,%d)-x-(%d,%d,%d) ; lf_Wn[%d] (%d,%d)\n",nb,M_Wn[nb]->rpop_j,M_Wn[nb]->rpop_b,M_Wn[nb]->nrows,M_Wn[nb]->cpop_j,M_Wn[nb]->cpop_b,M_Wn[nb]->ncols,nb,M_Wt[nb]->rpop_j,M_Wt[nb]->rpop_b,M_Wt[nb]->nrows,M_Wt[nb]->cpop_j,M_Wt[nb]->cpop_b,M_Wt[nb]->ncols,nb,lf_Wn[nb]->row_stride,lf_Wn[nb]->col_stride);
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (verbose && error_check){ } */}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
  	  length_AnZt_mm[nb] = wrap_AnZt_mm__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),lf_At[nb],lf_Zt[nb],&(lf_AnZt_mm[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
  	  length_YnWt_mm[nb] = wrap_AnZt_mm__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),lf_Yn[nb],lf_Wn[nb],&(lf_YnWt_mm[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){ }} */}}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZt_mm: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic();
	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_An[nb],A_ajdk,&(lf_At_ajdk[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Zn[nb],A_ajdk,&(lf_Zt_ajdk[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Yn[nb],Y_ajdk,&(lf_Yn_ajdk[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Wn[nb],Y_ajdk,&(lf_Wn_ajdk[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc();
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic(); 
	    length_AnZt_vv[nb] = wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,n_spacing_B,M_An[nb],M_Zn[nb],A_ajdk,lf_At_ajdk[nb],lf_Zt_ajdk[nb],&(lf_AnZt_vv[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic(); 
	    length_YnWt_vv[nb] = wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,n_spacing_B,M_Yn[nb],M_Wn[nb],Y_ajdk,lf_Yn_ajdk[nb],lf_Wn_ajdk[nb],&(lf_YnWt_vv[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){ }} */}}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZt_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  for (nb=0;nb<nbins;nb++){
	    if (verbose>2){ 
	      raprintf(lf_AnZt_mm[nb]->lf,"double",lf_AnZt_mm[nb]->row_stride,lf_AnZt_mm[nb]->col_stride," %% lf_AnZt_mm: ");
	      raprintf(lf_AnZt_vv[nb]->lf,"double",lf_AnZt_vv[nb]->row_stride,lf_AnZt_vv[nb]->col_stride," %% lf_AnZt_vv: ");
	      raprintf(lf_YnWt_mm[nb]->lf,"double",lf_YnWt_mm[nb]->row_stride,lf_YnWt_mm[nb]->col_stride," %% lf_YnWt_mm: ");
	      raprintf(lf_YnWt_vv[nb]->lf,"double",lf_YnWt_vv[nb]->row_stride,lf_YnWt_vv[nb]->col_stride," %% lf_YnWt_vv: ");
	      /* if (verbose>2){ } */}
	    printf(" %% lf_AnZt_mm[%d] error %0.16f\n",nb,dra_diff(lf_AnZt_mm[nb]->lf,lf_AnZt_vv[nb]->lf,length_AnZt_mm[nb],1));
	    printf(" %% lf_YnWt_mm[%d] error %0.16f\n",nb,dra_diff(lf_YnWt_mm[nb]->lf,lf_YnWt_vv[nb]->lf,length_YnWt_mm[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_AnZt_vv_test]\n");}
}
