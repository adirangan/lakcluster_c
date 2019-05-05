#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_AttAn_vv(void *vp)
{
  /* This function takes in M_At,M_Tt,M_Yt==M_At and ns_j and calculates :
     output_AttAn[na+ny*A_n_cols+ns*A_n_cols*Y_n_cols] = (At(na,:)-a_An*e_At)*diag(T(:,ns))*(Yn(:,ny)-e_An*a_Yn) = ;
     +  At(na,:) * diag(T(:,ns)) *  Yn(:,ny) ;
     - a_An*e_At * diag(T(:,ns)) *  Yn(:,ny) ;
     -  At(na,:) * diag(T(:,ns)) * e_An*a_Yn ;
     + a_An*e_At * diag(T(:,ns)) * e_An*a_Yn ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  int ns_j = *(int *)(vpra[ip++]);
  struct M_handle *M_Yt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *lf_AtTn = (struct L_handle *)(vpra[ip++]);
  struct L_handle *lf_YtTn = (struct L_handle *)(vpra[ip++]);
  struct L_handle *output_AttAn = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_y = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows)/* rup(M_At->nrows + M_At->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yt->nrows)/* rup(M_Yt->nrows + M_Yt->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int mx_j=0,mx_chunk=0,nt_j=0,nt_b=0,nt_a=0,tab_x=0,ta2_a=0,ta2_y=0,ta2_x=0;
  int /* na_j=0,na_b=0, */ma_a=0,na_j=0,na_b=0,na_a=0,tab_a_stride=0,tab_a=0;
  int /* ny_j=0,ny_b=0,ny_a=0, */ny_j=0,ny_b=0,ny_a=0,tab_y_stride=0,tab_y=0;
  long long int hhAY=0,hhAy=0,hhaY=0,hhay=0,hhay_[M_Tt->rpop_j];
  double output_tmp=0;
  unsigned char *A_tag=NULL,*T_tag=NULL,*Y_tag=NULL;
  __m128i *wAt_tag=NULL;
  __m128i *wTt_tag=NULL;
  __m128i *wYt_tag=NULL;
  __m128i *mcay_tag=NULL,*mcay_end=NULL;
  __m128i *mctx_tag=NULL,*mctx_end=NULL;
  unsigned int *ma_b_,*ma_a_;
  unsigned int *na_b_,*na_a_;
  unsigned int *ny_b_,*ny_a_;
  /* unsigned int *ny_b_,*ny_a_; */
  unsigned int *nt_b_,*nt_a_;
  int yper = maximum(1,minimum(M_Yt->rpop_j,maximum(1,/*262144*/ GLOBAL_2_cache_size / maximum(1,M_Yt->mc_length) * 3/4 )));
  int ny_block=0,ny_block_max=0;
  ny_block_max = rup(M_Yt->rpop_j,yper)/yper;
  if (verbose>1){ printf(" %% [entering get_AttAn_vv] tidx %d (yper %d,ny_block_max %d)\n",tidx,(int)yper,(int)(ny_block_max));}  
  if (verbose>1){ printf(" %% M_Yt->rpop_j %d, M_Yt->mc_length %d<--%d, yper = (l2 %d)/%d*3/4 = %d, ny_block_max %d\n",M_Yt->rpop_j,M_Yt->mc_length,M_Yt->ncols,(int)GLOBAL_2_cache_size,(int)M_Yt->mc_length,yper,ny_block_max);}
  if (verbose>2){ raprintf(D_An,"double",1,A_pcols," %% D_An: "); raprintf(a_An,"double",1,A_pcols," %% a_An: ");}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Yt->mr_b,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_b: "); bprintf(M_Yt->mr_j,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_j: ");}
  if (verbose>2){ bprintf(M_Yt->mc_b,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_b: "); bprintf(M_Yt->mc_j,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_j: ");}
  if (verbose>2){ raprintf(D_Yn,"double",1,Y_pcols," %% D_Yn: "); raprintf(a_Yn,"double",1,Y_pcols," %% a_Yn: ");}
  na_b_ = M_At->m_b_; na_a_ = M_At->m_a_;
  ma_b_ = M_At->n_b_; ma_a_ = M_At->n_a_;
  ny_b_ = M_Yt->m_b_; ny_a_ = M_Yt->m_a_;
  /* ny_b_ = M_Yt->n_b_; ny_a = M_Yt->n_a_; */
  nt_b_ = M_Tt->m_b_; nt_a_ = M_Tt->m_a_;
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_y){ case SPACING_j: tab_y_stride = M_Yt->rpop_j; break; case SPACING_b: tab_y_stride = M_Yt->rpop_b; break; case SPACING_a: tab_y_stride = M_Yt->nrows; break; default: break; /* switch (output_spacing_y){ } */}
  if (verbose>2){ printf(" %% tab_a_stride %d tab_y_stride %d\n",tab_a_stride,tab_y_stride);}
  output_AttAn->spacing_row = output_spacing_a; output_AttAn->row_stride = tab_a_stride;
  output_AttAn->spacing_col = output_spacing_y; output_AttAn->col_stride = tab_y_stride;
  fill_uchar_zero((unsigned char *)(output_AttAn->lf),tab_a_stride*tab_y_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"AttAn_vv")){ goto skip_AttAn_vv;}
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    nt_j=ns_j;
    nt_a = nt_a_[nt_j]; nt_b = nt_b_[nt_j];
    wTt_tag = (__m128i *)(&(M_Tt->wX[nt_b*M_Tt->mc_length]));
    mctx_tag = (__m128i *)(&(M_Tt->mc_j[0]));
    mctx_end = (__m128i *)(&(M_Tt->mc_j[M_Tt->mc_length]));
    hhay = -M_Tt->cpop_j + 2*popcount(&wTt_tag,&mctx_tag,&mctx_end); 
    GLOBAL_ops_count_one(tidx,1,M_Tt->mc_length*BIT8);
    ny_block=0;
    while (ny_block<rup(M_Yt->rpop_j,yper)/yper){
      na_j=0;
      while (na_j<M_At->rpop_j){
	na_a = na_a_[na_j]; na_b = na_b_[na_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	switch (output_spacing_y){ case SPACING_j: ta2_y=na_j; break; case SPACING_b: ta2_y=na_b; break; case SPACING_a: ta2_y=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	hhAy = *L2_get(lf_AtTn,na_j,na_b,na_a,nt_j,nt_b,nt_a);
	ny_j=ny_block*yper;
	while (ny_j<minimum(M_Yt->rpop_j,ny_block*yper+yper)){
	  ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j];
	  if (na_j<=ny_j){
	    switch (output_spacing_y){ case SPACING_j: tab_y=ny_j; break; case SPACING_b: tab_y=ny_b; break; case SPACING_a: tab_y=ny_a; break; default: break; /* switch (output_spacing_y){ } */}
	    switch (output_spacing_a){ case SPACING_j: ta2_a=ny_j; break; case SPACING_b: ta2_a=ny_b; break; case SPACING_a: ta2_a=ny_a; break; default: break; /* switch (output_spacing_y){ } */}
	    hhaY = *L2_get(lf_YtTn,ny_j,ny_b,ny_a,nt_j,nt_b,nt_a);
	    tab_x = tab_a + tab_y*tab_a_stride;
	    ta2_x = ta2_a + ta2_y*tab_a_stride;
	    hhAY = +(long long int)(M_At->cpop_j);
	    wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	    wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	    wYt_tag = (__m128i*)&(M_Yt->wX[ny_b*M_Yt->mc_length]);
	    mcay_tag = (__m128i*)&(M_At->mc_j[0]); mcay_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	    hhAY -= (long long int)(2*popcount_notxorxor(&wAt_tag,&wTt_tag,&wYt_tag,&mcay_tag,&mcay_end));
	    output_tmp = /* D_An[na_a/POPLENGTH]* */(hhAY - hhAy*a_Yn[ny_a/POPLENGTH] - a_An[na_a/POPLENGTH]*hhaY + a_An[na_a/POPLENGTH]*hhay*a_Yn[ny_a/POPLENGTH])/* *D_Yn[ny_a/POPLENGTH] */;
	    output_AttAn->lf[tab_x] = output_tmp;
	    output_AttAn->lf[ta2_x] = output_tmp;
	    /* if (na_j<=ny_j){ } */}
	  ny_j++; /* while (ny_j<M_Yt->rpop_j){ } */}
	na_j++; /* while (na_j<M_At->rpop_j){ } */}
      ny_block++; /* while (ny_block<rup(M_Yt->rpop_j,yper)/yper){ } */}
    GLOBAL_ops_count_one(tidx,M_Yt->rpop_j*M_At->rpop_j*3,M_Yt->rpop_j*M_At->rpop_j*M_At->mc_length*BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  if (GLOBAL_omp_type==GLOBAL_omp__on){
    nt_j = ns_j; nt_a = nt_a_[nt_j]; nt_b = nt_b_[nt_j];
    ny_block_max = rup(M_Yt->rpop_j,yper)/yper;
#pragma omp parallel private(mx_j,ny_block,na_j,na_a,na_b,ny_j,ny_a,ny_b,tab_y,tab_a,tab_x,ta2_y,ta2_a,ta2_x,hhay,hhaY,hhAy,hhAY,wAt_tag,wTt_tag,mctx_tag,mctx_end,wYt_tag,mcay_tag,mcay_end,output_tmp)
    { /* begin omp parallel */
      mx_j=0; tab_y=0; tab_a=0;
#pragma omp for schedule(dynamic)
      for (mx_j=0;mx_j<ny_block_max;mx_j++){
	ny_block = mx_j;	
	switch (output_spacing_y){ case SPACING_j: ta2_y=na_j; break; case SPACING_b: ta2_y=na_b; break; case SPACING_a: ta2_y=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	wTt_tag = (__m128i *)(&(M_Tt->wX[nt_b*M_Tt->mc_length]));
	mctx_tag = (__m128i *)(&(M_Tt->mc_j[0]));
	mctx_end = (__m128i *)(&(M_Tt->mc_j[M_Tt->mc_length]));
	hhay = -M_Tt->cpop_j + 2*popcount(&wTt_tag,&mctx_tag,&mctx_end); 
	GLOBAL_ops_count_one(tidx,1,M_Tt->mc_length*BIT8);
	na_j=0;
	while (na_j<M_At->rpop_j){
	  na_a = na_a_[na_j]; na_b = na_b_[na_j];
	  switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	  switch (output_spacing_y){ case SPACING_j: ta2_y=na_j; break; case SPACING_b: ta2_y=na_b; break; case SPACING_a: ta2_y=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	  hhAy = *L2_get(lf_AtTn,na_j,na_b,na_a,nt_j,nt_b,nt_a);
	  ny_j=ny_block*yper;
	  while (ny_j<minimum(M_Yt->rpop_j,ny_block*yper+yper)){
	    ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j];
	    if (na_j<=ny_j){
	      switch (output_spacing_y){ case SPACING_j: tab_y=ny_j; break; case SPACING_b: tab_y=ny_b; break; case SPACING_a: tab_y=ny_a; break; default: break; /* switch (output_spacing_y){ } */}
	      switch (output_spacing_a){ case SPACING_j: ta2_a=ny_j; break; case SPACING_b: ta2_a=ny_b; break; case SPACING_a: ta2_a=ny_a; break; default: break; /* switch (output_spacing_y){ } */}
	      hhaY = *L2_get(lf_YtTn,ny_j,ny_b,ny_a,nt_j,nt_b,nt_a);
	      tab_x = tab_a + tab_y*tab_a_stride;
	      ta2_x = ta2_a + ta2_y*tab_a_stride;
	      hhAY = +(long long int)(M_At->cpop_j);
	      wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	      wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	      wYt_tag = (__m128i*)&(M_Yt->wX[ny_b*M_Yt->mc_length]);
	      mcay_tag = (__m128i*)&(M_At->mc_j[0]); mcay_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	      hhAY -= (long long int)(2*popcount_notxorxor(&wAt_tag,&wTt_tag,&wYt_tag,&mcay_tag,&mcay_end));
	      output_tmp = /* D_An[na_a/POPLENGTH]* */(hhAY - hhAy*a_Yn[ny_a/POPLENGTH] - a_An[na_a/POPLENGTH]*hhaY + a_An[na_a/POPLENGTH]*hhay*a_Yn[ny_a/POPLENGTH])/* *D_Yn[ny_a/POPLENGTH] */;
	      output_AttAn->lf[tab_x] = output_tmp;
	      output_AttAn->lf[ta2_x] = output_tmp;
	      /* if (na_j<=ny_j){ } */}
	    ny_j++; /* while (ny_j<M_Yt->rpop_j){ } */}
	  na_j++; /* while (na_j<M_At->rpop_j){ } */}
	/* for (mx_j=0;mx_j<ny_block_max;mx_j++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,M_Yt->rpop_j*M_At->rpop_j*3,M_Yt->rpop_j*M_At->rpop_j*M_At->mc_length*BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
 skip_AttAn_vv:
  if (verbose>1){ printf(" %% [finished get_AttAn_vv] tidx %d\n",tidx);}
  return NULL;
}

int wrap_AttAn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_y,struct M_handle *M_At,struct M_handle *M_Tt,int ns_j,struct M_handle *M_Yt,double *A_ajdk,double *Y_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_YtTn,struct L_handle **output_AttAn_p)
{
  /* This function uses the M_handles M_At, M_Tt, M_Yt and ns_j to run a series of parallel calls to get_AttAn_vv ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 12)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_a=0,length_y=0,length_t=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_AttAn_vv__run] tidx %d \n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Yt,verbose," %% M_Yt: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_y){ case SPACING_j: length_y = M_Yt->rpop_j; break; case SPACING_b: length_y = M_Yt->rpop_b; break; case SPACING_a: length_y = M_Yt->nrows; break; default: break; /* switch (output_spacing_y){ } */}
  length = length_a*length_y; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_a,length_y,length);}
  length = length_a*length_y; if (*output_AttAn_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_AttAn_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: ");}
  if (verbose>2){ bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: ");}
  if (verbose>2){ bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Yt->mr_b,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_b: ");}
  if (verbose>2){ bprintf(M_Yt->mr_j,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_j: ");}
  if (verbose>2){ bprintf(M_Yt->mc_b,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_b: ");}
  if (verbose>2){ bprintf(M_Yt->mc_j,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_j: ");}
  length = length_a*length_y; if ((*output_AttAn_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_AttAn_vv__run\n",(*output_AttAn_p)->length,length);} memset((*output_AttAn_p)->lf,0,length*sizeof(double));
  if (ns_j<0 || ns_j>addressable_int_length-1){ printf(" %% Warning! ns_j %d in wrap_AttAn_vv_run\n",ns_j);}
  ip=0;
  vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = &(addressable_int[ns_j]); vpra[ip++] = M_Yt; vpra[ip++] = A_ajdk; vpra[ip++] = Y_ajdk; vpra[ip++] = lf_AtTn; vpra[ip++] = lf_YtTn; vpra[ip++] = *output_AttAn_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_y){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_y){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AttAn_vv,vpra)){ printf("Warning! cannot create thread %d in wrap_AttAn\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AttAn_vv(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_AttAn_p)->lf,"double",length_a,length_y*length_t," %% (*output_AttAn_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_AttAn_vv__run] tidx %d\n",*tidx);}
  return length;
}  
