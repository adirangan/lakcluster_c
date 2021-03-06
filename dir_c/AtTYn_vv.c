#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_AtTYn_vv(void *vp)
{
  /* This function takes in M_At,M_Tt,M_Yt and calculates :
     output_AtTYn[na+ny*A_n_cols+ns*A_n_cols*Y_n_cols] = (At(na,:)-a_An*e_At)*diag(T(:,ns))*(Yn(:,ny)-e_An*a_Yn) = ;
     +  At(na,:) * diag(T(:,ns)) *  Yn(:,ny) ;
     - a_An*e_At * diag(T(:,ns)) *  Yn(:,ny) ;
     -  At(na,:) * diag(T(:,ns)) * e_An*a_Yn ;
     + a_An*e_At * diag(T(:,ns)) * e_An*a_Yn ;
     For compatibility with earlier versions of AtTYn____WtSZn it may be more convenient to compute;
     output_AtTYn[na+ny*A_n_cols+ns*A_n_cols*Y_n_cols] = (At(na,:) - 0 )*diag(T(:,ns))*(Yn(:,ny)-e_An*a_Yn) = ;
     +  At(na,:) * diag(T(:,ns)) *  Yn(:,ny) ;
     -  At(na,:) * diag(T(:,ns)) * e_An*a_Yn ;
     however, we do not compute this here (see commented-out "//" computations below). ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Yt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *lf_AtTn = (struct L_handle *)(vpra[ip++]);
  struct L_handle *lf_YtTn = (struct L_handle *)(vpra[ip++]);
  struct L_handle *output_AtTYn = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_y = *(int *)(vpra[ip++]);
  int output_spacing_t = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows)/* rup(M_At->nrows + M_At->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yt->nrows)/* rup(M_Yt->nrows + M_Yt->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int mx_j=0,mx_chunk=0,nt_j=0,nt_b=0,nt_a=0,tab_t_stride=0,tab_t=0,tab_x=0;
  int /* na_j=0,na_b=0, */ma_a=0,na_j=0,na_b=0,na_a=0,tab_a_stride=0,tab_a=0;
  int /* ny_j=0,ny_b=0,ny_a=0, */ny_j=0,ny_b=0,ny_a=0,tab_y_stride=0,tab_y=0;
  long long int hhAY=0,hhAy=0,hhaY=0,hhay=0,hhay_[M_Tt->rpop_j];
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
  if (verbose>1){ printf(" %% [entering get_AtTYn_vv] tidx %d (yper %d,ny_block_max %d)\n",tidx,(int)yper,(int)(ny_block_max));}  
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
  switch (output_spacing_t){ case SPACING_j: tab_t_stride = M_Tt->rpop_j; break; case SPACING_b: tab_t_stride = M_Tt->rpop_b; break; case SPACING_a: tab_t_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_y){ case SPACING_j: tab_y_stride = M_Yt->rpop_j; break; case SPACING_b: tab_y_stride = M_Yt->rpop_b; break; case SPACING_a: tab_y_stride = M_Yt->nrows; break; default: break; /* switch (output_spacing_y){ } */}
  if (verbose>2){ printf(" %% tab_t_stride %d tab_a_stride %d tab_y_stride %d\n",tab_t_stride,tab_a_stride,tab_y_stride);}
  output_AtTYn->spacing_row = output_spacing_a; output_AtTYn->row_stride = tab_a_stride;
  output_AtTYn->spacing_col = output_spacing_y; output_AtTYn->col_stride = tab_y_stride;
  output_AtTYn->spacing_lyr = output_spacing_t; output_AtTYn->lyr_stride = tab_t_stride;
  fill_uchar_zero((unsigned char *)(output_AtTYn->lf),tab_a_stride*tab_y_stride*tab_t_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"AtTYn_vv")){ goto skip_AtTYn_vv;}
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    nt_j=0;
    while (nt_j<M_Tt->rpop_j){
      nt_a = nt_a_[nt_j]; nt_b = nt_b_[nt_j];
      switch (output_spacing_t){ case SPACING_j: tab_t=nt_j; break; case SPACING_b: tab_t=nt_b; break; case SPACING_a: tab_t=nt_a; break; default: break; /* switch (output_spacing_t){ } */}
      wTt_tag = (__m128i *)(&(M_Tt->wX[nt_b*M_Tt->mc_length]));
      mctx_tag = (__m128i *)(&(M_Tt->mc_j[0]));
      mctx_end = (__m128i *)(&(M_Tt->mc_j[M_Tt->mc_length]));
      hhay = -M_Tt->cpop_j + 2*popcount(&wTt_tag,&mctx_tag,&mctx_end); 
      GLOBAL_ops_count_one(tidx,1,(unsigned long long int)M_Tt->mc_length*BIT8);
      ny_block=0;
      while (ny_block<rup(M_Yt->rpop_j,yper)/yper){
      	na_j=0;
	while (na_j<M_At->rpop_j){
	  na_a = na_a_[na_j]; na_b = na_b_[na_j];
	  switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	  hhAy = *L2_get(lf_AtTn,na_j,na_b,na_a,nt_j,nt_b,nt_a);
	  ny_j=ny_block*yper;
	  while (ny_j<minimum(M_Yt->rpop_j,ny_block*yper+yper)){
	    ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j];
	    switch (output_spacing_y){ case SPACING_j: tab_y=ny_j; break; case SPACING_b: tab_y=ny_b; break; case SPACING_a: tab_y=ny_a; break; default: break; /* switch (output_spacing_y){ } */}
	    hhaY = *L2_get(lf_YtTn,ny_j,ny_b,ny_a,nt_j,nt_b,nt_a);
	    tab_x = tab_a + tab_y*tab_a_stride + tab_t*tab_a_stride*tab_y_stride;
	    hhAY = +(long long int)(M_At->cpop_j);
	    wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	    wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	    wYt_tag = (__m128i*)&(M_Yt->wX[ny_b*M_Yt->mc_length]);
	    mcay_tag = (__m128i*)&(M_At->mc_j[0]); mcay_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	    hhAY -= (long long int)(2*popcount_notxorxor(&wAt_tag,&wTt_tag,&wYt_tag,&mcay_tag,&mcay_end));
	    output_AtTYn->lf[tab_x] = /* D_An[na_a/POPLENGTH]* */(hhAY - hhAy*a_Yn[ny_a/POPLENGTH] - a_An[na_a/POPLENGTH]*hhaY + a_An[na_a/POPLENGTH]*hhay*a_Yn[ny_a/POPLENGTH])/* *D_Yn[ny_a/POPLENGTH] */;
	    //output_AtTYn->lf[tab_x] = /* D_An[na_a/POPLENGTH]* */(hhAY - hhAy*a_Yn[ny_a/POPLENGTH] /* - a_An[na_a/POPLENGTH]*hhaY + a_An[na_a/POPLENGTH]*hhay*a_Yn[ny_a/POPLENGTH] */)/* *D_Yn[ny_a/POPLENGTH] */;
	    ny_j++; /* while (ny_j<M_Yt->rpop_j){ } */}
	  na_j++; /* while (na_j<M_At->rpop_j){ } */}
	ny_block++; /* while (ny_block<rup(M_Yt->rpop_j,yper)/yper){ } */}
      GLOBAL_ops_count_one(tidx,(unsigned long long int)M_Yt->rpop_j*(unsigned long long int)M_At->rpop_j*3,(unsigned long long int)M_Yt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_At->mc_length*BIT8);
      nt_j++; /* while (nt_j<M_Tt->rpop_j){ } */}
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  if (GLOBAL_omp_type==GLOBAL_omp__on){
    ny_block_max = rup(M_Yt->rpop_j,yper)/yper;
#pragma omp parallel private(mx_j,ny_block,nt_j,nt_a,nt_b,na_j,na_a,na_b,ny_j,ny_a,ny_b,tab_t,tab_y,tab_a,tab_x,hhay,hhaY,hhAy,hhAY,wAt_tag,wTt_tag,mctx_tag,mctx_end,wYt_tag,mcay_tag,mcay_end)
    { /* begin omp parallel */
      mx_j=0; tab_t=0; tab_y=0; tab_a=0;
#pragma omp for schedule(dynamic)
      for (mx_j=0;mx_j<ny_block_max*M_Tt->rpop_j;mx_j++){
	nt_j = mx_j / ny_block_max;
	ny_block = mx_j % ny_block_max;
	nt_a = nt_a_[nt_j]; nt_b = nt_b_[nt_j];
	switch (output_spacing_t){ case SPACING_j: tab_t=nt_j; break; case SPACING_b: tab_t=nt_b; break; case SPACING_a: tab_t=nt_a; break; default: break; /* switch (output_spacing_t){ } */}
	wTt_tag = (__m128i *)(&(M_Tt->wX[nt_b*M_Tt->mc_length]));
	mctx_tag = (__m128i *)(&(M_Tt->mc_j[0]));
	mctx_end = (__m128i *)(&(M_Tt->mc_j[M_Tt->mc_length]));
	hhay = -M_Tt->cpop_j + 2*popcount(&wTt_tag,&mctx_tag,&mctx_end); 
	GLOBAL_ops_count_one(tidx,1,(unsigned long long int)M_Tt->mc_length*BIT8);
	na_j=0;
	while (na_j<M_At->rpop_j){
	  na_a = na_a_[na_j]; na_b = na_b_[na_j];
	  switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	  hhAy = *L2_get(lf_AtTn,na_j,na_b,na_a,nt_j,nt_b,nt_a);
	  ny_j=ny_block*yper;
	  while (ny_j<minimum(M_Yt->rpop_j,ny_block*yper+yper)){
	    ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j];
	    switch (output_spacing_y){ case SPACING_j: tab_y=ny_j; break; case SPACING_b: tab_y=ny_b; break; case SPACING_a: tab_y=ny_a; break; default: break; /* switch (output_spacing_y){ } */}
	    hhaY = *L2_get(lf_YtTn,ny_j,ny_b,ny_a,nt_j,nt_b,nt_a);
	    tab_x = tab_a + tab_y*tab_a_stride + tab_t*tab_a_stride*tab_y_stride;
	    hhAY = +(long long int)(M_At->cpop_j);
	    wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	    wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	    wYt_tag = (__m128i*)&(M_Yt->wX[ny_b*M_Yt->mc_length]);
	    mcay_tag = (__m128i*)&(M_At->mc_j[0]); mcay_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	    hhAY -= (long long int)(2*popcount_notxorxor(&wAt_tag,&wTt_tag,&wYt_tag,&mcay_tag,&mcay_end));
	    output_AtTYn->lf[tab_x] = /* D_An[na_a/POPLENGTH]* */(hhAY - hhAy*a_Yn[ny_a/POPLENGTH] - a_An[na_a/POPLENGTH]*hhaY + a_An[na_a/POPLENGTH]*hhay*a_Yn[ny_a/POPLENGTH])/* *D_Yn[ny_a/POPLENGTH] */;
	    //output_AtTYn->lf[tab_x] = /* D_An[na_a/POPLENGTH]* */(hhAY - hhAy*a_Yn[ny_a/POPLENGTH] /* - a_An[na_a/POPLENGTH]*hhaY + a_An[na_a/POPLENGTH]*hhay*a_Yn[ny_a/POPLENGTH] */)/* *D_Yn[ny_a/POPLENGTH] */;
	    ny_j++; /* while (ny_j<M_Yt->rpop_j){ } */}
	  na_j++; /* while (na_j<M_At->rpop_j){ } */}
	/* for (mx_j=0;mx_j<M_Tt->rpop_j*ny_block_max;mx_j++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_Yt->rpop_j*(unsigned long long int)M_At->rpop_j*3,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_Yt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_At->mc_length*BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp_unused){
    for (nt_j=0;nt_j<M_Tt->rpop_j;nt_j++){
      nt_a = nt_a_[nt_j]; nt_b = nt_b_[nt_j];
      switch (output_spacing_t){ case SPACING_j: tab_t=nt_j; break; case SPACING_b: tab_t=nt_b; break; case SPACING_a: tab_t=nt_a; break; default: break; /* switch (output_spacing_t){ } */}
      wTt_tag = (__m128i *)(&(M_Tt->wX[nt_b*M_Tt->mc_length]));
      mctx_tag = (__m128i *)(&(M_Tt->mc_j[0]));
      mctx_end = (__m128i *)(&(M_Tt->mc_j[M_Tt->mc_length]));
      hhay_[nt_j] = -M_Tt->cpop_j + 2*popcount(&wTt_tag,&mctx_tag,&mctx_end); 
      /* for (nt_j=0;nt_j<M_Tt->rpop_j;nt_j++){ } */}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)M_Tt->rpop_j,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_Tt->mc_length*BIT8);
    mx_chunk=1920;
#pragma omp parallel shared(mx_chunk) private(mx_j,nt_j,nt_a,nt_b,na_j,na_a,na_b,ny_j,ny_a,ny_b,tab_t,tab_y,tab_a,tab_x,hhaY,hhAy,hhAY,wAt_tag,wTt_tag,wYt_tag,mcay_tag,mcay_end)
    { /* begin omp parallel */
      mx_j=0; tab_t=0; tab_y=0; tab_a=0;
#pragma omp for schedule(dynamic,mx_chunk)
      for (mx_j=0;mx_j<M_Tt->rpop_j*M_Yt->rpop_j*M_At->rpop_j;mx_j++){
	nt_j = mx_j / (M_Yt->rpop_j*M_At->rpop_j); 
	ny_j = (mx_j % (M_Yt->rpop_j*M_At->rpop_j)) / M_At->rpop_j ;
	na_j = (mx_j % M_At->rpop_j) ;
	nt_a = nt_a_[nt_j]; nt_b = nt_b_[nt_j];
	switch (output_spacing_t){ case SPACING_j: tab_t=nt_j; break; case SPACING_b: tab_t=nt_b; break; case SPACING_a: tab_t=nt_a; break; default: break; /* switch (output_spacing_t){ } */}
	ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j];
	switch (output_spacing_y){ case SPACING_j: tab_y=ny_j; break; case SPACING_b: tab_y=ny_b; break; case SPACING_a: tab_y=ny_a; break; default: break; /* switch (output_spacing_y){ } */}
	hhaY = *L2_get(lf_YtTn,ny_j,ny_b,ny_a,nt_j,nt_b,nt_a);
	na_a = na_a_[na_j]; na_b = na_b_[na_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	hhAy = *L2_get(lf_AtTn,na_j,na_b,na_a,nt_j,nt_b,nt_a);
	tab_x = tab_a + tab_y*tab_a_stride + tab_t*tab_a_stride*tab_y_stride;
	hhAY = +(long long int)(M_At->cpop_j);
	wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	wYt_tag = (__m128i*)&(M_Yt->wX[ny_b*M_Yt->mc_length]);
	mcay_tag = (__m128i*)&(M_At->mc_j[0]); mcay_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	hhAY -= (long long int)(2*popcount_notxorxor(&wAt_tag,&wTt_tag,&wYt_tag,&mcay_tag,&mcay_end));
	output_AtTYn->lf[tab_x] = /* D_An[na_a/POPLENGTH]* */(hhAY - hhAy*a_Yn[ny_a/POPLENGTH] - a_An[na_a/POPLENGTH]*hhaY + a_An[na_a/POPLENGTH]*hhay_[nt_j]*a_Yn[ny_a/POPLENGTH])/* *D_Yn[ny_a/POPLENGTH] */;
	/* for (mx_j=0;mx_j<M_Tt->rpop_j*M_Yt->rpop_j*M_At->rpop_j;mx_j++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_Yt->rpop_j*(unsigned long long int)M_At->rpop_j*3,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_Yt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_At->mc_length*BIT8);
    /* else if (GLOBAL_omp_type==GLOBAL_omp_unused){ } */}
 skip_AtTYn_vv:
  if (verbose>1){ printf(" %% [finished get_AtTYn_vv] tidx %d\n",tidx);}
  return NULL;
}

int wrap_AtTYn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_y,int output_spacing_t,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yt,double *A_ajdk,double *Y_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_YtTn,struct L_handle **output_AtTYn_p)
{
  /* This function uses the M_handles M_At, M_Tt and M_Yt to run a series of parallel calls to get_AtTYn_vv ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 13)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_a=0,length_y=0,length_t=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_AtTYn_vv__run] tidx %d \n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Yt,verbose," %% M_Yt: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_y){ case SPACING_j: length_y = M_Yt->rpop_j; break; case SPACING_b: length_y = M_Yt->rpop_b; break; case SPACING_a: length_y = M_Yt->nrows; break; default: break; /* switch (output_spacing_y){ } */}
  switch (output_spacing_t){ case SPACING_j: length_t = M_Tt->rpop_j; break; case SPACING_b: length_t = M_Tt->rpop_b; break; case SPACING_a: length_t = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  length = length_a*length_y*length_t; if (verbose){ printf(" %% length %llu*%llu*%llu=%llu\n",length_a,length_y,length_t,length);}
  length = length_a*length_y*length_t; if (*output_AtTYn_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_AtTYn_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: ");}
  if (verbose>2){ bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: ");}
  if (verbose>2){ bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Yt->mr_b,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_b: ");}
  if (verbose>2){ bprintf(M_Yt->mr_j,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_j: ");}
  if (verbose>2){ bprintf(M_Yt->mc_b,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_b: ");}
  if (verbose>2){ bprintf(M_Yt->mc_j,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_j: ");}
  length = length_a*length_y*length_t; if ((*output_AtTYn_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_AtTYn_vv__run\n",(*output_AtTYn_p)->length,length);} memset((*output_AtTYn_p)->lf,0,length*sizeof(double));
  ip=0;
  vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Yt; vpra[ip++] = A_ajdk; vpra[ip++] = Y_ajdk; vpra[ip++] = lf_AtTn; vpra[ip++] = lf_YtTn; vpra[ip++] = *output_AtTYn_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_y){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_y){ } */}
  switch (output_spacing_t){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_t){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AtTYn_vv,vpra)){ printf("Warning! cannot create thread %d in wrap_AtTYn\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AtTYn_vv(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_AtTYn_p)->lf,"double",length_a,length_y*length_t," %% (*output_AtTYn_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_AtTYn_vv__run] tidx %d\n",*tidx);}
  return length;
}  

void *get_AtTYn_uu(void *vp)
{
  /* This function takes in M_At,M_Tt,M_Yt and calculates :
     output_AtTYn[na+ny*A_n_cols+ns*A_n_cols*Y_n_cols] = (At(na,:)-a_An*e_At)*diag(T(:,ns))*(Yn(:,ny)-e_An*a_Yn) = ;
     +  At(na,:) * diag(T(:,ns)) *  Yn(:,ny) ;
     - a_An*e_At * diag(T(:,ns)) *  Yn(:,ny) ;
     -  At(na,:) * diag(T(:,ns)) * e_An*a_Yn ;
     + a_An*e_At * diag(T(:,ns)) * e_An*a_Yn ;
     For compatibility with earlier versions of AtTYn____WtSZn it may be more convenient to compute;
     output_AtTYn[na+ny*A_n_cols+ns*A_n_cols*Y_n_cols] = (At(na,:) - 0 )*diag(T(:,ns))*(Yn(:,ny)-e_An*a_Yn) = ;
     +  At(na,:) * diag(T(:,ns)) *  Yn(:,ny) ;
     -  At(na,:) * diag(T(:,ns)) * e_An*a_Yn ;
     however, we do not compute this here (see commented-out "//" computations below). ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Yt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_AtTYn = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_y = *(int *)(vpra[ip++]);
  int output_spacing_t = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows)/* rup(M_At->nrows + M_At->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yt->nrows)/* rup(M_Yt->nrows + M_Yt->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int nt_j=0,nt_b=0,nt_a=0,tab_t_stride=0,tab_t=0,tab_x=0;
  int /* na_j=0,na_b=0, */ma_a=0,na_j=0,na_b=0,na_a=0,tab_a_stride=0,tab_a=0;
  int /* ny_j=0,ny_b=0,ny_a=0, */ny_j=0,ny_b=0,ny_a=0,tab_y_stride=0,tab_y=0;
  long long int hhAY=0,hhAy=0,hhaY=0,hhay=0,hhay_[M_Tt->rpop_j];
  unsigned char *A_tag=NULL,*T_tag=NULL,*Y_tag=NULL;
  unsigned int *ma_b_,*ma_a_;
  unsigned int *na_b_,*na_a_;
  unsigned int *ny_b_,*ny_a_;
  /* unsigned int *ny_b_,*ny_a_; */
  unsigned int *nt_b_,*nt_a_;
  double dtmp=0;
  if (verbose>1){ printf(" %% [entering get_AtTYn_uu] tidx %d\n",tidx);}
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
  switch (output_spacing_t){ case SPACING_j: tab_t_stride = M_Tt->rpop_j; break; case SPACING_b: tab_t_stride = M_Tt->rpop_b; break; case SPACING_a: tab_t_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_y){ case SPACING_j: tab_y_stride = M_Yt->rpop_j; break; case SPACING_b: tab_y_stride = M_Yt->rpop_b; break; case SPACING_a: tab_y_stride = M_Yt->nrows; break; default: break; /* switch (output_spacing_y){ } */}
  if (verbose>2){ printf(" %% tab_t_stride %d tab_a_stride %d tab_y_stride %d\n",tab_t_stride,tab_a_stride,tab_y_stride);}
  output_AtTYn->spacing_row = output_spacing_a; output_AtTYn->row_stride = tab_a_stride;
  output_AtTYn->spacing_col = output_spacing_y; output_AtTYn->col_stride = tab_y_stride;
  output_AtTYn->spacing_lyr = output_spacing_t; output_AtTYn->lyr_stride = tab_t_stride;
  fill_uchar_zero((unsigned char *)(output_AtTYn->lf),tab_a_stride*tab_y_stride*tab_t_stride*sizeof(double));
  nt_j=0;
  while (nt_j<M_Tt->rpop_j){
    nt_a = nt_a_[nt_j]; nt_b = nt_b_[nt_j];
    switch (output_spacing_t){ case SPACING_j: tab_t=nt_j; break; case SPACING_b: tab_t=nt_b; break; case SPACING_a: tab_t=nt_a; break; default: break; /* switch (output_spacing_t){ } */}
    T_tag = (unsigned char *)(&(M_Tt->wX[nt_b*M_Tt->mc_length]));
    ny_j=0;
    while (ny_j<M_Yt->rpop_j){
      ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j];
      switch (output_spacing_y){ case SPACING_j: tab_y=ny_j; break; case SPACING_b: tab_y=ny_b; break; case SPACING_a: tab_y=ny_a; break; default: break; /* switch (output_spacing_y){ } */}
      Y_tag = (unsigned char *)(&(M_Yt->wX[ny_b*M_Yt->mc_length]));
      na_j=0;
      while (na_j<M_At->rpop_j){
	na_a = na_a_[na_j]; na_b = na_b_[na_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	A_tag = (unsigned char *)(&(M_At->wX[na_b*M_At->mc_length]));
	tab_x = tab_a + tab_y*tab_a_stride + tab_t*tab_a_stride*tab_y_stride;
	dtmp=0;
	for (ma_a=0;ma_a<M_At->ncols;ma_a++){
	  dtmp += (bget____(A_tag,ma_a) - a_An[na_a/POPLENGTH])*bget____(T_tag,ma_a)*(bget____(Y_tag,ma_a) - a_Yn[ny_a/POPLENGTH])*bget__on(M_At->mc_j,ma_a);
	  /* for (ma_a=0;ma_a<M_At->ncols;ma_a++){ } */}
	output_AtTYn->lf[tab_x] = dtmp;
	na_j++;/* while (na_j<M_At->rpop_j){ } */}
      ny_j++;/* while (ny_j<M_Yt->rpop_j){ } */}
    nt_j++;/* while (nt_j<M_Tt->rpop_j){ } */}
  if (verbose>1){ printf(" %% [finished get_AtTYn_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_AtTYn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_y,int output_spacing_t,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yt,double *A_ajdk,double *Y_ajdk,struct L_handle **output_AtTYn_p)
{
  /* This function uses the M_handles M_At, M_Tt and M_Yt to run a series of parallel calls to get_AtTYn_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 10)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_a=0,length_y=0,length_t=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_AtTYn_uu__run] tidx %d \n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Yt,verbose," %% M_Yt: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_y){ case SPACING_j: length_y = M_Yt->rpop_j; break; case SPACING_b: length_y = M_Yt->rpop_b; break; case SPACING_a: length_y = M_Yt->nrows; break; default: break; /* switch (output_spacing_y){ } */}
  switch (output_spacing_t){ case SPACING_j: length_t = M_Tt->rpop_j; break; case SPACING_b: length_t = M_Tt->rpop_b; break; case SPACING_a: length_t = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  length = length_a*length_y*length_t; if (verbose){ printf(" %% length %llu*%llu*%llu=%llu\n",length_a,length_y,length_t,length);}
  length = length_a*length_y*length_t; if (*output_AtTYn_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_AtTYn_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: ");}
  if (verbose>2){ bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: ");}
  if (verbose>2){ bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Yt->mr_b,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_b: ");}
  if (verbose>2){ bprintf(M_Yt->mr_j,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_j: ");}
  if (verbose>2){ bprintf(M_Yt->mc_b,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_b: ");}
  if (verbose>2){ bprintf(M_Yt->mc_j,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_j: ");}
  length = length_a*length_y*length_t; if ((*output_AtTYn_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_AtTYn_uu__run\n",(*output_AtTYn_p)->length,length);} memset((*output_AtTYn_p)->lf,0,length*sizeof(double));
  ip=0;
  vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Yt; vpra[ip++] = A_ajdk; vpra[ip++] = Y_ajdk; vpra[ip++] = *output_AtTYn_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_y){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_y){ } */}
  switch (output_spacing_t){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_t){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AtTYn_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_AtTYn\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AtTYn_uu(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_AtTYn_p)->lf,"double",length_a,length_y*length_t," %% (*output_AtTYn_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_AtTYn_uu__run] tidx %d\n",*tidx);}
  return length;
}  

void wrap_AtTYn_vv_test()
{
  /* test for errors with input file: AtTYn_vv_error.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= AtTYn_vv;
     GLOBAL_TEST_TYP2= error;
     GLOBAL_NBINS= 4;
     GLOBAL_TEST_mrand= 0.015;
     GLOBAL_TEST_A_n_rows= 1930;
     GLOBAL_TEST_A_n_cols= 20;
     GLOBAL_TEST_Z_n_rows= 21;
     GLOBAL_TEST_Y_n_cols= 2122;
     GLOBAL_TEST_T_n_cols= 4;
     GLOBAL_TEST_niter= 1;
     END= 0;
  */
  /* test for speed with input file: AtTYn_vv_speed.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= AtTYn_vv;
     GLOBAL_TEST_TYP2= speed;
     GLOBAL_NBINS= 4;
     GLOBAL_TEST_mrand= 0.05;
     GLOBAL_TEST_A_n_rows= 19200;
     GLOBAL_TEST_A_n_cols= 1920;
     GLOBAL_TEST_Z_n_rows= 1920;
     GLOBAL_TEST_Y_n_cols= 19200;
     GLOBAL_TEST_T_n_cols= 2;
     GLOBAL_TEST_niter= 1;
     END= 0;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AtTn=NULL,**lf_YtTn=NULL,**lf_ZtSn=NULL,**lf_WtSn=NULL;
  struct L_handle **lf_AtTYn_vv=NULL,**lf_ZtSWn_vv=NULL; int *length_AtTYn_vv=NULL,*length_ZtSWn_vv=NULL;
  struct L_handle **lf_AtTYn_uu=NULL,**lf_ZtSWn_uu=NULL; int *length_AtTYn_uu=NULL,*length_ZtSWn_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_AtTYn_vv_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  lf_AtTn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_YtTn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_ZtSn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_WtSn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_AtTn[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    lf_YtTn[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    lf_ZtSn[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    lf_WtSn[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_AtTYn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AtTYn_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_ZtSWn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_ZtSWn_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_AtTYn_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Yn[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    lf_ZtSWn_vv[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_Wn[nb]->ncols*(unsigned long long int)M_Sn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (error_check){ 
    lf_AtTYn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AtTYn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    lf_ZtSWn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_ZtSWn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    for (nb=0;nb<nbins;nb++){
      lf_AtTYn_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Yn[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
      lf_ZtSWn_uu[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_Wn[nb]->ncols*(unsigned long long int)M_Sn[nb]->ncols);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (error_check){ } */}
  for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
	if (verbose){ printf(" %% %s; %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],NULL,NULL,NULL,&(lf_AtTn[nb]));
	  GLOBAL_pthread_toc(); 
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Yt[nb],M_Tt[nb],NULL,NULL,NULL,&(lf_YtTn[nb]));
	  GLOBAL_pthread_toc(); 
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],NULL,NULL,NULL,&(lf_ZtSn[nb]));
	  GLOBAL_pthread_toc(); 
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Wt[nb],M_St[nb],NULL,NULL,NULL,&(lf_WtSn[nb]));
	  GLOBAL_pthread_toc(); 
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc();
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
	for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic();
	  length_AtTYn_vv[nb] = wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_Yt[nb],A_ajdk,Y_ajdk,lf_AtTn[nb],lf_YtTn[nb],&(lf_AtTYn_vv[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  length_ZtSWn_vv[nb] = wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_Zt[nb],M_St[nb],M_Wt[nb],A_ajdk,Y_ajdk,lf_ZtSn[nb],lf_WtSn[nb],&(lf_ZtSWn_vv[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){ }} */}}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYn_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");  
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){	
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
	  for (nb=0;nb<nbins;nb++){
	      GLOBAL_pthread_tic();
	      length_AtTYn_uu[nb] = wrap_AtTYn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_Yt[nb],A_ajdk,Y_ajdk,&(lf_AtTYn_uu[nb]));
	      GLOBAL_pthread_toc();
	      GLOBAL_pthread_tic();
	      length_ZtSWn_uu[nb] = wrap_AtTYn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_Zt[nb],M_St[nb],M_Wt[nb],A_ajdk,Y_ajdk,&(lf_ZtSWn_uu[nb]));
	      GLOBAL_pthread_toc();
	      /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYn_uu: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");  
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_AtTYn_vv[%d] error %0.16f\n",nb,dra_diff(lf_AtTYn_vv[nb]->lf,lf_AtTYn_uu[nb]->lf,length_AtTYn_vv[nb],1));
	    printf(" %% lf_ZtSWn_vv[%d] error %0.16f\n",nb,dra_diff(lf_ZtSWn_vv[nb]->lf,lf_ZtSWn_uu[nb]->lf,length_ZtSWn_vv[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}	  
	  /* if (error_check){ } */}
	/* for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_AtTYn_vv_test]\n");}
}
