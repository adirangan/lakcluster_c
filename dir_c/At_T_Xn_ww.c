#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_At_T_Xn_ww(void *vp)
{
  /* Here M_Xt is only used to store the column mask for lf_Xn */
  int verbose=1;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Xn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Xt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_At_T_Xn_ww = (struct L_handle *)(vpra[ip++]);
  int trm_flag = *(int *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_t = *(int *)(vpra[ip++]);
  int output_spacing_w = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows)/* rup(M_At->nrows + M_At->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int ns_j=0,ns_a=0,ns_b=0,tab_t=0,tab_t_stride=0,na_j=0,na_a=0,na_b=0,tab_a=0,tab_a_stride=0,mw_j=0,mw_a=0,mw_b=0,tab_w=0,tab_w_stride=0,my_j=0,my_a=0,my_b=0,tab_y=0;
  double output_tmp=0;
  int vA=0,vS=0;
  unsigned char *Tt_tag=NULL,*At_tag=NULL;
  __m128i *wAt_tag=NULL,*wTt_tag=NULL,*wXn_tag=NULL,*mcat_tag=NULL,*mcat_end=NULL;
  long long int lltmp=0,n2=0;
  double output_At_T_Xn_base=0,output_at_T_Xn_base=0,output____T_Xn_base=0,output____T_Xn_base_[M_Tt->rpop_j],dtmp=0,output_At_T_Xn_tmp=0,output_at_T_Xn_tmp=0,output____T_Xn_tmp_[M_Xt->rpop_j*M_Tt->rpop_j];
  int tab_x=0,mr=0,mx=0,mx_j=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering get_At_T_Xn_ww] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% Calculating output_At_T_Xn_ww\n");}
  switch (output_spacing_t){ case SPACING_j: tab_t_stride = M_Tt->rpop_j; break; case SPACING_b: tab_t_stride = M_Tt->rpop_b; break; case SPACING_a: tab_t_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_w){ case SPACING_j: tab_w_stride = M_Xt->rpop_j; break; case SPACING_b: tab_w_stride = M_Xt->rpop_b; break; case SPACING_a: tab_w_stride = M_Xt->nrows; break; default: break; /* switch (output_spacing_w){ } */}
  output_At_T_Xn_ww->spacing_row = output_spacing_a; output_At_T_Xn_ww->row_stride = tab_a_stride;
  output_At_T_Xn_ww->spacing_col = output_spacing_w; output_At_T_Xn_ww->col_stride = tab_w_stride;
  output_At_T_Xn_ww->spacing_lyr = output_spacing_t; output_At_T_Xn_ww->lyr_stride = tab_t_stride;
  if (strstr(GLOBAL_skip,"At_T_Xn_ww")){ goto skip_At_T_Xn_ww;}
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    ns_j=0;
    while (ns_j<M_Tt->rpop_j){
      ns_a = M_Tt->m_a_[ns_j]; ns_b = M_Tt->m_b_[ns_j];
      switch (output_spacing_t){ case SPACING_j: tab_t=ns_j; break; case SPACING_b: tab_t=ns_b; break; case SPACING_a: tab_t=ns_a; break; default: break; /* switch (output_spacing_t){ } */}
      wAt_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
      wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
      wXn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
      mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
      lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wXn_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
      output____T_Xn_base = dtmp*M_Xn->min_d*M_Xn->mlt_d;
      memset(output____T_Xn_tmp_,0,M_Xt->rpop_j*sizeof(double));
      mw_j=0;
      while (mw_j<M_Xt->rpop_j){
	mw_a = M_Xt->m_a_[mw_j]; mw_b = M_Xt->m_b_[mw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=mw_j; break; case SPACING_b: tab_w=mw_b; break; case SPACING_a: tab_w=mw_a; break; default: break; /* switch (output_spacing_w){ } */}
	output_at_T_Xn_tmp = 0; 
	mr=M_Xn->ncols_per_z*mw_j/* spacing_j */; n2=1;
	for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){
	  wAt_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
	  wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
	  wXn_tag = (__m128i*)&(M_Xn->wX[(mr+mx)*M_Xn->mc_length]);
	  mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	  lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wXn_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
	  output_at_T_Xn_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){ } */}
	output____T_Xn_tmp_[mw_j] = output_at_T_Xn_tmp;
	mw_j++; /* while (mw_j<M_Xt->rpop_j){ } */}
      na_j=0;
      while (na_j<M_At->rpop_j){
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
	wXn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
	mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wXn_tag,&mcat_tag,&mcat_end);
	output_At_T_Xn_base = lltmp*M_Xn->min_d*M_Xn->mlt_d;
	mw_j=0;
	while (mw_j<M_Xt->rpop_j){
	  mw_a = M_Xt->m_a_[mw_j]; mw_b = M_Xt->m_b_[mw_j];
	  switch (output_spacing_w){ case SPACING_j: tab_w=mw_j; break; case SPACING_b: tab_w=mw_b; break; case SPACING_a: tab_w=mw_a; break; default: break; /* switch (output_spacing_w){ } */}
	  tab_x = tab_a + tab_w*tab_a_stride + tab_t*tab_a_stride*tab_w_stride;
	  output_At_T_Xn_tmp = 0;
	  mr=M_Xn->ncols_per_z*mw_j/* spacing_j */; n2=1;
	  for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){
	    wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	    wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
	    wXn_tag = (__m128i*)&(M_Xn->wX[(mr+mx)*M_Xn->mc_length]);
	    mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	    lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wXn_tag,&mcat_tag,&mcat_end);
	    output_At_T_Xn_tmp += n2*lltmp;
	    n2*=2;
	    /* for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){ } */}
	  output_at_T_Xn_base = output____T_Xn_base*a_An[na_a/POPLENGTH];
	  output_at_T_Xn_tmp = output____T_Xn_tmp_[mw_j]*a_An[na_a/POPLENGTH];
	  output_At_T_Xn_ww->lf[tab_x] = sqrt(D_An[na_a/POPLENGTH]) * ( (double)(output_At_T_Xn_base + output_At_T_Xn_tmp)/M_Xn->mlt_d - (double)(output_at_T_Xn_base + output_at_T_Xn_tmp)/M_Xn->mlt_d );
	  mw_j++; /* while (mw_j<M_Xt->rpop_j){ } */}
	na_j++; /* while (na_j<M_At->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>1){ raprintf(output_At_T_Xn_ww->lf,"double",tab_a_stride,tab_t_stride*tab_w_stride," %% output_At_T_Xn_ww->lf: ");}
    GLOBAL_ops_count_one(tidx,0,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Xt->rpop_j*(unsigned long long int)M_Xn->ncols_per_z*(unsigned long long int)M_Xn->mc_length*BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  if (GLOBAL_omp_type==GLOBAL_omp__on){
    memset(output____T_Xn_base_,0,M_Tt->rpop_j*sizeof(double));
    memset(output____T_Xn_tmp_,0,M_Xt->rpop_j*M_Tt->rpop_j*sizeof(double));
    ns_j=0;
    while (ns_j<M_Tt->rpop_j){
      ns_a = M_Tt->m_a_[ns_j]; ns_b = M_Tt->m_b_[ns_j];
      switch (output_spacing_t){ case SPACING_j: tab_t=ns_j; break; case SPACING_b: tab_t=ns_b; break; case SPACING_a: tab_t=ns_a; break; default: break; /* switch (output_spacing_t){ } */}
      wAt_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
      wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
      wXn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
      mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
      lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wXn_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
      output____T_Xn_base_[ns_j] = dtmp*M_Xn->min_d*M_Xn->mlt_d;
      mw_j=0;
      while (mw_j<M_Xt->rpop_j){
	mw_a = M_Xt->m_a_[mw_j]; mw_b = M_Xt->m_b_[mw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=mw_j; break; case SPACING_b: tab_w=mw_b; break; case SPACING_a: tab_w=mw_a; break; default: break; /* switch (output_spacing_w){ } */}
	output_at_T_Xn_tmp = 0; 
	mr=M_Xn->ncols_per_z*mw_j/* spacing_j */; n2=1;
	for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){
	  wAt_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
	  wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
	  wXn_tag = (__m128i*)&(M_Xn->wX[(mr+mx)*M_Xn->mc_length]);
	  mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	  lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wXn_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
	  output_at_T_Xn_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){ } */}
	output____T_Xn_tmp_[mw_j+ns_j*M_Xt->rpop_j] = output_at_T_Xn_tmp;
	mw_j++; /* while (mw_j<M_Xt->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    GLOBAL_ops_count_one(tidx,0,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_Xt->rpop_j*(unsigned long long int)M_Xn->mc_length*BIT8);
#pragma omp parallel private(mx_j,ns_j,ns_a,ns_b,na_j,na_b,na_a,mw_j,mw_b,mw_a,mr,mx,tab_t,tab_a,tab_x,tab_w,wAt_tag,wTt_tag,wXn_tag,mcat_tag,mcat_end,n2,lltmp,dtmp,output____T_Xn_base,output_At_T_Xn_base,output_at_T_Xn_base,output_At_T_Xn_tmp,output_at_T_Xn_tmp)
    { /* begin omp parallel */
      mx_j=0;
#pragma omp for schedule(dynamic)
      for (mx_j=0;mx_j<M_At->rpop_j*M_Tt->rpop_j;mx_j++){
	ns_j = mx_j / M_At->rpop_j; 
	na_j = mx_j % M_At->rpop_j; 
	ns_a = M_Tt->m_a_[ns_j]; ns_b = M_Tt->m_b_[ns_j];
	switch (output_spacing_t){ case SPACING_j: tab_t=ns_j; break; case SPACING_b: tab_t=ns_b; break; case SPACING_a: tab_t=ns_a; break; default: break; /* switch (output_spacing_t){ } */}
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
	wXn_tag = (__m128i*)&(M_Xn->wX[0/* start */*M_Xn->mc_length]);
	mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wXn_tag,&mcat_tag,&mcat_end);
	output_At_T_Xn_base = lltmp*M_Xn->min_d*M_Xn->mlt_d;
	mw_j=0;
	while (mw_j<M_Xt->rpop_j){
	  mw_a = M_Xt->m_a_[mw_j]; mw_b = M_Xt->m_b_[mw_j];
	  switch (output_spacing_w){ case SPACING_j: tab_w=mw_j; break; case SPACING_b: tab_w=mw_b; break; case SPACING_a: tab_w=mw_a; break; default: break; /* switch (output_spacing_w){ } */}
	  tab_x = tab_a + tab_w*tab_a_stride + tab_t*tab_a_stride*tab_w_stride;
	  output_At_T_Xn_tmp = 0;
	  mr=M_Xn->ncols_per_z*mw_j/* spacing_j */; n2=1;
	  for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){
	    wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	    wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
	    wXn_tag = (__m128i*)&(M_Xn->wX[(mr+mx)*M_Xn->mc_length]);
	    mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	    lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wXn_tag,&mcat_tag,&mcat_end);
	    output_At_T_Xn_tmp += n2*lltmp;
	    n2*=2;
	    /* for (mx=M_Xn->ncols_per_z-1;mx>0;mx--){ } */}
	  output_at_T_Xn_base = output____T_Xn_base_[ns_j]*a_An[na_a/POPLENGTH];
	  output_at_T_Xn_tmp = output____T_Xn_tmp_[mw_j+ns_j*M_Xt->rpop_j]*a_An[na_a/POPLENGTH];
	  output_At_T_Xn_ww->lf[tab_x] = sqrt(D_An[na_a/POPLENGTH]) * ( (double)(output_At_T_Xn_base + output_At_T_Xn_tmp)/M_Xn->mlt_d - (double)(output_at_T_Xn_base + output_at_T_Xn_tmp)/M_Xn->mlt_d );
	  mw_j++; /* while (mw_j<M_Xt->rpop_j){ } */}
	/* for (mx_j=0;mx_j<M_At->rpop_j*M_Tt->rpop_j;mx_j++){ } */}
      /* end omp parallel */}
    if (verbose>1){ raprintf(output_At_T_Xn_ww->lf,"double",tab_a_stride,tab_t_stride*tab_w_stride," %% output_At_T_Xn_ww->lf: ");}
    GLOBAL_ops_count_one(tidx,0,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Xt->rpop_j*(unsigned long long int)M_Xn->ncols_per_z*(unsigned long long int)M_Xn->mc_length*BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
 skip_At_T_Xn_ww:
  if (verbose>1){ printf(" %% [finished get_At_T_Xn_ww] tidx %d\n",tidx);}
  return NULL;
}

int wrap_At_T_Xn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_t,int spacing_w,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Xn,struct M_handle *M_Xt,double *A_ajdk,int trm_flag,struct L_handle **output_At_T_Xn_ww_p)
{
  /* calls get_At_T_Xn_ww;
     variable space in **vpra (should be at least size 11)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_t=0,length_w=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_At_T_Xn_ww__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Xt,verbose," %% M_Xt: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_t){ case SPACING_j: length_t = M_Tt->rpop_j; break; case SPACING_b: length_t = M_Tt->rpop_b; break; case SPACING_a: length_t = M_Tt->nrows; break; default: break; /* switch (spacing_t){ } */}
  switch (spacing_w){ case SPACING_j: length_w = M_Xt->rpop_j; break; case SPACING_b: length_w = M_Xt->rpop_b; break; case SPACING_a: length_w = M_Xt->nrows; break; default: break; /* switch (spacing_w){ } */}
  length = length_a*length_t*length_w; if (verbose){ printf(" %% length %llu*%llu*%llu=%llu\n",length_a,length_t,length_w,length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: ");}
  if (verbose>2){ bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: ");}
  if (verbose>2){ bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Xt->mr_b,M_Xt->bitj,1,M_Xt->nrows," %% M_Xt->mr_b: ");}
  if (verbose>2){ bprintf(M_Xt->mr_j,M_Xt->bitj,1,M_Xt->nrows," %% M_Xt->mr_j: ");}
  if (verbose>2){ bprintf(M_Xt->mc_b,M_Xt->bitj,1,M_Xt->ncols," %% M_Xt->mc_b: ");}
  if (verbose>2){ bprintf(M_Xt->mc_j,M_Xt->bitj,1,M_Xt->ncols," %% M_Xt->mc_j: ");}
  if (*output_At_T_Xn_ww_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_At_T_Xn_ww_p = L_handle_make(length);}
  if ((*output_At_T_Xn_ww_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_At_T_Xn_ww__run\n",(*output_At_T_Xn_ww_p)->length,length);}
  memset((*output_At_T_Xn_ww_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Xn; vpra[ip++] = M_Xt; vpra[ip++] = A_ajdk; vpra[ip++] = *output_At_T_Xn_ww_p;
  switch (trm_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break; default: break; /* switch (trm_flag){ } */}
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_t){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_t){ } */}
  switch (spacing_w){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_w){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_At_T_Xn_ww,vpra)){ printf("Warning! cannot create thread %d in wrap_At_T_Xn_ww__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_At_T_Xn_ww(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_At_T_Xn_ww__run] tidx %d\n",*tidx);}
  return length;
}

void *get_At_T_Xn_uu(void *vp)
{
  /* Here M_Xt is only used to store the column mask for lf_Xn */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct L_handle *lf_Xn = (struct L_handle *)(vpra[ip++]);
  struct M_handle *M_Xt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_At_T_Xn_uu = (struct L_handle *)(vpra[ip++]);
  int trm_flag = *(int *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_t = *(int *)(vpra[ip++]);
  int output_spacing_w = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows)/* rup(M_At->nrows + M_At->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int ns_j=0,ns_a=0,ns_b=0,tab_t=0,tab_t_stride=0,na_j=0,na_a=0,na_b=0,tab_a=0,tab_a_stride=0,mw_j=0,mw_a=0,mw_b=0,tab_w=0,tab_w_stride=0,my_j=0,my_a=0,my_b=0,tab_y=0;
  double output_tmp=0;
  int vA=0,vS=0;
  unsigned char *Tt_tag=NULL,*At_tag=NULL;
  __m128i *wAt_tag=NULL,*wTt_tag=NULL,*wXn_tag=NULL,*mcat_tag=NULL,*mcat_end=NULL;
  long long int lltmp=0,n2=0;
  double output_At_T_Xn_base=0,output_at_T_Xn_base=0,output____T_Xn_base=0,dtmp=0,output_At_T_Xn_tmp=0,output_at_T_Xn_tmp=0,output____T_Xn_tmp_[M_Xt->nrows];
  int tab_x=0,mr=0,mx=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering get_At_T_Xn_uu] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% Calculating output_At_T_Xn_uu\n");}
  switch (output_spacing_t){ case SPACING_j: tab_t_stride = M_Tt->rpop_j; break; case SPACING_b: tab_t_stride = M_Tt->rpop_b; break; case SPACING_a: tab_t_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_w){ case SPACING_j: tab_w_stride = M_Xt->rpop_j; break; case SPACING_b: tab_w_stride = M_Xt->rpop_b; break; case SPACING_a: tab_w_stride = M_Xt->nrows; break; default: break; /* switch (output_spacing_w){ } */}
  if (verbose>1){ printf(" %% tab_t_stride %d tab_a_stride %d tab_w_stride %d\n",tab_t_stride,tab_a_stride,tab_w_stride);}
  output_At_T_Xn_uu->spacing_row = output_spacing_a; output_At_T_Xn_uu->row_stride = tab_a_stride;
  output_At_T_Xn_uu->spacing_col = output_spacing_w; output_At_T_Xn_uu->col_stride = tab_w_stride;
  output_At_T_Xn_uu->spacing_lyr = output_spacing_t; output_At_T_Xn_uu->lyr_stride = tab_t_stride;
  if (strstr(GLOBAL_skip,"At_T_Xn_uu")){ goto skip_At_T_Xn_uu;}
  if (verbose>2){ M_handle_printf(M_At,1," %% M_At: ");}
  if (verbose>2){ M_handle_printf(M_Tt,1," %% M_Tt: ");}
  if (verbose>1){ raprintf(lf_Xn->lf,"double",lf_Xn->row_stride,lf_Xn->col_stride," %% lf_Xn->lf: ");}
  if (verbose>1){ raprintf(D_An,"double",1,A_pcols," %% D_An: "); raprintf(a_An,"double",1,A_pcols," %% a_An: ");}
  if (verbose>1){ printf(" %% POPLENGTH %d M_At->mc_length %d M_Tt->mc_length %d \n",POPLENGTH,M_At->mc_length,M_Tt->mc_length);}
  ns_j=0;
  while (ns_j<M_Tt->rpop_j){
    ns_a = M_Tt->m_a_[ns_j]; ns_b = M_Tt->m_b_[ns_j];
    switch (output_spacing_t){ case SPACING_j: tab_t=ns_j; break; case SPACING_b: tab_t=ns_b; break; case SPACING_a: tab_t=ns_a; break; default: break; /* switch (output_spacing_t){ } */}
    Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
    na_j=0;
    while (na_j<M_At->rpop_j){
      na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
      At_tag = (unsigned char *)(&(M_At->wX[na_b*M_At->mc_length]));
      mw_j=0;
      while (mw_j<M_Xt->rpop_j){
	mw_a = M_Xt->m_a_[mw_j]; mw_b = M_Xt->m_b_[mw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=mw_j; break; case SPACING_b: tab_w=mw_b; break; case SPACING_a: tab_w=mw_a; break; default: break; /* switch (output_spacing_w){ } */}
	output_tmp=0;
	my_j=0;
	while (my_j<M_At->cpop_j){
	  my_a = M_At->n_a_[my_j]; my_b = M_At->n_b_[my_j];
	  vA = bget____(At_tag,my_a); vS = bget____(Tt_tag,my_a);
	  output_tmp += (vA - a_An[na_a/POPLENGTH])*vS* (trm_flag ? (*L2_get(lf_Xn,mw_j,mw_b,mw_a,my_j,my_b,my_a)) : (*L2_get(lf_Xn,my_j,my_b,my_a,mw_j,mw_b,mw_a)));
	  my_j++; /* while (my_j<M_At->cpop_j){ } */}
	if (verbose>3){ printf("(%d,%d,%d) %lf ",tab_a,tab_w,tab_t,output_tmp);}
	output_At_T_Xn_uu->lf[tab_a + tab_w*tab_a_stride + tab_t*tab_a_stride*tab_w_stride] = sqrt(D_An[na_a/POPLENGTH]) * output_tmp;
	mw_j++; /* while (mw_j<M_Xt->rpop_j){ } */}
      if (verbose>3){ printf("\n");}
      na_j++; /* while (na_j<M_At->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
  if (verbose>1){ raprintf(output_At_T_Xn_uu->lf,"double",tab_a_stride,tab_t_stride*tab_w_stride," %% output_At_T_Xn_uu->lf: ");}
  GLOBAL_ops_count_one(tidx,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Xt->rpop_j*(unsigned long long int)M_At->cpop_j,0);
 skip_At_T_Xn_uu:
  if (verbose>1){ printf(" %% [finished get_At_T_Xn_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_At_T_Xn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_t,int spacing_w,struct M_handle *M_At,struct M_handle *M_Tt,struct L_handle *lf_Xn,struct M_handle *M_Xt,double *A_ajdk,int trm_flag,struct L_handle **output_At_T_Xn_uu_p)
{
  /* calls get_At_T_Xn_uu;
     variable space in **vpra (should be at least size 11)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_t=0,length_w=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_At_T_Xn_uu__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Xt,verbose," %% M_Xt: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_t){ case SPACING_j: length_t = M_Tt->rpop_j; break; case SPACING_b: length_t = M_Tt->rpop_b; break; case SPACING_a: length_t = M_Tt->nrows; break; default: break; /* switch (spacing_t){ } */}
  switch (spacing_w){ case SPACING_j: length_w = M_Xt->rpop_j; break; case SPACING_b: length_w = M_Xt->rpop_b; break; case SPACING_a: length_w = M_Xt->nrows; break; default: break; /* switch (spacing_w){ } */}
  length = length_a*length_t*length_w; if (verbose){ printf(" %% length %llu*%llu*%llu=%llu\n",length_a,length_t,length_w,length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: ");}
  if (verbose>2){ bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: ");}
  if (verbose>2){ bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Xt->mr_b,M_Xt->bitj,1,M_Xt->nrows," %% M_Xt->mr_b: ");}
  if (verbose>2){ bprintf(M_Xt->mr_j,M_Xt->bitj,1,M_Xt->nrows," %% M_Xt->mr_j: ");}
  if (verbose>2){ bprintf(M_Xt->mc_b,M_Xt->bitj,1,M_Xt->ncols," %% M_Xt->mc_b: ");}
  if (verbose>2){ bprintf(M_Xt->mc_j,M_Xt->bitj,1,M_Xt->ncols," %% M_Xt->mc_j: ");}
  if (*output_At_T_Xn_uu_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_At_T_Xn_uu_p = L_handle_make(length);}
  if ((*output_At_T_Xn_uu_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_At_T_Xn_uu__run\n",(*output_At_T_Xn_uu_p)->length,length);}
  memset((*output_At_T_Xn_uu_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = lf_Xn; vpra[ip++] = M_Xt; vpra[ip++] = A_ajdk; vpra[ip++] = *output_At_T_Xn_uu_p;
  switch (trm_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break; default: break; /* switch (trm_flag){ } */}
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_t){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_t){ } */}
  switch (spacing_w){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_w){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_At_T_Xn_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_At_T_Xn_uu__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_At_T_Xn_uu(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_At_T_Xn_uu__run] tidx %d\n",*tidx);}
  return length;
}

void wrap_At_T_Xn_ww_test()
{
  /* test for errors with input file: At_T_Xn_ww_error.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  unsigned char *bmc1_p=NULL;
  int X_n_cols=0,nr_j=0,nr_b=0,nr_a=0,nc_j=0,nc_b=0,nc_a=0;
  struct L_handle **lf_Xn=NULL;
  struct M_handle **M_Xn=NULL,**M_Xt=NULL;
  struct L_handle **lf_At_T_Xn_ww=NULL,**lf_At_T_AnZt_ww=NULL; int *length_At_T_Xn_ww=NULL;
  struct L_handle **lf_At_T_Xn_uu=NULL,**lf_At_T_AnZt_uu=NULL; int *length_At_T_Xn_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_At_T_Xn_ww_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  X_n_cols = 8;
  bmc1_p = wkspace_all0c(bsize(X_n_cols)); 
  bset_off(bmc1_p,0); bset__on(bmc1_p,1); bset_off(bmc1_p,2); bset__on(bmc1_p,3); bset_off(bmc1_p,4); bset_off(bmc1_p,5); bset__on(bmc1_p,6); bset__on(bmc1_p,7); 
  M_Xt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_Xt[nb] = M_handle_v_make(BITJ,X_n_cols,A_n_rows[nb],NULL,NULL,bmc1_p,M_An[nb]->mr_b);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_Xn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_Xn[nb] = L_handle_make((unsigned long long int)X_n_cols*(unsigned long long int)A_n_rows[nb]);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_Xn = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_Xn[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,A_n_rows[nb],X_n_cols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_At_T_Xn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_Xn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_At_T_Xn_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)X_n_cols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_At_T_Xn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_Xn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_At_T_Xn_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)X_n_cols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
  	if (verbose){ printf(" %% %s; %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_An[nb],A_ajdk,&(lf_An_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	for (nb=0;nb<nbins;nb++){
	  M_mxset(M_Xt[nb],bmc1_p,M_An[nb]->mr_j);
	  if (verbose>2){ M_handle_printf(M_Xt[nb],1," %% M_Xt[nb]: ");}
	  lf_Xn[nb]->spacing_row = n_spacing_A; 
	  switch(lf_Xn[nb]->spacing_row){ case SPACING_j: lf_Xn[nb]->row_stride = M_Xt[nb]->cpop_j; break; case SPACING_b: lf_Xn[nb]->row_stride = M_Xt[nb]->cpop_b; break; case SPACING_a: lf_Xn[nb]->row_stride = M_Xt[nb]->ncols; break; default: break; /* switch(lf_Xn[nb]->spacing_row){ } */}
	  lf_Xn[nb]->spacing_col = n_spacing_A;
	  switch(lf_Xn[nb]->spacing_col){ case SPACING_j: lf_Xn[nb]->col_stride = M_Xt[nb]->rpop_j; break; case SPACING_b: lf_Xn[nb]->col_stride = M_Xt[nb]->rpop_b; break; case SPACING_a: lf_Xn[nb]->col_stride = M_Xt[nb]->nrows; break; default: break; /* switch(lf_Xn[nb]->spacing_col){ } */}
	  lf_Xn[nb]->spacing_lyr = n_spacing_A;
	  switch(lf_Xn[nb]->spacing_lyr){ case SPACING_j: lf_Xn[nb]->lyr_stride = M_Tt[nb]->rpop_j; break; case SPACING_b: lf_Xn[nb]->lyr_stride = M_Tt[nb]->rpop_b; break; case SPACING_a: lf_Xn[nb]->lyr_stride = M_Tt[nb]->nrows; break; default: break; /* switch(lf_Xn[nb]->spacing_lyr){ } */}
	  for (nc_j=0;nc_j<M_Xt[nb]->rpop_j;nc_j++){ 
	    nc_a = M_Xt[nb]->m_a_[nc_j]; nc_b = M_Xt[nb]->m_b_[nc_j];
	    for (nr_j=0;nr_j<M_Xt[nb]->cpop_j;nr_j++){
	      nr_a = M_Xt[nb]->n_a_[nr_j]; nr_b = M_Xt[nb]->n_b_[nr_j];
	      if (verbose>2){ printf(" %% nb: %d --> nr (%d,%d,%d)/%d , nc (%d,%d,%d)/%d\n",nb,nr_j,nr_b,nr_a,M_Xt[nb]->cpop_j,nc_j,nc_b,nc_a,M_Xt[nb]->rpop_j);}
	      L2_set(lf_Xn[nb],nr_j,nr_b,nr_a,nc_j,nc_b,nc_a,rand01);
	      /* for (nr_j=0;nr_j<M_Xt[nb]->cpop_j;nr_j++){ } */}
	    /* for (nc_j=0;nc_j<M_Xt[nb]->rpop_j;nc_j++){ } */}
	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
	  wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Xt[nb]->mc_j,M_Xt[nb]->mc_b,M_Xt[nb]->mr_j,M_Xt[nb]->mr_b,lf_Xn[nb],lf_Xn[nb]->lf,&(M_Xt[nb]->ncols),&(M_Xt[nb]->nrows),&(M_Xn[nb]),&(GLOBAL_B_MLT),(addressable_0));
	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% xcalc trm_0: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
	  length_At_T_Xn_ww[nb] = wrap_At_T_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_Xn[nb],M_Xt[nb],A_ajdk,(addressable_0),&(lf_At_T_Xn_ww[nb]));
	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% At_T_Xn_ww : ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic(); 
	    length_At_T_Xn_uu[nb] = wrap_At_T_Xn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],lf_Xn[nb],M_Xt[nb],A_ajdk,(addressable_0),&(lf_At_T_Xn_uu[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% At_T_Xn_uu : ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_At_T_Xn_ww[%d] error %0.16f\n",nb,dra_diff(lf_At_T_Xn_ww[nb]->lf,lf_At_T_Xn_uu[nb]->lf,length_At_T_Xn_ww[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_At_T_Xn_ww_test]\n");}
}
