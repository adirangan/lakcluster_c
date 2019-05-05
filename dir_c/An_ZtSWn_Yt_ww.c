#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_An_ZtSWn_Yt_ww(void *vp)
{
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle **M_ZtSWn_ = (struct M_handle **)(vpra[ip++]); struct M_handle *M_ZtSWn=NULL;
  struct M_handle *M_Wt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Yn = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_An_ZtSWn_Yt_ww = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); int output_spacing_y = output_spacing_a;
  int output_spacing_s = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yn->ncols);
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int ns_j=0,ns_a=0,ns_b=0,tab_s=0,tab_s_stride=0,ma_j=0,ma_a=0,ma_b=0,tab_a=0,tab_a_stride=0,nz_j=0,nz_a=0,nz_b=0,tab_z=0,nw_j=0,nw_a=0,nw_b=0,my_j=0,my_a=0,my_b=0,tab_y=0;
  double output_tmp=0;
  int vY=0;
  unsigned char *An_tag=NULL,*Yn_tag=NULL;
  __m128i *wAn_tag=NULL,*wYn_tag=NULL,*wZtSWn_tag=NULL,*mcan_tag=NULL,*mcan_end=NULL;
  long long int n2=0;
  double output_An_ZtSWn_base=0,output_an_ZtSWn_base=0,output____ZtSWn_base=0,*dinp=NULL,dtmp=0,output_An_ZtSWn_tmp=0,output_an_ZtSWn_tmp=0,output____ZtSWn_tmp_[M_Wt->nrows],output_An_ZtSWn_Yt_tmp=0;
  int tab_x=0,mr=0,mx=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering get_An_ZtSWn_Yt_ww] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% Calculating output_An_ZtSWn_Yt_ww\n");}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  output_An_ZtSWn_Yt_ww->spacing_row = output_spacing_a; output_An_ZtSWn_Yt_ww->row_stride = tab_a_stride;
  output_An_ZtSWn_Yt_ww->spacing_col = output_spacing_s; output_An_ZtSWn_Yt_ww->col_stride = tab_s_stride;
  if (strstr(GLOBAL_skip,"An_ZtSWn_Yt_ww")){ goto skip_An_ZtSWn_Yt_ww;}
  ns_j=0;
  while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j]; M_ZtSWn = M_ZtSWn_[ns_j];
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    wAn_tag = (__m128i*)&(M_ZtSWn->wX[0/* start */*M_ZtSWn->mc_length]);
    wZtSWn_tag = (__m128i*)&(M_ZtSWn->wX[0/* start */*M_ZtSWn->mc_length]);
    mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
    dinp = &(A_ajdk[0+AJDK_1_1*A_pcols]);
    dtmp = popcount_pm0_lf(&wAn_tag,&wZtSWn_tag,&mcan_tag,&mcan_end,&dinp);
    output____ZtSWn_base = dtmp*M_ZtSWn->min_d*M_ZtSWn->mlt_d;
    memset(output____ZtSWn_tmp_,0,M_Wt->rpop_j*sizeof(double));
    nw_j=0;
    while (nw_j<M_Wt->rpop_j){
      nw_a = M_Wt->m_a_[nw_j]; nw_b = M_Wt->m_b_[nw_j];
      output_an_ZtSWn_tmp = 0; 
      mr=M_ZtSWn->ncols_per_z*nw_j/* spacing_j */; n2=1;
      for (mx=M_ZtSWn->ncols_per_z-1;mx>0;mx--){
	wAn_tag = (__m128i*)&(M_ZtSWn->wX[0/* start */*M_ZtSWn->mc_length]);
	wZtSWn_tag = (__m128i*)&(M_ZtSWn->wX[(mr+mx)*M_ZtSWn->mc_length]);
	mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	dinp = &(A_ajdk[0+AJDK_1_1*A_pcols]);
	dtmp = popcount_pm0_lf(&wAn_tag,&wZtSWn_tag,&mcan_tag,&mcan_end,&dinp);
	output_an_ZtSWn_tmp += n2*dtmp;
	n2*=2;
	/* for (mx=M_ZtSWn->ncols_per_z-1;mx>0;mx--){ } */}
      output____ZtSWn_tmp_[nw_j] = output_an_ZtSWn_tmp;
      nw_j++; /* while (nw_j<M_Wt->rpop_j){ } */}
    ma_j=0; my_j=0;
    while (ma_j<M_An->rpop_j && my_j<M_Yn->rpop_j){
      ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j]; my_a = M_Yn->m_a_[my_j]; my_b = M_Yn->m_b_[my_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=ma_j; break; case SPACING_b: tab_a=ma_b; break; case SPACING_a: tab_a=ma_a; break; default: break; /* switch (output_spacing_a){ } */}
      switch (output_spacing_y){ case SPACING_j: tab_y=my_j; break; case SPACING_b: tab_y=my_b; break; case SPACING_a: tab_y=my_a; break; default: break; /* switch (output_spacing_y){ } */}
      tab_x = tab_a + tab_s*tab_a_stride;
      output_An_ZtSWn_Yt_ww->lf[tab_x] = 0;
      wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
      wZtSWn_tag = (__m128i*)&(M_ZtSWn->wX[0/* start */*M_ZtSWn->mc_length]);
      mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
      dinp = &(A_ajdk[0+AJDK_0_1*A_pcols]);
      dtmp = popcount_pm0_lf(&wAn_tag,&wZtSWn_tag,&mcan_tag,&mcan_end,&dinp);
      output_An_ZtSWn_base = dtmp*M_ZtSWn->min_d*M_ZtSWn->mlt_d;
      wYn_tag = (__m128i*)&(M_Yn->wX[my_b*M_Yn->mc_length]);
      nw_j=0;
      while (nw_j<M_Wt->rpop_j){
	nw_a = M_Wt->m_a_[nw_j]; nw_b = M_Wt->m_b_[nw_j];
	output_An_ZtSWn_tmp = 0;
	mr=M_ZtSWn->ncols_per_z*nw_j/* spacing_j */; n2=1;
	for (mx=M_ZtSWn->ncols_per_z-1;mx>0;mx--){
	  wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	  wZtSWn_tag = (__m128i*)&(M_ZtSWn->wX[(mr+mx)*M_ZtSWn->mc_length]);
	  mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	  dinp = &(A_ajdk[0+AJDK_0_1*A_pcols]);
	  dtmp = popcount_pm0_lf(&wAn_tag,&wZtSWn_tag,&mcan_tag,&mcan_end,&dinp);
	  output_An_ZtSWn_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_ZtSWn->ncols_per_z-1;mx>0;mx--){ } */}
	output_an_ZtSWn_base = output____ZtSWn_base; output_an_ZtSWn_tmp = output____ZtSWn_tmp_[nw_j];
	output_An_ZtSWn_tmp = (double)(output_An_ZtSWn_base + output_An_ZtSWn_tmp)/M_ZtSWn->mlt_d - (double)(output_an_ZtSWn_base + output_an_ZtSWn_tmp)/M_ZtSWn->mlt_d;
	Yn_tag = (unsigned char *) wYn_tag; vY = bget____(Yn_tag,nw_a);
	output_An_ZtSWn_Yt_tmp = output_An_ZtSWn_tmp * D_Yn[nw_a/POPLENGTH]*(vY - a_Yn[nw_a/POPLENGTH]);
	output_An_ZtSWn_Yt_ww->lf[tab_x] += output_An_ZtSWn_Yt_tmp;
	nw_j++; /* while (nw_j<M_Wt->rpop_j){ } */}
      ma_j++; my_j++; /* while (ma_j<M_An->rpop_j && my_j<M_Yn->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  if (verbose>1){ raprintf(output_An_ZtSWn_Yt_ww->lf,"double",tab_a_stride,tab_s_stride," %% output_An_ZtSWn_Yt_ww->lf: ");}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_An->rpop_j*M_Wt->rpop_j,M_St->rpop_j*M_An->rpop_j*M_Wt->rpop_j*M_ZtSWn->ncols_per_z*M_ZtSWn->mc_length*BIT8);
 skip_An_ZtSWn_Yt_ww:
  if (verbose>1){ printf(" %% [finished get_An_ZtSWn_Yt_ww] tidx %d\n",tidx);}
  return NULL;
}

int wrap_An_ZtSWn_Yt_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,struct M_handle *M_An,struct M_handle *M_St,struct M_handle **M_ZtSWn_,struct M_handle *M_Wt,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle **output_An_ZtSWn_Yt_ww_p)
{
  /* calls get_An_ZtSWn_Yt_ww;
     variable space in **vpra (should be at least size 13)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_An_ZtSWn_Yt_ww__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Wt,verbose," %% M_Wt: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (spacing_s){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_a,length_s,length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Wt->mr_b,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_b: ");}
  if (verbose>2){ bprintf(M_Wt->mr_j,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_j: ");}
  if (verbose>2){ bprintf(M_Wt->mc_b,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_b: ");}
  if (verbose>2){ bprintf(M_Wt->mc_j,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_j: ");}
  if (*output_An_ZtSWn_Yt_ww_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_ZtSWn_Yt_ww_p = L_handle_make(length);}
  if ((*output_An_ZtSWn_Yt_ww_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_ZtSWn_Yt_ww__run\n",(*output_An_ZtSWn_Yt_ww_p)->length,length);}
  memset((*output_An_ZtSWn_Yt_ww_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_St; vpra[ip++] = M_ZtSWn_; vpra[ip++] = M_Wt; vpra[ip++] = M_Yn; vpra[ip++] = A_ajdk; vpra[ip++] = Y_ajdk; vpra[ip++] = *output_An_ZtSWn_Yt_ww_p;
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_ZtSWn_Yt_ww,vpra)){ printf("Warning! cannot create thread %d in wrap_An_ZtSWn_Yt_ww__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_ZtSWn_Yt_ww(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_An_ZtSWn_Yt_ww__run] tidx %d\n",*tidx);}
  return length;
}

void wrap_An_ZtSWn_Yt_ww_test()
{
  /* test for errors with input file: An_ZtSWn_Yt_ww_error.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AnZt=NULL,**lf_AnAt=NULL,**lf_YnWt=NULL,**lf_YnYt=NULL;
  struct L_handle **lf_AtTn=NULL,**lf_YtTn=NULL,**lf_ZtSn=NULL,**lf_WtSn=NULL;
  struct L_handle **lf_ZtSWn=NULL,**lf_AtTYn=NULL;
  int ns_j=0,ns_b=0,ns_a=0;
  struct M_handle ***M_ZtSWn_=NULL,***M_AtTYn_=NULL;
  struct L_handle **lf_An_ZtSWn_Yt_ww=NULL,**lf_An_AtTYn_Yt_ww=NULL; int *length_An_ZtSWn_Yt_ww=NULL,*length_An_AtTYn_Yt_ww=NULL;
  struct M_handle ***M_WtSZn_=NULL,***M_YtTAn_=NULL;
  struct L_handle **lf_Yn_WtSZn_At_ww=NULL,**lf_Yn_YtTAn_At_ww=NULL; int *length_Yn_WtSZn_At_ww=NULL,*length_Yn_YtTAn_At_ww=NULL;
  struct L_handle **lf_An_ZtSWn_Yt_vv=NULL,**lf_An_AtTYn_Yt_vv=NULL; int *length_An_ZtSWn_Yt_vv=NULL,*length_An_AtTYn_Yt_vv=NULL;
  struct L_handle **lf_AnZt_S_WnYt_vv=NULL,**lf_AnAt_T_YnYt_vv=NULL; int *length_AnZt_S_WnYt_vv=NULL,*length_AnAt_T_YnYt_vv=NULL;
  if (verbose){ printf(" %% [entering wrap_An_ZtSWn_Yt_ww_test]\n");}
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
  lf_AnZt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AnAt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_YnWt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_YnYt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_AnZt[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    lf_AnAt[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_An[nb]->nrows);
    lf_YnWt[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_Wn[nb]->nrows);
    lf_YnYt[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_Yn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_ZtSWn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AtTYn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_ZtSWn[nb] = L_handle_make((unsigned long long int)M_Zt[nb]->nrows*(unsigned long long int)M_Wt[nb]->nrows*(unsigned long long int)M_St[nb]->nrows);
    lf_AtTYn[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Yt[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_ZtSWn_ = (struct M_handle ***)wkspace_all0c(sizeof(struct M_handle **)*nbins);
  M_AtTYn_ = (struct M_handle ***)wkspace_all0c(sizeof(struct M_handle **)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_ZtSWn_[nb] = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*M_St[nb]->rpop_j);
    M_AtTYn_[nb] = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*M_Tt[nb]->rpop_j);
    for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){
      M_ZtSWn_[nb][ns_j] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Zt[nb]->nrows,M_Wt[nb]->nrows);
      M_AtTYn_[nb][ns_j] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_At[nb]->nrows,M_Yt[nb]->nrows);
      /* for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){ } */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_An_ZtSWn_Yt_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_ZtSWn_Yt_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_An_AtTYn_Yt_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_AtTYn_Yt_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_An_ZtSWn_Yt_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_St[nb]->nrows);
    lf_An_AtTYn_Yt_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_WtSZn_ = (struct M_handle ***)wkspace_all0c(sizeof(struct M_handle **)*nbins);
  M_YtTAn_ = (struct M_handle ***)wkspace_all0c(sizeof(struct M_handle **)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_WtSZn_[nb] = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*M_St[nb]->rpop_j);
    M_YtTAn_[nb] = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*M_Tt[nb]->rpop_j);
    for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){
      M_WtSZn_[nb][ns_j] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Wt[nb]->nrows,M_Zt[nb]->nrows);
      M_YtTAn_[nb][ns_j] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Yt[nb]->nrows,M_At[nb]->nrows);
      /* for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){ } */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_Yn_WtSZn_At_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Yn_WtSZn_At_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Yn_YtTAn_At_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Yn_YtTAn_At_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_Yn_WtSZn_At_ww[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_St[nb]->nrows);
    lf_Yn_YtTAn_At_ww[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_AnAt_T_YnYt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnAt_T_YnYt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_AnZt_S_WnYt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnZt_S_WnYt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_AnAt_T_YnYt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
    lf_AnZt_S_WnYt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_An_ZtSWn_Yt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_ZtSWn_Yt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_An_AtTYn_Yt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_AtTYn_Yt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_An_ZtSWn_Yt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
    lf_An_AtTYn_Yt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
  	if (verbose){ printf(" %% %s; %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_An[nb],A_ajdk,&(lf_An_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Zn[nb],A_ajdk,&(lf_Zn_ajdk[nb]));
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
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
  	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_An[nb],M_Zn[nb],A_ajdk,lf_An_ajdk[nb],lf_Zn_ajdk[nb],&(lf_AnZt[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
  	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_An[nb],M_An[nb],A_ajdk,lf_An_ajdk[nb],lf_An_ajdk[nb],&(lf_AnAt[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
  	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Yn[nb],M_Wn[nb],Y_ajdk,lf_Yn_ajdk[nb],lf_Wn_ajdk[nb],&(lf_YnWt[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
  	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Yn[nb],M_Yn[nb],Y_ajdk,lf_Yn_ajdk[nb],lf_Yn_ajdk[nb],&(lf_YnYt[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZt and YnWt: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
	  length_AnZt_S_WnYt_vv[nb] = wrap_AnZt_S_WnYt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_St[nb],M_Zn[nb],lf_AnZt[nb],lf_YnWt[nb],&(lf_AnZt_S_WnYt_vv[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
	  length_AnAt_T_YnYt_vv[nb] = wrap_AnZt_S_WnYt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_Tt[nb],M_An[nb],lf_AnAt[nb],lf_YnYt[nb],&(lf_AnAt_T_YnYt_vv[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZt_S_WnYt_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
	  wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],M_Yt[nb],A_ajdk,Y_ajdk,lf_AtTn[nb],lf_YtTn[nb],&(lf_AtTYn[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],M_Wt[nb],A_ajdk,Y_ajdk,lf_ZtSn[nb],lf_WtSn[nb],&(lf_ZtSWn[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYn and ZtSWn: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
	  for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){
	    ns_a = M_St[nb]->m_a_[ns_j]; ns_b = M_St[nb]->m_b_[ns_j];
	    GLOBAL_pthread_tic(); 
	    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Zt[nb]->mr_j,M_Zt[nb]->mr_b,M_Wt[nb]->mr_j,M_Wt[nb]->mr_b,lf_ZtSWn[nb],L3_lf_get(lf_ZtSWn[nb],lf_ZtSWn[nb]->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_Zt[nb]->nrows),&(M_Wt[nb]->nrows),&(M_ZtSWn_[nb][ns_j]),&(GLOBAL_B_MLT),(addressable_0));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Wt[nb]->mr_j,M_Wt[nb]->mr_b,M_Zt[nb]->mr_j,M_Zt[nb]->mr_b,lf_ZtSWn[nb],L3_lf_get(lf_ZtSWn[nb],lf_ZtSWn[nb]->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_Wt[nb]->nrows),&(M_Zt[nb]->nrows),&(M_WtSZn_[nb][ns_j]),&(GLOBAL_B_MLT),(addressable_1));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_At[nb]->mr_j,M_At[nb]->mr_b,M_Yt[nb]->mr_j,M_Yt[nb]->mr_b,lf_AtTYn[nb],L3_lf_get(lf_AtTYn[nb],lf_AtTYn[nb]->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_At[nb]->nrows),&(M_Yt[nb]->nrows),&(M_AtTYn_[nb][ns_j]),&(GLOBAL_B_MLT),(addressable_0));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Yt[nb]->mr_j,M_Yt[nb]->mr_b,M_At[nb]->mr_j,M_At[nb]->mr_b,lf_AtTYn[nb],L3_lf_get(lf_AtTYn[nb],lf_AtTYn[nb]->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_Yt[nb]->nrows),&(M_At[nb]->nrows),&(M_YtTAn_[nb][ns_j]),&(GLOBAL_B_MLT),(addressable_1));
	    GLOBAL_pthread_toc();
	    /* for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){ } */}
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% xcalc: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
	  length_An_ZtSWn_Yt_ww[nb] = wrap_An_ZtSWn_Yt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_St[nb],M_ZtSWn_[nb],M_Wt[nb],M_Yn[nb],A_ajdk,Y_ajdk,&(lf_An_ZtSWn_Yt_ww[nb]));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  length_Yn_WtSZn_At_ww[nb] = wrap_An_ZtSWn_Yt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_Yn[nb],M_St[nb],M_WtSZn_[nb],M_Zt[nb],M_An[nb],Y_ajdk,A_ajdk,&(lf_Yn_WtSZn_At_ww[nb]));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  length_An_AtTYn_Yt_ww[nb] = wrap_An_ZtSWn_Yt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_Tt[nb],M_AtTYn_[nb],M_Yt[nb],M_Yn[nb],A_ajdk,Y_ajdk,&(lf_An_AtTYn_Yt_ww[nb]));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  length_Yn_YtTAn_At_ww[nb] = wrap_An_ZtSWn_Yt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_Yn[nb],M_Tt[nb],M_YtTAn_[nb],M_At[nb],M_An[nb],Y_ajdk,A_ajdk,&(lf_Yn_YtTAn_At_ww[nb]));
	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_ZtSWn_Yt_ww Yn_WtSZn_At_ww: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic();
	  length_An_ZtSWn_Yt_vv[nb] = wrap_An_ZtSWn_Yt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_St[nb],M_Yn[nb],A_ajdk,Y_ajdk,lf_ZtSWn[nb],&(lf_An_ZtSWn_Yt_vv[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  length_An_AtTYn_Yt_vv[nb] = wrap_An_ZtSWn_Yt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_Tt[nb],M_Yn[nb],A_ajdk,Y_ajdk,lf_AtTYn[nb],&(lf_An_AtTYn_Yt_vv[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_ZtSWn_Yt_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_An_ZtSWn_Yt_ww[%d] - lf_An_ZtSWn_Yt_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_An_ZtSWn_Yt_ww[nb]->lf,lf_An_ZtSWn_Yt_vv[nb]->lf,length_An_ZtSWn_Yt_ww[nb],1));
	    printf(" %% lf_An_ZtSWn_Yt_ww[%d] - lf_Yn_WtSZn_At_ww[%d] %0.16f\n",nb,nb,dra_diff(lf_An_ZtSWn_Yt_ww[nb]->lf,lf_Yn_WtSZn_At_ww[nb]->lf,length_An_ZtSWn_Yt_ww[nb],1));
	    printf(" %% lf_An_ZtSWn_Yt_ww[%d] - lf_AnZt_S_WnYt_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_An_ZtSWn_Yt_ww[nb]->lf,lf_AnZt_S_WnYt_vv[nb]->lf,length_An_ZtSWn_Yt_ww[nb],1));
	    printf(" %% lf_An_ZtSWn_Yt_vv[%d] - lf_AnZt_S_WnYt_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_An_ZtSWn_Yt_vv[nb]->lf,lf_AnZt_S_WnYt_vv[nb]->lf,length_An_ZtSWn_Yt_ww[nb],1));
	    printf(" %% lf_An_AtTYn_Yt_ww[%d] - lf_An_AtTYn_Yt_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_An_AtTYn_Yt_ww[nb]->lf,lf_An_AtTYn_Yt_vv[nb]->lf,length_An_AtTYn_Yt_ww[nb],1));
	    printf(" %% lf_An_AtTYn_Yt_ww[%d] - lf_Yn_YtTAn_At_ww[%d] %0.16f\n",nb,nb,dra_diff(lf_An_AtTYn_Yt_ww[nb]->lf,lf_Yn_YtTAn_At_ww[nb]->lf,length_An_AtTYn_Yt_ww[nb],1));
	    printf(" %% lf_An_AtTYn_Yt_ww[%d] - lf_AnAt_T_YnYt_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_An_AtTYn_Yt_ww[nb]->lf,lf_AnAt_T_YnYt_vv[nb]->lf,length_An_AtTYn_Yt_ww[nb],1));
	    printf(" %% lf_An_AtTYn_Yt_vv[%d] - lf_AnAt_T_YnYt_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_An_AtTYn_Yt_vv[nb]->lf,lf_AnAt_T_YnYt_vv[nb]->lf,length_An_AtTYn_Yt_ww[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_An_ZtSWn_Yt_ww_test]\n");}
}
