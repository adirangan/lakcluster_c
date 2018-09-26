#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_AnZt_S_WnYt_vv(void *vp)
{
  /* This function takes in M_An, M_St, M_Zn, lf_AnZt, lf_YnWt and calculates output_AnZt_S_WnYt[mj+ns*A_n_rows] = AnZt(mj,:)*diag(S)*WnYt(:,mj); ;
     in this scenario, lf_YnWt = (Yn(:,:)-e_An*a_Yt)*diag(D_Yn)*(Wt(:,:)-a_Yn*e_Zt) ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zn = (struct M_handle *)(vpra[ip++]);
  struct L_handle *lf_AnZt = (struct L_handle *)(vpra[ip++]); /* used for direct calculation */
  struct L_handle *lf_YnWt = (struct L_handle *)(vpra[ip++]); /* used for direct calculation */
  struct L_handle *output_AnZt_S_WnYt = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); int output_spacing_An = output_spacing_a; int output_spacing_Zn = output_spacing_a;
  int output_spacing_s = *(int *)(vpra[ip++]); int output_spacing_Tt = output_spacing_s;
  int mx_j=0,mx_chunk=0,ns_j=0,ns_b=0,ns_a=0,tab_Tt_stride=0,tab_St=0/* ,na_j=0,na_b=0,na_a=0 */,ma_j=0,ma_b=0,ma_a=0,tab_An_stride=0,tab_An=0/* ,nz_j=0,nz_b=0,nz_a=0 */,mz_j=0,mz_b=0,mz_a=0,tab_Zn=0;
  double output_AnZt_S_WnYt_tmp=0; 
  unsigned int *ma_b_,*ma_a_;
  unsigned int *na_b_,*na_a_;
  unsigned int *mz_b_,*mz_a_;
  unsigned int *nz_b_,*nz_a_;
  unsigned int *ns_b_,*ns_a_;
  int vM=0,vS=0;
  unsigned char *St_tag=NULL,*mcs_tag=NULL,*mcs_end=NULL;
  if (verbose>1){ printf(" %% [entering get_AnZt_S_WnYt_vv] tidx %d \n",tidx);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Zn->mr_b,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_b: "); bprintf(M_Zn->mr_j,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_j: ");}
  if (verbose>2){ bprintf(M_Zn->mc_b,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_b: "); bprintf(M_Zn->mc_j,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_j: ");}
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  na_b_ = M_An->n_b_; na_a_ = M_An->n_a_;
  mz_b_ = M_Zn->m_b_; mz_a_ = M_Zn->m_a_;
  nz_b_ = M_Zn->n_b_; nz_a_ = M_Zn->n_a_;
  ns_b_ = M_St->m_b_; ns_a_ = M_St->m_a_;
  switch (output_spacing_Tt){ case SPACING_j: tab_Tt_stride = M_St->rpop_j; break; case SPACING_b: tab_Tt_stride = M_St->rpop_b; break; case SPACING_a: tab_Tt_stride = M_St->nrows; break; default: break; /* switch (output_spacing_Tt){ } */}
  switch (output_spacing_An){ case SPACING_j: tab_An_stride = M_An->rpop_j; break; case SPACING_b: tab_An_stride = M_An->rpop_b; break; case SPACING_a: tab_An_stride = M_An->nrows; break; default: break; /* switch (output_spacing_An){ } */}
  output_AnZt_S_WnYt->spacing_row = output_spacing_An; output_AnZt_S_WnYt->row_stride = tab_An_stride;
  output_AnZt_S_WnYt->spacing_col = output_spacing_Tt; output_AnZt_S_WnYt->col_stride = tab_Tt_stride;
  if (strstr(GLOBAL_skip,"AnZt_S_WnYt_vv")){ goto skip_AnZt_S_WnYt_vv;}
  if (lf_AnZt!=NULL && lf_YnWt!=NULL){
    if (verbose>2){ raprintf(lf_AnZt->lf,"double",M_An->rpop_b,M_Zn->rpop_b," %% lf_AnZt->lf: "); raprintf(lf_YnWt->lf,"double",M_An->rpop_b,M_Zn->rpop_b," %% lf_YnWt->lf: ");}
    if (GLOBAL_omp_type==GLOBAL_omp_off){ 
      ns_j=0; while (ns_j<M_St->rpop_j){
	ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	St_tag = (unsigned char *)(&(M_St->wX[ns_b*M_St->mc_length])); if (verbose>2){ printf(" %% ns_j %d ns_a %d ns_b %d \n",ns_j,ns_a,ns_b); bprintf((unsigned char *)St_tag,M_St->bitj,1,M_St->ncols," % % wSt_tag: ");}
	mcs_tag = (unsigned char *)(&(M_St->mc_j[0]));
	mcs_end = (unsigned char *)(&(M_St->mc_j[M_St->mc_length]));
	switch (output_spacing_Tt){ case SPACING_j: tab_St=ns_j; break; case SPACING_b: tab_St=ns_b; break; case SPACING_a: tab_St=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
	ma_j=0; while (ma_j<M_An->rpop_j){
	  ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
	  switch (output_spacing_An){ case SPACING_j: tab_An=ma_j; break; case SPACING_b: tab_An=ma_b; break; case SPACING_a: tab_An=ma_a; break; default: break; /* switch (output_spacing_An){ } */}
	  output_AnZt_S_WnYt_tmp=0;
	  mz_j=0; while (mz_j<M_Zn->rpop_j){
	    mz_a = M_Zn->m_a_[mz_j]; mz_b = M_Zn->m_b_[mz_j];
	    switch (output_spacing_Zn){ case SPACING_j: tab_Zn=mz_j; break; case SPACING_b: tab_Zn=mz_b; break; case SPACING_a: tab_Zn=mz_a; break; default: break; /* switch (output_spacing_Zn){ } */}
	    vS = bget____(St_tag,mz_a); vM = bget__on(mcs_tag,mz_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_St->mc_j\n",mz_j,mz_b,mz_a);}
	    output_AnZt_S_WnYt_tmp += (*L2_get(lf_AnZt,ma_j,ma_b,ma_a,mz_j,mz_b,mz_a)) * vS * (*L2_get(lf_YnWt,ma_j,ma_b,ma_a,mz_j,mz_b,mz_a));
	    mz_j++; /* while (mz_j<M_Zn->rpop_j){ } */}
	  output_AnZt_S_WnYt->lf[tab_An + tab_St*tab_An_stride] = output_AnZt_S_WnYt_tmp;
	  ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
	ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
      /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
    else if (GLOBAL_omp_type==GLOBAL_omp__on){ /* useful when T_n_cols is much less than number of processors (e.g., 1 to 6) */
      mx_chunk=128; 
#pragma omp parallel shared(mx_chunk) private(mx_j,ns_j,ns_a,ns_b,ma_j,ma_a,ma_b,mz_j,mz_a,mz_b,St_tag,mcs_tag,mcs_end,tab_St,tab_An,vS,vM,output_AnZt_S_WnYt_tmp)
      { /* begin omp parallel */
	mx_j=0; tab_St=0; tab_An=0; 
#pragma omp for schedule(dynamic,mx_chunk)
	for (mx_j=0;mx_j<M_St->rpop_j*M_An->rpop_j;mx_j++){
	  ns_j = mx_j / (M_An->rpop_j); ma_j = mx_j % (M_An->rpop_j);
	  ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	  St_tag = (unsigned char *)(&(M_St->wX[ns_b*M_St->mc_length]));
	  mcs_tag = (unsigned char *)(&(M_St->mc_j[0]));
	  mcs_end = (unsigned char *)(&(M_St->mc_j[M_St->mc_length]));
	  switch (output_spacing_Tt){ case SPACING_j: tab_St=ns_j; break; case SPACING_b: tab_St=ns_b; break; case SPACING_a: tab_St=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
	  ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
	  switch (output_spacing_An){ case SPACING_j: tab_An=ma_j; break; case SPACING_b: tab_An=ma_b; break; case SPACING_a: tab_An=ma_a; break; default: break; /* switch (output_spacing_An){ } */}
	  output_AnZt_S_WnYt_tmp=0;
	  mz_j=0; while (mz_j<M_Zn->rpop_j){
	    mz_a = M_Zn->m_a_[mz_j]; mz_b = M_Zn->m_b_[mz_j];
	    switch (output_spacing_Zn){ case SPACING_j: tab_Zn=mz_j; break; case SPACING_b: tab_Zn=mz_b; break; case SPACING_a: tab_Zn=mz_a; break; default: break; /* switch (output_spacing_Zn){ } */}
	    vS = bget____(St_tag,mz_a); vM = bget__on(mcs_tag,mz_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_St->mc_j\n",mz_j,mz_b,mz_a);}
	    output_AnZt_S_WnYt_tmp += (*L2_get(lf_AnZt,ma_j,ma_b,ma_a,mz_j,mz_b,mz_a)) * vS * (*L2_get(lf_YnWt,ma_j,ma_b,ma_a,mz_j,mz_b,mz_a));
	    mz_j++; /* while (mz_j<M_Zn->rpop_j){ } */}
	  output_AnZt_S_WnYt->lf[tab_An + tab_St*tab_An_stride] = output_AnZt_S_WnYt_tmp;
	  /* for (mx_j=0;mx_j<M_St->rpop_j*M_An->rpop_j;mx_j++){ } */}
	/* end omp parallel */}
      /* else if (GLOBAL_omp_type==GLOBAL_omp_unused){ } */}
    /* if (lf_AnZt!=NULL && lf_YnWt!=NULL){ } */}
  if (verbose>1){ printf(" %% [finished get_AnZt_S_WnYt_vv] tidx %d\n",tidx);}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_An->rpop_j*M_Zn->rpop_j*2,0);
 skip_AnZt_S_WnYt_vv:
  if (verbose>2){ raprintf(output_AnZt_S_WnYt->lf,"double",tab_An_stride,tab_Tt_stride," %% output_AnZt_S_WnYt->lf: ");}
  return NULL;
}

int wrap_AnZt_S_WnYt_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_St,struct M_handle *M_Zn,struct L_handle *lf_AnZt,struct L_handle *lf_YnWt,struct L_handle **output_p)
{
  /* This function calls get_AnZt_S_WnYt_vv ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 10)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_AnZt_S_WnYt_vv__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Zn,verbose," %% M_Zn: ");}
  switch (output_spacing_a){
  case SPACING_j: length_a = M_An->rpop_j; break;
  case SPACING_b: length_a = M_An->rpop_b; break;
  case SPACING_a: length_a = M_An->nrows; break;
  default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){
  case SPACING_j: length_s = M_St->rpop_j; break;
  case SPACING_b: length_s = M_St->rpop_b; break;
  case SPACING_a: length_s = M_St->nrows; break;
  default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_a,length_s,length);}
  length = length_a*length_s; if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zn->mr_b,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_b: "); bprintf(M_Zn->mr_j,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_j: ");}
  if (verbose>2){ bprintf(M_Zn->mc_b,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_b: "); bprintf(M_Zn->mc_j,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_AnZt_S_WnYt_v__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_St; vpra[ip++] = M_Zn; vpra[ip++] = lf_AnZt; vpra[ip++] = lf_YnWt; vpra[ip++] = *output_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AnZt_S_WnYt_vv,vpra)){ printf("Warning! cannot create thread %d in wrap_AnZt_S_WnYt_vv__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AnZt_S_WnYt_vv(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_AnZt_S_WnYt_vv__run] tidx %d\n",*tidx);}
  return length;
}  

void wrap_AnZt_S_WnYt_vv_test()
{
  /* test for errors with input file: AnZt_S_WnYt_vv_error.in ;
  */
  /* test for speed with input file: AnZt_S_WnYt_vv_speed.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AnZt=NULL,**lf_AnAt=NULL,**lf_YnWt=NULL,**lf_YnYt=NULL;
  struct L_handle **lf_AtTn=NULL,**lf_YtTn=NULL,**lf_ZtSn=NULL,**lf_WtSn=NULL,**lf_AtTYn=NULL,**lf_ZtSWn=NULL;
  struct L_handle **lf_AnZt_S_WnYt_vv=NULL,**lf_AnAt_T_YnYt_vv=NULL; int *length_AnZt_S_WnYt_vv=NULL,*length_AnAt_T_YnYt_vv=NULL;
  struct L_handle **lf_An_ZtSWn_Yt_vv=NULL,**lf_An_AtTYn_Yt_vv=NULL; int *length_An_ZtSWn_Yt_vv=NULL,*length_An_AtTYn_Yt_vv=NULL;
  struct L_handle **lf_AnZtSWnYt_uu=NULL,**lf_AnAtTYnYt_uu=NULL; int *length_AnZtSWnYt_uu=NULL,*length_AnAtTYnYt_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_AnZt_S_WnYt_vv_test]\n");}
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
    lf_ZtSn[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_Sn[nb]->ncols);
    lf_WtSn[nb] = L_handle_make((unsigned long long int)M_Wn[nb]->ncols*(unsigned long long int)M_Sn[nb]->ncols);
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
  lf_AnAt_T_YnYt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnAt_T_YnYt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_AnZt_S_WnYt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnZt_S_WnYt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_AnAt_T_YnYt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
    lf_AnZt_S_WnYt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (error_check){ 
    lf_AtTYn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
    lf_ZtSWn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
    for (nb=0;nb<nbins;nb++){
      lf_AtTYn[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Yn[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
      lf_ZtSWn[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_Wn[nb]->ncols*(unsigned long long int)M_Sn[nb]->ncols);
      /* for (nb=0;nb<nbins;nb++){ } */}
    lf_An_AtTYn_Yt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_AtTYn_Yt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
    lf_An_ZtSWn_Yt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_ZtSWn_Yt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
    for (nb=0;nb<nbins;nb++){ 
      lf_An_AtTYn_Yt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
      lf_An_ZtSWn_Yt_vv[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
      /* for (nb=0;nb<nbins;nb++){ } */}
    lf_AnAtTYnYt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnAtTYnYt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    lf_AnZtSWnYt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnZtSWnYt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    for (nb=0;nb<nbins;nb++){ 
      lf_AnAtTYnYt_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
      lf_AnZtSWnYt_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (error_check){ } */}
  for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
	if (verbose){ printf(" %% %s; %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
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
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_ajdk_v: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
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
	  length_AnAt_T_YnYt_vv[nb] = wrap_AnZt_S_WnYt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_Tt[nb],M_An[nb],lf_AnAt[nb],lf_YnYt[nb],&(lf_AnAt_T_YnYt_vv[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
	  length_AnZt_S_WnYt_vv[nb] = wrap_AnZt_S_WnYt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_St[nb],M_Zn[nb],lf_AnZt[nb],lf_YnWt[nb],&(lf_AnZt_S_WnYt_vv[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZt_S_WnYt_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
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
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYn_vv and ZtSWn_vv: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic();
	    length_An_AtTYn_Yt_vv[nb] = wrap_An_ZtSWn_Yt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_Tt[nb],M_Yn[nb],A_ajdk,Y_ajdk,lf_AtTYn[nb],&(lf_An_AtTYn_Yt_vv[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_An_ZtSWn_Yt_vv[nb] = wrap_An_ZtSWn_Yt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_St[nb],M_Yn[nb],A_ajdk,Y_ajdk,lf_ZtSWn[nb],&(lf_An_ZtSWn_Yt_vv[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_ZtSWn_Yt_vv: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic();
	    length_AnAtTYnYt_uu[nb] = wrap_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_At[nb],M_Tt[nb],M_Yt[nb],M_Yn[nb],A_ajdk,Y_ajdk,&(lf_AnAtTYnYt_uu[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_AnZtSWnYt_uu[nb] = wrap_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_Zt[nb],M_St[nb],M_Wt[nb],M_Yn[nb],A_ajdk,Y_ajdk,&(lf_AnZtSWnYt_uu[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZtSWnYt_uu: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_AnAt_T_YnYt_vv[%d] - lf_An_AtTYn_Yt_vv[%d] error %0.16f\n",nb,nb,dra_diff(lf_AnAt_T_YnYt_vv[nb]->lf,lf_An_AtTYn_Yt_vv[nb]->lf,length_AnAt_T_YnYt_vv[nb],1));
	    printf(" %% lf_AnZt_S_WnYt_vv[%d] - lf_An_ZtSWn_Yt_vv[%d] error %0.16f\n",nb,nb,dra_diff(lf_AnZt_S_WnYt_vv[nb]->lf,lf_An_ZtSWn_Yt_vv[nb]->lf,length_AnZt_S_WnYt_vv[nb],1));
	    printf(" %% lf_AnAt_T_YnYt_vv[%d] - lf_AnAtTYnYt_uu[%d] error %0.16f\n",nb,nb,dra_diff(lf_AnAt_T_YnYt_vv[nb]->lf,lf_AnAtTYnYt_uu[nb]->lf,length_AnAt_T_YnYt_vv[nb],1));
	    printf(" %% lf_AnZt_S_WnYt_vv[%d] - lf_AnZtSWnYt_uu[%d] error %0.16f\n",nb,nb,dra_diff(lf_AnZt_S_WnYt_vv[nb]->lf,lf_AnZtSWnYt_uu[nb]->lf,length_AnZt_S_WnYt_vv[nb],1));
	    printf(" %% lf_An_AtTYn_Yt_vv[%d] - lf_AnAtTYnYt_uu[%d] error %0.16f\n",nb,nb,dra_diff(lf_An_AtTYn_Yt_vv[nb]->lf,lf_AnAtTYnYt_uu[nb]->lf,length_An_AtTYn_Yt_vv[nb],1));
	    printf(" %% lf_An_ZtSWn_Yt_vv[%d] - lf_AnZtSWnYt_uu[%d] error %0.16f\n",nb,nb,dra_diff(lf_An_ZtSWn_Yt_vv[nb]->lf,lf_AnZtSWnYt_uu[nb]->lf,length_An_ZtSWn_Yt_vv[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  /* if (error_check){ } */}
	/* for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_An_ZtSWn_Yt_vv_test]\n");}
}
