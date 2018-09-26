#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_TAnZtS_ww_stage_0(void *vp)
{
  /* This function takes in a variety of inputs and calculates a variety of single-study terms */
  int verbose=0;
  int ip=0,length=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct dcc_single *E = (struct dcc_single *)(vpra[ip++]); struct dcc_ajdk *D=E->D; 
  struct M_handle *M_An = E->M_An;
  struct M_handle *M_At = E->M_At;
  struct M_handle *M_Yn = E->M_Yn;
  struct M_handle *M_Yt = E->M_Yt;
  struct M_handle *M_Tt = E->M_Tt;
  struct M_handle *M_St = E->M_St;
  struct M_handle *M_Zt = E->M_Zt;
  struct M_handle *M_Wt = E->M_Wt;
  struct L_handle *lf_An_ajdk = E->lf_An_ajdk;
  struct L_handle *lf_Yn_ajdk = E->lf_Yn_ajdk;
  struct L_handle *lf_Zn_ajdk = E->lf_Zn_ajdk;
  struct L_handle *lf_Wn_ajdk = E->lf_Wn_ajdk;
  double *A_ajdk = D->A_ajdk;
  int A_pcols = psize(M_At->nrows);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  double *Y_ajdk = D->Y_ajdk;
  int Y_pcols = psize(M_Yt->nrows);
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  struct L_handle *lf_tmp=NULL;
  int A_ncols = M_At->nrows, Y_ncols = M_Yt->nrows, T_ncols = M_Tt->nrows, A_nrows = M_At->ncols;
  length = 1; lf_tmp = E->lf_A_a0d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = 1; lf_tmp = E->lf_A_a2d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = 1; lf_tmp = E->lf_Y_a0d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = 1; lf_tmp = E->lf_Y_a2d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a1d1_WtSn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a1d1_ZtSn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a1d1_YtTn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a1d1_AtTn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_et_Tn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_et_Sn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = Y_ncols*T_ncols; lf_tmp = E->lf_a0d1_WtSn; lf_tmp->row_stride = Y_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  length = A_ncols*T_ncols; lf_tmp = E->lf_a0d1_ZtSn; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  length = Y_ncols*T_ncols; lf_tmp = E->lf_a0d1_YtTn; lf_tmp->row_stride = Y_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  length = A_ncols*T_ncols; lf_tmp = E->lf_a0d1_AtTn; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  int ns_j=0,ns_b=0,ns_a=0;
  int na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0;
  int ny_j=0,ny_b=0,ny_a=0,my_j=0,my_b=0,my_a=0;
  int nz_j=0,nz_b=0,nz_a=0,mz_j=0,mz_b=0,mz_a=0;
  int nw_j=0,nw_b=0,nw_a=0,mw_j=0,mw_b=0,mw_a=0;
  int vA=0,vY=0,vM=0,vT=0,vS=0;
  __m128i *wAt_tag=NULL;
  __m128i *wYt_tag=NULL;
  __m128i *wTt_tag=NULL;
  __m128i *wSt_tag=NULL;
  __m128i *mcax_tag=NULL,*mcax_end=NULL;
  __m128i *mcyx_tag=NULL,*mcyx_end=NULL;
  __m128i *mctx_tag=NULL,*mctx_end=NULL;
  __m128i *wA1_tag=NULL; __m128i *wA2_tag=NULL; __m128i *wY1_tag=NULL; __m128i *wY2_tag=NULL; __m128i *mc_tag=NULL; __m128i *mc_end=NULL;
  int nc=0; 
  double A_lftmp=0,A_lftmp_a0d1=0,A_lftmp_a2d1=0;
  double Y_lftmp=0,Y_lftmp_a0d1=0,Y_lftmp_a2d1=0;
  unsigned char *At_tag=NULL;
  unsigned char *An_tag=NULL;
  unsigned char *Yt_tag=NULL;
  unsigned char *Yn_tag=NULL;
  unsigned char *Tt_tag=NULL;
  /* unsigned char *Zt_tag=NULL; */
  /* unsigned char *Wt_tag=NULL; */
  unsigned char *St_tag=NULL;
  unsigned char *mca_tag=NULL,*mca_end=NULL;
  unsigned char *mct_tag=NULL,*mct_end=NULL;
  unsigned int *ma_b_,*ma_a_;
  unsigned int *na_b_,*na_a_;
  unsigned int *my_b_,*my_a_;
  unsigned int *ny_b_,*ny_a_;
  unsigned int *mz_b_,*mz_a_;
  unsigned int *nz_b_,*nz_a_;
  unsigned int *mw_b_,*mw_a_;
  unsigned int *nw_b_,*nw_a_;
  unsigned int *ns_b_,*ns_a_;
  int QR_TYnWtS_bother = E->A_rbother && D->Y_cbother && E->Z_rbother && D->Y_cbother;
  int QR_TAnZtS_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother;
  int QR_TYnYtT_bother = E->A_rbother && D->Y_cbother && E->A_rbother && D->Y_cbother;
  int QR_TAnAtT_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother;
  int QC_TYnWtS_bother = E->A_rbother && D->Y_cbother && E->Z_rbother && D->Y_cbother;
  int QC_TAnZtS_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother;
  int QC_TYnYtT_bother = E->A_rbother && D->Y_cbother && E->A_rbother && D->Y_cbother;
  int QC_TAnAtT_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother;
  if (verbose>1){ printf(" %% [entering get_TAnZtS_ww_stage_0] tidx %d\n",tidx);}  
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Yt->mr_b,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_b: "); bprintf(M_Yt->mr_j,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_j: ");}
  if (verbose>2){ bprintf(M_Yt->mc_b,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_b: "); bprintf(M_Yt->mc_j,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  if (verbose>2){ bprintf(M_Wt->mr_b,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_b: "); bprintf(M_Wt->mr_j,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_j: ");}
  if (verbose>2){ bprintf(M_Wt->mc_b,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_b: "); bprintf(M_Wt->mc_j,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_j: ");}
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  na_b_ = M_An->n_b_; na_a_ = M_An->n_a_;
  my_b_ = M_Yn->m_b_; my_a_ = M_Yn->m_a_;
  ny_b_ = M_Yn->n_b_; ny_a_ = M_Yn->n_a_;
  mz_b_ = M_Zt->n_b_; mz_a_ = M_Zt->n_a_;
  nz_b_ = M_Zt->m_b_; nz_a_ = M_Zt->m_a_;
  mw_b_ = M_Wt->n_b_; mw_a_ = M_Wt->n_a_;
  nw_b_ = M_Wt->m_b_; nw_a_ = M_Wt->m_a_;
  ns_b_ = M_St->m_b_; ns_a_ = M_St->m_a_;
  A_lftmp_a0d1=0; for (nc=0;nc<A_pcols;nc++){ A_lftmp_a0d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_0_1*A_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
  A_lftmp_a2d1=0; for (nc=0;nc<A_pcols;nc++){ A_lftmp_a2d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_2_1*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
  if (verbose>2){ printf(" %% A_lftmp_a0d1 %f A_lftmp_a2d1 %f\n",A_lftmp_a0d1,A_lftmp_a2d1);}
  L1_set(E->lf_A_a0d1, 0,0,0, A_lftmp_a0d1);
  L1_set(E->lf_A_a2d1, 0,0,0, A_lftmp_a2d1);
  Y_lftmp_a0d1=0; for (nc=0;nc<Y_pcols;nc++){ Y_lftmp_a0d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_0_1*Y_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_Yn->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
  Y_lftmp_a2d1=0; for (nc=0;nc<Y_pcols;nc++){ Y_lftmp_a2d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_2_1*Y_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_Yn->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
  if (verbose>2){ printf(" %% Y_lftmp_a0d1 %f Y_lftmp_a2d1 %f\n",Y_lftmp_a0d1,Y_lftmp_a2d1);}
  L1_set(E->lf_Y_a0d1, 0,0,0, Y_lftmp_a0d1);
  L1_set(E->lf_Y_a2d1, 0,0,0, Y_lftmp_a2d1);
  if (QR_TAnZtS_bother || QR_TAnAtT_bother){
    if (verbose>1){ printf(" %% a1d1_ZtSn, a1d1_AtTn\n");}
    ns_j=0; while (ns_j<M_Tt->rpop_j){
      ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
      Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
      mct_tag = (unsigned char *)(&(M_Tt->mc_j[0]));
      mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length]));
      if (QR_TAnAtT_bother){ 
	ma_j=0; while (ma_j<M_At->cpop_j){
	  ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	  vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",ma_j,ma_b,ma_a);}
	  L1_plusequals(E->lf_a1d1_AtTn, ns_j,ns_b,ns_a , vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));
	  GLOBAL_ops_count_one(tidx,3,0);
	  ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
	/* if (QR_TAnAtT_bother){ } */}
      St_tag = (unsigned char *)(&(M_St->wX[ns_b*M_St->mc_length]));
      mct_tag = (unsigned char *)(&(M_St->mc_j[0]));
      mct_end = (unsigned char *)(&(M_St->mc_j[M_St->mc_length]));
      if (QR_TAnZtS_bother){
	mz_j=0; while (mz_j<M_Zt->cpop_j){
	  mz_a = mz_a_[mz_j]; mz_b = mz_b_[mz_j];
	  vS = bget____(St_tag,mz_a); vM = bget__on(mct_tag,mz_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_St->mc_j\n",mz_j,mz_b,mz_a);}
	  L1_plusequals(E->lf_a1d1_ZtSn, ns_j,ns_b,ns_a , vS*(*L2_get(lf_Zn_ajdk , mz_j,mz_b,mz_a , 0,0,AJDK_1_1)));
	  GLOBAL_ops_count_one(tidx,3,0);
	  mz_j++; /* while (mz_j<M_Zt->cpop_j){ } */}
	/* if (QR_TAnZtS_bother){ } */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_a1d1_AtTn," %% a1d1_AtTn: ");}
    if (verbose>2){ lfprintf(E->lf_a1d1_ZtSn," %% a1d1_ZtSn: ");}
    /* if (QR_TAnZtS_bother || QR_TAnAtT_bother){ } */}
  if (QR_TYnWtS_bother || QR_TYnYtT_bother){
    if (verbose>1){ printf(" %% a1d1_WtSn, a1d1_YtTn\n");}
    ns_j=0; while (ns_j<M_Tt->rpop_j){
      ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
      Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
      mct_tag = (unsigned char *)(&(M_Tt->mc_j[0]));
      mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length]));
      if (QR_TYnYtT_bother){ 
	my_j=0; while (my_j<M_Yt->cpop_j){
	  my_a = my_a_[my_j]; my_b = my_b_[my_j];
	  vT = bget____(Tt_tag,my_a); vM = bget__on(mct_tag,my_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",my_j,my_b,my_a);}
	  L1_plusequals(E->lf_a1d1_YtTn, ns_j,ns_b,ns_a , vT*(*L2_get(lf_Yn_ajdk , my_j,my_b,my_a , 0,0,AJDK_1_1)));
	  GLOBAL_ops_count_one(tidx,3,0);
	  my_j++; /* while (my_j<M_Yt->cpop_j){ } */}
	/* if (QR_TYnYtT_bother){ } */}
      St_tag = (unsigned char *)(&(M_St->wX[ns_b*M_St->mc_length]));
      mct_tag = (unsigned char *)(&(M_St->mc_j[0]));
      mct_end = (unsigned char *)(&(M_St->mc_j[M_St->mc_length]));
      if (QR_TYnWtS_bother){
	mw_j=0; while (mw_j<M_Wt->cpop_j){
	  mw_a = mw_a_[mw_j]; mw_b = mw_b_[mw_j];
	  vS = bget____(St_tag,mw_a); vM = bget__on(mct_tag,mw_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_St->mc_j\n",mw_j,mw_b,mw_a);}
	  L1_plusequals(E->lf_a1d1_WtSn, ns_j,ns_b,ns_a , vS*(*L2_get(lf_Wn_ajdk , mw_j,mw_b,mw_a , 0,0,AJDK_1_1)));
	  GLOBAL_ops_count_one(tidx,3,0);
	  mw_j++; /* while (mw_j<M_Wt->cpop_j){ } */}
	/* if (QR_TYnWtT_bother){ } */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_a1d1_YtTn," %% a1d1_YtTn: ");}
    if (verbose>2){ lfprintf(E->lf_a1d1_WtSn," %% a1d1_WtSn: ");}
    /* if (QR_TYnWtS_bother || QR_TYnYtT_bother){ } */}
  if (QR_TYnWtS_bother || QR_TAnZtS_bother || QR_TYnYtT_bother || QR_TAnAtT_bother){
    if (verbose>1){ printf(" %% et_Tn, et_Sn\n");}
    ns_j=0; while (ns_j<M_Tt->rpop_j){
      ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
      if (QR_TAnAtT_bother || QR_TYnYtT_bother){
	wTt_tag = (__m128i *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
	mctx_tag = (__m128i *)(&(M_Tt->mc_j[0]));
	mctx_end = (__m128i *)(&(M_Tt->mc_j[M_Tt->mc_length]));
	L1_set(E->lf_et_Tn , ns_j,ns_b,ns_a , -M_Tt->cpop_j + 2*popcount(&wTt_tag,&mctx_tag,&mctx_end));
	GLOBAL_ops_count_one(tidx,2,M_Tt->mc_length*BIT8);
	/* if bother */}
      if (QR_TAnZtS_bother || QR_TYnWtS_bother){
	wSt_tag = (__m128i *)(&(M_St->wX[ns_b*M_St->mc_length]));
	mctx_tag = (__m128i *)(&(M_St->mc_j[0]));
	mctx_end = (__m128i *)(&(M_St->mc_j[M_St->mc_length]));
	L1_set(E->lf_et_Sn , ns_j,ns_b,ns_a , -M_St->cpop_j + 2*popcount(&wSt_tag,&mctx_tag,&mctx_end));
	GLOBAL_ops_count_one(tidx,2,M_St->mc_length*BIT8);
	/* if bother */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_et_Tn," %% et_Tn: ");}
    if (verbose>2){ lfprintf(E->lf_et_Sn," %% et_Sn: ");}
    /* if (QR_TYnWtS_bother || QR_TAnZtS_bother || QR_TYnYtT_bother || QR_TAnAtT_bother){ } */}
  if (verbose){ printf(" %% [finished get_TAnZtS_ww_stage_0] tidx %d\n",tidx);}
  return NULL;
}

void wrap_TAnZtS_ww_stage_0(int *tidx,void **vpra,pthread_t *thread_in,struct dcc_single *E)
{
  /* This function calls get_TAnZtS_ww_stage_0 ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 2)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int nb1=0,nbins=0;
  int ip=0;
  if (verbose){ printf(" %% [entering wrap_TAnZtS_ww_stage_0] tidx %d\n",*tidx);}
  ip=0; 
  vpra[ip++] = tidx; vpra[ip++] = E; 
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_TAnZtS_ww_stage_0,vpra)){ printf("Warning! cannot create thread %d in wrap_TAnZtS_ww_stage_0\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_TAnZtS_ww_stage_0(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_TAnZtS_ww_stage_0] tidx %d\n",*tidx);}
}

void dcc_TAnZtS_ww_stage_0(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb=0; struct dcc_single *E=NULL;
  if (verbose){ printf(" %% [entering dcc_TAnZtS_ww_stage_0]\n");}
  if (verbose){ printf(" %% calculating TAnZtS_ww_stage_0.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    GLOBAL_pthread_tic();
    wrap_TAnZtS_ww_stage_0(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E);
    GLOBAL_pthread_toc();
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% TAnZtS_ww_stage_0: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcc_TAnZtS_ww_stage_0]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void dcc_TAnZtS_ww_stage_1(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; struct dcc_single *E=NULL,*E_nb1=NULL,*E_nb2=NULL; struct dcc_double *F=NULL;
  int n_spacing_A = SPACING_a;
  struct L_handle *lf_tmp=NULL;
  if (verbose){ printf(" %% [entering dcc_TAnZtS_ww_stage_1]\n");}
  if (verbose){ printf(" %% calculating TAnZtS_ww_stage_1.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0);
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    if (E->Z_rbother && D->Y_cbother){ An_a0d1_ZtSn_excerpt(E->M_Wt,E->M_St,D->Y_ajdk,E->lf_WtSn,E->lf_a0d1_WtSn);}
    if (E->Z_rbother && D->A_cbother){ An_a0d1_ZtSn_excerpt(E->M_Zt,E->M_St,D->A_ajdk,E->lf_ZtSn,E->lf_a0d1_ZtSn);}
    if (E->A_rbother && D->Y_cbother){ An_a0d1_ZtSn_excerpt(E->M_Yt,E->M_Tt,D->Y_ajdk,E->lf_YtTn,E->lf_a0d1_YtTn);}
    if (E->A_rbother && D->A_cbother){ An_a0d1_ZtSn_excerpt(E->M_At,E->M_Tt,D->A_ajdk,E->lf_AtTn,E->lf_a0d1_AtTn);}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  GLOBAL_toc(0,verbose," %% a0d1_ZtSn a0d1_AtTn: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    if (E->Z_rbother && D->Y_cbother){
      GLOBAL_pthread_tic(); 
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->M_Wt->mr_j,E->M_Wt->mr_b,E->M_St->mr_j,E->M_St->mr_b,E->lf_a0d1_WtSn,E->lf_a0d1_WtSn->lf,&(E->M_Wt->nrows),&(E->M_St->nrows),&(E->M_a0d1_WtSn),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){
      GLOBAL_pthread_tic(); 
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->M_Zt->mr_j,E->M_Zt->mr_b,E->M_St->mr_j,E->M_St->mr_b,E->lf_a0d1_ZtSn,E->lf_a0d1_ZtSn->lf,&(E->M_Zt->nrows),&(E->M_St->nrows),&(E->M_a0d1_ZtSn),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->Y_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->M_Yt->mr_j,E->M_Yt->mr_b,E->M_Tt->mr_j,E->M_Tt->mr_b,E->lf_a0d1_YtTn,E->lf_a0d1_YtTn->lf,&(E->M_Yt->nrows),&(E->M_Tt->nrows),&(E->M_a0d1_YtTn),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->M_At->mr_j,E->M_At->mr_b,E->M_Tt->mr_j,E->M_Tt->mr_b,E->lf_a0d1_AtTn,E->lf_a0d1_AtTn->lf,&(E->M_At->nrows),&(E->M_Tt->nrows),&(E->M_a0d1_AtTn),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; E_nb1 = E_[nb1]; E_nb2 = E_[nb2]; F = F_[nbx];
      lf_tmp = F->lf_Yn_a0d1_WtSn; lf_tmp->row_stride = F->E_nb1->A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = D->T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
      lf_tmp = F->lf_An_a0d1_ZtSn; lf_tmp->row_stride = F->E_nb1->A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = D->T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
      lf_tmp = F->lf_Yn_a0d1_YtTn; lf_tmp->row_stride = F->E_nb1->A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = D->T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
      lf_tmp = F->lf_An_a0d1_AtTn; lf_tmp->row_stride = F->E_nb1->A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = D->T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
      if (E_nb1->A_rbother && D->Y_cbother && E_nb2->Z_rbother && D->Y_cbother){
	GLOBAL_pthread_tic();
	wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yn,E_nb2->lf_a0d1_WtSn,E_nb2->M_a0d1_WtSn,E_nb2->M_St,&(F->lf_Yn_a0d1_WtSn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->lf_a0d1_ZtSn,E_nb2->M_a0d1_ZtSn,E_nb2->M_St,&(F->lf_An_a0d1_ZtSn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->Y_cbother && E_nb2->A_rbother && D->Y_cbother){
	GLOBAL_pthread_tic();
	wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yn,E_nb2->lf_a0d1_YtTn,E_nb2->M_a0d1_YtTn,E_nb2->M_Tt,&(F->lf_Yn_a0d1_YtTn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->lf_a0d1_AtTn,E_nb2->M_a0d1_AtTn,E_nb2->M_Tt,&(F->lf_An_a0d1_AtTn));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% F->lf_Yn_a0d1_WtSn F->lf_An_a0d1_ZtSn F->lf_Yn_a0d1_YtTn F->lf_An_a0d1_AtTn: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcc_TAnZtS_ww_stage_1]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void *get_TAnZtS_ww_stage_2(void *vp)
{
  /* This function takes in a variety of inputs and calculates a variety of single-study terms */
  int verbose=0;
  int ip=0,length=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct dcc_double *F = (struct dcc_double *)(vpra[ip++]); struct dcc_ajdk *D=F->D; struct dcc_single *E_nb1 = F->E_nb1; struct dcc_single *E_nb2 = F->E_nb2;
  double *A_ajdk = D->A_ajdk;
  double *Y_ajdk = D->Y_ajdk;
  struct L_handle *lf_tmp=NULL;
  int A_ncols = D->A_ncols, Y_ncols = D->Y_ncols, T_ncols = D->T_ncols, A_nrows = E_nb1->A_nrows;
  length = A_nrows*T_ncols; lf_tmp = F->lf_TYnWtS_ww; lf_tmp->row_stride = A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  length = A_nrows*T_ncols; lf_tmp = F->lf_TAnZtS_ww; lf_tmp->row_stride = A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  length = A_nrows*T_ncols; lf_tmp = F->lf_TYnYtT_ww; lf_tmp->row_stride = A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  length = A_nrows*T_ncols; lf_tmp = F->lf_TAnAtT_ww; lf_tmp->row_stride = A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  int ns_j=0,ns_b=0,ns_a=0;
  int na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0;
  int ny_j=0,ny_b=0,ny_a=0,my_j=0,my_b=0,my_a=0;
  int nz_j=0,nz_b=0,nz_a=0,mz_j=0,mz_b=0,mz_a=0;
  int nw_j=0,nw_b=0,nw_a=0,mw_j=0,mw_b=0,mw_a=0;
  unsigned char *Tt_tag=NULL,*mct_tag=NULL,*mct_end=NULL;
  int vT=0,vM=0;
  double lftmp=0;
  int QR_TYnWtS_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QR_TAnZtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QR_TYnYtT_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QR_TAnAtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  int QC_TYnWtS_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QC_TAnZtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QC_TYnYtT_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QC_TAnAtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  if (verbose>1){ printf(" %% [entering get_TAnZtS_ww_stage_2] tidx %d\n",tidx);}  
  if (QR_TAnZtS_bother || QR_TAnAtT_bother){
    if (verbose>1){ printf(" %% TAnZtS and TAnAtT\n");}
    ns_j=0; while (ns_j<E_nb1->M_Tt->rpop_j){
      ns_a = E_nb1->M_Tt->m_a_[ns_j]; ns_b = E_nb1->M_Tt->m_b_[ns_j];
      if (QR_TAnZtS_bother){ 
	ma_j=0; while (ma_j<E_nb1->M_An->rpop_j){
	  ma_a = E_nb1->M_An->m_a_[ma_j]; ma_b = E_nb1->M_An->m_b_[ma_j];
	  lftmp = (*L2_get(F->lf_An_a0d1_ZtSn , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)) - (*L1_get(E_nb2->lf_a1d1_ZtSn, ns_j,ns_b,ns_a)) - (*L2_get(E_nb1->lf_An_ajdk, ma_j,ma_b,ma_a , 0,0,AJDK_1_1))*(*L1_get(E_nb2->lf_et_Sn, ns_j,ns_b,ns_a)) + (*L1_get(E_nb1->lf_A_a2d1, 0,0,0))*(*L1_get(E_nb2->lf_et_Sn, ns_j,ns_b,ns_a));
	  L2_plusequals(F->lf_TAnZtS_ww, ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , +lftmp);
	  GLOBAL_ops_count_one(tidx,1,0);
	  ma_j++; /* while (ma_j<E_nb1->M_An->rpop_j){ } */}
	/* if (QR_TAnZtS_bother){ } */}
      if (QR_TAnAtT_bother){
	ma_j=0; while (ma_j<E_nb1->M_An->rpop_j){
	  ma_a = E_nb1->M_An->m_a_[ma_j]; ma_b = E_nb1->M_An->m_b_[ma_j];
	  lftmp = (*L2_get(F->lf_An_a0d1_AtTn , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)) - (*L1_get(E_nb2->lf_a1d1_AtTn, ns_j,ns_b,ns_a)) - (*L2_get(E_nb1->lf_An_ajdk, ma_j,ma_b,ma_a , 0,0,AJDK_1_1))*(*L1_get(E_nb2->lf_et_Tn, ns_j,ns_b,ns_a)) + (*L1_get(E_nb1->lf_A_a2d1, 0,0,0))*(*L1_get(E_nb2->lf_et_Tn, ns_j,ns_b,ns_a));
	  L2_plusequals(F->lf_TAnAtT_ww, ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , +lftmp);
	  GLOBAL_ops_count_one(tidx,1,0);
	  ma_j++; /* while (ma_j<E_nb1->M_An->rpop_j){ } */}
	/* if (QR_TAnAtT_bother){ } */}
      ns_j++; /* while (ns_j<E_nb1->M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(F->lf_TAnZtS_ww," %% TAnZtS_ww: ");}
    if (verbose>2){ lfprintf(F->lf_TAnAtT_ww," %% TAnAtT_ww: ");}
    /* if (QR_TAnZtS_bother || QR_TAnAtT_bother){ } */}
  if (QR_TAnAtT_bother && F->nb1==F->nb2){
    ns_j=0; while (ns_j<E_nb1->M_Tt->rpop_j){
      ns_a = E_nb1->M_Tt->m_a_[ns_j]; ns_b = E_nb1->M_Tt->m_b_[ns_j];
      Tt_tag = (unsigned char *)(&(E_nb1->M_Tt->wX[ns_b*E_nb1->M_Tt->mc_length]));
      mct_tag = (unsigned char *)(&(E_nb1->M_Tt->mc_j[0]));
      mct_end = (unsigned char *)(&(E_nb1->M_Tt->mc_j[E_nb1->M_Tt->mc_length]));
      ma_j=0; while (ma_j<E_nb1->M_An->rpop_j){
  	ma_a = E_nb1->M_An->m_a_[ma_j]; ma_b = E_nb1->M_An->m_b_[ma_j];
  	vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",ma_j,ma_b,ma_a);}
  	lftmp = vT * ( (*L1_get(E_nb1->lf_A_a0d1,0,0,0)) - 2*(*L2_get(E_nb1->lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)) + (*L1_get(E_nb1->lf_A_a2d1,0,0,0)) );
  	L2_plusequals(F->lf_TAnAtT_ww, ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , -lftmp);
  	GLOBAL_ops_count_one(tidx,1,0);
  	ma_j++; /* while (ma_j<E_nb1->M_An->rpop_j){ } */}
      ns_j++; /* while (ns_j<E_nb1->M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(F->lf_TAnAtT_ww," %% TAnAtT_ww: ");}
    /* if bother */}
  if (QR_TYnWtS_bother || QR_TYnYtT_bother){
    if (verbose>1){ printf(" %% TYnWtS and TYnYtT\n");}
    ns_j=0; while (ns_j<E_nb1->M_Tt->rpop_j){
      ns_a = E_nb1->M_Tt->m_a_[ns_j]; ns_b = E_nb1->M_Tt->m_b_[ns_j];
      if (QR_TYnWtS_bother){ 
	my_j=0; while (my_j<E_nb1->M_Yn->rpop_j){
	  my_a = E_nb1->M_Yn->m_a_[my_j]; my_b = E_nb1->M_Yn->m_b_[my_j];
	  lftmp = (*L2_get(F->lf_Yn_a0d1_WtSn , my_j,my_b,my_a , ns_j,ns_b,ns_a)) - (*L1_get(E_nb2->lf_a1d1_WtSn, ns_j,ns_b,ns_a)) - (*L2_get(E_nb1->lf_Yn_ajdk, my_j,my_b,my_a , 0,0,AJDK_1_1))*(*L1_get(E_nb2->lf_et_Sn, ns_j,ns_b,ns_a)) + (*L1_get(E_nb1->lf_Y_a2d1, 0,0,0))*(*L1_get(E_nb2->lf_et_Sn, ns_j,ns_b,ns_a));
	  L2_plusequals(F->lf_TYnWtS_ww, my_j,my_b,my_a , ns_j,ns_b,ns_a , +lftmp);
	  GLOBAL_ops_count_one(tidx,1,0);
	  my_j++; /* while (my_j<E_nb1->M_Yn->rpop_j){ } */}
	/* if (QR_TYnWtS_bother){ } */}
      if (QR_TYnYtT_bother){
	my_j=0; while (my_j<E_nb1->M_Yn->rpop_j){
	  my_a = E_nb1->M_Yn->m_a_[my_j]; my_b = E_nb1->M_Yn->m_b_[my_j];
	  lftmp = (*L2_get(F->lf_Yn_a0d1_YtTn , my_j,my_b,my_a , ns_j,ns_b,ns_a)) - (*L1_get(E_nb2->lf_a1d1_YtTn, ns_j,ns_b,ns_a)) - (*L2_get(E_nb1->lf_Yn_ajdk, my_j,my_b,my_a , 0,0,AJDK_1_1))*(*L1_get(E_nb2->lf_et_Tn, ns_j,ns_b,ns_a)) + (*L1_get(E_nb1->lf_Y_a2d1, 0,0,0))*(*L1_get(E_nb2->lf_et_Tn, ns_j,ns_b,ns_a));
	  L2_plusequals(F->lf_TYnYtT_ww, my_j,my_b,my_a , ns_j,ns_b,ns_a , +lftmp);
	  GLOBAL_ops_count_one(tidx,1,0);
	  my_j++; /* while (my_j<E_nb1->M_Yn->rpop_j){ } */}
	/* if (QR_TYnYtT_bother){ } */}
      ns_j++; /* while (ns_j<E_nb1->M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(F->lf_TYnWtS_ww," %% TYnWtS_ww: ");}
    if (verbose>2){ lfprintf(F->lf_TYnYtT_ww," %% TYnYtT_ww: ");}
    /* if (QR_TYnWtS_bother || QR_TYnYtT_bother){ } */}
  if (QR_TYnYtT_bother && F->nb1==F->nb2){
    ns_j=0; while (ns_j<E_nb1->M_Tt->rpop_j){
      ns_a = E_nb1->M_Tt->m_a_[ns_j]; ns_b = E_nb1->M_Tt->m_b_[ns_j];
      Tt_tag = (unsigned char *)(&(E_nb1->M_Tt->wX[ns_b*E_nb1->M_Tt->mc_length]));
      mct_tag = (unsigned char *)(&(E_nb1->M_Tt->mc_j[0]));
      mct_end = (unsigned char *)(&(E_nb1->M_Tt->mc_j[E_nb1->M_Tt->mc_length]));
      my_j=0; while (my_j<E_nb1->M_Yn->rpop_j){
  	my_a = E_nb1->M_Yn->m_a_[my_j]; my_b = E_nb1->M_Yn->m_b_[my_j];
  	vT = bget____(Tt_tag,my_a); vM = bget__on(mct_tag,my_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",my_j,my_b,my_a);}
  	lftmp = vT * ( (*L1_get(E_nb1->lf_Y_a0d1,0,0,0)) - 2*(*L2_get(E_nb1->lf_Yn_ajdk , my_j,my_b,my_a , 0,0,AJDK_1_1)) + (*L1_get(E_nb1->lf_Y_a2d1,0,0,0)) );
  	L2_plusequals(F->lf_TYnYtT_ww, my_j,my_b,my_a , ns_j,ns_b,ns_a , -lftmp);
  	GLOBAL_ops_count_one(tidx,1,0);
  	my_j++; /* while (my_j<E_nb1->M_Yn->rpop_j){ } */}
      ns_j++; /* while (ns_j<E_nb1->M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(F->lf_TYnYtT_ww," %% TYnYtT_ww: ");}
    /* if bother */}
  if (verbose){ printf(" %% [finished get_TAnZtS_ww_stage_2] tidx %d\n",tidx);}
  return NULL;
}

void wrap_TAnZtS_ww_stage_2(int *tidx,void **vpra,pthread_t *thread_in,struct dcc_double *F)
{
  /* This function calls get_TAnZtS_ww_stage_2 ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 2)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int nb1=0,nbins=0;
  int ip=0;
  if (verbose){ printf(" %% [entering wrap_TAnZtS_ww_stage_2] tidx %d\n",*tidx);}
  ip=0; 
  vpra[ip++] = tidx; vpra[ip++] = F; 
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_TAnZtS_ww_stage_2,vpra)){ printf("Warning! cannot create thread %d in wrap_TAnZtS_ww_stage_2\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_TAnZtS_ww_stage_2(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_TAnZtS_ww_stage_2] tidx %d\n",*tidx);}
}

void dcc_TAnZtS_ww_stage_2(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; struct dcc_double *F=NULL;
  if (verbose){ printf(" %% [entering dcc_TAnZtS_ww_stage_2]\n");}
  if (verbose){ printf(" %% calculating TAnZtS_ww_stage_2.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx=nb1+nb2*nbins; F = F_[nbx];
    GLOBAL_pthread_tic();
    wrap_TAnZtS_ww_stage_2(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),F);
    GLOBAL_pthread_toc();
    /*  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% TAnZtS_ww_stage_2: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcc_TAnZtS_ww_stage_2]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void *get_TAnZtS_uu(void *vp)
{
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_TAnZtS_uu = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); int output_spacing_An = output_spacing_a; 
  int output_spacing_s = *(int *)(vpra[ip++]); int output_spacing_Tt = output_spacing_s;
  int ncols_A_p = psize(M_An->ncols) ;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*ncols_A_p]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*ncols_A_p]);
  int ns_j=0,ns_b=0,ns_a=0,tab_Tt_stride=0,tab_Tt=0,ma_j=0,ma_b=0,ma_a=0,tab_An_stride=0,tab_An=0,mz_j=0,mz_b=0,mz_a=0,na_j=0,na_b=0,na_a=0,nz_j=0,nz_b=0,nz_a=0;
  double output_tmp=0; 
  unsigned char *An_tag; unsigned char *Zn_tag; 
  unsigned char *Tt_tag; unsigned char *St_tag;
  int vA=0,vT=0,vZ=0,vS=0;
  if (verbose>1){ printf(" %% [entering get_TAnZtS_uu] tidx %d \n",tidx);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  switch (output_spacing_Tt){ case SPACING_j: tab_Tt_stride = M_Tt->rpop_j; break; case SPACING_b: tab_Tt_stride = M_Tt->rpop_b; break; case SPACING_a: tab_Tt_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_Tt){ } */}
  switch (output_spacing_An){ case SPACING_j: tab_An_stride = M_An->rpop_j; break; case SPACING_b: tab_An_stride = M_An->rpop_b; break; case SPACING_a: tab_An_stride = M_An->nrows; break; default: break; /* switch (output_spacing_An){ } */}
  output_TAnZtS_uu->spacing_row = output_spacing_An; output_TAnZtS_uu->row_stride = tab_An_stride;
  output_TAnZtS_uu->spacing_col = output_spacing_Tt; output_TAnZtS_uu->col_stride = tab_Tt_stride;
  fill_uchar_zero((unsigned char *)(output_TAnZtS_uu->lf),tab_An_stride*tab_Tt_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"TAnZtS_uu")){ goto skip_TAnZtS_uu;}
  ns_j=0; while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_Tt){ case SPACING_j: tab_Tt=ns_j; break; case SPACING_b: tab_Tt=ns_b; break; case SPACING_a: tab_Tt=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
    St_tag = (unsigned char*)(&(M_St->wX[ns_b*M_St->mc_length]));
    Tt_tag = (unsigned char*)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
    ma_j=0; while (ma_j<M_An->rpop_j){
      ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
      An_tag = (unsigned char*)(&(M_An->wX[ma_b*M_An->mc_length]));
      switch (output_spacing_An){ case SPACING_j: tab_An=ma_j; break; case SPACING_b: tab_An=ma_b; break; case SPACING_a: tab_An=ma_a; break; default: break; /* switch (output_spacing_An){ } */}
      output_tmp=0;
      mz_j=0; while (mz_j<M_Zn->rpop_j){
	mz_a = M_Zn->m_a_[mz_j]; mz_b = M_Zn->m_b_[mz_j];
	Zn_tag = (unsigned char*)(&(M_Zn->wX[mz_b*M_Zn->mc_length]));
	na_j=0; while (na_j<M_An->cpop_j){
	  na_a = M_An->n_a_[na_j]; na_b = M_An->n_b_[na_j];
	  vA = bget____(An_tag,na_a);
	  vT = bget____(Tt_tag,ma_a);
	  vZ = bget____(Zn_tag,na_a);
	  vS = bget____(St_tag,mz_a);
	  output_tmp += /* vT* */(vA - a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*(vZ-a_An[na_a/POPLENGTH])*vS;
	  na_j++; /* while (na_j<M_An->cpop_j){ } */}
	mz_j++; /* while (mz_j<M_Zn->rpop_j){ } */}
      output_TAnZtS_uu->lf[tab_An + tab_Tt*tab_An_stride] = output_tmp;
      ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_An->rpop_j*M_An->cpop_j*M_Zn->rpop_j*3,0);
  if (verbose>2){ lfprintf(output_TAnZtS_uu," %% output_TAnZtS_uu: ");}
 skip_TAnZtS_uu:
  if (verbose>1){ printf(" %% [finished get_TAnZtS_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_TAnZtS_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_Tt,struct M_handle *M_Zn,struct M_handle *M_St,double *A_ajdk,struct L_handle **output_p)
 {
   /* This function calls get_TAnZtS_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 9)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_TAnZtS_uu__run] tidx %d\n",*tidx);}
  if (verbose>2){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose>2){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose>2){ M_handle_printf(M_Zn,verbose," %% M_Zn: ");}
  if (verbose>2){ M_handle_printf(M_St,verbose," %% M_St: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Zn->mr_b,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_b: "); bprintf(M_Zn->mr_j,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_j: ");}
  if (verbose>2){ bprintf(M_Zn->mc_b,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_b: "); bprintf(M_Zn->mc_j,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_TAnZtS_uu__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_Tt; vpra[ip++] = M_Zn; vpra[ip++] = M_St; vpra[ip++] = A_ajdk; vpra[ip++] = *output_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_TAnZtS_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_TAnZtS_uu\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_TAnZtS_uu(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_TAnZtS_uu__run] tidx %d\n",*tidx);}
  return length;
}  

void *get_TAnAtT_uu(void *vp)
{
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_TAnAtT_uu = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); int output_spacing_An = output_spacing_a; 
  int output_spacing_s = *(int *)(vpra[ip++]); int output_spacing_Tt = output_spacing_s;
  int ncols_A_p = psize(M_An->ncols) ;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*ncols_A_p]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*ncols_A_p]);
  int ns_j=0,ns_b=0,ns_a=0,tab_Tt_stride=0,tab_Tt=0,ma_j=0,ma_b=0,ma_a=0,tab_An_stride=0,tab_An=0,mz_j=0,mz_b=0,mz_a=0,na_j=0,na_b=0,na_a=0,nz_j=0,nz_b=0,nz_a=0;
  double output_tmp=0; 
  unsigned char *An_tag; unsigned char *Zn_tag; 
  unsigned char *Tt_tag; unsigned char *St_tag;
  int vA=0,vT=0,vZ=0,vS=0;
  if (verbose>1){ printf(" %% [entering get_TAnAtT_uu] tidx %d \n",tidx);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  switch (output_spacing_Tt){ case SPACING_j: tab_Tt_stride = M_Tt->rpop_j; break; case SPACING_b: tab_Tt_stride = M_Tt->rpop_b; break; case SPACING_a: tab_Tt_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_Tt){ } */}
  switch (output_spacing_An){ case SPACING_j: tab_An_stride = M_An->rpop_j; break; case SPACING_b: tab_An_stride = M_An->rpop_b; break; case SPACING_a: tab_An_stride = M_An->nrows; break; default: break; /* switch (output_spacing_An){ } */}
  output_TAnAtT_uu->spacing_row = output_spacing_An; output_TAnAtT_uu->row_stride = tab_An_stride;
  output_TAnAtT_uu->spacing_col = output_spacing_Tt; output_TAnAtT_uu->col_stride = tab_Tt_stride;
  fill_uchar_zero((unsigned char *)(output_TAnAtT_uu->lf),tab_An_stride*tab_Tt_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"TAnAtT_uu")){ goto skip_TAnAtT_uu;}
  ns_j=0; while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_Tt){ case SPACING_j: tab_Tt=ns_j; break; case SPACING_b: tab_Tt=ns_b; break; case SPACING_a: tab_Tt=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
    St_tag = (unsigned char*)(&(M_St->wX[ns_b*M_St->mc_length]));
    Tt_tag = (unsigned char*)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
    ma_j=0; while (ma_j<M_An->rpop_j){
      ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
      An_tag = (unsigned char*)(&(M_An->wX[ma_b*M_An->mc_length]));
      switch (output_spacing_An){ case SPACING_j: tab_An=ma_j; break; case SPACING_b: tab_An=ma_b; break; case SPACING_a: tab_An=ma_a; break; default: break; /* switch (output_spacing_An){ } */}
      output_tmp=0;
      mz_j=0; while (mz_j<M_Zn->rpop_j){
	mz_a = M_Zn->m_a_[mz_j]; mz_b = M_Zn->m_b_[mz_j];
	Zn_tag = (unsigned char*)(&(M_Zn->wX[mz_b*M_Zn->mc_length]));
	na_j=0; while (na_j<M_An->cpop_j){
	  na_a = M_An->n_a_[na_j]; na_b = M_An->n_b_[na_j];
	  vA = bget____(An_tag,na_a);
	  vT = bget____(Tt_tag,ma_a);
	  vZ = bget____(Zn_tag,na_a);
	  vS = bget____(St_tag,mz_a);
	  if (ma_j!=mz_j){ output_tmp += /* vT* */(vA - a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*(vZ-a_An[na_a/POPLENGTH])*vS;}
	  na_j++; /* while (na_j<M_An->cpop_j){ } */}
	mz_j++; /* while (mz_j<M_Zn->rpop_j){ } */}
      output_TAnAtT_uu->lf[tab_An + tab_Tt*tab_An_stride] = output_tmp;
      ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_An->rpop_j*M_An->cpop_j*M_Zn->cpop_j*3,0);
  if (verbose>2){ lfprintf(output_TAnAtT_uu," %% output_TAnAtT_uu: ");}
 skip_TAnAtT_uu:
  if (verbose>1){ printf(" %% [finished get_TAnAtT_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_TAnAtT_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_Tt,struct M_handle *M_Zn,struct M_handle *M_St,double *A_ajdk,struct L_handle **output_p)
 {
   /* This function calls get_TAnAtT_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 9)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_TAnAtT_uu__run] tidx %d\n",*tidx);}
  if (verbose>2){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose>2){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose>2){ M_handle_printf(M_Zn,verbose," %% M_Zn: ");}
  if (verbose>2){ M_handle_printf(M_St,verbose," %% M_St: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Zn->mr_b,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_b: "); bprintf(M_Zn->mr_j,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_j: ");}
  if (verbose>2){ bprintf(M_Zn->mc_b,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_b: "); bprintf(M_Zn->mc_j,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_TAnAtT_uu__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_Tt; vpra[ip++] = M_Zn; vpra[ip++] = M_St; vpra[ip++] = A_ajdk; vpra[ip++] = *output_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_TAnAtT_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_TAnAtT_uu\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_TAnAtT_uu(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_TAnAtT_uu__run] tidx %d\n",*tidx);}
  return length;
}  

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

