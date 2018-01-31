void *get_singlestudy_uu(void *vp)
{
  /* This function takes in a variety of inputs and calculates a variety of single-study terms */
  int verbose=0;
  int ip=0,length=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct bcc_single *E = (struct bcc_single *)(vpra[ip++]); struct bcc_ajdk *D=E->D; 
  struct M_handle *M_An = E->M_An;
  struct M_handle *M_At = E->M_At;
  struct M_handle *M_Tt = E->M_Tt;
  struct M_handle *M_Yn = E->M_Yn;
  struct M_handle *M_Yt = E->M_Yt;
  /* struct M_handle *M_Wt = E->M_Wt; */
  struct M_handle *M_St = E->M_St;
  struct M_handle *M_Zt = E->M_Zt;
  struct L_handle *lf_An_ajdk = E->lf_An_ajdk;
  struct L_handle *lf_Yn_ajdk = E->lf_Yn_ajdk;
  /* struct L_handle *lf_Wn_ajdk = E->lf_Wn_ajdk; */
  struct L_handle *lf_Zn_ajdk = E->lf_Zn_ajdk;
  double *A_ajdk = D->A_ajdk;
  double *Y_ajdk = D->Y_ajdk;
  int A_pcols = psize(M_At->nrows);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yt->nrows);
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  struct L_handle *lf_tmp=NULL;
  int A_ncols = M_At->nrows, T_ncols = M_Tt->nrows, A_nrows = M_At->ncols;
  length = A_ncols; lf_tmp = E->lf_At_Yn_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = A_ncols; lf_tmp = E->lf_At_An_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a1d2_Zt_Sn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a2d2_Zt_Sn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a3d2_Zt_Sn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a1d2_At_Tn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a2d2_At_Tn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a3d2_At_Tn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_et_Tn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_et_Sn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = A_ncols; lf_tmp = E->lf_et_An; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = A_ncols; lf_tmp = E->lf_a1d1_At_en; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = A_ncols; lf_tmp = E->lf_a1d2_At_en; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = A_ncols; lf_tmp = E->lf_a3d2_At_en; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = 1; lf_tmp = E->lf_et_Yn_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_Tt_Yn_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  /* length = A_ncols*T_ncols; lf_tmp = E->lf_At_T_Yn_a1d1; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp); */
  length = 1; lf_tmp = E->lf_et_An_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_Tt_An_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  /* length = A_ncols*T_ncols; lf_tmp = E->lf_At_T_An_a1d1; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp); */
  length = T_ncols*A_nrows; lf_tmp = E->lf_T_AnAt_YnYt; lf_tmp->row_stride = T_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = A_nrows; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  length = T_ncols*A_nrows; lf_tmp = E->lf_T_AnAt_AnAt; lf_tmp->row_stride = T_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = A_nrows; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  int ns_j=0,ns_b=0,ns_a=0;
  int na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0;
  int ny_j=0,ny_b=0,ny_a=0,my_j=0,my_b=0,my_a=0;
  /* int nw_j=0,nw_b=0,nw_a=0,mw_j=0,mw_b=0,mw_a=0; */
  int nz_j=0,nz_b=0,nz_a=0,mz_j=0,mz_b=0,mz_a=0;
  int vA=0,vY=0,vM=0,vT=0,vS=0;
  __m128i *wAt_tag=NULL;
  __m128i *wTt_tag=NULL;
  __m128i *wSt_tag=NULL;
  __m128i *mcax_tag=NULL,*mcax_end=NULL;
  __m128i *mctx_tag=NULL,*mctx_end=NULL;
  __m128i *wA1_tag=NULL; __m128i *wA2_tag=NULL; __m128i *wY1_tag=NULL; __m128i *wY2_tag=NULL; __m128i *mc_tag=NULL; __m128i *mc_end=NULL;
  int nc=0; double A_dtmp=0,A_dtmp_a0d1=0,A_dtmp_a2d1=0,Y_dtmp=0,Y_dtmp_a0d1=0,Y_dtmp_a2d1=0,*dinp=NULL;
  unsigned char *At_tag=NULL;
  unsigned char *An_tag=NULL;
  unsigned char *Tt_tag=NULL;
  unsigned char *Yn_tag=NULL;
  /* unsigned char *Yt_tag=NULL; */
  /* unsigned char *Wt_tag=NULL; */
  /* unsigned char *Zt_tag=NULL; */
  unsigned char *St_tag=NULL;
  unsigned char *mca_tag=NULL,*mca_end=NULL;
  unsigned char *mcy_tag=NULL,*mcy_end=NULL;
  unsigned char *mct_tag=NULL,*mct_end=NULL;
  unsigned int *ma_b_,*ma_a_;
  unsigned int *na_b_,*na_a_;
  unsigned int *my_b_,*my_a_;
  unsigned int *ny_b_,*ny_a_;
  /* unsigned int *mw_b_,*mw_; */
  /* unsigned int *nw_b_,*nw_; */
  unsigned int *mz_b_,*mz_a_;
  unsigned int *nz_b_,*nz_a_;
  unsigned int *ns_b_,*ns_a_;
  int QR_AnZtSWnYt_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother;
  int QR_AnZtSZnAt_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother;
  int QR_AnAtTYnYt_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother;
  int QR_AnAtTAnAt_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother;
  int QC_AtTYnWtSZn_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother;
  int QC_AtTAnZtSZn_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother;
  int QC_AtTYnYtTAn_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother;
  int QC_AtTAnAtTAn_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother;
  if (verbose>1){ printf(" %% [entering get_singlestudy_uu] tidx %d\n",tidx);}  
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Yt->mr_b,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_b: "); bprintf(M_Yt->mr_j,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_j: ");}
  if (verbose>2){ bprintf(M_Yt->mc_b,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_b: "); bprintf(M_Yt->mc_j,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_j: ");}
  /* if (verbose>2){ bprintf(M_Wt->mr_b,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_b: "); bprintf(M_Wt->mr_j,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_j: ");} */
  /* if (verbose>2){ bprintf(M_Wt->mc_b,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_b: "); bprintf(M_Wt->mc_j,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_j: ");} */
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  na_b_ = M_An->n_b_; na_a_ = M_An->n_a_;
  my_b_ = M_Yn->m_b_; my_a_ = M_Yn->m_a_;
  ny_b_ = M_Yn->n_b_; ny_a_ = M_Yn->n_a_;
  /* mw_b_ = M_Wt->m_b_; mw_ = M_Wt->m_a_; */
  /* nw_b_ = M_Wt->n_b_; nw_ = M_Wt->n_a_; */
  mz_b_ = M_Zt->n_b_; mz_a_ = M_Zt->n_a_;
  nz_b_ = M_Zt->m_b_; nz_a_ = M_Zt->m_a_;
  ns_b_ = M_St->m_b_; ns_a_ = M_St->m_a_;
  if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){
    if (verbose>1){ printf(" %% Tt_Yn_a1d1, Tt_An_a1d1\n");}
    ns_j=0; while (ns_j<M_Tt->rpop_j){
      ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
      Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
      mct_tag = (unsigned char *)(&(M_Tt->mc_j[0]));
      mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length]));
      ma_j=0; while (ma_j<M_At->cpop_j){
	ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",ma_j,ma_b,ma_a);}
	if (QC_AtTYnYtTAn_bother){ L1_plusequals(E->lf_Tt_Yn_a1d1 , ns_j,ns_b,ns_a , vT*(*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
	if (QC_AtTAnAtTAn_bother){ L1_plusequals(E->lf_Tt_An_a1d1 , ns_j,ns_b,ns_a , vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
	  GLOBAL_ops_count_one(tidx,4,0);
	ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_Tt_Yn_a1d1," %% Tt_Yn_a1d1: ");}
    if (verbose>2){ lfprintf(E->lf_Tt_An_a1d1," %% Tt_An_a1d1: ");}
    /* if (verbose>1){ printf(" %% At_T_Yn_a1d1, At_T_An_a1d1\n");} */
    /* ns_j=0; while (ns_j<M_Tt->rpop_j){ */
    /*   ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j]; */
    /*   Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length])); */
    /*   mct_tag = (unsigned char *)(&(M_Tt->mc_j[0])); */
    /*   mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length])); */
    /*   na_j=0; while (na_j<M_At->rpop_j){ */
    /* 	na_a = na_a_[na_j]; na_b = na_b_[na_j]; */
    /* 	At_tag = (unsigned char *)(&(M_At->wX[na_b*M_At->mc_length])); */
    /* 	mca_tag = (unsigned char *)(&(M_At->mc_j[0])); */
    /* 	mca_end = (unsigned char *)(&(M_At->mc_j[M_At->mc_length])); */
    /* 	ma_j=0; while (ma_j<M_At->cpop_j){ */
    /* 	  ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j]; */
    /* 	  vA = bget____(At_tag,ma_a); vM = bget__on(mca_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_At->mc_j\n",ma_j,ma_b,ma_a);} */
    /* 	  vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_At->mc_j\n",ma_j,ma_b,ma_a);} */
    /* 	  if (QC_AtTYnYtTAn_bother){ L2_plusequals(E->lf_At_T_Yn_a1d1 , na_j,na_b,na_a , ns_j,ns_b,ns_a , vA*vT*(*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
    /* 	  if (QC_AtTAnAtTAn_bother){ L2_plusequals(E->lf_At_T_An_a1d1 , na_j,na_b,na_a , ns_j,ns_b,ns_a , vA*vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
    /* 	  GLOBAL_ops_count_one(tidx,4,0); */
    /* 	  ma_j++; /\* while (ma_j<M_At->cpop_j){ } *\/} */
    /* 	na_j++; /\* while (na_j<M_At->rpop_j){ } *\/} */
    /*   ns_j++; /\* while (ns_j<M_Tt->rpop_j){ } *\/} */
    /* if (verbose>2){ lfprintf(E->lf_At_T_Yn_a1d1," %% At_T_Yn_a1d1: ");} */
    /* if (verbose>2){ lfprintf(E->lf_At_T_An_a1d1," %% At_T_An_a1d1: ");} */
    /* if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){
    if (verbose>1){ printf(" %% At_Yn_a1d1, At_An_a1d1\n");}
    na_j=0; while (na_j<M_At->rpop_j){
      na_a = na_a_[na_j]; na_b = na_b_[na_j];
      At_tag = (unsigned char *)(&(M_At->wX[na_b*M_At->mc_length]));
      mca_tag = (unsigned char *)(&(M_At->mc_j[0]));
      mca_end = (unsigned char *)(&(M_At->mc_j[M_At->mc_length]));
      ma_j=0; while (ma_j<M_At->cpop_j){
	ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	vA = bget____(At_tag,ma_a); vM = bget__on(mca_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_At->mc_j\n",ma_j,ma_b,ma_a);}
	if (QC_AtTYnYtTAn_bother){ L1_plusequals(E->lf_At_Yn_a1d1 , na_j,na_b,na_a , vA*(*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
	if (QC_AtTAnAtTAn_bother){ L1_plusequals(E->lf_At_An_a1d1 , na_j,na_b,na_a , vA*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
	  GLOBAL_ops_count_one(tidx,2,0);
	ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
      na_j++; /* while (na_j<M_At->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_At_Yn_a1d1," %% At_Yn_a1d1: ");}
    if (verbose>2){ lfprintf(E->lf_At_An_a1d1," %% At_An_a1d1: ");}
    /* if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QR_AnZtSZnAt_bother || QR_AnAtTAnAt_bother){
    if (verbose>1){ printf(" %% a1d2_Zt_Sn, a2d2_Zt_Sn, a3d2_Zt_Sn, a1d2_At_Tn, a2d2_At_Tn, a3d2_At_Tn\n");}
    ns_j=0; while (ns_j<M_Tt->rpop_j){
      ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
      Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
      mct_tag = (unsigned char *)(&(M_Tt->mc_j[0]));
      mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length]));
      if (QR_AnAtTAnAt_bother){ 
	ma_j=0; while (ma_j<M_At->cpop_j){
	  ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	  vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",ma_j,ma_b,ma_a);}
	  L1_plusequals(E->lf_a1d2_At_Tn, ns_j,ns_b,ns_a , vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_2)));
	  L1_plusequals(E->lf_a2d2_At_Tn, ns_j,ns_b,ns_a , vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_2_2)));
	  L1_plusequals(E->lf_a3d2_At_Tn, ns_j,ns_b,ns_a , vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_3_2)));
	  GLOBAL_ops_count_one(tidx,3,0);
	  ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
	/* if (QR_AnAtTAnAt_bother){ } */}
      St_tag = (unsigned char *)(&(M_St->wX[ns_b*M_St->mc_length]));
      mct_tag = (unsigned char *)(&(M_St->mc_j[0]));
      mct_end = (unsigned char *)(&(M_St->mc_j[M_St->mc_length]));
      if (QR_AnZtSZnAt_bother){
	mz_j=0; while (mz_j<M_Zt->cpop_j){
	  mz_a = mz_a_[mz_j]; mz_b = mz_b_[mz_j];
	  vS = bget____(St_tag,mz_a); vM = bget__on(mct_tag,mz_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_St->mc_j\n",mz_j,mz_b,mz_a);}
	  L1_plusequals(E->lf_a1d2_Zt_Sn, ns_j,ns_b,ns_a , vS*(*L2_get(lf_Zn_ajdk , mz_j,mz_b,mz_a , 0,0,AJDK_1_2)));
	  L1_plusequals(E->lf_a2d2_Zt_Sn, ns_j,ns_b,ns_a , vS*(*L2_get(lf_Zn_ajdk , mz_j,mz_b,mz_a , 0,0,AJDK_2_2)));
	  L1_plusequals(E->lf_a3d2_Zt_Sn, ns_j,ns_b,ns_a , vS*(*L2_get(lf_Zn_ajdk , mz_j,mz_b,mz_a , 0,0,AJDK_3_2)));
	  GLOBAL_ops_count_one(tidx,3,0);
	  mz_j++; /* while (mz_j<M_Tt->cpop_j){ } */}
	/* if (QR_AnZtSZnAt_bother){ } */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_a1d2_At_Tn," %% a1d2_At_Tn: ");}
    if (verbose>2){ lfprintf(E->lf_a2d2_At_Tn," %% a2d2_At_Tn: ");}
    if (verbose>2){ lfprintf(E->lf_a3d2_At_Tn," %% a3d2_At_Tn: ");}
    if (verbose>2){ lfprintf(E->lf_a1d2_Zt_Sn," %% a1d2_Zt_Sn: ");}
    if (verbose>2){ lfprintf(E->lf_a2d2_Zt_Sn," %% a2d2_Zt_Sn: ");}
    if (verbose>2){ lfprintf(E->lf_a3d2_Zt_Sn," %% a3d2_Zt_Sn: ");}
    /* if (QR_AnZtSZnAt_bother || QR_AnAtTAnAt_bother){ } */}
  if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnAtTAnAt_bother){
    if (verbose>1){ printf(" %% et_An, a1d1_At_en, a1d2_At_en, a3d2_At_en\n");}
    na_j=0; while (na_j<M_At->rpop_j){
      na_a = na_a_[na_j]; na_b = na_b_[na_j];
      wAt_tag = (__m128i *)(&(M_At->wX[na_b*M_At->mc_length]));
      mcax_tag = (__m128i *)(&(M_At->mc_j[0]));
      mcax_end = (__m128i *)(&(M_At->mc_j[M_At->mc_length]));
      if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ L1_plusequals(E->lf_et_An , na_j,na_b,na_a , -M_At->cpop_j + 2*popcount(&wAt_tag,&mcax_tag,&mcax_end));}
      if (QR_AnAtTAnAt_bother){ L1_plusequals(E->lf_a1d1_At_en , na_j,na_b,na_a , (A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_1_1*A_pcols] : 0)*(*L1_get(E->lf_et_An , na_j,na_b,na_a)));}
      if (QR_AnAtTAnAt_bother){ L1_plusequals(E->lf_a1d2_At_en , na_j,na_b,na_a , (A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_1_2*A_pcols] : 0)*(*L1_get(E->lf_et_An , na_j,na_b,na_a)));}
      if (QR_AnAtTAnAt_bother){ L1_plusequals(E->lf_a3d2_At_en , na_j,na_b,na_a , (A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_3_2*A_pcols] : 0)*(*L1_get(E->lf_et_An , na_j,na_b,na_a)));}
      GLOBAL_ops_count_one(tidx,3,0);
      na_j++; /* while (na_j<M_At->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_et_An," %% et_An: ");}
    if (verbose>2){ lfprintf(E->lf_a1d1_At_en," %% a1d1_At_en: ");}
    if (verbose>2){ lfprintf(E->lf_a1d2_At_en," %% a1d2_At_en: ");}
    if (verbose>2){ lfprintf(E->lf_a3d2_At_en," %% a3d2_At_en: ");}
    /* if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnAtTAnAt_bother){ } */}
  if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnZtSZnAt_bother || QR_AnAtTAnAt_bother){  
    if (verbose>1){ printf(" %% et_Tn, et_Sn\n");}
    ns_j=0; while (ns_j<M_Tt->rpop_j){
      ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
      if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnAtTAnAt_bother){
	wTt_tag = (__m128i *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
	mctx_tag = (__m128i *)(&(M_Tt->mc_j[0]));
	mctx_end = (__m128i *)(&(M_Tt->mc_j[M_Tt->mc_length]));
	L1_set(E->lf_et_Tn , ns_j,ns_b,ns_a , -M_Tt->cpop_j + 2*popcount(&wTt_tag,&mctx_tag,&mctx_end));
	GLOBAL_ops_count_one(tidx,2,M_Tt->mc_length*BIT8);
	/* if bother */}
      if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QR_AnZtSZnAt_bother){ 
	wSt_tag = (__m128i *)(&(M_St->wX[ns_b*M_St->mc_length]));
	mctx_tag = (__m128i *)(&(M_St->mc_j[0]));
	mctx_end = (__m128i *)(&(M_St->mc_j[M_St->mc_length]));
	L1_set(E->lf_et_Sn , ns_j,ns_b,ns_a , -M_St->cpop_j + 2*popcount(&wSt_tag,&mctx_tag,&mctx_end));
	GLOBAL_ops_count_one(tidx,2,M_St->mc_length*BIT8);
	/* if bother */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_et_Tn," %% et_Tn: ");}
    if (verbose>2){ lfprintf(E->lf_et_Sn," %% et_Sn: ");}
    /* if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnZtSZnAt_bother ||  QR_AnAtTAnAt_bother){ } */}
  if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){
    if (verbose>1){ printf(" %% et_Yn_a1d1, et_An_a1d1\n");}
    ma_j=0; while (ma_j<M_At->cpop_j){
      ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
      if (QC_AtTYnYtTAn_bother){ L1_plusequals(E->lf_et_Yn_a1d1 , 0,0,0 , (*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
      if (QC_AtTAnAtTAn_bother){ L1_plusequals(E->lf_et_An_a1d1 , 0,0,0 , (*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
      GLOBAL_ops_count_one(tidx,2,0);
      ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_et_Yn_a1d1," %% et_Yn_a1d1: ");}
    if (verbose>2){ lfprintf(E->lf_et_An_a1d1," %% et_An_a1d1: ");}
    /* if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother){
    if (verbose>1){ printf(" %% T_AnAt_YnYt, T_AnAt_AnAt\n");}  
    A_dtmp_a0d1=0; for (nc=0;nc<A_pcols;nc++){ A_dtmp_a0d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_0_1*A_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    A_dtmp_a2d1=0; for (nc=0;nc<A_pcols;nc++){ A_dtmp_a2d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_2_1*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    if (verbose>2){ printf(" %% A_dtmp_a0d1 %f A_dtmp_a2d1 %f\n",A_dtmp_a0d1,A_dtmp_a2d1);}
    Y_dtmp_a0d1=0; for (nc=0;nc<Y_pcols;nc++){ Y_dtmp_a0d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_0_1*Y_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_Yn->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    Y_dtmp_a2d1=0; for (nc=0;nc<Y_pcols;nc++){ Y_dtmp_a2d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_2_1*Y_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_Yn->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    if (verbose>2){ printf(" %% Y_dtmp_a0d1 %f Y_dtmp_a2d1 %f\n",Y_dtmp_a0d1,Y_dtmp_a2d1);}
    ma_j=0; while (ma_j<M_At->cpop_j){
      ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
      A_dtmp = A_dtmp_a0d1;
      wA1_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
      wA2_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
      mc_tag = (__m128i*)&(M_An->mc_j[0]); mc_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
      dinp = &(A_ajdk[0+AJDK_0_1*A_pcols]);
      A_dtmp -= 2*popcount_xor_lf(&wA1_tag,&wA2_tag,&mc_tag,&mc_end,&dinp);
      A_dtmp -= (*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)) + (*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)) - A_dtmp_a2d1;
      GLOBAL_ops_count_one(tidx,1,M_An->mc_length*BIT8);
      if (D->Y_cbother){
	my_j = ma_j; my_b = ma_b ; my_a = ma_a;
	Y_dtmp = Y_dtmp_a0d1;
	wY1_tag = (__m128i*)&(M_Yn->wX[my_b*M_Yn->mc_length]);
	wY2_tag = (__m128i*)&(M_Yn->wX[my_b*M_Yn->mc_length]);
	mc_tag = (__m128i*)&(M_Yn->mc_j[0]); mc_end = (__m128i*)&(M_Yn->mc_j[M_Yn->mc_length]);
	dinp = &(Y_ajdk[0+AJDK_0_1*Y_pcols]);
	Y_dtmp -= 2*popcount_xor_lf(&wY1_tag,&wY2_tag,&mc_tag,&mc_end,&dinp);
	Y_dtmp -= (*L2_get(lf_Yn_ajdk , my_j,my_b,my_a , 0,0,AJDK_1_1)) + (*L2_get(lf_Yn_ajdk , my_j,my_b,my_a , 0,0,AJDK_1_1)) - Y_dtmp_a2d1;
	GLOBAL_ops_count_one(tidx,1,M_Yn->mc_length*BIT8);
	/* if bother */}
      ns_j=0; while (ns_j<M_Tt->rpop_j){
	ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
	Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
	mct_tag = (unsigned char *)(&(M_Tt->mc_j[0]));
	mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length]));
	vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",ma_j,ma_b,ma_a);}
	if (QR_AnAtTAnAt_bother){ L2_set(E->lf_T_AnAt_AnAt , ns_j,ns_b,ns_a , ma_j,ma_b,ma_a , vT*A_dtmp*A_dtmp);}
	if (QR_AnAtTYnYt_bother){ L2_set(E->lf_T_AnAt_YnYt , ns_j,ns_b,ns_a , ma_j,ma_b,ma_a , vT*A_dtmp*Y_dtmp);}
	GLOBAL_ops_count_one(tidx,4,0);
	ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
      ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_T_AnAt_YnYt," %% T_AnAt_YnYt: ");}
    if (verbose>2){ lfprintf(E->lf_T_AnAt_AnAt," %% T_AnAt_AnAt: ");}
    /* if (QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother){ } */}
  if (verbose){ printf(" %% [finished get_singlestudy_uu] tidx %d\n",tidx);}
  return NULL;
}

void wrap_singlestudy_uu(int *tidx,void **vpra,pthread_t *thread_in,struct bcc_single *E)
{
  /* This function calls get_singlestudy_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 2)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int nb1=0,nbins=0;
  int ip=0;
  double **output_singlestudy_uu_p=NULL;
  if (verbose){ printf(" %% [entering wrap_singlestudy_uu] tidx %d\n",*tidx);}
  ip=0; 
  vpra[ip++] = tidx; vpra[ip++] = E; 
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_singlestudy_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_singlestudy_uu\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_singlestudy_uu(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_singlestudy_uu] tidx %d\n",*tidx);}
}

void bcc_singlestudy_uu(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  if (verbose){ printf(" %% [entering bcc_singlestudy_uu]\n");}
  if (verbose){ printf(" %% calculating singlestudy_uu.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    GLOBAL_pthread_tic();
    wrap_singlestudy_uu(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E);
    GLOBAL_pthread_toc();
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% singlestudy_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished bcc_singlestudy_uu]\n");}
}

/******************************************************************/

void *get_singlestudy_vv(void *vp)
{
  /* This function takes in a variety of inputs and calculates a variety of single-study terms */
  int verbose=0;
  int ip=0,length=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct bcc_single *E = (struct bcc_single *)(vpra[ip++]); struct bcc_ajdk *D=E->D; 
  struct M_handle *M_An = E->M_An;
  struct M_handle *M_At = E->M_At;
  struct M_handle *M_Tt = E->M_Tt;
  struct M_handle *M_Yn = E->M_Yn;
  struct M_handle *M_Yt = E->M_Yt;
  /* struct M_handle *M_Wt = E->M_Wt; */
  struct M_handle *M_St = E->M_St;
  struct M_handle *M_Zt = E->M_Zt;
  struct L_handle *lf_An_ajdk = E->lf_An_ajdk;
  struct L_handle *lf_Yn_ajdk = E->lf_Yn_ajdk;
  /* struct L_handle *lf_Wn_ajdk = E->lf_Wn_ajdk; */
  struct L_handle *lf_Zn_ajdk = E->lf_Zn_ajdk;
  double *A_ajdk = D->A_ajdk;
  double *Y_ajdk = D->Y_ajdk;
  int A_pcols = psize(M_At->nrows);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yt->nrows);
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  struct L_handle *lf_tmp=NULL;
  int A_ncols = M_At->nrows, T_ncols = M_Tt->nrows, A_nrows = M_At->ncols;
  int tic_num = GLOBAL_NTICKS-1;
  int mx_j=0,mx_chunk=0;
  /* length = A_ncols; lf_tmp = E->lf_At_Yn_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp); */
  /* length = A_ncols; lf_tmp = E->lf_At_An_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp); */
  length = T_ncols; lf_tmp = E->lf_a1d2_Zt_Sn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a2d2_Zt_Sn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a3d2_Zt_Sn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a1d2_At_Tn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a2d2_At_Tn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_a3d2_At_Tn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_et_Tn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_et_Sn; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = A_ncols; lf_tmp = E->lf_et_An; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = A_ncols; lf_tmp = E->lf_a1d1_At_en; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = A_ncols; lf_tmp = E->lf_a1d2_At_en; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = A_ncols; lf_tmp = E->lf_a3d2_At_en; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = 1; lf_tmp = E->lf_et_Yn_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_Tt_Yn_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  /* length = A_ncols*T_ncols; lf_tmp = E->lf_At_T_Yn_a1d1; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp); */
  length = 1; lf_tmp = E->lf_et_An_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  length = T_ncols; lf_tmp = E->lf_Tt_An_a1d1; lf_tmp->row_stride = length; lf_tmp->spacing_row = SPACING_a; L_zero(lf_tmp);
  /* length = A_ncols*T_ncols; lf_tmp = E->lf_At_T_An_a1d1; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp); */
  length = T_ncols*A_nrows; lf_tmp = E->lf_T_AnAt_YnYt; lf_tmp->row_stride = T_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = A_nrows; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  length = T_ncols*A_nrows; lf_tmp = E->lf_T_AnAt_AnAt; lf_tmp->row_stride = T_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = A_nrows; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  int ns_j=0,ns_b=0,ns_a=0;
  int na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0;
  int ny_j=0,ny_b=0,ny_a=0,my_j=0,my_b=0,my_a=0;
  /* int nw_j=0,nw_b=0,nw_a=0,mw_j=0,mw_b=0,mw_a=0; */
  int nz_j=0,nz_b=0,nz_a=0,mz_j=0,mz_b=0,mz_a=0;
  int vA=0,vY=0,vM=0,vT=0,vS=0;
  __m128i *wAt_tag=NULL;
  __m128i *wTt_tag=NULL;
  __m128i *wSt_tag=NULL;
  __m128i *mcax_tag=NULL,*mcax_end=NULL;
  __m128i *mctx_tag=NULL,*mctx_end=NULL;
  __m128i *wA1_tag=NULL; __m128i *wA2_tag=NULL; __m128i *wY1_tag=NULL; __m128i *wY2_tag=NULL; __m128i *mc_tag=NULL; __m128i *mc_end=NULL;
  int nc=0; double A_dtmp=0,A_dtmp_a0d1=0,A_dtmp_a2d1=0,Y_dtmp=0,Y_dtmp_a0d1=0,Y_dtmp_a2d1=0,*dinp=NULL;
  unsigned char *At_tag=NULL;
  unsigned char *An_tag=NULL;
  unsigned char *Tt_tag=NULL;
  unsigned char *Yn_tag=NULL;
  /* unsigned char *Yt_tag=NULL; */
  /* unsigned char *Wt_tag=NULL; */
  /* unsigned char *Zt_tag=NULL; */
  unsigned char *St_tag=NULL;
  unsigned char *mca_tag=NULL,*mca_end=NULL;
  unsigned char *mcy_tag=NULL,*mcy_end=NULL;
  unsigned char *mct_tag=NULL,*mct_end=NULL;
  unsigned int *ma_b_,*ma_a_;
  unsigned int *na_b_,*na_a_;
  unsigned int *my_b_,*my_a_;
  unsigned int *ny_b_,*ny_a_;
  /* unsigned int *mw_b_,*mw_; */
  /* unsigned int *nw_b_,*nw_; */
  unsigned int *mz_b_,*mz_a_;
  unsigned int *nz_b_,*nz_a_;
  unsigned int *ns_b_,*ns_a_;
  int QR_AnZtSWnYt_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother;
  int QR_AnZtSZnAt_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother;
  int QR_AnAtTYnYt_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother;
  int QR_AnAtTAnAt_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother;
  int QC_AtTYnWtSZn_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother;
  int QC_AtTAnZtSZn_bother = E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother;
  int QC_AtTYnYtTAn_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother;
  int QC_AtTAnAtTAn_bother = E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother;
  if (verbose>1){ printf(" %% [entering get_singlestudy_vv] tidx %d\n",tidx);}  
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Yt->mr_b,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_b: "); bprintf(M_Yt->mr_j,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_j: ");}
  if (verbose>2){ bprintf(M_Yt->mc_b,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_b: "); bprintf(M_Yt->mc_j,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_j: ");}
  /* if (verbose>2){ bprintf(M_Wt->mr_b,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_b: "); bprintf(M_Wt->mr_j,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_j: ");} */
  /* if (verbose>2){ bprintf(M_Wt->mc_b,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_b: "); bprintf(M_Wt->mc_j,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_j: ");} */
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  na_b_ = M_An->n_b_; na_a_ = M_An->n_a_;
  my_b_ = M_Yn->m_b_; my_a_ = M_Yn->m_a_;
  ny_b_ = M_Yn->n_b_; ny_a_ = M_Yn->n_a_;
  /* mw_b_ = M_Wt->m_b_; mw_ = M_Wt->m_a_; */
  /* nw_b_ = M_Wt->n_b_; nw_ = M_Wt->n_a_; */
  mz_b_ = M_Zt->n_b_; mz_a_ = M_Zt->n_a_;
  nz_b_ = M_Zt->m_b_; nz_a_ = M_Zt->m_a_;
  ns_b_ = M_St->m_b_; ns_a_ = M_St->m_a_;
  if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){
    if (verbose>1){ printf(" %% Tt_Yn_a1d1, Tt_An_a1d1\n");}
    GLOBAL_tic(tic_num);
    ns_j=0; while (ns_j<M_Tt->rpop_j){
      ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
      Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
      mct_tag = (unsigned char *)(&(M_Tt->mc_j[0]));
      mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length]));
      ma_j=0; while (ma_j<M_At->cpop_j){
	ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",ma_j,ma_b,ma_a);}
	if (QC_AtTYnYtTAn_bother){ L1_plusequals(E->lf_Tt_Yn_a1d1 , ns_j,ns_b,ns_a , vT*(*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
	if (QC_AtTAnAtTAn_bother){ L1_plusequals(E->lf_Tt_An_a1d1 , ns_j,ns_b,ns_a , vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
	GLOBAL_ops_count_one(tidx,4,0);
	ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_Tt_Yn_a1d1," %% Tt_Yn_a1d1: ");}
    if (verbose>2){ lfprintf(E->lf_Tt_An_a1d1," %% Tt_An_a1d1: ");}
    GLOBAL_toc(tic_num,0+verbose," %% Tt_Yn_a1d1: ");
/*     if (verbose>1){ printf(" %% At_T_Yn_a1d1, At_T_An_a1d1\n");} */
/*     GLOBAL_tic(tic_num); */
/*     if (GLOBAL_omp_type==GLOBAL_omp_off){  */
/*       ns_j=0; while (ns_j<M_Tt->rpop_j){ */
/* 	ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j]; */
/* 	Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length])); */
/* 	mct_tag = (unsigned char *)(&(M_Tt->mc_j[0])); */
/* 	mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length])); */
/* 	na_j=0; while (na_j<M_At->rpop_j){ */
/* 	  na_a = na_a_[na_j]; na_b = na_b_[na_j]; */
/* 	  At_tag = (unsigned char *)(&(M_At->wX[na_b*M_At->mc_length])); */
/* 	  mca_tag = (unsigned char *)(&(M_At->mc_j[0])); */
/* 	  mca_end = (unsigned char *)(&(M_At->mc_j[M_At->mc_length])); */
/* 	  ma_j=0; while (ma_j<M_At->cpop_j){ */
/* 	    ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j]; */
/* 	    vA = bget____(At_tag,ma_a); vM = bget__on(mca_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_At->mc_j\n",ma_j,ma_b,ma_a);} */
/* 	    vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_At->mc_j\n",ma_j,ma_b,ma_a);} */
/* 	    if (QC_AtTYnYtTAn_bother){ L2_plusequals(E->lf_At_T_Yn_a1d1 , na_j,na_b,na_a , ns_j,ns_b,ns_a , vA*vT*(*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
/* 	    if (QC_AtTAnAtTAn_bother){ L2_plusequals(E->lf_At_T_An_a1d1 , na_j,na_b,na_a , ns_j,ns_b,ns_a , vA*vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
/* 	    GLOBAL_ops_count_one(tidx,4,0); */
/* 	    ma_j++; /\* while (ma_j<M_At->cpop_j){ } *\/} */
/* 	  na_j++; /\* while (na_j<M_At->rpop_j){ } *\/} */
/* 	ns_j++; /\* while (ns_j<M_Tt->rpop_j){ } *\/} */
/*       /\* if (GLOBAL_omp_type==GLOBAL_omp_off){ } *\/} */
/*     else if (GLOBAL_omp_type==GLOBAL_omp__on){  */
/*       mx_chunk=128; */
/* #pragma omp parallel private(mx_j,ns_j,ns_a,ns_b,Tt_tag,na_j,na_a,na_b,At_tag,ma_j,ma_a,ma_b,vA,vT) */
/*       { /\* begin omp parallel *\/ */
/* 	mx_j=0; */
/* #pragma omp for schedule(static) */
/* 	for (mx_j=0;mx_j<M_Tt->rpop_j*M_At->rpop_j;mx_j++){ */
/* 	  ns_j = mx_j/M_At->rpop_j; na_j = mx_j % M_At->rpop_j; */
/* 	  ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j]; */
/* 	  Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length])); */
/* 	  na_a = na_a_[na_j]; na_b = na_b_[na_j]; */
/* 	  At_tag = (unsigned char *)(&(M_At->wX[na_b*M_At->mc_length])); */
/* 	  ma_j=0; while (ma_j<M_At->cpop_j){ */
/* 	    ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j]; */
/* 	    vA = bget____(At_tag,ma_a); vT = bget____(Tt_tag,ma_a); */
/* 	    if (QC_AtTYnYtTAn_bother){ L2_plusequals(E->lf_At_T_Yn_a1d1 , na_j,na_b,na_a , ns_j,ns_b,ns_a , vA*vT*(*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
/* 	    if (QC_AtTAnAtTAn_bother){ L2_plusequals(E->lf_At_T_An_a1d1 , na_j,na_b,na_a , ns_j,ns_b,ns_a , vA*vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
/* 	    ma_j++; /\* while (ma_j<M_At->cpop_j){ } *\/} */
/* 	  /\* for (mx_j=0;mx_j<M_Tt->rpop_j*M_At->rpop_j;mx_j++){ } *\/} */
/* 	/\* end omp parallel *\/} */
/*       GLOBAL_ops_count_one(tidx,4*M_Tt->rpop_j*M_At->rpop_j*M_At->cpop_j,0); */
/*       /\* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } *\/} */
/*     if (verbose>2){ lfprintf(E->lf_At_T_Yn_a1d1," %% At_T_Yn_a1d1: ");} */
/*     if (verbose>2){ lfprintf(E->lf_At_T_An_a1d1," %% At_T_An_a1d1: ");} */
/*     GLOBAL_toc(tic_num,0+verbose," %% At_T_Yn_a1d1: "); */
    /* if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){
/*     if (verbose>1){ printf(" %% At_Yn_a1d1, At_An_a1d1\n");} */
/*     GLOBAL_tic(tic_num); */
/*     if (GLOBAL_omp_type==GLOBAL_omp_off){  */
/*       na_j=0; while (na_j<M_At->rpop_j){ */
/* 	na_a = na_a_[na_j]; na_b = na_b_[na_j]; */
/* 	At_tag = (unsigned char *)(&(M_At->wX[na_b*M_At->mc_length])); */
/* 	mca_tag = (unsigned char *)(&(M_At->mc_j[0])); */
/* 	mca_end = (unsigned char *)(&(M_At->mc_j[M_At->mc_length])); */
/* 	ma_j=0; while (ma_j<M_At->cpop_j){ */
/* 	  ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j]; */
/* 	  vA = bget____(At_tag,ma_a); vM = bget__on(mca_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_At->mc_j\n",ma_j,ma_b,ma_a);} */
/* 	  if (QC_AtTYnYtTAn_bother){ L1_plusequals(E->lf_At_Yn_a1d1 , na_j,na_b,na_a , vA*(*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
/* 	  if (QC_AtTAnAtTAn_bother){ L1_plusequals(E->lf_At_An_a1d1 , na_j,na_b,na_a , vA*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
/* 	  GLOBAL_ops_count_one(tidx,2,0); */
/* 	  ma_j++; /\* while (ma_j<M_At->cpop_j){ } *\/} */
/* 	na_j++; /\* while (na_j<M_At->rpop_j){ } *\/} */
/*       /\* if (GLOBAL_omp_type==GLOBAL_omp_off){ } *\/} */
/*     else if (GLOBAL_omp_type==GLOBAL_omp__on){  */
/*       mx_chunk=128; */
/* #pragma omp parallel private(mx_j,na_j,na_a,na_b,At_tag,ma_j,ma_a,ma_b,vA) */
/*       { /\* begin omp parallel *\/ */
/* 	mx_j=0; */
/* #pragma omp for schedule(static) */
/* 	for (mx_j=0;mx_j<M_At->rpop_j;mx_j++){ */
/* 	  na_j = mx_j; */
/* 	  na_a = na_a_[na_j]; na_b = na_b_[na_j]; */
/* 	  At_tag = (unsigned char *)(&(M_At->wX[na_b*M_At->mc_length])); */
/* 	  ma_j=0; while (ma_j<M_At->cpop_j){ */
/* 	    ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j]; */
/* 	    vA = bget____(At_tag,ma_a); */
/* 	    if (QC_AtTYnYtTAn_bother){ L1_plusequals(E->lf_At_Yn_a1d1 , na_j,na_b,na_a , vA*(*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
/* 	    if (QC_AtTAnAtTAn_bother){ L1_plusequals(E->lf_At_An_a1d1 , na_j,na_b,na_a , vA*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));} */
/* 	    ma_j++; /\* while (ma_j<M_At->cpop_j){ } *\/}	   */
/* 	  /\* for (mx_j=0;mx_j<M_At->rpop_j;mx_j++){ } *\/} */
/* 	/\* end omp parallel *\/} */
/*       GLOBAL_ops_count_one(tidx,2*M_At->rpop_j*M_At->cpop_j,0); */
/*       /\* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } *\/} */
/*     if (verbose>2){ lfprintf(E->lf_At_Yn_a1d1," %% At_Yn_a1d1: ");} */
/*     if (verbose>2){ lfprintf(E->lf_At_An_a1d1," %% At_An_a1d1: ");} */
/*     GLOBAL_toc(tic_num,0+verbose," %% At_Yn_a1d1: "); */
    /* if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QR_AnZtSZnAt_bother || QR_AnAtTAnAt_bother){
    if (verbose>1){ printf(" %% a1d2_Zt_Sn, a2d2_Zt_Sn, a3d2_Zt_Sn, a1d2_At_Tn, a2d2_At_Tn, a3d2_At_Tn\n");}
    GLOBAL_tic(tic_num);
    ns_j=0; while (ns_j<M_Tt->rpop_j){
      ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
      Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
      mct_tag = (unsigned char *)(&(M_Tt->mc_j[0]));
      mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length]));
      if (QR_AnAtTAnAt_bother){ 
	ma_j=0; while (ma_j<M_At->cpop_j){
	  ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	  vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",ma_j,ma_b,ma_a);}
	  L1_plusequals(E->lf_a1d2_At_Tn, ns_j,ns_b,ns_a , vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_2)));
	  L1_plusequals(E->lf_a2d2_At_Tn, ns_j,ns_b,ns_a , vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_2_2)));
	  L1_plusequals(E->lf_a3d2_At_Tn, ns_j,ns_b,ns_a , vT*(*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_3_2)));
	  GLOBAL_ops_count_one(tidx,3,0);
	  ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
	/* if (QR_AnAtTAnAt_bother){ } */}
      St_tag = (unsigned char *)(&(M_St->wX[ns_b*M_St->mc_length]));
      mct_tag = (unsigned char *)(&(M_St->mc_j[0]));
      mct_end = (unsigned char *)(&(M_St->mc_j[M_St->mc_length]));
      if (QR_AnZtSZnAt_bother){
	mz_j=0; while (mz_j<M_Zt->cpop_j){
	  mz_a = mz_a_[mz_j]; mz_b = mz_b_[mz_j];
	  vS = bget____(St_tag,mz_a); vM = bget__on(mct_tag,mz_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_St->mc_j\n",mz_j,mz_b,mz_a);}
	  L1_plusequals(E->lf_a1d2_Zt_Sn, ns_j,ns_b,ns_a , vS*(*L2_get(lf_Zn_ajdk , mz_j,mz_b,mz_a , 0,0,AJDK_1_2)));
	  L1_plusequals(E->lf_a2d2_Zt_Sn, ns_j,ns_b,ns_a , vS*(*L2_get(lf_Zn_ajdk , mz_j,mz_b,mz_a , 0,0,AJDK_2_2)));
	  L1_plusequals(E->lf_a3d2_Zt_Sn, ns_j,ns_b,ns_a , vS*(*L2_get(lf_Zn_ajdk , mz_j,mz_b,mz_a , 0,0,AJDK_3_2)));
	  GLOBAL_ops_count_one(tidx,3,0);
	  mz_j++; /* while (mz_j<M_Tt->cpop_j){ } */}
	/* if (QR_AnZtSZnAt_bother){ } */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_a1d2_At_Tn," %% a1d2_At_Tn: ");}
    if (verbose>2){ lfprintf(E->lf_a2d2_At_Tn," %% a2d2_At_Tn: ");}
    if (verbose>2){ lfprintf(E->lf_a3d2_At_Tn," %% a3d2_At_Tn: ");}
    if (verbose>2){ lfprintf(E->lf_a1d2_Zt_Sn," %% a1d2_Zt_Sn: ");}
    if (verbose>2){ lfprintf(E->lf_a2d2_Zt_Sn," %% a2d2_Zt_Sn: ");}
    if (verbose>2){ lfprintf(E->lf_a3d2_Zt_Sn," %% a3d2_Zt_Sn: ");}
    GLOBAL_toc(tic_num,0+verbose," %% a1d2_Zt_Sn: ");
    /* if (QR_AnZtSZnAt_bother || QR_AnAtTAnAt_bother){ } */}
  if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnAtTAnAt_bother){
    if (verbose>1){ printf(" %% et_An, a1d1_At_en, a1d2_At_en, a3d2_At_en\n");}
    GLOBAL_tic(tic_num);
    na_j=0; while (na_j<M_At->rpop_j){
      na_a = na_a_[na_j]; na_b = na_b_[na_j];
      wAt_tag = (__m128i *)(&(M_At->wX[na_b*M_At->mc_length]));
      mcax_tag = (__m128i *)(&(M_At->mc_j[0]));
      mcax_end = (__m128i *)(&(M_At->mc_j[M_At->mc_length]));
      if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ L1_plusequals(E->lf_et_An , na_j,na_b,na_a , -M_At->cpop_j + 2*popcount(&wAt_tag,&mcax_tag,&mcax_end));}
      if (QR_AnAtTAnAt_bother){ L1_plusequals(E->lf_a1d1_At_en , na_j,na_b,na_a , (A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_1_1*A_pcols] : 0)*(*L1_get(E->lf_et_An , na_j,na_b,na_a)));}
      if (QR_AnAtTAnAt_bother){ L1_plusequals(E->lf_a1d2_At_en , na_j,na_b,na_a , (A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_1_2*A_pcols] : 0)*(*L1_get(E->lf_et_An , na_j,na_b,na_a)));}
      if (QR_AnAtTAnAt_bother){ L1_plusequals(E->lf_a3d2_At_en , na_j,na_b,na_a , (A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_3_2*A_pcols] : 0)*(*L1_get(E->lf_et_An , na_j,na_b,na_a)));}
      GLOBAL_ops_count_one(tidx,3,0);
      na_j++; /* while (na_j<M_At->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_et_An," %% et_An: ");}
    if (verbose>2){ lfprintf(E->lf_a1d1_At_en," %% a1d1_At_en: ");}
    if (verbose>2){ lfprintf(E->lf_a1d2_At_en," %% a1d2_At_en: ");}
    if (verbose>2){ lfprintf(E->lf_a3d2_At_en," %% a3d2_At_en: ");}
    GLOBAL_toc(tic_num,0+verbose," %% a1d1_At_en: ");
    /* if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnAtTAnAt_bother){ } */}
  if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnZtSZnAt_bother || QR_AnAtTAnAt_bother){  
    if (verbose>1){ printf(" %% et_Tn, et_Sn\n");}
    GLOBAL_tic(tic_num);
    ns_j=0; while (ns_j<M_Tt->rpop_j){
      ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
      if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnAtTAnAt_bother){
	wTt_tag = (__m128i *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
	mctx_tag = (__m128i *)(&(M_Tt->mc_j[0]));
	mctx_end = (__m128i *)(&(M_Tt->mc_j[M_Tt->mc_length]));
	L1_set(E->lf_et_Tn , ns_j,ns_b,ns_a , -M_Tt->cpop_j + 2*popcount(&wTt_tag,&mctx_tag,&mctx_end));
	GLOBAL_ops_count_one(tidx,2,M_Tt->mc_length*BIT8);
	/* if bother */}
      if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QR_AnZtSZnAt_bother){ 
	wSt_tag = (__m128i *)(&(M_St->wX[ns_b*M_St->mc_length]));
	mctx_tag = (__m128i *)(&(M_St->mc_j[0]));
	mctx_end = (__m128i *)(&(M_St->mc_j[M_St->mc_length]));
	L1_set(E->lf_et_Sn , ns_j,ns_b,ns_a , -M_St->cpop_j + 2*popcount(&wSt_tag,&mctx_tag,&mctx_end));
	GLOBAL_ops_count_one(tidx,2,M_St->mc_length*BIT8);
	/* if bother */}
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_et_Tn," %% et_Tn: ");}
    if (verbose>2){ lfprintf(E->lf_et_Sn," %% et_Sn: ");}
    GLOBAL_toc(tic_num,0+verbose," %% et_Tn: ");
    /* if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother || QR_AnZtSZnAt_bother ||  QR_AnAtTAnAt_bother){ } */}
  if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){
    if (verbose>1){ printf(" %% et_Yn_a1d1, et_An_a1d1\n");}
    GLOBAL_tic(tic_num);
    ma_j=0; while (ma_j<M_At->cpop_j){
      ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
      if (QC_AtTYnYtTAn_bother){ L1_plusequals(E->lf_et_Yn_a1d1 , 0,0,0 , (*L2_get(lf_Yn_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
      if (QC_AtTAnAtTAn_bother){ L1_plusequals(E->lf_et_An_a1d1 , 0,0,0 , (*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)));}
      GLOBAL_ops_count_one(tidx,2,0);
      ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_et_Yn_a1d1," %% et_Yn_a1d1: ");}
    if (verbose>2){ lfprintf(E->lf_et_An_a1d1," %% et_An_a1d1: ");}
    GLOBAL_toc(tic_num,0+verbose," %% et_Yn_a1d1: ");
    /* if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother){
    if (verbose>1){ printf(" %% T_AnAt_YnYt, T_AnAt_AnAt\n");}  
    GLOBAL_tic(tic_num);
    A_dtmp_a0d1=0; for (nc=0;nc<A_pcols;nc++){ A_dtmp_a0d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_0_1*A_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    A_dtmp_a2d1=0; for (nc=0;nc<A_pcols;nc++){ A_dtmp_a2d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_2_1*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    if (verbose>2){ printf(" %% A_dtmp_a0d1 %f A_dtmp_a2d1 %f\n",A_dtmp_a0d1,A_dtmp_a2d1);}
    Y_dtmp_a0d1=0; for (nc=0;nc<Y_pcols;nc++){ Y_dtmp_a0d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_0_1*Y_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_Yn->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    Y_dtmp_a2d1=0; for (nc=0;nc<Y_pcols;nc++){ Y_dtmp_a2d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_2_1*Y_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_Yn->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    if (verbose>2){ printf(" %% Y_dtmp_a0d1 %f Y_dtmp_a2d1 %f\n",Y_dtmp_a0d1,Y_dtmp_a2d1);}
    ma_j=0; while (ma_j<M_At->cpop_j){
      ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
      A_dtmp = A_dtmp_a0d1;
      wA1_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
      wA2_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
      mc_tag = (__m128i*)&(M_An->mc_j[0]); mc_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
      dinp = &(A_ajdk[0+AJDK_0_1*A_pcols]);
      A_dtmp -= 2*popcount_xor_lf(&wA1_tag,&wA2_tag,&mc_tag,&mc_end,&dinp);
      A_dtmp -= (*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)) + (*L2_get(lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_1)) - A_dtmp_a2d1;
      GLOBAL_ops_count_one(tidx,1,M_An->mc_length*BIT8);
      if (D->Y_cbother){
	my_j = ma_j; my_b = ma_b ; my_a = ma_a;
	Y_dtmp = Y_dtmp_a0d1;
	wY1_tag = (__m128i*)&(M_Yn->wX[my_b*M_Yn->mc_length]);
	wY2_tag = (__m128i*)&(M_Yn->wX[my_b*M_Yn->mc_length]);
	mc_tag = (__m128i*)&(M_Yn->mc_j[0]); mc_end = (__m128i*)&(M_Yn->mc_j[M_Yn->mc_length]);
	dinp = &(Y_ajdk[0+AJDK_0_1*Y_pcols]);
	Y_dtmp -= 2*popcount_xor_lf(&wY1_tag,&wY2_tag,&mc_tag,&mc_end,&dinp);
	Y_dtmp -= (*L2_get(lf_Yn_ajdk , my_j,my_b,my_a , 0,0,AJDK_1_1)) + (*L2_get(lf_Yn_ajdk , my_j,my_b,my_a , 0,0,AJDK_1_1)) - Y_dtmp_a2d1;
	GLOBAL_ops_count_one(tidx,1,M_Yn->mc_length*BIT8);
	/* if bother */}
      ns_j=0; while (ns_j<M_Tt->rpop_j){
	ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
	Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
	mct_tag = (unsigned char *)(&(M_Tt->mc_j[0]));
	mct_end = (unsigned char *)(&(M_Tt->mc_j[M_Tt->mc_length]));
	vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",ma_j,ma_b,ma_a);}
	if (QR_AnAtTAnAt_bother){ L2_set(E->lf_T_AnAt_AnAt , ns_j,ns_b,ns_a , ma_j,ma_b,ma_a , vT*A_dtmp*A_dtmp);}
	if (QR_AnAtTYnYt_bother){ L2_set(E->lf_T_AnAt_YnYt , ns_j,ns_b,ns_a , ma_j,ma_b,ma_a , vT*A_dtmp*Y_dtmp);}
	GLOBAL_ops_count_one(tidx,4,0);
	ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
      ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
    if (verbose>2){ lfprintf(E->lf_T_AnAt_YnYt," %% T_AnAt_YnYt: ");}
    if (verbose>2){ lfprintf(E->lf_T_AnAt_AnAt," %% T_AnAt_AnAt: ");}
    GLOBAL_toc(tic_num,0+verbose," %% T_AnAt_YnYt: ");
    /* if (QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother){ } */}
  if (verbose){ printf(" %% [finished get_singlestudy_vv] tidx %d\n",tidx);}
  return NULL;
}

void wrap_singlestudy_vv(int *tidx,void **vpra,pthread_t *thread_in,struct bcc_single *E)
{
  /* This function calls get_singlestudy_vv ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 2)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int nb1=0,nbins=0;
  int ip=0;
  double **output_singlestudy_vv_p=NULL;
  if (verbose){ printf(" %% [entering wrap_singlestudy_vv] tidx %d\n",*tidx);}
  ip=0; 
  vpra[ip++] = tidx; vpra[ip++] = E; 
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_singlestudy_vv,vpra)){ printf("Warning! cannot create thread %d in wrap_singlestudy_vv\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_singlestudy_vv(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_singlestudy_vv] tidx %d\n",*tidx);}
}

void At_Yn_a1d1_excerpt(struct M_handle *M_Yn,double *Y_ajdk,struct L_handle *lf_Yn_ajdk,struct L_handle *lf_Yn_a1d1)
{
  int verbose=0;
  int A_pcols=0,my_j=0,my_a=0,my_b=0,ns_j=0,ns_a=0,ns_b=0;
  if (verbose){ printf(" %% [entering At_Yn_a1d1_excerpt]\n");}
  lf_Yn_a1d1->spacing_row = SPACING_b;
  lf_Yn_a1d1->spacing_col = SPACING_b;
  lf_Yn_a1d1->row_stride = M_Yn->rpop_b;
  lf_Yn_a1d1->col_stride = 1;
  L_zero(lf_Yn_a1d1); 
  my_j=0;
  while (my_j<M_Yn->rpop_j){
    my_a = M_Yn->m_a_[my_j]; my_b = M_Yn->m_b_[my_j];
    L1_set(lf_Yn_a1d1,my_j,my_b,my_a, (*L2_get(lf_Yn_ajdk , my_j,my_b,my_a , 0,0,AJDK_1_1)));
    my_j++; /* while (my_j<M_Yn->rpop_j){ } */}
  if (verbose){ lfprintf(lf_Yn_ajdk," %% lf_Yn_ajdk: ");}
  if (verbose){ lfprintf(lf_Yn_a1d1," %% lf_Yn_a1d1: ");}
  if (verbose){ printf(" %% [finished At_Yn_a1d1_excerpt]\n");}
}

void bcc_singlestudy_ww(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb1=0; struct bcc_single *E=NULL;
  int n_spacing_B = SPACING_a;
  if (verbose){ printf(" %% [entering bcc_singlestudy_ww]\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0);
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    if (E->A_rbother && D->Y_cbother){ At_Yn_a1d1_excerpt(E->M_Yn,D->Y_ajdk,E->lf_Yn_ajdk,E->lf_Yn_a1d1);}
    if (E->A_rbother && D->A_cbother){ At_Yn_a1d1_excerpt(E->M_An,D->A_ajdk,E->lf_An_ajdk,E->lf_An_a1d1);}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  GLOBAL_toc(0,verbose," %% Yn_a1d1 An_a1d1: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    if (E->A_rbother && D->Y_cbother){
      GLOBAL_pthread_tic(); 
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->M_Yn->mr_j,E->M_Yn->mr_b,E->M_1->mr_j,E->M_1->mr_b,E->lf_Yn_a1d1,E->lf_Yn_a1d1->lf,&(E->M_Yn->nrows),&(addressable_1),&(E->M_Yn_a1d1),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->M_An->mr_j,E->M_An->mr_b,E->M_1->mr_j,E->M_1->mr_b,E->lf_An_a1d1,E->lf_An_a1d1->lf,&(E->M_An->nrows),&(addressable_1),&(E->M_An_a1d1),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    if (E->A_rbother && D->Y_cbother){
      GLOBAL_pthread_tic();
      wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,E->M_At,E->lf_Yn_a1d1,E->M_Yn_a1d1,E->M_1,&(E->lf_At_Yn_a1d1));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,E->M_At,E->lf_An_a1d1,E->M_An_a1d1,E->M_1,&(E->lf_At_An_a1d1));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% E->lf_At_Yn_a1d1 E->lf_At_An_a1d1: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    GLOBAL_pthread_tic();
    wrap_singlestudy_vv(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E);
    GLOBAL_pthread_toc();
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% singlestudy_vv: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  if (verbose){ printf(" %% [finished bcc_singlestudy_ww]\n");}
}

/******************************************************************/

void *get_doublestudy_uu(void *vp)
{
  /* This function takes in a variety of inputs and calculates a variety of double-study terns */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct bcc_double *F = (struct bcc_double *)(vpra[ip++]); struct bcc_ajdk *D=F->D; struct bcc_single *E_nb1 = F->E_nb1, *E_nb2 = F->E_nb2;
  struct M_handle *M_At_nb1 = E_nb1->M_At;
  struct M_handle *M_At_nb2 = E_nb2->M_At;
  struct M_handle *M_Tt_nb2 = E_nb2->M_Tt;
  struct M_handle *M_St_nb2 = E_nb2->M_St;
  struct M_handle *M_Zt_nb2 = E_nb2->M_Zt;
  struct L_handle *lf_ZtSn_nb2 = E_nb2->lf_ZtSn;
  struct L_handle *lf_AtTn_nb2 = E_nb2->lf_AtTn;
  int A_pcols = psize(M_At_nb1->nrows);
  double *A_ajdk = D->A_ajdk;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  struct L_handle *lf_tmp=NULL;
  lf_tmp = F->lf_An_a2d2_Zt_Sn; lf_tmp->row_stride = F->E_nb1->A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = D->T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  lf_tmp = F->lf_An_a2d2_At_Tn; lf_tmp->row_stride = F->E_nb1->A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = D->T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
  int ns_j=0,ns_b=0,ns_a=0;
  int na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0;
  int nz_j=0,nz_b=0,nz_a=0;
  double ZtSn_dtmp_a0d0=0,ZtSn_dtmp=0,AtTn_dtmp_a0d0=0,AtTn_dtmp=0;
  int vA=0,vM=0;
  __m128i *wS_tag=NULL;
  __m128i *wZ_tag=NULL;
  __m128i *wA_tag=NULL;
  __m128i *wT_tag=NULL;
  __m128i *wc_tag=NULL;
  __m128i *wc_end=NULL;
  unsigned char *At_tag=NULL;
  unsigned char *mca_tag=NULL,*mca_end=NULL;
  int mx_j=0,mx_chunk=0;
  int QR_AnZtSWnYt_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QR_AnZtSZnAt_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QR_AnAtTYnYt_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QR_AnAtTAnAt_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  int QC_AtTYnWtSZn_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QC_AtTAnZtSZn_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QC_AtTYnYtTAn_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QC_AtTAnAtTAn_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  if (verbose>1){ printf(" %% [entering get_doublestudy_uu] tidx %d\n",tidx);}
  if (verbose>2){ bprintf(M_At_nb1->mr_b,M_At_nb1->bitj,1,M_At_nb1->nrows," %% M_At_nb1->mr_b: "); bprintf(M_At_nb1->mr_j,M_At_nb1->bitj,1,M_At_nb1->nrows," %% M_At_nb1->mr_j: ");}
  if (verbose>2){ bprintf(M_At_nb1->mc_b,M_At_nb1->bitj,1,M_At_nb1->ncols," %% M_At_nb1->mc_b: "); bprintf(M_At_nb1->mc_j,M_At_nb1->bitj,1,M_At_nb1->ncols," %% M_At_nb1->mc_j: ");}
  if (verbose>2){ bprintf(M_At_nb2->mr_b,M_At_nb2->bitj,1,M_At_nb2->nrows," %% M_At_nb2->mr_b: "); bprintf(M_At_nb2->mr_j,M_At_nb2->bitj,1,M_At_nb2->nrows," %% M_At_nb2->mr_j: ");}
  if (verbose>2){ bprintf(M_At_nb2->mc_b,M_At_nb2->bitj,1,M_At_nb2->ncols," %% M_At_nb2->mc_b: "); bprintf(M_At_nb2->mc_j,M_At_nb2->bitj,1,M_At_nb2->ncols," %% M_At_nb2->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt_nb2->mr_b,M_Zt_nb2->bitj,1,M_Zt_nb2->nrows," %% M_Zt_nb2->mr_b: "); bprintf(M_Zt_nb2->mr_j,M_Zt_nb2->bitj,1,M_Zt_nb2->nrows," %% M_Zt_nb2->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt_nb2->mc_b,M_Zt_nb2->bitj,1,M_Zt_nb2->ncols," %% M_Zt_nb2->mc_b: "); bprintf(M_Zt_nb2->mc_j,M_Zt_nb2->bitj,1,M_Zt_nb2->ncols," %% M_Zt_nb2->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt_nb2->mr_b,M_Tt_nb2->bitj,1,M_Tt_nb2->nrows," %% M_Tt_nb2->mr_b: "); bprintf(M_Tt_nb2->mr_j,M_Tt_nb2->bitj,1,M_Tt_nb2->nrows," %% M_Tt_nb2->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt_nb2->mc_b,M_Tt_nb2->bitj,1,M_Tt_nb2->ncols," %% M_Tt_nb2->mc_b: "); bprintf(M_Tt_nb2->mc_j,M_Tt_nb2->bitj,1,M_Tt_nb2->ncols," %% M_Tt_nb2->mc_j: ");}
  if (verbose>2){ bprintf(M_St_nb2->mr_b,M_St_nb2->bitj,1,M_St_nb2->nrows," %% M_St_nb2->mr_b: "); bprintf(M_St_nb2->mr_j,M_St_nb2->bitj,1,M_St_nb2->nrows," %% M_St_nb2->mr_j: ");}
  if (verbose>2){ bprintf(M_St_nb2->mc_b,M_St_nb2->bitj,1,M_St_nb2->ncols," %% M_St_nb2->mc_b: "); bprintf(M_St_nb2->mc_j,M_St_nb2->bitj,1,M_St_nb2->ncols," %% M_St_nb2->mc_j: ");}  
  if (GLOBAL_omp_type==GLOBAL_omp_off){ 
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_St_nb2->m_a_[ns_j]; ns_b = M_St_nb2->m_b_[ns_j];
      nz_j=0; na_j=0; while (nz_j<M_Zt_nb2->rpop_j && na_j<M_At_nb2->rpop_j){
	if (QR_AnZtSZnAt_bother){ nz_a = M_Zt_nb2->m_a_[nz_j]; nz_b = M_Zt_nb2->m_b_[nz_j];}
	if (QR_AnAtTAnAt_bother){ na_a = M_At_nb2->m_a_[na_j]; na_b = M_At_nb2->m_b_[na_j];}
	At_tag = (unsigned char *)(&(M_At_nb1->wX[na_b*M_At_nb1->mc_length]));
	mca_tag = (unsigned char *)(&(M_At_nb1->mc_j[0]));
	mca_end = (unsigned char *)(&(M_At_nb1->mc_j[M_At_nb1->mc_length]));
	ma_j=0; while (ma_j<M_At_nb1->cpop_j){
	  ma_a = M_At_nb1->n_a_[ma_j]; ma_b = M_At_nb1->n_b_[ma_j];
	  vA = bget____(At_tag,ma_a); vM = bget__on(mca_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_At->mc_j\n",ma_j,ma_b,ma_a);}
	  if (QR_AnZtSZnAt_bother){ L2_plusequals(F->lf_An_a2d2_Zt_Sn , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , vA*A_ajdk[nz_a/POPLENGTH + AJDK_2_2*A_pcols]*(*L2_get(lf_ZtSn_nb2 , nz_j,nz_b,nz_a , ns_j,ns_b,ns_a)));}
	  if (QR_AnAtTAnAt_bother){ L2_plusequals(F->lf_An_a2d2_At_Tn , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , vA*A_ajdk[na_a/POPLENGTH + AJDK_2_2*A_pcols]*(*L2_get(lf_AtTn_nb2 , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  GLOBAL_ops_count_one(tidx,4,0);
	  ma_j++; /* while (ma_j<M_At_nb1->cpop_j){ } */}
	nz_j++; na_j++; /* while (nz_j<M_Zt_nb2->rpop_j && na_j<M_At_nb2->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp__on){ 
    mx_chunk=128; 
#pragma omp parallel private(mx_j,ns_j,ns_a,ns_b,nz_j,nz_a,nz_b,na_j,na_a,na_b,At_tag,ma_j,ma_a,ma_b,vA)
    { /* begin omp parallel */
      mx_j=0;
//#pragma omp for schedule(dynamic,mx_chunk)
#pragma omp for schedule(static)
      for (mx_j=0;mx_j<M_Tt_nb2->rpop_j*M_At_nb1->cpop_j;mx_j++){
	ns_j = mx_j/M_At_nb1->cpop_j; ma_j = mx_j % M_At_nb1->cpop_j;
	ns_a = M_St_nb2->m_a_[ns_j]; ns_b = M_St_nb2->m_b_[ns_j];
	ma_a = M_At_nb1->n_a_[ma_j]; ma_b = M_At_nb1->n_b_[ma_j];
	nz_j=0; na_j=0; while (nz_j<M_Zt_nb2->rpop_j && na_j<M_At_nb2->rpop_j){
	  if (QR_AnZtSZnAt_bother){ nz_a = M_Zt_nb2->m_a_[nz_j]; nz_b = M_Zt_nb2->m_b_[nz_j];}
	  if (QR_AnAtTAnAt_bother){ na_a = M_At_nb2->m_a_[na_j]; na_b = M_At_nb2->m_b_[na_j];}
	  At_tag = (unsigned char *)(&(M_At_nb1->wX[na_b*M_At_nb1->mc_length]));
	  vA = bget____(At_tag,ma_a);
	  if (QR_AnZtSZnAt_bother){ L2_plusequals(F->lf_An_a2d2_Zt_Sn , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , vA*A_ajdk[nz_a/POPLENGTH + AJDK_2_2*A_pcols]*(*L2_get(lf_ZtSn_nb2 , nz_j,nz_b,nz_a , ns_j,ns_b,ns_a)));}
	  if (QR_AnAtTAnAt_bother){ L2_plusequals(F->lf_An_a2d2_At_Tn , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , vA*A_ajdk[na_a/POPLENGTH + AJDK_2_2*A_pcols]*(*L2_get(lf_AtTn_nb2 , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  nz_j++; na_j++; /* while (nz_j<M_Zt_nb2->rpop_j && na_j<M_At_nb2->rpop_j){ } */}
	/* for (mx_j=0;mx_j<M_Tt_nb2->rpop_j*M_At_nb1->cpop_j;mx_j++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,4*M_Tt_nb2->rpop_j*M_At_nb2->rpop_j*M_At_nb1->cpop_j,0);
    /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
  if (verbose>1){ lfprintf(F->lf_An_a2d2_Zt_Sn," %% An_a2d2_Zt_Sn: ");}
  if (verbose>1){ lfprintf(F->lf_An_a2d2_At_Tn," %% An_a2d2_At_Tn: ");}
  if (verbose){ printf(" %% [finished get_doublestudy_uu] tidx %d\n",tidx);}
  return NULL;
}

void wrap_doublestudy_uu(int *tidx,void **vpra,pthread_t *thread_in,struct bcc_double *F)
{
  /* This function calls get_doublestudy_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 2)
   */
  int verbose=0;
  int ip=0;
  if (verbose){ printf(" %% [entering wrap_doublestudy_uu] tidx %d\n",*tidx);}
  ip=0;
  vpra[ip++] = tidx; vpra[ip++] = F;
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_doublestudy_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_doublestudy_uu\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_doublestudy_uu(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_doublestudy_uu] tidx %d\n",*tidx);}
}

void bcc_doublestudy_uu(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_double **F_ = D->F_; 
  int nbx=0,nb1=0,nb2=0; struct bcc_double *F=NULL;
  if (verbose){ printf(" %% [entering bcc_doublestudy_uu]\n");}
  if (verbose){ printf(" %% calculating doublestudy_uu.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx];
    GLOBAL_pthread_tic();
    wrap_doublestudy_uu(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),F);
    GLOBAL_pthread_toc();
    /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% doublestudy_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished bcc_doublestudy_uu]\n");}
}

void bcc_doublestudy_ww(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_; 
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E=NULL,*E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  int n_spacing_B = SPACING_a;
  struct L_handle *lf_tmp=NULL;
  if (verbose){ printf(" %% [entering bcc_doublestudy_ww]\n");}
  if (verbose){ printf(" %% calculating doublestudy_ww.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0);
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    if (E->Z_rbother && D->A_cbother){ An_a2d2_ZtSn_excerpt(E->M_Zt,E->M_St,D->A_ajdk,E->lf_ZtSn,E->lf_a2d2_ZtSn);}
    if (E->A_rbother && D->A_cbother){ An_a2d2_ZtSn_excerpt(E->M_At,E->M_Tt,D->A_ajdk,E->lf_AtTn,E->lf_a2d2_AtTn);}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  GLOBAL_toc(0,verbose," %% a2d2_ZtSn a2d2_AtTn: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    if (E->Z_rbother && D->A_cbother){
      GLOBAL_pthread_tic(); 
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->M_Zt->mr_j,E->M_Zt->mr_b,E->M_St->mr_j,E->M_St->mr_b,E->lf_a2d2_ZtSn,E->lf_a2d2_ZtSn->lf,&(E->M_Zt->nrows),&(E->M_St->nrows),&(E->M_a2d2_ZtSn),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->M_At->mr_j,E->M_At->mr_b,E->M_Tt->mr_j,E->M_Tt->mr_b,E->lf_a2d2_AtTn,E->lf_a2d2_AtTn->lf,&(E->M_At->nrows),&(E->M_Tt->nrows),&(E->M_a2d2_AtTn),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; E_nb1 = E_[nb1]; E_nb2 = E_[nb2]; F = F_[nbx];
      lf_tmp = F->lf_An_a2d2_Zt_Sn; lf_tmp->row_stride = F->E_nb1->A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = D->T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
      lf_tmp = F->lf_An_a2d2_At_Tn; lf_tmp->row_stride = F->E_nb1->A_nrows; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = D->T_ncols; lf_tmp->spacing_col = SPACING_a; L_zero(lf_tmp);
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,E_nb1->M_An,E_nb2->lf_a2d2_ZtSn,E_nb2->M_a2d2_ZtSn,E_nb2->M_St,&(F->lf_An_a2d2_Zt_Sn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,E_nb1->M_An,E_nb2->lf_a2d2_AtTn,E_nb2->M_a2d2_AtTn,E_nb2->M_Tt,&(F->lf_An_a2d2_At_Tn));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% F->lf_An_a2d2_Zt_Sn F->lf_An_a2d2_At_Tn: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished bcc_doublestudy_ww]\n");}
}

/******************************************************************/

void *get_QR_AnZtSWnYt_uu(void *vp)
{
  /* This function takes in M_An, M_Zt, M_St, M_Wt, M_Yn and calculates the output_QR_AnZtSWnYt: ;
     output_QR_AnZtSWnYt[mj+ns*A_n_rowsn] = (An(mj,:) - e_An*a_At)*D*(Zt(:,:)-a_An*e_At)*diag(S)*(Wn(:,:)-e_Yn*a_Yt)*D*(Yt(:,mj) - a_Yt*e_Yn). ;
     For this calculation we remove (i.e., skip) any flattened loops, using AZ_flag and AY_flag to determine which loops to skip. ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Wt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Yn = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_QR_AnZtSWnYt = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); 
  int output_spacing_s = *(int *)(vpra[ip++]);
  int AZ_flag = *(int *)(vpra[ip++]);
  int AY_flag = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols)/* rup(M_An->ncols + M_An->ncols_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yn->ncols)/* rup(M_Yn->ncols + M_Yn->ncols_extend,POPLENGTH)/POPLENGTH */;
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int mx_j=0,mx_chunk=0,ns_j=0,ns_b=0,ns_a=0,tab_s_stride=0,tab_s=0,na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0,tab_a_stride=0,tab_a=0,ny_j=0,ny_b=0,ny_a=0,my_j=0,my_b=0,my_a=0;
  int vA=0,vY=0;
  unsigned int *ma_b_,*ma_a_;
  unsigned int *na_b_,*na_a_;
  unsigned int *my_b_,*my_a_;
  unsigned int *ny_b_,*ny_a_;
  unsigned int *ns_b_,*ns_a_;
  int ms_j=0,ms_b=0,ms_a=0,mz_j=0,mz_b=0,mz_a=0,mw_j=0,mw_b=0,mw_a=0;
  int ncols_Z_p = psize(M_Zt->nrows);
  double *D_Zn = (double *)&(A_ajdk[0 + AJDK_0_1*ncols_Z_p]);
  double *a_Zn = (double *)&(A_ajdk[0 + AJDK_1_0*ncols_Z_p]);
  int ncols_W_p = psize(M_Wt->nrows);
  double *D_Wn = (double *)&(Y_ajdk[0 + AJDK_0_1*ncols_W_p]);
  double *a_Wn = (double *)&(Y_ajdk[0 + AJDK_1_0*ncols_W_p]);
  int vZ=0,vS=0,vW=0;
  unsigned char *wAn_tag=NULL,*wZt_tag=NULL,*wSt_tag=NULL,*wWt_tag=NULL,*wYn_tag=NULL;
  int rtn_flag=0;
  if (verbose>1){ printf(" %% [entering get_QR_AnZtSWnYt_uu] tidx %d \n",tidx);}
  if (verbose>2){ raprintf(D_An,"double",1,A_pcols," %% D_An: "); raprintf(a_An,"double",1,A_pcols," %% a_An: ");}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  if (verbose>2){ raprintf(D_Yn,"double",1,Y_pcols," %% D_Yn: "); raprintf(a_Yn,"double",1,Y_pcols," %% a_Yn: ");}
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  na_b_ = M_An->n_b_; na_a_ = M_An->n_a_;
  my_b_ = M_Yn->m_b_; my_a_ = M_Yn->m_a_;
  ny_b_ = M_Yn->n_b_; ny_a_ = M_Yn->n_a_;
  ns_b_ = M_St->m_b_; ns_a_ = M_St->m_a_;
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  output_QR_AnZtSWnYt->spacing_row = output_spacing_a; output_QR_AnZtSWnYt->row_stride = tab_a_stride;
  output_QR_AnZtSWnYt->spacing_col = output_spacing_s; output_QR_AnZtSWnYt->col_stride = tab_s_stride;
  fill_uchar_zero((unsigned char *)(output_QR_AnZtSWnYt->lf),tab_a_stride*tab_s_stride*sizeof(double));
  ns_j=0;
  while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    wSt_tag = (unsigned char *)(&(M_St->wX[ns_b*M_St->mc_length]));
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    ma_j=0;my_j=0;
    while (ma_j<M_An->rpop_j && my_j<M_Yn->rpop_j){
      ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j]; my_a = my_a_[my_j]; my_b = my_b_[my_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=ma_j+tab_s*tab_a_stride; break; case SPACING_b: tab_a=ma_b+tab_s*tab_a_stride; break; case SPACING_a: tab_a=ma_a+tab_s*tab_a_stride; break; default: break; /* switch (output_spacing_a){ } */}
      wAn_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length])); wYn_tag = (unsigned char *)(&(M_Yn->wX[my_b*M_Yn->mc_length]));
      na_j=0; 
      while (na_j<M_An->cpop_j){
	na_a = na_a_[na_j]; na_b = na_b_[na_j];
	vA = bget____(wAn_tag,na_a);
	wZt_tag = (unsigned char *)(&(M_Zt->wX[na_b*M_Zt->mc_length]));
	ny_j=0; 
	while (ny_j<M_Yn->cpop_j){
	  ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j];
	  vY = bget____(wYn_tag,ny_a);
	  wWt_tag = (unsigned char *)(&(M_Wt->wX[ny_b*M_Wt->mc_length]));
	  mz_j=0; ms_j=0; mw_j=0;
	  while (mz_j<M_Zt->cpop_j && ms_j<M_St->cpop_j && mw_j<M_Wt->cpop_j){
	    mz_a = M_Zt->n_a_[mz_j]; mz_b = M_Zt->n_b_[mz_j];
	    ms_a = M_St->n_a_[ms_j]; ms_b = M_St->n_b_[ms_j];
	    mw_a = M_Wt->n_a_[mw_j]; mw_b = M_Wt->n_b_[mw_j];
	    vZ = bget____(wZt_tag,ms_a); vS = bget____(wSt_tag,ms_a); vW = bget____(wWt_tag,ms_a);
	    rtn_flag = 1; if (AZ_flag && ma_j==mz_j){ rtn_flag *= 0;} if (AY_flag && na_j==ny_j){ rtn_flag *= 0;}
	    output_QR_AnZtSWnYt->lf[tab_a] += rtn_flag*(vA - a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*(vZ - a_Zn[na_a/POPLENGTH])*vS*(vW - a_Wn[ny_a/POPLENGTH])*D_Yn[ny_a/POPLENGTH]*(vY - a_Yn[ny_a/POPLENGTH]);
	    mz_j++; ms_j++; mw_j++; /* while (mz_j<M_Zt->cpop_j){ } */}
	  ny_j++; /* while (ny_j<M_Yn->cpop_j){ } */}
	na_j++; /* while (na_j<M_An->cpop_j){ } */}
      ma_j++;my_j++;}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_An->rpop_j*M_An->cpop_j*M_Yn->cpop_j*M_Zt->cpop_j*2,0);
  if (verbose>1){ printf(" %% [finished get_QR_AnZtSWnYt_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_QR_AnZtSWnYt_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,int AZ_flag,int AY_flag,struct M_handle *M_An,struct M_handle *M_Zt,struct M_handle *M_St,struct M_handle *M_Wt,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle **output_p)
{
  /* This function calls get_QR_AnZtSWnYt_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 14)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_QR_AnZtSWnYt_uu__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Yn,verbose," %% M_Yn: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  length = length_a*length_s; if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_QR_AnZtSWnYt_uu__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_Zt; vpra[ip++] = M_St; vpra[ip++] = M_Wt; vpra[ip++] = M_Yn; vpra[ip++] = A_ajdk; vpra[ip++] = Y_ajdk; vpra[ip++] = *output_p;
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  switch (AZ_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break;}
  switch (AY_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break;}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_QR_AnZtSWnYt_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_QR_AnZtSWnYt_uu__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_QR_AnZtSWnYt_uu(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_QR_AnZtSWnYt_uu__run] tidx %d\n",*tidx);}
  return length;
}  

void bcc_QR_AnZtSWnYt_uu(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  int n_spacing_A = SPACING_a; int n_spacing_B = SPACING_b;
  int AZ_flag=0,AY_flag=0;
  if (verbose){ printf(" %% [entering bcc_QR_AnZtSWnYt_uu]\n");}
  if (verbose){ printf(" %% calculating F_[nbx]->QR_AnZtSWnYt_uu, etc.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	AZ_flag=0; AY_flag=0;
	GLOBAL_pthread_tic(); 
	wrap_QR_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_B,AZ_flag,AY_flag,E_nb1->M_An,E_nb2->M_Zt,E_nb2->M_St,E_nb2->M_Wt,E_nb1->M_Yn,D->A_ajdk,D->Y_ajdk,&(F->QR_AnZtSWnYt_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	AZ_flag=0; AY_flag=1;
	GLOBAL_pthread_tic(); 
	wrap_QR_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_B,AZ_flag,AY_flag,E_nb1->M_An,E_nb2->M_Zt,E_nb2->M_St,E_nb2->M_Zt,E_nb1->M_An,D->A_ajdk,D->A_ajdk,&(F->QR_AnZtSZnAt_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	AZ_flag=(nb1==nb2); AY_flag=0;
	GLOBAL_pthread_tic(); 
	wrap_QR_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_B,AZ_flag,AY_flag,E_nb1->M_An,E_nb2->M_At,E_nb2->M_Tt,E_nb2->M_Yt,E_nb1->M_Yn,D->A_ajdk,D->Y_ajdk,&(F->QR_AnAtTYnYt_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	AZ_flag=(nb1==nb2); AY_flag=1;
	GLOBAL_pthread_tic(); 
	wrap_QR_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_B,AZ_flag,AY_flag,E_nb1->M_An,E_nb2->M_At,E_nb2->M_Tt,E_nb2->M_At,E_nb1->M_An,D->A_ajdk,D->A_ajdk,&(F->QR_AnAtTAnAt_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% QR_AnZtSWnYt_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	bprintf(E_nb2->M_Sn->mc_j,D->bitj,1,D->T_ncols," %% S2_bmc_j: ");
	lfprintf(F->QR_AnZtSWnYt_uu," %% QR_AnZtSWnYt_uu: ");
	lfprintf(F->QR_AnZtSZnAt_uu," %% QR_AnZtSZnAt_uu: ");
	bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	lfprintf(F->QR_AnAtTYnYt_uu," %% QR_AnAtTYnYt_uu: ");
	lfprintf(F->QR_AnAtTAnAt_uu," %% QR_AnAtTAnAt_uu: ");
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished bcc_QR_AnZtSWnYt_uu]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void *get_QC_AtTYnWtSZn_uu(void *vp)
{
  /* This function calculates :
     output_QC_AtTYnWtSZn_uu[mj+ns*A_n_cols] = daig(D_An)*(At(mj,:)-a_An(mj)*e_At)*diag(T(:,ns))*(Yn(:,:)-e_An*a_Yt)*diag(D_Yn)*(Wt(:,:)-a_Yn*e_Zt)*diag(S(:,ns))*(Zn(:,mj)-e_Zn*a_At(mj));
     For this calculation we remove (i.e., skip) any flattened loops, using AZ_flag and AY_flag to determine which loops to skip. ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Yn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Wn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_QC_AtTYnWtSZn_uu = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_s = *(int *)(vpra[ip++]);
  int AZ_flag = *(int *)(vpra[ip++]);
  int AY_flag = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yn->ncols);
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int ns_j=0,ns_b=0,ns_a=0,tab_s_stride=0,tab_s=0;
  int na_j=0,na_b=0,na_a=0,tab_a_stride=0,tab_a=0,my_j=0,my_b=0,my_a=0,nz_j=0,nz_b=0,nz_a=0,mw_j=0,mw_b=0,mw_a=0,tab_x=0;
  int ny_j=0,ny_b=0,ny_a=0,nw_j=0,nw_b=0,nw_a=0;
  unsigned char *At_tag=NULL;
  unsigned char *Tt_tag=NULL;
  unsigned char *Yn_tag=NULL;
  unsigned char *Wn_tag=NULL;
  unsigned char *St_tag=NULL;
  unsigned char *Zt_tag=NULL;
  int vA=0,vT=0,vY=0,vW=0,vS=0,vZ=0;
  int rtn_flag=0;
  if (verbose>1){ printf(" %% [entering get_QC_AtTYnWtSZn_uu] tidx %d\n",tidx);}  
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  if (verbose>2){ bprintf(M_Wn->mr_b,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_b: "); bprintf(M_Wn->mr_j,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_j: ");}
  if (verbose>2){ bprintf(M_Wn->mc_b,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_b: "); bprintf(M_Wn->mc_j,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  if (verbose>2){ printf(" %% tab_s_stride %d tab_a_stride %d\n",tab_s_stride,tab_a_stride);}
  output_QC_AtTYnWtSZn_uu->spacing_row = output_spacing_a; output_QC_AtTYnWtSZn_uu->row_stride = tab_a_stride; 
  output_QC_AtTYnWtSZn_uu->spacing_col = output_spacing_s; output_QC_AtTYnWtSZn_uu->col_stride = tab_s_stride; 
  ns_j=0;
  while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    Tt_tag = (unsigned char *)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
    St_tag = (unsigned char *)&(M_St->wX[ns_b*M_St->mc_length]);
    na_j=0; nz_j=0;
    while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
      na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
      At_tag = (unsigned char *)&(M_At->wX[na_b*M_At->mc_length]);
      Zt_tag = (unsigned char *)&(M_Zt->wX[nz_b*M_Zt->mc_length]);
      switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
      tab_x = tab_a + tab_s*tab_a_stride;
      output_QC_AtTYnWtSZn_uu->lf[tab_x]=0;
      my_j=0;
      while (my_j<M_Yn->rpop_j){
	my_a = M_Yn->m_a_[my_j]; my_b = M_Yn->m_b_[my_j];
	Yn_tag = (unsigned char *)&(M_Yn->wX[my_b*M_Yn->mc_length]);
	vA = bget____(At_tag,my_a);
	vT = bget____(Tt_tag,my_a);
	mw_j=0;
	while (mw_j<M_Wn->rpop_j){
	  mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	  Wn_tag = (unsigned char *)&(M_Wn->wX[mw_b*M_Wn->mc_length]);
	  vS = bget____(St_tag,mw_a);
	  vZ = bget____(Zt_tag,mw_a);
	  ny_j=0; nw_j=0;
	  while (ny_j<M_Yn->cpop_j && nw_j<M_Wn->cpop_j){
	    ny_a = M_Yn->n_a_[ny_j]; ny_b = M_Yn->n_b_[ny_j];
	    nw_a = M_Wn->n_a_[nw_j]; nw_b = M_Wn->n_b_[nw_j];
	    vY = bget____(Yn_tag,ny_a);
	    vW = bget____(Wn_tag,nw_a);
	    rtn_flag=1; if (AZ_flag && my_j==mw_j){ rtn_flag *= 0;} if (AY_flag && na_j==ny_j){ rtn_flag *= 0;}
	    output_QC_AtTYnWtSZn_uu->lf[tab_x] += rtn_flag* D_An[na_a/POPLENGTH] * (vA - a_An[na_a/POPLENGTH]) * vT * (vY - a_Yn[ny_a/POPLENGTH]) * D_Yn[ny_a/POPLENGTH] * (vW - a_Yn[nw_a/POPLENGTH]) * vS * (vZ - a_An[nz_a/POPLENGTH]);
	    ny_j++; nw_j++; /* while (ny_j<M_Yn->cpop_j && nw_j<M_Wn->cpop_j){ } */}
	  mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
	my_j++; /* while (my_j<M_Yn->rpop_j){ } */}
      na_j++;nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_At->rpop_j*M_Yn->rpop_j*M_Wn->rpop_j*M_Yn->cpop_j*2,0);
  if (verbose>1){ raprintf(output_QC_AtTYnWtSZn_uu->lf,"double",tab_a_stride,tab_s_stride," %% output_QC_AtTYnWtSZn_uu->lf: ");}
  if (verbose>1){ printf(" %% [finished get_QC_AtTYnWtSZn_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_QC_AtTYnWtSZn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,int AZ_flag,int AY_flag,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yn,struct M_handle *M_Wn,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,double *Y_ajdk,struct L_handle **output_p)
{
  /* This function calls get_QC_AtTYnWtSZn_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 17)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_QC_AtTYnWtSZn_uu__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Yn,verbose," %% M_Yn: ");}
  if (verbose){ M_handle_printf(M_Wn,verbose," %% M_Wn: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Zt,verbose," %% M_Zt: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  length = length_a*length_s; if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  if (verbose>2){ bprintf(M_Wn->mr_b,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_b: "); bprintf(M_Wn->mr_j,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_j: ");}
  if (verbose>2){ bprintf(M_Wn->mc_b,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_b: "); bprintf(M_Wn->mc_j,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_QC_AtTYnWtSZn_uu__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Yn; vpra[ip++] = M_Wn; vpra[ip++] = M_St; vpra[ip++] = M_Zt; vpra[ip++] = A_ajdk; vpra[ip++] = Y_ajdk; vpra[ip++] = *output_p;
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  switch (AZ_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break;}
  switch (AY_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break;}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_QC_AtTYnWtSZn_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_QC_AtTYnWtSZn_uu__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_QC_AtTYnWtSZn_uu(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p),"double",length_a,length_s," %% (*output_p): ");}
  if (verbose){ printf(" %% [finished wrap_QC_AtTYnWtSZn_uu__run] tidx %d\n",*tidx);}
  return length;
}  

void bcc_QC_AtTYnWtSZn_uu(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  int n_spacing_A = SPACING_a,n_spacing_B = SPACING_b;
  int AZ_flag=0,AY_flag=0;
  if (verbose){ printf(" %% [entering bcc_QC_AtTYnWtSZn_uu]\n");}
  if (verbose){ printf(" %% calculating F_[nbx]->QC_AtTYnWtSZn_uu, etc.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	AZ_flag=0; AY_flag=0;
	GLOBAL_pthread_tic(); 
	wrap_QC_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_B,AZ_flag,AY_flag,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_Yn,E_nb2->M_Wn,E_nb2->M_St,E_nb2->M_Zt,D->A_ajdk,D->Y_ajdk,&(F->QC_AtTYnWtSZn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	AZ_flag=0; AY_flag=1;
	GLOBAL_pthread_tic(); 
	wrap_QC_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_B,AZ_flag,AY_flag,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_An,E_nb2->M_Zn,E_nb2->M_St,E_nb2->M_Zt,D->A_ajdk,D->A_ajdk,&(F->QC_AtTAnZtSZn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	AZ_flag=(nb1==nb2); AY_flag=0;
	GLOBAL_pthread_tic(); 
	wrap_QC_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_B,AZ_flag,AY_flag,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_Yn,E_nb2->M_Yn,E_nb2->M_Tt,E_nb2->M_At,D->A_ajdk,D->Y_ajdk,&(F->QC_AtTYnYtTAn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	AZ_flag=(nb1==nb2); AY_flag=1;
	GLOBAL_pthread_tic(); 
	wrap_QC_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_B,AZ_flag,AY_flag,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_An,E_nb2->M_An,E_nb2->M_Tt,E_nb2->M_At,D->A_ajdk,D->A_ajdk,&(F->QC_AtTAnAtTAn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% QC_AtTYnWtSZn_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	lfprintf(F->QC_AtTYnWtSZn_uu," %% QC_AtTYnWtSZn_uu: ");
	lfprintf(F->QC_AtTAnZtSZn_uu," %% QC_AtTAnZtSZn_uu: ");
	bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	lfprintf(F->QC_AtTYnYtTAn_uu," %% QC_AtTYnYtTAn_uu: ");
	lfprintf(F->QC_AtTAnAtTAn_uu," %% QC_AtTAnAtTAn_uu: ");
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished bcc_QC_AtTYnWtSZn_uu]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void *get_flattenloop(void *vp)
{
  /* This function corrects the row- and col-scores, removing flattened-loops.
     Along the way we also convert from SPACING_a to SPACING_b along the T_ncols dimension. */
  int verbose=0;
  int ip=0,nc=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct bcc_double *F = (struct bcc_double *)(vpra[ip++]); struct bcc_ajdk *D = F->D; struct bcc_single *E_nb1 = F->E_nb1, *E_nb2 = F->E_nb2;
  struct M_handle *M_An_nb1 = E_nb1->M_An;
  struct M_handle *M_At_nb1 = E_nb1->M_At;
  struct M_handle *M_Tt_nb1 = E_nb1->M_Tt;
  struct M_handle *M_Yt_nb1 = E_nb1->M_Yt;
  /* struct M_handle *M_Wt_nb1 = E_nb1->M_Wt; */
  /* struct M_handle *M_St_nb1 = E_nb1->M_St; */
  /* struct M_handle *M_Zn_nb1 = E_nb1->M_Zn; */
  /* struct M_handle *M_Zt_nb1 = E_nb1->M_Zt; */
  /* struct M_handle *M_An_nb2 = E_nb2->M_An; */
  struct M_handle *M_At_nb2 = E_nb2->M_At;
  struct M_handle *M_Tt_nb2 = E_nb2->M_Tt;
  /* struct M_handle *M_Yt_nb2 = E_nb2->M_Yt; */
  /* struct M_handle *M_Wt_nb2 = E_nb2->M_Wt; */
  struct M_handle *M_St_nb2 = E_nb2->M_St;
  /* struct M_handle *M_Zn_nb2 = E_nb2->M_Zn; */
  struct M_handle *M_Zt_nb2 = E_nb2->M_Zt;
  int use_At_T_XnXt_X_Xn_flag = (strstr(D->QC_strategy,"YnWt") ? 1 : 0);
  struct L_handle *lf_At_T_AnAt_T_An = F->lf_At_T_AnAt_T_An;
  struct L_handle *lf_At_T_AnZt_S_Zn = F->lf_At_T_AnZt_S_Zn;
  struct L_handle *lf_At_T_YnYt_T_An = F->lf_At_T_YnYt_T_An;
  struct L_handle *lf_At_T_YnWt_S_Zn = F->lf_At_T_YnWt_S_Zn;
  int use_AtTXn____XtXXn_flag = (strstr(D->QC_strategy,"ZtSWn") ? 1 : 0);
  struct L_handle *lf_AtTAn____AtTAn = F->lf_AtTAn____AtTAn;
  struct L_handle *lf_AtTAn____ZtSZn = F->lf_AtTAn____ZtSZn;
  struct L_handle *lf_AtTYn____YtTAn = F->lf_AtTYn____YtTAn;
  struct L_handle *lf_AtTYn____WtSZn = F->lf_AtTYn____WtSZn; 
  int use_An_XtXXn_Xt_flag = (strstr(D->QR_strategy,"ZtSWn") ? 1 : 0);
  struct L_handle *lf_An_ZtSWn_Yt = F->lf_An_ZtSWn_Yt;
  struct L_handle *lf_An_ZtSZn_At = F->lf_An_ZtSZn_At;
  struct L_handle *lf_An_AtTYn_Yt = F->lf_An_AtTYn_Yt;
  struct L_handle *lf_An_AtTAn_At = F->lf_An_AtTAn_At;
  int use_AnXt_X_XnXt_flag = (strstr(D->QR_strategy,"YnWt") ? 1 : 0);
  struct L_handle *lf_AnZt_S_WnYt = F->lf_AnZt_S_WnYt;
  struct L_handle *lf_AnZt_S_ZnAt = F->lf_AnZt_S_ZnAt;
  struct L_handle *lf_AnAt_T_YnYt = F->lf_AnAt_T_YnYt;
  struct L_handle *lf_AnAt_T_AnAt = F->lf_AnAt_T_AnAt;
  double *A_ajdk = D->A_ajdk; double *Y_ajdk = D->Y_ajdk;
  int A_pcols = psize(M_At_nb1->nrows);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yt_nb1->nrows);
  int A_ncols = M_At_nb1->nrows, T_ncols = M_Tt_nb1->nrows;
  int A_nrows_nb1 = M_At_nb1->ncols, A_nrows_nb2 = M_At_nb2->ncols;
  int /* Z_nrows_nb1 = M_Zt_nb1->ncols, */Z_nrows_nb2 = M_Zt_nb2->ncols;
  int QR_AnZtSWnYt_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QR_AnZtSZnAt_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QR_AnAtTYnYt_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QR_AnAtTAnAt_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  int QC_AtTYnWtSZn_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QC_AtTAnZtSZn_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QC_AtTYnYtTAn_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QC_AtTAnAtTAn_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;  
  struct L_handle *lf_tmp=NULL;
  lf_tmp = F->QR_AnZtSWnYt; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_AnZtSZnAt; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_AnAtTYnYt; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_AnAtTAnAt; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_AtTYnWtSZn; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_AtTAnZtSZn; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_AtTYnYtTAn; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_AtTAnAtTAn; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  struct L_handle *QR_AnZtSWnYt = F->QR_AnZtSWnYt;
  struct L_handle *QR_AnZtSZnAt = F->QR_AnZtSZnAt;
  struct L_handle *QR_AnAtTYnYt = F->QR_AnAtTYnYt;
  struct L_handle *QR_AnAtTAnAt = F->QR_AnAtTAnAt;
  struct L_handle *QC_AtTYnWtSZn = F->QC_AtTYnWtSZn;
  struct L_handle *QC_AtTAnZtSZn = F->QC_AtTAnZtSZn;
  struct L_handle *QC_AtTYnYtTAn = F->QC_AtTYnYtTAn;
  struct L_handle *QC_AtTAnAtTAn = F->QC_AtTAnAtTAn;
  int vM=0,vT=0/* ,vS=0 */;
  /* unsigned char *An_tag=NULL; */
  /* unsigned char *At_tag=NULL; */
  unsigned char *Tt_tag=NULL;
  /* unsigned char *Yt_tag=NULL; */
  /* unsigned char *Wt_tag=NULL; */
  /* unsigned char *Zn_tag=NULL; */
  /* unsigned char *Zt_tag=NULL; */
  /* unsigned char *St_tag=NULL; */
  /* unsigned char *mca_tag=NULL,*mca_end=NULL; */
  /* unsigned char *mcz_tag=NULL,*mcz_end=NULL; */
  unsigned char *mct_tag=NULL,*mct_end=NULL;
  /* unsigned char *mcs_tag=NULL,*mcs_end=NULL; */
  int ns_j=0,ns_b=0,ns_a=0,tab_Tt_stride=0,tab_Tt=0,tab_St_stride=0,tab_St=0;
  int na_j=0,na_b=0,na_a=0,tab_At_stride=0,tab_At=0,ma_j=0,ma_b=0,ma_a=0,tab_An_stride=0,tab_An=0,tab_YnWt=0;
  /* int ny_j=0,ny_b=0,ny_a=0,tab_Yt_stride=0,tab_Yt=0,my_j=0,my_b=0,my_a=0,tab_Yn_stride=0,tab_Yn=0; */
  /* int nw_j=0,nw_b=0,nw=0,tab_Wt_stride=0,tab_Wt=0,mw_j=0,mw_b=0,mw=0,tab_Wn_stride=0,tab_Wn=0; */
  /* int nz_j=0,nz_b=0,nz_a=0,tab_Zt_stride=0,tab_Zt=0,mz_j=0,mz_b=0,mz_a=0,tab_Zn_stride=0,tab_Zn=0; */
  double lf_A_a0d2=0,lf_A_a1d2=0,lf_A_a2d2=0,lf_A_a3d2=0,lf_A_a4d2=0;
  double lf_A_a0d1=0,lf_A_a1d1=0,lf_A_a2d1=0,lf_A_a3d1=0,lf_A_a4d1=0;
  double lf_Y_a0d1=0,lf_Y_a1d1=0,lf_Y_a2d1=0,lf_Y_a3d1=0,lf_Y_a4d1=0;
  double dtmp=0;
  if (verbose>1){ printf(" %% [entering get_flattenloop] tidx %d\n",tidx);} if (verbose>3){ printf(" %% /******************************************************************/\n");}
  if (QR_AnZtSWnYt_bother || QR_AnZtSZnAt_bother || QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother || QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){
    lf_A_a0d2=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a0d2 += (double)(A_ajdk ? A_ajdk[nc + AJDK_0_2*A_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_A_a1d2=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a1d2 += (double)(A_ajdk ? A_ajdk[nc + AJDK_1_2*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_A_a2d2=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a2d2 += (double)(A_ajdk ? A_ajdk[nc + AJDK_2_2*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_A_a3d2=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a3d2 += (double)(A_ajdk ? A_ajdk[nc + AJDK_3_2*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_A_a4d2=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a4d2 += (double)(A_ajdk ? A_ajdk[nc + AJDK_4_2*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_A_a0d1=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a0d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_0_1*A_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_A_a1d1=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a1d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_1_1*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_A_a2d1=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a2d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_2_1*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_A_a3d1=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a3d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_3_1*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_A_a4d1=0; for (nc=0;nc<A_pcols;nc++){ lf_A_a4d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_4_1*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An_nb1->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_Y_a0d1=0; for (nc=0;nc<Y_pcols;nc++){ lf_Y_a0d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_0_1*Y_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_Yt_nb1->mr_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_Y_a1d1=0; for (nc=0;nc<Y_pcols;nc++){ lf_Y_a1d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_1_1*Y_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_Yt_nb1->mr_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_Y_a2d1=0; for (nc=0;nc<Y_pcols;nc++){ lf_Y_a2d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_2_1*Y_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_Yt_nb1->mr_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_Y_a3d1=0; for (nc=0;nc<Y_pcols;nc++){ lf_Y_a3d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_3_1*Y_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_Yt_nb1->mr_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    lf_Y_a4d1=0; for (nc=0;nc<Y_pcols;nc++){ lf_Y_a4d1 += (double)(Y_ajdk ? Y_ajdk[nc + AJDK_4_1*Y_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_Yt_nb1->mr_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
    /*   if (QR_AnZtSWnYt_bother || QR_AnZtSZnAt_bother || QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother || QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QR_AnZtSWnYt_bother || QR_AnZtSZnAt_bother || QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother){
    if (verbose>3){ printf(" %% QR_AnZtSWnYt QR_AnZtSZnAt QR_AnAtTYnYt QR_AnAtTAnAt\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      ma_j=0; while (ma_j<M_An_nb1->rpop_j){
	ma_a = M_An_nb1->m_a_[ma_j]; ma_b = M_An_nb1->m_b_[ma_j];
	if (0){ /* do nothing */}
	else if (use_AnXt_X_XnXt_flag){ if (verbose>1){ printf(" %% adopting lf_AnZt_S_WnYt\n");}
	  if (QR_AnZtSWnYt_bother){ L2_set(QR_AnZtSWnYt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_AnZt_S_WnYt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	  if (QR_AnZtSZnAt_bother){ L2_set(QR_AnZtSZnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_AnZt_S_ZnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	  if (QR_AnAtTYnYt_bother){ L2_set(QR_AnAtTYnYt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_AnAt_T_YnYt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	  if (QR_AnAtTAnAt_bother){ L2_set(QR_AnAtTAnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_AnAt_T_AnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	  /* if (use_AnXt_X_XnXt_flag){ } */}
	else if (use_An_XtXXn_Xt_flag){ if (verbose>1){ printf(" %% adopting lf_An_ZtSWn_Yt\n");}
	  if (QR_AnZtSWnYt_bother){ L2_set(QR_AnZtSWnYt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_An_ZtSWn_Yt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	  if (QR_AnZtSZnAt_bother){ L2_set(QR_AnZtSZnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_An_ZtSZn_At , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	  if (QR_AnAtTYnYt_bother){ L2_set(QR_AnAtTYnYt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_An_AtTYn_Yt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	  if (QR_AnAtTAnAt_bother){ L2_set(QR_AnAtTAnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_An_AtTAn_At , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	  /* else if (use_An_XtXXn_Xt_flag){ } */}
	ma_j++; /* while (ma_j<M_An_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QR_AnZtSWnYt," %% QR_AnZtSWnYt: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QR_AnZtSZnAt," %% QR_AnZtSZnAt: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QR_AnAtTYnYt," %% QR_AnAtTYnYt: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QR_AnAtTAnAt," %% QR_AnAtTAnAt: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if (QR_AnZtSWnYt_bother || QR_AnZtSZnAt_bother || QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother){ } */}
  if (QR_AnZtSZnAt_bother || QR_AnAtTAnAt_bother){
    if (verbose>3){ printf(" %% corrections to QR_AnZtSZnAt, and QR_AnAtTAnAt (both always necessary)\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /*    QR_AnZtZnAt_1(:,nt0) = diag(AnZt_1*diag(S(:,nt0))*transpose(AnZt_1)) - ...
	  ( ... ;
	  + e_An*f_At*D_An.^2*f_An*(e_Zt*S(:,nt0)) ...
	  -2 * e_An*f_At*D_An.^2*diag(a_An)*Zt*S(:,nt0) ...
	  + e_An*f_At*D_An.^2*a_An.^2*(e_Zt*S(:,nt0)) ...
	  -2 * An * diag(a_An)*D_An.^2*f_An*(e_Zt*S(:,nt0)) ...
	  +4 * An * diag(a_An)*D_An.^2*diag(a_An)*Zt*S(:,nt0) ...
	  -2 * An * diag(a_An)*D_An.^2*a_An.^2*e_Zt*S(:,nt0) ...
	  + e_An*a_At.^2*D_An.^2*f_An*(e_Zt*S(:,nt0)) ...
	  -2 * e_An*a_At.^2*D_An.^2*diag(a_An)*Zt*S(:,nt0) ...
	  + e_An*a_At.^2*D_An.^2*a_An.^2*e_Zt*S(:,nt0) ...
	  ) ;
    */
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      ma_j=0; while (ma_j<M_An_nb1->rpop_j){
	ma_a = M_An_nb1->m_a_[ma_j]; ma_b = M_An_nb1->m_b_[ma_j];
	if (QR_AnZtSZnAt_bother){
	  dtmp = 
	    + 1*lf_A_a0d2*(*L1_get(E_nb2->lf_et_Sn , ns_j,ns_b,ns_a))
	    - 2*(*L1_get(E_nb2->lf_a1d2_Zt_Sn, ns_j,ns_b,ns_a))
	    + 1*lf_A_a2d2*(*L1_get(E_nb2->lf_et_Sn, ns_j,ns_b,ns_a))
	    - 2*(*L2_get(E_nb1->lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_2)) * (*L1_get(E_nb2->lf_et_Sn , ns_j,ns_b,ns_a))
	    + 4*(*L2_get(F->lf_An_a2d2_Zt_Sn , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a))
	    - 2*(*L2_get(E_nb1->lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_3_2)) * (*L1_get(E_nb2->lf_et_Sn , ns_j,ns_b,ns_a))
	    + 1*lf_A_a2d2*(*L1_get(E_nb2->lf_et_Sn , ns_j,ns_b,ns_a))
	    - 2*(*L1_get(E_nb2->lf_a3d2_Zt_Sn , ns_j,ns_b,ns_a))
	    + 1*lf_A_a4d2*(*L1_get(E_nb2->lf_et_Sn , ns_j,ns_b,ns_a));
	  L2_plusequals(QR_AnZtSZnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , -dtmp);
	  /* if (QR_AnZtSZnAt_bother){ } */}
	if (QR_AnAtTAnAt_bother){
	  dtmp = 
	    + 1*lf_A_a0d2*(*L1_get(E_nb2->lf_et_Tn , ns_j,ns_b,ns_a))
	    - 2*(*L1_get(E_nb2->lf_a1d2_At_Tn , ns_j,ns_b,ns_a))
	    + 1*lf_A_a2d2*(*L1_get(E_nb2->lf_et_Tn , ns_j,ns_b,ns_a))
	    - 2*(*L2_get(E_nb1->lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_2)) * (*L1_get(E_nb2->lf_et_Tn , ns_j,ns_b,ns_a))
	    + 4*(*L2_get(F->lf_An_a2d2_At_Tn , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a))
	    - 2*(*L2_get(E_nb1->lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_3_2)) * (*L1_get(E_nb2->lf_et_Tn , ns_j,ns_b,ns_a))
	    + 1*lf_A_a2d2*(*L1_get(E_nb2->lf_et_Tn , ns_j,ns_b,ns_a))
	    - 2*(*L1_get(E_nb2->lf_a3d2_At_Tn , ns_j,ns_b,ns_a))
	    + 1*lf_A_a4d2*(*L1_get(E_nb2->lf_et_Tn , ns_j,ns_b,ns_a));
	  L2_plusequals(QR_AnAtTAnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , -dtmp);
	  /* if (QR_AnAtTAnAt_bother){ } */}
	ma_j++; /* while (ma_j<M_An_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QR_AnZtSZnAt," %% QR_AnZtSZnAt: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QR_AnAtTAnAt," %% QR_AnAtTAnAt: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if (QR_AnZtSZnAt_bother || QR_AnAtTAnAt_bother){ } */}
  if (QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother){
    if (verbose>3){ printf(" %% corrections to QR_AnAtTYnYt, and QR_AnAtTAnAt (necessary if F->nb1==F->nb2)\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (F->nb1==F->nb2){
      /* QR_AnAtYnYt_1(:,nt0) = diag(AnAt_1*diag(T(:,nt0))*transpose(YnYt_1)) - T(:,nt0).*diag(AnAt_1).*diag(YnYt_1); */
      ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
	ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
	ma_j=0; while (ma_j<M_An_nb1->rpop_j){
	  ma_a = M_An_nb1->m_a_[ma_j]; ma_b = M_An_nb1->m_b_[ma_j];
	  if (QR_AnAtTYnYt_bother){ 
	    dtmp = (*L2_get(E_nb1->lf_T_AnAt_YnYt , ns_j,ns_b,ns_a , ma_j,ma_b,ma_a));
	    L2_plusequals(QR_AnAtTYnYt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , -dtmp);
	    /* if (QR_AnAtTYnYt_bother){ } */}
	  if (QR_AnAtTAnAt_bother){ 
	    dtmp = (*L2_get(E_nb1->lf_T_AnAt_AnAt , ns_j,ns_b,ns_a , ma_j,ma_b,ma_a));
	    L2_plusequals(QR_AnAtTAnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , -dtmp);
	    /* if (QR_AnAtTAnAt_bother){ } */}
	  ma_j++; /* while (ma_j<M_An_nb1->rpop_j){ } */}
	ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
      if (verbose>3){ lfprintf(QR_AnAtTYnYt," %% QR_AnAtTYnYt: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
      if (verbose>3){ lfprintf(QR_AnAtTAnAt," %% QR_AnAtTAnAt: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
      /* if (F->nb1==F->nb2){ } */}
    /* if (QR_AnAtTYnYt_bother || QR_AnAtTAnAt_bother){ } */}
  if (QR_AnAtTAnAt_bother){
    if (verbose>3){ printf(" %% corrections to QR_AnAtTAnAt (necessary if F->nb1==F->nb2)\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (F->nb1==F->nb2){
      /* 
	 QR_AnAtAnAt_1(:,nt0) = diag(AnAt_1*diag(T(:,nt0))*transpose(AnAt_1)) - T(:,nt0).*diag(AnAt_1).*diag(AnAt_1) - ...
	 ( ... ;
	 + e_An*f_At*D_An.^2*f_An*(e_At*T(:,nt0)) ...
	 -2 * e_An*f_At*D_An.^2*diag(a_An)*At*T(:,nt0) ...
	 + e_An*f_At*D_An.^2*a_An.^2*(e_At*T(:,nt0)) ...
	 -2 * An * diag(a_An)*D_An.^2*f_An*(e_At*T(:,nt0)) ...
	 +4 * An * diag(a_An)*D_An.^2*diag(a_An)*At*T(:,nt0) ...
	 -2 * An * diag(a_An)*D_An.^2*a_An.^2*e_At*T(:,nt0) ...
	 + e_An*a_At.^2*D_An.^2*f_An*(e_At*T(:,nt0)) ...
	 -2 * e_An*a_At.^2*D_An.^2*diag(a_An)*At*T(:,nt0) ...
	 + e_An*a_At.^2*D_An.^2*a_An.^2*e_At*T(:,nt0) ...
	 ) ...
	 + ...
	 T(:,nt0).*( ...
	 + e_An*f_At*D_An.^2*f_An ...
	 -4 * An*diag(a_An)*D_An.^2*f_An ...
	 +6 * e_An*a_At.^2*D_An.^2*f_An ...
	 -4 * An*diag(a_An.^3)*D_An.^2*f_An ...
	 +1 * e_An*a_At.^4*D_An.^2*f_An ...
	 )...
	 ; */
      ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
	ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
	Tt_tag = (unsigned char *)(&(M_Tt_nb2->wX[ns_b*M_Tt_nb2->mc_length]));
	mct_tag = (unsigned char *)(&(M_Tt_nb2->mc_j[0]));
	mct_end = (unsigned char *)(&(M_Tt_nb2->mc_j[M_Tt_nb2->mc_length]));    
	ma_j=0; while (ma_j<M_An_nb1->rpop_j){
	  ma_a = M_An_nb1->m_a_[ma_j]; ma_b = M_An_nb1->m_b_[ma_j];
	  vT = bget____(Tt_tag,ma_a); vM = bget__on(mct_tag,ma_a); if(!vM){ printf(" %% Warning! bit %d,%d,%d off in M_Tt->mc_j\n",ma_j,ma_b,ma_a);}
	  dtmp = 
	    vT*(
		+ 1*lf_A_a0d2
		- 4*(*L2_get(E_nb1->lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_1_2))
		+ 6*lf_A_a2d2
		- 4*(*L2_get(E_nb1->lf_An_ajdk , ma_j,ma_b,ma_a , 0,0,AJDK_3_2))
		+ 1*lf_A_a4d2
		);
	  L2_plusequals(QR_AnAtTAnAt , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , +dtmp);
	  ma_j++; /* while (ma_j<M_An_nb1->rpop_j){ } */}
	ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
      if (verbose>3){ lfprintf(QR_AnAtTAnAt," %% QR_AnAtTAnAt: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
      /* if (F->nb1==F->nb2){ } */}
    /* if (QR_AnAtTAnAt_bother){ } */}
  if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){
    if (verbose>3){ printf(" %% QC_AtTYnWtSZn, QC_AtTAnZtSZn, QC_AtTYnYtTAn, QC_AtTAnAtTAn\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      na_j=0; while (na_j<M_At_nb1->rpop_j){
	na_a = M_At_nb1->m_a_[na_j]; na_b = M_At_nb1->m_b_[na_j];
	if (0){ /* do nothing */}
	else if (use_AtTXn____XtXXn_flag){ if (verbose>1){ printf(" %% adopting lf_AtTYn____WtSZn\n");}
	  if (QC_AtTYnWtSZn_bother){ L2_set(QC_AtTYnWtSZn , na_j,na_b,na_a , ns_j,ns_b,ns_a , D_An[na_a/POPLENGTH] * (*L2_get(lf_AtTYn____WtSZn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  if (QC_AtTAnZtSZn_bother){ L2_set(QC_AtTAnZtSZn , na_j,na_b,na_a , ns_j,ns_b,ns_a , D_An[na_a/POPLENGTH] * (*L2_get(lf_AtTAn____ZtSZn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  if (QC_AtTYnYtTAn_bother){ L2_set(QC_AtTYnYtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a , D_An[na_a/POPLENGTH] * (*L2_get(lf_AtTYn____YtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  if (QC_AtTAnAtTAn_bother){ L2_set(QC_AtTAnAtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a , D_An[na_a/POPLENGTH] * (*L2_get(lf_AtTAn____AtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  /* if (use_At_T_XnXt_X_Xn_flag){ } */}
	else if (use_At_T_XnXt_X_Xn_flag){ if (verbose>1){ printf(" %% adopting lf_At_T_YnWt_S_Zn\n");}
	  if (QC_AtTYnWtSZn_bother){ L2_set(QC_AtTYnWtSZn , na_j,na_b,na_a , ns_j,ns_b,ns_a , D_An[na_a/POPLENGTH] * (*L2_get(lf_At_T_YnWt_S_Zn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  if (QC_AtTAnZtSZn_bother){ L2_set(QC_AtTAnZtSZn , na_j,na_b,na_a , ns_j,ns_b,ns_a , D_An[na_a/POPLENGTH] * (*L2_get(lf_At_T_AnZt_S_Zn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  if (QC_AtTYnYtTAn_bother){ L2_set(QC_AtTYnYtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a , D_An[na_a/POPLENGTH] * (*L2_get(lf_At_T_YnYt_T_An , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  if (QC_AtTAnAtTAn_bother){ L2_set(QC_AtTAnAtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a , D_An[na_a/POPLENGTH] * (*L2_get(lf_At_T_AnAt_T_An , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	  /* if (use_At_T_XnXt_X_Xn_flag){ } */}
	na_j++; /* while (na_j<M_At_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QC_AtTYnWtSZn," %% QC_AtTYnWtSZn: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QC_AtTAnZtSZn," %% QC_AtTAnZtSZn: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QC_AtTYnYtTAn," %% QC_AtTYnYtTAn: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QC_AtTAnAtTAn," %% QC_AtTAnAtTAn: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if (QC_AtTYnWtSZn_bother || QC_AtTAnZtSZn_bother || QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QC_AtTAnZtSZn_bother || QC_AtTAnAtTAn_bother){
    if (verbose>3){ printf(" %% corrections to QC_AtTAnZtSZn and QC_AtTAnAtTAn (always needed)\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /*
      QC_AtAnZtZn_1(:,nt0) = diag(D_An*(At-a_An*e_At)*diag(T(:,nt0))*AnZt_1*diag(S(:,nt0))*(Zn-e_Zn*a_At)) - ...
      (D_An * (...
      + f_An*e_At ...
      -2 * diag(a_An)*At ...
      + a_An.^2*e_At ...
      )*T(:,nt0)) .* ...
      (D_An * (...
      + f_An*e_Zt ...
      -2 * diag(a_An)*Zt ...
      + a_An.^2*e_Zt ...
      )*S(:,nt0));
    */  
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      na_j=0; while (na_j<M_At_nb1->rpop_j){
	na_a = M_At_nb1->m_a_[na_j]; na_b = M_At_nb1->m_b_[na_j];
	if (QC_AtTAnZtSZn_bother){
	  dtmp = 
	    A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*
	    (
	     +(1*(*L1_get(E_nb1->lf_et_Tn , ns_j,ns_b,ns_a)))
	     -(2*A_ajdk[na_a/POPLENGTH  + AJDK_1_0*A_pcols] * (*L2_get(E_nb1->lf_AtTn , na_j,na_b,na_a , ns_j,ns_b,ns_a)))
	     +(A_ajdk[na_a/POPLENGTH  + AJDK_2_0*A_pcols] * (*L1_get(E_nb1->lf_et_Tn , ns_j,ns_b,ns_a)))
	     )*
	    A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*
	    (
	     +(1*(*L1_get(E_nb2->lf_et_Sn , ns_j,ns_b,ns_a)))
	     -(2*A_ajdk[na_a/POPLENGTH  + AJDK_1_0*A_pcols] * (*L2_get(E_nb2->lf_ZtSn , na_j,na_b,na_a , ns_j,ns_b,ns_a)))
	     +(A_ajdk[na_a/POPLENGTH  + AJDK_2_0*A_pcols] * (*L1_get(E_nb2->lf_et_Sn , ns_j,ns_b,ns_a)))
	     );
	  L2_plusequals(QC_AtTAnZtSZn , na_j,na_b,na_a , ns_j,ns_b,ns_a , -dtmp);
	  /* if (QC_AtTAnZtSZn_bother){ } */}
	if (QC_AtTAnAtTAn_bother){
	  dtmp = 
	    A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*
	    (
	     +(1*(*L1_get(E_nb1->lf_et_Tn , ns_j,ns_b,ns_a)))
	     -(2*A_ajdk[na_a/POPLENGTH  + AJDK_1_0*A_pcols] * (*L2_get(E_nb1->lf_AtTn , na_j,na_b,na_a , ns_j,ns_b,ns_a)))
	     +(A_ajdk[na_a/POPLENGTH  + AJDK_2_0*A_pcols] * (*L1_get(E_nb1->lf_et_Tn , ns_j,ns_b,ns_a)))
	     )*
	    A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*
	    (
	     +(1*(*L1_get(E_nb2->lf_et_Tn , ns_j,ns_b,ns_a)))
	     -(2*A_ajdk[na_a/POPLENGTH  + AJDK_1_0*A_pcols] * (*L2_get(E_nb2->lf_AtTn , na_j,na_b,na_a , ns_j,ns_b,ns_a)))
	     +(A_ajdk[na_a/POPLENGTH  + AJDK_2_0*A_pcols] * (*L1_get(E_nb2->lf_et_Tn , ns_j,ns_b,ns_a)))
	     );
	  L2_plusequals(QC_AtTAnAtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a , -dtmp);
	  /* if (QC_AtTAnAtTAn_bother){ } */}
	na_j++; /* while (na_j<M_At_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QC_AtTAnZtSZn," %% QC_AtTAnZtSZn: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QC_AtTAnAtTAn," %% QC_AtTAnAtTAn: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if (QC_AtTAnZtSZn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){
    if (verbose>3){ printf(" %% correcting QC_AtTYnYtTAn and QC_AtTAnAtTAn (necessary if F->nb1==F->nb2)\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* QC_AtYnYtAn_1(:,nt0) = diag(D_An*(At-a_An*e_At)*diag(T(:,nt0))*YnYt_1*diag(T(:,nt0))*(An-e_An*a_At)) - ...
       (...
       +1 * D_An*f_An*e_At*mc_T(nt0)*e_An*f_Yt*D_Yn*f_Yn ... % + diag(D_An)*transpose(T(:,nt0).^2)*e_An*f_Yt*D_Yn*f_Yn ...
       -2 * D_An*diag(a_At)*At*mc_T(nt0)*e_An*f_Yt*D_Yn*f_Yn ... % -2 * D_An*transpose(An*diag(a_An))*diag(T(:,nt0).^2)*e_An*f_Yt*D_Yn*f_Yn ...
       +1 * D_An*a_An.^2*e_At*mc_T(nt0)*e_An*f_Yt*D_Yn*f_Yn ... % + (D_An*a_An.^2)*transpose(T(:,nt0).^2)*e_An*f_Yt*D_Yn*f_Yn ...
       -2 * D_An*f_An*e_At*mc_T(nt0)*Yn*diag(a_Yn)*D_Yn*f_Yn ... % -2 * diag(D_An)*transpose(T(:,nt0).^2)*Yn*diag(a_Yn)*D_Yn*f_Yn ...
       +4 * D_An*diag(a_At)*At*mc_T(nt0)*Yn*diag(a_Yn)*D_Yn*f_Yn ... % +4 * D_An*transpose(An*diag(a_An))*diag(T(:,nt0).^2)*Yn*diag(a_Yn)*D_Yn*f_Yn ...
       -2 * D_An*a_An.^2*e_At*mc_T(nt0)*Yn*diag(a_Yn)*D_Yn*f_Yn ... % -2 * (D_An*a_An.^2)*transpose(T(:,nt0).^2)*Yn*diag(a_Yn)*D_Yn*f_Yn ...
       +1 * D_An*f_An*e_At*mc_T(nt0)*e_An*a_Yt.^2*D_Yn*f_Yn ... % + diag(D_An)*transpose(T(:,nt0).^2)*e_An*a_Yt.^2*D_Yn*f_Yn ...
       -2 * D_An*diag(a_At)*At*mc_T(nt0)*e_An*a_Yt.^2*D_Yn*f_Yn ... % -2 * D_An*transpose(An*diag(a_An))*diag(T(:,nt0).^2)*e_An*a_Yt.^2*D_Yn*f_Yn ...
       +1 * D_An*a_An.^2*e_At*mc_T(nt0)*e_An*a_Yt.^2*D_Yn*f_Yn ...% + (D_An*a_An.^2)*transpose(T(:,nt0).^2)*e_An*a_Yt.^2*D_Yn*f_Yn ...
       );*/
    if (F->nb1==F->nb2){
      ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
	ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
	na_j=0; while (na_j<M_At_nb1->rpop_j){
	  na_a = M_At_nb1->m_a_[na_j]; na_b = M_At_nb1->m_b_[na_j];
	  if (QC_AtTYnYtTAn_bother){	  
	    dtmp = 
	      (
	       +(1*A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*(M_An_nb1->rpop_j)*lf_Y_a0d1)
	       -(2*A_ajdk[na_a/POPLENGTH + AJDK_1_1*A_pcols]*(*L1_get(E_nb1->lf_et_An , na_j,na_b,na_a))*lf_Y_a0d1)
	       +(1*A_ajdk[na_a/POPLENGTH + AJDK_2_1*A_pcols]*(M_An_nb1->rpop_j)*lf_Y_a0d1)
	       -(2*A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*(*L1_get(E_nb1->lf_et_Yn_a1d1 , 0,0,0)))
	       +(4*A_ajdk[na_a/POPLENGTH + AJDK_1_1*A_pcols]*(*L1_get(E_nb1->lf_At_Yn_a1d1 , na_j,na_b,na_a)))
	       -(2*A_ajdk[na_a/POPLENGTH + AJDK_2_1*A_pcols]*(*L1_get(E_nb1->lf_et_Yn_a1d1 , 0,0,0)))
	       +(1*A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*(M_An_nb1->rpop_j)*lf_Y_a2d1)
	       -(2*A_ajdk[na_a/POPLENGTH + AJDK_1_1*A_pcols]*(*L1_get(E_nb1->lf_et_An , na_j,na_b,na_a))*lf_Y_a2d1)
	       +(1*A_ajdk[na_a/POPLENGTH + AJDK_2_1*A_pcols]*(M_An_nb1->rpop_j)*lf_Y_a2d1)
	       );
	    L2_plusequals(QC_AtTYnYtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a , -dtmp);
	    /* if (QC_AtTYnYtTAn_bother){ } */}
	  if (QC_AtTAnAtTAn_bother){
	    dtmp = 
	      (
	       +(1*A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*(M_An_nb1->rpop_j)*lf_A_a0d1)
	       -(2*A_ajdk[na_a/POPLENGTH + AJDK_1_1*A_pcols]*(*L1_get(E_nb1->lf_et_An , na_j,na_b,na_a))*lf_A_a0d1)
	       +(1*A_ajdk[na_a/POPLENGTH + AJDK_2_1*A_pcols]*(M_An_nb1->rpop_j)*lf_A_a0d1)
	       -(2*A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*(*L1_get(E_nb1->lf_et_An_a1d1 , 0,0,0)))
	       +(4*A_ajdk[na_a/POPLENGTH + AJDK_1_1*A_pcols]*(*L1_get(E_nb1->lf_At_An_a1d1 , na_j,na_b,na_a)))
	       -(2*A_ajdk[na_a/POPLENGTH + AJDK_2_1*A_pcols]*(*L1_get(E_nb1->lf_et_An_a1d1 , 0,0,0)))
	       +(1*A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols]*(M_An_nb1->rpop_j)*lf_A_a2d1)
	       -(2*A_ajdk[na_a/POPLENGTH + AJDK_1_1*A_pcols]*(*L1_get(E_nb1->lf_et_An , na_j,na_b,na_a))*lf_A_a2d1)
	       +(1*A_ajdk[na_a/POPLENGTH + AJDK_2_1*A_pcols]*(M_An_nb1->rpop_j)*lf_A_a2d1)
	       );
	    L2_plusequals(QC_AtTAnAtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a , -dtmp);
	    /* if (QC_AtTAnAtTAn_bother){ } */}
	  na_j++; /* while (na_j<M_At_nb1->rpop_j){ } */}
	ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
      if (verbose>3){ lfprintf(QC_AtTYnYtTAn," %% QC_AtTYnYtTAn: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
      if (verbose>3){ lfprintf(QC_AtTAnAtTAn," %% QC_AtTAnAtTAn: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
      /* if (F->nb1==F->nb2){ } */}
    /* if (QC_AtTYnYtTAn_bother || QC_AtTAnAtTAn_bother){ } */}
  if (QC_AtTAnAtTAn_bother){
    if (verbose>3){ printf(" %% correcting QC_At_T_AnAt_T_An (necessary if F->nb1==F->nb2)\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /*
      QC_AtAnAtAn_1(:,nt0) = diag(D_An*(At-a_An*e_At)*diag(T(:,nt0))*AnAt_1*diag(T(:,nt0))*(An-e_An*a_At)) - ...
      (D_An * (...
      + f_An*e_At ...
      -2 * diag(a_An)*At ...
      + a_An.^2*e_At ...
      )*T(:,nt0)) .* ...
      (D_An * (...
      + f_An*e_At ...
      -2 * diag(a_An)*At ...
      + a_An.^2*e_At ...
      )*T(:,nt0)) ...
      - ...
      (...
      +1 * D_An*f_An*e_At*mc_T(nt0)*e_An*f_At*D_An*f_An ... % + diag(D_An)*transpose(T(:,nt0).^2)*e_An*f_At*D_An*f_An ...
      -2 * D_An*diag(a_At)*At*mc_T(nt0)*e_An*f_At*D_An*f_An ... % -2 * D_An*transpose(An*diag(a_An))*diag(T(:,nt0).^2)*e_An*f_At*D_An*f_An ...
      +1 * D_An*a_An.^2*e_At*mc_T(nt0)*e_An*f_At*D_An*f_An ... % + (D_An*a_An.^2)*transpose(T(:,nt0).^2)*e_An*f_At*D_An*f_An ...
      -2 * D_An*f_An*e_At*mc_T(nt0)*An*diag(a_An)*D_An*f_An ... % -2 * diag(D_An)*transpose(T(:,nt0).^2)*An*diag(a_An)*D_An*f_An ...
      +4 * D_An*diag(a_At)*At*mc_T(nt0)*An*diag(a_An)*D_An*f_An ... % +4 * D_An*transpose(An*diag(a_An))*diag(T(:,nt0).^2)*An*diag(a_An)*D_An*f_An ...
      -2 * D_An*a_An.^2*e_At*mc_T(nt0)*An*diag(a_An)*D_An*f_An ... % -2 * (D_An*a_An.^2)*transpose(T(:,nt0).^2)*An*diag(a_An)*D_An*f_An ...
      +1 * D_An*f_An*e_At*mc_T(nt0)*e_An*a_At.^2*D_An*f_An ... % + diag(D_An)*transpose(T(:,nt0).^2)*e_An*a_At.^2*D_An*f_An ...
      -2 * D_An*diag(a_At)*At*mc_T(nt0)*e_An*a_At.^2*D_An*f_An ... % -2 * D_An*transpose(An*diag(a_An))*diag(T(:,nt0).^2)*e_An*a_At.^2*D_An*f_An ...
      +1 * D_An*a_An.^2*e_At*mc_T(nt0)*e_An*a_At.^2*D_An*f_An ...% + (D_An*a_An.^2)*transpose(T(:,nt0).^2)*e_An*a_At.^2*D_An*f_An ...
      ) ...
      + ...
      (D_An.^2*( ...
      + f_An*e_At ...
      -4 * diag(a_At)*At ...
      +6 * a_An.^2*e_At ...
      -4 * diag(a_At.^3)*At ...
      + a_An.^4*e_At ...
      )*mc_T(nt0)*e_An);
    */
    if (F->nb1==F->nb2){
      ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
	ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
	na_j=0; while (na_j<M_At_nb1->rpop_j){
	  na_a = M_At_nb1->m_a_[na_j]; na_b = M_At_nb1->m_b_[na_j];
	  dtmp = 
	    (
	     +(1*A_ajdk[na_a/POPLENGTH + AJDK_0_2*A_pcols]*(M_An_nb1->rpop_j))
	     -(4*A_ajdk[na_a/POPLENGTH + AJDK_1_2*A_pcols]*(*L1_get(E_nb1->lf_et_An , na_j,na_b,na_a)))
	     +(6*A_ajdk[na_a/POPLENGTH + AJDK_2_2*A_pcols]*(M_An_nb1->rpop_j))
	     -(4*A_ajdk[na_a/POPLENGTH + AJDK_3_2*A_pcols]*(*L1_get(E_nb1->lf_et_An , na_j,na_b,na_a)))
	     +(1*A_ajdk[na_a/POPLENGTH + AJDK_4_2*A_pcols]*(M_An_nb1->rpop_j))
	     );
	  L2_plusequals(QC_AtTAnAtTAn , na_j,na_b,na_a , ns_j,ns_b,ns_a , +dtmp);
	  na_j++; /* while (na_j<M_At_nb1->rpop_j){ } */}
	ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
      if (verbose>3){ lfprintf(QC_AtTAnAtTAn," %% QC_AtTAnAtTAn: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
      /* if (F->nb1==F->nb2){ } */}
    /* if (QC_AtTAnAtTAn_bother){ } */}
  if (verbose>1){ printf(" %% [finished get_flattenloop] tidx %d\n",tidx);}
  return NULL;
}

void wrap_flattenloop(int *tidx,void **vpra,pthread_t *thread_in,struct bcc_double *F)
{
  /* This function calls get_flattenloop ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 2)
   */
  int verbose=0;
  int ip=0;
  if (verbose){ printf(" %% [entering wrap_flattenloop] tidx %d\n",*tidx);}
  ip=0;
  vpra[ip++] = tidx; vpra[ip++] = F;
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_flattenloop,vpra)){ printf("Warning! cannot create thread %d in wrap_flattenloop\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_flattenloop(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_flattenloop] tidx %d\n",*tidx);}
}

void bcc_flattenloop(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_double **F_ = D->F_; 
  int nbx=0,nb1=0,nb2=0; struct bcc_double *F=NULL;
  if (verbose){ printf(" %% [entering bcc_flattenloop]\n");}
  if (verbose){ printf(" %% calculating flattenloop.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx];
    GLOBAL_pthread_tic();
    wrap_flattenloop(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),F);
    GLOBAL_pthread_toc();
    /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% flattenloop: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished bcc_flattenloop]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_QR_AnZtSWnYt_error(int verbose,struct bcc_ajdk *D)
{
  int nbins=D->nbins; struct bcc_single **E_=D->E_; struct bcc_double **F_=D->F_;
  int nbx=0,nb1=0,nb2=0,length=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      length = E_nb1->A_nrows * D->T_ncols;
      printf(" %% nb1 %d nb2 %d length %d*%d=%d\n",nb1,nb2,E_nb1->A_nrows,D->T_ncols,length);
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% QR_AnZtSWnYt - QR_AnZtSWnYt_uu %0.16f\n",dra_diff(F->QR_AnZtSWnYt->lf,F->QR_AnZtSWnYt_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% QR_AnZtSZnAt - QR_AnZtSZnAt_uu %0.16f\n",dra_diff(F->QR_AnZtSZnAt->lf,F->QR_AnZtSZnAt_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% QR_AnAtTYnYt - QR_AnAtTYnYt_uu %0.16f\n",dra_diff(F->QR_AnAtTYnYt->lf,F->QR_AnAtTYnYt_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% QR_AnAtTAnAt - QR_AnAtTAnAt_uu %0.16f\n",dra_diff(F->QR_AnAtTAnAt->lf,F->QR_AnAtTAnAt_uu->lf,length,1));
	/* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
}

void bcc_QC_AtTYnWtSZn_error(int verbose,struct bcc_ajdk *D)
{
  int nbins=D->nbins; struct bcc_single **E_=D->E_; struct bcc_double **F_=D->F_;
  int nbx=0,nb1=0,nb2=0,length=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      length = D->A_ncols * D->T_ncols;
      printf(" %% nb1 %d nb2 %d length %d*%d=%d\n",nb1,nb2,D->A_ncols,D->T_ncols,length);
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% QC_AtTYnWtSZn - QC_AtTYnWtSZn_uu %0.16f\n",dra_diff(F->QC_AtTYnWtSZn->lf,F->QC_AtTYnWtSZn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% QC_AtTAnZtSZn - QC_AtTAnZtSZn_uu %0.16f\n",dra_diff(F->QC_AtTAnZtSZn->lf,F->QC_AtTAnZtSZn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% QC_AtTYnYtTAn - QC_AtTYnYtTAn_uu %0.16f\n",dra_diff(F->QC_AtTYnYtTAn->lf,F->QC_AtTYnYtTAn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% QC_AtTAnAtTAn - QC_AtTAnAtTAn_uu %0.16f\n",dra_diff(F->QC_AtTAnAtTAn->lf,F->QC_AtTAnAtTAn_uu->lf,length,1));
	/* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}} 
}

void bcc_flattenloop_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  bcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D);
  bcc_An_ajdk(D);
  bcc_lf_ZtSn(D);
  bcc_lf_ZtSWn(D);
  bcc_lf_YnWt(D);
  bcc_M_ZtSWn(D);
  bcc_M_YnWt(D);
  bcc_lf_AtTYn____WtSZn(D);
  bcc_lf_At_T_YnWt_S_Zn(D);
  if (error_check){ bcc_QC_AtTYnWtSZn_uu(D);}
  bcc_lf_AnZt_S_WnYt(D);
  bcc_lf_An_ZtSWn_Yt(D);
  if (error_check){ bcc_QR_AnZtSWnYt_uu(D);}
  bcc_singlestudy_ww(D);
  bcc_doublestudy_ww(D);
  bcc_flattenloop(D);
  if (error_check){ bcc_QC_AtTYnWtSZn_error(verbose,D);}
  if (error_check){ bcc_QR_AnZtSWnYt_error(verbose,D);}
  wkspace_printf();
}



