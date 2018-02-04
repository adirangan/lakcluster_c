void dcc_ajdk_copy(struct dcc_ajdk *D,struct dcc_ajdk *D_in)
{
  /* copy from D_in; assuming D has been preallocated. */
  int verbose=0;
  int nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0;
  if (verbose){ printf(" %% [entering dcc_ajdk_copy]\n");}
  D->A_ncols=D_in->A_ncols;D->T_ncols=D_in->T_ncols;
  /* copying A_bmc_b */
  D->A_cbother = (D->A_ncols>0);
  D->A_ncols_extend = (D->bitj - (D->A_ncols % D->bitj)) % D->bitj; 
  D->A_mc_length = bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% A_ncols %d, A_ncols_extend %d A_mc_length %d\n",D->A_ncols,D->A_ncols_extend,D->A_mc_length);}
  memcpy(D->A_bmc_b,D_in->A_bmc_b,D->A_mc_length);
  memcpy(D->A_bmc_j,D_in->A_bmc_j,D->A_mc_length);
  memcpy(D->A_bmc_j_rmv,D_in->A_bmc_j_rmv,D->A_mc_length);
  memcpy(D->A_bmc_j_rtn,D_in->A_bmc_j_rtn,D->A_mc_length);
  D->A_cpop_b= popcount_uchar_array(D->A_bmc_b,D->A_mc_length); D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->A_cpop_b %d D->A_cpop_j %d\n",D->A_cpop_b,D->A_cpop_j);}
  memcpy(D->A_umc_b,D_in->A_umc_b,D->A_ncols*sizeof(unsigned char));
  memcpy(D->A_umc_j,D_in->A_umc_j,D->A_ncols*sizeof(unsigned char));
  memcpy(D->A_umc_j_rmv,D_in->A_umc_j_rmv,D->A_ncols*sizeof(unsigned char));
  memcpy(D->A_umc_j_rtn,D_in->A_umc_j_rtn,D->A_ncols*sizeof(unsigned char));
  /* loading T_bmc_b */
  D->T_ncols_extend = (D->bitj - (D->T_ncols % D->bitj)) % D->bitj; 
  D->T_mc_length = bsize(D->T_ncols)/* rup(D->T_ncols+D->T_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% T_ncols %d, T_ncols_extend %d T_mc_length %d\n",D->T_ncols,D->T_ncols_extend,D->T_mc_length);}
  memcpy(D->T_bmc_b,D_in->T_bmc_b,D->T_mc_length);
  memcpy(D->T_bmc_j,D_in->T_bmc_j,D->T_mc_length);
  D->T_cpop_b= popcount_uchar_array(D->T_bmc_b,D->T_mc_length); D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length); if (verbose>1){ printf(" %%%% D->T_cpop_b %d D->T_cpop_j %d\n",D->T_cpop_b,D->T_cpop_j);}
  memcpy(D->T_umc_b,D_in->T_umc_b,D->T_ncols*sizeof(unsigned char));
  memcpy(D->T_umc_j,D_in->T_umc_j,D->T_ncols*sizeof(unsigned char));
  D->Irem=D_in->Irem;
  D->Ireq=D_in->Ireq;
  D->out_iteration=D_in->out_iteration;
  D->out_xdrop_ij=D_in->out_xdrop_ij;
  D->out_trace_length=D_in->out_trace_length;
  if (verbose){ printf(" %% [finished dcc_ajdk_copy] D->A_ncols=%d \n",D->A_ncols,D->T_ncols);}
}

void dcc_ajdk_load(struct dcc_ajdk *D)
{
  /* load from file */
  int verbose=0;
  int nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0;
  if (verbose){ printf(" %% [entering dcc_ajdk_load]\n");}
  D->A_ncols=GLOBAL_A_n_cols;D->T_ncols=GLOBAL_T_n_cols;
  /* loading A_bmc_b */
  D->A_cbother = (D->A_ncols>0);
  D->A_ncols_extend = (D->bitj - (D->A_ncols % D->bitj)) % D->bitj; 
  D->A_mc_length = bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% A_ncols %d, A_ncols_extend %d A_mc_length %d\n",D->A_ncols,D->A_ncols_extend,D->A_mc_length);}
  D->A_bmc_b = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_b){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} 
  D->A_bmc_j = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} 
  D->A_bmc_j_rmv = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rmv){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} 
  D->A_bmc_j_rtn = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rtn){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} 
  if (GLOBAL_A_n_cind==NULL || !strcmp(GLOBAL_A_n_cind,"\0")){ for (nc=0;nc<D->A_ncols;nc++){ bset__on(D->A_bmc_b,nc);}}
  else{ binary_read(GLOBAL_A_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->A_bmc_b)); if (nrows_tmp!=D->A_ncols){ printf(" %% Warning! A_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_A_n_name,GLOBAL_A_n_cind,nrows_tmp,D->A_ncols);}}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_b: "); if (verbose>2){ bprintf(D->A_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  for (nc=0;nc<D->A_mc_length;nc++){ D->A_bmc_j[nc] = D->A_bmc_b[nc];}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); if (verbose>2){ bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  D->A_cpop_b= popcount_uchar_array(D->A_bmc_b,D->A_mc_length); D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->A_cpop_b %d D->A_cpop_j %d\n",D->A_cpop_b,D->A_cpop_j);}
  D->A_umc_b = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_b){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_b[nc] = bget__on(D->A_bmc_b,nc);}
  D->A_umc_j = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);}
  D->A_umc_j_rmv = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rmv){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rmv[nc] = bget__on(D->A_bmc_j_rmv,nc);}
  D->A_umc_j_rtn = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rtn){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = bget__on(D->A_bmc_j_rtn,nc);}
  /* loading T_bmc_b */
  D->T_ncols_extend = (D->bitj - (D->T_ncols % D->bitj)) % D->bitj; 
  D->T_mc_length = bsize(D->T_ncols)/* rup(D->T_ncols+D->T_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% T_ncols %d, T_ncols_extend %d T_mc_length %d\n",D->T_ncols,D->T_ncols_extend,D->T_mc_length);}
  D->T_bmc_b = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_b){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} 
  D->T_bmc_j = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_j){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} 
  if (GLOBAL_T_n_cind==NULL || !strcmp(GLOBAL_T_n_cind,"\0")){ for (nc=0;nc<D->T_ncols;nc++){ bset__on(D->T_bmc_b,nc);}}
  else{ binary_read(GLOBAL_T_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->T_bmc_b)); if (nrows_tmp!=D->T_ncols){ printf(" %% Warning! T_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_T_n_name,GLOBAL_T_n_cind,nrows_tmp,D->T_ncols);}}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_b: "); if (verbose>2){ bprintf(D->T_bmc_b,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  for (nc=0;nc<D->T_mc_length;nc++){ D->T_bmc_j[nc] = D->T_bmc_b[nc];}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_j: "); if (verbose>2){ bprintf(D->T_bmc_j,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  D->T_cpop_b= popcount_uchar_array(D->T_bmc_b,D->T_mc_length); D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length);
  D->T_umc_b = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_b){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_b[nc] = bget__on(D->T_bmc_b,nc);}
  D->T_umc_j = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_j){ printf(" %% Warning! not enough memory in dcc_ajdk_load\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_j[nc] = bget__on(D->T_bmc_j,nc);}
  D->QC_TAnZtS_nrm = NULL; D->QC_TAnAtT_nrm = NULL;
  D->QC_sra=NULL; D->QC_lmc_a=NULL; D->QC_lmc_b=NULL; D->QC_lmc_j=NULL;
  D->QR_sra=NULL; D->QR_lmr_a=NULL; D->QR_lmr_b=NULL; D->QR_lmr_j=NULL; D->QR_lnb=NULL; D->QR_imr_a=NULL; D->QR_imr_b=NULL; D->QR_imr_j=NULL;
  D->Irem=0;D->Ireq=0; D->out_iteration=0; D->out_xdrop_ij=0; D->out_trace=NULL; D->out_xdrop_a=NULL; D->out_xdrop_b=NULL; D->out_trace_length = 0;
  if (verbose){ printf(" %% [finished dcc_ajdk_load] D->A_ncols=%d D->T_ncols=%d \n",D->A_ncols,D->T_ncols);}
}

void dcc_ajdk_init(double mrnd,struct dcc_ajdk *D)
{
  /* initialize based on test parameters */
  int verbose=0;
  int nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0;
  unsigned char *bb=NULL,*bj=NULL;
  double dtmp=0; int nx=0;
  if (verbose){ printf(" %% [entering dcc_ajdk_init]\n");}
  D->A_ncols=GLOBAL_TEST_A_n_cols;D->T_ncols=GLOBAL_TEST_T_n_cols;
  /* initializing A_bmc_b */
  D->A_cbother = (D->A_ncols>0);
  D->A_ncols_extend = (D->bitj - (D->A_ncols % D->bitj)) % D->bitj; 
  D->A_mc_length = bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% A_ncols %d, A_ncols_extend %d A_mc_length %d\n",D->A_ncols,D->A_ncols_extend,D->A_mc_length);}
  D->A_bmc_b = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_b){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} 
  D->A_bmc_j = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} 
  bb = D->A_bmc_b; bj = D->A_bmc_j; for (nx=0;nx<D->A_ncols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_b: "); if (verbose>2){ bprintf(D->A_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); if (verbose>2){ bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  D->A_bmc_j_rmv = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rmv){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} 
  D->A_bmc_j_rtn = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rtn){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} 
  D->A_cpop_b= popcount_uchar_array(D->A_bmc_b,D->A_mc_length); D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->A_cpop_b %d D->A_cpop_j %d\n",D->A_cpop_b,D->A_cpop_j);}
  D->A_umc_b = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_b){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_b[nc] = bget__on(D->A_bmc_b,nc);}
  D->A_umc_j = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);}
  D->A_umc_j_rmv = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rmv){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rmv[nc] = bget__on(D->A_bmc_j_rmv,nc);}
  D->A_umc_j_rtn = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rtn){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = bget__on(D->A_bmc_j_rtn,nc);}
  /* initializing T_bmc_b */
  D->T_ncols_extend = (D->bitj - (D->T_ncols % D->bitj)) % D->bitj; 
  D->T_mc_length = bsize(D->T_ncols)/* rup(D->T_ncols+D->T_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% T_ncols %d, T_ncols_extend %d T_mc_length %d\n",D->T_ncols,D->T_ncols_extend,D->T_mc_length);}
  D->T_bmc_b = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_b){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} 
  D->T_bmc_j = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_j){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} 
  bb = D->T_bmc_b; bj = D->T_bmc_j; for (nx=0;nx<D->T_ncols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}} bset__on(bb,0); bset__on(bj,0); bset__on(bb,D->T_ncols-1), bset__on(bj,D->T_ncols-1);
  sprintf(D->tmpTnchar," %%%% D->T_bmc_b: "); if (verbose>2){ bprintf(D->T_bmc_b,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_j: "); if (verbose>2){ bprintf(D->T_bmc_j,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  D->T_cpop_b= popcount_uchar_array(D->T_bmc_b,D->T_mc_length); D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length);
  D->T_umc_b = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_b){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_b[nc] = bget__on(D->T_bmc_b,nc);}
  D->T_umc_j = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_j){ printf(" %% Warning! not enough memory in dcc_ajdk_init\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_j[nc] = bget__on(D->T_bmc_j,nc);}
  D->QC_TAnZtS_nrm = NULL; D->QC_TAnAtT_nrm = NULL;
  D->QC_sra=NULL; D->QC_lmc_a=NULL; D->QC_lmc_b=NULL; D->QC_lmc_j=NULL;
  D->QR_sra=NULL; D->QR_lmr_a=NULL; D->QR_lmr_b=NULL; D->QR_lmr_j=NULL; D->QR_lnb=NULL; D->QR_imr_a=NULL; D->QR_imr_b=NULL; D->QR_imr_j=NULL;
  D->Irem=0;D->Ireq=0; D->out_iteration=0; D->out_xdrop_ij=0; D->out_trace=NULL; D->out_xdrop_a=NULL; D->out_xdrop_b=NULL; D->out_trace_length = 0;
  if (verbose){ printf(" %% [finished dcc_ajdk_init] D->A_ncols=%d D->T_ncols=%d \n",D->A_ncols,D->T_ncols);}
}

void dcc_single_copy_M_An(struct dcc_single *E,struct dcc_single *E_in)
{
  /* copy from file */
  int verbose=0;
  struct dcc_ajdk *D=E->D,*D_in=E_in->D;
  int nb_given = E->nb;
  int nr=0,nb=0,nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0,nrows_tmp_extend,brows_tmp;
  unsigned char *bXra_tmp=NULL;
  if (verbose){ printf(" %% [entering dcc_single_copy_M_An]\n");}
  /* copying A_bmr_b etc */
  nb = nb_given;
  /* copying A */
  D->bitj = D_in->bitj; D->A_ncols = D_in->A_ncols; E->A_nrows = E_in->A_nrows;
  E->A_rbother = (E->A_nrows>0);
  E->A_nrows_extend = (D->bitj - (E->A_nrows % D->bitj)) % D->bitj;
  E->A_mr_length = bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */;
  memcpy(E->A_bmr_b,E_in->A_bmr_b,E->A_mr_length);
  memcpy(E->A_bmr_j,E_in->A_bmr_j,E->A_mr_length);
  memcpy(E->A_bmr_j_rmv,E_in->A_bmr_j_rmv,E->A_mr_length);
  memcpy(E->A_bmr_j_rtn,E_in->A_bmr_j_rtn,E->A_mr_length);
  E->A_rpop_b = popcount_uchar_array(E->A_bmr_b,E->A_mr_length);
  E->A_rpop_j = popcount_uchar_array(E->A_bmr_j,E->A_mr_length);
  if (verbose>1){ printf(" %% A_nrows_extend_[%d] %d\n",nb,E->A_nrows_extend);}
  if (verbose>1){ printf(" %% A_mr_length_[%d] %d\n",nb,E->A_mr_length);}
  sprintf(D->tmpAnchar," %%%% E->A_bmr_b_[%.2d]:",nb); if (verbose>1){ bprintf(E->A_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% E->A_bmr_j_[%.2d]:",nb); if (verbose>1){ bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  if (verbose>1){ printf(" %% A_nrows_[%d]  = %d\n",nb,E->A_nrows);}
  if (verbose>1){ printf(" %% A_rpop_b_[%d] = %d\n",nb,E->A_rpop_b);}
  if (verbose>1){ printf(" %% A_rpop_j_[%d] = %d\n",nb,E->A_rpop_j);}
  memcpy(E->A_umr_b,E_in->A_umr_b,E->A_nrows*sizeof(unsigned char));
  memcpy(E->A_umr_j,E_in->A_umr_j,E->A_nrows*sizeof(unsigned char));
  memcpy(E->A_umr_j_rmv,E_in->A_umr_j_rmv,E->A_nrows*sizeof(unsigned char));
  memcpy(E->A_umr_j_rtn,E_in->A_umr_j_rtn,E->A_nrows*sizeof(unsigned char));
  if (verbose>1){ printf(" %% copying M_An_[%d]\n",nb);}
  M_handle_copy(E->M_An,E_in->M_An);
  M_handle_copy(E->M_At,E_in->M_At);
  /* copying Z */
  E->Z_nrows = E_in->Z_nrows;
  E->Z_rbother = (E->Z_nrows>0);
  E->Z_nrows_extend = (D->bitj - (E->Z_nrows % D->bitj)) % D->bitj;
  E->Z_mr_length = bsize(E->Z_nrows)/* rup(E->Z_nrows+E->Z_nrows_extend,POPLENGTH)/BIT8 */;
  memcpy(E->Z_bmr_b,E_in->Z_bmr_b,E->Z_mr_length);
  memcpy(E->Z_bmr_j,E_in->Z_bmr_j,E->Z_mr_length);
  E->Z_rpop_b = popcount_uchar_array(E->Z_bmr_b,E->Z_mr_length);
  E->Z_rpop_j = popcount_uchar_array(E->Z_bmr_j,E->Z_mr_length);
  if (verbose>1){ printf(" %% Z_nrows_extend_[%d] %d\n",nb,E->Z_nrows_extend);}
  if (verbose>1){ printf(" %% Z_mr_length_[%d] %d\n",nb,E->Z_mr_length);}
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_b_[%.2d]:",nb); if (verbose>1){ bprintf(E->Z_bmr_b,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_j_[%.2d]:",nb); if (verbose>1){ bprintf(E->Z_bmr_j,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  if (verbose>1){ printf(" %% Z_nrows_[%d]  = %d\n",nb,E->Z_nrows);}
  if (verbose>1){ printf(" %% Z_rpop_b_[%d] = %d\n",nb,E->Z_rpop_b);}
  if (verbose>1){ printf(" %% Z_rpop_j_[%d] = %d\n",nb,E->Z_rpop_j);}
  memcpy(E->Z_umr_b,E_in->Z_umr_b,E->Z_nrows*sizeof(unsigned char));
  memcpy(E->Z_umr_j,E_in->Z_umr_j,E->Z_nrows*sizeof(unsigned char));
  if (verbose>1){ printf(" %% copying M_Zn_[%d]\n",nb);}
  M_handle_copy(E->M_Zn,E_in->M_Zn);
  M_handle_copy(E->M_Zt,E_in->M_Zt);
  /* copying T */
  D->T_ncols = D_in->T_ncols;
  M_handle_copy(E->M_Tn,E_in->M_Tn);
  M_handle_copy(E->M_Tt,E_in->M_Tt);
  /* copying S */
  M_handle_copy(E->M_Sn,E_in->M_Sn);
  M_handle_copy(E->M_St,E_in->M_St);
  if (verbose){ printf(" %% [finished dcc_single_copy_M_An]\n");}
}

void dcc_single_load_M_An(struct dcc_single *E)
{
  /* load from file */
  int verbose=0;
  struct dcc_ajdk *D=E->D;
  int nb_given = E->nb;
  int nr=0,nb=0,nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0,nrows_tmp_extend,brows_tmp;
  unsigned char *bXra_tmp=NULL;
  if (verbose){ printf(" %% [entering dcc_single_load_M_An]\n");}
  /* loading A_bmr_b etc */
  nb = nb_given;
  /* loading A */
  binary_read_getsize(GLOBAL_A_t_name_[nb],&(D->bitj),&(D->A_ncols),&(E->A_nrows));
  E->A_rbother = (E->A_nrows>0);
  E->A_nrows_extend = (D->bitj - (E->A_nrows % D->bitj)) % D->bitj;
  E->A_mr_length = bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */;
  E->A_bmr_b = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_b){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");}
  E->A_bmr_j = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");}
  E->A_bmr_j_rmv = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rmv){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");}
  E->A_bmr_j_rtn = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rtn){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");}
  if (GLOBAL_A_n_rind_[nb]==NULL || !strcmp(GLOBAL_A_n_rind_[nb],"\0")){ for (nr=0;nr<E->A_nrows;nr++){ bset__on(E->A_bmr_b,nr);}}
  else{ binary_read(GLOBAL_A_n_rind_[nb],&bitj_tmp,&nrows_tmp,&ncols_tmp,&(E->A_bmr_b)); if (nrows_tmp!=E->A_nrows || bitj_tmp!=D->bitj){ printf(" %% Warning! A_: %s; improper mr_b_ %s, nrows %d instead of %d, bitj %d vs %d\n",GLOBAL_A_t_name_[nb],GLOBAL_A_n_rind_[nb],nrows_tmp,E->A_nrows,bitj_tmp,D->bitj);}}
  for (nr=0;nr<E->A_mr_length;nr++){ E->A_bmr_j[nr] = E->A_bmr_b[nr];}
  E->A_rpop_b = popcount_uchar_array(E->A_bmr_b,E->A_mr_length);
  E->A_rpop_j = popcount_uchar_array(E->A_bmr_j,E->A_mr_length);
  if (verbose>1){ printf(" %% reading GLOBAL_A_t_name_[%d] = %s\n",nb,GLOBAL_A_t_name_[nb]);}
  if (verbose>1){ printf(" %% read GLOBAL_A_t_name_[%d] %d-x-%d (bitj %d)\n",nb,D->A_ncols,E->A_nrows,D->bitj);}
  if (verbose>1){ printf(" %% A_nrows_extend_[%d] %d\n",nb,E->A_nrows_extend);}
  if (verbose>1){ printf(" %% A_mr_length_[%d] %d\n",nb,E->A_mr_length);}
  if (verbose>1){ printf(" %% reading GLOBAL_A_n_rind_[%d] = %s\n",nb,GLOBAL_A_n_rind_[nb]);} 
  sprintf(D->tmpAnchar," %%%% E->A_bmr_b_[%.2d]:",nb); if (verbose>1){ bprintf(E->A_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% E->A_bmr_j_[%.2d]:",nb); if (verbose>1){ bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  if (verbose>1){ printf(" %% A_nrows_[%d]  = %d\n",nb,E->A_nrows);}
  if (verbose>1){ printf(" %% A_rpop_b_[%d] = %d\n",nb,E->A_rpop_b);}
  if (verbose>1){ printf(" %% A_rpop_j_[%d] = %d\n",nb,E->A_rpop_j);}
  E->A_umr_b = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_b){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_b[nr] = bget__on(E->A_bmr_b,nr);}
  E->A_umr_j = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);}
  E->A_umr_j_rmv = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rmv){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rmv[nr] = bget__on(E->A_bmr_j_rmv,nr);}
  E->A_umr_j_rtn = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rtn){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = bget__on(E->A_bmr_j_rtn,nr);}
  if (verbose>1){ printf(" %% generating M_An_[%d]\n",nb);}
  E->M_An = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,GLOBAL_A_t_name_[nb],NULL,E->A_bmr_b,D->A_bmc_b); M_mxget(E->M_An);
  E->M_At = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,GLOBAL_A_n_name_[nb],NULL,D->A_bmc_b,E->A_bmr_b); M_mxget(E->M_At);
  /* loading Z */
  if (GLOBAL_Z_t_name_[nb]==NULL || !strcmp(GLOBAL_Z_t_name_[nb],"\0")){ bitj_tmp=D->bitj;ncols_tmp=D->A_ncols;E->Z_nrows=0;}
  else{ binary_read_getsize(GLOBAL_Z_t_name_[nb],&bitj_tmp,&ncols_tmp,&(E->Z_nrows));}
  if (bitj_tmp!=D->bitj || ncols_tmp!=D->A_ncols){ printf(" %% Warning! Z_: %s; improper bitj %d vs %d, ncols %d vs %d\n",GLOBAL_Z_t_name_[nb],bitj_tmp,D->bitj,ncols_tmp,D->A_ncols);} 
  E->Z_rbother = (E->Z_nrows>0);
  E->Z_nrows_extend = (D->bitj - (E->Z_nrows % D->bitj)) % D->bitj;
  E->Z_mr_length = bsize(E->Z_nrows)/* rup(E->Z_nrows+E->Z_nrows_extend,POPLENGTH)/BIT8 */;
  E->Z_bmr_b = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_b){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");}
  E->Z_bmr_j = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_j){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");}
  if (GLOBAL_Z_n_rind_[nb]==NULL || !strcmp(GLOBAL_Z_n_rind_[nb],"\0")){ for (nr=0;nr<E->Z_nrows;nr++){ bset__on(E->Z_bmr_b,nr);}}
  else{ binary_read(GLOBAL_Z_n_rind_[nb],&bitj_tmp,&nrows_tmp,&ncols_tmp,&(E->Z_bmr_b)); if (nrows_tmp!=E->Z_nrows || bitj_tmp!=D->bitj){ printf(" %% Warning! Z_: %s; improper mr_b_ %s, nrows %d instead of %d, bitj %d vs %d\n",GLOBAL_Z_t_name_[nb],GLOBAL_Z_n_rind_[nb],nrows_tmp,E->Z_nrows,bitj_tmp,D->bitj);}}
  for (nr=0;nr<E->Z_mr_length;nr++){ E->Z_bmr_j[nr] = E->Z_bmr_b[nr];}
  E->Z_rpop_b = popcount_uchar_array(E->Z_bmr_b,E->Z_mr_length);
  E->Z_rpop_j = popcount_uchar_array(E->Z_bmr_j,E->Z_mr_length);
  if (verbose>1){ printf(" %% reading GLOBAL_Z_t_name_[%d] = %s\n",nb,GLOBAL_Z_t_name_[nb]);}
  if (verbose>1){ printf(" %% read GLOBAL_Z_t_name_[%d] %d-x-%d (bitj %d)\n",nb,D->A_ncols,E->Z_nrows,D->bitj);}
  if (verbose>1){ printf(" %% Z_nrows_extend_[%d] %d\n",nb,E->Z_nrows_extend);}
  if (verbose>1){ printf(" %% Z_mr_length_[%d] %d\n",nb,E->Z_mr_length);}
  if (verbose>1){ printf(" %% reading Z_n_rind_[%d] = %s\n",nb,GLOBAL_Z_n_rind_[nb]);} 
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_b_[%.2d]:",nb); if (verbose>1){ bprintf(E->Z_bmr_b,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_j_[%.2d]:",nb); if (verbose>1){ bprintf(E->Z_bmr_j,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  if (verbose>1){ printf(" %% Z_nrows_[%d]  = %d\n",nb,E->Z_nrows);}
  if (verbose>1){ printf(" %% Z_rpop_b_[%d] = %d\n",nb,E->Z_rpop_b);}
  if (verbose>1){ printf(" %% Z_rpop_j_[%d] = %d\n",nb,E->Z_rpop_j);}
  E->Z_umr_b = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_b){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_b[nr] = bget__on(E->Z_bmr_b,nr);}
  E->Z_umr_j = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_j){ printf(" %% Warning! not enough memory in dcc_single_load_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_j[nr] = bget__on(E->Z_bmr_j,nr);}
  if (verbose>1){ printf(" %% generating M_Zn_[%d]\n",nb);}
  if (GLOBAL_Z_t_name_[nb]==NULL || !strcmp(GLOBAL_Z_t_name_[nb],"\0")){ E->M_Zn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,NULL,NULL,E->Z_bmr_b,D->A_bmc_b);}
  else{ E->M_Zn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,GLOBAL_Z_t_name_[nb],NULL,E->Z_bmr_b,D->A_bmc_b);} M_mxget(E->M_Zn);
  if (GLOBAL_Z_n_name_[nb]==NULL || !strcmp(GLOBAL_Z_n_name_[nb],"\0")){ E->M_Zt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,NULL,NULL,D->A_bmc_b,E->Z_bmr_b);}
  else{ E->M_Zt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,GLOBAL_Z_n_name_[nb],NULL,D->A_bmc_b,E->Z_bmr_b);} M_mxget(E->M_Zt);
  /* loading T */
  if (GLOBAL_T_t_name_[nb]==NULL || !strcmp(GLOBAL_T_t_name_[nb],"\0")){ bitj_tmp=D->bitj;D->T_ncols=1;nrows_tmp=E->A_nrows;}
  else{ binary_read_getsize(GLOBAL_T_t_name_[nb],&bitj_tmp,&(D->T_ncols),&nrows_tmp);}
  if (bitj_tmp!=D->bitj || nrows_tmp!=E->A_nrows){ printf(" %% Warning! T_: %s; improper bitj %d vs %d, nrows %d vs %d\n",GLOBAL_T_t_name_[nb],bitj_tmp,D->bitj,nrows_tmp,E->A_nrows);}
  if (GLOBAL_T_t_name_[nb]==NULL || !strcmp(GLOBAL_T_t_name_[nb],"\0")){ 
    bitj_tmp = BITJ; nrows_tmp = 1; ncols_tmp = E->A_nrows; nrows_tmp_extend = (bitj_tmp - (nrows_tmp % bitj_tmp)) % bitj_tmp; brows_tmp = bsize(nrows_tmp)/* (rup(nrows_tmp + nrows_tmp_extend,POPLENGTH))/BIT8 */;
    bXra_tmp = wkspace_all0c(ncols_tmp*brows_tmp); for (nc=0;nc<ncols_tmp;nc++){ bXra_tmp[nc*brows_tmp] |=  (1 << 7);}
    E->M_Tn = M_handle_v_make(D->bitj,E->A_nrows,D->T_ncols,NULL,bXra_tmp,E->A_bmr_b,D->T_bmc_b);}
  else{ E->M_Tn = M_handle_v_make(D->bitj,E->A_nrows,D->T_ncols,GLOBAL_T_t_name_[nb],NULL,E->A_bmr_b,D->T_bmc_b);} M_mxget(E->M_Tn);
  if (GLOBAL_T_n_name_[nb]==NULL || !strcmp(GLOBAL_T_n_name_[nb],"\0")){ 
    bitj_tmp = BITJ; nrows_tmp = E->A_nrows; ncols_tmp = 1; nrows_tmp_extend = (bitj_tmp - (nrows_tmp % bitj_tmp)) % bitj_tmp; brows_tmp = bsize(nrows_tmp)/* (rup(nrows_tmp + nrows_tmp_extend,POPLENGTH))/BIT8 */;
    bXra_tmp = wkspace_all0c(ncols_tmp*brows_tmp); for (nr=0;nr<nrows_tmp;nr++){ bset__on(bXra_tmp,nr);}
    E->M_Tt = M_handle_v_make(D->bitj,D->T_ncols,E->A_nrows,NULL,bXra_tmp,D->T_bmc_b,E->A_bmr_b);}
  else{ E->M_Tt = M_handle_v_make(D->bitj,D->T_ncols,E->A_nrows,GLOBAL_T_n_name_[nb],NULL,D->T_bmc_b,E->A_bmr_b);} M_mxget(E->M_Tt);
  /* loading S */
  if (GLOBAL_S_t_name_[nb]==NULL || !strcmp(GLOBAL_S_t_name_[nb],"\0")){ bitj_tmp=D->bitj;ncols_tmp=D->T_ncols;nrows_tmp=E->Z_nrows;}
  else{ binary_read_getsize(GLOBAL_S_t_name_[nb],&bitj_tmp,&ncols_tmp,&nrows_tmp);}
  if (bitj_tmp!=D->bitj || ncols_tmp!=D->T_ncols || nrows_tmp!=E->Z_nrows){ printf(" %% Warning! S_: %s; improper bitj %d vs %d, ncols %d vs %d, nrows %d vs %d\n",GLOBAL_S_t_name_[nb],bitj_tmp,D->bitj,ncols_tmp,D->T_ncols,nrows_tmp,E->Z_nrows);}
  if (GLOBAL_S_t_name_[nb]==NULL || !strcmp(GLOBAL_S_t_name_[nb],"\0")){ 
    bitj_tmp = BITJ; nrows_tmp = 1; ncols_tmp = E->Z_nrows; nrows_tmp_extend = (bitj_tmp - (nrows_tmp % bitj_tmp)) % bitj_tmp; brows_tmp = bsize(nrows_tmp)/* (rup(nrows_tmp + nrows_tmp_extend,POPLENGTH))/BIT8 */;
    bXra_tmp = wkspace_all0c(ncols_tmp*brows_tmp); for (nc=0;nc<ncols_tmp;nc++){ bXra_tmp[nc*brows_tmp] |=  (1 << 7);}
    E->M_Sn = M_handle_v_make(D->bitj,E->Z_nrows,D->T_ncols,NULL,bXra_tmp,E->Z_bmr_b,D->T_bmc_b);}
  else{ E->M_Sn = M_handle_v_make(D->bitj,E->Z_nrows,D->T_ncols,GLOBAL_S_t_name_[nb],NULL,E->Z_bmr_b,D->T_bmc_b);} M_mxget(E->M_Sn);
  if (GLOBAL_S_n_name_[nb]==NULL || !strcmp(GLOBAL_S_n_name_[nb],"\0")){ 
    bitj_tmp = BITJ; nrows_tmp = E->Z_nrows; ncols_tmp = 1; nrows_tmp_extend = (bitj_tmp - (nrows_tmp % bitj_tmp)) % bitj_tmp; brows_tmp = bsize(nrows_tmp)/* (rup(nrows_tmp + nrows_tmp_extend,POPLENGTH))/BIT8 */;
    bXra_tmp = wkspace_all0c(ncols_tmp*brows_tmp); for (nr=0;nr<nrows_tmp;nr++){ bset__on(bXra_tmp,nr);}
    E->M_St = M_handle_v_make(D->bitj,D->T_ncols,E->Z_nrows,NULL,bXra_tmp,D->T_bmc_b,E->Z_bmr_b);}
  else{ E->M_St = M_handle_v_make(D->bitj,D->T_ncols,E->Z_nrows,GLOBAL_S_n_name_[nb],NULL,D->T_bmc_b,E->Z_bmr_b);} M_mxget(E->M_St);
  E->QR_TAnZtS_nrm = NULL; E->QR_TAnAtT_nrm = NULL; E->QR_imr_a=NULL; E->QR_imr_b=NULL;
  if (verbose){ printf(" %% [finished dcc_single_load_M_An]\n");}
}

void dcc_single_init_M_An(char *error_vs_speed,double mrnd,struct dcc_single *E)
{
  /* initialize based on test parameters */
  int verbose=0;
  struct dcc_ajdk *D=E->D;
  int nb_given = E->nb;
  int nr=0,nb=0,nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0,nrows_tmp_extend,brows_tmp;
  unsigned char *bXra_tmp=NULL;
  unsigned char *bb=NULL,*bj=NULL,*bn=NULL,*bt=NULL,*bn2=NULL,*bt2=NULL;
  double dtmp=0; int nx=0;
  if (verbose){ printf(" %% [entering dcc_single_init_M_An]\n");}
  nb = nb_given;
  /* initializing A */
  E->A_nrows = maximum(0,GLOBAL_TEST_A_n_rows + (nb_given % 3) - 1);
  E->A_rbother = (E->A_nrows>0);
  E->A_nrows_extend = (D->bitj - (E->A_nrows % D->bitj)) % D->bitj;
  E->A_mr_length = bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */;
  E->A_bmr_b = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_b){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");}
  E->A_bmr_j = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");}
  bb = E->A_bmr_b; bj = E->A_bmr_j; for (nx=0;nx<E->A_nrows;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpAnchar," %%%% E->A_bmr_b: "); if (verbose>2){ bprintf(E->A_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% E->A_bmr_j: "); if (verbose>2){ bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  E->A_bmr_j_rmv = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rmv){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");}
  E->A_bmr_j_rtn = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rtn){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");}
  E->A_rpop_b = popcount_uchar_array(E->A_bmr_b,E->A_mr_length);
  E->A_rpop_j = popcount_uchar_array(E->A_bmr_j,E->A_mr_length);
  E->A_umr_b = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_b){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_b[nr] = bget__on(E->A_bmr_b,nr);}
  E->A_umr_j = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);}
  E->A_umr_j_rmv = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rmv){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rmv[nr] = bget__on(E->A_bmr_j_rmv,nr);}
  E->A_umr_j_rtn = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rtn){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = bget__on(E->A_bmr_j_rtn,nr);}
  if (verbose>1){ printf(" %% generating E->M_An\n");}
  bn = NULL; bt = NULL; E->M_An = NULL; E->M_At = NULL;
  wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->A_ncols,E->A_bmr_b,D->A_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->A_ncols,bn,bt,E->A_bmr_b,E->A_bmr_j,D->A_bmc_b,D->A_bmc_j,&(E->M_An),&(E->M_At));
  /* initializing Z */
  E->Z_nrows = maximum(0,GLOBAL_TEST_Z_n_rows + (nb_given % 5) - 2);
  E->Z_rbother = (E->Z_nrows>0);
  E->Z_nrows_extend = (D->bitj - (E->Z_nrows % D->bitj)) % D->bitj;
  E->Z_mr_length = bsize(E->Z_nrows)/* rup(E->Z_nrows+E->Z_nrows_extend,POPLENGTH)/BIT8 */;
  E->Z_bmr_b = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_b){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");}
  E->Z_bmr_j = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_j){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");}
  bb = E->Z_bmr_b; bj = E->Z_bmr_j; for (nx=0;nx<E->Z_nrows;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_b: "); if (verbose>2){ bprintf(E->Z_bmr_b,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_j: "); if (verbose>2){ bprintf(E->Z_bmr_j,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  E->Z_rpop_b = popcount_uchar_array(E->Z_bmr_b,E->Z_mr_length);
  E->Z_rpop_j = popcount_uchar_array(E->Z_bmr_j,E->Z_mr_length);
  E->Z_umr_b = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_b){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_b[nr] = bget__on(E->Z_bmr_b,nr);}
  E->Z_umr_j = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_j){ printf(" %% Warning! not enough memory in dcc_single_init_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_j[nr] = bget__on(E->Z_bmr_j,nr);}
  if (verbose>1){ printf(" %% generating E->M_Zn\n");}
  bn = NULL; bt = NULL; E->M_Zn = NULL; E->M_Zt = NULL;
  wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->Z_nrows,D->A_ncols,E->Z_bmr_b,D->A_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->Z_nrows,D->A_ncols,bn,bt,E->Z_bmr_b,E->Z_bmr_j,D->A_bmc_b,D->A_bmc_j,&(E->M_Zn),&(E->M_Zt));
  /* initializing T */
  bn = NULL; bt = NULL; E->M_Tn = NULL; E->M_Tt = NULL;
  wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->T_ncols,E->A_bmr_b,D->T_bmc_b,&(bn),&(bt));
  for (nx=0;nx<E->A_rpop_b;nx++){ bn2 = &(bn[nx*bsize(D->T_ncols)]); bset__on(bn2,0);}
  for (nx=0;nx<E->A_nrows;nx++){ bt2 = &(bt[0*bsize(E->A_nrows)]); bset__on(bt2,nx);}
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->T_ncols,bn,bt,E->A_bmr_b,E->A_bmr_j,D->T_bmc_b,D->T_bmc_j,&(E->M_Tn),&(E->M_Tt));
  /* initializing S */
  bn = NULL; bt = NULL; E->M_Sn = NULL; E->M_St = NULL;
  wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->Z_nrows,D->T_ncols,E->Z_bmr_b,D->T_bmc_b,&(bn),&(bt));
  for (nx=0;nx<E->Z_rpop_b;nx++){ bn2 = &(bn[nx*bsize(D->T_ncols)]); bset__on(bn2,0);}
  for (nx=0;nx<E->Z_nrows;nx++){ bt2 = &(bt[0*bsize(E->Z_nrows)]); bset__on(bt2,nx);}
  wrap_M_setup_test_excerpt_2(E->Z_nrows,D->T_ncols,bn,bt,E->Z_bmr_b,E->Z_bmr_j,D->T_bmc_b,D->T_bmc_j,&(E->M_Sn),&(E->M_St));
  /* initializing E->QR_TAnZtS_nrm */
  E->QR_TAnZtS_nrm = NULL; E->QR_TAnAtT_nrm = NULL; E->QR_imr_a=NULL; E->QR_imr_b=NULL;
  if (verbose){ printf(" %% [finished dcc_single_init_M_An]\n");}
}

void dcc_single_copy_lf(struct dcc_single *E,struct dcc_single *E_in)
{
  int verbose=0; /* double gamma = GLOBAL_gamma; */
  struct dcc_ajdk *D = E->D;
  int ns_b=0,n1a=0,n1b=0,n2a=0,n2b=0;
  char prefix[FNAMESIZE];
  int M_flag=0;
  int cdrop=0,ckeep=0;
  long long int mm=0;
  L_handle_copy(E->lf_At_rsum,E_in->lf_At_rsum);
  L_handle_copy(E->lf_Zt_rsum,E_in->lf_Zt_rsum);
  L_handle_copy(E->lf_AtTn,E_in->lf_AtTn);
  L_handle_copy(E->lf_ZtSn,E_in->lf_ZtSn);
  L_handle_copy(E->lf_An_ajdk,E_in->lf_An_ajdk);
  L_handle_copy(E->lf_Zn_ajdk,E_in->lf_Zn_ajdk);
  L_handle_copy(E->lf_a0d1,E_in->lf_a0d1);
  L_handle_copy(E->lf_a2d1,E_in->lf_a2d1);
  L_handle_copy(E->lf_a1d1_ZtSn,E_in->lf_a1d1_ZtSn);
  L_handle_copy(E->lf_a1d1_AtTn,E_in->lf_a1d1_AtTn);
  L_handle_copy(E->lf_et_Sn,E_in->lf_et_Sn);
  L_handle_copy(E->lf_et_Tn,E_in->lf_et_Tn);
  L_handle_copy(E->lf_a0d1_ZtSn,E_in->lf_a0d1_ZtSn); M_handle_copy(E->M_a0d1_ZtSn,E_in->M_a0d1_ZtSn);
  L_handle_copy(E->lf_a0d1_AtTn,E_in->lf_a0d1_AtTn); M_handle_copy(E->M_a0d1_AtTn,E_in->M_a0d1_AtTn);
}

void dcc_single_init_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,struct M_handle **M_p)
{
  int verbose=0;
  int length=0;
  length  = n1b*n2b; *L_p = L_handle_make(length); if (!(*L_p)){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  if (M_flag){ *M_p = M_handle_w_make(BITJ,GLOBAL_B_MLT,n1a,n2b);}
  if (verbose){ printf(" %% %s n1a %d n1b %d n2a %d n2b %d --> length %d \n",prefix,n1a,n1b,n2a,n2b,length);}
}

void dcc_single_init_lf(struct dcc_single *E)
{
  int verbose=0; /* double gamma = GLOBAL_gamma; */
  struct dcc_ajdk *D = E->D;
  int ns_b=0,n1a=0,n1b=0,n2a=0,n2b=0;
  char prefix[FNAMESIZE];
  int M_flag=0;
  int cdrop=0,ckeep=0;
  long long int mm=0;
  E->lf_At_rsum=NULL; E->lf_Zt_rsum=NULL;
  E->length = D->A_ncols; E->lf_At_rsum = L_handle_make(E->length); if (!E->lf_At_rsum){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  E->length = D->A_ncols; E->lf_Zt_rsum = L_handle_make(E->length); if (!E->lf_Zt_rsum){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  if (verbose){ printf(" %% initializing E->lf_AtTn, E->lf_ZtSn.\n");}
  E->lf_AtTn = NULL; E->lf_ZtSn = NULL;
  E->length = D->A_ncols*D->T_ncols; E->lf_AtTn = L_handle_make(E->length); if (!E->lf_AtTn){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  E->length = D->A_ncols*D->T_ncols; E->lf_ZtSn = L_handle_make(E->length); if (!E->lf_ZtSn){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing E->lf_An_ajdk, etc\n");}
  E->lf_An_ajdk = NULL; E->lf_Zn_ajdk = NULL;
  E->length = E->A_nrows*AJDK_TOT; E->lf_An_ajdk = L_handle_make(E->length); if (!E->lf_An_ajdk){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  E->length = E->Z_nrows*AJDK_TOT; E->lf_Zn_ajdk = L_handle_make(E->length); if (!E->lf_Zn_ajdk){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing E->lf_a0d1, etc\n");}
  E->length = 1; E->lf_a0d1 = L_handle_make(E->length);
  E->length = 1; E->lf_a2d1 = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_a1d1_ZtSn = L_handle_make(E->length); if (!E->lf_a1d1_ZtSn){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_a1d1_AtTn = L_handle_make(E->length); if (!E->lf_a1d1_AtTn){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_et_Sn = L_handle_make(E->length); if (!E->lf_et_Sn){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_et_Tn = L_handle_make(E->length); if (!E->lf_et_Tn){ printf(" %% Warning! not enough memory in dcc_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing storage for doublestudy\n");}
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_ZtSn");
  n1a = D->A_cbother*D->A_ncols; n1b = D->A_cbother*D->A_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcc_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_ZtSn),M_flag,&(E->M_a0d1_ZtSn));
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_AtTn");
  n1a = D->A_cbother*D->A_ncols; n1b = D->A_cbother*D->A_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcc_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_AtTn),M_flag,&(E->M_a0d1_AtTn));
}

void dcc_double_copy_lf(struct dcc_double *F,struct dcc_double *F_in)
{
  int verbose=0;
  struct dcc_ajdk *D=F->D;
  int n1a=0,n2a=0,n1b,n2b=0;
  int M_flag=0,rdrop=0,rdrop2=0;
  char prefix[FNAMESIZE];
  GLOBAL_tic(0);
  if (verbose>1){ printf(" %% dcc_double_copy_lf: %d,%d\n",F->nb1,F->nb2);}
  L_handle_copy(F->lf_An_a0d1_ZtSn,F_in->lf_An_a0d1_ZtSn);
  L_handle_copy(F->lf_An_a0d1_AtTn,F_in->lf_An_a0d1_AtTn);
  L_handle_copy(F->lf_TAnZtS_ww,F_in->lf_TAnZtS_ww);
  L_handle_copy(F->lf_TAnAtT_ww,F_in->lf_TAnAtT_ww);
  L_handle_copy(F->lf_TAnZtS_uu,F_in->lf_TAnZtS_uu);
  L_handle_copy(F->lf_TAnAtT_uu,F_in->lf_TAnAtT_uu);
  L_handle_copy(F->lf_D_AtTn_ZtSn_vv,F_in->lf_D_AtTn_ZtSn_vv);
  L_handle_copy(F->lf_D_AtTn_AtTn_vv,F_in->lf_D_AtTn_AtTn_vv);
  L_handle_copy(F->lf_D_AtTn_ZtSn_uu,F_in->lf_D_AtTn_ZtSn_uu);
  L_handle_copy(F->lf_D_AtTn_AtTn_uu,F_in->lf_D_AtTn_AtTn_uu);
  L_handle_copy(F->QR_TAnZtS,F_in->QR_TAnZtS);
  L_handle_copy(F->QR_TAnAtT,F_in->QR_TAnAtT);
  L_handle_copy(F->QC_TAnZtS,F_in->QC_TAnZtS);
  L_handle_copy(F->QC_TAnAtT,F_in->QC_TAnAtT);
  L_handle_copy(F->QR_TAnZtS_uu,F_in->QR_TAnZtS_uu);
  L_handle_copy(F->QR_TAnAtT_uu,F_in->QR_TAnAtT_uu);
  L_handle_copy(F->QC_TAnZtS_uu,F_in->QC_TAnZtS_uu);
  L_handle_copy(F->QC_TAnAtT_uu,F_in->QC_TAnAtT_uu);
  GLOBAL_toc(0,verbose," %% copying data_structures: ");
}

void dcc_double_init_lf(struct dcc_double *F)
{
  int verbose=0;
  struct dcc_ajdk *D=F->D;
  int n1a=0,n2a=0,n1b,n2b=0;
  int M_flag=0,rdrop=0,rdrop2=0;
  char prefix[FNAMESIZE];
  GLOBAL_tic(0);
  if (verbose>1){ printf(" %% dcc_double_init_lf: %d,%d\n",F->nb1,F->nb2);}
  if (verbose>1){ printf(" %% initializing lf_ \n");}
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_An_a0d1_ZtSn = L_handle_make(F->length); if (!F->lf_An_a0d1_ZtSn){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_An_a0d1_AtTn = L_handle_make(F->length); if (!F->lf_An_a0d1_AtTn){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_TAnZtS_ww = L_handle_make(F->length); if (!F->lf_TAnZtS_ww){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_TAnAtT_ww = L_handle_make(F->length); if (!F->lf_TAnAtT_ww){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_TAnZtS_uu = L_handle_make(F->length); if (!F->lf_TAnZtS_uu){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_TAnAtT_uu = L_handle_make(F->length); if (!F->lf_TAnAtT_uu){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  F->length = D->A_ncols*D->T_ncols; F->lf_D_AtTn_ZtSn_vv = L_handle_make(F->length); if (!F->lf_D_AtTn_ZtSn_vv){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  F->length = D->A_ncols*D->T_ncols; F->lf_D_AtTn_AtTn_vv = L_handle_make(F->length); if (!F->lf_D_AtTn_AtTn_vv){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  F->length = D->A_ncols*D->T_ncols; F->lf_D_AtTn_ZtSn_uu = L_handle_make(F->length); if (!F->lf_D_AtTn_ZtSn_uu){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  F->length = D->A_ncols*D->T_ncols; F->lf_D_AtTn_AtTn_uu = L_handle_make(F->length); if (!F->lf_D_AtTn_AtTn_uu){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  if (verbose>1){ printf(" %% initializing F->QX\n");}
  F->QR_TAnZtS=NULL; F->QR_TAnAtT=NULL; F->QC_TAnZtS=NULL; F->QC_TAnAtT=NULL;
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_TAnZtS = L_handle_make(F->length);
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_TAnAtT = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_TAnZtS = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_TAnAtT = L_handle_make(F->length); if (!F->QC_TAnAtT){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  if (verbose>1){ printf(" %% initializing F->QX_uu\n");}
  F->QR_TAnZtS_uu=NULL; F->QR_TAnAtT_uu=NULL; F->QC_TAnZtS_uu=NULL; F->QC_TAnAtT_uu=NULL;
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_TAnZtS_uu = L_handle_make(F->length);
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_TAnAtT_uu = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_TAnZtS_uu = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_TAnAtT_uu = L_handle_make(F->length); if (!F->QC_TAnAtT_uu){ printf(" %% Warning! not enough memory in dcc_double_init_lf\n");}
  GLOBAL_toc(0,verbose," %% initializing data_structures: ");
}

void dcc_copy_A_p(struct dcc_ajdk *D,struct dcc_ajdk *D_in)
{
  /* copies A_p for each column from D_in */
  int verbose=0; char tempchar[FNAMESIZE];
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_single **E_in_ = D_in->E_;
  int nb=0,nc=0,nc_p=0,nc_start=0,nc_final=0;
  int A_pcols=0,A_ncols=0,A_rpop_b_total=0,Z_rpop_b_total=0;
  double tmp_sum=0,tmp_num=0;
  if (verbose){ printf(" %% [entering dcc_copy_A_p]\n"); wkspace_printf();}
  if (verbose){ printf(" %% \n");}
  D->A_pcols = psize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/POPLENGTH */;
  if (verbose){ printf(" %% nbins %d; A_ncols %d A_pcols %d\n",nbins,D->A_ncols,D->A_pcols);}
  memcpy(D->AZ_rsum,D_in->AZ_rsum,D->A_ncols*sizeof(double)); 
  memcpy(D->A_p,D_in->A_p,D->A_pcols*sizeof(double)); 
  memcpy(D->A_ajdk,D_in->A_ajdk,AJDK_TOT*D->A_pcols*sizeof(double)); 
  for (nb=0;nb<nbins;nb++){
    L_handle_copy(E_[nb]->lf_At_rsum,E_in_[nb]->lf_At_rsum);
    L_handle_copy(E_[nb]->lf_Zt_rsum,E_in_[nb]->lf_Zt_rsum);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% [finished dcc_copy_A_p]\n");  wkspace_printf();}
}

void dcc_X_nrows_total(struct dcc_ajdk *D)
{
  /* initializes D->A_nrows_total and D->Z_nrows_total for each column using row-mask mr_b, not mr_j */
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb=0,nc=0;
  if (verbose){ printf(" %% [entering dcc_X_nrows_total]\n");}
  D->A_nrows_total=0; D->A_rpop_b_total=0; D->A_rpop_j_total=0;
  D->Z_nrows_total=0; D->Z_rpop_b_total=0; D->Z_rpop_j_total=0;
  for (nb=0;nb<nbins;nb++){
    D->A_nrows_total += E_[nb]->A_nrows; D->A_rpop_b_total += E_[nb]->A_rpop_b; D->A_rpop_j_total += E_[nb]->A_rpop_j;
    D->Z_nrows_total += E_[nb]->Z_nrows; D->Z_rpop_b_total += E_[nb]->Z_rpop_b; D->Z_rpop_j_total += E_[nb]->Z_rpop_j;
    /* for (nb=0;nb<nbins;nb++){ }*/}
  if (verbose){ printf(" %% D->A_nrows_total %d D->Z_nrows_total %d\n",D->A_nrows_total,D->Z_nrows_total);}
  if (verbose){ printf(" %% D->A_rpop_b_total %d D->Z_rpop_b_total %d\n",D->A_rpop_b_total,D->Z_rpop_b_total);}
  if (verbose){ printf(" %% D->A_rpop_j_total %d D->Z_rpop_j_total %d\n",D->A_rpop_j_total,D->Z_rpop_j_total);}
  if (verbose){ printf(" %% [finished dcc_X_nrows_total]\n");}
}

void dcc_load_A_p(struct dcc_ajdk *D)
{
  /* calculates A_p for each column using row-mask mr_b, not mr_j */
  int verbose=0; char tempchar[FNAMESIZE];
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb=0,nc=0,nc_p=0,nc_start=0,nc_final=0;
  int A_pcols=0,A_ncols=0,A_rpop_b_total=0,Z_rpop_b_total=0;
  double *AZ_rsum=NULL,*A_p=NULL,**A_ajdk_p=NULL;
  unsigned char *A_bmc_b=NULL;
  double tmp_sum=0,tmp_num=0;
  if (verbose){ printf(" %% [entering dcc_load_A_p]\n"); wkspace_printf();}
  if (verbose){ printf(" %% \n");}
  D->A_pcols = psize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/POPLENGTH */;
  if (verbose){ printf(" %% nbins %d; A_ncols %d A_pcols %d\n",nbins,D->A_ncols,D->A_pcols);}
  D->AZ_rsum = (double *) wkspace_all0c(D->A_ncols*sizeof(double)); 
  D->A_p = (double *) wkspace_all0c(D->A_pcols*sizeof(double)); 
  D->A_ajdk = (double *) wkspace_all0c(AJDK_TOT*D->A_pcols*sizeof(double)); 
  if (verbose){ printf(" %% D->A_ncols %d D->A_pcols %d \n",D->A_ncols,D->A_pcols); wkspace_printf();}
  A_ncols = D->A_ncols; A_pcols = D->A_pcols; A_rpop_b_total = D->A_rpop_b_total; Z_rpop_b_total = D->Z_rpop_b_total;
  AZ_rsum = D->AZ_rsum; A_p = D->A_p; A_ajdk_p = &(D->A_ajdk); A_bmc_b = D->A_bmc_b; 
  if (verbose){ printf(" %% calculating E_[nb]->lf_At_rsum etc.\n"); wkspace_printf();}
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){
    for (nc=0;nc<D->A_ncols;nc++){ E_[nb]->lf_At_rsum->lf[nc]=0;}
    GLOBAL_pthread_tic(); wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),TYPE_00,SPACING_a,E_[nb]->M_At,&(E_[nb]->lf_At_rsum)); GLOBAL_pthread_toc();
    for (nc=0;nc<D->A_ncols;nc++){ E_[nb]->lf_Zt_rsum->lf[nc]=0;}
    GLOBAL_pthread_tic(); wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),TYPE_00,SPACING_a,E_[nb]->M_Zt,&(E_[nb]->lf_Zt_rsum)); GLOBAL_pthread_toc();
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  if (verbose){ printf(" %% finished calculating E_[nb]->lf_At_rsum etc.\n"); wkspace_printf();}
  if (verbose){
    for (nb=0;nb<nbins;nb++){
      sprintf(D->tmpAnchar," %%%% E_[nb]->lf_At_[%d]_rsum->lf: ",nb); lfprintf(E_[nb]->lf_At_rsum,D->tmpAnchar);
      sprintf(D->tmpZnchar," %%%% E_[nb]->lf_Zt_[%d]_rsum->lf: ",nb); lfprintf(E_[nb]->lf_Zt_rsum,D->tmpZnchar);
      /* for (nb=0;nb<nbins;nb++){ } */}
    wkspace_printf(); /* if (verbose){ } */}
  if (verbose){ printf(" %% calculating AZ_rsum\n");}
  for (nc=0;nc<A_ncols;nc++){ AZ_rsum[nc]=0;}
  for (nb=0;nb<nbins;nb++){ 
    dra_plusequals(AZ_rsum,D->A_ncols,E_[nb]->lf_At_rsum->lf);
    dra_plusequals(AZ_rsum,D->A_ncols,E_[nb]->lf_Zt_rsum->lf);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ raprintf(AZ_rsum,"double",1,A_ncols," %% AZ_rsum: "); wkspace_printf();}
  if (verbose){ printf(" %% calculating A_p\n");}
  for (nc_p=0;nc_p<A_pcols;nc_p++){ A_p[nc_p]=0.5;}
  if (GLOBAL_TEST_sparse){
    for (nc_p=0;nc_p<A_pcols;nc_p++){ 
      nc_start = minimum(A_ncols-1,maximum(0,0 + (nc_p+0)*POPLENGTH - 0));
      nc_final = minimum(A_ncols-1,maximum(0,0 + (nc_p+1)*POPLENGTH - 1));
      if (verbose){ printf(" %% nc_p %d: [%d..%d]:",nc_p,nc_start,nc_final);}
      tmp_sum = 0; tmp_num=0;
      for (nc=nc_start;nc<=nc_final;nc++){
	if (bget__on(A_bmc_b,nc)){ tmp_sum += AZ_rsum[nc]; tmp_num+=1;}
	/* for (nc=nc_start;nc<=nc_final;nc++){ } */}
      if (verbose){ printf(" tmp_sum %f tmp_num %f*(%d+%d);",tmp_sum,tmp_num,A_rpop_b_total,Z_rpop_b_total);}
      if (tmp_num==0){ A_p[nc_p] = 0.5;}
      else /* if (tmp_num>0) */{ A_p[nc_p] = crop((tmp_sum)/(tmp_num*(A_rpop_b_total+Z_rpop_b_total)),0.01,0.99);}
      if (verbose){ printf(" A_p %f\n",A_p[nc_p]);}
      /*  for (nc_p=0;nc_p<A_pcols;nc_p++){ } */}
    /* if (GLOBAL_TEST_sparse){ } */}
  if (verbose){ raprintf(A_p,"double",1,A_pcols," %% A_p: "); wkspace_printf();}
  calc_A_ajdk(A_p,A_pcols,A_ajdk_p); if (verbose){ raprintf(*A_ajdk_p,"double",A_pcols,AJDK_TOT," %% A_ajdk: "); wkspace_printf();}
  sprintf(tempchar,"%s/AZ_rsum.mda",GLOBAL_DIR_NAME); mda_write_r8(tempchar,A_ncols,1,AZ_rsum);
  sprintf(tempchar,"%s/A_p.mda",GLOBAL_DIR_NAME); mda_write_r8(tempchar,A_pcols,1,A_p);
  if (verbose){ printf(" %% [finished dcc_load_A_p]\n");  wkspace_printf();}
}

void dcc_init_A_p(double mrnd,struct dcc_ajdk *D)
{
  /* initializes A_p for each column using row-mask mr_b, not mr_j */
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb=0,nc=0;
  if (verbose){ printf(" %% [entering dcc_init_A_p]\n");}
  D->A_pcols = psize(D->A_ncols); 
  D->A_p = (double *)wkspace_all0c(D->A_pcols*sizeof(double));
  for (nc=0;nc<D->A_pcols;nc++){ D->A_p[nc] = maximum(0.1,minimum(0.9,rand01));} if (GLOBAL_TEST_sparse==0){ printf(" %% turning off D and a\n"); for (nc=0;nc<D->A_pcols;nc++){ D->A_p[nc] = 0.5;}}
  D->A_ajdk = (double *)wkspace_all0c(AJDK_TOT*D->A_pcols*sizeof(double)); 
  calc_A_ajdk(D->A_p,D->A_pcols,&(D->A_ajdk));
  if (verbose){ printf(" %% [finished dcc_init_A_p]\n");}
}

void dcc_copy_QX(struct dcc_ajdk *D,struct dcc_ajdk *D_in)
{
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_single **E_in_ = D_in->E_;
  int nb1=0; struct dcc_single *E=NULL;
  D->out_iteration=0; D->out_trace_length = 6; memcpy(D->out_trace,D_in->out_trace,(D->A_ncols+D->A_nrows_total)*D->out_trace_length*sizeof(double));
  memcpy(D->out_xdrop_a,D_in->out_xdrop_a,(D->A_ncols+D->A_nrows_total)*2*sizeof(int));
  memcpy(D->out_xdrop_b,D_in->out_xdrop_b,(D->A_ncols+D->A_nrows_total)*2*sizeof(int));
  memcpy(D->QC_TAnZtS_nrm,D_in->QC_TAnZtS_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_TAnAtT_nrm,D_in->QC_TAnAtT_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_sra,D_in->QC_sra,D->A_ncols*sizeof(double));
  memcpy(D->QC_lmc_a,D_in->QC_lmc_a,D->A_ncols*sizeof(int));
  memcpy(D->QC_lmc_b,D_in->QC_lmc_b,D->A_ncols*sizeof(int));
  memcpy(D->QC_lmc_j,D_in->QC_lmc_j,D->A_ncols*sizeof(int));
  memcpy(D->QR_sra,D_in->QR_sra,D->A_nrows_total*sizeof(double));
  memcpy(D->QR_imr_a,D_in->QR_imr_a,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_imr_b,D_in->QR_imr_b,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_imr_j,D_in->QR_imr_j,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_lmr_a,D_in->QR_lmr_a,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_lmr_b,D_in->QR_lmr_b,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_lmr_j,D_in->QR_lmr_j,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_lnb,D_in->QR_lnb,D->A_nrows_total*sizeof(int));
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    memcpy(E->QR_imr_a,E_in_[nb1]->QR_imr_a,E->A_nrows*sizeof(int));
    memcpy(E->QR_imr_b,E_in_[nb1]->QR_imr_b,E->A_nrows*sizeof(int));
    memcpy(E->QR_imr_j,E_in_[nb1]->QR_imr_j,E->A_nrows*sizeof(int));
    memcpy(E->QR_TAnZtS_nrm,E_in_[nb1]->QR_TAnZtS_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_TAnAtT_nrm,E_in_[nb1]->QR_TAnAtT_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    /*  for (nb1=0;nb1<nbins;nb1++){ } */}
}

void dcc_init_QX(struct dcc_ajdk *D)
{
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL; int mx_a=0,mx_b=0,mx_j=0,mr_a=0,mr_b=0,mr_j=0;
  D->out_iteration=0; D->out_trace_length = 6; D->out_trace = (double *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*D->out_trace_length*sizeof(double)); if (!D->out_trace){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->out_xdrop_a = (int *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*2*sizeof(int)); if (!D->out_xdrop_a){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->out_xdrop_b = (int *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*2*sizeof(int)); if (!D->out_xdrop_b){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QC_TAnZtS_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_TAnZtS_nrm){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");}
  D->QC_TAnAtT_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_TAnAtT_nrm){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");}
  D->QC_sra=(double *)wkspace_all0c(D->A_ncols*sizeof(double)); if (!D->QC_sra){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QC_lmc_a=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_lmc_a){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QC_lmc_b=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_lmc_b){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QC_lmc_j=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_lmc_j){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QR_sra=(double *)wkspace_all0c(D->A_nrows_total*sizeof(double)); if (!D->QR_sra){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QR_imr_a=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_imr_a){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QR_imr_b=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_imr_b){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QR_imr_j=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_imr_j){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QR_lmr_a=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_lmr_a){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QR_lmr_b=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_lmr_b){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QR_lmr_j=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_lmr_j){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  D->QR_lnb=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_lnb){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
  mx_a=0; mx_b=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    E->QR_imr_a = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_imr_a){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
    E->QR_imr_b = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_imr_b){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
    E->QR_imr_j = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_imr_j){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");} 
    mr_a=0;mr_b=0;mr_j=0;
    while (mr_a<E->A_nrows){
      E->QR_imr_a[mr_a] = mx_a;
      if (bget__on(E->A_bmr_b,mr_a)){
	E->QR_imr_b[mr_a] = mx_b;
	if (bget__on(E->A_bmr_j,mr_a)){
	  E->QR_imr_j[mr_a] = mx_j;
	  mr_j++; mx_j++; /* if (bget__on(E->A_bmr_j,mr_a)){ } */}
	mr_b++; mx_b++; /* if (bget__on(E->A_bmr_b,mr_b)){ } */}
      mr_a++; mx_a++; /* while (mr_a<E->A_nrows){ } */}
    E->QR_TAnZtS_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_TAnZtS_nrm){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");}
    E->QR_TAnAtT_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_TAnAtT_nrm){ printf(" %% Warning! not enough memory in dcc_init_X_QX\n");}
    /*  for (nb1=0;nb1<nbins;nb1++){ } */}
}

void dcc_M_mxset(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb=0;
  if (verbose){ printf(" %% setting M_An->mr_j, etc.\n");}
  D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length);
  D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length);
  D->A_rpop_j_total=0;D->Z_rpop_j_total=0;
  for (nb=0;nb<nbins;nb++){ 
    E_[nb]->A_rpop_j = popcount_uchar_array(E_[nb]->A_bmr_j,E_[nb]->A_mr_length); D->A_rpop_j_total += E_[nb]->A_rpop_j;
    E_[nb]->Z_rpop_j = popcount_uchar_array(E_[nb]->Z_bmr_j,E_[nb]->Z_mr_length); D->Z_rpop_j_total += E_[nb]->Z_rpop_j;
    M_mxset(E_[nb]->M_An,E_[nb]->A_bmr_j,D->A_bmc_j); M_mxset(E_[nb]->M_At,D->A_bmc_j,E_[nb]->A_bmr_j);
    M_mxset(E_[nb]->M_Zn,E_[nb]->Z_bmr_j,D->A_bmc_j); M_mxset(E_[nb]->M_Zt,D->A_bmc_j,E_[nb]->Z_bmr_j);
    M_mxset(E_[nb]->M_Tn,E_[nb]->A_bmr_j,D->T_bmc_j); M_mxset(E_[nb]->M_Tt,D->T_bmc_j,E_[nb]->A_bmr_j);
    M_mxset(E_[nb]->M_Sn,E_[nb]->Z_bmr_j,D->T_bmc_j); M_mxset(E_[nb]->M_St,D->T_bmc_j,E_[nb]->Z_bmr_j);
    /* for (nb=0;nb<nbins;nb++){ } */}
}

void dcc_copy(struct dcc_ajdk *D,struct dcc_ajdk *D_in)
{
  int verbose=0;
  struct dcc_single **E_ = D->E_; struct dcc_single **E_in_ = D_in->E_;
  struct dcc_double **F_ = D->F_; struct dcc_double **F_in_ = D_in->F_;
  int nbins = GLOBAL_NBINS, nbx=0,nb1=0,nb2=0;
  if (verbose){ printf(" %% [entering dcc_copy] \n");}
  GLOBAL_tic(0);
  if (verbose){ printf(" %% calling dcc_ajdk_copy: \n");}
  dcc_ajdk_copy(D,D_in);
  if (verbose){ printf(" %% finished dcc_ajdk_copy: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ 
    if (verbose){ printf(" %% calling dcc_single_copy_M_An for nb1 %d: \n",nb1);} 
    dcc_single_copy_M_An(E_[nb1],E_in_[nb1]);
    if (verbose){ printf(" %% finished dcc_single_copy_M_An for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% calling dcc_X_nrows_total: \n");}
  dcc_X_nrows_total(D);
  if (verbose){ printf(" %% finished dcc_X_nrows_total: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ 
    if (verbose){ printf(" %% calling dcc_single_copy_lf for nb1 %d: \n",nb1);}
    dcc_single_copy_lf(E_[nb1],E_in_[nb1]);
    if (verbose){ printf(" %% finished dcc_single_copy_lf for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; 
      if (verbose){ printf(" %% calling dcc_double_copy_lf for nb1 %d nb2 %d nbx %d: \n",nb1,nb2,nbx);}
      dcc_double_copy_lf(F_[nbx],F_in_[nbx]); 
      if (verbose){ printf(" %% finished dcc_double_copy_lf for nb1 %d nb2 %d nbx %d: ",nb1,nb2,nbx); wkspace_printf();}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_toc(0,verbose," %% generating data matrices: ");
  if (verbose){ printf(" %% calling dcc_copy_A_p: \n");}
  dcc_copy_A_p(D,D_in);
  if (verbose){ printf(" %% finished dcc_copy_A_p: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcc_copy_QX\n");}
  dcc_copy_QX(D,D_in);
  if (verbose){ printf(" %% finished dcc_copy_QX: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcc_M_mxset\n");}
  dcc_M_mxset(D);
  if (verbose){ printf(" %% finished dcc_M_mxset: "); wkspace_printf();}
  if (verbose){ printf(" %% [finished dcc_copy] "); wkspace_printf();}
}

void dcc_load(struct dcc_ajdk **D_p,struct dcc_single ***E_p,struct dcc_double ***F_p)
{
  int verbose=GLOBAL_verbose;
  int nbins = GLOBAL_NBINS, nbx=0,nb1=0,nb2=0;
  struct dcc_ajdk *D=NULL;
  struct dcc_single *E=NULL,*E_nb1=NULL,*E_nb2=NULL;
  struct dcc_double *F=NULL;
  if (verbose){ printf(" %% [entering dcc_load] \n");}
  if (*D_p==NULL){
    (*D_p) = (struct dcc_ajdk *) wkspace_all0c(sizeof(struct dcc_ajdk)); D = (*D_p); D->bitj = BITJ; D->nbins = nbins; D->E_ = NULL; D->F_ = NULL;
    /* if (*D_p==NULL){ } */}
  if (*E_p==NULL){
    (*E_p) = (struct dcc_single **) wkspace_all0c(nbins*sizeof(struct dcc_single *));
    for (nb1=0;nb1<nbins;nb1++){ (*E_p)[nb1] = (struct dcc_single *) wkspace_all0c(sizeof(struct dcc_single)*1); E = (*E_p)[nb1]; E->D = (*D_p); E->nb=nb1;}
    (*D_p)->E_ = (*E_p);
    /* if (*E_p==NULL){ } */}
  if (*F_p==NULL){
    (*F_p) = (struct dcc_double **) wkspace_all0c(nbins*nbins*sizeof(struct dcc_double *));
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins;
	(*F_p)[nbx] = (struct dcc_double *) wkspace_all0c(sizeof(struct dcc_double)*1); F = (*F_p)[nbx]; F->D = (*D_p); F->nb1=nb1; F->nb2=nb2; F->E_nb1 = (*E_p)[nb1]; F->E_nb2 = (*E_p)[nb2];
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    (*D_p)->F_ = (*F_p);
    /* if (*F_p==NULL){ } */}
  D = (*D_p); 
  GLOBAL_tic(0);
  if (verbose){ printf(" %% calling dcc_ajdk_load: \n");}
  dcc_ajdk_load(D);
  if (verbose){ printf(" %% finished dcc_ajdk_load: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling dcc_single_load_M_An for nb1 %d: \n",nb1);} 
    dcc_single_load_M_An(E);
    if (verbose){ printf(" %% finished dcc_single_load_M_An for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% calling dcc_X_nrows_total: \n");}
  dcc_X_nrows_total(D);
  if (verbose){ printf(" %% finished dcc_X_nrows_total: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling dcc_single_init_lf for nb1 %d: \n",nb1);}
    dcc_single_init_lf(E);
    if (verbose){ printf(" %% finished dcc_single_init_lf for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; E_nb1 = (*E_p)[nb1]; E_nb2 = (*E_p)[nb2]; F = (*F_p)[nbx];
      if (verbose){ printf(" %% calling dcc_double_init_lf for nb1 %d nb2 %d nbx %d: \n",nb1,nb2,nbx);}
      dcc_double_init_lf(F); 
      if (verbose){ printf(" %% finished dcc_double_init_lf for nb1 %d nb2 %d nbx %d: ",nb1,nb2,nbx); wkspace_printf();}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_toc(0,verbose," %% generating data matrices: ");
  if (verbose){ printf(" %% calling dcc_load_A_p: \n");}
  dcc_load_A_p(D);
  if (verbose){ printf(" %% finished dcc_load_A_p: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcc_init_QX\n");}
  dcc_init_QX(D);
  if (verbose){ printf(" %% finished dcc_init_QX: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcc_M_mxset\n");}
  dcc_M_mxset(D);
  if (verbose){ printf(" %% finished dcc_M_mxset: "); wkspace_printf();}
  if (verbose){ printf(" %% [finished dcc_load] "); wkspace_printf();}
}

void dcc_init(char *error_vs_speed,double mrnd,struct dcc_ajdk **D_p,struct dcc_single ***E_p,struct dcc_double ***F_p)
{
  int verbose=GLOBAL_verbose;
  int nbins = GLOBAL_NBINS, nbx=0,nb1=0,nb2=0;
  struct dcc_ajdk *D=NULL;
  struct dcc_single *E=NULL,*E_nb1=NULL,*E_nb2=NULL;
  struct dcc_double *F=NULL;
  if (verbose){ printf(" %% [entering dcc_init] \n");}
  if (*D_p==NULL){
    (*D_p) = (struct dcc_ajdk *) wkspace_all0c(sizeof(struct dcc_ajdk)); D = (*D_p); D->bitj = BITJ; D->nbins = nbins; D->E_ = NULL; D->F_ = NULL;
    /* if (*D_p==NULL){ } */}
  if (*E_p==NULL){
    (*E_p) = (struct dcc_single **) wkspace_all0c(nbins*sizeof(struct dcc_single *));
    for (nb1=0;nb1<nbins;nb1++){ (*E_p)[nb1] = (struct dcc_single *) wkspace_all0c(sizeof(struct dcc_single)*1); E = (*E_p)[nb1]; E->D = (*D_p); E->nb=nb1;}
    (*D_p)->E_ = (*E_p);
    /* if (*E_p==NULL){ } */}
  if (*F_p==NULL){
    (*F_p) = (struct dcc_double **) wkspace_all0c(nbins*nbins*sizeof(struct dcc_double *));
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins;
	(*F_p)[nbx] = (struct dcc_double *) wkspace_all0c(sizeof(struct dcc_double)*1); F = (*F_p)[nbx]; F->D = (*D_p); F->nb1=nb1; F->nb2=nb2; F->E_nb1 = (*E_p)[nb1]; F->E_nb2 = (*E_p)[nb2];
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    (*D_p)->F_ = (*F_p);
    /* if (*F_p==NULL){ } */}
  D = (*D_p);
  GLOBAL_tic(0);
  if (verbose){ printf(" %% calling dcc_ajdk_init: \n");}
  dcc_ajdk_init(mrnd,D);
  if (verbose){ printf(" %% finished dcc_ajdk_init: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling dcc_single_init_M_An for nb1 %d: \n",nb1);}
    dcc_single_init_M_An(error_vs_speed,mrnd,E);
    if (verbose){ printf(" %% finished dcc_single_init_M_An for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% calling dcc_X_nrows_total: \n");}
  dcc_X_nrows_total(D);
  if (verbose){ printf(" %% finished dcc_X_nrows_total: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling dcc_single_init_lf for nb1 %d: \n",nb1);}
    dcc_single_init_lf(E);
    if (verbose){ printf(" %% finished dcc_single_init_lf for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; E_nb1 = (*E_p)[nb1]; E_nb2 = (*E_p)[nb2]; F = (*F_p)[nbx];
      if (verbose){ printf(" %% calling dcc_double_init_lf for nb1 %d nb2 %d nbx %d: \n",nb1,nb2,nbx);}
      dcc_double_init_lf(F);
      if (verbose){ printf(" %% finished dcc_double_init_lf for nb1 %d nb2 %d nbx %d: ",nb1,nb2,nbx); wkspace_printf();}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_toc(0,verbose," %% generating data matrices: ");
  if (verbose){ printf(" %% calling dcc_init_A_p: \n");}
  dcc_init_A_p(mrnd,D);
  if (verbose){ printf(" %% finished dcc_init_A_p: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcc_init_QX\n");}
  dcc_init_QX(D);
  if (verbose){ printf(" %% finished dcc_init_QX: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcc_M_mxset\n");}
  dcc_M_mxset(D);
  if (verbose){ printf(" %% finished dcc_M_mxset: "); wkspace_printf();}
  if (verbose){ printf(" %% [finished dcc_init] "); wkspace_printf();}
}

void dcc_init_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  dcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); dcc_init_QX(D);
  M_handle_printf(E_[0]->M_An,verbose," %% E_[0]->M_An[0]: "); M_handle_printf(E_[0]->M_At,verbose," %% E_[0]->M_At[0]: ");
  M_handle_printf(E_[1]->M_An,verbose," %% E_[1]->M_An[1]: "); M_handle_printf(E_[1]->M_At,verbose," %% E_[1]->M_At[1]: ");
  M_handle_printf(E_[0]->M_Zn,verbose," %% E_[0]->M_Zn[0]: "); M_handle_printf(E_[0]->M_Zt,verbose," %% E_[0]->M_Zt[0]: ");
  M_handle_printf(E_[1]->M_Zn,verbose," %% E_[1]->M_Zn[1]: "); M_handle_printf(E_[1]->M_Zt,verbose," %% E_[1]->M_Zt[1]: ");
  M_handle_printf(E_[0]->M_Tn,verbose," %% E_[0]->M_Tn[0]: "); M_handle_printf(E_[0]->M_Tt,verbose," %% E_[0]->M_Tt[0]: ");
  M_handle_printf(E_[1]->M_Tn,verbose," %% E_[1]->M_Tn[1]: "); M_handle_printf(E_[1]->M_Tt,verbose," %% E_[1]->M_Tt[1]: ");
  M_handle_printf(E_[0]->M_Sn,verbose," %% E_[0]->M_Sn[0]: "); M_handle_printf(E_[0]->M_St,verbose," %% E_[0]->M_St[0]: ");
  M_handle_printf(E_[1]->M_Sn,verbose," %% E_[1]->M_Sn[1]: "); M_handle_printf(E_[1]->M_St,verbose," %% E_[1]->M_St[1]: ");
  raprintf(D->A_p,"double",1,psize(D->A_ncols)," %% D->A_p: "); raprintf(D->A_ajdk,"double",psize(D->A_ncols),AJDK_TOT," %% D->A_ajdk: ");
  wkspace_printf();
}
