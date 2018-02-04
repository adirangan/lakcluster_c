void bcc_ajdk_copy(struct bcc_ajdk *D,struct bcc_ajdk *D_in)
{
  /* copy from D_in; assuming D has been preallocated. */
  int verbose=0;
  int nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0;
  if (verbose){ printf(" %% [entering bcc_ajdk_copy]\n");}
  D->A_ncols=D_in->A_ncols;D->Y_ncols=D_in->Y_ncols;D->T_ncols=D_in->T_ncols;
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
  /* loading Y_bmc_b */
  D->Y_cbother = (D->Y_ncols>0);
  D->Y_ncols_extend = (D->bitj - (D->Y_ncols % D->bitj)) % D->bitj; 
  D->Y_mc_length = bsize(D->Y_ncols)/* rup(D->Y_ncols+D->Y_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% Y_ncols %d, Y_ncols_extend %d Y_mc_length %d\n",D->Y_ncols,D->Y_ncols_extend,D->Y_mc_length);}
  memcpy(D->Y_bmc_b,D_in->Y_bmc_b,D->Y_mc_length);
  memcpy(D->Y_bmc_j,D_in->Y_bmc_j,D->Y_mc_length);
  D->Y_cpop_b= popcount_uchar_array(D->Y_bmc_b,D->Y_mc_length); D->Y_cpop_j = popcount_uchar_array(D->Y_bmc_j,D->Y_mc_length); if (verbose>1){ printf(" %%%% D->Y_cpop_b %d D->Y_cpop_j %d\n",D->Y_cpop_b,D->Y_cpop_j);}
  memcpy(D->Y_umc_b,D_in->Y_umc_b,D->Y_ncols*sizeof(unsigned char));
  memcpy(D->Y_umc_j,D_in->Y_umc_j,D->Y_ncols*sizeof(unsigned char));
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
  if (verbose){ printf(" %% [finished bcc_ajdk_copy] D->A_ncols=%d D->Y_ncols=%d D->T_ncols=%d \n",D->A_ncols,D->Y_ncols,D->T_ncols);}
}

void bcc_ajdk_load(struct bcc_ajdk *D)
{
  /* load from file */
  int verbose=0;
  int nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0;
  if (verbose){ printf(" %% [entering bcc_ajdk_load]\n");}
  D->A_ncols=GLOBAL_A_n_cols;D->Y_ncols=GLOBAL_Y_n_cols;D->T_ncols=GLOBAL_T_n_cols;
  /* loading A_bmc_b */
  D->A_cbother = (D->A_ncols>0);
  D->A_ncols_extend = (D->bitj - (D->A_ncols % D->bitj)) % D->bitj; 
  D->A_mc_length = bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% A_ncols %d, A_ncols_extend %d A_mc_length %d\n",D->A_ncols,D->A_ncols_extend,D->A_mc_length);}
  D->A_bmc_b = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} 
  D->A_bmc_j = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} 
  D->A_bmc_j_rmv = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rmv){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} 
  D->A_bmc_j_rtn = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rtn){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} 
  if (GLOBAL_A_n_cind==NULL || !strcmp(GLOBAL_A_n_cind,"\0")){ for (nc=0;nc<D->A_ncols;nc++){ bset__on(D->A_bmc_b,nc);}}
  else{ binary_read(GLOBAL_A_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->A_bmc_b)); if (nrows_tmp!=D->A_ncols){ printf(" %% Warning! A_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_A_n_name,GLOBAL_A_n_cind,nrows_tmp,D->A_ncols);}}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_b: "); if (verbose>2){ bprintf(D->A_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  for (nc=0;nc<D->A_mc_length;nc++){ D->A_bmc_j[nc] = D->A_bmc_b[nc];}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); if (verbose>2){ bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  D->A_cpop_b= popcount_uchar_array(D->A_bmc_b,D->A_mc_length); D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->A_cpop_b %d D->A_cpop_j %d\n",D->A_cpop_b,D->A_cpop_j);}
  D->A_umc_b = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_b[nc] = bget__on(D->A_bmc_b,nc);}
  D->A_umc_j = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);}
  D->A_umc_j_rmv = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rmv){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rmv[nc] = bget__on(D->A_bmc_j_rmv,nc);}
  D->A_umc_j_rtn = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rtn){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = bget__on(D->A_bmc_j_rtn,nc);}
  /* loading Y_bmc_b */
  D->Y_cbother = (D->Y_ncols>0);
  D->Y_ncols_extend = (D->bitj - (D->Y_ncols % D->bitj)) % D->bitj; 
  D->Y_mc_length = bsize(D->Y_ncols)/* rup(D->Y_ncols+D->Y_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% Y_ncols %d, Y_ncols_extend %d Y_mc_length %d\n",D->Y_ncols,D->Y_ncols_extend,D->Y_mc_length);} 
  D->Y_bmc_b = wkspace_all0c(D->Y_mc_length); if (!D->Y_bmc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} 
  D->Y_bmc_j = wkspace_all0c(D->Y_mc_length); if (!D->Y_bmc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} 
  if (GLOBAL_Y_n_cind==NULL || !strcmp(GLOBAL_Y_n_cind,"\0")){ for (nc=0;nc<D->Y_ncols;nc++){ bset__on(D->Y_bmc_b,nc);}}
  else{ binary_read(GLOBAL_Y_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->Y_bmc_b)); if (nrows_tmp!=D->Y_ncols){ printf(" %% Warning! Y_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_Y_n_name,GLOBAL_Y_n_cind,nrows_tmp,D->Y_ncols);}}
  sprintf(D->tmpYnchar," %%%% D->Y_bmc_b: "); if (verbose>2){ bprintf(D->Y_bmc_b,D->bitj,1,D->Y_ncols,D->tmpYnchar);}
  for (nc=0;nc<D->Y_mc_length;nc++){ D->Y_bmc_j[nc] = D->Y_bmc_b[nc];}
  sprintf(D->tmpYnchar," %%%% D->Y_bmc_j: "); if (verbose>2){ bprintf(D->Y_bmc_j,D->bitj,1,D->Y_ncols,D->tmpYnchar);}
  D->Y_cpop_b= popcount_uchar_array(D->Y_bmc_b,D->Y_mc_length); D->Y_cpop_j = popcount_uchar_array(D->Y_bmc_j,D->Y_mc_length); if (verbose>1){ printf(" %%%% D->Y_cpop_b %d D->Y_cpop_j %d\n",D->Y_cpop_b,D->Y_cpop_j);}
  D->Y_umc_b = (unsigned char *) wkspace_all0c(D->Y_ncols*sizeof(unsigned char)); if (!D->Y_umc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} for (nc=0;nc<D->Y_ncols;nc++){ D->Y_umc_b[nc] = bget__on(D->Y_bmc_b,nc);}
  D->Y_umc_j = (unsigned char *) wkspace_all0c(D->Y_ncols*sizeof(unsigned char)); if (!D->Y_umc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} for (nc=0;nc<D->Y_ncols;nc++){ D->Y_umc_j[nc] = bget__on(D->Y_bmc_j,nc);}
  /* loading T_bmc_b */
  D->T_ncols_extend = (D->bitj - (D->T_ncols % D->bitj)) % D->bitj; 
  D->T_mc_length = bsize(D->T_ncols)/* rup(D->T_ncols+D->T_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% T_ncols %d, T_ncols_extend %d T_mc_length %d\n",D->T_ncols,D->T_ncols_extend,D->T_mc_length);}
  D->T_bmc_b = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} 
  D->T_bmc_j = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} 
  if (GLOBAL_T_n_cind==NULL || !strcmp(GLOBAL_T_n_cind,"\0")){ for (nc=0;nc<D->T_ncols;nc++){ bset__on(D->T_bmc_b,nc);}}
  else{ binary_read(GLOBAL_T_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->T_bmc_b)); if (nrows_tmp!=D->T_ncols){ printf(" %% Warning! T_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_T_n_name,GLOBAL_T_n_cind,nrows_tmp,D->T_ncols);}}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_b: "); if (verbose>2){ bprintf(D->T_bmc_b,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  for (nc=0;nc<D->T_mc_length;nc++){ D->T_bmc_j[nc] = D->T_bmc_b[nc];}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_j: "); if (verbose>2){ bprintf(D->T_bmc_j,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  D->T_cpop_b= popcount_uchar_array(D->T_bmc_b,D->T_mc_length); D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length); if (verbose>1){ printf(" %%%% D->T_cpop_b %d D->T_cpop_j %d\n",D->T_cpop_b,D->T_cpop_j);}
  D->T_umc_b = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_b[nc] = bget__on(D->T_bmc_b,nc);}
  D->T_umc_j = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_load\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_j[nc] = bget__on(D->T_bmc_j,nc);}
  D->QC_AtTYnWtSZn_nrm = NULL; D->QC_AtTAnZtSZn_nrm = NULL; D->QC_AtTYnYtTAn_nrm = NULL; D->QC_AtTAnAtTAn_nrm = NULL;
  D->QC_sra=NULL; D->QC_lmc_a=NULL; D->QC_lmc_b=NULL; D->QC_lmc_j=NULL;
  D->QR_sra=NULL; D->QR_lmr_a=NULL; D->QR_lmr_b=NULL; D->QR_lmr_j=NULL; D->QR_lnb=NULL; D->QR_imr_a=NULL; D->QR_imr_b=NULL; D->QR_imr_j=NULL;
  D->Irem=0;D->Ireq=0; D->out_iteration=0; D->out_xdrop_ij=0; D->out_trace=NULL; D->out_xdrop_a=NULL; D->out_xdrop_b=NULL; D->out_trace_length = 0;
  if (verbose){ printf(" %% [finished bcc_ajdk_load] D->A_ncols=%d D->Y_ncols=%d D->T_ncols=%d \n",D->A_ncols,D->Y_ncols,D->T_ncols);}
}

void bcc_ajdk_init(double mrnd,struct bcc_ajdk *D)
{
  /* initialize based on test parameters */
  int verbose=0;
  int nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0;
  unsigned char *bb=NULL,*bj=NULL;
  double dtmp=0; int nx=0;
  if (verbose){ printf(" %% [entering bcc_ajdk_init]\n");}
  D->A_ncols=GLOBAL_TEST_A_n_cols;D->Y_ncols=GLOBAL_TEST_Y_n_cols;D->T_ncols=GLOBAL_TEST_T_n_cols;
  /* initializing A_bmc_b */
  D->A_cbother = (D->A_ncols>0);
  D->A_ncols_extend = (D->bitj - (D->A_ncols % D->bitj)) % D->bitj; 
  D->A_mc_length = bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% A_ncols %d, A_ncols_extend %d A_mc_length %d\n",D->A_ncols,D->A_ncols_extend,D->A_mc_length);}
  D->A_bmc_b = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} 
  D->A_bmc_j = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} 
  bb = D->A_bmc_b; bj = D->A_bmc_j; for (nx=0;nx<D->A_ncols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_b: "); if (verbose>2){ bprintf(D->A_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); if (verbose>2){ bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  D->A_bmc_j_rmv = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rmv){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} 
  D->A_bmc_j_rtn = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rtn){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} 
  D->A_cpop_b= popcount_uchar_array(D->A_bmc_b,D->A_mc_length); D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->A_cpop_b %d D->A_cpop_j %d\n",D->A_cpop_b,D->A_cpop_j);}
  D->A_umc_b = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_b[nc] = bget__on(D->A_bmc_b,nc);}
  D->A_umc_j = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);}
  D->A_umc_j_rmv = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rmv){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rmv[nc] = bget__on(D->A_bmc_j_rmv,nc);}
  D->A_umc_j_rtn = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rtn){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = bget__on(D->A_bmc_j_rtn,nc);}
  /* initializing Y_bmc_b */
  D->Y_cbother = (D->Y_ncols>0);
  D->Y_ncols_extend = (D->bitj - (D->Y_ncols % D->bitj)) % D->bitj; 
  D->Y_mc_length = bsize(D->Y_ncols)/* rup(D->Y_ncols+D->Y_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% Y_ncols %d, Y_ncols_extend %d Y_mc_length %d\n",D->Y_ncols,D->Y_ncols_extend,D->Y_mc_length);}
  D->Y_bmc_b = wkspace_all0c(D->Y_mc_length); if (!D->Y_bmc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} 
  D->Y_bmc_j = wkspace_all0c(D->Y_mc_length); if (!D->Y_bmc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} 
  bb = D->Y_bmc_b; bj = D->Y_bmc_j; for (nx=0;nx<D->Y_ncols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpYnchar," %%%% D->Y_bmc_b: "); if (verbose>2){ bprintf(D->Y_bmc_b,D->bitj,1,D->Y_ncols,D->tmpYnchar);}
  sprintf(D->tmpYnchar," %%%% D->Y_bmc_j: "); if (verbose>2){ bprintf(D->Y_bmc_j,D->bitj,1,D->Y_ncols,D->tmpYnchar);}
  D->Y_cpop_b= popcount_uchar_array(D->Y_bmc_b,D->Y_mc_length); D->Y_cpop_j = popcount_uchar_array(D->Y_bmc_j,D->Y_mc_length); if (verbose>1){ printf(" %%%% D->Y_cpop_b %d D->Y_cpop_j %d\n",D->Y_cpop_b,D->Y_cpop_j);}
  D->Y_umc_b = (unsigned char *) wkspace_all0c(D->Y_ncols*sizeof(unsigned char)); if (!D->Y_umc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} for (nc=0;nc<D->Y_ncols;nc++){ D->Y_umc_b[nc] = bget__on(D->Y_bmc_b,nc);}
  D->Y_umc_j = (unsigned char *) wkspace_all0c(D->Y_ncols*sizeof(unsigned char)); if (!D->Y_umc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} for (nc=0;nc<D->Y_ncols;nc++){ D->Y_umc_j[nc] = bget__on(D->Y_bmc_j,nc);}
  /* initializing T_bmc_b */
  D->T_ncols_extend = (D->bitj - (D->T_ncols % D->bitj)) % D->bitj; 
  D->T_mc_length = bsize(D->T_ncols)/* rup(D->T_ncols+D->T_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% T_ncols %d, T_ncols_extend %d T_mc_length %d\n",D->T_ncols,D->T_ncols_extend,D->T_mc_length);}
  D->T_bmc_b = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} 
  D->T_bmc_j = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} 
  bb = D->T_bmc_b; bj = D->T_bmc_j; for (nx=0;nx<D->T_ncols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}} bset__on(bb,0); bset__on(bj,0); bset__on(bb,D->T_ncols-1), bset__on(bj,D->T_ncols-1);
  sprintf(D->tmpTnchar," %%%% D->T_bmc_b: "); if (verbose>2){ bprintf(D->T_bmc_b,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_j: "); if (verbose>2){ bprintf(D->T_bmc_j,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  D->T_cpop_b= popcount_uchar_array(D->T_bmc_b,D->T_mc_length); D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length); if (verbose>1){ printf(" %%%% D->T_cpop_b %d D->T_cpop_j %d\n",D->T_cpop_b,D->T_cpop_j);}
  D->T_umc_b = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_b){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_b[nc] = bget__on(D->T_bmc_b,nc);}
  D->T_umc_j = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_j){ printf(" %% Warning! not enough memory in bcc_ajdk_init\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_j[nc] = bget__on(D->T_bmc_j,nc);}
  D->QC_AtTYnWtSZn_nrm = NULL; D->QC_AtTAnZtSZn_nrm = NULL; D->QC_AtTYnYtTAn_nrm = NULL; D->QC_AtTAnAtTAn_nrm = NULL;
  D->QC_sra=NULL; D->QC_lmc_a=NULL; D->QC_lmc_b=NULL; D->QC_lmc_j=NULL;
  D->QR_sra=NULL; D->QR_lmr_a=NULL; D->QR_lmr_b=NULL; D->QR_lmr_j=NULL; D->QR_lnb=NULL; D->QR_imr_a=NULL; D->QR_imr_b=NULL; D->QR_imr_j=NULL;
  D->Irem=0;D->Ireq=0; D->out_iteration=0; D->out_xdrop_ij=0; D->out_trace=NULL; D->out_xdrop_a=NULL; D->out_xdrop_b=NULL; D->out_trace_length = 0;
  if (verbose){ printf(" %% [finished bcc_ajdk_init] D->A_ncols=%d D->Y_ncols=%d D->T_ncols=%d \n",D->A_ncols,D->Y_ncols,D->T_ncols);}
}

void bcc_single_copy_M_An(struct bcc_single *E,struct bcc_single *E_in)
{
  /* copy from file */
  int verbose=0;
  struct bcc_ajdk *D=E->D,*D_in=E_in->D;
  int nb_given = E->nb;
  int nr=0,nb=0,nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0,nrows_tmp_extend,brows_tmp;
  unsigned char *bXra_tmp=NULL;
  if (verbose){ printf(" %% [entering bcc_single_copy_M_An]\n");}
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
  /* copying Y */
  D->Y_ncols = D_in->Y_ncols;
  M_handle_copy(E->M_Yn,E_in->M_Yn);
  M_handle_copy(E->M_Yt,E_in->M_Yt);
  /* copying W */
  M_handle_copy(E->M_Wn,E_in->M_Wn);
  M_handle_copy(E->M_Wt,E_in->M_Wt);
  /* copying T */
  D->T_ncols = D_in->T_ncols;
  M_handle_copy(E->M_Tn,E_in->M_Tn);
  M_handle_copy(E->M_Tt,E_in->M_Tt);
  /* copying S */
  M_handle_copy(E->M_Sn,E_in->M_Sn);
  M_handle_copy(E->M_St,E_in->M_St);
  /* copying a,j,y,t,k,h,u,r,z,v */
  M_handle_copy(E->M_an,E_in->M_an);
  M_handle_copy(E->M_jn,E_in->M_jn);
  M_handle_copy(E->M_kn,E_in->M_kn);
  M_handle_copy(E->M_hn,E_in->M_hn);
  M_handle_copy(E->M_at,E_in->M_at);
  M_handle_copy(E->M_jt,E_in->M_jt);
  M_handle_copy(E->M_kt,E_in->M_kt);
  M_handle_copy(E->M_ht,E_in->M_ht);
  M_handle_copy(E->M_yn,E_in->M_yn);
  M_handle_copy(E->M_un,E_in->M_un);
  M_handle_copy(E->M_yt,E_in->M_yt);
  M_handle_copy(E->M_ut,E_in->M_ut);
  M_handle_copy(E->M_zn,E_in->M_zn);
  M_handle_copy(E->M_vn,E_in->M_vn);
  M_handle_copy(E->M_zt,E_in->M_zt);
  M_handle_copy(E->M_vt,E_in->M_vt);
  M_handle_copy(E->M_tt,E_in->M_tt);
  M_handle_copy(E->M_rt,E_in->M_rt);
  if (verbose){ printf(" %% [finished bcc_single_copy_M_An]\n");}
}

void bcc_single_load_M_An(struct bcc_single *E)
{
  /* load from file */
  int verbose=0;
  struct bcc_ajdk *D=E->D;
  int nb_given = E->nb;
  int nr=0,nb=0,nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0,nrows_tmp_extend,brows_tmp;
  unsigned char *bXra_tmp=NULL;
  if (verbose){ printf(" %% [entering bcc_single_load_M_An]\n");}
  /* loading A_bmr_b etc */
  nb = nb_given;
  /* loading A */
  binary_read_getsize(GLOBAL_A_t_name_[nb],&(D->bitj),&(D->A_ncols),&(E->A_nrows));
  E->A_rbother = (E->A_nrows>0);
  E->A_nrows_extend = (D->bitj - (E->A_nrows % D->bitj)) % D->bitj;
  E->A_mr_length = bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */;
  E->A_bmr_b = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_b){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");}
  E->A_bmr_j = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");}
  E->A_bmr_j_rmv = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rmv){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");}
  E->A_bmr_j_rtn = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rtn){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");}
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
  E->A_umr_b = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_b){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_b[nr] = bget__on(E->A_bmr_b,nr);}
  E->A_umr_j = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);}
  E->A_umr_j_rmv = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rmv){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rmv[nr] = bget__on(E->A_bmr_j_rmv,nr);}
  E->A_umr_j_rtn = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rtn){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = bget__on(E->A_bmr_j_rtn,nr);}
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
  E->Z_bmr_b = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_b){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");}
  E->Z_bmr_j = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_j){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");}
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
  E->Z_umr_b = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_b){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_b[nr] = bget__on(E->Z_bmr_b,nr);}
  E->Z_umr_j = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_j){ printf(" %% Warning! not enough memory in bcc_single_load_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_j[nr] = bget__on(E->Z_bmr_j,nr);}
  if (verbose>1){ printf(" %% generating M_Zn_[%d]\n",nb);}
  if (GLOBAL_Z_t_name_[nb]==NULL || !strcmp(GLOBAL_Z_t_name_[nb],"\0")){ E->M_Zn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,NULL,NULL,E->Z_bmr_b,D->A_bmc_b);}
  else{ E->M_Zn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,GLOBAL_Z_t_name_[nb],NULL,E->Z_bmr_b,D->A_bmc_b);} M_mxget(E->M_Zn);
  if (GLOBAL_Z_n_name_[nb]==NULL || !strcmp(GLOBAL_Z_n_name_[nb],"\0")){ E->M_Zt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,NULL,NULL,D->A_bmc_b,E->Z_bmr_b);}
  else{ E->M_Zt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,GLOBAL_Z_n_name_[nb],NULL,D->A_bmc_b,E->Z_bmr_b);} M_mxget(E->M_Zt);
  /* loading Y */
  if (GLOBAL_Y_t_name_[nb]==NULL || !strcmp(GLOBAL_Y_t_name_[nb],"\0")){ bitj_tmp=D->bitj;D->Y_ncols=0;nrows_tmp=E->A_nrows;}
  else{ binary_read_getsize(GLOBAL_Y_t_name_[nb],&bitj_tmp,&(D->Y_ncols),&nrows_tmp);}
  if (bitj_tmp!=D->bitj || nrows_tmp!=E->A_nrows){ printf(" %% Warning! Y_: %s; improper bitj %d vs %d, nrows %d vs %d\n",GLOBAL_Y_t_name_[nb],bitj_tmp,D->bitj,nrows_tmp,E->A_nrows);}
  if (GLOBAL_Y_t_name_[nb]==NULL || !strcmp(GLOBAL_Y_t_name_[nb],"\0")){ E->M_Yn = M_handle_v_make(D->bitj,E->A_nrows,D->Y_ncols,NULL,NULL,E->A_bmr_b,D->Y_bmc_b);}
  else{ E->M_Yn = M_handle_v_make(D->bitj,E->A_nrows,D->Y_ncols,GLOBAL_Y_t_name_[nb],NULL,E->A_bmr_b,D->Y_bmc_b);} M_mxget(E->M_Yn);
  if (GLOBAL_Y_n_name_[nb]==NULL || !strcmp(GLOBAL_Y_n_name_[nb],"\0")){ E->M_Yt = M_handle_v_make(D->bitj,D->Y_ncols,E->A_nrows,NULL,NULL,D->Y_bmc_b,E->A_bmr_b);}
  else{ E->M_Yt = M_handle_v_make(D->bitj,D->Y_ncols,E->A_nrows,GLOBAL_Y_n_name_[nb],NULL,D->Y_bmc_b,E->A_bmr_b);} M_mxget(E->M_Yt);
  /* loading W */
  if (GLOBAL_W_t_name_[nb]==NULL || !strcmp(GLOBAL_W_t_name_[nb],"\0")){ bitj_tmp=D->bitj;ncols_tmp=D->Y_ncols;nrows_tmp=E->Z_nrows;}
  else{ binary_read_getsize(GLOBAL_W_t_name_[nb],&bitj_tmp,&ncols_tmp,&nrows_tmp);}
  if (bitj_tmp!=D->bitj || ncols_tmp!=D->Y_ncols || nrows_tmp!=E->Z_nrows){ printf(" %% Warning! W_: %s; improper bitj %d vs %d, ncols %d vs %d, nrows %d vs %d\n",GLOBAL_W_t_name_[nb],bitj_tmp,D->bitj,ncols_tmp,D->Y_ncols,nrows_tmp,E->Z_nrows);}
  if (GLOBAL_W_t_name_[nb]==NULL || !strcmp(GLOBAL_W_t_name_[nb],"\0")){ E->M_Wn = M_handle_v_make(D->bitj,E->Z_nrows,D->Y_ncols,NULL,NULL,E->Z_bmr_b,D->Y_bmc_b);}
  else{ E->M_Wn = M_handle_v_make(D->bitj,E->Z_nrows,D->Y_ncols,GLOBAL_W_t_name_[nb],NULL,E->Z_bmr_b,D->Y_bmc_b);} M_mxget(E->M_Wn);
  if (GLOBAL_W_n_name_[nb]==NULL || !strcmp(GLOBAL_W_n_name_[nb],"\0")){ E->M_Wt = M_handle_v_make(D->bitj,D->Y_ncols,E->Z_nrows,NULL,NULL,D->Y_bmc_b,E->Z_bmr_b);}
  else{ E->M_Wt = M_handle_v_make(D->bitj,D->Y_ncols,E->Z_nrows,GLOBAL_W_n_name_[nb],NULL,D->Y_bmc_b,E->Z_bmr_b);} M_mxget(E->M_Wt);
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
  E->QR_AnZtSWnYt_nrm = NULL; E->QR_AnZtSZnAt_nrm = NULL; E->QR_AnAtTYnYt_nrm = NULL; E->QR_AnAtTAnAt_nrm = NULL; E->QR_imr_a=NULL; E->QR_imr_b=NULL;
  /* initializing a,j,y,t,k,h,u,r,z,v */
  E->M_an = NULL; E->M_at = NULL; E->M_jn = NULL; E->M_jt = NULL; E->M_yn = NULL; E->M_yt; E->M_tt;
  E->M_kn = NULL; E->M_kt = NULL; E->M_hn = NULL; E->M_ht = NULL; E->M_un = NULL; E->M_ut; E->M_rt;
  E->M_zn = NULL; E->M_zt = NULL; E->M_vn = NULL; E->M_vt = NULL;
  E->M_an = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->M_An->mr_b,E->M_An->mc_b);
  E->M_jn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->M_An->mr_b,E->M_An->mc_b);
  E->M_kn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->M_An->mr_b,E->M_An->mc_b);
  E->M_hn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->M_An->mr_b,E->M_An->mc_b);
  E->M_at = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,E->M_At->mr_b,E->M_At->mc_b);
  E->M_jt = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,E->M_At->mr_b,E->M_At->mc_b);
  E->M_kt = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,E->M_At->mr_b,E->M_At->mc_b);
  E->M_ht = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,E->M_At->mr_b,E->M_At->mc_b);
  E->M_yn = M_handle_v_make(D->bitj,E->A_nrows,D->Y_ncols,NULL,E->M_Yn->wX,E->M_Yn->mr_b,E->M_Yn->mc_b);
  E->M_un = M_handle_v_make(D->bitj,E->A_nrows,D->Y_ncols,NULL,E->M_Yn->wX,E->M_Yn->mr_b,E->M_Yn->mc_b);
  E->M_yt = M_handle_v_make(D->bitj,D->Y_ncols,E->A_nrows,NULL,E->M_Yt->wX,E->M_Yt->mr_b,E->M_Yt->mc_b);
  E->M_ut = M_handle_v_make(D->bitj,D->Y_ncols,E->A_nrows,NULL,E->M_Yt->wX,E->M_Yt->mr_b,E->M_Yt->mc_b);
  E->M_zn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,NULL,E->M_Zn->wX,E->M_Zn->mr_b,E->M_Zn->mc_b);
  E->M_vn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,NULL,E->M_Zn->wX,E->M_Zn->mr_b,E->M_Zn->mc_b);
  E->M_zt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,NULL,E->M_Zt->wX,E->M_Zt->mr_b,E->M_Zt->mc_b);
  E->M_vt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,NULL,E->M_Zt->wX,E->M_Zt->mr_b,E->M_Zt->mc_b);
  E->M_tt = M_handle_v_make(D->bitj,D->T_ncols,E->A_nrows,NULL,E->M_Tt->wX,E->M_Tt->mr_b,E->M_Tt->mc_b);
  E->M_rt = M_handle_v_make(D->bitj,D->T_ncols,E->A_nrows,NULL,E->M_Tt->wX,E->M_Tt->mr_b,E->M_Tt->mc_b);
  if (verbose){ printf(" %% [finished bcc_single_load_M_An]\n");}
}

void bcc_single_init_M_An(char *error_vs_speed,double mrnd,struct bcc_single *E)
{
  /* initialize based on test parameters */
  int verbose=0;
  struct bcc_ajdk *D=E->D;
  int nb_given = E->nb;
  int nr=0,nb=0,nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0,nrows_tmp_extend,brows_tmp;
  unsigned char *bXra_tmp=NULL;
  unsigned char *bb=NULL,*bj=NULL,*bn=NULL,*bt=NULL,*bn2=NULL,*bt2=NULL;
  double dtmp=0; int nx=0;
  if (verbose){ printf(" %% [entering bcc_single_init_M_An]\n");}
  nb = nb_given;
  /* initializing A */
  E->A_nrows = maximum(0,GLOBAL_TEST_A_n_rows + (nb_given % 3) - 1);
  E->A_rbother = (E->A_nrows>0);
  E->A_nrows_extend = (D->bitj - (E->A_nrows % D->bitj)) % D->bitj;
  E->A_mr_length = bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */;
  E->A_bmr_b = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_b){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");}
  E->A_bmr_j = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");}
  bb = E->A_bmr_b; bj = E->A_bmr_j; for (nx=0;nx<E->A_nrows;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpAnchar," %%%% E->A_bmr_b: "); if (verbose>2){ bprintf(E->A_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% E->A_bmr_j: "); if (verbose>2){ bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  E->A_bmr_j_rmv = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rmv){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");}
  E->A_bmr_j_rtn = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rtn){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");}
  E->A_rpop_b = popcount_uchar_array(E->A_bmr_b,E->A_mr_length);
  E->A_rpop_j = popcount_uchar_array(E->A_bmr_j,E->A_mr_length);
  E->A_umr_b = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_b){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_b[nr] = bget__on(E->A_bmr_b,nr);}
  E->A_umr_j = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);}
  E->A_umr_j_rmv = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rmv){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rmv[nr] = bget__on(E->A_bmr_j_rmv,nr);}
  E->A_umr_j_rtn = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rtn){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = bget__on(E->A_bmr_j_rtn,nr);}
  if (verbose>1){ printf(" %% generating E->M_An\n");}
  bn = NULL; bt = NULL; E->M_An = NULL; E->M_At = NULL;
  wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->A_ncols,E->A_bmr_b,D->A_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->A_ncols,bn,bt,E->A_bmr_b,E->A_bmr_j,D->A_bmc_b,D->A_bmc_j,&(E->M_An),&(E->M_At));
  /* initializing Z */
  E->Z_nrows = maximum(0,GLOBAL_TEST_Z_n_rows + (nb_given % 5) - 2);
  E->Z_rbother = (E->Z_nrows>0);
  E->Z_nrows_extend = (D->bitj - (E->Z_nrows % D->bitj)) % D->bitj;
  E->Z_mr_length = bsize(E->Z_nrows)/* rup(E->Z_nrows+E->Z_nrows_extend,POPLENGTH)/BIT8 */;
  E->Z_bmr_b = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_b){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");}
  E->Z_bmr_j = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_j){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");}
  bb = E->Z_bmr_b; bj = E->Z_bmr_j; for (nx=0;nx<E->Z_nrows;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_b: "); if (verbose>2){ bprintf(E->Z_bmr_b,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_j: "); if (verbose>2){ bprintf(E->Z_bmr_j,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  E->Z_rpop_b = popcount_uchar_array(E->Z_bmr_b,E->Z_mr_length);
  E->Z_rpop_j = popcount_uchar_array(E->Z_bmr_j,E->Z_mr_length);
  E->Z_umr_b = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_b){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_b[nr] = bget__on(E->Z_bmr_b,nr);}
  E->Z_umr_j = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_j){ printf(" %% Warning! not enough memory in bcc_single_init_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_j[nr] = bget__on(E->Z_bmr_j,nr);}
  if (verbose>1){ printf(" %% generating E->M_Zn\n");}
  bn = NULL; bt = NULL; E->M_Zn = NULL; E->M_Zt = NULL;
  wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->Z_nrows,D->A_ncols,E->Z_bmr_b,D->A_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->Z_nrows,D->A_ncols,bn,bt,E->Z_bmr_b,E->Z_bmr_j,D->A_bmc_b,D->A_bmc_j,&(E->M_Zn),&(E->M_Zt));
  /* initializing Y */
  bn = NULL; bt = NULL; E->M_Yn = NULL; E->M_Yt = NULL;
  wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->Y_ncols,E->A_bmr_b,D->Y_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->Y_ncols,bn,bt,E->A_bmr_b,E->A_bmr_j,D->Y_bmc_b,D->Y_bmc_j,&(E->M_Yn),&(E->M_Yt));
  /* initializing W */
  bn = NULL; bt = NULL; E->M_Wn = NULL; E->M_Wt = NULL;
  wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->Z_nrows,D->Y_ncols,E->Z_bmr_b,D->Y_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->Z_nrows,D->Y_ncols,bn,bt,E->Z_bmr_b,E->Z_bmr_j,D->Y_bmc_b,D->Y_bmc_j,&(E->M_Wn),&(E->M_Wt));
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
  /* initializing E->QR_AnZtSWnYt_nrm */
  E->QR_AnZtSWnYt_nrm = NULL; E->QR_AnZtSZnAt_nrm = NULL; E->QR_AnAtTYnYt_nrm = NULL; E->QR_AnAtTAnAt_nrm = NULL; E->QR_imr_a=NULL; E->QR_imr_b=NULL;
  /* initializing a,j,y,t,k,h,u,r,z,v */
  E->M_an = NULL; E->M_at = NULL; E->M_jn = NULL; E->M_jt = NULL; E->M_yn = NULL; E->M_yt; E->M_tt;
  E->M_kn = NULL; E->M_kt = NULL; E->M_hn = NULL; E->M_ht = NULL; E->M_un = NULL; E->M_ut; E->M_rt;
  E->M_zn = NULL; E->M_zt = NULL; E->M_vn = NULL; E->M_vt = NULL;
  E->M_an = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->M_An->mr_b,E->M_An->mc_b);
  E->M_jn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->M_An->mr_b,E->M_An->mc_b);
  E->M_kn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->M_An->mr_b,E->M_An->mc_b);
  E->M_hn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->M_An->mr_b,E->M_An->mc_b);
  E->M_at = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,E->M_At->mr_b,E->M_At->mc_b);
  E->M_jt = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,E->M_At->mr_b,E->M_At->mc_b);
  E->M_kt = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,E->M_At->mr_b,E->M_At->mc_b);
  E->M_ht = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,E->M_At->mr_b,E->M_At->mc_b);
  E->M_yn = M_handle_v_make(D->bitj,E->A_nrows,D->Y_ncols,NULL,E->M_Yn->wX,E->M_Yn->mr_b,E->M_Yn->mc_b);
  E->M_un = M_handle_v_make(D->bitj,E->A_nrows,D->Y_ncols,NULL,E->M_Yn->wX,E->M_Yn->mr_b,E->M_Yn->mc_b);
  E->M_yt = M_handle_v_make(D->bitj,D->Y_ncols,E->A_nrows,NULL,E->M_Yt->wX,E->M_Yt->mr_b,E->M_Yt->mc_b);
  E->M_ut = M_handle_v_make(D->bitj,D->Y_ncols,E->A_nrows,NULL,E->M_Yt->wX,E->M_Yt->mr_b,E->M_Yt->mc_b);
  E->M_zn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,NULL,E->M_Zn->wX,E->M_Zn->mr_b,E->M_Zn->mc_b);
  E->M_vn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,NULL,E->M_Zn->wX,E->M_Zn->mr_b,E->M_Zn->mc_b);
  E->M_zt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,NULL,E->M_Zt->wX,E->M_Zt->mr_b,E->M_Zt->mc_b);
  E->M_vt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,NULL,E->M_Zt->wX,E->M_Zt->mr_b,E->M_Zt->mc_b);
  E->M_tt = M_handle_v_make(D->bitj,D->T_ncols,E->A_nrows,NULL,E->M_Tt->wX,E->M_Tt->mr_b,E->M_Tt->mc_b);
  E->M_rt = M_handle_v_make(D->bitj,D->T_ncols,E->A_nrows,NULL,E->M_Tt->wX,E->M_Tt->mr_b,E->M_Tt->mc_b);
  if (verbose){ printf(" %% [finished bcc_single_init_M_An]\n");}
}

void bcc_single_copy_lf(struct bcc_single *E,struct bcc_single *E_in)
{
  int verbose=0; /* double gamma = GLOBAL_gamma; */
  struct bcc_ajdk *D = E->D;
  int ns_b=0,n1a=0,n1b=0,n2a=0,n2b=0;
  char prefix[FNAMESIZE];
  int M_flag=0;
  int cdrop=0,ckeep=0;
  long long int mm=0;
  L_handle_copy(E->lf_At_rsum,E_in->lf_At_rsum);
  L_handle_copy(E->lf_Zt_rsum,E_in->lf_Zt_rsum);
  L_handle_copy(E->lf_Wt_rsum,E_in->lf_Wt_rsum);
  L_handle_copy(E->lf_Yt_rsum,E_in->lf_Yt_rsum);
  L_handle_copy(E->lf_AtTn,E_in->lf_AtTn);
  L_handle_copy(E->lf_ZtSn,E_in->lf_ZtSn);
  L_handle_copy(E->lf_YtTn,E_in->lf_YtTn);
  L_handle_copy(E->lf_WtSn,E_in->lf_WtSn);
  L_handle_copy(E->lf_An_ajdk,E_in->lf_An_ajdk);
  L_handle_copy(E->lf_Zn_ajdk,E_in->lf_Zn_ajdk);
  L_handle_copy(E->lf_Yn_ajdk,E_in->lf_Yn_ajdk);
  L_handle_copy(E->lf_Wn_ajdk,E_in->lf_Wn_ajdk);
  if (strstr(D->QR_strategy,"ZtSWn") || strstr(D->QC_strategy,"ZtSWn")){
    L_handle_copy(E->lf_ZtSWn,E_in->lf_ZtSWn);
    L_handle_copy(E->lf_ZtSZn,E_in->lf_ZtSZn);
    L_handle_copy(E->lf_AtTYn,E_in->lf_AtTYn);
    L_handle_copy(E->lf_AtTAn,E_in->lf_AtTAn);
    /* if strategy */}
  if (strstr(D->QR_strategy,"ZtSWn")){
    n1a = E->M_Zt->nrows; n1b = E->M_Zt->rpop_b; n2a = E->M_Wt->nrows; n2b = E->M_Wt->rpop_b;
    if (0){}
    else if (n2b<=n1b){ E->M_ZtSWn_trm = 0; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ M_handle_copy(E->M_ZtSWn_[ns_b],E_in->M_ZtSWn_[ns_b]);}}
    else if (n1b<=n2b){ E->M_ZtSWn_trm = 1; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ M_handle_copy(E->M_ZtSWn_[ns_b],E_in->M_ZtSWn_[ns_b]);}}
    n1a = E->M_Zt->nrows; n1b = E->M_Zt->rpop_b; n2a = E->M_Zt->nrows; n2b = E->M_Zt->rpop_b;
    if (0){}
    else if (n2b<=n1b){ E->M_ZtSZn_trm = 0; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ M_handle_copy(E->M_ZtSZn_[ns_b],E_in->M_ZtSZn_[ns_b]);}}
    else if (n1b<=n2b){ E->M_ZtSZn_trm = 1; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ M_handle_copy(E->M_ZtSZn_[ns_b],E_in->M_ZtSZn_[ns_b]);}}
    n1a = E->M_At->nrows; n1b = E->M_At->rpop_b; n2a = E->M_Yt->nrows; n2b = E->M_Yt->rpop_b;
    if (0){}
    else if (n2b<=n1b){ E->M_AtTYn_trm = 0; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ M_handle_copy(E->M_AtTYn_[ns_b],E_in->M_AtTYn_[ns_b]);}}
    else if (n1b<=n2b){ E->M_AtTYn_trm = 1; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ M_handle_copy(E->M_AtTYn_[ns_b],E_in->M_AtTYn_[ns_b]);}}
    n1a = E->M_At->nrows; n1b = E->M_At->rpop_b; n2a = E->M_At->nrows; n2b = E->M_At->rpop_b;
    if (0){}
    else if (n2b<=n1b){ E->M_AtTAn_trm = 0; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ M_handle_copy(E->M_AtTAn_[ns_b],E_in->M_AtTAn_[ns_b]);}}
    else if (n1b<=n2b){ E->M_AtTAn_trm = 1; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ M_handle_copy(E->M_AtTAn_[ns_b],E_in->M_AtTAn_[ns_b]);}}
    /* if strategy */}
  if (verbose>1){ printf(" %% copying E->lrup\n");}
  get_xdrop(D->A_rpop_j_total,D->A_cpop_j,NULL,&cdrop); ckeep=D->A_cpop_j-cdrop;
  /* QR YnWt lrup */
  if (strstr(D->QR_strategy,"YnWt")){
    L_handle_copy(E->lf_jn_ajdk,E_in->lf_jn_ajdk);
    L_handle_copy(E->lf_vn_ajdk,E_in->lf_vn_ajdk);
    /* if strategy */}
  /* QC ZtSWn lrup */
  if (strstr(D->QC_strategy,"ZtSWn")){
    L_handle_copy(E->lf_ktrn,E_in->lf_ktrn);
    L_handle_copy(E->lf_utrn,E_in->lf_utrn);
    /* if strategy */}
  if (strstr(D->QC_strategy,"ZtSWn")){
    L_handle_copy(E->lf_ktrun,E_in->lf_ktrun);
    L_handle_copy(E->lf_ktrkn,E_in->lf_ktrkn);
    /* if strategy */}  
  /* QC YnWt lrup */
  if (strstr(D->QC_strategy,"YnWt")){
    L_handle_copy(E->lf_kn_ajdk,E_in->lf_kn_ajdk);
    L_handle_copy(E->lf_zn_ajdk,E_in->lf_zn_ajdk);
    L_handle_copy(E->lf_hn_ajdk,E_in->lf_hn_ajdk);
    L_handle_copy(E->lf_vn_ajdk,E_in->lf_vn_ajdk);
    L_handle_copy(E->lf_an_ajdk,E_in->lf_an_ajdk);
    /* if strategy */}
  /* QC YnWt lrup */
  if (strstr(D->QC_strategy,"YnWt")){
    L_handle_copy(E->lf_attn,E_in->lf_attn);
    L_handle_copy(E->lf_jttn,E_in->lf_jttn);
    L_handle_copy(E->lf_ktrn,E_in->lf_ktrn);
    L_handle_copy(E->lf_htrn,E_in->lf_htrn);
    L_handle_copy(E->lf_ztsn,E_in->lf_ztsn);
    L_handle_copy(E->lf_vtsn,E_in->lf_vtsn);
    /* if strategy */}
  if (strstr(D->QC_strategy,"YnWt")){
    if (strstr(D->QC_strategy,"store one")){
      L_handle_copy(E->lf_ztsvn,E_in->lf_ztsvn);
      L_handle_copy(E->lf_attjn,E_in->lf_attjn);
      L_handle_copy(E->lf_ktrhn,E_in->lf_ktrhn);
      /* if strategy */}
    else /* store all */{
      L_handle_copy(E->lf_ztsvn,E_in->lf_ztsvn);
      L_handle_copy(E->lf_attjn,E_in->lf_attjn);
      L_handle_copy(E->lf_ktrhn,E_in->lf_ktrhn);
      /* if strategy */}
    /* if strategy */}
  if (verbose>1){ printf(" %% copying singlestudy\n");}
  L_handle_copy(E->lf_At_Yn_a1d1,E_in->lf_At_Yn_a1d1);
  L_handle_copy(E->lf_At_An_a1d1,E_in->lf_At_An_a1d1);
  L_handle_copy(E->lf_a1d2_Zt_Sn,E_in->lf_a1d2_Zt_Sn);
  L_handle_copy(E->lf_a2d2_Zt_Sn,E_in->lf_a2d2_Zt_Sn);
  L_handle_copy(E->lf_a3d2_Zt_Sn,E_in->lf_a3d2_Zt_Sn);
  L_handle_copy(E->lf_a1d2_At_Tn,E_in->lf_a1d2_At_Tn);
  L_handle_copy(E->lf_a2d2_At_Tn,E_in->lf_a2d2_At_Tn);
  L_handle_copy(E->lf_a3d2_At_Tn,E_in->lf_a3d2_At_Tn);
  L_handle_copy(E->lf_et_Tn,E_in->lf_et_Tn);
  L_handle_copy(E->lf_et_Sn,E_in->lf_et_Sn);
  L_handle_copy(E->lf_et_An,E_in->lf_et_An);
  L_handle_copy(E->lf_a1d1_At_en,E_in->lf_a1d1_At_en);
  L_handle_copy(E->lf_a1d2_At_en,E_in->lf_a1d2_At_en);
  L_handle_copy(E->lf_a3d2_At_en,E_in->lf_a3d2_At_en);
  L_handle_copy(E->lf_et_Yn_a1d1,E_in->lf_et_Yn_a1d1);
  L_handle_copy(E->lf_Tt_Yn_a1d1,E_in->lf_Tt_Yn_a1d1);
  L_handle_copy(E->lf_et_An_a1d1,E_in->lf_et_An_a1d1);
  L_handle_copy(E->lf_Tt_An_a1d1,E_in->lf_Tt_An_a1d1);
  L_handle_copy(E->lf_T_AnAt_YnYt,E_in->lf_T_AnAt_YnYt);
  L_handle_copy(E->lf_T_AnAt_AnAt,E_in->lf_T_AnAt_AnAt);
  if (verbose>1){ printf(" %% copying storage for singlestudy\n");}
  M_handle_copy(E->M_1,E_in->M_1);
  L_handle_copy(E->lf_Yn_a1d1,E_in->lf_Yn_a1d1); M_handle_copy(E->M_Yn_a1d1,E_in->M_Yn_a1d1);
  L_handle_copy(E->lf_An_a1d1,E_in->lf_An_a1d1); M_handle_copy(E->M_An_a1d1,E_in->M_An_a1d1);
  L_handle_copy(E->lf_a2d2_ZtSn,E_in->lf_a2d2_ZtSn); M_handle_copy(E->M_a2d2_ZtSn,E_in->M_a2d2_ZtSn);
  L_handle_copy(E->lf_a2d2_AtTn,E_in->lf_a2d2_AtTn); M_handle_copy(E->M_a2d2_AtTn,E_in->M_a2d2_AtTn);
}

void bcc_single_init_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,struct M_handle **M_p)
{
  int verbose=0;
  int length=0;
  length  = n1b*n2b; *L_p = L_handle_make(length); if (!(*L_p)){ printf(" %% Warning! not enough memory in bcc_double_init_lf\n");}
  if (M_flag){ *M_p = M_handle_w_make(BITJ,GLOBAL_B_MLT,n1a,n2b);}
  if (verbose){ printf(" %% %s n1a %d n1b %d n2a %d n2b %d --> length %d \n",prefix,n1a,n1b,n2a,n2b,length);}
}

void bcc_single_init_lf(struct bcc_single *E)
{
  int verbose=0; /* double gamma = GLOBAL_gamma; */
  struct bcc_ajdk *D = E->D;
  int ns_b=0,n1a=0,n1b=0,n2a=0,n2b=0;
  char prefix[FNAMESIZE];
  int M_flag=0;
  int cdrop=0,ckeep=0;
  long long int mm=0;
  E->lf_At_rsum=NULL; E->lf_Zt_rsum=NULL; E->lf_Wt_rsum=NULL; E->lf_Yt_rsum=NULL;
  E->length = D->A_ncols; E->lf_At_rsum = L_handle_make(E->length); if (!E->lf_At_rsum){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  E->length = D->A_ncols; E->lf_Zt_rsum = L_handle_make(E->length); if (!E->lf_Zt_rsum){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  E->length = D->Y_ncols; E->lf_Wt_rsum = L_handle_make(E->length); if (!E->lf_Wt_rsum){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  E->length = D->Y_ncols; E->lf_Yt_rsum = L_handle_make(E->length); if (!E->lf_Yt_rsum){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  if (verbose){ printf(" %% initializing E->lf_AtTn, E->lf_ZtSn, E->lf_YtTn, E->lf_WtSn.\n");}
  E->lf_AtTn = NULL; E->lf_ZtSn = NULL; E->lf_YtTn = NULL; E->lf_WtSn = NULL;
  E->length = D->A_ncols*D->T_ncols; E->lf_AtTn = L_handle_make(E->length); if (!E->lf_AtTn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  E->length = D->A_ncols*D->T_ncols; E->lf_ZtSn = L_handle_make(E->length); if (!E->lf_ZtSn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  E->length = D->Y_ncols*D->T_ncols; E->lf_YtTn = L_handle_make(E->length); if (!E->lf_YtTn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  E->length = D->Y_ncols*D->T_ncols; E->lf_WtSn = L_handle_make(E->length); if (!E->lf_WtSn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing E->lf_An_ajdk, etc\n");}
  E->lf_An_ajdk = NULL; E->lf_Zn_ajdk = NULL; E->lf_Yn_ajdk = NULL; E->lf_Wn_ajdk = NULL;
  E->length = E->A_nrows*AJDK_TOT; E->lf_An_ajdk = L_handle_make(E->length); if (!E->lf_An_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  E->length = E->Z_nrows*AJDK_TOT; E->lf_Zn_ajdk = L_handle_make(E->length); if (!E->lf_Zn_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  E->length = E->A_nrows*AJDK_TOT; E->lf_Yn_ajdk = L_handle_make(E->length); if (!E->lf_Yn_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  E->length = E->Z_nrows*AJDK_TOT; E->lf_Wn_ajdk = L_handle_make(E->length); if (!E->lf_Wn_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing E->lf_ZtSWn, etc\n");}
  E->lf_ZtSWn = NULL; E->lf_ZtSZn = NULL; E->lf_AtTYn = NULL; E->lf_AtTAn = NULL; 
  if (strstr(D->QR_strategy,"ZtSWn") || strstr(D->QC_strategy,"ZtSWn")){
    E->length = E->Z_rbother*E->M_Zn->cpop_b*E->M_Wn->cpop_b*D->T_cpop_b; E->lf_ZtSWn = L_handle_make(E->length); if (!E->lf_ZtSWn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = E->Z_rbother*E->M_Zn->cpop_b*E->M_Zn->cpop_b*D->T_cpop_b; E->lf_ZtSZn = L_handle_make(E->length); if (!E->lf_ZtSZn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = E->A_rbother*E->M_An->cpop_b*E->M_Yn->cpop_b*D->T_cpop_b; E->lf_AtTYn = L_handle_make(E->length); if (!E->lf_AtTYn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = E->A_rbother*E->M_An->cpop_b*E->M_An->cpop_b*D->T_cpop_b; E->lf_AtTAn = L_handle_make(E->length); if (!E->lf_AtTAn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    /* if strategy */}
  E->M_ZtSWn_ = NULL; E->M_ZtSZn_ = NULL; E->M_AtTYn_ = NULL; E->M_AtTAn_ = NULL;
  if (strstr(D->QR_strategy,"ZtSWn")){
    E->M_ZtSWn_ = (struct M_handle **) wkspace_all0c(sizeof(struct M_handle *)*D->T_cpop_b);
    n1a = E->M_Zt->nrows; n1b = E->M_Zt->rpop_b; n2a = E->M_Wt->nrows; n2b = E->M_Wt->rpop_b;
    if (0){}
    else if (n2b<=n1b){ E->M_ZtSWn_trm = 0; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ E->M_ZtSWn_[ns_b] = M_handle_w_make(BITJ,GLOBAL_B_MLT,n1a,n2b);}}
    else if (n1b<=n2b){ E->M_ZtSWn_trm = 1; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ E->M_ZtSWn_[ns_b] = M_handle_w_make(BITJ,GLOBAL_B_MLT,n2a,n1b);}}
    E->M_ZtSZn_ = (struct M_handle **) wkspace_all0c(sizeof(struct M_handle *)*D->T_cpop_b);
    n1a = E->M_Zt->nrows; n1b = E->M_Zt->rpop_b; n2a = E->M_Zt->nrows; n2b = E->M_Zt->rpop_b;
    if (0){}
    else if (n2b<=n1b){ E->M_ZtSZn_trm = 0; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ E->M_ZtSZn_[ns_b] = M_handle_w_make(BITJ,GLOBAL_B_MLT,n1a,n2b);}}
    else if (n1b<=n2b){ E->M_ZtSZn_trm = 1; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ E->M_ZtSZn_[ns_b] = M_handle_w_make(BITJ,GLOBAL_B_MLT,n2a,n1b);}}
    E->M_AtTYn_ = (struct M_handle **) wkspace_all0c(sizeof(struct M_handle *)*D->T_cpop_b);
    n1a = E->M_At->nrows; n1b = E->M_At->rpop_b; n2a = E->M_Yt->nrows; n2b = E->M_Yt->rpop_b;
    if (0){}
    else if (n2b<=n1b){ E->M_AtTYn_trm = 0; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ E->M_AtTYn_[ns_b] = M_handle_w_make(BITJ,GLOBAL_B_MLT,n1a,n2b);}}
    else if (n1b<=n2b){ E->M_AtTYn_trm = 1; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ E->M_AtTYn_[ns_b] = M_handle_w_make(BITJ,GLOBAL_B_MLT,n2a,n1b);}}
    E->M_AtTAn_ = (struct M_handle **) wkspace_all0c(sizeof(struct M_handle *)*D->T_cpop_b);
    n1a = E->M_At->nrows; n1b = E->M_At->rpop_b; n2a = E->M_At->nrows; n2b = E->M_At->rpop_b;
    if (0){}
    else if (n2b<=n1b){ E->M_AtTAn_trm = 0; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ E->M_AtTAn_[ns_b] = M_handle_w_make(BITJ,GLOBAL_B_MLT,n1a,n2b);}}
    else if (n1b<=n2b){ E->M_AtTAn_trm = 1; for (ns_b=0;ns_b<D->T_cpop_b;ns_b++){ E->M_AtTAn_[ns_b] = M_handle_w_make(BITJ,GLOBAL_B_MLT,n2a,n1b);}}
    /* if strategy */}
  if (verbose>1){ printf(" %% initializing E->lrup\n");}
  get_xdrop(D->A_rpop_j_total,D->A_cpop_j,NULL,&cdrop); ckeep=D->A_cpop_j-cdrop;
  /* QR YnWt lrup */
  E->lf_jn_ajdk=NULL;E->lf_vn_ajdk=NULL;
  if (strstr(D->QR_strategy,"YnWt")){
    E->length = E->A_nrows*AJDK_TOT; E->lf_jn_ajdk = L_handle_make(E->length); if (!E->lf_jn_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = E->Z_nrows*AJDK_TOT; E->lf_vn_ajdk = L_handle_make(E->length); if (!E->lf_vn_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    /* if strategy */}
  /* QC ZtSWn lrup */
  E->lf_ktrn=NULL;E->lf_utrn=NULL;
  if (strstr(D->QC_strategy,"ZtSWn")){
    E->length = D->A_ncols*D->T_ncols; E->lf_ktrn = L_handle_make(E->length); if (!E->lf_ktrn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = D->Y_ncols*D->T_ncols; E->lf_utrn = L_handle_make(E->length); if (!E->lf_utrn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    /* if strategy */}
  E->lf_ktrun=NULL;E->lf_ktrkn=NULL;
  if (strstr(D->QC_strategy,"ZtSWn")){
    E->length = E->A_rbother*E->M_An->cpop_b*E->M_Yn->cpop_b*D->T_cpop_b; E->lf_ktrun = L_handle_make(E->length); if (!E->lf_ktrun){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = E->A_rbother*E->M_An->cpop_b*E->M_An->cpop_b*D->T_cpop_b; E->lf_ktrkn = L_handle_make(E->length) ;if (!E->lf_ktrkn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    /* if strategy */}  
  /* QC YnWt lrup */
  E->lf_kn_ajdk=NULL; E->lf_zn_ajdk=NULL; E->lf_hn_ajdk=NULL; /* E->lf_vn_ajdk=NULL; */ E->lf_an_ajdk=NULL;
  if (strstr(D->QC_strategy,"YnWt")){
    E->length = E->A_nrows*AJDK_TOT; E->lf_kn_ajdk = L_handle_make(E->length); if (!E->lf_kn_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = E->Z_nrows*AJDK_TOT; E->lf_zn_ajdk = L_handle_make(E->length); if (!E->lf_zn_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = E->A_nrows*AJDK_TOT; E->lf_hn_ajdk = L_handle_make(E->length); if (!E->lf_hn_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    if (E->lf_vn_ajdk==NULL){ E->length = E->Z_nrows*AJDK_TOT; E->lf_vn_ajdk = L_handle_make(E->length); if (!E->lf_vn_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}}
    E->length = E->A_nrows*AJDK_TOT; E->lf_an_ajdk = L_handle_make(E->length); if (!E->lf_an_ajdk){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    /* if strategy */}
  /* QC YnWt lrup */
  E->lf_attn=NULL; E->lf_jttn=NULL; /* E->lf_ktrn=NULL; */ E->lf_htrn=NULL; E->lf_ztsn=NULL; E->lf_vtsn=NULL;
  if (strstr(D->QC_strategy,"YnWt")){
    E->length = D->A_ncols*D->T_ncols; E->lf_attn = L_handle_make(E->length); if (!E->lf_attn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = D->A_ncols*D->T_ncols; E->lf_jttn = L_handle_make(E->length); if (!E->lf_jttn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    if (E->lf_ktrn==NULL){ E->length = D->A_ncols*D->T_ncols; E->lf_ktrn = L_handle_make(E->length); if (!E->lf_ktrn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}}
    E->length = D->A_ncols*D->T_ncols; E->lf_htrn = L_handle_make(E->length); if (!E->lf_htrn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = D->A_ncols*D->T_ncols; E->lf_ztsn = L_handle_make(E->length); if (!E->lf_ztsn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    E->length = D->A_ncols*D->T_ncols; E->lf_vtsn = L_handle_make(E->length); if (!E->lf_vtsn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
    /* if strategy */}
  E->lf_ztsvn=NULL; E->lf_attjn=NULL; E->lf_ktrhn=NULL;
  if (strstr(D->QC_strategy,"YnWt")){
    if (strstr(D->QC_strategy,"store one")){
      E->length = E->Z_rbother*ckeep*cdrop; E->lf_ztsvn = L_handle_make(E->length); if (!E->lf_ztsvn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
      mm = 8*(double)E->length;
      if (verbose || mm>1024*1024*256){ printf(" %% E->length %d*%d --> %d; memory required %0.2fG\n",ckeep,cdrop,E->length,8*(double)E->length/1e9);}
      E->length = E->A_rbother*ckeep*cdrop; E->lf_attjn = L_handle_make(E->length); if (!E->lf_attjn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
      mm = 8*(double)E->length;
      if (verbose || mm>1024*1024*256){ printf(" %% E->length %d*%d --> %d; memory required %0.2fG\n",ckeep,cdrop,E->length,8*(double)E->length/1e9);}
      E->length = E->A_rbother*ckeep*cdrop; E->lf_ktrhn = L_handle_make(E->length); if (!E->lf_ktrhn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
      mm = 8*(double)E->length;
      if (verbose || mm>1024*1024*256){ printf(" %% E->length %d*%d --> %d; memory required %0.2fG\n",ckeep,cdrop,E->length,8*(double)E->length/1e9);}
      /* if strategy */}
    else /* store all */{
      E->length = E->Z_rbother*ckeep*cdrop*D->T_ncols; E->lf_ztsvn = L_handle_make(E->length); if (!E->lf_ztsvn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
      mm = 8*(double)E->length;
      if (verbose || mm>1024*1024*256){ printf(" %% E->length %d*%d*%d --> %d; memory required %0.2fG\n",ckeep,cdrop,D->T_ncols,E->length,8*(double)E->length/1e9);}
      E->length = E->A_rbother*ckeep*cdrop*D->T_ncols; E->lf_attjn = L_handle_make(E->length); if (!E->lf_attjn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
      mm = 8*(double)E->length;
      if (verbose || mm>1024*1024*256){ printf(" %% E->length %d*%d*%d --> %d; memory required %0.2fG\n",ckeep,cdrop,D->T_ncols,E->length,8*(double)E->length/1e9);}
      E->length = E->A_rbother*ckeep*cdrop*D->T_ncols; E->lf_ktrhn = L_handle_make(E->length); if (!E->lf_ktrhn){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
      mm = 8*(double)E->length;
      if (verbose || mm>1024*1024*256){ printf(" %% E->length %d*%d*%d --> %d; memory required %0.2fG\n",ckeep,cdrop,D->T_ncols,E->length,8*(double)E->length/1e9);}
      /* if strategy */}
    /* if strategy */}
  if (verbose>1){ printf(" %% initializing singlestudy\n");}
  E->length = D->A_ncols; E->lf_At_Yn_a1d1 = L_handle_make(E->length);
  E->length = D->A_ncols; E->lf_At_An_a1d1 = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_a1d2_Zt_Sn = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_a2d2_Zt_Sn = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_a3d2_Zt_Sn = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_a1d2_At_Tn = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_a2d2_At_Tn = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_a3d2_At_Tn = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_et_Tn = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_et_Sn = L_handle_make(E->length);
  E->length = D->A_ncols; E->lf_et_An = L_handle_make(E->length);
  E->length = D->A_ncols; E->lf_a1d1_At_en = L_handle_make(E->length);
  E->length = D->A_ncols; E->lf_a1d2_At_en = L_handle_make(E->length);
  E->length = D->A_ncols; E->lf_a3d2_At_en = L_handle_make(E->length);
  E->length = 1; E->lf_et_Yn_a1d1 = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_Tt_Yn_a1d1 = L_handle_make(E->length);
  /* E->length = D->A_ncols*D->T_ncols; E->lf_At_T_Yn_a1d1 = L_handle_make(E->length); */
  E->length = 1; E->lf_et_An_a1d1 = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_Tt_An_a1d1 = L_handle_make(E->length);
  /* E->length = D->A_ncols*D->T_ncols; E->lf_At_T_An_a1d1 = L_handle_make(E->length); */
  E->length = D->T_ncols*E->A_nrows; E->lf_T_AnAt_YnYt = L_handle_make(E->length);
  E->length = D->T_ncols*E->A_nrows; E->lf_T_AnAt_AnAt = L_handle_make(E->length); if (!E->lf_T_AnAt_AnAt){ printf(" %% Warning! not enough memory in bcc_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing storage for singlestudy\n");}
  E->M_1 = M_handle_v_make(D->bitj,1,1,NULL,NULL,NULL,NULL);
  M_flag = 1;
  sprintf(prefix,"E->M_Yn_a1d1");
  n1a = D->Y_cbother*E->A_nrows; n1b = D->Y_cbother*E->A_rpop_b;
  n2a = 1; n2b = 1; 
  bcc_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_Yn_a1d1),M_flag,&(E->M_Yn_a1d1));
  M_flag = 1;
  sprintf(prefix,"E->M_An_a1d1");
  n1a = D->A_cbother*E->A_nrows; n1b = D->A_cbother*E->A_rpop_b;
  n2a = 1; n2b = 1; 
  bcc_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_An_a1d1),M_flag,&(E->M_An_a1d1));
  if (verbose>1){ printf(" %% initializing storage for doublestudy\n");}
  M_flag = 1;
  sprintf(prefix,"E->M_a2d2_Zt_Sn");
  n1a = D->A_cbother*D->A_ncols; n1b = D->A_cbother*D->A_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  bcc_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a2d2_ZtSn),M_flag,&(E->M_a2d2_ZtSn));
  M_flag = 1;
  sprintf(prefix,"E->M_a2d2_At_Tn");
  n1a = D->A_cbother*D->A_ncols; n1b = D->A_cbother*D->A_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  bcc_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a2d2_AtTn),M_flag,&(E->M_a2d2_AtTn));
}

void bcc_double_copy_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,int *trm_flag_p,struct M_handle **M_p,struct L_handle **L_in_p,struct M_handle **M_in_p)
{
  int verbose=0;
  int length=0;
  length  = n1b*n2b; L_handle_copy(*L_p,*L_in_p);
  if (M_flag){
    if (0){}
    else if (n2b<=n1b){ *trm_flag_p = 0; M_handle_copy(*M_p,*M_in_p);}
    else if (n1b<=n2b){ *trm_flag_p = 1; M_handle_copy(*M_p,*M_in_p);}
    /* if (M_flag){ } */}
  if (verbose){ printf(" %% %s n1a %d n1b %d n2a %d n2b %d --> length %d trm_flag %d \n",prefix,n1a,n1b,n2a,n2b,length,*trm_flag_p);}
}

void bcc_double_copy_lf(struct bcc_double *F,struct bcc_double *F_in)
{
  int verbose=0;
  struct bcc_ajdk *D=F->D;
  int n1a=0,n2a=0,n1b,n2b=0;
  int M_flag=0,rdrop=0,rdrop2=0;
  char prefix[FNAMESIZE];
  GLOBAL_tic(0);
  if (verbose>1){ printf(" %% bcc_double_copy_lf: %d,%d\n",F->nb1,F->nb2);}
  if (verbose>1){ printf(" %% copying F->lf_YnWt, etc\n");}
  if (strstr(D->QR_strategy,"YnWt") || strstr(D->QC_strategy,"YnWt")){
    M_flag = (strstr(D->QC_strategy,"YnWt") ? 1 : 0);
    n1a = F->E_nb1->M_Yn->nrows; n1b = F->E_nb1->M_Yn->rpop_b; n2a = F->E_nb2->M_Wn->nrows; n2b = F->E_nb2->M_Wn->rpop_b; sprintf(prefix,"F->M_YnWt");
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_YnWt),M_flag,&(F->M_YnWt_trm),&(F->M_YnWt),&(F_in->lf_YnWt),&(F_in->M_YnWt));
    n1a = F->E_nb1->M_An->nrows; n1b = F->E_nb1->M_An->rpop_b; n2a = F->E_nb2->M_Zn->nrows; n2b = F->E_nb2->M_Zn->rpop_b; sprintf(prefix,"F->M_AnZt");
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_AnZt),M_flag,&(F->M_AnZt_trm),&(F->M_AnZt),&(F_in->lf_AnZt),&(F_in->M_AnZt));
    n1a = F->E_nb1->M_Yn->nrows; n1b = F->E_nb1->M_Yn->rpop_b; n2a = F->E_nb2->M_Yn->nrows; n2b = F->E_nb2->M_Yn->rpop_b; sprintf(prefix,"F->M_YnYt");
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_YnYt),M_flag,&(F->M_YnYt_trm),&(F->M_YnYt),&(F_in->lf_YnYt),&(F_in->M_YnYt));
    n1a = F->E_nb1->M_An->nrows; n1b = F->E_nb1->M_An->rpop_b; n2a = F->E_nb2->M_An->nrows; n2b = F->E_nb2->M_An->rpop_b; sprintf(prefix,"F->M_AnAt");
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_AnAt),M_flag,&(F->M_AnAt_trm),&(F->M_AnAt),&(F_in->lf_AnAt),&(F_in->M_AnAt));
    /* if strategy */}
  get_xdrop(D->A_rpop_j_total,D->A_cpop_j,&rdrop,NULL); rdrop2 = maximum(2,rdrop);
  /* QR YnWt lrup */
  if (strstr(D->QR_strategy,"YnWt")){
    L_handle_copy(F->lf_jnvt,F_in->lf_jnvt);
    L_handle_copy(F->lf_jnjt,F_in->lf_jnjt);
    /* if strategy */}
  /* QC YnWt lrup */
  if (strstr(D->QC_strategy,"YnWt")){
    M_flag = 1;
    sprintf(prefix,"F->M_unwt");
    n1a = D->Y_cbother*F->E_nb1->M_An->nrows; n1b = D->Y_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->Y_cbother*F->E_nb2->M_Zn->nrows; n2b = D->Y_cbother*F->E_nb2->M_Zn->rpop_b; 
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_unwt),M_flag,&(F->M_unwt_trm),&(F->M_unwt),&(F_in->lf_unwt),&(F_in->M_unwt));
    sprintf(prefix,"F->M_knzt");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->A_cbother*F->E_nb2->M_Zn->nrows; n2b = D->A_cbother*F->E_nb2->M_Zn->rpop_b; 
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_knzt),M_flag,&(F->M_knzt_trm),&(F->M_knzt),&(F_in->lf_knzt),&(F_in->M_knzt));
    sprintf(prefix,"F->M_hnvt");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->A_cbother*F->E_nb2->M_Zn->nrows; n2b = D->A_cbother*F->E_nb2->M_Zn->rpop_b; 
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_hnvt),M_flag,&(F->M_hnvt_trm),&(F->M_hnvt),&(F_in->lf_hnvt),&(F_in->M_hnvt));
    sprintf(prefix,"F->M_ynut");
    n1a = D->Y_cbother*F->E_nb1->M_An->nrows; n1b = D->Y_cbother*F->E_nb1->M_An->rpop_b;
    n2a = D->Y_cbother*F->E_nb2->M_An->nrows; n2b = D->Y_cbother*minimum(rdrop2,F->E_nb2->M_An->rpop_b); 
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_ynut),M_flag,&(F->M_ynut_trm),&(F->M_ynut),&(F_in->lf_ynut),&(F_in->M_ynut));
    sprintf(prefix,"F->M_unyt");
    n1a = D->Y_cbother*F->E_nb1->M_An->nrows; n1b = D->Y_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->Y_cbother*F->E_nb2->M_An->nrows; n2b = D->Y_cbother*F->E_nb2->M_An->rpop_b; 
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_unyt),M_flag,&(F->M_unyt_trm),&(F->M_unyt),&(F_in->lf_unyt),&(F_in->M_unyt));
    sprintf(prefix,"F->M_unut");
    n1a = D->Y_cbother*F->E_nb1->M_An->nrows; n1b = D->Y_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->Y_cbother*F->E_nb2->M_An->nrows; n2b = D->Y_cbother*minimum(rdrop2,F->E_nb2->M_An->rpop_b); 
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_unut),M_flag,&(F->M_unut_trm),&(F->M_unut),&(F_in->lf_unut),&(F_in->M_unut));
    sprintf(prefix,"F->M_ankt");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*F->E_nb1->M_An->rpop_b;
    n2a = D->A_cbother*F->E_nb2->M_An->nrows; n2b = D->A_cbother*minimum(rdrop2,F->E_nb2->M_An->rpop_b); 
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_ankt),M_flag,&(F->M_ankt_trm),&(F->M_ankt),&(F_in->lf_ankt),&(F_in->M_ankt));
    sprintf(prefix,"F->M_knat");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->A_cbother*F->E_nb2->M_An->nrows; n2b = D->A_cbother*F->E_nb2->M_An->rpop_b; 
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_knat),M_flag,&(F->M_knat_trm),&(F->M_knat),&(F_in->lf_knat),&(F_in->M_knat));
    sprintf(prefix,"F->M_knkt");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->A_cbother*F->E_nb2->M_An->nrows; n2b = D->A_cbother*minimum(rdrop2,F->E_nb2->M_An->rpop_b); 
    bcc_double_copy_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_knkt),M_flag,&(F->M_knkt_trm),&(F->M_knkt),&(F_in->lf_knkt),&(F_in->M_knkt));
    /* if strategy */}
  if (strstr(D->QC_strategy,"YnWt")){
    L_handle_copy(F->lf_kt_r_unwt_s_zn,F_in->lf_kt_r_unwt_s_zn);
    L_handle_copy(F->lf_kt_r_knzt_s_zn,F_in->lf_kt_r_knzt_s_zn);
    L_handle_copy(F->lf_kt_r_hnvt_s_zn,F_in->lf_kt_r_hnvt_s_zn);
    L_handle_copy(F->lf_at_t_ynut_r_kn,F_in->lf_at_t_ynut_r_kn);
    L_handle_copy(F->lf_kt_r_unyt_t_an,F_in->lf_kt_r_unyt_t_an);
    L_handle_copy(F->lf_kt_r_unut_r_kn,F_in->lf_kt_r_unut_r_kn);
    L_handle_copy(F->lf_at_t_ankt_r_kn,F_in->lf_at_t_ankt_r_kn);
    L_handle_copy(F->lf_kt_r_knat_t_an,F_in->lf_kt_r_knat_t_an);
    L_handle_copy(F->lf_kt_r_knkt_r_kn,F_in->lf_kt_r_knkt_r_kn);
    /* if strategy */}
  if (strstr(D->QC_strategy,"YnWt")){
    L_handle_copy(F->lf_attjn____vtszn,F_in->lf_attjn____vtszn);
    L_handle_copy(F->lf_attjn____jttan,F_in->lf_attjn____jttan);
    L_handle_copy(F->lf_attjn____htrkn,F_in->lf_attjn____htrkn);
    L_handle_copy(F->lf_ktrhn____jttan,F_in->lf_ktrhn____jttan);
    L_handle_copy(F->lf_ktrhn____htrkn,F_in->lf_ktrhn____htrkn);
    /* if strategy */}
  if (verbose>1){ printf(" %% copying lf_ \n");}
  L_handle_copy(F->lf_AnZtSWnYt,F_in->lf_AnZtSWnYt);
  L_handle_copy(F->lf_AnZtSZnAt,F_in->lf_AnZtSZnAt);
  L_handle_copy(F->lf_AnAtTYnYt,F_in->lf_AnAtTYnYt);
  L_handle_copy(F->lf_AnAtTAnAt,F_in->lf_AnAtTAnAt);
  L_handle_copy(F->lf_AnZt,F_in->lf_AnZt);
  L_handle_copy(F->lf_AnZt_S_ZnAt,F_in->lf_AnZt_S_ZnAt);
  L_handle_copy(F->lf_AnAt_T_YnYt,F_in->lf_AnAt_T_YnYt);
  L_handle_copy(F->lf_AnAt_T_AnAt,F_in->lf_AnAt_T_AnAt);
  L_handle_copy(F->lf_An_ZtSWn_Yt,F_in->lf_An_ZtSWn_Yt);
  L_handle_copy(F->lf_An_ZtSZn_At,F_in->lf_An_ZtSZn_At);
  L_handle_copy(F->lf_An_AtTYn_Yt,F_in->lf_An_AtTYn_Yt);
  L_handle_copy(F->lf_An_AtTAn_At,F_in->lf_An_AtTAn_At);
  L_handle_copy(F->lf_AtTYnWtSZn,F_in->lf_AtTYnWtSZn);
  L_handle_copy(F->lf_AtTAnZtSZn,F_in->lf_AtTAnZtSZn);
  L_handle_copy(F->lf_AtTYnYtTAn,F_in->lf_AtTYnYtTAn);
  L_handle_copy(F->lf_AtTAnAtTAn,F_in->lf_AtTAnAtTAn);
  L_handle_copy(F->lf_AtTYn____WtSZn,F_in->lf_AtTYn____WtSZn);
  L_handle_copy(F->lf_AtTAn____ZtSZn,F_in->lf_AtTAn____ZtSZn);
  L_handle_copy(F->lf_AtTYn____YtTAn,F_in->lf_AtTYn____YtTAn);
  L_handle_copy(F->lf_AtTAn____AtTAn,F_in->lf_AtTAn____AtTAn);
  L_handle_copy(F->lf_At_T_YnWt_S_Zn,F_in->lf_At_T_YnWt_S_Zn);
  L_handle_copy(F->lf_At_T_AnZt_S_Zn,F_in->lf_At_T_AnZt_S_Zn);
  L_handle_copy(F->lf_At_T_YnYt_T_An,F_in->lf_At_T_YnYt_T_An);
  L_handle_copy(F->lf_At_T_AnAt_T_An,F_in->lf_At_T_AnAt_T_An);
  if (verbose>1){ printf(" %% copying doublestudy\n");}
  L_handle_copy(F->lf_An_a2d2_Zt_Sn,F_in->lf_An_a2d2_Zt_Sn);
  L_handle_copy(F->lf_An_a2d2_At_Tn,F_in->lf_An_a2d2_At_Tn);
  if (verbose>1){ printf(" %% copying F->QX\n");}
  L_handle_copy(F->QR_AnZtSWnYt,F_in->QR_AnZtSWnYt);
  L_handle_copy(F->QR_AnZtSZnAt,F_in->QR_AnZtSZnAt);
  L_handle_copy(F->QR_AnAtTYnYt,F_in->QR_AnAtTYnYt);
  L_handle_copy(F->QR_AnAtTAnAt,F_in->QR_AnAtTAnAt);
  L_handle_copy(F->QC_AtTYnWtSZn,F_in->QC_AtTYnWtSZn);
  L_handle_copy(F->QC_AtTAnZtSZn,F_in->QC_AtTAnZtSZn);
  L_handle_copy(F->QC_AtTYnYtTAn,F_in->QC_AtTYnYtTAn);
  L_handle_copy(F->QC_AtTAnAtTAn,F_in->QC_AtTAnAtTAn);
  if (verbose>1){ printf(" %% copying F->QX_uu\n");}
  L_handle_copy(F->QR_AnZtSWnYt_uu,F_in->QR_AnZtSWnYt_uu);
  L_handle_copy(F->QR_AnZtSZnAt_uu,F_in->QR_AnZtSZnAt_uu);
  L_handle_copy(F->QR_AnAtTYnYt_uu,F_in->QR_AnAtTYnYt_uu);
  L_handle_copy(F->QR_AnAtTAnAt_uu,F_in->QR_AnAtTAnAt_uu);
  L_handle_copy(F->QC_AtTYnWtSZn_uu,F_in->QC_AtTYnWtSZn_uu);
  L_handle_copy(F->QC_AtTAnZtSZn_uu,F_in->QC_AtTAnZtSZn_uu);
  L_handle_copy(F->QC_AtTYnYtTAn_uu,F_in->QC_AtTYnYtTAn_uu);
  L_handle_copy(F->QC_AtTAnAtTAn_uu,F_in->QC_AtTAnAtTAn_uu);
  GLOBAL_toc(0,verbose," %% copying data_structures: ");
}

void bcc_double_init_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,int *trm_flag_p,struct M_handle **M_p)
{
  int verbose=0;
  int length=0;
  length  = n1b*n2b; *L_p = L_handle_make(length); if (!(*L_p)){ printf(" %% Warning! not enough memory in bcc_double_init_lf\n");}
  if (M_flag){
    if (0){}
    else if (n2b<=n1b){ *trm_flag_p = 0; *M_p = M_handle_w_make(BITJ,GLOBAL_B_MLT,n1a,n2b);}
    else if (n1b<=n2b){ *trm_flag_p = 1; *M_p = M_handle_w_make(BITJ,GLOBAL_B_MLT,n2a,n1b);}
    /* if (M_flag){ } */}
  if (verbose){ printf(" %% %s n1a %d n1b %d n2a %d n2b %d --> length %d trm_flag %d \n",prefix,n1a,n1b,n2a,n2b,length,*trm_flag_p);}
}

void bcc_double_init_lf(struct bcc_double *F)
{
  int verbose=0;
  struct bcc_ajdk *D=F->D;
  int n1a=0,n2a=0,n1b,n2b=0;
  int M_flag=0,rdrop=0,rdrop2=0;
  char prefix[FNAMESIZE];
  GLOBAL_tic(0);
  if (verbose>1){ printf(" %% bcc_double_init_lf: %d,%d\n",F->nb1,F->nb2);}
  if (verbose>1){ printf(" %% initializing F->lf_YnWt, etc\n");}
  F->lf_YnWt = NULL; F->lf_AnZt = NULL; F->lf_YnYt = NULL; F->lf_AnAt = NULL; F->M_YnWt = NULL; F->M_AnZt = NULL; F->M_YnYt = NULL; F->M_AnAt = NULL;
  if (strstr(D->QR_strategy,"YnWt") || strstr(D->QC_strategy,"YnWt")){
    M_flag = (strstr(D->QC_strategy,"YnWt") ? 1 : 0);
    n1a = F->E_nb1->M_Yn->nrows; n1b = F->E_nb1->M_Yn->rpop_b; n2a = F->E_nb2->M_Wn->nrows; n2b = F->E_nb2->M_Wn->rpop_b; sprintf(prefix,"F->M_YnWt");
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_YnWt),M_flag,&(F->M_YnWt_trm),&(F->M_YnWt));
    n1a = F->E_nb1->M_An->nrows; n1b = F->E_nb1->M_An->rpop_b; n2a = F->E_nb2->M_Zn->nrows; n2b = F->E_nb2->M_Zn->rpop_b; sprintf(prefix,"F->M_AnZt");
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_AnZt),M_flag,&(F->M_AnZt_trm),&(F->M_AnZt));
    n1a = F->E_nb1->M_Yn->nrows; n1b = F->E_nb1->M_Yn->rpop_b; n2a = F->E_nb2->M_Yn->nrows; n2b = F->E_nb2->M_Yn->rpop_b; sprintf(prefix,"F->M_YnYt");
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_YnYt),M_flag,&(F->M_YnYt_trm),&(F->M_YnYt));
    n1a = F->E_nb1->M_An->nrows; n1b = F->E_nb1->M_An->rpop_b; n2a = F->E_nb2->M_An->nrows; n2b = F->E_nb2->M_An->rpop_b; sprintf(prefix,"F->M_AnAt");
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_AnAt),M_flag,&(F->M_AnAt_trm),&(F->M_AnAt));
    /* if strategy */}
  get_xdrop(D->A_rpop_j_total,D->A_cpop_j,&rdrop,NULL); rdrop2 = maximum(2,rdrop);
  /* QR YnWt lrup */
  F->lf_jnvt=NULL; F->lf_jnjt=NULL;
  if (strstr(D->QR_strategy,"YnWt")){
    n1a = F->E_nb1->M_An->nrows; n1b = F->E_nb1->M_An->rpop_b; 
    n2a = F->E_nb2->M_Zn->nrows; n2b = F->E_nb2->M_Zn->rpop_b; 
    F->length = n1b*n2b; F->lf_jnvt = L_handle_make(F->length); if (!F->lf_jnvt){ printf(" %% Warning! not enough memory in bcc_double_init_lf\n");}
    n1a = F->E_nb1->M_An->nrows; n1b = F->E_nb1->M_An->rpop_b; 
    n2a = F->E_nb2->M_An->nrows; n2b = F->E_nb2->M_An->rpop_b; 
    F->length = n1b*n2b; F->lf_jnjt = L_handle_make(F->length); if (!F->lf_jnjt){ printf(" %% Warning! not enough memory in bcc_double_init_lf\n");}
    /* if strategy */}
  /* QC YnWt lrup */
  F->lf_unwt=NULL; F->lf_knzt=NULL; F->lf_hnvt=NULL; F->lf_ynut=NULL; F->lf_unyt=NULL; F->lf_unut=NULL; F->lf_ankt=NULL; F->lf_knat=NULL; F->lf_knkt=NULL;
  F->M_unwt=NULL; F->M_knzt=NULL; F->M_hnvt=NULL; F->M_ynut=NULL; F->M_unyt=NULL; F->M_unut=NULL; F->M_ankt=NULL; F->M_knat=NULL; F->M_knkt=NULL;
  if (strstr(D->QC_strategy,"YnWt")){
    M_flag = 1;
    sprintf(prefix,"F->M_unwt");
    n1a = D->Y_cbother*F->E_nb1->M_An->nrows; n1b = D->Y_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->Y_cbother*F->E_nb2->M_Zn->nrows; n2b = D->Y_cbother*F->E_nb2->M_Zn->rpop_b; 
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_unwt),M_flag,&(F->M_unwt_trm),&(F->M_unwt));
    sprintf(prefix,"F->M_knzt");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->A_cbother*F->E_nb2->M_Zn->nrows; n2b = D->A_cbother*F->E_nb2->M_Zn->rpop_b; 
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_knzt),M_flag,&(F->M_knzt_trm),&(F->M_knzt));
    sprintf(prefix,"F->M_hnvt");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->A_cbother*F->E_nb2->M_Zn->nrows; n2b = D->A_cbother*F->E_nb2->M_Zn->rpop_b; 
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_hnvt),M_flag,&(F->M_hnvt_trm),&(F->M_hnvt));
    sprintf(prefix,"F->M_ynut");
    n1a = D->Y_cbother*F->E_nb1->M_An->nrows; n1b = D->Y_cbother*F->E_nb1->M_An->rpop_b;
    n2a = D->Y_cbother*F->E_nb2->M_An->nrows; n2b = D->Y_cbother*minimum(rdrop2,F->E_nb2->M_An->rpop_b); 
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_ynut),M_flag,&(F->M_ynut_trm),&(F->M_ynut));
    sprintf(prefix,"F->M_unyt");
    n1a = D->Y_cbother*F->E_nb1->M_An->nrows; n1b = D->Y_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->Y_cbother*F->E_nb2->M_An->nrows; n2b = D->Y_cbother*F->E_nb2->M_An->rpop_b; 
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_unyt),M_flag,&(F->M_unyt_trm),&(F->M_unyt));
    sprintf(prefix,"F->M_unut");
    n1a = D->Y_cbother*F->E_nb1->M_An->nrows; n1b = D->Y_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->Y_cbother*F->E_nb2->M_An->nrows; n2b = D->Y_cbother*minimum(rdrop2,F->E_nb2->M_An->rpop_b); 
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_unut),M_flag,&(F->M_unut_trm),&(F->M_unut));
    sprintf(prefix,"F->M_ankt");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*F->E_nb1->M_An->rpop_b;
    n2a = D->A_cbother*F->E_nb2->M_An->nrows; n2b = D->A_cbother*minimum(rdrop2,F->E_nb2->M_An->rpop_b); 
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_ankt),M_flag,&(F->M_ankt_trm),&(F->M_ankt));
    sprintf(prefix,"F->M_knat");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->A_cbother*F->E_nb2->M_An->nrows; n2b = D->A_cbother*F->E_nb2->M_An->rpop_b; 
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_knat),M_flag,&(F->M_knat_trm),&(F->M_knat));
    sprintf(prefix,"F->M_knkt");
    n1a = D->A_cbother*F->E_nb1->M_An->nrows; n1b = D->A_cbother*minimum(rdrop2,F->E_nb1->M_An->rpop_b);
    n2a = D->A_cbother*F->E_nb2->M_An->nrows; n2b = D->A_cbother*minimum(rdrop2,F->E_nb2->M_An->rpop_b); 
    bcc_double_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(F->lf_knkt),M_flag,&(F->M_knkt_trm),&(F->M_knkt));
    /* if strategy */}
  F->lf_kt_r_unwt_s_zn=NULL; F->lf_kt_r_knzt_s_zn=NULL; F->lf_kt_r_hnvt_s_zn=NULL; F->lf_at_t_ynut_r_kn=NULL; F->lf_kt_r_unyt_t_an=NULL; F->lf_kt_r_unut_r_kn=NULL; F->lf_at_t_ankt_r_kn=NULL; F->lf_kt_r_knat_t_an=NULL; F->lf_kt_r_knkt_r_kn=NULL;
  if (strstr(D->QC_strategy,"YnWt")){
    F->length = D->A_ncols*D->T_ncols; F->lf_kt_r_unwt_s_zn = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_kt_r_knzt_s_zn = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_kt_r_hnvt_s_zn = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_at_t_ynut_r_kn = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_kt_r_unyt_t_an = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_kt_r_unut_r_kn = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_at_t_ankt_r_kn = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_kt_r_knat_t_an = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_kt_r_knkt_r_kn = L_handle_make(F->length);
    /* if strategy */}
  F->lf_attjn____vtszn=NULL; F->lf_attjn____jttan=NULL; F->lf_attjn____htrkn=NULL; F->lf_ktrhn____jttan=NULL; F->lf_ktrhn____htrkn=NULL;
  if (strstr(D->QC_strategy,"YnWt")){
    F->length = D->A_ncols*D->T_ncols; F->lf_attjn____vtszn = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_attjn____jttan = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_attjn____htrkn = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_ktrhn____jttan = L_handle_make(F->length);
    F->length = D->A_ncols*D->T_ncols; F->lf_ktrhn____htrkn = L_handle_make(F->length);
    /* if strategy */}
  if (verbose>1){ printf(" %% initializing lf_ \n");}
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_AnZtSWnYt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_AnZtSZnAt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_AnAtTYnYt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_AnAtTAnAt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_AnZt_S_WnYt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_AnZt_S_ZnAt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_AnAt_T_YnYt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_AnAt_T_AnAt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_An_ZtSWn_Yt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_An_ZtSZn_At = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_An_AtTYn_Yt = L_handle_make(F->length);
  F->length = F->E_nb1->M_An->nrows*D->T_ncols; F->lf_An_AtTAn_At = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_AtTYnWtSZn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_AtTAnZtSZn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_AtTYnYtTAn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_AtTAnAtTAn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_AtTYn____WtSZn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_AtTAn____ZtSZn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_AtTYn____YtTAn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_AtTAn____AtTAn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_At_T_YnWt_S_Zn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_At_T_AnZt_S_Zn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_At_T_YnYt_T_An = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->lf_At_T_AnAt_T_An = L_handle_make(F->length); if (!F->lf_At_T_AnAt_T_An){ printf(" %% Warning! not enough memory in bcc_double_init_lf\n");}
  if (verbose>1){ printf(" %% initializing doublestudy\n");}
  F->length = F->E_nb1->A_nrows*F->D->T_ncols; F->lf_An_a2d2_Zt_Sn = L_handle_make(F->length); if (!F->lf_An_a2d2_Zt_Sn){ printf(" %% Warning! not enough memory in bcc_double_init_lf\n");}
  F->length = F->E_nb1->A_nrows*F->D->T_ncols; F->lf_An_a2d2_At_Tn = L_handle_make(F->length); if (!F->lf_An_a2d2_At_Tn){ printf(" %% Warning! not enough memory in bcc_double_init_lf\n");}
  if (verbose>1){ printf(" %% initializing F->QX\n");}
  F->QR_AnZtSWnYt=NULL; F->QR_AnZtSZnAt=NULL; F->QR_AnAtTYnYt=NULL; F->QR_AnAtTAnAt=NULL; F->QC_AtTYnWtSZn=NULL; F->QC_AtTAnZtSZn=NULL; F->QC_AtTYnYtTAn=NULL; F->QC_AtTAnAtTAn=NULL;
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_AnZtSWnYt = L_handle_make(F->length);
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_AnZtSZnAt = L_handle_make(F->length);
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_AnAtTYnYt = L_handle_make(F->length);
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_AnAtTAnAt = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_AtTYnWtSZn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_AtTAnZtSZn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_AtTYnYtTAn = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_AtTAnAtTAn = L_handle_make(F->length); if (!F->QC_AtTAnAtTAn){ printf(" %% Warning! not enough memory in bcc_double_init_lf\n");}
  if (verbose>1){ printf(" %% initializing F->QX_uu\n");}
  F->QR_AnZtSWnYt_uu=NULL; F->QR_AnZtSZnAt_uu=NULL; F->QR_AnAtTYnYt_uu=NULL; F->QR_AnAtTAnAt_uu=NULL; F->QC_AtTYnWtSZn_uu=NULL; F->QC_AtTAnZtSZn_uu=NULL; F->QC_AtTYnYtTAn_uu=NULL; F->QC_AtTAnAtTAn_uu=NULL;
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_AnZtSWnYt_uu = L_handle_make(F->length);
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_AnZtSZnAt_uu = L_handle_make(F->length);
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_AnAtTYnYt_uu = L_handle_make(F->length);
  F->length = F->E_nb1->A_nrows*D->T_ncols; F->QR_AnAtTAnAt_uu = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_AtTYnWtSZn_uu = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_AtTAnZtSZn_uu = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_AtTYnYtTAn_uu = L_handle_make(F->length);
  F->length = D->A_ncols*D->T_ncols; F->QC_AtTAnAtTAn_uu = L_handle_make(F->length); if (!F->QC_AtTAnAtTAn_uu){ printf(" %% Warning! not enough memory in bcc_double_init_lf\n");}
  GLOBAL_toc(0,verbose," %% initializing data_structures: ");
}

void bcc_copy_A_p(struct bcc_ajdk *D,struct bcc_ajdk *D_in)
{
  /* copies A_p and Y_p for each column from D_in */
  int verbose=0; char tempchar[FNAMESIZE];
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_single **E_in_ = D_in->E_;
  int nb=0,nc=0,nc_p=0,nc_start=0,nc_final=0;
  int A_pcols=0,Y_pcols=0,A_ncols=0,Y_ncols=0,A_rpop_b_total=0,Z_rpop_b_total=0;
  double tmp_sum=0,tmp_num=0;
  if (verbose){ printf(" %% [entering bcc_copy_A_p]\n"); wkspace_printf();}
  if (verbose){ printf(" %% \n");}
  D->A_pcols = psize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/POPLENGTH */; D->Y_pcols = psize(D->Y_ncols)/* rup(D->Y_ncols+D->Y_ncols_extend,POPLENGTH)/POPLENGTH */;
  if (verbose){ printf(" %% nbins %d; A_ncols %d A_pcols %d, Y_ncols %d Y_pcols %d\n",nbins,D->A_ncols,D->A_pcols,D->Y_ncols,D->Y_pcols);}
  memcpy(D->AZ_rsum,D_in->AZ_rsum,D->A_ncols*sizeof(double)); 
  memcpy(D->A_p,D_in->A_p,D->A_pcols*sizeof(double)); 
  memcpy(D->A_ajdk,D_in->A_ajdk,AJDK_TOT*D->A_pcols*sizeof(double)); 
  memcpy(D->YW_rsum,D_in->YW_rsum,D->Y_ncols*sizeof(double)); 
  memcpy(D->Y_p,D_in->Y_p,D->Y_pcols*sizeof(double)); 
  memcpy(D->Y_ajdk,D_in->Y_ajdk,AJDK_TOT*D->Y_pcols*sizeof(double)); 
  for (nb=0;nb<nbins;nb++){
    L_handle_copy(E_[nb]->lf_At_rsum,E_in_[nb]->lf_At_rsum);
    L_handle_copy(E_[nb]->lf_Zt_rsum,E_in_[nb]->lf_Zt_rsum);
    L_handle_copy(E_[nb]->lf_Yt_rsum,E_in_[nb]->lf_Yt_rsum);
    L_handle_copy(E_[nb]->lf_Wt_rsum,E_in_[nb]->lf_Wt_rsum);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% [finished bcc_copy_A_p]\n");  wkspace_printf();}
}

void bcc_X_nrows_total(struct bcc_ajdk *D)
{
  /* initializes D->A_nrows_total and D->Z_nrows_total for each column using row-mask mr_b, not mr_j */
  int verbose=0;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0,nc=0;
  if (verbose){ printf(" %% [entering bcc_X_nrows_total]\n");}
  D->A_nrows_total=0; D->A_rpop_b_total=0; D->A_rpop_j_total=0;
  D->Z_nrows_total=0; D->Z_rpop_b_total=0; D->Z_rpop_j_total=0;
  for (nb=0;nb<nbins;nb++){
    D->A_nrows_total += E_[nb]->A_nrows; D->A_rpop_b_total += E_[nb]->A_rpop_b; D->A_rpop_j_total += E_[nb]->A_rpop_j;
    D->Z_nrows_total += E_[nb]->Z_nrows; D->Z_rpop_b_total += E_[nb]->Z_rpop_b; D->Z_rpop_j_total += E_[nb]->Z_rpop_j;
    /* for (nb=0;nb<nbins;nb++){ }*/}
  if (verbose){ printf(" %% D->A_nrows_total %d D->Z_nrows_total %d\n",D->A_nrows_total,D->Z_nrows_total);}
  if (verbose){ printf(" %% D->A_rpop_b_total %d D->Z_rpop_b_total %d\n",D->A_rpop_b_total,D->Z_rpop_b_total);}
  if (verbose){ printf(" %% D->A_rpop_j_total %d D->Z_rpop_j_total %d\n",D->A_rpop_j_total,D->Z_rpop_j_total);}
  if (verbose){ printf(" %% [finished bcc_X_nrows_total]\n");}
}

void bcc_load_A_p(struct bcc_ajdk *D)
{
  /* calculates A_p and Y_p for each column using row-mask mr_b, not mr_j */
  int verbose=0; char tempchar[FNAMESIZE];
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0,nc=0,nc_p=0,nc_start=0,nc_final=0;
  int A_pcols=0,Y_pcols=0,A_ncols=0,Y_ncols=0,A_rpop_b_total=0,Z_rpop_b_total=0;
  double *AZ_rsum=NULL,*A_p=NULL,**A_ajdk_p=NULL,*YW_rsum=NULL,*Y_p=NULL,**Y_ajdk_p=NULL;
  unsigned char *A_bmc_b=NULL,*Y_bmc_b=NULL;
  double tmp_sum=0,tmp_num=0;
  if (verbose){ printf(" %% [entering bcc_load_A_p]\n"); wkspace_printf();}
  if (verbose){ printf(" %% \n");}
  D->A_pcols = psize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/POPLENGTH */; D->Y_pcols = psize(D->Y_ncols)/* rup(D->Y_ncols+D->Y_ncols_extend,POPLENGTH)/POPLENGTH */;
  if (verbose){ printf(" %% nbins %d; A_ncols %d A_pcols %d, Y_ncols %d Y_pcols %d\n",nbins,D->A_ncols,D->A_pcols,D->Y_ncols,D->Y_pcols);}
  D->AZ_rsum = (double *) wkspace_all0c(D->A_ncols*sizeof(double)); 
  D->A_p = (double *) wkspace_all0c(D->A_pcols*sizeof(double)); 
  D->A_ajdk = (double *) wkspace_all0c(AJDK_TOT*D->A_pcols*sizeof(double)); 
  D->YW_rsum = (double *) wkspace_all0c(D->Y_ncols*sizeof(double)); 
  D->Y_p = (double *) wkspace_all0c(D->Y_pcols*sizeof(double)); 
  D->Y_ajdk = (double *) wkspace_all0c(AJDK_TOT*D->Y_pcols*sizeof(double)); 
  if (verbose){ printf(" %% D->A_ncols %d D->Y_ncols %d D->A_pcols %d D->Y_pcols %d\n",D->A_ncols,D->Y_ncols,D->A_pcols,D->Y_pcols); wkspace_printf();}
  A_ncols = D->A_ncols; Y_ncols = D->Y_ncols; A_pcols = D->A_pcols; Y_pcols = D->Y_pcols; A_rpop_b_total = D->A_rpop_b_total; Z_rpop_b_total = D->Z_rpop_b_total;
  AZ_rsum = D->AZ_rsum; A_p = D->A_p; A_ajdk_p = &(D->A_ajdk); A_bmc_b = D->A_bmc_b; 
  YW_rsum = D->YW_rsum; Y_p = D->Y_p; Y_ajdk_p = &(D->Y_ajdk); Y_bmc_b = D->Y_bmc_b;
  if (verbose){ printf(" %% calculating E_[nb]->lf_At_rsum etc.\n"); wkspace_printf();}
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){
    for (nc=0;nc<D->A_ncols;nc++){ E_[nb]->lf_At_rsum->lf[nc]=0;}
    GLOBAL_pthread_tic(); wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),TYPE_00,SPACING_a,E_[nb]->M_At,&(E_[nb]->lf_At_rsum)); GLOBAL_pthread_toc();
    for (nc=0;nc<D->A_ncols;nc++){ E_[nb]->lf_Zt_rsum->lf[nc]=0;}
    GLOBAL_pthread_tic(); wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),TYPE_00,SPACING_a,E_[nb]->M_Zt,&(E_[nb]->lf_Zt_rsum)); GLOBAL_pthread_toc();
    for (nc=0;nc<D->Y_ncols;nc++){ E_[nb]->lf_Yt_rsum->lf[nc]=0;}
    GLOBAL_pthread_tic(); wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),TYPE_00,SPACING_a,E_[nb]->M_Yt,&(E_[nb]->lf_Yt_rsum)); GLOBAL_pthread_toc();
    for (nc=0;nc<D->Y_ncols;nc++){ E_[nb]->lf_Wt_rsum->lf[nc]=0;}
    GLOBAL_pthread_tic(); wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),TYPE_00,SPACING_a,E_[nb]->M_Wt,&(E_[nb]->lf_Wt_rsum)); GLOBAL_pthread_toc();
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  if (verbose){ printf(" %% finished calculating E_[nb]->lf_At_rsum etc.\n"); wkspace_printf();}
  if (verbose){
    for (nb=0;nb<nbins;nb++){
      sprintf(D->tmpAnchar," %%%% E_[nb]->lf_At_[%d]_rsum->lf: ",nb); lfprintf(E_[nb]->lf_At_rsum,D->tmpAnchar);
      sprintf(D->tmpZnchar," %%%% E_[nb]->lf_Zt_[%d]_rsum->lf: ",nb); lfprintf(E_[nb]->lf_Zt_rsum,D->tmpZnchar);
      sprintf(D->tmpYnchar," %%%% E_[nb]->lf_Yt_[%d]_rsum->lf: ",nb); lfprintf(E_[nb]->lf_Yt_rsum,D->tmpYnchar);
      sprintf(D->tmpWnchar," %%%% E_[nb]->lf_Wt_[%d]_rsum->lf: ",nb); lfprintf(E_[nb]->lf_Wt_rsum,D->tmpWnchar);
      /* for (nb=0;nb<nbins;nb++){ } */}
    wkspace_printf(); /* if (verbose){ } */}
  if (verbose){ printf(" %% calculating AZ_rsum\n");}
  for (nc=0;nc<A_ncols;nc++){ AZ_rsum[nc]=0;}
  for (nb=0;nb<nbins;nb++){ 
    dra_plusequals(AZ_rsum,D->A_ncols,E_[nb]->lf_At_rsum->lf);
    dra_plusequals(AZ_rsum,D->A_ncols,E_[nb]->lf_Zt_rsum->lf);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ raprintf(AZ_rsum,"double",1,A_ncols," %% AZ_rsum: "); wkspace_printf();}
  if (verbose){ printf(" %% calculating YW_rsum\n");}
  for (nc=0;nc<Y_ncols;nc++){ YW_rsum[nc]=0;}
  for (nb=0;nb<nbins;nb++){ 
    dra_plusequals(YW_rsum,D->Y_ncols,E_[nb]->lf_Yt_rsum->lf);
    dra_plusequals(YW_rsum,D->Y_ncols,E_[nb]->lf_Wt_rsum->lf);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ raprintf(YW_rsum,"double",1,Y_ncols," %% YW_rsum: "); wkspace_printf();}
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
  if (verbose){ printf(" %% calculating Y_p\n");}
  for (nc_p=0;nc_p<Y_pcols;nc_p++){ Y_p[nc_p]=0.5;}
  if (GLOBAL_TEST_sparse){
    for (nc_p=0;nc_p<Y_pcols;nc_p++){ 
      nc_start = minimum(Y_ncols-1,maximum(0,0 + (nc_p+0)*POPLENGTH - 0));
      nc_final = minimum(Y_ncols-1,maximum(0,0 + (nc_p+1)*POPLENGTH - 1));
      if (verbose){ printf(" %% nc_p %d: [%d..%d]:",nc_p,nc_start,nc_final);}
      tmp_sum = 0; tmp_num=0;
      for (nc=nc_start;nc<=nc_final;nc++){
	if (bget__on(Y_bmc_b,nc)){ tmp_sum += YW_rsum[nc]; tmp_num+=1;}
	/* for (nc=nc_start;nc<=nc_final;nc++){ } */}
      if (verbose){ printf(" tmp_sum %f tmp_num %f*(%d+%d);",tmp_sum,tmp_num,A_rpop_b_total,Z_rpop_b_total);}
      if (tmp_num==0){ Y_p[nc_p] = 0.5;}
      else /* if (tmp_num>0) */{ Y_p[nc_p] = crop((tmp_sum)/(tmp_num*(A_rpop_b_total+Z_rpop_b_total)),0.01,0.99);}
      if (verbose){ printf(" Y_p %f\n",Y_p[nc_p]);}
      /*  for (nc_p=0;nc_p<Y_pcols;nc_p++){ } */}
    /* if (GLOBAL_TEST_sparse){ } */}
  if (verbose){ raprintf(Y_p,"double",1,Y_pcols," %% Y_p: "); wkspace_printf();}
  calc_A_ajdk(Y_p,Y_pcols,Y_ajdk_p); if (verbose){ raprintf(*Y_ajdk_p,"double",Y_pcols,AJDK_TOT," %% Y_ajdk: "); wkspace_printf();}
  sprintf(tempchar,"%s/AZ_rsum.mda",GLOBAL_DIR_NAME); mda_write_r8(tempchar,A_ncols,1,AZ_rsum);
  sprintf(tempchar,"%s/YW_rsum.mda",GLOBAL_DIR_NAME); mda_write_r8(tempchar,Y_ncols,1,YW_rsum);
  sprintf(tempchar,"%s/A_p.mda",GLOBAL_DIR_NAME); mda_write_r8(tempchar,A_pcols,1,A_p);
  sprintf(tempchar,"%s/Y_p.mda",GLOBAL_DIR_NAME); mda_write_r8(tempchar,Y_pcols,1,Y_p);
  if (verbose){ printf(" %% [finished bcc_load_A_p]\n");  wkspace_printf();}
}

void bcc_init_A_p(double mrnd,struct bcc_ajdk *D)
{
  /* initializes A_p and Y_p for each column using row-mask mr_b, not mr_j */
  int verbose=0;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0,nc=0;
  if (verbose){ printf(" %% [entering bcc_init_A_p]\n");}
  D->A_pcols = psize(D->A_ncols); 
  D->A_p = (double *)wkspace_all0c(D->A_pcols*sizeof(double));
  for (nc=0;nc<D->A_pcols;nc++){ D->A_p[nc] = maximum(0.1,minimum(0.9,rand01));} if (GLOBAL_TEST_sparse==0){ printf(" %% turning off D and a\n"); for (nc=0;nc<D->A_pcols;nc++){ D->A_p[nc] = 0.5;}}
  D->A_ajdk = (double *)wkspace_all0c(AJDK_TOT*D->A_pcols*sizeof(double)); 
  calc_A_ajdk(D->A_p,D->A_pcols,&(D->A_ajdk));
  D->Y_pcols = psize(D->Y_ncols);
  D->Y_p = (double *)wkspace_all0c(D->Y_pcols*sizeof(double));
  for (nc=0;nc<D->Y_pcols;nc++){ D->Y_p[nc] = maximum(0.1,minimum(0.9,rand01));} if (GLOBAL_TEST_sparse==0){ printf(" %% turning off D and a\n"); for (nc=0;nc<D->Y_pcols;nc++){ D->Y_p[nc] = 0.5;}}
  D->Y_ajdk = (double *)wkspace_all0c(AJDK_TOT*D->Y_pcols*sizeof(double));
  calc_A_ajdk(D->Y_p,D->Y_pcols,&(D->Y_ajdk));
  if (verbose){ printf(" %% [finished bcc_init_A_p]\n");}
}

void bcc_copy_QX(struct bcc_ajdk *D,struct bcc_ajdk *D_in)
{
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_single **E_in_ = D_in->E_;
  int nb1=0; struct bcc_single *E=NULL;
  D->out_iteration=0; D->out_trace_length = 6; memcpy(D->out_trace,D_in->out_trace,(D->A_ncols+D->A_nrows_total)*D->out_trace_length*sizeof(double));
  memcpy(D->out_xdrop_a,D_in->out_xdrop_a,(D->A_ncols+D->A_nrows_total)*2*sizeof(int));
  memcpy(D->out_xdrop_b,D_in->out_xdrop_b,(D->A_ncols+D->A_nrows_total)*2*sizeof(int));
  memcpy(D->QC_AtTYnWtSZn_nrm,D_in->QC_AtTYnWtSZn_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_AtTAnZtSZn_nrm,D_in->QC_AtTAnZtSZn_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_AtTYnYtTAn_nrm,D_in->QC_AtTYnYtTAn_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_AtTAnAtTAn_nrm,D_in->QC_AtTAnAtTAn_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
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
    memcpy(E->QR_AnZtSWnYt_nrm,E_in_[nb1]->QR_AnZtSWnYt_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_AnZtSZnAt_nrm,E_in_[nb1]->QR_AnZtSZnAt_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_AnAtTYnYt_nrm,E_in_[nb1]->QR_AnAtTYnYt_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_AnAtTAnAt_nrm,E_in_[nb1]->QR_AnAtTAnAt_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    /*  for (nb1=0;nb1<nbins;nb1++){ } */}
}

void bcc_init_QX(struct bcc_ajdk *D)
{
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb1=0; struct bcc_single *E=NULL; int mx_a=0,mx_b=0,mx_j=0,mr_a=0,mr_b=0,mr_j=0;
  D->out_iteration=0; D->out_trace_length = 6; D->out_trace = (double *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*D->out_trace_length*sizeof(double)); if (!D->out_trace){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->out_xdrop_a = (int *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*2*sizeof(int)); if (!D->out_xdrop_a){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->out_xdrop_b = (int *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*2*sizeof(int)); if (!D->out_xdrop_b){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QC_AtTYnWtSZn_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_AtTYnWtSZn_nrm){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");}
  D->QC_AtTAnZtSZn_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_AtTAnZtSZn_nrm){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");}
  D->QC_AtTYnYtTAn_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_AtTYnYtTAn_nrm){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");}
  D->QC_AtTAnAtTAn_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_AtTAnAtTAn_nrm){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");}
  D->QC_sra=(double *)wkspace_all0c(D->A_ncols*sizeof(double)); if (!D->QC_sra){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QC_lmc_a=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_lmc_a){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QC_lmc_b=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_lmc_b){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QC_lmc_j=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_lmc_j){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QR_sra=(double *)wkspace_all0c(D->A_nrows_total*sizeof(double)); if (!D->QR_sra){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QR_imr_a=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_imr_a){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QR_imr_b=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_imr_b){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QR_imr_j=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_imr_j){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QR_lmr_a=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_lmr_a){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QR_lmr_b=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_lmr_b){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QR_lmr_j=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_lmr_j){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  D->QR_lnb=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_lnb){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
  mx_a=0; mx_b=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    E->QR_imr_a = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_imr_a){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
    E->QR_imr_b = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_imr_b){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
    E->QR_imr_j = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_imr_j){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");} 
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
    E->QR_AnZtSWnYt_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_AnZtSWnYt_nrm){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");}
    E->QR_AnZtSZnAt_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_AnZtSZnAt_nrm){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");}
    E->QR_AnAtTYnYt_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_AnAtTYnYt_nrm){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");}
    E->QR_AnAtTAnAt_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_AnAtTAnAt_nrm){ printf(" %% Warning! not enough memory in bcc_init_X_QX\n");}
    /*  for (nb1=0;nb1<nbins;nb1++){ } */}
}

void bcc_M_mxset(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0;
  if (verbose){ printf(" %% setting M_An->mr_j, etc.\n");}
  D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length);
  D->Y_cpop_j = popcount_uchar_array(D->Y_bmc_j,D->Y_mc_length);
  D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length);
  D->A_rpop_j_total=0;D->Z_rpop_j_total=0;
  for (nb=0;nb<nbins;nb++){ 
    E_[nb]->A_rpop_j = popcount_uchar_array(E_[nb]->A_bmr_j,E_[nb]->A_mr_length); D->A_rpop_j_total += E_[nb]->A_rpop_j;
    E_[nb]->Z_rpop_j = popcount_uchar_array(E_[nb]->Z_bmr_j,E_[nb]->Z_mr_length); D->Z_rpop_j_total += E_[nb]->Z_rpop_j;
    M_mxset(E_[nb]->M_An,E_[nb]->A_bmr_j,D->A_bmc_j); M_mxset(E_[nb]->M_At,D->A_bmc_j,E_[nb]->A_bmr_j);
    M_mxset(E_[nb]->M_Zn,E_[nb]->Z_bmr_j,D->A_bmc_j); M_mxset(E_[nb]->M_Zt,D->A_bmc_j,E_[nb]->Z_bmr_j);
    M_mxset(E_[nb]->M_Yn,E_[nb]->A_bmr_j,D->Y_bmc_j); M_mxset(E_[nb]->M_Yt,D->Y_bmc_j,E_[nb]->A_bmr_j);
    M_mxset(E_[nb]->M_Wn,E_[nb]->Z_bmr_j,D->Y_bmc_j); M_mxset(E_[nb]->M_Wt,D->Y_bmc_j,E_[nb]->Z_bmr_j);
    M_mxset(E_[nb]->M_Tn,E_[nb]->A_bmr_j,D->T_bmc_j); M_mxset(E_[nb]->M_Tt,D->T_bmc_j,E_[nb]->A_bmr_j);
    M_mxset(E_[nb]->M_Sn,E_[nb]->Z_bmr_j,D->T_bmc_j); M_mxset(E_[nb]->M_St,D->T_bmc_j,E_[nb]->Z_bmr_j);
    /* for (nb=0;nb<nbins;nb++){ } */}
}

void bcc_copy(struct bcc_ajdk *D,struct bcc_ajdk *D_in)
{
  int verbose=0;
  struct bcc_single **E_ = D->E_; struct bcc_single **E_in_ = D_in->E_;
  struct bcc_double **F_ = D->F_; struct bcc_double **F_in_ = D_in->F_;
  int nbins = GLOBAL_NBINS, nbx=0,nb1=0,nb2=0;
  if (verbose){ printf(" %% [entering bcc_copy] \n");}
  GLOBAL_tic(0);
  if (verbose){ printf(" %% calling bcc_ajdk_copy: \n");}
  bcc_ajdk_copy(D,D_in);
  if (verbose){ printf(" %% finished bcc_ajdk_copy: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ 
    if (verbose){ printf(" %% calling bcc_single_copy_M_An for nb1 %d: \n",nb1);} 
    bcc_single_copy_M_An(E_[nb1],E_in_[nb1]);
    if (verbose){ printf(" %% finished bcc_single_copy_M_An for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% calling bcc_X_nrows_total: \n");}
  bcc_X_nrows_total(D);
  if (verbose){ printf(" %% finished bcc_X_nrows_total: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ 
    if (verbose){ printf(" %% calling bcc_single_copy_lf for nb1 %d: \n",nb1);}
    bcc_single_copy_lf(E_[nb1],E_in_[nb1]);
    if (verbose){ printf(" %% finished bcc_single_copy_lf for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; 
      if (verbose){ printf(" %% calling bcc_double_copy_lf for nb1 %d nb2 %d nbx %d: \n",nb1,nb2,nbx);}
      bcc_double_copy_lf(F_[nbx],F_in_[nbx]); 
      if (verbose){ printf(" %% finished bcc_double_copy_lf for nb1 %d nb2 %d nbx %d: ",nb1,nb2,nbx); wkspace_printf();}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_toc(0,verbose," %% generating data matrices: ");
  if (verbose){ printf(" %% calling bcc_copy_A_p: \n");}
  bcc_copy_A_p(D,D_in);
  if (verbose){ printf(" %% finished bcc_copy_A_p: "); wkspace_printf();}
  if (verbose){ printf(" %% calling bcc_copy_QX\n");}
  bcc_copy_QX(D,D_in);
  if (verbose){ printf(" %% finished bcc_copy_QX: "); wkspace_printf();}
  if (verbose){ printf(" %% calling bcc_M_mxset\n");}
  bcc_M_mxset(D);
  if (verbose){ printf(" %% finished bcc_M_mxset: "); wkspace_printf();}
  if (verbose){ printf(" %% [finished bcc_copy] "); wkspace_printf();}
}

void bcc_load(struct bcc_ajdk **D_p,struct bcc_single ***E_p,struct bcc_double ***F_p,char *QR_strategy,char *QC_strategy)
{
  int verbose=GLOBAL_verbose;
  int nbins = GLOBAL_NBINS, nbx=0,nb1=0,nb2=0;
  struct bcc_ajdk *D=NULL;
  struct bcc_single *E=NULL,*E_nb1=NULL,*E_nb2=NULL;
  struct bcc_double *F=NULL;
  if (verbose){ printf(" %% [entering bcc_load] \n");}
  if (*D_p==NULL){
    (*D_p) = (struct bcc_ajdk *) wkspace_all0c(sizeof(struct bcc_ajdk)); D = (*D_p); D->bitj = BITJ; D->nbins = nbins; D->E_ = NULL; D->F_ = NULL;
    sprintf(D->QR_strategy,"%s",QR_strategy); sprintf(D->QC_strategy,"%s",QC_strategy);
    /* if (*D_p==NULL){ } */}
  if (*E_p==NULL){
    (*E_p) = (struct bcc_single **) wkspace_all0c(nbins*sizeof(struct bcc_single *));
    for (nb1=0;nb1<nbins;nb1++){ (*E_p)[nb1] = (struct bcc_single *) wkspace_all0c(sizeof(struct bcc_single)*1); E = (*E_p)[nb1]; E->D = (*D_p); E->nb=nb1;}
    (*D_p)->E_ = (*E_p);
    /* if (*E_p==NULL){ } */}
  if (*F_p==NULL){
    (*F_p) = (struct bcc_double **) wkspace_all0c(nbins*nbins*sizeof(struct bcc_double *));
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins;
	(*F_p)[nbx] = (struct bcc_double *) wkspace_all0c(sizeof(struct bcc_double)*1); F = (*F_p)[nbx]; F->D = (*D_p); F->nb1=nb1; F->nb2=nb2; F->E_nb1 = (*E_p)[nb1]; F->E_nb2 = (*E_p)[nb2];
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    (*D_p)->F_ = (*F_p);
    /* if (*F_p==NULL){ } */}
  D = (*D_p); 
  GLOBAL_tic(0);
  if (verbose){ printf(" %% calling bcc_ajdk_load: \n");}
  bcc_ajdk_load(D);
  if (verbose){ printf(" %% finished bcc_ajdk_load: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling bcc_single_load_M_An for nb1 %d: \n",nb1);} 
    bcc_single_load_M_An(E);
    if (verbose){ printf(" %% finished bcc_single_load_M_An for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% calling bcc_X_nrows_total: \n");}
  bcc_X_nrows_total(D);
  if (verbose){ printf(" %% finished bcc_X_nrows_total: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling bcc_single_init_lf for nb1 %d: \n",nb1);}
    bcc_single_init_lf(E);
    if (verbose){ printf(" %% finished bcc_single_init_lf for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; E_nb1 = (*E_p)[nb1]; E_nb2 = (*E_p)[nb2]; F = (*F_p)[nbx];
      if (verbose){ printf(" %% calling bcc_double_init_lf for nb1 %d nb2 %d nbx %d: \n",nb1,nb2,nbx);}
      bcc_double_init_lf(F); 
      if (verbose){ printf(" %% finished bcc_double_init_lf for nb1 %d nb2 %d nbx %d: ",nb1,nb2,nbx); wkspace_printf();}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_toc(0,verbose," %% generating data matrices: ");
  if (verbose){ printf(" %% calling bcc_load_A_p: \n");}
  bcc_load_A_p(D);
  if (verbose){ printf(" %% finished bcc_load_A_p: "); wkspace_printf();}
  if (verbose){ printf(" %% calling bcc_init_QX\n");}
  bcc_init_QX(D);
  if (verbose){ printf(" %% finished bcc_init_QX: "); wkspace_printf();}
  if (verbose){ printf(" %% calling bcc_M_mxset\n");}
  bcc_M_mxset(D);
  if (verbose){ printf(" %% finished bcc_M_mxset: "); wkspace_printf();}
  if (verbose){ printf(" %% [finished bcc_load] "); wkspace_printf();}
}

void bcc_init(char *error_vs_speed,double mrnd,struct bcc_ajdk **D_p,struct bcc_single ***E_p,struct bcc_double ***F_p,char *QR_strategy,char *QC_strategy)
{
  int verbose=GLOBAL_verbose;
  int nbins = GLOBAL_NBINS, nbx=0,nb1=0,nb2=0;
  struct bcc_ajdk *D=NULL;
  struct bcc_single *E=NULL,*E_nb1=NULL,*E_nb2=NULL;
  struct bcc_double *F=NULL;
  if (verbose){ printf(" %% [entering bcc_init] \n");}
  if (*D_p==NULL){
    (*D_p) = (struct bcc_ajdk *) wkspace_all0c(sizeof(struct bcc_ajdk)); D = (*D_p); D->bitj = BITJ; D->nbins = nbins; D->E_ = NULL; D->F_ = NULL;
    sprintf(D->QR_strategy,"%s",QR_strategy); sprintf(D->QC_strategy,"%s",QC_strategy);
    /* if (*D_p==NULL){ } */}
  if (*E_p==NULL){
    (*E_p) = (struct bcc_single **) wkspace_all0c(nbins*sizeof(struct bcc_single *));
    for (nb1=0;nb1<nbins;nb1++){ (*E_p)[nb1] = (struct bcc_single *) wkspace_all0c(sizeof(struct bcc_single)*1); E = (*E_p)[nb1]; E->D = (*D_p); E->nb=nb1;}
    (*D_p)->E_ = (*E_p);
    /* if (*E_p==NULL){ } */}
  if (*F_p==NULL){
    (*F_p) = (struct bcc_double **) wkspace_all0c(nbins*nbins*sizeof(struct bcc_double *));
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins;
	(*F_p)[nbx] = (struct bcc_double *) wkspace_all0c(sizeof(struct bcc_double)*1); F = (*F_p)[nbx]; F->D = (*D_p); F->nb1=nb1; F->nb2=nb2; F->E_nb1 = (*E_p)[nb1]; F->E_nb2 = (*E_p)[nb2];
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    (*D_p)->F_ = (*F_p);
    /* if (*F_p==NULL){ } */}
  D = (*D_p);
  GLOBAL_tic(0);
  if (verbose){ printf(" %% calling bcc_ajdk_init: \n");}
  bcc_ajdk_init(mrnd,D);
  if (verbose){ printf(" %% finished bcc_ajdk_init: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling bcc_single_init_M_An for nb1 %d: \n",nb1);}
    bcc_single_init_M_An(error_vs_speed,mrnd,E);
    if (verbose){ printf(" %% finished bcc_single_init_M_An for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% calling bcc_X_nrows_total: \n");}
  bcc_X_nrows_total(D);
  if (verbose){ printf(" %% finished bcc_X_nrows_total: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling bcc_single_init_lf for nb1 %d: \n",nb1);}
    bcc_single_init_lf(E);
    if (verbose){ printf(" %% finished bcc_single_init_lf for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; E_nb1 = (*E_p)[nb1]; E_nb2 = (*E_p)[nb2]; F = (*F_p)[nbx];
      if (verbose){ printf(" %% calling bcc_double_init_lf for nb1 %d nb2 %d nbx %d: \n",nb1,nb2,nbx);}
      bcc_double_init_lf(F);
      if (verbose){ printf(" %% finished bcc_double_init_lf for nb1 %d nb2 %d nbx %d: ",nb1,nb2,nbx); wkspace_printf();}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_toc(0,verbose," %% generating data matrices: ");
  if (verbose){ printf(" %% calling bcc_init_A_p: \n");}
  bcc_init_A_p(mrnd,D);
  if (verbose){ printf(" %% finished bcc_init_A_p: "); wkspace_printf();}
  if (verbose){ printf(" %% calling bcc_init_QX\n");}
  bcc_init_QX(D);
  if (verbose){ printf(" %% finished bcc_init_QX: "); wkspace_printf();}
  if (verbose){ printf(" %% calling bcc_M_mxset\n");}
  bcc_M_mxset(D);
  if (verbose){ printf(" %% finished bcc_M_mxset: "); wkspace_printf();}
  if (verbose){ printf(" %% [finished bcc_init] "); wkspace_printf();}
}

void bcc_init_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  bcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D);
  M_handle_printf(E_[0]->M_An,verbose," %% E_[0]->M_An[0]: "); M_handle_printf(E_[0]->M_At,verbose," %% E_[0]->M_At[0]: ");
  M_handle_printf(E_[1]->M_An,verbose," %% E_[1]->M_An[1]: "); M_handle_printf(E_[1]->M_At,verbose," %% E_[1]->M_At[1]: ");
  M_handle_printf(E_[0]->M_Zn,verbose," %% E_[0]->M_Zn[0]: "); M_handle_printf(E_[0]->M_Zt,verbose," %% E_[0]->M_Zt[0]: ");
  M_handle_printf(E_[1]->M_Zn,verbose," %% E_[1]->M_Zn[1]: "); M_handle_printf(E_[1]->M_Zt,verbose," %% E_[1]->M_Zt[1]: ");
  M_handle_printf(E_[0]->M_Wn,verbose," %% E_[0]->M_Wn[0]: "); M_handle_printf(E_[0]->M_Wt,verbose," %% E_[0]->M_Wt[0]: ");
  M_handle_printf(E_[1]->M_Wn,verbose," %% E_[1]->M_Wn[1]: "); M_handle_printf(E_[1]->M_Wt,verbose," %% E_[1]->M_Wt[1]: ");
  M_handle_printf(E_[0]->M_Yn,verbose," %% E_[0]->M_Yn[0]: "); M_handle_printf(E_[0]->M_Yt,verbose," %% E_[0]->M_Yt[0]: ");
  M_handle_printf(E_[1]->M_Yn,verbose," %% E_[1]->M_Yn[1]: "); M_handle_printf(E_[1]->M_Yt,verbose," %% E_[1]->M_Yt[1]: ");
  M_handle_printf(E_[0]->M_Tn,verbose," %% E_[0]->M_Tn[0]: "); M_handle_printf(E_[0]->M_Tt,verbose," %% E_[0]->M_Tt[0]: ");
  M_handle_printf(E_[1]->M_Tn,verbose," %% E_[1]->M_Tn[1]: "); M_handle_printf(E_[1]->M_Tt,verbose," %% E_[1]->M_Tt[1]: ");
  M_handle_printf(E_[0]->M_Sn,verbose," %% E_[0]->M_Sn[0]: "); M_handle_printf(E_[0]->M_St,verbose," %% E_[0]->M_St[0]: ");
  M_handle_printf(E_[1]->M_Sn,verbose," %% E_[1]->M_Sn[1]: "); M_handle_printf(E_[1]->M_St,verbose," %% E_[1]->M_St[1]: ");
  raprintf(D->A_p,"double",1,psize(D->A_ncols)," %% D->A_p: "); raprintf(D->A_ajdk,"double",psize(D->A_ncols),AJDK_TOT," %% D->A_ajdk: ");
  raprintf(D->Y_p,"double",1,psize(D->Y_ncols)," %% D->Y_p: "); raprintf(D->Y_ajdk,"double",psize(D->Y_ncols),AJDK_TOT," %% D->Y_ajdk: ");
  wkspace_printf();
}
