#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void dcg_ajdk_copy(struct dcg_ajdk *D,struct dcg_ajdk *D_in)
{
  /* copy from D_in; assuming D has been preallocated. */
  int verbose=0;
  int nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0;
  if (verbose){ printf(" %% [entering dcg_ajdk_copy]\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
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
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  /* copying K_bmc_b */
  memcpy(D->K_bmc_b,D_in->K_bmc_b,D->A_mc_length);
  memcpy(D->K_bmc_j,D_in->K_bmc_j,D->A_mc_length);
  memcpy(D->K_bmc_j_rmv,D_in->K_bmc_j_rmv,D->A_mc_length);
  memcpy(D->K_bmc_j_rtn,D_in->K_bmc_j_rtn,D->A_mc_length);
  D->K_cpop_b= popcount_uchar_array(D->K_bmc_b,D->A_mc_length); D->K_cpop_j = popcount_uchar_array(D->K_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->K_cpop_b %d D->K_cpop_j %d\n",D->K_cpop_b,D->K_cpop_j);}
  memcpy(D->K_umc_b,D_in->K_umc_b,D->A_ncols*sizeof(unsigned char));
  memcpy(D->K_umc_j,D_in->K_umc_j,D->A_ncols*sizeof(unsigned char));
  memcpy(D->K_umc_j_rmv,D_in->K_umc_j_rmv,D->A_ncols*sizeof(unsigned char));
  memcpy(D->K_umc_j_rtn,D_in->K_umc_j_rtn,D->A_ncols*sizeof(unsigned char));
  D->out_xgrow_ij=D_in->out_xgrow_ij;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcg_ajdk_copy] D->A_ncols=%d D->Y_ncols=%d D->T_ncols=%d; \n",D->A_ncols,D->Y_ncols,D->T_ncols);}
}

void dcg_ajdk_load(struct dcg_ajdk *D)
{
  /* load from file */
  int verbose=0;
  int nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0;
  unsigned char *bAb=NULL,*bAj=NULL,*bKj=NULL;
  if (verbose){ printf(" %% [entering dcg_ajdk_load]\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  D->A_ncols=GLOBAL_A_n_cols;D->Y_ncols=GLOBAL_Y_n_cols;D->T_ncols=GLOBAL_T_n_cols;
  /* loading A_bmc_b */
  D->A_cbother = (D->A_ncols>0);
  D->A_ncols_extend = (D->bitj - (D->A_ncols % D->bitj)) % D->bitj; 
  D->A_mc_length = bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% A_ncols %d, A_ncols_extend %d A_mc_length %d\n",D->A_ncols,D->A_ncols_extend,D->A_mc_length);}
  D->A_bmc_b = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  D->A_bmc_j = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  D->A_bmc_j_rmv = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rmv){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  D->A_bmc_j_rtn = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rtn){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  if (GLOBAL_A_n_cind==NULL || !strcmp(GLOBAL_A_n_cind,"\0")){ for (nc=0;nc<D->A_ncols;nc++){ bset__on(D->A_bmc_b,nc);}}
  else{ binary_read(GLOBAL_A_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->A_bmc_b)); if (nrows_tmp!=D->A_ncols){ printf(" %% Warning! A_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_A_n_name,GLOBAL_A_n_cind,nrows_tmp,D->A_ncols);}}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_b: "); if (verbose>2){ bprintf(D->A_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  for (nc=0;nc<D->A_mc_length;nc++){ D->A_bmc_j[nc] = D->A_bmc_b[nc];}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); if (verbose>2){ bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  D->A_cpop_b= popcount_uchar_array(D->A_bmc_b,D->A_mc_length); D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->A_cpop_b %d D->A_cpop_j %d\n",D->A_cpop_b,D->A_cpop_j);}
  D->A_umc_b = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_b[nc] = bget__on(D->A_bmc_b,nc);}
  D->A_umc_j = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);}
  D->A_umc_j_rmv = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rmv){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rmv[nc] = bget__on(D->A_bmc_j_rmv,nc);}
  D->A_umc_j_rtn = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rtn){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = bget__on(D->A_bmc_j_rtn,nc);}
  /* loading Y_bmc_b */
  D->Y_cbother = (D->Y_ncols>0);
  D->Y_ncols_extend = (D->bitj - (D->Y_ncols % D->bitj)) % D->bitj; 
  D->Y_mc_length = bsize(D->Y_ncols)/* rup(D->Y_ncols+D->Y_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% Y_ncols %d, Y_ncols_extend %d Y_mc_length %d\n",D->Y_ncols,D->Y_ncols_extend,D->Y_mc_length);} 
  D->Y_bmc_b = wkspace_all0c(D->Y_mc_length); if (!D->Y_bmc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  D->Y_bmc_j = wkspace_all0c(D->Y_mc_length); if (!D->Y_bmc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  if (GLOBAL_Y_n_cind==NULL || !strcmp(GLOBAL_Y_n_cind,"\0")){ for (nc=0;nc<D->Y_ncols;nc++){ bset__on(D->Y_bmc_b,nc);}}
  else{ binary_read(GLOBAL_Y_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->Y_bmc_b)); if (nrows_tmp!=D->Y_ncols){ printf(" %% Warning! Y_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_Y_n_name,GLOBAL_Y_n_cind,nrows_tmp,D->Y_ncols);}}
  sprintf(D->tmpYnchar," %%%% D->Y_bmc_b: "); if (verbose>2){ bprintf(D->Y_bmc_b,D->bitj,1,D->Y_ncols,D->tmpYnchar);}
  for (nc=0;nc<D->Y_mc_length;nc++){ D->Y_bmc_j[nc] = D->Y_bmc_b[nc];}
  sprintf(D->tmpYnchar," %%%% D->Y_bmc_j: "); if (verbose>2){ bprintf(D->Y_bmc_j,D->bitj,1,D->Y_ncols,D->tmpYnchar);}
  D->Y_cpop_b= popcount_uchar_array(D->Y_bmc_b,D->Y_mc_length); D->Y_cpop_j = popcount_uchar_array(D->Y_bmc_j,D->Y_mc_length); if (verbose>1){ printf(" %%%% D->Y_cpop_b %d D->Y_cpop_j %d\n",D->Y_cpop_b,D->Y_cpop_j);}
  D->Y_umc_b = (unsigned char *) wkspace_all0c(D->Y_ncols*sizeof(unsigned char)); if (!D->Y_umc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->Y_ncols;nc++){ D->Y_umc_b[nc] = bget__on(D->Y_bmc_b,nc);}
  D->Y_umc_j = (unsigned char *) wkspace_all0c(D->Y_ncols*sizeof(unsigned char)); if (!D->Y_umc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->Y_ncols;nc++){ D->Y_umc_j[nc] = bget__on(D->Y_bmc_j,nc);}
  /* loading T_bmc_b */
  D->T_ncols_extend = (D->bitj - (D->T_ncols % D->bitj)) % D->bitj; 
  D->T_mc_length = bsize(D->T_ncols)/* rup(D->T_ncols+D->T_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% T_ncols %d, T_ncols_extend %d T_mc_length %d\n",D->T_ncols,D->T_ncols_extend,D->T_mc_length);}
  D->T_bmc_b = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  D->T_bmc_j = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  if (GLOBAL_T_n_cind==NULL || !strcmp(GLOBAL_T_n_cind,"\0")){ for (nc=0;nc<D->T_ncols;nc++){ bset__on(D->T_bmc_b,nc);}}
  else{ binary_read(GLOBAL_T_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->T_bmc_b)); if (nrows_tmp!=D->T_ncols){ printf(" %% Warning! T_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_T_n_name,GLOBAL_T_n_cind,nrows_tmp,D->T_ncols);}}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_b: "); if (verbose>2){ bprintf(D->T_bmc_b,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  for (nc=0;nc<D->T_mc_length;nc++){ D->T_bmc_j[nc] = D->T_bmc_b[nc];}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_j: "); if (verbose>2){ bprintf(D->T_bmc_j,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  D->T_cpop_b= popcount_uchar_array(D->T_bmc_b,D->T_mc_length); D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length); if (verbose>1){ printf(" %%%% D->T_cpop_b %d D->T_cpop_j %d\n",D->T_cpop_b,D->T_cpop_j);}
  D->T_umc_b = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_b[nc] = bget__on(D->T_bmc_b,nc);}
  D->T_umc_j = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_j[nc] = bget__on(D->T_bmc_j,nc);}
  D->QC_TYnWtS_nrm = NULL; D->QC_TAnZtS_nrm = NULL; D->QC_TYnYtT_nrm = NULL; D->QC_TAnAtT_nrm = NULL;
  D->QC_A_svalue=NULL; D->QC_A_index_local_mc_a=NULL; D->QC_A_index_local_mc_b=NULL; D->QC_A_index_local_mc_j=NULL;
  D->QR_A_svalue=NULL; D->QR_A_index_local_mr_a=NULL; D->QR_A_index_local_mr_b=NULL; D->QR_A_index_local_mr_j=NULL; D->QR_A_index_local_nb=NULL; D->QR_A_index_global_mr_a=NULL; D->QR_A_index_global_mr_b=NULL; D->QR_A_index_global_mr_j=NULL;
  D->Irem=0;D->Ireq=0; D->out_iteration=0; D->out_xdrop_ij=0; D->out_trace=NULL; D->out_xdrop_a=NULL; D->out_xdrop_b=NULL; D->out_trace_length = 0;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  /* loading K_bmc_b */
  D->K_bmc_b = wkspace_all0c(D->A_mc_length); if (!D->K_bmc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  D->K_bmc_j = wkspace_all0c(D->A_mc_length); if (!D->K_bmc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  D->K_bmc_j_rmv = wkspace_all0c(D->A_mc_length); if (!D->K_bmc_j_rmv){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  D->K_bmc_j_rtn = wkspace_all0c(D->A_mc_length); if (!D->K_bmc_j_rtn){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} 
  if (GLOBAL_A_n_cind==NULL || !strcmp(GLOBAL_A_n_cind,"\0")){ for (nc=0;nc<D->A_ncols;nc++){ bset__on(D->K_bmc_b,nc);}}
  else{ binary_read(GLOBAL_A_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->K_bmc_b)); if (nrows_tmp!=D->A_ncols){ printf(" %% Warning! A_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_A_n_name,GLOBAL_A_n_cind,nrows_tmp,D->A_ncols);}}
  sprintf(D->tmpAnchar," %%%% D->K_bmc_b: "); if (verbose>2){ bprintf(D->K_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  if (GLOBAL_K_n_cind==NULL || !strcmp(GLOBAL_K_n_cind,"\0")){ ura_dup(D->K_bmc_j,D->A_mc_length,D->K_bmc_b);}
  else{ binary_read(GLOBAL_K_n_cind,&bitj_tmp,&nrows_tmp,&ncols_tmp,&(D->K_bmc_j)); if (nrows_tmp!=D->A_ncols){ printf(" %% Warning! A_: %s; improper mc_b_ %s, nrows %d instead of %d\n",GLOBAL_A_n_name,GLOBAL_K_n_cind,nrows_tmp,D->A_ncols);}}
  sprintf(D->tmpAnchar," %%%% D->K_bmc_j: "); if (verbose>2){ bprintf(D->K_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  D->K_cpop_b= popcount_uchar_array(D->K_bmc_b,D->A_mc_length); D->K_cpop_j = popcount_uchar_array(D->K_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->K_cpop_b %d D->K_cpop_j %d\n",D->K_cpop_b,D->K_cpop_j);}
  /* redefining D->A_bmc_j */
  bAb = D->A_bmc_b; bAj = D->A_bmc_j; bKj = D->K_bmc_j;
  for (nc=0;nc<D->A_ncols;nc++){ if (bget__on(bAb,nc) && bget__on(bKj,nc)){ bset_off(bAj,nc);} if (bget__on(bAb,nc) && bget_off(bKj,nc)){ bset__on(bAj,nc);}} 
  D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length);
  D->K_umc_b = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->K_umc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->K_umc_b[nc] = bget__on(D->K_bmc_b,nc);}
  D->K_umc_j = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->K_umc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->K_umc_j[nc] = bget__on(D->K_bmc_j,nc);}
  D->K_umc_j_rmv = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->K_umc_j_rmv){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->K_umc_j_rmv[nc] = bget__on(D->K_bmc_j_rmv,nc);}
  D->K_umc_j_rtn = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->K_umc_j_rtn){ printf(" %% Warning! not enough memory in dcg_ajdk_load\n");} for (nc=0;nc<D->A_ncols;nc++){ D->K_umc_j_rtn[nc] = bget__on(D->K_bmc_j_rtn,nc);}
  /* initializing E->QC_RJnAtT_nrm */
  D->QC_RUnWtS_nrm = NULL; D->QC_RUnYtT_nrm = NULL; D->QC_RJnZtS_nrm = NULL; D->QC_RJnAtT_nrm = NULL; D->QC_TKnVtS_nrm = NULL; D->QC_TKnKtT_nrm = NULL;
  D->QC_K_svalue=NULL; D->QC_K_index_local_mc_a=NULL; D->QC_K_index_local_mc_b=NULL; D->QC_K_index_local_mc_j=NULL;
  D->QR_J_svalue=NULL; D->QR_J_index_local_mr_a=NULL; D->QR_J_index_local_mr_b=NULL; D->QR_J_index_local_mr_j=NULL; D->QR_J_index_local_nb=NULL; D->QR_J_index_global_mr_a=NULL; D->QR_J_index_global_mr_b=NULL; D->QR_J_index_global_mr_j=NULL;
  D->out_xgrow_ij=0; D->out_xgrow_a=NULL; D->out_xgrow_b=NULL;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcg_ajdk_load] D->A_ncols=%d D->Y_ncols=%d D->T_ncols=%d \n",D->A_ncols,D->Y_ncols,D->T_ncols);}
}

void dcg_ajdk_init(double mrnd,struct dcg_ajdk *D)
{
  /* initialize based on test parameters */
  int verbose=0;
  int nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0;
  unsigned char *bb=NULL,*bj=NULL;
  double dtmp=0; int nx=0;
  if (verbose){ printf(" %% [entering dcg_ajdk_init]\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  D->A_ncols=GLOBAL_TEST_A_n_cols;D->Y_ncols=GLOBAL_TEST_Y_n_cols;D->T_ncols=GLOBAL_TEST_T_n_cols;
  /* initializing A_bmc_b */
  D->A_cbother = (D->A_ncols>0);
  D->A_ncols_extend = (D->bitj - (D->A_ncols % D->bitj)) % D->bitj; 
  D->A_mc_length = bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% A_ncols %d, A_ncols_extend %d A_mc_length %d\n",D->A_ncols,D->A_ncols_extend,D->A_mc_length);}
  D->A_bmc_b = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  D->A_bmc_j = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  bb = D->A_bmc_b; bj = D->A_bmc_j; for (nx=0;nx<D->A_ncols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_b: "); if (verbose>2){ bprintf(D->A_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); if (verbose>2){ bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  D->A_bmc_j_rmv = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rmv){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  D->A_bmc_j_rtn = wkspace_all0c(D->A_mc_length); if (!D->A_bmc_j_rtn){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  D->A_cpop_b= popcount_uchar_array(D->A_bmc_b,D->A_mc_length); D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->A_cpop_b %d D->A_cpop_j %d\n",D->A_cpop_b,D->A_cpop_j);}
  D->A_umc_b = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_b[nc] = bget__on(D->A_bmc_b,nc);}
  D->A_umc_j = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);}
  D->A_umc_j_rmv = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rmv){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rmv[nc] = bget__on(D->A_bmc_j_rmv,nc);}
  D->A_umc_j_rtn = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->A_umc_j_rtn){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = bget__on(D->A_bmc_j_rtn,nc);}
  /* initializing Y_bmc_b */
  D->Y_cbother = (D->Y_ncols>0);
  D->Y_ncols_extend = (D->bitj - (D->Y_ncols % D->bitj)) % D->bitj; 
  D->Y_mc_length = bsize(D->Y_ncols)/* rup(D->Y_ncols+D->Y_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% Y_ncols %d, Y_ncols_extend %d Y_mc_length %d\n",D->Y_ncols,D->Y_ncols_extend,D->Y_mc_length);}
  D->Y_bmc_b = wkspace_all0c(D->Y_mc_length); if (!D->Y_bmc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  D->Y_bmc_j = wkspace_all0c(D->Y_mc_length); if (!D->Y_bmc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  bb = D->Y_bmc_b; bj = D->Y_bmc_j; for (nx=0;nx<D->Y_ncols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpYnchar," %%%% D->Y_bmc_b: "); if (verbose>2){ bprintf(D->Y_bmc_b,D->bitj,1,D->Y_ncols,D->tmpYnchar);}
  sprintf(D->tmpYnchar," %%%% D->Y_bmc_j: "); if (verbose>2){ bprintf(D->Y_bmc_j,D->bitj,1,D->Y_ncols,D->tmpYnchar);}
  D->Y_cpop_b= popcount_uchar_array(D->Y_bmc_b,D->Y_mc_length); D->Y_cpop_j = popcount_uchar_array(D->Y_bmc_j,D->Y_mc_length); if (verbose>1){ printf(" %%%% D->Y_cpop_b %d D->Y_cpop_j %d\n",D->Y_cpop_b,D->Y_cpop_j);}
  D->Y_umc_b = (unsigned char *) wkspace_all0c(D->Y_ncols*sizeof(unsigned char)); if (!D->Y_umc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->Y_ncols;nc++){ D->Y_umc_b[nc] = bget__on(D->Y_bmc_b,nc);}
  D->Y_umc_j = (unsigned char *) wkspace_all0c(D->Y_ncols*sizeof(unsigned char)); if (!D->Y_umc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->Y_ncols;nc++){ D->Y_umc_j[nc] = bget__on(D->Y_bmc_j,nc);}
  /* initializing T_bmc_b */
  D->T_ncols_extend = (D->bitj - (D->T_ncols % D->bitj)) % D->bitj; 
  D->T_mc_length = bsize(D->T_ncols)/* rup(D->T_ncols+D->T_ncols_extend,POPLENGTH)/BIT8 */;
  if (verbose>1){ printf(" %%%% T_ncols %d, T_ncols_extend %d T_mc_length %d\n",D->T_ncols,D->T_ncols_extend,D->T_mc_length);}
  D->T_bmc_b = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  D->T_bmc_j = wkspace_all0c(D->T_mc_length); if (!D->T_bmc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  bb = D->T_bmc_b; bj = D->T_bmc_j; for (nx=0;nx<D->T_ncols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}} bset__on(bb,0); bset__on(bj,0); bset__on(bb,D->T_ncols-1), bset__on(bj,D->T_ncols-1);
  sprintf(D->tmpTnchar," %%%% D->T_bmc_b: "); if (verbose>2){ bprintf(D->T_bmc_b,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_j: "); if (verbose>2){ bprintf(D->T_bmc_j,D->bitj,1,D->T_ncols,D->tmpTnchar);}
  D->T_cpop_b= popcount_uchar_array(D->T_bmc_b,D->T_mc_length); D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length); if (verbose>1){ printf(" %%%% D->T_cpop_b %d D->T_cpop_j %d\n",D->T_cpop_b,D->T_cpop_j);}
  D->T_umc_b = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_b[nc] = bget__on(D->T_bmc_b,nc);}
  D->T_umc_j = (unsigned char *) wkspace_all0c(D->T_ncols*sizeof(unsigned char)); if (!D->T_umc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_j[nc] = bget__on(D->T_bmc_j,nc);}
  D->QC_TYnWtS_nrm = NULL; D->QC_TAnZtS_nrm = NULL; D->QC_TYnYtT_nrm = NULL; D->QC_TAnAtT_nrm = NULL;
  D->QC_A_svalue=NULL; D->QC_A_index_local_mc_a=NULL; D->QC_A_index_local_mc_b=NULL; D->QC_A_index_local_mc_j=NULL;
  D->QR_A_svalue=NULL; D->QR_A_index_local_mr_a=NULL; D->QR_A_index_local_mr_b=NULL; D->QR_A_index_local_mr_j=NULL; D->QR_A_index_local_nb=NULL; D->QR_A_index_global_mr_a=NULL; D->QR_A_index_global_mr_b=NULL; D->QR_A_index_global_mr_j=NULL;
  D->Irem=0;D->Ireq=0; D->out_iteration=0; D->out_xdrop_ij=0; D->out_trace=NULL; D->out_xdrop_a=NULL; D->out_xdrop_b=NULL; D->out_trace_length = 0;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  D->K_bmc_b = wkspace_all0c(D->A_mc_length); if (!D->K_bmc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  D->K_bmc_j = wkspace_all0c(D->A_mc_length); if (!D->K_bmc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  ura_dup(D->K_bmc_b,D->A_mc_length,D->A_bmc_b);
  bj = D->K_bmc_j; for (nx=0;nx<D->A_ncols;nx++){ if (bget__on(D->A_bmc_b,nx) && bget_off(D->A_bmc_j,nx)){ bset__on(bj,nx);}}
  sprintf(D->tmpAnchar," %%%% D->K_bmc_b: "); if (verbose>2){ bprintf(D->K_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% D->K_bmc_j: "); if (verbose>2){ bprintf(D->K_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  D->K_bmc_j_rmv = wkspace_all0c(D->A_mc_length); if (!D->K_bmc_j_rmv){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  D->K_bmc_j_rtn = wkspace_all0c(D->A_mc_length); if (!D->K_bmc_j_rtn){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} 
  D->K_cpop_b= popcount_uchar_array(D->K_bmc_b,D->A_mc_length); D->K_cpop_j = popcount_uchar_array(D->K_bmc_j,D->A_mc_length); if (verbose>1){ printf(" %%%% D->K_cpop_b %d D->K_cpop_j %d\n",D->K_cpop_b,D->K_cpop_j);}
  D->K_umc_b = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->K_umc_b){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->K_umc_b[nc] = bget__on(D->K_bmc_b,nc);}
  D->K_umc_j = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->K_umc_j){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->K_umc_j[nc] = bget__on(D->K_bmc_j,nc);}
  D->K_umc_j_rmv = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->K_umc_j_rmv){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->K_umc_j_rmv[nc] = bget__on(D->K_bmc_j_rmv,nc);}
  D->K_umc_j_rtn = (unsigned char *) wkspace_all0c(D->A_ncols*sizeof(unsigned char)); if (!D->K_umc_j_rtn){ printf(" %% Warning! not enough memory in dcg_ajdk_init\n");} for (nc=0;nc<D->A_ncols;nc++){ D->K_umc_j_rtn[nc] = bget__on(D->K_bmc_j_rtn,nc);}
  /* initializing E->QC_RJnAtT_nrm */
  D->QC_RUnWtS_nrm = NULL; D->QC_RUnYtT_nrm = NULL; D->QC_RJnZtS_nrm = NULL; D->QC_RJnAtT_nrm = NULL; D->QC_TKnVtS_nrm = NULL; D->QC_TKnKtT_nrm = NULL;
  D->QC_K_svalue=NULL; D->QC_K_index_local_mc_a=NULL; D->QC_K_index_local_mc_b=NULL; D->QC_K_index_local_mc_j=NULL;
  D->QR_J_svalue=NULL; D->QR_J_index_local_mr_a=NULL; D->QR_J_index_local_mr_b=NULL; D->QR_J_index_local_mr_j=NULL; D->QR_J_index_local_nb=NULL; D->QR_J_index_global_mr_a=NULL; D->QR_J_index_global_mr_b=NULL; D->QR_J_index_global_mr_j=NULL;
  D->out_xgrow_ij=0; D->out_xgrow_a=NULL; D->out_xgrow_b=NULL;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcg_ajdk_init] D->A_ncols=%d D->Y_ncols=%d D->T_ncols=%d \n",D->A_ncols,D->Y_ncols,D->T_ncols);}
}

void dcg_single_copy_M_An(struct dcg_single *E,struct dcg_single *E_in)
{
  /* copy from variable */
  int verbose=0;
  struct dcg_ajdk *D=E->D,*D_in=E_in->D;
  int nb_given = E->nb;
  int nr=0,nb=0,nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0,nrows_tmp_extend,brows_tmp;
  unsigned char *bXra_tmp=NULL;
  if (verbose){ printf(" %% [entering dcg_single_copy_M_An]\n");}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
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
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* copying J */
  memcpy(E->J_bmr_b,E_in->J_bmr_b,E->A_mr_length);
  memcpy(E->J_bmr_j,E_in->J_bmr_j,E->A_mr_length);
  memcpy(E->J_bmr_j_rmv,E_in->J_bmr_j_rmv,E->A_mr_length);
  memcpy(E->J_bmr_j_rtn,E_in->J_bmr_j_rtn,E->A_mr_length);
  E->J_rpop_b = popcount_uchar_array(E->J_bmr_b,E->A_mr_length);
  E->J_rpop_j = popcount_uchar_array(E->J_bmr_j,E->A_mr_length);
  if (verbose>1){ printf(" %% A_nrows_extend_[%d] %d\n",nb,E->A_nrows_extend);}
  if (verbose>1){ printf(" %% A_mr_length_[%d] %d\n",nb,E->A_mr_length);}
  sprintf(D->tmpAnchar," %%%% E->J_bmr_b_[%.2d]:",nb); if (verbose>1){ bprintf(E->J_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% E->J_bmr_j_[%.2d]:",nb); if (verbose>1){ bprintf(E->J_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  if (verbose>1){ printf(" %% A_nrows_[%d]  = %d\n",nb,E->A_nrows);}
  if (verbose>1){ printf(" %% J_rpop_b_[%d] = %d\n",nb,E->J_rpop_b);}
  if (verbose>1){ printf(" %% J_rpop_j_[%d] = %d\n",nb,E->J_rpop_j);}
  memcpy(E->J_umr_b,E_in->J_umr_b,E->A_nrows*sizeof(unsigned char));
  memcpy(E->J_umr_j,E_in->J_umr_j,E->A_nrows*sizeof(unsigned char));
  memcpy(E->J_umr_j_rmv,E_in->J_umr_j_rmv,E->A_nrows*sizeof(unsigned char));
  memcpy(E->J_umr_j_rtn,E_in->J_umr_j_rtn,E->A_nrows*sizeof(unsigned char));
  if (verbose>1){ printf(" %% copying M_Jn_[%d]\n",nb);}
  M_handle_copy(E->M_Jn,E_in->M_Jn);
  M_handle_copy(E->M_Jt,E_in->M_Jt);
  if (verbose>1){ printf(" %% copying M_Kn_[%d]\n",nb);}
  M_handle_copy(E->M_Kn,E_in->M_Kn);
  M_handle_copy(E->M_Kt,E_in->M_Kt);
  if (verbose>1){ printf(" %% copying M_Vn_[%d]\n",nb);}
  M_handle_copy(E->M_Vn,E_in->M_Vn);
  M_handle_copy(E->M_Vt,E_in->M_Vt);
  /* copying U */
  M_handle_copy(E->M_Un,E_in->M_Un);
  M_handle_copy(E->M_Ut,E_in->M_Ut);
  /* copying R */
  M_handle_copy(E->M_Rn,E_in->M_Rn);
  M_handle_copy(E->M_Rt,E_in->M_Rt);
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished dcg_single_copy_M_An]\n");}
}

void dcg_single_load_M_An(struct dcg_single *E)
{
  /* load from file */
  int verbose=0;
  struct dcg_ajdk *D=E->D;
  int nb_given = E->nb;
  int nr=0,nb=0,nc=0,nx=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0,nrows_tmp_extend,brows_tmp;
  unsigned char *bXra_tmp=NULL;
  unsigned char *bAb=NULL,*bAj=NULL,*bJj=NULL;
  unsigned char *bn=NULL,*bt=NULL,*bn2=NULL,*bt2=NULL;
  if (verbose){ printf(" %% [entering dcg_single_load_M_An]\n");}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */  
  /* loading A_bmr_b etc */
  nb = nb_given;
  /* loading A */
  binary_read_getsize(GLOBAL_A_t_name_[nb],&(D->bitj),&(D->A_ncols),&(E->A_nrows));
  E->A_rbother = (E->A_nrows>0);
  E->A_nrows_extend = (D->bitj - (E->A_nrows % D->bitj)) % D->bitj;
  E->A_mr_length = bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */;
  E->A_bmr_b = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_b){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
  E->A_bmr_j = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
  E->A_bmr_j_rmv = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rmv){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
  E->A_bmr_j_rtn = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rtn){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
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
  E->A_umr_b = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_b){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_b[nr] = bget__on(E->A_bmr_b,nr);}
  E->A_umr_j = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);}
  E->A_umr_j_rmv = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rmv){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rmv[nr] = bget__on(E->A_bmr_j_rmv,nr);}
  E->A_umr_j_rtn = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rtn){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = bget__on(E->A_bmr_j_rtn,nr);}
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
  E->Z_bmr_b = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_b){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
  E->Z_bmr_j = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_j){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
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
  E->Z_umr_b = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_b){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_b[nr] = bget__on(E->Z_bmr_b,nr);}
  E->Z_umr_j = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_j){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_j[nr] = bget__on(E->Z_bmr_j,nr);}
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
  E->QR_TYnWtS_nrm = NULL; E->QR_TAnZtS_nrm = NULL; E->QR_TYnYtT_nrm = NULL; E->QR_TAnAtT_nrm = NULL; E->QR_A_index_global_mr_a=NULL; E->QR_A_index_global_mr_b=NULL;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */  
  E->J_bmr_b = wkspace_all0c(E->A_mr_length); if (!E->J_bmr_b){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
  E->J_bmr_j = wkspace_all0c(E->A_mr_length); if (!E->J_bmr_j){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
  E->J_bmr_j_rmv = wkspace_all0c(E->A_mr_length); if (!E->J_bmr_j_rmv){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
  E->J_bmr_j_rtn = wkspace_all0c(E->A_mr_length); if (!E->J_bmr_j_rtn){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");}
  if (GLOBAL_A_n_rind_[nb]==NULL || !strcmp(GLOBAL_A_n_rind_[nb],"\0")){ ura_dup(E->J_bmr_b,E->A_mr_length,E->A_bmr_b);}
  else{ binary_read(GLOBAL_A_n_rind_[nb],&bitj_tmp,&nrows_tmp,&ncols_tmp,&(E->J_bmr_b)); if (nrows_tmp!=E->A_nrows || bitj_tmp!=D->bitj){ printf(" %% Warning! J_: %s; improper mr_b_ %s, nrows %d instead of %d, bitj %d vs %d\n",GLOBAL_A_t_name_[nb],GLOBAL_A_n_rind_[nb],nrows_tmp,E->A_nrows,bitj_tmp,D->bitj);}}
  if (GLOBAL_J_n_rind_[nb]==NULL || !strcmp(GLOBAL_J_n_rind_[nb],"\0")){ ura_dup(E->J_bmr_j,E->A_mr_length,E->J_bmr_b);}
  else{ binary_read(GLOBAL_J_n_rind_[nb],&bitj_tmp,&nrows_tmp,&ncols_tmp,&(E->J_bmr_j)); if (nrows_tmp!=E->A_nrows || bitj_tmp!=D->bitj){ printf(" %% Warning! J_: %s; improper mr_b_ %s, nrows %d instead of %d, bitj %d vs %d\n",GLOBAL_A_t_name_[nb],GLOBAL_J_n_rind_[nb],nrows_tmp,E->A_nrows,bitj_tmp,D->bitj);}}
  E->J_rpop_b = popcount_uchar_array(E->J_bmr_b,E->A_mr_length);
  E->J_rpop_j = popcount_uchar_array(E->J_bmr_j,E->A_mr_length);
  /* redefining E->A_bmr_j */
  bAb = E->A_bmr_b; bAj = E->A_bmr_j; bJj = E->J_bmr_j;
  for (nr=0;nr<E->A_nrows;nr++){ if (bget__on(bAb,nr) && bget__on(bJj,nr)){ bset_off(bAj,nr);} if (bget__on(bAb,nr) && bget_off(bJj,nr)){ bset__on(bAj,nr);}} 
  E->A_rpop_j = popcount_uchar_array(E->A_bmr_j,E->A_mr_length);
  E->J_umr_b = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->J_umr_b){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->J_umr_b[nr] = bget__on(E->J_bmr_b,nr);}
  E->J_umr_j = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->J_umr_j){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->J_umr_j[nr] = bget__on(E->J_bmr_j,nr);}
  E->J_umr_j_rmv = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->J_umr_j_rmv){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->J_umr_j_rmv[nr] = bget__on(E->J_bmr_j_rmv,nr);}
  E->J_umr_j_rtn = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->J_umr_j_rtn){ printf(" %% Warning! not enough memory in dcg_single_load_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->J_umr_j_rtn[nr] = bget__on(E->J_bmr_j_rtn,nr);}
  E->J_rpop_b = popcount_uchar_array(E->J_bmr_b,E->A_mr_length);
  E->J_rpop_j = popcount_uchar_array(E->J_bmr_j,E->A_mr_length);
  /* initializing J,K,U,R,V */
  E->M_Jn = NULL; E->M_Jt = NULL;
  E->M_Kn = NULL; E->M_Kt = NULL;
  E->M_Un = NULL; E->M_Ut = NULL; E->M_Rn = NULL; E->M_Rt = NULL;
  E->M_Vn = NULL; E->M_Vt = NULL; 
  /* E->M_Jn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->J_bmr_b,D->A_bmc_b); ura_dup(E->M_Jn->mr_j,E->M_Jn->mr_length,E->J_bmr_j); ura_dup(E->M_Jn->mc_j,E->M_Jn->mc_length,D->A_bmc_j); M_mxget(E->M_Jn); */
  /* E->M_Jt = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,D->A_bmc_b,E->J_bmr_b); ura_dup(E->M_Jt->mr_j,E->M_Jt->mr_length,D->A_bmc_j); ura_dup(E->M_Jt->mc_j,E->M_Jt->mc_length,E->J_bmr_j); M_mxget(E->M_Jt); */
  /* E->M_Kn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->A_bmr_b,D->K_bmc_b); ura_dup(E->M_Kn->mr_j,E->M_Kn->mr_length,E->A_bmr_j); ura_dup(E->M_Kn->mc_j,E->M_Kn->mc_length,D->K_bmc_j); M_mxget(E->M_Kn); */
  /* E->M_Kt = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,D->K_bmc_b,E->A_bmr_b); ura_dup(E->M_Kt->mr_j,E->M_Kt->mr_length,D->K_bmc_j); ura_dup(E->M_Kt->mc_j,E->M_Kt->mc_length,E->A_bmr_j); M_mxget(E->M_Kt); */
  /* E->M_Un = M_handle_v_make(D->bitj,E->A_nrows,D->Y_ncols,NULL,E->M_Yn->wX,E->J_bmr_b,D->Y_bmc_b); ura_dup(E->M_Un->mr_j,E->M_Un->mr_length,E->J_bmr_j); ura_dup(E->M_Un->mc_j,E->M_Un->mc_length,D->Y_bmc_j); M_mxget(E->M_Un); */
  /* E->M_Ut = M_handle_v_make(D->bitj,D->Y_ncols,E->A_nrows,NULL,E->M_Yt->wX,D->Y_bmc_b,E->J_bmr_b); ura_dup(E->M_Ut->mr_j,E->M_Ut->mr_length,D->Y_bmc_j); ura_dup(E->M_Ut->mc_j,E->M_Ut->mc_length,E->J_bmr_j); M_mxget(E->M_Ut); */
  /* E->M_Rn = M_handle_v_make(D->bitj,E->A_nrows,D->T_ncols,NULL,E->M_Tn->wX,E->J_bmr_b,D->T_bmc_b); ura_dup(E->M_Rn->mr_j,E->M_Rn->mr_length,E->J_bmr_j); ura_dup(E->M_Rn->mc_j,E->M_Rn->mc_length,D->T_bmc_j); M_mxget(E->M_Rn); */
  /* E->M_Rt = M_handle_v_make(D->bitj,D->T_ncols,E->A_nrows,NULL,E->M_Tt->wX,D->T_bmc_b,E->J_bmr_b); ura_dup(E->M_Rt->mr_j,E->M_Rt->mr_length,D->T_bmc_j); ura_dup(E->M_Rt->mc_j,E->M_Rt->mc_length,E->J_bmr_j); M_mxget(E->M_Rt); */
  /* E->M_Vn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,NULL,E->M_Zn->wX,E->Z_bmr_b,D->K_bmc_b); ura_dup(E->M_Vn->mr_j,E->M_Vn->mr_length,E->Z_bmr_j); ura_dup(E->M_Vn->mc_j,E->M_Vn->mc_length,D->K_bmc_j); M_mxget(E->M_Vn); */
  /* E->M_Vt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,NULL,E->M_Zt->wX,D->K_bmc_b,E->Z_bmr_b); ura_dup(E->M_Vt->mr_j,E->M_Vt->mr_length,D->K_bmc_j); ura_dup(E->M_Vt->mc_j,E->M_Vt->mc_length,E->Z_bmr_j); M_mxget(E->M_Vt); */
  bn = E->M_An->wX; bt = E->M_At->wX; E->M_Jn = NULL; E->M_Jt = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->A_ncols,E->J_bmr_b,D->A_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->A_ncols,bn,bt,E->J_bmr_b,E->J_bmr_j,D->A_bmc_b,D->A_bmc_j,&(E->M_Jn),&(E->M_Jt));
  bn  = E->M_An->wX; bt  = E->M_At->wX; E->M_Kn = NULL; E->M_Kt = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->A_ncols,E->A_bmr_b,D->K_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->A_ncols,bn,bt,E->A_bmr_b,E->A_bmr_j,D->K_bmc_b,D->K_bmc_j,&(E->M_Kn),&(E->M_Kt));
  bn  = E->M_Yn->wX; bt  = E->M_Yt->wX; E->M_Un = NULL; E->M_Ut = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->Y_ncols,E->J_bmr_b,D->Y_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->Y_ncols,bn,bt,E->J_bmr_b,E->J_bmr_j,D->Y_bmc_b,D->Y_bmc_j,&(E->M_Un),&(E->M_Ut));
  bn  = E->M_Zn->wX; bt  = E->M_Zt->wX; E->M_Vn = NULL; E->M_Vt = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->Z_nrows,D->A_ncols,E->Z_bmr_b,D->K_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->Z_nrows,D->A_ncols,bn,bt,E->Z_bmr_b,E->Z_bmr_j,D->K_bmc_b,D->K_bmc_j,&(E->M_Vn),&(E->M_Vt));
  bn  = E->M_Tn->wX; bt  = E->M_Tt->wX; E->M_Rn = NULL; E->M_Rt = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->T_ncols,E->J_bmr_b,D->T_bmc_b,&(bn),&(bt));
  for (nx=0;nx<E->A_rpop_b;nx++){ bn2 = &(bn[nx*bsize(D->T_ncols)]); bset__on(bn2,0);}
  for (nx=0;nx<E->A_nrows;nx++){ bt2 = &(bt[0*bsize(E->A_nrows)]); bset__on(bt2,nx);}
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->T_ncols,bn,bt,E->J_bmr_b,E->J_bmr_j,D->T_bmc_b,D->T_bmc_j,&(E->M_Rn),&(E->M_Rt));
  /* initializing E->QR_RJnAtT_nrm */
  E->QR_RUnWtS_nrm = NULL; E->QR_RUnYtT_nrm = NULL; E->QR_RJnZtS_nrm = NULL; E->QR_RJnAtT_nrm = NULL; E->QR_J_index_global_mr_a=NULL; E->QR_J_index_global_mr_b=NULL;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */  
  if (verbose){ printf(" %% [finished dcg_single_load_M_An]\n");}
}

void dcg_single_init_M_An(char *error_vs_speed,double mrnd,struct dcg_single *E)
{
  /* initialize based on test parameters */
  int verbose=0;
  struct dcg_ajdk *D=E->D;
  int nb_given = E->nb;
  int nr=0,nb=0,nc=0;
  int bitj_tmp=0,nrows_tmp=0,ncols_tmp=0,nrows_tmp_extend,brows_tmp;
  unsigned char *bXra_tmp=NULL;
  unsigned char *bb=NULL,*bj=NULL,*bn=NULL,*bt=NULL,*bn2=NULL,*bt2=NULL;
  double dtmp=0; int nx=0;
  if (verbose){ printf(" %% [entering dcg_single_init_M_An]\n");}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  nb = nb_given;
  /* initializing A */
  E->A_nrows = maximum(0,GLOBAL_TEST_A_n_rows + (nb_given % 3) - 1);
  E->A_rbother = (E->A_nrows>0);
  E->A_nrows_extend = (D->bitj - (E->A_nrows % D->bitj)) % D->bitj;
  E->A_mr_length = bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */;
  E->A_bmr_b = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_b){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  E->A_bmr_j = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  bb = E->A_bmr_b; bj = E->A_bmr_j; for (nx=0;nx<E->A_nrows;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpAnchar," %%%% E->A_bmr_b: "); if (verbose>2){ bprintf(E->A_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% E->A_bmr_j: "); if (verbose>2){ bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  E->A_bmr_j_rmv = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rmv){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  E->A_bmr_j_rtn = wkspace_all0c(E->A_mr_length); if (!E->A_bmr_j_rtn){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  E->A_rpop_b = popcount_uchar_array(E->A_bmr_b,E->A_mr_length);
  E->A_rpop_j = popcount_uchar_array(E->A_bmr_j,E->A_mr_length);
  E->A_umr_b = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_b){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_b[nr] = bget__on(E->A_bmr_b,nr);}
  E->A_umr_j = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);}
  E->A_umr_j_rmv = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rmv){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rmv[nr] = bget__on(E->A_bmr_j_rmv,nr);}
  E->A_umr_j_rtn = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->A_umr_j_rtn){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = bget__on(E->A_bmr_j_rtn,nr);}
  if (verbose>1){ printf(" %% generating E->M_An\n");}
  bn = NULL; bt = NULL; E->M_An = NULL; E->M_At = NULL;
  wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->A_ncols,E->A_bmr_b,D->A_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->A_ncols,bn,bt,E->A_bmr_b,E->A_bmr_j,D->A_bmc_b,D->A_bmc_j,&(E->M_An),&(E->M_At));
  /* initializing Z */
  E->Z_nrows = maximum(0,GLOBAL_TEST_Z_n_rows + (nb_given % 5) - 2);
  E->Z_rbother = (E->Z_nrows>0);
  E->Z_nrows_extend = (D->bitj - (E->Z_nrows % D->bitj)) % D->bitj;
  E->Z_mr_length = bsize(E->Z_nrows)/* rup(E->Z_nrows+E->Z_nrows_extend,POPLENGTH)/BIT8 */;
  E->Z_bmr_b = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_b){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  E->Z_bmr_j = wkspace_all0c(E->Z_mr_length); if (!E->Z_bmr_j){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  bb = E->Z_bmr_b; bj = E->Z_bmr_j; for (nx=0;nx<E->Z_nrows;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_b: "); if (verbose>2){ bprintf(E->Z_bmr_b,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  sprintf(D->tmpZnchar," %%%% E->Z_bmr_j: "); if (verbose>2){ bprintf(E->Z_bmr_j,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
  E->Z_rpop_b = popcount_uchar_array(E->Z_bmr_b,E->Z_mr_length);
  E->Z_rpop_j = popcount_uchar_array(E->Z_bmr_j,E->Z_mr_length);
  E->Z_umr_b = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_b){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_b[nr] = bget__on(E->Z_bmr_b,nr);}
  E->Z_umr_j = (unsigned char *) wkspace_all0c(E->Z_nrows*sizeof(unsigned char)); if (!E->Z_umr_j){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_j[nr] = bget__on(E->Z_bmr_j,nr);}
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
  /* initializing E->QR_TAnZtS_nrm */
  E->QR_TYnWtS_nrm = NULL; E->QR_TAnZtS_nrm = NULL; E->QR_TYnYtT_nrm = NULL; E->QR_TAnAtT_nrm = NULL; E->QR_A_index_global_mr_a=NULL; E->QR_A_index_global_mr_b=NULL;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  E->J_bmr_b = wkspace_all0c(E->A_mr_length); if (!E->J_bmr_b){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  E->J_bmr_j = wkspace_all0c(E->A_mr_length); if (!E->J_bmr_j){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  ura_dup(E->J_bmr_b,E->A_mr_length,E->A_bmr_b);
  bj = E->J_bmr_j; for (nx=0;nx<E->A_nrows;nx++){ if (bget__on(E->A_bmr_b,nx) && bget_off(E->A_bmr_j,nx)){ bset__on(bj,nx);}}
  sprintf(D->tmpAnchar," %%%% E->J_bmr_b: "); if (verbose>2){ bprintf(E->J_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  sprintf(D->tmpAnchar," %%%% E->J_bmr_j: "); if (verbose>2){ bprintf(E->J_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
  E->J_bmr_j_rmv = wkspace_all0c(E->A_mr_length); if (!E->J_bmr_j_rmv){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  E->J_bmr_j_rtn = wkspace_all0c(E->A_mr_length); if (!E->J_bmr_j_rtn){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");}
  E->J_rpop_b = popcount_uchar_array(E->J_bmr_b,E->A_mr_length);
  E->J_rpop_j = popcount_uchar_array(E->J_bmr_j,E->A_mr_length);
  E->J_umr_b = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->J_umr_b){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->J_umr_b[nr] = bget__on(E->J_bmr_b,nr);}
  E->J_umr_j = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->J_umr_j){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->J_umr_j[nr] = bget__on(E->J_bmr_j,nr);}
  E->J_umr_j_rmv = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->J_umr_j_rmv){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->J_umr_j_rmv[nr] = bget__on(E->J_bmr_j_rmv,nr);}
  E->J_umr_j_rtn = (unsigned char *) wkspace_all0c(E->A_nrows*sizeof(unsigned char)); if (!E->J_umr_j_rtn){ printf(" %% Warning! not enough memory in dcg_single_init_M_An\n");} for (nr=0;nr<E->A_nrows;nr++){ E->J_umr_j_rtn[nr] = bget__on(E->J_bmr_j_rtn,nr);}
  if (verbose>1){ printf(" %% generating E->M_An\n");}
  E->M_Jn = NULL; E->M_Jt = NULL;
  E->M_Kn = NULL; E->M_Kt = NULL;
  E->M_Un = NULL; E->M_Ut = NULL; E->M_Rn = NULL; E->M_Rt = NULL;
  E->M_Vn = NULL; E->M_Vt = NULL; 
  /* E->M_Jn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->J_bmr_b,D->A_bmc_b); ura_dup(E->M_Jn->mr_j,E->M_Jn->mr_length,E->J_bmr_j); ura_dup(E->M_Jn->mc_j,E->M_Jn->mc_length,D->A_bmc_j); M_mxget(E->M_Jn); */
  /* E->M_Jt = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,D->A_bmc_b,E->J_bmr_b); ura_dup(E->M_Jt->mr_j,E->M_Jt->mr_length,D->A_bmc_j); ura_dup(E->M_Jt->mc_j,E->M_Jt->mc_length,E->J_bmr_j); M_mxget(E->M_Jt); */
  /* E->M_Kn = M_handle_v_make(D->bitj,E->A_nrows,D->A_ncols,NULL,E->M_An->wX,E->A_bmr_b,D->K_bmc_b); ura_dup(E->M_Kn->mr_j,E->M_Kn->mr_length,E->A_bmr_j); ura_dup(E->M_Kn->mc_j,E->M_Kn->mc_length,D->K_bmc_j); M_mxget(E->M_Kn); */
  /* E->M_Kt = M_handle_v_make(D->bitj,D->A_ncols,E->A_nrows,NULL,E->M_At->wX,D->K_bmc_b,E->A_bmr_b); ura_dup(E->M_Kt->mr_j,E->M_Kt->mr_length,D->K_bmc_j); ura_dup(E->M_Kt->mc_j,E->M_Kt->mc_length,E->A_bmr_j); M_mxget(E->M_Kt); */
  /* E->M_Un = M_handle_v_make(D->bitj,E->A_nrows,D->Y_ncols,NULL,E->M_Yn->wX,E->J_bmr_b,D->Y_bmc_b); ura_dup(E->M_Un->mr_j,E->M_Un->mr_length,E->J_bmr_j); ura_dup(E->M_Un->mc_j,E->M_Un->mc_length,D->Y_bmc_j); M_mxget(E->M_Un); */
  /* E->M_Ut = M_handle_v_make(D->bitj,D->Y_ncols,E->A_nrows,NULL,E->M_Yt->wX,D->Y_bmc_b,E->J_bmr_b); ura_dup(E->M_Ut->mr_j,E->M_Ut->mr_length,D->Y_bmc_j); ura_dup(E->M_Ut->mc_j,E->M_Ut->mc_length,E->J_bmr_j); M_mxget(E->M_Ut); */
  /* E->M_Rn = M_handle_v_make(D->bitj,E->A_nrows,D->T_ncols,NULL,E->M_Tn->wX,E->J_bmr_b,D->T_bmc_b); ura_dup(E->M_Rn->mr_j,E->M_Rn->mr_length,E->J_bmr_j); ura_dup(E->M_Rn->mc_j,E->M_Rn->mc_length,D->T_bmc_j); M_mxget(E->M_Rn); */
  /* E->M_Rt = M_handle_v_make(D->bitj,D->T_ncols,E->A_nrows,NULL,E->M_Tt->wX,D->T_bmc_b,E->J_bmr_b); ura_dup(E->M_Rt->mr_j,E->M_Rt->mr_length,D->T_bmc_j); ura_dup(E->M_Rt->mc_j,E->M_Rt->mc_length,E->J_bmr_j); M_mxget(E->M_Rt); */
  /* E->M_Vn = M_handle_v_make(D->bitj,E->Z_nrows,D->A_ncols,NULL,E->M_Zn->wX,E->Z_bmr_b,D->K_bmc_b); ura_dup(E->M_Vn->mr_j,E->M_Vn->mr_length,E->Z_bmr_j); ura_dup(E->M_Vn->mc_j,E->M_Vn->mc_length,D->K_bmc_j); M_mxget(E->M_Vn); */
  /* E->M_Vt = M_handle_v_make(D->bitj,D->A_ncols,E->Z_nrows,NULL,E->M_Zt->wX,D->K_bmc_b,E->Z_bmr_b); ura_dup(E->M_Vt->mr_j,E->M_Vt->mr_length,D->K_bmc_j); ura_dup(E->M_Vt->mc_j,E->M_Vt->mc_length,E->Z_bmr_j); M_mxget(E->M_Vt); */
  bn = E->M_An->wX; bt = E->M_At->wX; E->M_Jn = NULL; E->M_Jt = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->A_ncols,E->J_bmr_b,D->A_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->A_ncols,bn,bt,E->J_bmr_b,E->J_bmr_j,D->A_bmc_b,D->A_bmc_j,&(E->M_Jn),&(E->M_Jt));
  bn  = E->M_An->wX; bt  = E->M_At->wX; E->M_Kn = NULL; E->M_Kt = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->A_ncols,E->A_bmr_b,D->K_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->A_ncols,bn,bt,E->A_bmr_b,E->A_bmr_j,D->K_bmc_b,D->K_bmc_j,&(E->M_Kn),&(E->M_Kt));
  bn  = E->M_Yn->wX; bt  = E->M_Yt->wX; E->M_Un = NULL; E->M_Ut = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->Y_ncols,E->J_bmr_b,D->Y_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->Y_ncols,bn,bt,E->J_bmr_b,E->J_bmr_j,D->Y_bmc_b,D->Y_bmc_j,&(E->M_Un),&(E->M_Ut));
  bn  = E->M_Zn->wX; bt  = E->M_Zt->wX; E->M_Vn = NULL; E->M_Vt = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->Z_nrows,D->A_ncols,E->Z_bmr_b,D->K_bmc_b,&(bn),&(bt));
  wrap_M_setup_test_excerpt_2(E->Z_nrows,D->A_ncols,bn,bt,E->Z_bmr_b,E->Z_bmr_j,D->K_bmc_b,D->K_bmc_j,&(E->M_Vn),&(E->M_Vt));
  bn  = E->M_Tn->wX; bt  = E->M_Tt->wX; E->M_Rn = NULL; E->M_Rt = NULL;
  //wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,E->A_nrows,D->T_ncols,E->J_bmr_b,D->T_bmc_b,&(bn),&(bt));
  for (nx=0;nx<E->A_rpop_b;nx++){ bn2 = &(bn[nx*bsize(D->T_ncols)]); bset__on(bn2,0);}
  for (nx=0;nx<E->A_nrows;nx++){ bt2 = &(bt[0*bsize(E->A_nrows)]); bset__on(bt2,nx);}
  wrap_M_setup_test_excerpt_2(E->A_nrows,D->T_ncols,bn,bt,E->J_bmr_b,E->J_bmr_j,D->T_bmc_b,D->T_bmc_j,&(E->M_Rn),&(E->M_Rt));
  /* initializing E->QR_RJnAtT_nrm */
  E->QR_RUnWtS_nrm = NULL; E->QR_RUnYtT_nrm = NULL; E->QR_RJnZtS_nrm = NULL; E->QR_RJnAtT_nrm = NULL; E->QR_J_index_global_mr_a=NULL; E->QR_J_index_global_mr_b=NULL;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose){ printf(" %% [finished dcg_single_init_M_An]\n");}
}

void dcg_single_copy_lf(struct dcg_single *E,struct dcg_single *E_in)
{
  int verbose=0; /* double gamma = GLOBAL_gamma; */
  struct dcg_ajdk *D = E->D;
  int ns_b=0,n1a=0,n1b=0,n2a=0,n2b=0;
  char prefix[FNAMESIZE];
  int M_flag=0;
  int cdrop=0,ckeep=0;
  long long int mm=0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
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
  L_handle_copy(E->lf_A_a0d1,E_in->lf_A_a0d1);
  L_handle_copy(E->lf_A_a2d1,E_in->lf_A_a2d1);
  L_handle_copy(E->lf_Y_a0d1,E_in->lf_Y_a0d1);
  L_handle_copy(E->lf_Y_a2d1,E_in->lf_Y_a2d1);
  L_handle_copy(E->lf_a1d1_WtSn,E_in->lf_a1d1_WtSn);
  L_handle_copy(E->lf_a1d1_ZtSn,E_in->lf_a1d1_ZtSn);
  L_handle_copy(E->lf_a1d1_YtTn,E_in->lf_a1d1_YtTn);
  L_handle_copy(E->lf_a1d1_AtTn,E_in->lf_a1d1_AtTn);
  L_handle_copy(E->lf_et_Sn,E_in->lf_et_Sn);
  L_handle_copy(E->lf_et_Tn,E_in->lf_et_Tn);
  L_handle_copy(E->lf_a0d1_WtSn,E_in->lf_a0d1_WtSn); M_handle_copy(E->M_a0d1_WtSn,E_in->M_a0d1_WtSn);
  L_handle_copy(E->lf_a0d1_ZtSn,E_in->lf_a0d1_ZtSn); M_handle_copy(E->M_a0d1_ZtSn,E_in->M_a0d1_ZtSn);
  L_handle_copy(E->lf_a0d1_YtTn,E_in->lf_a0d1_YtTn); M_handle_copy(E->M_a0d1_YtTn,E_in->M_a0d1_YtTn);
  L_handle_copy(E->lf_a0d1_AtTn,E_in->lf_a0d1_AtTn); M_handle_copy(E->M_a0d1_AtTn,E_in->M_a0d1_AtTn);
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  L_handle_copy(E->lf_Jt_rsum,E_in->lf_Jt_rsum);
  L_handle_copy(E->lf_Kt_rsum,E_in->lf_Kt_rsum);
  L_handle_copy(E->lf_Wt_rsum,E_in->lf_Wt_rsum);
  L_handle_copy(E->lf_Ut_rsum,E_in->lf_Ut_rsum);
  L_handle_copy(E->lf_Vt_rsum,E_in->lf_Vt_rsum);
  L_handle_copy(E->lf_JtRn,E_in->lf_JtRn);
  L_handle_copy(E->lf_KtTn,E_in->lf_KtTn);
  L_handle_copy(E->lf_UtRn,E_in->lf_UtRn);
  L_handle_copy(E->lf_VtSn,E_in->lf_VtSn);
  L_handle_copy(E->lf_Jn_ajdk,E_in->lf_Jn_ajdk);
  L_handle_copy(E->lf_Kn_ajdk,E_in->lf_Kn_ajdk);
  L_handle_copy(E->lf_Un_ajdk,E_in->lf_Un_ajdk);
  L_handle_copy(E->lf_Vn_ajdk,E_in->lf_Vn_ajdk);
  L_handle_copy(E->lf_J_a0d1,E_in->lf_J_a0d1);
  L_handle_copy(E->lf_J_a2d1,E_in->lf_J_a2d1);
  L_handle_copy(E->lf_K_a0d1,E_in->lf_J_a0d1);
  L_handle_copy(E->lf_K_a2d1,E_in->lf_J_a2d1);
  L_handle_copy(E->lf_U_a0d1,E_in->lf_U_a0d1);
  L_handle_copy(E->lf_U_a2d1,E_in->lf_U_a2d1);
  L_handle_copy(E->lf_V_a0d1,E_in->lf_V_a0d1);
  L_handle_copy(E->lf_V_a2d1,E_in->lf_V_a2d1);
  L_handle_copy(E->lf_a1d1_JtRn,E_in->lf_a1d1_JtRn);
  L_handle_copy(E->lf_a1d1_KtTn,E_in->lf_a1d1_KtTn);
  L_handle_copy(E->lf_a1d1_UtRn,E_in->lf_a1d1_UtRn);
  L_handle_copy(E->lf_a1d1_VtSn,E_in->lf_a1d1_VtSn);
  L_handle_copy(E->lf_et_Rn,E_in->lf_et_Rn);
  L_handle_copy(E->lf_a0d1_JtRn,E_in->lf_a0d1_JtRn); M_handle_copy(E->M_a0d1_JtRn,E_in->M_a0d1_JtRn);
  L_handle_copy(E->lf_a0d1_KtTn,E_in->lf_a0d1_KtTn); M_handle_copy(E->M_a0d1_KtTn,E_in->M_a0d1_KtTn);
  L_handle_copy(E->lf_a0d1_UtRn,E_in->lf_a0d1_UtRn); M_handle_copy(E->M_a0d1_UtRn,E_in->M_a0d1_UtRn);
  L_handle_copy(E->lf_a0d1_VtSn,E_in->lf_a0d1_VtSn); M_handle_copy(E->M_a0d1_VtSn,E_in->M_a0d1_VtSn);
}

void dcg_single_init_lf_excerpt(char *prefix,int n1a,int n1b,int n2a,int n2b,struct L_handle **L_p,int M_flag,struct M_handle **M_p)
{
  int verbose=0;
  unsigned long long int length=0;
  length  = (unsigned long long int)n1b*(unsigned long long int)n2b; *L_p = L_handle_make(length); if (!(*L_p)){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  if (M_flag){ *M_p = M_handle_w_make(BITJ,GLOBAL_B_MLT,n1a,n2b);}
  if (verbose){ printf(" %% %s n1a %d n1b %d n2a %d n2b %d --> length %llu \n",prefix,n1a,n1b,n2a,n2b,length);}
}

void dcg_single_init_lf(struct dcg_single *E)
{
  int verbose=2; /* double gamma = GLOBAL_gamma; */
  struct dcg_ajdk *D = E->D;
  int ns_b=0,n1a=0,n1b=0,n2a=0,n2b=0;
  char prefix[FNAMESIZE];
  int M_flag=0;
  int cdrop=0,ckeep=0;
  long long int mm=0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% initializing E->lf_At_rsum, etc.\n");}
  E->lf_At_rsum=NULL; E->lf_Zt_rsum=NULL; E->lf_Wt_rsum=NULL; E->lf_Yt_rsum=NULL;
  E->length = D->A_ncols; E->lf_At_rsum = L_handle_make(E->length); if (!E->lf_At_rsum){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->A_ncols; E->lf_Zt_rsum = L_handle_make(E->length); if (!E->lf_Zt_rsum){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->Y_ncols; E->lf_Wt_rsum = L_handle_make(E->length); if (!E->lf_Wt_rsum){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->Y_ncols; E->lf_Yt_rsum = L_handle_make(E->length); if (!E->lf_Yt_rsum){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing E->lf_AtTn, E->lf_ZtSn, E->lf_YtTn, E->lf_WtSn.\n");}
  E->lf_AtTn = NULL; E->lf_ZtSn = NULL; E->lf_YtTn = NULL; E->lf_WtSn = NULL;
  E->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; E->lf_AtTn = L_handle_make(E->length); if (!E->lf_AtTn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; E->lf_ZtSn = L_handle_make(E->length); if (!E->lf_ZtSn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; E->lf_YtTn = L_handle_make(E->length); if (!E->lf_YtTn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; E->lf_WtSn = L_handle_make(E->length); if (!E->lf_WtSn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing E->lf_An_ajdk, etc\n");}
  E->lf_An_ajdk = NULL; E->lf_Zn_ajdk = NULL; E->lf_Yn_ajdk = NULL; E->lf_Wn_ajdk = NULL;
  E->length = (unsigned long long int)E->A_nrows*(unsigned long long int)AJDK_TOT; E->lf_An_ajdk = L_handle_make(E->length); if (!E->lf_An_ajdk){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)E->Z_nrows*(unsigned long long int)AJDK_TOT; E->lf_Zn_ajdk = L_handle_make(E->length); if (!E->lf_Zn_ajdk){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)E->A_nrows*(unsigned long long int)AJDK_TOT; E->lf_Yn_ajdk = L_handle_make(E->length); if (!E->lf_Yn_ajdk){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)E->Z_nrows*(unsigned long long int)AJDK_TOT; E->lf_Wn_ajdk = L_handle_make(E->length); if (!E->lf_Wn_ajdk){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing E->lf_A_a0d1, etc\n");}
  E->length = 1; E->lf_A_a0d1 = L_handle_make(E->length);
  E->length = 1; E->lf_A_a2d1 = L_handle_make(E->length);
  E->length = 1; E->lf_Y_a0d1 = L_handle_make(E->length);
  E->length = 1; E->lf_Y_a2d1 = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_a1d1_WtSn = L_handle_make(E->length); if (!E->lf_a1d1_WtSn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_a1d1_ZtSn = L_handle_make(E->length); if (!E->lf_a1d1_ZtSn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_a1d1_YtTn = L_handle_make(E->length); if (!E->lf_a1d1_YtTn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_a1d1_AtTn = L_handle_make(E->length); if (!E->lf_a1d1_AtTn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_et_Sn = L_handle_make(E->length); if (!E->lf_et_Sn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_et_Tn = L_handle_make(E->length); if (!E->lf_et_Tn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing storage for doublestudy\n");}
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_WtSn");
  n1a = D->Y_cbother*D->Y_ncols; n1b = D->Y_cbother*D->Y_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcg_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_WtSn),M_flag,&(E->M_a0d1_WtSn));
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_ZtSn");
  n1a = D->A_cbother*D->A_ncols; n1b = D->A_cbother*D->A_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcg_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_ZtSn),M_flag,&(E->M_a0d1_ZtSn));
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_YtTn");
  n1a = D->Y_cbother*D->Y_ncols; n1b = D->Y_cbother*D->Y_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcg_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_YtTn),M_flag,&(E->M_a0d1_YtTn));
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_AtTn");
  n1a = D->A_cbother*D->A_ncols; n1b = D->A_cbother*D->A_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcg_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_AtTn),M_flag,&(E->M_a0d1_AtTn));
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% initializing E->lf_Jt_rsum, etc.\n");}
  E->lf_Jt_rsum=NULL; E->lf_Kt_rsum=NULL; E->lf_Wt_rsum=NULL; E->lf_Yt_rsum=NULL;
  E->length = D->A_ncols; E->lf_Jt_rsum = L_handle_make(E->length); if (!E->lf_Jt_rsum){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->A_ncols; E->lf_Kt_rsum = L_handle_make(E->length); if (!E->lf_Kt_rsum){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->A_ncols; E->lf_Vt_rsum = L_handle_make(E->length); if (!E->lf_Vt_rsum){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->Y_ncols; E->lf_Ut_rsum = L_handle_make(E->length); if (!E->lf_Ut_rsum){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  if (verbose){ printf(" %% initializing E->lf_JtRn, E->lf_KtTn, E->lf_UtRn, E->lf_VtSn.\n");}
  E->lf_JtRn = NULL; E->lf_KtTn = NULL; E->lf_VtSn = NULL; E->lf_VtSn = NULL;
  E->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; E->lf_JtRn = L_handle_make(E->length); if (!E->lf_JtRn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; E->lf_KtTn = L_handle_make(E->length); if (!E->lf_KtTn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; E->lf_VtSn = L_handle_make(E->length); if (!E->lf_VtSn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; E->lf_UtRn = L_handle_make(E->length); if (!E->lf_UtRn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing E->lf_An_ajdk, etc\n");}
  E->lf_Jn_ajdk = NULL; E->lf_Kn_ajdk = NULL; E->lf_Un_ajdk = NULL; E->lf_Vn_ajdk = NULL;
  E->length = (unsigned long long int)E->A_nrows*(unsigned long long int)AJDK_TOT; E->lf_Jn_ajdk = L_handle_make(E->length); if (!E->lf_Jn_ajdk){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)E->A_nrows*(unsigned long long int)AJDK_TOT; E->lf_Kn_ajdk = L_handle_make(E->length); if (!E->lf_Kn_ajdk){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)E->Z_nrows*(unsigned long long int)AJDK_TOT; E->lf_Vn_ajdk = L_handle_make(E->length); if (!E->lf_Vn_ajdk){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = (unsigned long long int)E->A_nrows*(unsigned long long int)AJDK_TOT; E->lf_Un_ajdk = L_handle_make(E->length); if (!E->lf_Un_ajdk){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing E->lf_A_a0d1, etc\n");}
  E->length = 1; E->lf_J_a0d1 = L_handle_make(E->length);
  E->length = 1; E->lf_J_a2d1 = L_handle_make(E->length);
  E->length = 1; E->lf_K_a0d1 = L_handle_make(E->length);
  E->length = 1; E->lf_K_a2d1 = L_handle_make(E->length);
  E->length = 1; E->lf_U_a0d1 = L_handle_make(E->length);
  E->length = 1; E->lf_U_a2d1 = L_handle_make(E->length);
  E->length = 1; E->lf_V_a0d1 = L_handle_make(E->length);
  E->length = 1; E->lf_V_a2d1 = L_handle_make(E->length);
  E->length = D->T_ncols; E->lf_a1d1_JtRn = L_handle_make(E->length); if (!E->lf_a1d1_JtRn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_a1d1_KtTn = L_handle_make(E->length); if (!E->lf_a1d1_KtTn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_a1d1_VtSn = L_handle_make(E->length); if (!E->lf_a1d1_VtSn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_a1d1_UtRn = L_handle_make(E->length); if (!E->lf_a1d1_UtRn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  E->length = D->T_ncols; E->lf_et_Rn = L_handle_make(E->length); if (!E->lf_et_Rn){ printf(" %% Warning! not enough memory in dcg_single_init_lf\n");}
  if (verbose>1){ printf(" %% initializing storage for doublestudy\n");}
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_VtSn");
  n1a = D->A_cbother*D->A_ncols; n1b = D->A_cbother*D->A_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcg_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_VtSn),M_flag,&(E->M_a0d1_VtSn));
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_UtRn");
  n1a = D->Y_cbother*D->Y_ncols; n1b = D->Y_cbother*D->Y_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcg_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_UtRn),M_flag,&(E->M_a0d1_UtRn));
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_JtRn");
  n1a = D->A_cbother*D->A_ncols; n1b = D->A_cbother*D->A_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcg_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_JtRn),M_flag,&(E->M_a0d1_JtRn));
  M_flag = 1;
  sprintf(prefix,"E->M_a0d1_KtTn");
  n1a = D->A_cbother*D->A_ncols; n1b = D->A_cbother*D->A_cpop_b;
  n2a = D->T_ncols; n2b = D->T_cpop_b; 
  dcg_single_init_lf_excerpt(prefix,n1a,n1b,n2a,n2b,&(E->lf_a0d1_KtTn),M_flag,&(E->M_a0d1_KtTn));
}

void dcg_double_copy_lf(struct dcg_double *F,struct dcg_double *F_in)
{
  int verbose=0;
  struct dcg_ajdk *D=F->D;
  int n1a=0,n2a=0,n1b,n2b=0;
  int M_flag=0,rdrop=0,rdrop2=0;
  char prefix[FNAMESIZE];
  GLOBAL_tic(0);
  if (verbose>1){ printf(" %% dcg_double_copy_lf: %d,%d\n",F->nb1,F->nb2);}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  L_handle_copy(F->lf_Yn_a0d1_WtSn,F_in->lf_Yn_a0d1_WtSn);
  L_handle_copy(F->lf_An_a0d1_ZtSn,F_in->lf_An_a0d1_ZtSn);
  L_handle_copy(F->lf_Yn_a0d1_YtTn,F_in->lf_Yn_a0d1_YtTn);
  L_handle_copy(F->lf_An_a0d1_AtTn,F_in->lf_An_a0d1_AtTn);
  L_handle_copy(F->lf_TYnWtS_ww,F_in->lf_TYnWtS_ww);
  L_handle_copy(F->lf_TAnZtS_ww,F_in->lf_TAnZtS_ww);
  L_handle_copy(F->lf_TYnYtT_ww,F_in->lf_TYnYtT_ww);
  L_handle_copy(F->lf_TAnAtT_ww,F_in->lf_TAnAtT_ww);
  L_handle_copy(F->lf_TYnWtS_uu,F_in->lf_TYnWtS_uu);
  L_handle_copy(F->lf_TAnZtS_uu,F_in->lf_TAnZtS_uu);
  L_handle_copy(F->lf_TYnYtT_uu,F_in->lf_TYnYtT_uu);
  L_handle_copy(F->lf_TAnAtT_uu,F_in->lf_TAnAtT_uu);
  L_handle_copy(F->lf_D_YtTn_WtSn_vv,F_in->lf_D_YtTn_WtSn_vv);
  L_handle_copy(F->lf_D_AtTn_ZtSn_vv,F_in->lf_D_AtTn_ZtSn_vv);
  L_handle_copy(F->lf_D_YtTn_YtTn_vv,F_in->lf_D_YtTn_YtTn_vv);
  L_handle_copy(F->lf_D_AtTn_AtTn_vv,F_in->lf_D_AtTn_AtTn_vv);
  L_handle_copy(F->lf_D_YtTn_WtSn_uu,F_in->lf_D_YtTn_WtSn_uu);
  L_handle_copy(F->lf_D_AtTn_ZtSn_uu,F_in->lf_D_AtTn_ZtSn_uu);
  L_handle_copy(F->lf_D_YtTn_YtTn_uu,F_in->lf_D_YtTn_YtTn_uu);
  L_handle_copy(F->lf_D_AtTn_AtTn_uu,F_in->lf_D_AtTn_AtTn_uu);
  L_handle_copy(F->QR_TYnWtS,F_in->QR_TYnWtS);
  L_handle_copy(F->QR_TAnZtS,F_in->QR_TAnZtS);
  L_handle_copy(F->QR_TYnYtT,F_in->QR_TYnYtT);
  L_handle_copy(F->QR_TAnAtT,F_in->QR_TAnAtT);
  L_handle_copy(F->QC_TYnWtS,F_in->QC_TYnWtS);
  L_handle_copy(F->QC_TAnZtS,F_in->QC_TAnZtS);
  L_handle_copy(F->QC_TYnYtT,F_in->QC_TYnYtT);
  L_handle_copy(F->QC_TAnAtT,F_in->QC_TAnAtT);
  L_handle_copy(F->QR_TYnWtS_uu,F_in->QR_TYnWtS_uu);
  L_handle_copy(F->QR_TAnZtS_uu,F_in->QR_TAnZtS_uu);
  L_handle_copy(F->QR_TYnYtT_uu,F_in->QR_TYnYtT_uu);
  L_handle_copy(F->QR_TAnAtT_uu,F_in->QR_TAnAtT_uu);
  L_handle_copy(F->QC_TYnWtS_uu,F_in->QC_TYnWtS_uu);
  L_handle_copy(F->QC_TAnZtS_uu,F_in->QC_TAnZtS_uu);
  L_handle_copy(F->QC_TYnYtT_uu,F_in->QC_TYnYtT_uu);
  L_handle_copy(F->QC_TAnAtT_uu,F_in->QC_TAnAtT_uu);
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  L_handle_copy(F->lf_Un_a0d1_WtSn,F_in->lf_Un_a0d1_WtSn);
  L_handle_copy(F->lf_Un_a0d1_YtTn,F_in->lf_Un_a0d1_YtTn);
  L_handle_copy(F->lf_Jn_a0d1_ZtSn,F_in->lf_Jn_a0d1_ZtSn);
  L_handle_copy(F->lf_Jn_a0d1_AtTn,F_in->lf_Jn_a0d1_AtTn);
  L_handle_copy(F->lf_RUnWtS_ww,F_in->lf_RUnWtS_ww);
  L_handle_copy(F->lf_RUnYtT_ww,F_in->lf_RUnYtT_ww);
  L_handle_copy(F->lf_RJnZtS_ww,F_in->lf_RJnZtS_ww);
  L_handle_copy(F->lf_RJnAtT_ww,F_in->lf_RJnAtT_ww);
  L_handle_copy(F->lf_RUnWtS_uu,F_in->lf_RUnWtS_uu);
  L_handle_copy(F->lf_RUnYtT_uu,F_in->lf_RUnYtT_uu);
  L_handle_copy(F->lf_RJnZtS_uu,F_in->lf_RJnZtS_uu);
  L_handle_copy(F->lf_RJnAtT_uu,F_in->lf_RJnAtT_uu);
  L_handle_copy(F->lf_D_UtRn_WtSn_vv,F_in->lf_D_UtRn_WtSn_vv);
  L_handle_copy(F->lf_D_UtRn_YtTn_vv,F_in->lf_D_UtRn_YtTn_vv);
  L_handle_copy(F->lf_D_JtRn_ZtSn_vv,F_in->lf_D_JtRn_ZtSn_vv);
  L_handle_copy(F->lf_D_JtRn_AtTn_vv,F_in->lf_D_JtRn_AtTn_vv);
  L_handle_copy(F->lf_D_KtTn_VtSn_vv,F_in->lf_D_KtTn_VtSn_vv);
  L_handle_copy(F->lf_D_KtTn_KtTn_vv,F_in->lf_D_KtTn_KtTn_vv);
  L_handle_copy(F->lf_D_UtRn_WtSn_uu,F_in->lf_D_UtRn_WtSn_uu);
  L_handle_copy(F->lf_D_UtRn_YtTn_uu,F_in->lf_D_UtRn_YtTn_uu);
  L_handle_copy(F->lf_D_JtRn_ZtSn_uu,F_in->lf_D_JtRn_ZtSn_uu);
  L_handle_copy(F->lf_D_JtRn_AtTn_uu,F_in->lf_D_JtRn_AtTn_uu);
  L_handle_copy(F->lf_D_KtTn_VtSn_uu,F_in->lf_D_KtTn_VtSn_uu);
  L_handle_copy(F->lf_D_KtTn_KtTn_uu,F_in->lf_D_KtTn_KtTn_uu);
  L_handle_copy(F->QR_RUnWtS,F_in->QR_RUnWtS);
  L_handle_copy(F->QR_RUnYtT,F_in->QR_RUnYtT);
  L_handle_copy(F->QR_RJnZtS,F_in->QR_RJnZtS);
  L_handle_copy(F->QR_RJnAtT,F_in->QR_RJnAtT);
  L_handle_copy(F->QR_RUnWtS_uu,F_in->QR_RUnWtS_uu);
  L_handle_copy(F->QR_RUnYtT_uu,F_in->QR_RUnYtT_uu);
  L_handle_copy(F->QR_RJnZtS_uu,F_in->QR_RJnZtS_uu);
  L_handle_copy(F->QR_RJnAtT_uu,F_in->QR_RJnAtT_uu);
  L_handle_copy(F->QC_RUnWtS,F_in->QC_RUnWtS);
  L_handle_copy(F->QC_RUnYtT,F_in->QC_RUnYtT);
  L_handle_copy(F->QC_RJnZtS,F_in->QC_RJnZtS);
  L_handle_copy(F->QC_RJnAtT,F_in->QC_RJnAtT);
  L_handle_copy(F->QC_TKnVtS,F_in->QC_TKnVtS);
  L_handle_copy(F->QC_TKnKtT,F_in->QC_TKnKtT);
  L_handle_copy(F->QC_RUnWtS_uu,F_in->QC_RUnWtS_uu);
  L_handle_copy(F->QC_RUnYtT_uu,F_in->QC_RUnYtT_uu);
  L_handle_copy(F->QC_RJnZtS_uu,F_in->QC_RJnZtS_uu);
  L_handle_copy(F->QC_RJnAtT_uu,F_in->QC_RJnAtT_uu);
  L_handle_copy(F->QC_TKnVtS_uu,F_in->QC_TKnVtS_uu);
  L_handle_copy(F->QC_TKnKtT_uu,F_in->QC_TKnKtT_uu);
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  GLOBAL_toc(0,verbose," %% copying data_structures: ");
}

void dcg_double_init_lf(struct dcg_double *F)
{
  int verbose=0;
  struct dcg_ajdk *D=F->D;
  int n1a=0,n2a=0,n1b,n2b=0;
  int M_flag=0,rdrop=0,rdrop2=0;
  char prefix[FNAMESIZE];
  GLOBAL_tic(0);
  if (verbose>1){ printf(" %% dcg_double_init_lf: %d,%d\n",F->nb1,F->nb2);}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% initializing lf_ \n");}
  F->length = (unsigned long long int)F->E_nb1->M_Yn->nrows*(unsigned long long int)D->T_ncols; F->lf_Yn_a0d1_WtSn = L_handle_make(F->length); if (!F->lf_Yn_a0d1_WtSn){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_An->nrows*(unsigned long long int)D->T_ncols; F->lf_An_a0d1_ZtSn = L_handle_make(F->length); if (!F->lf_An_a0d1_ZtSn){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Yn->nrows*(unsigned long long int)D->T_ncols; F->lf_Yn_a0d1_YtTn = L_handle_make(F->length); if (!F->lf_Yn_a0d1_YtTn){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_An->nrows*(unsigned long long int)D->T_ncols; F->lf_An_a0d1_AtTn = L_handle_make(F->length); if (!F->lf_An_a0d1_AtTn){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Yn->nrows*(unsigned long long int)D->T_ncols; F->lf_TYnWtS_ww = L_handle_make(F->length); if (!F->lf_TYnWtS_ww){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_An->nrows*(unsigned long long int)D->T_ncols; F->lf_TAnZtS_ww = L_handle_make(F->length); if (!F->lf_TAnZtS_ww){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Yn->nrows*(unsigned long long int)D->T_ncols; F->lf_TYnYtT_ww = L_handle_make(F->length); if (!F->lf_TYnYtT_ww){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_An->nrows*(unsigned long long int)D->T_ncols; F->lf_TAnAtT_ww = L_handle_make(F->length); if (!F->lf_TAnAtT_ww){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Yn->nrows*(unsigned long long int)D->T_ncols; F->lf_TYnWtS_uu = L_handle_make(F->length); if (!F->lf_TYnWtS_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_An->nrows*(unsigned long long int)D->T_ncols; F->lf_TAnZtS_uu = L_handle_make(F->length); if (!F->lf_TAnZtS_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Yn->nrows*(unsigned long long int)D->T_ncols; F->lf_TYnYtT_uu = L_handle_make(F->length); if (!F->lf_TYnYtT_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_An->nrows*(unsigned long long int)D->T_ncols; F->lf_TAnAtT_uu = L_handle_make(F->length); if (!F->lf_TAnAtT_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->lf_D_YtTn_WtSn_vv = L_handle_make(F->length); if (!F->lf_D_YtTn_WtSn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_AtTn_ZtSn_vv = L_handle_make(F->length); if (!F->lf_D_AtTn_ZtSn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->lf_D_YtTn_YtTn_vv = L_handle_make(F->length); if (!F->lf_D_YtTn_YtTn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_AtTn_AtTn_vv = L_handle_make(F->length); if (!F->lf_D_AtTn_AtTn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->lf_D_YtTn_WtSn_uu = L_handle_make(F->length); if (!F->lf_D_YtTn_WtSn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_AtTn_ZtSn_uu = L_handle_make(F->length); if (!F->lf_D_AtTn_ZtSn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->lf_D_YtTn_YtTn_uu = L_handle_make(F->length); if (!F->lf_D_YtTn_YtTn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_AtTn_AtTn_uu = L_handle_make(F->length); if (!F->lf_D_AtTn_AtTn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  if (verbose>1){ printf(" %% initializing F->QX\n");}
  F->QR_TYnWtS=NULL; F->QR_TYnYtT=NULL; F->QC_TYnWtS=NULL; F->QC_TYnYtT=NULL;
  F->QR_TAnZtS=NULL; F->QR_TAnAtT=NULL; F->QC_TAnZtS=NULL; F->QC_TAnAtT=NULL;
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_TYnWtS = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_TAnZtS = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_TYnYtT = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_TAnAtT = L_handle_make(F->length);
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->QC_TYnWtS = L_handle_make(F->length);
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_TAnZtS = L_handle_make(F->length);
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->QC_TYnYtT = L_handle_make(F->length); if (!F->QC_TYnYtT){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_TAnAtT = L_handle_make(F->length); if (!F->QC_TAnAtT){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  if (verbose>1){ printf(" %% initializing F->QX_uu\n");}
  F->QR_TYnWtS_uu=NULL; F->QR_TYnYtT_uu=NULL; F->QC_TYnWtS_uu=NULL; F->QC_TYnYtT_uu=NULL;
  F->QR_TAnZtS_uu=NULL; F->QR_TAnAtT_uu=NULL; F->QC_TAnZtS_uu=NULL; F->QC_TAnAtT_uu=NULL;
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_TYnWtS_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_TAnZtS_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_TYnYtT_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_TAnAtT_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->QC_TYnWtS_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_TAnZtS_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->QC_TYnYtT_uu = L_handle_make(F->length); if (!F->QC_TYnYtT_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_TAnAtT_uu = L_handle_make(F->length); if (!F->QC_TAnAtT_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% initializing lf_ \n");}
  F->length = (unsigned long long int)F->E_nb1->M_Un->nrows*(unsigned long long int)D->T_ncols; F->lf_Un_a0d1_WtSn = L_handle_make(F->length); if (!F->lf_Un_a0d1_WtSn){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Un->nrows*(unsigned long long int)D->T_ncols; F->lf_Un_a0d1_YtTn = L_handle_make(F->length); if (!F->lf_Un_a0d1_YtTn){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Jn->nrows*(unsigned long long int)D->T_ncols; F->lf_Jn_a0d1_ZtSn = L_handle_make(F->length); if (!F->lf_Jn_a0d1_ZtSn){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Jn->nrows*(unsigned long long int)D->T_ncols; F->lf_Jn_a0d1_AtTn = L_handle_make(F->length); if (!F->lf_Jn_a0d1_AtTn){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Un->nrows*(unsigned long long int)D->T_ncols; F->lf_RUnWtS_ww = L_handle_make(F->length); if (!F->lf_RUnWtS_ww){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Un->nrows*(unsigned long long int)D->T_ncols; F->lf_RUnYtT_ww = L_handle_make(F->length); if (!F->lf_RUnYtT_ww){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Jn->nrows*(unsigned long long int)D->T_ncols; F->lf_RJnZtS_ww = L_handle_make(F->length); if (!F->lf_RJnZtS_ww){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Jn->nrows*(unsigned long long int)D->T_ncols; F->lf_RJnAtT_ww = L_handle_make(F->length); if (!F->lf_RJnAtT_ww){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Un->nrows*(unsigned long long int)D->T_ncols; F->lf_RUnWtS_uu = L_handle_make(F->length); if (!F->lf_RUnWtS_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Un->nrows*(unsigned long long int)D->T_ncols; F->lf_RUnYtT_uu = L_handle_make(F->length); if (!F->lf_RUnYtT_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Jn->nrows*(unsigned long long int)D->T_ncols; F->lf_RJnZtS_uu = L_handle_make(F->length); if (!F->lf_RJnZtS_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)F->E_nb1->M_Jn->nrows*(unsigned long long int)D->T_ncols; F->lf_RJnAtT_uu = L_handle_make(F->length); if (!F->lf_RJnAtT_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->lf_D_UtRn_WtSn_vv = L_handle_make(F->length); if (!F->lf_D_UtRn_WtSn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->lf_D_UtRn_YtTn_vv = L_handle_make(F->length); if (!F->lf_D_UtRn_YtTn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_JtRn_ZtSn_vv = L_handle_make(F->length); if (!F->lf_D_JtRn_ZtSn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_JtRn_AtTn_vv = L_handle_make(F->length); if (!F->lf_D_JtRn_AtTn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_KtTn_VtSn_vv = L_handle_make(F->length); if (!F->lf_D_KtTn_VtSn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_KtTn_KtTn_vv = L_handle_make(F->length); if (!F->lf_D_KtTn_KtTn_vv){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->lf_D_UtRn_WtSn_uu = L_handle_make(F->length); if (!F->lf_D_UtRn_WtSn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->lf_D_UtRn_YtTn_uu = L_handle_make(F->length); if (!F->lf_D_UtRn_YtTn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_JtRn_ZtSn_uu = L_handle_make(F->length); if (!F->lf_D_JtRn_ZtSn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_JtRn_AtTn_uu = L_handle_make(F->length); if (!F->lf_D_JtRn_AtTn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_KtTn_VtSn_uu = L_handle_make(F->length); if (!F->lf_D_KtTn_VtSn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->lf_D_KtTn_KtTn_uu = L_handle_make(F->length); if (!F->lf_D_KtTn_KtTn_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  if (verbose>1){ printf(" %% initializing F->QX\n");}
  F->QR_RUnWtS=NULL; F->QR_RUnYtT=NULL; F->QC_RUnWtS=NULL; F->QC_RUnYtT=NULL;
  F->QR_RJnZtS=NULL; F->QR_RJnAtT=NULL; F->QC_RJnZtS=NULL; F->QC_RJnAtT=NULL;
  F->QC_TKnVtS=NULL; F->QC_TKnKtT=NULL;
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_RUnWtS = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_RUnYtT = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_RJnZtS = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_RJnAtT = L_handle_make(F->length);
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->QC_RUnWtS = L_handle_make(F->length);
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->QC_RUnYtT = L_handle_make(F->length); if (!F->QC_RUnYtT){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_RJnZtS = L_handle_make(F->length);
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_RJnAtT = L_handle_make(F->length); if (!F->QC_RJnAtT){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_TKnVtS = L_handle_make(F->length);
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_TKnKtT = L_handle_make(F->length); if (!F->QC_TKnKtT){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  if (verbose>1){ printf(" %% initializing F->QX_uu\n");}
  F->QR_RUnWtS_uu=NULL; F->QR_RUnYtT_uu=NULL; F->QC_RUnWtS_uu=NULL; F->QC_RUnYtT_uu=NULL;
  F->QR_RJnZtS_uu=NULL; F->QR_RJnAtT_uu=NULL; F->QC_RJnZtS_uu=NULL; F->QC_RJnAtT_uu=NULL;
  F->QC_TKnVtS_uu=NULL; F->QC_TKnKtT_uu=NULL;
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_RUnWtS_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_RUnYtT_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_RJnZtS_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)F->E_nb1->A_nrows*(unsigned long long int)D->T_ncols; F->QR_RJnAtT_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->QC_RUnWtS_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)D->Y_ncols*(unsigned long long int)D->T_ncols; F->QC_RUnYtT_uu = L_handle_make(F->length); if (!F->QC_RUnYtT_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_RJnZtS_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_RJnAtT_uu = L_handle_make(F->length); if (!F->QC_RJnAtT_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_TKnVtS_uu = L_handle_make(F->length);
  F->length = (unsigned long long int)D->A_ncols*(unsigned long long int)D->T_ncols; F->QC_TKnKtT_uu = L_handle_make(F->length); if (!F->QC_TKnKtT_uu){ printf(" %% Warning! not enough memory in dcg_double_init_lf\n");}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  GLOBAL_toc(0,verbose," %% initializing data_structures: ");
}

void dcg_copy_A_p(struct dcg_ajdk *D,struct dcg_ajdk *D_in)
{
  /* copies A_p and Y_p for each column from D_in */
  int verbose=0; char tempchar[FNAMESIZE];
  int nbins = D->nbins; struct dcg_single **E_ = D->E_; struct dcg_single **E_in_ = D_in->E_;
  int nb=0,nc=0,nc_p=0,nc_start=0,nc_final=0;
  int A_pcols=0,Y_pcols=0,A_ncols=0,Y_ncols=0,A_rpop_b_total=0,Z_rpop_b_total=0;
  double tmp_sum=0,tmp_num=0;
  if (verbose){ printf(" %% [entering dcg_copy_A_p]\n"); wkspace_printf();}
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
  if (verbose){ printf(" %% [finished dcg_copy_A_p]\n");  wkspace_printf();}
}

void dcg_X_nrows_total(struct dcg_ajdk *D)
{
  int verbose=0;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_;
  int nb=0,nc=0;
  if (verbose){ printf(" %% [entering dcg_X_nrows_total]\n");}
  D->A_nrows_total=0; D->A_rpop_b_total=0; D->A_rpop_j_total=0;
  D->Z_nrows_total=0; D->Z_rpop_b_total=0; D->Z_rpop_j_total=0;
  D->J_rpop_b_total=0; D->J_rpop_j_total=0;
  for (nb=0;nb<nbins;nb++){
    D->A_nrows_total += E_[nb]->A_nrows; D->A_rpop_b_total += E_[nb]->A_rpop_b; D->A_rpop_j_total += E_[nb]->A_rpop_j;
    D->Z_nrows_total += E_[nb]->Z_nrows; D->Z_rpop_b_total += E_[nb]->Z_rpop_b; D->Z_rpop_j_total += E_[nb]->Z_rpop_j;
    D->J_rpop_b_total += E_[nb]->J_rpop_b; D->J_rpop_j_total += E_[nb]->J_rpop_j;
    /* for (nb=0;nb<nbins;nb++){ }*/}
  if (verbose){ printf(" %% D->A_nrows_total %d D->Z_nrows_total %d\n",D->A_nrows_total,D->Z_nrows_total);}
  if (verbose){ printf(" %% D->A_rpop_b_total %d D->Z_rpop_b_total %d\n",D->A_rpop_b_total,D->Z_rpop_b_total);}
  if (verbose){ printf(" %% D->A_rpop_j_total %d D->Z_rpop_j_total %d\n",D->A_rpop_j_total,D->Z_rpop_j_total);}
  if (verbose){ printf(" %% D->J_rpop_b_total %d D->J_rpop_j_total %d\n",D->J_rpop_b_total,D->J_rpop_j_total);}
  if (verbose){ printf(" %% [finished dcg_X_nrows_total]\n");}
}

void dcg_load_A_p(struct dcg_ajdk *D)
{
  /* calculates A_p and Y_p for each column using row-mask mr_b, not mr_j */
  int verbose=0; char tempchar[FNAMESIZE];
  int nbins = D->nbins; struct dcg_single **E_ = D->E_;
  int nb=0,nc=0,nc_p=0,nc_start=0,nc_final=0;
  int A_pcols=0,Y_pcols=0,A_ncols=0,Y_ncols=0,A_rpop_b_total=0,Z_rpop_b_total=0;
  double *AZ_rsum=NULL,*A_p=NULL,**A_ajdk_p=NULL,*YW_rsum=NULL,*Y_p=NULL,**Y_ajdk_p=NULL;
  double **J_ajdk_p=NULL,**K_ajdk_p=NULL,**U_ajdk_p=NULL;
  unsigned char *A_bmc_b=NULL,*Y_bmc_b=NULL;
  double tmp_sum=0,tmp_num=0;
  int tmp_d1=0,tmp_d2=0;
  if (verbose){ printf(" %% [entering dcg_load_A_p]\n"); wkspace_printf();}
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
  if (GLOBAL_A_p_name!=NULL &&  strcmp(GLOBAL_A_p_name,"\0")){
    mda_read_d2_r8(GLOBAL_A_p_name,&tmp_d1,&tmp_d2,NULL); if (tmp_d1!=A_pcols){ printf(" %% Warning! tmp_d1 %d A_pcols %d in dcg.c\n",tmp_d1,A_pcols);}
    mda_read_d2_r8(GLOBAL_A_p_name,NULL,NULL,A_p);
    /* if (GLOBAL_A_p_name!=NULL &&  strcmp(GLOBAL_A_p_name,"\0")){ } */}
  else /* if (GLOBAL_A_p_name==NULL || !strcmp(GLOBAL_A_p_name,"\0")) */{
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
    /* if (GLOBAL_A_p_name==NULL || !strcmp(GLOBAL_A_p_name,"\0")){ } */}
  if (verbose){ raprintf(A_p,"double",1,A_pcols," %% A_p: "); wkspace_printf();}
  calc_A_ajdk(A_p,A_pcols,A_ajdk_p); if (verbose){ raprintf(*A_ajdk_p,"double",A_pcols,AJDK_TOT," %% A_ajdk: "); wkspace_printf();}
  if (GLOBAL_Y_p_name!=NULL &&  strcmp(GLOBAL_Y_p_name,"\0")){
    mda_read_d2_r8(GLOBAL_Y_p_name,&tmp_d1,&tmp_d2,NULL); if (tmp_d1!=Y_pcols){ printf(" %% Warning! tmp_d1 %d Y_pcols %d in dcg.c\n",tmp_d1,Y_pcols);}
    mda_read_d2_r8(GLOBAL_Y_p_name,NULL,NULL,Y_p);
    /* if (GLOBAL_Y_p_name!=NULL &&  strcmp(GLOBAL_Y_p_name,"\0")){ } */}
  else /* if (GLOBAL_Y_p_name==NULL || !strcmp(GLOBAL_Y_p_name,"\0")) */{
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
    /* if (GLOBAL_Y_p_name==NULL || !strcmp(GLOBAL_Y_p_name,"\0")){ } */}
  if (verbose){ raprintf(Y_p,"double",1,Y_pcols," %% Y_p: "); wkspace_printf();}
  calc_A_ajdk(Y_p,Y_pcols,Y_ajdk_p); if (verbose){ raprintf(*Y_ajdk_p,"double",Y_pcols,AJDK_TOT," %% Y_ajdk: "); wkspace_printf();}
  sprintf(tempchar,"%s/AZ_rsum.mda",GLOBAL_DIR_NAME); mda_write_d2_r8(tempchar,A_ncols,1,AZ_rsum);
  sprintf(tempchar,"%s/YW_rsum.mda",GLOBAL_DIR_NAME); mda_write_d2_r8(tempchar,Y_ncols,1,YW_rsum);
  sprintf(tempchar,"%s/A_p.mda",GLOBAL_DIR_NAME); mda_write_d2_r8(tempchar,A_pcols,1,A_p);
  sprintf(tempchar,"%s/Y_p.mda",GLOBAL_DIR_NAME); mda_write_d2_r8(tempchar,Y_pcols,1,Y_p);
  if (verbose){ printf(" %% [finished dcg_load_A_p]\n");  wkspace_printf();}
}

void dcg_init_A_p(double mrnd,struct dcg_ajdk *D)
{
  /* initializes A_p and Y_p for each column using row-mask mr_b, not mr_j */
  int verbose=0;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_;
  int nb=0,nc=0;
  if (verbose){ printf(" %% [entering dcg_init_A_p]\n");}
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
  if (verbose){ printf(" %% [finished dcg_init_A_p]\n");}
}

void dcg_copy_QX(struct dcg_ajdk *D,struct dcg_ajdk *D_in)
{
  int nbins = D->nbins; struct dcg_single **E_ = D->E_; struct dcg_single **E_in_ = D_in->E_;
  int nb1=0; struct dcg_single *E=NULL;
  D->out_iteration=0; D->out_trace_length = 6; memcpy(D->out_trace,D_in->out_trace,(D->A_ncols+D->A_nrows_total)*D->out_trace_length*sizeof(double));
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  memcpy(D->out_xdrop_a,D_in->out_xdrop_a,(D->A_ncols+D->A_nrows_total)*2*sizeof(int));
  memcpy(D->out_xdrop_b,D_in->out_xdrop_b,(D->A_ncols+D->A_nrows_total)*2*sizeof(int));
  memcpy(D->QC_TYnWtS_nrm,D_in->QC_TYnWtS_nrm,D->Y_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_TAnZtS_nrm,D_in->QC_TAnZtS_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_TYnYtT_nrm,D_in->QC_TYnYtT_nrm,D->Y_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_TAnAtT_nrm,D_in->QC_TAnAtT_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_A_svalue,D_in->QC_A_svalue,D->A_ncols*sizeof(double));
  memcpy(D->QC_A_index_local_mc_a,D_in->QC_A_index_local_mc_a,D->A_ncols*sizeof(int));
  memcpy(D->QC_A_index_local_mc_b,D_in->QC_A_index_local_mc_b,D->A_ncols*sizeof(int));
  memcpy(D->QC_A_index_local_mc_j,D_in->QC_A_index_local_mc_j,D->A_ncols*sizeof(int));
  memcpy(D->QR_A_svalue,D_in->QR_A_svalue,D->A_nrows_total*sizeof(double));
  memcpy(D->QR_A_index_global_mr_a,D_in->QR_A_index_global_mr_a,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_A_index_global_mr_b,D_in->QR_A_index_global_mr_b,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_A_index_global_mr_j,D_in->QR_A_index_global_mr_j,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_A_index_local_mr_a,D_in->QR_A_index_local_mr_a,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_A_index_local_mr_b,D_in->QR_A_index_local_mr_b,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_A_index_local_mr_j,D_in->QR_A_index_local_mr_j,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_A_index_local_nb,D_in->QR_A_index_local_nb,D->A_nrows_total*sizeof(int));
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    memcpy(E->QR_A_index_global_mr_a,E_in_[nb1]->QR_A_index_global_mr_a,E->A_nrows*sizeof(int));
    memcpy(E->QR_A_index_global_mr_b,E_in_[nb1]->QR_A_index_global_mr_b,E->A_nrows*sizeof(int));
    memcpy(E->QR_A_index_global_mr_j,E_in_[nb1]->QR_A_index_global_mr_j,E->A_nrows*sizeof(int));
    memcpy(E->QR_TYnWtS_nrm,E_in_[nb1]->QR_TYnWtS_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_TAnZtS_nrm,E_in_[nb1]->QR_TAnZtS_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_TYnYtT_nrm,E_in_[nb1]->QR_TYnYtT_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_TAnAtT_nrm,E_in_[nb1]->QR_TAnAtT_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    /*  for (nb1=0;nb1<nbins;nb1++){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  memcpy(D->out_xgrow_a,D_in->out_xgrow_a,(D->A_ncols+D->A_nrows_total)*2*sizeof(int));
  memcpy(D->out_xgrow_b,D_in->out_xgrow_b,(D->A_ncols+D->A_nrows_total)*2*sizeof(int));
  memcpy(D->QC_RUnWtS_nrm,D_in->QC_RUnWtS_nrm,D->Y_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_RUnYtT_nrm,D_in->QC_RUnYtT_nrm,D->Y_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_RJnZtS_nrm,D_in->QC_RJnZtS_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_RJnAtT_nrm,D_in->QC_RJnAtT_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_TKnVtS_nrm,D_in->QC_TKnVtS_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_TKnKtT_nrm,D_in->QC_TKnKtT_nrm,D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double));
  memcpy(D->QC_K_svalue,D_in->QC_K_svalue,D->A_ncols*sizeof(double));
  memcpy(D->QC_K_index_local_mc_a,D_in->QC_K_index_local_mc_a,D->A_ncols*sizeof(int));
  memcpy(D->QC_K_index_local_mc_b,D_in->QC_K_index_local_mc_b,D->A_ncols*sizeof(int));
  memcpy(D->QC_K_index_local_mc_j,D_in->QC_K_index_local_mc_j,D->A_ncols*sizeof(int));
  memcpy(D->QR_J_svalue,D_in->QR_J_svalue,D->A_nrows_total*sizeof(double));
  memcpy(D->QR_J_index_global_mr_a,D_in->QR_J_index_global_mr_a,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_J_index_global_mr_b,D_in->QR_J_index_global_mr_b,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_J_index_global_mr_j,D_in->QR_J_index_global_mr_j,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_J_index_local_mr_a,D_in->QR_J_index_local_mr_a,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_J_index_local_mr_b,D_in->QR_J_index_local_mr_b,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_J_index_local_mr_j,D_in->QR_J_index_local_mr_j,D->A_nrows_total*sizeof(int));
  memcpy(D->QR_J_index_local_nb,D_in->QR_J_index_local_nb,D->A_nrows_total*sizeof(int));
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    memcpy(E->QR_RUnWtS_nrm,E_in_[nb1]->QR_RUnWtS_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_RUnYtT_nrm,E_in_[nb1]->QR_RUnYtT_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_RJnZtS_nrm,E_in_[nb1]->QR_RJnZtS_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    memcpy(E->QR_RJnAtT_nrm,E_in_[nb1]->QR_RJnAtT_nrm,E->A_nrows*D->T_ncols*nbins*sizeof(double));
    /*  for (nb1=0;nb1<nbins;nb1++){ } */}
}

void dcg_init_QX(struct dcg_ajdk *D)
{
  int verbose=0;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_;
  int nb1=0; struct dcg_single *E=NULL; int mx_a=0,mx_b=0,mx_j=0,mr_a=0,mr_b=0,mr_j=0;
  int mc_a=0,mc_b=0,mc_j=0;
  D->out_iteration=0; D->out_trace_length = 6; D->out_trace = (double *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*D->out_trace_length*sizeof(double)); if (!D->out_trace){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  D->out_xdrop_a = (int *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*2*sizeof(int)); if (!D->out_xdrop_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->out_xdrop_b = (int *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*2*sizeof(int)); if (!D->out_xdrop_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QC_TYnWtS_nrm = (double *) wkspace_all0c(D->Y_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_TYnWtS_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_TAnZtS_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_TAnZtS_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_TYnYtT_nrm = (double *) wkspace_all0c(D->Y_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_TYnYtT_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_TAnAtT_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_TAnAtT_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_A_svalue=(double *)wkspace_all0c(D->A_ncols*sizeof(double)); if (!D->QC_A_svalue){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QC_A_index_local_mc_a=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_A_index_local_mc_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QC_A_index_local_mc_b=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_A_index_local_mc_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QC_A_index_local_mc_j=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_A_index_local_mc_j){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_A_svalue=(double *)wkspace_all0c(D->A_nrows_total*sizeof(double)); if (!D->QR_A_svalue){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_A_index_global_mr_a=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_A_index_global_mr_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_A_index_global_mr_b=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_A_index_global_mr_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_A_index_global_mr_j=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_A_index_global_mr_j){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_A_index_local_mr_a=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_A_index_local_mr_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_A_index_local_mr_b=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_A_index_local_mr_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_A_index_local_mr_j=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_A_index_local_mr_j){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_A_index_local_nb=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_A_index_local_nb){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  mx_a=0; mx_b=0; mx_j=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    E->QR_A_index_global_mr_a = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_A_index_global_mr_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
    E->QR_A_index_global_mr_b = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_A_index_global_mr_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
    E->QR_A_index_global_mr_j = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_A_index_global_mr_j){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
    mr_a=0;mr_b=0;mr_j=0;
    while (mr_a<E->A_nrows){
      E->QR_A_index_global_mr_a[mr_a] = mx_a;
      if (bget__on(E->A_bmr_b,mr_a)){
	E->QR_A_index_global_mr_b[mr_a] = mx_b;
	if (bget__on(E->A_bmr_j,mr_a)){
	  E->QR_A_index_global_mr_j[mr_a] = mx_j;
	  mr_j++; mx_j++; /* if (bget__on(E->A_bmr_j,mr_a)){ } */}
	mr_b++; mx_b++; /* if (bget__on(E->A_bmr_b,mr_b)){ } */}
      mr_a++; mx_a++; /* while (mr_a<E->A_nrows){ } */}
    E->QR_TYnWtS_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_TYnWtS_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
    E->QR_TAnZtS_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_TAnZtS_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
    E->QR_TYnYtT_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_TYnYtT_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
    E->QR_TAnAtT_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_TAnAtT_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
    /*  for (nb1=0;nb1<nbins;nb1++){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  D->out_xgrow_a = (int *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*2*sizeof(int)); if (!D->out_xgrow_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->out_xgrow_b = (int *)wkspace_all0c((D->A_ncols+D->A_nrows_total)*2*sizeof(int)); if (!D->out_xgrow_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QC_RUnWtS_nrm = (double *) wkspace_all0c(D->Y_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_RUnWtS_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_RUnYtT_nrm = (double *) wkspace_all0c(D->Y_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_RUnYtT_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_RJnZtS_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_RJnZtS_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_RJnAtT_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_RJnAtT_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_TKnVtS_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_TKnVtS_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_TKnKtT_nrm = (double *) wkspace_all0c(D->A_ncols*D->T_ncols*nbins*nbins*sizeof(double)); if (!D->QC_TKnKtT_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
  D->QC_K_svalue=(double *)wkspace_all0c(D->A_ncols*sizeof(double)); if (!D->QC_K_svalue){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QC_K_index_local_mc_a=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_K_index_local_mc_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QC_K_index_local_mc_b=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_K_index_local_mc_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QC_K_index_local_mc_j=(int *)wkspace_all0c(D->A_ncols*sizeof(int)); if (!D->QC_K_index_local_mc_j){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_J_svalue=(double *)wkspace_all0c(D->A_nrows_total*sizeof(double)); if (!D->QR_J_svalue){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_J_index_global_mr_a=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_J_index_global_mr_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_J_index_global_mr_b=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_J_index_global_mr_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_J_index_global_mr_j=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_J_index_global_mr_j){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_J_index_local_mr_a=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_J_index_local_mr_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_J_index_local_mr_b=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_J_index_local_mr_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_J_index_local_mr_j=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_J_index_local_mr_j){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  D->QR_J_index_local_nb=(int *)wkspace_all0c(D->A_nrows_total*sizeof(int)); if (!D->QR_J_index_local_nb){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
  mx_a=0; mx_b=0; mx_j=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    E->QR_J_index_global_mr_a = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_J_index_global_mr_a){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
    E->QR_J_index_global_mr_b = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_J_index_global_mr_b){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
    E->QR_J_index_global_mr_j = (int *)wkspace_all0c(E->A_nrows*sizeof(int)); if (!E->QR_J_index_global_mr_j){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");} 
    mr_a=0;mr_b=0;mr_j=0;
    while (mr_a<E->A_nrows){
      E->QR_J_index_global_mr_a[mr_a] = mx_a;
      if (bget__on(E->J_bmr_b,mr_a)){
	E->QR_J_index_global_mr_b[mr_a] = mx_b;
	if (bget__on(E->J_bmr_j,mr_a)){
	  E->QR_J_index_global_mr_j[mr_a] = mx_j;
	  mr_j++; mx_j++; /* if (bget__on(E->J_bmr_j,mr_a)){ } */}
	mr_b++; mx_b++; /* if (bget__on(E->J_bmr_b,mr_b)){ } */}
      mr_a++; mx_a++; /* while (mr_a<E->A_nrows){ } */}
    E->QR_RUnWtS_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_RUnWtS_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
    E->QR_RUnYtT_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_RUnYtT_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
    E->QR_RJnZtS_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_RJnZtS_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
    E->QR_RJnAtT_nrm = (double *) wkspace_all0c(E->A_nrows*D->T_ncols*nbins*sizeof(double)); if (!E->QR_RJnAtT_nrm){ printf(" %% Warning! not enough memory in dcg_init_X_QX\n");}
    /*  for (nb1=0;nb1<nbins;nb1++){ } */}
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  D->out_xgrow_ij=0;
  mx_a=0; mx_b=0; mx_j=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    mr_a=0;mr_b=0;mr_j=0;
    while (mr_a<E->A_nrows){
      if (bget__on(E->A_bmr_b,mr_a)){
	if (bget__on(E->A_bmr_j,mr_a)){
	  if (E->QR_A_index_global_mr_a[mr_a]!=mx_a){ printf(" %% Warning! index mismatch in mx_a in dcg_init_QX\n");}
	  D->out_xgrow_a[0 + D->out_xgrow_ij*2] = mx_a; D->out_xgrow_a[1 + D->out_xgrow_ij*2] = -1;
	  if (E->QR_A_index_global_mr_b[mr_a]!=mx_b){ printf(" %% Warning! index mismatch in mx_b in dcg_init_QX\n");}
	  D->out_xgrow_b[0 + D->out_xgrow_ij*2] = mx_b; D->out_xgrow_b[1 + D->out_xgrow_ij*2] = -1;
	  D->out_xgrow_ij++;
	  mr_j++; mx_j++; /* if (bget__on(E->A_bmr_j,mr_a)){ } */}
	mr_b++; mx_b++; /* if (bget__on(E->A_bmr_b,mr_b)){ } */}
      mr_a++; mx_a++; /* while (mr_a<E->A_nrows){ } */}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose>-1){ printf(" %% writing %d/(A:%d,%d,%d vs J:%d,%d,%d) rows to out_xgrow\n",mx_j,D->A_nrows_total,D->A_rpop_b_total,D->A_rpop_j_total,D->A_nrows_total,D->J_rpop_b_total,D->J_rpop_j_total);}
  mc_a=0;mc_b=0;mc_j=0;
  while (mc_a<D->A_ncols){
    if (bget__on(D->A_bmc_b,mc_a)){
      if (bget__on(D->A_bmc_j,mc_a)){
	D->out_xgrow_a[1 + D->out_xgrow_ij*2] = mc_a; D->out_xgrow_a[0 + D->out_xgrow_ij*2] = -1; 
	D->out_xgrow_b[1 + D->out_xgrow_ij*2] = mc_b; D->out_xgrow_b[0 + D->out_xgrow_ij*2] = -1; 
	D->out_xgrow_ij++;
	mc_j++; /* if (bget__on(D->A_bmc_j,mc_a)){ } */}
      mc_b++; /* if (bget__on(D->A_bmc_b,mc_b)){ } */}
    mc_a++; /* while (mc_a<D->A_ncols){ } */}
  if (verbose>-1){ printf(" %% writing %d/(A:%d,%d,%d vs K:%d,%d,%d) cols to out_xgrow\n",mc_j,D->A_ncols,D->A_cpop_b,D->A_cpop_j,D->A_ncols,D->K_cpop_b,D->K_cpop_j);}
}

void dcg_M_mxset(struct dcg_ajdk *D)
{
  /* updates cpop_j and rpop_j, copies D->A_bmc_j into M_An->mc_j, and copies E_[nb]->A_bmr_j into M_An->mr_j. */
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_;
  int nb=0;
  if (verbose){ printf(" %% setting M_An->mr_j, etc.\n");}
  D->A_cpop_j = popcount_uchar_array(D->A_bmc_j,D->A_mc_length);
  D->K_cpop_j = popcount_uchar_array(D->K_bmc_j,D->A_mc_length);
  D->Y_cpop_j = popcount_uchar_array(D->Y_bmc_j,D->Y_mc_length);
  D->T_cpop_j = popcount_uchar_array(D->T_bmc_j,D->T_mc_length);
  D->A_rpop_j_total=0;D->Z_rpop_j_total=0;
  for (nb=0;nb<nbins;nb++){ 
    E_[nb]->A_rpop_j = popcount_uchar_array(E_[nb]->A_bmr_j,E_[nb]->A_mr_length); D->A_rpop_j_total += E_[nb]->A_rpop_j;
    E_[nb]->J_rpop_j = popcount_uchar_array(E_[nb]->J_bmr_j,E_[nb]->A_mr_length); 
    E_[nb]->Z_rpop_j = popcount_uchar_array(E_[nb]->Z_bmr_j,E_[nb]->Z_mr_length); D->Z_rpop_j_total += E_[nb]->Z_rpop_j;
    M_mxset(E_[nb]->M_An,E_[nb]->A_bmr_j,D->A_bmc_j); M_mxset(E_[nb]->M_At,D->A_bmc_j,E_[nb]->A_bmr_j);
    M_mxset(E_[nb]->M_Zn,E_[nb]->Z_bmr_j,D->A_bmc_j); M_mxset(E_[nb]->M_Zt,D->A_bmc_j,E_[nb]->Z_bmr_j);
    M_mxset(E_[nb]->M_Yn,E_[nb]->A_bmr_j,D->Y_bmc_j); M_mxset(E_[nb]->M_Yt,D->Y_bmc_j,E_[nb]->A_bmr_j);
    M_mxset(E_[nb]->M_Wn,E_[nb]->Z_bmr_j,D->Y_bmc_j); M_mxset(E_[nb]->M_Wt,D->Y_bmc_j,E_[nb]->Z_bmr_j);
    M_mxset(E_[nb]->M_Tn,E_[nb]->A_bmr_j,D->T_bmc_j); M_mxset(E_[nb]->M_Tt,D->T_bmc_j,E_[nb]->A_bmr_j);
    M_mxset(E_[nb]->M_Sn,E_[nb]->Z_bmr_j,D->T_bmc_j); M_mxset(E_[nb]->M_St,D->T_bmc_j,E_[nb]->Z_bmr_j);
    M_mxset(E_[nb]->M_Jn,E_[nb]->J_bmr_j,D->A_bmc_j); M_mxset(E_[nb]->M_Jt,D->A_bmc_j,E_[nb]->J_bmr_j);
    M_mxset(E_[nb]->M_Kn,E_[nb]->A_bmr_j,D->K_bmc_j); M_mxset(E_[nb]->M_Kt,D->K_bmc_j,E_[nb]->A_bmr_j);
    M_mxset(E_[nb]->M_Vn,E_[nb]->Z_bmr_j,D->K_bmc_j); M_mxset(E_[nb]->M_Vt,D->K_bmc_j,E_[nb]->Z_bmr_j);
    M_mxset(E_[nb]->M_Un,E_[nb]->J_bmr_j,D->Y_bmc_j); M_mxset(E_[nb]->M_Ut,D->Y_bmc_j,E_[nb]->J_bmr_j);
    M_mxset(E_[nb]->M_Rn,E_[nb]->J_bmr_j,D->T_bmc_j); M_mxset(E_[nb]->M_Rt,D->T_bmc_j,E_[nb]->J_bmr_j);
    /* for (nb=0;nb<nbins;nb++){ } */}
}

void dcg_copy(struct dcg_ajdk *D,struct dcg_ajdk *D_in)
{
  int verbose=0;
  struct dcg_single **E_ = D->E_; struct dcg_single **E_in_ = D_in->E_;
  struct dcg_double **F_ = D->F_; struct dcg_double **F_in_ = D_in->F_;
  int nbins = GLOBAL_NBINS, nbx=0,nb1=0,nb2=0;
  if (verbose){ printf(" %% [entering dcg_copy] \n");}
  GLOBAL_tic(0);
  if (verbose){ printf(" %% calling dcg_ajdk_copy: \n");}
  dcg_ajdk_copy(D,D_in);
  if (verbose){ printf(" %% finished dcg_ajdk_copy: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ 
    if (verbose){ printf(" %% calling dcg_single_copy_M_An for nb1 %d: \n",nb1);} 
    dcg_single_copy_M_An(E_[nb1],E_in_[nb1]);
    if (verbose){ printf(" %% finished dcg_single_copy_M_An for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% calling dcg_X_nrows_total: \n");}
  dcg_X_nrows_total(D);
  if (verbose){ printf(" %% finished dcg_X_nrows_total: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ 
    if (verbose){ printf(" %% calling dcg_single_copy_lf for nb1 %d: \n",nb1);}
    dcg_single_copy_lf(E_[nb1],E_in_[nb1]);
    if (verbose){ printf(" %% finished dcg_single_copy_lf for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; 
      if (verbose){ printf(" %% calling dcg_double_copy_lf for nb1 %d nb2 %d nbx %d: \n",nb1,nb2,nbx);}
      dcg_double_copy_lf(F_[nbx],F_in_[nbx]); 
      if (verbose){ printf(" %% finished dcg_double_copy_lf for nb1 %d nb2 %d nbx %d: ",nb1,nb2,nbx); wkspace_printf();}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_toc(0,verbose," %% generating data matrices: ");
  if (verbose){ printf(" %% calling dcg_copy_A_p: \n");}
  dcg_copy_A_p(D,D_in);
  if (verbose){ printf(" %% finished dcg_copy_A_p: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcg_copy_QX\n");}
  dcg_copy_QX(D,D_in);
  if (verbose){ printf(" %% finished dcg_copy_QX: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcg_M_mxset\n");}
  dcg_M_mxset(D);
  if (verbose){ printf(" %% finished dcg_M_mxset: "); wkspace_printf();}
  if (verbose){ printf(" %% [finished dcg_copy] "); wkspace_printf();}
}

void dcg_load(struct dcg_ajdk **D_p,struct dcg_single ***E_p,struct dcg_double ***F_p)
{
  int verbose=GLOBAL_verbose+2;
  int nbins = GLOBAL_NBINS, nbx=0,nb1=0,nb2=0;
  struct dcg_ajdk *D=NULL;
  struct dcg_single *E=NULL,*E_nb1=NULL,*E_nb2=NULL;
  struct dcg_double *F=NULL;
  if (verbose){ printf(" %% [entering dcg_load] \n");}
  if (*D_p==NULL){
    (*D_p) = (struct dcg_ajdk *) wkspace_all0c(sizeof(struct dcg_ajdk)); D = (*D_p); D->bitj = BITJ; D->nbins = nbins; D->E_ = NULL; D->F_ = NULL;
    /* sprintf(D->QR_strategy,"%s",QR_strategy); sprintf(D->QC_strategy,"%s",QC_strategy); */
    /* if (*D_p==NULL){ } */}
  if (*E_p==NULL){
    (*E_p) = (struct dcg_single **) wkspace_all0c(nbins*sizeof(struct dcg_single *));
    for (nb1=0;nb1<nbins;nb1++){ (*E_p)[nb1] = (struct dcg_single *) wkspace_all0c(sizeof(struct dcg_single)*1); E = (*E_p)[nb1]; E->D = (*D_p); E->nb=nb1;}
    (*D_p)->E_ = (*E_p);
    /* if (*E_p==NULL){ } */}
  if (*F_p==NULL){
    (*F_p) = (struct dcg_double **) wkspace_all0c(nbins*nbins*sizeof(struct dcg_double *));
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins;
	(*F_p)[nbx] = (struct dcg_double *) wkspace_all0c(sizeof(struct dcg_double)*1); F = (*F_p)[nbx]; F->D = (*D_p); F->nb1=nb1; F->nb2=nb2; F->E_nb1 = (*E_p)[nb1]; F->E_nb2 = (*E_p)[nb2];
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    (*D_p)->F_ = (*F_p);
    /* if (*F_p==NULL){ } */}
  D = (*D_p); 
  GLOBAL_tic(0);
  if (verbose){ printf(" %% calling dcg_ajdk_load: \n");}
  dcg_ajdk_load(D);
  if (verbose){ printf(" %% finished dcg_ajdk_load: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling dcg_single_load_M_An for nb1 %d: \n",nb1);} 
    dcg_single_load_M_An(E);
    if (verbose){ printf(" %% finished dcg_single_load_M_An for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% calling dcg_X_nrows_total: \n");}
  dcg_X_nrows_total(D);
  if (verbose){ printf(" %% finished dcg_X_nrows_total: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling dcg_single_init_lf for nb1 %d: \n",nb1);}
    dcg_single_init_lf(E);
    if (verbose){ printf(" %% finished dcg_single_init_lf for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; E_nb1 = (*E_p)[nb1]; E_nb2 = (*E_p)[nb2]; F = (*F_p)[nbx];
      if (verbose){ printf(" %% calling dcg_double_init_lf for nb1 %d nb2 %d nbx %d: \n",nb1,nb2,nbx);}
      dcg_double_init_lf(F); 
      if (verbose){ printf(" %% finished dcg_double_init_lf for nb1 %d nb2 %d nbx %d: ",nb1,nb2,nbx); wkspace_printf();}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_toc(0,verbose," %% generating data matrices: ");
  if (verbose){ printf(" %% calling dcg_load_A_p: \n");}
  dcg_load_A_p(D);
  if (verbose){ printf(" %% finished dcg_load_A_p: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcg_init_QX\n");}
  dcg_init_QX(D);
  if (verbose){ printf(" %% finished dcg_init_QX: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcg_M_mxset\n");}
  dcg_M_mxset(D);
  if (verbose){ printf(" %% finished dcg_M_mxset: "); wkspace_printf();}
  if (verbose){ printf(" %% [finished dcg_load] "); wkspace_printf();}
}

void dcg_init(char *error_vs_speed,double mrnd,struct dcg_ajdk **D_p,struct dcg_single ***E_p,struct dcg_double ***F_p)
{
  int verbose=GLOBAL_verbose;
  int nbins = GLOBAL_NBINS, nbx=0,nb1=0,nb2=0;
  struct dcg_ajdk *D=NULL;
  struct dcg_single *E=NULL,*E_nb1=NULL,*E_nb2=NULL;
  struct dcg_double *F=NULL;
  if (verbose){ printf(" %% [entering dcg_init] \n");}
  if (*D_p==NULL){
    (*D_p) = (struct dcg_ajdk *) wkspace_all0c(sizeof(struct dcg_ajdk)); D = (*D_p); D->bitj = BITJ; D->nbins = nbins; D->E_ = NULL; D->F_ = NULL;
    /* sprintf(D->QR_strategy,"%s",QR_strategy); sprintf(D->QC_strategy,"%s",QC_strategy); */
    /* if (*D_p==NULL){ } */}
  if (*E_p==NULL){
    (*E_p) = (struct dcg_single **) wkspace_all0c(nbins*sizeof(struct dcg_single *));
    for (nb1=0;nb1<nbins;nb1++){ (*E_p)[nb1] = (struct dcg_single *) wkspace_all0c(sizeof(struct dcg_single)*1); E = (*E_p)[nb1]; E->D = (*D_p); E->nb=nb1;}
    (*D_p)->E_ = (*E_p);
    /* if (*E_p==NULL){ } */}
  if (*F_p==NULL){
    (*F_p) = (struct dcg_double **) wkspace_all0c(nbins*nbins*sizeof(struct dcg_double *));
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins;
	(*F_p)[nbx] = (struct dcg_double *) wkspace_all0c(sizeof(struct dcg_double)*1); F = (*F_p)[nbx]; F->D = (*D_p); F->nb1=nb1; F->nb2=nb2; F->E_nb1 = (*E_p)[nb1]; F->E_nb2 = (*E_p)[nb2];
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    (*D_p)->F_ = (*F_p);
    /* if (*F_p==NULL){ } */}
  D = (*D_p);
  GLOBAL_tic(0);
  if (verbose){ printf(" %% calling dcg_ajdk_init: \n");}
  dcg_ajdk_init(mrnd,D);
  if (verbose){ printf(" %% finished dcg_ajdk_init: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling dcg_single_init_M_An for nb1 %d: \n",nb1);}
    dcg_single_init_M_An(error_vs_speed,mrnd,E);
    if (verbose){ printf(" %% finished dcg_single_init_M_An for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% calling dcg_X_nrows_total: \n");}
  dcg_X_nrows_total(D);
  if (verbose){ printf(" %% finished dcg_X_nrows_total: "); wkspace_printf();}
  for (nb1=0;nb1<nbins;nb1++){ E = (*E_p)[nb1];
    if (verbose){ printf(" %% calling dcg_single_init_lf for nb1 %d: \n",nb1);}
    dcg_single_init_lf(E);
    if (verbose){ printf(" %% finished dcg_single_init_lf for nb1 %d: ",nb1); wkspace_printf();}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; E_nb1 = (*E_p)[nb1]; E_nb2 = (*E_p)[nb2]; F = (*F_p)[nbx];
      if (verbose){ printf(" %% calling dcg_double_init_lf for nb1 %d nb2 %d nbx %d: \n",nb1,nb2,nbx);}
      dcg_double_init_lf(F);
      if (verbose){ printf(" %% finished dcg_double_init_lf for nb1 %d nb2 %d nbx %d: ",nb1,nb2,nbx); wkspace_printf();}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_toc(0,verbose," %% generating data matrices: ");
  if (verbose){ printf(" %% calling dcg_init_A_p: \n");}
  dcg_init_A_p(mrnd,D);
  if (verbose){ printf(" %% finished dcg_init_A_p: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcg_init_QX\n");}
  dcg_init_QX(D);
  if (verbose){ printf(" %% finished dcg_init_QX: "); wkspace_printf();}
  if (verbose){ printf(" %% calling dcg_M_mxset\n");}
  dcg_M_mxset(D);
  if (verbose){ printf(" %% finished dcg_M_mxset: "); wkspace_printf();}
  if (verbose){ printf(" %% [finished dcg_init] "); wkspace_printf();}
}

void dcg_init_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  struct dcg_ajdk *D=NULL;struct dcg_single **E_=NULL; struct dcg_double **F_=NULL;
  dcg_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); dcg_init_QX(D);
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
  M_handle_printf(E_[0]->M_Jn,verbose," %% E_[0]->M_Jn[0]: "); M_handle_printf(E_[0]->M_Jt,verbose," %% E_[0]->M_Jt[0]: ");
  M_handle_printf(E_[1]->M_Jn,verbose," %% E_[1]->M_Jn[1]: "); M_handle_printf(E_[1]->M_Jt,verbose," %% E_[1]->M_Jt[1]: ");
  M_handle_printf(E_[0]->M_Kn,verbose," %% E_[0]->M_Kn[0]: "); M_handle_printf(E_[0]->M_Kt,verbose," %% E_[0]->M_Kt[0]: ");
  M_handle_printf(E_[1]->M_Kn,verbose," %% E_[1]->M_Kn[1]: "); M_handle_printf(E_[1]->M_Kt,verbose," %% E_[1]->M_Kt[1]: ");
  M_handle_printf(E_[0]->M_Vn,verbose," %% E_[0]->M_Vn[0]: "); M_handle_printf(E_[0]->M_Vt,verbose," %% E_[0]->M_Vt[0]: ");
  M_handle_printf(E_[1]->M_Vn,verbose," %% E_[1]->M_Vn[1]: "); M_handle_printf(E_[1]->M_Vt,verbose," %% E_[1]->M_Vt[1]: ");
  M_handle_printf(E_[0]->M_Un,verbose," %% E_[0]->M_Un[0]: "); M_handle_printf(E_[0]->M_Ut,verbose," %% E_[0]->M_Ut[0]: ");
  M_handle_printf(E_[1]->M_Un,verbose," %% E_[1]->M_Un[1]: "); M_handle_printf(E_[1]->M_Ut,verbose," %% E_[1]->M_Ut[1]: ");
  M_handle_printf(E_[0]->M_Rn,verbose," %% E_[0]->M_Rn[0]: "); M_handle_printf(E_[0]->M_Rt,verbose," %% E_[0]->M_Rt[0]: ");
  M_handle_printf(E_[1]->M_Rn,verbose," %% E_[1]->M_Rn[1]: "); M_handle_printf(E_[1]->M_Rt,verbose," %% E_[1]->M_Rt[1]: ");
  wkspace_printf();
}
