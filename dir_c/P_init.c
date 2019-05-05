#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void P_handle_make_excerpt_A_nrows_(struct P_handle *P,struct L_handle ***U__,struct M_handle ***M__)
{
  int nbins = P->D->nbins; int nb=0;
  if (U__!=NULL){
    *(U__) = (struct L_handle **) wkspace_all0c(nbins*sizeof(struct L_handle *));
    for (nb=0;nb<nbins;nb++){ 
      (*U__)[nb] = L_handle_make((unsigned long long int)P->D->E_[nb]->A_nrows/* *(unsigned long long int)P->D->T_ncols */*(unsigned long long int)P->rank);
      (*U__)[nb]->spacing_row = SPACING_a; (*U__)[nb]->row_stride = P->D->E_[nb]->A_nrows;
      (*U__)[nb]->spacing_col = SPACING_a; (*U__)[nb]->col_stride = P->rank;
      (*U__)[nb]->spacing_lyr = SPACING_a; (*U__)[nb]->lyr_stride = 1 /* P->D->T_ncols */;
      L_zero((*U__)[nb]);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (U__!=NULL){ } */}    
  if (M__!=NULL){
    *(M__) = (struct M_handle **) wkspace_all0c(nbins*sizeof(struct M_handle *));
    for (nb=0;nb<nbins;nb++){ (*M__)[nb] = M_handle_w_make(BITJ,P->b_mlt,P->D->E_[nb]->A_nrows/* *P->D->T_ncols */,P->rank);}
    /* if (M__!=NULL){ } */}
}

void P_handle_make_excerpt_Z_nrows_(struct P_handle *P,struct L_handle ***U__,struct M_handle ***M__)
{
  int nbins = P->D->nbins; int nb=0;
  if (U__!=NULL){
    *(U__) = (struct L_handle **) wkspace_all0c(nbins*sizeof(struct L_handle *));
    for (nb=0;nb<nbins;nb++){ 
      (*U__)[nb] = L_handle_make((unsigned long long int)P->D->E_[nb]->Z_nrows/* *(unsigned long long int)P->D->T_ncols */*(unsigned long long int)P->rank);
      (*U__)[nb]->spacing_row = SPACING_a; (*U__)[nb]->row_stride = P->D->E_[nb]->Z_nrows;
      (*U__)[nb]->spacing_col = SPACING_a; (*U__)[nb]->col_stride = P->rank;
      (*U__)[nb]->spacing_lyr = SPACING_a; (*U__)[nb]->lyr_stride = 1 /* P->D->T_ncols */;
      L_zero((*U__)[nb]);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (U__!=NULL){ } */}
  if (M__!=NULL){
    *(M__) = (struct M_handle **) wkspace_all0c(nbins*sizeof(struct M_handle *));
    for (nb=0;nb<nbins;nb++){ (*M__)[nb] = M_handle_w_make(BITJ,P->b_mlt,P->D->E_[nb]->Z_nrows/* *P->D->T_ncols */,P->rank);}
    /* if (M__!=NULL){ } */}
}

void P_handle_make_excerpt_A_ncols_(struct P_handle *P,struct L_handle ***V__,struct M_handle ***M__)
{
  int nbins = P->D->nbins; int nb=0;
  if (V__!=NULL){
    *(V__) = (struct L_handle **) wkspace_all0c(nbins*sizeof(struct L_handle *));
    for (nb=0;nb<nbins;nb++){ 
      (*V__)[nb] = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank*(unsigned long long int)P->D->T_ncols);
      (*V__)[nb]->spacing_row = SPACING_a; (*V__)[nb]->row_stride = P->D->A_ncols;
      (*V__)[nb]->spacing_col = SPACING_a; (*V__)[nb]->col_stride = P->rank;
      (*V__)[nb]->spacing_lyr = SPACING_a; (*V__)[nb]->lyr_stride = P->D->T_ncols;
      L_zero((*V__)[nb]);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (V__!=NULL){ } */}
  if (M__!=NULL){
    *(M__) = (struct M_handle **) wkspace_all0c(nbins*sizeof(struct M_handle *));
    for (nb=0;nb<nbins;nb++){ (*M__)[nb] = M_handle_w_make(BITJ,P->b_mlt,P->D->A_ncols*P->D->T_ncols,P->rank);}
    /* if (M__!=NULL){ } */}
}

void P_handle_make_excerpt_Y_ncols_(struct P_handle *P,struct L_handle ***V__,struct M_handle ***M__)
{
  int nbins = P->D->nbins; int nb=0;
  if (V__!=NULL){
    *(V__) = (struct L_handle **) wkspace_all0c(nbins*sizeof(struct L_handle *));
    for (nb=0;nb<nbins;nb++){ 
      (*V__)[nb] = L_handle_make((unsigned long long int)P->D->Y_ncols*(unsigned long long int)P->rank*(unsigned long long int)P->D->T_ncols);
      (*V__)[nb]->spacing_row = SPACING_a; (*V__)[nb]->row_stride = P->D->Y_ncols;
      (*V__)[nb]->spacing_col = SPACING_a; (*V__)[nb]->col_stride = P->rank;
      (*V__)[nb]->spacing_lyr = SPACING_a; (*V__)[nb]->lyr_stride = P->D->T_ncols;
      L_zero((*V__)[nb]);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (V__!=NULL){ } */}
  if (M__!=NULL){
    *(M__) = (struct M_handle **) wkspace_all0c(nbins*sizeof(struct M_handle *));
    for (nb=0;nb<nbins;nb++){ (*M__)[nb] = M_handle_w_make(BITJ,P->b_mlt,P->D->Y_ncols*P->D->T_ncols,P->rank);}
    /* if (M__!=NULL){ } */}
}

void P_handle_make_excerpt_A_ncols(struct P_handle *P,struct L_handle **V_,struct M_handle **M_)
{
  int nbins = P->D->nbins; int nb=0;
  if (V_!=NULL){
    *V_ = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank);
    (*V_)->spacing_row = SPACING_a; (*V_)->row_stride = P->D->A_ncols;
    (*V_)->spacing_col = SPACING_a; (*V_)->col_stride = P->rank;
    (*V_)->spacing_lyr = SPACING_a; (*V_)->lyr_stride = 1;
    L_zero(*V_);
    /* if (V_!=NULL){ } */}
  if (M_!=NULL){ *M_ = M_handle_w_make(BITJ,P->b_mlt,P->D->A_ncols,P->rank);}
}

void P_handle_make_excerpt_Y_ncols(struct P_handle *P,struct L_handle **V_,struct M_handle **M_)
{
  int nbins = P->D->nbins; int nb=0;
  *V_ = L_handle_make((unsigned long long int)P->D->Y_ncols*(unsigned long long int)P->rank);
  (*V_)->spacing_row = SPACING_a; (*V_)->row_stride = P->D->Y_ncols;
  (*V_)->spacing_col = SPACING_a; (*V_)->col_stride = P->rank;
  (*V_)->spacing_lyr = SPACING_a; (*V_)->lyr_stride = 1;
  L_zero(*V_);
  *M_ = M_handle_w_make(BITJ,P->b_mlt,P->D->Y_ncols,P->rank);
}

void dcc_ajdk_set_mx_j_to_mx_b(struct dcc_ajdk *D)
{
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL;
  int nl=0;
  for (nl=0;nl<D->A_mc_length;nl++){ D->A_bmc_j[nl] = D->A_bmc_b[nl];}
  for (nl=0;nl<D->Y_mc_length;nl++){ D->Y_bmc_j[nl] = D->Y_bmc_b[nl];}
  for (nl=0;nl<D->T_mc_length;nl++){ D->T_bmc_j[nl] = D->T_bmc_b[nl];}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    for (nl=0;nl<E->A_mr_length;nl++){ E->A_bmr_j[nl] = E->A_bmr_b[nl];}
    for (nl=0;nl<E->Z_mr_length;nl++){ E->Z_bmr_j[nl] = E->Z_bmr_b[nl];}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  dcc_M_mxset(D);
}

void dcc_ajdk_set_QC_index_local(struct dcc_ajdk *D)
{
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL;
  int na_a=0,na_b=0,na_j=0;
  if (verbose){ printf(" %% [entering dcc_ajdk_set_QC_index_local]\n");}
  na_a=0;na_b=0;na_j=0;
  while (na_a<D->A_ncols){
    if (bget__on(D->A_bmc_b,na_a)){
      if (bget__on(D->A_bmc_j,na_a)){
	D->QC_index_local_mc_a[na_j] = na_a; 
	D->QC_index_local_mc_b[na_j] = na_b; 
	D->QC_index_local_mc_j[na_j] = na_j;
	na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
      na_b++; /* if (bget__on(D->A_bmc_b,na_a)){ } */}
    na_a++;/* while (na_a<D->A_ncols){ } */}
  if (verbose>1){
    raprintf(D->QC_index_local_mc_a,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_a pre :");
    raprintf(D->QC_index_local_mc_b,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_b pre :");
    raprintf(D->QC_index_local_mc_j,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_j pre :");
    /* if (verbose){ } */}
  if (verbose){ printf(" %% [finished dcc_ajdk_set_QC_index_local]\n");}
}

void dcc_ajdk_copy_QR_index_global(struct dcc_ajdk *D)
{
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL;
  int ma_a=0,ma_b=0,ma_j=0;
  if (verbose){ printf(" %% [entering dcc_ajdk_copy_QR_index_global]\n");}
  D->A_nrows_total=0;D->A_rpop_b_total=0;D->A_rpop_j_total=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; 
    if (verbose){
      printf(" %% nb1 %d E->A_nrows %d E->A_rpop_b %d E->A_rpop_j %d\n",nb1,E->A_nrows,E->A_rpop_b,E->A_rpop_j);
      if (verbose>1){ sprintf(D->tmpAnchar," %%%% E->A_bmr_b: "); bprintf(E->A_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);}
      if (verbose>1){ sprintf(D->tmpAnchar," %%%% E->A_bmr_j: "); bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
      /* if (verbose){ } */}
    ma_a=0;ma_b=0;ma_j=0;
    while (ma_a<E->A_nrows){
      if (bget__on(E->A_bmr_b,ma_a)){
	if (bget__on(E->A_bmr_j,ma_a)){
	  D->QR_index_global_mr_a[D->A_rpop_j_total] = E->QR_index_global_mr_a[ma_a]; 
	  D->QR_index_global_mr_b[D->A_rpop_j_total] = E->QR_index_global_mr_b[ma_a];
	  ma_j++; D->A_rpop_j_total++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	ma_b++; D->A_rpop_b_total++; /* while (ma_a<E->A_nrows){ } */}
      ma_a++; D->A_nrows_total++; /* while (ma_a<E->A_nrows){ } */}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
    if (verbose>1){
      for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
	sprintf(D->tmpAnchar," %%%% E_[%d]->A_bmr_b: ",nb1); bprintf(E->A_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);
	/* for (nb1=0;nb1<nbins;nb1++){ } */}
      for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
	sprintf(D->tmpAnchar," %%%% E_[%d]->A_bmr_j: ",nb1); bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);
	/* for (nb1=0;nb1<nbins;nb1++){ } */}
      raprintf( D->QR_index_global_mr_a,   "int",1,D->A_rpop_j_total," %%  D->QR_index_global_mr_a :");
      raprintf( D->QR_index_global_mr_b,   "int",1,D->A_rpop_j_total," %%  D->QR_index_global_mr_b :");
      /* if (verbose){ } */}
  if (verbose){ printf(" %% [finished dcc_ajdk_copy_QR_index_global]\n");}
}

struct P_handle *P_handle_make(char *pca_infix,char *out_xdrop_fname,char *V_fname,struct dcc_ajdk *D,int iteration_num,int iteration_max,int iteration_min,int rank,double tolerance,int b_mlt)
{
  int verbose=0;
  int nbins = D->nbins,nb=0,nr=0;
  struct P_handle *P=NULL;
  P = (struct P_handle *)wkspace_all0c(sizeof(struct P_handle)*1);
  unsigned char *bmc1_p=NULL;
  if (verbose){ printf(" %% [entering P_handle_make]\n");}
  if (pca_infix!=NULL){ sprintf(P->infix,pca_infix);} else{ sprintf(P->infix,"\0");}
  if (out_xdrop_fname!=NULL){ sprintf(P->out_xdrop_name,out_xdrop_fname);} else{ sprintf(P->out_xdrop_name,"\0");}
  if (verbose>1){ printf(" %% P->infix %s P->out_xdrop_name %s, initializing structure\n",P->infix,P->out_xdrop_name);}
  P->D = D;
  P->out_xdrop_nlines=0; /* number of lines in out_xdrop */
  P->out_xdrop_nrows=0; /* number of row-indices stored in out_xdrop */
  P->out_xdrop_ncols=0; /* number of col-indices stored in out_xdrop */
  P->out_xdrop=NULL; /* actual array for out_xdrop */
  P->mr_index_sort=NULL; /* row-indices from out_xdrop */
  P->mc_index_sort=NULL; /* col-indices from out_xdrop */
  P->mr_index_local_nb=NULL; /* local nb from row-indices */
  P->mr_index_local_mr=NULL; /* local mr from row-indices */
  P->iteration_num=iteration_num; /* number of iterations */
  P->iteration_min=iteration_min; /* maximum iteration */
  P->iteration_max=iteration_max; /* minimum iteration */
  P->rkeep=NULL; /* array of r-values (rows remaining) */
  P->rdrop=NULL; /* array of rdrop-values */
  P->ckeep=NULL; /* array of c-values (cols remaining) */
  P->cdrop=NULL; /* array of cdrop-values */
  P->rank=rank; /* rank of pca */
  P->tolerance=tolerance; /* tolerance used to assess relative-error in power-iteration */
  P->A_rpop_j_total=NULL; /* output array (size iteration_num) */
  P->A_cpop_j=NULL; /* output array (size iteration_num) */
  P->Irem=NULL; /* output array (size iteration_num) */
  P->U_=NULL; /* left singular-vectors (size P->D->A_nrows_total-x-P->rank-x-P->iteration_num) */
  P->V_=NULL; /* right singular_vectors (size P->D->A_ncols-x-P->rank-x-P->iteration_num) */
  if (V_fname!=NULL){ sprintf(P->V_name,V_fname);} else{ sprintf(P->V_name,"\0");}
  P->AnV_=NULL; /* projection of V_ onto An */
  P->ZnV_=NULL; /* projection of V_ onto Zn */
  P->b_mlt=b_mlt; /* precision used for xcalc */
  if (P->out_xdrop_name!=NULL &&  strcmp(P->out_xdrop_name,"\0")){ if (verbose>1){ printf(" %% calling xcc_out_xdrop_load\n");} xcc_out_xdrop_load(P->out_xdrop_name,&(P->out_xdrop_nlines),&(P->out_xdrop));}
  if (P->out_xdrop_name==NULL || !strcmp(P->out_xdrop_name,"\0")){ 
    if (P->iteration_num>1 || P->iteration_min>0 || P->iteration_max>0){
      if (verbose>1){ printf(" %% calling xcc_out_xdrop_init\n");} 
      xcc_out_xdrop_init(P->D->A_rpop_b_total,P->D->QR_index_global_mr_a,P->D->A_cpop_b,P->D->QC_index_local_mc_a,&(P->out_xdrop_nlines),&(P->out_xdrop));
      /* if (P->iteration_num>1 || P->iteration_min>0 || P->iteration_max>0){ } */}
    if (P->iteration_num==1 && P->iteration_min==0 && P->iteration_max==0){
      if (verbose>1){ printf(" %% calling xcc_out_xdrop_ini2\n");} 
      xcc_out_xdrop_ini2(P->D->A_rpop_b_total,P->D->QR_index_global_mr_a,P->D->A_cpop_b,P->D->QC_index_local_mc_a,&(P->out_xdrop_nlines),&(P->out_xdrop));
      /* if (P->iteration_num==1 && P->iteration_min==0 && P->iteration_max==0){ } */}
    /* if (P->out_xdrop_name==NULL || !strcmp(P->out_xdrop_name,"\0")){ } */}
  if (verbose>1){ printf(" %% calling xcc_out_xdrop_sort\n");}
  xcc_out_xdrop_sort(P->out_xdrop_nlines,P->out_xdrop,&(P->out_xdrop_nrows),&(P->mr_index_sort),&(P->out_xdrop_ncols),&(P->mc_index_sort));
  /* P->U_ = (double *) wkspace_all0c(P->D->A_nrows_total*P->rank*P->iteration_num*sizeof(double)); */
  P->U_ = NULL;
  if (P->V_name==NULL || !strcmp(P->V_name,"\0")){ 
    if (verbose>1){ printf(" %% allocating P->V_\n");}
    P->A_rpop_j_total = (int *) wkspace_all0c(P->iteration_num*sizeof(int));
    P->A_cpop_j = (int *) wkspace_all0c(P->iteration_num*sizeof(int));
    P->Irem = (int *) wkspace_all0c(P->iteration_num*sizeof(int));
    P->V_ = (double *) wkspace_all0c(P->D->A_ncols*P->rank*P->iteration_num*sizeof(double));
    /* if (P->V_name==NULL || !strcmp(P->V_name,"\0")){ } */}
  if (P->V_name!=NULL &&  strcmp(P->V_name,"\0")){ 
    if (verbose>1){ printf(" %% reading P->V_\n");}
    mda_read_d3_r8(P->V_name,NULL,&(P->rank),&(P->iteration_num),NULL);
    if (verbose>2){ printf(" %% P->V_name %s P->rank %d P->iteration_num %d\n",P->V_name,P->rank,P->iteration_num);}
    P->A_rpop_j_total = (int *) wkspace_all0c(P->iteration_num*sizeof(int));
    P->A_cpop_j = (int *) wkspace_all0c(P->iteration_num*sizeof(int));
    P->Irem = (int *) wkspace_all0c(P->iteration_num*sizeof(int));
    P->V_ = (double *) wkspace_all0c(P->D->A_ncols*P->rank*P->iteration_num*sizeof(double));
    mda_read_d3_r8(P->V_name,NULL,NULL,NULL,P->V_);
    P->AnV_ = (double *) wkspace_all0c(P->D->A_nrows_total*P->rank*P->iteration_num*sizeof(double));
    P->ZnV_ = (double *) wkspace_all0c(P->D->Z_nrows_total*P->rank*P->iteration_num*sizeof(double));
    /* if (P->V_name!=NULL &&  strcmp(P->V_name,"\0")){ } */}
  /* Left singular vector U minimizes:
     (1/A_nrows) * (AnAtU/A_ncols - YnYtU/Y_ncols).^2 - (1/Z_nrows) * (ZnAtU/A_nrows - WnYtU/Y_ncols).^2,
     implying that U is an eigenvector of:
     + (1/A_nrows) * (AnAtAnAt/A_ncols^2 - AnAtYnYt/A_ncols/Y_ncols - YnYtAnAt/Y_ncols/A_ncols + YnYtYnYt/Y_ncols^2)
     - (1/Z_nrows) * (AnZtZnAt/A_ncols^2 - AnZtWnYt/A_ncols/Y_ncols - YnWtZnAt/Y_ncols/A_ncols + YnWtWnYt/Y_ncols^2)
   */
  /* P_handle_make_excerpt_A_nrows_(P,&(P->lf_U_),&(P->M_U_)); */
  P->lf_U_ = NULL; P->M_U_ = NULL;
  /* Right singular vector V minimizes:
     (1/A_ncols) * (AtAnV/A_nrows - ZtZnV/Z_nrows).^2 - (1/Y_ncols) * (YtAnV/A_nrows - WtZnV/Z_nrows).^2,
     implying that V is an eigenvector of:
     + (1/A_ncols) * (AtAnAtAn/A_nrows^2 - AtAnZtZn/Z_nrows/A_nrows - ZtZnAtAn/A_nrows/Z_nrows + ZtZnZtZn/Z_nrows^2)
     - (1/Y_ncols) * (AtYnYtAn/A_nrows^2 - AtYnWtZn/Z_nrows/A_nrows - ZtWnYtAn/A_nrows/Z_nrows + ZtWnWtZn/Z_nrows^2)
   */
  if (verbose>1){ printf(" %% allocating P->lf_V\n");}
  P_handle_make_excerpt_A_ncols(P,&(P->lf_V),&(P->M_V));
  P->lf_R1 = L_handle_make((unsigned long long int)P->rank*(unsigned long long int)P->rank);
  P->lf_R1->spacing_row = SPACING_a; P->lf_R1->row_stride = P->rank;
  P->lf_R1->spacing_col = SPACING_a; P->lf_R1->col_stride = P->rank;
  P->lf_R1->spacing_lyr = SPACING_a; P->lf_R1->lyr_stride = 0;
  if (verbose>1){ printf(" %% allocating P->lf_R1\n");}
  L_zero(P->lf_R1);
  P->lf_R2 = L_handle_make((unsigned long long int)P->rank*(unsigned long long int)P->rank);
  P->lf_R2->spacing_row = SPACING_a; P->lf_R2->row_stride = P->rank;
  P->lf_R2->spacing_col = SPACING_a; P->lf_R2->col_stride = P->rank;
  P->lf_R2->spacing_lyr = SPACING_a; P->lf_R2->lyr_stride = 0;
  L_zero(P->lf_R2);
  bmc1_p = wkspace_all0c(bsize(P->rank)); for (nr=0;nr<P->rank;nr++){ bset__on(bmc1_p,nr);}
  P->M_rank = M_handle_v_make(BITJ,P->rank,1,NULL,NULL,bmc1_p,NULL);
  M_mxset(P->M_rank,bmc1_p,NULL);
  /*
    Calculating V requires:
    AnV_ (one for each bin), as well as AtAnV (summed over all bins) and YtAnV (summed over all bins);
    ZnV_ (one for each bin), as well as ZtZnV (summed over all bins) and WtZnV (summed over all bins);
    Once we have AtAnV and ZtZnV, we can form AtAnV/A_nrows - ZtZnV/Z_nrows =: VA;
    then form AnVA_ (one for each bin) as well as AtAnVA (summed over all bins);
    also form ZnVA_ (one for each bin) as well as ZtZnVA (summed over all bins);
    then finally form AtAnVA/A_nrows - ZtZnVA/Z_nrows =: VAZ;
    Once we have YtAnV and WtZnV, we can form YtAnV/A_nrows - WtZnV/Z_nrows =: VY;
    then form YnVY_ (one for each bin) as well as AtYnVY (summed over all bins);
    also form WnVY_ (one for each bin) as well as ZtWnVY (summed over all bins);
    then finally form AtYnVY/A_nrows - ZtWnVY/Z_nrows =: VYW;
    Once we have VAZ and VYW we can form VAZ/A_ncols - VYW/Y_ncols =: V;
   */
  if (verbose>1){ printf(" %% allocating P->lf_AnV_, etc.\n");}
  P_handle_make_excerpt_A_nrows_(P,&(P->lf_AnV_),&(P->M_AnV_)); P_handle_make_excerpt_Z_nrows_(P,&(P->lf_ZnV_),&(P->M_ZnV_));
  P_handle_make_excerpt_A_ncols_(P,&(P->lf_AtAnV_),&(P->M_AtAnV_)); P_handle_make_excerpt_A_ncols(P,&(P->lf_AtAnV),&(P->M_AtAnV)); 
  P_handle_make_excerpt_Y_ncols_(P,&(P->lf_YtAnV_),&(P->M_YtAnV_)); P_handle_make_excerpt_Y_ncols(P,&(P->lf_YtAnV),&(P->M_YtAnV));
  P_handle_make_excerpt_A_ncols_(P,&(P->lf_ZtZnV_),&(P->M_ZtZnV_)); P_handle_make_excerpt_A_ncols(P,&(P->lf_ZtZnV),&(P->M_ZtZnV)); 
  P_handle_make_excerpt_Y_ncols_(P,&(P->lf_WtZnV_),&(P->M_WtZnV_)); P_handle_make_excerpt_Y_ncols(P,&(P->lf_WtZnV),&(P->M_WtZnV));
  P_handle_make_excerpt_A_ncols(P,&(P->lf_VA),&(P->M_VA));
  P_handle_make_excerpt_A_nrows_(P,&(P->lf_AnVA_),&(P->M_AnVA_)); P_handle_make_excerpt_Z_nrows_(P,&(P->lf_ZnVA_),&(P->M_ZnVA_));
  P_handle_make_excerpt_A_ncols_(P,&(P->lf_AtAnVA_),&(P->M_AtAnVA_)); P_handle_make_excerpt_A_ncols(P,&(P->lf_AtAnVA),&(P->M_AtAnVA));
  P_handle_make_excerpt_A_ncols_(P,&(P->lf_ZtZnVA_),&(P->M_ZtZnVA_)); P_handle_make_excerpt_A_ncols(P,&(P->lf_ZtZnVA),&(P->M_ZtZnVA));
  P_handle_make_excerpt_Y_ncols(P,&(P->lf_VY),&(P->M_VY));
  P_handle_make_excerpt_A_nrows_(P,&(P->lf_YnVY_),&(P->M_YnVY_)); P_handle_make_excerpt_Z_nrows_(P,&(P->lf_WnVY_),&(P->M_WnVY_));
  P_handle_make_excerpt_A_ncols_(P,&(P->lf_AtYnVY_),&(P->M_AtYnVY_)); P_handle_make_excerpt_A_ncols(P,&(P->lf_AtYnVY),&(P->M_AtYnVY));
  P_handle_make_excerpt_A_ncols_(P,&(P->lf_ZtWnVY_),&(P->M_ZtWnVY_)); P_handle_make_excerpt_A_ncols(P,&(P->lf_ZtWnVY),&(P->M_ZtWnVY));
  P_handle_make_excerpt_A_ncols(P,&(P->lf_VAZ),&(P->M_VAZ));
  P_handle_make_excerpt_A_ncols(P,&(P->lf_VYW),&(P->M_VYW));
  P->nx=0;
  /* used for error checking */
  if (verbose>1){ printf(" %% initializing error checking structures\n");}
  P->lf_At = NULL; P->lf_Zt = NULL; P->lf_Yt = NULL; P->lf_Wt = NULL;
  P->lf_An = NULL; P->lf_Zn = NULL; P->lf_Yn = NULL; P->lf_Wn = NULL;
  P->lf_AtAn = NULL; P->lf_ZtZn = NULL; P->lf_AtYn = NULL; P->lf_ZtWn = NULL;
  P->lf_AtAnAtAn = NULL; P->lf_AtAnZtZn = NULL; P->lf_ZtZnAtAn = NULL; P->lf_ZtZnZtZn = NULL;
  P->lf_AtYnYtAn = NULL; P->lf_AtYnWtZn = NULL; P->lf_ZtWnYtAn = NULL; P->lf_ZtWnWtZn = NULL;
  P->lf_S = NULL; 
  P_handle_make_excerpt_A_ncols(P,&(P->lf_Vn),NULL); 
  P_handle_make_excerpt_A_ncols(P,&(P->lf_Vt),NULL);
  P_handle_make_excerpt_A_ncols(P,&(P->lf_Vx),NULL);
  P->lf_VtAt = NULL; P->lf_VtZt = NULL;
  P->lf_VtAtAn = NULL; P->lf_VtZtZn = NULL; P->lf_VtAtYn = NULL; P->lf_VtZtWn = NULL;
  P->lf_VAt = NULL;
  P->lf_VAtAt = NULL; P->lf_VAtZt = NULL;
  P->lf_VAtAtAn = NULL; P->lf_VAtZtZn = NULL; 
  P->lf_VAZt = NULL;
  P->lf_VYt = NULL;
  P->lf_VYtYt = NULL; P->lf_VYtWt = NULL;
  P->lf_VYtYtAn = NULL; P->lf_VYtWtZn = NULL; 
  P->lf_VYWt = NULL;
  if (verbose){ printf(" %% [finished P_handle_make]\n");}
  return P;
}

void P_xdrop_init(struct dcc_ajdk *D,struct P_handle *P)
{
  dcc_out_xdrop_lkp(D,P->out_xdrop_nrows,P->mr_index_sort,&(P->mr_index_local_nb),&(P->mr_index_local_mr));
  get_xdrop_array_sub(D->A_rpop_j_total,D->A_cpop_j,P->iteration_num,P->iteration_max,P->iteration_min,NULL,&(P->rdrop),&(P->cdrop),&(P->rkeep),&(P->ckeep));
}

void P_handle_printf(int verbose,struct P_handle *P,char *prefix)
{
  int nl=0,nx=0;
  char tmpchar[FNAMESIZE];
  if (verbose>0){ printf("%s iteration_num %d iteration_max %d iteration_min %d rank %d b_mlt %d\n",prefix,P->iteration_num,P->iteration_max,P->iteration_min,P->rank,P->b_mlt);}
  if (verbose>0){ printf("%s nx %d\n",prefix,P->nx);}
  if (verbose>0){ printf("%s out_xdrop_name: %s\n",prefix,P->out_xdrop_name);}
  if (verbose>0){ printf("%s nlines %d, nrows %d, ncols %d\n",prefix,P->out_xdrop_nlines,P->out_xdrop_nrows,P->out_xdrop_ncols);}
  if (verbose>0){
    if (P->rkeep!=NULL){ sprintf(tmpchar,"%s rkeep:",prefix); raprintf(P->rkeep,"int",1,P->iteration_num,tmpchar);}
    if (P->rdrop!=NULL){ sprintf(tmpchar,"%s rdrop:",prefix); raprintf(P->rdrop,"int",1,P->iteration_num,tmpchar);}
    if (P->ckeep!=NULL){ sprintf(tmpchar,"%s ckeep:",prefix); raprintf(P->ckeep,"int",1,P->iteration_num,tmpchar);}
    if (P->cdrop!=NULL){ sprintf(tmpchar,"%s cdrop:",prefix); raprintf(P->cdrop,"int",1,P->iteration_num,tmpchar);}
    /* if (verbose>0){ } */}
  if (verbose>1){
    if (P->A_rpop_j_total!=NULL){ sprintf(tmpchar,"%s A_rpop_j_total:",prefix); raprintf(P->A_rpop_j_total,"int",1,P->iteration_num,tmpchar);}
    if (P->A_cpop_j!=NULL){ sprintf(tmpchar,"%s A_cpop_j:",prefix); raprintf(P->A_cpop_j,"int",1,P->iteration_num,tmpchar);}
    if (P->Irem!=NULL){ sprintf(tmpchar,"%s Irem:",prefix); raprintf(P->Irem,"int",1,P->iteration_num,tmpchar);}
    for (nx=0;nx<P->iteration_num;nx++){
      printf(" %% nx %d\n",nx);
      for (nl=0;nl<P->rank;nl++){
    	/* if (P->U_!=NULL){ sprintf(tmpchar,"%s U_[%d]:",prefix,nl); raprintf(&(P->U_[0 + nl*P->D->A_nrows_total + nx*P->D->A_nrows_total*P->rank]),"double",1,P->D->A_nrows_total,tmpchar);} */
    	if (P->V_!=NULL){ sprintf(tmpchar,"%s V_[%d]:",prefix,nl); raprintf(&(P->V_[0 + nl*P->D->A_ncols + nx*P->D->A_ncols*P->rank]),"double",1,P->D->A_ncols,tmpchar);}
    	/* for (nl=0;nl<P->rank;nl++){ } */}
      /* for (nx=0;nx<P->iteration_num;nx++){ } */}
    /* if (verbose>2){ } */}
}

void P_handle_dmp(struct P_handle *P)
{
  char prefix[FNAMESIZE];
  char fname[FNAMESIZE];
  if (P->infix==NULL || !strcmp(P->infix,"\0")){ sprintf(prefix,"pca_k%d_B%d",P->rank,P->b_mlt);}
  if (P->infix!=NULL &&  strcmp(P->infix,"\0")){ sprintf(prefix,"%s_k%d_B%d",P->infix,P->rank,P->b_mlt);}
  sprintf(fname,"%s/%s_A_rpop_j_total.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_i4(fname,P->iteration_num,1,P->A_rpop_j_total);
  sprintf(fname,"%s/%s_A_cpop_j.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_i4(fname,P->iteration_num,1,P->A_cpop_j);
  sprintf(fname,"%s/%s_Irem.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_i4(fname,P->iteration_num,1,P->Irem);
  /* sprintf(fname,"%s/%s_U_.mda",GLOBAL_DIR_NAME,prefix); mda_write_d3_r8(fname,P->D->A_nrows_total,P->rank,P->iteration_num,P->U_); */
  sprintf(fname,"%s/%s_V_.mda",GLOBAL_DIR_NAME,prefix); mda_write_d3_r8(fname,P->D->A_ncols,P->rank,P->iteration_num,P->V_);
}

void P_handle_proj_dmp(struct P_handle *P)
{
  char prefix[FNAMESIZE];
  char fname[FNAMESIZE];
  if (P->infix==NULL || !strcmp(P->infix,"\0")){ sprintf(prefix,"pca_proj_k%d_B%d",P->rank,P->b_mlt);}
  if (P->infix!=NULL &&  strcmp(P->infix,"\0")){ sprintf(prefix,"%s_k%d_B%d",P->infix,P->rank,P->b_mlt);}
   sprintf(fname,"%s/%s_A_rpop_j_total.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_i4(fname,P->iteration_num,1,P->A_rpop_j_total);
  sprintf(fname,"%s/%s_A_cpop_j.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_i4(fname,P->iteration_num,1,P->A_cpop_j);
  sprintf(fname,"%s/%s_Irem.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_i4(fname,P->iteration_num,1,P->Irem);
  sprintf(fname,"%s/%s_AnV_.mda",GLOBAL_DIR_NAME,prefix); mda_write_d3_r8(fname,P->D->A_nrows_total,P->rank,P->iteration_num,P->AnV_);
  sprintf(fname,"%s/%s_ZnV_.mda",GLOBAL_DIR_NAME,prefix); mda_write_d3_r8(fname,P->D->Z_nrows_total,P->rank,P->iteration_num,P->ZnV_);
}
