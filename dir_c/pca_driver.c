#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void pca_srt(struct P_handle *P)
{
  /* Copies all D->A_rpop_b_total==P->out_xdrop_nrows entries from P->mr_index_sort (i.e., P->mr_index_local_nb and P->mr_index_local_mr) ;
     into D->QR_index_local_nb and D->QR_index_local_mr_a, respectively. ;
     Copies all D->A_cpop_b==P->out_xdrop_ncols entries from P->mc_index_sort ;
     into D->QC_index_local_mc_a. ;
  */
  int verbose=0;
  struct dcc_ajdk *D=P->D;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL;
  int nc=0,nr=0;
  if (verbose>1){ printf(" %% [entering pca_srt]\n");}
  if (D->A_rpop_b_total!=P->out_xdrop_nrows){ printf(" %% Note: D->A_rpop_b_total %d != P->out_xdrop_nrows %d\n",D->A_rpop_b_total,P->out_xdrop_nrows);}
  if (P->mc_index_sort!=NULL){
    for (nc=0;nc<P->out_xdrop_ncols;nc++){ D->QC_index_local_mc_a[nc] = P->mc_index_sort[nc];}
    /* if (P->mc_index_sort!=NULL){ } */}
  if (P->mr_index_local_nb!=NULL && P->mr_index_local_mr!=NULL){
    for (nr=0;nr<P->out_xdrop_nrows;nr++){
      D->QR_index_local_mr_a[nr] = P->mr_index_local_mr[nr];
      D->QR_index_local_nb[nr] = P->mr_index_local_nb[nr];
      /* for (nr=0;nr<D->P->out_xdrop_nrows;nr++){ } */}
    /* if (P->mc_index_sort!=NULL){ } */};
  if (verbose>1){ printf(" %% finished copying mr_index_local_nb,mr_index_local_mr and mc_index_sort\n");}
  if (verbose>1){
    raprintf(D->QR_index_local_mr_a,   "int",1,P->out_xdrop_nrows," %% D->QR_index_local_mr_a pos :");
    raprintf(  D->QR_index_local_nb,   "int",1,P->out_xdrop_nrows," %%   D->QR_index_local_nb pos :");
    raprintf(D->QC_index_local_mc_a,   "int",1,P->out_xdrop_ncols," %% D->QC_index_local_mc_a pos :");
    /* if (verbose>1){ } */}
  if (verbose>1){ printf(" %% [finished pca_srt]\n");}
}

void pca_mxA(struct P_handle *P,int rdrop,int cdrop)
{
  /* reading rdrop and cdrop from input: ;
     copies D->A_bmc_b to D->A_umc_j ; 
     define D->A_umc_j_rmv and D->A_umc_j_rtn ; removing first cdrop entries of D->QC_index_local_mc_a ;
     copies E->A_bmr_b to E->A_umr_j ; 
     define E->A_umr_j_rmv and E->A_umr_j_rtn ; removing first rdrop entries of D->QR_index_local_mr_a ; using D->QR_index_local_nb to index nb ;
     call dcc_sumscores_mxB ; copying all instances of umc and umr to associated bmc and bmr respectively ;
   */
  int verbose=0; 
  struct dcc_ajdk *D=P->D;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_single *E=NULL;
  int nb1=0,nr=0,nc=0,nl=0;
  for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_b,nc);} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = D->A_umc_j[nc]; D->A_umc_j_rmv[nc] = 0;}
  for (nc=0;nc<minimum(cdrop,D->A_cpop_b);nc++){ D->A_umc_j_rmv[D->QC_index_local_mc_a[nc]]=1; D->A_umc_j_rtn[D->QC_index_local_mc_a[nc]]=0;}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_b,nr);} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = E->A_umr_j[nr]; E->A_umr_j_rmv[nr] = 0;}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nr=0;nr<minimum(rdrop,D->A_rpop_b_total);nr++){ E_[D->QR_index_local_nb[nr]]->A_umr_j_rmv[D->QR_index_local_mr_a[nr]]=1; E_[D->QR_index_local_nb[nr]]->A_umr_j_rtn[D->QR_index_local_mr_a[nr]]=0;}
  if (verbose>1){
    if (verbose>1){ raprintf(D->A_umc_j    ,"char",1,D->A_ncols," %% D->A_umc_j    : ");}
    if (verbose>1){ raprintf(D->A_umc_j_rmv,"char",1,D->A_ncols," %% D->A_umc_j_rmv: ");}
    if (verbose>1){ raprintf(D->A_umc_j_rtn,"char",1,D->A_ncols," %% D->A_umc_j_rtn: ");}
    for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j    : ",nb1); raprintf(E->A_umr_j    ,"char",1,E->A_nrows,D->tmpAnchar);
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j_rmv: ",nb1); raprintf(E->A_umr_j_rmv,"char",1,E->A_nrows,D->tmpAnchar);
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j_rtn: ",nb1); raprintf(E->A_umr_j_rtn,"char",1,E->A_nrows,D->tmpAnchar);
      /* for (nb1=0;nb1<nbins;nb1++){ } */}
    /* if (verbose>1){ } */}
  if (verbose>1){ printf(" %% calling dcc_sumscores_mxB\n");}
  dcc_sumscores_mxB(D); 
  if (verbose>1){ printf(" %% finished dcc_sumscores_mxB\n");}
  if (verbose>1){
    sprintf(D->tmpAnchar," %%%% D->A_bmc_j    : "); bprintf(D->A_bmc_j    ,D->bitj,1,D->A_ncols,D->tmpAnchar);
    sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rmv: "); bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols,D->tmpAnchar);
    sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rtn: "); bprintf(D->A_bmc_j_rtn,D->bitj,1,D->A_ncols,D->tmpAnchar);
    for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j    : ",nb1); bprintf(E->A_bmr_j    ,D->bitj,1,E->A_nrows,D->tmpAnchar);
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j_rmv: ",nb1); bprintf(E->A_bmr_j_rmv,D->bitj,1,E->A_nrows,D->tmpAnchar);
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j_rtn: ",nb1); bprintf(E->A_bmr_j_rtn,D->bitj,1,E->A_nrows,D->tmpAnchar);
      /* for (nb1=0;nb1<nbins;nb1++){ } */}
    /* if (verbose>1){ } */}
  if (verbose>1){ printf(" %% [finished pca_mxA]\n");}
}

void pca_printf(int verbose,struct P_handle *P)
{
  int nb=0,nl=0,na_a=0;
  char strtmp[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering pca_printf]\n");}
  for (nb=0;nb<P->D->nbins;nb++){
    for (na_a=0;na_a<P->D->A_ncols;na_a++){
      if (bget____(P->D->E_[nb]->M_An->mc_j,na_a)!=bget____(P->D->E_[0]->M_An->mc_j,na_a)){
	printf(" %% Warning! mc_j mismatch %d vs %d at na_a %d nb %d in pca_printf\n",bget____(P->D->E_[nb]->M_An->mc_j,na_a),bget____(P->D->E_[0]->M_An->mc_j,na_a),na_a,nb);
	/* if (bget____(P->D->E_[nb]->M_An->mc_j,na_a)!=bget____(P->D->E_[0]->M_An->mc_j,na_a)){ } */}
      /* for (na_a=0;na_a<P->D->A_ncols;na_a++){ } */}
    /*  for (nb=0;nb<P->D->nbins;nb++){ } */}
  if (verbose>2){ 
    for (nb=0;nb<P->D->nbins;nb++){ 
      sprintf(strtmp," %%%% [rpop_j: %d < %d] P->D->E_[%d]->M_An->mr_j: ",P->D->E_[nb]->M_An->rpop_j,P->D->A_rpop_j_total,nb);
      bprintf(P->D->E_[nb]->M_An->mr_j,P->D->E_[nb]->M_An->bitj,1,P->D->E_[nb]->M_An->nrows,strtmp); 
      /*  for (nb=0;nb<P->D->nbins;nb++){ } */}
    nb=0;
    sprintf(strtmp," %%%% [cpop_j: %d < %d] P->D->E_[%d]->M_An->mc_j: ",P->D->E_[nb]->M_An->cpop_j,P->D->A_cpop_j,nb);
    bprintf(P->D->E_[nb]->M_An->mc_j,P->D->E_[nb]->M_An->bitj,1,P->D->E_[nb]->M_An->ncols,strtmp); 
    /* if (verbose>2){ } */}
  if (verbose>0){
    lfprintf(P->lf_R1," %% P->lf_R1: ");
    /* if (verbose>0){ } */}
  if (verbose>1){
    printf(" %% P->lf_V length %d row_stride %d col_stride %d lyr_stride %d\n",P->lf_V->length,P->lf_V->row_stride,P->lf_V->col_stride,P->lf_V->lyr_stride);
    lfprintf(P->lf_V," %% P->lf_V: ");
    /* if (verbose>1){ } */}
  if (verbose>1){ printf(" %% [finished pca_printf]\n");}
}

void dra_dot(int nrows,double *a,double *b,double *r)
{
  int nr=0;
  if (r!=NULL){
    *r = 0;
    for (nr=0;nr<nrows;nr++){ *r += a[nr]*b[nr];} 
     /* if (r!=NULL){ } */}
}

void dra_norm(int nrows,double *a,double *r)
{
  if (r!=NULL){ dra_dot(nrows,a,a,r); *r = sqrt(*r);}
}

void dra_graham_schmidt(int nrows,int ncols,double *A,double *R)
{
  /* orthonormalizes columns of A in place, returning matrix R */
  int verbose=0;
  int nr=0,nc1=0,nc2=0;
  double r=0;
  if (verbose>1){ printf(" %% [entering dra_graham_schmidt] nrows %d ncols %d\n",nrows,ncols);}
  for (nc1=0;nc1<ncols;nc1++){
    dra_norm(nrows,&(A[0+nc1*nrows]),&r);
    if (R!=NULL){ R[nc1+nc1*ncols] = r;}
    if (r<1e-12){ printf(" %% Warning! small norm in dra_graham_schmidt\n");}
    dra_times(&(A[0+nc1*nrows]),nrows,1.0/maximum(1e-12,r));
    for (nc2=nc1+1;nc2<ncols;nc2++){
      dra_dot(nrows,&(A[0+nc1*nrows]),&(A[0+nc2*nrows]),&r);
      if (R!=NULL){ R[nc1+nc2*ncols]=r;}
      dra_plustimesequals(&(A[0+nc2*nrows]),nrows,&(A[0+nc1*nrows]),-r);
      /* for (nc2=nc1+1;nc2<ncols;nc2++){ } */}
    /* for (nc1=0;nc1<ncols;nc1++){ } */}
  if (verbose>1){
    for (nc1=0;nc1<ncols;nc1++){
      dra_norm(nrows,&(A[0+nc1*nrows]),&r);
      printf(" %% column %d norm %f\n",nc1,r);
      /* for (nc1=0;nc1<ncols;nc1++){ } */}
    /* if (verbose>1){ } */}
  if (verbose>1){ printf(" %% [finished dra_graham_schmidt] nrows %d ncols %d\n",nrows,ncols);}
}

void pca_stage_a0(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  int na_a=0,na_b=0,na_j=0,nr=0,nb=0;
  struct M_handle *M_An=NULL;
  double dtmp=0;
  if (verbose>1){ printf(" %% [entering pca_stage_a0]\n");}
  L_zero(P->lf_V);
  M_An = P->D->E_[nb]->M_An;
  na_j=0;
  while (na_j<M_An->cpop_j){
    na_a = M_An->n_a_[na_j]; na_b = M_An->n_b_[na_j];
    for (nr=0;nr<P->rank;nr++){
      /* dtmp = (nr==0?1:0) + (nr>0?2*rand01-1:0); */
      dtmp = (nr==0?1:0) + (nr>0? pow(-1,na_a%nr)*na_a:0);
      L2_set(P->lf_V,na_j,na_b,na_a,0,0,nr,dtmp);
      /* for (nr=0;nr<P->rank;nr++){ } */}
    na_j++; /* while (na_j<M_An->cpop_j){ } */}
  if (verbose>1){ printf(" %% [finished pca_stage_a0]\n");}
}

void pca_stage_a1(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  if (verbose>1){ printf(" %% [entering pca_stage_a1]\n");}
  dra_graham_schmidt(P->lf_V->row_stride,P->lf_V->col_stride,P->lf_V->lf,P->lf_R1->lf);
  if (verbose>1){ printf(" %% [finished pca_stage_a1]\n");}
}

void pca_stage_a0b(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  int na_a=0,na_b=0,na_j=0,nr=0,nb=0;
  struct M_handle *M_An=NULL;
  double dtmp=0;
  if (verbose>1){ printf(" %% [entering pca_stage_a0b]\n");}
  L_zero(P->lf_Vt);
  M_An = P->D->E_[nb]->M_An;
  na_j=0;
  while (na_j<M_An->cpop_j){
    na_a = M_An->n_a_[na_j]; na_b = M_An->n_b_[na_j];
    for (nr=0;nr<P->rank;nr++){
      /* dtmp = *(L2_get(P->lf_V,na_j,na_b,na_a,0,0,nr)); */
      dtmp = (nr==0?1:0) + (nr>0? pow(-1,na_a%nr)*na_a:0);
      L2_set(P->lf_Vt,na_j,na_b,na_a,0,0,nr,dtmp);
      /* for (nr=0;nr<P->rank;nr++){ } */}
    na_j++; /* while (na_j<M_An->cpop_j){ } */}
  if (verbose>1){ printf(" %% [finished pca_stage_a0b]\n");}
}

void pca_stage_a1b(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  if (verbose>1){ printf(" %% [entering pca_stage_a1b]\n");}
  dra_graham_schmidt(P->lf_Vt->row_stride,P->lf_Vt->col_stride,P->lf_Vt->lf,P->lf_R2->lf);
  if (verbose>1){ printf(" %% [finished pca_stage_a1b]\n");}
}

void pca_stage_b(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  struct dcc_ajdk *D = P->D; 
  struct dcc_single **E_ = D->E_;
  int nbins = D->nbins; 
  int nb=0; struct dcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  int ns_j=0,ns_b=0,ns_a=0;
  double dtmp=0;
  if (verbose>1){ printf(" %% [entering pca_stage_b]\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  if (D->A_cbother){ 
    GLOBAL_pthread_tic(); 
    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),P->D->A_bmc_j,P->D->A_bmc_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_V,P->lf_V->lf,&(P->D->A_ncols),&(P->M_rank->nrows),&(P->M_V),&(P->b_mlt),(addressable_0));
    GLOBAL_pthread_toc(); /* if bother */}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_An,P->M_V,P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_AnV_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Zn,P->M_V,P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_ZnV_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% An_Xn_ww : ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      printf(" %% nb %d: P->lf_AnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_AnV_[nb]->length,P->lf_AnV_[nb]->row_stride,P->lf_AnV_[nb]->col_stride,P->lf_AnV_[nb]->lyr_stride);
      lfprintf(P->lf_AnV_[nb]," %% P->lf_AnV_[nb]->lf: ");
      printf(" %% nb %d: P->lf_ZnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_ZnV_[nb]->length,P->lf_ZnV_[nb]->row_stride,P->lf_ZnV_[nb]->col_stride,P->lf_ZnV_[nb]->lyr_stride);
      lfprintf(P->lf_ZnV_[nb]," %% P->lf_ZnV_[nb]->lf: ");
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->A_bmr_j,E->A_bmr_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_AnV_[nb],P->lf_AnV_[nb]->lf,&(E->A_nrows),&(P->M_rank->nrows),&(P->M_AnV_[nb]),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->Z_bmr_j,E->Z_bmr_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_ZnV_[nb],P->lf_ZnV_[nb]->lf,&(E->Z_nrows),&(P->M_rank->nrows),&(P->M_ZnV_[nb]),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc trm_0: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_At_T_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_At,E->M_Tt,P->M_AnV_[nb],P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_AtAnV_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic();
      wrap_At_T_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Yt,E->M_Tt,P->M_AnV_[nb],P->M_rank,P->D->Y_ajdk,(addressable_0),&(P->lf_YtAnV_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_At_T_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Zt,E->M_St,P->M_ZnV_[nb],P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_ZtZnV_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->Y_cbother){
      GLOBAL_pthread_tic();
      wrap_At_T_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Wt,E->M_St,P->M_ZnV_[nb],P->M_rank,P->D->Y_ajdk,(addressable_0),&(P->lf_WtZnV_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% At_T_Xn_ww : ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->A_rbother && D->A_cbother){ 
	printf(" %% nb %d: P->lf_AtAnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_AtAnV_[nb]->length,P->lf_AtAnV_[nb]->row_stride,P->lf_AtAnV_[nb]->col_stride,P->lf_AtAnV_[nb]->lyr_stride);
	lfprintf(P->lf_AtAnV_[nb]," %% P->lf_AtAnV_[nb]->lf: "); /* if bother */}
      if (E->A_rbother && D->Y_cbother){ 
	printf(" %% nb %d: P->lf_YtAnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_YtAnV_[nb]->length,P->lf_YtAnV_[nb]->row_stride,P->lf_YtAnV_[nb]->col_stride,P->lf_YtAnV_[nb]->lyr_stride);
	lfprintf(P->lf_YtAnV_[nb]," %% P->lf_YtAnV_[nb]->lf: "); /* if bother */}
      if (E->Z_rbother && D->A_cbother){ 
	printf(" %% nb %d: P->lf_ZtZnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_ZtZnV_[nb]->length,P->lf_ZtZnV_[nb]->row_stride,P->lf_ZtZnV_[nb]->col_stride,P->lf_ZtZnV_[nb]->lyr_stride);
	lfprintf(P->lf_ZtZnV_[nb]," %% P->lf_ZtZnV_[nb]->lf: "); /* if bother */}
      if (E->Z_rbother && D->Y_cbother){ 
	printf(" %% nb %d: P->lf_WtZnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_WtZnV_[nb]->length,P->lf_WtZnV_[nb]->row_stride,P->lf_WtZnV_[nb]->col_stride,P->lf_WtZnV_[nb]->lyr_stride);
	lfprintf(P->lf_WtZnV_[nb]," %% P->lf_WtZnV_[nb]->lf: "); /* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  L_zero(P->lf_VA); L_zero(P->lf_VY);
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    ns_j = 0; ns_a = 0; ns_b = 0; 
    if (E->A_rbother && D->A_cbother){ 
      dtmp = maximum(1,P->D->A_rpop_j_total);
      dra_plustimesequals(L2_get(P->lf_VA,0,0,0,0,0,0),P->lf_VA->length,L3_get(P->lf_AtAnV_[nb],0,0,0,0,0,0,ns_j,ns_b,ns_a),+1.0/dtmp); /* if bother */}
    if (E->A_rbother && D->Y_cbother){ 
      dtmp = maximum(1,P->D->A_rpop_j_total);
      dra_plustimesequals(L2_get(P->lf_VY,0,0,0,0,0,0),P->lf_VY->length,L3_get(P->lf_YtAnV_[nb],0,0,0,0,0,0,ns_j,ns_b,ns_a),+1.0/dtmp); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      dtmp = maximum(1,P->D->Z_rpop_j_total);
      dra_plustimesequals(L2_get(P->lf_VA,0,0,0,0,0,0),P->lf_VA->length,L3_get(P->lf_ZtZnV_[nb],0,0,0,0,0,0,ns_j,ns_b,ns_a),-1.0/dtmp); /* if bother */}
    if (E->Z_rbother && D->Y_cbother){ 
      maximum(1,P->D->Z_rpop_j_total);
      dra_plustimesequals(L2_get(P->lf_VY,0,0,0,0,0,0),P->lf_VY->length,L3_get(P->lf_WtZnV_[nb],0,0,0,0,0,0,ns_j,ns_b,ns_a),-1.0/dtmp); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
      if (E->A_rbother && D->A_cbother){ 
	printf(" %% P->lf_VA length %d row_stride %d col_stride %d lyr_stride %d\n",P->lf_VA->length,P->lf_VA->row_stride,P->lf_VA->col_stride,P->lf_VA->lyr_stride);
	lfprintf(P->lf_VA," %% P->lf_VA->lf: "); /* if bother */}
      if (E->A_rbother && D->Y_cbother){ 
	printf(" %% P->lf_VY length %d row_stride %d col_stride %d lyr_stride %d\n",P->lf_VY->length,P->lf_VY->row_stride,P->lf_VY->col_stride,P->lf_VY->lyr_stride);
	lfprintf(P->lf_VY," %% P->lf_VY->lf: "); /* if bother */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/  
  if (verbose>1){ printf(" %% [finished pca_stage_b]\n");}
}

void pca_stage_c1(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  struct dcc_ajdk *D = P->D; 
  struct dcc_single **E_ = D->E_;
  int nbins = D->nbins; 
  int nb=0; struct dcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  int ns_j=0,ns_b=0,ns_a=0;
  double dtmp=0;
  if (verbose>1){ printf(" %% [entering pca_stage_c1]\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  if (D->A_cbother){ 
    GLOBAL_pthread_tic(); 
    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),P->D->A_bmc_j,P->D->A_bmc_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_VA,P->lf_VA->lf,&(P->D->A_ncols),&(P->M_rank->nrows),&(P->M_VA),&(P->b_mlt),(addressable_0));
    GLOBAL_pthread_toc(); /* if bother */}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_An,P->M_VA,P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_AnVA_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Zn,P->M_VA,P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_ZnVA_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% An_Xn_ww : ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      printf(" %% nb %d: P->lf_AnVA_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_AnVA_[nb]->length,P->lf_AnVA_[nb]->row_stride,P->lf_AnVA_[nb]->col_stride,P->lf_AnVA_[nb]->lyr_stride);
      lfprintf(P->lf_AnVA_[nb]," %% P->lf_AnVA_[nb]->lf: ");
      printf(" %% nb %d: P->lf_ZnVA_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_ZnVA_[nb]->length,P->lf_ZnVA_[nb]->row_stride,P->lf_ZnVA_[nb]->col_stride,P->lf_ZnVA_[nb]->lyr_stride);
      lfprintf(P->lf_ZnVA_[nb]," %% P->lf_ZnVA_[nb]->lf: ");
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->A_bmr_j,E->A_bmr_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_AnVA_[nb],P->lf_AnVA_[nb]->lf,&(E->A_nrows),&(P->M_rank->nrows),&(P->M_AnVA_[nb]),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->Z_bmr_j,E->Z_bmr_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_ZnVA_[nb],P->lf_ZnVA_[nb]->lf,&(E->Z_nrows),&(P->M_rank->nrows),&(P->M_ZnVA_[nb]),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc trm_0: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_At_T_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_At,E->M_Tt,P->M_AnVA_[nb],P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_AtAnVA_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){
      GLOBAL_pthread_tic();
      wrap_At_T_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Zt,E->M_St,P->M_ZnVA_[nb],P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_ZtZnVA_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% At_T_Xn_ww : ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->A_rbother && D->A_cbother){ 
	printf(" %% nb %d: P->lf_AtAnVA_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_AtAnVA_[nb]->length,P->lf_AtAnVA_[nb]->row_stride,P->lf_AtAnVA_[nb]->col_stride,P->lf_AtAnVA_[nb]->lyr_stride);
	lfprintf(P->lf_AtAnVA_[nb]," %% P->lf_AtAnVA_[nb]->lf: "); /* if bother */}
      if (E->Z_rbother && D->A_cbother){ 
	printf(" %% nb %d: P->lf_ZtZnVA_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_ZtZnVA_[nb]->length,P->lf_ZtZnVA_[nb]->row_stride,P->lf_ZtZnVA_[nb]->col_stride,P->lf_ZtZnVA_[nb]->lyr_stride);
	lfprintf(P->lf_ZtZnVA_[nb]," %% P->lf_ZtZnVA_[nb]->lf: "); /* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  L_zero(P->lf_VAZ); 
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    ns_j = 0; ns_a = 0; ns_b = 0; 
    if (E->A_rbother && D->A_cbother){ 
      dtmp = maximum(1,P->D->A_rpop_j_total);
      dra_plustimesequals(L2_get(P->lf_VAZ,0,0,0,0,0,0),P->lf_VAZ->length,L3_get(P->lf_AtAnVA_[nb],0,0,0,0,0,0,ns_j,ns_b,ns_a),+1.0/dtmp); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      dtmp = maximum(1,P->D->Z_rpop_j_total);
      dra_plustimesequals(L2_get(P->lf_VAZ,0,0,0,0,0,0),P->lf_VAZ->length,L3_get(P->lf_ZtZnVA_[nb],0,0,0,0,0,0,ns_j,ns_b,ns_a),-1.0/dtmp); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
      if (E->A_rbother && D->A_cbother){ 
	printf(" %% P->lf_VAZ length %d row_stride %d col_stride %d lyr_stride %d\n",P->lf_VAZ->length,P->lf_VAZ->row_stride,P->lf_VAZ->col_stride,P->lf_VAZ->lyr_stride);
	lfprintf(P->lf_VAZ," %% P->lf_VAZ->lf: "); /* if bother */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/  
  if (verbose>1){ printf(" %% [finished pca_stage_c1]\n");}
}

void pca_stage_c2(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  struct dcc_ajdk *D = P->D; 
  struct dcc_single **E_ = D->E_;
  int nbins = D->nbins; 
  int nb=0; struct dcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  int ns_j=0,ns_b=0,ns_a=0;
  double dtmp=0;
  if (verbose>1){ printf(" %% [entering pca_stage_c2]\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  if (D->Y_cbother){ 
    GLOBAL_pthread_tic(); 
    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),P->D->Y_bmc_j,P->D->Y_bmc_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_VY,P->lf_VY->lf,&(P->D->Y_ncols),&(P->M_rank->nrows),&(P->M_VY),&(P->b_mlt),(addressable_0));
    GLOBAL_pthread_toc(); /* if bother */}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Yn,P->M_VY,P->M_rank,P->D->Y_ajdk,(addressable_0),&(P->lf_YnVY_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Wn,P->M_VY,P->M_rank,P->D->Y_ajdk,(addressable_0),&(P->lf_WnVY_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% An_Xn_ww : ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      printf(" %% nb %d: P->lf_YnVY_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_YnVY_[nb]->length,P->lf_YnVY_[nb]->row_stride,P->lf_YnVY_[nb]->col_stride,P->lf_YnVY_[nb]->lyr_stride);
      lfprintf(P->lf_YnVY_[nb]," %% P->lf_YnVY_[nb]->lf: ");
      printf(" %% nb %d: P->lf_WnVY_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_WnVY_[nb]->length,P->lf_WnVY_[nb]->row_stride,P->lf_WnVY_[nb]->col_stride,P->lf_WnVY_[nb]->lyr_stride);
      lfprintf(P->lf_WnVY_[nb]," %% P->lf_WnVY_[nb]->lf: ");
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->Y_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->A_bmr_j,E->A_bmr_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_YnVY_[nb],P->lf_YnVY_[nb]->lf,&(E->A_nrows),&(P->M_rank->nrows),&(P->M_YnVY_[nb]),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->Y_cbother){
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E->Z_bmr_j,E->Z_bmr_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_WnVY_[nb],P->lf_WnVY_[nb]->lf,&(E->Z_nrows),&(P->M_rank->nrows),&(P->M_WnVY_[nb]),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc trm_0: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother && D->Y_cbother){ 
      GLOBAL_pthread_tic();
      wrap_At_T_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_At,E->M_Tt,P->M_YnVY_[nb],P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_AtYnVY_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother && D->Y_cbother){
      GLOBAL_pthread_tic();
      wrap_At_T_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Zt,E->M_St,P->M_WnVY_[nb],P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_ZtWnVY_[nb]));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% At_T_Xn_ww : ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->A_rbother && D->A_cbother && D->Y_cbother){ 
	printf(" %% nb %d: P->lf_AtYnVY_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_AtYnVY_[nb]->length,P->lf_AtYnVY_[nb]->row_stride,P->lf_AtYnVY_[nb]->col_stride,P->lf_AtYnVY_[nb]->lyr_stride);
	lfprintf(P->lf_AtYnVY_[nb]," %% P->lf_AtYnVY_[nb]->lf: "); /* if bother */}
      if (E->Z_rbother && D->A_cbother && D->Y_cbother){ 
	printf(" %% nb %d: P->lf_ZtWnVY_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_ZtWnVY_[nb]->length,P->lf_ZtWnVY_[nb]->row_stride,P->lf_ZtWnVY_[nb]->col_stride,P->lf_ZtWnVY_[nb]->lyr_stride);
	lfprintf(P->lf_ZtWnVY_[nb]," %% P->lf_ZtWnVY_[nb]->lf: "); /* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  L_zero(P->lf_VYW);
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    ns_j = 0; ns_a = 0; ns_b = 0; 
    if (E->A_rbother && D->A_cbother && D->Y_cbother){ 
      dtmp = maximum(1,P->D->A_rpop_j_total);
      dra_plustimesequals(L2_get(P->lf_VYW,0,0,0,0,0,0),P->lf_VYW->length,L3_get(P->lf_AtYnVY_[nb],0,0,0,0,0,0,ns_j,ns_b,ns_a),+1.0/dtmp); /* if bother */}
    if (E->Z_rbother && D->A_cbother && D->Y_cbother){ 
      dtmp = maximum(1,P->D->Z_rpop_j_total);
      dra_plustimesequals(L2_get(P->lf_VYW,0,0,0,0,0,0),P->lf_VYW->length,L3_get(P->lf_ZtWnVY_[nb],0,0,0,0,0,0,ns_j,ns_b,ns_a),-1.0/dtmp); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
      if (E->A_rbother && D->A_cbother && D->Y_cbother){ 
	printf(" %% P->lf_VYW length %d row_stride %d col_stride %d lyr_stride %d\n",P->lf_VYW->length,P->lf_VYW->row_stride,P->lf_VYW->col_stride,P->lf_VYW->lyr_stride);
	lfprintf(P->lf_VYW," %% P->lf_VYW->lf: "); /* if bother */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/  
  if (verbose>1){ printf(" %% [finished pca_stage_c2]\n");}
}

void pca_stage_d(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  struct dcc_ajdk *D = P->D; 
  struct dcc_single **E_ = D->E_;
  int nbins = D->nbins; 
  int nb=0; struct dcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  int ns_j=0,ns_b=0,ns_a=0;
  double dtmp=0;
  if (verbose>1){ printf(" %% [entering pca_stage_d]\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  L_zero(P->lf_V);
  if (D->A_cbother){ 
    dtmp = maximum(1,P->D->A_cpop_j);
    dra_plustimesequals(L2_get(P->lf_V,0,0,0,0,0,0),P->lf_V->length,L2_get(P->lf_VAZ,0,0,0,0,0,0),+1.0/dtmp); /* if bother */}
  if (D->Y_cbother){ 
    dtmp = maximum(1,P->D->Y_cpop_j);
    dra_plustimesequals(L2_get(P->lf_V,0,0,0,0,0,0),P->lf_V->length,L2_get(P->lf_VYW,0,0,0,0,0,0),-1.0/dtmp); /* if bother */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
      if (D->A_cbother){ 
	printf(" %% P->lf_V length %d row_stride %d col_stride %d lyr_stride %d\n",P->lf_V->length,P->lf_V->row_stride,P->lf_V->col_stride,P->lf_V->lyr_stride);
	lfprintf(P->lf_V," %% P->lf_V->lf: "); /* if bother */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/  
  if (verbose>1){ printf(" %% [finished pca_stage_d]\n");}
}

void pca_stage_e(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  int nl=0,nb1=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering pca_stage_e]\n");}
  P->A_rpop_j_total[P->nx] = P->D->A_rpop_j_total;
  P->A_cpop_j[P->nx] = P->D->A_cpop_j;
  P->Irem[P->nx] = 0; for (nb1=0;nb1<P->D->nbins;nb1++){ P->Irem[P->nx] += (P->D->E_[nb1]->A_rpop_j>0?1:0);}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  for (nl=0;nl<P->lf_V->row_stride*P->lf_V->col_stride;nl++){
    P->V_[nl + P->nx*P->lf_V->row_stride*P->lf_V->col_stride] = P->lf_V->lf[nl];
    /* for (nl=0;nl<P->lf_V->row_stride*P->lf_V->col_stride;nl++){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
    sprintf(tmpchar," %%%% P->V_[%d]: ",P->nx);
    raprintf(&(P->V_[P->nx*P->lf_V->row_stride*P->lf_V->col_stride]),"double",P->lf_V->row_stride,P->lf_V->col_stride,tmpchar);
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/  
  if (verbose>1){ printf(" %% [finished pca_stage_e]\n");}
}

void pca_stage_a9(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  int nl=0,nb1=0;
  char tmpchar[FNAMESIZE];
  int itmp=0; double dtmp_row=0.0,dtmp_col=0.0;
  if (verbose>1){ printf(" %% [entering pca_stage_a9]\n");}
  P->A_rpop_j_total[P->nx] = P->D->A_rpop_j_total;
  P->A_cpop_j[P->nx] = P->D->A_cpop_j;
  P->Irem[P->nx] = 0; for (nb1=0;nb1<P->D->nbins;nb1++){ P->Irem[P->nx] += (P->D->E_[nb1]->A_rpop_j>0?1:0);}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  for (nl=0;nl<P->lf_V->row_stride*P->lf_V->col_stride;nl++){
    P->lf_V->lf[nl] = P->V_[nl + P->nx*P->lf_V->row_stride*P->lf_V->col_stride];
    /* for (nl=0;nl<P->lf_V->row_stride*P->lf_V->col_stride;nl++){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>4){
    sprintf(tmpchar," %%%% P->lf_V: ");
    lfprintf(P->lf_V," %% P->lf_V->lf: ");
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/  
  if (verbose>2){
    itmp=0;
    for (nl=0;nl<P->lf_V->row_stride;nl++){
      dtmp_row = P->lf_V->lf[nl + 0*P->lf_V->row_stride];
      dtmp_col = P->lf_V->lf[nl + 1*P->lf_V->row_stride];
      if ((itmp<1024) && ((dtmp_row!=0) || (dtmp_col!=0))){ 
	printf(" %% itmp %4d: P->lf_V->lf[%d] = %+16.8f , %+16.8f\n",itmp,dtmp_row,dtmp_col);
	itmp++; /* if ((itmp<1024) && ((dtmp_row!=0) || (dtmp_col!=0))){ } */}
      /* for (nl=0;nl<P->lf_V->length;nl++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/  
  if (verbose>1){ printf(" %% [finished pca_stage_a9]\n");}
}

void pca_stage_b9(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  int flag_ww_vs_uu=0;
  struct dcc_ajdk *D = P->D; 
  struct dcc_single **E_ = D->E_;
  int nbins = D->nbins; 
  int nb=0; struct dcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  int ns_j=0,ns_b=0,ns_a=0;
  int nl=0,itmp=0; double dtmp_row=0.0,dtmp_col=0.0;
  if (verbose>1){ printf(" %% [entering pca_stage_b9]\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  if (D->A_cbother){ 
    GLOBAL_pthread_tic(); 
    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),P->D->A_bmc_j,P->D->A_bmc_b,P->M_rank->mr_j,P->M_rank->mr_b,P->lf_V,P->lf_V->lf,&(P->D->A_ncols),&(P->M_rank->nrows),&(P->M_V),&(P->b_mlt),(addressable_0));
    GLOBAL_pthread_toc(); /* if bother */}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (flag_ww_vs_uu==1){
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->A_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_An_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_An,P->M_V,P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_AnV_[nb]));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E->Z_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_An_Xn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Zn,P->M_V,P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_ZnV_[nb]));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% An_Xn_ww : ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /* if (flag_ww_vs_uu==1){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (flag_ww_vs_uu==0){
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->A_rbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_An_Xn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_An,P->lf_V,P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_AnV_[nb]));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E->Z_rbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_An_Xn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Zn,P->lf_V,P->M_rank,P->D->A_ajdk,(addressable_0),&(P->lf_ZnV_[nb]));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% An_Xn_uu : ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /* if (flag_ww_vs_uu==0){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>4){
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      printf(" %% nb %d: P->lf_AnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_AnV_[nb]->length,P->lf_AnV_[nb]->row_stride,P->lf_AnV_[nb]->col_stride,P->lf_AnV_[nb]->lyr_stride);
      lfprintf(P->lf_AnV_[nb]," %% P->lf_AnV_[nb]->lf: ");
      printf(" %% nb %d: P->lf_ZnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_ZnV_[nb]->length,P->lf_ZnV_[nb]->row_stride,P->lf_ZnV_[nb]->col_stride,P->lf_ZnV_[nb]->lyr_stride);
      lfprintf(P->lf_ZnV_[nb]," %% P->lf_ZnV_[nb]->lf: ");
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      printf(" %% nb %d: P->lf_AnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_AnV_[nb]->length,P->lf_AnV_[nb]->row_stride,P->lf_AnV_[nb]->col_stride,P->lf_AnV_[nb]->lyr_stride);
      itmp=0;
      for (nl=0;nl<P->lf_AnV_[nb]->row_stride;nl++){
	dtmp_row = P->lf_AnV_[nb]->lf[nl + 0*P->lf_AnV_[nb]->row_stride];
	dtmp_col = P->lf_AnV_[nb]->lf[nl + 1*P->lf_AnV_[nb]->row_stride];
	if ((itmp<32) && ((dtmp_row!=0) || (dtmp_col!=0))){
	  printf(" %% itmp %2d: P->lf_AnV_[%d]->lf[%d] = %+16.8f , %+16.8f\n",itmp,nb,dtmp_row,dtmp_col);
	  itmp++; /* if ((itmp<32) && ((dtmp_row!=0) || (dtmp_col!=0))){ } */}
	/* for (nl=0;nl<P->lf_AnV_[nb]->row_stride;nl++){ } */}
      printf(" %% nb %d: P->lf_ZnV_[nb] length %d row_stride %d col_stride %d lyr_stride %d\n",nb,P->lf_ZnV_[nb]->length,P->lf_ZnV_[nb]->row_stride,P->lf_ZnV_[nb]->col_stride,P->lf_ZnV_[nb]->lyr_stride);
      itmp=0;
      for (nl=0;nl<P->lf_ZnV_[nb]->row_stride;nl++){
	dtmp_row = P->lf_ZnV_[nb]->lf[nl + 0*P->lf_ZnV_[nb]->row_stride];
	dtmp_col = P->lf_ZnV_[nb]->lf[nl + 1*P->lf_ZnV_[nb]->row_stride];
	if ((itmp<32) && ((dtmp_row!=0) || (dtmp_col!=0))){
	  printf(" %% itmp %2d: P->lf_ZnV_[%d]->lf[%d] = %+16.8f , %+16.8f\n",itmp,nb,dtmp_row,dtmp_col);
	  itmp++; /* if ((itmp<32) && ((dtmp_row!=0) || (dtmp_col!=0))){ } */}
	/* for (nl=0;nl<P->lf_ZnV_[nb]->row_stride;nl++){ } */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>1){ printf(" %% [finished pca_stage_b9]\n");}
}

void pca_stage_e9(struct P_handle *P)
{
  int verbose=GLOBAL_verbose;
  struct dcc_ajdk *D = P->D; 
  struct dcc_single **E_ = D->E_;
  int nbins = D->nbins; 
  int nb=0; struct dcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  int ns_j=0,ns_b=0,ns_a=0;
  int na_a=0,ma_a=0,ma_x=0,mz_a=0,mz_x=0,nb1=0;
  int tab_A=0,tab_Z=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering pca_stage_e9]\n");}
  P->A_rpop_j_total[P->nx] = P->D->A_rpop_j_total;
  P->A_cpop_j[P->nx] = P->D->A_cpop_j;
  P->Irem[P->nx] = 0; for (nb1=0;nb1<P->D->nbins;nb1++){ P->Irem[P->nx] += (P->D->E_[nb1]->A_rpop_j>0?1:0);}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  ma_x=0; mz_x=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    for (ma_a=0;ma_a<E->A_nrows;ma_a++){
      for (na_a=0;na_a<P->rank;na_a++){
	tab_A = (ma_x + ma_a) + na_a*D->A_nrows_total + P->nx*(D->A_nrows_total*P->rank);
	P->AnV_[tab_A] = *L2_get(P->lf_AnV_[nb],0,0,ma_a,0,0,na_a);
	/* for (na_a=0;na_a<P->rank;na_a++){ } */}
      /* for (ma_a=0;ma_a<E->A_nrows;ma_a++){ } */}
    ma_x += E->A_nrows;
    for (mz_a=0;mz_a<E->Z_nrows;mz_a++){
      for (na_a=0;na_a<P->rank;na_a++){
	tab_Z = (mz_x + mz_a) + na_a*D->Z_nrows_total + P->nx*(D->Z_nrows_total*P->rank);
	P->ZnV_[tab_Z] = *L2_get(P->lf_ZnV_[nb],0,0,mz_a,0,0,na_a);
	/* for (na_a=0;na_a<P->rank;na_a++){ } */}
      /* for (mz_a=0;mz_a<E->A_nrows;mz_a++){ } */}
    mz_x += E->Z_nrows;
    /* for (nb=0;nb<nbins;nb++){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      sprintf(tmpchar," %%%% P->AnV_[%d]: ",P->nx);
      raprintf(&(P->AnV_[P->nx*D->A_nrows_total*P->rank]),"double",D->A_nrows_total,P->rank,tmpchar);
      sprintf(tmpchar," %%%% P->ZnV_[%d]: ",P->nx);
      raprintf(&(P->ZnV_[P->nx*D->Z_nrows_total*P->rank]),"double",D->Z_nrows_total,P->rank,tmpchar);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/  
  if (verbose>1){ printf(" %% [finished pca_stage_e9]\n");}
}

void L2_AtZn_uu(struct L_handle *lf_At,struct L_handle *lf_Zn,struct L_handle *lf_ZnAt)
{
  /* calculates AnZt = An*Zt */
  int ma=0,mz=0,na=0;
  double dtmp = 0;
  if (lf_At->col_stride != lf_Zn->row_stride){ printf(" %% Warning! lf_At->col_stride != lf_Zn->col_stride: lf_At (%d,%d) lf_Zn (%d,%d) lf_ZnAt (%d,%d) in L2_AtZt_uu\n",lf_At->row_stride,lf_At->col_stride,lf_Zn->row_stride,lf_Zn->col_stride,lf_ZnAt->row_stride,lf_ZnAt->col_stride);}
  if (lf_ZnAt->row_stride != lf_At->row_stride){ printf(" %% Warning! lf_ZnAt->row_stride != lf_At->row_stride: lf_At (%d,%d) lf_Zn (%d,%d) lf_ZnAt (%d,%d) in L2_AtZt_uu\n",lf_At->row_stride,lf_At->col_stride,lf_Zn->row_stride,lf_Zn->col_stride,lf_ZnAt->row_stride,lf_ZnAt->col_stride);}
  if (lf_ZnAt->col_stride != lf_Zn->col_stride){ printf(" %% Warning! lf_ZnAt->col_stride != lf_Zn->row_stride: lf_At (%d,%d) lf_Zn (%d,%d) lf_ZnAt (%d,%d) in L2_AtZt_uu\n",lf_At->row_stride,lf_At->col_stride,lf_Zn->row_stride,lf_Zn->col_stride,lf_ZnAt->row_stride,lf_ZnAt->col_stride);}
  if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zn->row_stride*lf_Zn->col_stride>0){
    for (ma=0;ma<lf_At->row_stride;ma++){
      for (mz=0;mz<lf_Zn->col_stride;mz++){
	L2_set(lf_ZnAt,0,0,ma,0,0,mz,0);
	dtmp = 0;
	for (na=0;na<lf_At->col_stride;na++){
	  dtmp += (*L2_get(lf_At,0,0,ma,0,0,na))*(*L2_get(lf_Zn,0,0,na,0,0,mz));
	  /* for (na=0;na<lf_At->col_stride;na++){ } */}
	L2_set(lf_ZnAt,0,0,ma,0,0,mz,dtmp);
	/* for (mz=0;mz<lf_Zn->col_stride;mz++){ } */}
      /* for (ma=0;ma<lf_At->row_stride;ma++){ } */}
    /* if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zn->row_stride*lf_Zn->col_stride>0){ } */}
}

void L2_AtZt_uu(struct L_handle *lf_At,struct L_handle *lf_Zt,struct L_handle *lf_ZnAt)
{
  /* calculates AnZt = An*Zt */
  int ma=0,mz=0,na=0;
  double dtmp = 0;
  if (lf_At->col_stride != lf_Zt->col_stride){ printf(" %% Warning! lf_At->col_stride != lf_Zt->col_stride: lf_At (%d,%d) lf_Zt (%d,%d) lf_ZnAt (%d,%d) in L2_AtZt_uu\n",lf_At->row_stride,lf_At->col_stride,lf_Zt->row_stride,lf_Zt->col_stride,lf_ZnAt->row_stride,lf_ZnAt->col_stride);}
  if (lf_ZnAt->row_stride != lf_At->row_stride){ printf(" %% Warning! lf_ZnAt->row_stride != lf_At->row_stride: lf_At (%d,%d) lf_Zt (%d,%d) lf_ZnAt (%d,%d) in L2_AtZt_uu\n",lf_At->row_stride,lf_At->col_stride,lf_Zt->row_stride,lf_Zt->col_stride,lf_ZnAt->row_stride,lf_ZnAt->col_stride);}
  if (lf_ZnAt->col_stride != lf_Zt->row_stride){ printf(" %% Warning! lf_ZnAt->col_stride != lf_Zt->row_stride: lf_At (%d,%d) lf_Zt (%d,%d) lf_ZnAt (%d,%d) in L2_AtZt_uu\n",lf_At->row_stride,lf_At->col_stride,lf_Zt->row_stride,lf_Zt->col_stride,lf_ZnAt->row_stride,lf_ZnAt->col_stride);}
  if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zt->row_stride*lf_Zt->col_stride>0){
    for (ma=0;ma<lf_At->row_stride;ma++){
      for (mz=0;mz<lf_Zt->row_stride;mz++){
	L2_set(lf_ZnAt,0,0,ma,0,0,mz,0);
	dtmp = 0;
	for (na=0;na<lf_At->col_stride;na++){
	  dtmp += (*L2_get(lf_At,0,0,ma,0,0,na))*(*L2_get(lf_Zt,0,0,mz,0,0,na));
	  /* for (na=0;na<lf_At->col_stride;na++){ } */}
	L2_set(lf_ZnAt,0,0,ma,0,0,mz,dtmp);
	/* for (mz=0;mz<lf_Zt->row_stride;mz++){ } */}
      /* for (ma=0;ma<lf_At->row_stride;ma++){ } */}
    /* if (lf_At->row_stride*lf_At->col_stride>0 && lf_Zt->row_stride*lf_Zt->col_stride>0){ } */}
}

void L2_AnZt_mm(struct L_handle *lf_An,struct L_handle *lf_Zn,struct L_handle *lf_AnZt)
{
  /* calculates AnZt = An*Zt */
  if (lf_An->row_stride != lf_Zn->row_stride){ printf(" %% Warning! lf_An (%d,%d) lf_Zn (%d,%d) in L2_AnZt_mm\n",lf_An->row_stride,lf_An->col_stride,lf_Zn->row_stride,lf_Zn->col_stride);}
  if (lf_AnZt->row_stride != lf_An->col_stride){ printf(" %% Warning! lf_AnZt (%d,%d) lf_An (%d,%d) in L2_AnZt_mm\n",lf_AnZt->row_stride,lf_AnZt->col_stride,lf_An->row_stride,lf_An->col_stride);}
  if (lf_AnZt->col_stride != lf_Zn->col_stride){ printf(" %% Warning! lf_AnZt (%d,%d) lf_Zn (%d,%d) in L2_AnZt_mm\n",lf_AnZt->row_stride,lf_AnZt->col_stride,lf_Zn->row_stride,lf_Zn->col_stride);}
  if (lf_An->row_stride*lf_An->col_stride>0 && lf_Zn->row_stride*lf_Zn->col_stride>0){
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,lf_An->col_stride,lf_Zn->col_stride,lf_An->row_stride,1.0,lf_An->lf,lf_An->row_stride,lf_Zn->lf,lf_Zn->row_stride,0.0,lf_AnZt->lf,lf_AnZt->row_stride);
    /* if (lf_An->row_stride*lf_An->col_stride>0 && lf_Zn->row_stride*lf_Zn->col_stride>0){ } */}
}

void pca_uu(struct P_handle *P)
{
  int verbose=GLOBAL_verbose; int iteration_max = GLOBAL_TEST_niter;
  struct dcc_ajdk *D=P->D;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL;
  int na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0; int ma_x=0;
  int ny_j=0,ny_b=0,ny_a=0,mz_j=0,mz_b=0,mz_a=0; int mz_x=0;
  unsigned char *A_tag=NULL,*Z_tag=NULL,*Y_tag=NULL,*W_tag=NULL;
  struct M_handle *M_An=NULL,*M_Zn=NULL,*M_Yn=NULL,*M_Wn=NULL;
  double *A_ajdk=NULL,*Y_ajdk=NULL;
  int A_pcols=0,Y_pcols=0;
  double dtmp=0;
  int nl=0;
  if (verbose>1){ printf(" %% [entering pca_uu]\n");}
  if (P->lf_At==NULL){ 
    P->lf_At = L_handle_make((unsigned long long int)D->A_nrows_total*(unsigned long long int)D->A_ncols); 
    P->lf_At->spacing_row = SPACING_a; P->lf_At->row_stride = D->A_nrows_total;
    P->lf_At->spacing_col = SPACING_a; P->lf_At->col_stride = D->A_ncols;
    /* if (P->lf_At==NULL){ } */}
  if (P->lf_Zt==NULL){ 
    P->lf_Zt = L_handle_make((unsigned long long int)D->Z_nrows_total*(unsigned long long int)D->A_ncols); 
    P->lf_Zt->spacing_row = SPACING_a; P->lf_Zt->row_stride = D->Z_nrows_total;
    P->lf_Zt->spacing_col = SPACING_a; P->lf_Zt->col_stride = D->A_ncols;
    /* if (P->lf_Zt==NULL){ } */}
  if (P->lf_Yt==NULL){ 
    P->lf_Yt = L_handle_make((unsigned long long int)D->A_nrows_total*(unsigned long long int)D->Y_ncols); 
    P->lf_Yt->spacing_row = SPACING_a; P->lf_Yt->row_stride = D->A_nrows_total;
    P->lf_Yt->spacing_col = SPACING_a; P->lf_Yt->col_stride = D->Y_ncols;
    /* if (P->lf_Yt==NULL){ } */}
  if (P->lf_Wt==NULL){ 
    P->lf_Wt = L_handle_make((unsigned long long int)D->Z_nrows_total*(unsigned long long int)D->Y_ncols); 
    P->lf_Wt->spacing_row = SPACING_a; P->lf_Wt->row_stride = D->Z_nrows_total;
    P->lf_Wt->spacing_col = SPACING_a; P->lf_Wt->col_stride = D->Y_ncols;
    /* if (P->lf_Wt==NULL){ } */}
  L_zero(P->lf_At); L_zero(P->lf_Zt); L_zero(P->lf_Yt); L_zero(P->lf_Wt);
  if (P->lf_An==NULL){ 
    P->lf_An = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_nrows_total); 
    P->lf_An->spacing_row = SPACING_a; P->lf_An->row_stride = D->A_ncols;
    P->lf_An->spacing_col = SPACING_a; P->lf_An->col_stride = D->A_nrows_total;
    /* if (P->lf_An==NULL){ } */}
  if (P->lf_Zn==NULL){ 
    P->lf_Zn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->Z_nrows_total); 
    P->lf_Zn->spacing_row = SPACING_a; P->lf_Zn->row_stride = D->A_ncols;
    P->lf_Zn->spacing_col = SPACING_a; P->lf_Zn->col_stride = D->Z_nrows_total;
    /* if (P->lf_Zn==NULL){ } */}
  if (P->lf_Yn==NULL){ 
    P->lf_Yn = L_handle_make((unsigned long long int)D->Y_ncols*(unsigned long long int)D->A_nrows_total); 
    P->lf_Yn->spacing_row = SPACING_a; P->lf_Yn->row_stride = D->Y_ncols;
    P->lf_Yn->spacing_col = SPACING_a; P->lf_Yn->col_stride = D->A_nrows_total;
    /* if (P->lf_Yn==NULL){ } */}
  if (P->lf_Wn==NULL){ 
    P->lf_Wn = L_handle_make((unsigned long long int)D->Y_ncols*(unsigned long long int)D->Z_nrows_total); 
    P->lf_Wn->spacing_row = SPACING_a; P->lf_Wn->row_stride = D->Y_ncols;
    P->lf_Wn->spacing_col = SPACING_a; P->lf_Wn->col_stride = D->Z_nrows_total;
    /* if (P->lf_Wn==NULL){ } */}
  L_zero(P->lf_An); L_zero(P->lf_Zn); L_zero(P->lf_Yt); L_zero(P->lf_Wn);
  /* %%%%%%%%%%%%%%%% */
  if (P->lf_AtAn==NULL){
    P->lf_AtAn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_AtAn->spacing_row = SPACING_a; P->lf_AtAn->row_stride = D->A_ncols;
    P->lf_AtAn->spacing_col = SPACING_a; P->lf_AtAn->col_stride = D->A_ncols;
    /* if (P->lf_AtAn==NULL){ } */}
  if (P->lf_ZtZn==NULL){
    P->lf_ZtZn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_ZtZn->spacing_row = SPACING_a; P->lf_ZtZn->row_stride = D->A_ncols;
    P->lf_ZtZn->spacing_col = SPACING_a; P->lf_ZtZn->col_stride = D->A_ncols;
    /* if (P->lf_AtAn==NULL){ } */}
  if (P->lf_AtYn==NULL){
    P->lf_AtYn = L_handle_make((unsigned long long int)D->Y_ncols*(unsigned long long int)D->A_ncols);
    P->lf_AtYn->spacing_row = SPACING_a; P->lf_AtYn->row_stride = D->Y_ncols;
    P->lf_AtYn->spacing_col = SPACING_a; P->lf_AtYn->col_stride = D->A_ncols;
    /* if (P->lf_AtYn==NULL){ } */}
  if (P->lf_ZtWn==NULL){
    P->lf_ZtWn = L_handle_make((unsigned long long int)D->Y_ncols*(unsigned long long int)D->A_ncols);
    P->lf_ZtWn->spacing_row = SPACING_a; P->lf_ZtWn->row_stride = D->Y_ncols;
    P->lf_ZtWn->spacing_col = SPACING_a; P->lf_ZtWn->col_stride = D->A_ncols;
    /* if (P->lf_ZtWn==NULL){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (P->lf_AtAnAtAn==NULL){
    P->lf_AtAnAtAn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_AtAnAtAn->spacing_row = SPACING_a; P->lf_AtAnAtAn->row_stride = D->A_ncols;
    P->lf_AtAnAtAn->spacing_col = SPACING_a; P->lf_AtAnAtAn->col_stride = D->A_ncols;
    /* if (P->lf_AtAnAtAn==NULL){ } */}
  if (P->lf_AtAnZtZn==NULL){
    P->lf_AtAnZtZn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_AtAnZtZn->spacing_row = SPACING_a; P->lf_AtAnZtZn->row_stride = D->A_ncols;
    P->lf_AtAnZtZn->spacing_col = SPACING_a; P->lf_AtAnZtZn->col_stride = D->A_ncols;
    /* if (P->lf_AtAnZtZn==NULL){ } */}
  if (P->lf_ZtZnAtAn==NULL){
    P->lf_ZtZnAtAn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_ZtZnAtAn->spacing_row = SPACING_a; P->lf_ZtZnAtAn->row_stride = D->A_ncols;
    P->lf_ZtZnAtAn->spacing_col = SPACING_a; P->lf_ZtZnAtAn->col_stride = D->A_ncols;
    /* if (P->lf_ZtZnAtAn==NULL){ } */}
  if (P->lf_ZtZnZtZn==NULL){
    P->lf_ZtZnZtZn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_ZtZnZtZn->spacing_row = SPACING_a; P->lf_ZtZnZtZn->row_stride = D->A_ncols;
    P->lf_ZtZnZtZn->spacing_col = SPACING_a; P->lf_ZtZnZtZn->col_stride = D->A_ncols;
    /* if (P->lf_ZtZnZtZn==NULL){ } */}
  if (P->lf_AtYnYtAn==NULL){
    P->lf_AtYnYtAn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_AtYnYtAn->spacing_row = SPACING_a; P->lf_AtYnYtAn->row_stride = D->A_ncols;
    P->lf_AtYnYtAn->spacing_col = SPACING_a; P->lf_AtYnYtAn->col_stride = D->A_ncols;
    /* if (P->lf_AtYnYtAn==NULL){ } */}
  if (P->lf_AtYnWtZn==NULL){
    P->lf_AtYnWtZn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_AtYnWtZn->spacing_row = SPACING_a; P->lf_AtYnWtZn->row_stride = D->A_ncols;
    P->lf_AtYnWtZn->spacing_col = SPACING_a; P->lf_AtYnWtZn->col_stride = D->A_ncols;
    /* if (P->lf_AtYnWtZn==NULL){ } */}
  if (P->lf_ZtWnYtAn==NULL){
    P->lf_ZtWnYtAn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_ZtWnYtAn->spacing_row = SPACING_a; P->lf_ZtWnYtAn->row_stride = D->A_ncols;
    P->lf_ZtWnYtAn->spacing_col = SPACING_a; P->lf_ZtWnYtAn->col_stride = D->A_ncols;
    /* if (P->lf_ZtWnYtAn==NULL){ } */}
  if (P->lf_ZtWnWtZn==NULL){
    P->lf_ZtWnWtZn = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_ZtWnWtZn->spacing_row = SPACING_a; P->lf_ZtWnWtZn->row_stride = D->A_ncols;
    P->lf_ZtWnWtZn->spacing_col = SPACING_a; P->lf_ZtWnWtZn->col_stride = D->A_ncols;
    /* if (P->lf_ZtWnWtZn==NULL){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (P->lf_S==NULL){
    P->lf_S = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)D->A_ncols);
    P->lf_S->spacing_row = SPACING_a; P->lf_S->row_stride = D->A_ncols;
    P->lf_S->spacing_col = SPACING_a; P->lf_S->col_stride = D->A_ncols;
    /* if (P->lf_S==NULL){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (P->lf_Vt==NULL){
    P->lf_Vt = L_handle_make((unsigned long long int)D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_Vt->spacing_row = SPACING_a; P->lf_Vt->row_stride = D->A_ncols;
    P->lf_Vt->spacing_col = SPACING_a; P->lf_Vt->col_stride = P->rank;
    /* if (P->lf_Vt==NULL){ } */}
  if (P->lf_R2==NULL){
    P->lf_R2 = L_handle_make((unsigned long long int)P->rank*(unsigned long long int)P->rank); 
    P->lf_R2->spacing_row = SPACING_a; P->lf_R2->row_stride = P->rank;
    P->lf_R2->spacing_col = SPACING_a; P->lf_R2->col_stride = P->rank;
    /* if (P->lf_R2==NULL){ } */}
  /* %%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% extracting lf_At\n");}
  ma_x=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; M_An = E->M_An; A_ajdk = D->A_ajdk; A_pcols = psize(M_An->ncols);
    ma_j=0; 
    while (ma_j<M_An->rpop_j){
      ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
      A_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length]));
      dtmp=0;
      na_j=0;
      while (na_j<M_An->cpop_j){
	na_a = M_An->n_a_[na_j]; na_b = M_An->n_b_[na_j];
	dtmp = (bget____(A_tag,na_a) - (A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*sqrt(A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols] : 1);
	L2_set(P->lf_At,0,0,ma_x+ma_a,0,0,na_a,dtmp);
	na_j++; /* while (na_j<M_An->cpop_j){ } */}
      ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
    ma_x += M_An->nrows;
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  /* %%%%%%%%%%%%%%%% */
  L2_transpose(P->lf_An,P->lf_At);
  /* %%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% extracting lf_Zt\n");}
  mz_x=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; M_Zn = E->M_Zn; A_ajdk = D->A_ajdk; A_pcols = psize(M_Zn->ncols);
    mz_j=0; 
    while (mz_j<M_Zn->rpop_j){
      mz_a = M_Zn->m_a_[mz_j]; mz_b = M_Zn->m_b_[mz_j];
      Z_tag = (unsigned char *)(&(M_Zn->wX[mz_b*M_Zn->mc_length]));
      dtmp=0;
      na_j=0;
      while (na_j<M_Zn->cpop_j){
	na_a = M_Zn->n_a_[na_j]; na_b = M_Zn->n_b_[na_j];
	dtmp = (bget____(Z_tag,na_a) - (A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*sqrt(A_ajdk ? A_ajdk[na_a/POPLENGTH + AJDK_0_1*A_pcols] : 1);
	L2_set(P->lf_Zt,0,0,mz_x+mz_a,0,0,na_a,dtmp);
	na_j++; /* while (na_j<M_Zn->cpop_j){ } */}
      mz_j++; /* while (mz_j<M_Zn->rpop_j){ } */}
    mz_x += M_Zn->nrows;
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  /* %%%%%%%%%%%%%%%% */
  L2_transpose(P->lf_Zn,P->lf_Zt);
  /* %%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% extracting lf_Yt\n");}
  ma_x=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; M_Yn = E->M_Yn; Y_ajdk = D->Y_ajdk; Y_pcols = psize(M_Yn->ncols);
    ma_j=0; 
    while (ma_j<M_Yn->rpop_j){
      ma_a = M_Yn->m_a_[ma_j]; ma_b = M_Yn->m_b_[ma_j];
      Y_tag = (unsigned char *)(&(M_Yn->wX[ma_b*M_Yn->mc_length]));
      dtmp=0;
      ny_j=0;
      while (ny_j<M_Yn->cpop_j){
	ny_a = M_Yn->n_a_[ny_j]; ny_b = M_Yn->n_b_[ny_j];
	dtmp = (bget____(Y_tag,ny_a) - (Y_ajdk ? Y_ajdk[ny_a/POPLENGTH + AJDK_1_0*Y_pcols] : 0))*sqrt(Y_ajdk ? Y_ajdk[ny_a/POPLENGTH + AJDK_0_1*Y_pcols] : 1);
	L2_set(P->lf_Yt,0,0,ma_x+ma_a,0,0,ny_a,dtmp);
	ny_j++; /* while (ny_j<M_Yn->cpop_j){ } */}
      ma_j++; /* while (ma_j<M_Yn->rpop_j){ } */}
    ma_x += M_Yn->nrows;
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  /* %%%%%%%%%%%%%%%% */
  L2_transpose(P->lf_Yn,P->lf_Yt);
  /* %%%%%%%%%%%%%%%% */
  if (verbose>1){ printf(" %% extracting lf_Wt\n");}
  mz_x=0;
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; M_Wn = E->M_Wn; Y_ajdk = D->Y_ajdk; Y_pcols = psize(M_Wn->ncols);
    mz_j=0; 
    while (mz_j<M_Wn->rpop_j){
      mz_a = M_Wn->m_a_[mz_j]; mz_b = M_Wn->m_b_[mz_j];
      W_tag = (unsigned char *)(&(M_Wn->wX[mz_b*M_Wn->mc_length]));
      dtmp=0;
      ny_j=0;
      while (ny_j<M_Wn->cpop_j){
	ny_a = M_Wn->n_a_[ny_j]; ny_b = M_Wn->n_b_[ny_j];
	dtmp = (bget____(W_tag,ny_a) - (Y_ajdk ? Y_ajdk[ny_a/POPLENGTH + AJDK_1_0*Y_pcols] : 0))*sqrt(Y_ajdk ? Y_ajdk[ny_a/POPLENGTH + AJDK_0_1*Y_pcols] : 1);
	L2_set(P->lf_Wt,0,0,mz_x+mz_a,0,0,ny_a,dtmp);
	ny_j++; /* while (ny_j<M_Wn->cpop_j){ } */}
      mz_j++; /* while (mz_j<M_Wn->rpop_j){ } */}
    mz_x += M_Wn->nrows;
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  /* %%%%%%%%%%%%%%%% */
  L2_transpose(P->lf_Wn,P->lf_Wt);
  /* %%%%%%%%%%%%%%%% */
  if (verbose>0){
    printf(" %% P->lf_At (%d,%d)\n",P->lf_At->row_stride,P->lf_At->col_stride);
    printf(" %% P->lf_An (%d,%d)\n",P->lf_An->row_stride,P->lf_An->col_stride);
    printf(" %% P->lf_Zt (%d,%d)\n",P->lf_Zt->row_stride,P->lf_Zt->col_stride);
    printf(" %% P->lf_Zn (%d,%d)\n",P->lf_Zn->row_stride,P->lf_Zn->col_stride);
    printf(" %% P->lf_Yt (%d,%d)\n",P->lf_Yt->row_stride,P->lf_Yt->col_stride);
    printf(" %% P->lf_Yn (%d,%d)\n",P->lf_Yn->row_stride,P->lf_Yn->col_stride);
    printf(" %% P->lf_Wt (%d,%d)\n",P->lf_Wt->row_stride,P->lf_Wt->col_stride);
    printf(" %% P->lf_Wn (%d,%d)\n",P->lf_Wn->row_stride,P->lf_Wn->col_stride);
    /* if (verbose>0){ } */}
  if (verbose>1){
    lfprintf(P->lf_At," P->lf_At: ");
    lfprintf(P->lf_An," P->lf_An: ");
    lfprintf(P->lf_Zt," P->lf_Zt: ");
    lfprintf(P->lf_Zn," P->lf_Zn: ");
    lfprintf(P->lf_Yt," P->lf_Yt: ");
    lfprintf(P->lf_Yn," P->lf_Yn: ");
    lfprintf(P->lf_Wt," P->lf_Wt: ");
    lfprintf(P->lf_Wn," P->lf_Wn: ");
    /* if (verbose>1){ } */}
  /* %%%%%%%%%%%%%%%% */
  /* Construct S */
  /* %%%%%%%%%%%%%%%% */
  L2_AnZt_mm(P->lf_At,P->lf_At,P->lf_AtAn);
  L2_AnZt_mm(P->lf_Zt,P->lf_Zt,P->lf_ZtZn);
  L2_AnZt_mm(P->lf_Yt,P->lf_At,P->lf_AtYn);
  L2_AnZt_mm(P->lf_Wt,P->lf_Zt,P->lf_ZtWn);
  L2_AnZt_mm(P->lf_AtAn,P->lf_AtAn,P->lf_AtAnAtAn);
  L2_AnZt_mm(P->lf_AtAn,P->lf_ZtZn,P->lf_AtAnZtZn);
  L2_AnZt_mm(P->lf_ZtZn,P->lf_AtAn,P->lf_ZtZnAtAn);
  L2_AnZt_mm(P->lf_ZtZn,P->lf_ZtZn,P->lf_ZtZnZtZn);
  L2_AnZt_mm(P->lf_AtYn,P->lf_AtYn,P->lf_AtYnYtAn);
  L2_AnZt_mm(P->lf_AtYn,P->lf_ZtWn,P->lf_AtYnWtZn);
  L2_AnZt_mm(P->lf_ZtWn,P->lf_AtYn,P->lf_ZtWnYtAn);
  L2_AnZt_mm(P->lf_ZtWn,P->lf_ZtWn,P->lf_ZtWnWtZn);
  L_zero(P->lf_S);
  dtmp = maximum(1,D->A_rpop_j_total*D->A_rpop_j_total*D->A_cpop_j); dra_plustimesequals(P->lf_S->lf,P->lf_S->length,P->lf_AtAnAtAn->lf,+1.0/dtmp);
  dtmp = maximum(1,D->A_rpop_j_total*D->Z_rpop_j_total*D->A_cpop_j); dra_plustimesequals(P->lf_S->lf,P->lf_S->length,P->lf_AtAnZtZn->lf,-1.0/dtmp);
  dtmp = maximum(1,D->Z_rpop_j_total*D->A_rpop_j_total*D->A_cpop_j); dra_plustimesequals(P->lf_S->lf,P->lf_S->length,P->lf_ZtZnAtAn->lf,-1.0/dtmp);
  dtmp = maximum(1,D->Z_rpop_j_total*D->Z_rpop_j_total*D->A_cpop_j); dra_plustimesequals(P->lf_S->lf,P->lf_S->length,P->lf_ZtZnZtZn->lf,+1.0/dtmp);
  dtmp = maximum(1,D->A_rpop_j_total*D->A_rpop_j_total*D->Y_cpop_j); dra_plustimesequals(P->lf_S->lf,P->lf_S->length,P->lf_AtYnYtAn->lf,-1.0/dtmp);
  dtmp = maximum(1,D->A_rpop_j_total*D->Z_rpop_j_total*D->Y_cpop_j); dra_plustimesequals(P->lf_S->lf,P->lf_S->length,P->lf_AtYnWtZn->lf,+1.0/dtmp);
  dtmp = maximum(1,D->Z_rpop_j_total*D->A_rpop_j_total*D->Y_cpop_j); dra_plustimesequals(P->lf_S->lf,P->lf_S->length,P->lf_ZtWnYtAn->lf,+1.0/dtmp);
  dtmp = maximum(1,D->Z_rpop_j_total*D->Z_rpop_j_total*D->Y_cpop_j); dra_plustimesequals(P->lf_S->lf,P->lf_S->length,P->lf_ZtWnWtZn->lf,-1.0/dtmp);
  if (verbose>1){
    lfprintf(P->lf_S," P->lf_S: ");
    /* if (verbose>1){ } */}
  /* %%%%%%%%%%%%%%%% */
  /* Single iteration using S */
  /* %%%%%%%%%%%%%%%% */
  if (iteration_max<=1){
  pca_stage_a0b(P);
  pca_stage_a1b(P);
  //lfprintf(P->lf_Vt," P->lf_Vt: ");
  L2_transpose(P->lf_Vn,P->lf_Vt);
  L2_AtZt_uu(P->lf_S,P->lf_Vn,P->lf_Vt);
  //lfprintf(P->lf_Vt," %% P->lf_Vt->lf: ");
  printf(" %% Single iteration using S: error %0.16f\n",dra_diff(P->lf_Vt->lf,P->lf_V->lf,P->lf_V->length,1));
  /* if (iteration_max<=1){ } */}
  /* %%%%%%%%%%%%%%%% */
  /* Single iteration using D */
  /* %%%%%%%%%%%%%%%% */
  if (iteration_max<=1){
  pca_stage_a0b(P);
  pca_stage_a1b(P);
  //lfprintf(P->lf_Vt," P->lf_Vt: ");
  L2_transpose(P->lf_Vn,P->lf_Vt);
  //lfprintf(P->lf_Vn," P->lf_Vn: ");
  if (P->lf_VtAt==NULL){
    P->lf_VtAt = L_handle_make((unsigned long long int)P->D->A_nrows_total*(unsigned long long int)P->rank); 
    P->lf_VtAt->spacing_row = SPACING_a; P->lf_VtAt->row_stride = P->D->A_nrows_total;
    P->lf_VtAt->spacing_col = SPACING_a; P->lf_VtAt->col_stride = P->rank;
    /* if (P->lf_VtAt==NULL){ } */}
  L2_AtZt_uu(P->lf_At,P->lf_Vn,P->lf_VtAt);
  //lfprintf(P->lf_VtAt," %% P->lf_VtAt->lf: ");
  if (P->lf_VtZt==NULL){
    P->lf_VtZt = L_handle_make((unsigned long long int)P->D->Z_nrows_total*(unsigned long long int)P->rank); 
    P->lf_VtZt->spacing_row = SPACING_a; P->lf_VtZt->row_stride = P->D->Z_nrows_total;
    P->lf_VtZt->spacing_col = SPACING_a; P->lf_VtZt->col_stride = P->rank;
    /* if (P->lf_VtZt==NULL){ } */}
  L2_AtZt_uu(P->lf_Zt,P->lf_Vn,P->lf_VtZt);
  //lfprintf(P->lf_VtZt," %% P->lf_VtZt->lf: ");
  if (P->lf_VtZtZn==NULL){
    P->lf_VtZtZn = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_VtZtZn->spacing_row = SPACING_a; P->lf_VtZtZn->row_stride = P->D->A_ncols;
    P->lf_VtZtZn->spacing_col = SPACING_a; P->lf_VtZtZn->col_stride = P->rank;
    /* if (P->lf_VtZtZn==NULL){ } */}
  if (P->lf_VtAtAn==NULL){
    P->lf_VtAtAn = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_VtAtAn->spacing_row = SPACING_a; P->lf_VtAtAn->row_stride = P->D->A_ncols;
    P->lf_VtAtAn->spacing_col = SPACING_a; P->lf_VtAtAn->col_stride = P->rank;
    /* if (P->lf_VtAtAn==NULL){ } */}
  if (P->lf_VtZtWn==NULL){
    P->lf_VtZtWn = L_handle_make((unsigned long long int)P->D->Y_ncols*(unsigned long long int)P->rank); 
    P->lf_VtZtWn->spacing_row = SPACING_a; P->lf_VtZtWn->row_stride = P->D->Y_ncols;
    P->lf_VtZtWn->spacing_col = SPACING_a; P->lf_VtZtWn->col_stride = P->rank;
    /* if (P->lf_VtZtWn==NULL){ } */}
  if (P->lf_VtAtYn==NULL){
    P->lf_VtAtYn = L_handle_make((unsigned long long int)P->D->Y_ncols*(unsigned long long int)P->rank); 
    P->lf_VtAtYn->spacing_row = SPACING_a; P->lf_VtAtYn->row_stride = P->D->Y_ncols;
    P->lf_VtAtYn->spacing_col = SPACING_a; P->lf_VtAtYn->col_stride = P->rank;
    /* if (P->lf_VtAtYn==NULL){ } */}
  L2_AtZn_uu(P->lf_An,P->lf_VtAt,P->lf_VtAtAn);
  L2_AtZn_uu(P->lf_Yn,P->lf_VtAt,P->lf_VtAtYn);
  L2_AtZn_uu(P->lf_Zn,P->lf_VtZt,P->lf_VtZtZn);
  L2_AtZn_uu(P->lf_Wn,P->lf_VtZt,P->lf_VtZtWn);
  //lfprintf(P->lf_VtAtAn," %% P->lf_VtAtAn->lf: ");
  //lfprintf(P->lf_VtAtYn," %% P->lf_VtAtYn->lf: ");
  //lfprintf(P->lf_VtZtZn," %% P->lf_VtZtZn->lf: ");
  //lfprintf(P->lf_VtZtWn," %% P->lf_VtZtWn->lf: ");
  if (P->lf_VAt==NULL){
    P->lf_VAt = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_VAt->spacing_row = SPACING_a; P->lf_VAt->row_stride = P->D->A_ncols;
    P->lf_VAt->spacing_col = SPACING_a; P->lf_VAt->col_stride = P->rank;
    /* if (P->lf_VAt==NULL){ } */}
  if (P->lf_VYt==NULL){
    P->lf_VYt = L_handle_make((unsigned long long int)P->D->Y_ncols*(unsigned long long int)P->rank); 
    P->lf_VYt->spacing_row = SPACING_a; P->lf_VYt->row_stride = P->D->Y_ncols;
    P->lf_VYt->spacing_col = SPACING_a; P->lf_VYt->col_stride = P->rank;
    /* if (P->lf_VYt==NULL){ } */}
  L_zero(P->lf_VAt);
  dtmp = P->D->A_nrows_total; if (dtmp>0){ dra_plustimesequals(P->lf_VAt->lf,P->lf_VAt->length,P->lf_VtAtAn->lf,+1.0/dtmp);}
  dtmp = P->D->Z_nrows_total; if (dtmp>0){ dra_plustimesequals(P->lf_VAt->lf,P->lf_VAt->length,P->lf_VtZtZn->lf,-1.0/dtmp);}
  //lfprintf(P->lf_VAt," %% P->lf_VAt->lf: ");
  L_zero(P->lf_VYt);
  dtmp = P->D->A_nrows_total; if (dtmp>0){ dra_plustimesequals(P->lf_VYt->lf,P->lf_VYt->length,P->lf_VtAtYn->lf,+1.0/dtmp);}
  dtmp = P->D->Z_nrows_total; if (dtmp>0){ dra_plustimesequals(P->lf_VYt->lf,P->lf_VYt->length,P->lf_VtZtWn->lf,-1.0/dtmp);}
  //lfprintf(P->lf_VYt," %% P->lf_VYt->lf: ");
  if (P->lf_VAtAt==NULL){
    P->lf_VAtAt = L_handle_make((unsigned long long int)P->D->A_nrows_total*(unsigned long long int)P->rank); 
    P->lf_VAtAt->spacing_row = SPACING_a; P->lf_VAtAt->row_stride = P->D->A_nrows_total;
    P->lf_VAtAt->spacing_col = SPACING_a; P->lf_VAtAt->col_stride = P->rank;
    /* if (P->lf_VAtAt==NULL){ } */}
  L2_AtZn_uu(P->lf_At,P->lf_VAt,P->lf_VAtAt);
  //lfprintf(P->lf_VAtAt," %% P->lf_VAtAt->lf: ");
  if (P->lf_VAtZt==NULL){
    P->lf_VAtZt = L_handle_make((unsigned long long int)P->D->Z_nrows_total*(unsigned long long int)P->rank); 
    P->lf_VAtZt->spacing_row = SPACING_a; P->lf_VAtZt->row_stride = P->D->Z_nrows_total;
    P->lf_VAtZt->spacing_col = SPACING_a; P->lf_VAtZt->col_stride = P->rank;
    /* if (P->lf_VAtZt==NULL){ } */}
  L2_AtZn_uu(P->lf_Zt,P->lf_VAt,P->lf_VAtZt);
  //lfprintf(P->lf_VAtZt," %% P->lf_VAtZt->lf: ");
  if (P->lf_VAtZtZn==NULL){
    P->lf_VAtZtZn = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_VAtZtZn->spacing_row = SPACING_a; P->lf_VAtZtZn->row_stride = P->D->A_ncols;
    P->lf_VAtZtZn->spacing_col = SPACING_a; P->lf_VAtZtZn->col_stride = P->rank;
    /* if (P->lf_VAtZtZn==NULL){ } */}
  if (P->lf_VAtAtAn==NULL){
    P->lf_VAtAtAn = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_VAtAtAn->spacing_row = SPACING_a; P->lf_VAtAtAn->row_stride = P->D->A_ncols;
    P->lf_VAtAtAn->spacing_col = SPACING_a; P->lf_VAtAtAn->col_stride = P->rank;
    /* if (P->lf_VAtAtAn==NULL){ } */}
  L2_AtZn_uu(P->lf_An,P->lf_VAtAt,P->lf_VAtAtAn);
  L2_AtZn_uu(P->lf_Zn,P->lf_VAtZt,P->lf_VAtZtZn);
  //lfprintf(P->lf_VAtAtAn," %% P->lf_VAtAtAn->lf: ");
  //lfprintf(P->lf_VAtZtZn," %% P->lf_VAtZtZn->lf: ");
  if (P->lf_VAZt==NULL){
    P->lf_VAZt = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_VAZt->spacing_row = SPACING_a; P->lf_VAZt->row_stride = P->D->A_ncols;
    P->lf_VAZt->spacing_col = SPACING_a; P->lf_VAZt->col_stride = P->rank;
    /* if (P->lf_VAZt==NULL){ } */}
  L_zero(P->lf_VAZt);
  dtmp = P->D->A_nrows_total; if (dtmp>0){ dra_plustimesequals(P->lf_VAZt->lf,P->lf_VAZt->length,P->lf_VAtAtAn->lf,+1.0/dtmp);}
  dtmp = P->D->Z_nrows_total; if (dtmp>0){ dra_plustimesequals(P->lf_VAZt->lf,P->lf_VAZt->length,P->lf_VAtZtZn->lf,-1.0/dtmp);}
  //lfprintf(P->lf_VAZt," %% P->lf_VAZt->lf: ");
  if (P->lf_VYtYt==NULL){
    P->lf_VYtYt = L_handle_make((unsigned long long int)P->D->A_nrows_total*(unsigned long long int)P->rank); 
    P->lf_VYtYt->spacing_row = SPACING_a; P->lf_VYtYt->row_stride = P->D->A_nrows_total;
    P->lf_VYtYt->spacing_col = SPACING_a; P->lf_VYtYt->col_stride = P->rank;
    /* if (P->lf_VYtYt==NULL){ } */}
  L2_AtZn_uu(P->lf_Yt,P->lf_VYt,P->lf_VYtYt);
  //lfprintf(P->lf_VYtYt," %% P->lf_VYtYt->lf: ");
  if (P->lf_VYtWt==NULL){
    P->lf_VYtWt = L_handle_make((unsigned long long int)P->D->Z_nrows_total*(unsigned long long int)P->rank); 
    P->lf_VYtWt->spacing_row = SPACING_a; P->lf_VYtWt->row_stride = P->D->Z_nrows_total;
    P->lf_VYtWt->spacing_col = SPACING_a; P->lf_VYtWt->col_stride = P->rank;
    /* if (P->lf_VYtWt==NULL){ } */}
  L2_AtZn_uu(P->lf_Wt,P->lf_VYt,P->lf_VYtWt);
  //lfprintf(P->lf_VYtWt," %% P->lf_VYtWt->lf: ");
  if (P->lf_VYtWtZn==NULL){
    P->lf_VYtWtZn = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_VYtWtZn->spacing_row = SPACING_a; P->lf_VYtWtZn->row_stride = P->D->A_ncols;
    P->lf_VYtWtZn->spacing_col = SPACING_a; P->lf_VYtWtZn->col_stride = P->rank;
    /* if (P->lf_VYtWtZn==NULL){ } */}
  if (P->lf_VYtYtAn==NULL){
    P->lf_VYtYtAn = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_VYtYtAn->spacing_row = SPACING_a; P->lf_VYtYtAn->row_stride = P->D->A_ncols;
    P->lf_VYtYtAn->spacing_col = SPACING_a; P->lf_VYtYtAn->col_stride = P->rank;
    /* if (P->lf_VYtYtAn==NULL){ } */}
  L2_AtZn_uu(P->lf_An,P->lf_VYtYt,P->lf_VYtYtAn);
  L2_AtZn_uu(P->lf_Zn,P->lf_VYtWt,P->lf_VYtWtZn);
  //lfprintf(P->lf_VYtYtAn," %% P->lf_VYtYtAn->lf: ");
  //lfprintf(P->lf_VYtWtZn," %% P->lf_VYtWtZn->lf: ");
  if (P->lf_VYWt==NULL){
    P->lf_VYWt = L_handle_make((unsigned long long int)P->D->A_ncols*(unsigned long long int)P->rank); 
    P->lf_VYWt->spacing_row = SPACING_a; P->lf_VYWt->row_stride = P->D->A_ncols;
    P->lf_VYWt->spacing_col = SPACING_a; P->lf_VYWt->col_stride = P->rank;
    /* if (P->lf_VYWt==NULL){ } */}
  L_zero(P->lf_VYWt);
  dtmp = P->D->A_nrows_total; if (dtmp>0){ dra_plustimesequals(P->lf_VYWt->lf,P->lf_VYWt->length,P->lf_VYtYtAn->lf,+1.0/dtmp);}
  dtmp = P->D->Z_nrows_total; if (dtmp>0){ dra_plustimesequals(P->lf_VYWt->lf,P->lf_VYWt->length,P->lf_VYtWtZn->lf,-1.0/dtmp);}
  //lfprintf(P->lf_VYWt," %% P->lf_VYWt->lf: ");
  L_zero(P->lf_Vt);
  dtmp = P->D->A_ncols; if (dtmp>0){ dra_plustimesequals(P->lf_Vt->lf,P->lf_Vt->length,P->lf_VAZt->lf,+1.0/dtmp);}
  dtmp = P->D->Y_ncols; if (dtmp>0){ dra_plustimesequals(P->lf_Vt->lf,P->lf_Vt->length,P->lf_VYWt->lf,-1.0/dtmp);}
  //lfprintf(P->lf_Vt," %% P->lf_Vt->lf: ");
  printf(" %% Single iteration using D: error %0.16f\n",dra_diff(P->lf_Vt->lf,P->lf_V->lf,P->lf_V->length,1));
  /* if (iteration_max<=1){ } */}
  /* %%%%%%%%%%%%%%%% */
  /* Multiple iterations using S */
  /* %%%%%%%%%%%%%%%% */
  if (iteration_max>1){
    pca_stage_a0b(P);
    pca_stage_a1b(P);
    nl=0;
    while (nl<iteration_max){
      L2_AnZt_mm(P->lf_S,P->lf_Vt,P->lf_Vx); L_handle_copy(P->lf_Vt,P->lf_Vx); pca_stage_a1b(P);
      nl++; /* while (nl<iteration_max){ } */}
    if (verbose>1){ lfprintf(P->lf_Vt," %% P->lf_Vt->lf: ");}
    if (verbose>1){ lfprintf(P->lf_R2," %% P->lf_R2->lf: ");}
    printf(" %% Multiple iterations using S: lf_V error %0.16f\n",dra_diff(P->lf_Vt->lf,P->lf_V->lf,P->lf_V->length,1));
    printf(" %% Multiple iterations using S: lf_R1 error %0.16f\n",dra_diff(P->lf_R2->lf,P->lf_R1->lf,P->lf_R1->length,1));
    /* if (iteration_max>1){ } */}
  if (verbose>1){ printf(" %% [finished pca_uu]\n");}
}

void pca_init_test()
{
  /* runs pca_init test */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0); int iteration_max = GLOBAL_TEST_niter;
  int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  struct P_handle *P=NULL;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering pca_init_test]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); 
  dcc_ajdk_set_mx_j_to_mx_b(D); dcc_X_nrows_total(D); dcc_M_mxset(D); dcc_init_QX(D); dcc_ajdk_set_QC_index_local(D); dcc_ajdk_copy_QR_index_global(D); 
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making P_handle.\n");}
  P = P_handle_make(GLOBAL_pca_infix,GLOBAL_pca_out_xdrop,NULL,D,GLOBAL_pca_iteration_num,GLOBAL_pca_iteration_max,GLOBAL_pca_iteration_min,GLOBAL_pca_rank,GLOBAL_pca_tolerance,GLOBAL_B_MLT);
  P_xdrop_init(D,P);  
  if (verbose>1){ printf(" %% P->D->A_rpop_b_total %d P->out_xdrop_nrows %d\n",P->D->A_rpop_b_total,P->out_xdrop_nrows);}
  for (P->nx=0;P->nx<P->iteration_num;P->nx++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*P->iteration_num)/maximum(1,P->nx);
    if (verbose>1){ printf(" %% pre: iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    if (verbose>1){ printf(" %% reset D->QR_index_local and D->QC_index_local to match P->mr_index_local and P->mc_index_local, respectively\n");}
    pca_srt(P);
    if (verbose>1){ printf(" %% Define D->A_umc_j_rmv and D->A_umc_j_rtn, retaining only P->ckeep[P->iteration_num-1-P->nx] entries from D->QC_index_local\n");}
    if (verbose>1){ printf(" %% Define E_[nb]->A_umr_j_rmv and E_[nb]->A_umr_j_rtn, retaining only P->rkeep[P->iteration_num-1-P->nx] entries from D->QR_index_local\n");}
    pca_mxA(P,P->D->A_rpop_b_total - P->rkeep[P->iteration_num-1-P->nx],P->D->A_cpop_b - P->ckeep[P->iteration_num-1-P->nx]);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j_rtn into D->A_bmc_j, and E_[nb]->A_bmr_j_rtn into E_[nb]->A_bmr_j\n");}
    dcc_lrup_mxdup(D);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j into E_[nb]->M_An->mc_j, and E_[nb]->A_bmr_j into E_[nb]->M_An->mr_j, etc.\n");} 
    dcc_M_mxset(D);
    if (verbose>1){ printf(" %% mid: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    if (iteration_max<=1){
      pca_stage_a0(P);
      pca_stage_a1(P);
      //lfprintf(P->lf_V ," P->lf_V1: ");
      pca_stage_b(P);
      pca_stage_c1(P);
      pca_stage_c2(P);
      pca_stage_d(P);
      /* if (iteration_max<=1){ } */}
    if (iteration_max>1){
      pca_stage_a0(P);
      pca_stage_a1(P);
      nl=0;
      while (nl<iteration_max){
        if (verbose>2){ pca_printf(verbose,P);}
        pca_stage_b(P); pca_stage_c1(P); pca_stage_c2(P); pca_stage_d(P); pca_stage_a1(P);
	if (verbose>2){ 
	  sprintf(tmpchar," %% P->nx %0.3d nl %0.3d P->lf_R1->lf: ",P->nx,nl); 
	  lfprintf(P->lf_R1,tmpchar);
	  /* if (verbose>2){ } */}
        nl++; /* while (nl<iteration_max){ } */}
      pca_stage_e(P);
      /* if (iteration_max>1){ } */}
    if (error_check){ pca_uu(P);}
    if (verbose>1){ printf(" %% pos: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    /* for (P->nx=0;P->nx<P->iteration_num;P->nx++){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  wkspace_printf();
  if (verbose>1){ printf(" %% [finished pca_init_test]\n");}
}

double pca_R_rfro(struct L_handle *lf_R1,struct L_handle *lf_R2)
{
  int verbose=0;
  double ddif=0,dfro=0,drfro=0;
  int nr=0,nc=0;
  if (verbose){ printf(" %% [entering pca_R_rfro]\n");}
  if (verbose){ lfprintf(lf_R1," %% lf_R1: "); lfprintf(lf_R2," %% lf_R2: ");}
  for (nr=0;nr<lf_R1->row_stride;nr++){ for (nc=0;nc<lf_R1->col_stride;nc++){
    ddif += pow((*L2_get(lf_R1,0,0,nr,0,0,nc))-(*L2_get(lf_R2,0,0,nr,0,0,nc)),2);
    dfro += pow((*L2_get(lf_R1,0,0,nr,0,0,nc)),2);
    /* for (nr=0;nr<lf_R1->row_stride;nr++){ for (nc=0;nc<lf_R1->col_stride;nc++){ }} */}}
  ddif = sqrt(ddif); if (dfro<=0){ dfro = 1.0;} dfro = sqrt(dfro);
  drfro = ddif/dfro; if (drfro>1e12){ drfro=1e12;}
  if (verbose){ printf(" %% [finished pca_R_rfro]: drfro %0.16f\n",drfro);}
  return drfro;
}

void pca_iter_test()
{
  /* runs pca_iter test */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0); int iteration_max = GLOBAL_TEST_niter;
  int continue_flag=0; double drfro=0;
  int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  struct P_handle *P=NULL;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering pca_iter_test]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); 
  dcc_ajdk_set_mx_j_to_mx_b(D); dcc_X_nrows_total(D); dcc_M_mxset(D); dcc_init_QX(D); dcc_ajdk_set_QC_index_local(D); dcc_ajdk_copy_QR_index_global(D); 
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making P_handle.\n");}
  P = P_handle_make(GLOBAL_pca_infix,GLOBAL_pca_out_xdrop,NULL,D,GLOBAL_pca_iteration_num,GLOBAL_pca_iteration_max,GLOBAL_pca_iteration_min,GLOBAL_pca_rank,GLOBAL_pca_tolerance,GLOBAL_B_MLT);
  P_xdrop_init(D,P);  
  if (verbose>1){ printf(" %% P->D->A_rpop_b_total %d P->out_xdrop_nrows %d\n",P->D->A_rpop_b_total,P->out_xdrop_nrows);}
  for (P->nx=0;P->nx<P->iteration_num;P->nx++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*P->iteration_num)/maximum(1,P->nx);
    if (verbose>1){ printf(" %% pre: iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    if (verbose>1){ printf(" %% reset D->QR_index_local and D->QC_index_local to match P->mr_index_local and P->mc_index_local, respectively\n");}
    pca_srt(P);
    if (verbose>1){ printf(" %% Define D->A_umc_j_rmv and D->A_umc_j_rtn, retaining only P->ckeep[P->iteration_num-1-P->nx] entries from D->QC_index_local\n");}
    if (verbose>1){ printf(" %% Define E_[nb]->A_umr_j_rmv and E_[nb]->A_umr_j_rtn, retaining only P->rkeep[P->iteration_num-1-P->nx] entries from D->QR_index_local\n");}
    pca_mxA(P,P->D->A_rpop_b_total - P->rkeep[P->iteration_num-1-P->nx],P->D->A_cpop_b - P->ckeep[P->iteration_num-1-P->nx]);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j_rtn into D->A_bmc_j, and E_[nb]->A_bmr_j_rtn into E_[nb]->A_bmr_j\n");}
    dcc_lrup_mxdup(D);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j into E_[nb]->M_An->mc_j, and E_[nb]->A_bmr_j into E_[nb]->M_An->mr_j, etc.\n");} 
    dcc_M_mxset(D);
    if (verbose>-1){ printf(" %% mid: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    if (P->nx==0){ if (verbose>1){ printf(" %% initializing P->lf_V\n");} pca_stage_a0(P);} 
    else /* if (P->nx>0) */{ /* do nothing */}
    pca_stage_a1(P);
    continue_flag=1; nl=0;
    while (continue_flag){
      pca_stage_b(P); pca_stage_c1(P); pca_stage_c2(P); pca_stage_d(P); pca_stage_a1(P);
      drfro = pca_R_rfro(P->lf_R1,P->lf_R2); L_handle_copy(P->lf_R2,P->lf_R1);
      if (verbose>2){ 
	sprintf(tmpchar," %% P->nx %0.3d nl %0.3d drfro %0.16f P->lf_R1->lf: ",P->nx,nl,drfro); 
	lfprintf(P->lf_R1,tmpchar);
	/* if (verbose>2){ } */}
      nl++; continue_flag = (iteration_max>1) && (nl<=iteration_max) && (drfro>P->tolerance);
      /* while (continue_flag){ } */}
    if (verbose>-1){ printf(" %% P->nx %d: nl %d drfro %0.16f exiting inner loop\n",P->nx,nl,drfro);}
    pca_stage_e(P);
    if (verbose>1){ printf(" %% pos: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    /* for (P->nx=0;P->nx<P->iteration_num;P->nx++){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  P_handle_printf(verbose,P," %% P: ");
  P_handle_dmp(P);
  wkspace_printf();
  if (verbose>1){ printf(" %% [finished pca_iter_test]\n");}
}

void pca_driver()
{
  /* runs main pca driver */
  int verbose=GLOBAL_verbose; int nl=0,nr=0; double ct=0,rt=0,r=0,it=0,et=0; int iteration_max = GLOBAL_TEST_niter;
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  struct P_handle *P=NULL;
  int continue_flag=0; double drfro=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering pca_driver]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_load(&D,&E_,&F_); 
  dcc_ajdk_set_mx_j_to_mx_b(D); dcc_X_nrows_total(D); dcc_M_mxset(D); dcc_init_QX(D); dcc_ajdk_set_QC_index_local(D); dcc_ajdk_copy_QR_index_global(D); 
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making P_handle.\n");}
  P = P_handle_make(GLOBAL_pca_infix,GLOBAL_pca_out_xdrop,NULL,D,GLOBAL_pca_iteration_num,GLOBAL_pca_iteration_max,GLOBAL_pca_iteration_min,GLOBAL_pca_rank,GLOBAL_pca_tolerance,GLOBAL_B_MLT);
  P_xdrop_init(D,P);  
  if (verbose>1){ printf(" %% P->D->A_rpop_b_total %d P->out_xdrop_nrows %d\n",P->D->A_rpop_b_total,P->out_xdrop_nrows);}
  for (P->nx=0;P->nx<P->iteration_num;P->nx++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*P->iteration_num)/maximum(1,P->nx);
    if (verbose>1){ printf(" %% pre: iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    if (verbose>1){ printf(" %% reset D->QR_index_local and D->QC_index_local to match P->mr_index_local and P->mc_index_local, respectively\n");}
    pca_srt(P);
    if (verbose>1){ printf(" %% Define D->A_umc_j_rmv and D->A_umc_j_rtn, retaining only P->ckeep[P->iteration_num-1-P->nx] entries from D->QC_index_local\n");}
    if (verbose>1){ printf(" %% Define E_[nb]->A_umr_j_rmv and E_[nb]->A_umr_j_rtn, retaining only P->rkeep[P->iteration_num-1-P->nx] entries from D->QR_index_local\n");}
    pca_mxA(P,P->D->A_rpop_b_total - P->rkeep[P->iteration_num-1-P->nx],P->D->A_cpop_b - P->ckeep[P->iteration_num-1-P->nx]);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j_rtn into D->A_bmc_j, and E_[nb]->A_bmr_j_rtn into E_[nb]->A_bmr_j\n");}
    dcc_lrup_mxdup(D);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j into E_[nb]->M_An->mc_j, and E_[nb]->A_bmr_j into E_[nb]->M_An->mr_j, etc.\n");} 
    dcc_M_mxset(D);
    if (verbose>-1){ printf(" %% mid: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    if (P->nx==0){ if (verbose>1){ printf(" %% initializing P->lf_V\n");} pca_stage_a0(P);} 
    else /* if (P->nx>0) */{ /* do nothing */}
    pca_stage_a1(P);
    continue_flag=1; nl=0;
    while (continue_flag){
      pca_stage_b(P); pca_stage_c1(P); pca_stage_c2(P); pca_stage_d(P); pca_stage_a1(P);
      drfro = pca_R_rfro(P->lf_R1,P->lf_R2); L_handle_copy(P->lf_R2,P->lf_R1);
      if (verbose>2){ 
	sprintf(tmpchar," %% P->nx %0.3d nl %0.3d drfro %0.16f P->lf_R1->lf diagonal: ",P->nx,nl,drfro); 
	for (nr=0;nr<P->rank;nr++){ sprintf(tmpchar,"%s %0.16f",tmpchar,*L2_get(P->lf_R1,0,0,nr,0,0,nr));} 
	sprintf(tmpchar,"%s\n",tmpchar); printf("%s",tmpchar);
	/* if (verbose>2){ } */}
      if (verbose>2){ 
	sprintf(tmpchar," %% P->nx %0.3d nl %0.3d drfro %0.16f P->lf_R1->lf: ",P->nx,nl,drfro); 
	lfprintf(P->lf_R1,tmpchar);
	/* if (verbose>2){ } */}
      nl++; continue_flag = (iteration_max>1) && (nl<=iteration_max) && (drfro>P->tolerance);
      /* while (continue_flag){ } */}
    if (verbose>-1){ printf(" %% P->nx %d: nl %d drfro %0.16f exiting inner loop\n",P->nx,nl,drfro);}
    pca_stage_e(P);
    if (verbose>1){ printf(" %% pos: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    /* for (P->nx=0;P->nx<P->iteration_num;P->nx++){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  P_handle_printf(verbose,P," %% P: ");
  P_handle_dmp(P);
  wkspace_printf();
  if (verbose>1){ printf(" %% [finished pca_driver]\n");}
}

void pca_proj_test()
{
  /* runs pca projection test */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0); int iteration_max = GLOBAL_TEST_niter;
  int continue_flag=0; double drfro=0;
  int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  struct P_handle *P=NULL;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering pca_proj_test]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); 
  dcc_ajdk_set_mx_j_to_mx_b(D); dcc_X_nrows_total(D); dcc_M_mxset(D); dcc_init_QX(D); dcc_ajdk_set_QC_index_local(D); dcc_ajdk_copy_QR_index_global(D); 
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making P_handle.\n");}
  P = P_handle_make(GLOBAL_pca_infix,GLOBAL_pca_out_xdrop,GLOBAL_pca_V_,D,GLOBAL_pca_iteration_num,GLOBAL_pca_iteration_max,GLOBAL_pca_iteration_min,GLOBAL_pca_rank,GLOBAL_pca_tolerance,GLOBAL_B_MLT);
  P_xdrop_init(D,P); for (nl=0;nl<P->iteration_num;nl++){ P->rkeep[nl] = P->D->A_rpop_b_total; P->rdrop[nl] = 0;}
  if (verbose>1){ printf(" %% P->D->A_rpop_b_total %d P->out_xdrop_nrows %d\n",P->D->A_rpop_b_total,P->out_xdrop_nrows);}
  for (P->nx=0;P->nx<P->iteration_num;P->nx++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*P->iteration_num)/maximum(1,P->nx);
    if (verbose>1){ printf(" %% pre: iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    if (verbose>1){ printf(" %% reset D->QR_index_local and D->QC_index_local to match P->mr_index_local and P->mc_index_local, respectively\n");}
    pca_srt(P);
    if (verbose>1){ printf(" %% Define D->A_umc_j_rmv and D->A_umc_j_rtn, retaining only P->ckeep[P->iteration_num-1-P->nx] entries from D->QC_index_local\n");}
    if (verbose>1){ printf(" %% Define E_[nb]->A_umr_j_rmv and E_[nb]->A_umr_j_rtn, retaining only P->rkeep[P->iteration_num-1-P->nx] entries from D->QR_index_local\n");}
    pca_mxA(P,P->D->A_rpop_b_total - P->rkeep[P->iteration_num-1-P->nx],P->D->A_cpop_b - P->ckeep[P->iteration_num-1-P->nx]);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j_rtn into D->A_bmc_j, and E_[nb]->A_bmr_j_rtn into E_[nb]->A_bmr_j\n");}
    dcc_lrup_mxdup(D);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j into E_[nb]->M_An->mc_j, and E_[nb]->A_bmr_j into E_[nb]->M_An->mr_j, etc.\n");} 
    dcc_M_mxset(D);
    if (verbose>-1){ printf(" %% mid: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    pca_stage_a9(P);
    pca_stage_b9(P);
    pca_stage_e9(P);
    if (verbose>1){ printf(" %% pos: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    /* for (P->nx=0;P->nx<P->iteration_num;P->nx++){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  P_handle_printf(verbose,P," %% P: ");
  P_handle_proj_dmp(P);
  wkspace_printf();
  if (verbose>1){ printf(" %% [finished pca_proj_test]\n");}
}

void pca_proj_driver()
{
  /* runs main pca proj_driver */
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0; int iteration_max = GLOBAL_TEST_niter;
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  struct P_handle *P=NULL;
  int continue_flag=0; double drfro=0;
  if (verbose>1){ printf(" %% [entering pca_proj_driver]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_load(&D,&E_,&F_); 
  dcc_ajdk_set_mx_j_to_mx_b(D); dcc_X_nrows_total(D); dcc_M_mxset(D); dcc_init_QX(D); dcc_ajdk_set_QC_index_local(D); dcc_ajdk_copy_QR_index_global(D); 
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making P_handle.\n");}
  P = P_handle_make(GLOBAL_pca_infix,GLOBAL_pca_out_xdrop,GLOBAL_pca_V_,D,GLOBAL_pca_iteration_num,GLOBAL_pca_iteration_max,GLOBAL_pca_iteration_min,GLOBAL_pca_rank,GLOBAL_pca_tolerance,GLOBAL_B_MLT);
  P_xdrop_init(D,P); for (nl=0;nl<P->iteration_num;nl++){ P->rkeep[nl] = P->D->A_rpop_b_total; P->rdrop[nl] = 0;}
  if (verbose>1){ printf(" %% P->D->A_rpop_b_total %d P->out_xdrop_nrows %d\n",P->D->A_rpop_b_total,P->out_xdrop_nrows);}
  /* if (verbose>1){ printf(" %% P->rank %d P->iteration_num %d\n",P->rank,P->iteration_num);} */
  for (P->nx=0;P->nx<P->iteration_num;P->nx++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*P->iteration_num)/maximum(1,P->nx);
    if (verbose>1){ printf(" %% pre: iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    if (verbose>1){ printf(" %% reset D->QR_index_local and D->QC_index_local to match P->mr_index_local and P->mc_index_local, respectively\n");}
    pca_srt(P);
    if (verbose>1){ printf(" %% Define D->A_umc_j_rmv and D->A_umc_j_rtn, retaining only P->ckeep[P->iteration_num-1-P->nx] entries from D->QC_index_local\n");}
    if (verbose>1){ printf(" %% Define E_[nb]->A_umr_j_rmv and E_[nb]->A_umr_j_rtn, retaining only P->rkeep[P->iteration_num-1-P->nx] entries from D->QR_index_local\n");}
    pca_mxA(P,P->D->A_rpop_b_total - P->rkeep[P->iteration_num-1-P->nx],P->D->A_cpop_b - P->ckeep[P->iteration_num-1-P->nx]);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j_rtn into D->A_bmc_j, and E_[nb]->A_bmr_j_rtn into E_[nb]->A_bmr_j\n");}
    dcc_lrup_mxdup(D);
    if (verbose>1){ printf(" %% Copy D->A_bmc_j into E_[nb]->M_An->mc_j, and E_[nb]->A_bmr_j into E_[nb]->M_An->mr_j, etc.\n");} 
    dcc_M_mxset(D);
    if (verbose>-1){ printf(" %% mid: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    pca_stage_a9(P);
    pca_stage_b9(P);
    pca_stage_e9(P);
    if (verbose>1){ printf(" %% pos: iteration %.5d/%.5d, D->A_rpop_j_total %d / D->A_rpop_b_total %d / D->A_nrows_total %d ; D->A_cpop_j %d / D->A_cpop_b %d / D->A_ncols %d \n",P->nx,P->iteration_num,D->A_rpop_j_total,D->A_rpop_b_total,D->A_nrows_total,D->A_cpop_j,D->A_cpop_b,D->A_ncols);}
    /* for (P->nx=0;P->nx<P->iteration_num;P->nx++){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  P_handle_printf(verbose,P," %% P: ");
  P_handle_proj_dmp(P);
  wkspace_printf();
  if (verbose>1){ printf(" %% [finished pca_proj_driver]\n");}
}
