#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void dcg_An_ajdk(struct dcg_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_;
  int nb=0; struct dcg_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcg_An_ajdk]\n");}
  if (verbose){ printf(" %% calculating E_[nb]->lf_An_ajdk, E_[nb]->lf_Zn_ajdk.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_An,D->A_ajdk,&(E->lf_An_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Yn,D->Y_ajdk,&(E->lf_Yn_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Zn,D->A_ajdk,&(E->lf_Zn_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Wn,D->Y_ajdk,&(E->lf_Wn_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Jn,D->A_ajdk,&(E->lf_Jn_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Kn,D->A_ajdk,&(E->lf_Kn_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Un,D->Y_ajdk,&(E->lf_Un_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Vn,D->A_ajdk,&(E->lf_Vn_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% An_ajdk_v: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){ 
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      printf(" %% nb %d\n",nb);
      lfprintf(E->lf_An_ajdk," %% lf_An_ajdk->lf: ");
      lfprintf(E->lf_Yn_ajdk," %% lf_Yn_ajdk->lf: ");
      lfprintf(E->lf_Zn_ajdk," %% lf_Zn_ajdk->lf: ");
      lfprintf(E->lf_Wn_ajdk," %% lf_Wn_ajdk->lf: ");
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      lfprintf(E->lf_Jn_ajdk," %% lf_Jn_ajdk->lf: ");
      lfprintf(E->lf_Kn_ajdk," %% lf_Kn_ajdk->lf: ");
      lfprintf(E->lf_Un_ajdk," %% lf_Un_ajdk->lf: ");
      lfprintf(E->lf_Vn_ajdk," %% lf_Vn_ajdk->lf: ");
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcg_An_ajdk]\n");}
}

void dcg_lf_ZtSn(struct dcg_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_;
  int nb=0; struct dcg_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcg_lf_ZtSn]\n");}
  if (verbose){ printf(" %% calculating E_[nb]->lf_AtTn, E_[nb]->lf_ZtSn.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_At,E->M_Tt,NULL,NULL,NULL,&(E->lf_AtTn));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Yt,E->M_Tt,NULL,NULL,NULL,&(E->lf_YtTn));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Zt,E->M_St,NULL,NULL,NULL,&(E->lf_ZtSn));
      GLOBAL_pthread_toc(); /* if bother */} 
    if (E->Z_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Wt,E->M_St,NULL,NULL,NULL,&(E->lf_WtSn));
      GLOBAL_pthread_toc(); /* if bother */} 
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Jt,E->M_Rt,NULL,NULL,NULL,&(E->lf_JtRn));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Kt,E->M_Tt,NULL,NULL,NULL,&(E->lf_KtTn));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->A_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Ut,E->M_Rt,NULL,NULL,NULL,&(E->lf_UtRn));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Vt,E->M_St,NULL,NULL,NULL,&(E->lf_VtSn));
      GLOBAL_pthread_toc(); /* if bother */} 
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% ZtSn_vv: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      printf(" %% nb %d\n",nb);
      bprintf(E->M_An->mr_j,D->bitj,1,E->A_nrows," %% A_bmr_j: ");
      bprintf(E->M_An->mc_j,D->bitj,1,D->A_ncols," %% A_bmc_j: ");
      bprintf(E->M_Jn->mr_j,D->bitj,1,E->A_nrows," %% J_bmr_j: ");
      bprintf(E->M_Jn->mc_j,D->bitj,1,D->A_ncols," %% J_bmr_j: ");
      bprintf(E->M_Kn->mr_j,D->bitj,1,E->A_nrows," %% K_bmr_j: ");
      bprintf(E->M_Kn->mc_j,D->bitj,1,D->A_ncols," %% K_bmc_j: ");
      lfprintf(E->lf_AtTn," %% lf_AtTn: ");
      lfprintf(E->lf_YtTn," %% lf_YtTn: ");
      lfprintf(E->lf_ZtSn," %% lf_ZtSn: ");
      lfprintf(E->lf_WtSn," %% lf_WtSn: ");
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      lfprintf(E->lf_JtRn," %% lf_JtRn: ");
      lfprintf(E->lf_KtTn," %% lf_KtTn: ");
      lfprintf(E->lf_UtRn," %% lf_UtRn: ");
      lfprintf(E->lf_VtSn," %% lf_VtSn: ");
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcg_lf_ZtSn]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void dcg_lf_D_AtTn_ZtSn_vv(struct dcg_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_; struct dcg_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; struct dcg_single *E_nb1=NULL,*E_nb2=NULL; struct dcg_double *F=NULL;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcg_lf_D_AtTn_ZtSn_vv]\n");}
  if (verbose){ printf(" %% calculating F_[nbx]->lf_D_AtTn_ZtSn, F_[nbx]->lf_D_AtTn_AtTn.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx=nb1+nb2*nbins; E_nb1 = E_[nb1]; E_nb2 = E_[nb2]; F = F_[nbx];
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->Y_cbother){
      	GLOBAL_pthread_tic();
      	wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yt,E_nb1->M_Tt,E_nb2->M_Wt,E_nb2->M_St,D->Y_ajdk,E_nb1->lf_YtTn,E_nb2->lf_WtSn,&(F->lf_D_YtTn_WtSn_vv));
      	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->Y_cbother){
      	if (nb1==nb2){
      	  GLOBAL_pthread_tic();
      	  wrap_D_AtTn_AtTn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yt,E_nb1->M_Tt,E_nb2->M_Yt,E_nb2->M_Tt,D->Y_ajdk,E_nb1->lf_YtTn,E_nb2->lf_YtTn,&(F->lf_D_YtTn_YtTn_vv));
      	  GLOBAL_pthread_toc();
      	  /* if (nb1==nb2){ } */}
      	if (nb1!=nb2){
      	  GLOBAL_pthread_tic();
      	  wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yt,E_nb1->M_Tt,E_nb2->M_Yt,E_nb2->M_Tt,D->Y_ajdk,E_nb1->lf_YtTn,E_nb2->lf_YtTn,&(F->lf_D_YtTn_YtTn_vv));
      	  GLOBAL_pthread_toc();
      	  /* if (nb1!=nb2){ } */}
      	/* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
      	GLOBAL_pthread_tic();
      	wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb2->M_Zt,E_nb2->M_St,D->A_ajdk,E_nb1->lf_AtTn,E_nb2->lf_ZtSn,&(F->lf_D_AtTn_ZtSn_vv));
      	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){
      	if (nb1==nb2){
      	  GLOBAL_pthread_tic();
      	  wrap_D_AtTn_AtTn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb2->M_At,E_nb2->M_Tt,D->A_ajdk,E_nb1->lf_AtTn,E_nb2->lf_AtTn,&(F->lf_D_AtTn_AtTn_vv));
      	  GLOBAL_pthread_toc();
      	  /* if (nb1==nb2){ } */}
      	if (nb1!=nb2){
      	  GLOBAL_pthread_tic();
      	  wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb2->M_At,E_nb2->M_Tt,D->A_ajdk,E_nb1->lf_AtTn,E_nb2->lf_AtTn,&(F->lf_D_AtTn_AtTn_vv));
      	  GLOBAL_pthread_toc();
      	  /* if (nb1!=nb2){ } */}
      	/* if bother */}
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->Y_cbother){
      	GLOBAL_pthread_tic();
      	wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Ut,E_nb1->M_Rt,E_nb2->M_Wt,E_nb2->M_St,D->Y_ajdk,E_nb1->lf_UtRn,E_nb2->lf_WtSn,&(F->lf_D_UtRn_WtSn_vv));
      	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->Y_cbother){
      	GLOBAL_pthread_tic();
      	wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Ut,E_nb1->M_Rt,E_nb2->M_Yt,E_nb2->M_Tt,D->Y_ajdk,E_nb1->lf_UtRn,E_nb2->lf_YtTn,&(F->lf_D_UtRn_YtTn_vv));
      	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Jt,E_nb1->M_Rt,E_nb2->M_Zt,E_nb2->M_St,D->A_ajdk,E_nb1->lf_JtRn,E_nb2->lf_ZtSn,&(F->lf_D_JtRn_ZtSn_vv));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){
      	GLOBAL_pthread_tic();
      	wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Jt,E_nb1->M_Rt,E_nb2->M_At,E_nb2->M_Tt,D->A_ajdk,E_nb1->lf_JtRn,E_nb2->lf_AtTn,&(F->lf_D_JtRn_AtTn_vv));
      	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
      	GLOBAL_pthread_tic();
      	wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Kt,E_nb1->M_Tt,E_nb2->M_Vt,E_nb2->M_St,D->A_ajdk,E_nb1->lf_KtTn,E_nb2->lf_VtSn,&(F->lf_D_KtTn_VtSn_vv));
      	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){
      	if (nb1==nb2){
      	  GLOBAL_pthread_tic();
      	  wrap_D_AtTn_AtTn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Kt,E_nb1->M_Tt,E_nb2->M_Kt,E_nb2->M_Tt,D->A_ajdk,E_nb1->lf_KtTn,E_nb2->lf_KtTn,&(F->lf_D_KtTn_KtTn_vv));
      	  GLOBAL_pthread_toc();
      	  /* if (nb1==nb2){ } */}
      	if (nb1!=nb2){
      	  GLOBAL_pthread_tic();
      	  wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Kt,E_nb1->M_Tt,E_nb2->M_Kt,E_nb2->M_Tt,D->A_ajdk,E_nb1->lf_KtTn,E_nb2->lf_KtTn,&(F->lf_D_KtTn_KtTn_vv));
      	  GLOBAL_pthread_toc();
      	  /* if (nb1!=nb2){ } */}
      	/* if bother */}
      /*  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% D_AtTn_ZtSn_vv: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A_bmr_j: ");
	bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A_bmc_j: ");
	bprintf(E_nb1->M_Jn->mr_j,D->bitj,1,E_nb1->A_nrows," %% J_bmr_j: ");
	bprintf(E_nb1->M_Jn->mc_j,D->bitj,1,D->A_ncols," %% J_bmc_j: ");
	bprintf(E_nb1->M_Kn->mr_j,D->bitj,1,E_nb1->A_nrows," %% K_bmr_j: ");
	bprintf(E_nb1->M_Kn->mc_j,D->bitj,1,D->A_ncols," %% K_bmc_j: ");
	bprintf(E_nb1->M_Tn->mr_j,D->bitj,1,E_nb1->A_nrows," %% T_bmr_j: ");
	bprintf(E_nb1->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T_bmc_j: ");
	bprintf(E_nb1->M_Rn->mr_j,D->bitj,1,E_nb1->A_nrows," %% R_bmr_j: ");
	bprintf(E_nb1->M_Rn->mc_j,D->bitj,1,D->T_ncols," %% R_bmc_j: ");
	lfprintf(F->lf_D_YtTn_WtSn_vv," %% lf_D_YtTn_WtSn_vv: ");
	lfprintf(F->lf_D_YtTn_YtTn_vv," %% lf_D_YtTn_YtTn_vv: ");
	lfprintf(F->lf_D_AtTn_ZtSn_vv," %% lf_D_AtTn_ZtSn_vv: ");
	lfprintf(F->lf_D_AtTn_AtTn_vv," %% lf_D_AtTn_AtTn_vv: ");
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	lfprintf(F->lf_D_UtRn_YtTn_vv," %% lf_D_UtRn_YtTn_vv: ");
	lfprintf(F->lf_D_UtRn_WtSn_vv," %% lf_D_UtRn_WtSn_vv: ");
	lfprintf(F->lf_D_KtTn_KtTn_vv," %% lf_D_KtTn_KtTn_vv: ");
	lfprintf(F->lf_D_KtTn_VtSn_vv," %% lf_D_KtTn_VtSn_vv: ");
	lfprintf(F->lf_D_JtRn_AtTn_vv," %% lf_D_JtRn_AtTn_vv: ");
	lfprintf(F->lf_D_JtRn_ZtSn_vv," %% lf_D_JtRn_ZtSn_vv: ");
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcg_lf_D_AtTn_ZtSn_vv]\n");}
}

void dcg_lf_D_AtTn_ZtSn_uu(struct dcg_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_; struct dcg_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct dcg_single *E_nb1=NULL,*E_nb2=NULL; struct dcg_double *F=NULL,*F_trn=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcg_lf_D_AtTn_ZtSn_uu]\n");}
  if (verbose){ printf(" %% calculating F_[nbx]->lf_D_AtTn_ZtSn_uu, F_[nbx]->lf_D_AtTn_AtTn_uu.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
      if (D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yt,E_nb1->M_Tt,E_nb2->M_Wt,E_nb2->M_St,D->Y_ajdk,&(F->lf_D_YtTn_WtSn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	if (nb1==nb2){
	  GLOBAL_pthread_tic(); 
	  wrap_D_AtTn_AtTn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yt,E_nb1->M_Tt,E_nb2->M_Yt,E_nb2->M_Tt,D->Y_ajdk,&(F->lf_D_YtTn_YtTn_uu));
	  GLOBAL_pthread_toc(); 
	  /* if (nb1==nb2){ } */}
	if (nb1!=nb2){
	  GLOBAL_pthread_tic(); 
	  wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yt,E_nb1->M_Tt,E_nb2->M_Yt,E_nb2->M_Tt,D->Y_ajdk,&(F->lf_D_YtTn_YtTn_uu));
	  GLOBAL_pthread_toc(); 
	  /* if (nb1!=nb2){ } */}
	/* if bother */}
      if (D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb2->M_Zt,E_nb2->M_St,D->A_ajdk,&(F->lf_D_AtTn_ZtSn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	if (nb1==nb2){
	  GLOBAL_pthread_tic(); 
	  wrap_D_AtTn_AtTn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb2->M_At,E_nb2->M_Tt,D->A_ajdk,&(F->lf_D_AtTn_AtTn_uu));
	  GLOBAL_pthread_toc(); 
	  /* if (nb1==nb2){ } */}
	if (nb1!=nb2){
	  GLOBAL_pthread_tic(); 
	  wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb2->M_At,E_nb2->M_Tt,D->A_ajdk,&(F->lf_D_AtTn_AtTn_uu));
	  GLOBAL_pthread_toc(); 
	  /* if (nb1!=nb2){ } */}
	/* if bother */}
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->Y_cbother){
	GLOBAL_pthread_tic();
	wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Ut,E_nb1->M_Rt,E_nb2->M_Wt,E_nb2->M_St,D->Y_ajdk,&(F->lf_D_UtRn_WtSn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->Y_cbother){
	GLOBAL_pthread_tic();
	wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Ut,E_nb1->M_Rt,E_nb2->M_Yt,E_nb2->M_Tt,D->Y_ajdk,&(F->lf_D_UtRn_YtTn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Jt,E_nb1->M_Rt,E_nb2->M_Zt,E_nb2->M_St,D->A_ajdk,&(F->lf_D_JtRn_ZtSn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Jt,E_nb1->M_Rt,E_nb2->M_At,E_nb2->M_Tt,D->A_ajdk,&(F->lf_D_JtRn_AtTn_uu));
	GLOBAL_pthread_toc(); /* if bother */}      
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Kt,E_nb1->M_Tt,E_nb2->M_Vt,E_nb2->M_St,D->A_ajdk,&(F->lf_D_KtTn_VtSn_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){
	if (nb1==nb2){
	  GLOBAL_pthread_tic();
	  wrap_D_AtTn_AtTn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Kt,E_nb1->M_Tt,E_nb2->M_Kt,E_nb2->M_Tt,D->A_ajdk,&(F->lf_D_KtTn_KtTn_uu));
	  GLOBAL_pthread_toc();
	  /* if (nb1==nb2){ } */}
	if (nb1!=nb2){
	  GLOBAL_pthread_tic();
	  wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Kt,E_nb1->M_Tt,E_nb2->M_Kt,E_nb2->M_Tt,D->A_ajdk,&(F->lf_D_KtTn_KtTn_uu));
	  GLOBAL_pthread_toc();
	  /* if (nb1!=nb2){ } */}
	/* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% D_AtTn_ZtSn_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A_bmr_j: ");
	bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A_bmc_j: ");
	bprintf(E_nb1->M_Jn->mr_j,D->bitj,1,E_nb1->A_nrows," %% J_bmr_j: ");
	bprintf(E_nb1->M_Jn->mc_j,D->bitj,1,D->A_ncols," %% J_bmc_j: ");
	bprintf(E_nb1->M_Kn->mr_j,D->bitj,1,E_nb1->A_nrows," %% K_bmr_j: ");
	bprintf(E_nb1->M_Kn->mc_j,D->bitj,1,D->A_ncols," %% K_bmc_j: ");
	bprintf(E_nb1->M_Tn->mr_j,D->bitj,1,E_nb1->A_nrows," %% T_bmr_j: ");
	bprintf(E_nb1->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T_bmc_j: ");
	bprintf(E_nb1->M_Rn->mr_j,D->bitj,1,E_nb1->A_nrows," %% R_bmr_j: ");
	bprintf(E_nb1->M_Rn->mc_j,D->bitj,1,D->T_ncols," %% R_bmc_j: ");
	lfprintf(F->lf_D_YtTn_WtSn_uu," %% lf_D_YtTn_WtSn_uu: ");
	lfprintf(F->lf_D_YtTn_YtTn_uu," %% lf_D_YtTn_YtTn_uu: ");
	lfprintf(F->lf_D_AtTn_ZtSn_uu," %% lf_D_AtTn_ZtSn_uu: ");
	lfprintf(F->lf_D_AtTn_AtTn_uu," %% lf_D_AtTn_AtTn_uu: ");
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	lfprintf(F->lf_D_UtRn_YtTn_uu," %% lf_D_UtRn_YtTn_uu: ");
	lfprintf(F->lf_D_UtRn_WtSn_uu," %% lf_D_UtRn_WtSn_uu: ");
	lfprintf(F->lf_D_KtTn_KtTn_uu," %% lf_D_KtTn_KtTn_uu: ");
	lfprintf(F->lf_D_KtTn_VtSn_uu," %% lf_D_KtTn_VtSn_uu: ");
	lfprintf(F->lf_D_JtRn_AtTn_uu," %% lf_D_JtRn_AtTn_uu: ");
	lfprintf(F->lf_D_JtRn_ZtSn_uu," %% lf_D_JtRn_ZtSn_uu: ");
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcg_lf_D_AtTn_ZtSn_uu]\n");}
}

void dcg_lf_D_AtTn_ZtSn_error(int verbose,struct dcg_ajdk *D)
{
  int nbins=D->nbins; struct dcg_single **E_=D->E_; struct dcg_double **F_=D->F_;
  int nbx=0,nb1=0,nb2=0,length=0; struct dcg_single *E_nb1=NULL,*E_nb2=NULL; struct dcg_double *F=NULL;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      length = D->A_ncols * D->T_ncols;
      printf(" %% nb1 %d nb2 %d length %d*%d=%d\n",nb1,nb2,D->A_ncols,D->T_ncols,length);
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_D_YtTn_WtSn_vv - lf_D_YtTn_WtSn_uu %0.16f\n",dra_diff(F->lf_D_YtTn_WtSn_vv->lf,F->lf_D_YtTn_WtSn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_D_YtTn_YtTn_vv - lf_D_YtTn_YtTn_uu %0.16f\n",dra_diff(F->lf_D_YtTn_YtTn_vv->lf,F->lf_D_YtTn_YtTn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_D_AtTn_ZtSn_vv - lf_D_AtTn_ZtSn_uu %0.16f\n",dra_diff(F->lf_D_AtTn_ZtSn_vv->lf,F->lf_D_AtTn_ZtSn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_D_AtTn_AtTn_vv - lf_D_AtTn_AtTn_uu %0.16f\n",dra_diff(F->lf_D_AtTn_AtTn_vv->lf,F->lf_D_AtTn_AtTn_uu->lf,length,1));
	/* if bother */}
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_D_UtRn_WtSn_vv - lf_D_UtRn_WtSn_uu %0.16f\n",dra_diff(F->lf_D_UtRn_WtSn_vv->lf,F->lf_D_UtRn_WtSn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_D_UtRn_UtRn_vv - lf_D_UtRn_UtRn_uu %0.16f\n",dra_diff(F->lf_D_UtRn_YtTn_vv->lf,F->lf_D_UtRn_YtTn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_D_KtTn_VtSn_vv - lf_D_KtTn_VtSn_uu %0.16f\n",dra_diff(F->lf_D_KtTn_VtSn_vv->lf,F->lf_D_KtTn_VtSn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_D_KtTn_KtTn_vv - lf_D_KtTn_KtTn_uu %0.16f\n",dra_diff(F->lf_D_KtTn_KtTn_vv->lf,F->lf_D_KtTn_KtTn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_D_JtRn_ZtSn_vv - lf_D_JtRn_ZtSn_uu %0.16f\n",dra_diff(F->lf_D_JtRn_ZtSn_vv->lf,F->lf_D_JtRn_ZtSn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_D_JtRn_AtTn_vv - lf_D_JtRn_AtTn_uu %0.16f\n",dra_diff(F->lf_D_JtRn_AtTn_vv->lf,F->lf_D_JtRn_AtTn_uu->lf,length,1));
	/* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
}

void dcg_lf_D_AtTn_ZtSn_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  struct dcg_ajdk *D=NULL;struct dcg_single **E_=NULL; struct dcg_double **F_=NULL;
  dcg_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); dcg_init_QX(D);
  dcg_An_ajdk(D);
  dcg_lf_ZtSn(D);
  dcg_lf_D_AtTn_ZtSn_vv(D);
  dcg_lf_D_AtTn_ZtSn_uu(D);
  if (error_check){ dcg_lf_D_AtTn_ZtSn_error(verbose,D);}
  wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */


void dcg_lf_TAnZtS_ww(struct dcg_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_; struct dcg_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; struct dcg_single *E_nb1=NULL,*E_nb2=NULL; struct dcg_double *F=NULL;
  if (verbose){ printf(" %% [entering dcg_TAnZtS_ww]\n");}
  dcg_TAnZtS_ww_stage_0(D);
  dcg_TAnZtS_ww_stage_1(D);
  dcg_TAnZtS_ww_stage_2(D);
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A_bmr_j: ");
	bprintf(E_nb1->M_Jn->mr_j,D->bitj,1,E_nb1->A_nrows," %% J_bmr_j: ");
	bprintf(E_nb1->M_Kn->mc_j,D->bitj,1,D->A_ncols," %% K_bmc_j: ");
	bprintf(E_nb1->M_Tt->mr_j,D->bitj,1,D->T_ncols," %% T_bmc_j: ");
	/* %%%%%%%% */
	lfprintf(F->lf_TYnWtS_ww," %% lf_TYnWtS_ww: ");
	lfprintf(F->lf_TAnZtS_ww," %% lf_TAnZtS_ww: ");
	lfprintf(F->lf_TYnYtT_ww," %% lf_TYnYtT_ww: ");
	lfprintf(F->lf_TAnAtT_ww," %% lf_TAnAtT_ww: ");
	/* %%%%%%%% */
	lfprintf(F->lf_RUnWtS_ww," %% lf_RUnWtS_ww: ");
	lfprintf(F->lf_RUnYtT_ww," %% lf_RUnYtT_ww: ");
	lfprintf(F->lf_RJnZtS_ww," %% lf_RJnZtS_ww: ");
	lfprintf(F->lf_RJnAtT_ww," %% lf_RJnAtT_ww: ");
	/* lfprintf(F->lf_TKnVtS_ww," %% lf_TKnVtS_ww: "); */
	/* lfprintf(F->lf_TKnKtT_ww," %% lf_TKnKtT_ww: "); */
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcg_TAnZtS_ww]\n");}
}

void dcg_lf_TAnZtS_uu(struct dcg_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_; struct dcg_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; struct dcg_single *E_nb1=NULL,*E_nb2=NULL; struct dcg_double *F=NULL;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcg_TAnZtS_uu]\n");}
  if (verbose){ printf(" %% calculating TAnZtS_uu.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx=nb1+nb2*nbins; E_nb1 = E_[nb1]; E_nb2 = E_[nb2]; F = F_[nbx];
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (E_nb1->A_rbother && D->Y_cbother && E_nb2->Z_rbother && D->Y_cbother){
	GLOBAL_pthread_tic();
	wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yn,E_nb1->M_Tt,E_nb2->M_Wn,E_nb2->M_St,D->Y_ajdk,&(F->lf_TYnWtS_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb1->M_Tt,E_nb2->M_Zn,E_nb2->M_St,D->A_ajdk,&(F->lf_TAnZtS_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->Y_cbother && E_nb2->A_rbother && D->Y_cbother){
	if (nb1==nb2){
	  GLOBAL_pthread_tic();
	  wrap_TAnAtT_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yn,E_nb1->M_Tt,E_nb2->M_Yn,E_nb2->M_Tt,D->Y_ajdk,&(F->lf_TYnYtT_uu));
	  GLOBAL_pthread_toc();
	  /* if (nb1==nb2){ } */}
	if (nb1!=nb2){
	  GLOBAL_pthread_tic();
	  wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Yn,E_nb1->M_Tt,E_nb2->M_Yn,E_nb2->M_Tt,D->Y_ajdk,&(F->lf_TYnYtT_uu));
	  GLOBAL_pthread_toc();
	  /* if (nb1!=nb2){ } */}
	/* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){
	if (nb1==nb2){
	  GLOBAL_pthread_tic();
	  wrap_TAnAtT_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb1->M_Tt,E_nb2->M_An,E_nb2->M_Tt,D->A_ajdk,&(F->lf_TAnAtT_uu));
	  GLOBAL_pthread_toc();
	  /* if (nb1==nb2){ } */}
	if (nb1!=nb2){
	  GLOBAL_pthread_tic();
	  wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb1->M_Tt,E_nb2->M_An,E_nb2->M_Tt,D->A_ajdk,&(F->lf_TAnAtT_uu));
	  GLOBAL_pthread_toc();
	  /* if (nb1!=nb2){ } */}
	/* if bother */}
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (E_nb1->A_rbother && D->Y_cbother && E_nb2->Z_rbother && D->Y_cbother){
	GLOBAL_pthread_tic();
	wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Un,E_nb1->M_Rt,E_nb2->M_Wn,E_nb2->M_St,D->Y_ajdk,&(F->lf_RUnWtS_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->Y_cbother && E_nb2->A_rbother && D->Y_cbother){
	GLOBAL_pthread_tic();
	wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Un,E_nb1->M_Rt,E_nb2->M_Yn,E_nb2->M_Tt,D->Y_ajdk,&(F->lf_RUnYtT_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Jn,E_nb1->M_Rt,E_nb2->M_Zn,E_nb2->M_St,D->A_ajdk,&(F->lf_RJnZtS_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Jn,E_nb1->M_Rt,E_nb2->M_An,E_nb2->M_Tt,D->A_ajdk,&(F->lf_RJnAtT_uu));
	GLOBAL_pthread_toc(); /* if bother */}
      /* if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){ */
      /* 	GLOBAL_pthread_tic(); */
      /* 	wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Kn,E_nb1->M_Tt,E_nb2->M_Vn,E_nb2->M_St,D->A_ajdk,&(F->lf_TKnVtS_uu)); */
      /* 	GLOBAL_pthread_toc(); /\* if bother *\/} */
      /* if (E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother){ */
      /* 	if (nb1==nb2){ */
      /* 	  GLOBAL_pthread_tic(); */
      /* 	  wrap_TAnAtT_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Kn,E_nb1->M_Tt,E_nb2->M_Kn,E_nb2->M_Tt,D->A_ajdk,&(F->lf_TKnKtT_uu)); */
      /* 	  GLOBAL_pthread_toc(); */
      /* 	  /\* if (nb1==nb2){ } *\/} */
      /* 	if (nb1!=nb2){ */
      /* 	  GLOBAL_pthread_tic(); */
      /* 	  wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_Kn,E_nb1->M_Tt,E_nb2->M_Kn,E_nb2->M_Tt,D->A_ajdk,&(F->lf_TKnKtT_uu)); */
      /* 	  GLOBAL_pthread_toc(); */
      /* 	  /\* if (nb1!=nb2){ } *\/} */
      /* 	/\* if bother *\/} */
      /*  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% TAnZtS_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcg_TAnZtS_uu]\n");}
}

void dcg_lf_TAnZtS_error(int verbose,struct dcg_ajdk *D)
{
  int nbins=D->nbins; struct dcg_single **E_=D->E_; struct dcg_double **F_=D->F_;
  int nbx=0,nb1=0,nb2=0,length=0; struct dcg_single *E_nb1=NULL,*E_nb2=NULL; struct dcg_double *F=NULL;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      length = E_nb1->A_nrows * D->T_ncols;
      printf(" %% nb1 %d nb2 %d length %d*%d=%d\n",nb1,nb2,E_nb1->A_nrows,D->T_ncols,length);
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (D->Y_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_TYnWtS_ww - lf_TYnWtS_uu %0.16f\n",dra_diff(F->lf_TYnWtS_ww->lf,F->lf_TYnWtS_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_TAnZtS_ww - lf_TAnZtS_uu %0.16f\n",dra_diff(F->lf_TAnZtS_ww->lf,F->lf_TAnZtS_uu->lf,length,1));
	/* if bother */}
      if (D->Y_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_TYnYtT_ww - lf_TYnYtT_uu %0.16f\n",dra_diff(F->lf_TYnYtT_ww->lf,F->lf_TYnYtT_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_TAnAtT_ww - lf_TAnAtT_uu %0.16f\n",dra_diff(F->lf_TAnAtT_ww->lf,F->lf_TAnAtT_uu->lf,length,1));
	/* if bother */}
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (D->Y_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_RUnWtS_ww - lf_RUnWtS_uu %0.16f\n",dra_diff(F->lf_RUnWtS_ww->lf,F->lf_RUnWtS_uu->lf,length,1));
	/* if bother */}
      if (D->Y_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_RUnYtT_ww - lf_RUnYtT_uu %0.16f\n",dra_diff(F->lf_RUnYtT_ww->lf,F->lf_RUnYtT_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_RJnZtS_ww - lf_RJnZtS_uu %0.16f\n",dra_diff(F->lf_RJnZtS_ww->lf,F->lf_RJnZtS_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_RJnAtT_ww - lf_RJnAtT_uu %0.16f\n",dra_diff(F->lf_RJnAtT_ww->lf,F->lf_RJnAtT_uu->lf,length,1));
	/* if bother */}
      /* if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){  */
      /* 	printf(" %% lf_TKnVtS_ww - lf_TKnVtS_uu %0.16f\n",dra_diff(F->lf_TKnVtS_ww->lf,F->lf_TKnVtS_uu->lf,length,1)); */
      /* 	/\* if bother *\/} */
      /* if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){  */
      /* 	printf(" %% lf_TKnKtT_ww - lf_TKnKtT_uu %0.16f\n",dra_diff(F->lf_TKnKtT_ww->lf,F->lf_TKnKtT_uu->lf,length,1)); */
      /* 	/\* if bother *\/} */
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
}

void dcg_lf_TAnZtS_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  struct dcg_ajdk *D=NULL;struct dcg_single **E_=NULL; struct dcg_double **F_=NULL;
  dcg_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); dcg_init_QX(D);
  dcg_An_ajdk(D);
  dcg_lf_ZtSn(D);
  dcg_lf_TAnZtS_ww(D);
  dcg_lf_TAnZtS_uu(D);
  if (error_check){ dcg_lf_TAnZtS_error(verbose,D);}
  wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void *get_dcg_halfloop(void *vp)
{
  /* This function converts from SPACING_a to SPACING_b along the T_ncols dimension.
     Warning! Later on we expect D->T_bmc_b to include only contiguous bits.
     This ensures that D->T_bmc_j will serve as a mask for QR_TAnAtT etc. after switching to SPACING_b. */
  int verbose=0;
  int ip=0,nc=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct dcg_double *F = (struct dcg_double *)(vpra[ip++]); struct dcg_ajdk *D = F->D; struct dcg_single *E_nb1 = F->E_nb1, *E_nb2 = F->E_nb2;
  /* %%%%%%%% */
  struct M_handle *M_An_nb1 = E_nb1->M_An;
  struct M_handle *M_At_nb1 = E_nb1->M_At;
  struct M_handle *M_Yn_nb1 = E_nb1->M_Yn;
  struct M_handle *M_Yt_nb1 = E_nb1->M_Yt;
  struct M_handle *M_Tt_nb1 = E_nb1->M_Tt;
  struct M_handle *M_At_nb2 = E_nb2->M_At;
  struct M_handle *M_Tt_nb2 = E_nb2->M_Tt;
  /* %%%%%%%% */
  struct M_handle *M_Jn_nb1 = E_nb1->M_Jn;
  struct M_handle *M_Jt_nb1 = E_nb1->M_Jt;
  struct M_handle *M_Kn_nb1 = E_nb1->M_Kn;
  struct M_handle *M_Kt_nb1 = E_nb1->M_Kt;
  struct M_handle *M_Un_nb1 = E_nb1->M_Un;
  struct M_handle *M_Ut_nb1 = E_nb1->M_Ut;
  struct M_handle *M_Rt_nb1 = E_nb1->M_Rt;
  /* %%%%%%%% */
  int A_ncols = M_At_nb1->nrows, Y_ncols = M_Yt_nb1->nrows, T_ncols = M_Tt_nb1->nrows;
  int A_nrows_nb1 = M_At_nb1->ncols, A_nrows_nb2 = M_At_nb2->ncols;
  /* %%%%%%%% */
  struct L_handle *lf_TYnWtS = F->lf_TYnWtS_ww;
  struct L_handle *lf_TAnZtS = F->lf_TAnZtS_ww;
  struct L_handle *lf_TYnYtT = F->lf_TYnYtT_ww;
  struct L_handle *lf_TAnAtT = F->lf_TAnAtT_ww;
  /* %%%%%%%% */
  struct L_handle *lf_RUnWtS = F->lf_RUnWtS_ww;
  struct L_handle *lf_RUnYtT = F->lf_RUnYtT_ww;
  struct L_handle *lf_RJnZtS = F->lf_RJnZtS_ww;
  struct L_handle *lf_RJnAtT = F->lf_RJnAtT_ww;
  /* struct L_handle *lf_TKnVtS = F->lf_TKnVtS_ww; */
  /* struct L_handle *lf_TKnKtT = F->lf_TKnKtT_ww; */
  /* %%%%%%%% */
  struct L_handle *lf_D_AtTn_ZtSn = F->lf_D_AtTn_ZtSn_vv;
  struct L_handle *lf_D_AtTn_AtTn = F->lf_D_AtTn_AtTn_vv;
  struct L_handle *lf_D_YtTn_WtSn = F->lf_D_YtTn_WtSn_vv;
  struct L_handle *lf_D_YtTn_YtTn = F->lf_D_YtTn_YtTn_vv;
  /* %%%%%%%% */
  struct L_handle *lf_D_JtRn_ZtSn = F->lf_D_JtRn_ZtSn_vv;
  struct L_handle *lf_D_JtRn_AtTn = F->lf_D_JtRn_AtTn_vv;
  struct L_handle *lf_D_KtTn_VtSn = F->lf_D_KtTn_VtSn_vv;
  struct L_handle *lf_D_KtTn_KtTn = F->lf_D_KtTn_KtTn_vv;
  struct L_handle *lf_D_UtRn_WtSn = F->lf_D_UtRn_WtSn_vv;
  struct L_handle *lf_D_UtRn_YtTn = F->lf_D_UtRn_YtTn_vv;
  /* %%%%%%%% */
  int QR_TYnWtS_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QR_TAnZtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QR_TYnYtT_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QR_TAnAtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  int QC_TYnWtS_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QC_TAnZtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QC_TYnYtT_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QC_TAnAtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  /* %%%%%%%% */
  int QR_RUnWtS_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QR_RUnYtT_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QR_RJnZtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QR_RJnAtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  /* int QR_TKnVtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother; */
  /* int QR_TKnKtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother; */
  int QC_RUnWtS_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->Z_rbother && D->Y_cbother;
  int QC_RUnYtT_bother = E_nb1->A_rbother && D->Y_cbother && E_nb2->A_rbother && D->Y_cbother;
  int QC_RJnZtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QC_RJnAtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  int QC_TKnVtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QC_TKnKtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  /* %%%%%%%% */
  struct L_handle *lf_tmp=NULL;
  lf_tmp = F->QR_TYnWtS; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_TAnZtS; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_TYnYtT; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_TAnAtT; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_TYnWtS; lf_tmp->row_stride = Y_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_TAnZtS; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_TYnYtT; lf_tmp->row_stride = Y_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_TAnAtT; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  /* %%%%%%%% */
  lf_tmp = F->QR_RUnWtS; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_RUnYtT; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_RJnZtS; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_RJnAtT; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  /* lf_tmp = F->QR_TKnVtS; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp); */
  /* lf_tmp = F->QR_TKnKtT; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp); */
  lf_tmp = F->QC_RUnWtS; lf_tmp->row_stride = Y_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_RUnYtT; lf_tmp->row_stride = Y_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_RJnZtS; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_RJnAtT; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_TKnVtS; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_TKnKtT; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  /* %%%%%%%% */
  struct L_handle *QR_TYnWtS = F->QR_TYnWtS;
  struct L_handle *QR_TAnZtS = F->QR_TAnZtS;
  struct L_handle *QR_TYnYtT = F->QR_TYnYtT;
  struct L_handle *QR_TAnAtT = F->QR_TAnAtT;
  struct L_handle *QC_TYnWtS = F->QC_TYnWtS;
  struct L_handle *QC_TAnZtS = F->QC_TAnZtS;
  struct L_handle *QC_TYnYtT = F->QC_TYnYtT;
  struct L_handle *QC_TAnAtT = F->QC_TAnAtT;
  /* %%%%%%%% */
  struct L_handle *QR_RUnWtS = F->QR_RUnWtS;
  struct L_handle *QR_RUnYtT = F->QR_RUnYtT;
  struct L_handle *QR_RJnZtS = F->QR_RJnZtS;
  struct L_handle *QR_RJnAtT = F->QR_RJnAtT;
  /* struct L_handle *QR_TKnVtS = F->QR_TKnVtS; */
  /* struct L_handle *QR_TKnKtT = F->QR_TKnKtT; */
  struct L_handle *QC_RUnWtS = F->QC_RUnWtS;
  struct L_handle *QC_RUnYtT = F->QC_RUnYtT;
  struct L_handle *QC_RJnZtS = F->QC_RJnZtS;
  struct L_handle *QC_RJnAtT = F->QC_RJnAtT;
  struct L_handle *QC_TKnVtS = F->QC_TKnVtS;
  struct L_handle *QC_TKnKtT = F->QC_TKnKtT;
  /* %%%%%%%% */
  int ns_j=0,ns_b=0,ns_a=0;
  int na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0;
  int ny_j=0,ny_b=0,ny_a=0,my_j=0,my_b=0,my_a=0;
  /* %%%%%%%% */
  if (QR_TYnWtS_bother || QR_TYnYtT_bother){
    if (verbose>3){ printf(" %% QR_TYnWtS QR_TYnYtT\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      my_j=0; while (my_j<M_Yn_nb1->rpop_j){
	my_a = M_Yn_nb1->m_a_[my_j]; my_b = M_Yn_nb1->m_b_[my_j];
	if (QR_TYnWtS_bother){ L2_set(QR_TYnWtS , my_j,my_b,my_a , ns_j,ns_b,ns_a , (*L2_get(lf_TYnWtS , my_j,my_b,my_a , ns_j,ns_b,ns_a)));}
	if (QR_TYnYtT_bother){ L2_set(QR_TYnYtT , my_j,my_b,my_a , ns_j,ns_b,ns_a , (*L2_get(lf_TYnYtT , my_j,my_b,my_a , ns_j,ns_b,ns_a)));}
	my_j++; /* while (my_j<M_Yn_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QR_TYnWtS," %% QR_TYnWtS: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QR_TYnYtT," %% QR_TYnYtT: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if bother */}
  if (QR_TAnZtS_bother || QR_TAnAtT_bother){
    if (verbose>3){ printf(" %% QR_TAnZtS QR_TAnAtT\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      ma_j=0; while (ma_j<M_An_nb1->rpop_j){
	ma_a = M_An_nb1->m_a_[ma_j]; ma_b = M_An_nb1->m_b_[ma_j];
	if (QR_TAnZtS_bother){ L2_set(QR_TAnZtS , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_TAnZtS , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	if (QR_TAnAtT_bother){ L2_set(QR_TAnAtT , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_TAnAtT , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	ma_j++; /* while (ma_j<M_An_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QR_TAnZtS," %% QR_TAnZtS: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QR_TAnAtT," %% QR_TAnAtT: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if bother */}
  if (QC_TYnWtS_bother || QC_TYnYtT_bother){
    if (verbose>3){ printf(" %% QC_TYnWtS, QC_TYnYtT\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      ny_j=0; while (ny_j<M_Yt_nb1->rpop_j){
	ny_a = M_Yt_nb1->m_a_[ny_j]; ny_b = M_Yt_nb1->m_b_[ny_j];
	if (QC_TYnWtS_bother){ L2_set(QC_TYnWtS , ny_j,ny_b,ny_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_YtTn_WtSn , ny_j,ny_b,ny_a , ns_j,ns_b,ns_a)));}
	if (QC_TYnYtT_bother){ L2_set(QC_TYnYtT , ny_j,ny_b,ny_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_YtTn_YtTn , ny_j,ny_b,ny_a , ns_j,ns_b,ns_a)));}
	ny_j++; /* while (ny_j<M_Yt_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QC_TYnWtS," %% QC_TYnWtS: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QC_TYnYtT," %% QC_TYnYtT: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if (QC_TYnWtS_bother || QC_TAnAtT_bother){ } */}
  if (QC_TAnZtS_bother || QC_TAnAtT_bother){
    if (verbose>3){ printf(" %% QC_TAnZtS, QC_TAnAtT\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      na_j=0; while (na_j<M_At_nb1->rpop_j){
	na_a = M_At_nb1->m_a_[na_j]; na_b = M_At_nb1->m_b_[na_j];
	if (QC_TAnZtS_bother){ L2_set(QC_TAnZtS , na_j,na_b,na_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_AtTn_ZtSn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	if (QC_TAnAtT_bother){ L2_set(QC_TAnAtT , na_j,na_b,na_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_AtTn_AtTn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	na_j++; /* while (na_j<M_At_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QC_TAnZtS," %% QC_TAnZtS: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QC_TAnAtT," %% QC_TAnAtT: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if (QC_TAnZtS_bother || QC_TAnAtT_bother){ } */}
  /* %%%%%%%% */

  if (QR_RUnWtS_bother || QR_RUnYtT_bother){
    if (verbose>3){ printf(" %% QR_RUnWtS QR_RUnYtT\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      my_j=0; while (my_j<M_Un_nb1->rpop_j){
	my_a = M_Un_nb1->m_a_[my_j]; my_b = M_Un_nb1->m_b_[my_j];
	if (QR_RUnWtS_bother){ L2_set(QR_RUnWtS , my_j,my_b,my_a , ns_j,ns_b,ns_a , (*L2_get(lf_RUnWtS , my_j,my_b,my_a , ns_j,ns_b,ns_a)));}
	if (QR_RUnYtT_bother){ L2_set(QR_RUnYtT , my_j,my_b,my_a , ns_j,ns_b,ns_a , (*L2_get(lf_RUnYtT , my_j,my_b,my_a , ns_j,ns_b,ns_a)));}
	my_j++; /* while (my_j<M_Un_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QR_RUnWtS," %% QR_RUnWtS: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QR_RUnYtT," %% QR_RUnYtT: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if bother */}
  if (QR_RJnZtS_bother || QR_RJnAtT_bother){
    if (verbose>3){ printf(" %% QR_RJnZtS QR_RJnAtT\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      ma_j=0; while (ma_j<M_Jn_nb1->rpop_j){
	ma_a = M_Jn_nb1->m_a_[ma_j]; ma_b = M_Jn_nb1->m_b_[ma_j];
	if (QR_RJnZtS_bother){ L2_set(QR_RJnZtS , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_RJnZtS , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	if (QR_RJnAtT_bother){ L2_set(QR_RJnAtT , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_RJnAtT , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));}
	ma_j++; /* while (ma_j<M_Jn_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QR_RJnZtS," %% QR_RJnZtS: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QR_RJnAtT," %% QR_RJnAtT: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if bother */}
  /* if (QR_TKnVtS_bother || QR_TKnKtT_bother){ */
  /*   if (verbose>3){ printf(" %% QR_TKnVtS QR_TKnKtT\n");} if (verbose>3){ printf(" %% /\******************************************************************\/\n");} */
  /*   ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){ */
  /*     ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j]; */
  /*     ma_j=0; while (ma_j<M_Kn_nb1->rpop_j){ */
  /* 	ma_a = M_Kn_nb1->m_a_[ma_j]; ma_b = M_Kn_nb1->m_b_[ma_j]; */
  /* 	if (QR_TKnVtS_bother){ L2_set(QR_TKnVtS , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_TKnVtS , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));} */
  /* 	if (QR_TKnKtT_bother){ L2_set(QR_TKnKtT , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a , (*L2_get(lf_TKnKtT , ma_j,ma_b,ma_a , ns_j,ns_b,ns_a)));} */
  /* 	ma_j++; /\* while (ma_j<M_Kn_nb1->rpop_j){ } *\/} */
  /*     ns_j++; /\* while (ns_j<M_Tt_nb2->rpop_j){ } *\/} */
  /*   if (verbose>3){ lfprintf(QR_TKnVtS," %% QR_TKnVtS: ");} if (verbose>3){ printf(" %% /\******************************************************************\/\n");} */
  /*   if (verbose>3){ lfprintf(QR_TKnKtT," %% QR_TKnKtT: ");} if (verbose>3){ printf(" %% /\******************************************************************\/\n");} */
  /*   /\* if bother *\/} */
  if (QC_RUnWtS_bother || QC_RUnYtT_bother){
    if (verbose>3){ printf(" %% QC_RUnWtS, QC_RUnYtT\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      ny_j=0; while (ny_j<M_Ut_nb1->rpop_j){
	ny_a = M_Ut_nb1->m_a_[ny_j]; ny_b = M_Ut_nb1->m_b_[ny_j];
	if (QC_RUnWtS_bother){ L2_set(QC_RUnWtS , ny_j,ny_b,ny_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_UtRn_WtSn , ny_j,ny_b,ny_a , ns_j,ns_b,ns_a)));}
	if (QC_RUnYtT_bother){ L2_set(QC_RUnYtT , ny_j,ny_b,ny_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_UtRn_YtTn , ny_j,ny_b,ny_a , ns_j,ns_b,ns_a)));}
	ny_j++; /* while (ny_j<M_Ut_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QC_RUnWtS," %% QC_RUnWtS: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QC_RUnYtT," %% QC_RUnYtT: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if (QC_RUnWtS_bother || QC_TAnAtT_bother){ } */}
  if (QC_RJnZtS_bother || QC_RJnAtT_bother){
    if (verbose>3){ printf(" %% QC_RJnZtS, QC_RJnAtT\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      na_j=0; while (na_j<M_Jt_nb1->rpop_j){
	na_a = M_Jt_nb1->m_a_[na_j]; na_b = M_Jt_nb1->m_b_[na_j];
	if (QC_RJnZtS_bother){ L2_set(QC_RJnZtS , na_j,na_b,na_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_JtRn_ZtSn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	if (QC_RJnAtT_bother){ L2_set(QC_RJnAtT , na_j,na_b,na_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_JtRn_AtTn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	na_j++; /* while (na_j<M_Jt_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QC_RJnZtS," %% QC_RJnZtS: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QC_RJnAtT," %% QC_RJnAtT: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if (QC_RJnZtS_bother || QC_RJnAtT_bother){ } */}
  if (QC_TKnVtS_bother || QC_TKnKtT_bother){
    if (verbose>3){ printf(" %% QC_TKnVtS, QC_TKnKtT\n");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    ns_j=0; while (ns_j<M_Tt_nb2->rpop_j){
      ns_a = M_Tt_nb2->m_a_[ns_j]; ns_b = M_Tt_nb2->m_b_[ns_j];
      na_j=0; while (na_j<M_Kt_nb1->rpop_j){
	na_a = M_Kt_nb1->m_a_[na_j]; na_b = M_Kt_nb1->m_b_[na_j];
	if (QC_TKnVtS_bother){ L2_set(QC_TKnVtS , na_j,na_b,na_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_KtTn_VtSn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	if (QC_TKnKtT_bother){ L2_set(QC_TKnKtT , na_j,na_b,na_a , ns_j,ns_b,ns_a , (*L2_get(lf_D_KtTn_KtTn , na_j,na_b,na_a , ns_j,ns_b,ns_a)));}
	na_j++; /* while (na_j<M_Kt_nb1->rpop_j){ } */}
      ns_j++; /* while (ns_j<M_Tt_nb2->rpop_j){ } */}
    if (verbose>3){ lfprintf(QC_TKnVtS," %% QC_TKnVtS: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    if (verbose>3){ lfprintf(QC_TKnKtT," %% QC_TKnKtT: ");} if (verbose>3){ printf(" %% /******************************************************************/\n");}
    /* if (QC_TKnVtS_bother || QC_TKnKtT_bother){ } */}
  /* %%%%%%%% */
  if (verbose>1){ printf(" %% [finished get_dcg_halfloop] tidx %d\n",tidx);}
  return NULL;
}

void wrap_dcg_halfloop(int *tidx,void **vpra,pthread_t *thread_in,struct dcg_double *F)
{
  /* This function calls get_dcg_halfloop ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 2)
   */
  int verbose=0;
  int ip=0;
  if (verbose){ printf(" %% [entering wrap_dcg_halfloop] tidx %d\n",*tidx);}
  ip=0;
  vpra[ip++] = tidx; vpra[ip++] = F;
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_dcg_halfloop,vpra)){ printf("Warning! cannot create thread %d in wrap_dcg_halfloop\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_dcg_halfloop(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_dcg_halfloop] tidx %d\n",*tidx);}
}

void dcg_wrap_dcg_halfloop(struct dcg_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcg_double **F_ = D->F_; 
  int nbx=0,nb1=0,nb2=0; struct dcg_double *F=NULL;
  if (verbose){ printf(" %% [entering dcg_wrap_dcg_halfloop]\n");}
  if (verbose){ printf(" %% calculating halfloop.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx];
    GLOBAL_pthread_tic();
    wrap_dcg_halfloop(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),F);
    GLOBAL_pthread_toc();
    /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% halfloop: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcg_wrap_dcg_halfloop]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void dcg_lrup_mxdup(struct dcg_ajdk *D)
{
  /* Copies bmc_j_rtn into bmc_j, and bmr_j_rtn into bmr_j */
  /* Also copies umc_j_rtn into umc_j, and umr_j_rtn into umr_j */
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_; struct dcg_single *E=NULL;
  int nb=0;
  ura_dup(D->A_bmc_j,D->A_mc_length,D->A_bmc_j_rtn); ura_dup(D->A_umc_j,D->A_ncols,D->A_umc_j_rtn);
  for (nb=0;nb<nbins;nb++){ E = E_[nb]; ura_dup(E->A_bmr_j,E->A_mr_length,E->A_bmr_j_rtn); ura_dup(E->A_umr_j,E->A_nrows,E->A_umr_j_rtn); /* for (nb=0;nb<nbins;nb++){ } */}
  ura_dup(D->K_bmc_j,D->A_mc_length,D->K_bmc_j_rtn); ura_dup(D->K_umc_j,D->A_ncols,D->K_umc_j_rtn);
  for (nb=0;nb<nbins;nb++){ E = E_[nb]; ura_dup(E->J_bmr_j,E->A_mr_length,E->J_bmr_j_rtn); ura_dup(E->J_umr_j,E->A_nrows,E->J_umr_j_rtn); /* for (nb=0;nb<nbins;nb++){ } */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
