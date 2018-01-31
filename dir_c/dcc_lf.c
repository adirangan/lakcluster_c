/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void dcc_An_ajdk(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb=0; struct dcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcc_An_ajdk]\n");}
  if (verbose){ printf(" %% calculating E_[nb]->lf_An_ajdk, E_[nb]->lf_Zn_ajdk.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_An,D->A_ajdk,&(E->lf_An_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Zn,D->A_ajdk,&(E->lf_Zn_ajdk));
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
      lfprintf(E->lf_Zn_ajdk," %% lf_Zn_ajdk->lf: ");
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcc_An_ajdk]\n");}
}

void dcc_lf_ZtSn(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb=0; struct dcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcc_lf_ZtSn]\n");}
  if (verbose){ printf(" %% calculating E_[nb]->lf_AtTn, E_[nb]->lf_ZtSn.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    if (E->A_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_At,E->M_Tt,NULL,NULL,NULL,&(E->lf_AtTn));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->A_cbother){ 
      GLOBAL_pthread_tic(); 
      wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Zt,E->M_St,NULL,NULL,NULL,&(E->lf_ZtSn));
      GLOBAL_pthread_toc(); /* if bother */} 
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% ZtSn_vv: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      printf(" %% nb %d\n",nb);
      bprintf(E->M_An->mc_j,D->bitj,1,D->A_ncols," %% A_bmc_j: ");
      lfprintf(E->lf_AtTn," %% lf_AtTn: ");
      lfprintf(E->lf_ZtSn," %% lf_ZtSn: ");
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcc_lf_ZtSn]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void dcc_lf_D_AtTn_ZtSn_vv(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; struct dcc_single *E_nb1=NULL,*E_nb2=NULL; struct dcc_double *F=NULL;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcc_lf_D_AtTn_ZtSn_vv]\n");}
  if (verbose){ printf(" %% calculating F_[nbx]->lf_D_AtTn_ZtSn, F_[nbx]->lf_D_AtTn_AtTn.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx=nb1+nb2*nbins; E_nb1 = E_[nb1]; E_nb2 = E_[nb2]; F = F_[nbx];
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
      /*  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% D_AtTn_ZtSn_vv: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_At->mr_j,D->bitj,1,D->A_ncols," %% A_bmc_j: ");
	bprintf(E_nb1->M_Tt->mr_j,D->bitj,1,D->T_ncols," %% T_bmc_j: ");
	lfprintf(F->lf_D_AtTn_ZtSn_vv," %% lf_D_AtTn_ZtSn_vv: ");
	lfprintf(F->lf_D_AtTn_AtTn_vv," %% lf_D_AtTn_AtTn_vv: ");
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcc_lf_D_AtTn_ZtSn_vv]\n");}
}

void dcc_lf_D_AtTn_ZtSn_uu(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct dcc_single *E_nb1=NULL,*E_nb2=NULL; struct dcc_double *F=NULL,*F_trn=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcc_lf_D_AtTn_ZtSn_uu]\n");}
  if (verbose){ printf(" %% calculating F_[nbx]->lf_D_AtTn_ZtSn_uu, F_[nbx]->lf_D_AtTn_AtTn_uu.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
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
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% D_AtTn_ZtSn_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_At->mr_j,D->bitj,1,D->A_ncols," %% A_bmc_j: ");
	bprintf(E_nb1->M_Tt->mr_j,D->bitj,1,D->T_ncols," %% T_bmc_j: ");
	lfprintf(F->lf_D_AtTn_ZtSn_uu," %% lf_D_AtTn_ZtSn_uu: ");
	lfprintf(F->lf_D_AtTn_AtTn_uu," %% lf_D_AtTn_AtTn_uu: ");
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcc_lf_D_AtTn_ZtSn_uu]\n");}
}


void dcc_lf_D_AtTn_ZtSn_error(int verbose,struct dcc_ajdk *D)
{
  int nbins=D->nbins; struct dcc_single **E_=D->E_; struct dcc_double **F_=D->F_;
  int nbx=0,nb1=0,nb2=0,length=0; struct dcc_single *E_nb1=NULL,*E_nb2=NULL; struct dcc_double *F=NULL;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      length = D->A_ncols * D->T_ncols;
      printf(" %% nb1 %d nb2 %d length %d*%d=%d\n",nb1,nb2,D->A_ncols,D->T_ncols,length);
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_D_AtTn_ZtSn_vv - lf_D_AtTn_ZtSn_uu %0.16f\n",dra_diff(F->lf_D_AtTn_ZtSn_vv->lf,F->lf_D_AtTn_ZtSn_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_D_AtTn_AtTn_vv - lf_D_AtTn_AtTn_uu %0.16f\n",dra_diff(F->lf_D_AtTn_AtTn_vv->lf,F->lf_D_AtTn_AtTn_uu->lf,length,1));
	/* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
}

void dcc_lf_D_AtTn_ZtSn_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  dcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); dcc_init_QX(D);
  dcc_An_ajdk(D);
  dcc_lf_ZtSn(D);
  dcc_lf_D_AtTn_ZtSn_vv(D);
  dcc_lf_D_AtTn_ZtSn_uu(D);
  if (error_check){ dcc_lf_D_AtTn_ZtSn_error(verbose,D);}
  wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */


void dcc_lf_TAnZtS_ww(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; struct dcc_single *E_nb1=NULL,*E_nb2=NULL; struct dcc_double *F=NULL;
  if (verbose){ printf(" %% [entering dcc_TAnZtS_ww]\n");}
  dcc_TAnZtS_ww_stage_0(D);
  dcc_TAnZtS_ww_stage_1(D);
  dcc_TAnZtS_ww_stage_2(D);
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A_bmr_j: ");
	bprintf(E_nb1->M_Tt->mr_j,D->bitj,1,D->T_ncols," %% T_bmc_j: ");
	lfprintf(F->lf_TAnZtS_ww," %% lf_TAnZtS_ww: ");
	lfprintf(F->lf_TAnAtT_ww," %% lf_TAnAtT_ww: ");
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished dcc_TAnZtS_ww]\n");}
}

void dcc_lf_TAnZtS_uu(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; struct dcc_single *E_nb1=NULL,*E_nb2=NULL; struct dcc_double *F=NULL;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering dcc_TAnZtS_uu]\n");}
  if (verbose){ printf(" %% calculating TAnZtS_uu.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx=nb1+nb2*nbins; E_nb1 = E_[nb1]; E_nb2 = E_[nb2]; F = F_[nbx];
      if (E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother){
	GLOBAL_pthread_tic();
	wrap_TAnZtS_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb1->M_Tt,E_nb2->M_Zn,E_nb2->M_St,D->A_ajdk,&(F->lf_TAnZtS_uu));
	GLOBAL_pthread_toc(); /* if bother */}
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
      /*  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% TAnZtS_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcc_TAnZtS_uu]\n");}
}


void dcc_lf_TAnZtS_error(int verbose,struct dcc_ajdk *D)
{
  int nbins=D->nbins; struct dcc_single **E_=D->E_; struct dcc_double **F_=D->F_;
  int nbx=0,nb1=0,nb2=0,length=0; struct dcc_single *E_nb1=NULL,*E_nb2=NULL; struct dcc_double *F=NULL;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      length = E_nb1->A_nrows * D->T_ncols;
      printf(" %% nb1 %d nb2 %d length %d*%d=%d\n",nb1,nb2,E_nb1->A_nrows,D->T_ncols,length);
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	printf(" %% lf_TAnZtS_ww - lf_TAnZtS_uu %0.16f\n",dra_diff(F->lf_TAnZtS_ww->lf,F->lf_TAnZtS_uu->lf,length,1));
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	printf(" %% lf_TAnAtT_ww - lf_TAnAtT_uu %0.16f\n",dra_diff(F->lf_TAnAtT_ww->lf,F->lf_TAnAtT_uu->lf,length,1));
	/* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
}

void dcc_lf_TAnZtS_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  dcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); dcc_init_QX(D);
  dcc_An_ajdk(D);
  dcc_lf_ZtSn(D);
  dcc_lf_TAnZtS_ww(D);
  dcc_lf_TAnZtS_uu(D);
  if (error_check){ dcc_lf_TAnZtS_error(verbose,D);}
  wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void *get_halfloop(void *vp)
{
  /* This function converts from SPACING_a to SPACING_b along the T_ncols dimension. */
  int verbose=0;
  int ip=0,nc=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct dcc_double *F = (struct dcc_double *)(vpra[ip++]); struct dcc_ajdk *D = F->D; struct dcc_single *E_nb1 = F->E_nb1, *E_nb2 = F->E_nb2;
  struct M_handle *M_An_nb1 = E_nb1->M_An;
  struct M_handle *M_At_nb1 = E_nb1->M_At;
  struct M_handle *M_Tt_nb1 = E_nb1->M_Tt;
  struct M_handle *M_At_nb2 = E_nb2->M_At;
  struct M_handle *M_Tt_nb2 = E_nb2->M_Tt;
  int A_ncols = M_At_nb1->nrows, T_ncols = M_Tt_nb1->nrows;
  int A_nrows_nb1 = M_At_nb1->ncols, A_nrows_nb2 = M_At_nb2->ncols;
  struct L_handle *lf_TAnZtS = F->lf_TAnZtS_ww;
  struct L_handle *lf_TAnAtT = F->lf_TAnAtT_ww;
  struct L_handle *lf_D_AtTn_ZtSn = F->lf_D_AtTn_ZtSn_vv;
  struct L_handle *lf_D_AtTn_AtTn = F->lf_D_AtTn_AtTn_vv;
  int QR_TAnZtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QR_TAnAtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;
  int QC_TAnZtS_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->Z_rbother && D->A_cbother;
  int QC_TAnAtT_bother = E_nb1->A_rbother && D->A_cbother && E_nb2->A_rbother && D->A_cbother;  
  struct L_handle *lf_tmp=NULL;
  lf_tmp = F->QR_TAnZtS; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QR_TAnAtT; lf_tmp->row_stride = A_nrows_nb1; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_TAnZtS; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  lf_tmp = F->QC_TAnAtT; lf_tmp->row_stride = A_ncols; lf_tmp->spacing_row = SPACING_a; lf_tmp->col_stride = T_ncols; lf_tmp->spacing_col = SPACING_b; L_zero(lf_tmp);
  struct L_handle *QR_TAnZtS = F->QR_TAnZtS;
  struct L_handle *QR_TAnAtT = F->QR_TAnAtT;
  struct L_handle *QC_TAnZtS = F->QC_TAnZtS;
  struct L_handle *QC_TAnAtT = F->QC_TAnAtT;
  int ns_j=0,ns_b=0,ns_a=0;
  int na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0;
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
  if (verbose>1){ printf(" %% [finished get_halfloop] tidx %d\n",tidx);}
  return NULL;
}

void wrap_halfloop(int *tidx,void **vpra,pthread_t *thread_in,struct dcc_double *F)
{
  /* This function calls get_halfloop ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 2)
   */
  int verbose=0;
  int ip=0;
  if (verbose){ printf(" %% [entering wrap_halfloop] tidx %d\n",*tidx);}
  ip=0;
  vpra[ip++] = tidx; vpra[ip++] = F;
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_halfloop,vpra)){ printf("Warning! cannot create thread %d in wrap_halfloop\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_halfloop(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_halfloop] tidx %d\n",*tidx);}
}

void dcc_halfloop(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_double **F_ = D->F_; 
  int nbx=0,nb1=0,nb2=0; struct dcc_double *F=NULL;
  if (verbose){ printf(" %% [entering dcc_halfloop]\n");}
  if (verbose){ printf(" %% calculating halfloop.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx];
    GLOBAL_pthread_tic();
    wrap_halfloop(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),F);
    GLOBAL_pthread_toc();
    /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc(); 
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% halfloop: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose){ printf(" %% [finished dcc_halfloop]\n");}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void dcc_lrup_mxdup(struct dcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_single *E=NULL;
  int nb=0;
  ura_dup(D->A_bmc_j,D->A_mc_length,D->A_bmc_j_rtn); ura_dup(D->A_umc_j,D->A_ncols,D->A_umc_j_rtn);
  for (nb=0;nb<nbins;nb++){ E = E_[nb]; ura_dup(E->A_bmr_j,E->A_mr_length,E->A_bmr_j_rtn); ura_dup(E->A_umr_j,E->A_nrows,E->A_umr_j_rtn); /* for (nb=0;nb<nbins;nb++){ } */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
