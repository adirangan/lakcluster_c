/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_An_ajdk(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering bcc_An_ajdk]\n");}
  if (verbose){ printf(" %% calculating E_[nb]->lf_An_ajdk, E_[nb]->lf_Zn_ajdk, E_[nb]->lf_Yn_ajdk, E_[nb]->lf_Wn_ajdk.\n");}
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
    if (E->A_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Yn,D->Y_ajdk,&(E->lf_Yn_ajdk));
      GLOBAL_pthread_toc(); /* if bother */}
    if (E->Z_rbother && D->Y_cbother){ 
      GLOBAL_pthread_tic();
      wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_Wn,D->Y_ajdk,&(E->lf_Wn_ajdk));
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
      lfprintf(E->lf_Yn_ajdk," %% lf_Yn_ajdk->lf: ");
      lfprintf(E->lf_Wn_ajdk," %% lf_Wn_ajdk->lf: ");
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished bcc_An_ajdk]\n");}
}

void bcc_lf_ZtSn(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (1 || strstr(D->QR_strategy,"ZtSWn") || strstr(D->QC_strategy,"ZtSWn")){ /* always recalculate for use in doublestudy and getscores */
    if (verbose){ printf(" %% [entering bcc_lf_ZtSn]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->lf_AtTn, E_[nb]->lf_ZtSn, E_[nb]->lf_YtTn, E_[nb]->lf_WtSn.\n");}
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
      if (E->A_rbother && D->Y_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Yt,E->M_Tt,NULL,NULL,NULL,&(E->lf_YtTn));
	GLOBAL_pthread_toc(); /* if bother */} 
      if (E->Z_rbother && D->Y_cbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_Wt,E->M_St,NULL,NULL,NULL,&(E->lf_WtSn));
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
	bprintf(E->M_Yn->mc_j,D->bitj,1,D->A_ncols," %% Y_bmc_j: ");
	lfprintf(E->lf_YtTn," %% lf_YtTn: ");
	lfprintf(E->lf_WtSn," %% lf_WtSn: ");
	/* for (nb=0;nb<nbins;nb++){ } */}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lf_ZtSn]\n");}
    /* if strategy */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_lf_ZtSWn(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  int n_spacing_A = SPACING_b;
  if (strstr(D->QR_strategy,"ZtSWn") || strstr(D->QC_strategy,"ZtSWn")){
    if (verbose){ printf(" %% [entering bcc_lf_ZtSWn]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->lf_AtTAn, E_[nb]->lf_AtTYn, E_[nb]->lf_ZtSZn, E_[nb]->lf_ZtSWn.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->Z_rbother && D->A_cbother && D->Y_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Zt,E->M_St,E->M_Wt,D->A_ajdk,D->Y_ajdk,E->lf_ZtSn,E->lf_WtSn,&(E->lf_ZtSWn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E->Z_rbother && D->A_cbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AtTAn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_Zt,E->M_St,E->M_Zt,D->A_ajdk,D->A_ajdk,E->lf_ZtSn,E->lf_ZtSn,&(E->lf_ZtSZn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E->A_rbother && D->A_cbother && D->Y_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_At,E->M_Tt,E->M_Yt,D->A_ajdk,D->Y_ajdk,E->lf_AtTn,E->lf_YtTn,&(E->lf_AtTYn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E->A_rbother && D->A_cbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AtTAn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_At,E->M_Tt,E->M_At,D->A_ajdk,D->A_ajdk,E->lf_AtTn,E->lf_AtTn,&(E->lf_AtTAn));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% ZtSWn_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	printf(" %% nb %d\n",nb);
	bprintf(E->M_An->mc_j,D->bitj,1,D->A_ncols," %% A_bmc_j: ");
	bprintf(E->M_Yn->mc_j,D->bitj,1,D->Y_ncols," %% Y_bmc_j: ");
	bprintf(E->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T_bmc_j: ");
	if (E->Z_rbother && D->A_cbother && D->Y_cbother){ lfprintf(E->lf_ZtSWn," %% lf_ZtSWn: "); /* if bother */} else{ printf(" %% no lf_ZtSWn\n");}
	if (E->Z_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_ZtSZn," %% lf_ZtSZn: "); /* if bother */} else{ printf(" %% no lf_ZtSZn\n");}
	if (E->A_rbother && D->A_cbother && D->Y_cbother){ lfprintf(E->lf_AtTYn," %% lf_AtTYn: "); /* if bother */} else{ printf(" %% no lf_AtTYn\n");}
	if (E->A_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_AtTAn," %% lf_AtTAn: "); /* if bother */} else{ printf(" %% no lf_AtTAn\n");}
	/* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lf_ZtSWn]\n");}
    /* if strategy */}
}

void bcc_M_ZtSWn_excerpt(int trm_flag,struct M_handle *M_St,struct M_handle *M_in1,struct M_handle *M_in2,struct L_handle *L_in,struct M_handle **M_out_)
{
  int ns_j=0,ns_b=0,ns_a=0;
  if (trm_flag==0){
    for (ns_j=0;ns_j<M_St->rpop_j;ns_j++){
      ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      GLOBAL_pthread_tic(); 
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_in1->mr_j,M_in1->mr_b,M_in2->mr_j,M_in2->mr_b,L_in,L3_lf_get(L_in,L_in->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_in1->nrows),&(M_in2->nrows),&(M_out_[ns_j]),&(GLOBAL_B_MLT),(addressable_0));
      GLOBAL_pthread_toc(); 
      /* for (ns_j=0;ns_j<M_St->rpop_j;ns_j++){ } */}
    /* if trm==0 */}
  if (trm_flag==1){
    for (ns_j=0;ns_j<M_St->rpop_j;ns_j++){
      ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      GLOBAL_pthread_tic();
      wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_in2->mr_j,M_in2->mr_b,M_in1->mr_j,M_in1->mr_b,L_in,L3_lf_get(L_in,L_in->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_in2->nrows),&(M_in1->nrows),&(M_out_[ns_j]),&(GLOBAL_B_MLT),(addressable_1));
      GLOBAL_pthread_toc(); 
      /* for (ns_j=0;ns_j<M_St->rpop_j;ns_j++){ } */}
    /* trm==1 */}
}

void bcc_M_ZtSWn(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0,ns_j=0,ns_b=0,ns_a=0; struct bcc_single *E=NULL;
  struct M_handle *M_Xt=NULL,*M_in1=NULL,*M_in2=NULL,**M_out_=NULL;
  struct L_handle *L_in=NULL;
  int trm_flag=0;
  if (strstr(D->QR_strategy,"ZtSWn")){
    if (verbose){ printf(" %% [entering bcc_M_ZtSWn]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->M_AtTAn, E_[nb]->M_AtTYn, E_[nb]->M_ZtSZn, E_[nb]->M_ZtSWn.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
      if (E->Z_rbother && D->A_cbother && D->Y_cbother){ 
	M_Xt = E->M_St; M_in1= E->M_Zt; M_in2 = E->M_Wt; M_out_ = E->M_ZtSWn_; L_in = E->lf_ZtSWn; trm_flag = E->M_ZtSWn_trm; if (verbose>0){ printf(" %% bcc_M_ZtSWn: trm_flag %d\n",trm_flag);}
	bcc_M_ZtSWn_excerpt(trm_flag,M_Xt,M_in1,M_in2,L_in,M_out_);
	/* if bother */}
      /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
      if (E->Z_rbother && D->A_cbother && D->A_cbother){ 
	M_Xt = E->M_St; M_in1= E->M_Zt; M_in2 = E->M_Zt; M_out_ = E->M_ZtSZn_; L_in = E->lf_ZtSZn; trm_flag = E->M_ZtSZn_trm; if (verbose>0){ printf(" %% bcc_M_ZtSZn: trm_flag %d\n",trm_flag);}
	bcc_M_ZtSWn_excerpt(trm_flag,M_Xt,M_in1,M_in2,L_in,M_out_);
	/* if bother */}
      /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
      if (E->A_rbother && D->A_cbother && D->Y_cbother){ 
	M_Xt = E->M_Tt; M_in1= E->M_At; M_in2 = E->M_Yt; M_out_ = E->M_AtTYn_; L_in = E->lf_AtTYn; trm_flag = E->M_AtTYn_trm; if (verbose>0){ printf(" %% bcc_M_AtTYn: trm_flag %d\n",trm_flag);}
	bcc_M_ZtSWn_excerpt(trm_flag,M_Xt,M_in1,M_in2,L_in,M_out_);
	/* if bother */}
      /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
      if (E->A_rbother && D->A_cbother && D->A_cbother){ 
	M_Xt = E->M_Tt; M_in1= E->M_At; M_in2 = E->M_At; M_out_ = E->M_AtTAn_; L_in = E->lf_AtTAn; trm_flag = E->M_AtTAn_trm; if (verbose>0){ printf(" %% bcc_M_AtTAn: trm_flag %d\n",trm_flag);}
	bcc_M_ZtSWn_excerpt(trm_flag,M_Xt,M_in1,M_in2,L_in,M_out_);
	/* if bother */}
      /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc ZtSWn: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose){ printf(" %% [finished bcc_M_ZtSWn]\n");}
    /* if strategy */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_lf_YnWt(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_b;
  if (strstr(D->QR_strategy,"YnWt") || strstr(D->QC_strategy,"YnWt")){
    if (verbose){ printf(" %% [entering bcc_lf_YnWt]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_YnWt, F_[nbx]->lf_YnYt, F_[nbx]->lf_AnZt, F_[nbx]->lf_AnAt.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_Yn,E_nb2->M_Wn,D->Y_ajdk,E_nb1->lf_Yn_ajdk,E_nb2->lf_Wn_ajdk,&(F->lf_YnWt));
	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_Zn,D->A_ajdk,E_nb1->lf_An_ajdk,E_nb2->lf_Zn_ajdk,&(F->lf_AnZt));
	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1==nb2){
	    GLOBAL_pthread_tic(); 
	    wrap_AnAt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_Yn,E_nb2->M_Yn,D->Y_ajdk,E_nb1->lf_Yn_ajdk,E_nb2->lf_Yn_ajdk,&(F->lf_YnYt));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1==nb2){ } */}
	  if (nb1<nb2){
	    GLOBAL_pthread_tic(); 
	    wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_Yn,E_nb2->M_Yn,D->Y_ajdk,E_nb1->lf_Yn_ajdk,E_nb2->lf_Yn_ajdk,&(F->lf_YnYt));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1!=nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1==nb2){
	    GLOBAL_pthread_tic(); 
	    wrap_AnAt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_An,D->A_ajdk,E_nb1->lf_An_ajdk,E_nb2->lf_An_ajdk,&(F->lf_AnAt));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1==nb2){ } */}
	  if (nb1<nb2){
	    GLOBAL_pthread_tic(); 
	    wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_An,D->A_ajdk,E_nb1->lf_An_ajdk,E_nb2->lf_An_ajdk,&(F->lf_AnAt));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1!=nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% YnWt_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_transpose(F->lf_YnYt,F_trn->lf_YnYt);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_transpose(F->lf_AnAt,F_trn->lf_AnAt);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	  bprintf(E_nb2->M_Zn->mr_j,D->bitj,1,E_nb2->Z_nrows," %% Z2_bmr_j: ");
	  lfprintf(F->lf_YnWt," %% lf_YnWt: ");
	  lfprintf(F->lf_AnZt," %% lf_AnZt: ");
	  bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	  bprintf(E_nb2->M_An->mr_j,D->bitj,1,E_nb2->A_nrows," %% A2_bmr_j: ");
	  lfprintf(F->lf_YnYt," %% lf_YnYt: ");
	  lfprintf(F->lf_AnAt," %% lf_AnAt: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lf_YnWt]\n");}
    /* if strategy */}
}

void bcc_M_YnWt_excerpt(int trm_flag,struct M_handle *M_in1,struct M_handle *M_in2,struct L_handle *L_in,struct M_handle **M_out_p)
{
  if (trm_flag==0){
    GLOBAL_pthread_tic(); 
    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_in1->mr_j,M_in1->mr_b,M_in2->mr_j,M_in2->mr_b,L_in,L_in->lf,&(M_in1->nrows),&(M_in2->nrows),&(*M_out_p),&(GLOBAL_B_MLT),(addressable_0));
    GLOBAL_pthread_toc(); 
    /* if trm==0 */}
  if (trm_flag==1){
    GLOBAL_pthread_tic();
    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_in2->mr_j,M_in2->mr_b,M_in1->mr_j,M_in1->mr_b,L_in,L_in->lf,&(M_in2->nrows),&(M_in1->nrows),&(*M_out_p),&(GLOBAL_B_MLT),(addressable_1));
    GLOBAL_pthread_toc(); 
    /* trm==1 */}
}

void bcc_M_YnWt(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  struct M_handle *M_in1=NULL,*M_in2=NULL,**M_out_p=NULL;
  struct L_handle *L_in=NULL;
  int trm_flag=0;
  if (strstr(D->QC_strategy,"YnWt")){
    if (verbose){ printf(" %% [entering bcc_M_YnWt]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->M_YnWt, F_[nbx]->M_YnYt, F_[nbx]->M_AnZt, F_[nbx]->M_AnAt.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  M_in1 = E_nb1->M_Yn; M_in2 = E_nb2->M_Wn; L_in = F->lf_YnWt; M_out_p = &(F->M_YnWt); trm_flag = F->M_YnWt_trm; if (verbose>0){ printf(" %% bcc_M_YnWt: trm_flag %d\n",trm_flag);}
	  bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p);
	  /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  M_in1 = E_nb1->M_An; M_in2 = E_nb2->M_Zn; L_in = F->lf_AnZt; M_out_p = &(F->M_AnZt); trm_flag = F->M_AnZt_trm; if (verbose>0){ printf(" %% bcc_M_AnZt: trm_flag %d\n",trm_flag);}
	  bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p);
	  /* if bother */}
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  M_in1 = E_nb1->M_Yn; M_in2 = E_nb2->M_Yn; L_in = F->lf_YnYt; M_out_p = &(F->M_YnYt); trm_flag = F->M_YnYt_trm; if (verbose>0){ printf(" %% bcc_M_YnYt: trm_flag %d\n",trm_flag);}
	  bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p);
	  /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  M_in1 = E_nb1->M_An; M_in2 = E_nb2->M_An; L_in = F->lf_AnAt; M_out_p = &(F->M_AnAt); trm_flag = F->M_AnAt_trm; if (verbose>0){ printf(" %% bcc_M_AnAt: trm_flag %d\n",trm_flag);}
	  bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p);
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc YnWt: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose){ printf(" %% [finished bcc_M_YnWt]\n");}
    /* if strategy */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_lf_AnZt_S_WnYt(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QR_strategy,"YnWt")){
    if (verbose){ printf(" %% [entering bcc_lf_AnZt_S_WnYt]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_AnZt_S_WnYt, etc.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
  	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_S_WnYt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_St,E_nb2->M_Zn,F->lf_AnZt,F->lf_YnWt,&(F->lf_AnZt_S_WnYt));
  	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
  	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_S_WnYt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_St,E_nb2->M_Zn,F->lf_AnZt,F->lf_AnZt,&(F->lf_AnZt_S_ZnAt));
  	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
  	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_S_WnYt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_Tt,E_nb2->M_An,F->lf_AnAt,F->lf_YnYt,&(F->lf_AnAt_T_YnYt));
  	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
  	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_S_WnYt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_Tt,E_nb2->M_An,F->lf_AnAt,F->lf_AnAt,&(F->lf_AnAt_T_AnAt));
  	  GLOBAL_pthread_toc(); /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% AnZt_S_WnYt_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	  bprintf(E_nb2->M_Sn->mc_j,D->bitj,1,D->T_ncols," %% S2_bmc_j: ");
	  lfprintf(F->lf_AnZt_S_WnYt," %% lf_AnZt_S_WnYt: ");
	  lfprintf(F->lf_AnZt_S_ZnAt," %% lf_AnZt_S_ZnAt: ");
	  bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	  bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	  lfprintf(F->lf_AnAt_T_YnYt," %% lf_AnAt_T_YnYt: ");
	  lfprintf(F->lf_AnAt_T_AnAt," %% lf_AnAt_T_AnAt: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lf_AnZt_S_WnYt]\n");}
    /* if strategy */}
}

void bcc_lf_An_ZtSWn_Yt_excerpt(int trm_flag,int n_spacing_A,struct M_handle *M_An,struct M_handle *M_St,struct M_handle **M_ZtSWn_,struct M_handle *M_Zt,struct M_handle *M_Wt,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle **L_out_p)
{
  struct M_handle **M_WtSZn_ = M_ZtSWn_;
  if (trm_flag==0){
    GLOBAL_pthread_tic();
    wrap_An_ZtSWn_Yt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,M_An,M_St,M_ZtSWn_,M_Wt,M_Yn,A_ajdk,Y_ajdk,&(*L_out_p));
    GLOBAL_pthread_toc(); /* if trm==0 */}
  if (trm_flag==1){
    GLOBAL_pthread_tic();
    wrap_An_ZtSWn_Yt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,M_Yn,M_St,M_WtSZn_,M_Zt,M_An,Y_ajdk,A_ajdk,&(*L_out_p));
    GLOBAL_pthread_toc(); /* if trm==1 */}
}

void bcc_lf_An_ZtSWn_Yt(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QR_strategy,"ZtSWn")){
    if (verbose){ printf(" %% [entering bcc_lf_An_ZtSWn_Yt]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_An_ZtSWn_Yt, etc.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  bcc_lf_An_ZtSWn_Yt_excerpt(E_nb2->M_ZtSWn_trm,n_spacing_A,E_nb1->M_An,E_nb2->M_St,E_nb2->M_ZtSWn_,E_nb2->M_Zt,E_nb2->M_Wt,E_nb1->M_Yn,D->A_ajdk,D->Y_ajdk,&(F->lf_An_ZtSWn_Yt)); /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  bcc_lf_An_ZtSWn_Yt_excerpt(E_nb2->M_ZtSZn_trm,n_spacing_A,E_nb1->M_An,E_nb2->M_St,E_nb2->M_ZtSZn_,E_nb2->M_Zt,E_nb2->M_Zt,E_nb1->M_An,D->A_ajdk,D->A_ajdk,&(F->lf_An_ZtSZn_At)); /* if bother */}
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  bcc_lf_An_ZtSWn_Yt_excerpt(E_nb2->M_ZtSWn_trm,n_spacing_A,E_nb1->M_An,E_nb2->M_Tt,E_nb2->M_AtTYn_,E_nb2->M_At,E_nb2->M_Yt,E_nb1->M_Yn,D->A_ajdk,D->Y_ajdk,&(F->lf_An_AtTYn_Yt)); /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  bcc_lf_An_ZtSWn_Yt_excerpt(E_nb2->M_ZtSZn_trm,n_spacing_A,E_nb1->M_An,E_nb2->M_Tt,E_nb2->M_AtTAn_,E_nb2->M_At,E_nb2->M_At,E_nb1->M_An,D->A_ajdk,D->A_ajdk,&(F->lf_An_AtTAn_At)); /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% An_ZtSWn_Yt_ww: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	  bprintf(E_nb2->M_Sn->mc_j,D->bitj,1,D->T_ncols," %% S2_bmc_j: ");
	  lfprintf(F->lf_An_ZtSWn_Yt," %% lf_An_ZtSWn_Yt: ");
	  lfprintf(F->lf_An_ZtSZn_At," %% lf_An_ZtSZn_At: ");
	  bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	  bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	  lfprintf(F->lf_An_AtTYn_Yt," %% lf_An_AtTYn_Yt: ");
	  lfprintf(F->lf_An_AtTAn_At," %% lf_An_AtTAn_At: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lf_An_ZtSWn_Yt]\n");}
    /* if strategy */}
}

void bcc_lf_AnZtSWnYt(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering bcc_lf_AnZtSWnYt]\n");}
  if (verbose){ printf(" %% calculating F_[nbx]->lf_AnZtSWnYt, etc.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_Zt,E_nb2->M_St,E_nb2->M_Wt,E_nb1->M_Yn,D->A_ajdk,D->Y_ajdk,&(F->lf_AnZtSWnYt));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_Zt,E_nb2->M_St,E_nb2->M_Zt,E_nb1->M_An,D->A_ajdk,D->A_ajdk,&(F->lf_AnZtSZnAt));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_At,E_nb2->M_Tt,E_nb2->M_Yt,E_nb1->M_Yn,D->A_ajdk,D->Y_ajdk,&(F->lf_AnAtTYnYt));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_An,E_nb2->M_At,E_nb2->M_Tt,E_nb2->M_At,E_nb1->M_An,D->A_ajdk,D->A_ajdk,&(F->lf_AnAtTAnAt));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% AnZtSWnYt_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	bprintf(E_nb2->M_Sn->mc_j,D->bitj,1,D->T_ncols," %% S2_bmc_j: ");
	lfprintf(F->lf_AnZtSWnYt," %% lf_AnZtSWnYt: ");
	lfprintf(F->lf_AnZtSZnAt," %% lf_AnZtSZnAt: ");
	bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	lfprintf(F->lf_AnAtTYnYt," %% lf_AnAtTYnYt: ");
	lfprintf(F->lf_AnAtTAnAt," %% lf_AnAtTAnAt: ");
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished bcc_lf_AnZtSWnYt]\n");}
}

void bcc_lf_AnZtSWnYt_error(int verbose,struct bcc_ajdk *D)
{
  int nbins=D->nbins; struct bcc_single **E_=D->E_; struct bcc_double **F_=D->F_;
  int nbx=0,nb1=0,nb2=0,length=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      length = E_nb1->A_nrows * D->T_ncols;
      printf(" %% nb1 %d nb2 %d length %d*%d=%d\n",nb1,nb2,E_nb1->A_nrows,D->T_ncols,length);
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	if (strstr(D->QR_strategy,"ZtSWn") && strstr(D->QR_strategy,"YnWt")){ printf(" %% lf_An_ZtSWn_Yt - lf_AnZt_S_WnYt %0.16f\n",dra_diff(F->lf_An_ZtSWn_Yt->lf,F->lf_AnZt_S_WnYt->lf,length,1));}
	if (strstr(D->QR_strategy,"ZtSWn")){ printf(" %% lf_An_ZtSWn_Yt - lf_AnZtSWnYt   %0.16f\n",dra_diff(F->lf_An_ZtSWn_Yt->lf,F->lf_AnZtSWnYt->lf,length,1));}
	if (strstr(D->QR_strategy,"YnWt")){  printf(" %% lf_AnZt_S_WnYt - lf_AnZtSWnYt   %0.16f\n",dra_diff(F->lf_AnZt_S_WnYt->lf,F->lf_AnZtSWnYt->lf,length,1));}
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	if (strstr(D->QR_strategy,"ZtSWn") && strstr(D->QR_strategy,"YnWt")){ printf(" %% lf_An_ZtSZn_At - lf_AnZt_S_ZnAt %0.16f\n",dra_diff(F->lf_An_ZtSZn_At->lf,F->lf_AnZt_S_ZnAt->lf,length,1));}
	if (strstr(D->QR_strategy,"ZtSWn")){ printf(" %% lf_An_ZtSZn_At - lf_AnZtSZnAt   %0.16f\n",dra_diff(F->lf_An_ZtSZn_At->lf,F->lf_AnZtSZnAt->lf,length,1));}
	if (strstr(D->QR_strategy,"YnWt")){  printf(" %% lf_AnZt_S_ZnAt - lf_AnZtSZnAt   %0.16f\n",dra_diff(F->lf_AnZt_S_ZnAt->lf,F->lf_AnZtSZnAt->lf,length,1));}
	/* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	if (strstr(D->QR_strategy,"ZtSWn") && strstr(D->QR_strategy,"YnWt")){ printf(" %% lf_An_AtTYn_Yt - lf_AnAt_T_YnYt %0.16f\n",dra_diff(F->lf_An_AtTYn_Yt->lf,F->lf_AnAt_T_YnYt->lf,length,1));}
	if (strstr(D->QR_strategy,"ZtSWn")){ printf(" %% lf_An_AtTYn_Yt - lf_AnAtTYnYt   %0.16f\n",dra_diff(F->lf_An_AtTYn_Yt->lf,F->lf_AnAtTYnYt->lf,length,1));}
	if (strstr(D->QR_strategy,"YnWt")){  printf(" %% lf_AnAt_T_YnYt - lf_AnAtTYnYt   %0.16f\n",dra_diff(F->lf_AnAt_T_YnYt->lf,F->lf_AnAtTYnYt->lf,length,1));}
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	if (strstr(D->QR_strategy,"ZtSWn") && strstr(D->QR_strategy,"YnWt")){ printf(" %% lf_An_AtTAn_At - lf_AnAt_T_AnAt %0.16f\n",dra_diff(F->lf_An_AtTAn_At->lf,F->lf_AnAt_T_AnAt->lf,length,1));}
	if (strstr(D->QR_strategy,"ZtSWn")){ printf(" %% lf_An_AtTAn_At - lf_AnAtTAnAt   %0.16f\n",dra_diff(F->lf_An_AtTAn_At->lf,F->lf_AnAtTAnAt->lf,length,1));}
	if (strstr(D->QR_strategy,"YnWt")){  printf(" %% lf_AnAt_T_AnAt - lf_AnAtTAnAt   %0.16f\n",dra_diff(F->lf_AnAt_T_AnAt->lf,F->lf_AnAtTAnAt->lf,length,1));}
	/* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
}

void bcc_lf_AnZtSWnYt_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  bcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D);
  bcc_An_ajdk(D);
  bcc_lf_ZtSn(D);
  bcc_lf_ZtSWn(D);
  bcc_M_ZtSWn(D);
  bcc_lf_YnWt(D);
  /* bcc_M_YnWt(D); */
  bcc_lf_AnZt_S_WnYt(D);
  bcc_lf_An_ZtSWn_Yt(D);
  bcc_lf_AnZtSWnYt(D);
  if (error_check){ bcc_lf_AnZtSWnYt_error(verbose,D);}
  wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_lf_AtTYn____WtSZn(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QC_strategy,"ZtSWn")){
    if (verbose){ printf(" %% [entering bcc_lf_AtTYn____WtSZn]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_AtTYn____WtSZn, etc.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
  	  GLOBAL_pthread_tic(); 
	  wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_Yt,E_nb2->M_Wt,E_nb2->M_St,E_nb2->M_Zt,D->Y_ajdk,E_nb1->lf_AtTYn,E_nb2->lf_ZtSWn,&(F->lf_AtTYn____WtSZn));
  	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
  	  GLOBAL_pthread_tic(); 
	  wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_At,E_nb2->M_Zt,E_nb2->M_St,E_nb2->M_Zt,D->A_ajdk,E_nb1->lf_AtTAn,E_nb2->lf_ZtSZn,&(F->lf_AtTAn____ZtSZn));
  	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1<=nb2){
	    GLOBAL_pthread_tic(); 
	    wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_Yt,E_nb2->M_Yt,E_nb2->M_Tt,E_nb2->M_At,D->Y_ajdk,E_nb1->lf_AtTYn,E_nb2->lf_AtTYn,&(F->lf_AtTYn____YtTAn));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1<=nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1<=nb2){
	    GLOBAL_pthread_tic(); 
	    wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_At,E_nb2->M_At,E_nb2->M_Tt,E_nb2->M_At,D->A_ajdk,E_nb1->lf_AtTAn,E_nb2->lf_AtTAn,&(F->lf_AtTAn____AtTAn));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1<=nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% AtTYn____WtSZn_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_duplicate(F->lf_AtTYn____YtTAn,F_trn->lf_AtTYn____YtTAn);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_duplicate(F->lf_AtTAn____AtTAn,F_trn->lf_AtTAn____AtTAn);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	  bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	  lfprintf(F->lf_AtTYn____WtSZn," %% lf_AtTYn____WtSZn: ");
	  lfprintf(F->lf_AtTAn____ZtSZn," %% lf_AtTAn____ZtSZn: ");
	  bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	  bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	  lfprintf(F->lf_AtTYn____YtTAn," %% lf_AtTYn____YtTAn: ");
	  lfprintf(F->lf_AtTAn____AtTAn," %% lf_AtTAn____AtTAn: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lf_AtTYn____WtSZn]\n");}
    /* if strategy */}
}

void bcc_lf_At_T_YnWt_S_Zn_excerpt(int trm_flag,int n_spacing_A,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yn,struct M_handle *M_Wn,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,struct M_handle *M_YnWt,struct L_handle **L_out_p)
{
  struct M_handle *M_WnYt = M_YnWt;
  if (trm_flag==0){
    GLOBAL_pthread_tic();
    wrap_At_T_YnWt_S_Zn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,M_At,M_Tt,M_Yn,M_Wn,M_St,M_Zt,A_ajdk,M_YnWt,&(*L_out_p));
    GLOBAL_pthread_toc(); /* if trm==0 */}
  if (trm_flag==1){
    GLOBAL_pthread_tic();
    wrap_At_T_YnWt_S_Zn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,M_Zt,M_St,M_Wn,M_Yn,M_Tt,M_At,A_ajdk,M_WnYt,&(*L_out_p));
    GLOBAL_pthread_toc(); /* if trm==1 */}
}

void bcc_lf_At_T_YnWt_S_Zn(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QC_strategy,"YnWt")){
    if (verbose){ printf(" %% [entering bcc_lf_At_T_YnWt_S_Zn]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_At_T_YnWt_S_Zn, F_[nbx]->lf_At_T_AnZt_S_Zn, F_[nbx]->lf_At_T_YnYt_T_An, F_[nbx]->lf_At_T_AnAt_T_An.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_YnWt_trm,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_Yn,E_nb2->M_Wn,E_nb2->M_St,E_nb2->M_Zt,D->A_ajdk,F->M_YnWt,&(F->lf_At_T_YnWt_S_Zn)); /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_AnZt_trm,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_An,E_nb2->M_Zn,E_nb2->M_St,E_nb2->M_Zt,D->A_ajdk,F->M_AnZt,&(F->lf_At_T_AnZt_S_Zn)); /* if bother */}
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1<=nb2){
	    bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_YnYt_trm,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_Yn,E_nb2->M_Yn,E_nb2->M_Tt,E_nb2->M_At,D->A_ajdk,F->M_YnYt,&(F->lf_At_T_YnYt_T_An)); 
	    /* if (nb1<=nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1<=nb2){
	    bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_AnAt_trm,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_An,E_nb2->M_An,E_nb2->M_Tt,E_nb2->M_At,D->A_ajdk,F->M_AnAt,&(F->lf_At_T_AnAt_T_An)); 
	    /* if (nb1<=nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% At_T_YnWt_S_Zn_ww: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_duplicate(F->lf_At_T_YnYt_T_An,F_trn->lf_At_T_YnYt_T_An);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_duplicate(F->lf_At_T_AnAt_T_An,F_trn->lf_At_T_AnAt_T_An);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	  bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	  lfprintf(F->lf_At_T_YnWt_S_Zn," %% lf_At_T_YnWt_S_Zn: ");
	  lfprintf(F->lf_At_T_AnZt_S_Zn," %% lf_At_T_AnZt_S_Zn: ");
	  bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	  bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	  lfprintf(F->lf_At_T_YnYt_T_An," %% lf_At_T_YnYt_T_An: ");
	  lfprintf(F->lf_At_T_AnAt_T_An," %% lf_At_T_AnAt_T_An: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lf_At_T_YnWt_S_Zn]\n");}
    /* if strategy */}
}

void bcc_lf_AtTYnWtSZn(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  int n_spacing_A = SPACING_a;
  if (verbose){ printf(" %% [entering bcc_lf_AtTYnWtSZn]\n");}
  if (verbose){ printf(" %% calculating F_[nbx]->lf_AtTYnWtSZn, etc.\n");}
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_Yn,E_nb2->M_Wn,E_nb2->M_St,E_nb2->M_Zt,D->A_ajdk,D->Y_ajdk,&(F->lf_AtTYnWtSZn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_An,E_nb2->M_Zn,E_nb2->M_St,E_nb2->M_Zt,D->A_ajdk,D->A_ajdk,&(F->lf_AtTAnZtSZn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_Yn,E_nb2->M_Yn,E_nb2->M_Tt,E_nb2->M_At,D->A_ajdk,D->Y_ajdk,&(F->lf_AtTYnYtTAn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_At,E_nb1->M_Tt,E_nb1->M_An,E_nb2->M_An,E_nb2->M_Tt,E_nb2->M_At,D->A_ajdk,D->A_ajdk,&(F->lf_AtTAnAtTAn));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  GLOBAL_pthread_tuc();
  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% AtTYnWtSZn_uu: ");
  GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if (verbose>2){  
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	lfprintf(F->lf_AtTYnWtSZn," %% lf_AtTYnWtSZn: ");
	lfprintf(F->lf_AtTAnZtSZn," %% lf_AtTAnZtSZn: ");
	bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	lfprintf(F->lf_AtTYnYtTAn," %% lf_AtTYnYtTAn: ");
	lfprintf(F->lf_AtTAnAtTAn," %% lf_AtTAnAtTAn: ");
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (verbose>2){ } */}
  if (verbose){ printf(" %% [finished bcc_lf_AtTYnWtSZn]\n");}
}

void bcc_lf_AtTYnWtSZn_error(int verbose,struct bcc_ajdk *D)
{
  int nbins=D->nbins; struct bcc_single **E_=D->E_; struct bcc_double **F_=D->F_;
  int nbx=0,nb1=0,nb2=0,length=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
      length = D->A_ncols * D->T_ncols;
      printf(" %% nb1 %d nb2 %d length %d*%d=%d\n",nb1,nb2,D->A_ncols,D->T_ncols,length);
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	if (strstr(D->QC_strategy,"YnWt") && strstr(D->QC_strategy,"ZtSWn")){ printf(" %% lf_At_T_YnWt_S_Zn - lf_AtTYn____WtSZn %0.16f\n",dra_diff(F->lf_At_T_YnWt_S_Zn->lf,F->lf_AtTYn____WtSZn->lf,length,1));}
	if (strstr(D->QC_strategy,"YnWt")){  printf(" %% lf_At_T_YnWt_S_Zn - lf_AtTYnWtSZn     %0.16f\n",dra_diff(F->lf_At_T_YnWt_S_Zn->lf,F->lf_AtTYnWtSZn->lf,length,1));}
	if (strstr(D->QC_strategy,"ZtSWn")){ printf(" %% lf_AtTYn____WtSZn - lf_AtTYnWtSZn     %0.16f\n",dra_diff(F->lf_AtTYn____WtSZn->lf,F->lf_AtTYnWtSZn->lf,length,1));}
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	if (strstr(D->QC_strategy,"YnWt") && strstr(D->QC_strategy,"ZtSWn")){ printf(" %% lf_At_T_AnZt_S_Zn - lf_AtTAn____ZtSZn %0.16f\n",dra_diff(F->lf_At_T_AnZt_S_Zn->lf,F->lf_AtTAn____ZtSZn->lf,length,1));}
	if (strstr(D->QC_strategy,"YnWt")){  printf(" %% lf_At_T_AnZt_S_Zn - lf_AtTAnZtSZn     %0.16f\n",dra_diff(F->lf_At_T_AnZt_S_Zn->lf,F->lf_AtTAnZtSZn->lf,length,1));}
	if (strstr(D->QC_strategy,"ZtSWn")){ printf(" %% lf_AtTAn____ZtSZn - lf_AtTAnZtSZn     %0.16f\n",dra_diff(F->lf_AtTAn____ZtSZn->lf,F->lf_AtTAnZtSZn->lf,length,1));}
	/* if bother */}
      if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	if (strstr(D->QC_strategy,"YnWt") && strstr(D->QC_strategy,"ZtSWn")){ printf(" %% lf_At_T_YnYt_T_An - lf_AtTYn____YtTAn %0.16f\n",dra_diff(F->lf_At_T_YnYt_T_An->lf,F->lf_AtTYn____YtTAn->lf,length,1));}
	if (strstr(D->QC_strategy,"YnWt")){  printf(" %% lf_At_T_YnYt_T_An - lf_AtTYnYtTAn     %0.16f\n",dra_diff(F->lf_At_T_YnYt_T_An->lf,F->lf_AtTYnYtTAn->lf,length,1));}
	if (strstr(D->QC_strategy,"ZtSWn")){ printf(" %% lf_AtTYn____YtTAn - lf_AtTYnYtTAn     %0.16f\n",dra_diff(F->lf_AtTYn____YtTAn->lf,F->lf_AtTYnYtTAn->lf,length,1));}
	/* if bother */}
      if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	if (strstr(D->QC_strategy,"YnWt") && strstr(D->QC_strategy,"ZtSWn")){ printf(" %% lf_At_T_AnAt_T_An - lf_AtTAn____AtTAn %0.16f\n",dra_diff(F->lf_At_T_AnAt_T_An->lf,F->lf_AtTAn____AtTAn->lf,length,1));}
	if (strstr(D->QC_strategy,"YnWt")){  printf(" %% lf_At_T_AnAt_T_An - lf_AtTAnAtTAn     %0.16f\n",dra_diff(F->lf_At_T_AnAt_T_An->lf,F->lf_AtTAnAtTAn->lf,length,1));}
	if (strstr(D->QC_strategy,"ZtSWn")){ printf(" %% lf_AtTAn____AtTAn - lf_AtTAnAtTAn     %0.16f\n",dra_diff(F->lf_AtTAn____AtTAn->lf,F->lf_AtTAnAtTAn->lf,length,1));}
	/* if bother */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}} 
}

void bcc_lf_AtTYnWtSZn_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  bcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D);
  bcc_An_ajdk(D);
  bcc_lf_ZtSn(D);
  bcc_lf_ZtSWn(D);
  /* bcc_M_ZtSWn(D); */
  bcc_lf_YnWt(D);
  bcc_M_YnWt(D);
  bcc_lf_AtTYn____WtSZn(D);
  bcc_lf_At_T_YnWt_S_Zn(D);
  bcc_lf_AtTYnWtSZn(D);
  if (error_check){ bcc_lf_AtTYnWtSZn_error(verbose,D);}
  wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

