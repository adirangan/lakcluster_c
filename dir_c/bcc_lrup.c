void bcc_lrup_QC_YnWt_stage_a0(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QC_strategy,"YnWt")){
    if (verbose){ printf(" %% [entering bcc_lrup_QC_YnWt_stage_a0]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->lf_kn_ajdk, E_[nb]->lf_zn_ajdk, E_[nb]->lf_hn_ajdk, E_[nb]->lf_vn_ajdk, E_[nb]->lf_an_ajdk.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->A_rbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_kn,D->A_ajdk,&(E->lf_kn_ajdk));
	GLOBAL_pthread_toc(); 
	GLOBAL_pthread_tic();
	wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_hn,D->A_ajdk,&(E->lf_hn_ajdk));
	GLOBAL_pthread_toc(); 
	GLOBAL_pthread_tic();
	wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_an,D->A_ajdk,&(E->lf_an_ajdk));
	GLOBAL_pthread_toc(); 
	/* if bother */}
      if (E->Z_rbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_zn,D->A_ajdk,&(E->lf_zn_ajdk));
	GLOBAL_pthread_toc(); 
	GLOBAL_pthread_tic();
	wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_vn,D->A_ajdk,&(E->lf_vn_ajdk));
	GLOBAL_pthread_toc(); 
	/* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc(); 
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% An_ajdk_v: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){ 
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	printf(" %% nb %d\n",nb);
	bprintf(E->A_bmr_j_rtn,D->bitj,1,E->A_nrows," %% A_bmr_j_rtn: ");
	bprintf(E->M_an->mr_j ,D->bitj,1,E->A_nrows," %% M_an->mr_j : ");
	bprintf(E->M_jn->mr_j ,D->bitj,1,E->A_nrows," %% M_jn->mr_j : ");
	bprintf(E->M_hn->mr_j ,D->bitj,1,E->A_nrows," %% M_hn->mr_j : ");
	bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols," %% A_bmr_j_rtn: ");
	bprintf(E->M_an->mc_j ,D->bitj,1,D->A_ncols," %% M_an->mc_j : ");
	bprintf(E->M_jn->mc_j ,D->bitj,1,D->A_ncols," %% M_jn->mc_j : ");
	bprintf(E->M_hn->mc_j ,D->bitj,1,D->A_ncols," %% M_hn->mc_j : ");
	lfprintf(E->lf_an_ajdk," %% lf_an_ajdk->lf: ");
	lfprintf(E->lf_jn_ajdk," %% lf_jn_ajdk->lf: ");
	lfprintf(E->lf_hn_ajdk," %% lf_hn_ajdk->lf: ");
	bprintf(E->M_zn->mc_j ,D->bitj,1,D->A_ncols," %% M_zn->mc_j : ");
	bprintf(E->M_vn->mc_j ,D->bitj,1,D->A_ncols," %% M_vn->mc_j : ");
	lfprintf(E->lf_zn_ajdk," %% lf_zn_ajdk->lf: ");
	lfprintf(E->lf_vn_ajdk," %% lf_vn_ajdk->lf: ");
	/* for (nb=0;nb<nbins;nb++){ } */}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_YnWt_stage_a0]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_YnWt_stage_a1(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_j;
  if (strstr(D->QC_strategy,"YnWt")){
    if (verbose){ printf(" %% [entering bcc_lrup_QC_YnWt_stage_a1]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_unwt, F_[nbx]->lf_knzt, F_[nbx]->lf_hnvt, F_[nbx]->lf_ynut, F_[nbx]->lf_unyt, F_[nbx]->lf_unut, F_[nbx]->lf_ankt, F_[nbx]->lf_knat, F_[nbx]->lf_knkt.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  if (verbose){ printf(" %% unwt %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_un->rpop_j,E_nb1->M_un->cpop_j,E_nb2->M_Wn->rpop_j,E_nb2->M_Wn->cpop_j,F->lf_unwt->length);}
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_un,E_nb2->M_Wn,D->Y_ajdk,E_nb1->lf_Yn_ajdk,E_nb2->lf_Wn_ajdk,&(F->lf_unwt));
	  GLOBAL_pthread_toc(); 
	  /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  if (verbose){ printf(" %% knzt %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_kn->rpop_j,E_nb1->M_kn->cpop_j,E_nb2->M_zn->rpop_j,E_nb2->M_zn->cpop_j,F->lf_knzt->length);}
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_kn,E_nb2->M_zn,D->A_ajdk,E_nb1->lf_kn_ajdk,E_nb2->lf_zn_ajdk,&(F->lf_knzt));
	  GLOBAL_pthread_toc(); 
	  if (verbose){ printf(" %% hnvt %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_hn->rpop_j,E_nb1->M_hn->cpop_j,E_nb2->M_vn->rpop_j,E_nb2->M_vn->cpop_j,F->lf_hnvt->length);}
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_hn,E_nb2->M_vn,D->A_ajdk,E_nb1->lf_hn_ajdk,E_nb2->lf_vn_ajdk,&(F->lf_hnvt));
	  GLOBAL_pthread_toc(); 
	  /* if bother */}
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (verbose){ printf(" %% ynut %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_yn->rpop_j,E_nb1->M_yn->cpop_j,E_nb2->M_un->rpop_j,E_nb2->M_un->cpop_j,F->lf_ynut->length);}
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_yn,E_nb2->M_un,D->Y_ajdk,E_nb1->lf_Yn_ajdk,E_nb2->lf_Yn_ajdk,&(F->lf_ynut));
	  GLOBAL_pthread_toc(); 
	  if (0){
	    if (verbose){ printf(" %% unyt %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_un->rpop_j,E_nb1->M_un->cpop_j,E_nb2->M_yn->rpop_j,E_nb2->M_yn->cpop_j,F->lf_unyt->length);}
	    GLOBAL_pthread_tic(); 
	    wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_un,E_nb2->M_yn,D->Y_ajdk,E_nb1->lf_Yn_ajdk,E_nb2->lf_Yn_ajdk,&(F->lf_unyt));
	    GLOBAL_pthread_toc(); 
	    /* skip */}
	  if (nb1==nb2){
	    if (verbose){ printf(" %% unut %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_un->rpop_j,E_nb1->M_un->cpop_j,E_nb2->M_un->rpop_j,E_nb2->M_un->cpop_j,F->lf_unut->length);}
	    GLOBAL_pthread_tic(); 
	    wrap_AnAt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_un,E_nb2->M_un,D->Y_ajdk,E_nb1->lf_Yn_ajdk,E_nb2->lf_Yn_ajdk,&(F->lf_unut));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1==nb2){ } */}
	  if (nb1<nb2){
	    if (verbose){ printf(" %% unut %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_un->rpop_j,E_nb1->M_un->cpop_j,E_nb2->M_un->rpop_j,E_nb2->M_un->cpop_j,F->lf_unut->length);}
	    GLOBAL_pthread_tic(); 
	    wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_un,E_nb2->M_un,D->Y_ajdk,E_nb1->lf_Yn_ajdk,E_nb2->lf_Yn_ajdk,&(F->lf_unut));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1<nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (verbose){ printf(" %% ankt %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_an->rpop_j,E_nb1->M_an->cpop_j,E_nb2->M_kn->rpop_j,E_nb2->M_kn->cpop_j,F->lf_ankt->length);}
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_an,E_nb2->M_kn,D->A_ajdk,E_nb1->lf_an_ajdk,E_nb2->lf_kn_ajdk,&(F->lf_ankt));
	  GLOBAL_pthread_toc(); 
	  if (0){
	    if (verbose){ printf(" %% knat %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_kn->rpop_j,E_nb1->M_kn->cpop_j,E_nb2->M_an->rpop_j,E_nb2->M_an->cpop_j,F->lf_knat->length);}
	    GLOBAL_pthread_tic(); 
	    wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_kn,E_nb2->M_an,D->A_ajdk,E_nb1->lf_kn_ajdk,E_nb2->lf_an_ajdk,&(F->lf_knat));
	    GLOBAL_pthread_toc(); 
	    /* skip */}
	  if (nb1==nb2){
	    if (verbose){ printf(" %% knkt %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_kn->rpop_j,E_nb1->M_kn->cpop_j,E_nb2->M_kn->rpop_j,E_nb2->M_kn->cpop_j,F->lf_knkt->length);}
	    GLOBAL_pthread_tic(); 
	    wrap_AnAt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_kn,E_nb2->M_kn,D->A_ajdk,E_nb1->lf_kn_ajdk,E_nb2->lf_kn_ajdk,&(F->lf_knkt));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1==nb2){ } */}
	  if (nb1<nb2){
	    if (verbose){ printf(" %% knkt %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_kn->rpop_j,E_nb1->M_kn->cpop_j,E_nb2->M_kn->rpop_j,E_nb2->M_kn->cpop_j,F->lf_knkt->length);}
	    GLOBAL_pthread_tic(); 
	    wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_kn,E_nb2->M_kn,D->A_ajdk,E_nb1->lf_kn_ajdk,E_nb2->lf_kn_ajdk,&(F->lf_knkt));
	    GLOBAL_pthread_toc(); 
	    /* if (nb1<nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% YnWt_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (1){
	    if (verbose){ printf(" %% unyt %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_un->rpop_j,E_nb1->M_un->cpop_j,E_nb2->M_yn->rpop_j,E_nb2->M_yn->cpop_j,F->lf_unyt->length);}
	    L2_transpose(F->lf_unyt,F_trn->lf_ynut);
	    /* fill */}
	  if (nb1>nb2){
	    if (verbose){ printf(" %% unut %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_un->rpop_j,E_nb1->M_un->cpop_j,E_nb2->M_un->rpop_j,E_nb2->M_un->cpop_j,F->lf_unut->length);}
	    L2_transpose(F->lf_unut,F_trn->lf_unut);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (1){
	    if (verbose){ printf(" %% knat %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_kn->rpop_j,E_nb1->M_kn->cpop_j,E_nb2->M_an->rpop_j,E_nb2->M_an->cpop_j,F->lf_knat->length);}
	    L2_transpose(F->lf_knat,F_trn->lf_ankt);
	    /* fill */}
	  if (nb1>nb2){
	    if (verbose){ printf(" %% knkt %d-x-%d and %d-x-%d--> %d \n",E_nb1->M_kn->rpop_j,E_nb1->M_kn->cpop_j,E_nb2->M_kn->rpop_j,E_nb2->M_kn->cpop_j,F->lf_knkt->length);}
	    L2_transpose(F->lf_knkt,F_trn->lf_knkt);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% YnWt_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(E_nb1->A_bmr_j_rtn,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j_rtn: ");
	  bprintf(E_nb1->A_bmr_j_rmv,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j_rmv: ");
	  bprintf(E_nb2->M_Zn->mr_j ,D->bitj,1,E_nb2->Z_nrows," %% Z2_bmr_j    : ");
	  lfprintf(F->lf_unwt," %% lf_unwt: ");
	  lfprintf(F->lf_knzt," %% lf_knzt: ");
	  lfprintf(F->lf_hnvt," %% lf_hnvt: ");
	  lfprintf(F->lf_ynut," %% lf_ynut: ");
	  lfprintf(F->lf_unyt," %% lf_unyt: ");
	  lfprintf(F->lf_unut," %% lf_unut: ");
	  lfprintf(F->lf_ankt," %% lf_ankt: ");
	  lfprintf(F->lf_knat," %% lf_knat: ");
	  lfprintf(F->lf_knkt," %% lf_knkt: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_YnWt_stage_a1]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_YnWt_stage_a2(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  struct M_handle *M_in1=NULL,*M_in2=NULL,**M_out_p=NULL;
  struct L_handle *L_in=NULL;
  int trm_flag=0;
  if (strstr(D->QC_strategy,"YnWt")){
    if (verbose){ printf(" %% [entering bcc_lrup_QC_YnWt_stage_a2]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->M_unwt, F_[nbx]->M_knzt, F_[nbx]->M_hnvt, F_[nbx]->M_ynut, F_[nbx]->M_unyt, F_[nbx]->M_unut, F_[nbx]->M_ankt, F_[nbx]->M_knat, F_[nbx]->M_knkt.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  M_in1 = E_nb1->M_un; M_in2 = E_nb2->M_Wn; L_in = F->lf_unwt; M_out_p = &(F->M_unwt); trm_flag = F->M_unwt_trm; if (verbose>0){ printf(" %% bcc_M_unwt: trm_flag %d\n",trm_flag);} bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p);
	  /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  M_in1 = E_nb1->M_kn; M_in2 = E_nb2->M_zn; L_in = F->lf_knzt; M_out_p = &(F->M_knzt); trm_flag = F->M_knzt_trm; if (verbose>0){ printf(" %% bcc_M_knzt: trm_flag %d\n",trm_flag);} bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p);
	  M_in1 = E_nb1->M_hn; M_in2 = E_nb2->M_vn; L_in = F->lf_hnvt; M_out_p = &(F->M_hnvt); trm_flag = F->M_hnvt_trm; if (verbose>0){ printf(" %% bcc_M_hnvt: trm_flag %d\n",trm_flag);} bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p);
	  /* if bother */}
	if (D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  M_in1 = E_nb1->M_yn; M_in2 = E_nb2->M_un; L_in = F->lf_ynut; M_out_p = &(F->M_ynut); trm_flag = F->M_ynut_trm; if (verbose>0){ printf(" %% bcc_M_ynut: trm_flag %d\n",trm_flag);} bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p);
	  if (0){ M_in1 = E_nb1->M_un; M_in2 = E_nb2->M_yn; L_in = F->lf_unyt; M_out_p = &(F->M_unyt); trm_flag = F->M_unyt_trm; if (verbose>0){ printf(" %% bcc_M_unyt: trm_flag %d\n",trm_flag);} bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p); /* skip */}
	  if (nb1<=nb2){ M_in1 = E_nb1->M_un; M_in2 = E_nb2->M_un; L_in = F->lf_unut; M_out_p = &(F->M_unut); trm_flag = F->M_unut_trm; if (verbose>0){ printf(" %% bcc_M_unut: trm_flag %d\n",trm_flag);} bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p); /* if (nb1<=nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  M_in1 = E_nb1->M_an; M_in2 = E_nb2->M_kn; L_in = F->lf_ankt; M_out_p = &(F->M_ankt); trm_flag = F->M_ankt_trm; if (verbose>0){ printf(" %% bcc_M_ankt: trm_flag %d\n",trm_flag);} bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p);
	  if (0){ M_in1 = E_nb1->M_kn; M_in2 = E_nb2->M_an; L_in = F->lf_knat; M_out_p = &(F->M_knat); trm_flag = F->M_knat_trm; if (verbose>0){ printf(" %% bcc_M_knat: trm_flag %d\n",trm_flag);} bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p); /* skip */}
	  if (nb1<=nb2){ M_in1 = E_nb1->M_kn; M_in2 = E_nb2->M_kn; L_in = F->lf_knkt; M_out_p = &(F->M_knkt); trm_flag = F->M_knkt_trm; if (verbose>0){ printf(" %% bcc_M_knkt: trm_flag %d\n",trm_flag);} bcc_M_YnWt_excerpt(trm_flag,M_in1,M_in2,L_in,M_out_p); /* if (nb1<=nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% xcalc YnWt: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose){ printf(" %% [finished bcc_lrup_QC_YnWt_stage_a2]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_YnWt_stage_a3(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QC_strategy,"YnWt")){
    if (verbose){ printf(" %% [entering bcc_lrup_QC_YnWt_stage_a3]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_kt_r_unwt_s_zn, F_[nbx]->lf_kt_r_knzt_s_zn, F_[nbx]->lf_kt_r_hnvt_s_zn, F_[nbx]->lf_at_t_ynut_r_kn, F_[nbx]->lf_kt_r_unyt_t_an, F_[nbx]->lf_kt_r_unut_r_kn, F_[nbx]->lf_at_t_ankt_r_kn, F_[nbx]->lf_kt_r_knat_t_an, F_[nbx]->lf_kt_r_knkt_r_kn\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_unwt_trm,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,E_nb1->M_un,E_nb2->M_Wn,E_nb2->M_St,E_nb2->M_zt,D->A_ajdk,F->M_unwt,&(F->lf_kt_r_unwt_s_zn));
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_knzt_trm,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,E_nb1->M_kn,E_nb2->M_zn,E_nb2->M_St,E_nb2->M_zt,D->A_ajdk,F->M_knzt,&(F->lf_kt_r_knzt_s_zn));
	  bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_hnvt_trm,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,E_nb1->M_hn,E_nb2->M_vn,E_nb2->M_St,E_nb2->M_zt,D->A_ajdk,F->M_hnvt,&(F->lf_kt_r_hnvt_s_zn));
	  /* if bother */}
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_ynut_trm,n_spacing_A,E_nb1->M_at,E_nb1->M_tt,E_nb1->M_yn,E_nb2->M_un,E_nb2->M_rt,E_nb2->M_kt,D->A_ajdk,F->M_ynut,&(F->lf_at_t_ynut_r_kn));
	  if (0){ bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_unyt_trm,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,E_nb1->M_un,E_nb2->M_yn,E_nb2->M_tt,E_nb2->M_at,D->A_ajdk,F->M_unyt,&(F->lf_kt_r_unyt_t_an)); /* skip */}
	  if (nb1<=nb2){bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_unut_trm,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,E_nb1->M_un,E_nb2->M_un,E_nb2->M_rt,E_nb2->M_kt,D->A_ajdk,F->M_unut,&(F->lf_kt_r_unut_r_kn)); /* if (nb1<=nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_ankt_trm,n_spacing_A,E_nb1->M_at,E_nb1->M_tt,E_nb1->M_an,E_nb2->M_kn,E_nb2->M_rt,E_nb2->M_kt,D->A_ajdk,F->M_ankt,&(F->lf_at_t_ankt_r_kn));
	  if (0){ bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_knat_trm,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,E_nb1->M_kn,E_nb2->M_an,E_nb2->M_tt,E_nb2->M_at,D->A_ajdk,F->M_knat,&(F->lf_kt_r_knat_t_an)); /* skip */}
	  if (nb1<=nb2){ bcc_lf_At_T_YnWt_S_Zn_excerpt(F->M_knkt_trm,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,E_nb1->M_kn,E_nb2->M_kn,E_nb2->M_rt,E_nb2->M_kt,D->A_ajdk,F->M_knkt,&(F->lf_kt_r_knkt_r_kn)); /* if (nb1<=nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% At_T_YnWt_S_Zn_ww: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (1){ L2_duplicate(F->lf_kt_r_unyt_t_an,F_trn->lf_at_t_ynut_r_kn); /* fill */}
	  if (nb1>nb2){ L2_duplicate(F->lf_kt_r_unut_r_kn,F_trn->lf_kt_r_unut_r_kn); /* if (nb1>nb2){ } */}
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (1){ L2_duplicate(F->lf_kt_r_knat_t_an,F_trn->lf_at_t_ankt_r_kn); /* fill */}
	  if (nb1>nb2){ L2_duplicate(F->lf_kt_r_knkt_r_kn,F_trn->lf_kt_r_knkt_r_kn); /* if (nb1>nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(D->A_bmc_j_rtn   ,D->bitj,1,D->A_ncols," %% D_bmc_j_rtn: ");
	  bprintf(E_nb1->M_An->mc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j   : ");
	  bprintf(E_nb2->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j   : ");
	  lfprintf(F->lf_kt_r_unwt_s_zn," %% lf_kt_r_unwt_s_zn: ");
	  lfprintf(F->lf_kt_r_knzt_s_zn," %% lf_kt_r_knzt_s_zn: ");
	  lfprintf(F->lf_kt_r_hnvt_s_zn," %% lf_kt_r_hnvt_s_zn: ");
	  lfprintf(F->lf_at_t_ynut_r_kn," %% lf_at_t_ynut_r_kn: ");
	  lfprintf(F->lf_kt_r_unyt_t_an," %% lf_kt_r_unyt_t_an: ");
	  lfprintf(F->lf_kt_r_unut_r_kn," %% lf_kt_r_unut_r_kn: ");
	  lfprintf(F->lf_at_t_ankt_r_kn," %% lf_at_t_ankt_r_kn: ");
	  lfprintf(F->lf_kt_r_knat_t_an," %% lf_kt_r_knat_t_an: ");
	  lfprintf(F->lf_kt_r_knkt_r_kn," %% lf_kt_r_knkt_r_kn: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_YnWt_stage_a3]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_YnWt_stage_b0(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QC_strategy,"YnWt")){ 
    if (verbose){ printf(" %% [entering bcc_lrup_QC_YnWt_stage_b0]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->lf_attn, E_[nb]->lf_jttn, E_[nb]->lf_ktrn, E_[nb]->lf_htrn, E_[nb]->lf_ztsn, E_[nb]->lf_vtsn.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->A_rbother && D->A_cbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_at,E->M_tt,NULL,NULL,NULL,&(E->lf_attn));
	GLOBAL_pthread_toc(); 
	GLOBAL_pthread_tic(); 
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_jt,E->M_tt,NULL,NULL,NULL,&(E->lf_jttn));
	GLOBAL_pthread_toc(); 
	GLOBAL_pthread_tic(); 
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_kt,E->M_rt,NULL,NULL,NULL,&(E->lf_ktrn));
	GLOBAL_pthread_toc(); 
	GLOBAL_pthread_tic(); 
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_ht,E->M_rt,NULL,NULL,NULL,&(E->lf_htrn));
	GLOBAL_pthread_toc(); 
	/* if bother */}
      if (E->Z_rbother && D->A_cbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_zt,E->M_St,NULL,NULL,NULL,&(E->lf_ztsn));
	GLOBAL_pthread_toc();
	GLOBAL_pthread_tic(); 
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_vt,E->M_St,NULL,NULL,NULL,&(E->lf_vtsn));
	GLOBAL_pthread_toc();
	/* if bother */} 
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% ZtSn_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	printf(" %% nb %d\n",nb);
	bprintf(D->A_bmc_j_rtn,D->bitj,1,D->A_ncols," %% A_bmc_j_rtn: ");
	bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols," %% A_bmc_j_rmv: ");
	lfprintf(E->lf_attn," %% lf_attn: ");
	lfprintf(E->lf_jttn," %% lf_jttn: ");
	lfprintf(E->lf_ktrn," %% lf_ktrn: ");
	lfprintf(E->lf_htrn," %% lf_htrn: ");
	lfprintf(E->lf_ztsn," %% lf_ztsn: ");
	lfprintf(E->lf_vtsn," %% lf_vtsn: ");
	/* for (nb=0;nb<nbins;nb++){ } */}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_YnWt_stage_b0]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_YnWt_stage_b1(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  int n_spacing_A = SPACING_j;
  if (strstr(D->QC_strategy,"YnWt") && !strstr(D->QC_strategy,"store one")){ 
    if (verbose){ printf(" %% [entering bcc_lrup_QC_YnWt_stage_b1]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->lf_ztsvn, E_[nb]->lf_attjn, E_[nb]->lf_ktrhn.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->Z_rbother && D->A_cbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_zt,E->M_St,E->M_vt,D->A_ajdk,D->A_ajdk,E->lf_ztsn,E->lf_vtsn,&(E->lf_ztsvn));
	GLOBAL_pthread_toc();
	/* if bother */}
      if (E->A_rbother && D->A_cbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_at,E->M_tt,E->M_jt,D->A_ajdk,D->A_ajdk,E->lf_attn,E->lf_jttn,&(E->lf_attjn));
	GLOBAL_pthread_toc();
	GLOBAL_pthread_tic();
	wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_kt,E->M_rt,E->M_ht,D->A_ajdk,D->A_ajdk,E->lf_ktrn,E->lf_htrn,&(E->lf_ktrhn));
	GLOBAL_pthread_toc();
	/* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% ZtSWn_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	printf(" %% nb %d\n",nb);
	bprintf(D->A_bmc_j_rtn,D->bitj,1,D->A_ncols," %% A_bmc_j_rtn: ");
	bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols," %% A_bmc_j_rmv: ");
	bprintf(E->M_Tn->mc_j ,D->bitj,1,D->T_ncols," %% T_bmc_j    : ");
	if (E->Z_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_ztsvn," %% lf_ztsvn: "); /* if bother */} else{ printf(" %% no lf_ztsvn\n");}
	if (E->A_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_attjn," %% lf_attjn: "); /* if bother */} else{ printf(" %% no lf_attjn\n");}
	if (E->A_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_ktrhn," %% lf_ktrhn: "); /* if bother */} else{ printf(" %% no lf_ktrhn\n");}
	/* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_YnWt_stage_b1]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_YnWt_stage_b2(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QC_strategy,"YnWt") && !strstr(D->QC_strategy,"store one")){ 
    if (verbose){ printf(" %% [entering bcc_lrup_QC_YnWt_stage_b2]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_attjn____vtszn, F_[nbx]->lf_attjn____jttan, F_[nbx]->lf_attjn____htrkn, F_[nbx]->lf_ktrhn____jttan, F_[nbx]->lf_ktrhn____htrkn.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
  	  GLOBAL_pthread_tic(); 
	  wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_at,E_nb1->M_tt,E_nb1->M_jt,E_nb2->M_vt,E_nb2->M_St,E_nb2->M_zt,D->A_ajdk,E_nb1->lf_attjn,E_nb2->lf_ztsvn,&(F->lf_attjn____vtszn));
  	  GLOBAL_pthread_toc();
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1<=nb2){
	    GLOBAL_pthread_tic(); 
	    wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_at,E_nb1->M_tt,E_nb1->M_jt,E_nb2->M_jt,E_nb2->M_tt,E_nb2->M_at,D->A_ajdk,E_nb1->lf_attjn,E_nb2->lf_attjn,&(F->lf_attjn____jttan));
	    GLOBAL_pthread_toc();
	    /* if (nb1<=nb2){ } */}
  	  GLOBAL_pthread_tic(); 
	  wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_at,E_nb1->M_tt,E_nb1->M_jt,E_nb2->M_ht,E_nb2->M_rt,E_nb2->M_kt,D->A_ajdk,E_nb1->lf_attjn,E_nb2->lf_ktrhn,&(F->lf_attjn____htrkn));
  	  GLOBAL_pthread_toc();
	  if (0){
	    GLOBAL_pthread_tic(); 
	    wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,E_nb1->M_ht,E_nb2->M_jt,E_nb2->M_tt,E_nb2->M_at,D->A_ajdk,E_nb1->lf_ktrhn,E_nb2->lf_attjn,&(F->lf_ktrhn____jttan));
	    GLOBAL_pthread_toc();
	    /* skip */}
	  if (nb1<=nb2){
	    GLOBAL_pthread_tic(); 
	    wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,E_nb1->M_ht,E_nb2->M_ht,E_nb2->M_rt,E_nb2->M_kt,D->A_ajdk,E_nb1->lf_ktrhn,E_nb2->lf_ktrhn,&(F->lf_ktrhn____htrkn));
	    GLOBAL_pthread_toc();
	    /* if (nb1<=nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% AtTYn____WtSZn_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_duplicate(F->lf_attjn____jttan,F_trn->lf_attjn____jttan);
	    /* if (nb1>nb2){ } */}
	  if (1){
	    L2_duplicate(F->lf_ktrhn____jttan,F_trn->lf_attjn____htrkn);
	    /* fill */}
	  if (nb1>nb2){
	    L2_duplicate(F->lf_ktrhn____htrkn,F_trn->lf_ktrhn____htrkn);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(D->A_bmc_j_rtn   ,D->bitj,1,D->A_ncols," %% A_bmc_j_rtn: ");
	  bprintf(D->A_bmc_j_rmv   ,D->bitj,1,D->A_ncols," %% A_bmc_j_rmv: ");
	  bprintf(E_nb1->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T_bmc_j    : ");
	  lfprintf(F->lf_attjn____vtszn," %% lf_attjn____vtszn: ");
	  lfprintf(F->lf_attjn____jttan," %% lf_attjn____jttan: ");
	  lfprintf(F->lf_attjn____htrkn," %% lf_attjn____htrkn: ");
	  lfprintf(F->lf_ktrhn____jttan," %% lf_ktrhn____jttan: ");
	  lfprintf(F->lf_ktrhn____htrkn," %% lf_ktrhn____htrkn: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_YnWt_stage_b2]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_YnWt_stage_b3(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nb=0; struct bcc_single *E=NULL;
  int n_spacing_A = SPACING_j;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  int ns_j=0,ns_a=0,ns_b=0;
  if (nbins>0 && strstr(D->QC_strategy,"YnWt") && strstr(D->QC_strategy,"store one")){ 
    if (verbose){ printf(" %% [entering bcc_lrup_QC_YnWt_stage_b3]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->lf_ztsvn, E_[nb]->lf_attjn, E_[nb]->lf_ktrhn for each ns_j individually.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    ns_j=0;
    while (ns_j<E_[0]->M_Tt->rpop_j){
      ns_a = E_[0]->M_Tt->m_a_[ns_j]; ns_b = E_[0]->M_Tt->m_b_[ns_j];
      if (verbose){ printf(" %% ns_j %d, ns_b %d, ns_a %d\n",ns_j,ns_b,ns_a);}
      GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	if (E->Z_rbother && D->A_cbother && D->A_cbother){ 
	  GLOBAL_pthread_tic();
	  wrap_AttYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E->M_zt,E->M_St,ns_j,E->M_vt,D->A_ajdk,D->A_ajdk,E->lf_ztsn,E->lf_vtsn,&(E->lf_ztsvn));
	  GLOBAL_pthread_toc();
	  /* if bother */}
	if (E->A_rbother && D->A_cbother && D->A_cbother){ 
	  GLOBAL_pthread_tic();
	  wrap_AttYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E->M_at,E->M_tt,ns_j,E->M_jt,D->A_ajdk,D->A_ajdk,E->lf_attn,E->lf_jttn,&(E->lf_attjn));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  wrap_AttYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E->M_kt,E->M_rt,ns_j,E->M_ht,D->A_ajdk,D->A_ajdk,E->lf_ktrn,E->lf_htrn,&(E->lf_ktrhn));
	  GLOBAL_pthread_toc();
	  /* if bother */}
	/* for (nb=0;nb<nbins;nb++){ } */}
      GLOBAL_pthread_tuc();
      GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	  if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	    GLOBAL_pthread_tic(); 
	    wrap_AttYn____WtsZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_at,E_nb1->M_tt,ns_j,E_nb1->M_jt,E_nb2->M_vt,E_nb2->M_St,E_nb2->M_zt,D->A_ajdk,E_nb1->lf_attjn,E_nb2->lf_ztsvn,&(F->lf_attjn____vtszn));
	    GLOBAL_pthread_toc();
	    /* if bother */}
	  if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	    if (nb1<=nb2){
	      GLOBAL_pthread_tic(); 
	      wrap_AttYn____WtsZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_at,E_nb1->M_tt,ns_j,E_nb1->M_jt,E_nb2->M_jt,E_nb2->M_tt,E_nb2->M_at,D->A_ajdk,E_nb1->lf_attjn,E_nb2->lf_attjn,&(F->lf_attjn____jttan));
	      GLOBAL_pthread_toc();
	      /* if (nb1<=nb2){ } */}
	    GLOBAL_pthread_tic(); 
	    wrap_AttYn____WtsZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_at,E_nb1->M_tt,ns_j,E_nb1->M_jt,E_nb2->M_ht,E_nb2->M_rt,E_nb2->M_kt,D->A_ajdk,E_nb1->lf_attjn,E_nb2->lf_ktrhn,&(F->lf_attjn____htrkn));
	    GLOBAL_pthread_toc();
	    if (0){
	      GLOBAL_pthread_tic(); 
	      wrap_AttYn____WtsZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,ns_j,E_nb1->M_ht,E_nb2->M_jt,E_nb2->M_tt,E_nb2->M_at,D->A_ajdk,E_nb1->lf_ktrhn,E_nb2->lf_attjn,&(F->lf_ktrhn____jttan));
	      GLOBAL_pthread_toc();
	      /* skip */}
	    if (nb1<=nb2){
	      GLOBAL_pthread_tic(); 
	      wrap_AttYn____WtsZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,E_nb1->M_kt,E_nb1->M_rt,ns_j,E_nb1->M_ht,E_nb2->M_ht,E_nb2->M_rt,E_nb2->M_kt,D->A_ajdk,E_nb1->lf_ktrhn,E_nb2->lf_ktrhn,&(F->lf_ktrhn____htrkn));
	      GLOBAL_pthread_toc();
	      /* if (nb1<=nb2){ } */}
	    /* if bother */}
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      GLOBAL_pthread_tuc();
      ns_j++; /* while (ns_j<E_[0]->M_Tt->rpop_j){ } */}
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% ztsvn, attjn, ktrhn, attjn____vtszn, etc.: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_duplicate(F->lf_attjn____jttan,F_trn->lf_attjn____jttan);
	    /* if (nb1>nb2){ } */}
	  if (1){
	    L2_duplicate(F->lf_ktrhn____jttan,F_trn->lf_attjn____htrkn);
	    /* fill */}
	  if (nb1>nb2){
	    L2_duplicate(F->lf_ktrhn____htrkn,F_trn->lf_ktrhn____htrkn);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	printf(" %% nb %d\n",nb);
	bprintf(D->A_bmc_j_rtn,D->bitj,1,D->A_ncols," %% A_bmc_j_rtn: ");
	bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols," %% A_bmc_j_rmv: ");
	bprintf(E->M_Tn->mc_j ,D->bitj,1,D->T_ncols," %% T_bmc_j    : ");
	if (E->Z_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_ztsvn," %% lf_ztsvn: "); /* if bother */} else{ printf(" %% no lf_ztsvn\n");}
	if (E->A_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_attjn," %% lf_attjn: "); /* if bother */} else{ printf(" %% no lf_attjn\n");}
	if (E->A_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_ktrhn," %% lf_ktrhn: "); /* if bother */} else{ printf(" %% no lf_ktrhn\n");}
	/* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(D->A_bmc_j_rtn   ,D->bitj,1,D->A_ncols," %% A_bmc_j_rtn: ");
	  bprintf(D->A_bmc_j_rmv   ,D->bitj,1,D->A_ncols," %% A_bmc_j_rmv: ");
	  bprintf(E_nb1->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% T_bmc_j    : ");
	  lfprintf(F->lf_attjn____vtszn," %% lf_attjn____vtszn: ");
	  lfprintf(F->lf_attjn____jttan," %% lf_attjn____jttan: ");
	  lfprintf(F->lf_attjn____htrkn," %% lf_attjn____htrkn: ");
	  lfprintf(F->lf_ktrhn____jttan," %% lf_ktrhn____jttan: ");
	  lfprintf(F->lf_ktrhn____htrkn," %% lf_ktrhn____htrkn: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_YnWt_stage_b3]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_YnWt_stage_c(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  if (strstr(D->QC_strategy,"YnWt")){ 
    if (verbose){ printf(" %% [entering bcc_lrup_QC_YnWt_stage_c]\n");}
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	  bprintf(D->T_bmc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	  lfprintf(F->lf_At_T_YnWt_S_Zn," %% pre lf_At_T_YnWt_S_Zn: ");
	  lfprintf(F->lf_At_T_AnZt_S_Zn," %% pre lf_At_T_AnZt_S_Zn: ");
	  lfprintf(F->lf_At_T_YnYt_T_An," %% pre lf_At_T_YnYt_T_An: ");
	  lfprintf(F->lf_At_T_AnAt_T_An," %% pre lf_At_T_AnAt_T_An: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% updating F_[nbx]->lf_At_T_YnWt_S_Zn, etc..\n");}
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  L2_clean(F->lf_At_T_YnWt_S_Zn,D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_YnWt_S_Zn,F->lf_kt_r_unwt_s_zn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  L2_clean(F->lf_At_T_AnZt_S_Zn,D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnZt_S_Zn,F->lf_attjn____vtszn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnZt_S_Zn,F->lf_kt_r_knzt_s_zn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnZt_S_Zn,F->lf_kt_r_hnvt_s_zn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  /* if bother */}
	if (D->A_cbother && D->Y_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  L2_clean(F->lf_At_T_YnYt_T_An,D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_YnYt_T_An,F->lf_at_t_ynut_r_kn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_YnYt_T_An,F->lf_kt_r_unyt_t_an,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_YnYt_T_An,F->lf_kt_r_unut_r_kn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  /* if bother */}
	if (D->A_cbother && D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  L2_clean(F->lf_At_T_AnAt_T_An,D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnAt_T_An,F->lf_attjn____jttan,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnAt_T_An,F->lf_attjn____htrkn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnAt_T_An,F->lf_ktrhn____jttan,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnAt_T_An,F->lf_ktrhn____htrkn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnAt_T_An,F->lf_at_t_ankt_r_kn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnAt_T_An,F->lf_kt_r_knat_t_an,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  LL2_plustimesequals(F->lf_At_T_AnAt_T_An,F->lf_kt_r_knkt_r_kn,-1.0 , D->A_ncols,D->A_bmc_b,D->A_bmc_j_rtn,D->T_ncols,D->T_bmc_b,D->T_bmc_j);
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols," %% A1_bmc_j: ");
	  bprintf(D->T_bmc_j,D->bitj,1,D->T_ncols," %% T2_bmc_j: ");
	  lfprintf(F->lf_At_T_YnWt_S_Zn," %% pos lf_At_T_YnWt_S_Zn: ");
	  lfprintf(F->lf_At_T_AnZt_S_Zn," %% pos lf_At_T_AnZt_S_Zn: ");
	  lfprintf(F->lf_At_T_YnYt_T_An," %% pos lf_At_T_YnYt_T_An: ");
	  lfprintf(F->lf_At_T_AnAt_T_An," %% pos lf_At_T_AnAt_T_An: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_YnWt_stage_c]\n");}
    /* if strategy */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_lrup_QR_YnWt_stage_0(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QR_strategy,"YnWt") && !strstr(D->QR_strategy,"condense")){
    if (verbose){ printf(" %% [entering bcc_lrup_QR_YnWt_stage_0]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->lf_jn_ajdk, E_[nb]->lf_vn_ajdk.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->A_rbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_jn,D->A_ajdk,&(E->lf_jn_ajdk));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E->Z_rbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,E->M_vn,D->A_ajdk,&(E->lf_vn_ajdk));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc(); 
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% An_ajdk_v: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");  
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){ 
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	printf(" %% nb %d\n",nb);
	bprintf(E->A_bmr_j_rtn,D->bitj,1,E->A_nrows," %% A_bmr_j_rtn: ");
	bprintf(E->M_jn->mr_j ,D->bitj,1,E->A_nrows," %% M_jn->mr_j : ");
	bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols," %% A_bmr_j_rtn: ");
	bprintf(E->M_jn->mc_j ,D->bitj,1,D->A_ncols," %% M_jn->mc_j : ");
	lfprintf(E->lf_jn_ajdk," %% lf_jn_ajdk->lf: ");
	bprintf(E->M_vn->mc_j ,D->bitj,1,D->A_ncols," %% M_vn->mc_j : ");
	lfprintf(E->lf_vn_ajdk," %% lf_vn_ajdk->lf: ");
	/* for (nb=0;nb<nbins;nb++){ } */}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QR_YnWt_stage_0]\n");}
    /* if strategy */}
}

void bcc_lrup_QR_YnWt_stage_1(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_b;
  if (strstr(D->QR_strategy,"YnWt") && !strstr(D->QR_strategy,"condense")){
    if (verbose){ printf(" %% [entering bcc_lrup_QR_YnWt_stage_1]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_jnvt, F_[nbx]->lf_jnjt using AnZt_vv.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  GLOBAL_pthread_tic(); 
	  if (verbose>1){
	    printf(" %% M_jn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb1->M_jn->rpop_j,E_nb1->M_jn->rpop_b,E_nb1->M_jn->nrows,E_nb1->M_jn->cpop_j,E_nb1->M_jn->cpop_b,E_nb1->M_jn->ncols);
	    printf(" %% M_vn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb2->M_vn->rpop_j,E_nb2->M_vn->rpop_b,E_nb2->M_vn->nrows,E_nb2->M_vn->cpop_j,E_nb2->M_vn->cpop_b,E_nb2->M_vn->ncols);
	    /* if (verbose>1){ } */}
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_jn,E_nb2->M_vn,D->A_ajdk,E_nb1->lf_jn_ajdk,E_nb2->lf_vn_ajdk,&(F->lf_jnvt));
	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1==nb2){
	    GLOBAL_pthread_tic(); 
	    if (verbose>1){
	      printf(" %% M_jn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb1->M_jn->rpop_j,E_nb1->M_jn->rpop_b,E_nb1->M_jn->nrows,E_nb1->M_jn->cpop_j,E_nb1->M_jn->cpop_b,E_nb1->M_jn->ncols);
	      printf(" %% M_jn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb2->M_jn->rpop_j,E_nb2->M_jn->rpop_b,E_nb2->M_jn->nrows,E_nb2->M_jn->cpop_j,E_nb2->M_jn->cpop_b,E_nb2->M_jn->ncols);
	      /* if (verbose>1){ } */}
	    wrap_AnAt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_jn,E_nb2->M_jn,D->A_ajdk,E_nb1->lf_jn_ajdk,E_nb2->lf_jn_ajdk,&(F->lf_jnjt));
	    GLOBAL_pthread_toc();
	    /* if (nb1==nb2){ } */}
	  if (nb1<nb2){
	    GLOBAL_pthread_tic(); 
	    if (verbose>1){
	      printf(" %% M_jn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb1->M_jn->rpop_j,E_nb1->M_jn->rpop_b,E_nb1->M_jn->nrows,E_nb1->M_jn->cpop_j,E_nb1->M_jn->cpop_b,E_nb1->M_jn->ncols);
	      printf(" %% M_jn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb2->M_jn->rpop_j,E_nb2->M_jn->rpop_b,E_nb2->M_jn->nrows,E_nb2->M_jn->cpop_j,E_nb2->M_jn->cpop_b,E_nb2->M_jn->ncols);
	      /* if (verbose>1){ } */}
	    wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E_nb1->M_jn,E_nb2->M_jn,D->A_ajdk,E_nb1->lf_jn_ajdk,E_nb2->lf_jn_ajdk,&(F->lf_jnjt));
	    GLOBAL_pthread_toc();
	    /* if (nb1<nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% YnWt_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_transpose(F->lf_jnjt,F_trn->lf_jnjt);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(E_nb1->A_bmr_j_rtn,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j_rtn: ");
	  bprintf(E_nb2->M_Zn->mr_j ,D->bitj,1,E_nb2->Z_nrows," %% Z2_bmr_j    : ");
	  lfprintf(F->lf_jnvt," %% lf_jnvt: ");
	  bprintf(E_nb1->A_bmr_j_rtn,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j_rtn: ");
	  bprintf(E_nb2->A_bmr_j_rtn,D->bitj,1,E_nb2->A_nrows," %% A2_bmr_j_rtn: ");
	  lfprintf(F->lf_jnjt," %% lf_jnjt: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QR_YnWt_stage_1]\n");}
    /* if strategy */}
}

void bcc_lrup_QR_YnWt_stage_2(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL,*F_trn=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_b; int n_spacing_B = SPACING_j;
  if (strstr(D->QR_strategy,"YnWt") && strstr(D->QR_strategy,"condense")){
    if (verbose){ printf(" %% [entering bcc_lrup_QR_YnWt_stage_2]\n");}
    if (verbose){ printf(" %% calculating F_[nbx]->lf_jnvt, F_[nbx]->lf_jnjt using M_At_to_L2 and AnZt_mm.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ E_nb1 = E_[nb1];
	if (D->A_cbother && E_nb1->A_rbother){ 
	  GLOBAL_pthread_tic(); 
	  wrap_M_At_to_L2_run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,n_spacing_A,E_nb1->M_jt,D->A_ajdk,&(E_nb1->lf_jn));
	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->A_cbother && E_nb1->Z_rbother){ 
	  GLOBAL_pthread_tic(); 
	  wrap_M_At_to_L2_run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,n_spacing_A,E_nb1->M_vt,D->A_ajdk,&(E_nb1->lf_vn));
	  GLOBAL_pthread_toc();
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ } */}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% M_At_to_L2: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    if (verbose>2){
      for (nb1=0;nb1<nbins;nb1++){ E_nb1 = E_[nb1];
	printf(" %% E_[%d]->M_jn (%d,%d,%d)-x-(%d,%d,%d)\n",nb1,E_nb1->M_jn->rpop_j,E_nb1->M_jn->rpop_b,E_nb1->M_jn->nrows,E_nb1->M_jn->cpop_j,E_nb1->M_jn->cpop_b,E_nb1->M_jn->ncols);
	printf(" %% E_[%d]->lf_jn (%d,%d) <-- (%s,%s)\n",nb1,E_nb1->lf_jn->row_stride,E_nb1->lf_jn->col_stride,SPACING_name[E_nb1->lf_jn->spacing_row],SPACING_name[E_nb1->lf_jn->spacing_col]);
	printf(" %% E_[%d]->M_vn (%d,%d,%d)-x-(%d,%d,%d)\n",nb1,E_nb1->M_vn->rpop_j,E_nb1->M_vn->rpop_b,E_nb1->M_vn->nrows,E_nb1->M_vn->cpop_j,E_nb1->M_vn->cpop_b,E_nb1->M_vn->ncols);
	printf(" %% E_[%d]->lf_vn (%d,%d) <-- (%s,%s)\n",nb1,E_nb1->lf_vn->row_stride,E_nb1->lf_vn->col_stride,SPACING_name[E_nb1->lf_vn->spacing_row],SPACING_name[E_nb1->lf_vn->spacing_col]);
	/* for (nb1=0;nb1<nbins;nb1++){ } */}
      /* if (verbose>2){ } */}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ 
	  GLOBAL_pthread_tic(); 
	  if (verbose>1){
	    printf(" %% M_jn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb1->M_jn->rpop_j,E_nb1->M_jn->rpop_b,E_nb1->M_jn->nrows,E_nb1->M_jn->cpop_j,E_nb1->M_jn->cpop_b,E_nb1->M_jn->ncols);
	    printf(" %% M_vn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb2->M_vn->rpop_j,E_nb2->M_vn->rpop_b,E_nb2->M_vn->nrows,E_nb2->M_vn->cpop_j,E_nb2->M_vn->cpop_b,E_nb2->M_vn->ncols);
	    /* if (verbose>1){ } */}
	  wrap_AnZt_mm__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E_nb1->lf_jn,E_nb2->lf_vn,&(F->lf_jnvt));
	  GLOBAL_pthread_toc(); /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1<=nb2){
	    GLOBAL_pthread_tic(); 
	    if (verbose>1){
	      printf(" %% M_jn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb1->M_jn->rpop_j,E_nb1->M_jn->rpop_b,E_nb1->M_jn->nrows,E_nb1->M_jn->cpop_j,E_nb1->M_jn->cpop_b,E_nb1->M_jn->ncols);
	      printf(" %% M_jn (%d,%d,%d)-x-(%d,%d,%d)\n",E_nb2->M_jn->rpop_j,E_nb2->M_jn->rpop_b,E_nb2->M_jn->nrows,E_nb2->M_jn->cpop_j,E_nb2->M_jn->cpop_b,E_nb2->M_jn->ncols);
	      /* if (verbose>1){ } */}
	    wrap_AnZt_mm__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),E_nb1->lf_jn,E_nb2->lf_jn,&(F->lf_jnjt));
	    GLOBAL_pthread_toc();
	    /* if (nb1<=nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% AnZt_mm: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    if (verbose>2){
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx];
	  printf(" %% F_[%d,%d]->lf_jnjt (%d,%d) <-- (%s,%s)\n",nb1,nb2,F->lf_jnjt->row_stride,F->lf_jnjt->col_stride,SPACING_name[F->lf_jnjt->spacing_row],SPACING_name[F->lf_jnjt->spacing_col]);
	  printf(" %% F_[%d,%d]->lf_jnvt (%d,%d) <-- (%s,%s)\n",nb1,nb2,F->lf_jnvt->row_stride,F->lf_jnvt->col_stride,SPACING_name[F->lf_jnvt->spacing_row],SPACING_name[F->lf_jnvt->spacing_col]);
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2; F_trn = F_[nb2+nb1*nbins];
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ 
	  if (nb1>nb2){
	    L2_transpose(F->lf_jnjt,F_trn->lf_jnjt);
	    /* if (nb1>nb2){ } */}
	  /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(E_nb1->A_bmr_j_rtn,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j_rtn: ");
	  bprintf(E_nb2->M_Zn->mr_j ,D->bitj,1,E_nb2->Z_nrows," %% Z2_bmr_j    : ");
	  lfprintf(F->lf_jnvt," %% lf_jnvt: ");
	  bprintf(E_nb1->A_bmr_j_rtn,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j_rtn: ");
	  bprintf(E_nb2->A_bmr_j_rtn,D->bitj,1,E_nb2->A_nrows," %% A2_bmr_j_rtn: ");
	  lfprintf(F->lf_jnjt," %% lf_jnjt: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QR_YnWt_stage_2]\n");}
    /* if strategy */}
}

void bcc_lrup_QR_YnWt_stage_3(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nbx=0,nb1=0,nb2=0; struct bcc_single *E_nb1=NULL,*E_nb2=NULL; struct bcc_double *F=NULL;
  if (strstr(D->QR_strategy,"YnWt")){
    if (verbose){ printf(" %% [entering bcc_lrup_QR_YnWt_stage_3]\n");}
    if (verbose){ printf(" %% updating F_[nbx]->lf_AnZt, F_[nbx]->lf_AnAt.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->Z_rbother){ dra_subtequals(F->lf_AnZt->lf,F->lf_AnZt->length,F->lf_jnvt->lf); /* if bother */}
	if (D->A_cbother && E_nb1->A_rbother && E_nb2->A_rbother){ dra_subtequals(F->lf_AnAt->lf,F->lf_AnAt->length,F->lf_jnjt->lf); /* if bother */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; F = F_[nbx]; E_nb1 = F->E_nb1; E_nb2 = F->E_nb2;
	  printf(" %% nb1 %d nb2 %d\n",nb1,nb2);
	  bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	  bprintf(E_nb2->M_Zn->mr_j,D->bitj,1,E_nb2->Z_nrows," %% Z2_bmr_j: ");
	  lfprintf(F->lf_AnZt," %% lf_AnZt: ");
	  bprintf(E_nb1->M_An->mr_j,D->bitj,1,E_nb1->A_nrows," %% A1_bmr_j: ");
	  bprintf(E_nb2->M_An->mr_j,D->bitj,1,E_nb2->A_nrows," %% A2_bmr_j: ");
	  lfprintf(F->lf_AnAt," %% lf_AnAt: ");
	  /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QR_YnWt_stage_3]\n");}
    /* if strategy */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_lrup_QC_ZtSWn_stage_0(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  int n_type = TYPE_pm;
  int n_spacing_A = SPACING_a;
  if (strstr(D->QC_strategy,"ZtSWn")){
    if (verbose){ printf(" %% [entering bcc_lrup_QC_ZtSWn_stage_0]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->lf_ktrn, E_[nb]->lf_utrn.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
      if (E->A_rbother && D->A_cbother){ 
	GLOBAL_pthread_tic(); 
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_kt,E->M_rt,NULL,NULL,NULL,&(E->lf_ktrn));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E->A_rbother && D->Y_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,E->M_ut,E->M_rt,NULL,NULL,NULL,&(E->lf_utrn));
	GLOBAL_pthread_toc(); /* if bother */} 
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% ZtSn_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){  
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	printf(" %% nb %d\n",nb);
	bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols," %% A_bmr_j_rtn: ");
	bprintf(E->M_kt->mr_j ,D->bitj,1,D->A_ncols," %% M_kt->mr_j : ");
	bprintf(E->A_bmr_j_rtn,D->bitj,1,E->A_nrows," %% A_bmr_j_rtn: ");
	bprintf(E->M_kt->mc_j ,D->bitj,1,E->A_nrows," %% M_kt->mc_j : ");
	lfprintf(E->lf_ktrn," %% lf_ktrn: ");
	bprintf(E->M_ut->mc_j ,D->bitj,1,E->A_nrows," %% M_ut->mc_j : ");
	lfprintf(E->lf_utrn," %% lf_utrn: ");
	/* for (nb=0;nb<nbins;nb++){ } */}
      /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_ZtSWn_stage_0]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_ZtSWn_stage_1(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  int n_spacing_A = SPACING_b;
  if (strstr(D->QC_strategy,"ZtSWn")){
    if (verbose){ printf(" %% [entering bcc_lrup_QC_ZtSWn_stage_1]\n");}
    if (verbose){ printf(" %% calculating E_[nb]->lf_ktrkn, E_[nb]->lf_ktrun.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
     if (E->A_rbother && D->A_cbother && D->Y_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_kt,E->M_rt,E->M_ut,D->A_ajdk,D->Y_ajdk,E->lf_ktrn,E->lf_utrn,&(E->lf_ktrun));
	GLOBAL_pthread_toc(); /* if bother */}
      if (E->A_rbother && D->A_cbother && D->A_cbother){ 
	GLOBAL_pthread_tic();
	wrap_AtTAn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,E->M_kt,E->M_rt,E->M_kt,D->A_ajdk,D->A_ajdk,E->lf_ktrn,E->lf_ktrn,&(E->lf_ktrkn));
	GLOBAL_pthread_toc(); /* if bother */}
      /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_pthread_tuc();
    GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose," %% ZtSWn_vv: ");
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	printf(" %% nb %d\n",nb);
	bprintf(D->A_bmc_j_rtn,D->bitj,1,D->A_ncols," %% A_bmc_j_rtn: ");
	bprintf(E->M_Yn->mc_j ,D->bitj,1,D->Y_ncols," %% bmc_Y_j    : ");
	bprintf(E->M_Tn->mc_j ,D->bitj,1,D->T_ncols," %% bmc_T_j    : ");
	if (E->A_rbother && D->A_cbother && D->Y_cbother){ lfprintf(E->lf_ktrun," %% lf_ktrun: "); /* if bother */} else{ printf(" %% no lf_ktrun\n");}
	if (E->A_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_ktrkn," %% lf_ktrkn: "); /* if bother */} else{ printf(" %% no lf_ktrkn\n");}
	/* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_ZtSWn_stage_1]\n");}
    /* if strategy */}
}

void bcc_lrup_QC_ZtSWn_stage_2(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0; struct bcc_single *E=NULL;
  if (strstr(D->QC_strategy,"ZtSWn")){
    if (verbose){ printf(" %% [entering bcc_lrup_QC_ZtSWn_stage_2]\n");}
    if (verbose){ printf(" %% updating E_[nb]->lf_AtTAn, E_[nb]->lf_AtTYn.\n");}
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
    for (nb=0;nb<nbins;nb++){ E = E_[nb];
     if (E->A_rbother && D->A_cbother && D->Y_cbother){ dra_subtequals(E->lf_AtTYn->lf,E->lf_AtTYn->length,E->lf_ktrun->lf); /* if bother */}
     if (E->A_rbother && D->A_cbother && D->A_cbother){ dra_subtequals(E->lf_AtTAn->lf,E->lf_AtTAn->length,E->lf_ktrkn->lf); /* if bother */}
     /* for (nb=0;nb<nbins;nb++){ } */}
    GLOBAL_ops_toc(-1,0,verbose," %% total time: ");
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    if (verbose>2){
      for (nb=0;nb<nbins;nb++){ E = E_[nb];
	printf(" %% nb %d\n",nb);
	bprintf(E->M_An->mc_j,D->bitj,1,D->A_ncols," %% bmc_A_j: ");
	bprintf(E->M_Yn->mc_j,D->bitj,1,D->Y_ncols," %% bmc_Y_j: ");
	bprintf(E->M_Tn->mc_j,D->bitj,1,D->T_ncols," %% bmc_T_j: ");
	if (E->A_rbother && D->A_cbother && D->Y_cbother){ lfprintf(E->lf_AtTYn," %% lf_AtTYn: "); /* if bother */} else{ printf(" %% no lf_AtTYn\n");}
	if (E->A_rbother && D->A_cbother && D->A_cbother){ lfprintf(E->lf_AtTAn," %% lf_AtTAn: "); /* if bother */} else{ printf(" %% no lf_AtTAn\n");}
	/* for (nb=0;nb<nbins;nb++){ } */}
    /* if (verbose>2){ } */}
    if (verbose){ printf(" %% [finished bcc_lrup_QC_ZtSWn_stage_2]\n");}
    /* if strategy */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_test_mxcut(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_single *E=NULL;
  int nb=0,nr=0,nc=0,cdrop=0,rdrop=0;
  double mrnd = GLOBAL_gamma;
  if (verbose){ printf(" %% [entering bcc_test_mxcut]\n");}
  for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc]=D->A_umc_j[nc]*(rand01>mrnd);}
  if (verbose>2){ raprintf(D->A_umc_j,"char",1,D->A_ncols," %% D->A_umc_j: ");}
  fill_uchar_zero(D->A_bmc_j,bsize(D->A_ncols));
  for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j,D->A_umc_j,nc);}
  cdrop = D->A_cpop_j - popcount_uchar_array(D->A_bmc_j,D->A_mc_length);
  if (verbose>0){ printf(" %% %d columns cut\n",cdrop);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); if (verbose>2){ bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  for (nc=0;nc<D->Y_ncols;nc++){ D->Y_umc_j[nc]=D->Y_umc_j[nc]*(1 || rand01>mrnd);} /* do not remove any columns */
  if (verbose>2){ raprintf(D->Y_umc_j,"char",1,D->Y_ncols," %% D->Y_umc_j: ");}
  fill_uchar_zero(D->Y_bmc_j,bsize(D->Y_ncols));
  for (nc=0;nc<D->Y_ncols;nc++){ b_copy_u(D->Y_bmc_j,D->Y_umc_j,nc);}
  sprintf(D->tmpYnchar," %%%% D->Y_bmc_j: "); if (verbose>2){ bprintf(D->Y_bmc_j,D->bitj,1,D->Y_ncols,D->tmpYnchar);}
  for (nc=0;nc<D->T_ncols;nc++){ D->T_umc_j[nc]=D->T_umc_j[nc]*(1 || nc==0 || rand01>mrnd);} /* do not remove any columns ; otherwise ensure that the first column of T,S is retained */
  if (verbose>2){ raprintf(D->T_umc_j,"char",1,D->T_ncols," %% D->T_umc_j: ");}
  fill_uchar_zero(D->T_bmc_j,bsize(D->T_ncols));
  for (nc=0;nc<D->T_ncols;nc++){ b_copy_u(D->T_bmc_j,D->T_umc_j,nc);}
  sprintf(D->tmpTnchar," %%%% D->T_bmc_j: "); if (verbose>2){ bprintf(D->T_bmc_j,D->bitj,1,D->T_ncols,D->tmpTnchar);}      
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr]=E->A_umr_j[nr]*(rand01>mrnd);}
    sprintf(D->tmpAnchar," %%%% E_[nb]->A_umr_j_[%.2d]: ",nb); if (verbose>2){ raprintf(E->A_umr_j,"char",1,E->A_nrows,D->tmpAnchar);}
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    fill_uchar_zero(E->A_bmr_j,bsize(E->A_nrows));
    for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j,E->A_umr_j,nr);}
    rdrop = E->A_rpop_j - popcount_uchar_array(E->A_bmr_j,E->A_mr_length);
    if (verbose>0){ printf(" %% %d rows cut\n",rdrop);}
    sprintf(D->tmpAnchar," %%%% E_[nb]->A_bmr_j_[%.2d]: ",nb); if (verbose>2){ bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    for (nr=0;nr<E->Z_nrows;nr++){ E->Z_umr_j[nr]=E->Z_umr_j[nr]*(1 || rand01>mrnd);} /* do not remove any rows */
    sprintf(D->tmpZnchar," %%%% E_[nb]->Z_umr_j_[%.2d]: ",nb); if (verbose>2){ raprintf(E->Z_umr_j,"char",1,E->Z_nrows,D->tmpZnchar);}
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    fill_uchar_zero(E->Z_bmr_j,bsize(E->Z_nrows));
    for (nr=0;nr<E->Z_nrows;nr++){ b_copy_u(E->Z_bmr_j,E->Z_umr_j,nr);}
    sprintf(D->tmpZnchar," %%%% E_[nb]->Z_bmr_j_[%.2d]: ",nb); if (verbose>2){ bprintf(E->Z_bmr_j,D->bitj,1,E->Z_nrows,D->tmpZnchar);}
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% [finished bcc_test_mxcut]\n");}
}

void bcc_lrup_mxcut(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_single *E=NULL;
  int nb=0,nr=0,nc=0,rdrop=0,rkeep=0,rdrop_max=0,cdrop=0,ckeep=0,cdrop_max=0,drop_flag=0;
  double mrnd = GLOBAL_gamma;
  if (verbose){ printf(" %% [entering bcc_lrup_mxcut]\n");}
  get_xdrop(D->A_rpop_j_total,D->A_cpop_j,&rdrop_max,&cdrop_max);
  for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);}
  cdrop=0; ckeep=0; for (nc=0;nc<D->A_ncols;nc++){ drop_flag = (cdrop>=cdrop_max || rand01>mrnd); D->A_umc_j_rtn[nc]=D->A_umc_j[nc]*drop_flag; D->A_umc_j_rmv[nc]=D->A_umc_j[nc]-D->A_umc_j_rtn[nc]; cdrop += D->A_umc_j_rmv[nc]; ckeep += D->A_umc_j_rtn[nc];}
  if (verbose>2){ raprintf(D->A_umc_j    ,"char",1,D->A_ncols," %% D->A_umc_j    : ");}
  if (verbose>2){ raprintf(D->A_umc_j_rmv,"char",1,D->A_ncols," %% D->A_umc_j_rmv: ");}
  if (verbose>2){ raprintf(D->A_umc_j_rtn,"char",1,D->A_ncols," %% D->A_umc_j_rtn: ");}
  if (verbose>0){ printf(" %% D->A_cpop_j %d ckeep %d cdrop %d\n",D->A_cpop_j,ckeep,cdrop);}
  fill_uchar_zero(D->A_bmc_j,bsize(D->A_ncols)); for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j,D->A_umc_j,nc);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j    : "); if (verbose>2){ bprintf(D->A_bmc_j    ,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  fill_uchar_zero(D->A_bmc_j_rmv,bsize(D->A_ncols)); for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j_rmv,D->A_umc_j_rmv,nc);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rmv: "); if (verbose>2){ bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  fill_uchar_zero(D->A_bmc_j_rtn,bsize(D->A_ncols)); for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j_rtn,D->A_umc_j_rtn,nc);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rtn: "); if (verbose>2){ bprintf(D->A_bmc_j_rtn,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);}
    rdrop=0; rkeep=0; for (nr=0;nr<E->A_nrows;nr++){ drop_flag = (rdrop>=rdrop_max || rand01>mrnd); E->A_umr_j_rtn[nr] = E->A_umr_j[nr]*drop_flag; E->A_umr_j_rmv[nr] = E->A_umr_j[nr]-E->A_umr_j_rtn[nr]; rdrop += E->A_umr_j_rmv[nr]; rkeep += E->A_umr_j_rtn[nr];}
    sprintf(D->tmpAnchar," %%%% E_[nb]->A_umr_j_[%.2d]    : ",nb); if (verbose>2){ raprintf(E->A_umr_j    ,"char",1,E->A_nrows,D->tmpAnchar);}
    sprintf(D->tmpAnchar," %%%% E_[nb]->A_umr_j_rmv_[%.2d]: ",nb); if (verbose>2){ raprintf(E->A_umr_j_rmv,"char",1,E->A_nrows,D->tmpAnchar);}
    sprintf(D->tmpAnchar," %%%% E_[nb]->A_umr_j_rtn_[%.2d]: ",nb); if (verbose>2){ raprintf(E->A_umr_j_rtn,"char",1,E->A_nrows,D->tmpAnchar);}
    if (verbose>0){ printf(" %% E->A_rpop_j %d rkeep %d rdrop %d\n",E->A_rpop_j,rkeep,rdrop);}
    fill_uchar_zero(E->A_bmr_j    ,bsize(E->A_nrows)); for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j    ,E->A_umr_j,nr);}
    sprintf(D->tmpAnchar," %%%% E_[nb]->A_bmr_j_    [%.2d]: ",nb); if (verbose>2){ bprintf(E->A_bmr_j    ,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    fill_uchar_zero(E->A_bmr_j_rmv,bsize(E->A_nrows)); for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j_rmv,E->A_umr_j_rmv,nr);}
    sprintf(D->tmpAnchar," %%%% E_[nb]->A_bmr_j_rmv_[%.2d]: ",nb); if (verbose>2){ bprintf(E->A_bmr_j_rmv,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    fill_uchar_zero(E->A_bmr_j_rtn,bsize(E->A_nrows)); for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j_rtn,E->A_umr_j_rtn,nr);}
    sprintf(D->tmpAnchar," %%%% E_[nb]->A_bmr_j_rtn_[%.2d]: ",nb); if (verbose>2){ bprintf(E->A_bmr_j_rtn,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% [finished bcc_lrup_mxcut]\n");}
}

void bcc_lrup_mxset(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_single *E=NULL;
  int nb=0;
  if (verbose){ printf(" %% [entering bcc_jvkur_mxset]\n");}
  for (nb=0;nb<nbins;nb++){ E = E_[nb];
    M_mxset(E->M_an,E->A_bmr_j_rtn,D->A_bmc_j_rtn); M_mxset(E->M_at,D->A_bmc_j_rtn,E->A_bmr_j_rtn);
    M_mxset(E->M_jn,E->A_bmr_j_rtn,D->A_bmc_j_rmv); M_mxset(E->M_jt,D->A_bmc_j_rmv,E->A_bmr_j_rtn);
    M_mxset(E->M_kn,E->A_bmr_j_rmv,D->A_bmc_j_rtn); M_mxset(E->M_kt,D->A_bmc_j_rtn,E->A_bmr_j_rmv);
    M_mxset(E->M_hn,E->A_bmr_j_rmv,D->A_bmc_j_rmv); M_mxset(E->M_ht,D->A_bmc_j_rmv,E->A_bmr_j_rmv);
    M_mxset(E->M_yn,E->A_bmr_j_rtn,D->Y_bmc_j    ); M_mxset(E->M_yt,D->Y_bmc_j    ,E->A_bmr_j_rtn);
    M_mxset(E->M_un,E->A_bmr_j_rmv,D->Y_bmc_j    ); M_mxset(E->M_ut,D->Y_bmc_j    ,E->A_bmr_j_rmv);
    M_mxset(E->M_zn,E->Z_bmr_j    ,D->A_bmc_j_rtn); M_mxset(E->M_zt,D->A_bmc_j_rtn,E->Z_bmr_j    );
    M_mxset(E->M_vn,E->Z_bmr_j    ,D->A_bmc_j_rmv); M_mxset(E->M_vt,D->A_bmc_j_rmv,E->Z_bmr_j    );
    M_mxset(E->M_tt,D->T_bmc_j    ,E->A_bmr_j_rtn);
    M_mxset(E->M_rt,D->T_bmc_j    ,E->A_bmr_j_rmv);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% [finished bcc_jvkur_mxset]\n");}
}

void bcc_lrup_mxdup(struct bcc_ajdk *D)
{
  int verbose=GLOBAL_verbose;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_single *E=NULL;
  int nb=0;
  ura_dup(D->A_bmc_j,D->A_mc_length,D->A_bmc_j_rtn); ura_dup(D->A_umc_j,D->A_ncols,D->A_umc_j_rtn);
  for (nb=0;nb<nbins;nb++){ E = E_[nb]; ura_dup(E->A_bmr_j,E->A_mr_length,E->A_bmr_j_rtn); ura_dup(E->A_umr_j,E->A_nrows,E->A_umr_j_rtn); /* for (nb=0;nb<nbins;nb++){ } */}
}

void bcc_lrup_test()
{
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0); int iteration_max = GLOBAL_TEST_niter;
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  int nl=0;
  bcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D);
  bcc_An_ajdk(D);
  bcc_lf_ZtSn(D);
  bcc_lf_ZtSWn(D);
  bcc_lf_YnWt(D);
  bcc_M_ZtSWn(D);
  bcc_M_YnWt(D);
  bcc_lf_AtTYn____WtSZn(D);
  bcc_lf_At_T_YnWt_S_Zn(D);
  if (error_check){ bcc_lf_AtTYnWtSZn(D);}
  bcc_lf_AnZt_S_WnYt(D);
  bcc_lf_An_ZtSWn_Yt(D);
  if (error_check){ bcc_lf_AnZtSWnYt(D);}
  if (error_check){ bcc_lf_AtTYnWtSZn_error(verbose,D);}
  if (error_check){ bcc_lf_AnZtSWnYt_error(verbose,D);}
  for (nl=0;nl<iteration_max;nl++){
    if (verbose>0){ printf(" %% iteration %d/%d\n",nl,iteration_max);}
    bcc_lrup_mxcut(D);
    bcc_lrup_mxset(D);
    if (strstr(D->QR_strategy,"condense")){ bcc_lrup_QR_YnWt_stage_2(D); /* if strategy */}
    else{ /* use AnZt_vv */ bcc_lrup_QR_YnWt_stage_0(D); bcc_lrup_QR_YnWt_stage_1(D); /* if strategy */}
    bcc_lrup_QR_YnWt_stage_3(D);
    bcc_lrup_QC_ZtSWn_stage_0(D);
    bcc_lrup_QC_ZtSWn_stage_1(D);
    bcc_lrup_QC_ZtSWn_stage_2(D);
    bcc_lrup_QC_YnWt_stage_a0(D);
    bcc_lrup_QC_YnWt_stage_a1(D);
    bcc_lrup_QC_YnWt_stage_a2(D);
    bcc_lrup_QC_YnWt_stage_a3(D);
    bcc_lrup_QC_YnWt_stage_b0(D);
    if (strstr(D->QC_strategy,"store one")){ bcc_lrup_QC_YnWt_stage_b3(D); /* if strategy */}
    else{ /* store all */ bcc_lrup_QC_YnWt_stage_b1(D); bcc_lrup_QC_YnWt_stage_b2(D); /* if strategy */}
    bcc_lrup_mxdup(D);
    bcc_M_mxset(D);
    bcc_M_ZtSWn(D);
    bcc_M_YnWt(D);
    bcc_lf_AtTYn____WtSZn(D);
    /* bcc_lf_At_T_YnWt_S_Zn(D); */
    bcc_lrup_QC_YnWt_stage_c(D);
    if (error_check){ bcc_lf_AtTYnWtSZn(D);}
    bcc_lf_AnZt_S_WnYt(D);
    /* fix later */ /* bcc_lf_An_ZtSWn_Yt(D); */
    if (error_check){ bcc_lf_AnZtSWnYt(D);}
    if (error_check){ bcc_lf_AtTYnWtSZn_error(verbose,D);}
    if (error_check){ bcc_lf_AnZtSWnYt_error(verbose,D);}
    /* for (nl=0;nl<iteration_max;nl++){ } */}  
  wkspace_printf();
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void bcc_lrup_flattenloop_test()
{
  /* tests a combination of lrup and flattenloop */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0); int iteration_max = GLOBAL_TEST_niter;
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  int nl=0;
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
  for (nl=0;nl<iteration_max;nl++){
    if (verbose>0){ printf(" %% iteration %d/%d\n",nl,iteration_max);}
    bcc_lrup_mxcut(D);
    bcc_lrup_mxset(D);
    if (strstr(D->QR_strategy,"condense")){ bcc_lrup_QR_YnWt_stage_2(D); /* if strategy */}
    else{ /* use AnZt_vv */ bcc_lrup_QR_YnWt_stage_0(D); bcc_lrup_QR_YnWt_stage_1(D); /* if strategy */}
    bcc_lrup_QR_YnWt_stage_3(D);
    bcc_lrup_QC_ZtSWn_stage_0(D);
    bcc_lrup_QC_ZtSWn_stage_1(D);
    bcc_lrup_QC_ZtSWn_stage_2(D);
    bcc_lrup_QC_YnWt_stage_a0(D);
    bcc_lrup_QC_YnWt_stage_a1(D);
    bcc_lrup_QC_YnWt_stage_a2(D);
    bcc_lrup_QC_YnWt_stage_a3(D);
    bcc_lrup_QC_YnWt_stage_b0(D);
    if (strstr(D->QC_strategy,"store one")){ bcc_lrup_QC_YnWt_stage_b3(D); /* if strategy */}
    else{ /* store all */ bcc_lrup_QC_YnWt_stage_b1(D); bcc_lrup_QC_YnWt_stage_b2(D); /* if strategy */}
    bcc_lrup_mxdup(D);
    bcc_M_mxset(D);
    bcc_M_ZtSWn(D);
    bcc_M_YnWt(D);
    bcc_lf_AtTYn____WtSZn(D);
    /* bcc_lf_At_T_YnWt_S_Zn(D); */
    bcc_lrup_QC_YnWt_stage_c(D);
    if (error_check){ bcc_QC_AtTYnWtSZn_uu(D);}
    bcc_lf_AnZt_S_WnYt(D);
    /* fix later */ /* bcc_lf_An_ZtSWn_Yt(D); */
    if (error_check){ bcc_QR_AnZtSWnYt_uu(D);}
    bcc_An_ajdk(D);
    bcc_lf_ZtSn(D);
    bcc_singlestudy_ww(D);
    bcc_doublestudy_ww(D);
    bcc_flattenloop(D);
    if (error_check){ bcc_QC_AtTYnWtSZn_error(verbose,D);}
    if (error_check){ bcc_QR_AnZtSWnYt_error(verbose,D);}
    /* for (nl=0;nl<iteration_max;nl++){ } */}  
  wkspace_printf();
}
