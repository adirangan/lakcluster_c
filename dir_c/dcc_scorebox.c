
void dcc_scorebox_mxA(struct dcc_ajdk *D,int rdrop,int cdrop)
{
  /* reading rdrop and cdrop from input: ;
     copies D->A_bmc_j to D->A_umc_j ; 
     define D->A_umc_j_rmv and D->A_umc_j_rtn ; removing first cdrop entries of D->QC_lmc_a ;
     copies E->A_bmr_j to E->A_umr_j ; 
     define E->A_umr_j_rmv and E->A_umr_j_rtn ; removing first rdrop entries of D->QR_lmr_a ; using D->QR_lnb to index nb ;
     call dcc_sumscores_mxB ; copying all instances of umc and umr to associated bmc and bmr respectively ;
   */
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_single *E=NULL;
  int nb1=0,nr=0,nc=0;
  if (verbose){ printf(" %% [entering dcc_scorebox_mxA]\n");}
  for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = D->A_umc_j[nc]; D->A_umc_j_rmv[nc] = 0;}
  for (nc=0;nc<minimum(cdrop,D->A_cpop_j);nc++){ D->A_umc_j_rmv[D->QC_lmc_a[nc]]=1; D->A_umc_j_rtn[D->QC_lmc_a[nc]]=0;}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = E->A_umr_j[nr]; E->A_umr_j_rmv[nr] = 0;}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nr=0;nr<minimum(rdrop,D->A_rpop_j_total);nr++){ E_[D->QR_lnb[nr]]->A_umr_j_rmv[D->QR_lmr_a[nr]]=1; E_[D->QR_lnb[nr]]->A_umr_j_rtn[D->QR_lmr_a[nr]]=0;}
  if (verbose>1){
    if (verbose){ raprintf(D->A_umc_j    ,"char",1,D->A_ncols," %% D->A_umc_j    : ");}
    if (verbose){ raprintf(D->A_umc_j_rmv,"char",1,D->A_ncols," %% D->A_umc_j_rmv: ");}
    if (verbose){ raprintf(D->A_umc_j_rtn,"char",1,D->A_ncols," %% D->A_umc_j_rtn: ");}
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
  if (verbose){ printf(" %% [finished dcc_scorebox_mxA]\n");}
}

void dcc_scorebox_sra(struct dcc_ajdk *D)
{
  /* Calculates D->QR_sra and D->QC_sra from E->QR_TAnAtT_nrm and D->QC_TAnAtT_nrm, respectively. */
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL;
  int na_a=0,na_b=0,na_j=0,ma_a=0,ma_b=0,ma_j=0;
  int Irow=0,Icol=0;
  if (verbose){ printf(" %% [entering dcc_scorebox_sra]\n");}
  D->Irem=0; for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; D->Irem += (E->A_rpop_j>0?1:0);} 
  D->Ireq = minimum(nbins,minimum(D->Irem,GLOBAL_Ireq));
  if (GLOBAL_Ireq<=0){ Irow=0; Icol=0;} else if (GLOBAL_Ireq>0){ Irow = D->Irem - D->Ireq; Icol = D->Irem*D->Irem - D->Ireq*D->Ireq;}
  if (verbose){ printf(" %% found D->Irem %d/%d, using D->Ireq %d; Irow %d Icol %d\n",D->Irem,nbins,D->Ireq,Irow,Icol);}
  if (verbose){
    printf(" %% D->A_ncols %d D->A_cpop_b %d D->A_cpop_j %d\n",D->A_ncols,D->A_cpop_b,D->A_cpop_j);
    if (verbose>1){ sprintf(D->tmpAnchar," %%%% D->A_bmc_b: "); bprintf(D->A_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
    if (verbose>1){ sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
    /* if (verbose){ } */}
  if (D->A_cbother && D->A_rpop_b_total){
    na_a=0;na_b=0;na_j=0;
    while (na_a<D->A_ncols){
      if (bget__on(D->A_bmc_b,na_a)){
	if (bget__on(D->A_bmc_j,na_a)){
	  D->QC_sra[na_j] = D->QC_TAnAtT_nrm[na_a + 0*D->A_ncols + Icol*D->A_ncols*D->T_ncols];
	  na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	na_b++; /* if (bget__on(D->A_bmc_b,na_a)){ } */}
      na_a++;/* while (na_a<D->A_ncols){ } */}
    if (verbose>1){
      raprintf(  D->QC_sra,"double",1,D->A_cpop_j," %%   D->QC_sra :");
      /* if (verbose){ } */}
    /* if (D->A_cbother && D->A_rpop_b_total){ } */}
  if (D->A_rpop_b_total && D->A_cbother){
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
	    D->QR_sra[D->A_rpop_j_total] = E->QR_TAnAtT_nrm[ma_a + 0*E->A_nrows + Irow*E->A_nrows*D->T_ncols];
	    ma_j++; D->A_rpop_j_total++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	  ma_b++; D->A_rpop_b_total++; /* while (ma_a<E->A_nrows){ } */}
	ma_a++; D->A_nrows_total++; /* while (ma_a<E->A_nrows){ } */}
      /* for (nb1=0;nb1<nbins;nb1++){ } */}
    if (verbose>1){
      raprintf(  D->QR_sra,"double",1,D->A_rpop_j_total," %% D->QR_sra :");
      /* if (verbose){ } */}
    /* if (D->A_rpop_b_total && D->A_cbother){ } */}
  if (verbose){ printf(" %% [finished dcc_scorebox_sra]\n");}
}

void dcc_scorebox_srt(struct dcc_ajdk *D,int nrows,int *mr_lnb,int *mr_lmr,int ncols,int *mc_srt)
{
  /* Copies last D->A_rpop_j_total entries from mr_srt (i.e., mr_lnb and mr_lmr) ;
     into D->QR_lnb and D->QR_lmr_a, respectively. ;
     Copies last D->A_cpop_j entries from mc_srt ;
     into D->QC_lmc_a. ;
  */
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL;
  int na_j=0,nr=0;
  if (verbose){ printf(" %% [entering dcc_scorebox_srt]\n");}
  if (mc_srt!=NULL && ncols>=D->A_cpop_j){ 
    for (na_j=0;na_j<D->A_cpop_j;na_j++){ D->QC_lmc_a[na_j] = mc_srt[ncols-D->A_cpop_j+na_j];} 
    /* if (mc_srt!=NULL){ } */}
  if (mr_lnb!=NULL && mr_lmr!=NULL && nrows>=D->A_rpop_j_total){ 
    for (nr=0;nr<D->A_rpop_j_total;nr++){ 
      D->QR_lmr_a[nr] = mr_lmr[nrows-D->A_rpop_j_total+nr];
      D->QR_lnb[nr] = mr_lnb[nrows-D->A_rpop_j_total+nr];
      /* for (nr=0;nr<D->A_rpop_j_total;nr++){ } */}
    /* if (mc_srt!=NULL){ } */};
  if (verbose>1){ printf(" %% finished copying mr_lnb,mr_lmr and mc_srt\n");}    
  if (verbose>1){
    raprintf(D->QR_lmr_a,   "int",1,D->A_rpop_j_total," %% D->QR_lmr_a pos :");
    raprintf(  D->QR_lnb,   "int",1,D->A_rpop_j_total," %%   D->QR_lnb pos :");
    raprintf(D->QC_lmc_a,   "int",1,D->A_cpop_j," %% D->QC_lmc_a pos :");
    /* if (verbose){ } */}
  if (verbose){ printf(" %% [finished dcc_scorebox_srt]\n");}
}

void dcc_out_xdrop_lkp(struct dcc_ajdk *D,int nrows,int *mr_srt,int **mr_lnb_,int **mr_lmr_)
{
  /* Converts mr_srt to lnb (local nb) and lmr (local mr) indices. ;
     Note that all indices are of 'a' type (i.e., referring to mr_a). */
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nl=0,tmpd=0,nb=0,mx_a_tot=0;
  struct dcc_single *E=NULL;
  if (verbose){ printf(" %% [entering dcc_out_xdrop_lkp]\n");}
  mx_a_tot = 0; for (nb=0;nb<nbins;nb++){ E = E_[nb]; mx_a_tot += E->A_nrows; /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% nbins %d, mx_a_tot %d, nrows %d\n",nbins,mx_a_tot,nrows);}
  if (*mr_lnb_==NULL){ (*mr_lnb_) = (int *) wkspace_all0c(nrows*sizeof(int));}
  if (*mr_lmr_==NULL){ (*mr_lmr_) = (int *) wkspace_all0c(nrows*sizeof(int));}
  for (nl=0;nl<nrows;nl++){  
    tmpd = mr_srt[nl];
    nb=0; while (tmpd>=E_[nb]->A_nrows){ tmpd-=E_[nb]->A_nrows; nb++;}
    (*mr_lnb_)[nl] = nb; (*mr_lmr_)[nl] = tmpd;
    /* for (nl=0;nl<nrows;nl++){ } */}
  if (verbose>1){ raprintf(*mr_lnb_,"int",1,nrows," %% mr_lnb: ");}
  if (verbose>1){ raprintf(*mr_lmr_,"int",1,nrows," %% mr_lmr: ");}
  if (verbose){ printf(" %% [finished dcc_out_xdrop_lkp]\n");}
}

