void dcc_sumscores_mxB(struct dcc_ajdk *D)
{
  /* writes D->A_umc_j into D->A_bmc_j  */
  /* writes D->A_umc_j_rmv into D->A_bmc_j_rmv  */
  /* writes D->A_umc_j_rtn into D->A_bmc_j_rtn */
  /* writes E->A_umr_j into E->A_bmr_j  */
  /* writes E->A_umr_j_rmv into E->A_bmr_j_rmv  */
  /* writes E->A_umr_j_rtn into E->A_bmr_j_rtn */
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0,nr=0,nc=0; struct dcc_single *E=NULL;
  if (verbose){ raprintf(D->A_umc_j,"char",1,D->A_ncols," %% D->A_umc_j: ");}
  if (verbose){ raprintf(D->A_umc_j_rmv,"char",1,D->A_ncols," %% D->A_umc_j_rmv: ");}
  if (verbose){ raprintf(D->A_umc_j_rtn,"char",1,D->A_ncols," %% D->A_umc_j_rtn: ");}
  fill_uchar_zero(D->A_bmc_j,bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */); for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j,D->A_umc_j,nc);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); if (verbose){ bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  fill_uchar_zero(D->A_bmc_j_rmv,bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */); for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j_rmv,D->A_umc_j_rmv,nc);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rmv: "); if (verbose){ bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  fill_uchar_zero(D->A_bmc_j_rtn,bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */); for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j_rtn,D->A_umc_j_rtn,nc);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rtn: "); if (verbose){ bprintf(D->A_bmc_j_rtn,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j: ",nb1); if (verbose){ raprintf(E->A_umr_j,"char",1,E->A_nrows,D->tmpAnchar);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j_rmv: ",nb1); if (verbose){ raprintf(E->A_umr_j_rmv,"char",1,E->A_nrows,D->tmpAnchar);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j_rtn: ",nb1); if (verbose){ raprintf(E->A_umr_j_rtn,"char",1,E->A_nrows,D->tmpAnchar);}
    fill_uchar_zero(E->A_bmr_j,bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */); for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j,E->A_umr_j,nr);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j: ",nb1); if (verbose){ bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    fill_uchar_zero(E->A_bmr_j_rmv,bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */); for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j_rmv,E->A_umr_j_rmv,nr);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j_rmv: ",nb1); if (verbose){ bprintf(E->A_bmr_j_rmv,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    fill_uchar_zero(E->A_bmr_j_rtn,bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */); for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j_rtn,E->A_umr_j_rtn,nr);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j_rtn: ",nb1); if (verbose){ bprintf(E->A_bmr_j_rtn,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  D->A_cpop_j_rmv = popcount_uchar_array(D->A_bmc_j_rmv,D->A_mc_length);
  D->A_cpop_j_rtn = popcount_uchar_array(D->A_bmc_j_rtn,D->A_mc_length);
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    E->A_rpop_j_rmv = popcount_uchar_array(E->A_bmr_j_rmv,E->A_mr_length);
    E->A_rpop_j_rtn = popcount_uchar_array(E->A_bmr_j_rtn,E->A_mr_length);
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% [finished dcc_sumscores_mxB]\n");}
}

void dcc_sumscores_mxC(struct dcc_ajdk *D)
{
  /* writes data to out_trace ; 
     0: out_iteration ;
     1: A_rpop_j_total ; rows remaining
     2: A_cpop_j ; columns remaining
     3: average of QR_sra ; average row-score
     4: average of QC_sra ; average col-score
     5: nb_rem ; covariate categories remaining
  */
  int verbose=0;
  double tmp_mean=0;
  /* 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MA = nrows_rem; NA = ncols_rem; MZ = size(Z_orig,1); NZ = ncols_rem;
    tmp_R = sum(QR); R_d = log(MA) + log(NA) + log(MA) + log(NA);
    tmp_C = sum(QC); C_d = log(NA) + log(MA) + log(NA) + log(MA);
    out_trace(iteration,:) = [iteration , MA , NA , sign(tmp_R)*exp(log(abs(tmp_R))-R_d) , sign(tmp_C)*exp(log(abs(tmp_C))-C_d) sum(cov_bother_flag_AAAA)+sum(cov_bother_flag_AZZA)];
    for nr=1:length(rdrop);
    out_xdrop(out_xdrop_ij,:) = [A_n_rind_vals_lookup(rdrop(end+1-nr))-1 , -1];
    out_xdrop_ij = out_xdrop_ij+1;
    end;%for nr=1:length(rdrop);
    for nc=1:length(cdrop);
    out_xdrop(out_xdrop_ij,:) = [-1 , A_n_cind_vals_lookup(cdrop(end+1-nc))-1];
    out_xdrop_ij = out_xdrop_ij+1;
    end;%for nc=1:length(cdrop);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */
  D->out_trace[0 + D->out_iteration*D->out_trace_length] = D->out_iteration;
  D->out_trace[1 + D->out_iteration*D->out_trace_length] = D->A_rpop_j_total;
  D->out_trace[2 + D->out_iteration*D->out_trace_length] = D->A_cpop_j;
  ra_stats(D->QR_sra,"double",D->A_rpop_j_total,NULL,NULL,&tmp_mean,NULL); if (!isfinite(tmp_mean)){ tmp_mean = 1.0;}
  D->out_trace[3 + D->out_iteration*D->out_trace_length] = tmp_mean; 
  ra_stats(D->QC_sra,"double",D->A_cpop_j,NULL,NULL,&tmp_mean,NULL); if (!isfinite(tmp_mean)){ tmp_mean = 1.0;}
  D->out_trace[4 + D->out_iteration*D->out_trace_length] = tmp_mean;
  D->out_trace[5 + D->out_iteration*D->out_trace_length] = D->Irem;
  D->out_iteration++;
  if (verbose){ printf(" %% [finished dcc_sumscores_mxC]\n");}
}

void dcc_sumscores_dmp(struct dcc_ajdk *D)
{
  /* dump temporary variables used within dcc to disc */
  /* Note that out_trace holds:
     0: out_iteration ;
     1: A_rpop_j_total ; rows remaining
     2: A_cpop_j ; columns remaining
     3: average of QR_sra ; average row-score
     4: average of QC_sra ; average col-score
     5: nb_rem ; covariate categories remaining
  */
  int verbose=0;
  int nl=0; FILE *fp=NULL; char tempchar[FNAMESIZE];
  FILE *fp_a=NULL; char tempchar_a[FNAMESIZE];
  FILE *fp_b=NULL; char tempchar_b[FNAMESIZE];
  if (verbose){ printf(" %% [entering dcc_sumscores_dmp]\n");}
  sprintf(tempchar,"%s/out_trace.txt",GLOBAL_DIR_NAME); if ((fp=fopen(tempchar,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",tempchar); exit(RET_READ_FAIL);}
  for (nl=0;nl<D->out_iteration;nl++){
    fprintf(fp,"%d %d %d %0.16f %0.16f %d\n",(int)D->out_trace[0 + nl*D->out_trace_length],(int)D->out_trace[1 + nl*D->out_trace_length],(int)D->out_trace[2 + nl*D->out_trace_length],D->out_trace[3 + nl*D->out_trace_length],D->out_trace[4 + nl*D->out_trace_length],(int)D->out_trace[5 + nl*D->out_trace_length]);
    /* for (nl=0;nl<D->out_iteration;nl++){ } */}
  fclose(fp);fp=NULL; 
  sprintf(tempchar_a,"%s/out_xdrop_a.txt",GLOBAL_DIR_NAME); if ((fp_a=fopen(tempchar_a,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",tempchar_a); exit(RET_READ_FAIL);}
  sprintf(tempchar_b,"%s/out_xdrop_b.txt",GLOBAL_DIR_NAME); if ((fp_b=fopen(tempchar_b,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",tempchar_b); exit(RET_READ_FAIL);}
  for (nl=0;nl<D->out_xdrop_ij;nl++){
    fprintf(fp_a,"%d %d\n",(int)D->out_xdrop_a[0 + nl*2],(int)D->out_xdrop_a[1 + nl*2]);
    fprintf(fp_b,"%d %d\n",(int)D->out_xdrop_b[0 + nl*2],(int)D->out_xdrop_b[1 + nl*2]);
    /* for (nl=0;nl<D->out_xdrop_ij;nl++){ } */}
  fclose(fp_a);fp_a=NULL; 
  fclose(fp_b);fp_b=NULL; 
  if (verbose){ printf(" %% [finished dcc_sumscores_dmp]\n");}
}

void dcc_sumscores_mxA(struct dcc_ajdk *D)
{
  /* obtains rdrop and cdrop from get_xdrop ;
     copies D->A_bmc_j to D->A_umc_j ; 
     define D->A_umc_j_rmv and D->A_umc_j_rtn ; removing first cdrop entries of D->QC_lmc_a ;
     copies E->A_bmr_j to E->A_umr_j ; 
     define E->A_umr_j_rmv and E->A_umr_j_rtn ; removing first rdrop entries of D->QR_lmr_a ; using D->QR_lnb to index nb ;
     call dcc_sumscores_mxB ; copying all instances of umc and umr to associated bmc and bmr respectively ;
     call dcc_sumscores_mxC ; writing data to out_trace ;
     define D->out_xdrop_a and D->out_xdrop_b ; dump row-indices first, then col-indices ;
   */
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_single *E=NULL;
  double lrij=D->A_rpop_j_total,lcij=D->A_cpop_j; 
  int rdrop=0;int cdrop=0;
  int nb1=0,nr=0,nc=0;
  if (verbose){ printf(" %% [entering dcc_sumscores_mxA]\n");}
  get_xdrop(D->A_rpop_j_total,D->A_cpop_j,&rdrop,&cdrop);
  for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = D->A_umc_j[nc]; D->A_umc_j_rmv[nc] = 0;}
  for (nc=0;nc<minimum(cdrop,lcij);nc++){ D->A_umc_j_rmv[D->QC_lmc_a[nc]]=1; D->A_umc_j_rtn[D->QC_lmc_a[nc]]=0;}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = E->A_umr_j[nr]; E->A_umr_j_rmv[nr] = 0;}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nr=0;nr<minimum(rdrop,lrij);nr++){ E_[D->QR_lnb[nr]]->A_umr_j_rmv[D->QR_lmr_a[nr]]=1; E_[D->QR_lnb[nr]]->A_umr_j_rtn[D->QR_lmr_a[nr]]=0;}
  dcc_sumscores_mxB(D); dcc_sumscores_mxC(D);
  if (verbose){ raprintf(D->out_trace,"double_trn",D->out_trace_length,D->out_iteration," %% D->out_trace: ");}
  for (nr=0;nr<minimum(rdrop,lrij);nr++){ 
    D->out_xdrop_a[0 + D->out_xdrop_ij*2] = D->QR_imr_a[nr]; D->out_xdrop_a[1 + D->out_xdrop_ij*2] = -1;
    D->out_xdrop_b[0 + D->out_xdrop_ij*2] = D->QR_imr_b[nr]; D->out_xdrop_b[1 + D->out_xdrop_ij*2] = -1;
    D->out_xdrop_ij++; /* for (nr=0;nr<minimum(rdrop,lrij);nr++){ } */}
  for (nc=0;nc<minimum(cdrop,lcij);nc++){ 
    D->out_xdrop_a[1 + D->out_xdrop_ij*2] = D->QC_lmc_a[nc]; D->out_xdrop_a[0 + D->out_xdrop_ij*2] = -1; 
    D->out_xdrop_b[1 + D->out_xdrop_ij*2] = D->QC_lmc_b[nc]; D->out_xdrop_b[0 + D->out_xdrop_ij*2] = -1; 
    D->out_xdrop_ij++; /* for (nc=0;nc<minimum(cdrop,lcij);nc++){ } */}
  if (verbose){ raprintf(D->out_xdrop_a,"int",2,D->out_xdrop_ij," %% D->out_xdrop_a: ");}
  if (verbose){ raprintf(D->out_xdrop_b,"int",2,D->out_xdrop_ij," %% D->out_xdrop_b: ");}
  if (verbose){ printf(" %% [finished dcc_sumscores_mxA]\n");}
}

void dcc_sumscores_xij(struct dcc_ajdk *D)
{
  /* sorts scores in increasing order */
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL;
  int na_a=0,na_b=0,na_j=0,ma_a=0,ma_b=0,ma_j=0;
  int Irow=0,Icol=0;
  unsigned int nn=0;
  if (verbose){ printf(" %% [entering dcc_sumscores_xij]\n");}
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
	  D->QC_lmc_a[na_j] = na_a; D->QC_lmc_b[na_j] = na_b; D->QC_lmc_j[na_j] = na_j;
	  na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	na_b++; /* if (bget__on(D->A_bmc_b,na_a)){ } */}
      na_a++;/* while (na_a<D->A_ncols){ } */}
    if (verbose>1){
      raprintf(  D->QC_sra,"double",1,D->A_cpop_j," %%   D->QC_sra pre :");
      raprintf(D->QC_lmc_a,   "int",1,D->A_cpop_j," %% D->QC_lmc_a pre :");
      raprintf(D->QC_lmc_b,   "int",1,D->A_cpop_j," %% D->QC_lmc_b pre :");
      raprintf(D->QC_lmc_j,   "int",1,D->A_cpop_j," %% D->QC_lmc_j pre :");
      /* if (verbose){ } */}
    nn = dQuickSort_xij(0,D->QC_sra,1,D->QC_lmc_a,D->QC_lmc_b,D->QC_lmc_j,NULL,NULL,NULL,0,D->A_cpop_j-1);
    if (verbose>1){ printf(" %% finished sorting QC_sra, maximum recursion_level %d\n",nn);}
    if (verbose>1){
      raprintf(  D->QC_sra,"double",1,D->A_cpop_j," %%   D->QC_sra pos :");
      raprintf(D->QC_lmc_a,   "int",1,D->A_cpop_j," %% D->QC_lmc_a pos :");
      raprintf(D->QC_lmc_b,   "int",1,D->A_cpop_j," %% D->QC_lmc_b pos :");
      raprintf(D->QC_lmc_j,   "int",1,D->A_cpop_j," %% D->QC_lmc_j pos :");
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
	    D->QR_lmr_a[D->A_rpop_j_total] = ma_a; D->QR_lmr_b[D->A_rpop_j_total] = ma_b; D->QR_lmr_j[D->A_rpop_j_total] = ma_j; D->QR_lnb[D->A_rpop_j_total] = nb1; 
	    D->QR_imr_a[D->A_rpop_j_total] = E->QR_imr_a[ma_a]; D->QR_imr_b[D->A_rpop_j_total] = E->QR_imr_b[ma_a];
	    ma_j++; D->A_rpop_j_total++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	  ma_b++; D->A_rpop_b_total++; /* while (ma_a<E->A_nrows){ } */}
	ma_a++; D->A_nrows_total++; /* while (ma_a<E->A_nrows){ } */}
      /* for (nb1=0;nb1<nbins;nb1++){ } */}
    if (verbose>1){
      raprintf(  D->QR_sra,"double",1,D->A_rpop_j_total," %%   D->QR_sra pre :");
      raprintf(D->QR_lmr_a,   "int",1,D->A_rpop_j_total," %% D->QR_lmr_a pre :");
      raprintf(D->QR_lmr_b,   "int",1,D->A_rpop_j_total," %% D->QR_lmr_b pre :");
      raprintf(D->QR_lmr_j,   "int",1,D->A_rpop_j_total," %% D->QR_lmr_j pre :");
      raprintf( D->QR_lnb,   "int",1,D->A_rpop_j_total," %%  D->QR_lnb pre :");
      raprintf( D->QR_imr_a,   "int",1,D->A_rpop_j_total," %%  D->QR_imr_a pre :");
      raprintf( D->QR_imr_b,   "int",1,D->A_rpop_j_total," %%  D->QR_imr_b pre :");
      /* if (verbose){ } */}
    nn = dQuickSort_xij(0,D->QR_sra,1,D->QR_lmr_a,D->QR_lmr_b,D->QR_lmr_j,D->QR_lnb,D->QR_imr_a,D->QR_imr_b,0,D->A_rpop_j_total-1);
    if (verbose>1){ printf(" %% finished sorting QR_sra, maximum recursion_level %d\n",nn);}
    if (verbose>1){
      raprintf(  D->QR_sra,"double",1,D->A_rpop_j_total," %%   D->QR_sra pos :");
      raprintf(D->QR_lmr_a,   "int",1,D->A_rpop_j_total," %% D->QR_lmr_a pos :");
      raprintf(D->QR_lmr_b,   "int",1,D->A_rpop_j_total," %% D->QR_lmr_b pos :");
      raprintf(D->QR_lmr_j,   "int",1,D->A_rpop_j_total," %% D->QR_lmr_j pos :");
      raprintf( D->QR_lnb,   "int",1,D->A_rpop_j_total," %%  D->QR_lnb pos :");
      raprintf( D->QR_imr_a,   "int",1,D->A_rpop_j_total," %%  D->QR_imr_a pos :");
      raprintf( D->QR_imr_b,   "int",1,D->A_rpop_j_total," %%  D->QR_imr_b pos :");
      /* if (verbose){ } */}
    /* if (D->A_rpop_b_total && D->A_cbother){ } */}
  if (verbose){ printf(" %% [finished dcc_sumscores_xij]\n");}
}

void dcc_sumscores_cmb(struct dcc_ajdk *D)
{
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0,nb2=0,nbx=0; struct dcc_single *E=NULL; unsigned char bmc1=255;
  if (verbose){ printf(" %% [entering dcc_sumscores_cmb]\n");}
  /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    QC = QC_AAAA_min - QC_AAZZ_min;
    for nb1=0:nbins-1;
    QR{1+nb1} = QR_AAAA{1+nb1} - QR_AZZA{1+nb1};
    end;%for nb1=0:nbins-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */
  if (verbose){ printf(" %% [entering dcc_sumscores_cmb]\n");}
  if (D->A_cbother && D->A_rpop_b_total){
    if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ for (nbx=0;nbx<nbins*nbins;nbx++){ dra_plustimesequals_s___m_m(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,1/*D->T_ncols*/,&(D->QC_TAnZtS_nrm[nbx*D->A_ncols*D->T_ncols]),-1.0,D->A_bmc_j,&bmc1/*D->T_bmc_j*/);}}
    if (verbose){ for (nbx=0;nbx<nbins*nbins;nbx++){ sprintf(D->tmpAnchar," %%%% D->QC_TAnAtT_nrm[%03d]: ",nbx); raprintf(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,1,D->tmpAnchar);}}
    /* if (D->A_cbother && D->A_rpop_b_total){ } */}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; 
    if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){
      if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ for (nb2=0;nb2<nbins;nb2++){ dra_plustimesequals_s___m_m(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,1/*D->T_ncols*/,&(E->QR_TAnZtS_nrm[nb2*E->A_nrows*D->T_ncols]),-1.0,E->A_bmr_j,&bmc1/*D->T_bmc_j*/);}}
      if (verbose){ for (nb2=0;nb2<nbins;nb2++){ sprintf(D->tmpAnchar," %%%% E_[%d]->QR_TAnAtT_nrm[%d]: ",nb1,nb2); raprintf(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,1,D->tmpAnchar);}}
      /* if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ } */}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% [finished dcc_sumscores_cmb]\n");}
}

void dcc_sumscores_srt(struct dcc_ajdk *D)
{
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0,nb2=0,nbx=0; struct dcc_single *E=NULL; 
  /* unsigned char bmc1=255; */
  int ma_j=0,ma_b=0,ma_a=0,na_j=0,na_b=0,na_a=0;
  if (verbose){ printf(" %% [entering dcc_sumscores_srt]\n");}
  /* 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     QC_AAAA_nrm = min(QC_AAAA_nrm,2);
     QC_AAZZ_nrm = min(QC_AAZZ_nrm,2);
     for nb1=0:nbins-1;
     QR_AAAA_nrm{1+nb1} = sort(QR_AAAA_nrm{1+nb1},2);
     QR_AZZA_nrm{1+nb1} = sort(QR_AZZA_nrm{1+nb1},2);
     end;%for nb1=0:nbins-1;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];     
    if (verbose){ printf(" %% nb1 %d\n",nb1);}
    if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ 
      ma_a=0; ma_b=0; ma_j=0; 
      while (ma_a<E->A_nrows){
	if (bget__on(E->A_bmr_b,ma_a)){
	  if (bget__on(E->A_bmr_j,ma_a)){
	    if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_TAnZtS_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	    dQuickSort(0,&(E->QR_TAnZtS_nrm[ma_a+0*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,0,nbins-1);   
	    if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_TAnZtS_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	    ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	  ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	ma_a++; /* while (ma_a<E->A_nrows){ } */}
      /* if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ } */}
    if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ 
      ma_a=0; ma_b=0; ma_j=0; 
      while (ma_a<E->A_nrows){
	if (bget__on(E->A_bmr_b,ma_a)){
	  if (bget__on(E->A_bmr_j,ma_a)){
	    if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_TAnAtT_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	    dQuickSort(0,&(E->QR_TAnAtT_nrm[ma_a+0*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,0,nbins-1);   
	    if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_TAnAtT_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	    ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	  ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	ma_a++; /* while (ma_a<E->A_nrows){ } */}
      /* if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ } */}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ 
    na_a=0; na_b=0; na_j=0; 
    while (na_a<D->A_ncols){
      if (bget__on(D->A_bmc_b,na_a)){
	if (bget__on(D->A_bmc_j,na_a)){
	  if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_TAnZtS_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	  dQuickSort(0,&(D->QC_TAnZtS_nrm[na_a+0*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,0,nbins*nbins-1);   
	  if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_TAnZtS_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	  na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
      na_a++; /* while (na_a<D->A_ncols){ } */}
    /* if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ } */}
  if (D->A_cbother && D->A_rpop_b_total){ 
    na_a=0; na_b=0; na_j=0; 
    while (na_a<D->A_ncols){
      if (bget__on(D->A_bmc_b,na_a)){
	if (bget__on(D->A_bmc_j,na_a)){
	  if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_TAnAtT_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	  dQuickSort(0,&(D->QC_TAnAtT_nrm[na_a+0*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,0,nbins*nbins-1);   
	  if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_TAnAtT_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	  na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
      na_a++; /* while (na_a<D->A_ncols){ } */}
    /* if (D->A_cbother && D->A_rpop_b_total){ } */}
  if (verbose){ printf(" %% [finished dcc_sumscores_srt]\n");}
}

void dcc_sumscores_ifT(struct dcc_ajdk *D)
{
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0,nb2=0,nbx=0; struct dcc_single *E=NULL; 
  int n_mds = maximum(0,D->T_cpop_j-1); int n_mds_max = 6; double mds_scale_factor=1.0;
  if (verbose){ printf(" %% [entering dcc_sumscores_ifT], n_mds %d\n",n_mds);}
  /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_mds = 5; mds_scale_factor = [1.00,0.30,0.12,0.10,0.09];
    QC_all = temp_D_At_T_AnAt_T_An{1+0}/((A_nrows_(1+0) - 1)*A_ncols*(A_ncols-1));
    QR_all = temp_AnAt_T_AnAt{1+0}/(A_nrows_(1+0)*(A_nrows_(1+0) - 1)*A_ncols);
    QC = QC_all(:,1) - sum(QC_all(:,2:end),2)/n_mds/mds_scale_factor(n_mds);
    QR = (QR_all(:,1).^2 - sum(QR_all(:,2:end).^2,2)/n_mds/mds_scale_factor(n_mds));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */
  if (n_mds>n_mds_max){ printf(" %% Warning! n_mds>n_mds_max %d\n",n_mds_max); n_mds = n_mds_max;} 
  if (n_mds<=0){
    if (verbose){ printf(" %% n_mds==0, do nothing\n");}
    /* if (n_mds<=0){ } */}
  else /* if (n_mds>0) */{
    mds_scale_factor = GLOBAL_kappa_squared_loop_scale_factor_[n_mds-1];
    if (GLOBAL_kappa_squared>0);{ mds_scale_factor = GLOBAL_kappa_squared;}
    if (verbose){ printf(" %% n_mds==%d, mds_scale_factor=%0.3f\n",n_mds,mds_scale_factor);}
    for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; 
      for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins; 
	if (verbose){ printf(" %% nb1 %d nb2 %d nbx %d\n",nb1,nb2,nbx);}
	if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ 
	  dra_mds_pow_s___m_m(&(E->QR_TAnZtS_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,n_mds*mds_scale_factor,E->A_bmr_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(E->QR_TAnZtS_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_TAnZtS_nrm: ");}
	  /* if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ } */}
	if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ 
	  dra_mds_pow_s___m_m(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,n_mds*mds_scale_factor,E->A_bmr_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_TAnAtT_nrm: ");}
	  /* if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ } */}
	if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ 
	  dra_mds_nrm_s___m_m(&(D->QC_TAnZtS_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,n_mds*mds_scale_factor,D->A_bmc_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(D->QC_TAnZtS_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_TAnZtS_nrm: ");}
	  /* if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ } */}
	if (D->A_cbother && D->A_rpop_b_total){ 
	  dra_mds_nrm_s___m_m(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,n_mds*mds_scale_factor,D->A_bmc_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_TAnAtT_nrm: ");}
	  /* if (D->A_cbother && D->A_rpop_b_total){ } */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}    
    /* if (n_mds>0){ } */}
  if (verbose){ printf(" %% [finished dcc_sumscores_ifT]\n");}
}

void dcc_sumscores_nrm(struct dcc_ajdk *D)
{
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; long long int lld=0; struct dcc_single *E=NULL; struct dcc_double *F=NULL;
  if (verbose){ printf(" %% [entering dcc_sumscores_nrm]\n");}
  for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins;
      if (verbose){ printf(" %% nb1 %d nb2 %d nbx %d\n",nb1,nb2,nbx);}
      E = E_[nb1]; F = F_[nbx];
      if (verbose>3){
	lfprintf(F->QR_TAnZtS," %% F->QR_TAnZtS: ");
	lfprintf(F->QR_TAnAtT," %% F->QR_TAnAtT: ");
	lfprintf(F->QC_TAnZtS," %% F->QC_TAnZtS: ");
	lfprintf(F->QC_TAnAtT," %% F->QC_TAnAtT: ");
	/* if (verbose>3){ } */}
      if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->A_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_TAnZtS_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	lld = (long long int)(E_[nb2]->M_Zn->rpop_j - 0*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_Zn->cpop_j);
	dra_plusdivequals_s___m_m(&(E->QR_TAnZtS_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_TAnZtS->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	if (verbose){ raprintf(&(E->QR_TAnZtS_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_TAnZtS_nrm: ");}
	/* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->A_cbother){ } */}
      if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->A_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	lld = (long long int)(E_[nb2]->M_An->rpop_j - 1*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_An->cpop_j);
	dra_plusdivequals_s___m_m(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_TAnAtT->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	if (verbose){ raprintf(&(E->QR_TAnAtT_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_TAnAtT_nrm: ");}
	/* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->A_cbother){ } */}
      if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->Z_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_TAnZtS_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	lld = (long long int)(E_[nb1]->M_An->rpop_j) * (long long int)(E_[nb2]->M_Zn->rpop_j);
	dra_plusdivequals_s___m_m(&(D->QC_TAnZtS_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_TAnZtS->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	if (verbose){ raprintf(&(D->QC_TAnZtS_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_TAnZtS_nrm: ");}
	/* if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->Z_rbother){ } */}
      if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->A_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	lld = (long long int)(E_[nb1]->M_An->rpop_j) * (long long int)(E_[nb2]->M_An->rpop_j - 1*(nb1==nb2?1:0));
	dra_plusdivequals_s___m_m(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_TAnAtT->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	if (verbose){ raprintf(&(D->QC_TAnAtT_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_TAnAtT_nrm: ");}
	/* if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->A_rbother){ } */}
      /* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
  if (verbose){ printf(" %% [finished dcc_sumscores_nrm]\n");}
}

void dcc_sumscores_test()
{
  /* tests a combination of lrup and sumscores */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0); int iteration_max = GLOBAL_TEST_niter;
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  int nl=0; double ct=0,rt=0,r=0;
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); dcc_init_QX(D);
  GLOBAL_toc(1,1+verbose," %% loading time: ");
  if (verbose>1){ printf(" %% calculating initial half-loop-subscores.\n");}
  GLOBAL_tic(1);
  dcc_An_ajdk(D);
  dcc_lf_ZtSn(D);
  dcc_lf_D_AtTn_ZtSn_vv(D);
  if (error_check){ dcc_lf_D_AtTn_ZtSn_uu(D);}
  dcc_lf_TAnZtS_ww(D);
  if (error_check){ dcc_lf_TAnZtS_uu(D);}
  GLOBAL_toc(1,1+verbose," %% initial subscore time: ");
  if (error_check){ if (verbose>1){ printf(" %% checking errors.\n");}}
  if (error_check){ dcc_lf_D_AtTn_ZtSn_error(verbose,D);}
  if (error_check){ dcc_lf_TAnZtS_error(verbose,D);}
  dcc_halfloop(D);
  if (verbose>1){ printf(" %% beginning iteration.\n");}
  nl=0;
  while ((iteration_max<=0 || nl<iteration_max) && (D->A_rpop_j_total>0 || D->A_cpop_j>0)){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %0.5fs(%2.1fh)/%0.5fs(%2.1fh) = %.1f\n",nl,iteration_max,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r);}
    if (verbose>1){ printf(" %% combining subscores to form initial half-loop-scores.\n");}
    dcc_sumscores_nrm(D);
    dcc_sumscores_ifT(D);
    dcc_sumscores_srt(D);
    dcc_sumscores_cmb(D);
    if (verbose>1){ printf(" %% finding rows and columns with low scores.\n");}
    dcc_sumscores_xij(D);
    dcc_sumscores_mxA(D);
    dcc_sumscores_dmp(D);
    dcc_lrup_mxdup(D);
    dcc_M_mxset(D);
    if (verbose>1){ printf(" %% recalculating subscores.\n");}
    dcc_An_ajdk(D);
    dcc_lf_ZtSn(D);
    dcc_lf_D_AtTn_ZtSn_vv(D);
    if (error_check){ dcc_lf_D_AtTn_ZtSn_uu(D);}
    dcc_lf_TAnZtS_ww(D);
    if (error_check){ dcc_lf_TAnZtS_uu(D);}
    if (error_check){ if (verbose>1){ printf(" %% checking errors.\n");}}
    if (error_check){ dcc_lf_D_AtTn_ZtSn_error(verbose,D);}
    if (error_check){ dcc_lf_TAnZtS_error(verbose,D);}
    dcc_halfloop(D);
    nl++; /* while (nl<iteration_max && D->A_cpop_j>2 && D->A_rpop_j_total>2){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  wkspace_printf();
}

void dcc_time_sumscores_test()
{
  /* tests a combination of lrup and sumscores */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0); int iteration_max = GLOBAL_TEST_niter; int xdrop_total = 0;
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_); dcc_init_QX(D);
  GLOBAL_toc(1,1+verbose," %% loading time: ");
  xdrop_total = get_xdrop_total(D->A_rpop_j_total,D->A_cpop_j);
  GLOBAL_tic(1);
  GLOBAL_tic(4); dcc_An_ajdk(D); GLOBAL_toc(4,1+verbose," %% lf_An_ajdk: ");
  GLOBAL_tic(4); dcc_lf_ZtSn(D); GLOBAL_toc(4,1+verbose," %% lf_ZtSn: ");
  GLOBAL_tic(4); dcc_lf_D_AtTn_ZtSn_vv(D); GLOBAL_toc(4,1+verbose," %% lf_D_AtTn_ZtSn: ");
  GLOBAL_tic(4); dcc_lf_TAnZtS_ww(D); GLOBAL_toc(4,1+verbose," %% lf_TAnZtS: ");
  GLOBAL_toc(1,1+verbose," %% initial subscore : ");
  GLOBAL_tic(4); dcc_halfloop(D); GLOBAL_toc(4,1+verbose," %% halfloop: ");
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); it = rt;
  if (verbose>1){ printf(" %% elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; beginning iteration.\n",ct,ct/3600,rt,rt/3600,r);}
  GLOBAL_tic(3);
  nl=0;
  while ((iteration_max<=0 || nl<iteration_max) && (D->A_rpop_j_total>0 || D->A_cpop_j>0)){
    GLOBAL_toc(3,0,""); ct = GLOBAL_elct[3]; rt = GLOBAL_elrt[3]; r=ct/maximum(1,rt); et = rt/maximum(1,nl)*xdrop_total;
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh) + %6.1fs(%2.1fh)\n",nl,xdrop_total,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,it,it/3600,et,et/3600);}
    if (verbose>1){ printf(" %% combining subscores to form initial loop-scores.\n");}
    GLOBAL_tic(4); dcc_sumscores_nrm(D); GLOBAL_toc(4,1+verbose," %% sumscores_nrm: ");
    GLOBAL_tic(4); dcc_sumscores_ifT(D); GLOBAL_toc(4,1+verbose," %% sumscores_ifT: ");
    GLOBAL_tic(4); dcc_sumscores_srt(D); GLOBAL_toc(4,1+verbose," %% sumscores_srt: ");
    GLOBAL_tic(4); dcc_sumscores_cmb(D); GLOBAL_toc(4,1+verbose," %% sumscores_cmb: ");
    if (verbose>1){ printf(" %% finding rows and columns with low scores.\n");}
    GLOBAL_tic(4); dcc_sumscores_xij(D); GLOBAL_toc(4,1+verbose," %% sumscores_xij: ");
    GLOBAL_tic(4); dcc_sumscores_mxA(D); GLOBAL_toc(4,1+verbose," %% sumscores_mxA: ");
    GLOBAL_tic(4); dcc_sumscores_dmp(D); GLOBAL_toc(4,1+verbose," %% sumscores_dmp: ");
    GLOBAL_tic(4); dcc_lrup_mxdup(D); GLOBAL_toc(4,1+verbose," %% lrup_mxdup: ");
    GLOBAL_tic(4); dcc_M_mxset(D); GLOBAL_toc(4,1+verbose," %% M_mxset: ");
    GLOBAL_tic(4); dcc_An_ajdk(D); GLOBAL_toc(4,1+verbose," %% An_ajdk: ");
    GLOBAL_tic(4); dcc_lf_ZtSn(D); GLOBAL_toc(4,1+verbose," %% lf_ZtSn: ");
    GLOBAL_tic(4); dcc_lf_D_AtTn_ZtSn_vv(D); GLOBAL_toc(4,1+verbose," %% lf_D_AtTn_ZtSn: ");
    GLOBAL_tic(4); dcc_lf_TAnZtS_ww(D); GLOBAL_toc(4,1+verbose," %% lf_TAnZtS: ");
    GLOBAL_tic(4); dcc_halfloop(D); GLOBAL_toc(4,1+verbose," %% halfloop: ");
    nl++; /* while (nl<iteration_max && D->A_cpop_j>2 && D->A_rpop_j_total>2){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  wkspace_printf();
}



