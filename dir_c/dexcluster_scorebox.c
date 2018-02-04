
void dexcluster_scorebox_excerpt_0(int verbose,struct dcc_ajdk *D,struct S_handle *S)
{
  int tab=0;
  double tmp_mean=0;
  if (verbose>1){ printf(" %% updating scorebox.\n");}
  tab = S->nr + S->nc*S->row_num;
  S->A_rpop_j_total[tab] = D->A_rpop_j_total;
  S->A_cpop_j[tab] = D->A_cpop_j;
  S->Irem[tab] = D->Irem;
  ra_stats(D->QR_sra,"double",D->A_rpop_j_total,NULL,NULL,&tmp_mean,NULL); if (!isfinite(tmp_mean)){ tmp_mean = 1.0;}
  S->QR_avg[tab] = tmp_mean; 
  ra_stats(D->QC_sra,"double",D->A_cpop_j,NULL,NULL,&tmp_mean,NULL); if (!isfinite(tmp_mean)){ tmp_mean = 1.0;}
  S->QC_avg[tab] = tmp_mean; 
}

void dexcluster_scorebox_excerpt_1(int verbose,struct dcc_ajdk *D,struct S_handle *S)
{
  if (verbose>1){ printf(" %% directly calculating halfloop-subscores.\n");}
  GLOBAL_tic(1);
  dcc_An_ajdk(D); 
  dcc_lf_ZtSn(D); 
  dcc_lf_D_AtTn_ZtSn_vv(D); 
  dcc_lf_TAnZtS_ww(D); 
  GLOBAL_toc(1,verbose," %% subscore calculation time: ");
  GLOBAL_tic(1);
  if (verbose>1){ printf(" %% correcting for collapsed-halfloops.\n");}
  dcc_halfloop(D); 
  GLOBAL_toc(1,verbose," %% collapsed-halfloop correction time: ");
  if (verbose>1){ printf(" %% combining subscores.\n");}
  dcc_sumscores_nrm(D);
  dcc_sumscores_ifT(D);
  dcc_sumscores_srt(D);
  dcc_sumscores_cmb(D);
  if (verbose>1){ printf(" %% calculating QX_sra.\n");}
  dcc_scorebox_sra(D);
  dexcluster_scorebox_excerpt_0(verbose,D,S);
  S_handle_printf(verbose,S," %% scorebox: ");
}

void dexcluster_scorebox_excerpt_2(int verbose,struct dcc_ajdk *D,struct S_handle *S)
{
  if (verbose>1){ printf(" %% use low-rank-update to recalculate halfloop-subscores after removing some rows and/or columns.\n");}
  GLOBAL_tic(1);
  dcc_lrup_mxdup(D);
  dcc_M_mxset(D); 
  dcc_An_ajdk(D); 
  dcc_lf_ZtSn(D); 
  dcc_lf_D_AtTn_ZtSn_vv(D); 
  dcc_lf_TAnZtS_ww(D); 
  dcc_halfloop(D); 
  GLOBAL_toc(1,verbose," %% low-rank update time: ");
  if (verbose>1){ printf(" %% combining subscores.\n");}
  dcc_sumscores_nrm(D);
  dcc_sumscores_ifT(D);
  dcc_sumscores_srt(D);
  dcc_sumscores_cmb(D);
  if (verbose>1){ printf(" %% calculating QX_sra.\n");}
  dcc_scorebox_sra(D);
  dexcluster_scorebox_excerpt_0(verbose,D,S);
}

void dexcluster_scorebox_rc()
{
  /* runs main dexcluster scorebox, varying rows first */
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  struct dcc_ajdk *D_ori=NULL;struct dcc_single **E_ori_=NULL; struct dcc_double **F_ori_=NULL;
  struct dcc_ajdk *D_sub=NULL;struct dcc_single **E_sub_=NULL; struct dcc_double **F_sub_=NULL;
  struct S_handle *S=NULL;
  if (verbose){ printf(" %% [entering dexcluster_scorebox]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_load(&D_ori,&E_ori_,&F_ori_); dcc_init_QX(D_ori);
  dcc_load(&D_sub,&E_sub_,&F_sub_); dcc_init_QX(D_sub);
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making S_handle.\n");}
  S = S_handle_make(GLOBAL_scorebox_out_xdrop,GLOBAL_scorebox_row_max,GLOBAL_scorebox_row_num,GLOBAL_scorebox_col_max,GLOBAL_scorebox_col_num);
  S_init_dcc(D_ori,S);  
  for (S->nc=0;S->nc<S->col_num;S->nc++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*S->col_num)/maximum(1,S->nc);
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",S->nc,S->col_num,D_ori->A_rpop_j_total,D_ori->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    S->nr=0;
    if (S->nc==0){
      if (verbose>1){ printf(" %% initial drop: S->rdrop[%d] %d S->cdrop[%d] %d\n",S->nr,S->rdrop[S->nr],S->nc,S->cdrop[S->nc]);}
      dcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_lnb,S->mr_lmr,S->out_xdrop_ncols,S->mc_srt);
      dcc_scorebox_mxA(D_ori,S->rdrop[S->nr],S->cdrop[S->nc]);
      if (verbose>1){ printf(" %% calling dcc_lrup_mxdup\n");} dcc_lrup_mxdup(D_ori);
      if (verbose>1){ printf(" %% calling dcc_M_mxset\n");} dcc_M_mxset(D_ori);
      dexcluster_scorebox_excerpt_1(verbose,D_ori,S);
      /* if (S->nc==0){ } */}
    if (S->nc>0){
      if (verbose>1){ printf(" %% cdrop: S->cdrop[%d] %d\n",S->nc,S->cdrop[S->nc]);}
      dcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_lnb,S->mr_lmr,S->out_xdrop_ncols,S->mc_srt);
      dcc_scorebox_mxA(D_ori,0,S->cdrop[S->nc]);
      dexcluster_scorebox_excerpt_2(verbose,D_ori,S);
      /* if (S->nc>0){ } */}
    if (verbose>1){ printf(" %% copying D_ori to D_sub.\n");} dcc_copy(D_sub,D_ori);    
    for (S->nr=1;S->nr<S->row_num;S->nr++){
      if (verbose>1){ printf(" %% S->nc %d; rdrop: S->rdrop[%d] %d\n",S->nc,S->nr,S->rdrop[S->nr]);}
      dcc_scorebox_srt(D_sub,S->out_xdrop_nrows,S->mr_lnb,S->mr_lmr,S->out_xdrop_ncols,S->mc_srt);
      dcc_scorebox_mxA(D_sub,S->rdrop[S->nr],0);
      dexcluster_scorebox_excerpt_2(verbose,D_sub,S);
      /* for (S->nr=1;S->nr<S->row_num;S->nr++){ } */}
    /* for (S->nc=0;S->nc<S->col_num;S->nc++){ } */}
  S_handle_printf(verbose,S," %% scorebox: ");    
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  S_handle_dmp(S);
  wkspace_printf();
}

void dexcluster_scorebox_cr()
{
  /* runs main dexcluster scorebox, varying cols first */
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  struct dcc_ajdk *D_ori=NULL;struct dcc_single **E_ori_=NULL; struct dcc_double **F_ori_=NULL;
  struct dcc_ajdk *D_sub=NULL;struct dcc_single **E_sub_=NULL; struct dcc_double **F_sub_=NULL;
  struct S_handle *S=NULL;
  if (verbose){ printf(" %% [entering dexcluster_scorebox]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_load(&D_ori,&E_ori_,&F_ori_); dcc_init_QX(D_ori);
  dcc_load(&D_sub,&E_sub_,&F_sub_); dcc_init_QX(D_sub);
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making S_handle.\n");}
  S = S_handle_make(GLOBAL_scorebox_out_xdrop,GLOBAL_scorebox_row_max,GLOBAL_scorebox_row_num,GLOBAL_scorebox_col_max,GLOBAL_scorebox_col_num);
  S_init_dcc(D_ori,S);  
  for (S->nr=0;S->nr<S->row_num;S->nr++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*S->row_num)/maximum(1,S->nr);
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",S->nr,S->row_num,D_ori->A_rpop_j_total,D_ori->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    S->nc=0;
    if (S->nr==0){
      if (verbose>1){ printf(" %% initial drop: S->rdrop[%d] %d S->cdrop[%d] %d\n",S->nr,S->rdrop[S->nr],S->nc,S->cdrop[S->nc]);}
      dcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_lnb,S->mr_lmr,S->out_xdrop_ncols,S->mc_srt);
      dcc_scorebox_mxA(D_ori,S->rdrop[S->nr],S->cdrop[S->nc]);
      if (verbose>1){ printf(" %% calling dcc_lrup_mxdup\n");} dcc_lrup_mxdup(D_ori);
      if (verbose>1){ printf(" %% calling dcc_M_mxset\n");} dcc_M_mxset(D_ori);
      dexcluster_scorebox_excerpt_1(verbose,D_ori,S);
      /* if (S->nr==0){ } */}
    if (S->nr>0){
      if (verbose>1){ printf(" %% rdrop: S->rdrop[%d] %d\n",S->nr,S->rdrop[S->nr]);}
      dcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_lnb,S->mr_lmr,S->out_xdrop_ncols,S->mc_srt);
      dcc_scorebox_mxA(D_ori,S->rdrop[S->nr],0);
      dexcluster_scorebox_excerpt_2(verbose,D_ori,S);
      /* if (S->nr>0){ } */}
    if (verbose>1){ printf(" %% copying D_ori to D_sub.\n");} dcc_copy(D_sub,D_ori);
    for (S->nc=1;S->nc<S->col_num;S->nc++){
      if (verbose>1){ printf(" %% S->nr %d; cdrop: S->cdrop[%d] %d\n",S->nr,S->nc,S->cdrop[S->nc]);}
      dcc_scorebox_srt(D_sub,S->out_xdrop_nrows,S->mr_lnb,S->mr_lmr,S->out_xdrop_ncols,S->mc_srt);
      dcc_scorebox_mxA(D_sub,0,S->cdrop[S->nc]);
      dexcluster_scorebox_excerpt_2(verbose,D_sub,S);
      /* for (S->nc=0;S->nc<S->col_num;S->nc++){ } */}
    /* for (S->nr=1;S->nr<S->row_num;S->nr++){ } */}
  S_handle_printf(verbose,S," %% scorebox: ");
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  S_handle_dmp(S);
  wkspace_printf();
}

void dexcluster_scorebox_xx()
{
  /* runs main dexcluster scorebox, calculating full score at each stage */
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  struct dcc_ajdk *D_ori=NULL;struct dcc_single **E_ori_=NULL; struct dcc_double **F_ori_=NULL;
  struct dcc_ajdk *D_sub=NULL;struct dcc_single **E_sub_=NULL; struct dcc_double **F_sub_=NULL;
  struct S_handle *S=NULL;
  if (verbose){ printf(" %% [entering dexcluster_scorebox]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_load(&D_ori,&E_ori_,&F_ori_); dcc_init_QX(D_ori);
  dcc_load(&D_sub,&E_sub_,&F_sub_); dcc_init_QX(D_sub);
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making S_handle.\n");}
  S = S_handle_make(GLOBAL_scorebox_out_xdrop,GLOBAL_scorebox_row_max,GLOBAL_scorebox_row_num,GLOBAL_scorebox_col_max,GLOBAL_scorebox_col_num);
  S_init_dcc(D_ori,S);  
  for (S->nr=0;S->nr<S->row_num;S->nr++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*S->row_num)/maximum(1,S->nr);
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",S->nr,S->row_num,D_ori->A_rpop_j_total,D_ori->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    S->nc=0;
    if (S->nr==0){
      if (verbose>1){ printf(" %% initial drop: S->rdrop[%d] %d S->cdrop[%d] %d\n",S->nr,S->rdrop[S->nr],S->nc,S->cdrop[S->nc]);}
      dcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_lnb,S->mr_lmr,S->out_xdrop_ncols,S->mc_srt);
      dcc_scorebox_mxA(D_ori,S->rdrop[S->nr],S->cdrop[S->nc]);
      if (verbose>1){ printf(" %% calling dcc_lrup_mxdup\n");} dcc_lrup_mxdup(D_ori);
      if (verbose>1){ printf(" %% calling dcc_M_mxset\n");} dcc_M_mxset(D_ori);
      dexcluster_scorebox_excerpt_1(verbose,D_ori,S);
      /* if (S->nr==0){ } */}
    if (S->nr>0){
      if (verbose>1){ printf(" %% rdrop: S->rdrop[%d] %d\n",S->nr,S->rdrop[S->nr]);}
      dcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_lnb,S->mr_lmr,S->out_xdrop_ncols,S->mc_srt);
      dcc_scorebox_mxA(D_ori,S->rdrop[S->nr],0);
      if (verbose>1){ printf(" %% calling dcc_lrup_mxdup\n");} dcc_lrup_mxdup(D_ori);
      if (verbose>1){ printf(" %% calling dcc_M_mxset\n");} dcc_M_mxset(D_ori);
      dexcluster_scorebox_excerpt_1(verbose,D_ori,S);
      /* if (S->nr>0){ } */}
    if (verbose>1){ printf(" %% copying D_ori to D_sub.\n");} dcc_copy(D_sub,D_ori);
    for (S->nc=1;S->nc<S->col_num;S->nc++){
      if (verbose>1){ printf(" %% S->nr %d; cdrop: S->cdrop[%d] %d\n",S->nr,S->nc,S->cdrop[S->nc]);}
      dcc_scorebox_srt(D_sub,S->out_xdrop_nrows,S->mr_lnb,S->mr_lmr,S->out_xdrop_ncols,S->mc_srt);
      dcc_scorebox_mxA(D_sub,0,S->cdrop[S->nc]);
      if (verbose>1){ printf(" %% calling dcc_lrup_mxdup\n");} dcc_lrup_mxdup(D_sub);
      if (verbose>1){ printf(" %% calling dcc_M_mxset\n");} dcc_M_mxset(D_sub);
      dexcluster_scorebox_excerpt_1(verbose,D_sub,S);
      /* for (S->nc=0;S->nc<S->col_num;S->nc++){ } */}
    /* for (S->nr=1;S->nr<S->row_num;S->nr++){ } */}
  S_handle_printf(verbose,S," %% scorebox: ");
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  S_handle_dmp(S);
  wkspace_printf();
}

void dexcluster_scorebox()
{
  if (0){}
  else if (GLOBAL_scorebox_row_max> GLOBAL_scorebox_col_max){ dexcluster_scorebox_rc();}
  else if (GLOBAL_scorebox_row_max<=GLOBAL_scorebox_col_max){ dexcluster_scorebox_cr();}
  else{ /* full calculation at each step: */ dexcluster_scorebox_xx();}
}
