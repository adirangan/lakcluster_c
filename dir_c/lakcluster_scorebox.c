#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void lakcluster_scorebox_excerpt_0(int verbose,struct bcc_ajdk *D,struct S_handle *S)
{
  int tab=0;
  double tmp_mean=0;
  if (verbose>1){ printf(" %% updating scorebox.\n");}
  tab = S->nr + S->nc*S->row_num;
  S->A_rpop_j_total[tab] = D->A_rpop_j_total;
  S->A_cpop_j[tab] = D->A_cpop_j;
  S->Irem[tab] = D->Irem;
  ra_stats(D->QR_svalue,"double",D->A_rpop_j_total,NULL,NULL,&tmp_mean,NULL); if (!isfinite(tmp_mean)){ tmp_mean = 1.0;}
  S->QR_avg[tab] = tmp_mean; 
  ra_stats(D->QC_svalue,"double",D->A_cpop_j,NULL,NULL,&tmp_mean,NULL); if (!isfinite(tmp_mean)){ tmp_mean = 1.0;}
  S->QC_avg[tab] = tmp_mean; 
}

void lakcluster_scorebox_excerpt_1(int verbose,struct bcc_ajdk *D,struct S_handle *S)
{
  if (verbose>1){ printf(" %% directly calculating loop-subscores.\n");}
  GLOBAL_tic(1);
  bcc_An_ajdk(D);
  bcc_lf_ZtSn(D);
  bcc_lf_ZtSWn(D);
  bcc_lf_YnWt(D);
  bcc_M_ZtSWn(D);
  bcc_M_YnWt(D);
  bcc_lf_AtTYn____WtSZn(D);
  bcc_lf_At_T_YnWt_S_Zn(D);
  bcc_lf_AnZt_S_WnYt(D);
  /* fix later */ /* bcc_lf_An_ZtSWn_Yt(D); */
  GLOBAL_toc(1,verbose," %% subscore calculation time: ");
  GLOBAL_tic(1);
  if (verbose>1){ printf(" %% correcting for collapsed-loops.\n");}
  bcc_singlestudy_ww(D);
  bcc_doublestudy_ww(D);
  bcc_flattenloop(D);
  GLOBAL_toc(1,verbose," %% collapsed-loop correction time: ");
  if (verbose>1){ printf(" %% combining subscores.\n");}
  bcc_sumscores_nrm(D);
  bcc_sumscores_ifT(D);
  bcc_sumscores_srt(D);
  bcc_sumscores_cmb(D);
  if (verbose>1){ printf(" %% calculating QX_svalue.\n");}
  bcc_scorebox_svalue(D);
  lakcluster_scorebox_excerpt_0(verbose,D,S);
  S_handle_printf(verbose,S," %% scorebox: ");
}

void lakcluster_scorebox_excerpt_2(int verbose,struct bcc_ajdk *D,struct S_handle *S)
{
  if (verbose>1){ printf(" %% use low-rank-update to recalculate loop-subscores after removing some rows and/or columns.\n");}
  GLOBAL_tic(1);
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
  bcc_lrup_QC_YnWt_stage_c(D);
  bcc_lf_AnZt_S_WnYt(D);
  /* fix later */ /* bcc_lf_An_ZtSWn_Yt(D); */
  GLOBAL_toc(1,verbose," %% low-rank update time: ");
  if (verbose>1){ printf(" %% correcting for collapsed-loops.\n");}
  GLOBAL_tic(1);
  bcc_An_ajdk(D);
  bcc_lf_ZtSn(D);
  bcc_singlestudy_ww(D);
  bcc_doublestudy_ww(D);
  bcc_flattenloop(D);
  GLOBAL_toc(1,verbose," %% collapsed-loop correction time: ");
  if (verbose>1){ printf(" %% combining subscores.\n");}
  bcc_sumscores_nrm(D);
  bcc_sumscores_ifT(D);
  bcc_sumscores_srt(D);
  bcc_sumscores_cmb(D);
  if (verbose>1){ printf(" %% calculating QX_svalue.\n");}
  bcc_scorebox_svalue(D);
  lakcluster_scorebox_excerpt_0(verbose,D,S);
}

void lakcluster_scorebox_rc()
{
  /* runs main lakcluster scorebox, varying rows first */
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  struct bcc_ajdk *D_ori=NULL;struct bcc_single **E_ori_=NULL; struct bcc_double **F_ori_=NULL;
  struct bcc_ajdk *D_sub=NULL;struct bcc_single **E_sub_=NULL; struct bcc_double **F_sub_=NULL;
  struct S_handle *S=NULL;
  if (verbose){ printf(" %% [entering lakcluster_scorebox]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  bcc_load(&D_ori,&E_ori_,&F_ori_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D_ori);
  bcc_load(&D_sub,&E_sub_,&F_sub_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D_sub);
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making S_handle.\n");}
  S = S_handle_make(GLOBAL_scorebox_out_xdrop,GLOBAL_scorebox_row_max,GLOBAL_scorebox_row_num,GLOBAL_scorebox_col_max,GLOBAL_scorebox_col_num);
  S_init_bcc(D_ori,S);  
  for (S->nc=0;S->nc<S->col_num;S->nc++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*S->col_num)/maximum(1,S->nc);
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",S->nc,S->col_num,D_ori->A_rpop_j_total,D_ori->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    S->nr=0;
    if (S->nc==0){
      if (verbose>1){ printf(" %% initial drop: S->rdrop[%d] %d S->cdrop[%d] %d\n",S->nr,S->rdrop[S->nr],S->nc,S->cdrop[S->nc]);}
      bcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_index_local_nb,S->mr_index_local_mr,S->out_xdrop_ncols,S->mc_index_sort);
      bcc_scorebox_mxA(D_ori,S->rdrop[S->nr],S->cdrop[S->nc]);
      if (verbose>1){ printf(" %% calling bcc_lrup_mxdup\n");} bcc_lrup_mxdup(D_ori);
      if (verbose>1){ printf(" %% calling bcc_M_mxset\n");} bcc_M_mxset(D_ori);
      lakcluster_scorebox_excerpt_1(verbose,D_ori,S);
      /* if (S->nc==0){ } */}
    if (S->nc>0){
      if (verbose>1){ printf(" %% cdrop: S->cdrop[%d] %d\n",S->nc,S->cdrop[S->nc]);}
      bcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_index_local_nb,S->mr_index_local_mr,S->out_xdrop_ncols,S->mc_index_sort);
      bcc_scorebox_mxA(D_ori,0,S->cdrop[S->nc]);
      lakcluster_scorebox_excerpt_2(verbose,D_ori,S);
      /* if (S->nc>0){ } */}
    if (verbose>1){ printf(" %% copying D_ori to D_sub.\n");} bcc_copy(D_sub,D_ori);    
    for (S->nr=1;S->nr<S->row_num;S->nr++){
      if (verbose>1){ printf(" %% S->nc %d; rdrop: S->rdrop[%d] %d\n",S->nc,S->nr,S->rdrop[S->nr]);}
      bcc_scorebox_srt(D_sub,S->out_xdrop_nrows,S->mr_index_local_nb,S->mr_index_local_mr,S->out_xdrop_ncols,S->mc_index_sort);
      bcc_scorebox_mxA(D_sub,S->rdrop[S->nr],0);
      lakcluster_scorebox_excerpt_2(verbose,D_sub,S);
      /* for (S->nr=1;S->nr<S->row_num;S->nr++){ } */}
    /* for (S->nc=0;S->nc<S->col_num;S->nc++){ } */}
  S_handle_printf(verbose,S," %% scorebox: ");    
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  S_handle_dmp(S);
  wkspace_printf();
}

void lakcluster_scorebox_cr()
{
  /* runs main lakcluster scorebox, varying cols first */
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  struct bcc_ajdk *D_ori=NULL;struct bcc_single **E_ori_=NULL; struct bcc_double **F_ori_=NULL;
  struct bcc_ajdk *D_sub=NULL;struct bcc_single **E_sub_=NULL; struct bcc_double **F_sub_=NULL;
  struct S_handle *S=NULL;
  if (verbose){ printf(" %% [entering lakcluster_scorebox]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  bcc_load(&D_ori,&E_ori_,&F_ori_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D_ori);
  bcc_load(&D_sub,&E_sub_,&F_sub_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D_sub);
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making S_handle.\n");}
  S = S_handle_make(GLOBAL_scorebox_out_xdrop,GLOBAL_scorebox_row_max,GLOBAL_scorebox_row_num,GLOBAL_scorebox_col_max,GLOBAL_scorebox_col_num);
  S_init_bcc(D_ori,S);  
  for (S->nr=0;S->nr<S->row_num;S->nr++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*S->row_num)/maximum(1,S->nr);
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",S->nr,S->row_num,D_ori->A_rpop_j_total,D_ori->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    S->nc=0;
    if (S->nr==0){
      if (verbose>1){ printf(" %% initial drop: S->rdrop[%d] %d S->cdrop[%d] %d\n",S->nr,S->rdrop[S->nr],S->nc,S->cdrop[S->nc]);}
      bcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_index_local_nb,S->mr_index_local_mr,S->out_xdrop_ncols,S->mc_index_sort);
      bcc_scorebox_mxA(D_ori,S->rdrop[S->nr],S->cdrop[S->nc]);
      if (verbose>1){ printf(" %% calling bcc_lrup_mxdup\n");} bcc_lrup_mxdup(D_ori);
      if (verbose>1){ printf(" %% calling bcc_M_mxset\n");} bcc_M_mxset(D_ori);
      lakcluster_scorebox_excerpt_1(verbose,D_ori,S);
      /* if (S->nr==0){ } */}
    if (S->nr>0){
      if (verbose>1){ printf(" %% rdrop: S->rdrop[%d] %d\n",S->nr,S->rdrop[S->nr]);}
      bcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_index_local_nb,S->mr_index_local_mr,S->out_xdrop_ncols,S->mc_index_sort);
      bcc_scorebox_mxA(D_ori,S->rdrop[S->nr],0);
      lakcluster_scorebox_excerpt_2(verbose,D_ori,S);
      /* if (S->nr>0){ } */}
    if (verbose>1){ printf(" %% copying D_ori to D_sub.\n");} bcc_copy(D_sub,D_ori);
    for (S->nc=1;S->nc<S->col_num;S->nc++){
      if (verbose>1){ printf(" %% S->nr %d; cdrop: S->cdrop[%d] %d\n",S->nr,S->nc,S->cdrop[S->nc]);}
      bcc_scorebox_srt(D_sub,S->out_xdrop_nrows,S->mr_index_local_nb,S->mr_index_local_mr,S->out_xdrop_ncols,S->mc_index_sort);
      bcc_scorebox_mxA(D_sub,0,S->cdrop[S->nc]);
      lakcluster_scorebox_excerpt_2(verbose,D_sub,S);
      /* for (S->nc=0;S->nc<S->col_num;S->nc++){ } */}
    /* for (S->nr=1;S->nr<S->row_num;S->nr++){ } */}
  S_handle_printf(verbose,S," %% scorebox: ");
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  S_handle_dmp(S);
  wkspace_printf();
}

void lakcluster_scorebox_xx()
{
  /* runs main lakcluster scorebox, calculating full score at each stage */
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  struct bcc_ajdk *D_ori=NULL;struct bcc_single **E_ori_=NULL; struct bcc_double **F_ori_=NULL;
  struct bcc_ajdk *D_sub=NULL;struct bcc_single **E_sub_=NULL; struct bcc_double **F_sub_=NULL;
  struct S_handle *S=NULL;
  if (verbose){ printf(" %% [entering lakcluster_scorebox]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  bcc_load(&D_ori,&E_ori_,&F_ori_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D_ori);
  bcc_load(&D_sub,&E_sub_,&F_sub_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D_sub);
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (verbose>1){ printf(" %% making S_handle.\n");}
  S = S_handle_make(GLOBAL_scorebox_out_xdrop,GLOBAL_scorebox_row_max,GLOBAL_scorebox_row_num,GLOBAL_scorebox_col_max,GLOBAL_scorebox_col_num);
  S_init_bcc(D_ori,S);  
  for (S->nr=0;S->nr<S->row_num;S->nr++){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); et = (rt*S->row_num)/maximum(1,S->nr);
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",S->nr,S->row_num,D_ori->A_rpop_j_total,D_ori->A_cpop_j,ct,ct/3600,rt,rt/3600,r,et,et/3600);}
    S->nc=0;
    if (S->nr==0){
      if (verbose>1){ printf(" %% initial drop: S->rdrop[%d] %d S->cdrop[%d] %d\n",S->nr,S->rdrop[S->nr],S->nc,S->cdrop[S->nc]);}
      if (verbose>2){ printf(" %% This step sets D->QC_index_local_mc_a to be the final D->A_cpop_j entries of S->mc_index_sort.\n");}
      if (verbose>2){ printf(" %% This step also sets D->QR_index_local_nb and D->QR_index_local_mr_a to be the final D->A_rpop_j entries of mr_index_sort (i.e., mr_index_local_nb and mr_index_local_mr).\n");}
      bcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_index_local_nb,S->mr_index_local_mr,S->out_xdrop_ncols,S->mc_index_sort);
      if (verbose>2){ printf(" %% This step does the following in order: \n");}
      if (verbose>2){ printf(" %% copies D->A_bmc_j to D->A_umc_j ; \n");}
      if (verbose>2){ printf(" %% define D->A_umc_j_rmv and D->A_umc_j_rtn ; removing first cdrop entries of D->QC_index_local_mc_a ; \n");}
      if (verbose>2){ printf(" %% copies E->A_bmr_j to E->A_umr_j ; \n");}
      if (verbose>2){ printf(" %% define E->A_umr_j_rmv and E->A_umr_j_rtn ; removing first rdrop entries of D->QR_index_local_mr_a ; using D->QR_index_local_nb to index nb ; \n");}
      if (verbose>2){ printf(" %% call bcc_sumscores_mxB ; copying all instances of umc and umr to associated bmc and bmr respectively ; \n");}
      if (verbose>2){ printf(" %% Note that bcc_sumscores_mxB also calls bcc_lrup_mxset, which copies all of the A_bmr_j_rtn etc (i.e., retained) and A_bmr_j_rmv etc (i.e., removed) bitmasks into the temporary variables used for lrup.\n");}
      bcc_scorebox_mxA(D_ori,S->rdrop[S->nr],S->cdrop[S->nc]);
      if (verbose>1){ printf(" %% calling bcc_lrup_mxdup\n");} 
      if (verbose>2){ printf(" %% This step copies D->A_bmc_j_rtn into D->A_bmc_j, and E_[nb]->A_bmr_j_rtn into E_[nb]->A_bmr_j.\n");} 
      bcc_lrup_mxdup(D_ori);
      if (verbose>1){ printf(" %% calling bcc_M_mxset\n");} 
      if (verbose>2){ printf(" %% This step updates cpop_j and rpop_j, copies D->A_bmc_j into M_An->mc_j, and copies E_[nb]->A_bmr_j into M_An->mr_j.\n");}
      bcc_M_mxset(D_ori);
      if (verbose>1){ printf(" %% calling lakcluster_scorebox_excerpt_1\n");} 
      if (verbose>2){ printf(" %% This step calculates the original scores (without using lrup)\n");}
      lakcluster_scorebox_excerpt_1(verbose,D_ori,S);
      /* if (S->nr==0){ } */}
    if (S->nr>0){
      if (verbose>1){ printf(" %% rdrop: S->rdrop[%d] %d\n",S->nr,S->rdrop[S->nr]);}
      bcc_scorebox_srt(D_ori,S->out_xdrop_nrows,S->mr_index_local_nb,S->mr_index_local_mr,S->out_xdrop_ncols,S->mc_index_sort);
      bcc_scorebox_mxA(D_ori,S->rdrop[S->nr],0);
      if (verbose>1){ printf(" %% calling bcc_lrup_mxdup\n");} bcc_lrup_mxdup(D_ori);
      if (verbose>1){ printf(" %% calling bcc_M_mxset\n");} bcc_M_mxset(D_ori);
      lakcluster_scorebox_excerpt_1(verbose,D_ori,S);
      /* if (S->nr>0){ } */}
    if (verbose>1){ printf(" %% copying D_ori to D_sub.\n");} bcc_copy(D_sub,D_ori);
    for (S->nc=1;S->nc<S->col_num;S->nc++){
      if (verbose>1){ printf(" %% S->nr %d; cdrop: S->cdrop[%d] %d\n",S->nr,S->nc,S->cdrop[S->nc]);}
      bcc_scorebox_srt(D_sub,S->out_xdrop_nrows,S->mr_index_local_nb,S->mr_index_local_mr,S->out_xdrop_ncols,S->mc_index_sort);
      bcc_scorebox_mxA(D_sub,0,S->cdrop[S->nc]);
      if (verbose>1){ printf(" %% calling bcc_lrup_mxdup\n");} bcc_lrup_mxdup(D_sub);
      if (verbose>1){ printf(" %% calling bcc_M_mxset\n");} bcc_M_mxset(D_sub);
      lakcluster_scorebox_excerpt_1(verbose,D_sub,S);
      /* for (S->nc=0;S->nc<S->col_num;S->nc++){ } */}
    /* for (S->nr=1;S->nr<S->row_num;S->nr++){ } */}
  S_handle_printf(verbose,S," %% scorebox: ");
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  S_handle_dmp(S);
  wkspace_printf();
}

void lakcluster_scorebox()
{
  if (0){}
  else if (GLOBAL_scorebox_row_max> GLOBAL_scorebox_col_max){ lakcluster_scorebox_rc();}
  else if (GLOBAL_scorebox_row_max<=GLOBAL_scorebox_col_max){ lakcluster_scorebox_cr();}
  else{ /* full calculation at each step: */ lakcluster_scorebox_xx();}
}
