#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void timing_dump(int nt)
{
  FILE *fp=NULL; char tempchar[FNAMESIZE];
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  sprintf(tempchar,"%s/timing.m",GLOBAL_DIR_NAME); if ((fp=fopen(tempchar,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",tempchar); exit(RET_READ_FAIL);}
  fprintf(fp,"elct = %0.2f; elrt = %0.2f;\n",GLOBAL_elct[nt],GLOBAL_elrt[nt]);
  fclose(fp);fp=NULL; 
}

void lakcluster_driver()
{
  /* runs main lakcluster driver */
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0; int xdrop_length;
  int *rdrop=NULL,*cdrop=NULL,*rkeep=NULL,*ckeep=NULL,*rcsum=NULL,*ccsum=NULL;
  double tau_c_est=0,tau_r_est=0;
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  if (verbose){ printf(" %% [entering lakcluster_driver]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  /* loads matrices from disc and sets masks. ;
     bcc_ajdk_load:
      Load binary mask D->A_bmc_b with size D->A_ncols from disc. ; Similarly for D->Y_bmc_b and D->Y_ncols, etc.;
      Copy unsigned char mask D->A_umc_b from D->A_bmc_b. ; Similary for D->Y_umc_b, etc. ;
      Set D->A_bmc_j to equal D->A_bmc_b, etc.; also calculate D->A_cpop_b, D->A_cpop_j, etc. ;
      Set D->Y_cbother. ;
     bcc_single_load_M_An:
      Load binary mask E->A_bmr_b withsize E->A_nrows from disc. ; Similarly for E->Z_bmr_b and E->Z_nrows, etc. ;
      Copy unsigned char mask E->A_umr_b from E->A_bmr_b. ; Similarly for E->Z_umr_b, etc. ;
      Set E->A_bmr_j to equal E->A_bmr_b, etc.;  also calculate E->A_rpop_b, E->A_rpop_j, etc. ;
      Set E->Z_rbother. ;
      Load E->M_An, E->M_At, E->M_Zn, etc. from disc. ;
      Initialize (but do not load) lrup (e.g., E->M_an, E->M_jn, etc.). ;
     bcc_X_nrows_total:
      Uses E->A_nrows, E->A_rpop_b, E->A_rpop_j to set D->A_nrows_total, D->A_rpop_b_total, D->A_rpop_j_total;
      USes E->Z_nrows, E->Z_rpop_b, E->Z_rpop_j to set D->Z_nrows_total, D->Z_rpop_b_total, D->Z_rpop_j_total;
     bcc_single_init_lf:
      Initialize L_handle and M_handle within E. ;
     bcc_double_init_lf:
      Initialize L_handle and M_handle within F. ;
     bcc_load_A_p:
      Calculates D->A_p and D->Y_p for each column using row-mask mr_b, not mr_j. ;
     bcc_init_QX:
      Initialize scores D->QC_AtTYnWtSZn_nrm, etc. ;
      Also define E->QR_index_global_mr_a, E->QR_index_global_mr_b, E->QR_index_global_mr_j. ;
     bcc_M_mxset:
      Define D->A_cpop_j using D->A_bmc_j. ; Similarly for D->Y_cpop_j, D->T_cpop_j.
      Define E->A_rpop_j using E->A_bmr_j. ; Similarly for E->Z_rpop_j. ;
      Also define D->A_rpop_j_total using E->A_rpop_j. ;
      Along the way set E->M_An->m_a_, etc. using M_mxset. ;
   */
  bcc_load(&D,&E_,&F_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D); 
  GLOBAL_toc(1,verbose," %% loading time: ");
  if (GLOBAL_scramble_num>0){
    GLOBAL_tic(1);
    wrap_bcc_scramble(D);
    GLOBAL_toc(1,verbose," %% shuffling time: ");
    /* if (GLOBAL_scramble_num>0){ } */}
  /*
    get_xdrop_array: Set xdrop_length to be the total number of iterations 
    required to reduce D->A_rpop_j_total,D->A_cpop_j to zero. ;
    Note that gamma is read from GLOBAL_gamma. ;
   */
  get_xdrop_array(D->A_rpop_j_total,D->A_cpop_j,&xdrop_length,&rdrop,&cdrop,&rkeep,&ckeep);
  ira_cumulative_sum(rkeep,xdrop_length,&rcsum); ira_cumulative_sum(ckeep,xdrop_length,&ccsum);
  if (verbose>1){ printf(" %% calculating initial loop-subscores.\n");}
  GLOBAL_tic(1);
  /* bcc_An_ajdk: calculating E_[nb]->lf_An_ajdk, E_[nb]->lf_Zn_ajdk, E_[nb]->lf_Yn_ajdk, E_[nb]->lf_Wn_ajdk. */
  bcc_An_ajdk(D);
  /* bcc_lf_ZtSn: calculating E_[nb]->lf_AtTn, E_[nb]->lf_ZtSn, E_[nb]->lf_YtTn, E_[nb]->lf_WtSn. */
  bcc_lf_ZtSn(D);
  /* bcc_lf_ZtSWn: calculating E_[nb]->lf_AtTAn, E_[nb]->lf_AtTYn, E_[nb]->lf_ZtSZn, E_[nb]->lf_ZtSWn. */
  bcc_lf_ZtSWn(D);
  /* bcc_lf_YnWt: calculating F_[nbx]->lf_YnWt, F_[nbx]->lf_YnYt, F_[nbx]->lf_AnZt, F_[nbx]->lf_AnAt. */
  bcc_lf_YnWt(D);
  /* bcc_M_ZtSWn: calculating E_[nb]->M_AtTAn, E_[nb]->M_AtTYn, E_[nb]->M_ZtSZn, E_[nb]->M_ZtSWn. */
  bcc_M_ZtSWn(D);
  /* bcc_YnWt: calculating F_[nbx]->M_YnWt, F_[nbx]->M_YnYt, F_[nbx]->M_AnZt, F_[nbx]->M_AnAt. */
  bcc_M_YnWt(D);
  /* bcc_lf_AtTYn____WtSZn: calculating F_[nbx]->lf_AtTYn____WtSZn, etc. */
  bcc_lf_AtTYn____WtSZn(D);
  /* bcc_lf_At_T_YnWt_S_Zn: calculating F_[nbx]->lf_At_T_YnWt_S_Zn, F_[nbx]->lf_At_T_AnZt_S_Zn, F_[nbx]->lf_At_T_YnYt_T_An, F_[nbx]->lf_At_T_AnAt_T_An. */
  bcc_lf_At_T_YnWt_S_Zn(D);
  /* bcc_lf_AnZt_S_WnYt: calculating F_[nbx]->lf_AnZt_S_WnYt, etc. */
  bcc_lf_AnZt_S_WnYt(D);
  /* fix later */ /* bcc_lf_An_ZtSWn_Yt(D); */
  GLOBAL_toc(1,verbose," %% initial subscore time: ");
  if (verbose>1){ printf(" %% correcting for collapsed-loops.\n");}
  GLOBAL_tic(1);
  /* bcc_singlestudy_ww: This function takes in a variety of inputs and calculates a variety of single-study terms. */
  bcc_singlestudy_ww(D);
  /* bcc_doublestudy_ww: This function takes in a variety of inputs and calculates a variety of double-study terms. */
  bcc_doublestudy_ww(D);
  /* bcc_flattenloop: 
      This function corrects the row- and col-scores, removing flattened-loops.
      Along the way we also convert from SPACING_a to SPACING_b along the T_ncols dimension.
      Warning! Later on we expect D->T_bmc_b to include only contiguous bits.
      This ensures that D->T_bmc_j will serve as a mask for QR_AnAtTAnAt etc. after switching to SPACING_b. 
  */
  bcc_flattenloop(D);
  GLOBAL_toc(1,verbose," %% collapsed-loop correction time: ");
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); it = rt;
  if (verbose>1){ printf(" %% elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; beginning iteration.\n",ct,ct/3600,rt,rt/3600,r);}
  GLOBAL_tic(3);
  nl=0;
  while (D->A_rpop_j_total>0 || D->A_cpop_j>0){
    GLOBAL_toc(3,0,""); ct = GLOBAL_elct[3]; rt = GLOBAL_elrt[3]; r=ct/maximum(1,rt); tau_c_est = rt/ccsum[nl]; tau_r_est = rt/rcsum[nl]; 
    //et = maximum(rcsum[xdrop_length-1]*tau_r_est,ccsum[xdrop_length-1]*tau_c_est);
    et = ccsum[xdrop_length-1]*tau_c_est;
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh) + %6.1fs(%2.1fh)\n",nl,xdrop_length,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,it,it/3600,et,et/3600);}
    if (verbose>1){ printf(" %% combining subscores to form initial loop-scores.\n");}
    /* bcc_sumscores_nrm: 
        If GLOBAL_Ireq>0, we normalize by the number of summands in each bin. ;
	Otherwise, if GLOBAL_Ireq<=0 we normalize by the number of summands across all bins, 
	and then sum the results. ; 
    */
    bcc_sumscores_nrm(D);
    /* bcc_sumscores_ifT: 
        If GLOBAL_Ireq>0, we correct each E_[nb1]->QR_AnAtTAnAt_nrm[0+nb2*A_nrows*T_ncols] 
	for T, across all values of nb1 and nb2. ;
	Otherwise, if GLOBAL_Ireq<=0 also correct each E_[nb1]->QR_AnAtTAnAt_nrm[0+nb2*A_nrows*T_ncols] 
	for T, but only use the nb2==0 term for each nb1.
    */
    bcc_sumscores_ifT(D);
    /* bcc_sumscores_srt:
        If GLOBAL_Ireq>0, then we sort the scores by bin. ;
        Otherwise, if GLOBAL_Ireq<=0, we use the summed scores over all bins (already stored in nb==0).
    */
    bcc_sumscores_srt(D);
    /* bcc_sumscores_cmb: Combines scores such as AtTAnAtTAn and AtTYnYtTAn to form final score. */
    bcc_sumscores_cmb(D); 
    if (verbose>1){ printf(" %% finding rows and columns with low scores.\n");}
    if (verbose>2){ printf(" %%  This step sets D->QC_index_local_mc, as well as D->QR_index_local_mr, D->QR_index_local_nb and D->QR_index_global_mr.\n");}
    /* bcc_sumscores_xij:
	Sets D->QC_index_local_mc_x.
	These arrays store the indices (i.e., index locations) of the masks 
	M_An->mc_x==D->A_bmc_x.
        Then sorts col-scores in increasing order across columns.
	At this point the first few entries of D->QC_index_local_mc_x correspond to the columns that are least likely to participate in a bicluster. ;
	Sets D->QR_index_local_mr_x and D->QR_index_local_nb.
	These arrays store the indices (i.e., index locations) of the masks 
	M_An[nb]->mr_x==E->A_bmr_x.
	Also copies D->QR_index_global_mr_x from E->QR_index_global_mr_x.
	These arrays store the overall index locations of the masks M_An->mr_x accumulated across bins. ;
        Then sorts row-scores in increasing order across rows.
	At this point the first few entries of D->QR_index_global_mr_x correspond to the rows that are least likely to participate in a bicluster. ;
	Also redefines D->A_rpop_j_total using E->A_bmr_b and E->A_bmr_j. 
     */
    bcc_sumscores_xij(D);
    if (verbose>1){ printf(" %% finding rows and columns with low scores.\n obtains rdrop and cdrop from get_xdrop ;");}
    if (verbose>2){ printf(" %% This step copies D->A_bmc_j to D->A_umc_j, defines D->A_umc_j_rmv and D->A_umc_j_rtn, and then removes the first cdrop entries of D->QC_index_local_mc_a.\n");}
    if (verbose>2){ printf(" %% This step also copies E->A_bmr_j to E->A_umr_j, defines E->A_umr_j_rmv and E->A_umr_j_rtn, and then removes the first rdrop entries of D->QR_index_local_mr_a, using D->QR_index_local_nb to index nb.\n");}
    if (verbose>2){ printf(" %% In addition, this step calls bcc_sumscores_mxB, which copies all instances of umc and umr to the associated bmc and bmr respectively.\n");}
    if (verbose>2){ printf(" %% bcc_sumscores_mxB also calls bcc_lrup_mxset, which copies all of the A_bmr_j_rtn etc (i.e., retained) and A_bmr_j_rmv etc (i.e., removed) bitmasks into the temporary variables used for lrup.\n");}
    /* bcc_sumscores_mxA:
        Obtains rdrop and cdrop from get_xdrop (based on GLOBAL_gamma). ;
	Copies D->A_bmc_j to D->A_umc_j. ; 
	Defines D->A_umc_j_rmv and D->A_umc_j_rtn ; removing first cdrop entries of D->QC_index_local_mc_a. ;
	Copies E->A_bmr_j to E->A_umr_j. ; 
	Defines E->A_umr_j_rmv and E->A_umr_j_rtn ; removing first rdrop entries of D->QR_index_local_mr_a ; using D->QR_index_local_nb to index nb. ;
	Call bcc_sumscores_mxB ; 
	 Copies all instances of umc and umr to associated bmc and bmr respectively ;
	 Note that bcc_sumscores_mxB also calls 
	 bcc_lrup_mxset:
	  Copies all of the 
	  A_bmr_j_rtn etc (i.e., retained) and 
	  A_bmr_j_rmv etc (i.e., removed) bitmasks 
	  into the temporary variables used for lrup 
	  (e.g., E->M_an used in low-rank update).
	Call bcc_sumscores_mxC:
	 Writing data to out_trace ;
	Finally: add to D->out_xdrop_a and D->out_xdrop_b (using D->QR_index_global_mr_x for rows, and D->QC_index_local_mc for cols). ; 
	Note: D->out_xdrop_x lists row-indices first, then col-indices ;
    */
    bcc_sumscores_mxA(D);
    if (verbose>1){ printf(" %% writing scores to disk.\n");}
    /* bcc_sumscores_dmp: dump temporary variables used within bcc to disc. */
    bcc_sumscores_dmp(D);
    if (verbose>1){ printf(" %% use low-rank-update to calculate new loop-subscores after removing rows and columns with low scores.\n");}
    /* bcc_lrup_QX_: These functions use the temporary-variables E->M_an, E->M_kn, E->M_jn, E->M_hn, etc. to update the scores QR and QC. */
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
    if (verbose>1){ printf(" %% copy D->A_bmc_j_rtn into D->A_bmc_j, and E_[nb]->A_bmr_j_rtn into E_[nb]->A_bmr_j.\n");}
    /* bcc_lrup_mxdup: 
        Copies D->A_bmc_j_rtn into D->A_bmc_j, and E->A_bmr_j_rtn into E->A_bmr_j.
        Copies D->A_umc_j_rtn into D->A_umc_j, and E->A_umr_j_rtn into E->A_umr_j.
    */
    bcc_lrup_mxdup(D);
    if (verbose>1){ printf(" %% Updates cpop_j and rpop_j, copies D->A_bmc_j into M_An->mc_j, and copies E_[nb]->A_bmr_j into M_An->mr_j.\n");}
    /* bcc_M_mxset:
        Updates D->A_cpop_j, etc. and E->A_rpop_j, etc. ;
        Copies D->A_bmc_j into E->M_An->mc_j, and updates E->M_An->n_a_, etc. ;
        Copies E_[nb]->A_bmr_j into E->M_An->mr_j and updates E->M_An->m_a_, etc.
    */
    bcc_M_mxset(D);
    if (verbose>1){ printf(" %% Convert each matrix of form lf_ZtSWn into binary M_ZtSWn\n");}
    /* bcc_M_ZtSWn: calculating E_[nb]->M_AtTAn, E_[nb]->M_AtTYn, E_[nb]->M_ZtSZn, E_[nb]->M_ZtSWn. */
    bcc_M_ZtSWn(D);
    if (verbose>1){ printf(" %% Convert each matrix of form lf_YnWt into binary M_YnWt\n");}
    /* bcc_M_YnWt: calculating F_[nbx]->M_YnWt, F_[nbx]->M_YnYt, F_[nbx]->M_AnZt, F_[nbx]->M_AnAt. */
    bcc_M_YnWt(D);
    if (verbose>1){ printf(" %% update column-scores if QC_strategy ZtSWn requires it.\n");}
    /* bcc_lf_AtTYn____WtSZn: calculating F_[nbx]->lf_AtTYn____WtSZn, etc. */
    bcc_lf_AtTYn____WtSZn(D);
    if (verbose>1){ printf(" %% update column-scores if QC_strategy YnWt requires it.\n");}
    /* bcc_lrup_QC_YnWt_stage_c: 
        Updates terms for QC such as F->lf_At_T_YnWt_S_Zn 
        using low-rank terms such as F->lf_kt_r_unwt_s_zn, etc. 
    */
    bcc_lrup_QC_YnWt_stage_c(D);
    if (verbose>1){ printf(" %% update row-scores if QR_strategy YnWt requires it.\n");}
    /* bcc_lf_AnZt_S_WnYt: calculating F_[nbx]->lf_AnZt_S_WnYt, etc. */
    bcc_lf_AnZt_S_WnYt(D);
    /* fix later */ /* bcc_lf_An_ZtSWn_Yt(D); */
    if (verbose>1){ printf(" %% correcting for collapsed-loops once again.\n");}
  /* bcc_An_ajdk: calculating E_[nb]->lf_An_ajdk, E_[nb]->lf_Zn_ajdk, E_[nb]->lf_Yn_ajdk, E_[nb]->lf_Wn_ajdk. */
    bcc_An_ajdk(D);
  /* bcc_lf_ZtSn: calculating E_[nb]->lf_AtTn, E_[nb]->lf_ZtSn, E_[nb]->lf_YtTn, E_[nb]->lf_WtSn. */
    bcc_lf_ZtSn(D);
  /* bcc_singlestudy_ww: This function takes in a variety of inputs and calculates a variety of single-study terms. */
    bcc_singlestudy_ww(D);
  /* bcc_doublestudy_ww: This function takes in a variety of inputs and calculates a variety of double-study terms. */
    bcc_doublestudy_ww(D);
  /* bcc_flattenloop: 
      This function corrects the row- and col-scores, removing flattened-loops.
      Along the way we also convert from SPACING_a to SPACING_b along the T_ncols dimension.
      Warning! Later on we expect D->T_bmc_b to include only contiguous bits.
      This ensures that D->T_bmc_j will serve as a mask for QR_AnAtTAnAt etc. after switching to SPACING_b. 
  */
    bcc_flattenloop(D);
    nl++; /* while (D->A_rpop_j_total>0 || D->A_cpop_j>0){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  timing_dump(2);
  wkspace_printf();
}

