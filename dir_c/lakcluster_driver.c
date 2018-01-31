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
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0; int xdrop_total;
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  if (verbose){ printf(" %% [entering lakcluster_driver]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  bcc_load(&D,&E_,&F_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D);
  GLOBAL_toc(1,verbose," %% loading time: ");
  xdrop_total = get_xdrop_total(D->A_rpop_j_total,D->A_cpop_j);
  if (verbose>1){ printf(" %% calculating initial loop-subscores.\n");}
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
  GLOBAL_toc(1,verbose," %% initial subscore time: ");
  if (verbose>1){ printf(" %% correcting for collapsed-loops.\n");}
  GLOBAL_tic(1);
  bcc_singlestudy_ww(D);
  bcc_doublestudy_ww(D);
  bcc_flattenloop(D);
  GLOBAL_toc(1,verbose," %% collapsed-loop correction time: ");
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); it = rt;
  if (verbose>1){ printf(" %% elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; beginning iteration.\n",ct,ct/3600,rt,rt/3600,r);}
  GLOBAL_tic(3);
  nl=0;
  while (D->A_rpop_j_total>0 || D->A_cpop_j>0){
    GLOBAL_toc(3,0,""); ct = GLOBAL_elct[3]; rt = GLOBAL_elrt[3]; r=ct/maximum(1,rt); et = rt/maximum(1,nl)*xdrop_total;
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh) + %6.1fs(%2.1fh)\n",nl,xdrop_total,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,it,it/3600,et,et/3600);}
    if (verbose>1){ printf(" %% combining subscores to form initial loop-scores.\n");}
    bcc_sumscores_nrm(D);
    bcc_sumscores_ifT(D);
    bcc_sumscores_srt(D);
    bcc_sumscores_cmb(D); 
    if (verbose>1){ printf(" %% finding rows and columns with low scores.\n");}
    bcc_sumscores_xij(D);
    bcc_sumscores_mxA(D);
    bcc_sumscores_dmp(D);
    if (verbose>1){ printf(" %% use low-rank-update to calculate new loop-subscores after removing rows and columns with low scores.\n");}
    bcc_lrup_QR_YnWt_stage_0(D);
    bcc_lrup_QR_YnWt_stage_1(D);
    bcc_lrup_QR_YnWt_stage_2(D);
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
    if (verbose>1){ printf(" %% correcting for collapsed-loops once again.\n");}
    bcc_An_ajdk(D);
    bcc_lf_ZtSn(D);
    bcc_singlestudy_ww(D);
    bcc_doublestudy_ww(D);
    bcc_flattenloop(D);
    nl++; /* while (D->A_rpop_j_total>0 || D->A_cpop_j>0){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  timing_dump(2);
  wkspace_printf();
}
