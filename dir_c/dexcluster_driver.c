#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void dexcluster_driver()
{
  /* runs main dexcluster driver */
  int verbose=GLOBAL_verbose; int nl=0; double ct=0,rt=0,r=0,it=0,et=0; int xdrop_length;
  struct dcc_ajdk *D=NULL;struct dcc_single **E_=NULL; struct dcc_double **F_=NULL;
  if (verbose){ printf(" %% [entering dexcluster_driver]\n");}
  if (verbose>1){ printf(" %% loading data from disk.\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  dcc_load(&D,&E_,&F_); dcc_init_QX(D);
  GLOBAL_toc(1,1+verbose," %% loading time: ");
  if (GLOBAL_scramble_num>0){
    GLOBAL_tic(1);
    wrap_dcc_scramble(D);
    GLOBAL_toc(1,verbose," %% shuffling time: ");
    /* if (GLOBAL_scramble_num>0){ } */}
  xdrop_length = get_xdrop_length((double)(D->A_rpop_j_total),(double)(D->A_cpop_j));
  GLOBAL_tic(1);
  dcc_An_ajdk(D); 
  dcc_lf_ZtSn(D); 
  dcc_lf_D_AtTn_ZtSn_vv(D); 
  dcc_lf_TAnZtS_ww(D); 
  GLOBAL_toc(1,1+verbose," %% initial subscore : ");
  dcc_halfloop(D); 
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); it = rt;
  if (verbose>1){ printf(" %% elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; beginning iteration.\n",ct,ct/3600,rt,rt/3600,r);}
  GLOBAL_tic(3);
  nl=0;
  while (D->A_rpop_j_total>0 || D->A_cpop_j>0){
    GLOBAL_toc(3,0,""); ct = GLOBAL_elct[3]; rt = GLOBAL_elrt[3]; r=ct/maximum(1,rt); et = rt/maximum(1,nl)*xdrop_length;
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh) + %6.1fs(%2.1fh)\n",nl,xdrop_length,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,it,it/3600,et,et/3600);}
    if (verbose>1){ printf(" %% combining subscores to form initial loop-scores.\n");}
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
    dcc_An_ajdk(D); 
    dcc_lf_ZtSn(D); 
    dcc_lf_D_AtTn_ZtSn_vv(D); 
    dcc_lf_TAnZtS_ww(D); 
    dcc_halfloop(D); 
    nl++; /* while (nl<iteration_max && D->A_cpop_j>2 && D->A_rpop_j_total>2){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  timing_dump(2);
  wkspace_printf();
}
