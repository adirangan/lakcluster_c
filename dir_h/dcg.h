
  /*     full        -->      split       
    +-----+-----++-+     +-----+-----++-+
    |AAAAA|YYYYY||T| --> |AAAAJ|YYYYY||T|
    |AAAAA|YYYYY||T| --> |AAAAJ|YYYYY||T|
    |AAAAA|YYYYY||T| --> |AAAAJ|YYYYY||T|
    |AAAAA|YYYYY||T| --> |KKKKH|UUUUU||R|
    +-----+-----++-+ --> +-----+-----++-+
    |ZZZZZ|WWWWW||S| --> |ZZZZV|WWWWW||S|
    |ZZZZZ|WWWWW||S| --> |ZZZZV|WWWWW||S|
    |ZZZZZ|WWWWW||S| --> |ZZZZV|WWWWW||S|
    |ZZZZZ|WWWWW||S| --> |ZZZZV|WWWWW||S|
    +-----+-----++-+     +-----+-----++-+
  */

struct dcg_ajdk /* for use in passing inputs */
{
  int bitj;
  int nbins; struct dcg_single **E_; struct dcg_double **F_;
  /* char QR_strategy[FNAMESIZE]; char QC_strategy[FNAMESIZE]; */
  char tmpAnchar[FNAMESIZE]; char tmpZnchar[FNAMESIZE]; char tmpAtchar[FNAMESIZE]; char tmpZtchar[FNAMESIZE];
  char tmpYnchar[FNAMESIZE]; char tmpWnchar[FNAMESIZE]; char tmpYtchar[FNAMESIZE]; char tmpWtchar[FNAMESIZE];
  char tmpTnchar[FNAMESIZE]; char tmpSnchar[FNAMESIZE]; char tmpTtchar[FNAMESIZE]; char tmpStchar[FNAMESIZE];
  int A_pcols; int Y_pcols; int A_ncols; int Y_ncols;
  double *AZ_rsum;double *YW_rsum;double *A_p;double *Y_p; double *A_ajdk;double *Y_ajdk; 
  int A_nrows_total; unsigned int A_rpop_b_total; unsigned int A_rpop_j_total;
  int Z_nrows_total; unsigned int Z_rpop_b_total; unsigned int Z_rpop_j_total;
  int A_ncols_extend;int A_cbother;
  unsigned char *A_umc_b;unsigned char *A_umc_j;unsigned char *A_umc_j_rmv;unsigned char *A_umc_j_rtn;
  int Y_ncols_extend;int Y_cbother;
  unsigned char *Y_umc_b;unsigned char *Y_umc_j;
  int T_ncols;int T_ncols_extend;
  unsigned char *T_umc_b;unsigned char *T_umc_j;
  int A_mc_length;int Y_mc_length;int T_mc_length;
  unsigned char *A_bmc_b;unsigned char *A_bmc_j;unsigned char *A_bmc_j_rmv;unsigned char *A_bmc_j_rtn;
  unsigned char *Y_bmc_b;unsigned char *Y_bmc_j;
  unsigned char *T_bmc_b;unsigned char *T_bmc_j;
  unsigned int A_cpop_b;unsigned int A_cpop_j; unsigned int A_cpop_j_rmv; unsigned int A_cpop_j_rtn;
  unsigned int Y_cpop_b;unsigned int Y_cpop_j;
  unsigned int T_cpop_b;unsigned int T_cpop_j;
  double *QC_TYnWtS_nrm;
  double *QC_TAnZtS_nrm;
  double *QC_TYnYtT_nrm;
  double *QC_TAnAtT_nrm;
  double *QC_A_svalue; int *QC_A_index_local_mc_a; int *QC_A_index_local_mc_b; int *QC_A_index_local_mc_j;
  double *QR_A_svalue; int *QR_A_index_local_mr_a; int *QR_A_index_local_mr_b; int *QR_A_index_local_mr_j; int *QR_A_index_local_nb; int *QR_A_index_global_mr_a; int *QR_A_index_global_mr_b; int *QR_A_index_global_mr_j;
  int Irem;int Ireq;
  int out_iteration; int out_xdrop_ij; int out_trace_length;
  double *out_trace; int *out_xdrop_a; int *out_xdrop_b;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  unsigned int J_rpop_b_total; unsigned int J_rpop_j_total;
  unsigned char *K_umc_b;unsigned char *K_umc_j;unsigned char *K_umc_j_rmv;unsigned char *K_umc_j_rtn;
  unsigned char *K_bmc_b;unsigned char *K_bmc_j;unsigned char *K_bmc_j_rmv;unsigned char *K_bmc_j_rtn;
  unsigned int K_cpop_b;unsigned int K_cpop_j; unsigned int K_cpop_j_rmv; unsigned int K_cpop_j_rtn;
  double *QC_RUnWtS_nrm;
  double *QC_RUnYtT_nrm;
  double *QC_RJnZtS_nrm;
  double *QC_RJnAtT_nrm;
  double *QC_TKnVtS_nrm;
  double *QC_TKnKtT_nrm;
  double *QC_K_svalue; int *QC_K_index_local_mc_a; int *QC_K_index_local_mc_b; int *QC_K_index_local_mc_j;
  double *QR_J_svalue; int *QR_J_index_local_mr_a; int *QR_J_index_local_mr_b; int *QR_J_index_local_mr_j; int *QR_J_index_local_nb; int *QR_J_index_global_mr_a; int *QR_J_index_global_mr_b; int *QR_J_index_global_mr_j;
  int out_xgrow_ij; 
  int *out_xgrow_a; int *out_xgrow_b;
};

struct dcg_single /* for use in passing inputs */
{
  int bitj; int nb; struct dcg_ajdk *D;
  long long int length;long long int length_tmp;long long int length_max;
  int A_nrows;int A_nrows_max;int A_nrows_extend; int A_rbother;
  unsigned char *A_umr_b;unsigned char *A_umr_j;unsigned char *A_umr_j_rmv;unsigned char *A_umr_j_rtn;
  int Z_nrows;int Z_nrows_max;int Z_nrows_extend; int Z_rbother;
  unsigned char *Z_umr_b;unsigned char *Z_umr_j;
  int A_mr_length;int Z_mr_length;
  unsigned char *A_bmr_b;unsigned char *A_bmr_j;unsigned char *A_bmr_j_rmv;unsigned char *A_bmr_j_rtn;
  unsigned char *Z_bmr_b;unsigned char *Z_bmr_j;
  unsigned int A_rpop_b;unsigned int A_rpop_j; unsigned int A_rpop_j_rmv; unsigned int A_rpop_j_rtn;
  unsigned int Z_rpop_b;unsigned int Z_rpop_j;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct M_handle *M_An;struct M_handle *M_At; struct M_handle *M_Yn;struct M_handle *M_Yt; struct M_handle *M_Tn;struct M_handle *M_Tt; 
  struct M_handle *M_Zn;struct M_handle *M_Zt; struct M_handle *M_Wn;struct M_handle *M_Wt; struct M_handle *M_Sn;struct M_handle *M_St;
  struct L_handle *lf_At_rsum; struct L_handle *lf_Zt_rsum; struct L_handle *lf_Yt_rsum; struct L_handle *lf_Wt_rsum; 
  struct L_handle *lf_An_ajdk;struct L_handle *lf_Zn_ajdk;struct L_handle *lf_Yn_ajdk;struct L_handle *lf_Wn_ajdk;
  struct L_handle *lf_AtTn;struct L_handle *lf_ZtSn;struct L_handle *lf_YtTn;struct L_handle *lf_WtSn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_A_a0d1;
  struct L_handle *lf_A_a2d1;
  struct L_handle *lf_Y_a0d1;
  struct L_handle *lf_Y_a2d1;
  struct L_handle *lf_a1d1_WtSn;
  struct L_handle *lf_a1d1_ZtSn;
  struct L_handle *lf_a1d1_YtTn;
  struct L_handle *lf_a1d1_AtTn;
  struct L_handle *lf_et_Tn;
  struct L_handle *lf_et_Sn;
  struct L_handle *lf_a0d1_WtSn;
  struct L_handle *lf_a0d1_ZtSn;
  struct L_handle *lf_a0d1_YtTn;
  struct L_handle *lf_a0d1_AtTn;
  struct M_handle *M_a0d1_WtSn;
  struct M_handle *M_a0d1_ZtSn;
  struct M_handle *M_a0d1_YtTn;
  struct M_handle *M_a0d1_AtTn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  double *QR_TYnWtS_nrm;
  double *QR_TAnZtS_nrm;
  double *QR_TYnYtT_nrm;
  double *QR_TAnAtT_nrm;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  int *QR_A_index_global_mr_a; int *QR_A_index_global_mr_b; int *QR_A_index_global_mr_j;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  unsigned char *J_umr_b;unsigned char *J_umr_j;unsigned char *J_umr_j_rmv;unsigned char *J_umr_j_rtn;
  unsigned char *J_bmr_b;unsigned char *J_bmr_j;unsigned char *J_bmr_j_rmv;unsigned char *J_bmr_j_rtn;
  unsigned int J_rpop_b;unsigned int J_rpop_j; unsigned int J_rpop_j_rmv; unsigned int J_rpop_j_rtn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct M_handle *M_Jn;struct M_handle *M_Jt; struct M_handle *M_Kn;struct M_handle *M_Kt; 
  struct M_handle *M_Un;struct M_handle *M_Ut; struct M_handle *M_Vn;struct M_handle *M_Vt; 
  struct M_handle *M_Rn;struct M_handle *M_Rt; 
  struct L_handle *lf_Jt_rsum; struct L_handle *lf_Kt_rsum; struct L_handle *lf_Vt_rsum; struct L_handle *lf_Ut_rsum;
  struct L_handle *lf_Jn_ajdk; struct L_handle *lf_Kn_ajdk; struct L_handle *lf_Vn_ajdk; struct L_handle *lf_Un_ajdk;
  struct L_handle *lf_JtRn;struct L_handle *lf_KtTn;struct L_handle *lf_VtSn;struct L_handle *lf_UtRn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_J_a0d1;
  struct L_handle *lf_J_a2d1;
  struct L_handle *lf_K_a0d1;
  struct L_handle *lf_K_a2d1;
  struct L_handle *lf_U_a0d1;
  struct L_handle *lf_U_a2d1;
  struct L_handle *lf_V_a0d1;
  struct L_handle *lf_V_a2d1;
  struct L_handle *lf_a1d1_JtRn;
  struct L_handle *lf_a1d1_KtTn;
  struct L_handle *lf_a1d1_VtSn;
  struct L_handle *lf_a1d1_UtRn;
  struct L_handle *lf_et_Rn;
  struct L_handle *lf_a0d1_VtSn;
  struct L_handle *lf_a0d1_UtRn;
  struct L_handle *lf_a0d1_JtRn;
  struct L_handle *lf_a0d1_KtTn;
  struct M_handle *M_a0d1_VtSn;
  struct M_handle *M_a0d1_UtRn;
  struct M_handle *M_a0d1_JtRn;
  struct M_handle *M_a0d1_KtTn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  double *QR_RUnWtS_nrm;
  double *QR_RUnYtT_nrm;
  double *QR_RJnZtS_nrm;
  double *QR_RJnAtT_nrm;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  int *QR_J_index_global_mr_a; int *QR_J_index_global_mr_b; int *QR_J_index_global_mr_j;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
};

struct dcg_double /* for use in passing inputs */
{
  int bitj; int nb1; int nb2; struct dcg_ajdk *D; struct dcg_single *E_nb1; struct dcg_single *E_nb2;
  long long int length;long long int length_tmp;long long int length_max;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_Yn_a0d1_WtSn; struct L_handle *lf_Yn_a0d1_YtTn;
  struct L_handle *lf_An_a0d1_ZtSn; struct L_handle *lf_An_a0d1_AtTn;
  struct L_handle *lf_TYnWtS_ww; struct L_handle *lf_TYnYtT_ww;
  struct L_handle *lf_TAnZtS_ww; struct L_handle *lf_TAnAtT_ww;
  struct L_handle *lf_TYnWtS_uu; struct L_handle *lf_TYnYtT_uu;
  struct L_handle *lf_TAnZtS_uu; struct L_handle *lf_TAnAtT_uu;
  struct L_handle *lf_D_YtTn_WtSn_vv; struct L_handle *lf_D_YtTn_YtTn_vv;
  struct L_handle *lf_D_AtTn_ZtSn_vv; struct L_handle *lf_D_AtTn_AtTn_vv;
  struct L_handle *lf_D_YtTn_WtSn_uu; struct L_handle *lf_D_YtTn_YtTn_uu;
  struct L_handle *lf_D_AtTn_ZtSn_uu; struct L_handle *lf_D_AtTn_AtTn_uu;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *QR_TYnWtS; struct L_handle *QR_TAnZtS; struct L_handle *QR_TYnYtT; struct L_handle *QR_TAnAtT;
  struct L_handle *QC_TYnWtS; struct L_handle *QC_TAnZtS; struct L_handle *QC_TYnYtT; struct L_handle *QC_TAnAtT;
  struct L_handle *QR_TYnWtS_uu; struct L_handle *QR_TAnZtS_uu; struct L_handle *QR_TYnYtT_uu; struct L_handle *QR_TAnAtT_uu;
  struct L_handle *QC_TYnWtS_uu; struct L_handle *QC_TAnZtS_uu; struct L_handle *QC_TYnYtT_uu; struct L_handle *QC_TAnAtT_uu;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_Un_a0d1_WtSn; struct L_handle *lf_Un_a0d1_YtTn;
  struct L_handle *lf_Jn_a0d1_ZtSn; struct L_handle *lf_Jn_a0d1_AtTn;
  struct L_handle *lf_RUnWtS_ww; struct L_handle *lf_RUnYtT_ww;
  struct L_handle *lf_RJnZtS_ww; struct L_handle *lf_RJnAtT_ww;
  struct L_handle *lf_RUnWtS_uu; struct L_handle *lf_RUnYtT_uu;
  struct L_handle *lf_RJnZtS_uu; struct L_handle *lf_RJnAtT_uu;
  struct L_handle *lf_D_UtRn_WtSn_vv; struct L_handle *lf_D_UtRn_YtTn_vv;
  struct L_handle *lf_D_JtRn_ZtSn_vv; struct L_handle *lf_D_JtRn_AtTn_vv;
  struct L_handle *lf_D_KtTn_VtSn_vv; struct L_handle *lf_D_KtTn_KtTn_vv;
  struct L_handle *lf_D_UtRn_WtSn_uu; struct L_handle *lf_D_UtRn_YtTn_uu;
  struct L_handle *lf_D_JtRn_ZtSn_uu; struct L_handle *lf_D_JtRn_AtTn_uu;
  struct L_handle *lf_D_KtTn_VtSn_uu; struct L_handle *lf_D_KtTn_KtTn_uu;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *QR_RUnWtS; struct L_handle *QR_RJnZtS; struct L_handle *QR_RUnYtT; struct L_handle *QR_RJnAtT;
  struct L_handle *QC_RUnWtS; struct L_handle *QC_RJnZtS; struct L_handle *QC_RUnYtT; struct L_handle *QC_RJnAtT;
  struct L_handle *QC_TKnVtS; struct L_handle *QC_TKnKtT; 
  struct L_handle *QR_RUnWtS_uu; struct L_handle *QR_RJnZtS_uu; struct L_handle *QR_RUnYtT_uu; struct L_handle *QR_RJnAtT_uu;
  struct L_handle *QC_RUnWtS_uu; struct L_handle *QC_RJnZtS_uu; struct L_handle *QC_RUnYtT_uu; struct L_handle *QC_RJnAtT_uu;
  struct L_handle *QC_TKnVtS_uu; struct L_handle *QC_TKnKtT_uu;
};
