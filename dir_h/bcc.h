struct bcc_ajdk /* for use in passing inputs */
{
  int bitj;
  int nbins; struct bcc_single **E_; struct bcc_double **F_;
  char QR_strategy[FNAMESIZE]; char QC_strategy[FNAMESIZE];
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
  double *QC_AtTYnWtSZn_nrm;
  double *QC_AtTAnZtSZn_nrm;
  double *QC_AtTYnYtTAn_nrm;
  double *QC_AtTAnAtTAn_nrm;
  double *QC_svalue; int *QC_index_local_mc_a; int *QC_index_local_mc_b; int *QC_index_local_mc_j;
  double *QR_svalue; int *QR_index_local_mr_a; int *QR_index_local_mr_b; int *QR_index_local_mr_j; int *QR_index_local_nb; int *QR_index_global_mr_a; int *QR_index_global_mr_b; int *QR_index_global_mr_j;
  int Irem;int Ireq;
  int out_iteration; int out_xdrop_ij; int out_trace_length;
  double *out_trace; int *out_xdrop_a; int *out_xdrop_b;
};

struct bcc_single /* for use in passing inputs */
{
  int bitj; int nb; struct bcc_ajdk *D;
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
  struct L_handle *lf_AtTAn;struct L_handle *lf_AtTYn;struct L_handle *lf_ZtSZn;struct L_handle *lf_ZtSWn;
  struct M_handle **M_AtTAn_;struct M_handle **M_AtTYn_;struct M_handle **M_ZtSZn_;struct M_handle **M_ZtSWn_;
  int M_AtTAn_trm;int M_AtTYn_trm;int M_ZtSZn_trm;int M_ZtSWn_trm;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  /*     old         -->      new        
    +-----+-----++-+     +-----+-----++-+
    |AAAAA|YYYYY||T| --> |aaaaj|yyyyy||t|
    |AAAAA|YYYYY||T| --> |aaaaj|yyyyy||t|
    |AAAAA|YYYYY||T| --> |aaaaj|yyyyy||t|
    |AAAAA|YYYYY||T| --> |kkkkh|uuuuu||r|
    +-----+-----++-+ --> +-----+-----++-+
    |ZZZZZ|WWWWW||S| --> |zzzzv|WWWWW||S|
    |ZZZZZ|WWWWW||S| --> |zzzzv|WWWWW||S|
    |ZZZZZ|WWWWW||S| --> |zzzzv|WWWWW||S|
    |ZZZZZ|WWWWW||S| --> |zzzzv|WWWWW||S|
    +-----+-----++-+     +-----+-----++-+
  */
  struct M_handle *M_an; struct M_handle *M_at; struct M_handle *M_jn; struct M_handle *M_jt; struct M_handle *M_yn; struct M_handle *M_yt; struct M_handle *M_tt;  struct M_handle *M_kn; struct M_handle *M_kt; struct M_handle *M_hn; struct M_handle *M_ht; struct M_handle *M_un; struct M_handle *M_ut; struct M_handle *M_rt;
  struct M_handle *M_zn; struct M_handle *M_zt; struct M_handle *M_vn; struct M_handle *M_vt;
  /* QR YnWt lrup */
  struct L_handle *lf_jn_ajdk;struct L_handle *lf_vn_ajdk;
  struct L_handle *lf_jn;struct L_handle *lf_vn;
  /* QC ZtSWn lrup */
  struct L_handle *lf_ktrn;struct L_handle *lf_utrn;
  struct L_handle *lf_ktrun;struct L_handle *lf_ktrkn;
  /* QC YnWt lrup */
  struct L_handle *lf_kn_ajdk; struct L_handle *lf_zn_ajdk; struct L_handle *lf_hn_ajdk; /* struct L_handle *lf_vn_ajdk; */ struct L_handle *lf_an_ajdk;
  /* QC YnWt lrup */
  struct L_handle *lf_attn; struct L_handle *lf_jttn; /* struct L_handle *lf_ktrn; */ struct L_handle *lf_htrn; struct L_handle *lf_ztsn; struct L_handle *lf_vtsn; 
  struct L_handle *lf_ztsvn; struct L_handle *lf_attjn; struct L_handle *lf_ktrhn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_At_Yn_a1d1;
  struct L_handle *lf_At_An_a1d1;
  struct L_handle *lf_a1d2_Zt_Sn;
  struct L_handle *lf_a2d2_Zt_Sn;
  struct L_handle *lf_a3d2_Zt_Sn;
  struct L_handle *lf_a1d2_At_Tn;
  struct L_handle *lf_a2d2_At_Tn;
  struct L_handle *lf_a3d2_At_Tn;
  struct L_handle *lf_et_Tn;
  struct L_handle *lf_et_Sn;
  struct L_handle *lf_et_An;
  struct L_handle *lf_a1d1_At_en;
  struct L_handle *lf_a1d2_At_en;
  struct L_handle *lf_a3d2_At_en;
  struct L_handle *lf_et_Yn_a1d1;
  struct L_handle *lf_Tt_Yn_a1d1;
  /* struct L_handle *lf_At_T_Yn_a1d1; */
  struct L_handle *lf_et_An_a1d1;
  struct L_handle *lf_Tt_An_a1d1;
  /* struct L_handle *lf_At_T_An_a1d1; */
  struct L_handle *lf_T_AnAt_YnYt;
  struct L_handle *lf_T_AnAt_AnAt;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct M_handle *M_1;
  struct L_handle *lf_Yn_a1d1;
  struct L_handle *lf_An_a1d1;
  struct M_handle *M_Yn_a1d1;
  struct M_handle *M_An_a1d1;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_a2d2_ZtSn;
  struct L_handle *lf_a2d2_AtTn;
  struct M_handle *M_a2d2_ZtSn;
  struct M_handle *M_a2d2_AtTn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  double *QR_AnZtSWnYt_nrm;
  double *QR_AnZtSZnAt_nrm;
  double *QR_AnAtTYnYt_nrm;
  double *QR_AnAtTAnAt_nrm;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  int *QR_index_global_mr_a; int *QR_index_global_mr_b; int *QR_index_global_mr_j;
};

struct bcc_double /* for use in passing inputs */
{
  int bitj; int nb1; int nb2; struct bcc_ajdk *D; struct bcc_single *E_nb1; struct bcc_single *E_nb2;
  long long int length;long long int length_tmp;long long int length_max;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_AnAt;struct L_handle *lf_AnZt;struct L_handle *lf_YnYt;struct L_handle *lf_YnWt;
  struct M_handle *M_YnWt;struct M_handle *M_AnZt;struct M_handle *M_YnYt;struct M_handle *M_AnAt;
  int M_YnWt_trm; int M_AnZt_trm; int M_YnYt_trm; int M_AnAt_trm;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  /* QR YnWt lrup */
  struct L_handle *lf_jnvt;struct L_handle *lf_jnjt;
  /* QC YnWt lrup */
  struct L_handle *lf_unwt; struct L_handle *lf_knzt; struct L_handle *lf_hnvt; struct L_handle *lf_ynut; struct L_handle *lf_unyt; struct L_handle *lf_unut; struct L_handle *lf_ankt; struct L_handle *lf_knat; struct L_handle *lf_knkt;
  struct M_handle *M_unwt; struct M_handle *M_knzt; struct M_handle *M_hnvt; struct M_handle *M_ynut; struct M_handle *M_unyt; struct M_handle *M_unut; struct M_handle *M_ankt; struct M_handle *M_knat; struct M_handle *M_knkt;
  int M_unwt_trm; int M_knzt_trm; int M_hnvt_trm; int M_ynut_trm; int M_unyt_trm; int M_unut_trm; int M_ankt_trm; int M_knat_trm; int M_knkt_trm;
  struct L_handle *lf_kt_r_unwt_s_zn; struct L_handle *lf_kt_r_knzt_s_zn; struct L_handle *lf_kt_r_hnvt_s_zn; struct L_handle *lf_at_t_ynut_r_kn; struct L_handle *lf_kt_r_unyt_t_an; struct L_handle *lf_kt_r_unut_r_kn; struct L_handle *lf_at_t_ankt_r_kn; struct L_handle *lf_kt_r_knat_t_an; struct L_handle *lf_kt_r_knkt_r_kn;
  struct L_handle *lf_attjn____vtszn; struct L_handle *lf_attjn____jttan; struct L_handle *lf_attjn____htrkn; struct L_handle *lf_ktrhn____jttan; struct L_handle *lf_ktrhn____htrkn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_An_a2d2_Zt_Sn;
  struct L_handle *lf_An_a2d2_At_Tn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_AnZtSWnYt; struct L_handle *lf_AnZtSZnAt; struct L_handle *lf_AnAtTYnYt; struct L_handle *lf_AnAtTAnAt;
  struct L_handle *lf_AnZt_S_WnYt; struct L_handle *lf_AnZt_S_ZnAt; struct L_handle *lf_AnAt_T_YnYt; struct L_handle *lf_AnAt_T_AnAt;
  struct L_handle *lf_An_ZtSWn_Yt; struct L_handle *lf_An_ZtSZn_At; struct L_handle *lf_An_AtTYn_Yt; struct L_handle *lf_An_AtTAn_At;
  struct L_handle *lf_AtTYnWtSZn; struct L_handle *lf_AtTAnZtSZn; struct L_handle *lf_AtTYnYtTAn; struct L_handle *lf_AtTAnAtTAn;
  struct L_handle *lf_AtTYn____WtSZn; struct L_handle *lf_AtTAn____ZtSZn; struct L_handle *lf_AtTYn____YtTAn; struct L_handle *lf_AtTAn____AtTAn;
  struct L_handle *lf_At_T_YnWt_S_Zn; struct L_handle *lf_At_T_AnZt_S_Zn; struct L_handle *lf_At_T_YnYt_T_An; struct L_handle *lf_At_T_AnAt_T_An;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *QR_AnZtSWnYt; struct L_handle *QR_AnZtSZnAt; struct L_handle *QR_AnAtTYnYt; struct L_handle *QR_AnAtTAnAt;
  struct L_handle *QC_AtTYnWtSZn; struct L_handle *QC_AtTAnZtSZn; struct L_handle *QC_AtTYnYtTAn; struct L_handle *QC_AtTAnAtTAn;
  struct L_handle *QR_AnZtSWnYt_uu; struct L_handle *QR_AnZtSZnAt_uu; struct L_handle *QR_AnAtTYnYt_uu; struct L_handle *QR_AnAtTAnAt_uu;
  struct L_handle *QC_AtTYnWtSZn_uu; struct L_handle *QC_AtTAnZtSZn_uu; struct L_handle *QC_AtTYnYtTAn_uu; struct L_handle *QC_AtTAnAtTAn_uu;
};
