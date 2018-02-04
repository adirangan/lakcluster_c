struct dcc_ajdk /* for use in passing inputs */
{
  int bitj;
  int nbins; struct dcc_single **E_; struct dcc_double **F_;
  char tmpAnchar[FNAMESIZE]; char tmpZnchar[FNAMESIZE]; char tmpAtchar[FNAMESIZE]; char tmpZtchar[FNAMESIZE];
  char tmpTnchar[FNAMESIZE]; char tmpSnchar[FNAMESIZE]; char tmpTtchar[FNAMESIZE]; char tmpStchar[FNAMESIZE];
  int A_pcols; int A_ncols;
  double *AZ_rsum;double *A_p;double *A_ajdk; 
  int A_nrows_total; unsigned int A_rpop_b_total; unsigned int A_rpop_j_total;
  int Z_nrows_total; unsigned int Z_rpop_b_total; unsigned int Z_rpop_j_total;
  int A_ncols_extend;int A_cbother;
  unsigned char *A_umc_b;unsigned char *A_umc_j;unsigned char *A_umc_j_rmv;unsigned char *A_umc_j_rtn;
  int T_ncols;int T_ncols_extend;
  unsigned char *T_umc_b;unsigned char *T_umc_j;
  int A_mc_length;int T_mc_length;
  unsigned char *A_bmc_b;unsigned char *A_bmc_j;unsigned char *A_bmc_j_rmv;unsigned char *A_bmc_j_rtn;
  unsigned char *T_bmc_b;unsigned char *T_bmc_j;
  unsigned int A_cpop_b;unsigned int A_cpop_j; unsigned int A_cpop_j_rmv; unsigned int A_cpop_j_rtn;
  unsigned int T_cpop_b;unsigned int T_cpop_j;
  double *QC_TAnZtS_nrm;
  double *QC_TAnAtT_nrm;
  double *QC_sra; int *QC_lmc_a; int *QC_lmc_b; int *QC_lmc_j;
  double *QR_sra; int *QR_lmr_a; int *QR_lmr_b; int *QR_lmr_j; int *QR_lnb; int *QR_imr_a; int *QR_imr_b; int *QR_imr_j;
  int Irem;int Ireq;
  int out_iteration; int out_xdrop_ij; int out_trace_length;
  double *out_trace; int *out_xdrop_a; int *out_xdrop_b;
};

struct dcc_single /* for use in passing inputs */
{
  int bitj; int nb; struct dcc_ajdk *D;
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
  struct M_handle *M_An;struct M_handle *M_At; struct M_handle *M_Tn;struct M_handle *M_Tt; 
  struct M_handle *M_Zn;struct M_handle *M_Zt; struct M_handle *M_Sn;struct M_handle *M_St;
  struct L_handle *lf_At_rsum; struct L_handle *lf_Zt_rsum;
  struct L_handle *lf_An_ajdk;struct L_handle *lf_Zn_ajdk;
  struct L_handle *lf_AtTn;struct L_handle *lf_ZtSn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_a0d1;
  struct L_handle *lf_a2d1;
  struct L_handle *lf_a1d1_ZtSn;
  struct L_handle *lf_a1d1_AtTn;
  struct L_handle *lf_et_Tn;
  struct L_handle *lf_et_Sn;
  struct L_handle *lf_a0d1_ZtSn;
  struct L_handle *lf_a0d1_AtTn;
  struct M_handle *M_a0d1_ZtSn;
  struct M_handle *M_a0d1_AtTn;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  double *QR_TAnZtS_nrm;
  double *QR_TAnAtT_nrm;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  int *QR_imr_a; int *QR_imr_b; int *QR_imr_j;
};

struct dcc_double /* for use in passing inputs */
{
  int bitj; int nb1; int nb2; struct dcc_ajdk *D; struct dcc_single *E_nb1; struct dcc_single *E_nb2;
  long long int length;long long int length_tmp;long long int length_max;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *lf_An_a0d1_ZtSn; struct L_handle *lf_An_a0d1_AtTn;
  struct L_handle *lf_TAnZtS_ww; struct L_handle *lf_TAnAtT_ww;
  struct L_handle *lf_TAnZtS_uu; struct L_handle *lf_TAnAtT_uu;
  struct L_handle *lf_D_AtTn_ZtSn_vv; struct L_handle *lf_D_AtTn_AtTn_vv;
  struct L_handle *lf_D_AtTn_ZtSn_uu; struct L_handle *lf_D_AtTn_AtTn_uu;
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  struct L_handle *QR_TAnZtS; struct L_handle *QR_TAnAtT;
  struct L_handle *QC_TAnZtS; struct L_handle *QC_TAnAtT;
  struct L_handle *QR_TAnZtS_uu; struct L_handle *QR_TAnAtT_uu;
  struct L_handle *QC_TAnZtS_uu; struct L_handle *QC_TAnAtT_uu;
};
