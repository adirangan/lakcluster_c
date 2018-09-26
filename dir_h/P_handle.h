struct P_handle
{
  char infix[FNAMESIZE]; /* label */
  char out_xdrop_name[FNAMESIZE]; /* name of out_xdrop file */
  int out_xdrop_nlines; /* number of lines in out_xdrop */
  int out_xdrop_nrows; /* number of row-indices stored in out_xdrop */
  int out_xdrop_ncols; /* number of col-indices stored in out_xdrop */
  int *out_xdrop; /* actual array for out_xdrop */
  int *mr_index_sort; /* row-indices from out_xdrop */
  int *mc_index_sort; /* col-indices from out_xdrop */
  int *mr_index_local_nb; /* local nb from row-indices */
  int *mr_index_local_mr; /* local mr from row-indices */
  int rank; /* rank of pca */
  double tolerance; /* tolerance used to assess relative-error in power-iteration */
  int iteration_num; /* number of iterations */
  int iteration_min; /* maximum iteration */
  int iteration_max; /* minimum iteration */
  int *rkeep; /* array of r-values (rows remaining) */
  int *rdrop; /* array of rdrop-values */
  int *ckeep; /* array of c-values (cols remaining) */
  int *cdrop; /* array of cdrop-values */
  int *A_rpop_j_total; /* output array (size row_num-x-col_num) */
  int *A_cpop_j; /* output array (size row_num-x-col_num) */
  int *Irem; /* output array (size row_num-x-col_num) */
  struct dcc_ajdk *D;
  double *U_; /* left singular-vectors (size P->D->A_nrows_total-x-P->rank-x-P->iteration_num) */
  double *V_; /* right singular_vectors (size P->D->A_ncols-x-P->rank-x-P->iteration_num) */
  char V_name[FNAMESIZE]; /* name of V_ file */
  double *AnV_; /* projection of V_ onto An */
  double *ZnV_; /* projection of V_ onto Zn */
  struct L_handle **lf_U_; struct M_handle **M_U_;
  struct L_handle *lf_V; struct M_handle *M_V; struct L_handle *lf_R1; struct L_handle *lf_R2; struct M_handle *M_rank;
  struct L_handle **lf_AnV_; struct M_handle **M_AnV_;
  struct L_handle **lf_ZnV_; struct M_handle **M_ZnV_;
  struct L_handle **lf_AtAnV_; struct M_handle **M_AtAnV_; struct L_handle *lf_AtAnV; struct M_handle *M_AtAnV;
  struct L_handle **lf_YtAnV_; struct M_handle **M_YtAnV_; struct L_handle *lf_YtAnV; struct M_handle *M_YtAnV;
  struct L_handle **lf_ZtZnV_; struct M_handle **M_ZtZnV_; struct L_handle *lf_ZtZnV; struct M_handle *M_ZtZnV;
  struct L_handle **lf_WtZnV_; struct M_handle **M_WtZnV_; struct L_handle *lf_WtZnV; struct M_handle *M_WtZnV;
  struct L_handle *lf_VA; struct M_handle *M_VA;
  struct L_handle **lf_AnVA_; struct M_handle **M_AnVA_;
  struct L_handle **lf_ZnVA_; struct M_handle **M_ZnVA_;
  struct L_handle **lf_AtAnVA_; struct M_handle **M_AtAnVA_; struct L_handle *lf_AtAnVA; struct M_handle *M_AtAnVA;
  struct L_handle **lf_ZtZnVA_; struct M_handle **M_ZtZnVA_; struct L_handle *lf_ZtZnVA; struct M_handle *M_ZtZnVA;
  struct L_handle *lf_VY; struct M_handle *M_VY;
  struct L_handle **lf_YnVY_; struct M_handle **M_YnVY_;
  struct L_handle **lf_WnVY_; struct M_handle **M_WnVY_;
  struct L_handle **lf_AtYnVY_; struct M_handle **M_AtYnVY_; struct L_handle *lf_AtYnVY; struct M_handle *M_AtYnVY;
  struct L_handle **lf_ZtWnVY_; struct M_handle **M_ZtWnVY_; struct L_handle *lf_ZtWnVY; struct M_handle *M_ZtWnVY;  
  struct L_handle *lf_VAZ; struct M_handle *M_VAZ;
  struct L_handle *lf_VYW; struct M_handle *M_VYW;
  int b_mlt; /* precision used for xcalc */
  int nx; /* temporary row/col index */
  /* used for error checking */
  struct L_handle *lf_At; struct L_handle *lf_Zt; struct L_handle *lf_Yt; struct L_handle *lf_Wt;
  struct L_handle *lf_An; struct L_handle *lf_Zn; struct L_handle *lf_Yn; struct L_handle *lf_Wn;
  struct L_handle *lf_AtAn; struct L_handle *lf_ZtZn; struct L_handle *lf_AtYn; struct L_handle *lf_ZtWn;
  struct L_handle *lf_AtAnAtAn; struct L_handle *lf_AtAnZtZn; struct L_handle *lf_ZtZnAtAn; struct L_handle *lf_ZtZnZtZn;
  struct L_handle *lf_AtYnYtAn; struct L_handle *lf_AtYnWtZn; struct L_handle *lf_ZtWnYtAn; struct L_handle *lf_ZtWnWtZn;
  struct L_handle *lf_S; struct L_handle *lf_Vn; struct L_handle *lf_Vt; struct L_handle *lf_Vx; 
  struct L_handle *lf_VtAt; struct L_handle *lf_VtZt;
  struct L_handle *lf_VtAtAn; struct L_handle *lf_VtZtZn; struct L_handle *lf_VtAtYn; struct L_handle *lf_VtZtWn;
  struct L_handle *lf_VAt;
  struct L_handle *lf_VAtAt; struct L_handle *lf_VAtZt;
  struct L_handle *lf_VAtAtAn; struct L_handle *lf_VAtZtZn;
  struct L_handle *lf_VAZt;
  struct L_handle *lf_VYt;
  struct L_handle *lf_VYtYt; struct L_handle *lf_VYtWt;
  struct L_handle *lf_VYtYtAn; struct L_handle *lf_VYtWtZn;
  struct L_handle *lf_VYWt;
};
