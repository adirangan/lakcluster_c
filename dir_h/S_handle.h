struct S_handle
{
  char out_xdrop_name[FNAMESIZE]; /* name of out_xdrop file */
  int out_xdrop_nlines; /* number of lines in out_xdrop */
  int out_xdrop_nrows; /* number of row-indices stored in out_xdrop */
  int out_xdrop_ncols; /* number of col-indices stored in out_xdrop */
  int *out_xdrop; /* actual array for out_xdrop */
  int *mr_index_sort; /* row-indices from out_xdrop */
  int *mc_index_sort; /* col-indices from out_xdrop */
  int *mr_index_local_nb; /* local nb from row-indices */
  int *mr_index_local_mr; /* local mr from row-indices */
  int row_max; /* maximum row index for score calculation */
  int row_num; /* number of row indices for score calculation */
  int col_max; /* maximum col index for score calculation */
  int col_num; /* number of col indices for score calculation */
  int *rkeep; /* array of r-values (rows remaining) */
  int *rdrop; /* array of rdrop-values */
  int *ckeep; /* array of c-values (cols remaining) */
  int *cdrop; /* array of cdrop-values */
  int *A_rpop_j_total; /* output array (size row_num-x-col_num) */
  int *A_cpop_j; /* output array (size row_num-x-col_num) */
  int *Irem; /* output array (size row_num-x-col_num) */
  double *QR_avg; /* output array (size row_num-x-col_num) */
  double *QC_avg; /* output array (size row_num-x-col_num) */
  int nr; /* temporary row index */
  int nc; /* temporary row index */
};
