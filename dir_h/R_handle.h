struct R_handle
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
  unsigned long int rseed; /* seed to use when scrambling */
  int nr; /* temporary row index */
  int nc; /* temporary row index */
};
