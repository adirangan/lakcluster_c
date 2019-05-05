struct M_handle
{
  unsigned char *mark;
  int header_length; /* header length in bytes expected for array data file */
  char A_name[FNAMESIZE]; /* name of array data file */ 
  FILE *A_fp; /* stores array data file pointer */ unsigned char *Ara; /* stores array data in memory */
  int bitj; /* set to 16, this is the chunk-size used when storing rows and columns */
  int nrows; int ncols; /* dimensions of A */
  int nrows_extend; int ncols_extend; /* additional rows and columns required to reach a multiple of bitj */
  long long int rpop_b; long long int rpop_j; /* popcounts for mr_b and mr_j, respectively */
  long long int cpop_b; long long int cpop_j; /* popcounts for mc_b and mc_j, respectively */
  unsigned char *wX; /* holds data in memory, spaced by POPLENGTH */
  unsigned char *mr_b; unsigned char *mr_j; /* bitmasks for the rows; original and current */
  unsigned char *mc_b; unsigned char *mc_j; /* bitmasks for the cols; original and current */
  double *rsum; /* holds the row-sums of A */
  unsigned int *m_a_; /* m_a_[nr_j] = nr_a */ /* position nr_a of mr_b,mr_j corresponds to nr_j^th positive entry of mr_j*/
  unsigned int *m_b_; /* m_b_[nr_j] = nr_b */ /* and nr_b^th positive entry of mr_b */  
  unsigned int *n_a_; /* n_a_[nc_j] = nc_a */ /* position nc_a of mc_b,mc_j corresponds to nc_j^th positive entry of mc_j*/
  unsigned int *n_b_; /* n_b_[nc_j] = nc_b */ /* and nc_b^th positive entry of mc_b */  
  unsigned long mr_length; unsigned long mc_length; /* bsize(nrows) and bsize(ncols), respectively */
  unsigned long cc_length; /* (ncols + ncols_extend)/BIT8 */
  unsigned long wt_length; /* length of tail required to cap lengthof M_handle to a multiple of POPLENGTH */
  unsigned char *wt; /* tail (not used) */
  unsigned long long int length_total; /* total length of M_handle */
  int ncols_per_z; /* for use with xcalc_dra_to_M_x_bitbuffer */
  double max_d; double min_d; double mlt_d; /* for use with xcalc_dra_to_M_x_bitbuffer */
};
