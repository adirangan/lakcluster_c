struct L_handle
{
  unsigned long long int length; /* total length of lf -- often of the form row_stride*col_stride*lyr_stride */
  double *lf; /* holds data in memory -- often of the form lf[nr + nc*row_stride + nl*row_stride*col_stride] */
  int spacing_row; /* spacing of row-index */
  unsigned long long int row_stride; /* associated stride */
  int spacing_col; /* spacing of column-index */
  unsigned long long int col_stride; /* associated stride */
  int spacing_lyr; /* spacing of layer-index */
  unsigned long long int lyr_stride; /* associated stride -- often unused */
};
