#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

struct S_handle *S_handle_make(char *fname,int row_max,int row_num,int col_max,int col_num)
{
  struct S_handle *S=NULL;
  S = (struct S_handle *)wkspace_all0c(sizeof(struct S_handle)*1);
  sprintf(S->out_xdrop_name,fname);
  S->out_xdrop_nlines=0; /* number of lines in out_xdrop */
  S->out_xdrop_nrows=0; /* number of row-indices stored in out_xdrop */
  S->out_xdrop_ncols=0; /* number of col-indices stored in out_xdrop */
  S->out_xdrop=NULL; /* actual array for out_xdrop */
  S->mr_index_sort=NULL; /* row-indices from out_xdrop */
  S->mc_index_sort=NULL; /* col-indices from out_xdrop */
  S->mr_index_local_nb=NULL; /* local nb from row-indices */
  S->mr_index_local_mr=NULL; /* local mr from row-indices */
  S->row_max=row_max; /* maximum row index for score calculation */
  S->row_num=row_num; /* number of row indices for score calculation */
  S->col_max=col_max; /* maximum col index for score calculation */
  S->col_num=col_num; /* number of col indices for score calculation */
  S->rkeep=NULL; /* array of r-values (rows remaining) */
  S->rdrop=NULL; /* array of rdrop-values */
  S->ckeep=NULL; /* array of c-values (cols remaining) */
  S->cdrop=NULL; /* array of cdrop-values */
  S->A_rpop_j_total=NULL; /* output array (size row_num-x-col_num) */
  S->A_cpop_j=NULL; /* output array (size row_num-x-col_num) */
  S->Irem=NULL; /* output array (size row_num-x-col_num) */
  S->QR_avg=NULL; /* output array (size row_num-x-col_num) */
  S->QC_avg=NULL; /* output array (size row_num-x-col_num) */
  xcc_out_xdrop_load(S->out_xdrop_name,&(S->out_xdrop_nlines),&(S->out_xdrop));
  xcc_out_xdrop_sort(S->out_xdrop_nlines,S->out_xdrop,&(S->out_xdrop_nrows),&(S->mr_index_sort),&(S->out_xdrop_ncols),&(S->mc_index_sort));
  S->A_rpop_j_total = (int *) wkspace_all0c(S->row_num*S->col_num*sizeof(int));
  S->A_cpop_j = (int *) wkspace_all0c(S->row_num*S->col_num*sizeof(int));
  S->Irem = (int *) wkspace_all0c(S->row_num*S->col_num*sizeof(int));
  S->QR_avg = (double *) wkspace_all0c(S->row_num*S->col_num*sizeof(double));
  S->QC_avg = (double *) wkspace_all0c(S->row_num*S->col_num*sizeof(double));
  S->nr=0;
  S->nc=0;
  return S;
}


void S_init_bcc(struct bcc_ajdk *D,struct S_handle *S)
{
  bcc_out_xdrop_lkp(D,S->out_xdrop_nrows,S->mr_index_sort,&(S->mr_index_local_nb),&(S->mr_index_local_mr));
  xcc_scorebox_xdrop_array(D->A_cpop_j,S->col_max,S->col_num,&(S->ckeep),&(S->cdrop));
  xcc_scorebox_xdrop_array(D->A_rpop_j_total,S->row_max,S->row_num,&(S->rkeep),&(S->rdrop));
}

void S_init_dcc(struct dcc_ajdk *D,struct S_handle *S)
{
  dcc_out_xdrop_lkp(D,S->out_xdrop_nrows,S->mr_index_sort,&(S->mr_index_local_nb),&(S->mr_index_local_mr));
  xcc_scorebox_xdrop_array(D->A_cpop_j,S->col_max,S->col_num,&(S->ckeep),&(S->cdrop));
  xcc_scorebox_xdrop_array(D->A_rpop_j_total,S->row_max,S->row_num,&(S->rkeep),&(S->rdrop));
}

void S_handle_printf(int verbose,struct S_handle *S,char *prefix)
{
  int nl=0;
  char tmpchar[FNAMESIZE];
  if (verbose>0){ printf("%s nr %d nc %d\n",prefix,S->nr,S->nc);}
  if (verbose>1){ printf("%s out_xdrop_name: %s\n",prefix,S->out_xdrop_name);}
  if (verbose>1){ printf("%s nlines %d, nrows %d, ncols %d\n",prefix,S->out_xdrop_nlines,S->out_xdrop_nrows,S->out_xdrop_ncols);}
  if (verbose>1){ printf("%s row_max %d row_num %d col_max %d col_num %d\n",prefix,S->row_max,S->row_num,S->col_max,S->col_num);}
  if (verbose>2){
    if (S->rkeep!=NULL){ sprintf(tmpchar,"%s rkeep:",prefix); raprintf(S->rkeep,"int",1,S->row_num,tmpchar);}
    if (S->rdrop!=NULL){ sprintf(tmpchar,"%s rdrop:",prefix); raprintf(S->rdrop,"int",1,S->row_num,tmpchar);}
    if (S->ckeep!=NULL){ sprintf(tmpchar,"%s ckeep:",prefix); raprintf(S->ckeep,"int",1,S->col_num,tmpchar);}
    if (S->cdrop!=NULL){ sprintf(tmpchar,"%s cdrop:",prefix); raprintf(S->cdrop,"int",1,S->col_num,tmpchar);}
    if (S->A_rpop_j_total!=NULL){ sprintf(tmpchar,"%s A_rpop_j_total:",prefix); raprintf(S->A_rpop_j_total,"int",S->row_num,S->col_num,tmpchar);}
    if (S->A_cpop_j!=NULL){ sprintf(tmpchar,"%s A_cpop_j:",prefix); raprintf(S->A_cpop_j,"int",S->row_num,S->col_num,tmpchar);}
    if (S->Irem!=NULL){ sprintf(tmpchar,"%s Irem:",prefix); raprintf(S->Irem,"int",S->row_num,S->col_num,tmpchar);}
    if (S->QR_avg!=NULL){ sprintf(tmpchar,"%s QR_avg:",prefix); raprintf(S->QR_avg,"double",S->row_num,S->col_num,tmpchar);}
    if (S->QC_avg!=NULL){ sprintf(tmpchar,"%s QC_avg:",prefix); raprintf(S->QC_avg,"double",S->row_num,S->col_num,tmpchar);}
    /* if (verbose>1){ } */}
}

void S_handle_dmp(struct S_handle *S)
{
  char prefix[FNAMESIZE];
  char fname[FNAMESIZE];
  sprintf(prefix,"scorebox_r%dn%d_c%dn%d",S->row_max,S->row_num,S->col_max,S->col_num);
  sprintf(fname,"%s/%s_A_rpop_j_total.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_i4(fname,S->row_num,S->col_num,S->A_rpop_j_total);
  sprintf(fname,"%s/%s_A_cpop_j.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_i4(fname,S->row_num,S->col_num,S->A_cpop_j);
  sprintf(fname,"%s/%s_Irem.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_i4(fname,S->row_num,S->col_num,S->Irem);
  sprintf(fname,"%s/%s_QR_avg.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_r8(fname,S->row_num,S->col_num,S->QR_avg);
  sprintf(fname,"%s/%s_QC_avg.mda",GLOBAL_DIR_NAME,prefix); mda_write_d2_r8(fname,S->row_num,S->col_num,S->QC_avg);
}
