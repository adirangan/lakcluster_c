void mda_write_r8(char *fname,int nrows,int ncols,double *ra)
{
  FILE *fp=NULL; int tmp_d=0;
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",fname); exit(RET_READ_FAIL);}
  tmp_d = 2; fwrite(&tmp_d,sizeof(int),1,fp); tmp_d = nrows; fwrite(&tmp_d,sizeof(int),1,fp); tmp_d = ncols; fwrite(&tmp_d,sizeof(int),1,fp); fwrite(ra,sizeof(double),nrows*ncols,fp);
  fclose(fp);fp=NULL;
}

void mda_write_i4(char *fname,int nrows,int ncols,int *ra)
{
  FILE *fp=NULL; int tmp_d=0;
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",fname); exit(RET_READ_FAIL);}
  tmp_d = 2; fwrite(&tmp_d,sizeof(int),1,fp); tmp_d = nrows; fwrite(&tmp_d,sizeof(int),1,fp); tmp_d = ncols; fwrite(&tmp_d,sizeof(int),1,fp); fwrite(ra,sizeof(int),nrows*ncols,fp);
  fclose(fp);fp=NULL;
}
