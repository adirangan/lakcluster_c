#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void mda_write_d2_r8(char *fname,int nrows,int ncols,double *dra)
{
  FILE *fp=NULL; int tmp_d=0;
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",fname); exit(RET_READ_FAIL);}
  tmp_d = 2; fwrite(&tmp_d,sizeof(int),1,fp); tmp_d = nrows; fwrite(&tmp_d,sizeof(int),1,fp); tmp_d = ncols; fwrite(&tmp_d,sizeof(int),1,fp); fwrite(dra,sizeof(double),nrows*ncols,fp);
  fclose(fp);fp=NULL;
}

void mda_write_d3_r8(char *fname,int nrows,int ncols,int nlyrs,double *dra)
{
  FILE *fp=NULL; int tmp_d=0;
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",fname); exit(RET_READ_FAIL);}
  tmp_d = 3; fwrite(&tmp_d,sizeof(int),1,fp); 
  tmp_d = nrows; fwrite(&tmp_d,sizeof(int),1,fp); 
  tmp_d = ncols; fwrite(&tmp_d,sizeof(int),1,fp); 
  tmp_d = nlyrs; fwrite(&tmp_d,sizeof(int),1,fp); 
  fwrite(dra,sizeof(double),nrows*ncols*nlyrs,fp);
  fclose(fp);fp=NULL;
}

void mda_write_d2_i4(char *fname,int nrows,int ncols,int *ira)
{
  FILE *fp=NULL; int tmp_d=0;
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",fname); exit(RET_READ_FAIL);}
  tmp_d = 2; fwrite(&tmp_d,sizeof(int),1,fp); tmp_d = nrows; fwrite(&tmp_d,sizeof(int),1,fp); tmp_d = ncols; fwrite(&tmp_d,sizeof(int),1,fp); fwrite(ira,sizeof(int),nrows*ncols,fp);
  fclose(fp);fp=NULL;
}

void mda_read_d1_r8(char *fname,int *nrows_p,double *dra)
{
  FILE *fp=NULL; int tmp_d=0; int nrows=0;
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s when reading from disc.\n",fname); exit(RET_READ_FAIL);}
  fread(&tmp_d,sizeof(int),1,fp); if (tmp_d!=1){ printf(" %% Warning! tmp_d %d in mda_read_d1_r8\n",tmp_d);}
  fread(&tmp_d,sizeof(int),1,fp); nrows = tmp_d;
  if (nrows_p!=NULL){ *nrows_p = nrows;}
  if (dra!=NULL){ fread(dra,sizeof(double),nrows,fp);}
  fclose(fp);fp=NULL;
}

void mda_read_d2_r8(char *fname,int *nrows_p,int *ncols_p,double *dra)
{
  FILE *fp=NULL; int tmp_d=0; int nrows=0,ncols=0;
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s when reading from disc.\n",fname); exit(RET_READ_FAIL);}
  fread(&tmp_d,sizeof(int),1,fp); if (tmp_d!=2){ printf(" %% Warning! tmp_d %d in mda_read_d2_r8\n",tmp_d);}
  fread(&tmp_d,sizeof(int),1,fp); nrows = tmp_d;
  fread(&tmp_d,sizeof(int),1,fp); ncols = tmp_d;
  if (nrows_p!=NULL){ *nrows_p = nrows;}
  if (ncols_p!=NULL){ *ncols_p = ncols;}
  if (dra!=NULL){ fread(dra,sizeof(double),nrows*ncols,fp);}
  fclose(fp);fp=NULL;
}

void mda_read_d3_r8(char *fname,int *nrows_p,int *ncols_p,int *nlyrs_p,double *dra)
{
  FILE *fp=NULL; int tmp_d=0; int nrows=0,ncols=0,nlyrs=0;
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s when reading from disc.\n",fname); exit(RET_READ_FAIL);}
  fread(&tmp_d,sizeof(int),1,fp); if (tmp_d!=3){ printf(" %% Warning! tmp_d %d in mda_read_d3_r8\n",tmp_d);}
  fread(&tmp_d,sizeof(int),1,fp); nrows = tmp_d;
  fread(&tmp_d,sizeof(int),1,fp); ncols = tmp_d;
  fread(&tmp_d,sizeof(int),1,fp); nlyrs = tmp_d;
  if (nrows_p!=NULL){ *nrows_p = nrows;}
  if (ncols_p!=NULL){ *ncols_p = ncols;}
  if (nlyrs_p!=NULL){ *nlyrs_p = nlyrs;}
  if (dra!=NULL){ fread(dra,sizeof(double),nrows*ncols*nlyrs,fp);}
  fclose(fp);fp=NULL;
}
