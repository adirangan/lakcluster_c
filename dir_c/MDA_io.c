#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void MDA_io_test()
{
  int verbose=1;
  int n_d_0 = 2;
  int d_0_[n_d_0];
  d_0_[0] = 2; d_0_[1] = 3;
  int n_i=0;
  int nd=0;
  n_i = 1; for (nd=0;nd<n_d_0;nd++){ n_i*=d_0_[nd];}
  int i4_0_[n_i];
  i4_0_[0] = 0; i4_0_[2] = 10; i4_0_[4] = 100;
  i4_0_[1] = 1; i4_0_[3] = 20; i4_0_[5] = 200;
  int *i4_1_=NULL;
  int *d_1_=NULL;
  int n_d_1=NULL;
  char fname[32];
  sprintf(fname,"MDA_io.test");
  raprintf(i4_0_,"int",d_0_[0],d_0_[1]," %% i4_0_: ");
  MDA_write_i4(n_d_0,d_0_,i4_0_,fname);
  MDA_read_i4(&n_d_1,&d_1_,&i4_1_,fname);
  raprintf(i4_1_,"int",d_1_[0],d_1_[1]," %% i4_1_: ");
  double r8_0_[n_i];
  r8_0_[0] = 0.1; r8_0_[2] = 10.1; r8_0_[4] = 100.1;
  r8_0_[1] = 1.1; r8_0_[3] = 20.1; r8_0_[5] = 200.1;
  int n_d_2=NULL;
  int *r8_2_=NULL;
  int *d_2_=NULL;
  raprintf(r8_0_,"double",d_0_[0],d_0_[1]," %% r8_0_: ");
  MDA_write_r8(n_d_0,d_0_,r8_0_,fname);
  MDA_read_r8(&n_d_2,&d_2_,&r8_2_,fname);
  raprintf(r8_2_,"double",d_2_[0],d_2_[1]," %% r8_2_: ");
  unsigned long long int ulli_0_[n_i];
  ulli_0_[0] = 100l; ulli_0_[2] = 1010l; ulli_0_[4] = 10100l;
  ulli_0_[1] = 101l; ulli_0_[3] = 1020l; ulli_0_[5] = 10200l;
  int n_d_3=NULL;
  int *ulli_3_=NULL;
  int *d_3_=NULL;
  raprintf(ulli_0_,"unsigned long long int",d_0_[0],d_0_[1]," %% ulli_0_: ");
  MDA_write_ulli(n_d_0,d_0_,ulli_0_,fname);
  MDA_read_ulli(&n_d_3,&d_3_,&ulli_3_,fname);
  raprintf(ulli_3_,"unsigned long long int",d_3_[0],d_3_[1]," %% ulli_3_: ");
  wkspace_printf();
}

void MDA_write_i4(int n_d,int *d_,int *i4_,char *fname)
{
  int verbose=1;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int nd=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(&n_d,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(d_,sizeof(int),n_d,fp); if (s!=n_d){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  for (nd=0;nd<n_d;nd++){ n_i*=d_[nd];}
  s=fwrite(i4_,sizeof(int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}

void MDA_read_i4(int *n_d_p,int **d_p,int **i4_p,char *fname)
{
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int nd=0;
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fread(n_d_p,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  *d_p = (int *)wkspace_all0c(sizeof(int)*(*n_d_p));
  s=fread(*d_p,sizeof(int),*n_d_p,fp); if (s!=*n_d_p){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  for (nd=0;nd<*n_d_p;nd++){ n_i*=(*d_p)[nd];}
  *i4_p = (int *)wkspace_all0c(sizeof(int)*n_i);
  s=fread(*i4_p,sizeof(int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}

void MDA_write_r8(int n_d,int *d_,double *r8_,char *fname)
{
  int verbose=1;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int nd=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(&n_d,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(d_,sizeof(int),n_d,fp); if (s!=n_d){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  for (nd=0;nd<n_d;nd++){ n_i*=d_[nd];}
  s=fwrite(r8_,sizeof(double),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}

void MDA_read_r8(int *n_d_p,int **d_p,double **r8_p,char *fname)
{
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int nd=0;
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fread(n_d_p,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  *d_p = (int *)wkspace_all0c(sizeof(int)*(*n_d_p));
  s=fread(*d_p,sizeof(int),*n_d_p,fp); if (s!=*n_d_p){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  for (nd=0;nd<*n_d_p;nd++){ n_i*=(*d_p)[nd];}
  *r8_p = (double *)wkspace_all0c(sizeof(double)*n_i);
  s=fread(*r8_p,sizeof(double),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}

void MDA_write_ulli(int n_d,int *d_,unsigned long long int *ulli_,char *fname)
{
  int verbose=1;
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int nd=0;
  if (verbose){ printf(" %% writing %s\n",fname);}
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(&n_d,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  s=fwrite(d_,sizeof(int),n_d,fp); if (s!=n_d){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  for (nd=0;nd<n_d;nd++){ n_i*=d_[nd];}
  s=fwrite(ulli_,sizeof(unsigned long long int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not write to %s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}

void MDA_read_ulli(int *n_d_p,int **d_p,unsigned long long int **ulli_p,char *fname)
{
  FILE *fp=NULL;
  size_t n_i=1,s=0;
  int nd=0;
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s.\n",fname); exit(RET_READ_FAIL);}
  s=fread(n_d_p,sizeof(int),1,fp); if (s!=1){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  *d_p = (int *)wkspace_all0c(sizeof(int)*(*n_d_p));
  s=fread(*d_p,sizeof(int),*n_d_p,fp); if (s!=*n_d_p){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  for (nd=0;nd<*n_d_p;nd++){ n_i*=(*d_p)[nd];}
  *ulli_p = (unsigned long long int *)wkspace_all0c(sizeof(unsigned long long int)*n_i);
  s=fread(*ulli_p,sizeof(unsigned long long int),n_i,fp); if (s!=n_i){ printf(" %% Warning, could not read from%s.\n",fname); exit(RET_READ_FAIL);}
  fclose(fp);fp=NULL;
}
