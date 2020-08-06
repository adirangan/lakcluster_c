#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void raprintf(void *v_,char *type,int rows,int cols,char *prefix)
{
  /* prints out ar_ys of varying types */
  int nr=0,nc=0,tmp=0;
  float *f_=NULL;
  double *d_=NULL; int d_z_flag=0;
  int *i_=NULL;
  unsigned int *ui_=NULL;
  long *l_=NULL;
  long long *ll_=NULL;
  unsigned long int *ul_=NULL;
  unsigned long long int *ull_=NULL;
  char *c_=NULL;
  unsigned char *uc_=NULL;
  /* fftw_complex *ff_=NULL; */
  /* double dtmp=0; */
  double printftol=0.000000001;
  if (strcmp(type,"float")==0){
    f_ = (float *) v_;
    for (nr=0;nr<rows;nr++){ 
      printf("%s",prefix); 
      for (nc=0;nc<cols;nc++){ 
	if (fabs(f_[nr+nc*rows]-(int)f_[nr+nc*rows])<printftol){ printf(" %s%d",(int)f_[nr+nc*rows]>0 ? "+" : ((int)f_[nr+nc*rows]<0 ? "" : " "),(int)f_[nr+nc*rows]);}
	else{ printf(" %s%f",f_[nr+nc*rows]>0 ? "+" : (f_[nr+nc*rows]<0 ? "" : " "),f_[nr+nc*rows]);}}
      printf("\n");}}
  else if (strcmp(type,"double")==0){
    if (rows>1 && cols==1){ tmp=rows; rows=cols; cols=tmp;}
    d_ = (double *) v_;
    if (0){
      d_z_flag=1; for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ if (fabs(d_[nr+nc*rows]-(int)d_[nr+nc*rows])>printftol){ d_z_flag=0;}}}
      for (nr=0;nr<rows;nr++){ 
	printf("%s",prefix); 
	for (nc=0;nc<cols;nc++){ 
	  if (d_z_flag && fabs(d_[nr+nc*rows]-(int)d_[nr+nc*rows])<printftol){ printf(" %s%d",(int)d_[nr+nc*rows]>0 ? "+" : ((int)d_[nr+nc*rows]<0 ? "" : " "),(int)d_[nr+nc*rows]);}
	  else{ /* printf(" %s%.3f",d_[nr+nc*rows]>0 ? "+" : (d_[nr+nc*rows]<0 ? "" : " "),d_[nr+nc*rows]); */ 
	    if (fabs(d_[nr+nc*rows])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*rows]);}
	    /* if !d_z_flag */}
	  /* for (nc=0;nc<cols;nc++){ } */}
	printf("\n");} /* if 0 */}
    if (1){
      for (nr=0;nr<rows;nr++){ printf("%s",prefix); 
	for (nc=0;nc<cols;nc++){ if (fabs(d_[nr+nc*rows])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*rows]);} /* for (nc=0;nc<cols;nc++){ } */}
	printf("\n"); /* for (nr=0;nr<rows;nr++){ } */} /* if 1 */}
    /* if double */}
  else if (strcmp(type,"double_trn")==0 && rows>1){
    d_ = (double *) v_;
    for (nc=0;nc<cols;nc++){ printf("%s",prefix); 
      for (nr=0;nr<rows;nr++){ if (fabs(d_[nr+nc*rows])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",d_[nr+nc*rows]);} /* for (nr=0;nr<rows;nr++){ } */}
      printf("\n"); /* for (nc=0;nc<cols;nc++){ } */}
    /* if double_trn */}
  else if (strcmp(type,"long")==0){
    l_ = (long *) v_;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %ld",l_[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"long long int")==0){
    ll_ = (long long *) v_;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %lld",ll_[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"int")==0){
    i_ = (int *) v_;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %d",i_[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"unsigned int")==0){
    ui_ = (unsigned int *) v_;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %d",(int)ui_[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"unsigned long int")==0){
    ul_ = (unsigned long int *) v_;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %ld",ul_[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"unsigned long long int")==0){
    ull_ = (unsigned long long int *) v_;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %lld",ull_[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"alpha")==0){
    c_ = (char *) v_;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %c",(int)c_[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"char")==0){
    c_ = (char *) v_;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %d",(int)c_[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"unsigned char")==0){
    uc_ = (unsigned char *) v_;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %d",(int)uc_[nr+nc*rows]);} printf("\n");}}
  else{ printf(" warning, poor type %s in raprintf\n",type);}
}

void ra_fprintf(char *fname,void *v_,char *type,int rows,int cols,char *prefix)
{
  /* prints out arrays of varying types */
  FILE *fp=NULL;
  int nr=0,nc=0,tmp=0;
  float *f_=NULL;
  double *d_=NULL; int d_z_flag=1;
  int *i_=NULL;
  unsigned int *ui_=NULL;
  long *l_=NULL;
  long long *ll_=NULL;
  char *c_=NULL;
  unsigned char *uc_=NULL;
  /* fftw_complex *ff_=NULL; */
  /* double dtmp=0; */
  double printftol=0.000000001;
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! cannot open %s in ra_fprintf\n",fname); fp=stdout;}
  if (strcmp(type,"float")==0){
    f_ = (float *) v_;
    for (nr=0;nr<rows;nr++){ 
      fprintf(fp,"%s",prefix); 
      for (nc=0;nc<cols;nc++){ 
	if (fabs(f_[nr+nc*rows]-(int)f_[nr+nc*rows])<printftol){ fprintf(fp," %s%d",(int)f_[nr+nc*rows]>0 ? "+" : ((int)f_[nr+nc*rows]<0 ? "" : " "),(int)f_[nr+nc*rows]);}
	else{ fprintf(fp," %s%f",f_[nr+nc*rows]>0 ? "+" : (f_[nr+nc*rows]<0 ? "" : " "),f_[nr+nc*rows]);}}
      fprintf(fp,"\n");}}
  else if (strcmp(type,"double")==0){
    if (rows>1 && cols==1){ tmp=rows; rows=cols; cols=tmp;}
    d_ = (double *) v_;
    if (0){
      d_z_flag=1; for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ if (fabs(d_[nr+nc*rows]-(int)d_[nr+nc*rows])>printftol){ d_z_flag=0;}}}
      for (nr=0;nr<rows;nr++){ 
	fprintf(fp,"%s",prefix); 
	for (nc=0;nc<cols;nc++){ 
	  if (d_z_flag && fabs(d_[nr+nc*rows]-(int)d_[nr+nc*rows])<printftol){ fprintf(fp," %s%d",(int)d_[nr+nc*rows]>0 ? "+" : ((int)d_[nr+nc*rows]<0 ? "" : " "),(int)d_[nr+nc*rows]);}
	  else{ /* fprintf(fp," %s%.3f",d_[nr+nc*rows]>0 ? "+" : (d_[nr+nc*rows]<0 ? "" : " "),d_[nr+nc*rows]); */ 
	    if (fabs(d_[nr+nc*rows])<printftol){ fprintf(fp,"   0.0  ");} else{ fprintf(fp," %0.16f",d_[nr+nc*rows]);}
	    /* if !d_z_flag */}
	  /* for (nc=0;nc<cols;nc++){ } */}
	fprintf(fp,"\n");} /* if 0 */}
    if (1){
      for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); 
	for (nc=0;nc<cols;nc++){ if (fabs(d_[nr+nc*rows])<printftol){ fprintf(fp,"   0.0  ");} else{ fprintf(fp," %+07.3f",d_[nr+nc*rows]);} /* for (nc=0;nc<cols;nc++){ } */}
	fprintf(fp,"\n"); /* for (nr=0;nr<rows;nr++){ } */} /* if 1 */}
    /* if double */}
  else if (strcmp(type,"double_trn")==0 && rows>1){
    d_ = (double *) v_;
    for (nc=0;nc<cols;nc++){ fprintf(fp,"%s",prefix); 
      for (nr=0;nr<rows;nr++){ 
	if (d_z_flag && fabs(d_[nr+nc*rows]-(int)d_[nr+nc*rows])<printftol){ fprintf(fp," %s%d",(int)d_[nr+nc*rows]>0 ? "+" : ((int)d_[nr+nc*rows]<0 ? "" : " "),(int)d_[nr+nc*rows]);}
	else{ if (fabs(d_[nr+nc*rows])<printftol){ fprintf(fp,"   0.0  ");} else{ fprintf(fp," %0.16f",d_[nr+nc*rows]);} /* if !d_z_flag */}
	/* for (nr=0;nr<rows;nr++){ } */}
      fprintf(fp,"\n"); /* for (nc=0;nc<cols;nc++){ } */}
    /* if double_trn */}
  else if (strcmp(type,"long")==0 || strcmp(type,"long int")==0){
    l_ = (long *) v_;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %ld",l_[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"long long int")==0){
    ll_ = (long long *) v_;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %lld",ll_[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"int")==0){
    i_ = (int *) v_;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %d",i_[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"unsigned int")==0){
    ui_ = (unsigned int *) v_;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %d",(int)ui_[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"alpha")==0){
    c_ = (char *) v_;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %c",(int)c_[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"char")==0){
    c_ = (char *) v_;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %d",(int)c_[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"unsigned char")==0){
    uc_ = (unsigned char *) v_;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %d",(int)uc_[nr+nc*rows]);} fprintf(fp,"\n");}}
  else{ fprintf(fp," warning, poor type %s in raprintf\n",type);}
  if (fp!=stdout){ fclose(fp); fp=NULL;}
}

void getBinW(unsigned char *w_, char *str_, int k)
{
  /* extracts the first k bits of unsigned char *w_ and stores as output *str_ of zeros and ones */
  /* test with: b1=0;do{ getBinW(&b1,bstr_1,8); printf(" %% (%u --> %s)",(int)b1,bstr_1);} while (++b1!=0); printf("\n"); */
  unsigned char wc = 0;
  int i=0,j=0,l=0;
  *(str_+k) = '\0';
  while (i<k){
    j=i/BIT8;l=7-(i%BIT8); wc=w_[j];
    str_[i] = ((wc >> l) & 1) + '0';
    i++; /* while (i<k){ } */}
}

void bprintf(unsigned char *w_,int bitj,int nrows,int ncols,char *prefix)
{
  /* prints out binary array stored in *w_, assuming w_ is stored in column-major order (i.e., columns varying quickly, rows varying slowly)
     This is useful for printing row-vectors on a single line. */
  char kstr_[ncols+2];
  int ncols_extend = (bitj - (ncols % bitj)) % bitj,nr=0;
  for (nr=0;nr<nrows;nr++){ getBinW(&(w_[nr*((ncols+ncols_extend)/BIT8)]),kstr_,ncols); printf("%s%s -- %d,%d\n",prefix,kstr_,nr,ncols);}
}
