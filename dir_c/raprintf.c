void raprintf(void *vra,char *type,int rows,int cols,char *prefix)
{
  /* prints out arrays of varying types */
  int nr=0,nc=0,tmp=0;
  float *fra=NULL;
  double *dra=NULL; int draz_flag=0;
  int *ira=NULL;
  unsigned int *uira=NULL;
  long *lra=NULL;
  long long *llra=NULL;
  char *cra=NULL;
  unsigned char *ucra=NULL;
  /* fftw_complex *ffra=NULL; */
  /* double dtmp=0; */
  double printftol=0.000000001;
  if (strcmp(type,"float")==0){
    fra = (float *) vra;
    for (nr=0;nr<rows;nr++){ 
      printf("%s",prefix); 
      for (nc=0;nc<cols;nc++){ 
	if (fabs(fra[nr+nc*rows]-(int)fra[nr+nc*rows])<printftol){ printf(" %s%d",(int)fra[nr+nc*rows]>0 ? "+" : ((int)fra[nr+nc*rows]<0 ? "" : " "),(int)fra[nr+nc*rows]);}
	else{ printf(" %s%f",fra[nr+nc*rows]>0 ? "+" : (fra[nr+nc*rows]<0 ? "" : " "),fra[nr+nc*rows]);}}
      printf("\n");}}
  else if (strcmp(type,"double")==0){
    if (rows>1 && cols==1){ tmp=rows; rows=cols; cols=tmp;}
    dra = (double *) vra;
    if (0){
      draz_flag=1; for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ if (fabs(dra[nr+nc*rows]-(int)dra[nr+nc*rows])>printftol){ draz_flag=0;}}}
      for (nr=0;nr<rows;nr++){ 
	printf("%s",prefix); 
	for (nc=0;nc<cols;nc++){ 
	  if (draz_flag && fabs(dra[nr+nc*rows]-(int)dra[nr+nc*rows])<printftol){ printf(" %s%d",(int)dra[nr+nc*rows]>0 ? "+" : ((int)dra[nr+nc*rows]<0 ? "" : " "),(int)dra[nr+nc*rows]);}
	  else{ /* printf(" %s%.3f",dra[nr+nc*rows]>0 ? "+" : (dra[nr+nc*rows]<0 ? "" : " "),dra[nr+nc*rows]); */ 
	    if (fabs(dra[nr+nc*rows])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",dra[nr+nc*rows]);}
	    /* if !draz_flag */}
	  /* for (nc=0;nc<cols;nc++){ } */}
	printf("\n");} /* if 0 */}
    if (1){
      for (nr=0;nr<rows;nr++){ printf("%s",prefix); 
	for (nc=0;nc<cols;nc++){ if (fabs(dra[nr+nc*rows])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",dra[nr+nc*rows]);} /* for (nc=0;nc<cols;nc++){ } */}
	printf("\n"); /* for (nr=0;nr<rows;nr++){ } */} /* if 1 */}
    /* if double */}
  else if (strcmp(type,"double_trn")==0 && rows>1){
    dra = (double *) vra;
    for (nc=0;nc<cols;nc++){ printf("%s",prefix); 
      for (nr=0;nr<rows;nr++){ if (fabs(dra[nr+nc*rows])<printftol){ printf("    .   ");} else{ printf(" %+07.3f",dra[nr+nc*rows]);} /* for (nr=0;nr<rows;nr++){ } */}
      printf("\n"); /* for (nc=0;nc<cols;nc++){ } */}
    /* if double_trn */}
  else if (strcmp(type,"long")==0){
    lra = (long *) vra;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %ld",lra[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"long long int")==0){
    llra = (long long *) vra;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %lld",llra[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"int")==0){
    ira = (int *) vra;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %d",ira[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"unsigned int")==0){
    uira = (unsigned int *) vra;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %d",(int)uira[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"alpha")==0){
    cra = (char *) vra;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %c",(int)cra[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"char")==0){
    cra = (char *) vra;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %d",(int)cra[nr+nc*rows]);} printf("\n");}}
  else if (strcmp(type,"unsigned char")==0){
    ucra = (unsigned char *) vra;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %d",(int)ucra[nr+nc*rows]);} printf("\n");}}
  else{ printf(" warning, poor type %s in raprintf\n",type);}
}

void ra_fprintf(char *fname,void *vra,char *type,int rows,int cols,char *prefix)
{
  /* prints out arrays of varying types */
  FILE *fp=NULL;
  int nr=0,nc=0,tmp=0;
  float *fra=NULL;
  double *dra=NULL; int draz_flag=1;
  int *ira=NULL;
  unsigned int *uira=NULL;
  long *lra=NULL;
  long long *llra=NULL;
  char *cra=NULL;
  unsigned char *ucra=NULL;
  /* fftw_complex *ffra=NULL; */
  /* double dtmp=0; */
  double printftol=0.000000001;
  if ((fp=fopen(fname,"w"))==NULL){ printf(" %% Warning! cannot open %s in ra_fprintf\n",fname); fp=stdout;}
  if (strcmp(type,"float")==0){
    fra = (float *) vra;
    for (nr=0;nr<rows;nr++){ 
      fprintf(fp,"%s",prefix); 
      for (nc=0;nc<cols;nc++){ 
	if (fabs(fra[nr+nc*rows]-(int)fra[nr+nc*rows])<printftol){ fprintf(fp," %s%d",(int)fra[nr+nc*rows]>0 ? "+" : ((int)fra[nr+nc*rows]<0 ? "" : " "),(int)fra[nr+nc*rows]);}
	else{ fprintf(fp," %s%f",fra[nr+nc*rows]>0 ? "+" : (fra[nr+nc*rows]<0 ? "" : " "),fra[nr+nc*rows]);}}
      fprintf(fp,"\n");}}
  else if (strcmp(type,"double")==0){
    if (rows>1 && cols==1){ tmp=rows; rows=cols; cols=tmp;}
    dra = (double *) vra;
    if (0){
      draz_flag=1; for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ if (fabs(dra[nr+nc*rows]-(int)dra[nr+nc*rows])>printftol){ draz_flag=0;}}}
      for (nr=0;nr<rows;nr++){ 
	fprintf(fp,"%s",prefix); 
	for (nc=0;nc<cols;nc++){ 
	  if (draz_flag && fabs(dra[nr+nc*rows]-(int)dra[nr+nc*rows])<printftol){ fprintf(fp," %s%d",(int)dra[nr+nc*rows]>0 ? "+" : ((int)dra[nr+nc*rows]<0 ? "" : " "),(int)dra[nr+nc*rows]);}
	  else{ /* fprintf(fp," %s%.3f",dra[nr+nc*rows]>0 ? "+" : (dra[nr+nc*rows]<0 ? "" : " "),dra[nr+nc*rows]); */ 
	    if (fabs(dra[nr+nc*rows])<printftol){ fprintf(fp,"   0.0  ");} else{ fprintf(fp," %0.16f",dra[nr+nc*rows]);}
	    /* if !draz_flag */}
	  /* for (nc=0;nc<cols;nc++){ } */}
	fprintf(fp,"\n");} /* if 0 */}
    if (1){
      for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); 
	for (nc=0;nc<cols;nc++){ if (fabs(dra[nr+nc*rows])<printftol){ fprintf(fp,"   0.0  ");} else{ fprintf(fp," %+07.3f",dra[nr+nc*rows]);} /* for (nc=0;nc<cols;nc++){ } */}
	fprintf(fp,"\n"); /* for (nr=0;nr<rows;nr++){ } */} /* if 1 */}
    /* if double */}
  else if (strcmp(type,"double_trn")==0 && rows>1){
    dra = (double *) vra;
    for (nc=0;nc<cols;nc++){ fprintf(fp,"%s",prefix); 
      for (nr=0;nr<rows;nr++){ 
	if (draz_flag && fabs(dra[nr+nc*rows]-(int)dra[nr+nc*rows])<printftol){ fprintf(fp," %s%d",(int)dra[nr+nc*rows]>0 ? "+" : ((int)dra[nr+nc*rows]<0 ? "" : " "),(int)dra[nr+nc*rows]);}
	else{ if (fabs(dra[nr+nc*rows])<printftol){ fprintf(fp,"   0.0  ");} else{ fprintf(fp," %0.16f",dra[nr+nc*rows]);} /* if !draz_flag */}
	/* for (nr=0;nr<rows;nr++){ } */}
      fprintf(fp,"\n"); /* for (nc=0;nc<cols;nc++){ } */}
    /* if double_trn */}
  else if (strcmp(type,"long")==0 || strcmp(type,"long int")==0){
    lra = (long *) vra;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %ld",lra[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"long long int")==0){
    llra = (long long *) vra;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %lld",llra[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"int")==0){
    ira = (int *) vra;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %d",ira[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"unsigned int")==0){
    uira = (unsigned int *) vra;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %d",(int)uira[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"alpha")==0){
    cra = (char *) vra;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %c",(int)cra[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"char")==0){
    cra = (char *) vra;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %d",(int)cra[nr+nc*rows]);} fprintf(fp,"\n");}}
  else if (strcmp(type,"unsigned char")==0){
    ucra = (unsigned char *) vra;
    for (nr=0;nr<rows;nr++){ fprintf(fp,"%s",prefix); for (nc=0;nc<cols;nc++){ fprintf(fp," %d",(int)ucra[nr+nc*rows]);} fprintf(fp,"\n");}}
  else{ fprintf(fp," warning, poor type %s in raprintf\n",type);}
  if (fp!=stdout){ fclose(fp); fp=NULL;}
}

void getBinW(unsigned char *w, char *str, int k)
{
  /* extracts the first k bits of unsigned char *w and stores as output *str of zeros and ones */
  /* test with: b1=0;do{ getBinW(&b1,bstr1,5); printf(" %% (%u --> %s)",(int)b1,bstr1);} while (++b1!=0); printf("\n"); */
  unsigned char wc = 0;
  int ii=0,jj=0,kk=0;
  *(str+k) = '\0';
  while (ii<k){
    jj=ii/BIT8;kk=7-(ii%BIT8); wc=w[jj];
    str[ii] = ((wc >> kk) & 1) + '0';
    ii++; /* while (ii<k){ } */}
}

void bprintf(unsigned char *w,int bitj,int nrows,int ncols,char *prefix)
{
  /* prints out binary array stored in *w, assuming w is stored in column-major order (first column compressed first). 
     This is useful for printing row-vectors on a single line. */
  char kstr[ncols+2];
  int ncols_extend = (bitj - (ncols % bitj)) % bitj,nr=0;
  for (nr=0;nr<nrows;nr++){ getBinW(&(w[nr*((ncols+ncols_extend)/BIT8)]),kstr,ncols); printf("%s%s -- %d,%d\n",prefix,kstr,nr,ncols);}
}
