#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void binary_read(char *filename,int *bitj_p,int *nrows_p,int *ncols_p,unsigned char **A_p)
{
  /* reads a binary array stored in *filename; stores ouput in **A_p, allocated as necessary */
  int verbose=0;
  int alloc_flag = (*A_p==NULL);
  /*uint8_t*/ int bitj_x=0;
  /* uint32_t */ int nrows_x=0;
  /* uint32_t */ int ncols_x=0;
  int bitj=0,nrows=0,ncols=0,nrows_extend=0,wp_length=0;
  FILE *fp=NULL;
  if (verbose){ printf(" %% [entering binary_read] filename %s\n",filename); }
  if ((fp=fopen(filename,"r"))==NULL){ printf(" Warning! cannot open %s in binary_read\n",filename);}
  fread(&bitj_x,sizeof(/*uint8_t*/ int),1,fp); bitj = (int)bitj_x;
  if (verbose){ printf(" %% bitj %d\n",bitj);}
  fread(&nrows_x,sizeof(/* uint32_t */ int),1,fp); nrows = (int)nrows_x;
  if (verbose){ printf(" %% nrows %d\n",nrows);}
  fread(&ncols_x,sizeof(/* uint32_t */ int),1,fp); ncols = (int)ncols_x;
  if (verbose){ printf(" %% ncols %d\n",ncols);}
  if (bitj_p!=NULL){ *bitj_p = bitj;}
  if (nrows_p!=NULL){ *nrows_p = nrows;}
  if (ncols_p!=NULL){ *ncols_p = ncols;}
  nrows_extend = (bitj - (nrows % bitj)) % bitj;
  if (verbose){ printf(" %% nrows_extend %d\n",nrows_extend);}
  wp_length = (nrows+nrows_extend)/BIT8;
  if (verbose){ printf(" %% wp_length %d\n",wp_length);}
  if (A_p!=NULL){
    if (alloc_flag){ if (verbose){ printf(" %% allocating %d bytes\n",ncols*wp_length);} *A_p = wkspace_all0c(ncols*wp_length); if (!(*A_p)){ printf(" Warning! not enough memory in binary_read\n");}}
    if (verbose){ printf(" %% reading from disc\n");}
    fread(*A_p,sizeof(unsigned char),ncols*wp_length,fp);
    /* if (A_p!=NULL){ } */}
  if (fp!=NULL){ fclose(fp);fp=NULL;}
  if (verbose){ printf(" %% [finished binary_read] filename %s\n",filename); }
}

void binary_read_getsize(char *filename,int *bitj_p,int *nrows_p,int *ncols_p)
{
  /* reads a binary array stored in *filename, and returns the dimensions */
  /*uint8_t*/ int bitj_x=0;
  /* uint32_t */ int nrows_x=0;
  /* uint32_t */ int ncols_x=0;
  int bitj=0,nrows=0,ncols=0;
  FILE *fp=NULL;
  if ((fp=fopen(filename,"r"))==NULL){ printf(" Warning! cannot open %s in binary_read_getsize\n",filename);}
  fread(&bitj_x,sizeof(/*uint8_t*/ int),1,fp); bitj = (int)bitj_x;
  fread(&nrows_x,sizeof(/* uint32_t */ int),1,fp); nrows = (int)nrows_x;
  fread(&ncols_x,sizeof(/* uint32_t */ int),1,fp); ncols = (int)ncols_x;
  if (bitj_p!=NULL){ *bitj_p = bitj;}
  if (nrows_p!=NULL){ *nrows_p = nrows;}
  if (ncols_p!=NULL){ *ncols_p = ncols;}
  if (fp!=NULL){ fclose(fp);fp=NULL;}
}

void binary_write(char *filename,int bitj,int nrows,int ncols,unsigned char *A)
{
  /* Writes a binary array *A to filename, typically named ./filename.b16 ; 
     Assumes that &(A[nc*bsize(nrows)]) is the starting address of column nc. ;
     This means that for a typical M_handle M_An, one should call:
     binary_write(filename,bitj,M_An->ncols,M_An->rpop_b,M_An->wX);
  */
  int verbose=0;
  FILE *fpout=NULL;
  int nc=0,nrows_extend=0,wp_length;
  if (verbose){ printf(" %% [entering binary_write] bitj %d nrows %d ncols %d filename %s\n",bitj,nrows,ncols,filename);}
  if (verbose>1){ bprintf(A,bitj,nrows,ncols," %% A: ");}
  if ((fpout=fopen(filename,"w"))==NULL){ printf(" Warning! cannot open %s in binary_write\n",filename);}
  fwrite(&bitj,sizeof(/*uint8_t*/ int),1,fpout);fwrite(&nrows,sizeof(/* uint32_t */ int),1,fpout);fwrite(&ncols,sizeof(/* uint32_t */ int),1,fpout);
  if (verbose){ printf(" %% wrote bitj %d nrows %d ncols %d\n",bitj,nrows,ncols);}
  nrows_extend = (bitj - (nrows % bitj)) % bitj;
  if (verbose){ printf(" %% nrows_extend %d\n",nrows_extend);}
  wp_length = (nrows+nrows_extend)/BIT8;
  if (verbose){ printf(" %% wp_length %d\n",wp_length);}
  for (nc=0;nc<ncols;nc++){
    fwrite(&(A[nc*bsize(nrows)]),sizeof(unsigned char),wp_length,fpout);
    /* for (nc=0;nc<ncols;nc++){ } */}
  if (fpout!=NULL){ fclose(fpout);fpout=NULL;}
  if (verbose){ printf(" %% [finished binary_write] filename %s\n",filename);}
}

void uchar_write_as_binary(char *filename,int bitj,int nrows,int ncols,unsigned char *u_)
{
  /* Writes a binary array *u_ to filename, typically named ./filename.b16 ;
     We assume u_[nr + nc*nrows] is the (nr,nc) entry of u_ */
  int verbose=0;
  FILE *fpout=NULL;
  /*uint8_t*/ int bitj_x=0;
  /* uint32_t */ int nrows_x=0;
  /* uint32_t */ int ncols_x=0;
  int nrows_extend=0,ncols_extend=0;
  unsigned long wp_length=0;
  int nc=0,nr=0;
  unsigned char *wkspace_mark=NULL,*wp=NULL;
  unsigned char b1=0,b2=0;
  char bstr1[9],bstr2[9];
  if (verbose){ printf(" %% [entering uchar_write_as_binary] bitj %d nrows %d ncols %d filename %s\n",bitj,nrows,ncols,filename);}
  if (verbose){ raprintf(u_,"char",nrows,ncols," %% u_: ");}
  if ((fpout=fopen(filename,"w"))==NULL){ printf(" Warning! cannot open %s in uchar_write_as_binary\n",filename);}
  bitj_x = (/*uint8_t*/ int)bitj; nrows_x = (/* uint32_t */ int)nrows; ncols_x = (/* uint32_t */ int)ncols;
  fwrite(&bitj_x,sizeof(/*uint8_t*/ int),1,fpout);fwrite(&nrows_x,sizeof(/* uint32_t */ int),1,fpout);fwrite(&ncols_x,sizeof(/* uint32_t */ int),1,fpout);
  if (verbose){ printf(" %% wrote bitj %d nrows %d ncols %d\n",bitj,nrows,ncols);}
  nrows_extend = (bitj - (nrows % bitj)) % bitj;
  ncols_extend = (bitj - (ncols % bitj)) % bitj;
  wp_length = sizeof(long)*ceil((double)(nrows+nrows_extend)/(double)sizeof(long)); wp_length = wp_length/BIT8;
  if (verbose){ printf(" %% nrows %d+%d --> wp_length %d/BIT8=%d\n",nrows,nrows_extend,(int)wp_length*BIT8,(int)wp_length);}
  GLOBAL_wkspace_point = 0; wkspace_mark = wkspace_base; wp = wkspace_all0c(wp_length); if (!wp){ printf(" Warning! not enough memory in uchar_write_as_binary\n");}
  nc=0;
  while (nc<ncols){
    fill_uchar_zero(wp,wp_length);
    nr=0; 
    while (nr<nrows){ 
      if (verbose>1){ 
	printf("nr %d nc %d: ",nr,nc);
	b1 = (u_[nr + nc*nrows]>0) ? (1 << (7-(nr%BIT8))) : 0;getBinW(&b1,bstr1,8);
	printf("bitpattern %s set, ",bstr1);
	b1 = wp[nr/BIT8];getBinW(&b1,bstr1,8);
	b2 = (wp[nr/BIT8] | ((u_[nr + nc*nrows]>0) ? (1 << (7-(nr%BIT8))) : 0));getBinW(&b2,bstr2,8);
	printf("updating %s into %s.\n",bstr1,bstr2);
	/* if (verbose){ } */}
      wp[nr/BIT8] |= ((u_[nr + nc*nrows]>0) ? (1 << (7-(nr%BIT8))) : 0) ; nr++; /* while (nr<nrows){ } */}
    fwrite(wp,sizeof(unsigned char),(nrows + nrows_extend)/BIT8,fpout);
    nc ++;
    if (verbose){ printf("\r%d cols complete.", nc); fflush(stdout);}
    /* while (nc<ncols){ } */}
  if (verbose){ printf("\r finished writing %s\n",filename); fflush(stdout);}
  if (fpout!=NULL){ fclose(fpout);fpout=NULL;}
  if (verbose){ printf(" %% [finished uchar_write_as_binary] filename %s\n",filename);}
  wkspace_reset(wkspace_mark); GLOBAL_wkspace_point = 1;
}

void wrap_rand(char *filename,int bitj,int nrows,int ncols,int type_set)
{
  /* Here we generate a random/arbitrary unsigned-char array *u_ ;
     We assume that u_[nr + nc*nrows] is the (nr,nc) entry of u_ ;
     if type_set==0, then we don't actually fill this array at all. ;
     otherwise, if type_set!=0, we fill the array randomly ; */
  int verbose=0;
  FILE *fpout=NULL;
  /*uint8_t*/ int bitj_x=0;
  /* uint32_t */ int nrows_x=0;
  /* uint32_t */ int ncols_x=0;
  int nrows_extend=0,ncols_extend=0;
  unsigned long wp_length=0;
  int nc=0,nr=0;
  unsigned char *wkspace_mark=NULL,*wp=NULL;
  if (verbose){ printf(" %% [entering wrap_rand] bitj %d nrows %d ncols %d filename %s\n",bitj,nrows,ncols,filename);}
  if ((fpout=fopen(filename,"w"))==NULL){ printf(" Warning! cannot open %s in wrap_rand\n",filename);}
  bitj_x = (/*uint8_t*/ int)bitj; nrows_x = (/* uint32_t */ int)nrows; ncols_x = (/* uint32_t */ int)ncols;
  fwrite(&bitj_x,sizeof(/*uint8_t*/ int),1,fpout);fwrite(&nrows_x,sizeof(/* uint32_t */ int),1,fpout);fwrite(&ncols_x,sizeof(/* uint32_t */ int),1,fpout);
  if (verbose){ printf(" %% wrote bitj %d nrows %d ncols %d\n",bitj,nrows,ncols);}
  nrows_extend = (bitj - (nrows % bitj)) % bitj;
  ncols_extend = (bitj - (ncols % bitj)) % bitj;
  wp_length = sizeof(long)*ceil((double)(nrows+nrows_extend)/(double)sizeof(long)); wp_length = wp_length/BIT8;
  if (verbose){ printf(" %% nrows %d+%d --> wp_length %d/BIT8=%d\n",nrows,nrows_extend,(int)wp_length*BIT8,(int)wp_length);}
  GLOBAL_wkspace_point = 0; wkspace_mark = wkspace_base; wp = wkspace_all0c(wp_length); if (!wp){ printf(" Warning! not enough memory in wrap_rand\n");}
  nc=0;
  while (nc<ncols){
    fill_uchar_zero(wp,wp_length);
    nr=0; 
    if (type_set==0){/* do nothing */}
    else{ while (nr<nrows){ wp[nr/BIT8] = (unsigned char)rand(); nr+=BIT8; /* while (nr<nrows){ } */} /* else */}
    fwrite(wp,sizeof(unsigned char),(nrows + nrows_extend)/BIT8,fpout);
    nc ++;
    if (verbose){ printf("\r%d cols complete.", nc); fflush(stdout);}
    /* while (nc<ncols){ } */}
  if (verbose){ printf("\r finished writing %s\n",filename); fflush(stdout);}
  if (fpout!=NULL){ fclose(fpout);fpout=NULL;}
  if (verbose){ printf(" %% [finished wrap_rand] filename %s\n",filename);}
  wkspace_reset(wkspace_mark); GLOBAL_wkspace_point = 1;
}

void wrap_transpose(char *fname_0in,char *fname_out,int n_bytes_per_read)
{
  /* Slowly converts a binary array (stored in *fname_0in) to a transpose of that array (stored in *fname_out) ;
     This is predicated on : ;
     1. a header of size 3*sizeof(int) (containing (int) versions of bitj, n_row and n_col) ;
     2. grabbing n_bytes_per_read from rows of *fname_0in and converting them into values ranging across n_bytes_per_read separate columns of *fname_out ;
     3. storage of (n_col+n_col_extend)*n_bytes_per_read is allocated to wp (workspace). */
  int verbose=GLOBAL_verbose;
  FILE *fp_0in=NULL;
  /*uint8_t*/ int bitj_x=0;
  /* uint32_t */ int n_row_x=0;
  /* uint32_t */ int n_col_x=0;
  int bitj=0,n_row=0,n_col=0,n_row_extend=0,n_col_extend=0;
  unsigned long long int wp_length=0;
  int current_row=0,nf=0,nc=0,nr=0,nrx=0;
  FILE *fp_out=NULL;
  unsigned char *wkspace_mark=NULL;
  unsigned char *wp=NULL;
  unsigned char *ip=NULL;
  unsigned char b1=0,b2=0;
  unsigned long long int ulli_tab=0;
  char bstr1[9],bstr2[9],kstr[FNAMESIZE];
  double ct=0,rt=0,tratio=0,tau_est=0;
  if (verbose){ printf(" %% [entering wrap_transpose] fname_0in %s --> %s n_bytes_per_read %d\n",fname_0in,fname_out,n_bytes_per_read);}
  if (access(fname_0in,F_OK)!=0){ printf(" Warning! cannot access %s in wrap_transpose\n",fname_0in);}
  if ((fp_0in=fopen(fname_0in,"r"))==NULL){ printf(" Warning! cannot open %s in wrap_transpose\n",fname_0in);}
  fread(&bitj_x,sizeof(/*uint8_t*/ int),1,fp_0in);fread(&n_row_x,sizeof(/* uint32_t */ int),1,fp_0in);fread(&n_col_x,sizeof(/* uint32_t */ int),1,fp_0in); bitj=(int)bitj_x; n_row=(int)n_row_x; n_col=(int)n_col_x;
  if (fp_0in!=NULL){ fclose(fp_0in);fp_0in=NULL;}
  n_bytes_per_read = BIT8*(n_bytes_per_read/BIT8); if (n_bytes_per_read==0){ n_bytes_per_read=BIT8;}
  n_row_extend = (bitj - (n_row % bitj)) % bitj;
  n_col_extend = (bitj - (n_col % bitj)) % bitj;
  if (verbose){ printf(" %% read bitj %d n_row %d+%d n_col %d+%d\n",bitj,n_row,n_row_extend,n_col,n_col_extend);}
  if ((fp_0in=fopen(fname_0in,"r"))==NULL){ printf(" Warning! cannot open %s in wrap_transpose\n",fname_0in);}
  if ((fp_out=fopen(fname_out,"w"))==NULL){ printf(" Warning! cannot open %s in wrap_transpose\n",fname_out);}
  bitj_x = (/*uint8_t*/ int)bitj; n_row_x = (/* uint32_t */ int)n_row; n_col_x = (/* uint32_t */ int)n_col;
  fwrite(&bitj_x,sizeof(/*uint8_t*/ int),1,fp_out);fwrite(&n_col_x,sizeof(/* uint32_t */ int),1,fp_out);fwrite(&n_row_x,sizeof(/* uint32_t */ int),1,fp_out);
  wp_length = (unsigned long long int)(n_col+n_col_extend) * (unsigned long long int)n_bytes_per_read;
  wp_length = rup(wp_length,sizeof(long)); wp_length /= BIT8;
  if (verbose){ printf(" %% wp_length %lldUB --> %0.2fKB --> %0.2fMB --> %0.2fGB\n",wp_length,(double)wp_length/(double)pow(1024,1),(double)wp_length/(double)pow(1024,2),(double)wp_length/(double)pow(1024,3));}
  GLOBAL_wkspace_point = 0; wkspace_mark = wkspace_base; 
  wp = wkspace_all0c(wp_length); if (!wp){ printf(" Warning! not enough memory, reduce n_bytes_per_read %d in wrap_transpose\n",n_bytes_per_read);}
  fill_uchar_zero(wp,wp_length);
  ip = wkspace_all0c(n_bytes_per_read / BIT8 + 1); if (!ip){ printf(" Warning! not enough memory, reduce n_bytes_per_read %d in wrap_transpose\n",n_bytes_per_read);}
  GLOBAL_tic(3);
  current_row=0;
  while (current_row<n_row){
    GLOBAL_toc(3,0,""); ct = GLOBAL_elct[3]; rt = GLOBAL_elrt[3]; tratio=ct/maximum(1,rt); tau_est = rt*(double)n_row/(double)maximum(1,current_row);
    if ((verbose>0) && (current_row%1024==0)){ printf(" %% current_row %.5d/%.5d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh)\n",current_row,n_row,ct,ct/3600,rt,rt/3600,tratio,tau_est,tau_est/3600);}
    fseek(fp_0in,(size_t)3*sizeof(int) + (size_t)(current_row/BIT8),SEEK_SET);
    fill_uchar_zero(wp,wp_length);
    nr = minimum(n_bytes_per_read,n_row /* + n_row_extend */ - current_row);
    nrx = rup(nr,BIT8);
    nc=0;
    while (nc<n_col){ 
      fread(ip,sizeof(unsigned char),nrx/BIT8,fp_0in);
      if (verbose>2){
	getBinW(ip,kstr,nrx);
	printf("ip %s\n",kstr);
	/* if (verbose){ } */}
      nf=0;
      while (nf<nr){
	ulli_tab = (unsigned long long int)(nc/8) + (unsigned long long int)nf*(unsigned long long int)((n_col+n_col_extend)/BIT8);
	if (verbose>2){
	  printf(" nc %.6d current_row %.6d nf %.6d (nr %.6d), ",nc,current_row,nf,nr);
	  b1 = ip[nf/BIT8];getBinW(&b1,bstr1,8);
	  b2 = ((((ip[nf/BIT8] >> (7-(nf%BIT8))) & 1 ) << (7 - (nc%BIT8))));getBinW(&b2,bstr2,8);
	  printf("reading %s --> %s, ",bstr1,bstr2);
	  b1 = wp[ulli_tab];getBinW(&b1,bstr1,8);
	  b2 = (wp[ulli_tab] | ((((ip[nf/BIT8] >> (7-(nf%BIT8))) & 1 ) << (7 - (nc%BIT8)))));getBinW(&b2,bstr2,8);
	  printf("updating %s --> %s.\n",bstr1,bstr2);
	  /* if (verbose){ } */}
	wp[ulli_tab] |= ((((ip[nf/BIT8] >> (7-(nf%BIT8))) & 1 ) << (7 - (nc%BIT8))));
	nf++;
	/* while (nf<n_bytes_per_read){ } */}
      nc++;
      fseek(fp_0in,(n_row+n_row_extend)/BIT8 - nrx/BIT8,SEEK_CUR);
      /* while (nc<n_col){ } */}
    while (nc<n_col+n_col_extend){
      nf=0;
      while (nf<nr){
	ulli_tab = (unsigned long long int)(nc/8) + (unsigned long long int)nf*(unsigned long long int)((n_col+n_col_extend)/BIT8);
    	wp[ulli_tab] &= (~0 ^ (1 << (7-(nc%BIT8))));
    	nf++;
    	/* while (nf<n_bytes_per_read){ } */}
      nc++;
      /* while (nc<n_col){ } */}
    fwrite(wp,sizeof(unsigned char),(unsigned long long int)((n_col+n_col_extend)/BIT8) * (unsigned long long int)nr,fp_out);
    current_row += n_bytes_per_read;
    if (verbose>2){ printf("\r%d rows complete.", current_row); fflush(stdout);}
    /* while (current_row<n_row){ } */}
  if (verbose>2){ printf("\r finished writing %s\n",fname_out); fflush(stdout);}
  if (fp_out!=NULL){ fclose(fp_out);fp_out=NULL;}
  if (fp_0in!=NULL){ fclose(fp_0in);fp_0in=NULL;}
  if (verbose){ printf(" %% [finished wrap_transpose] fname_0in %s\n",fname_0in);}
  wkspace_reset(wkspace_mark); GLOBAL_wkspace_point = 1;
}

void wrap_transpose_test()
{
  int verbose=1;
  char fname_0in[1024];
  char fname_1in[1024];
  char fname_out[1024];
  unsigned char *u_=NULL,*v_=NULL,*w_=NULL,*x_=NULL;
  int x_n_row=0,x_n_col=0;
  int u_n_row=0,u_n_col=0;
  int v_n_row=0,v_n_col=0;
  int nrow=0,ncol=0;
  int bitj=16;
  int e=0,nl=0;
  if (verbose){ printf(" %% [entering wrap_transpose_test]\n");}
  u_n_row = 137; u_n_col = 42;
  u_ = (unsigned char *)wkspace_all0c(sizeof(unsigned char)*u_n_row*u_n_col);
  for (ncol=0;ncol<u_n_col;ncol++){
    for (nrow=0;nrow<u_n_row;nrow++){
      u_[nrow+ncol*u_n_row] = (rand01>0.5 ? 1 : 0);
      /* for (nrow=0;nrow<u_n_row;nrow++){ } */}
    /* for (ncol=0;ncol<u_n_col;ncol++){ } */}
  sprintf(fname_0in,"./dir_test/wrap_transpose_test_0in.b16");
  sprintf(fname_1in,"./dir_test/wrap_transpose_test_1in.b16");
  sprintf(fname_out,"./dir_test/wrap_transpose_test_out.b16");
  uchar_write_as_binary(fname_0in,bitj,u_n_row,u_n_col,u_);
  wrap_transpose(fname_0in,fname_out,1*BIT8);
  wrap_transpose(fname_out,fname_1in,1*BIT8);
  binary_read(fname_0in,&bitj,&u_n_row,&u_n_col,&w_);
  binary_read(fname_1in,&bitj,&x_n_row,&x_n_col,&x_);
  binary_read(fname_out,&bitj,&v_n_row,&v_n_col,&v_);
  for (nl=0;nl<(rup(u_n_row,BIT8)/BIT8)*u_n_col;nl++){ e+=abs((int)w_[nl]-(int)x_[nl]);}
  if (verbose>0){ printf(" %% error %d\n",e);}
  if (verbose>1){
    printf("%% w_: \n");
    bprintf(w_,bitj,u_n_col,u_n_row," %% w_: ");
    printf("%% \n");
    printf("%% x_: \n");
    bprintf(x_,bitj,x_n_col,x_n_row," %% x_: ");
    //printf("%% v_: \n");
    //bprintf(v_,bitj,v_n_col,v_n_row," %% v_: ");
    /* if (verbose>1){ } */}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_transpose_test]\n");}  
}

void wrap_transpose_driver()
{
  wrap_transpose(GLOBAL_fname_b16_0in,GLOBAL_fname_b16_out,GLOBAL_n_bytes_per_read);
}
