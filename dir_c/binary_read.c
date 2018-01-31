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

void binary_write(int bitj,unsigned char *bra,int nrows,int ncols,char *filename)
{
  /* Writes a binary array *bra to filename, typically named ./filename.b16 ;
     We assume bra[nr + nc*nrows] is the (nr,nc) entry of bra */
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
  if (verbose){ printf(" %% [entering binary_write] bitj %d nrows %d ncols %d filename %s\n",bitj,nrows,ncols,filename);}
  if (verbose){ raprintf(bra,"char",nrows,ncols," %% bra: ");}
  if ((fpout=fopen(filename,"w"))==NULL){ printf(" Warning! cannot open %s in binary_write\n",filename);}
  bitj_x = (/*uint8_t*/ int)bitj; nrows_x = (/* uint32_t */ int)nrows; ncols_x = (/* uint32_t */ int)ncols;
  fwrite(&bitj_x,sizeof(/*uint8_t*/ int),1,fpout);fwrite(&nrows_x,sizeof(/* uint32_t */ int),1,fpout);fwrite(&ncols_x,sizeof(/* uint32_t */ int),1,fpout);
  if (verbose){ printf(" %% wrote bitj %d nrows %d ncols %d\n",bitj,nrows,ncols);}
  nrows_extend = (bitj - (nrows % bitj)) % bitj;
  ncols_extend = (bitj - (ncols % bitj)) % bitj;
  wp_length = sizeof(long)*ceil((double)(nrows+nrows_extend)/(double)sizeof(long)); wp_length = wp_length/BIT8;
  if (verbose){ printf(" %% nrows %d+%d --> wp_length %d/BIT8=%d\n",nrows,nrows_extend,(int)wp_length*BIT8,(int)wp_length);}
  GLOBAL_wkspace_point = 0; wkspace_mark = wkspace_base; wp = wkspace_all0c(wp_length); if (!wp){ printf(" Warning! not enough memory in binary_write\n");}
  nc=0;
  while (nc<ncols){
    fill_uchar_zero(wp,wp_length);
    nr=0; 
    while (nr<nrows){ 
      if (verbose>1){ 
	printf("nr %d nc %d: ",nr,nc);
	b1 = (bra[nr + nc*nrows]>0) ? (1 << (7-(nr%BIT8))) : 0;getBinW(&b1,bstr1,8);
	printf("bitpattern %s set, ",bstr1);
	b1 = wp[nr/BIT8];getBinW(&b1,bstr1,8);
	b2 = (wp[nr/BIT8] | ((bra[nr + nc*nrows]>0) ? (1 << (7-(nr%BIT8))) : 0));getBinW(&b2,bstr2,8);
	printf("updating %s into %s.\n",bstr1,bstr2);
	/* if (verbose){ } */}
      wp[nr/BIT8] |= ((bra[nr + nc*nrows]>0) ? (1 << (7-(nr%BIT8))) : 0) ; nr++; /* while (nr<nrows){ } */}
    fwrite(wp,sizeof(unsigned char),(nrows + nrows_extend)/BIT8,fpout);
    nc ++;
    if (verbose){ printf("\r%d cols complete.", nc); fflush(stdout);}
    /* while (nc<ncols){ } */}
  if (verbose){ printf("\r finished writing %s\n",filename); fflush(stdout);}
  if (fpout!=NULL){ fclose(fpout);fpout=NULL;}
  if (verbose){ printf(" %% [finished binary_write] filename %s\n",filename);}
  wkspace_reset(wkspace_mark); GLOBAL_wkspace_point = 1;
}

void wrap_rand(int bitj,int nrows,int ncols,int type_set,char *filename)
{
  /* Here we generate a random/arbitrary binary array *bra ;
     We assume that bra[nr + nc*nrows] is the (nr,nc) entry of bra ;
     if type_set==0, then we don't actually fill this array at all. ;
     otherwise, if type_set!-0, we fill the array randomly ; */
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

void wrap_transpose(char *filename,char *outname,int nbytes_per_read)
{
  /* Slowly converts a binary array (stored in *filename) to a transpose of that array (stored in *outname) ;
     This is predicated on : ;
     1. a header of size 3*sizeof(int) (containing (int) versions of bitj, nrows and ncols) ;
     2. grabbing nbytes_per_read from rows of *filename and converting them into values ranging across nbytes_per_read separate columns of *outname ;
     3. storage of (ncols+ncols_extend)*nbyptes_per_read is allocated to wp (workspace). */
  int verbose=0;
  FILE *fpin=NULL;
  /*uint8_t*/ int bitj_x=0;
  /* uint32_t */ int nrows_x=0;
  /* uint32_t */ int ncols_x=0;
  int bitj=0,nrows=0,ncols=0,nrows_extend=0,ncols_extend=0;
  unsigned long wp_length=0;
  int current_row=0,nf=0,nc=0,nr=0,nrx=0;
  FILE *fpout=NULL;
  unsigned char *wkspace_mark=NULL;
  unsigned char *wp=NULL;
  unsigned char *ip=NULL;
  unsigned char b1=0,b2=0;
  char bstr1[9],bstr2[9],kstr[FNAMESIZE];
  if (verbose){ printf(" %% [finished wrap_transpose] filename %s --> %s nbytes_per_read %d\n",filename,outname,nbytes_per_read);}
  if (access(filename,F_OK)!=0){ printf(" Warning! cannot access %s in wrap_transpose\n",filename);}
  if ((fpin=fopen(filename,"r"))==NULL){ printf(" Warning! cannot open %s in wrap_transpose\n",filename);}
  fread(&bitj_x,sizeof(/*uint8_t*/ int),1,fpin);fread(&nrows_x,sizeof(/* uint32_t */ int),1,fpin);fread(&ncols_x,sizeof(/* uint32_t */ int),1,fpin); bitj=(int)bitj_x; nrows=(int)nrows_x; ncols=(int)ncols_x;
  if (fpin!=NULL){ fclose(fpin);fpin=NULL;}
  nbytes_per_read = BIT8*(nbytes_per_read/BIT8); if (nbytes_per_read==0){ nbytes_per_read=BIT8;}
  nrows_extend = (bitj - (nrows % bitj)) % bitj;
  ncols_extend = (bitj - (ncols % bitj)) % bitj;
  if (verbose){ printf(" %% read bitj %d nrows %d+%d ncols %d+%d\n",bitj,nrows,nrows_extend,ncols,ncols_extend);}
  if ((fpin=fopen(filename,"r"))==NULL){ printf(" Warning! cannot open %s in wrap_transpose\n",filename);}
  if ((fpout=fopen(outname,"w"))==NULL){ printf(" Warning! cannot open %s in wrap_transpose\n",outname);}
  bitj_x = (/*uint8_t*/ int)bitj; nrows_x = (/* uint32_t */ int)nrows; ncols_x = (/* uint32_t */ int)ncols;
  fwrite(&bitj_x,sizeof(/*uint8_t*/ int),1,fpout);fwrite(&ncols_x,sizeof(/* uint32_t */ int),1,fpout);fwrite(&nrows_x,sizeof(/* uint32_t */ int),1,fpout);
  wp_length = (ncols+ncols_extend) * nbytes_per_read; wp_length = rup(wp_length,sizeof(long)); wp_length /= BIT8;
  GLOBAL_wkspace_point = 0; wkspace_mark = wkspace_base; 
  wp = wkspace_all0c(wp_length); if (!wp){ printf(" Warning! not enough memory, reduce nbytes_per_read %d in wrap_transpose\n",nbytes_per_read);}
  fill_uchar_zero(wp,wp_length);
  ip = wkspace_all0c(nbytes_per_read / BIT8 + 1); if (!ip){ printf(" Warning! not enough memory, reduce nbytes_per_read %d in wrap_transpose\n",nbytes_per_read);}
  current_row=0;
  while (current_row<nrows){
    fseek(fpin,3*sizeof(int) + (current_row)/BIT8,SEEK_SET);
    fill_uchar_zero(wp,wp_length);
    nr = minimum(nbytes_per_read,nrows /* + nrows_extend */ - current_row);
    nrx = rup(nr,BIT8);
    nc=0;
    while (nc<ncols){ 
      fread(ip,sizeof(unsigned char),nrx/BIT8,fpin);
      if (verbose){
	getBinW(ip,kstr,nrx);
	printf("ip %s\n",kstr);
	/* if (verbose){ } */}
      nf=0;
      while (nf<nr){
	if (verbose){
	  printf(" nc %.6d current_row %.6d nf %.6d (nr %.6d), ",nc,current_row,nf,nr);
	  b1 = ip[nf/BIT8];getBinW(&b1,bstr1,8);
	  b2 = ((((ip[nf/BIT8] >> (7-(nf%BIT8))) & 1 ) << (7 - (nc%BIT8))));getBinW(&b2,bstr2,8);
	  printf("reading %s --> %s, ",bstr1,bstr2);
	  b1 = wp[nc/8 + nf*((ncols+ncols_extend)/BIT8)];getBinW(&b1,bstr1,8);
	  b2 = (wp[nc/8 + nf*((ncols+ncols_extend)/BIT8)] | ((((ip[nf/BIT8] >> (7-(nf%BIT8))) & 1 ) << (7 - (nc%BIT8)))));getBinW(&b2,bstr2,8);
	  printf("updating %s --> %s.\n",bstr1,bstr2);
	  /* if (verbose){ } */}
	wp[nc/8 + nf*((ncols+ncols_extend)/BIT8)] |= ((((ip[nf/BIT8] >> (7-(nf%BIT8))) & 1 ) << (7 - (nc%BIT8))));
	nf++;
	/* while (nf<nbytes_per_read){ } */}
      nc++;
      fseek(fpin,(nrows+nrows_extend)/BIT8 - nrx/BIT8,SEEK_CUR);
      /* while (nc<ncols){ } */}
    while (nc<ncols+ncols_extend){
      nf=0;
      while (nf<nr){
    	wp[nc/8 + nf*((ncols+ncols_extend)/BIT8)] &= (~0 ^ (1 << (7-(nc%BIT8))));
    	nf++;
    	/* while (nf<nbytes_per_read){ } */}
      nc++;
      /* while (nc<ncols){ } */}
    fwrite(wp,sizeof(unsigned char),(ncols+ncols_extend)/BIT8 * nr,fpout);
    current_row += nbytes_per_read;
    if (verbose){ printf("\r%d rows complete.", current_row); fflush(stdout);}
    /* while (current_row<nrows){ } */}
  if (verbose){ printf("\r finished writing %s\n",outname); fflush(stdout);}
  if (fpout!=NULL){ fclose(fpout);fpout=NULL;}
  if (fpin!=NULL){ fclose(fpin);fpin=NULL;}
  if (verbose){ printf(" %% [finished wrap_transpose] filename %s\n",filename);}
  wkspace_reset(wkspace_mark); GLOBAL_wkspace_point = 1;
}
