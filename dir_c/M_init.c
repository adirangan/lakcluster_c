/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

unsigned long long int M_wkspace_copy(struct M_handle *M,struct M_handle *M_in)
{
  /* copies memory to M from M_in. */
  int verbose=0;
  unsigned long long int length=0,length_total=0;
  int stp=0;  
  if (verbose){ printf(" %% [entering M_wkspace_copy] M->mc_length %d \n",(int)M->mc_length);}
  if (M->A_fp!=NULL && M->Ara==NULL){ 
    length = M->rpop_b*M->mc_length;
    memcpy(M->wX,M_in->wX,length); 
    length_total += length; 
    /* if A_fp */} 
  else if (M->A_fp==NULL && M->Ara!=NULL){ 
    length = 0; 
    M->wX=M->Ara; 
    length_total += length; 
    /* if Ara */} 
  else /* empty */ { length = 0; M->wX=NULL; length_total += length; } 
  stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->mr_length; memcpy(M->mr_b,M_in->mr_b,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->mr_length; memcpy(M->mr_j,M_in->mr_j,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->mc_length; memcpy(M->mc_b,M_in->mc_b,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->mc_length; memcpy(M->mc_j,M_in->mc_j,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->nrows*sizeof(double); memcpy(M->rsum,M_in->rsum,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->nrows*sizeof(unsigned int); memcpy(M->m_a_,M_in->m_a_,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->nrows*sizeof(unsigned int); memcpy(M->m_b_,M_in->m_b_,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->ncols*sizeof(unsigned int); memcpy(M->n_a_,M_in->n_a_,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->ncols*sizeof(unsigned int); memcpy(M->n_b_,M_in->n_b_,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  length = M->wt_length; memcpy(M->wt,M_in->wt,length); length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld <-- %lld\n",stp,length,length_total);}
  M->length_total = length_total;
  if (verbose){ printf(" %% [finished M_wkspace_copy] M->mc_length %d M->length_total \n",(int)M->mc_length,length_total); wkspace_printf();}
  return M->length_total;
}

unsigned long long int M_wkspace_alloc(struct M_handle *M,int estim_flag)
{
  /* Allocates memory to M, ;
     If estim_flag==1, no allocation is carried out but estimated space is computed. ;
     Otherwise, space is allocated (and filled with 0). ;
     The M->mark holds the wkspace_base at the start of the allocation (cannot free without disrupting wkspace_point) .; */
  int verbose=0;
  unsigned long long int length=0,length_total=0;
  int stp=0;  
  if (verbose){ printf(" %% [entering M_wkspace_alloc] e%d M->mc_length %d \n",estim_flag,(int)M->mc_length);}
  M->mark = wkspace_base;
  if (M->A_fp!=NULL && M->Ara==NULL){ 
    length = M->rpop_b*M->mc_length;
    if (!estim_flag){ 
      M->wX = wkspace_all0c(length); 
      if (!(M->wX)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}
      /* !estim_flag */} 
    length_total += length; 
    /* if A_fp */} 
  else if (M->A_fp==NULL && M->Ara!=NULL){ 
    length = 0; 
    M->wX=M->Ara; 
    length_total += length; 
    /* if Ara */} 
  else /* empty */ { length = 0; M->wX=NULL; length_total += length; } 
  stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  length = M->mr_length; if (!estim_flag){ M->mr_b = wkspace_all0c(length); if (!(M->mr_b)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  length = M->mr_length; if (!estim_flag){ M->mr_j = wkspace_all0c(length); if (!(M->mr_j)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  length = M->mc_length; if (!estim_flag){ M->mc_b = wkspace_all0c(length); if (!(M->mc_b)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  length = M->mc_length; if (!estim_flag){ M->mc_j = wkspace_all0c(length); if (!(M->mc_j)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  length = M->nrows*sizeof(double); if (!estim_flag){ M->rsum = (double *)wkspace_all0c(length); if (!(M->rsum)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  length = M->nrows*sizeof(unsigned int); if (!estim_flag){ M->m_a_ = (unsigned int *)wkspace_all0c(length); if (!(M->m_a_)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  length = M->nrows*sizeof(unsigned int); if (!estim_flag){ M->m_b_ = (unsigned int *)wkspace_all0c(length); if (!(M->m_b_)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  length = M->ncols*sizeof(unsigned int); if (!estim_flag){ M->n_a_ = (unsigned int *)wkspace_all0c(length); if (!(M->n_a_)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  length = M->ncols*sizeof(unsigned int); if (!estim_flag){ M->n_b_ = (unsigned int *)wkspace_all0c(length); if (!(M->n_b_)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  if (verbose){ printf(" %% stp %d length %lld, wkspace_left %lld (%lld MB)\n",stp,length,wkspace_left,wkspace_left/1048576);}
  M->wt_length = rup(length_total,POPLENGTH)-length_total; 
  length = M->wt_length; if (!estim_flag){ M->wt = wkspace_all0c(length); if (!(M->wt)){ printf(" %% Warning! stp %d out of memory in M_wkspace_alloc (wkspace_left %lld -- %lld MB)\n",stp,wkspace_left,wkspace_left/1048567);}} length_total += length; stp++;
  M->length_total=length_total;
  if (verbose){ printf(" %% [finished M_wkspace_alloc] e%d M->mc_length %d M->length_total %lld wkspace_left %lld (%lld MB)\n",estim_flag,(int)M->mc_length,length_total,wkspace_left,wkspace_left/1048576); wkspace_printf();}
  return M->length_total;
}

void M_handle_printf(struct M_handle *M,int verbose,char *prefix)
{
  /* prints out basic diagnostics for M_handle */
  char tmpchar[FNAMESIZE];
  unsigned char *bAra=NULL,*bX=NULL;
  int bitj=0,A_n_rows=0,A_n_cols=0;
  int nr_j=0,nr_b=0,nr_a=0,nc_j=0,nc_b=0,nc_a=0;
  if (M!=NULL){
    if (verbose){ printf("%sheader_length %d nr_a %d nc_a %d rpop_b %d(%d) cpop_b %d(%d)\n",prefix,M->header_length,M->nrows,M->ncols,(int)M->rpop_b,(int)M->rpop_j,(int)M->cpop_b,(int)M->cpop_j);}
    if (verbose){ printf("%smr_length %d mc_length %d\n",prefix,(int)M->mr_length,(int)M->mc_length);}
    if (verbose){ printf("%s wt_length %d length_total %lld\n",prefix,(int)M->wt_length,M->length_total);}
    if (verbose){ sprintf(tmpchar,"%s mr_b: ",prefix); bprintf(M->mr_b,M->bitj,1,M->nrows,tmpchar);}
    if (verbose){ sprintf(tmpchar,"%s mr_j: ",prefix); bprintf(M->mr_j,M->bitj,1,M->nrows,tmpchar);}
    if (verbose){ sprintf(tmpchar,"%s mc_b: ",prefix); bprintf(M->mc_b,M->bitj,1,M->ncols,tmpchar);}
    if (verbose){ sprintf(tmpchar,"%s mc_j: ",prefix); bprintf(M->mc_j,M->bitj,1,M->ncols,tmpchar);}
    if (verbose>1 && M->A_name!=NULL && M->A_fp!=NULL){ printf("%sA_name %s\n",prefix,M->A_name); binary_read(M->A_name,&bitj,&A_n_cols,&A_n_rows,&bAra); sprintf(tmpchar,"%s Afp: ",prefix); bprintf(bAra,bitj,A_n_rows,A_n_cols,tmpchar); }
    if (verbose>1 && M->Ara!=NULL){ sprintf(tmpchar,"%s Ara : ",prefix); bprintf(M->Ara,POPLENGTH,M->nrows,M->ncols,tmpchar);}
    if (verbose>1 && M->wX!=NULL){
      nr_j=0;nr_b=0;nr_a=0;
      while (nr_j<M->rpop_j && nr_b<M->rpop_b && nr_a<M->nrows){
	if (0){ /* do nothing */}
	else if (bget__on(M->mr_b,nr_a) && bget__on(M->mr_j,nr_a)){
	  printf("%s wX:",prefix);
	  nc_j=0;nc_b=0;nc_a=0;
	  while (nc_j<M->cpop_j && nc_b<M->cpop_b && nc_a<M->ncols){
	    if (0){ /* do nothing */}
	    else if (bget__on(M->mc_b,nc_a) && bget__on(M->mc_j,nc_a)){
	      //printf("%c",((int)(((M->wX[nc_a/BIT8 + nr_b*M->mc_length] >> (7-(nc_a%BIT8))) & 1) << 1)-(int)1)>0 ? '1':'0');
	      //printf("%c",bget__on(&(M->wX[0+nr_b*M->mc_length]),nc_a) ? '1' : '0');
	      bX = &(M->wX[0+nr_b*M->mc_length]); printf("%c",bget__on(bX,nc_a) ? 'X' : '0');
	      nc_j++;nc_b++;nc_a++;
	      /*  if (bget__on(M->mc_b,nc_a) && bget__on(M->mc_j,nc_a)){ } */}
	    else if (bget__on(M->mc_b,nc_a) && bget_off(M->mc_j,nc_a)){
	      nc_b++;nc_a++;
	      /*  if (bget__on(M->mc_b,nc_a) && bget_off(M->mc_j,nc_a)){ } */}
	    else if (bget_off(M->mc_b,nc_a)){
	      nc_a++;
	      /*  if (bget_off(M->mc_b,nc_a)){ } */}
	    /* while (nc_j<M->cpop_j && nc_b<M->cpop_b && nc_a<M->ncols){ } */}
	  printf("\n");
	  nr_j++;nr_b++;nr_a++;
	  /*  if (bget__on(M->mr_b,nr_a) && bget__on(M->mr_j,nr_a)){ } */}
	else if (bget__on(M->mr_b,nr_a) && bget_off(M->mr_j,nr_a)){
	  nr_b++;nr_a++;
	  /*  if (bget__on(M->mr_b,nr_a) && bget_off(M->mr_j,nr_a)){ } */}
	else if (bget_off(M->mr_b,nr_a)){
	  nr_a++;
	  /*  if (bget_off(M->mr_b,nr_a)){ } */}
	/* while (nr_j<M->rpop_j && nr_b<M->rpop_b && nr_a<M->nrows){ } */}
      /* if (verbose && M->wX!=NULL){ } */}
    /* if (M!=NULL){ } */}
}

struct M_handle * M_handle_make(int bitj,int nrows,int ncols,char *A_filename,unsigned char *Ara)
{
  /* designed for At*X; A stores a transpose of the first matrix in the product; 
     presumably X is stored in another M_handle ; We also assume the bitj of these two M_handles is the same ; */
  int verbose=0;
  struct M_handle *M = (struct M_handle *)wkspace_all0c(1*sizeof(struct M_handle));
  M->header_length = sizeof(/*uint8_t*/ int) + sizeof(/* uint32_t */ int) + sizeof(/* uint32_t */ int) /* 3*sizeof(int) */ /* bytes */;
  M->bitj=bitj;M->nrows=nrows;M->ncols=ncols;
  if ((A_filename!=NULL && strcmp(A_filename,"\0")) && Ara==NULL){
    sprintf(M->A_name,"%s",A_filename);
    if ((M->A_fp=fopen(A_filename,"r"))==NULL){ printf(" %% Warning! can't open %s in M_handle_make\n",A_filename);}
    fseeko(M->A_fp,(off_t)(M->header_length),SEEK_SET);
    binary_read_getsize(M->A_name,&(M->bitj),&(M->ncols),&(M->nrows));
    M->Ara=NULL;
    /* if ((A_filename!=NULL && strcmp(A_filename,"\0")) && Ara=NULL){ } */}
  else if ((A_filename==NULL || !strcmp(A_filename,"\0")) && Ara!=NULL){
    sprintf(M->A_name,"%s","\0"); M->A_fp=NULL;
    M->bitj = bitj; M->ncols=ncols; M->nrows=nrows; M->Ara=Ara; 
    /* else if ((A_filename==NULL || !strcmp(A_filename,"\0")) && Ara!=NULL){ } */}
  else if ((A_filename==NULL || !strcmp(A_filename,"\0")) && Ara==NULL){
    sprintf(M->A_name,"%s","\0"); M->A_fp=NULL;
    M->bitj = bitj; M->ncols=ncols; M->nrows=nrows; M->Ara=wkspace_all0c(1);
    /* else if ((A_filename==NULL || !strcmp(A_filename,"\0")) && Ara!=NULL){ } */}
  M->nrows_extend = (M->bitj - (M->nrows % M->bitj)) % M->bitj; M->mr_length = bsize(M->nrows)/* (rup(M->nrows+M->nrows_extend,POPLENGTH))/BIT8 */; 
  M->ncols_extend = (M->bitj - (M->ncols % M->bitj)) % M->bitj; M->mc_length = bsize(M->ncols)/* (rup(M->ncols+M->ncols_extend,POPLENGTH))/BIT8 */; 
  M->cc_length = (M->ncols + M->ncols_extend)/BIT8; /* used for loading columns from M_>A_fp */
  if (verbose){ printf(" bitj %d rows %d+%d cols %d+%d mr_length %d mc_length %d \n",(int)M->bitj,(int)M->nrows,(int)M->nrows_extend,(int)M->ncols,(int)M->ncols_extend,(int)M->mr_length,(int)M->mc_length); wkspace_printf(); /* if (verbose){ } */}
  M->rsum=NULL;
  M->m_a_=NULL;
  M->m_b_=NULL;
  M->m_b_=NULL;
  M->n_b_=NULL;
  return M;
}

void M_handle_copy(struct M_handle *M,struct M_handle *M_in)
{
  int verbose=0;
  /* copies M from M_in ; assumes M preallocated. */
  if (M!=NULL && M_in!=NULL){
    M->header_length = M_in->header_length; M->bitj = M_in->bitj; M->nrows = M_in->nrows; M->ncols = M_in->ncols;
    if ((M_in->A_name!=NULL && strcmp(M_in->A_name,"\0")) && M_in->Ara==NULL){
      sprintf(M->A_name,"%s",M_in->A_name);
      if ((M->A_fp=fopen(M->A_name,"r"))==NULL){ printf(" %% Warning! can't open %s in M_handle_copy\n",M->A_name);}
      fseeko(M->A_fp,(off_t)(M->header_length),SEEK_SET);
      binary_read_getsize(M->A_name,&(M->bitj),&(M->ncols),&(M->nrows));
      M->Ara=NULL;
      /* if M_in->A_name */}
    else if ((M_in->A_name==NULL || !strcmp(M_in->A_name,"\0")) && M_in->Ara!=NULL){
      sprintf(M->A_name,"%s","\0"); M->A_fp=NULL;
      M->bitj = M_in->bitj; M->ncols=M_in->ncols; M->nrows=M_in->nrows; M->Ara=M_in->Ara; 
      /* if M_in->Ara */}
    else if ((M_in->A_name==NULL || !strcmp(M_in->A_name,"\0")) && M_in->Ara==NULL){
      sprintf(M->A_name,"%s","\0"); M->A_fp=NULL;
      M->bitj = M_in->bitj; M->ncols=M_in->ncols; M->nrows=M_in->nrows; M->Ara=wkspace_all0c(1);
      /* if neither */}
    M->nrows_extend = (M->bitj - (M->nrows % M->bitj)) % M->bitj; M->mr_length = bsize(M->nrows)/* (rup(M->nrows+M->nrows_extend,POPLENGTH))/BIT8 */; 
    M->ncols_extend = (M->bitj - (M->ncols % M->bitj)) % M->bitj; M->mc_length = bsize(M->ncols)/* (rup(M->ncols+M->ncols_extend,POPLENGTH))/BIT8 */; 
    M->cc_length = (M->ncols + M->ncols_extend)/BIT8; /* used for loading columns from M_>A_fp */
    if (verbose){ printf(" bitj %d rows %d+%d cols %d+%d mr_length %d mc_length %d \n",(int)M->bitj,(int)M->nrows,(int)M->nrows_extend,(int)M->ncols,(int)M->ncols_extend,(int)M->mr_length,(int)M->mc_length); wkspace_printf(); /* if (verbose){ } */}
    M_wkspace_copy(M,M_in);
    /* if (M!=NULL && M_in!=NULL){ } */}
}

void M_fp_free(struct M_handle *M){ /* frees fp within M */ if (M->A_fp!=NULL){ fclose(M->A_fp);M->A_fp=NULL;}}

void M_load(struct M_handle *M)
{
  /* This function copies M_A_fp into M_wX */
  int verbose=0;
  int nr=0,nr_b=0;
  if (verbose>1){ printf(" %% [entering M_load]\n");}
  nr=0; nr_b=0; if (M->A_fp!=NULL){ fseeko(M->A_fp,(off_t)(M->header_length),SEEK_SET);} 
  while (nr<M->nrows && nr_b<M->rpop_b){
    while (bget_off(M->mr_b,nr) && nr<M->nrows && nr_b<M->rpop_b){ if (verbose>2){ printf(" %% jumping row %d\n",nr);} if (M->A_fp!=NULL){ fseeko(M->A_fp,M->cc_length,SEEK_CUR);} nr++;}
    if (nr<M->nrows && nr_b<M->rpop_b){ if (verbose>2){ printf(" %% reading row %d (%d)\n",nr,nr_b);} if (M->A_fp!=NULL){ fread(&(M->wX[nr_b*M->mc_length]),1,M->cc_length,M->A_fp);} nr_b++; nr++;}
    /* while (nr<M->nrows){ } */}
  if (verbose>3){ bprintf(M->wX,POPLENGTH,M->rpop_b,M->ncols," %% M->wX: ");} 
  if (verbose>1){ printf(" %% [finished M_load]\n"); wkspace_printf();}
}

void M_mxget_excerpt(int verbose,int nrows,unsigned char *mr_b,unsigned char *mr_j,long long int *rpop_b_p,long long int *rpop_j_p,unsigned int *m_a_,unsigned int *m_b_)
{
  int mr_length = bsize(nrows);
  long long int rpop_b=0,rpop_j=0;
  int m_a=0,m_b=0,m_j=0;
  rpop_b = popcount_uchar_array(mr_b,mr_length); rpop_j = popcount_uchar_array(mr_j,mr_length);
  if (verbose>1){ printf(" %% rpop_b %d rpop_j %d nrows %d\n",rpop_b,rpop_j,nrows);}
  m_a=0;m_b=0;m_j=0;
  while (m_a<nrows){
    while (m_a<nrows && bget_off(mr_b,m_a)){ /* stop when mr_b[m_a]==1 */ m_a++;}
    if (m_a<nrows && bget__on(mr_b,m_a) && bget__on(mr_j,m_a)){ 
      m_b_[m_j] = m_b;
      m_a_[m_j] = m_a;
      m_j++; m_b++; m_a++;
      /* if mr_j[m_a]==1 */}
    else if (m_a<nrows && bget__on(mr_b,m_a) && bget_off(mr_j,m_a)){ 
      if (verbose>2){ printf(" %% m_a %d m_b %d not processing...\n",m_a,m_b);} 
      m_b++; m_a++; 
      /* else */}
    else if (m_a<nrows && bget_off(mr_b,m_a)){ 
      if (verbose){ printf(" %% Warning! mr_b changed in M_mxget\n");} 
      m_a++; 
      /* else */}
    else if (m_a>=nrows){ /* finished */}
    /* while (m_a<nrows){ } */}
  if (rpop_b_p!=NULL){ *rpop_b_p = rpop_b;}
  if (rpop_j_p!=NULL){ *rpop_j_p = rpop_j;}
}

void M_mxget(struct M_handle *M)
{
  /* 
     M->m_b_[mr_j] = position within M->Xra  of mr_j^th 1 of M->mr_j & M->mr_b;
     M->m_a_[mr_j] = position within M->mr_b of mr_j^th 1 of M->mr_j & M->mr_b;
     M->n_b_[nr_j] = position within M->Xra  of nr_j^th 1 of M->mc_j & M->mc_b;
     M->n_a_[nr_j] = position within M->mc_b of nr_j^th 1 of M->mc_j & M->mc_b;
  */
  int verbose=0;
  int nw_a=0,nw_b=0,nw_j=0;
  if (verbose){ printf(" %% [entering M_mxget] \n");}
  if (verbose>1){ bprintf(M->mr_b,M->bitj,1,M->nrows," %% M->mr_b: ");}
  if (verbose>1){ bprintf(M->mr_j,M->bitj,1,M->nrows," %% M->mr_j: ");}
  M_mxget_excerpt(verbose,M->nrows,M->mr_b,M->mr_j,&(M->rpop_b),&(M->rpop_j),M->m_a_,M->m_b_);
  if (verbose>1){ raprintf(M->m_b_,"unsigned int",1,M->rpop_j," %% M->m_b_ : ");}
  if (verbose>1){ raprintf(M->m_a_,"unsigned int",1,M->rpop_j," %% M->m_a_ : ");}
  if (verbose>1){ bprintf(M->mc_b,M->bitj,1,M->ncols," %% M->mc_b: ");}
  if (verbose>1){ bprintf(M->mc_j,M->bitj,1,M->ncols," %% M->mc_j: ");}
  M_mxget_excerpt(verbose,M->ncols,M->mc_b,M->mc_j,&(M->cpop_b),&(M->cpop_j),M->n_a_,M->n_b_);
  if (verbose>1){ raprintf(M->n_b_,"unsigned int",1,M->cpop_j," %% M->n_b_ : ");}
  if (verbose>1){ raprintf(M->n_a_,"unsigned int",1,M->cpop_j," %% M->n_a_ : ");}
  if (verbose){ printf(" %% [finished M_mxget]\n"); wkspace_printf();}
}

void M_mxget_bkp(struct M_handle *M)
{
  /* 
     M->m_b_[mr_j] = position within M->Xra  of mr_j^th 1 of M->mr_j & M->mr_b;
     M->m_a_[mr_j] = position within M->mr_b of mr_j^th 1 of M->mr_j & M->mr_b;
     M->n_b_[nr_j] = position within M->Xra  of nr_j^th 1 of M->mc_j & M->mc_b;
     M->n_a_[nr_j] = position within M->mc_b of nr_j^th 1 of M->mc_j & M->mc_b;
  */
  int verbose=0;
  int nw_a=0,nw_b=0,nw_j=0;
  if (verbose){ printf(" %% [entering M_mxget] \n");}
  if (verbose>1){ bprintf(M->mr_b,M->bitj,1,M->nrows," %% M->mr_b: ");}
  if (verbose>1){ bprintf(M->mr_j,M->bitj,1,M->nrows," %% M->mr_j: ");}
  M->rpop_b = popcount_uchar_array(M->mr_b,M->mr_length); M->rpop_j = popcount_uchar_array(M->mr_j,M->mr_length);
  if (verbose>1){ printf(" %% rpop_b %d rpop_j %d nrows %d\n",M->rpop_b,M->rpop_j,M->nrows);}
  M->cpop_b = popcount_uchar_array(M->mc_b,M->mc_length); M->cpop_j = popcount_uchar_array(M->mc_j,M->mc_length);
  if (verbose>1){ printf(" %% cpop_b %d cpop_j %d ncols %d\n",M->cpop_b,M->cpop_j,M->ncols);}
  nw_a=0;nw_b=0;nw_j=0;
  while (nw_a<M->nrows){
    while (nw_a<M->nrows && bget_off(M->mr_b,nw_a)){ /* stop when M->mr_b[nw_a]==1 */ nw_a++;}
    if (nw_a<M->nrows && bget__on(M->mr_b,nw_a) && bget__on(M->mr_j,nw_a)){ 
      M->m_b_[nw_j] = nw_b;
      M->m_a_[nw_j] = nw_a;
      nw_j++; nw_b++; nw_a++;
      /* if M->mr_j[nw_a]==1 */}
    else if (nw_a<M->nrows && bget__on(M->mr_b,nw_a) && bget_off(M->mr_j,nw_a)){ 
      if (verbose>2){ printf(" %% nw_a %d nw_b %d not processing...\n",nw_a,nw_b);} 
      nw_b++; nw_a++; 
      /* else */}
    else if (nw_a<M->nrows && bget_off(M->mr_b,nw_a)){ 
      if (verbose){ printf(" %% Warning! M->mr_b changed in M_mxget\n");} 
      nw_a++; 
      /* else */}
    else if (nw_a>=M->nrows){ /* finished */}
    /* while (nw_a<M->nrows){ } */}
  if (verbose>1){ raprintf(M->m_b_,"unsigned int",1,M->rpop_j," %% M->m_b_ : ");}
  if (verbose>1){ raprintf(M->m_a_,"unsigned int",1,M->rpop_j," %% M->m_a_ : ");}
  if (verbose>1){ bprintf(M->mc_b,M->bitj,1,M->ncols," %% M->mc_b: ");}
  if (verbose>1){ bprintf(M->mc_j,M->bitj,1,M->ncols," %% M->mc_j: ");}
  nw_a=0;nw_b=0;nw_j=0;
  while (nw_a<M->ncols){
    while (nw_a<M->ncols && bget_off(M->mc_b,nw_a)){ /* stop when M->mc_b[nw_a]==1 */ nw_a++;}
    if (nw_a<M->ncols && bget__on(M->mc_b,nw_a) && bget__on(M->mc_j,nw_a)){ 
      M->n_b_[nw_j] = nw_b;
      M->n_a_[nw_j] = nw_a;
      nw_j++; nw_b++; nw_a++;
      /* if M->mc_j[nw_a]==1 */}
    else if (nw_a<M->ncols && bget__on(M->mc_b,nw_a) && bget_off(M->mc_j,nw_a)){ 
      if (verbose>2){ printf(" %% nw_a %d nw_b %d not processing...\n",nw_a,nw_b);} 
      nw_b++; nw_a++; 
      /* else */}
    else if (nw_a<M->ncols && bget_off(M->mc_b,nw_a)){ 
      if (verbose){ printf(" %% Warning! M->mc_b changed in M_mxget\n");} 
      nw_a++; 
      /* else */}
    else if (nw_a>=M->ncols){ /* finished */}
    /* while (nw_a<M->ncols){ } */}
  if (verbose>1){ raprintf(M->n_b_,"unsigned int",1,M->cpop_j," %% M->n_b_ : ");}
  if (verbose>1){ raprintf(M->n_a_,"unsigned int",1,M->cpop_j," %% M->n_a_ : ");}
  if (verbose){ printf(" %% [finished M_mxget]\n"); wkspace_printf();}
}

struct M_handle * M_handle_v_make(int bitj,int nrows,int ncols,char *A_filename,unsigned char *Ara,unsigned char *mr_b,unsigned char *mc_b)
{
  /* This function generates an M_handle for eventually running a series of parallel calls to get_AX_vv;
     We assume that A is a fixed binary array of +1/-1;
     The rows of A are stored in succession in either A_filename or Ara;
     The relevant masks are stored in mr_b, mc_b; we assume that mr_j and mc_j (to be put in later) can change arbitrarily. ;
     However, we assume that all future masks mr_j and mc_j will be subsets of the masks mr_b,mc_b. ;
     As a result, we only load the rows of A_filename corresponding to the nonzero entries of mr_b (of quantity rpop_b<=nrows). ;
     Ultimately, we expect M->wX to point to contiguous memory of the form:
     [<-- first row indicated in M->mr_b -->],[<-- second row indicated in M->mr_b -->],[...],[<-- last row indicated in M->mr_b -->]
     with each row of length M->mc_length (equal to M->ncols rounded up to the nearest POPLENGTH)
     Note that if A_filename is provided, it is carefully loaded into M->wX row by row.
     However, if Ara is provided it is simply linked to by M->wX; no loading is performed.
     This means that Ara should be prepared with the mask mr_b in mind.
   */
  int verbose=0;
  unsigned char *wkspace_mark=NULL;
  int nf=0,nr=0,nc=0,ira[MAX_THREADS];
  struct M_handle *M=NULL;
  /* int ipop = sizeof(int)*BIT8; wX_mult *= ipop; */
  if (verbose){ printf(" %% [entering M_handle_v_make] bitj %d nrows %d ncols %d\n",bitj,nrows,ncols);}
  wkspace_mark=wkspace_base;
  for (nf=0;nf<MAX_THREADS;nf++){ ira[nf]=nf;}
  M = M_handle_make(bitj,nrows,ncols,A_filename,Ara); if (mr_b){ M->rpop_b = popcount_uchar_array(mr_b,M->mr_length);} else{ M->rpop_b = M->nrows;} M_wkspace_alloc(M,0);
  if (mr_b!=NULL){ for (nr=0;nr<M->mr_length;nr++){ M->mr_b[nr] = mr_b[nr];}} else{ for (nr=0;nr<M->nrows;nr++){ bset__on(M->mr_b,nr);}}
  if (mr_b!=NULL){ for (nr=0;nr<M->mr_length;nr++){ M->mr_j[nr] = mr_b[nr];}} else{ for (nr=0;nr<M->nrows;nr++){ bset__on(M->mr_j,nr);}}
  if (mc_b!=NULL){ for (nc=0;nc<M->mc_length;nc++){ M->mc_b[nc] = mc_b[nc];}} else{ for (nc=0;nc<M->ncols;nc++){ bset__on(M->mc_b,nc);}}
  if (mc_b!=NULL){ for (nc=0;nc<M->mc_length;nc++){ M->mc_j[nc] = mc_b[nc];}} else{ for (nc=0;nc<M->ncols;nc++){ bset__on(M->mc_j,nc);}}
  M->rpop_b = popcount_uchar_array(M->mr_b,M->mr_length);
  M->rpop_j = popcount_uchar_array(M->mr_j,M->mr_length);
  M->cpop_b = popcount_uchar_array(M->mc_b,M->mc_length);
  M->cpop_j = popcount_uchar_array(M->mc_j,M->mc_length);
  if (verbose>1){ printf(" %% (%d,%d,%d)-x-(%d,%d,%d)\n",M->nrows,M->rpop_b,M->rpop_j,M->ncols,M->cpop_b,M->cpop_j);}
  if (verbose>2){ bprintf(M->mr_b,M->bitj,1,M->nrows," %% M->mr_b: ");}
  if (verbose>2){ bprintf(M->mr_j,M->bitj,1,M->nrows," %% M->mr_j: ");}
  if (verbose>2){ bprintf(M->mc_b,M->bitj,1,M->ncols," %% M->mc_b: ");}
  if (verbose>2){ bprintf(M->mc_j,M->bitj,1,M->ncols," %% M->mc_j: ");}
  M_load(M);
  if (verbose){ printf(" %% [finished M_handle_v_make]\n"); wkspace_printf();}
  return M;
}

struct M_handle * M_handle_w_make(int bitj,int b_mlt,int nrows_bZ,int ncols_bZ_tmp)
{
  /* This function generates an M_handle for use with xcalc operations ;
   */
  int ncols_bZ_per=0,brows_bZ=0,ncols_bZ,length_bZ;
  unsigned char **bZ_p=NULL;
  struct M_handle *M=NULL;
  ncols_bZ_per = 1+b_mlt; brows_bZ = bsize(nrows_bZ); ncols_bZ=(ncols_bZ_per*ncols_bZ_tmp); length_bZ = brows_bZ*ncols_bZ;
  bZ_p = (unsigned char **)wkspace_all0c(1*sizeof(unsigned char *)); *bZ_p = wkspace_all0c(length_bZ);
  M = M_handle_v_make(BITJ,ncols_bZ,nrows_bZ,NULL,(*bZ_p),NULL,NULL);
  return M;
}
