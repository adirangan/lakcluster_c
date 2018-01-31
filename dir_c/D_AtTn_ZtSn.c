void *get_D_AtTn_ZtSn_vv(void *vp)
{
  /* This function takes in M_At, M_St, M_Zt, lf_AtTn, lf_ZtSn, and calculates output[nk+ns*nrows_An] = (AtTn)(nk,ns).*D(nk).*(ZtSn)(nk,ns) ;
     The output should look like: ;
     ((At(:,:)-a_An*e_At)*T(:,ns)).*diag(D_An(:)).*((Zt(:,:)-a_An*e_Zt)*S(:,ns)) ;
     = ;
     + (At(:,:)*T(:,ns)).*diag(D_An(:)).*(Zt(:,:)*S(:,ns)) ;
     - (a_An*e_At*T(:,ns)).*diag(D_An(:)).*(Zt(:,:)*S(:,ns)) ;
     - (At(:,:)*T(:,ns)).*diag(D_An(:)).*(a_An*e_Zt*S(:,ns)) ;
     + (a_An*e_At*T(:,ns)).*diag(D_An(:)).*(a_An*e_Zt*S(:,ns)) ;
     In the situation where M_An==M_Zn, then we can also correct for repeated indices within An:
     D_AtTn_AtTn = ((At(:,:)-a_An*e_At)*T(:,ns)).*diag(D_An(:)).*((Zt(:,:)-a_An*e_Zt)*S(:,ns)) ;
     - D_An(nj)*( At->cpop_j - 2*popcount(A(:,nj))*a_An(nj) + At->cpop_j*a_An(j)^2) ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *lf_AtTn = (struct L_handle *)(vpra[ip++]); 
  struct L_handle *lf_ZtSn = (struct L_handle *)(vpra[ip++]); 
  struct L_handle *output_D_AtTn_ZtSn_vv = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); int output_spacing_At = output_spacing_a; 
  int output_spacing_s = *(int *)(vpra[ip++]); int output_spacing_Tt = output_spacing_s;
  int ncols_A_p = psize(M_At->nrows) ;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*ncols_A_p]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*ncols_A_p]);
  int mx_j=0,mx_chunk=0,ns_j=0,ns_b=0,ns_a=0,tab_Tt_stride=0,tab_Tt=0,na_j=0,na_b=0,na_a=0,tab_At_stride=0,tab_At=0;
  double output_tmp=0; 
  __m128i *wA_tag; __m128i *wS_tag; __m128i *wT_tag; __m128i *mc_tag,*mc_end;
  double dtmp_etTn=0,dtmp_etSn=0;
  if (verbose>1){ printf(" %% [entering get_D_AtTn_ZtSn_vv] tidx %d \n",tidx);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  switch (output_spacing_Tt){ case SPACING_j: tab_Tt_stride = M_Tt->rpop_j; break; case SPACING_b: tab_Tt_stride = M_Tt->rpop_b; break; case SPACING_a: tab_Tt_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_Tt){ } */}
  switch (output_spacing_At){ case SPACING_j: tab_At_stride = M_At->rpop_j; break; case SPACING_b: tab_At_stride = M_At->rpop_b; break; case SPACING_a: tab_At_stride = M_At->nrows; break; default: break; /* switch (output_spacing_At){ } */}
  output_D_AtTn_ZtSn_vv->spacing_row = output_spacing_At; output_D_AtTn_ZtSn_vv->row_stride = tab_At_stride;
  output_D_AtTn_ZtSn_vv->spacing_col = output_spacing_Tt; output_D_AtTn_ZtSn_vv->col_stride = tab_Tt_stride;
  fill_uchar_zero((unsigned char *)(output_D_AtTn_ZtSn_vv->lf),tab_At_stride*tab_Tt_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"D_AtTn_ZtSn_vv")){ goto skip_D_AtTn_ZtSn_vv;}
  if (lf_AtTn!=NULL && lf_ZtSn!=NULL){
    if (verbose>2){ lfprintf(lf_AtTn," %% lf_AtTn: ");} 
    if (verbose>2){ lfprintf(lf_ZtSn," %% lf_ZtSn: ");} 
    if (GLOBAL_omp_type==GLOBAL_omp_off){
      ns_j=0; while (ns_j<M_St->rpop_j){
	ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	switch (output_spacing_Tt){ case SPACING_j: tab_Tt=ns_j; break; case SPACING_b: tab_Tt=ns_b; break; case SPACING_a: tab_Tt=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
	dtmp_etSn=0;
	wS_tag = (__m128i*)((unsigned long *)(&(M_St->wX[ns_b*M_St->mc_length])));
	mc_tag = (__m128i*)((unsigned long *)(&(M_St->mc_j[0])));
	mc_end = (__m128i*)((unsigned long *)(&(M_St->mc_j[M_St->mc_length])));
	dtmp_etSn = -M_St->cpop_j + 2*popcount(&wS_tag,&mc_tag,&mc_end); 
	if (verbose>2){ printf(" %% ns %d,%d,%d, dtmp_etSn %0.2f\n",ns_j,ns_b,ns_a,dtmp_etSn);}
	dtmp_etTn=0;
	wT_tag = (__m128i*)((unsigned long *)(&(M_Tt->wX[ns_b*M_Tt->mc_length])));
	mc_tag = (__m128i*)((unsigned long *)(&(M_Tt->mc_j[0])));
	mc_end = (__m128i*)((unsigned long *)(&(M_Tt->mc_j[M_Tt->mc_length])));
	dtmp_etTn = -M_Tt->cpop_j + 2*popcount(&wT_tag,&mc_tag,&mc_end); 
	if (verbose>2){ printf(" %% ns %d,%d,%d, dtmp_etTn %0.2f\n",ns_j,ns_b,ns_a,dtmp_etTn);}
	na_j=0; while (na_j<M_At->rpop_j){
	  na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
	  switch (output_spacing_At){ case SPACING_j: tab_At=na_j; break; case SPACING_b: tab_At=na_b; break; case SPACING_a: tab_At=na_a; break; default: break; /* switch (output_spacing_At){ } */}
	  output_tmp=0;
	  output_tmp = (*L2_get(lf_AtTn,na_j,na_b,na_a,ns_j,ns_b,ns_a) - dtmp_etTn*a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*(*L2_get(lf_ZtSn,na_j,na_b,na_a,ns_j,ns_b,ns_a) - a_An[na_a/POPLENGTH]*dtmp_etSn);
	  output_D_AtTn_ZtSn_vv->lf[tab_At + tab_Tt*tab_At_stride] = output_tmp;
	  na_j++; /* while (na_j<M_At->rpop_j){ } */}
	ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
      GLOBAL_ops_count_one(tidx,M_At->rpop_j*M_Tt->rpop_j,M_At->rpop_j*M_Tt->rpop_j*M_At->mc_length*BIT8);
      /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
    else if (GLOBAL_omp_type==GLOBAL_omp__on){
      mx_chunk=1; 
#pragma omp parallel private(mx_j,ns_j,ns_a,ns_b,tab_Tt,wS_tag,wT_tag,mc_tag,mc_end,dtmp_etSn,dtmp_etTn,na_j,na_a,na_b,tab_At,wA_tag,output_tmp)
      { /* begin omp parallel */
	mx_j=0;ns_j=0;ns_a=0;ns_b=0;tab_Tt=0;wS_tag=NULL;wT_tag=NULL;mc_tag=NULL;mc_end=NULL;dtmp_etSn=0;na_j=0;na_a=0;na_b=0;tab_At=0;wA_tag=NULL;output_tmp=0;
#pragma omp for schedule(static)
	for (mx_j=0;mx_j<M_St->rpop_j*M_At->rpop_j;mx_j++){
	  ns_j = mx_j/M_At->rpop_j; na_j = mx_j % M_At->rpop_j;
	  ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	  switch (output_spacing_Tt){ case SPACING_j: tab_Tt=ns_j; break; case SPACING_b: tab_Tt=ns_b; break; case SPACING_a: tab_Tt=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
	  dtmp_etSn=0;
	  wS_tag = (__m128i*)((unsigned long *)(&(M_St->wX[ns_b*M_St->mc_length])));
	  mc_tag = (__m128i*)((unsigned long *)(&(M_St->mc_j[0])));
	  mc_end = (__m128i*)((unsigned long *)(&(M_St->mc_j[M_St->mc_length])));
	  dtmp_etSn = -M_St->cpop_j + 2*popcount(&wS_tag,&mc_tag,&mc_end); 
	  dtmp_etTn=0;
	  wT_tag = (__m128i*)((unsigned long *)(&(M_Tt->wX[ns_b*M_Tt->mc_length])));
	  mc_tag = (__m128i*)((unsigned long *)(&(M_Tt->mc_j[0])));
	  mc_end = (__m128i*)((unsigned long *)(&(M_Tt->mc_j[M_Tt->mc_length])));
	  dtmp_etTn = -M_Tt->cpop_j + 2*popcount(&wT_tag,&mc_tag,&mc_end); 
	  na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
	  switch (output_spacing_At){ case SPACING_j: tab_At=na_j; break; case SPACING_b: tab_At=na_b; break; case SPACING_a: tab_At=na_a; break; default: break; /* switch (output_spacing_At){ } */}
	  output_tmp = (*L2_get(lf_AtTn,na_j,na_b,na_a,ns_j,ns_b,ns_a) - dtmp_etTn*a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*(*L2_get(lf_ZtSn,na_j,na_b,na_a,ns_j,ns_b,ns_a) - a_An[na_a/POPLENGTH]*dtmp_etSn);
	  output_D_AtTn_ZtSn_vv->lf[tab_At + tab_Tt*tab_At_stride] = output_tmp;
	  /* for (mx_j=0;mx_j<M_St->rpop_j*M_At->rpop_j;mx_j++){ } */}
	/* end omp parallel */}
      GLOBAL_ops_count_one(tidx,M_At->rpop_j*M_Tt->rpop_j,M_At->rpop_j*M_Tt->rpop_j*M_At->mc_length*BIT8);
      /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
    /* if (structures exist){ } */}
  else{ printf(" %% Warning! structures lf_AtTn and lf_ZtSn does not exist in call to get_D_AtTn_ZtSn_vv\n");}
  if (verbose>2){ lfprintf(output_D_AtTn_ZtSn_vv," %% output_D_AtTn_ZtSn_vv: ");}
 skip_D_AtTn_ZtSn_vv:
  if (verbose>1){ printf(" %% [finished get_D_AtTn_ZtSn_vv] tidx %d\n",tidx);}
  return NULL;
}

int wrap_D_AtTn_ZtSn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_ZtSn,struct L_handle **output_p)
 {
   /* This function calls get_D_AtTn_ZtSn_vv ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 11)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_D_AtTn_ZtSn_vv__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Zt,verbose," %% M_Zt: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_D_AtTn_ZtSn_vv__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Zt; vpra[ip++] = M_St; vpra[ip++] = A_ajdk; vpra[ip++] = lf_AtTn; vpra[ip++] = lf_ZtSn; vpra[ip++] = *output_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_D_AtTn_ZtSn_vv,vpra)){ printf("Warning! cannot create thread %d in wrap_D_AtTn_ZtSn_vv\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_D_AtTn_ZtSn_vv(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_D_AtTn_ZtSn_vv__run] tidx %d\n",*tidx);}
  return length;
}  

void *get_D_AtTn_AtTn_vv(void *vp)
{
  /* This function takes in M_At, M_St==M_Tt, M_Zt==M_At, lf_AtTn, lf_ZtSn==lf_AtTn, and calculates output[nk+ns*nrows_An] = (AtTn)(nk,ns).*D(nk).*(ZtSn)(nk,ns) ;
     The output should look like: ;
     ((At(:,:)-a_An*e_At)*T(:,ns)).*diag(D_An(:)).*((Zt(:,:)-a_An*e_Zt)*S(:,ns)) ;
     = ;
     + (At(:,:)*T(:,ns)).*diag(D_An(:)).*(Zt(:,:)*S(:,ns)) ;
     - (a_An*e_At*T(:,ns)).*diag(D_An(:)).*(Zt(:,:)*S(:,ns)) ;
     - (At(:,:)*T(:,ns)).*diag(D_An(:)).*(a_An*e_Zt*S(:,ns)) ;
     + (a_An*e_At*T(:,ns)).*diag(D_An(:)).*(a_An*e_Zt*S(:,ns)) ;
     In the situation where M_An==M_Zn, then we can also correct for repeated indices within An:
     D_AtTn_AtTn = ((At(:,:)-a_An*e_At)*T(:,ns)).*diag(D_An(:)).*((Zt(:,:)-a_An*e_Zt)*S(:,ns)) ;
     - D_An(nj)*( At->cpop_j - 2*popcount(A(:,nj))*a_An(nj) + At->cpop_j*a_An(j)^2) ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *lf_AtTn = (struct L_handle *)(vpra[ip++]); 
  struct L_handle *lf_ZtSn = (struct L_handle *)(vpra[ip++]); 
  struct L_handle *output_D_AtTn_AtTn_vv = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); int output_spacing_At = output_spacing_a; 
  int output_spacing_s = *(int *)(vpra[ip++]); int output_spacing_Tt = output_spacing_s;
  int ncols_A_p = psize(M_At->nrows) ;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*ncols_A_p]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*ncols_A_p]);
  int mx_j=0,mx_chunk=0,ns_j=0,ns_b=0,ns_a=0,tab_Tt_stride=0,tab_Tt=0,na_j=0,na_b=0,na_a=0,tab_At_stride=0,tab_At=0;
  double output_tmp=0; 
  __m128i *wA_tag; __m128i *wS_tag; __m128i *wT_tag; __m128i *mc_tag,*mc_end;
  double dtmp_etTn=0,dtmp_etSn=0,dtmp_etAn=0;
  if (verbose>1){ printf(" %% [entering get_D_AtTn_AtTn_vv] tidx %d \n",tidx);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  switch (output_spacing_Tt){ case SPACING_j: tab_Tt_stride = M_Tt->rpop_j; break; case SPACING_b: tab_Tt_stride = M_Tt->rpop_b; break; case SPACING_a: tab_Tt_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_Tt){ } */}
  switch (output_spacing_At){ case SPACING_j: tab_At_stride = M_At->rpop_j; break; case SPACING_b: tab_At_stride = M_At->rpop_b; break; case SPACING_a: tab_At_stride = M_At->nrows; break; default: break; /* switch (output_spacing_At){ } */}
  output_D_AtTn_AtTn_vv->spacing_row = output_spacing_At; output_D_AtTn_AtTn_vv->row_stride = tab_At_stride;
  output_D_AtTn_AtTn_vv->spacing_col = output_spacing_Tt; output_D_AtTn_AtTn_vv->col_stride = tab_Tt_stride;
  fill_uchar_zero((unsigned char *)(output_D_AtTn_AtTn_vv->lf),tab_At_stride*tab_Tt_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"D_AtTn_AtTn_vv")){ goto skip_D_AtTn_AtTn_vv;}
  if (lf_AtTn!=NULL && lf_ZtSn!=NULL){
    if (verbose>2){ lfprintf(lf_AtTn," %% lf_AtTn: ");} 
    if (verbose>2){ lfprintf(lf_ZtSn," %% lf_ZtSn: ");} 
    if (GLOBAL_omp_type==GLOBAL_omp_off){
      ns_j=0; while (ns_j<M_St->rpop_j){
	ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	switch (output_spacing_Tt){ case SPACING_j: tab_Tt=ns_j; break; case SPACING_b: tab_Tt=ns_b; break; case SPACING_a: tab_Tt=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
	dtmp_etSn=0;
	wS_tag = (__m128i*)((unsigned long *)(&(M_St->wX[ns_b*M_St->mc_length])));
	mc_tag = (__m128i*)((unsigned long *)(&(M_St->mc_j[0])));
	mc_end = (__m128i*)((unsigned long *)(&(M_St->mc_j[M_St->mc_length])));
	dtmp_etSn = -M_St->cpop_j + 2*popcount(&wS_tag,&mc_tag,&mc_end); 
	if (verbose>2){ printf(" %% ns %d,%d,%d, dtmp_etSn %0.2f\n",ns_j,ns_b,ns_a,dtmp_etSn);}
	dtmp_etTn=0;
	wT_tag = (__m128i*)((unsigned long *)(&(M_Tt->wX[ns_b*M_Tt->mc_length])));
	mc_tag = (__m128i*)((unsigned long *)(&(M_Tt->mc_j[0])));
	mc_end = (__m128i*)((unsigned long *)(&(M_Tt->mc_j[M_Tt->mc_length])));
	dtmp_etTn = -M_Tt->cpop_j + 2*popcount(&wT_tag,&mc_tag,&mc_end); 
	if (verbose>2){ printf(" %% ns %d,%d,%d, dtmp_etTn %0.2f\n",ns_j,ns_b,ns_a,dtmp_etTn);}
	na_j=0; while (na_j<M_At->rpop_j){
	  na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
	  switch (output_spacing_At){ case SPACING_j: tab_At=na_j; break; case SPACING_b: tab_At=na_b; break; case SPACING_a: tab_At=na_a; break; default: break; /* switch (output_spacing_At){ } */}
	  output_tmp=0;
	  output_tmp = (*L2_get(lf_AtTn,na_j,na_b,na_a,ns_j,ns_b,ns_a) - dtmp_etTn*a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*(*L2_get(lf_ZtSn,na_j,na_b,na_a,ns_j,ns_b,ns_a) - a_An[na_a/POPLENGTH]*dtmp_etSn);
	  wA_tag = (__m128i *)(&(M_At->wX[na_b*M_At->mc_length]));
	  mc_tag = (__m128i *)(&(M_At->mc_j[0]));
	  mc_end = (__m128i *)(&(M_At->mc_j[M_At->mc_length]));
	  dtmp_etAn = -M_At->cpop_j + 2*popcount(&wA_tag,&mc_tag,&mc_end);
	  output_tmp -= D_An[na_a/POPLENGTH]*(M_At->cpop_j - 2*dtmp_etAn*a_An[na_a/POPLENGTH] + M_At->cpop_j*a_An[na_a/POPLENGTH]*a_An[na_a/POPLENGTH]);
	  output_D_AtTn_AtTn_vv->lf[tab_At + tab_Tt*tab_At_stride] = output_tmp;
	  na_j++; /* while (na_j<M_At->rpop_j){ } */}
	ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
      GLOBAL_ops_count_one(tidx,M_At->rpop_j*M_Tt->rpop_j,M_At->rpop_j*M_Tt->rpop_j*M_At->mc_length*BIT8);
      /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
    else if (GLOBAL_omp_type==GLOBAL_omp__on){
      mx_chunk=1; 
#pragma omp parallel private(mx_j,ns_j,ns_a,ns_b,tab_Tt,wS_tag,wT_tag,mc_tag,mc_end,dtmp_etSn,dtmp_etTn,na_j,na_a,na_b,tab_At,wA_tag,dtmp_etAn,output_tmp)
      { /* begin omp parallel */
	mx_j=0;ns_j=0;ns_a=0;ns_b=0;tab_Tt=0;wS_tag=NULL;wT_tag=NULL;mc_tag=NULL;mc_end=NULL;dtmp_etSn=0;na_j=0;na_a=0;na_b=0;tab_At=0;wA_tag=NULL;dtmp_etAn=0;output_tmp=0;
#pragma omp for schedule(static)
	for (mx_j=0;mx_j<M_St->rpop_j*M_At->rpop_j;mx_j++){
	  ns_j = mx_j/M_At->rpop_j; na_j = mx_j % M_At->rpop_j;
	  ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	  switch (output_spacing_Tt){ case SPACING_j: tab_Tt=ns_j; break; case SPACING_b: tab_Tt=ns_b; break; case SPACING_a: tab_Tt=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
	  dtmp_etSn=0;
	  wS_tag = (__m128i*)((unsigned long *)(&(M_St->wX[ns_b*M_St->mc_length])));
	  mc_tag = (__m128i*)((unsigned long *)(&(M_St->mc_j[0])));
	  mc_end = (__m128i*)((unsigned long *)(&(M_St->mc_j[M_St->mc_length])));
	  dtmp_etSn = -M_St->cpop_j + 2*popcount(&wS_tag,&mc_tag,&mc_end); 
	  dtmp_etTn=0;
	  wT_tag = (__m128i*)((unsigned long *)(&(M_Tt->wX[ns_b*M_Tt->mc_length])));
	  mc_tag = (__m128i*)((unsigned long *)(&(M_Tt->mc_j[0])));
	  mc_end = (__m128i*)((unsigned long *)(&(M_Tt->mc_j[M_Tt->mc_length])));
	  dtmp_etTn = -M_Tt->cpop_j + 2*popcount(&wT_tag,&mc_tag,&mc_end); 
	  na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
	  switch (output_spacing_At){ case SPACING_j: tab_At=na_j; break; case SPACING_b: tab_At=na_b; break; case SPACING_a: tab_At=na_a; break; default: break; /* switch (output_spacing_At){ } */}
	  output_tmp = (*L2_get(lf_AtTn,na_j,na_b,na_a,ns_j,ns_b,ns_a) - dtmp_etTn*a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*(*L2_get(lf_ZtSn,na_j,na_b,na_a,ns_j,ns_b,ns_a) - a_An[na_a/POPLENGTH]*dtmp_etSn);
	  wA_tag = (__m128i *)(&(M_At->wX[na_b*M_At->mc_length]));
	  mc_tag = (__m128i *)(&(M_At->mc_j[0]));
	  mc_end = (__m128i *)(&(M_At->mc_j[M_At->mc_length]));
	  dtmp_etAn = -M_At->cpop_j + 2*popcount(&wA_tag,&mc_tag,&mc_end);
	  output_tmp -= D_An[na_a/POPLENGTH]*(M_At->cpop_j - 2*dtmp_etAn*a_An[na_a/POPLENGTH] + M_At->cpop_j*a_An[na_a/POPLENGTH]*a_An[na_a/POPLENGTH]);
	  output_D_AtTn_AtTn_vv->lf[tab_At + tab_Tt*tab_At_stride] = output_tmp;
	  /* for (mx_j=0;mx_j<M_St->rpop_j*M_At->rpop_j;mx_j++){ } */}
	/* end omp parallel */}
      GLOBAL_ops_count_one(tidx,M_At->rpop_j*M_Tt->rpop_j,M_At->rpop_j*M_Tt->rpop_j*M_At->mc_length*BIT8);
      /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
    /* if (structures exist){ } */}
  else{ printf(" %% Warning! structures lf_AtTn and lf_ZtSn does not exist in call to get_D_AtTn_AtTn_vv\n");}
  if (verbose>2){ lfprintf(output_D_AtTn_AtTn_vv," %% output_D_AtTn_AtTn_vv: ");}
 skip_D_AtTn_AtTn_vv:
  if (verbose>1){ printf(" %% [finished get_D_AtTn_AtTn_vv] tidx %d\n",tidx);}
  return NULL;
}

int wrap_D_AtTn_AtTn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle *lf_AtTn,struct L_handle *lf_ZtSn,struct L_handle **output_p)
 {
   /* This function calls get_D_AtTn_AtTn_vv ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 11)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_D_AtTn_AtTn_vv__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Zt,verbose," %% M_Zt: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_D_AtTn_AtTn_vv__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Zt; vpra[ip++] = M_St; vpra[ip++] = A_ajdk; vpra[ip++] = lf_AtTn; vpra[ip++] = lf_ZtSn; vpra[ip++] = *output_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_D_AtTn_AtTn_vv,vpra)){ printf("Warning! cannot create thread %d in wrap_D_AtTn_AtTn_vv\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_D_AtTn_AtTn_vv(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_D_AtTn_AtTn_vv__run] tidx %d\n",*tidx);}
  return length;
}  

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void *get_D_AtTn_ZtSn_uu(void *vp)
{
  /* This function takes in M_At, M_St, M_Zt, lf_AtTn, lf_ZtSn, and calculates output[nk+ns*nrows_An] = (AtTn)(nk,ns).*D(nk).*(ZtSn)(nk,ns) ;
     The output should look like: ;
     ((At(:,:)-a_An*e_At)*T(:,ns)).*diag(D_An(:)).*((Zt(:,:)-a_An*e_Zt)*S(:,ns)) ;
     = ;
     + (At(:,:)*T(:,ns)).*diag(D_An(:)).*(Zt(:,:)*S(:,ns)) ;
     - (a_An*e_At*T(:,ns)).*diag(D_An(:)).*(Zt(:,:)*S(:,ns)) ;
     - (At(:,:)*T(:,ns)).*diag(D_An(:)).*(a_An*e_Zt*S(:,ns)) ;
     + (a_An*e_At*T(:,ns)).*diag(D_An(:)).*(a_An*e_Zt*S(:,ns)) ;
     In the situation where M_An==M_Zn, then we can also correct for repeated indices within An:
     D_AtTn_AtTn = ((At(:,:)-a_An*e_At)*T(:,ns)).*diag(D_An(:)).*((Zt(:,:)-a_An*e_Zt)*S(:,ns)) ;
     - D_An(nj)*( At->cpop_j - 2*popcount(A(:,nj))*a_An(nj) + At->cpop_j*a_An(j)^2) ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_D_AtTn_ZtSn_uu = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); int output_spacing_At = output_spacing_a; 
  int output_spacing_s = *(int *)(vpra[ip++]); int output_spacing_Tt = output_spacing_s;
  int ncols_A_p = psize(M_At->nrows) ;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*ncols_A_p]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*ncols_A_p]);
  int ns_j=0,ns_b=0,ns_a=0,tab_Tt_stride=0,tab_Tt=0,na_j=0,na_b=0,na_a=0,tab_At_stride=0,tab_At=0,nz_j=0,nz_b=0,nz_a=0,ma_j=0,ma_b=0,ma_a=0,mz_j=0,mz_b=0,mz_a=0;
  double output_tmp=0; 
  unsigned char *At_tag; unsigned char *Zt_tag; 
  unsigned char *Tt_tag; unsigned char *St_tag;
  int vA=0,vT=0,vZ=0,vS=0;
  if (verbose>1){ printf(" %% [entering get_D_AtTn_ZtSn_uu] tidx %d \n",tidx);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  switch (output_spacing_Tt){ case SPACING_j: tab_Tt_stride = M_Tt->rpop_j; break; case SPACING_b: tab_Tt_stride = M_Tt->rpop_b; break; case SPACING_a: tab_Tt_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_Tt){ } */}
  switch (output_spacing_At){ case SPACING_j: tab_At_stride = M_At->rpop_j; break; case SPACING_b: tab_At_stride = M_At->rpop_b; break; case SPACING_a: tab_At_stride = M_At->nrows; break; default: break; /* switch (output_spacing_At){ } */}
  output_D_AtTn_ZtSn_uu->spacing_row = output_spacing_At; output_D_AtTn_ZtSn_uu->row_stride = tab_At_stride;
  output_D_AtTn_ZtSn_uu->spacing_col = output_spacing_Tt; output_D_AtTn_ZtSn_uu->col_stride = tab_Tt_stride;
  fill_uchar_zero((unsigned char *)(output_D_AtTn_ZtSn_uu->lf),tab_At_stride*tab_Tt_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"D_AtTn_ZtSn_uu")){ goto skip_D_AtTn_ZtSn_uu;}
  ns_j=0; while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_Tt){ case SPACING_j: tab_Tt=ns_j; break; case SPACING_b: tab_Tt=ns_b; break; case SPACING_a: tab_Tt=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
    St_tag = (unsigned char*)(&(M_St->wX[ns_b*M_St->mc_length]));
    Tt_tag = (unsigned char*)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
    na_j=0; nz_j=0; while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
      na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
      nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
      At_tag = (unsigned char*)(&(M_At->wX[na_b*M_At->mc_length]));
      Zt_tag = (unsigned char*)(&(M_Zt->wX[nz_b*M_Zt->mc_length]));
      switch (output_spacing_At){ case SPACING_j: tab_At=na_j; break; case SPACING_b: tab_At=na_b; break; case SPACING_a: tab_At=na_a; break; default: break; /* switch (output_spacing_At){ } */}
      output_tmp=0;
      ma_j=0; while (ma_j<M_At->cpop_j){
	ma_a = M_At->n_a_[ma_j]; ma_b = M_At->n_b_[ma_j];
	vA = bget____(At_tag,ma_a);
	vT = bget____(Tt_tag,ma_a);
	mz_j=0; while (mz_j<M_Zt->cpop_j){
	  mz_a = M_Zt->n_a_[mz_j]; mz_b = M_Zt->n_b_[mz_j];
	  vZ = bget____(Zt_tag,mz_a);
	  vS = bget____(St_tag,mz_a);
	  output_tmp += vT*(vA - a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*(vZ-a_An[na_a/POPLENGTH])*vS;
	  mz_j++; /* while (mz_j<M_Zt->cpop_j){ } */}
	ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
      output_D_AtTn_ZtSn_uu->lf[tab_At + tab_Tt*tab_At_stride] = output_tmp;
      na_j++; nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_At->rpop_j*M_At->cpop_j*M_Zt->cpop_j*3,0);
  if (verbose>2){ lfprintf(output_D_AtTn_ZtSn_uu," %% output_D_AtTn_ZtSn_uu: ");}
 skip_D_AtTn_ZtSn_uu:
  if (verbose>1){ printf(" %% [finished get_D_AtTn_ZtSn_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_D_AtTn_ZtSn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle **output_p)
 {
   /* This function calls get_D_AtTn_ZtSn_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 9)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_D_AtTn_ZtSn_uu__run] tidx %d\n",*tidx);}
  if (verbose>2){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose>2){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose>2){ M_handle_printf(M_Zt,verbose," %% M_Zt: ");}
  if (verbose>2){ M_handle_printf(M_St,verbose," %% M_St: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_D_AtTn_ZtSn_uu__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Zt; vpra[ip++] = M_St; vpra[ip++] = A_ajdk; vpra[ip++] = *output_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_D_AtTn_ZtSn_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_D_AtTn_ZtSn_uu\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_D_AtTn_ZtSn_uu(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_D_AtTn_ZtSn_uu__run] tidx %d\n",*tidx);}
  return length;
}  

void *get_D_AtTn_AtTn_uu(void *vp)
{
  /* This function takes in M_At, M_St, M_Zt, lf_AtTn, lf_ZtSn, and calculates output[nk+ns*nrows_An] = (AtTn)(nk,ns).*D(nk).*(ZtSn)(nk,ns) ;
     The output should look like: ;
     ((At(:,:)-a_An*e_At)*T(:,ns)).*diag(D_An(:)).*((Zt(:,:)-a_An*e_Zt)*S(:,ns)) ;
     = ;
     + (At(:,:)*T(:,ns)).*diag(D_An(:)).*(Zt(:,:)*S(:,ns)) ;
     - (a_An*e_At*T(:,ns)).*diag(D_An(:)).*(Zt(:,:)*S(:,ns)) ;
     - (At(:,:)*T(:,ns)).*diag(D_An(:)).*(a_An*e_Zt*S(:,ns)) ;
     + (a_An*e_At*T(:,ns)).*diag(D_An(:)).*(a_An*e_Zt*S(:,ns)) ;
     In the situation where M_An==M_Zn, then we can also correct for repeated indices within An:
     D_AtTn_AtTn = ((At(:,:)-a_An*e_At)*T(:,ns)).*diag(D_An(:)).*((Zt(:,:)-a_An*e_Zt)*S(:,ns)) ;
     - D_An(nj)*( At->cpop_j - 2*popcount(A(:,nj))*a_An(nj) + At->cpop_j*a_An(j)^2) ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_D_AtTn_AtTn_uu = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); int output_spacing_At = output_spacing_a; 
  int output_spacing_s = *(int *)(vpra[ip++]); int output_spacing_Tt = output_spacing_s;
  int ncols_A_p = psize(M_At->nrows) ;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*ncols_A_p]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*ncols_A_p]);
  int ns_j=0,ns_b=0,ns_a=0,tab_Tt_stride=0,tab_Tt=0,na_j=0,na_b=0,na_a=0,tab_At_stride=0,tab_At=0,nz_j=0,nz_b=0,nz_a=0,ma_j=0,ma_b=0,ma_a=0,mz_j=0,mz_b=0,mz_a=0;
  double output_tmp=0; 
  unsigned char *At_tag; unsigned char *Zt_tag; 
  unsigned char *Tt_tag; unsigned char *St_tag;
  int vA=0,vT=0,vZ=0,vS=0;
  if (verbose>1){ printf(" %% [entering get_D_AtTn_AtTn_uu] tidx %d \n",tidx);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  switch (output_spacing_Tt){ case SPACING_j: tab_Tt_stride = M_Tt->rpop_j; break; case SPACING_b: tab_Tt_stride = M_Tt->rpop_b; break; case SPACING_a: tab_Tt_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_Tt){ } */}
  switch (output_spacing_At){ case SPACING_j: tab_At_stride = M_At->rpop_j; break; case SPACING_b: tab_At_stride = M_At->rpop_b; break; case SPACING_a: tab_At_stride = M_At->nrows; break; default: break; /* switch (output_spacing_At){ } */}
  output_D_AtTn_AtTn_uu->spacing_row = output_spacing_At; output_D_AtTn_AtTn_uu->row_stride = tab_At_stride;
  output_D_AtTn_AtTn_uu->spacing_col = output_spacing_Tt; output_D_AtTn_AtTn_uu->col_stride = tab_Tt_stride;
  fill_uchar_zero((unsigned char *)(output_D_AtTn_AtTn_uu->lf),tab_At_stride*tab_Tt_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"D_AtTn_AtTn_uu")){ goto skip_D_AtTn_AtTn_uu;}
  ns_j=0; while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_Tt){ case SPACING_j: tab_Tt=ns_j; break; case SPACING_b: tab_Tt=ns_b; break; case SPACING_a: tab_Tt=ns_a; break; default: break; /* switch (output_spacing_Tt){ } */}
    St_tag = (unsigned char*)(&(M_St->wX[ns_b*M_St->mc_length]));
    Tt_tag = (unsigned char*)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
    na_j=0; nz_j=0; while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
      na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
      nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
      At_tag = (unsigned char*)(&(M_At->wX[na_b*M_At->mc_length]));
      Zt_tag = (unsigned char*)(&(M_Zt->wX[nz_b*M_Zt->mc_length]));
      switch (output_spacing_At){ case SPACING_j: tab_At=na_j; break; case SPACING_b: tab_At=na_b; break; case SPACING_a: tab_At=na_a; break; default: break; /* switch (output_spacing_At){ } */}
      output_tmp=0;
      ma_j=0; while (ma_j<M_At->cpop_j){
	ma_a = M_At->n_a_[ma_j]; ma_b = M_At->n_b_[ma_j];
	vA = bget____(At_tag,ma_a);
	vT = bget____(Tt_tag,ma_a);
	mz_j=0; while (mz_j<M_Zt->cpop_j){
	  mz_a = M_Zt->n_a_[mz_j]; mz_b = M_Zt->n_b_[mz_j];
	  vZ = bget____(Zt_tag,mz_a);
	  vS = bget____(St_tag,mz_a);
	  if (ma_j!=mz_j){ output_tmp += vT*(vA - a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*(vZ-a_An[na_a/POPLENGTH])*vS;}
	  mz_j++; /* while (mz_j<M_Zt->cpop_j){ } */}
	ma_j++; /* while (ma_j<M_At->cpop_j){ } */}
      output_D_AtTn_AtTn_uu->lf[tab_At + tab_Tt*tab_At_stride] = output_tmp;
      na_j++; nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_At->rpop_j*M_At->cpop_j*M_Zt->cpop_j*3,0);
  if (verbose>2){ lfprintf(output_D_AtTn_AtTn_uu," %% output_D_AtTn_AtTn_uu: ");}
 skip_D_AtTn_AtTn_uu:
  if (verbose>1){ printf(" %% [finished get_D_AtTn_AtTn_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_D_AtTn_AtTn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle **output_p)
 {
   /* This function calls get_D_AtTn_AtTn_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 9)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_D_AtTn_AtTn_uu__run] tidx %d\n",*tidx);}
  if (verbose>2){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose>2){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose>2){ M_handle_printf(M_Zt,verbose," %% M_Zt: ");}
  if (verbose>2){ M_handle_printf(M_St,verbose," %% M_St: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_D_AtTn_AtTn_uu__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Zt; vpra[ip++] = M_St; vpra[ip++] = A_ajdk; vpra[ip++] = *output_p; 
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_D_AtTn_AtTn_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_D_AtTn_AtTn_uu\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_D_AtTn_AtTn_uu(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_D_AtTn_AtTn_uu__run] tidx %d\n",*tidx);}
  return length;
}  

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void wrap_D_AtTn_ZtSn_vv_test()
{
  /* test for errors with input file: D_AtTn_ZtSn_vv_error.in ;
  */
  /* test for speed with input file: D_AtTn_ZtSn_vv_speed.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AtTn=NULL,**lf_YtTn=NULL,**lf_ZtSn=NULL,**lf_WtSn=NULL,**lf_AtTYn=NULL,**lf_ZtSWn=NULL;
  struct L_handle **lf_D_AtTn_ZtSn_vv=NULL,**lf_D_AtTn_AtTn_vv=NULL; int *length_D_AtTn_ZtSn_vv=NULL,*length_D_AtTn_AtTn_vv=NULL;
  struct L_handle **lf_D_AtTn_ZtSn_uu=NULL,**lf_D_AtTn_AtTn_uu=NULL; int *length_D_AtTn_ZtSn_uu=NULL,*length_D_AtTn_AtTn_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_D_AtTn_ZtSn_vv_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  lf_AtTn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_YtTn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_ZtSn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_WtSn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_AtTn[nb] = L_handle_make(M_An[nb]->ncols*M_Tn[nb]->ncols);
    lf_YtTn[nb] = L_handle_make(M_Yn[nb]->ncols*M_Tn[nb]->ncols);
    lf_ZtSn[nb] = L_handle_make(M_Zn[nb]->ncols*M_Sn[nb]->ncols);
    lf_WtSn[nb] = L_handle_make(M_Wn[nb]->ncols*M_Sn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_D_AtTn_ZtSn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_D_AtTn_ZtSn_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_D_AtTn_AtTn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_D_AtTn_AtTn_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_D_AtTn_ZtSn_vv[nb] = L_handle_make(M_An[nb]->ncols*M_Tn[nb]->ncols);
    lf_D_AtTn_AtTn_vv[nb] = L_handle_make(M_An[nb]->ncols*M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (error_check){ 
    lf_D_AtTn_ZtSn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_D_AtTn_ZtSn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    lf_D_AtTn_AtTn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_D_AtTn_AtTn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    for (nb=0;nb<nbins;nb++){ 
      lf_D_AtTn_ZtSn_uu[nb] = L_handle_make(M_An[nb]->ncols*M_Tn[nb]->ncols);
      lf_D_AtTn_AtTn_uu[nb] = L_handle_make(M_An[nb]->ncols*M_Tn[nb]->ncols);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (error_check){ } */}
  for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
	if (verbose){ printf(" %% %s; %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_An[nb],A_ajdk,&(lf_An_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Zn[nb],A_ajdk,&(lf_Zn_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Yn[nb],Y_ajdk,&(lf_Yn_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Wn[nb],Y_ajdk,&(lf_Wn_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_ajdk_v: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],NULL,NULL,NULL,&(lf_AtTn[nb]));
	  GLOBAL_pthread_toc(); 
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Yt[nb],M_Tt[nb],NULL,NULL,NULL,&(lf_YtTn[nb]));
	  GLOBAL_pthread_toc(); 
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],NULL,NULL,NULL,&(lf_ZtSn[nb]));
	  GLOBAL_pthread_toc(); 
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Wt[nb],M_St[nb],NULL,NULL,NULL,&(lf_WtSn[nb]));
	  GLOBAL_pthread_toc(); 
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc();
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
	  length_D_AtTn_ZtSn_vv[nb] = wrap_D_AtTn_ZtSn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_Zt[nb],M_St[nb],A_ajdk,lf_AtTn[nb],lf_ZtSn[nb],&(lf_D_AtTn_ZtSn_vv[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
	  length_D_AtTn_AtTn_vv[nb] = wrap_D_AtTn_AtTn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_At[nb],M_Tt[nb],A_ajdk,lf_AtTn[nb],lf_AtTn[nb],&(lf_D_AtTn_AtTn_vv[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% D_AtTn_ZtSn_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic(); 
	    length_D_AtTn_ZtSn_uu[nb] = wrap_D_AtTn_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_Zt[nb],M_St[nb],A_ajdk,&(lf_D_AtTn_ZtSn_uu[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic(); 
	    length_D_AtTn_AtTn_uu[nb] = wrap_D_AtTn_AtTn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_At[nb],M_Tt[nb],A_ajdk,&(lf_D_AtTn_AtTn_uu[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc();
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYn_uu and ZtSWn_uu: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_D_AtTn_ZtSn_vv[%d] - lf_D_AtTn_ZtSn_uu[%d] error %0.16f\n",nb,nb,dra_diff(lf_D_AtTn_ZtSn_vv[nb]->lf,lf_D_AtTn_ZtSn_uu[nb]->lf,length_D_AtTn_ZtSn_vv[nb],1));
	    printf(" %% lf_D_AtTn_AtTn_vv[%d] - lf_D_AtTn_AtTn_uu[%d] error %0.16f\n",nb,nb,dra_diff(lf_D_AtTn_AtTn_vv[nb]->lf,lf_D_AtTn_AtTn_uu[nb]->lf,length_D_AtTn_AtTn_vv[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  /* if (error_check){ } */}
	/* for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_An_ZtSWn_Yt_vv_test]\n");}
}
