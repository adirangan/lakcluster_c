void *get_An_ZtSWn_Yt_vv(void *vp)
{
  /* This function takes in M_An, M_Yn, lf_ZtSWn and calculates the output_An_ZtSWn_Yt: ;
     output_An_ZtSWn_Yt[mj+ns*A_n_rowsn] = (An(mj,:) - e_An*a_At)*D*ZtSWn*D*(Yt(:,mj) - a_Yt*e_Yn). ;
     In this scenario we expect ZtSWn to be: ;
     ZtSWn = (Zt(:,:)-a_An*e_At)*diag(S)*(Wn(:,:)-e_Yn*a_Yt). ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Yn = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *lf_ZtSWn = (struct L_handle *)(vpra[ip++]); /* used for direct calculation */
  struct L_handle *output_An_ZtSWn_Yt = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); 
  int output_spacing_s = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols)/* rup(M_An->ncols + M_An->ncols_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yn->ncols)/* rup(M_Yn->ncols + M_Yn->ncols_extend,POPLENGTH)/POPLENGTH */;
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int mx_j=0,mx_chunk=0,ns_j=0,ns_b=0,ns_a=0,tab_s_stride=0,tab_s=0,na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0,tab_a_stride=0,tab_a=0,ny_j=0,ny_b=0,ny_a=0,my_j=0,my_b=0,my_a=0;
  int vA=0,vY=0;
  double output_An_ZtSWn_Yt_tmp=0; 
  unsigned int *ma_b_,*ma_a_;
  unsigned int *na_b_,*na_a_;
  unsigned int *my_b_,*my_a_;
  unsigned int *ny_b_,*ny_a_;
  unsigned int *ns_b_,*ns_a_;
  int ms_j=0,ms_b=0,ms_a=0,mz_j=0,mz_b=0,mz_a=0,mw_j=0,mw_b=0,mw_a=0;
  __m128i *wAn_tag=NULL;
  __m128i *wYn_tag=NULL;
  if (verbose>1){ printf(" %% [entering get_An_ZtSWn_Yt_vv] tidx %d \n",tidx);}
  if (verbose>2){ raprintf(D_An,"double",1,A_pcols," %% D_An: "); raprintf(a_An,"double",1,A_pcols," %% a_An: ");}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  if (verbose>2){ raprintf(D_Yn,"double",1,Y_pcols," %% D_Yn: "); raprintf(a_Yn,"double",1,Y_pcols," %% a_Yn: ");}
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  na_b_ = M_An->n_b_; na_a_ = M_An->n_a_;
  my_b_ = M_Yn->m_b_; my_a_ = M_Yn->m_a_;
  ny_b_ = M_Yn->n_b_; ny_a_ = M_Yn->n_a_;
  ns_b_ = M_St->m_b_; ns_a_ = M_St->m_a_;
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  output_An_ZtSWn_Yt->spacing_row = output_spacing_a; output_An_ZtSWn_Yt->row_stride = tab_a_stride;
  output_An_ZtSWn_Yt->spacing_col = output_spacing_s; output_An_ZtSWn_Yt->col_stride = tab_s_stride;
  if (strstr(GLOBAL_skip,"An_ZtSWn_Yt_vv")){ goto skip_An_ZtSWn_Yt_vv;}
  if (lf_ZtSWn!=NULL){
    if (GLOBAL_omp_type==GLOBAL_omp_off){
      ns_j=0;
      while (ns_j<M_St->rpop_j){
	ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
	switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
	ma_j=0;my_j=0;
	while (ma_j<M_An->rpop_j && my_j<M_Yn->rpop_j){
	  ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j]; my_a = my_a_[my_j]; my_b = my_b_[my_j];
	  switch (output_spacing_a){ case SPACING_j: tab_a=ma_j+tab_s*tab_a_stride; break; case SPACING_b: tab_a=ma_b+tab_s*tab_a_stride; break; case SPACING_a: tab_a=ma_a+tab_s*tab_a_stride; break; default: break; /* switch (output_spacing_a){ } */}
	  wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]); wYn_tag = (__m128i*)&(M_Yn->wX[my_b*M_Yn->mc_length]);
	  /* na_j=0;  */
	  /* while (na_j<M_An->cpop_j){ */
	  /*   na_a = na_a_[na_j]; na_b = na_b_[na_j]; */
	  /*   vA = ((int)(((((unsigned char *)wAn_tag)[na_a/BIT8] >> (7-(na_a%BIT8))) & 1) << 1)-(int)1); /\* vA = bget____(wAn_tag,na_a); *\/ */
	  /*   ny_j=0; output_An_ZtSWn_Yt_tmp=0; */
	  /*   while (ny_j<M_Yn->cpop_j){ */
	  /*     ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j]; */
	  /*     vY = ((int)(((((unsigned char *)wYn_tag)[ny_a/BIT8] >> (7-(ny_a%BIT8))) & 1) << 1)-(int)1); /\* vA = bget____(wAn_tag,ny_a); *\/ */
	  /*     output_An_ZtSWn_Yt_tmp += (*L3_get(lf_ZtSWn,na_j,na_b,na_a,ny_j,ny_b,ny_a,ns_j,ns_b,ns_a)) * D_Yn[ny_a/POPLENGTH]*(vY - a_Yn[ny_a/POPLENGTH]); */
	  /*     ny_j++; /\* while (ny_j<M_Yn->cpop_j){ } *\/} */
	  /*   output_An_ZtSWn_Yt_tmp *= (vA - a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]; */
	  /*   output_An_ZtSWn_Yt->lf[tab_a] += output_An_ZtSWn_Yt_tmp; */
	  /*   na_j++; /\* while (na_j<M_An->cpop_j){ } *\/} */
	  ny_j=0; while (ny_j<M_Yn->cpop_j){ ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j]; vY = ((int)(((((unsigned char *)wYn_tag)[ny_a/BIT8] >> (7-(ny_a%BIT8))) & 1) << 1)-(int)1); /* vA = bget____(wAn_tag,ny_a); */ 
	    output_An_ZtSWn_Yt_tmp=0;
	    na_j=0; while (na_j<M_An->cpop_j){ na_a = na_a_[na_j]; na_b = na_b_[na_j]; vA = ((int)(((((unsigned char *)wAn_tag)[na_a/BIT8] >> (7-(na_a%BIT8))) & 1) << 1)-(int)1); /* vA = bget____(wAn_tag,na_a); */
	      output_An_ZtSWn_Yt_tmp += (vA - a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH] * (*L3_get(lf_ZtSWn,na_j,na_b,na_a,ny_j,ny_b,ny_a,ns_j,ns_b,ns_a));
	      na_j++; /* while (na_j<M_An->cpop_j){ } */}
	    output_An_ZtSWn_Yt_tmp *= D_Yn[ny_a/POPLENGTH]*(vY - a_Yn[ny_a/POPLENGTH]);
	    output_An_ZtSWn_Yt->lf[tab_a] += output_An_ZtSWn_Yt_tmp;
	    ny_j++; /* while (ny_j<M_Yn->cpop_j){ } */}
	  ma_j++;my_j++;}
	ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
      GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_An->rpop_j*M_Yn->cpop_j*(M_An->cpop_j*2 + 1),0);
      /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
    else if (GLOBAL_omp_type==GLOBAL_omp__on){
      mx_chunk=128; 
#pragma omp parallel shared(mx_chunk) private(mx_j,ns_j,ns_a,ns_b,ma_j,ma_a,ma_b,my_j,my_a,my_b,tab_s,tab_a,wAn_tag,wYn_tag,output_An_ZtSWn_Yt_tmp,ny_j,ny_a,ny_b,vY,na_j,na_a,na_b,vA)
      { /* begin omp parallel */
	mx_j=0; tab_s=0; tab_a=0; 
#pragma omp for schedule(dynamic,mx_chunk)
	for (mx_j=0;mx_j<M_St->rpop_j*M_An->rpop_j;mx_j++){
	  ns_j = mx_j / M_An->rpop_j; ma_j = mx_j % M_An->rpop_j; my_j = ma_j;
	  ns_a = ns_a_[ns_j]; ns_b = ns_b_[ns_j];
	  switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
	  ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j]; my_a = my_a_[my_j]; my_b = my_b_[my_j];
	  switch (output_spacing_a){ case SPACING_j: tab_a=ma_j+tab_s*tab_a_stride; break; case SPACING_b: tab_a=ma_b+tab_s*tab_a_stride; break; case SPACING_a: tab_a=ma_a+tab_s*tab_a_stride; break; default: break; /* switch (output_spacing_a){ } */}
	  wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]); wYn_tag = (__m128i*)&(M_Yn->wX[my_b*M_Yn->mc_length]);
	  ny_j=0; while (ny_j<M_Yn->cpop_j){ ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j]; vY = ((int)(((((unsigned char *)wYn_tag)[ny_a/BIT8] >> (7-(ny_a%BIT8))) & 1) << 1)-(int)1); /* vA = bget____(wAn_tag,ny_a); */ 
	    output_An_ZtSWn_Yt_tmp=0;
	    na_j=0; while (na_j<M_An->cpop_j){ na_a = na_a_[na_j]; na_b = na_b_[na_j]; vA = ((int)(((((unsigned char *)wAn_tag)[na_a/BIT8] >> (7-(na_a%BIT8))) & 1) << 1)-(int)1); /* vA = bget____(wAn_tag,na_a); */
	      output_An_ZtSWn_Yt_tmp += (vA - a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH] * (*L3_get(lf_ZtSWn,na_j,na_b,na_a,ny_j,ny_b,ny_a,ns_j,ns_b,ns_a));
	      na_j++; /* while (na_j<M_An->cpop_j){ } */}
	    output_An_ZtSWn_Yt_tmp *= D_Yn[ny_a/POPLENGTH]*(vY - a_Yn[ny_a/POPLENGTH]);
#pragma omp critical
	    output_An_ZtSWn_Yt->lf[tab_a] += output_An_ZtSWn_Yt_tmp;
	    ny_j++; /* while (ny_j<M_Yn->cpop_j){ } */}
	  /* for (mx_j=0;mx_j<M_St->rpop_j*M_An->rpop_j;mx_j++){ } */}
	/* end omp parallel */}
      GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_An->rpop_j*M_Yn->cpop_j*(M_An->cpop_j*2 + 1),0);
      /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
    /* if (lf_ZtSWn!=NULL){ } */}
  else /* if (lf_ZtSWn==NULL){ } */{
    printf(" %% Warning! lf_ZtSWn==NULL in get_An_ZtSWn_Yt_vv\n");
    /* if (lf_ZtSWn==NULL){ } */}
 skip_An_ZtSWn_Yt_vv:
  if (verbose>1){ printf(" %% [finished get_An_ZtSWn_Yt_vv] tidx %d\n",tidx);}
  if (verbose>2){ raprintf(output_An_ZtSWn_Yt->lf,"double",tab_a_stride,tab_s_stride," %% output_An_ZtSWn_Yt->lf: ");}
  return NULL;
}

int wrap_An_ZtSWn_Yt_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_St,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle *lf_ZtSWn,struct L_handle **output_p)
{
  /* This function calls get_An_ZtSWn_Yt_vv ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 12)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_An_ZtSWn_Yt_vv__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Yn,verbose," %% M_Yn: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  length = length_a*length_s; if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_An_ZtSWn_Yt_vv__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_St; vpra[ip++] = M_Yn; vpra[ip++] = A_ajdk; vpra[ip++] = Y_ajdk; vpra[ip++] = lf_ZtSWn; vpra[ip++] = *output_p;
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_ZtSWn_Yt_vv,vpra)){ printf("Warning! cannot create thread %d in wrap_An_ZtSWn_Yt_vv__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_ZtSWn_Yt_vv(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_An_ZtSWn_Yt_vv__run] tidx %d\n",*tidx);}
  return length;
}  

void *get_AnZtSWnYt_uu(void *vp)
{
  /* This function takes in M_An, M_Zt, M_St, M_Wt, M_Yn and calculates the output_AnZtSWnYt: ;
     output_AnZtSWnYt[mj+ns*A_n_rowsn] = (An(mj,:) - e_An*a_At)*D*ZtSWn*D*(Yt(:,mj) - a_Yt*e_Yn). ;
     In this scenario we expect ZtSWn to be: ;
     ZtSWn = (Zt(:,:)-a_An*e_At)*diag(S)*(Wn(:,:)-e_Yn*a_Yt). ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Wt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Yn = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_AnZtSWnYt = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]); 
  int output_spacing_s = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols)/* rup(M_An->ncols + M_An->ncols_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yn->ncols)/* rup(M_Yn->ncols + M_Yn->ncols_extend,POPLENGTH)/POPLENGTH */;
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int mx_j=0,mx_chunk=0,ns_j=0,ns_b=0,ns_a=0,tab_s_stride=0,tab_s=0,na_j=0,na_b=0,na_a=0,ma_j=0,ma_b=0,ma_a=0,tab_a_stride=0,tab_a=0,ny_j=0,ny_b=0,ny_a=0,my_j=0,my_b=0,my_a=0;
  int vA=0,vY=0;
  double output_AnZtSWnYt_tmp=0; 
  unsigned int *ma_b_,*ma_a_;
  unsigned int *na_b_,*na_a_;
  unsigned int *my_b_,*my_a_;
  unsigned int *ny_b_,*ny_a_;
  unsigned int *ns_b_,*ns_a_;
  int ms_j=0,ms_b=0,ms_a=0,mz_j=0,mz_b=0,mz_a=0,mw_j=0,mw_b=0,mw_a=0;
  int ncols_Z_p = psize(M_Zt->nrows);
  double *D_Zn = (double *)&(A_ajdk[0 + AJDK_0_1*ncols_Z_p]);
  double *a_Zn = (double *)&(A_ajdk[0 + AJDK_1_0*ncols_Z_p]);
  int ncols_W_p = psize(M_Wt->nrows);
  double *D_Wn = (double *)&(Y_ajdk[0 + AJDK_0_1*ncols_W_p]);
  double *a_Wn = (double *)&(Y_ajdk[0 + AJDK_1_0*ncols_W_p]);
  int vZ=0,vS=0,vW=0;
  unsigned char *wAn_tag=NULL,*wZt_tag=NULL,*wSt_tag=NULL,*wWt_tag=NULL,*wYn_tag=NULL;
  if (verbose>1){ printf(" %% [entering get_AnZtSWnYt_uu] tidx %d \n",tidx);}
  if (verbose>2){ raprintf(D_An,"double",1,A_pcols," %% D_An: "); raprintf(a_An,"double",1,A_pcols," %% a_An: ");}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  if (verbose>2){ raprintf(D_Yn,"double",1,Y_pcols," %% D_Yn: "); raprintf(a_Yn,"double",1,Y_pcols," %% a_Yn: ");}
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  na_b_ = M_An->n_b_; na_a_ = M_An->n_a_;
  my_b_ = M_Yn->m_b_; my_a_ = M_Yn->m_a_;
  ny_b_ = M_Yn->n_b_; ny_a_ = M_Yn->n_a_;
  ns_b_ = M_St->m_b_; ns_a_ = M_St->m_a_;
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  output_AnZtSWnYt->spacing_row = output_spacing_a; output_AnZtSWnYt->row_stride = tab_a_stride;
  output_AnZtSWnYt->spacing_col = output_spacing_s; output_AnZtSWnYt->col_stride = tab_s_stride;
  fill_uchar_zero((unsigned char *)(output_AnZtSWnYt->lf),tab_a_stride*tab_s_stride*sizeof(double));
  ns_j=0;
  while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    wSt_tag = (unsigned char *)(&(M_St->wX[ns_b*M_St->mc_length]));
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    ma_j=0;my_j=0;
    while (ma_j<M_An->rpop_j && my_j<M_Yn->rpop_j){
      ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j]; my_a = my_a_[my_j]; my_b = my_b_[my_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=ma_j+tab_s*tab_a_stride; break; case SPACING_b: tab_a=ma_b+tab_s*tab_a_stride; break; case SPACING_a: tab_a=ma_a+tab_s*tab_a_stride; break; default: break; /* switch (output_spacing_a){ } */}
      wAn_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length])); wYn_tag = (unsigned char *)(&(M_Yn->wX[my_b*M_Yn->mc_length]));
      na_j=0; 
      while (na_j<M_An->cpop_j){
	na_a = na_a_[na_j]; na_b = na_b_[na_j];
	vA = bget____(wAn_tag,na_a);
	wZt_tag = (unsigned char *)(&(M_Zt->wX[na_b*M_Zt->mc_length]));
	ny_j=0; 
	while (ny_j<M_Yn->cpop_j){
	  ny_a = ny_a_[ny_j]; ny_b = ny_b_[ny_j];
	  vY = bget____(wYn_tag,ny_a);
	  wWt_tag = (unsigned char *)(&(M_Wt->wX[ny_b*M_Wt->mc_length]));
	  output_AnZtSWnYt_tmp=0;
	  mz_j=0; ms_j=0; mw_j=0;
	  while (mz_j<M_Zt->cpop_j && ms_j<M_St->cpop_j && mw_j<M_Wt->cpop_j){
	    mz_a = M_Zt->n_a_[mz_j]; mz_b = M_Zt->n_b_[mz_j];
	    ms_a = M_St->n_a_[ms_j]; ms_b = M_St->n_b_[ms_j];
	    mw_a = M_Wt->n_a_[mw_j]; mw_b = M_Wt->n_b_[mw_j];
	    vZ = bget____(wZt_tag,ms_a); vS = bget____(wSt_tag,ms_a); vW = bget____(wWt_tag,ms_a);
	    output_AnZtSWnYt_tmp += (vZ - a_Zn[na_a/POPLENGTH])*vS*(vW - a_Wn[ny_a/POPLENGTH]);
	    mz_j++; ms_j++; mw_j++; /* while (mz_j<M_Zt->cpop_j){ } */}
	  output_AnZtSWnYt->lf[tab_a] += (vA - a_An[na_a/POPLENGTH])*D_An[na_a/POPLENGTH]*output_AnZtSWnYt_tmp*D_Yn[ny_a/POPLENGTH]*(vY - a_Yn[ny_a/POPLENGTH]);
	  ny_j++; /* while (ny_j<M_Yn->cpop_j){ } */}
	na_j++; /* while (na_j<M_An->cpop_j){ } */}
      ma_j++;my_j++;}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_An->rpop_j*M_An->cpop_j*M_Yn->cpop_j*M_Zt->cpop_j*2,0);
  if (verbose>1){ printf(" %% [finished get_AnZtSWnYt_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_AnZtSWnYt_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_An,struct M_handle *M_Zt,struct M_handle *M_St,struct M_handle *M_Wt,struct M_handle *M_Yn,double *A_ajdk,double *Y_ajdk,struct L_handle **output_p)
{
  /* This function calls get_AnZtSWnYt_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 12)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_AnZtSWnYt_uu__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Yn,verbose," %% M_Yn: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %d*%d=%d\n",length_a,length_s,length);}
  length = length_a*length_s; if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: "); bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: "); bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_AnZtSWnYt_uu__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_Zt; vpra[ip++] = M_St; vpra[ip++] = M_Wt; vpra[ip++] = M_Yn; vpra[ip++] = A_ajdk; vpra[ip++] = Y_ajdk; vpra[ip++] = *output_p;
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AnZtSWnYt_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_AnZtSWnYt_uu__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AnZtSWnYt_uu(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p)->lf,"double",length_a,length_s," %% (*output_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_AnZtSWnYt_uu__run] tidx %d\n",*tidx);}
  return length;
}  

void wrap_An_ZtSWn_Yt_vv_test()
{
  /* test for errors with input file: An_ZtSWn_Yt_vv_error.in ;
  */
  /* test for speed with input file: An_ZtSWn_Yt_vv_speed.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AtTn=NULL,**lf_YtTn=NULL,**lf_ZtSn=NULL,**lf_WtSn=NULL,**lf_AtTYn=NULL,**lf_ZtSWn=NULL;
  struct L_handle **lf_An_ZtSWn_Yt_vv=NULL,**lf_An_AtTYn_Yt_vv=NULL; int *length_An_ZtSWn_Yt_vv=NULL,*length_An_AtTYn_Yt_vv=NULL;
  struct L_handle **lf_AnZtSWnYt_uu=NULL,**lf_AnAtTYnYt_uu=NULL; int *length_AnZtSWnYt_uu=NULL,*length_AnAtTYnYt_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_An_ZtSWn_Yt_vv_test]\n");}
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
  lf_ZtSWn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AtTYn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_ZtSWn[nb] = L_handle_make(M_Zn[nb]->ncols*M_Wn[nb]->ncols*M_Sn[nb]->ncols);
    lf_AtTYn[nb] = L_handle_make(M_An[nb]->ncols*M_Yn[nb]->ncols*M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_An_ZtSWn_Yt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_ZtSWn_Yt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_An_AtTYn_Yt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_AtTYn_Yt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_An_ZtSWn_Yt_vv[nb] = L_handle_make(M_An[nb]->nrows*M_Sn[nb]->ncols);
    lf_An_AtTYn_Yt_vv[nb] = L_handle_make(M_An[nb]->nrows*M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (error_check){ 
    lf_AnZtSWnYt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnZtSWnYt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    lf_AnAtTYnYt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnAtTYnYt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    for (nb=0;nb<nbins;nb++){ 
      lf_AnZtSWnYt_uu[nb] = L_handle_make(M_An[nb]->nrows*M_Sn[nb]->ncols);
      lf_AnAtTYnYt_uu[nb] = L_handle_make(M_An[nb]->nrows*M_Tn[nb]->ncols);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (error_check){ } */}
  for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
	if (verbose){ printf(" %% %s; %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
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
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic();
	  wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],M_Wt[nb],A_ajdk,Y_ajdk,lf_ZtSn[nb],lf_WtSn[nb],&(lf_ZtSWn[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic(); 
	  wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],M_Yt[nb],A_ajdk,Y_ajdk,lf_AtTn[nb],lf_YtTn[nb],&(lf_AtTYn[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYn_vv and ZtSWn_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic();
	  length_An_ZtSWn_Yt_vv[nb] = wrap_An_ZtSWn_Yt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_St[nb],M_Yn[nb],A_ajdk,Y_ajdk,lf_ZtSWn[nb],&(lf_An_ZtSWn_Yt_vv[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  length_An_AtTYn_Yt_vv[nb] = wrap_An_ZtSWn_Yt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_Tt[nb],M_Yn[nb],A_ajdk,Y_ajdk,lf_AtTYn[nb],&(lf_An_AtTYn_Yt_vv[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_ZtSWn_Yt_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	if (error_check){
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic();
	    length_AnZtSWnYt_uu[nb] = wrap_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_Zt[nb],M_St[nb],M_Wt[nb],M_Yn[nb],A_ajdk,Y_ajdk,&(lf_AnZtSWnYt_uu[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_AnAtTYnYt_uu[nb] = wrap_AnZtSWnYt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],M_At[nb],M_Tt[nb],M_Yt[nb],M_Yn[nb],A_ajdk,Y_ajdk,&(lf_AnAtTYnYt_uu[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZtSWnYt_uu: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_An_ZtSWn_Yt_vv[%d] error %0.16f\n",nb,dra_diff(lf_An_ZtSWn_Yt_vv[nb]->lf,lf_AnZtSWnYt_uu[nb]->lf,length_An_ZtSWn_Yt_vv[nb],1));
	    printf(" %% lf_An_AtTYn_Yt_vv[%d] error %0.16f\n",nb,dra_diff(lf_An_AtTYn_Yt_vv[nb]->lf,lf_AnAtTYnYt_uu[nb]->lf,length_An_AtTYn_Yt_vv[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
	/* for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_An_ZtSWn_Yt_vv_test]\n");}
}
