void *get_AnZt_vv__run(void *vp)
{
  /* This function serves as a module which calculates the matrix-matrix product AnZt, where both An and Zt are treated as +1/-1 values ;
     If the terms A_ajdk = (alpha=(p-q)).^j.*(D=1/(4*p*q)).^k are given, as well as lf_An_ajdk and lf_Zn_ajdk (recording the row-vectors storing An*diag(alpha.^j.*D.^k) and Zn*diag(alpha.^j.*D.^k), respectively), ;
     then we calculate (An-e_An*a_At)*diag(D)*(Zt-a_An*e_Zt). ;
     Otherwise, we simply calculate An*Zt. ;
     The calculation is performed using vector-operations. ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zn = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *lf_An_ajdk = (struct L_handle *)(vpra[ip++]);
  struct L_handle *lf_Zn_ajdk = (struct L_handle *)(vpra[ip++]);
  struct L_handle *output_AnZt = (struct L_handle *)(vpra[ip++]);
  int type_flag = *(int *)(vpra[ip++]); 
  int output_spacing_r = *(int *)(vpra[ip++]);
  int output_spacing_c = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols)/* rup(M_An->ncols + M_An->ncols_extend,POPLENGTH)/POPLENGTH */;
  int type_flag_single = 0;
  int mx_j=0,mx_chunk=0,ma_a=0,ma_b=0,ma_j=0,mz_a=0,mz_b=0,mz_j=0,nc=0,tab_An_ajdk=0,tab_r=0,tab_r_stride=0,tab_Zn_ajdk=0,tab_c=0,tab_c_stride=0;
  double dtmp=0,dtmp_a0d1=0,dtmp_a2d1=0,output_tmp=0;
  __m128i *wA_tag;
  __m128i *wZ_tag;
  __m128i *mc_tag,*mc_end;
  double *dinp=NULL;
  char tempchar[FNAMESIZE];
  int zper = maximum(1,minimum(M_Zn->rpop_j,maximum(1,/*262144*/ GLOBAL_2_cache_size / maximum(1,M_Zn->mc_length) * 3/4 )));
  int mz_block=0,mz_block_max=0;
  mz_block_max = rup(M_Zn->rpop_j,zper)/zper;
  if (type_flag==TYPE_p0){ printf(" %% Warning! TYPE_p0 not implemented within get_AnZt_vv__run\n");}
  if      (type_flag==TYPE_pm){ type_flag_single = TYPE_p_;}
  else if (type_flag==TYPE_00){ type_flag_single = TYPE_0_;}
  if (verbose>1){ printf(" %% [entering get_AnZt_vv__run] tidx %d rpop_b %d(%d) cpop_b %d(%d) xpop_b %d(%d) type_flag %d(%d) (%d,%d)\n",tidx,(int)(M_An->rpop_b),(int)(M_An->rpop_j),(int)(M_An->cpop_b),(int)(M_An->cpop_j),(int)(M_Zn->rpop_b),(int)(M_Zn->rpop_j),type_flag,type_flag_single,(int)zper,(int)(mz_block_max));}  
  if (verbose>1){ printf(" %% M_Zn->rpop_j %d, M_Zn->mc_length %d<--%d, zper = (l2 %d)/%d*3/4 = %d, mz_block_max %d\n",M_Zn->rpop_j,M_Zn->mc_length,M_Zn->ncols,(int)GLOBAL_2_cache_size,(int)M_Zn->mc_length,zper,mz_block_max);}
  if (verbose>3){ bprintf((unsigned char *)M_An->wX,POPLENGTH,M_An->rpop_b,minimum(M_An->mc_length*BIT8,M_An->ncols)," %% M_An->wX: ");}
  if (verbose>3){ bprintf((unsigned char *)M_Zn->wX,POPLENGTH,M_Zn->rpop_b,minimum(M_Zn->mc_length*BIT8,M_Zn->ncols)," %% M_Zn->wX: ");}
  unsigned int *ma_b_,*ma_a_;
  unsigned int *mz_b_,*mz_a_;
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  mz_b_ = M_Zn->m_b_; mz_a_ = M_Zn->m_a_;
  if (verbose>1){ sprintf(tempchar," %% ma_b_[%d]: ",tidx); raprintf(ma_b_,"unsigned int",1,M_An->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% ma_a_[%d]: ",tidx); raprintf(ma_a_,"unsigned int",1,M_An->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% mz_b_[%d]: ",tidx); raprintf(mz_b_,"unsigned int",1,M_Zn->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% mz_a_[%d]: ",tidx); raprintf(mz_a_,"unsigned int",1,M_Zn->rpop_j,tempchar);}
  switch (output_spacing_r){ case SPACING_j: tab_r_stride = M_An->rpop_j; break; case SPACING_b: tab_r_stride = M_An->rpop_b; break; case SPACING_a: tab_r_stride = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: tab_c_stride = M_Zn->rpop_j; break; case SPACING_b: tab_c_stride = M_Zn->rpop_b; break; case SPACING_a: tab_c_stride = M_Zn->nrows; break; default: break; /* switch (output_spacing_c){ } */}
  if (verbose>1 && lf_An_ajdk){ raprintf(lf_An_ajdk->lf,"double",lf_An_ajdk->row_stride,AJDK_TOT," %% lf_An_ajdk->lf: ");}
  if (verbose>1 && lf_Zn_ajdk){ raprintf(lf_Zn_ajdk->lf,"double",lf_Zn_ajdk->row_stride,AJDK_TOT," %% lf_Zn_ajdk->lf: ");}
  if (verbose>1){ printf(" %% calculating output_AnZt\n");}
  output_AnZt->spacing_row = output_spacing_r; output_AnZt->row_stride = tab_r_stride;
  output_AnZt->spacing_col = output_spacing_c; output_AnZt->col_stride = tab_c_stride;
  fill_uchar_zero((unsigned char *)(output_AnZt->lf),tab_r_stride*tab_c_stride*sizeof(double));
  if (strstr(GLOBAL_skip,"AnZt_vv")){ goto skip_AnZt_vv;}
  dtmp_a0d1=0; for (nc=0;nc<A_pcols;nc++){ dtmp_a0d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_0_1*A_pcols] : 1) * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
  dtmp_a2d1=0; for (nc=0;nc<A_pcols;nc++){ dtmp_a2d1 += (double)(A_ajdk ? A_ajdk[nc + AJDK_2_1*A_pcols] : 0) * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8);}
  if (verbose>2){ printf(" %% dtmp_a0d1 %f dtmp_a2d1 %f\n",dtmp_a0d1,dtmp_a2d1);}
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    mz_block=0;
    while (mz_block<mz_block_max){
      ma_j=0;
      while (ma_j<M_An->rpop_j){
	ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	if (lf_An_ajdk){ switch (lf_An_ajdk->spacing_row){ case SPACING_j: tab_An_ajdk=ma_j; break; case SPACING_b: tab_An_ajdk=ma_b; break; case SPACING_a: tab_An_ajdk=ma_a; break; default: break; /* switch (lf_An_ajdk->spacing_row){ } */}}
	switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
	mz_j=mz_block*zper;
	while (mz_j<minimum(M_Zn->rpop_j,mz_block*zper+zper)){
	  mz_a = mz_a_[mz_j]; mz_b = mz_b_[mz_j];
	  if (lf_Zn_ajdk){ switch (lf_Zn_ajdk->spacing_row){ case SPACING_j: tab_Zn_ajdk=mz_j; break; case SPACING_b: tab_Zn_ajdk=mz_b; break; case SPACING_a: tab_Zn_ajdk=mz_a; break; default: break; /* switch (lf_Zn_ajdk->spacing_row){ } */}}
	  switch (output_spacing_c){ case SPACING_j: tab_c=mz_j; break; case SPACING_b: tab_c=mz_b; break; case SPACING_a: tab_c=mz_a; break; default: break; /* switch (output_spacing_c){ } */}
	  if      (type_flag==TYPE_pm){ output_tmp = dtmp_a0d1;}
	  else if (type_flag==TYPE_00){ output_tmp = 0;}
	  wA_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	  wZ_tag = (__m128i*)&(M_Zn->wX[mz_b*M_Zn->mc_length]);
	  mc_tag = (__m128i*)&(M_An->mc_j[0]); mc_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	  if (A_ajdk){
	    dinp = &(A_ajdk[0+AJDK_0_1*A_pcols]);
	    if      (type_flag==TYPE_pm){ dtmp = popcount_xor_lf(&wA_tag,&wZ_tag,&mc_tag,&mc_end,&dinp);}
	    else if (type_flag==TYPE_00){ dtmp = popcount_and_lf(&wA_tag,&wZ_tag,&mc_tag,&mc_end,&dinp);}
	    /* if (A_ajdk){ } */}
	  else /* if (A_ajdk==NULL){ } */{
	    if      (type_flag==TYPE_pm){ dtmp = popcount_xor(&wA_tag,&wZ_tag,&mc_tag,&mc_end);}
	    else if (type_flag==TYPE_00){ dtmp = popcount_and(&wA_tag,&wZ_tag,&mc_tag,&mc_end);}
	    /* if (A_ajdk==NULL){ } */}
	  if      (type_flag==TYPE_pm){ output_tmp -= (2*dtmp);}
	  else if (type_flag==TYPE_00){ output_tmp += (dtmp);}
	  if (lf_An_ajdk && lf_Zn_ajdk){ output_tmp -= lf_An_ajdk->lf[tab_An_ajdk + AJDK_1_1*lf_An_ajdk->row_stride] + lf_Zn_ajdk->lf[tab_Zn_ajdk + AJDK_1_1*lf_Zn_ajdk->row_stride] - dtmp_a2d1;}
	  output_AnZt->lf[tab_r + tab_c*tab_r_stride] = output_tmp;
	  mz_j++; /* while (mz_j<M_Zn->rpop_j){ } */}
	ma_j++;}
      mz_block++; /* while (mz_block<mz_block_max){ } */}
    GLOBAL_ops_count_one(tidx,M_An->rpop_j*M_Zn->rpop_j,M_An->rpop_j*M_Zn->rpop_j*M_An->mc_length*BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp__on){
#pragma omp parallel private(mz_block,ma_j,ma_a,ma_b,mz_j,mz_a,mz_b,tab_An_ajdk,tab_Zn_ajdk,tab_r,tab_c,wA_tag,wZ_tag,mc_tag,mc_end,dinp,dtmp,output_tmp)
    { /* begin omp parallel */
      mz_block=0; tab_An_ajdk=0; tab_Zn_ajdk=0; tab_r=0; tab_c=0; dtmp=0;
      /* #pragma omp for schedule(dynamic,mx_chunk) */
#pragma omp for schedule(dynamic)
      for (mz_block=0;mz_block<mz_block_max;mz_block++){
	ma_j=0;
	while (ma_j<M_An->rpop_j){
	  ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	  if (lf_An_ajdk){ switch (lf_An_ajdk->spacing_row){ case SPACING_j: tab_An_ajdk=ma_j; break; case SPACING_b: tab_An_ajdk=ma_b; break; case SPACING_a: tab_An_ajdk=ma_a; break; default: break; /* switch (lf_An_ajdk->spacing_row){ } */}}
	  switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
	  mz_j=mz_block*zper;
	  while (mz_j<minimum(M_Zn->rpop_j,mz_block*zper+zper)){
	    mz_a = mz_a_[mz_j]; mz_b = mz_b_[mz_j];
	    if (lf_Zn_ajdk){ switch (lf_Zn_ajdk->spacing_row){ case SPACING_j: tab_Zn_ajdk=mz_j; break; case SPACING_b: tab_Zn_ajdk=mz_b; break; case SPACING_a: tab_Zn_ajdk=mz_a; break; default: break; /* switch (lf_Zn_ajdk->spacing_row){ } */}}
	    switch (output_spacing_c){ case SPACING_j: tab_c=mz_j; break; case SPACING_b: tab_c=mz_b; break; case SPACING_a: tab_c=mz_a; break; default: break; /* switch (output_spacing_c){ } */}
	    if      (type_flag==TYPE_pm){ output_tmp = dtmp_a0d1;}
	    else if (type_flag==TYPE_00){ output_tmp = 0;}
	    wA_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	    wZ_tag = (__m128i*)&(M_Zn->wX[mz_b*M_Zn->mc_length]);
	    mc_tag = (__m128i*)&(M_An->mc_j[0]); mc_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	    if (A_ajdk){
	      dinp = &(A_ajdk[0+AJDK_0_1*A_pcols]);
	      if      (type_flag==TYPE_pm){ dtmp = popcount_xor_lf(&wA_tag,&wZ_tag,&mc_tag,&mc_end,&dinp);}
	      else if (type_flag==TYPE_00){ dtmp = popcount_and_lf(&wA_tag,&wZ_tag,&mc_tag,&mc_end,&dinp);}
	      /* if (A_ajdk){ } */}
	    else /* if (A_ajdk==NULL){ } */{
	      if      (type_flag==TYPE_pm){ dtmp = popcount_xor(&wA_tag,&wZ_tag,&mc_tag,&mc_end);}
	      else if (type_flag==TYPE_00){ dtmp = popcount_and(&wA_tag,&wZ_tag,&mc_tag,&mc_end);}
	      /* if (A_ajdk==NULL){ } */}
	    if      (type_flag==TYPE_pm){ output_tmp -= (2*dtmp);}
	    else if (type_flag==TYPE_00){ output_tmp += (dtmp);}
	    if (lf_An_ajdk && lf_Zn_ajdk){ output_tmp -= lf_An_ajdk->lf[tab_An_ajdk + AJDK_1_1*lf_An_ajdk->row_stride] + lf_Zn_ajdk->lf[tab_Zn_ajdk + AJDK_1_1*lf_Zn_ajdk->row_stride] - dtmp_a2d1;}
	    output_AnZt->lf[tab_r + tab_c*tab_r_stride] = output_tmp;
	    mz_j++; /* while (mz_j<M_Zn->rpop_j){ } */}
	  ma_j++;}
	/* for (mz_block=0;mz_block<mz_block_max;mz_block++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,M_An->rpop_j*M_Zn->rpop_j,M_An->rpop_j*M_Zn->rpop_j*M_An->mc_length*BIT8);
    /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp_unused){
    mx_chunk=1920;
/* #pragma omp parallel shared(mx_chunk) private(mx_j,ma_j,ma_a,ma_b,mz_j,mz_a,mz_b,tab_An_ajdk,tab_Zn_ajdk,tab_r,tab_c,wA_tag,wZ_tag,mc_tag,mc_end,dinp,dtmp) */
#pragma omp parallel private(mx_j,ma_j,ma_a,ma_b,mz_j,mz_a,mz_b,tab_An_ajdk,tab_Zn_ajdk,tab_r,tab_c,wA_tag,wZ_tag,mc_tag,mc_end,dinp,dtmp,output_tmp)
    { /* begin omp parallel */
      mx_j=0; tab_An_ajdk=0; tab_Zn_ajdk=0; tab_r=0; tab_c=0; dtmp=0;
/* #pragma omp for schedule(dynamic,mx_chunk) */
#pragma omp for schedule(dynamic)
      for (mx_j=0;mx_j<M_An->rpop_j*M_Zn->rpop_j;mx_j++){
	ma_j = mx_j/M_Zn->rpop_j; mz_j = mx_j%M_Zn->rpop_j;
	ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	if (lf_An_ajdk){ switch (lf_An_ajdk->spacing_row){ case SPACING_j: tab_An_ajdk=ma_j; break; case SPACING_b: tab_An_ajdk=ma_b; break; case SPACING_a: tab_An_ajdk=ma_a; break; default: break; /* switch (lf_An_ajdk->spacing_row){ } */}}
	switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
	mz_a = mz_a_[mz_j]; mz_b = mz_b_[mz_j];
	if (lf_Zn_ajdk){ switch (lf_Zn_ajdk->spacing_row){ case SPACING_j: tab_Zn_ajdk=mz_j; break; case SPACING_b: tab_Zn_ajdk=mz_b; break; case SPACING_a: tab_Zn_ajdk=mz_a; break; default: break; /* switch (lf_Zn_ajdk->spacing_row){ } */}}
	switch (output_spacing_c){ case SPACING_j: tab_c=mz_j; break; case SPACING_b: tab_c=mz_b; break; case SPACING_a: tab_c=mz_a; break; default: break; /* switch (output_spacing_c){ } */}
	if      (type_flag==TYPE_pm){ output_tmp = dtmp_a0d1;}
	else if (type_flag==TYPE_00){ output_tmp = 0;}
	wA_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	wZ_tag = (__m128i*)&(M_Zn->wX[mz_b*M_Zn->mc_length]);
	mc_tag = (__m128i*)&(M_An->mc_j[0]); mc_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	if (A_ajdk){
	  dinp = &(A_ajdk[0+AJDK_0_1*A_pcols]);
	  if      (type_flag==TYPE_pm){ dtmp = popcount_xor_lf(&wA_tag,&wZ_tag,&mc_tag,&mc_end,&dinp);}
	  else if (type_flag==TYPE_00){ dtmp = popcount_and_lf(&wA_tag,&wZ_tag,&mc_tag,&mc_end,&dinp);}
	  /* if (A_ajdk){ } */}
	else /* if (A_ajdk==NULL){ } */{
	  if      (type_flag==TYPE_pm){ dtmp = popcount_xor(&wA_tag,&wZ_tag,&mc_tag,&mc_end);}
	  else if (type_flag==TYPE_00){ dtmp = popcount_and(&wA_tag,&wZ_tag,&mc_tag,&mc_end);}
	  /* if (A_ajdk==NULL){ } */}
	if      (type_flag==TYPE_pm){ output_tmp -= (2*dtmp);}
	else if (type_flag==TYPE_00){ output_tmp += (dtmp);}
	if (lf_An_ajdk && lf_Zn_ajdk){ output_tmp -= lf_An_ajdk->lf[tab_An_ajdk + AJDK_1_1*lf_An_ajdk->row_stride] + lf_Zn_ajdk->lf[tab_Zn_ajdk + AJDK_1_1*lf_Zn_ajdk->row_stride] - dtmp_a2d1;}
	output_AnZt->lf[tab_r + tab_c*tab_r_stride] = output_tmp;
	/* for (mx_j=0;mx_j<M_An->rpop_j*M_Zn->rpop_j;mx_j++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,M_An->rpop_j*M_Zn->rpop_j,M_An->rpop_j*M_Zn->rpop_j*M_An->mc_length*BIT8);
    /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
 skip_AnZt_vv:
  if (verbose>1){ printf(" %% [finished get_AnZt_vv__run] tidx %d \n",tidx);}
  return NULL;
}

int wrap_AnZt_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,int output_spacing_c,struct M_handle *M_An,struct M_handle *M_Zn,double *A_ajdk,struct L_handle *lf_An_ajdk,struct L_handle *lf_Zn_ajdk,struct L_handle **output_AnZt_p)
{
  /* This function uses the M_handles M_A and M_Zn to run a series of parallel calls to get_AnZt_vv__run ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 14)
     the type of multiplication is determined by type_flag: ;
     type_flag==TYPE_pm: get_AX_pm ;
     type_flag==TYPE_00: get_AX_00 ;
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_r=0,length_c=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_AnZt_vv__run] tidx %d type_flag %d\n",*tidx,type_flag);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_Zn,verbose," %% M_Zn: ");}
  switch (output_spacing_r){ case SPACING_j: length_r = M_An->rpop_j; break; case SPACING_b: length_r = M_An->rpop_b; break; case SPACING_a: length_r = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: length_c = M_Zn->rpop_j; break; case SPACING_b: length_c = M_Zn->rpop_b; break; case SPACING_a: length_c = M_Zn->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  length = length_r*length_c; if (verbose){ printf(" %% length %d*%d=%d\n",length_r,length_c,length);}
  if (*output_AnZt_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_AnZt_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Zn->mr_b,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_b: ");}
  if (verbose>2){ bprintf(M_Zn->mr_j,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_j: ");}
  if (verbose>2){ bprintf(M_Zn->mc_b,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_b: ");}
  if (verbose>2){ bprintf(M_Zn->mc_j,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_j: ");}
  length = length_r*length_c; if ((*output_AnZt_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_AnZt_vv__run\n",(*output_AnZt_p)->length,length);} memset((*output_AnZt_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_Zn; vpra[ip++] = A_ajdk; vpra[ip++] = lf_An_ajdk; vpra[ip++] = lf_Zn_ajdk; vpra[ip++] = *output_AnZt_p;
  switch (type_flag){ case TYPE_p0: vpra[ip++] = &addressable_type_p0; break; case TYPE_pm: vpra[ip++] = &addressable_type_pm; break; case TYPE_00: vpra[ip++] = &addressable_type_00; break; default: break; /* switch (type_flag){ } */}
  switch (output_spacing_r){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_c){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AnZt_vv__run,vpra)){ printf("Warning! cannot create thread %d in wrap_AnZt_vv__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AnZt_vv__run(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_AnZt_p)->lf,"double",length_r,length_c," %% (*output_AnZt_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_AnZt_vv__run] tidx %d\n",*tidx);}
  return length;
}  

void *get_AnZt_uu__run(void *vp)
{
  /* This function serves as a module which calculates the matrix-matrix product AnZt, where both An and Zt are treated as +1/-1 values ;
     If the terms A_ajdk = (alpha=(p-q)).^j.*(D=1/(4*p*q)).^k are given, as well as lf_An_ajdk and lf_Zn_ajdk (recording the row-vectors storing An*diag(alpha.^j.*D.^k) and Zn*diag(alpha.^j.*D.^k), respectively), ;
     then we calculate (An-e_An*a_At)*diag(D)*(Zt-a_An*e_Zt). ;
     Otherwise, we simply calculate An*Zt. ;
     The calculation is performed using element-wise-operations. ;  
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zn = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_AnZt = (struct L_handle *)(vpra[ip++]);
  int type_flag = *(int *)(vpra[ip++]); 
  int output_spacing_r = *(int *)(vpra[ip++]);
  int output_spacing_c = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols)/* rup(M_An->ncols + M_An->ncols_extend,POPLENGTH)/POPLENGTH */;
  int type_flag_single = 0;
  int mx_j=0,mx_chunk=0,ma_a=0,ma_b=0,ma_j=0,mz_a=0,mz_b=0,mz_j=0,nc_a=0,nc_b=0,nc_j=0,tab_r=0,tab_r_stride=0,tab_c=0,tab_c_stride=0;
  double dtmp=0,output_tmp=0;
  unsigned char *A_tag=NULL,*Z_tag=NULL;
  double *dinp=NULL;
  char tempchar[FNAMESIZE];
  //int zper = maximum(1,minimum(M_Zn->rpop_j,maximum(1,GLOBAL_2_cache_size / maximum(1,M_Zn->mc_length)*3/4))),mz_block=0;//zper = M_Zn->rpop_j;
  int zper = maximum(1,minimum(M_Zn->rpop_j,maximum(1,/*262144*/ GLOBAL_2_cache_size / maximum(1,M_Zn->mc_length) * 3/4 ))),mz_block=0;//zper = M_Zn->rpop_j;
  //int zper = 1, mz_block=0;
  if (type_flag==TYPE_p0){ printf(" %% Warning! TYPE_p0 not implemented within get_AnZt_uu__run\n");}
  if      (type_flag==TYPE_pm){ type_flag_single = TYPE_p_;}
  else if (type_flag==TYPE_00){ type_flag_single = TYPE_0_;}
  if (verbose>1){ printf(" %% [entering get_AnZt_uu__run] tidx %d rpop_b %d(%d) cpop_b %d(%d) xpop_b %d(%d) type_flag %d(%d) (%d,%d)\n",tidx,(int)(M_An->rpop_b),(int)(M_An->rpop_j),(int)(M_An->cpop_b),(int)(M_An->cpop_j),(int)(M_Zn->rpop_b),(int)(M_Zn->rpop_j),type_flag,type_flag_single,(int)zper,(int)(rup(M_Zn->rpop_j,zper)/zper));}  
  if (verbose>1){ printf(" %% zper = %d/%d*3/4 = %d\n",(int)GLOBAL_2_cache_size,(int)M_Zn->mc_length,zper);}
  if (verbose>3){ bprintf((unsigned char *)M_An->wX,POPLENGTH,M_An->rpop_b,minimum(M_An->mc_length*BIT8,M_An->ncols)," %% M_An->wX: ");}
  if (verbose>3){ bprintf((unsigned char *)M_Zn->wX,POPLENGTH,M_Zn->rpop_b,minimum(M_Zn->mc_length*BIT8,M_Zn->ncols)," %% M_Zn->wX: ");}
  unsigned int *ma_b_,*ma_a_;
  unsigned int *mz_b_,*mz_a_;
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  mz_b_ = M_Zn->m_b_; mz_a_ = M_Zn->m_a_;
  if (verbose>1){ sprintf(tempchar," %% ma_b_[%d]: ",tidx); raprintf(ma_b_,"unsigned int",1,M_An->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% ma_a_[%d]: ",tidx); raprintf(ma_a_,"unsigned int",1,M_An->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% mz_b_[%d]: ",tidx); raprintf(mz_b_,"unsigned int",1,M_Zn->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% mz_a_[%d]: ",tidx); raprintf(mz_a_,"unsigned int",1,M_Zn->rpop_j,tempchar);}
  switch (output_spacing_r){ case SPACING_j: tab_r_stride = M_An->rpop_j; break; case SPACING_b: tab_r_stride = M_An->rpop_b; break; case SPACING_a: tab_r_stride = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: tab_c_stride = M_Zn->rpop_j; break; case SPACING_b: tab_c_stride = M_Zn->rpop_b; break; case SPACING_a: tab_c_stride = M_Zn->nrows; break; default: break; /* switch (output_spacing_c){ } */}
  if (verbose>1){ printf(" %% calculating output_AnZt\n");}
  output_AnZt->spacing_row = output_spacing_r; output_AnZt->row_stride = tab_r_stride;
  output_AnZt->spacing_col = output_spacing_c; output_AnZt->col_stride = tab_c_stride;
  memset(output_AnZt->lf,0,tab_r_stride*tab_c_stride*sizeof(double));
  ma_j=0;
  while (ma_j<M_An->rpop_j){
    ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
    switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
    mz_j=0;
    while (mz_j<M_Zn->rpop_j){
      mz_a = mz_a_[mz_j]; mz_b = mz_b_[mz_j];
      switch (output_spacing_c){ case SPACING_j: tab_c=mz_j; break; case SPACING_b: tab_c=mz_b; break; case SPACING_a: tab_c=mz_a; break; default: break; /* switch (output_spacing_c){ } */}
      A_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length]));
      Z_tag = (unsigned char *)(&(M_Zn->wX[mz_b*M_Zn->mc_length]));
      dtmp=0;
      nc_j=0;
      while (nc_j<M_An->cpop_j){
	nc_a = M_An->n_a_[nc_j]; nc_b = M_An->n_b_[nc_j];
      	if      (type_flag==TYPE_pm){ dtmp += (bget____(A_tag,nc_a) - (A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*(bget____(Z_tag,nc_a) - (A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*(A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_0_1*A_pcols] : 1);}
      	else if (type_flag==TYPE_00){ dtmp += (bget__on(A_tag,nc_a) - (A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*(bget__on(Z_tag,nc_a) - (A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*(A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_0_1*A_pcols] : 1);}
	nc_j++; /* while (nc_j<M_An->cpop_j){ } */}
      /* for (nc_a=0;nc_a<M_An->ncols;nc_a++){ */
      /* 	if      (type_flag==TYPE_pm){ dtmp += (bget____(A_tag,nc_a) - (A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*(bget____(Z_tag,nc_a) - (A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*bget__on(M_An->mc_j,nc_a)*(A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_0_1*A_pcols] : 1);} */
      /* 	else if (type_flag==TYPE_00){ dtmp += (bget__on(A_tag,nc_a) - (A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*(bget__on(Z_tag,nc_a) - (A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_1_0*A_pcols] : 0))*bget__on(M_An->mc_j,nc_a)*(A_ajdk ? A_ajdk[nc_a/POPLENGTH + AJDK_0_1*A_pcols] : 1);} */
      /* 	/\* for (nc_a=0;nc_a<M_An->ncols;nc_a++){ } *\/} */
      output_tmp = dtmp;
      output_AnZt->lf[tab_r + tab_c*tab_r_stride] = output_tmp;
      mz_j++;	/* while (mz_j<M_Zn->rpop_j){ } */}
    ma_j++;}
  GLOBAL_ops_count_one(tidx,M_An->rpop_j*M_Zn->rpop_j*M_An->cpop_j*2,0);
  if (verbose>1){ raprintf(output_AnZt->lf,"double",tab_r_stride,tab_c_stride," %% output_AnZt->lf: ");}
  if (verbose>1){ printf(" %% [finished get_AnZt_uu__run] tidx %d \n",tidx);}  
  return NULL;
}

int wrap_AnZt_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,int output_spacing_c,struct M_handle *M_An,struct M_handle *M_Zn,double *A_ajdk,struct L_handle **output_AnZt_p)
{
  /* This function uses the M_handles M_A and M_Zn to run a series of parallel calls to get_AnZt_uu__run ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 8)
     the type of multiplication is determined by type_flag: ;
     type_flag==TYPE_pm: get_AX_pm ;
     type_flag==TYPE_00: get_AX_00 ;
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int length_r=0,length_c=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_AnZt_uu__run] tidx %d type_flag %d\n",*tidx,type_flag);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_Zn,verbose," %% M_Zn: ");}
  switch (output_spacing_r){ case SPACING_j: length_r = M_An->rpop_j; break; case SPACING_b: length_r = M_An->rpop_b; break; case SPACING_a: length_r = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: length_c = M_Zn->rpop_j; break; case SPACING_b: length_c = M_Zn->rpop_b; break; case SPACING_a: length_c = M_Zn->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  length = length_r*length_c; if (verbose){ printf(" %% length %d*%d=%d\n",length_r,length_c,length);}
  if (*output_AnZt_p==NULL){ if (verbose){ printf(" %% allocating output size %d*%d\n",length,(int)sizeof(double));} *output_AnZt_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Zn->mr_b,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_b: ");}
  if (verbose>2){ bprintf(M_Zn->mr_j,M_Zn->bitj,1,M_Zn->nrows," %% M_Zn->mr_j: ");}
  if (verbose>2){ bprintf(M_Zn->mc_b,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_b: ");}
  if (verbose>2){ bprintf(M_Zn->mc_j,M_Zn->bitj,1,M_Zn->ncols," %% M_Zn->mc_j: ");}
  length = length_r*length_c; if ((*output_AnZt_p)->length<length){ printf(" %% Warning! length %d<%d in wrap_AnZt__run\n",(*output_AnZt_p)->length,length);} memset((*output_AnZt_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_Zn; vpra[ip++] = A_ajdk; vpra[ip++] = *output_AnZt_p;
  switch (type_flag){ case TYPE_p0: vpra[ip++] = &addressable_type_p0; break; case TYPE_pm: vpra[ip++] = &addressable_type_pm; break; case TYPE_00: vpra[ip++] = &addressable_type_00; break; default: break; /* switch (type_flag){ } */}
  switch (output_spacing_r){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_c){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AnZt_uu__run,vpra)){ printf("Warning! cannot create thread %d in wrap_AnZt_uu__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AnZt_uu__run(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_AnZt_p)->lf,"double",length_r,length_c," %% (*output_AnZt_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_AnZt_uu__run] tidx %d\n",*tidx);}
  return length;
}  

void wrap_AnZt_vv_test()
{
  /* test for errors with input file: AnZt_vv_error.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= AnZt_vv;
     GLOBAL_TEST_TYP2= error;
     GLOBAL_NBINS= 4;
     GLOBAL_TEST_mrand= 0.5;
     GLOBAL_TEST_A_n_rows= 19;
     GLOBAL_TEST_A_n_cols= 2002;
     GLOBAL_TEST_Z_n_rows= 21;
     GLOBAL_TEST_Y_n_cols= 1950;
     GLOBAL_TEST_T_n_cols= 1;
     GLOBAL_TEST_niter= 1;
     END= 0;
  */
  /* test for speed with input file: AnZt_vv_speed.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= AnZt_vv;
     GLOBAL_TEST_TYP2= speed;
     GLOBAL_NBINS= 4;
     GLOBAL_TEST_mrand= 0.05;
     GLOBAL_TEST_A_n_rows= 1920;
     GLOBAL_TEST_A_n_cols= 19200;
     GLOBAL_TEST_Z_n_rows= 1920;
     GLOBAL_TEST_Y_n_cols= 19200;
     GLOBAL_TEST_T_n_cols= 1;
     GLOBAL_TEST_niter= 1;
     END= 0;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AnZt_vv=NULL,**lf_YnWt_vv=NULL; int *length_AnZt_vv=NULL,*length_YnWt_vv=NULL;
  struct L_handle **lf_AnZt_uu=NULL,**lf_YnWt_uu=NULL; int *length_AnZt_uu=NULL,*length_YnWt_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_AnZt_vv_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  lf_AnZt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnZt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_YnWt_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_YnWt_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_AnZt_vv[nb] = L_handle_make(M_An[nb]->nrows*M_Zn[nb]->nrows);
    lf_YnWt_vv[nb] = L_handle_make(M_An[nb]->nrows*M_Zn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (error_check){  
    lf_AnZt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AnZt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    lf_YnWt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_YnWt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
    for (nb=0;nb<nbins;nb++){
      lf_AnZt_uu[nb] = L_handle_make(M_An[nb]->nrows*M_Zn[nb]->nrows);
      lf_YnWt_uu[nb] = L_handle_make(M_An[nb]->nrows*M_Zn[nb]->nrows);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (error_check){ } */}
  for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
  	if (verbose){ printf(" %% %s; B: %s; A: %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
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
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
  	  length_AnZt_vv[nb] = wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,n_spacing_B,M_An[nb],M_Zn[nb],A_ajdk,lf_An_ajdk[nb],lf_Zn_ajdk[nb],&(lf_AnZt_vv[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
  	  length_YnWt_vv[nb] = wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,n_spacing_B,M_Yn[nb],M_Wn[nb],Y_ajdk,lf_Yn_ajdk[nb],lf_Wn_ajdk[nb],&(lf_YnWt_vv[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){ }} */}}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZt_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic(); 
	    length_AnZt_uu[nb] = wrap_AnZt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,n_spacing_B,M_An[nb],M_Zn[nb],A_ajdk,&(lf_AnZt_uu[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic(); 
	    length_YnWt_uu[nb] = wrap_AnZt_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,n_spacing_B,M_Yn[nb],M_Wn[nb],Y_ajdk,&(lf_YnWt_uu[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZt_uu: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_AnZt_vv[%d] error %0.16f\n",nb,dra_diff(lf_AnZt_vv[nb]->lf,lf_AnZt_uu[nb]->lf,length_AnZt_vv[nb],1));
	    printf(" %% lf_YnWt_vv[%d] error %0.16f\n",nb,dra_diff(lf_YnWt_vv[nb]->lf,lf_YnWt_uu[nb]->lf,length_YnWt_vv[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_AnZt_vv_test]\n");}
}
