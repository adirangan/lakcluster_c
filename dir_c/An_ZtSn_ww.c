#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_An_ZtSn_ww(void *vp)
{
  int verbose=1;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct L_handle *lf_ZtSn = (struct L_handle *)(vpra[ip++]);
  struct M_handle *M_ZtSn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct L_handle *output_An_ZtSn_ww = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_s = *(int *)(vpra[ip++]);
  int ns_j=0,ns_a=0,ns_b=0,tab_s=0,tab_s_stride=0,ma_j=0,ma_a=0,ma_b=0,tab_a=0,tab_a_stride=0,nz_j=0,nz_a=0,nz_b=0,tab_z=0;
  double output_tmp=0;
  int vA=0;
  unsigned char *An_tag=NULL;
  __m128i *wAn_tag=NULL,*wZtSn_tag=NULL,*mcan_tag=NULL,*mcan_end=NULL;
  long long int n2=0;
  double output_An_ZtSn_base=0,lltmp=0,output_An_ZtSn_tmp=0;
  int tab_x=0,mr=0,mx=0;
  int mx_j=0,mx_chunk=0;
  if (verbose>1){ printf(" %% [entering get_An_ZtSn_ww] tidx %d\n",tidx);}
  if (verbose>2){ M_handle_printf(M_An,1," %% M_An: ");}
  if (verbose>2){ M_handle_printf(M_St,1," %% M_St: ");}
  if (verbose>1){ raprintf(lf_ZtSn->lf,"double",lf_ZtSn->row_stride,lf_ZtSn->col_stride*lf_ZtSn->lyr_stride," %% lf_ZtSn->lf: ");}
  if (verbose>1){ printf(" %% POPLENGTH %d M_An->mc_length %d M_ZtSn->mc_length %d\n",POPLENGTH,M_An->mc_length,M_ZtSn->mc_length);}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  if (verbose>1){ printf(" %% tab_s_stride %d tab_a_stride %d\n",tab_s_stride,tab_a_stride);}
  output_An_ZtSn_ww->spacing_row = output_spacing_a; output_An_ZtSn_ww->row_stride = tab_a_stride;
  output_An_ZtSn_ww->spacing_col = output_spacing_s; output_An_ZtSn_ww->col_stride = tab_s_stride;
  if (strstr(GLOBAL_skip,"An_ZtSn_ww")){ goto skip_An_ZtSn_ww;}
  if (GLOBAL_omp_type==GLOBAL_omp_off){ 
    ma_j=0;
    while (ma_j<M_An->rpop_j){
      ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=ma_j; break; case SPACING_b: tab_a=ma_b; break; case SPACING_a: tab_a=ma_a; break; default: break; /* switch (output_spacing_a){ } */}
      wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
      wZtSn_tag = (__m128i*)&(M_ZtSn->wX[0/* start */*M_ZtSn->mc_length]);
      mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
      lltmp = popcount_pm0(&wAn_tag,&wZtSn_tag,&mcan_tag,&mcan_end);
      output_An_ZtSn_base = lltmp*M_ZtSn->min_d*M_ZtSn->mlt_d;
      ns_j=0;
      while (ns_j<M_St->rpop_j){
	ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
	tab_x = tab_a + tab_s*tab_a_stride;
	output_An_ZtSn_tmp = 0;
	mr=M_ZtSn->ncols_per_z*ns_j/* spacing_j */; n2=1;
	for (mx=M_ZtSn->ncols_per_z-1;mx>0;mx--){
	  wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	  wZtSn_tag = (__m128i*)&(M_ZtSn->wX[(mr+mx)*M_ZtSn->mc_length]);
	  mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	  lltmp = popcount_pm0(&wAn_tag,&wZtSn_tag,&mcan_tag,&mcan_end);
	  output_An_ZtSn_tmp += n2*lltmp;
	  n2*=2;
	  /* for (mx=M_ZtSn->ncols_per_z-1;mx>0;mx--){ } */}
	output_An_ZtSn_ww->lf[tab_x] = (double)(output_An_ZtSn_base + output_An_ZtSn_tmp)/M_ZtSn->mlt_d;
	ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
      ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
    GLOBAL_ops_count_one(tidx,M_An->rpop_j*M_St->rpop_j,M_An->rpop_j*M_St->rpop_j*M_ZtSn->ncols_per_z*M_ZtSn->mc_length*BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp__on){ 
    mx_chunk=128;
#pragma omp parallel private(mx_j,ma_j,ma_a,ma_b,tab_a,wAn_tag,wZtSn_tag,mcan_tag,mcan_end,lltmp,output_An_ZtSn_base,ns_j,ns_a,ns_b,tab_s,tab_x,output_An_ZtSn_tmp,mr,n2,mx)
    { /* begin omp parallel */
      mx_j=0;
#pragma omp for schedule(static)
      for (mx_j=0;mx_j<M_An->rpop_j;mx_j++){
	ma_j = mx_j;
	ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=ma_j; break; case SPACING_b: tab_a=ma_b; break; case SPACING_a: tab_a=ma_a; break; default: break; /* switch (output_spacing_a){ } */}
	wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	wZtSn_tag = (__m128i*)&(M_ZtSn->wX[0/* start */*M_ZtSn->mc_length]);
	mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	lltmp = popcount_pm0(&wAn_tag,&wZtSn_tag,&mcan_tag,&mcan_end);
	output_An_ZtSn_base = lltmp*M_ZtSn->min_d*M_ZtSn->mlt_d;
	ns_j=0;
	while (ns_j<M_St->rpop_j){
	  ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	  switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
	  tab_x = tab_a + tab_s*tab_a_stride;
	  output_An_ZtSn_tmp = 0;
	  mr=M_ZtSn->ncols_per_z*ns_j/* spacing_j */; n2=1;
	  for (mx=M_ZtSn->ncols_per_z-1;mx>0;mx--){
	    wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	    wZtSn_tag = (__m128i*)&(M_ZtSn->wX[(mr+mx)*M_ZtSn->mc_length]);
	    mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	    lltmp = popcount_pm0(&wAn_tag,&wZtSn_tag,&mcan_tag,&mcan_end);
	    output_An_ZtSn_tmp += n2*lltmp;
	    n2*=2;
	    /* for (mx=M_ZtSn->ncols_per_z-1;mx>0;mx--){ } */}
	  output_An_ZtSn_ww->lf[tab_x] = (double)(output_An_ZtSn_base + output_An_ZtSn_tmp)/M_ZtSn->mlt_d;
	  ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
	/* for (mx_j=0;mx_j<M_An->rpop_j;mx_j++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,M_An->rpop_j*M_St->rpop_j,M_An->rpop_j*M_St->rpop_j*M_ZtSn->ncols_per_z*M_ZtSn->mc_length*BIT8);    
    /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
 skip_An_ZtSn_ww:
  if (verbose>1){ raprintf(output_An_ZtSn_ww->lf,"double",tab_a_stride,tab_s_stride," %% output_An_ZtSn_ww->lf: ");}
  if (verbose>1){ printf(" %% [finished get_An_ZtSn_ww] tidx %d\n",tidx);}
  return NULL;
}

int wrap_An_ZtSn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,struct M_handle *M_An,struct L_handle *lf_ZtSn,struct M_handle *M_ZtSn,struct M_handle *M_St,struct L_handle **output_An_ZtSn_ww_p)
{
  /* calls get_An_ZtSn_ww;
     variable space in **vpra (should be at least size 8)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_s=0,length=0,ip=0;
  struct M_handle *tmp_M=NULL;
  struct L_handle *tmp_L=NULL;
  char tmpchar[FNAMESIZE];
  if (verbose){ printf(" %% [entering wrap_An_ZtSn_ww__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (spacing_s){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_a,length_s,length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: ");}
  if (verbose>2){ bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: ");}
  if (verbose>2){ bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (*output_An_ZtSn_ww_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_ZtSn_ww_p = L_handle_make(length);}
  if ((*output_An_ZtSn_ww_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_ZtSn_ww__run\n",(*output_An_ZtSn_ww_p)->length,length);}
  memset((*output_An_ZtSn_ww_p)->lf,0,length*sizeof(double));
  if (verbose>2){
    tmp_M = M_An; sprintf(tmpchar,"M_An");
    if (tmp_M!=NULL){ printf(" %% %s: nr_a %d nc_a %d rpop_b %d(%d) cpop_b %d(%d)\n",tmpchar,tmp_M->nrows,tmp_M->ncols,(int)tmp_M->rpop_b,(int)tmp_M->rpop_j,(int)tmp_M->cpop_b,(int)tmp_M->cpop_j);}
    tmp_M = M_St; sprintf(tmpchar,"M_St");
    if (tmp_M!=NULL){ printf(" %% %s: nr_a %d nc_a %d rpop_b %d(%d) cpop_b %d(%d)\n",tmpchar,tmp_M->nrows,tmp_M->ncols,(int)tmp_M->rpop_b,(int)tmp_M->rpop_j,(int)tmp_M->cpop_b,(int)tmp_M->cpop_j);}
    tmp_M = M_ZtSn; sprintf(tmpchar,"M_ZtSn");
    if (tmp_M!=NULL){ printf(" %% %s: nr_a %d nc_a %d rpop_b %d(%d) cpop_b %d(%d) (ncols_per_z %d)\n",tmpchar,tmp_M->nrows,tmp_M->ncols,(int)tmp_M->rpop_b,(int)tmp_M->rpop_j,(int)tmp_M->cpop_b,(int)tmp_M->cpop_j,(int)tmp_M->ncols_per_z);}
    tmp_L = lf_ZtSn; sprintf(tmpchar,"lf_ZtSn");
    if (tmp_L!=NULL){ printf(" %% %s: length %lld spacing_row %d row_stride %d spacing_col %d col_stride %d spacing_lyr %d lyr_stride %d\n",tmpchar,tmp_L->length,(int)tmp_L->spacing_row,(int)tmp_L->row_stride,(int)tmp_L->spacing_col,(int)tmp_L->col_stride,(int)tmp_L->spacing_lyr,(int)tmp_L->lyr_stride);}
    /* if (verbose>2){ } */}
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = lf_ZtSn; vpra[ip++] = M_ZtSn; vpra[ip++] = M_St; vpra[ip++] = *output_An_ZtSn_ww_p;
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_ZtSn_ww,vpra)){ printf("Warning! cannot create thread %d in wrap_An_ZtSn_ww__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_ZtSn_ww(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_An_ZtSn_ww__run] tidx %d\n",*tidx);}
  return length;
}

void *get_An_ZtSn_uu(void *vp)
{
  int verbose=1;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct L_handle *lf_ZtSn = (struct L_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct L_handle *output_An_ZtSn_uu = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_s = *(int *)(vpra[ip++]);
  int ns_j=0,ns_a=0,ns_b=0,tab_s=0,tab_s_stride=0,ma_j=0,ma_a=0,ma_b=0,tab_a=0,tab_a_stride=0,nz_j=0,nz_a=0,nz_b=0,tab_z=0;
  double output_tmp=0;
  int vA=0;
  unsigned char *An_tag=NULL;
  __m128i *wAn_tag=NULL,*wZtSn_tag=NULL,*mcan_tag=NULL,*mcan_end=NULL;
  long long int n2=0;
  double output_An_ZtSn_base=0,lltmp=0,output_An_ZtSn_tmp=0;
  int tab_x=0,mr=0,mx=0;
  if (verbose>1){ printf(" %% [entering get_An_ZtSn_uu] tidx %d\n",tidx);}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  if (verbose>1){ printf(" %% tab_s_stride %d tab_a_stride %d\n",tab_s_stride,tab_a_stride);}
  output_An_ZtSn_uu->spacing_row = output_spacing_a; output_An_ZtSn_uu->row_stride = tab_a_stride;
  output_An_ZtSn_uu->spacing_col = output_spacing_s; output_An_ZtSn_uu->col_stride = tab_s_stride;
  if (verbose>2){ M_handle_printf(M_An,1," %% M_An: ");}
  if (verbose>2){ M_handle_printf(M_St,1," %% M_St: ");}
  if (verbose>1){ raprintf(lf_ZtSn->lf,"double",lf_ZtSn->row_stride,lf_ZtSn->col_stride*lf_ZtSn->lyr_stride," %% lf_ZtSn->lf: ");}
  ma_j=0;
  while (ma_j<M_An->rpop_j){
    ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
    switch (output_spacing_a){ case SPACING_j: tab_a=ma_j; break; case SPACING_b: tab_a=ma_b; break; case SPACING_a: tab_a=ma_a; break; default: break; /* switch (output_spacing_a){ } */}
    An_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length]));
    ns_j=0;
    while (ns_j<M_St->rpop_j){
      ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
      output_tmp=0;
      nz_j=0;
      while (nz_j<M_An->cpop_j){
	nz_a = M_An->n_a_[nz_j]; nz_b = M_An->n_b_[nz_j];
	vA = bget____(An_tag,nz_a);
	output_tmp += vA * (*L2_get(lf_ZtSn,nz_j,nz_b,nz_a,ns_j,ns_b,ns_a));
	nz_j++; /* while (nz_j<M_An->cpop_j){ } */}
      if (verbose>3){ printf("(%d,%d) %lf ",tab_a,tab_s,output_tmp);}
      output_An_ZtSn_uu->lf[tab_a + tab_s*tab_a_stride] = output_tmp;
      ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
    ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
  if (verbose>1){ raprintf(output_An_ZtSn_uu->lf,"double",tab_a_stride,tab_s_stride," %% output_An_ZtSn_uu->lf: ");}
  if (verbose>1){ printf(" %% [finished get_An_ZtSn_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_An_ZtSn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,struct M_handle *M_An,struct L_handle *lf_ZtSn,struct M_handle *M_St,struct L_handle **output_An_ZtSn_uu_p)
{
  /* calls get_An_ZtSn_uu;
     variable space in **vpra (should be at least size 7)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_An_ZtSn_uu__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (spacing_s){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_a,length_s,length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: ");}
  if (verbose>2){ bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: ");}
  if (verbose>2){ bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (*output_An_ZtSn_uu_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_ZtSn_uu_p = L_handle_make(length);}
  if ((*output_An_ZtSn_uu_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_ZtSn_uu__run\n",(*output_An_ZtSn_uu_p)->length,length);}
  memset((*output_An_ZtSn_uu_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = lf_ZtSn; vpra[ip++] = M_St; vpra[ip++] = *output_An_ZtSn_uu_p;
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_ZtSn_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_An_ZtSn_uu__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_ZtSn_uu(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_An_ZtSn_uu__run] tidx %d\n",*tidx);}
  return length;
}

void An_a0d1_ZtSn_excerpt(struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle *lf_ZtSn,struct L_handle *lf_a0d1_ZtSn)
{
  int verbose=0;
  int A_pcols=0,nz_j=0,nz_a=0,nz_b=0,ns_j=0,ns_a=0,ns_b=0;
  if (verbose){ printf(" %% [entering An_a0d1_ZtSn_excerpt]\n");}
  lf_a0d1_ZtSn->spacing_row = SPACING_b;
  lf_a0d1_ZtSn->spacing_col = SPACING_b;
  lf_a0d1_ZtSn->row_stride = M_Zt->rpop_b;
  lf_a0d1_ZtSn->col_stride = M_St->rpop_b;
  L_zero(lf_a0d1_ZtSn); 
  A_pcols = psize(M_Zt->nrows);
  nz_j=0;
  while (nz_j<M_Zt->rpop_j){
    nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
    ns_j=0;
    while (ns_j<M_St->rpop_j){
      ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      if (verbose>1){ printf(" %% %d,%d,%d , %d,%d,%d --> %d\n",nz_j,nz_b,nz_a,ns_j,ns_b,ns_a,nz_a + ns_a*M_Zt->nrows);}
      L2_set(lf_a0d1_ZtSn,nz_j,nz_b,nz_a,ns_j,ns_b,ns_a, A_ajdk[nz_a/POPLENGTH + AJDK_0_1*A_pcols] * (*L2_get(lf_ZtSn,nz_j,nz_b,nz_a,ns_j,ns_b,ns_a)));
      ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
    nz_j++; /* while (nz_j<M_Zt->rpop_j){ } */}
  if (verbose){ printf(" %% [finished An_a0d1_ZtSn_excerpt]\n");}
}

void An_a2d2_ZtSn_excerpt(struct M_handle *M_Zt,struct M_handle *M_St,double *A_ajdk,struct L_handle *lf_ZtSn,struct L_handle *lf_a2d2_ZtSn)
{
  int verbose=0;
  int A_pcols=0,nz_j=0,nz_a=0,nz_b=0,ns_j=0,ns_a=0,ns_b=0;
  if (verbose){ printf(" %% [entering An_a2d2_ZtSn_excerpt]\n");}
  lf_a2d2_ZtSn->spacing_row = SPACING_b;
  lf_a2d2_ZtSn->spacing_col = SPACING_b;
  lf_a2d2_ZtSn->row_stride = M_Zt->rpop_b;
  lf_a2d2_ZtSn->col_stride = M_St->rpop_b;
  L_zero(lf_a2d2_ZtSn); 
  A_pcols = psize(M_Zt->nrows);
  nz_j=0;
  while (nz_j<M_Zt->rpop_j){
    nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
    ns_j=0;
    while (ns_j<M_St->rpop_j){
      ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      if (verbose>1){ printf(" %% %d,%d,%d , %d,%d,%d --> %d\n",nz_j,nz_b,nz_a,ns_j,ns_b,ns_a,nz_a + ns_a*M_Zt->nrows);}
      L2_set(lf_a2d2_ZtSn,nz_j,nz_b,nz_a,ns_j,ns_b,ns_a, A_ajdk[nz_a/POPLENGTH + AJDK_2_2*A_pcols] * (*L2_get(lf_ZtSn,nz_j,nz_b,nz_a,ns_j,ns_b,ns_a)));
      ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
    nz_j++; /* while (nz_j<M_Zt->rpop_j){ } */}
  if (verbose){ printf(" %% [finished An_a2d2_ZtSn_excerpt]\n");}
}

void wrap_An_ZtSn_ww_test()
{
  /* test for errors with input file: An_ZtSn_ww_error.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AtTn=NULL,**lf_ZtSn=NULL;
  struct L_handle **lf_a2d2_AtTn=NULL,**lf_a2d2_ZtSn=NULL;
  int ns_j=0,ns_b=0,ns_a=0;
  struct M_handle **M_a2d2_ZtSn=NULL,**M_a2d2_AtTn=NULL;
  struct L_handle **lf_An_a2d2_ZtSn_ww=NULL,**lf_An_a2d2_AtTn_ww=NULL; int *length_An_a2d2_ZtSn_ww=NULL,*length_An_a2d2_AtTn_ww=NULL;
  struct L_handle **lf_An_a2d2_ZtSn_uu=NULL,**lf_An_a2d2_AtTn_uu=NULL; int *length_An_a2d2_ZtSn_uu=NULL,*length_An_a2d2_AtTn_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_An_ZtSn_ww_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  lf_ZtSn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AtTn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_ZtSn[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_Sn[nb]->ncols);
    lf_AtTn[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_a2d2_ZtSn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_a2d2_AtTn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_a2d2_ZtSn[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_Sn[nb]->ncols);
    lf_a2d2_AtTn[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_a2d2_ZtSn = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  M_a2d2_AtTn = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_a2d2_ZtSn[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Zt[nb]->nrows,M_St[nb]->nrows);
    M_a2d2_AtTn[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_At[nb]->nrows,M_At[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_An_a2d2_ZtSn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_a2d2_ZtSn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_An_a2d2_AtTn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_a2d2_AtTn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_An_a2d2_ZtSn_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_St[nb]->nrows);
    lf_An_a2d2_AtTn_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_An_a2d2_ZtSn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_a2d2_ZtSn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_An_a2d2_AtTn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_a2d2_AtTn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_An_a2d2_ZtSn_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows);
    lf_An_a2d2_AtTn_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
  	if (verbose){ printf(" %% %s; %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_An[nb],A_ajdk,&(lf_An_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Zn[nb],A_ajdk,&(lf_Zn_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; 
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],NULL,NULL,NULL,&(lf_AtTn[nb]));
	  GLOBAL_pthread_toc(); 
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],NULL,NULL,NULL,&(lf_ZtSn[nb]));
	  GLOBAL_pthread_toc(); 
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc();
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	for (nb=0;nb<nbins;nb++){ 
	  An_a2d2_ZtSn_excerpt(M_Zt[nb],M_St[nb],A_ajdk,lf_ZtSn[nb],lf_a2d2_ZtSn[nb]);
	  An_a2d2_ZtSn_excerpt(M_At[nb],M_Tt[nb],A_ajdk,lf_AtTn[nb],lf_a2d2_AtTn[nb]);
	  /* for (nb=0;nb<nbins;nb++){ } */}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic(); 
	  wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Zt[nb]->mr_j,M_Zt[nb]->mr_b,M_St[nb]->mr_j,M_St[nb]->mr_b,lf_a2d2_ZtSn[nb],lf_a2d2_ZtSn[nb]->lf,&(M_Zt[nb]->nrows),&(M_St[nb]->nrows),&(M_a2d2_ZtSn[nb]),&(GLOBAL_B_MLT),(addressable_0));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_At[nb]->mr_j,M_At[nb]->mr_b,M_Tt[nb]->mr_j,M_Tt[nb]->mr_b,lf_a2d2_AtTn[nb],lf_a2d2_AtTn[nb]->lf,&(M_At[nb]->nrows),&(M_Tt[nb]->nrows),&(M_a2d2_AtTn[nb]),&(GLOBAL_B_MLT),(addressable_0));
	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% xcalc: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
	  length_An_a2d2_ZtSn_ww[nb] = wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],lf_a2d2_ZtSn[nb],M_a2d2_ZtSn[nb],M_St[nb],&(lf_An_a2d2_ZtSn_ww[nb]));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  length_An_a2d2_AtTn_ww[nb] = wrap_An_ZtSn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],lf_a2d2_AtTn[nb],M_a2d2_AtTn[nb],M_Tt[nb],&(lf_An_a2d2_AtTn_ww[nb]));
	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_a2d2_ZtSn_ww An_a2d2_AtTn_ww: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){	
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic();
	    length_An_a2d2_ZtSn_uu[nb] = wrap_An_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],lf_a2d2_ZtSn[nb],M_St[nb],&(lf_An_a2d2_ZtSn_uu[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_An_a2d2_AtTn_uu[nb] = wrap_An_ZtSn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_An[nb],lf_a2d2_AtTn[nb],M_Tt[nb],&(lf_An_a2d2_AtTn_uu[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc();
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_a2d2_ZtSn_uu An_a2d2_AtTn_uu: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_An_a2d2_ZtSn_ww[%d] error %0.16f\n",nb,dra_diff(lf_An_a2d2_ZtSn_ww[nb]->lf,lf_An_a2d2_ZtSn_uu[nb]->lf,length_An_a2d2_ZtSn_ww[nb],1));
	    printf(" %% lf_An_a2d2_AtTn_ww[%d] error %0.16f\n",nb,dra_diff(lf_An_a2d2_AtTn_ww[nb]->lf,lf_An_a2d2_AtTn_uu[nb]->lf,length_An_a2d2_AtTn_ww[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_An_ZtSn_ww_test]\n");}
}
