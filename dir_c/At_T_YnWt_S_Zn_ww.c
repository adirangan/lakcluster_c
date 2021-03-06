#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_At_T_YnWt_S_Zn_ww(void *vp)
{
  /* no preallocation, using inline computation */
  int verbose=1;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Yn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Wn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct M_handle *M_YnWt = (struct M_handle *)(vpra[ip++]);
  struct L_handle *output_At_T_YnWt_S_Zn_ww = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_t = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows)/* rup(M_At->nrows + M_At->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int nt_j=0,nt_a=0,nt_b=0,ns_j=0,ns_a=0,ns_b=0,tab_t=0,tab_t_stride=0,na_j=0,na_a=0,na_b=0,nz_j=0,nz_a=0,nz_b=0,tab_a=0,tab_a_stride=0,mw_j=0,mw_a=0,mw_b=0;
  double output_tmp=0;
  int vZ=0,vS=0;
  unsigned char *St_tag=NULL,*Zt_tag=NULL;
  __m128i *wAt_tag=NULL,*wTt_tag=NULL,*wYnWt_tag=NULL,*mcat_tag=NULL,*mcat_end=NULL;
  long long int lltmp=0,n2=0;
  double dtmp=0;
  double output_At_T_YnWt_base=0,output_at_T_YnWt_base=0,output____T_YnWt_base=0,output_At_T_YnWt_tmp=0,output_at_T_YnWt_tmp=0,output____T_YnWt_tmp_[M_Wn->rpop_j*M_St->rpop_j],output_At_T_YnWt______tmp=0;
  int tab_x=0,mr=0,mx=0 /* ncols_per */,mx_j=0 /* omp */;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering get_At_T_YnWt_S_Zn_ww] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% Calculating output_At_T_YnWt_uu\n");}
  switch (output_spacing_t){ case SPACING_j: tab_t_stride = M_Tt->rpop_j; break; case SPACING_b: tab_t_stride = M_Tt->rpop_b; break; case SPACING_a: tab_t_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  if (verbose>1){ printf(" %% tab_t_stride %d tab_a_stride %d\n",tab_t_stride,tab_a_stride);}
  output_At_T_YnWt_S_Zn_ww->spacing_row = output_spacing_a; output_At_T_YnWt_S_Zn_ww->row_stride = tab_a_stride;
  output_At_T_YnWt_S_Zn_ww->spacing_col = output_spacing_t; output_At_T_YnWt_S_Zn_ww->col_stride = tab_t_stride;
  if (strstr(GLOBAL_skip,"At_T_YnWt_S_Zn_ww")){ goto skip_At_T_YnWt_S_Zn_ww;}
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    nt_j=0; ns_j=0;
    while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){
      nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      switch (output_spacing_t){ case SPACING_j: tab_t=nt_j; break; case SPACING_b: tab_t=nt_b; break; case SPACING_a: tab_t=nt_a; break; default: break; /* switch (output_spacing_t){ } */}
      St_tag = (unsigned char *)&(M_St->wX[ns_b*M_St->mc_length]);
      wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
      wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
      wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
      mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
      lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
      output____T_YnWt_base = dtmp*M_YnWt->min_d*M_YnWt->mlt_d;
      memset(output____T_YnWt_tmp_,0,M_Wn->rpop_j*sizeof(double));
      mw_j=0;
      while (mw_j<M_Wn->rpop_j){
	mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	output_at_T_YnWt_tmp = 0; 
	mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
	for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	  wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	  wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	  wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	  mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	  lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
	  output_at_T_YnWt_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
	output____T_YnWt_tmp_[mw_j] = output_at_T_YnWt_tmp;
	mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
      na_j=0; nz_j=0;
      while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	tab_x = tab_a + tab_t*tab_a_stride;
	output_At_T_YnWt_S_Zn_ww->lf[tab_x] = 0;
	Zt_tag = (unsigned char *)&(M_Zt->wX[nz_b*M_Zt->mc_length]);
	wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
	output_At_T_YnWt_base = lltmp*M_YnWt->min_d*M_YnWt->mlt_d;
	mw_j=0;
	while (mw_j<M_Wn->rpop_j){
	  mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	  vS = bget____(St_tag,mw_a);
	  vZ = bget____(Zt_tag,mw_a);
	  output_At_T_YnWt_tmp = 0;
	  mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
	  for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	    wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	    wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	    wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	    mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	    lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
	    output_At_T_YnWt_tmp += n2*lltmp;
	    n2*=2;
	    /* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
	  output_at_T_YnWt_base = output____T_YnWt_base*a_An[na_a/POPLENGTH];
	  output_at_T_YnWt_tmp = output____T_YnWt_tmp_[mw_j]*a_An[na_a/POPLENGTH];
	  output_At_T_YnWt______tmp = (double)(output_At_T_YnWt_base + output_At_T_YnWt_tmp)/M_YnWt->mlt_d - (double)(output_at_T_YnWt_base + output_at_T_YnWt_tmp)/M_YnWt->mlt_d;
	  output_At_T_YnWt_S_Zn_ww->lf[tab_x] += output_At_T_YnWt______tmp*vS*(vZ - a_An[nz_a/POPLENGTH]);
	  mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
	na_j++; nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
      nt_j++; ns_j++; /* while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){ } */}
    if (verbose>1){ raprintf(output_At_T_YnWt_S_Zn_ww->lf,"double",tab_a_stride,tab_t_stride," %% output_At_T_YnWt_S_Zn_ww->lf: ");}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Wn->rpop_j,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Wn->rpop_j*(unsigned long long int)M_YnWt->ncols_per_z*(unsigned long long int)M_YnWt->mc_length*(unsigned long long int)BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  if (GLOBAL_omp_type==GLOBAL_omp__on){
    memset(output____T_YnWt_tmp_,0,M_Wn->rpop_j*M_St->rpop_j*sizeof(double));
    nt_j=0; ns_j=0;
    while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){
      nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      mw_j=0;
      while (mw_j<M_Wn->rpop_j){
	mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	output_at_T_YnWt_tmp = 0; 
	mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
	for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	  wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	  wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	  wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	  mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	  lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
	  output_at_T_YnWt_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
	output____T_YnWt_tmp_[mw_j + ns_j*M_Wn->rpop_j] = output_at_T_YnWt_tmp;
	mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
      nt_j++; ns_j++; /* while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){ } */}
#pragma omp parallel private(mx_j,nt_j,nt_a,nt_b,ns_j,ns_a,ns_b,na_j,na_b,na_a,nz_j,nz_a,nz_b,mw_j,mw_b,mw_a,mr,mx,tab_t,tab_a,tab_x,St_tag,Zt_tag,wAt_tag,wTt_tag,wYnWt_tag,mcat_tag,mcat_end,vS,vZ,n2,lltmp,dtmp,output____T_YnWt_base,output_At_T_YnWt_base,output_at_T_YnWt_base,output_At_T_YnWt_tmp,output_at_T_YnWt_tmp,output_At_T_YnWt______tmp)
    { /* begin omp parallel */
      mx_j=0;
#pragma omp for schedule(dynamic)
      for (mx_j=0;mx_j<M_At->rpop_j*M_Tt->rpop_j;mx_j++){
	nt_j = mx_j / M_At->rpop_j; ns_j = nt_j;
	na_j = mx_j % M_At->rpop_j; nz_j = na_j;
	nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	switch (output_spacing_t){ case SPACING_j: tab_t=nt_j; break; case SPACING_b: tab_t=nt_b; break; case SPACING_a: tab_t=nt_a; break; default: break; /* switch (output_spacing_t){ } */}
	St_tag = (unsigned char *)&(M_St->wX[ns_b*M_St->mc_length]);
	wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
	output____T_YnWt_base = dtmp*M_YnWt->min_d*M_YnWt->mlt_d;
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	tab_x = tab_a + tab_t*tab_a_stride;
	output_At_T_YnWt_S_Zn_ww->lf[tab_x] = 0;
	Zt_tag = (unsigned char *)&(M_Zt->wX[nz_b*M_Zt->mc_length]);
	wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
	output_At_T_YnWt_base = lltmp*M_YnWt->min_d*M_YnWt->mlt_d;
	mw_j=0;
	while (mw_j<M_Wn->rpop_j){
	  mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	  vS = bget____(St_tag,mw_a);
	  vZ = bget____(Zt_tag,mw_a);
	  output_At_T_YnWt_tmp = 0;
	  mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
	  for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	    wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	    wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	    wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	    mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	    lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
	    output_At_T_YnWt_tmp += n2*lltmp;
	    n2*=2;
	    /* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
	  output_at_T_YnWt_base = output____T_YnWt_base*a_An[na_a/POPLENGTH];
	  output_at_T_YnWt_tmp = output____T_YnWt_tmp_[mw_j + ns_j*M_Wn->rpop_j]*a_An[na_a/POPLENGTH];
	  output_At_T_YnWt______tmp = (double)(output_At_T_YnWt_base + output_At_T_YnWt_tmp)/M_YnWt->mlt_d - (double)(output_at_T_YnWt_base + output_at_T_YnWt_tmp)/M_YnWt->mlt_d;
	  output_At_T_YnWt_S_Zn_ww->lf[tab_x] += output_At_T_YnWt______tmp*vS*(vZ - a_An[nz_a/POPLENGTH]);
	  mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
	/* for (mx_j=0;mx_j<M_At->rpop_j*M_Tt->rpop_j;mx_j++){ } */}
      /* end omp parallel */}
    if (verbose>1){ raprintf(output_At_T_YnWt_S_Zn_ww->lf,"double",tab_a_stride,tab_t_stride," %% output_At_T_YnWt_S_Zn_ww->lf: ");}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Wn->rpop_j,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Wn->rpop_j*(unsigned long long int)M_YnWt->ncols_per_z*(unsigned long long int)M_YnWt->mc_length*(unsigned long long int)BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
 skip_At_T_YnWt_S_Zn_ww:
  if (verbose>1){ printf(" %% [finished get_At_T_YnWt_S_Zn_ww] tidx %d\n",tidx);}
  return NULL;
}

void *get_At_T_YnWt_S_Zn_ww_precompute(void *vp)
{
  /* Preallocation and precomputation, rather than inline computation. ;
     This is slightly faster than the inline computation, but does not ;
     parallelize as well due to the additional memory requirements. ; */
  int verbose=1;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Yn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Wn = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct M_handle *M_YnWt = (struct M_handle *)(vpra[ip++]);
  struct L_handle *output_At_T_YnWt_S_Zn_ww = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_t = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows)/* rup(M_At->nrows + M_At->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int nt_j=0,nt_a=0,nt_b=0,ns_j=0,ns_a=0,ns_b=0,tab_t=0,tab_t_stride=0,na_j=0,na_a=0,na_b=0,nz_j=0,nz_a=0,nz_b=0,tab_a=0,tab_a_stride=0,mw_j=0,mw_a=0,mw_b=0;
  double output_tmp=0;
  int vZ=0,vS=0;
  unsigned char *St_tag=NULL,*Zt_tag=NULL;
  __m128i *wAt_tag=NULL,*wTt_tag=NULL,*wYnWt_tag=NULL,*mcat_tag=NULL,*mcat_end=NULL;
  long long int lltmp=0,n2=0;
  double output_At_T_YnWt_base=0,output_At_T_YnWt_base_[M_Zt->rpop_j*M_St->rpop_j];
  double output_at_T_YnWt_base=0,output_at_T_YnWt_base_[M_At->rpop_j*M_St->rpop_j];
  double output____T_YnWt_base=0,output____T_YnWt_base_[M_St->rpop_j];
  double output____T_YnWt_tmp,output____T_YnWt_tmp_[M_Wn->rpop_j*M_St->rpop_j];
  double output_At_T_YnWt_tmp=0,output_at_T_YnWt_tmp=0,output_At_T_YnWt______tmp=0;
  double dtmp=0;
  int tab_x=0,mr=0,mx=0 /*per_z*/,mx_j=0 /*omp*/;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering get_At_T_YnWt_S_Zn_ww] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% Calculating output_At_T_YnWt_uu\n");}
  switch (output_spacing_t){ case SPACING_j: tab_t_stride = M_Tt->rpop_j; break; case SPACING_b: tab_t_stride = M_Tt->rpop_b; break; case SPACING_a: tab_t_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  if (verbose>1){ printf(" %% tab_t_stride %d tab_a_stride %d\n",tab_t_stride,tab_a_stride);}
  output_At_T_YnWt_S_Zn_ww->spacing_row = output_spacing_a; output_At_T_YnWt_S_Zn_ww->row_stride = tab_a_stride;
  output_At_T_YnWt_S_Zn_ww->spacing_col = output_spacing_t; output_At_T_YnWt_S_Zn_ww->col_stride = tab_t_stride;
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    memset(output____T_YnWt_base_,0,M_St->rpop_j*sizeof(double));
    memset(output____T_YnWt_tmp_,0,M_Wn->rpop_j*M_St->rpop_j*sizeof(double));
    nt_j=0; ns_j=0;
    while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){
      nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
      wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
      wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
      mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
      lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
      output____T_YnWt_base_[ns_j] = dtmp*M_YnWt->min_d*M_YnWt->mlt_d;
      mw_j=0;
      while (mw_j<M_Wn->rpop_j){
	mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	output_at_T_YnWt_tmp = 0; 
	mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
	for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	  wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	  wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	  wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	  mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	  lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
	  output_at_T_YnWt_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
	output____T_YnWt_tmp_[mw_j + ns_j*M_Wn->rpop_j] = output_at_T_YnWt_tmp;
	mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
      nt_j++; ns_j++; /* while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){ } */}
    memset(output_at_T_YnWt_base_,0,M_At->rpop_j*M_Tt->rpop_j*sizeof(double));
    nt_j=0; ns_j=0;
    while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){
      nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      na_j=0; nz_j=0;
      while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
	output_at_T_YnWt_base_[na_j + ns_j*M_Zt->rpop_j] = output____T_YnWt_base_[ns_j]*a_An[na_a/POPLENGTH];
	na_j++; nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
      nt_j++; ns_j++; /* while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){ } */}
    memset(output_At_T_YnWt_base_,0,M_Zt->rpop_j*M_St->rpop_j*sizeof(double));
    nt_j=0; ns_j=0;
    while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){
      nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      na_j=0; nz_j=0;
      while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
	Zt_tag = (unsigned char *)&(M_Zt->wX[nz_b*M_Zt->mc_length]);
	wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
	output_At_T_YnWt_base_[nz_j + ns_j*M_Zt->rpop_j] = lltmp*M_YnWt->min_d*M_YnWt->mlt_d;
	na_j++; nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
      nt_j++; ns_j++; /* while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){ } */}
    nt_j=0; ns_j=0;
    while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){
      nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      switch (output_spacing_t){ case SPACING_j: tab_t=nt_j; break; case SPACING_b: tab_t=nt_b; break; case SPACING_a: tab_t=nt_a; break; default: break; /* switch (output_spacing_t){ } */}
      St_tag = (unsigned char *)&(M_St->wX[ns_b*M_St->mc_length]);
      output____T_YnWt_base = output____T_YnWt_base_[ns_j];
      na_j=0; nz_j=0;
      while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	output_At_T_YnWt_base = output_At_T_YnWt_base_[nz_j + ns_j*M_Zt->rpop_j];
	output_at_T_YnWt_base = output____T_YnWt_base*a_An[na_a/POPLENGTH];
	tab_x = tab_a + tab_t*tab_a_stride;
	output_At_T_YnWt_S_Zn_ww->lf[tab_x] = 0;
	Zt_tag = (unsigned char *)&(M_Zt->wX[nz_b*M_Zt->mc_length]);
	mw_j=0;
	while (mw_j<M_Wn->rpop_j){
	  mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	  output____T_YnWt_tmp = output____T_YnWt_tmp_[mw_j + ns_j*M_Wn->rpop_j];
	  vS = bget____(St_tag,mw_a);
	  vZ = bget____(Zt_tag,mw_a);
	  output_At_T_YnWt_tmp = 0;
	  mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
	  for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	    wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	    wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	    wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	    mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	    lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
	    output_At_T_YnWt_tmp += n2*lltmp;
	    n2*=2;
	    /* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
	  output_at_T_YnWt_tmp = output____T_YnWt_tmp*a_An[na_a/POPLENGTH];
	  output_At_T_YnWt______tmp = (double)(output_At_T_YnWt_base + output_At_T_YnWt_tmp)/M_YnWt->mlt_d - (double)(output_at_T_YnWt_base + output_at_T_YnWt_tmp)/M_YnWt->mlt_d;
	  output_At_T_YnWt_S_Zn_ww->lf[tab_x] += output_At_T_YnWt______tmp*vS*(vZ - a_An[nz_a/POPLENGTH]);
	  mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
	na_j++; nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
      nt_j++; ns_j++; /* while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){ } */}
    if (verbose>1){ raprintf(output_At_T_YnWt_S_Zn_ww->lf,"double",tab_a_stride,tab_t_stride," %% output_At_T_YnWt_S_Zn_ww->lf: ");}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Wn->rpop_j,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Wn->rpop_j*(unsigned long long int)M_YnWt->ncols_per_z*(unsigned long long int)M_YnWt->mc_length*(unsigned long long int)BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  if (GLOBAL_omp_type==GLOBAL_omp__on){
    memset(output____T_YnWt_base_,0,M_St->rpop_j*sizeof(double));
    memset(output____T_YnWt_tmp_,0,M_Wn->rpop_j*M_St->rpop_j*sizeof(double));
    nt_j=0; ns_j=0;
    while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){
      nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
      wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
      wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
      mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
      lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
      output____T_YnWt_base_[ns_j] = dtmp*M_YnWt->min_d*M_YnWt->mlt_d;
      mw_j=0;
      while (mw_j<M_Wn->rpop_j){
	mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	output_at_T_YnWt_tmp = 0; 
	mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
	for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	  wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	  wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	  wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	  mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	  lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
	  output_at_T_YnWt_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
	output____T_YnWt_tmp_[mw_j + ns_j*M_Wn->rpop_j] = output_at_T_YnWt_tmp;
	mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
      nt_j++; ns_j++; /* while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){ } */}
    memset(output_at_T_YnWt_base_,0,M_At->rpop_j*M_Tt->rpop_j*sizeof(double));
    nt_j=0; ns_j=0;
    while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){
      nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      na_j=0; nz_j=0;
      while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
	output_at_T_YnWt_base_[na_j + ns_j*M_Zt->rpop_j] = output____T_YnWt_base_[ns_j]*a_An[na_a/POPLENGTH];
	na_j++; nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
      nt_j++; ns_j++; /* while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){ } */}
    memset(output_At_T_YnWt_base_,0,M_Zt->rpop_j*M_St->rpop_j*sizeof(double));
    nt_j=0; ns_j=0;
    while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){
      nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      na_j=0; nz_j=0;
      while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
	Zt_tag = (unsigned char *)&(M_Zt->wX[nz_b*M_Zt->mc_length]);
	wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
	output_At_T_YnWt_base_[nz_j + ns_j*M_Zt->rpop_j] = lltmp*M_YnWt->min_d*M_YnWt->mlt_d;
	na_j++; nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
      nt_j++; ns_j++; /* while (nt_j<M_Tt->rpop_j && ns_j<M_St->rpop_j){ } */}
#pragma omp parallel private(mx_j,nt_j,nt_a,nt_b,ns_j,ns_a,ns_b,na_j,na_b,na_a,nz_j,nz_a,nz_b,mw_j,mw_b,mw_a,mr,mx,tab_t,tab_a,tab_x,St_tag,Zt_tag,wAt_tag,wTt_tag,wYnWt_tag,mcat_tag,mcat_end,vS,vZ,n2,lltmp,output____T_YnWt_base,output_At_T_YnWt_base,output_at_T_YnWt_base,output____T_YnWt_tmp,output_At_T_YnWt_tmp,output_at_T_YnWt_tmp,output_At_T_YnWt______tmp)
    { /* begin omp parallel */
      mx_j=0;
#pragma omp for schedule(dynamic)
      for (mx_j=0;mx_j<M_At->rpop_j*M_Tt->rpop_j;mx_j++){
	nt_j = mx_j / M_At->rpop_j; ns_j = nt_j;
	na_j = mx_j % M_At->rpop_j; nz_j = na_j;
	nt_a = M_Tt->m_a_[nt_j]; nt_b = M_Tt->m_b_[nt_j]; ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
	switch (output_spacing_t){ case SPACING_j: tab_t=nt_j; break; case SPACING_b: tab_t=nt_b; break; case SPACING_a: tab_t=nt_a; break; default: break; /* switch (output_spacing_t){ } */}
	St_tag = (unsigned char *)&(M_St->wX[ns_b*M_St->mc_length]);
	output____T_YnWt_base = output____T_YnWt_base_[ns_j];
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	output_At_T_YnWt_base = output_At_T_YnWt_base_[nz_j + ns_j*M_Zt->rpop_j];
	output_at_T_YnWt_base = output____T_YnWt_base*a_An[na_a/POPLENGTH];
	tab_x = tab_a + tab_t*tab_a_stride;
	output_At_T_YnWt_S_Zn_ww->lf[tab_x] = 0;
	Zt_tag = (unsigned char *)&(M_Zt->wX[nz_b*M_Zt->mc_length]);
	mw_j=0;
	while (mw_j<M_Wn->rpop_j){
	  mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	  output____T_YnWt_tmp = output____T_YnWt_tmp_[mw_j + ns_j*M_Wn->rpop_j];
	  vS = bget____(St_tag,mw_a);
	  vZ = bget____(Zt_tag,mw_a);
	  output_At_T_YnWt_tmp = 0;
	  mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
	  for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	    wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	    wTt_tag = (__m128i*)&(M_Tt->wX[nt_b*M_Tt->mc_length]);
	    wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	    mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	    lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
	    output_At_T_YnWt_tmp += n2*lltmp;
	    n2*=2;
	    /* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
	  output_at_T_YnWt_tmp = output____T_YnWt_tmp*a_An[na_a/POPLENGTH];
	  output_At_T_YnWt______tmp = (double)(output_At_T_YnWt_base + output_At_T_YnWt_tmp)/M_YnWt->mlt_d - (double)(output_at_T_YnWt_base + output_at_T_YnWt_tmp)/M_YnWt->mlt_d;
	  output_At_T_YnWt_S_Zn_ww->lf[tab_x] += output_At_T_YnWt______tmp*vS*(vZ - a_An[nz_a/POPLENGTH]);
	  mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
	/* for (mx_j=0;mx_j<M_At->rpop_j*M_Tt->rpop_j;mx_j++){ } */}
      /* end omp parallel */}
    if (verbose>1){ raprintf(output_At_T_YnWt_S_Zn_ww->lf,"double",tab_a_stride,tab_t_stride," %% output_At_T_YnWt_S_Zn_ww->lf: ");}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Wn->rpop_j,(unsigned long long int)M_Tt->rpop_j*(unsigned long long int)M_At->rpop_j*(unsigned long long int)M_Wn->rpop_j*(unsigned long long int)M_YnWt->ncols_per_z*(unsigned long long int)M_YnWt->mc_length*(unsigned long long int)BIT8);
    /* if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
  if (verbose>1){ printf(" %% [finished get_At_T_YnWt_S_Zn_ww] tidx %d\n",tidx);}
  return NULL;
}

int wrap_At_T_YnWt_S_Zn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_t,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yn,struct M_handle *M_Wn,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,struct M_handle *M_YnWt,struct L_handle **output_At_T_YnWt_S_Zn_ww_p)
{
  /* calls get_At_T_YnWt_S_Zn_ww;
     variable space in **vpra (should be at least size 13)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_t=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_At_T_YnWt_S_Zn_ww__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_t){ case SPACING_j: length_t = M_Tt->rpop_j; break; case SPACING_b: length_t = M_Tt->rpop_b; break; case SPACING_a: length_t = M_Tt->nrows; break; default: break; /* switch (spacing_t){ } */}
  length = length_a*length_t; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_a,length_t,length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: ");}
  if (verbose>2){ bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: ");}
  if (verbose>2){ bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: ");}
  if (verbose>2){ bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: ");}
  if (verbose>2){ bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (*output_At_T_YnWt_S_Zn_ww_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_At_T_YnWt_S_Zn_ww_p = L_handle_make(length);}
  if ((*output_At_T_YnWt_S_Zn_ww_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_At_T_YnWt_S_Zn_ww__run\n",(*output_At_T_YnWt_S_Zn_ww_p)->length,length);}
  memset((*output_At_T_YnWt_S_Zn_ww_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Yn; vpra[ip++] = M_Wn; vpra[ip++] = M_St; vpra[ip++] = M_Zt; vpra[ip++] = A_ajdk; vpra[ip++] = M_YnWt; vpra[ip++] = *output_At_T_YnWt_S_Zn_ww_p;
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_t){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_t){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_At_T_YnWt_S_Zn_ww,vpra)){ printf("Warning! cannot create thread %d in wrap_At_T_YnWt_S_Zn_ww__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_At_T_YnWt_S_Zn_ww(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_At_T_YnWt_S_Zn_ww__run] tidx %d\n",*tidx);}
  return length;
}

void wrap_At_T_YnWt_S_Zn_ww_test()
{
  /* test for errors with input file: At_T_YnWt_S_Zn_ww_error.in ;
  */
  /* test for speed with input file: At_T_YnWt_S_Zn_ww_speed.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AnZt=NULL,**lf_AnAt=NULL,**lf_YnWt=NULL,**lf_YnYt=NULL;
  struct L_handle **lf_AtTn=NULL,**lf_YtTn=NULL,**lf_ZtSn=NULL,**lf_WtSn=NULL;
  struct L_handle **lf_ZtSWn=NULL,**lf_AtTYn=NULL,**lf_ZtSZn=NULL,**lf_AtTAn=NULL;
  struct M_handle **M_YnWt=NULL,**M_AnZt=NULL,**M_WnYt=NULL,**M_ZnAt=NULL;
  struct L_handle **lf_At_T_YnWt_S_Zn_ww=NULL,**lf_At_T_AnZt_S_Zn_ww=NULL; int *length_At_T_YnWt_S_Zn_ww=NULL,*length_At_T_AnZt_S_Zn_ww=NULL;
  struct L_handle **lf_Zt_S_WnYt_T_An_ww=NULL,**lf_Zt_S_ZnAt_T_An_ww=NULL; int *length_Zt_S_WnYt_T_An_ww=NULL,*length_Zt_S_ZnAt_T_An_ww=NULL;
  struct L_handle **lf_AtTYn____WtSZn_vv=NULL,**lf_AtTAn____ZtSZn_vv=NULL; int *length_AtTYn____WtSZn_vv=NULL,*length_AtTAn____ZtSZn_vv=NULL;
  struct L_handle **lf_At_T_YnWt_S_Zn_vv=NULL,**lf_At_T_AnZt_S_Zn_vv=NULL; int *length_At_T_YnWt_S_Zn_vv=NULL,*length_At_T_AnZt_S_Zn_vv=NULL;
  struct L_handle **lf_AtTYnWtSZn_uu=NULL,**lf_AtTAnZtSZn_uu=NULL; int *length_AtTYnWtSZn_uu=NULL,*length_AtTAnZtSZn_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_At_T_YnWt_S_Zn_ww_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  lf_AtTn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_YtTn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_ZtSn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_WtSn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_AtTn[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    lf_YtTn[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    lf_ZtSn[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    lf_WtSn[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->ncols*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_YnWt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  M_AnZt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  M_WnYt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  M_ZnAt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_YnWt[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Yn[nb]->nrows,M_Wn[nb]->nrows);
    M_AnZt[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_An[nb]->nrows,M_Zn[nb]->nrows);
    M_WnYt[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Wn[nb]->nrows,M_Yn[nb]->nrows);
    M_ZnAt[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Zn[nb]->nrows,M_An[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_At_T_YnWt_S_Zn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_YnWt_S_Zn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_At_T_AnZt_S_Zn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_AnZt_S_Zn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Zt_S_WnYt_T_An_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Zt_S_WnYt_T_An_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Zt_S_ZnAt_T_An_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Zt_S_ZnAt_T_An_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_At_T_YnWt_S_Zn_ww[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
    lf_At_T_AnZt_S_Zn_ww[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
    lf_Zt_S_WnYt_T_An_ww[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
    lf_Zt_S_ZnAt_T_An_ww[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_ZtSWn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AtTYn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_ZtSZn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AtTAn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_ZtSWn[nb] = L_handle_make((unsigned long long int)M_Zt[nb]->nrows*(unsigned long long int)M_Wt[nb]->nrows*(unsigned long long int)M_St[nb]->nrows);
    lf_AtTYn[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Yt[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows);
    lf_ZtSZn[nb] = L_handle_make((unsigned long long int)M_Zt[nb]->nrows*(unsigned long long int)M_Zt[nb]->nrows*(unsigned long long int)M_St[nb]->nrows);
    lf_AtTAn[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_AtTYn____WtSZn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AtTYn____WtSZn_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_AtTAn____ZtSZn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AtTAn____ZtSZn_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_AtTYn____WtSZn_vv[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
    lf_AtTAn____ZtSZn_vv[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_AnZt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AnAt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_YnWt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_YnYt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_AnZt[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    lf_AnAt[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_An[nb]->nrows);
    lf_YnWt[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_Wn[nb]->nrows);
    lf_YnYt[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_Yn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_At_T_YnWt_S_Zn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_YnWt_S_Zn_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_At_T_AnZt_S_Zn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_AnZt_S_Zn_vv = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_At_T_YnWt_S_Zn_vv[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
    lf_At_T_AnZt_S_Zn_vv[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_AtTYnWtSZn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AtTYnWtSZn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_AtTAnZtSZn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_AtTAnZtSZn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_AtTYnWtSZn_uu[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
    lf_AtTAnZtSZn_uu[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
	if (verbose){ printf(" %% %s; %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
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
	if (error_check){
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic();
	    wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],M_Yt[nb],A_ajdk,Y_ajdk,lf_AtTn[nb],lf_YtTn[nb],&(lf_AtTYn[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],M_At[nb],A_ajdk,A_ajdk,lf_AtTn[nb],lf_AtTn[nb],&(lf_AtTAn[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],M_Wt[nb],A_ajdk,Y_ajdk,lf_ZtSn[nb],lf_WtSn[nb],&(lf_ZtSWn[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],M_Zt[nb],A_ajdk,A_ajdk,lf_ZtSn[nb],lf_ZtSn[nb],&(lf_ZtSZn[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc();
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYn and ZtSWn: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic();
	    length_AtTYn____WtSZn_vv[nb] = wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_Yt[nb],M_Wt[nb],M_St[nb],M_Zt[nb],Y_ajdk,lf_AtTYn[nb],lf_ZtSWn[nb],&(lf_AtTYn____WtSZn_vv[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_AtTAn____ZtSZn_vv[nb] = wrap_AtTYn____WtSZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_At[nb],M_Zt[nb],M_St[nb],M_Zt[nb],A_ajdk,lf_AtTAn[nb],lf_ZtSZn[nb],&(lf_AtTAn____ZtSZn_vv[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYn____WtSZn_vv: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  /* if (error_check){ } */}
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
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_An[nb],M_Zn[nb],A_ajdk,lf_An_ajdk[nb],lf_Zn_ajdk[nb],&(lf_AnZt[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_An[nb],M_An[nb],A_ajdk,lf_An_ajdk[nb],lf_An_ajdk[nb],&(lf_AnAt[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Yn[nb],M_Wn[nb],Y_ajdk,lf_Yn_ajdk[nb],lf_Wn_ajdk[nb],&(lf_YnWt[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic(); 
	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Yn[nb],M_Yn[nb],Y_ajdk,lf_Yn_ajdk[nb],lf_Yn_ajdk[nb],&(lf_YnYt[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZt and YnWt: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
	  wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Yn[nb]->mr_j,M_Yn[nb]->mr_b,M_Wn[nb]->mr_j,M_Wn[nb]->mr_b,lf_YnWt[nb],lf_YnWt[nb]->lf,&(M_Yn[nb]->nrows),&(M_Wn[nb]->nrows),&(M_YnWt[nb]),&(GLOBAL_B_MLT),(addressable_0));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_An[nb]->mr_j,M_An[nb]->mr_b,M_Zn[nb]->mr_j,M_Zn[nb]->mr_b,lf_AnZt[nb],lf_AnZt[nb]->lf,&(M_An[nb]->nrows),&(M_Zn[nb]->nrows),&(M_AnZt[nb]),&(GLOBAL_B_MLT),(addressable_0));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
	  wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Wn[nb]->mr_j,M_Wn[nb]->mr_b,M_Yn[nb]->mr_j,M_Yn[nb]->mr_b,lf_YnWt[nb],lf_YnWt[nb]->lf,&(M_Wn[nb]->nrows),&(M_Yn[nb]->nrows),&(M_WnYt[nb]),&(GLOBAL_B_MLT),(addressable_1));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Zn[nb]->mr_j,M_Zn[nb]->mr_b,M_An[nb]->mr_j,M_An[nb]->mr_b,lf_AnZt[nb],lf_AnZt[nb]->lf,&(M_Zn[nb]->nrows),&(M_An[nb]->nrows),&(M_ZnAt[nb]),&(GLOBAL_B_MLT),(addressable_1));
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
	  length_At_T_YnWt_S_Zn_ww[nb] = wrap_At_T_YnWt_S_Zn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_Yn[nb],M_Wn[nb],M_St[nb],M_Zt[nb],A_ajdk,M_YnWt[nb],&(lf_At_T_YnWt_S_Zn_ww[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  length_At_T_AnZt_S_Zn_ww[nb] = wrap_At_T_YnWt_S_Zn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_An[nb],M_Zn[nb],M_St[nb],M_Zt[nb],A_ajdk,M_AnZt[nb],&(lf_At_T_AnZt_S_Zn_ww[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% At_T_YnWt_S_Zn_ww: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic();
	  length_Zt_S_WnYt_T_An_ww[nb] = wrap_At_T_YnWt_S_Zn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_Zt[nb],M_St[nb],M_Wn[nb],M_Yn[nb],M_Tt[nb],M_At[nb],A_ajdk,M_WnYt[nb],&(lf_Zt_S_WnYt_T_An_ww[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  length_Zt_S_ZnAt_T_An_ww[nb] = wrap_At_T_YnWt_S_Zn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_Zt[nb],M_St[nb],M_Zn[nb],M_An[nb],M_Tt[nb],M_At[nb],A_ajdk,M_ZnAt[nb],&(lf_Zt_S_ZnAt_T_An_ww[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% Zt_S_WnYt_T_An_ww: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic();
	  length_At_T_YnWt_S_Zn_vv[nb] = wrap_At_T_YnWt_S_Zn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_St[nb],M_Zt[nb],A_ajdk,lf_YnWt[nb],&(lf_At_T_YnWt_S_Zn_vv[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  length_At_T_AnZt_S_Zn_vv[nb] = wrap_At_T_YnWt_S_Zn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_St[nb],M_Zt[nb],A_ajdk,lf_AnZt[nb],&(lf_At_T_AnZt_S_Zn_vv[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% At_T_YnWt_S_Zn_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/* /\*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*\/ */
	/* GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0; */
	/* GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; */
	/* for (nb=0;nb<nbins;nb++){ */
	/*   GLOBAL_pthread_tic(); */
	/*   length_AtTYnWtSZn_uu[nb] = wrap_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_Yn[nb],M_Wn[nb],M_St[nb],M_Zt[nb],A_ajdk,Y_ajdk,&(lf_AtTYnWtSZn_uu[nb])); */
	/*   GLOBAL_pthread_toc(); */
	/*   GLOBAL_pthread_tic(); */
	/*   length_AtTAnZtSZn_uu[nb] = wrap_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_An[nb],M_Zn[nb],M_St[nb],M_Zt[nb],A_ajdk,A_ajdk,&(lf_AtTAnZtSZn_uu[nb])); */
	/*   GLOBAL_pthread_toc(); */
	/*   /\* for (nb=0;nb<nbins;nb++){ } *\/} */
	/* GLOBAL_pthread_tuc();  */
	/* GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYnWtSZn_uu: "); */
	/* GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: "); */
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_AtTYn____WtSZn_vv[%d] - lf_At_T_YnWt_S_Zn_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_AtTYn____WtSZn_vv[nb]->lf,lf_At_T_YnWt_S_Zn_vv[nb]->lf,length_AtTYn____WtSZn_vv[nb],1));
	    printf(" %% lf_AtTYn____WtSZn_vv[%d] - lf_At_T_YnWt_S_Zn_ww[%d] %0.16f\n",nb,nb,dra_diff(lf_AtTYn____WtSZn_vv[nb]->lf,lf_At_T_YnWt_S_Zn_ww[nb]->lf,length_AtTYn____WtSZn_vv[nb],1));
	    printf(" %% lf_Zt_S_WnYt_T_An_ww[%d] - lf_At_T_YnWt_S_Zn_ww[%d] %0.16f\n",nb,nb,dra_diff(lf_Zt_S_WnYt_T_An_ww[nb]->lf,lf_At_T_YnWt_S_Zn_ww[nb]->lf,length_Zt_S_WnYt_T_An_ww[nb],1));
	    printf(" %% lf_AtTAn____ZtSZn_vv[%d] - lf_At_T_AnZt_S_Zn_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_AtTAn____ZtSZn_vv[nb]->lf,lf_At_T_AnZt_S_Zn_vv[nb]->lf,length_AtTAn____ZtSZn_vv[nb],1));
	    printf(" %% lf_AtTAn____ZtSZn_vv[%d] - lf_At_T_AnZt_S_Zn_ww[%d] %0.16f\n",nb,nb,dra_diff(lf_AtTAn____ZtSZn_vv[nb]->lf,lf_At_T_AnZt_S_Zn_ww[nb]->lf,length_AtTAn____ZtSZn_vv[nb],1));
	    printf(" %% lf_Zt_S_ZnAt_T_An_ww[%d] - lf_At_T_AnZt_S_Zn_ww[%d] %0.16f\n",nb,nb,dra_diff(lf_Zt_S_ZnAt_T_An_ww[nb]->lf,lf_At_T_AnZt_S_Zn_ww[nb]->lf,length_Zt_S_ZnAt_T_An_ww[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
	/* for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_At_T_YnWt_S_Zn_ww_test]\n");}
}
