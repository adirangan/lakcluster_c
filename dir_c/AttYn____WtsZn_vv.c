#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_AttYn____WtsZn_vv(void *vp)
{
  /* Given a particular ns_j, This function calculates :
     output_AttYn____WtsZn_vv[mj+ns*A_n_cols] = (At(mj,:)-a_An(mj)*e_At)*diag(T(:,ns))*(Yn(:,:)-e_An*a_Yt)*diag(D_Yn)*(Wt(:,:)-a_Yn*e_Zt)*diag(S(:,ns))*(Zn(:,mj)-e_Zn*a_At(mj)) = ;
     + At(mj,:)*diag(T(:,ns))*  Yn(:,:)  *diag(D_Yn)*  Wt(:,:)  *diag(S(:,ns))*Zn(:,mj) ;
     - At(mj,:)*diag(T(:,ns))* e_An*a_Yt *diag(D_Yn)*  Wt(:,:)  *diag(S(:,ns))*Zn(:,mj) ;
     - At(mj,:)*diag(T(:,ns))*  Yn(:,:)  *diag(D_Yn)* a_Yn*e_Zt *diag(S(:,ns))*Zn(:,mj) ;
     + At(mj,:)*diag(T(:,ns))* e_An*a_Yt *diag(D_Yn)* a_Yn*e_Zt *diag(S(:,ns))*Zn(:,mj) ;
     - a_An(mj)*diag(T(:,ns))*  Yn(:,:)  *diag(D_Yn)*  Wt(:,:)  *diag(S(:,ns))*Zn(:,mj) ;
     + a_An(mj)*diag(T(:,ns))* e_An*a_Yt *diag(D_Yn)*  Wt(:,:)  *diag(S(:,ns))*Zn(:,mj) ;
     + a_An(mj)*diag(T(:,ns))*  Yn(:,:)  *diag(D_Yn)* a_Yn*e_Zt *diag(S(:,ns))*Zn(:,mj) ;
     - a_An(mj)*diag(T(:,ns))* e_An*a_Yt *diag(D_Yn)* a_Yn*e_Zt *diag(S(:,ns))*Zn(:,mj) ;
     + At(mj,:)*diag(T(:,ns))*  Yn(:,:)  *diag(D_Yn)*  Wt(:,:)  *diag(S(:,ns))*a_At(mj) ;
     - At(mj,:)*diag(T(:,ns))* e_An*a_Yt *diag(D_Yn)*  Wt(:,:)  *diag(S(:,ns))*a_At(mj) ;
     - At(mj,:)*diag(T(:,ns))*  Yn(:,:)  *diag(D_Yn)* a_Yn*e_Zt *diag(S(:,ns))*a_At(mj) ;
     + At(mj,:)*diag(T(:,ns))* e_An*a_Yt *diag(D_Yn)* a_Yn*e_Zt *diag(S(:,ns))*a_At(mj) ;
     - a_An(mj)*diag(T(:,ns))*  Yn(:,:)  *diag(D_Yn)*  Wt(:,:)  *diag(S(:,ns))*a_At(mj) ;
     + a_An(mj)*diag(T(:,ns))* e_An*a_Yt *diag(D_Yn)*  Wt(:,:)  *diag(S(:,ns))*a_At(mj) ;
     + a_An(mj)*diag(T(:,ns))*  Yn(:,:)  *diag(D_Yn)* a_Yn*e_Zt *diag(S(:,ns))*a_At(mj) ;
     - a_An(mj)*diag(T(:,ns))* e_An*a_Yt *diag(D_Yn)* a_Yn*e_Zt *diag(S(:,ns))*a_At(mj) ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  int ns_j = *(int *)(vpra[ip++]);
  struct M_handle *M_Yt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Wt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *lf_AttYn = (struct L_handle *)(vpra[ip++]);
  struct L_handle *lf_ZtsWn = (struct L_handle *)(vpra[ip++]);
  struct L_handle *output_AttYn____WtsZn_vv = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_s = *(int *)(vpra[ip++]);
  int Y_pcols = psize(M_Yt->nrows)/* rup(M_Yt->nrows + M_Yt->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int ns_b=0,ns_a=0,tab_s_stride=0,tab_s=0;
  int na_j=0,na_b=0,na_a=0,tab_a_stride=0,tab_a=0,ny_j=0,ny_b=0,ny_a=0,nz_j=0,nz_b=0,nz_a=0,nw_j=0,nw_b=0,nw_a=0,tab_x=0;
  int mx_j=0,mx_chunk=0;
  if (verbose>1){ printf(" %% [entering get_AttYn____WtsZn_vv] tidx %d\n",tidx);}  
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Yt->mr_b,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_b: "); bprintf(M_Yt->mr_j,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_j: ");}
  if (verbose>2){ bprintf(M_Yt->mc_b,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_b: "); bprintf(M_Yt->mc_j,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_j: ");}
  if (verbose>2){ bprintf(M_Wt->mr_b,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_b: "); bprintf(M_Wt->mr_j,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_j: ");}
  if (verbose>2){ bprintf(M_Wt->mc_b,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_b: "); bprintf(M_Wt->mc_j,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  if (verbose>2){ raprintf(D_Yn,"double",1,Y_pcols," %% D_Yn: "); raprintf(a_Yn,"double",1,Y_pcols," %% a_Yn: ");}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  if (verbose>2){ printf(" %% tab_s_stride %d tab_a_stride %d\n",tab_s_stride,tab_a_stride);}
  output_AttYn____WtsZn_vv->spacing_row = output_spacing_a; output_AttYn____WtsZn_vv->row_stride = tab_a_stride; 
  output_AttYn____WtsZn_vv->spacing_col = output_spacing_s; output_AttYn____WtsZn_vv->col_stride = tab_s_stride; 
  if (strstr(GLOBAL_skip,"AttYn____WtsZn_vv")){ goto skip_AttYn____WtsZn_vv;}
  if (GLOBAL_omp_type==GLOBAL_omp_off){ 
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    na_j=0; nz_j=0;
    while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
      na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
      tab_x = tab_a + tab_s*tab_a_stride;
      output_AttYn____WtsZn_vv->lf[tab_x]=0;
      ny_j=0; nw_j=0;
      while (ny_j<M_Yt->rpop_j && nw_j<M_Wt->rpop_j){
	ny_a = M_Yt->m_a_[ny_j]; ny_b = M_Yt->m_b_[ny_j];
	nw_a = M_Wt->m_a_[nw_j]; nw_b = M_Wt->m_b_[nw_j];
	output_AttYn____WtsZn_vv->lf[tab_x] += (*L2_get(lf_AttYn,na_j,na_b,na_a,ny_j,ny_b,ny_a)) *D_Yn[ny_a/POPLENGTH]* (*L2_get(lf_ZtsWn,nz_j,nz_b,nz_a,nw_j,nw_b,nw_a));
	ny_j++;nw_j++; /* while (ny_j<M_Yt->rpop_j && nw_j<M_Wt->rpop_j){ } */}
      na_j++;nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp__on){
    mx_chunk=128; 
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
#pragma omp parallel shared(mx_chunk) private(mx_j,na_j,nz_j,na_a,na_b,nz_a,nz_b,tab_a,tab_x,ny_j,nw_j,ny_a,ny_b,nw_a,nw_b)
    { /* begin omp parallel */
      mx_j=0;
#pragma omp for schedule(dynamic,mx_chunk)
      for (mx_j=0;mx_j<M_At->rpop_j;mx_j++){
	na_j = mx_j; nz_j = na_j;
	na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
	switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
	tab_x = tab_a + tab_s*tab_a_stride;
	output_AttYn____WtsZn_vv->lf[tab_x]=0;
	ny_j=0; nw_j=0;
	while (ny_j<M_Yt->rpop_j && nw_j<M_Wt->rpop_j){
	  ny_a = M_Yt->m_a_[ny_j]; ny_b = M_Yt->m_b_[ny_j];
	  nw_a = M_Wt->m_a_[nw_j]; nw_b = M_Wt->m_b_[nw_j];
	  output_AttYn____WtsZn_vv->lf[tab_x] += (*L2_get(lf_AttYn,na_j,na_b,na_a,ny_j,ny_b,ny_a)) *D_Yn[ny_a/POPLENGTH]* (*L2_get(lf_ZtsWn,nz_j,nz_b,nz_a,nw_j,nw_b,nw_a));
	  ny_j++;nw_j++; /* while (ny_j<M_Yt->rpop_j && nw_j<M_Wt->rpop_j){ } */}
	/* for (mx_j=0;mx_j<M_At->rpop_j;mx_j++){ } */}	
      /* end omp parallel */}
    /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
  GLOBAL_ops_count_one(tidx,M_At->rpop_j*M_Yt->rpop_j,0);
 skip_AttYn____WtsZn_vv:
  if (verbose>1){ printf(" %% [finished get_AttYn____WtsZn_vv] tidx %d\n",tidx);}
  return NULL;
}

int wrap_AttYn____WtsZn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,int ns_j,struct M_handle *M_Yt,struct M_handle *M_Wt,struct M_handle *M_St,struct M_handle *M_Zt,double *Y_ajdk,struct L_handle *lf_AttYn,struct L_handle *lf_ZtsWn,struct L_handle **output_p)
{
  /* This function calls get_AttYn____WtsZn_vv ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 16)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_a=0,length_s=0,length=0,ip=0;
  int ns_a=0,ns_b=0,tab_s=0;
  if (verbose){ printf(" %% [entering wrap_AttYn____WtsZn_vv__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Yt,verbose," %% M_Yt: ");}
  if (verbose){ M_handle_printf(M_Wt,verbose," %% M_Wt: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Zt,verbose," %% M_Zt: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  ns_a = M_Tt->m_a_[ns_j]; ns_b = M_Tt->m_b_[ns_j];
  switch (output_spacing_s){ case SPACING_j: tab_s = ns_j; break; case SPACING_b: tab_s = ns_b; break; case SPACING_a: tab_s = ns_a; break; default: break; /* switch (output_spacing_s){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_a,length_s,length);}
  length = length_a*length_s; if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Yt->mr_b,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_b: "); bprintf(M_Yt->mr_j,M_Yt->bitj,1,M_Yt->nrows," %% M_Yt->mr_j: ");}
  if (verbose>2){ bprintf(M_Yt->mc_b,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_b: "); bprintf(M_Yt->mc_j,M_Yt->bitj,1,M_Yt->ncols," %% M_Yt->mc_j: ");}
  if (verbose>2){ bprintf(M_Wt->mr_b,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_b: "); bprintf(M_Wt->mr_j,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_j: ");}
  if (verbose>2){ bprintf(M_Wt->mc_b,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_b: "); bprintf(M_Wt->mc_j,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_AttYn____WtsZn_vv__run\n",(*output_p)->length,length);} 
  memset(&((*output_p)->lf[tab_s*length_a]),0,length_a*sizeof(double));
  if (ns_j<0 || ns_j>addressable_int_length-1){ printf(" %% Warning! ns_j %d in wrap_AttYn____WtsZn_vv_run\n",ns_j);}
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = &(addressable_int[ns_j]); vpra[ip++] = M_Yt; vpra[ip++] = M_Wt; vpra[ip++] = M_St; vpra[ip++] = M_Zt; vpra[ip++] = Y_ajdk; vpra[ip++] = lf_AttYn; vpra[ip++] = lf_ZtsWn; vpra[ip++] = *output_p;
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AttYn____WtsZn_vv,vpra)){ printf("Warning! cannot create thread %d in wrap_AttYn____WtsZn_vv__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AttYn____WtsZn_vv(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p),"double",length_a,length_s," %% (*output_p): ");}
  if (verbose){ printf(" %% [finished wrap_AttYn____WtsZn_vv__run] tidx %d\n",*tidx);}
  return length;
}  

void wrap_AttYn____WtsZn_vv_test()
{
  /* test for errors with input file: AttYn____WtsZn_vv_error.in ;
  */
  /* test for speed with input file: AttYn____WtsZn_vv_speed.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  int ns_j=0,ns_a=0,ns_b=0;
  struct L_handle **lf_AnZt=NULL,**lf_AnAt=NULL,**lf_YnWt=NULL,**lf_YnYt=NULL;
  struct L_handle **lf_AtTn=NULL,**lf_YtTn=NULL,**lf_ZtSn=NULL,**lf_WtSn=NULL;
  struct L_handle **lf_ZtsWn=NULL,**lf_AttYn=NULL,**lf_ZtsZn=NULL,**lf_AttAn=NULL;
  struct L_handle **lf_AttYn____WtsZn_vv=NULL,**lf_AttAn____ZtsZn_vv=NULL;
  struct L_handle **lf_At_T_YnWt_S_Zn_vv=NULL,**lf_At_T_AnZt_S_Zn_vv=NULL; int *length_At_T_YnWt_S_Zn_vv=NULL,*length_At_T_AnZt_S_Zn_vv=NULL;
  struct L_handle **lf_AtTYnWtSZn_uu=NULL,**lf_AtTAnZtSZn_uu=NULL; int *length_AtTYnWtSZn_uu=NULL,*length_AtTAnZtSZn_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_AttYn____WtsZn_vv_test]\n");}
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
  lf_ZtsWn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AttYn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_ZtsZn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AttAn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_ZtsWn[nb] = L_handle_make((unsigned long long int)M_Zt[nb]->nrows*(unsigned long long int)M_Wt[nb]->nrows);
    lf_AttYn[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Yt[nb]->nrows);
    lf_ZtsZn[nb] = L_handle_make((unsigned long long int)M_Zt[nb]->nrows*(unsigned long long int)M_Zt[nb]->nrows);
    lf_AttAn[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_At[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_AttYn____WtsZn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AttAn____ZtsZn_vv = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    lf_AttYn____WtsZn_vv[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Sn[nb]->ncols);
    lf_AttAn____ZtsZn_vv[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Tn[nb]->ncols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (error_check){
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
    /* if (error_check){ } */}
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
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  	for (nb=0;nb<nbins;nb++){
	  ns_j=0;
	  while (ns_j<M_Tt[nb]->rpop_j){
	    ns_a = M_Tt[nb]->m_a_[ns_j]; ns_b = M_Tt[nb]->m_b_[ns_j];
	    printf(" %% nb %d; ns %d,%d,%d\n",nb,ns_j,ns_b,ns_a);
	    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; GLOBAL_pthread_tic();
	    wrap_AttYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],ns_j,M_Yt[nb],A_ajdk,Y_ajdk,lf_AtTn[nb],lf_YtTn[nb],&(lf_AttYn[nb]));
	    GLOBAL_pthread_toc(); GLOBAL_pthread_tuc();
	    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; GLOBAL_pthread_tic();
	    wrap_AttAn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],ns_j,M_At[nb],A_ajdk,A_ajdk,lf_AtTn[nb],lf_AtTn[nb],&(lf_AttAn[nb]));
	    GLOBAL_pthread_toc(); GLOBAL_pthread_tuc();
	    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; GLOBAL_pthread_tic();
	    wrap_AttYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],ns_j,M_Wt[nb],A_ajdk,Y_ajdk,lf_ZtSn[nb],lf_WtSn[nb],&(lf_ZtsWn[nb]));
	    GLOBAL_pthread_toc(); GLOBAL_pthread_tuc();
	    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; GLOBAL_pthread_tic();
	    wrap_AttAn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],ns_j,M_Zt[nb],A_ajdk,A_ajdk,lf_ZtSn[nb],lf_ZtSn[nb],&(lf_ZtsZn[nb]));
	    GLOBAL_pthread_toc(); GLOBAL_pthread_tuc();	    
	    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; GLOBAL_pthread_tic();
	    wrap_AttYn____WtsZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],ns_j,M_Yt[nb],M_Wt[nb],M_St[nb],M_Zt[nb],Y_ajdk,lf_AttYn[nb],lf_ZtsWn[nb],&(lf_AttYn____WtsZn_vv[nb]));
	    GLOBAL_pthread_toc(); GLOBAL_pthread_tuc();	    
	    GLOBAL_nf_cur=0; GLOBAL_nf_opn=0; GLOBAL_pthread_tic();
	    wrap_AttYn____WtsZn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],ns_j,M_At[nb],M_Zt[nb],M_St[nb],M_Zt[nb],A_ajdk,lf_AttAn[nb],lf_ZtsZn[nb],&(lf_AttAn____ZtsZn_vv[nb]));
	    GLOBAL_pthread_toc(); GLOBAL_pthread_tuc();	    
	    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	    ns_j++; /* while (ns_j<M_Tt[nb]->rpop_j){ } */}
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AttYn and ZtsWn and AttYn____WtsZn_vv and AttAn____ZtsZn_vv: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
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
	    length_At_T_YnWt_S_Zn_vv[nb] = wrap_At_T_YnWt_S_Zn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_St[nb],M_Zt[nb],A_ajdk,lf_YnWt[nb],&(lf_At_T_YnWt_S_Zn_vv[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_At_T_AnZt_S_Zn_vv[nb] = wrap_At_T_YnWt_S_Zn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_St[nb],M_Zt[nb],A_ajdk,lf_AnZt[nb],&(lf_At_T_AnZt_S_Zn_vv[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% At_T_YnWt_S_Zn_vv: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic();
	    length_AtTYnWtSZn_uu[nb] = wrap_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_Yn[nb],M_Wn[nb],M_St[nb],M_Zt[nb],A_ajdk,Y_ajdk,&(lf_AtTYnWtSZn_uu[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_AtTAnZtSZn_uu[nb] = wrap_AtTYnWtSZn_uu__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],M_An[nb],M_Zn[nb],M_St[nb],M_Zt[nb],A_ajdk,A_ajdk,&(lf_AtTAnZtSZn_uu[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc(); 
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AtTYnWtSZn_uu: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_At_T_YnWt_S_Zn_vv[%d] -     lf_AtTYnWtSZn_uu[%d] %0.16f\n",nb,nb,dra_diff(lf_At_T_YnWt_S_Zn_vv[nb]->lf,lf_AtTYnWtSZn_uu[nb]->lf,length_At_T_YnWt_S_Zn_vv[nb],1));
	    printf(" %% lf_AttYn____WtsZn_vv[%d] -     lf_AtTYnWtSZn_uu[%d] %0.16f\n",nb,nb,dra_diff(lf_AttYn____WtsZn_vv[nb]->lf,lf_AtTYnWtSZn_uu[nb]->lf,length_At_T_YnWt_S_Zn_vv[nb],1));
	    printf(" %% lf_AttYn____WtsZn_vv[%d] - lf_At_T_YnWt_S_Zn_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_AttYn____WtsZn_vv[nb]->lf,lf_At_T_YnWt_S_Zn_vv[nb]->lf,length_At_T_YnWt_S_Zn_vv[nb],1));
	    printf(" %% lf_At_T_AnZt_S_Zn_vv[%d] -     lf_AtTAnZtSZn_uu[%d] %0.16f\n",nb,nb,dra_diff(lf_At_T_AnZt_S_Zn_vv[nb]->lf,lf_AtTAnZtSZn_uu[nb]->lf,length_At_T_AnZt_S_Zn_vv[nb],1));
	    printf(" %% lf_AttAn____ZtsZn_vv[%d] -     lf_AtTAnZtSZn_uu[%d] %0.16f\n",nb,nb,dra_diff(lf_AttAn____ZtsZn_vv[nb]->lf,lf_AtTAnZtSZn_uu[nb]->lf,length_At_T_AnZt_S_Zn_vv[nb],1));
	    printf(" %% lf_AttAn____ZtsZn_vv[%d] - lf_At_T_AnZt_S_Zn_vv[%d] %0.16f\n",nb,nb,dra_diff(lf_AttAn____ZtsZn_vv[nb]->lf,lf_At_T_AnZt_S_Zn_vv[nb]->lf,length_At_T_AnZt_S_Zn_vv[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
	/* for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_AttYn____WtsZn_vv_test]\n");}
}
