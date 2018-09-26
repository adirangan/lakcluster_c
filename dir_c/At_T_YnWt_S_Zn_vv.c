#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_At_T_YnWt_S_Zn_vv(void *vp)
{
  /* This function calculates :
     output_At_T_YnWt_S_Zn_vv[mj+ns*A_n_cols] = (At(mj,:)-a_An(mj)*e_At)*diag(T(:,ns))*(Yn(:,:)-e_An*a_Yt)*diag(D_Yn)*(Wt(:,:)-a_Yn*e_Zt)*diag(S(:,ns))*(Zn(:,mj)-e_Zn*a_At(mj));
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Zt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *lf_YnWt = (struct L_handle *)(vpra[ip++]);
  struct L_handle *output_At_T_YnWt_S_Zn_vv = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_s = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows)/* rup(M_At->nrows + M_At->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int ns_j=0,ns_b=0,ns_a=0,tab_s_stride=0,tab_s=0;
  int na_j=0,na_b=0,na_a=0,tab_a_stride=0,tab_a=0,my_j=0,my_b=0,my_a=0,nz_j=0,nz_b=0,nz_a=0,mw_j=0,mw_b=0,mw_a=0,tab_x=0;
  unsigned char *At_tag=NULL;
  unsigned char *Tt_tag=NULL;
  unsigned char *St_tag=NULL;
  unsigned char *Zt_tag=NULL;
  int vA=0,vT=0,vS=0,vZ=0;
  double output______YnWt_S_Zn_tmp=0;
  if (verbose>1){ printf(" %% [entering get_At_T_YnWt_S_Zn_vv] tidx %d\n",tidx);}  
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  if (verbose>2){ printf(" %% tab_s_stride %d tab_a_stride %d\n",tab_s_stride,tab_a_stride);}
  output_At_T_YnWt_S_Zn_vv->spacing_row = output_spacing_a; output_At_T_YnWt_S_Zn_vv->row_stride = tab_a_stride; 
  output_At_T_YnWt_S_Zn_vv->spacing_col = output_spacing_s; output_At_T_YnWt_S_Zn_vv->col_stride = tab_s_stride; 
  if (strstr(GLOBAL_skip,"At_T_YnWt_S_Zn_vv")){ goto skip_At_T_YnWt_S_Zn_vv;}
  ns_j=0;
  while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    Tt_tag = (unsigned char *)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
    St_tag = (unsigned char *)&(M_St->wX[ns_b*M_St->mc_length]);
    na_j=0; nz_j=0;
    while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
      na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
      At_tag = (unsigned char *)&(M_At->wX[na_b*M_At->mc_length]);
      Zt_tag = (unsigned char *)&(M_Zt->wX[nz_b*M_Zt->mc_length]);
      switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
      tab_x = tab_a + tab_s*tab_a_stride;
      output_At_T_YnWt_S_Zn_vv->lf[tab_x]=0;
      my_j=0;
      while (my_j<M_At->cpop_j){
	my_a = M_At->n_a_[my_j]; my_b = M_At->n_b_[my_j];
	vA = bget____(At_tag,my_a);
	vT = bget____(Tt_tag,my_a);
	output______YnWt_S_Zn_tmp = 0;
	mw_j=0;
	while (mw_j<M_Zt->cpop_j){
	  mw_a = M_Zt->n_a_[mw_j]; mw_b = M_Zt->n_b_[mw_j];
	  vS = bget____(St_tag,mw_a);
	  vZ = bget____(Zt_tag,mw_a);
	  output______YnWt_S_Zn_tmp += (*L2_get(lf_YnWt,my_j,my_b,my_a,mw_j,mw_b,mw_a)) *vS*(vZ - a_An[nz_a/POPLENGTH]);
	  mw_j++; /* while (mw_j<M_Zt->cpop_j){ } */}
	output_At_T_YnWt_S_Zn_vv->lf[tab_x] += (vA - a_An[na_a/POPLENGTH])*vT*output______YnWt_S_Zn_tmp;
	my_j++; /* while (my_j<M_At->cpop_j){ } */}
      na_j++;nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_At->rpop_j*M_At->cpop_j*M_Zt->cpop_j*2,0);
  if (verbose>1){ raprintf(output_At_T_YnWt_S_Zn_vv->lf,"double",tab_a_stride,tab_s_stride," %% output_At_T_YnWt_S_Zn_vv->lf: ");}
 skip_At_T_YnWt_S_Zn_vv:
  if (verbose>1){ printf(" %% [finished get_At_T_YnWt_S_Zn_vv] tidx %d\n",tidx);}
  return NULL;
}

int wrap_At_T_YnWt_S_Zn_vv__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,struct L_handle *lf_YnWt,struct L_handle **output_p)
{
  /* This function calls get_At_T_YnWt_S_Zn_vv ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 15)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_At_T_YnWt_S_Zn_vv__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Zt,verbose," %% M_Zt: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_a,length_s,length);}
  length = length_a*length_s; if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_At_T_YnWt_S_Zn_vv__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_St; vpra[ip++] = M_Zt; vpra[ip++] = A_ajdk; vpra[ip++] = lf_YnWt; vpra[ip++] = *output_p;
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_At_T_YnWt_S_Zn_vv,vpra)){ printf("Warning! cannot create thread %d in wrap_At_T_YnWt_S_Zn_vv__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_At_T_YnWt_S_Zn_vv(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p),"double",length_a,length_s," %% (*output_p): ");}
  if (verbose){ printf(" %% [finished wrap_At_T_YnWt_S_Zn_vv__run] tidx %d\n",*tidx);}
  return length;
}  

void *get_AtTYnWtSZn_uu(void *vp)
{
  /* This function calculates :
     output_AtTYnWtSZn_uu[mj+ns*A_n_cols] = (At(mj,:)-a_An(mj)*e_At)*diag(T(:,ns))*(Yn(:,:)-e_An*a_Yt)*diag(D_Yn)*(Wt(:,:)-a_Yn*e_Zt)*diag(S(:,ns))*(Zn(:,mj)-e_Zn*a_At(mj));
  */
  int verbose=0;
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
  double *Y_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_AtTYnWtSZn_uu = (struct L_handle *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_s = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int Y_pcols = psize(M_Yn->ncols);
  double *D_Yn = (double *)&(Y_ajdk[0 + AJDK_0_1*Y_pcols]);
  double *a_Yn = (double *)&(Y_ajdk[0 + AJDK_1_0*Y_pcols]);
  int ns_j=0,ns_b=0,ns_a=0,tab_s_stride=0,tab_s=0;
  int na_j=0,na_b=0,na_a=0,tab_a_stride=0,tab_a=0,my_j=0,my_b=0,my_a=0,nz_j=0,nz_b=0,nz_a=0,mw_j=0,mw_b=0,mw_a=0,tab_x=0;
  int ny_j=0,ny_b=0,ny_a=0,nw_j=0,nw_b=0,nw_a=0;
  unsigned char *At_tag=NULL;
  unsigned char *Tt_tag=NULL;
  unsigned char *Yn_tag=NULL;
  unsigned char *Wn_tag=NULL;
  unsigned char *St_tag=NULL;
  unsigned char *Zt_tag=NULL;
  int vA=0,vT=0,vY=0,vW=0,vS=0,vZ=0;
  double output______YnWt_S_Zn_tmp=0,output______YnWt______tmp=0;
  if (verbose>1){ printf(" %% [entering get_AtTYnWtSZn_uu] tidx %d\n",tidx);}  
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  if (verbose>2){ bprintf(M_Wn->mr_b,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_b: "); bprintf(M_Wn->mr_j,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_j: ");}
  if (verbose>2){ bprintf(M_Wn->mc_b,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_b: "); bprintf(M_Wn->mc_j,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  if (verbose>2){ printf(" %% tab_s_stride %d tab_a_stride %d\n",tab_s_stride,tab_a_stride);}
  output_AtTYnWtSZn_uu->spacing_row = output_spacing_a; output_AtTYnWtSZn_uu->row_stride = tab_a_stride; 
  output_AtTYnWtSZn_uu->spacing_col = output_spacing_s; output_AtTYnWtSZn_uu->col_stride = tab_s_stride; 
  ns_j=0;
  while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    Tt_tag = (unsigned char *)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
    St_tag = (unsigned char *)&(M_St->wX[ns_b*M_St->mc_length]);
    na_j=0; nz_j=0;
    while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){
      na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j]; nz_a = M_Zt->m_a_[nz_j]; nz_b = M_Zt->m_b_[nz_j];
      At_tag = (unsigned char *)&(M_At->wX[na_b*M_At->mc_length]);
      Zt_tag = (unsigned char *)&(M_Zt->wX[nz_b*M_Zt->mc_length]);
      switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
      tab_x = tab_a + tab_s*tab_a_stride;
      output_AtTYnWtSZn_uu->lf[tab_x]=0;
      my_j=0;
      while (my_j<M_Yn->rpop_j){
	my_a = M_Yn->m_a_[my_j]; my_b = M_Yn->m_b_[my_j];
	Yn_tag = (unsigned char *)&(M_Yn->wX[my_b*M_Yn->mc_length]);
	vA = bget____(At_tag,my_a);
	vT = bget____(Tt_tag,my_a);
	output______YnWt_S_Zn_tmp = 0;
	mw_j=0;
	while (mw_j<M_Wn->rpop_j){
	  mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	  Wn_tag = (unsigned char *)&(M_Wn->wX[mw_b*M_Wn->mc_length]);
	  vS = bget____(St_tag,mw_a);
	  vZ = bget____(Zt_tag,mw_a);
	  output______YnWt______tmp = 0;
	  ny_j=0; nw_j=0;
	  while (ny_j<M_Yn->cpop_j && nw_j<M_Wn->cpop_j){
	    ny_a = M_Yn->n_a_[ny_j]; ny_b = M_Yn->n_b_[ny_j];
	    nw_a = M_Wn->n_a_[nw_j]; nw_b = M_Wn->n_b_[nw_j];
	    vY = bget____(Yn_tag,ny_a);
	    vW = bget____(Wn_tag,nw_a);
	    output______YnWt______tmp += (vY - a_Yn[ny_a/POPLENGTH]) * D_Yn[ny_a/POPLENGTH] * (vW - a_Yn[nw_a/POPLENGTH]);
	    ny_j++; nw_j++; /* while (ny_j<M_Yn->cpop_j && nw_j<M_Wn->cpop_j){ } */}
	  output______YnWt_S_Zn_tmp += output______YnWt______tmp * vS * (vZ - a_An[nz_a/POPLENGTH]);
	  mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
	output_AtTYnWtSZn_uu->lf[tab_x] += (vA - a_An[na_a/POPLENGTH]) * vT * output______YnWt_S_Zn_tmp;
	my_j++; /* while (my_j<M_Yn->rpop_j){ } */}
      na_j++;nz_j++; /* while (na_j<M_At->rpop_j && nz_j<M_Zt->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_At->rpop_j*M_Yn->rpop_j*M_Wn->rpop_j*M_Yn->cpop_j*2,0);
  if (verbose>1){ raprintf(output_AtTYnWtSZn_uu->lf,"double",tab_a_stride,tab_s_stride," %% output_AtTYnWtSZn_uu->lf: ");}
  if (verbose>1){ printf(" %% [finished get_AtTYnWtSZn_uu] tidx %d\n",tidx);}
  return NULL;
}

int wrap_AtTYnWtSZn_uu__run(int *tidx,void **vpra,pthread_t *thread_in,int output_spacing_a,int output_spacing_s,struct M_handle *M_At,struct M_handle *M_Tt,struct M_handle *M_Yn,struct M_handle *M_Wn,struct M_handle *M_St,struct M_handle *M_Zt,double *A_ajdk,double *Y_ajdk,struct L_handle **output_p)
{
  /* This function calls get_AtTYnWtSZn_uu ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 15)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_a=0,length_s=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_AtTYnWtSZn_uu__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Yn,verbose," %% M_Yn: ");}
  if (verbose){ M_handle_printf(M_Wn,verbose," %% M_Wn: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Zt,verbose," %% M_Zt: ");}
  switch (output_spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  length = length_a*length_s; if (verbose){ printf(" %% length %llu*%llu=%llu\n",length_a,length_s,length);}
  length = length_a*length_s; if (*output_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: "); bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: "); bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Tt->mr_b,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_b: "); bprintf(M_Tt->mr_j,M_Tt->bitj,1,M_Tt->nrows," %% M_Tt->mr_j: ");}
  if (verbose>2){ bprintf(M_Tt->mc_b,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_b: "); bprintf(M_Tt->mc_j,M_Tt->bitj,1,M_Tt->ncols," %% M_Tt->mc_j: ");}
  if (verbose>2){ bprintf(M_Yn->mr_b,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_b: "); bprintf(M_Yn->mr_j,M_Yn->bitj,1,M_Yn->nrows," %% M_Yn->mr_j: ");}
  if (verbose>2){ bprintf(M_Yn->mc_b,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_b: "); bprintf(M_Yn->mc_j,M_Yn->bitj,1,M_Yn->ncols," %% M_Yn->mc_j: ");}
  if (verbose>2){ bprintf(M_Wn->mr_b,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_b: "); bprintf(M_Wn->mr_j,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_j: ");}
  if (verbose>2){ bprintf(M_Wn->mc_b,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_b: "); bprintf(M_Wn->mc_j,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_j: ");}
  if (verbose>2){ bprintf(M_St->mr_b,M_St->bitj,1,M_St->nrows," %% M_St->mr_b: "); bprintf(M_St->mr_j,M_St->bitj,1,M_St->nrows," %% M_St->mr_j: ");}
  if (verbose>2){ bprintf(M_St->mc_b,M_St->bitj,1,M_St->ncols," %% M_St->mc_b: "); bprintf(M_St->mc_j,M_St->bitj,1,M_St->ncols," %% M_St->mc_j: ");}
  if (verbose>2){ bprintf(M_Zt->mr_b,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_b: "); bprintf(M_Zt->mr_j,M_Zt->bitj,1,M_Zt->nrows," %% M_Zt->mr_j: ");}
  if (verbose>2){ bprintf(M_Zt->mc_b,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_b: "); bprintf(M_Zt->mc_j,M_Zt->bitj,1,M_Zt->ncols," %% M_Zt->mc_j: ");}
  length = length_a*length_s; if ((*output_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_AtTYnWtSZn_uu__run\n",(*output_p)->length,length);} memset((*output_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = M_Yn; vpra[ip++] = M_Wn; vpra[ip++] = M_St; vpra[ip++] = M_Zt; vpra[ip++] = A_ajdk; vpra[ip++] = Y_ajdk; vpra[ip++] = *output_p;
  switch (output_spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_s){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_AtTYnWtSZn_uu,vpra)){ printf("Warning! cannot create thread %d in wrap_AtTYnWtSZn_uu__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_AtTYnWtSZn_uu(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_p),"double",length_a,length_s," %% (*output_p): ");}
  if (verbose){ printf(" %% [finished wrap_AtTYnWtSZn_uu__run] tidx %d\n",*tidx);}
  return length;
}  

void wrap_At_T_YnWt_S_Zn_vv_test()
{
  /* test for errors with input file: At_T_YnWt_S_Zn_vv_error.in ;
  */
  /* test for speed with input file: At_T_YnWt_S_Zn_vv_speed.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AnZt=NULL,**lf_AnAt=NULL,**lf_YnWt=NULL,**lf_YnYt=NULL;
  struct L_handle **lf_At_T_YnWt_S_Zn_vv=NULL,**lf_At_T_AnZt_S_Zn_vv=NULL; int *length_At_T_YnWt_S_Zn_vv=NULL,*length_At_T_AnZt_S_Zn_vv=NULL;
  struct L_handle **lf_AtTYnWtSZn_uu=NULL,**lf_AtTAnZtSZn_uu=NULL; int *length_AtTYnWtSZn_uu=NULL,*length_AtTAnZtSZn_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_At_T_YnWt_S_Zn_vv_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
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
  if (error_check){ 
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
	if (error_check){
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
	    printf(" %% lf_At_T_YnWt_S_Zn_vv[%d] error %0.16f\n",nb,dra_diff(lf_At_T_YnWt_S_Zn_vv[nb]->lf,lf_AtTYnWtSZn_uu[nb]->lf,length_At_T_YnWt_S_Zn_vv[nb],1));
	    printf(" %% lf_At_T_AnZt_S_Zn_vv[%d] error %0.16f\n",nb,dra_diff(lf_At_T_AnZt_S_Zn_vv[nb]->lf,lf_AtTAnZtSZn_uu[nb]->lf,length_At_T_AnZt_S_Zn_vv[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
	/* for (n_type=1;n_type<=1;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_At_T_YnWt_S_Zn_vv_test]\n");}
}
