#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_An_v__run(void *vp)
{
  /* This function serves as a module which calculates the matrix-vector product An*v, where v=ones(size(An,2),1);
     Here An is treated as a binary representation of a matrix with 
     type_flag==TYPE_p_: An taking on +1/-1 values ;
     type_flag==TYPE_0_: An taking on +1/+0 values ;
     This function takes into account M_A->mr_b and M_A->mr_j, as well as M_A->mc_b and M_A->mc_j, skipping rows and columns if necessary ; 
     We set up ma_b and ma_j to count from the beginning of M_A->mr_b and M_A->mr_j ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct L_handle *output_An_v = (struct L_handle *)(vpra[ip++]);
  int type_flag = *(int *)(vpra[ip++]); 
  int output_spacing_r = *(int *)(vpra[ip++]);
  int type_flag_single=0;
  int mx_chunk=0,ma_a=0,ma_b=0,ma_j=0,nc_a=0,tab_r=0,tab_r_stride=0;
  double dtmp=0;
  __m128i *wA_tag;
  __m128i *mc_tag,*mc_end;
  unsigned char *A_tag=NULL;
  char tempchar[FNAMESIZE];
  if (type_flag==TYPE_p0){ printf(" %% Warning! TYPE_p0 not implemented within get_An_v__run\n");}
  if      (type_flag==TYPE_pm){ type_flag_single = TYPE_p_;}
  else if (type_flag==TYPE_00){ type_flag_single = TYPE_0_;}
  if (verbose>1){ printf(" %% [entering get_An_v__run] tidx %d rpop_b %d(%d) cpop_b %d(%d) type_flag %d\n",tidx,(int)(M_An->rpop_b),(int)(M_An->rpop_j),(int)(M_An->cpop_b),(int)(M_An->cpop_j),type_flag);}  
  if (verbose>2){ bprintf((unsigned char *)M_An->wX,POPLENGTH,M_An->rpop_b,minimum(M_An->mc_length*BIT8,M_An->ncols)," %% M_An->wX: ");}
  unsigned int *ma_b_,*ma_a_;
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  if (verbose>1){ sprintf(tempchar," %% ma_b_[%d]: ",tidx); raprintf(ma_b_,"unsigned int",1,M_An->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% ma_a_[%d]: ",tidx); raprintf(ma_a_,"unsigned int",1,M_An->rpop_j,tempchar);}
  switch (output_spacing_r){ case SPACING_j: tab_r_stride = M_An->rpop_j; break; case SPACING_b: tab_r_stride = M_An->rpop_b; break; case SPACING_a: tab_r_stride = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  if (verbose>1){ printf(" %% calculating output_An_v\n");}
  output_An_v->spacing_row = output_spacing_r; output_An_v->row_stride = tab_r_stride;
  ma_j=0;
  while (ma_j<M_An->rpop_j){
    ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
    switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
    if      (type_flag_single==TYPE_p_){ output_An_v->lf[tab_r] = -M_An->cpop_j;}
    else if (type_flag_single==TYPE_0_){ output_An_v->lf[tab_r] = 0;}
    wA_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
    mc_tag = (__m128i*)&(M_An->mc_j[0]); mc_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
    if      (type_flag_single==TYPE_p_){ dtmp = popcount(&wA_tag,&mc_tag,&mc_end);}
    else if (type_flag_single==TYPE_0_){ dtmp = popcount(&wA_tag,&mc_tag,&mc_end);}
    if      (type_flag_single==TYPE_p_){ output_An_v->lf[tab_r] += (2*dtmp);}
    else if (type_flag_single==TYPE_0_){ output_An_v->lf[tab_r] += (dtmp);}
    ma_j++;}
  GLOBAL_ops_count_one(tidx,M_An->rpop_j,M_An->rpop_j*M_An->mc_length*BIT8);
  if (verbose>1){ printf(" %% [finished get_An_v__run] tidx %d \n",tidx);}  
  return NULL;
}

int wrap_An_v__run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,struct M_handle *M_An,struct L_handle **output_An_v_p)
{
  /* This function uses the M_handle M_A to run a series of parallel calls to get_An_v__run ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 5)
     the type of multiplication is determined by type_flag: ;
     type_flag==TYPE_pm: ;
     type_flag==TYPE_00: ;
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_r=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_An_v__run] tidx %d type_flag %d\n",*tidx,type_flag);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  switch (output_spacing_r){ case SPACING_j: length_r = M_An->rpop_j; break; case SPACING_b: length_r = M_An->rpop_b; break; case SPACING_a: length_r = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  length = length_r; if (verbose){ printf(" %% length %llu\n",length_r);}
  length = length_r; if (*output_An_v_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_v_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  length = length_r*1; if ((*output_An_v_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_v__run\n",(*output_An_v_p)->length,length);} memset((*output_An_v_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = *output_An_v_p; 
  switch (type_flag){ case TYPE_p0: vpra[ip++] = &addressable_type_p0; break; case TYPE_pm: vpra[ip++] = &addressable_type_pm; break; case TYPE_00: vpra[ip++] = &addressable_type_00; break; default: break; /* switch (type_flag){ } */}
  switch (output_spacing_r){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_r){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_v__run,vpra)){ printf("Warning! cannot create thread %d in wrap_An_v__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_v__run(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_An_v_p)->lf,"double",1,length_r," %% (*output_An_v_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_An_v__run] tidx %d\n",*tidx);}
  return length;
}  

void *get_An_u__run(void *vp)
{
  /* This function serves as a module which calculates the matrix-vector product An*v, where v=ones(size(An,2),1);
     Here An is treated as a binary representation of a matrix with 
     type_flag==TYPE_p_: An taking on +1/-1 values ;
     type_flag==TYPE_0_: An taking on +1/+0 values ;
     This function takes into account M_A->mr_b and M_A->mr_j, as well as M_A->mc_b and M_A->mc_j, skipping rows and columns if necessary ; 
     We set up ma_b and ma_j to count from the beginning of M_A->mr_b and M_A->mr_j ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct L_handle *output_An_u = (struct L_handle *)(vpra[ip++]);
  int type_flag = *(int *)(vpra[ip++]); 
  int output_spacing_r = *(int *)(vpra[ip++]);
  int type_flag_single=0;
  int ma_a=0,ma_b=0,ma_j=0,nc_a=0,tab_r=0,tab_r_stride=0;
  double dtmp=0;
  unsigned char *A_tag=NULL;
  char tempchar[FNAMESIZE];
  if (type_flag==TYPE_p0){ printf(" %% Warning! TYPE_p0 not implemented within get_An_u__run\n");}
  if      (type_flag==TYPE_pm){ type_flag_single = TYPE_p_;}
  else if (type_flag==TYPE_00){ type_flag_single = TYPE_0_;}
  if (verbose>1){ printf(" %% [entering get_An_u__run] tidx %d rpop_b %d(%d) cpop_b %d(%d) type_flag %d\n",tidx,(int)(M_An->rpop_b),(int)(M_An->rpop_j),(int)(M_An->cpop_b),(int)(M_An->cpop_j),type_flag);}  
  if (verbose>2){ bprintf((unsigned char *)M_An->wX,POPLENGTH,M_An->rpop_b,minimum(M_An->mc_length*BIT8,M_An->ncols)," %% M_An->wX: ");}
  unsigned int *ma_b_,*ma_a_;
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  if (verbose>1){ sprintf(tempchar," %% ma_b_[%d]: ",tidx); raprintf(ma_b_,"unsigned int",1,M_An->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% ma_a_[%d]: ",tidx); raprintf(ma_a_,"unsigned int",1,M_An->rpop_j,tempchar);}
  switch (output_spacing_r){ case SPACING_j: tab_r_stride = M_An->rpop_j; break; case SPACING_b: tab_r_stride = M_An->rpop_b; break; case SPACING_a: tab_r_stride = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  if (verbose>1){ printf(" %% calculating output_An_u\n");}
  output_An_u->spacing_row = output_spacing_r; output_An_u->row_stride = tab_r_stride;
  memset(output_An_u->lf,0,tab_r_stride*1*sizeof(double));
  ma_j=0;
  while (ma_j<M_An->rpop_j){
    ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
    switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
    A_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length]));
    dtmp=0;
    for (nc_a=0;nc_a<M_An->ncols;nc_a++){
      if      (type_flag_single==TYPE_p_){ dtmp += bget____(A_tag,nc_a)*bget__on(M_An->mc_j,nc_a);}
      else if (type_flag_single==TYPE_0_){ dtmp += bget__on(A_tag,nc_a)*bget__on(M_An->mc_j,nc_a);}
      /* for (nc_a=0;nc_a<M_An->ncols;nc_a++){ } */}
    output_An_u->lf[tab_r] = dtmp;
    ma_j++;}
  if (verbose>1){ printf(" %% [finished get_An_u__run] tidx %d\n",tidx);}  
  return NULL;
}

int wrap_An_u__run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,struct M_handle *M_An,struct L_handle **output_An_u_p)
{
  /* This function uses the M_handle M_A to run a series of parallel calls to get_An_u__run ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 5)
     the type of multiplication is determined by type_flag: ;
     type_flag==TYPE_pm: ;
     type_flag==TYPE_00: ;
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_r=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_An_u__run] tidx %d type_flag %d\n",*tidx,type_flag);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  switch (output_spacing_r){ case SPACING_j: length_r = M_An->rpop_j; break; case SPACING_b: length_r = M_An->rpop_b; break; case SPACING_a: length_r = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  length = length_r; if (verbose){ printf(" %% length %llu\n",length_r);}
  length = length_r; if (*output_An_u_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_u_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  length = length_r*1; if ((*output_An_u_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_u__run\n",(*output_An_u_p)->length,length);} memset((*output_An_u_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = *output_An_u_p; 
  switch (type_flag){ case TYPE_p0: vpra[ip++] = &addressable_type_p0; break; case TYPE_pm: vpra[ip++] = &addressable_type_pm; break; case TYPE_00: vpra[ip++] = &addressable_type_00; break; default: break; /* switch (type_flag){ } */}
  switch (output_spacing_r){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_r){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_u__run,vpra)){ printf("Warning! cannot create thread %d in wrap_An_u__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_u__run(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_An_u_p),"double",1,length_r," %% (*output_An_u_p): ");}
  if (verbose){ printf(" %% [finished wrap_An_u__run] tidx %d\n",*tidx);}
  return length;
}

void wrap_An_v_test()
{
  /* test for errors with input file: An_v_error.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= An_v;
     GLOBAL_TEST_TYP2= error;
     GLOBAL_NBINS= 4;
     GLOBAL_TEST_mrand= 0.5;
     GLOBAL_TEST_A_n_rows= 21;
     GLOBAL_TEST_A_n_cols= 2050;
     GLOBAL_TEST_Z_n_rows= 22;
     GLOBAL_TEST_Y_n_cols= 2150;
     GLOBAL_TEST_T_n_cols= 1;
     GLOBAL_TEST_niter= 1;
     END= 0;
  */
  /* test for speed with input file: An_v_speed.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= An_v;
     GLOBAL_TEST_TYP2= speed;
     GLOBAL_NBINS= 4;
     GLOBAL_TEST_mrand= 0.05;
     GLOBAL_TEST_A_n_rows= 2130;
     GLOBAL_TEST_A_n_cols= 192000;
     GLOBAL_TEST_Z_n_rows= 2120;
     GLOBAL_TEST_Y_n_cols= 192000;
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
  struct L_handle **lf_An_v=NULL,**lf_Zn_v=NULL,**lf_Yn_v=NULL,**lf_Wn_v=NULL; int *length_An_v=NULL,*length_Zn_v=NULL,*length_Yn_v=NULL,*length_Wn_v=NULL;
  struct L_handle **lf_An_u=NULL,**lf_Zn_u=NULL,**lf_Yn_u=NULL,**lf_Wn_u=NULL; int *length_An_u=NULL,*length_Zn_u=NULL,*length_Yn_u=NULL,*length_Wn_u=NULL;
  if (verbose){ printf(" %% [entering wrap_An_v_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  lf_An_v = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_v = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Zn_v = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Zn_v = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Yn_v = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Yn_v = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Wn_v = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Wn_v = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_An_v[nb] = L_handle_make(M_An[nb]->nrows);
    lf_Zn_v[nb] = L_handle_make(M_Zn[nb]->nrows);
    lf_Yn_v[nb] = L_handle_make(M_Yn[nb]->nrows);
    lf_Wn_v[nb] = L_handle_make(M_Wn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (error_check){  
    lf_An_u = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_u = (int *)wkspace_all0c(sizeof(int)*nbins);
    lf_Zn_u = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Zn_u = (int *)wkspace_all0c(sizeof(int)*nbins);
    lf_Yn_u = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Yn_u = (int *)wkspace_all0c(sizeof(int)*nbins);
    lf_Wn_u = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Wn_u = (int *)wkspace_all0c(sizeof(int)*nbins);
    for (nb=0;nb<nbins;nb++){
      lf_An_u[nb] = L_handle_make(M_An[nb]->nrows);
      lf_Zn_u[nb] = L_handle_make(M_Zn[nb]->nrows);
      lf_Yn_u[nb] = L_handle_make(M_Yn[nb]->nrows);
      lf_Wn_u[nb] = L_handle_make(M_Wn[nb]->nrows);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (error_check){ } */}
  for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){
      if (verbose){ printf(" %% %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B]);}
      GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
      GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
      for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){
	GLOBAL_pthread_tic(); 
	length_An_v[nb] = wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,M_An[nb],&(lf_An_v[nb]));
	GLOBAL_pthread_toc();
	GLOBAL_pthread_tic(); 
	length_Zn_v[nb] = wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,M_An[nb],&(lf_An_v[nb]));
	GLOBAL_pthread_toc();
	GLOBAL_pthread_tic(); 
	length_Yn_v[nb] = wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,M_An[nb],&(lf_An_v[nb]));
	GLOBAL_pthread_toc();
	GLOBAL_pthread_tic(); 
	length_Wn_v[nb] = wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,M_An[nb],&(lf_An_v[nb]));
	GLOBAL_pthread_toc();
	/* for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){ }} */}}
      GLOBAL_pthread_tuc(); 
      GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_v: ");
      GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");      
      if (error_check){
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	for (nb=0;nb<nbins;nb++){
	  GLOBAL_pthread_tic(); 
	  length_An_u[nb] = wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,M_An[nb],&(lf_An_u[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic(); 
	  length_Zn_u[nb] = wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,M_An[nb],&(lf_An_u[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic(); 
	  length_Yn_u[nb] = wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,M_An[nb],&(lf_An_u[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic(); 
	  length_Wn_u[nb] = wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,M_An[nb],&(lf_An_u[nb]));
	  GLOBAL_pthread_toc();
	  /* for (nb=0;nb<nbins;nb++){ } */}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_u: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	for (nb=0;nb<nbins;nb++){
	  printf(" %% lf_An_v[%d] error %0.16f\n",nb,dra_diff(lf_An_v[nb]->lf,lf_An_u[nb]->lf,length_An_v[nb],1));
	  printf(" %% lf_Zn_v[%d] error %0.16f\n",nb,dra_diff(lf_Zn_v[nb]->lf,lf_Zn_u[nb]->lf,length_Zn_v[nb],1));
	  printf(" %% lf_Yn_v[%d] error %0.16f\n",nb,dra_diff(lf_Yn_v[nb]->lf,lf_Yn_u[nb]->lf,length_Yn_v[nb],1));
	  printf(" %% lf_Wn_v[%d] error %0.16f\n",nb,dra_diff(lf_Wn_v[nb]->lf,lf_Wn_u[nb]->lf,length_Wn_v[nb],1));
	  /* for (nb=0;nb<nbins;nb++){ } */}
	/* if (error_check){ } */}
      /* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ }} */}}
  if (verbose){ printf(" %% [finished wrap_An_v_test]\n");}
}  
