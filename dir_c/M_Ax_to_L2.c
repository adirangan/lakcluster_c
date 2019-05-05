#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void *get_M_An_to_L2_run(void *vp)
{
  /* This function compactifies M_An into an L_handle. ;
     Here An is treated as a binary representation of a matrix with 
     type_flag==TYPE_p_: An taking on +1/-1 values ;
     type_flag==TYPE_0_: An taking on +1/+0 values ;
     We use type_flag, output_spacing_r and output_spacing_c. ;
   */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle **L_p = (struct L_handle **)(vpra[ip++]);
  int type_flag = *(int *)(vpra[ip++]); 
  int output_spacing_r = *(int *)(vpra[ip++]);
  int output_spacing_c = *(int *)(vpra[ip++]);
  int type_flag_single=0;
  if (type_flag==TYPE_p0){ printf(" %% Warning! TYPE_p0 not implemented within get_M_An_to_L2_run\n");}
  if      (type_flag==TYPE_pm){ type_flag_single = TYPE_p_;}
  else if (type_flag==TYPE_00){ type_flag_single = TYPE_0_;}
  if (verbose>1){ printf(" %% [entering get_M_An_to_L2_run] tidx %d rpop_b %d(%d) cpop_b %d(%d) type_flag %s spacing (%s,%s)\n",tidx,(int)(M_An->rpop_b),(int)(M_An->rpop_j),(int)(M_An->cpop_b),(int)(M_An->cpop_j),TYPE_name[type_flag],SPACING_name[output_spacing_r],SPACING_name[output_spacing_c]);}  
  int tab_r_stride=0,tab_c_stride=0;
  int tab_r=0,tab_c=0;
  int ma_a=0,ma_b=0,ma_j=0;
  int na_a=0,na_b=0,na_j=0;
  unsigned char *An_tag=NULL;
  double *dinp=NULL,*ainp=NULL,sdinp=0,dtmp=0;
  int A_pcols = psize(M_An->ncols)/* rup(M_An->ncols + M_An->ncols_extend,POPLENGTH)/POPLENGTH */;
  struct L_handle *L=NULL;
  switch (output_spacing_r){ case SPACING_j: tab_r_stride = M_An->rpop_j; break; case SPACING_b: tab_r_stride = M_An->rpop_b; break; case SPACING_a: tab_r_stride = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: tab_c_stride = M_An->cpop_j; break; case SPACING_b: tab_c_stride = M_An->cpop_b; break; case SPACING_a: tab_c_stride = M_An->ncols; break; default: break; /* switch (output_spacing_c){ } */}
  if (*L_p==NULL){ *L_p = L_handle_make((unsigned long long int)tab_r_stride*(unsigned long long int)tab_c_stride);} L = *L_p;
  if (L->length < tab_r_stride * tab_c_stride){ printf(" %% Warning! L->length %d < %d*%d in get_M_An_to_L2_run\n",L->length,tab_r_stride,tab_c_stride);}
  L_zero(L);
  L->spacing_row = output_spacing_r; L->row_stride = tab_r_stride;
  L->spacing_col = output_spacing_c; L->col_stride = tab_c_stride;
  L->spacing_lyr = SPACING_b; L->lyr_stride = 0;
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    ma_j=0;
    while (ma_j<M_An->rpop_j){
      ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
      switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
      An_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length]));
      na_j=0;
      while (na_j<M_An->cpop_j){
	na_a = M_An->n_a_[na_j]; na_b = M_An->n_b_[na_j];
	switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	if (A_ajdk){
	  dinp = &(A_ajdk[na_a/POPLENGTH+AJDK_0_1*A_pcols]); ainp = &(A_ajdk[na_a/POPLENGTH+AJDK_1_0*A_pcols]);
	  sdinp = sqrt(*dinp);
	  if      (type_flag_single==TYPE_p_){ dtmp = (bget____(An_tag,na_a) - *ainp)*sdinp;}
	  else if (type_flag_single==TYPE_0_){ dtmp = (bget__on(An_tag,na_a) - *ainp)*sdinp;}
	  /* if (A_ajdk){ } */}
	else /* if (A_ajdk==NULL){ } */{
	  if      (type_flag_single==TYPE_p_){ dtmp = bget____(An_tag,na_a);}
	  else if (type_flag_single==TYPE_0_){ dtmp = bget__on(An_tag,na_a);}
	  /* if (A_ajdk==NULL){ } */}
	L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	na_j++; /* while (na_j<M_An->cpop_j){ } */}
      ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)tab_r_stride*(unsigned long long int)tab_c_stride,0);
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp__on){
#pragma omp parallel private(ma_j,ma_a,ma_b,na_j,na_a,na_b,tab_r,tab_c,An_tag,dinp,ainp,sdinp,dtmp)
    { /* begin omp parallel */
      tab_r=0; tab_c=0; dinp=NULL; ainp=NULL; sdinp=0; dtmp=0;
#pragma omp for schedule(dynamic)
      for (ma_j=0;ma_j<M_An->rpop_j;ma_j++){
	ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
	switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
	An_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length]));
	na_j=0;
	while (na_j<M_An->cpop_j){
	  na_a = M_An->n_a_[na_j]; na_b = M_An->n_b_[na_j];
	  switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	  if (A_ajdk){
	    dinp = &(A_ajdk[na_a/POPLENGTH+AJDK_0_1*A_pcols]); ainp = &(A_ajdk[na_a/POPLENGTH+AJDK_1_0*A_pcols]);
	    sdinp = sqrt(*dinp);
	    if      (type_flag_single==TYPE_p_){ dtmp = (bget____(An_tag,na_a) - *ainp)*sdinp;}
	    else if (type_flag_single==TYPE_0_){ dtmp = (bget__on(An_tag,na_a) - *ainp)*sdinp;}
	    /* if (A_ajdk){ } */}
	  else /* if (A_ajdk==NULL){ } */{
	    if      (type_flag_single==TYPE_p_){ dtmp = bget____(An_tag,na_a);}
	    else if (type_flag_single==TYPE_0_){ dtmp = bget__on(An_tag,na_a);}
	    /* if (A_ajdk==NULL){ } */}
	  L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	  na_j++; /* while (na_j<M_An->cpop_j){ } */}
	/* for (ma_j=0;ma_j<M_An->rpop_j;ma_j++){ } */}
      /* end omp parallel */}    
    GLOBAL_ops_count_one(tidx,(unsigned long long int)tab_r_stride*(unsigned long long int)tab_c_stride,0);
    /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
  if (verbose>1){ printf(" %% [finished get_M_An_to_L2_run] tidx %d \n",tidx);}  
  return NULL;
}

int wrap_M_An_to_L2_run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,int output_spacing_c,struct M_handle *M_An,double *A_ajdk,struct L_handle **output_L_p)
{
  /* This function uses the M_handle M_An to run a series of parallel calls to get_M_An_to_L2_run ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 7)
     the type of multiplication is determined by type_flag: ;
     type_flag==TYPE_pm: ;
     type_flag==TYPE_00: ;
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_r=0,length_c=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_M_An_to_L2_run] tidx %d type_flag %d\n",*tidx,type_flag);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  switch (output_spacing_r){ case SPACING_j: length_r = M_An->rpop_j; break; case SPACING_b: length_r = M_An->rpop_b; break; case SPACING_a: length_r = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: length_c = M_An->cpop_j; break; case SPACING_b: length_c = M_An->cpop_b; break; case SPACING_a: length_c = M_An->ncols; break; default: break; /* switch (output_spacing_c){ } */}
  length = length_r*length_c; if (verbose){ printf(" %% length %llu\n",length);}
  length = length_r*length_c; if (*output_L_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_L_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  length = length_r*length_c; if ((*output_L_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_M_An_to_L2_run\n",(*output_L_p)->length,length);} memset((*output_L_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = A_ajdk; vpra[ip++] = output_L_p; 
  switch (type_flag){ case TYPE_p0: vpra[ip++] = &addressable_type_p0; break; case TYPE_pm: vpra[ip++] = &addressable_type_pm; break; case TYPE_00: vpra[ip++] = &addressable_type_00; break; default: break; /* switch (type_flag){ } */}
  switch (output_spacing_r){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_c){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_M_An_to_L2_run,vpra)){ printf("Warning! cannot create thread %d in wrap_M_An_to_L2_run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_M_An_to_L2_run(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_L_p)->lf,"double",1,length_r," %% (*output_L_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_M_An_to_L2_run] tidx %d\n",*tidx);}
  return length;
}  

void *get_M_At_to_L2_run(void *vp)
{
  /* This function compactifies M_At into an L_handle. ;
     Here At is treated as a binary representation of a matrix with 
     type_flag==TYPE_p_: At taking on +1/-1 values ;
     type_flag==TYPE_0_: At taking on +1/+0 values ;
     We use type_flag, output_spacing_r and output_spacing_c. ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle **L_p = (struct L_handle **)(vpra[ip++]);
  int type_flag = *(int *)(vpra[ip++]); 
  int output_spacing_r = *(int *)(vpra[ip++]);
  int output_spacing_c = *(int *)(vpra[ip++]);
  int type_flag_single=0;
  if (type_flag==TYPE_p0){ printf(" %% Warning! TYPE_p0 not implemented within get_M_At_to_L2_run\n");}
  if      (type_flag==TYPE_pm){ type_flag_single = TYPE_p_;}
  else if (type_flag==TYPE_00){ type_flag_single = TYPE_0_;}
  if (verbose>1){ printf(" %% [entering get_M_At_to_L2_run] tidx %d rpop_b %d(%d) cpop_b %d(%d) type_flag %s spacing (%s,%s)\n",tidx,(int)(M_At->rpop_b),(int)(M_At->rpop_j),(int)(M_At->cpop_b),(int)(M_At->cpop_j),TYPE_name[type_flag],SPACING_name[output_spacing_r],SPACING_name[output_spacing_c]);}  
  int tab_r_stride=0,tab_c_stride=0;
  int tab_r=0,tab_c=0;
  int ma_a=0,ma_b=0,ma_j=0;
  int na_a=0,na_b=0,na_j=0;
  unsigned char *At_tag=NULL;
  double *dinp=NULL,*ainp=NULL,sdinp=0,dtmp=0;
  int A_pcols = psize(M_At->nrows)/* rup(M_An->ncols + M_An->ncols_extend,POPLENGTH)/POPLENGTH */;
  struct L_handle *L=NULL;
  switch (output_spacing_r){ case SPACING_j: tab_r_stride = M_At->rpop_j; break; case SPACING_b: tab_r_stride = M_At->rpop_b; break; case SPACING_a: tab_r_stride = M_At->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: tab_c_stride = M_At->cpop_j; break; case SPACING_b: tab_c_stride = M_At->cpop_b; break; case SPACING_a: tab_c_stride = M_At->ncols; break; default: break; /* switch (output_spacing_c){ } */}
  if (*L_p==NULL){ *L_p = L_handle_make((unsigned long long int)tab_r_stride*(unsigned long long int)tab_c_stride);} L = *L_p;
  if (L->length < (unsigned long long int)tab_r_stride * (unsigned long long int)tab_c_stride){ printf(" %% Warning! L->length %llu < %d*%d in get_M_At_to_L2_run\n",L->length,tab_r_stride,tab_c_stride);}
  L_zero(L);
  L->spacing_row = output_spacing_r; L->row_stride = tab_r_stride;
  L->spacing_col = output_spacing_c; L->col_stride = tab_c_stride;
  L->spacing_lyr = SPACING_b; L->lyr_stride = 0;
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    ma_j=0;
    while (ma_j<M_At->rpop_j){
      ma_a = M_At->m_a_[ma_j]; ma_b = M_At->m_b_[ma_j];
      switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
      if (A_ajdk){
	dinp = &(A_ajdk[ma_a/POPLENGTH+AJDK_0_1*A_pcols]); sdinp = sqrt(*dinp);
	ainp = &(A_ajdk[ma_a/POPLENGTH+AJDK_1_0*A_pcols]);
	/* if (A_ajdk){ } */}
      At_tag = (unsigned char *)(&(M_At->wX[ma_b*M_At->mc_length]));
      if (A_ajdk && type_flag_single==TYPE_p_){
	na_j=0;
	while (na_j<M_At->cpop_j){
	  na_a = M_At->n_a_[na_j]; na_b = M_At->n_b_[na_j];
	  switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	  dtmp = (bget____(At_tag,na_a) - *ainp)*sdinp;
	  L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	  na_j++; /* while (na_j<M_At->cpop_j){ } */}
	/* if (A_ajdk && type_flag_single==TYPE_p_){ } */}
      if (A_ajdk && type_flag_single==TYPE_0_){
	na_j=0;
	while (na_j<M_At->cpop_j){
	  na_a = M_At->n_a_[na_j]; na_b = M_At->n_b_[na_j];
	  switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	  dtmp = (bget__on(At_tag,na_a) - *ainp)*sdinp;
	  L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	  na_j++; /* while (na_j<M_At->cpop_j){ } */}
	/* if (A_ajdk && type_flag_single==TYPE_0_){ } */}
      if (A_ajdk==NULL && type_flag_single==TYPE_p_){
	na_j=0;
	while (na_j<M_At->cpop_j){
	  na_a = M_At->n_a_[na_j]; na_b = M_At->n_b_[na_j];
	  switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	  dtmp = bget____(At_tag,na_a);
	  L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	  na_j++; /* while (na_j<M_At->cpop_j){ } */}
	/* if (A_ajdk && type_flag_single==TYPE_p_){ } */}
      if (A_ajdk==NULL && type_flag_single==TYPE_0_){
	na_j=0;
	while (na_j<M_At->cpop_j){
	  na_a = M_At->n_a_[na_j]; na_b = M_At->n_b_[na_j];
	  switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	  dtmp = bget__on(At_tag,na_a);
	  L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	  na_j++; /* while (na_j<M_At->cpop_j){ } */}
	/* if (A_ajdk && type_flag_single==TYPE_0_){ } */}
      ma_j++; /* while (ma_j<M_At->rpop_j){ } */}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)tab_r_stride*(unsigned long long int)tab_c_stride,0);
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp__on){
#pragma omp parallel private(ma_j,ma_a,ma_b,na_j,na_a,na_b,tab_r,tab_c,At_tag,dinp,ainp,sdinp,dtmp)
    { /* begin omp parallel */
      tab_r=0; tab_c=0; dinp=NULL; ainp=NULL; sdinp=0; dtmp=0;
#pragma omp for schedule(dynamic)
      for (ma_j=0;ma_j<M_At->rpop_j;ma_j++){
	ma_a = M_At->m_a_[ma_j]; ma_b = M_At->m_b_[ma_j];
	switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
	if (A_ajdk){
	  dinp = &(A_ajdk[ma_a/POPLENGTH+AJDK_0_1*A_pcols]); sdinp = sqrt(*dinp);
	  ainp = &(A_ajdk[ma_a/POPLENGTH+AJDK_1_0*A_pcols]);
	  /* if (A_ajdk){ } */}
	At_tag = (unsigned char *)(&(M_At->wX[ma_b*M_At->mc_length]));
	if (A_ajdk && type_flag_single==TYPE_p_){
	  na_j=0;
	  while (na_j<M_At->cpop_j){
	    na_a = M_At->n_a_[na_j]; na_b = M_At->n_b_[na_j];
	    switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	    dtmp = (bget____(At_tag,na_a) - *ainp)*sdinp;
	    L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	    na_j++; /* while (na_j<M_At->cpop_j){ } */}
	  /* if (A_ajdk && type_flag_single==TYPE_p_){ } */}
	if (A_ajdk && type_flag_single==TYPE_0_){
	  na_j=0;
	  while (na_j<M_At->cpop_j){
	    na_a = M_At->n_a_[na_j]; na_b = M_At->n_b_[na_j];
	    switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	    dtmp = (bget__on(At_tag,na_a) - *ainp)*sdinp;
	    L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	    na_j++; /* while (na_j<M_At->cpop_j){ } */}
	  /* if (A_ajdk && type_flag_single==TYPE_0_){ } */}
	if (A_ajdk==NULL && type_flag_single==TYPE_p_){
	  na_j=0;
	  while (na_j<M_At->cpop_j){
	    na_a = M_At->n_a_[na_j]; na_b = M_At->n_b_[na_j];
	    switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	    dtmp = bget____(At_tag,na_a);
	    L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	    na_j++; /* while (na_j<M_At->cpop_j){ } */}
	  /* if (A_ajdk && type_flag_single==TYPE_p_){ } */}
	if (A_ajdk==NULL && type_flag_single==TYPE_0_){
	  na_j=0;
	  while (na_j<M_At->cpop_j){
	    na_a = M_At->n_a_[na_j]; na_b = M_At->n_b_[na_j];
	    switch (output_spacing_c){ case SPACING_j: tab_c=na_j; break; case SPACING_b: tab_c=na_b; break; case SPACING_a: tab_c=na_a; break; default: break; /* switch (output_spacing_c){ } */}      
	    dtmp = bget__on(At_tag,na_a);
	    L2_set(L,ma_j,ma_b,ma_a,na_j,na_b,na_a,dtmp);
	    na_j++; /* while (na_j<M_At->cpop_j){ } */}
	  /* if (A_ajdk && type_flag_single==TYPE_0_){ } */}
	/* for (ma_j=0;ma_j<M_At->rpop_j;ma_j++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,(unsigned long long int)tab_r_stride*(unsigned long long int)tab_c_stride,0);
    /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
  if (verbose>1){ printf(" %% [finished get_M_At_to_L2_run] tidx %d \n",tidx);}  
  return NULL;
}

int wrap_M_At_to_L2_run(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,int output_spacing_c,struct M_handle *M_At,double *A_ajdk,struct L_handle **output_L_p)
{
  /* This function uses the M_handle M_At to run a series of parallel calls to get_M_At_to_L2_run ; 
     No reloading of data is performed ; we assume all data is preloaded ;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 7)
     the type of multiplication is determined by type_flag: ;
     type_flag==TYPE_pm: ;
     type_flag==TYPE_00: ;
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length_r=0,length_c=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_M_At_to_L2_run] tidx %d type_flag %d\n",*tidx,type_flag);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  switch (output_spacing_r){ case SPACING_j: length_r = M_At->rpop_j; break; case SPACING_b: length_r = M_At->rpop_b; break; case SPACING_a: length_r = M_At->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: length_c = M_At->cpop_j; break; case SPACING_b: length_c = M_At->cpop_b; break; case SPACING_a: length_c = M_At->ncols; break; default: break; /* switch (output_spacing_c){ } */}
  length = length_r*length_c; if (verbose){ printf(" %% length %llu\n",length);}
  length = length_r*length_c; if (*output_L_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_L_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: ");}
  if (verbose>2){ bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: ");}
  if (verbose>2){ bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  length = length_r*length_c; if ((*output_L_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_M_At_to_L2_run\n",(*output_L_p)->length,length);} memset((*output_L_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = A_ajdk; vpra[ip++] = output_L_p; 
  switch (type_flag){ case TYPE_p0: vpra[ip++] = &addressable_type_p0; break; case TYPE_pm: vpra[ip++] = &addressable_type_pm; break; case TYPE_00: vpra[ip++] = &addressable_type_00; break; default: break; /* switch (type_flag){ } */}
  switch (output_spacing_r){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_r){ } */}
  switch (output_spacing_c){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing_c){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_M_At_to_L2_run,vpra)){ printf("Warning! cannot create thread %d in wrap_M_At_to_L2_run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_M_At_to_L2_run(vpra);} /* must join threads later */;
  if (verbose>2){ raprintf((*output_L_p)->lf,"double",1,length_r," %% (*output_L_p)->lf: ");}
  if (verbose){ printf(" %% [finished wrap_M_At_to_L2_run] tidx %d\n",*tidx);}
  return length;
}  

void wrap_M_Ax_to_L2_test()
{
  /* test for errors with input file: M_Ax_to_L2_error.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 1;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= M_Ax_to_L2;
     GLOBAL_TEST_TYP2= error;
     GLOBAL_NBINS= 1;
     GLOBAL_TEST_mrand= 0.5;
     GLOBAL_TEST_A_n_rows= 32;
     GLOBAL_TEST_A_n_cols= 64;
     GLOBAL_TEST_Z_n_rows= 48;
     GLOBAL_TEST_Y_n_cols= 96;
     GLOBAL_TEST_T_n_cols= 1;
     GLOBAL_TEST_niter= 1;
     END= 0;
  */
  /* test for speed with input file: M_Ax_to_L2_speed.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= M_Ax_to_L2;
     GLOBAL_TEST_TYP2= speed;
     GLOBAL_NBINS= 1;
     GLOBAL_TEST_mrand= 0.05;
     GLOBAL_TEST_A_n_rows= 9800;
     GLOBAL_TEST_A_n_cols= 350;
     GLOBAL_TEST_Z_n_rows= 9800;
     GLOBAL_TEST_Y_n_cols= 0;
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
  struct L_handle **lf_An=NULL,**lf_At=NULL,**lf_Bn=NULL,**lf_Bt=NULL;
  int *length_An=NULL,*length_At=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  if (verbose){ printf(" %% [entering wrap_M_Ax_to_L2_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  lf_An = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_At = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Bn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); 
  lf_Bt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); 
  for (nb=0;nb<nbins;nb++){
    lf_An[nb] = L_handle_make((unsigned long long int)A_n_rows[nb]*(unsigned long long int)A_n_cols);
    lf_At[nb] = L_handle_make((unsigned long long int)A_n_rows[nb]*(unsigned long long int)A_n_cols);
    lf_Bn[nb] = L_handle_make((unsigned long long int)A_n_rows[nb]*(unsigned long long int)A_n_cols);
    lf_Bt[nb] = L_handle_make((unsigned long long int)A_n_rows[nb]*(unsigned long long int)A_n_cols);
    /* for (nb=0;nb<nbins;nb++){ } */}
  for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
  	if (verbose){ printf(" %% %s; B: %s; A: %s\n",TYPE_name[n_type],SPACING_name[n_spacing_B],SPACING_name[n_spacing_A]);}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic(); 
	    length_An[nb] = wrap_M_An_to_L2_run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_B,M_An[nb],A_ajdk,&(lf_At[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){ }} */}}
	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% M_An_to_L2: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic();
	    length_At[nb] = wrap_M_At_to_L2_run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_B,n_spacing_A,M_At[nb],A_ajdk,&(lf_An[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){ }} */}}
	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% M_At_to_L2: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	if (error_check){
	  for (nb=0;nb<nbins;nb++){
	    L2_transpose(lf_Bt[nb],lf_At[nb]);
	    L2_transpose(lf_Bn[nb],lf_An[nb]);
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% M_An[%d] (%d,%d,%d)-x-(%d,%d,%d) lf_At[%d] (%d,%d) error %0.16f\n",nb,M_An[nb]->rpop_j,M_An[nb]->rpop_b,M_An[nb]->nrows,M_An[nb]->cpop_j,M_An[nb]->cpop_b,M_An[nb]->ncols,nb,lf_At[nb]->row_stride,lf_At[nb]->col_stride,dra_diff(lf_At[nb]->lf,lf_Bn[nb]->lf,length_At[nb],1));
	    printf(" %% M_At[%d] (%d,%d,%d)-x-(%d,%d,%d) lf_An[%d] (%d,%d) error %0.16f\n",nb,M_At[nb]->rpop_j,M_At[nb]->rpop_b,M_At[nb]->nrows,M_At[nb]->cpop_j,M_At[nb]->cpop_b,M_At[nb]->ncols,nb,lf_An[nb]->row_stride,lf_An[nb]->col_stride,dra_diff(lf_An[nb]->lf,lf_Bt[nb]->lf,length_An[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_M_Ax_to_L2_test]\n");}
}
