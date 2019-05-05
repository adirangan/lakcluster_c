#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void xcalc_setmaxmin(struct L_handle *lf_YnWt,double *YnWt,int rpop_j,int rpop_b,int nrows,int cpop_j,int cpop_b,int ncols,unsigned char *mr_j,unsigned char *mr_b,unsigned char *mc_j,unsigned char *mc_b,int trm_flag,double *max_data_p,double *min_data_p)
{
  /* trm_flag==0 --> rpop_j, rpop_b, nrows, mr_j, mr_b all associated with Yn->mr_j and Yn->mr_b ; 
     trm_flag==1 --> rpop_j, rpop_b, nrows, mr_j, mr_b all associated with Wn->mr_j and Wn->mr_b ; 
     consequently, when trm_flag==1 we read lf_YnWt->lf using nc_x + nr_x*lf_YnWt->row_stride ;  */
  int verbose=0;
  int nc_a=0,nc_b=0,nc_j=0,nr_a=0,nr_b=0,nr_j=0,m_set=0;
  double lftmp=0,min_data=0,max_data=0,sum1_data=0,sum2_data=0,avg_data=0,var_data=0,std_data=0;
  if (verbose){ printf(" %% [entering xcalc_setmaxmin] rpop_j %d rpop_b %d nrows %d cpop_j %d cpop_b %d ncols %d\n",rpop_j,rpop_b,nrows,cpop_j,cpop_b,ncols);}
  m_set=0; nc_a=0;nc_b=0;nc_j=0;
  while (nc_a<ncols){
    if (mc_b==NULL || bget__on(mc_b,nc_a)){
      if (mc_j==NULL || bget__on(mc_j,nc_a)){
	if (verbose>1){ printf(" %% nc_a %d nc_b %d nc_j %d\n",nc_a,nc_b,nc_j);}
	nr_a=0;nr_b=0;nr_j=0;
	while (nr_a<nrows){
	  if (mr_b==NULL || bget__on(mr_b,nr_a)){
	    if (mr_j==NULL || bget__on(mr_j,nr_a)){
	      if (verbose>1){ printf(" %% nr_a %d nr_b %d nr_j %d \n",nr_a,nr_b,nr_j);}
	      lftmp = trm_flag ? (*L2_lf_get(lf_YnWt,YnWt,nc_j,nc_b,nc_a,nr_j,nr_b,nr_a)) : (*L2_lf_get(lf_YnWt,YnWt,nr_j,nr_b,nr_a,nc_j,nc_b,nc_a)) ;
	      if (!m_set){ min_data = lftmp; max_data = lftmp; m_set=1;} else /* if m_set */{ min_data = minimum(lftmp,min_data); max_data = maximum(lftmp,max_data);}
	      if (verbose>1){ printf(" %% lftmp %f --> min_data %f max_data %f\n",lftmp,min_data,max_data);}
	      sum1_data += lftmp; sum2_data += lftmp*lftmp;
	      nr_j++; /* if (mr_j==NULL || bget__on(mr_j,nr_a)){ } */}
	    nr_b++; /* if (mr_b==NULL || bget__on(mr_b,nr_a)){ } */}
	  nr_a++; /* while (nr_a<nrows){ } */}
	nc_j++; /* if (mc_j==NULL || bget__on(mc_j,nc_a)){ } */}
      nc_b++; /* if (mc_b==NULL || bget__on(mc_b,nc_a)){ } */}
    nc_a++; /* while (nc_a<ncols){ } */}
  if (min_data>=max_data){ max_data = min_data+1;}
  avg_data = sum1_data / maximum(1,rpop_j*cpop_j) ; var_data = sum2_data / maximum(1,rpop_j*cpop_j) - avg_data*avg_data ; std_data = sqrt(maximum(0,var_data)) ;
  if (verbose){ printf(" %% min %f max %f; sum1 %f avg %f sum2 %f var %f std %f\n",min_data,max_data,sum1_data,avg_data,sum2_data,var_data,std_data);}
  if (min_data_p!=NULL){ *min_data_p = min_data;} if (max_data_p!=NULL){ *max_data_p = max_data;}
  if (verbose){ printf(" %% [finished xcalc_setmaxmin]\n");}  
}

void *get_xcalc(void *vp)
{
  /* This function is designed for use with TYPE_pZ and TYPE_0Z computation, with simple parallelization. ;
     This can allow an extra bit to be used as a buffer for future additive updates (not implemented now);
     Given double array *YnWt, we create bit array bZ, with first column 1's.
     For example given *YnWt = 
     [ -5 ,  1 ] ;
     [  3 ,  0 ] ;
     [  2 , -1 ] ;
     [  8 , -2 ] ;
     we note min=-5 and max = 8, for a total of 13 in range.
     with a minbuffer of 0.25, we would set min = -5 - 0.25*13 = -9;
     and then subtract to obtain
     [  4 , 10 ] ;
     [ 12 ,  9 ] ;
     [ 11 ,  8 ] ;
     [ 17 ,  7 ] ;
     and then multiply by mlt = D_MLT/(max-min) (e.g., mlt = 128), to obtain
     [  4*mlt , 10*mlt ] ;
     [ 12*mlt ,  9*mlt ] ;
     [ 11*mlt ,  8*mlt ] ;
     [ 17*mlt ,  7*mlt ] ;
     at this point ceil(log2(max+1)) = 12
     so we need at least 12 columns per row
     since:
     4*mlt   = 001000000000
     12*mlt  = 011000000000
     11*mlt  = 010110000000
     17*mlt  = 100010000000
     10*mlt  = 010100000000
     9*mlt   = 010010000000
     8*mlt   = 010000000000
     7*mlt   = 001110000000
     the final bZ is something like:
     [ 1 001000000000 010100000000 ]
     [ 1 011000000000 010010000000 ]
     [ 1 010110000000 010000000000 ]
     [ 1 100010000000 001110000000 ]
     [ 0 000000000000 000000000000 ]
     [ 0 000000000000 000000000000 ]
     [ 0 000000000000 000000000000 ]
     [ 0 000000000000 000000000000 ]
     With this compressed notation, we can multiply bZ by a vector "A" of length 4 producing: ;
     A*bZ = [ xx y_{11} y_{10} ... y_{0} z_{11} z_{10} ... z_{0} ],
     with A*YnWt = B, where
     B = [ min*mlt*xx + 2^11*y_{11} + 2^10*y_{10} + ... + 2^0*y_{0} , min*mlt*xx + 2^11*z_{11} + 2^10*z_{10} + ... 2^0*z_{0} ] / mlt . ;
     The number of bits of precision retained will be roughly log2(D_MLT). ;
     Note that we pad bZ to account for POPLENGTH (i.e., for direct use as M_X->wX). ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  unsigned char *mr_j = (unsigned char *)(vpra[ip++]); /* original M_Yn->mr_j used to construct YnWt */
  unsigned char *mr_b = (unsigned char *)(vpra[ip++]); /* original M_Yn->mr_b used to construct YnWt */
  unsigned char *mc_j = (unsigned char *)(vpra[ip++]); /* original M_Wn->mr_j used to construct YnWt */
  unsigned char *mc_b = (unsigned char *)(vpra[ip++]); /* original M_Wn->mr_b used to construct YnWt */
  struct L_handle *lf_YnWt = (struct L_handle *)(vpra[ip++]);
  double *YnWt = (double *)(vpra[ip++]);
  int Yn_nrows = *(int *)(vpra[ip++]); /* original M_Yn->nrows used to construct YnWt */
  int Wn_nrows = *(int *)(vpra[ip++]); /* original M_Wn->nrows used to construct YnWt */
  struct M_handle **M_X_p = (struct M_handle **)(vpra[ip++]);
  int b_mlt = *(int *)(vpra[ip++]); /* number of bits to store */
  int trm_flag = *(int *)(vpra[ip++]); /* 0 --> YnWt, 1 --> WnYt (must switch masks and rows at input) */
  int rpop_b = 0,cpop_b=0; int rpop_j = 0,cpop_j=0;
  struct M_handle *M_X = NULL;
  int length = Yn_nrows*Wn_nrows;
  int ncols_bZ=0;
  int nrows_bZ=0,brows_bZ=0,ncols_per=0;
  double lftmp=0,min_data=0,max_data=0,min=0,max=0,mlt=0;
  int nb=0,nr_a=0,nc_a=0,nr_m8=0,nr_m7=0,nr_d8=0,tab_nb=0,nr_j=0,nc_j=0,nr_b=0,nc_b=0;
  unsigned long long int cz=0,cz_tmp=0,n2=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned char **bZ_p=NULL,*bZ;
  int mx_j=0,mx_chunk=0;
  long long int Wn_rpop_b=0,Wn_rpop_j=0; unsigned int Wn_m_a_[Wn_nrows],Wn_m_b_[Wn_nrows];
  /* int Yn_rpop_b=0,Yn_rpop_j=0; unsigned int Yn_m_a_[Wn_nrows],Yn_m_b_[Wn_nrows]; */
  if (verbose){ printf(" %% [entering get_xcalc] tidx %d trm_flag %d\n",tidx,trm_flag);}
  if (mr_j==NULL || mr_b==NULL || mc_j==NULL || mc_b==NULL){ printf(" %% Warning! null mask in get_xcalc\n");}
  if (verbose>2){ bprintf(mr_b,BITJ,1,Yn_nrows," %% mr_b: ");}
  if (verbose>2){ bprintf(mr_j,BITJ,1,Yn_nrows," %% mr_j: ");}
  if (verbose>2){ bprintf(mc_b,BITJ,1,Wn_nrows," %% mc_b: ");}
  if (verbose>2){ bprintf(mc_j,BITJ,1,Wn_nrows," %% mc_j: ");}
  rpop_j = popcount_uchar_array(mr_j,bsize(Yn_nrows)); cpop_j = popcount_uchar_array(mc_j,bsize(Wn_nrows));
  rpop_b = popcount_uchar_array(mr_b,bsize(Yn_nrows)); cpop_b = popcount_uchar_array(mc_b,bsize(Wn_nrows));
  if (verbose>2){ printf(" %% rpop_j %d rpop_b %d Yn_nrows %d; cpop_j %d cpop_b %d Wn_nrows %d\n",rpop_j,rpop_b,Yn_nrows,cpop_j,cpop_b,Wn_nrows);}
  if (verbose>2){ printf(" %% row_stride %d col_stride %d\n",lf_YnWt->row_stride,lf_YnWt->col_stride);}
  if (verbose>2){ raprintf(lf_YnWt->lf,"double",lf_YnWt->row_stride,lf_YnWt->col_stride," %% lf_YnWt->lf: ");}
  xcalc_setmaxmin(lf_YnWt,YnWt,rpop_j,rpop_b,Yn_nrows,cpop_j,cpop_b,Wn_nrows,mr_j,mr_b,mc_j,mc_b,trm_flag,&max_data,&min_data);
  min = min_data - 0.0625*(max_data - min_data); max = max_data + 0.0625*(max_data - min_data); mlt = (double)(GLOBAL_D_MLT)/(max - min); 
  if (verbose>2){ printf(" %% min_data %0.2f max_data %0.2f min %0.2f max %0.2f\n",min_data,max_data,min,max);}
  if (*M_X_p==NULL){ if (GLOBAL_thread_count>1){ printf(" %% Warning! *M_X_p undefined in get_xcalc; must pre-allocate to avoid thread conflict!\n");} bZ_p = (unsigned char **)wkspace_all0c(1*sizeof(unsigned char *)); *bZ_p = NULL;}
  else if (*M_X_p!=NULL){ M_X = *M_X_p; bZ_p = &(M_X->Ara);}
  if (length>0 && lf_YnWt!=NULL){
    ncols_per = b_mlt + 1; ncols_bZ = ncols_per*cpop_j; if (ncols_per>63){ printf(" %% Warning! overflow in get_xcalc!\n");} if (verbose>2){ printf(" %% ncols_per %d, ncols_bZ %d\n",ncols_per,ncols_bZ);}
    nrows_bZ = Yn_nrows; brows_bZ = bsize(nrows_bZ);
    if (*bZ_p==NULL){ if (GLOBAL_thread_count>1){ printf(" %% Warning! *bZ_p undefined in get_xcalc; must pre-allocate to avoid thread conflict!\n");} *bZ_p = wkspace_all0c(ncols_bZ*brows_bZ);} 
    fill_uchar_zero(*bZ_p,ncols_bZ*brows_bZ);
    if (strstr(GLOBAL_skip,"xcalc")){ goto skip_xcalc;}
    if (verbose>2){ printf(" %% nrows_bZ %d, brows_bZ %d\n",nrows_bZ,brows_bZ);}
    nr_m8=0;nr_m7=7;nr_d8=0; for (nr_a=0;nr_a<Yn_nrows;nr_a++){ if (bget__on(mr_j,nr_a)){ (*bZ_p)[nr_d8 + 0*brows_bZ] |= (1 << (nr_m7));} nr_m8++; nr_m7--; if (nr_m8>=BIT8){ nr_d8++; nr_m8=0; nr_m7=7;}}
    if (verbose>2){ bprintf(&((*bZ_p)[0*brows_bZ]),POPLENGTH,1,Yn_nrows," %% (*bZ_p)[0*brows_bZ]: ");}
    if (GLOBAL_omp_type==GLOBAL_omp_off){
      nc_j=0; nc_b=0; 
      for (nc_a=0;nc_a<Wn_nrows;nc_a++){
	if (bget__on(mc_b,nc_a)){
	  if (bget__on(mc_j,nc_a)){
	    if (verbose>2){ printf(" %% nc_a %d nc_b %d nc_j %d\n",nc_a,nc_b,nc_j);}
	    nr_m8=0;nr_m7=7;nr_d8=0;nr_j=0;nr_b=0;
	    for (nr_a=0;nr_a<Yn_nrows;nr_a++){
	      if (bget__on(mr_b,nr_a)){
		if (bget__on(mr_j,nr_a)){
		  if (verbose>2){ printf(" %% nr_a %d nr_b %d nr_j %d\n",nr_a,nr_b,nr_j);}
		  lftmp = trm_flag ? (*L2_lf_get(lf_YnWt,YnWt,nc_j,nc_b,nc_a,nr_j,nr_b,nr_a)) : (*L2_lf_get(lf_YnWt,YnWt,nr_j,nr_b,nr_a,nc_j,nc_b,nc_a)) ;
		  cz = (unsigned long long int)(mlt*(double)(lftmp-min)); cz = crop(cz,0+1,GLOBAL_D_MLT-1);
		  if (verbose>2){ printf(" %% cz %f --> %lld\n",lftmp,cz);}
		  if (verbose>2){ printf(" %% nr_d8 %d ncols_per %d nc_j %d brows_bZ %d tab_nb [%d,..,%d]\n",nr_d8,ncols_per,nc_j,brows_bZ,nr_d8 + (ncols_per*nc_j + 1)*brows_bZ,nr_d8 + (ncols_per*nc_j + ncols_per-1)*brows_bZ);}
		  nb=0; tab_nb = nr_d8 + (ncols_per*nc_j + ncols_per-1-nb)*brows_bZ; cz_tmp = cz;
		  for (nb=0;nb<ncols_per-1;nb++){ (*bZ_p)[tab_nb] |= ((unsigned char)(cz_tmp & 1) << (nr_m7)); cz_tmp >>= 1; tab_nb -= brows_bZ;}
		  nr_j++; /* if (bget__on(mr_j,nr_a)){ } */}
		nr_b++; /* if (bget__on(mr_b,nr_a)){ } */}
	      nr_m8++; nr_m7--; if (nr_m8>=BIT8){ nr_d8++; nr_m8=0; nr_m7=7;}
	      /* for (nr_a=0;nr_a<Yn_nrows;nr_a++){ } */}
	    nc_j++; /* if (bget__on(mc_j,nc_a)){ } */}
	  nc_b++; /* if (bget__on(mc_b,nc_a)){ } */}
	/* for (nc_a=0;nc_a<Wn_nrows;nc_a++){ } */}
      /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
    else if (GLOBAL_omp_type==GLOBAL_omp__on){
      M_mxget_excerpt(verbose,Wn_nrows,mc_b,mc_j,&Wn_rpop_b,&Wn_rpop_j,Wn_m_a_,Wn_m_b_);
      if (verbose>1){
	raprintf(Wn_m_a_,"unsigned int",1,Wn_rpop_j," %% Wn_m_a_: ");
	raprintf(Wn_m_b_,"unsigned int",1,Wn_rpop_j," %% Wn_m_b_: ");
	/* if (verbose){ } */}
      /* M_mxget_excerpt(verbose,Yn_nrows,mr_b,mr_j,&Yn_rpop_b,&Yn_rpop_j,Yn_m_a_,Yn_m_b_); */
      mx_chunk = 128;
#pragma omp parallel private(mx_j,nc_j,nc_a,nc_b,nr_m8,nr_m7,nr_d8,nr_j,nr_a,nr_b,lftmp,cz,nb,tab_nb,cz_tmp)
      { /* begin omp parallel */
	mx_j=0;
#pragma omp for schedule(dynamic,mx_chunk)
	for (mx_j=0;mx_j<Wn_rpop_j;mx_j++){
	  nc_j=mx_j; nc_a = Wn_m_a_[nc_j]; nc_b = Wn_m_b_[nc_j];
	  if (verbose>2){ printf(" %% nc_a %d nc_b %d nc_j %d\n",nc_a,nc_b,nc_j);}
	  nr_m8=0;nr_m7=7;nr_d8=0;nr_j=0;nr_b=0;
	  for (nr_a=0;nr_a<Yn_nrows;nr_a++){
	    if (bget__on(mr_b,nr_a)){
	      if (bget__on(mr_j,nr_a)){
		if (verbose>2){ printf(" %% nr_a %d nr_b %d nr_j %d\n",nr_a,nr_b,nr_j);}
		lftmp = trm_flag ? (*L2_lf_get(lf_YnWt,YnWt,nc_j,nc_b,nc_a,nr_j,nr_b,nr_a)) : (*L2_lf_get(lf_YnWt,YnWt,nr_j,nr_b,nr_a,nc_j,nc_b,nc_a)) ;
		cz = (unsigned long long int)(mlt*(double)(lftmp-min)); cz = crop(cz,0+1,GLOBAL_D_MLT-1);
		if (verbose>2){ printf(" %% cz %f --> %lld\n",lftmp,cz);}
		if (verbose>2){ printf(" %% nr_d8 %d ncols_per %d nc_j %d brows_bZ %d tab_nb [%d,..,%d]\n",nr_d8,ncols_per,nc_j,brows_bZ,nr_d8 + (ncols_per*nc_j + 1)*brows_bZ,nr_d8 + (ncols_per*nc_j + ncols_per-1)*brows_bZ);}
		nb=0; tab_nb = nr_d8 + (ncols_per*nc_j + ncols_per-1-nb)*brows_bZ; cz_tmp = cz;
		for (nb=0;nb<ncols_per-1;nb++){ (*bZ_p)[tab_nb] |= ((unsigned char)(cz_tmp & 1) << (nr_m7)); cz_tmp >>= 1; tab_nb -= brows_bZ;}
		nr_j++; /* if (bget__on(mr_j,nr_a)){ } */}
	      nr_b++; /* if (bget__on(mr_b,nr_a)){ } */}
	    nr_m8++; nr_m7--; if (nr_m8>=BIT8){ nr_d8++; nr_m8=0; nr_m7=7;}
	    /* for (nr_a=0;nr_a<Yn_nrows;nr_a++){ } */}
	  /* for (mx_j=0;mx_j<Wn_rpop_j;mx_j++){ } */}
	/* end omp parallel */}
      /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
    /* if (length>0 && lf_YnWt!=NULL){ } */}
  //if (verbose>3 && (*bZ_p)!=NULL){ bprintf((*bZ_p),POPLENGTH,ncols_bZ,brows_bZ*BIT8," %% (*bZ_p)[ ]: ");}
  //if (verbose>2){ for (nc_a=0;nc_a<ncols_per*cpop_j;nc_a++){ bprintf(&((*bZ_p)[nc_a*brows_bZ]),POPLENGTH,1,brows_bZ*BIT8," %% (*bZ_p)[ ]: ");}}
  if (verbose>2){ 
    nr_m8=0;nr_m7=7;nr_d8=0;nr_j=0;nr_b=0;
    for (nr_a=0;nr_a<Yn_nrows;nr_a++){
      if (bget__on(mr_b,nr_a)){
	if (bget__on(mr_j,nr_a)){
	  printf(" %% nr(%.3d,%.3d.%.3d): ",nr_a,nr_b,nr_j);
	  nc_j=0; nc_b=0; 
	  for (nc_a=0;nc_a<Wn_nrows;nc_a++){
	    if (bget__on(mc_b,nc_a)){
	      if (bget__on(mc_j,nc_a)){
		cz = (unsigned long long int)(mlt*(double)(lftmp-min)); cz = crop(cz,0+1,GLOBAL_D_MLT-1);
		lftmp = trm_flag ? (*L2_lf_get(lf_YnWt,YnWt,nc_j,nc_b,nc_a,nr_j,nr_b,nr_a)) : (*L2_lf_get(lf_YnWt,YnWt,nr_j,nr_b,nr_a,nc_j,nc_b,nc_a)) ;		           
		printf(" %+3.2f (%.3lld): nc(%.3d,%.3d,%.3d) ",lftmp,cz,nc_a,nc_b,nc_j);
		n2=1; cz_tmp=0;
		nb=0; tab_nb = nr_d8 + (ncols_per*nc_j + ncols_per-1-nb)*brows_bZ;
		for (nb=0;nb<ncols_per-1;nb++){ bZ = &((*bZ_p)[tab_nb]); printf("+%d*%d",n2,bget__on(bZ,nr_m8)); cz_tmp = cz_tmp += n2*bget__on(bZ,nr_m8); tab_nb -= brows_bZ; n2*=2;}
		printf(" %.3lld]",cz_tmp);
		nc_j++; /* if (bget__on(mc_j,nc_a)){ } */}
	      nc_b++; /* if (bget__on(mc_b,nc_a)){ } */}
	    /* for (nc_a=0;nc_a<Wn_nrows;nc_a++){ } */}
	  printf("\n");
	  nr_j++; /* if (bget__on(mr_j,nr_a)){ } */}
	nr_b++; /* if (bget__on(mr_b,nr_a)){ } */}    
      nr_m8++; nr_m7--; if (nr_m8>=BIT8){ nr_d8++; nr_m8=0; nr_m7=7;}
      /* for (nr_a=0;nr_a<Yn_nrows;nr_a++){ } */}
    /* if (verbose>2){ } */}
 skip_xcalc:
  if (*M_X_p==NULL){ if (GLOBAL_thread_count>1){ printf(" %% Warning! *M_X_p undefined in get_xcalc; must pre-allocate to avoid thread conflict!\n");} M_X = M_handle_v_make(BITJ,ncols_bZ,nrows_bZ,NULL,(*bZ_p),NULL,NULL); *M_X_p = M_X;}
  else if (*M_X_p!=NULL){
    M_X->nrows = ncols_bZ; M_X->ncols = nrows_bZ; 
    M_X->nrows_extend = (BITJ - (M_X->nrows % BITJ)) % BITJ; M_X->mr_length = bsize(M_X->nrows);
    M_X->ncols_extend = (BITJ - (M_X->ncols % BITJ)) % BITJ; M_X->mc_length = bsize(M_X->ncols);
    for (nc_j=0;nc_j<ncols_bZ;nc_j++){ bset__on(M_X->mr_b,nc_j); bset__on(M_X->mr_j,nc_j);} M_X->rpop_b = ncols_bZ; M_X->rpop_j = ncols_bZ;
    for (nr_a=0;nr_a<M_X->mc_length;nr_a++){ M_X->mc_b[nr_a] = mr_b[nr_a]; M_X->mc_j[nr_a] = mr_j[nr_a];} M_X->cpop_b = rpop_b; M_X->cpop_j = rpop_j;
    M_X->wX = (*bZ_p);
    /* else if (M_X_p!=NULL){ } */}
  M_X->ncols_per_z = ncols_per; M_X->max_d = max; M_X->min_d = min; M_X->mlt_d = mlt; if (!isfinite(M_X->mlt_d) || M_X->mlt_d<=0){ M_X->mlt_d=1;}
  M_mxget(M_X);  
  if (verbose>2){ M_handle_printf(M_X,verbose," %% M_X: ");}
  GLOBAL_ops_count_one(tidx,(unsigned long long int)lf_YnWt->row_stride*(unsigned long long int)lf_YnWt->col_stride*(unsigned long long int)ncols_per,(unsigned long long int)lf_YnWt->row_stride*(unsigned long long int)lf_YnWt->col_stride*(unsigned long long int)ncols_per*(unsigned long long int)BIT8);
  if (verbose){ printf(" %% [finished get_xcalc] tidx %d trm_flag %d\n",tidx,trm_flag);}
  return NULL;
}

void wrap_xcalc(int *tidx,void **vpra,pthread_t *thread_in,unsigned char *mr_j,unsigned char *mr_b,unsigned char *mc_j,unsigned char *mc_b,struct L_handle *lf_YnWt,double *YnWt,int *Yn_nrows_p,int *Wn_nrows_p,struct M_handle **M_X_p,int *b_mlt_p,int trm_flag)
{
  /* calls get_xcalc;
     variable space in **vpra (should be at least size 12)
   */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  int ip=0;
  if (verbose){ printf(" %% [entering wrap_xcalc] tidx %d trm_flag %d\n",*tidx,trm_flag);}
  ip=0; vpra[ip++] = tidx; vpra[ip++] = mr_j; vpra[ip++] = mr_b; vpra[ip++] = mc_j; vpra[ip++] = mc_b; vpra[ip++] = lf_YnWt; vpra[ip++] = YnWt; vpra[ip++] = Yn_nrows_p; vpra[ip++] = Wn_nrows_p; vpra[ip++] = M_X_p; vpra[ip++] = b_mlt_p;
  switch (trm_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break; default: break; /* switch (trm_flag){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_xcalc,vpra)){ printf("Warning! cannot create thread %d in wrap_xcalc\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_xcalc(vpra);} /* must join threads later */;
  if (verbose>2){ M_handle_printf(*M_X_p,verbose," %% M_X: ");}
  if (verbose){ printf(" %% [finished wrap_xcalc] tidx %d trm_flag %d\n",*tidx,trm_flag);}
}

void *get_At_T_YnWt_ww(void *vp)
{
  int verbose=1;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_At = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Tt = (struct M_handle *)(vpra[ip++]);
  struct L_handle *lf_YnWt = (struct L_handle *)(vpra[ip++]);
  struct M_handle *M_YnWt = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_Wn = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_At_T_YnWt_uu = (struct L_handle *)(vpra[ip++]);
  struct L_handle *output_At_T_YnWt_ww = (struct L_handle *)(vpra[ip++]);
  int trm_flag = *(int *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_t = *(int *)(vpra[ip++]);
  int output_spacing_w = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_At->nrows)/* rup(M_At->nrows + M_At->nrows_extend,POPLENGTH)/POPLENGTH */;
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int ns_j=0,ns_a=0,ns_b=0,tab_t=0,tab_t_stride=0,na_j=0,na_a=0,na_b=0,tab_a=0,tab_a_stride=0,mw_j=0,mw_a=0,mw_b=0,tab_w=0,tab_w_stride=0,my_j=0,my_a=0,my_b=0,tab_y=0;
  double output_tmp=0;
  int vA=0,vS=0;
  unsigned char *Tt_tag=NULL,*At_tag=NULL;
  __m128i *wAt_tag=NULL,*wTt_tag=NULL,*wYnWt_tag=NULL,*mcat_tag=NULL,*mcat_end=NULL;
  long long int lltmp=0,n2=0;
  double output_At_T_YnWt_base=0,output_at_T_YnWt_base=0,output____T_YnWt_base=0,dtmp=0,output_At_T_YnWt_tmp=0,output_at_T_YnWt_tmp=0,output____T_YnWt_tmp_[M_Wn->nrows];
  int tab_x=0,mr=0,mx=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering get_At_T_YnWt_ww] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% Calculating output_At_T_YnWt_uu\n");}
  switch (output_spacing_t){ case SPACING_j: tab_t_stride = M_Tt->rpop_j; break; case SPACING_b: tab_t_stride = M_Tt->rpop_b; break; case SPACING_a: tab_t_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_w){ case SPACING_j: tab_w_stride = M_Wn->rpop_j; break; case SPACING_b: tab_w_stride = M_Wn->rpop_b; break; case SPACING_a: tab_w_stride = M_Wn->nrows; break; default: break; /* switch (output_spacing_w){ } */}
  if (verbose>1){ printf(" %% tab_t_stride %d tab_a_stride %d tab_w_stride %d\n",tab_t_stride,tab_a_stride,tab_w_stride);}
  output_At_T_YnWt_uu->spacing_row = output_spacing_a; output_At_T_YnWt_uu->row_stride = tab_a_stride;
  output_At_T_YnWt_uu->spacing_col = output_spacing_w; output_At_T_YnWt_uu->col_stride = tab_w_stride;
  output_At_T_YnWt_uu->spacing_lyr = output_spacing_t; output_At_T_YnWt_uu->lyr_stride = tab_t_stride;
  if (strstr(GLOBAL_skip,"At_T_YnWt_ww")){ goto skip_At_T_YnWt_ww;}
  if (verbose>2){ M_handle_printf(M_At,1," %% M_At: ");}
  if (verbose>2){ M_handle_printf(M_Tt,1," %% M_Tt: ");}
  if (verbose>1){ raprintf(lf_YnWt->lf,"double",lf_YnWt->row_stride,lf_YnWt->col_stride," %% lf_YnWt->lf: ");}
  if (verbose>1){ raprintf(D_An,"double",1,A_pcols," %% D_An: "); raprintf(a_An,"double",1,A_pcols," %% a_An: ");}
  if (verbose>1){ printf(" %% POPLENGTH %d M_At->mc_length %d M_Tt->mc_length %d M_YnWt->mc_length %d \n",POPLENGTH,M_At->mc_length,M_Tt->mc_length,M_YnWt->mc_length);}
  ns_j=0;
  while (ns_j<M_Tt->rpop_j){
    ns_a = M_Tt->m_a_[ns_j]; ns_b = M_Tt->m_b_[ns_j];
    switch (output_spacing_t){ case SPACING_j: tab_t=ns_j; break; case SPACING_b: tab_t=ns_b; break; case SPACING_a: tab_t=ns_a; break; default: break; /* switch (output_spacing_t){ } */}
    Tt_tag = (unsigned char *)(&(M_Tt->wX[ns_b*M_Tt->mc_length]));
    na_j=0;
    while (na_j<M_At->rpop_j){
      na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
      At_tag = (unsigned char *)(&(M_At->wX[na_b*M_At->mc_length]));
      mw_j=0;
      while (mw_j<M_Wn->rpop_j){
	mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=mw_j; break; case SPACING_b: tab_w=mw_b; break; case SPACING_a: tab_w=mw_a; break; default: break; /* switch (output_spacing_w){ } */}
	output_tmp=0;
	my_j=0;
	while (my_j<M_At->cpop_j){
	  my_a = M_At->n_a_[my_j]; my_b = M_At->n_b_[my_j];
	  vA = bget____(At_tag,my_a); vS = bget____(Tt_tag,my_a);
	  output_tmp += (vA - a_An[na_a/POPLENGTH])*vS* (trm_flag ? (*L2_get(lf_YnWt,mw_j,mw_b,mw_a,my_j,my_b,my_a)) : (*L2_get(lf_YnWt,my_j,my_b,my_a,mw_j,mw_b,mw_a)));
	  my_j++; /* while (my_j<M_At->cpop_j){ } */}
	if (verbose>3){ printf("(%d,%d,%d) %lf ",tab_a,tab_w,tab_t,output_tmp);}
	output_At_T_YnWt_uu->lf[tab_a + tab_w*tab_a_stride + tab_t*tab_a_stride*tab_w_stride] = output_tmp;
	mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
      if (verbose>3){ printf("\n");}
      na_j++; /* while (na_j<M_At->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
  if (verbose>1){ raprintf(output_At_T_YnWt_uu->lf,"double",tab_a_stride,tab_t_stride*tab_w_stride," %% output_At_T_YnWt_uu->lf: ");}
  if (verbose>1){ printf(" %% Calculating output_At_T_YnWt_ww\n");}
  switch (output_spacing_t){ case SPACING_j: tab_t_stride = M_Tt->rpop_j; break; case SPACING_b: tab_t_stride = M_Tt->rpop_b; break; case SPACING_a: tab_t_stride = M_Tt->nrows; break; default: break; /* switch (output_spacing_t){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_At->rpop_j; break; case SPACING_b: tab_a_stride = M_At->rpop_b; break; case SPACING_a: tab_a_stride = M_At->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_w){ case SPACING_j: tab_w_stride = M_Wn->rpop_j; break; case SPACING_b: tab_w_stride = M_Wn->rpop_b; break; case SPACING_a: tab_w_stride = M_Wn->nrows; break; default: break; /* switch (output_spacing_w){ } */}
  output_At_T_YnWt_ww->spacing_row = output_spacing_a; output_At_T_YnWt_ww->row_stride = tab_a_stride;
  output_At_T_YnWt_ww->spacing_col = output_spacing_w; output_At_T_YnWt_ww->col_stride = tab_w_stride;
  output_At_T_YnWt_ww->spacing_lyr = output_spacing_t; output_At_T_YnWt_ww->lyr_stride = tab_t_stride;
  ns_j=0;
  while (ns_j<M_Tt->rpop_j){
    ns_a = M_Tt->m_a_[ns_j]; ns_b = M_Tt->m_b_[ns_j];
    switch (output_spacing_t){ case SPACING_j: tab_t=ns_j; break; case SPACING_b: tab_t=ns_b; break; case SPACING_a: tab_t=ns_a; break; default: break; /* switch (output_spacing_t){ } */}
    wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
    wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
    wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
    mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
    lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
    output____T_YnWt_base = dtmp*M_YnWt->min_d*M_YnWt->mlt_d;
    memset(output____T_YnWt_tmp_,0,M_Wn->rpop_j*sizeof(double));
    mw_j=0;
    while (mw_j<M_Wn->rpop_j){
      mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
      switch (output_spacing_w){ case SPACING_j: tab_w=mw_j; break; case SPACING_b: tab_w=mw_b; break; case SPACING_a: tab_w=mw_a; break; default: break; /* switch (output_spacing_w){ } */}
      output_at_T_YnWt_tmp = 0; 
      mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
      for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	wAt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
	wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
	wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end); dtmp = lltmp*1.0;
	output_at_T_YnWt_tmp += n2*dtmp;
	n2*=2;
	/* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
      output____T_YnWt_tmp_[mw_j] = output_at_T_YnWt_tmp;
      mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
    na_j=0;
    while (na_j<M_At->rpop_j){
      na_a = M_At->m_a_[na_j]; na_b = M_At->m_b_[na_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=na_j; break; case SPACING_b: tab_a=na_b; break; case SPACING_a: tab_a=na_a; break; default: break; /* switch (output_spacing_a){ } */}
      wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
      wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
      wYnWt_tag = (__m128i*)&(M_YnWt->wX[0/* start */*M_YnWt->mc_length]);
      mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
      lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
      output_At_T_YnWt_base = lltmp*M_YnWt->min_d*M_YnWt->mlt_d;
      mw_j=0;
      while (mw_j<M_Wn->rpop_j){
	mw_a = M_Wn->m_a_[mw_j]; mw_b = M_Wn->m_b_[mw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=mw_j; break; case SPACING_b: tab_w=mw_b; break; case SPACING_a: tab_w=mw_a; break; default: break; /* switch (output_spacing_w){ } */}
	tab_x = tab_a + tab_w*tab_a_stride + tab_t*tab_a_stride*tab_w_stride;
	output_At_T_YnWt_tmp = 0;
	mr=M_YnWt->ncols_per_z*mw_j/* spacing_j */; n2=1;
	for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){
	  wAt_tag = (__m128i*)&(M_At->wX[na_b*M_At->mc_length]);
	  wTt_tag = (__m128i*)&(M_Tt->wX[ns_b*M_Tt->mc_length]);
	  wYnWt_tag = (__m128i*)&(M_YnWt->wX[(mr+mx)*M_YnWt->mc_length]);
	  mcat_tag = (__m128i*)&(M_At->mc_j[0]); mcat_end = (__m128i*)&(M_At->mc_j[M_At->mc_length]);
	  lltmp = popcount_pmpm0(&wAt_tag,&wTt_tag,&wYnWt_tag,&mcat_tag,&mcat_end);
	  output_At_T_YnWt_tmp += n2*lltmp;
	  n2*=2;
	  /* for (mx=M_YnWt->ncols_per_z-1;mx>0;mx--){ } */}
	output_at_T_YnWt_base = output____T_YnWt_base*a_An[na_a/POPLENGTH];
	output_at_T_YnWt_tmp = output____T_YnWt_tmp_[mw_j]*a_An[na_a/POPLENGTH];
	output_At_T_YnWt_ww->lf[tab_x] = (double)(output_At_T_YnWt_base + output_At_T_YnWt_tmp)/M_YnWt->mlt_d - (double)(output_at_T_YnWt_base + output_at_T_YnWt_tmp)/M_YnWt->mlt_d;
	mw_j++; /* while (mw_j<M_Wn->rpop_j){ } */}
      na_j++; /* while (na_j<M_At->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
  if (verbose>1){ raprintf(output_At_T_YnWt_ww->lf,"double",tab_a_stride,tab_t_stride*tab_w_stride," %% output_At_T_YnWt_ww->lf: ");}
  if (verbose){
    ns_j=0;
    while (ns_j<M_Tt->rpop_j){
      ns_a = M_Tt->m_a_[ns_j]; ns_b = M_Tt->m_b_[ns_j];
      switch (output_spacing_t){ case SPACING_j: tab_t=ns_j; break; case SPACING_b: tab_t=ns_b; break; case SPACING_a: tab_t=ns_a; break; default: break; /* switch (output_spacing_t){ } */}
      if (verbose>1){
	sprintf(tmpchar," %% output_At_T_YnWt_uu->lf[%d]: ",tab_t);
	raprintf(&(output_At_T_YnWt_uu->lf[tab_t*tab_a_stride*tab_w_stride]),"double",tab_a_stride,tab_w_stride,tmpchar);
	sprintf(tmpchar," %% output_At_T_YnWt_ww->lf[%d]: ",tab_t);
	raprintf(&(output_At_T_YnWt_ww->lf[tab_t*tab_a_stride*tab_w_stride]),"double",tab_a_stride,tab_w_stride,tmpchar);
	/* if (verbose>1){ } */}
      printf(" %% Difference: %0.16f\n",dra_diff(&(output_At_T_YnWt_uu->lf[tab_t*tab_a_stride*tab_w_stride]),&(output_At_T_YnWt_ww->lf[tab_t*tab_a_stride*tab_w_stride]),tab_a_stride*tab_w_stride,1));
      ns_j++; /* while (ns_j<M_Tt->rpop_j){ } */}
    /* if (verbose){ } */}
  GLOBAL_ops_count_one(tidx,M_Tt->rpop_j*M_At->rpop_j*M_Wn->rpop_j,M_Tt->rpop_j*M_At->rpop_j*M_Wn->rpop_j*M_YnWt->ncols_per_z*M_YnWt->mc_length*BIT8);
 skip_At_T_YnWt_ww:
  if (verbose>1){ printf(" %% [finished get_At_T_YnWt_ww] tidx %d\n",tidx);}
  return NULL;
}

int wrap_At_T_YnWt_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_t,int spacing_w,struct M_handle *M_At,struct M_handle *M_Tt,struct L_handle *lf_YnWt,struct M_handle *M_YnWt,struct M_handle *M_Wn,double *A_ajdk,int trm_flag,struct L_handle **output_At_T_YnWt_uu_p,struct L_handle **output_At_T_YnWt_ww_p)
{
  /* calls get_At_T_YnWt_ww;
     variable space in **vpra (should be at least size 13)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_t=0,length_w=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_At_T_YnWt_ww__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_At,verbose," %% M_At: ");}
  if (verbose){ M_handle_printf(M_Tt,verbose," %% M_Tt: ");}
  if (verbose){ M_handle_printf(M_Wn,verbose," %% M_Wn: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_At->rpop_j; break; case SPACING_b: length_a = M_At->rpop_b; break; case SPACING_a: length_a = M_At->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_t){ case SPACING_j: length_t = M_Tt->rpop_j; break; case SPACING_b: length_t = M_Tt->rpop_b; break; case SPACING_a: length_t = M_Tt->nrows; break; default: break; /* switch (spacing_t){ } */}
  switch (spacing_w){ case SPACING_j: length_w = M_Wn->rpop_j; break; case SPACING_b: length_w = M_Wn->rpop_b; break; case SPACING_a: length_w = M_Wn->nrows; break; default: break; /* switch (spacing_w){ } */}
  length = length_a*length_t*length_w; if (verbose){ printf(" %% length %llu*%llu*%llu=%llu\n",length_a,length_t,length_w,length);}
  if (verbose>2){ bprintf(M_At->mr_b,M_At->bitj,1,M_At->nrows," %% M_At->mr_b: ");}
  if (verbose>2){ bprintf(M_At->mr_j,M_At->bitj,1,M_At->nrows," %% M_At->mr_j: ");}
  if (verbose>2){ bprintf(M_At->mc_b,M_At->bitj,1,M_At->ncols," %% M_At->mc_b: ");}
  if (verbose>2){ bprintf(M_At->mc_j,M_At->bitj,1,M_At->ncols," %% M_At->mc_j: ");}
  if (verbose>2){ bprintf(M_Wn->mr_b,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_b: ");}
  if (verbose>2){ bprintf(M_Wn->mr_j,M_Wn->bitj,1,M_Wn->nrows," %% M_Wn->mr_j: ");}
  if (verbose>2){ bprintf(M_Wn->mc_b,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_b: ");}
  if (verbose>2){ bprintf(M_Wn->mc_j,M_Wn->bitj,1,M_Wn->ncols," %% M_Wn->mc_j: ");}
  if (*output_At_T_YnWt_uu_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_At_T_YnWt_uu_p = L_handle_make(length);}
  if ((*output_At_T_YnWt_uu_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_At_T_YnWt_ww__run\n",(*output_At_T_YnWt_uu_p)->length,length);}
  memset((*output_At_T_YnWt_uu_p)->lf,0,length*sizeof(double));
  if (*output_At_T_YnWt_ww_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_At_T_YnWt_ww_p = L_handle_make(length);}
  if ((*output_At_T_YnWt_ww_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_At_T_YnWt_ww__run\n",(*output_At_T_YnWt_ww_p)->length,length);}
  memset((*output_At_T_YnWt_ww_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_At; vpra[ip++] = M_Tt; vpra[ip++] = lf_YnWt; vpra[ip++] = M_YnWt; vpra[ip++] = M_Wn; vpra[ip++] = A_ajdk; vpra[ip++] = *output_At_T_YnWt_uu_p; vpra[ip++] = *output_At_T_YnWt_ww_p;
  switch (trm_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break; default: break; /* switch (trm_flag){ } */}
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_t){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_t){ } */}
  switch (spacing_w){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_w){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_At_T_YnWt_ww,vpra)){ printf("Warning! cannot create thread %d in wrap_At_T_YnWt_ww__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_At_T_YnWt_ww(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_At_T_YnWt_ww__run] tidx %d\n",*tidx);}
  return length;
}

void wrap_At_T_YnWt_ww_test()
{
  /* test for errors with input file: At_T_YnWt_ww_error.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AnZt=NULL,**lf_YnWt=NULL;
  struct M_handle **M_YnWt=NULL,**M_AnZt=NULL;
  struct L_handle **lf_At_T_YnWt_ww=NULL,**lf_At_T_AnZt_ww=NULL; int *length_At_T_YnWt_ww=NULL,*length_At_T_AnZt_ww=NULL;
  struct L_handle **lf_At_T_YnWt_uu=NULL,**lf_At_T_AnZt_uu=NULL; int *length_At_T_YnWt_uu=NULL,*length_At_T_AnZt_uu=NULL;
  struct M_handle **M_WnYt=NULL,**M_ZnAt=NULL;
  struct L_handle **lf_Zt_S_WnYt_ww=NULL,**lf_Zt_S_ZnAt_ww=NULL; int *length_Zt_S_WnYt_ww=NULL,*length_Zt_S_ZnAt_ww=NULL;
  struct L_handle **lf_Zt_S_WnYt_uu=NULL,**lf_Zt_S_ZnAt_uu=NULL; int *length_Zt_S_WnYt_uu=NULL,*length_Zt_S_ZnAt_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_At_T_YnWt_ww_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  lf_AnZt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_YnWt = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_AnZt[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    lf_YnWt[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_YnWt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  M_AnZt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_YnWt[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Yn[nb]->nrows,M_Wn[nb]->nrows);
    M_AnZt[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_An[nb]->nrows,M_Zn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_At_T_YnWt_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_YnWt_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_At_T_AnZt_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_AnZt_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_At_T_YnWt_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)M_Wn[nb]->nrows);
    lf_At_T_AnZt_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_At_T_YnWt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_YnWt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_At_T_AnZt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_At_T_AnZt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_At_T_YnWt_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)M_Wn[nb]->nrows);
    lf_At_T_AnZt_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->ncols*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)M_Zn[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_WnYt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  M_ZnAt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_WnYt[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Wn[nb]->nrows,M_Yn[nb]->nrows);
    M_ZnAt[nb] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Zn[nb]->nrows,M_An[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_Zt_S_WnYt_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Zt_S_WnYt_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Zt_S_ZnAt_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Zt_S_ZnAt_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_Zt_S_WnYt_ww[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_St[nb]->nrows*(unsigned long long int)M_Yn[nb]->nrows);
    lf_Zt_S_ZnAt_ww[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_St[nb]->nrows*(unsigned long long int)M_An[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_Zt_S_WnYt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Zt_S_WnYt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Zt_S_ZnAt_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Zt_S_ZnAt_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_Zt_S_WnYt_uu[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_St[nb]->nrows*(unsigned long long int)M_Yn[nb]->nrows);
    lf_Zt_S_ZnAt_uu[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->ncols*(unsigned long long int)M_St[nb]->nrows*(unsigned long long int)M_An[nb]->nrows);
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
  	  wrap_AnZt_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,n_spacing_A,M_Yn[nb],M_Wn[nb],Y_ajdk,lf_Yn_ajdk[nb],lf_Wn_ajdk[nb],&(lf_YnWt[nb]));
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
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% xcalc trm_0: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
	  wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Wn[nb]->mr_j,M_Wn[nb]->mr_b,M_Yn[nb]->mr_j,M_Yn[nb]->mr_b,lf_YnWt[nb],lf_YnWt[nb]->lf,&(M_Wn[nb]->nrows),&(M_Yn[nb]->nrows),&(M_WnYt[nb]),&(GLOBAL_B_MLT),(addressable_1));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Zn[nb]->mr_j,M_Zn[nb]->mr_b,M_An[nb]->mr_j,M_An[nb]->mr_b,lf_AnZt[nb],lf_AnZt[nb]->lf,&(M_Zn[nb]->nrows),&(M_An[nb]->nrows),&(M_ZnAt[nb]),&(GLOBAL_B_MLT),(addressable_1));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% xcalc trm_1: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
	  length_At_T_YnWt_ww[nb] = wrap_At_T_YnWt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],lf_YnWt[nb],M_YnWt[nb],M_Wn[nb],A_ajdk,(addressable_0),&(lf_At_T_YnWt_uu[nb]),&(lf_At_T_YnWt_ww[nb]));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
	  length_At_T_AnZt_ww[nb] = wrap_At_T_YnWt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_At[nb],M_Tt[nb],lf_AnZt[nb],M_AnZt[nb],M_Zn[nb],A_ajdk,(addressable_0),&(lf_At_T_AnZt_uu[nb]),&(lf_At_T_AnZt_ww[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% At_T_YnWt At_T_AnZt: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic(); 
	  length_Zt_S_WnYt_ww[nb] = wrap_At_T_YnWt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_Zt[nb],M_St[nb],lf_YnWt[nb],M_WnYt[nb],M_Yn[nb],A_ajdk,(addressable_1),&(lf_Zt_S_WnYt_uu[nb]),&(lf_Zt_S_WnYt_ww[nb]));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic(); 
	  length_Zt_S_ZnAt_ww[nb] = wrap_At_T_YnWt_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_Zt[nb],M_St[nb],lf_AnZt[nb],M_ZnAt[nb],M_An[nb],A_ajdk,(addressable_1),&(lf_Zt_S_ZnAt_uu[nb]),&(lf_Zt_S_ZnAt_ww[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc(); 
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% Zt_S_WnYt Zt_S_ZnAt: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_At_T_YnWt_ww_test]\n");}
}

void *get_An_ZtSWn_ww(void *vp)
{
  int verbose=1;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  struct M_handle *M_St = (struct M_handle *)(vpra[ip++]);
  struct L_handle *lf_ZtSWn = (struct L_handle *)(vpra[ip++]);
  struct M_handle **M_ZtSWn_ = (struct M_handle **)(vpra[ip++]); struct M_handle *M_ZtSWn=NULL;
  struct M_handle *M_Wt = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_An_ZtSWn_uu = (struct L_handle *)(vpra[ip++]);
  struct L_handle *output_An_ZtSWn_ww = (struct L_handle *)(vpra[ip++]);
  int trm_flag = *(int *)(vpra[ip++]);
  int output_spacing_a = *(int *)(vpra[ip++]);
  int output_spacing_s = *(int *)(vpra[ip++]);
  int output_spacing_w = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols);
  double *D_An = (double *)&(A_ajdk[0 + AJDK_0_1*A_pcols]);
  double *a_An = (double *)&(A_ajdk[0 + AJDK_1_0*A_pcols]);
  int ns_j=0,ns_a=0,ns_b=0,tab_s=0,tab_s_stride=0,ma_j=0,ma_a=0,ma_b=0,tab_a=0,tab_a_stride=0,nz_j=0,nz_a=0,nz_b=0,tab_z=0,nw_j=0,nw_a=0,nw_b=0,tab_w=0,tab_w_stride=0;
  double output_tmp=0;
  int vA=0;
  unsigned char *An_tag=NULL;
  __m128i *wAn_tag=NULL,*wZtSWn_tag=NULL,*mcan_tag=NULL,*mcan_end=NULL;
  long long int n2=0;
  double output_An_ZtSWn_base=0,output_an_ZtSWn_base=0,output____ZtSWn_base=0,*dinp=NULL,dtmp=0,output_An_ZtSWn_tmp=0,output_an_ZtSWn_tmp=0,output____ZtSWn_tmp_[M_Wt->nrows];
  int tab_x=0,mr=0,mx=0;
  char tmpchar[FNAMESIZE];
  if (verbose>1){ printf(" %% [entering get_An_ZtSWn_ww] tidx %d\n",tidx);}
  if (verbose>1){ printf(" %% Calculating output_An_ZtSWn_uu\n");}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_w){ case SPACING_j: tab_w_stride = M_Wt->rpop_j; break; case SPACING_b: tab_w_stride = M_Wt->rpop_b; break; case SPACING_a: tab_w_stride = M_Wt->nrows; break; default: break; /* switch (output_spacing_w){ } */}
  if (verbose>1){ printf(" %% tab_s_stride %d tab_a_stride %d tab_w_stride %d\n",tab_s_stride,tab_a_stride,tab_w_stride);}
  output_An_ZtSWn_uu->spacing_row = output_spacing_a; output_An_ZtSWn_uu->row_stride = tab_a_stride;
  output_An_ZtSWn_uu->spacing_col = output_spacing_w; output_An_ZtSWn_uu->col_stride = tab_w_stride;
  output_An_ZtSWn_uu->spacing_lyr = output_spacing_s; output_An_ZtSWn_uu->lyr_stride = tab_s_stride;
  if (strstr(GLOBAL_skip,"An_ZtSWn_ww")){ goto skip_An_ZtSWn_ww;}
  if (verbose>2){ M_handle_printf(M_An,1," %% M_An: ");}
  if (verbose>2){ M_handle_printf(M_St,1," %% M_St: ");}
  if (verbose>1){ raprintf(lf_ZtSWn->lf,"double",lf_ZtSWn->row_stride,lf_ZtSWn->col_stride*lf_ZtSWn->lyr_stride," %% lf_ZtSWn->lf: ");}
  if (verbose>1){ raprintf(D_An,"double",1,A_pcols," %% D_An: "); raprintf(a_An,"double",1,A_pcols," %% a_An: ");}
  if (verbose>1){ printf(" %% POPLENGTH %d M_An->mc_length %d M_St->mc_length %d\n",POPLENGTH,M_An->mc_length,M_St->mc_length);}
  if (verbose>1){ printf(" %% "); for (ns_j=0;ns_j<M_St->rpop_j;ns_j++){ printf("M_ZtSWn_[ns_j==%d]->mc_length = %d ",ns_j,M_ZtSWn_[ns_j]->mc_length);} printf("\n");}
  ns_j=0;
  while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    ma_j=0;
    while (ma_j<M_An->rpop_j){
      ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=ma_j; break; case SPACING_b: tab_a=ma_b; break; case SPACING_a: tab_a=ma_a; break; default: break; /* switch (output_spacing_a){ } */}
      An_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length]));
      nw_j=0;
      while (nw_j<M_Wt->rpop_j){
	nw_a = M_Wt->m_a_[nw_j]; nw_b = M_Wt->m_b_[nw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=nw_j; break; case SPACING_b: tab_w=nw_b; break; case SPACING_a: tab_w=nw_a; break; default: break; /* switch (output_spacing_w){ } */}
	output_tmp=0;
	nz_j=0;
	while (nz_j<M_An->cpop_j){
	  nz_a = M_An->n_a_[nz_j]; nz_b = M_An->n_b_[nz_j];
	  vA = bget____(An_tag,nz_a);
	  output_tmp += (vA - a_An[nz_a/POPLENGTH])*D_An[nz_a/POPLENGTH] * (trm_flag ? (*L3_get(lf_ZtSWn,nw_j,nw_b,nw_a,nz_j,nz_b,nz_a,ns_j,ns_b,ns_a)) : (*L3_get(lf_ZtSWn,nz_j,nz_b,nz_a,nw_j,nw_b,nw_a,ns_j,ns_b,ns_a)));
	  nz_j++; /* while (nz_j<M_An->cpop_j){ } */}
	if (verbose>3){ printf("(%d,%d,%d) %lf ",tab_a,tab_w,tab_s,output_tmp);}
	output_An_ZtSWn_uu->lf[tab_a + tab_w*tab_a_stride + tab_s*tab_a_stride*tab_w_stride] = output_tmp;
	nw_j++; /* while (nw_j<M_Wt->rpop_j){ } */}
      if (verbose>3){ printf("\n");}
      ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  if (verbose>1){ raprintf(output_An_ZtSWn_uu->lf,"double",tab_a_stride,tab_w_stride*tab_s_stride," %% output_An_ZtSWn_uu->lf: ");}
  if (verbose>1){ printf(" %% Calculating output_An_ZtSWn_ww\n");}
  switch (output_spacing_s){ case SPACING_j: tab_s_stride = M_St->rpop_j; break; case SPACING_b: tab_s_stride = M_St->rpop_b; break; case SPACING_a: tab_s_stride = M_St->nrows; break; default: break; /* switch (output_spacing_s){ } */}
  switch (output_spacing_a){ case SPACING_j: tab_a_stride = M_An->rpop_j; break; case SPACING_b: tab_a_stride = M_An->rpop_b; break; case SPACING_a: tab_a_stride = M_An->nrows; break; default: break; /* switch (output_spacing_a){ } */}
  switch (output_spacing_w){ case SPACING_j: tab_w_stride = M_Wt->rpop_j; break; case SPACING_b: tab_w_stride = M_Wt->rpop_b; break; case SPACING_a: tab_w_stride = M_Wt->nrows; break; default: break; /* switch (output_spacing_w){ } */}
  output_An_ZtSWn_ww->spacing_row = output_spacing_a; output_An_ZtSWn_ww->row_stride = tab_a_stride;
  output_An_ZtSWn_ww->spacing_col = output_spacing_w; output_An_ZtSWn_ww->col_stride = tab_w_stride;
  output_An_ZtSWn_ww->spacing_lyr = output_spacing_s; output_An_ZtSWn_ww->lyr_stride = tab_s_stride;
  ns_j=0;
  while (ns_j<M_St->rpop_j){
    ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j]; M_ZtSWn = M_ZtSWn_[ns_j];
    switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
    wAn_tag = (__m128i*)&(M_ZtSWn->wX[0/* start */*M_ZtSWn->mc_length]);
    wZtSWn_tag = (__m128i*)&(M_ZtSWn->wX[0/* start */*M_ZtSWn->mc_length]);
    mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
    dinp = &(A_ajdk[0+AJDK_1_1*A_pcols]);
    dtmp = popcount_pm0_lf(&wAn_tag,&wZtSWn_tag,&mcan_tag,&mcan_end,&dinp);
    output____ZtSWn_base = dtmp*M_ZtSWn->min_d*M_ZtSWn->mlt_d;
    memset(output____ZtSWn_tmp_,0,M_Wt->rpop_j*sizeof(double));
    nw_j=0;
    while (nw_j<M_Wt->rpop_j){
      nw_a = M_Wt->m_a_[nw_j]; nw_b = M_Wt->m_b_[nw_j];
      switch (output_spacing_w){ case SPACING_j: tab_w=nw_j; break; case SPACING_b: tab_w=nw_b; break; case SPACING_a: tab_w=nw_a; break; default: break; /* switch (output_spacing_w){ } */}
      output_an_ZtSWn_tmp = 0; 
      mr=M_ZtSWn->ncols_per_z*nw_j/* spacing_j */; n2=1;
      for (mx=M_ZtSWn->ncols_per_z-1;mx>0;mx--){
	wAn_tag = (__m128i*)&(M_ZtSWn->wX[0/* start */*M_ZtSWn->mc_length]);
	wZtSWn_tag = (__m128i*)&(M_ZtSWn->wX[(mr+mx)*M_ZtSWn->mc_length]);
	mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	dinp = &(A_ajdk[0+AJDK_1_1*A_pcols]);
	dtmp = popcount_pm0_lf(&wAn_tag,&wZtSWn_tag,&mcan_tag,&mcan_end,&dinp);
	output_an_ZtSWn_tmp += n2*dtmp;
	n2*=2;
	/* for (mx=M_ZtSWn->ncols_per_z-1;mx>0;mx--){ } */}
      output____ZtSWn_tmp_[nw_j] = output_an_ZtSWn_tmp;
      nw_j++; /* while (nw_j<M_Wt->rpop_j){ } */}
    ma_j=0;
    while (ma_j<M_An->rpop_j){
      ma_a = M_An->m_a_[ma_j]; ma_b = M_An->m_b_[ma_j];
      switch (output_spacing_a){ case SPACING_j: tab_a=ma_j; break; case SPACING_b: tab_a=ma_b; break; case SPACING_a: tab_a=ma_a; break; default: break; /* switch (output_spacing_a){ } */}
      wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
      wZtSWn_tag = (__m128i*)&(M_ZtSWn->wX[0/* start */*M_ZtSWn->mc_length]);
      mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
      dinp = &(A_ajdk[0+AJDK_0_1*A_pcols]);
      dtmp = popcount_pm0_lf(&wAn_tag,&wZtSWn_tag,&mcan_tag,&mcan_end,&dinp);
      output_An_ZtSWn_base = dtmp*M_ZtSWn->min_d*M_ZtSWn->mlt_d;
      nw_j=0;
      while (nw_j<M_Wt->rpop_j){
	nw_a = M_Wt->m_a_[nw_j]; nw_b = M_Wt->m_b_[nw_j];
	switch (output_spacing_w){ case SPACING_j: tab_w=nw_j; break; case SPACING_b: tab_w=nw_b; break; case SPACING_a: tab_w=nw_a; break; default: break; /* switch (output_spacing_w){ } */}
	tab_x = tab_a + tab_w*tab_a_stride + tab_s*tab_a_stride*tab_w_stride;
	output_An_ZtSWn_tmp = 0;
	mr=M_ZtSWn->ncols_per_z*nw_j/* spacing_j */; n2=1;
	for (mx=M_ZtSWn->ncols_per_z-1;mx>0;mx--){
	  wAn_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	  wZtSWn_tag = (__m128i*)&(M_ZtSWn->wX[(mr+mx)*M_ZtSWn->mc_length]);
	  mcan_tag = (__m128i*)&(M_An->mc_j[0]); mcan_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	  dinp = &(A_ajdk[0+AJDK_0_1*A_pcols]);
	  dtmp = popcount_pm0_lf(&wAn_tag,&wZtSWn_tag,&mcan_tag,&mcan_end,&dinp);
	  output_An_ZtSWn_tmp += n2*dtmp;
	  n2*=2;
	  /* for (mx=M_ZtSWn->ncols_per_z-1;mx>0;mx--){ } */}
	output_an_ZtSWn_base = output____ZtSWn_base; output_an_ZtSWn_tmp = output____ZtSWn_tmp_[nw_j];
	output_An_ZtSWn_ww->lf[tab_x] = (double)(output_An_ZtSWn_base + output_An_ZtSWn_tmp)/M_ZtSWn->mlt_d - (double)(output_an_ZtSWn_base + output_an_ZtSWn_tmp)/M_ZtSWn->mlt_d;
	nw_j++; /* while (nw_j<M_Wt->rpop_j){ } */}
      ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
    ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
  if (verbose>1){ raprintf(output_An_ZtSWn_ww->lf,"double",tab_a_stride,tab_s_stride*tab_w_stride," %% output_An_ZtSWn_ww->lf: ");}
  if (verbose){
    ns_j=0;
    while (ns_j<M_St->rpop_j){
      ns_a = M_St->m_a_[ns_j]; ns_b = M_St->m_b_[ns_j];
      switch (output_spacing_s){ case SPACING_j: tab_s=ns_j; break; case SPACING_b: tab_s=ns_b; break; case SPACING_a: tab_s=ns_a; break; default: break; /* switch (output_spacing_s){ } */}
      if (verbose>1){
	sprintf(tmpchar," %% output_An_ZtSWn_uu->lf[%d]: ",tab_s);
	raprintf(&(output_An_ZtSWn_uu->lf[tab_s*tab_a_stride*tab_w_stride]),"double",tab_a_stride,tab_w_stride,tmpchar);
	sprintf(tmpchar," %% output_An_ZtSWn_ww->lf[%d]: ",tab_s);
	raprintf(&(output_An_ZtSWn_ww->lf[tab_s*tab_a_stride*tab_w_stride]),"double",tab_a_stride,tab_w_stride,tmpchar);
	/* if (verbose>1){ } */}
      printf(" %% Difference: %0.16f\n",dra_diff(&(output_An_ZtSWn_uu->lf[tab_s*tab_a_stride*tab_w_stride]),&(output_An_ZtSWn_ww->lf[tab_s*tab_a_stride*tab_w_stride]),tab_a_stride*tab_w_stride,1));
      ns_j++; /* while (ns_j<M_St->rpop_j){ } */}
    /* if (verbose){ } */}
  GLOBAL_ops_count_one(tidx,M_St->rpop_j*M_An->rpop_j*M_Wt->rpop_j,M_St->rpop_j*M_An->rpop_j*M_Wt->rpop_j*M_ZtSWn->ncols_per_z*M_ZtSWn->mc_length*BIT8);
 skip_An_ZtSWn_ww:
  if (verbose>1){ printf(" %% [finished get_An_ZtSWn_ww] tidx %d\n",tidx);}
  return NULL;
}

int wrap_An_ZtSWn_ww__run(int *tidx,void **vpra,pthread_t *thread_in,int spacing_a,int spacing_s,int spacing_w,struct M_handle *M_An,struct M_handle *M_St,struct L_handle *lf_ZtSWn,struct M_handle **M_ZtSWn_,struct M_handle *M_Wt,double *A_ajdk,int trm_flag,struct L_handle **output_An_ZtSWn_uu_p,struct L_handle **output_An_ZtSWn_ww_p)
{
  /* calls get_An_ZtSWn_ww;
     variable space in **vpra (should be at least size 13)
   */
  int verbose=0;
  unsigned long long int length_a=0,length_s=0,length_w=0,length=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_An_ZtSWn_ww__run] tidx %d\n",*tidx);}
  if (verbose){ M_handle_printf(M_An,verbose," %% M_An: ");}
  if (verbose){ M_handle_printf(M_St,verbose," %% M_St: ");}
  if (verbose){ M_handle_printf(M_Wt,verbose," %% M_Wt: ");}
  switch (spacing_a){ case SPACING_j: length_a = M_An->rpop_j; break; case SPACING_b: length_a = M_An->rpop_b; break; case SPACING_a: length_a = M_An->nrows; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: length_s = M_St->rpop_j; break; case SPACING_b: length_s = M_St->rpop_b; break; case SPACING_a: length_s = M_St->nrows; break; default: break; /* switch (spacing_s){ } */}
  switch (spacing_w){ case SPACING_j: length_w = M_Wt->rpop_j; break; case SPACING_b: length_w = M_Wt->rpop_b; break; case SPACING_a: length_w = M_Wt->nrows; break; default: break; /* switch (spacing_w){ } */}
  length = length_a*length_s*length_w; if (verbose){ printf(" %% length %llu*%llu*%llu=%llu\n",length_a,length_s,length_w,length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  if (verbose>2){ bprintf(M_Wt->mr_b,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_b: ");}
  if (verbose>2){ bprintf(M_Wt->mr_j,M_Wt->bitj,1,M_Wt->nrows," %% M_Wt->mr_j: ");}
  if (verbose>2){ bprintf(M_Wt->mc_b,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_b: ");}
  if (verbose>2){ bprintf(M_Wt->mc_j,M_Wt->bitj,1,M_Wt->ncols," %% M_Wt->mc_j: ");}
  if (*output_An_ZtSWn_uu_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_ZtSWn_uu_p = L_handle_make(length);}
  if ((*output_An_ZtSWn_uu_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_ZtSWn_ww__run\n",(*output_An_ZtSWn_uu_p)->length,length);}
  memset((*output_An_ZtSWn_uu_p)->lf,0,length*sizeof(double));
  if (*output_An_ZtSWn_ww_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_ZtSWn_ww_p = L_handle_make(length);}
  if ((*output_An_ZtSWn_ww_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_ZtSWn_ww__run\n",(*output_An_ZtSWn_ww_p)->length,length);}
  memset((*output_An_ZtSWn_ww_p)->lf,0,length*sizeof(double));
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = M_St; vpra[ip++] = lf_ZtSWn; vpra[ip++] = M_ZtSWn_; vpra[ip++] = M_Wt; vpra[ip++] = A_ajdk; vpra[ip++] = *output_An_ZtSWn_uu_p; vpra[ip++] = *output_An_ZtSWn_ww_p;
  switch (trm_flag){ case 0: vpra[ip++] = &addressable_0; break; case 1: vpra[ip++] = &addressable_1; break; default: break; /* switch (trm_flag){ } */}
  switch (spacing_a){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_a){ } */}
  switch (spacing_s){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_s){ } */}
  switch (spacing_w){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (spacing_w){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_ZtSWn_ww,vpra)){ printf("Warning! cannot create thread %d in wrap_An_ZtSWn_ww__run\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_ZtSWn_ww(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_An_ZtSWn_ww__run] tidx %d\n",*tidx);}
  return length;
}

void wrap_An_ZtSWn_ww_test()
{
  /* test for errors with input file: An_ZtSWn_ww_error.in ;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  struct L_handle **lf_AtTn=NULL,**lf_YtTn=NULL,**lf_ZtSn=NULL,**lf_WtSn=NULL;
  struct L_handle **lf_ZtSWn=NULL,**lf_AtTYn=NULL;
  int ns_j=0,ns_b=0,ns_a=0;
  struct M_handle ***M_ZtSWn_=NULL,***M_AtTYn_=NULL;
  struct L_handle **lf_An_ZtSWn_ww=NULL,**lf_An_AtTYn_ww=NULL; int *length_An_ZtSWn_ww=NULL,*length_An_AtTYn_ww=NULL;
  struct L_handle **lf_An_ZtSWn_uu=NULL,**lf_An_AtTYn_uu=NULL; int *length_An_ZtSWn_uu=NULL,*length_An_AtTYn_uu=NULL;
  struct M_handle ***M_WtSZn_=NULL,***M_YtTAn_=NULL;
  struct L_handle **lf_Yn_WtSZn_ww=NULL,**lf_Yn_YtTAn_ww=NULL; int *length_Yn_WtSZn_ww=NULL,*length_Yn_YtTAn_ww=NULL;
  struct L_handle **lf_Yn_WtSZn_uu=NULL,**lf_Yn_YtTAn_uu=NULL; int *length_Yn_WtSZn_uu=NULL,*length_Yn_YtTAn_uu=NULL;
  if (verbose){ printf(" %% [entering wrap_An_ZtSWn_ww_test]\n");}
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
  lf_ZtSWn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_AtTYn = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_ZtSWn[nb] = L_handle_make((unsigned long long int)M_Zt[nb]->nrows*(unsigned long long int)M_Wt[nb]->nrows*(unsigned long long int)M_St[nb]->nrows);
    lf_AtTYn[nb] = L_handle_make((unsigned long long int)M_At[nb]->nrows*(unsigned long long int)M_Yt[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_ZtSWn_ = (struct M_handle ***)wkspace_all0c(sizeof(struct M_handle **)*nbins);
  M_AtTYn_ = (struct M_handle ***)wkspace_all0c(sizeof(struct M_handle **)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_ZtSWn_[nb] = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*M_St[nb]->rpop_j);
    M_AtTYn_[nb] = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*M_Tt[nb]->rpop_j);
    for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){
      M_ZtSWn_[nb][ns_j] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Zt[nb]->nrows,M_Wt[nb]->nrows);
      M_AtTYn_[nb][ns_j] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_At[nb]->nrows,M_Yt[nb]->nrows);
      /* for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){ } */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_An_ZtSWn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_ZtSWn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_An_AtTYn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_AtTYn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_An_ZtSWn_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_St[nb]->nrows*(unsigned long long int)M_Wt[nb]->nrows);
    lf_An_AtTYn_ww[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)M_Yt[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_An_ZtSWn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_ZtSWn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_An_AtTYn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_An_AtTYn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_An_ZtSWn_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)M_Wt[nb]->nrows);
    lf_An_AtTYn_uu[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)M_Yt[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  M_WtSZn_ = (struct M_handle ***)wkspace_all0c(sizeof(struct M_handle **)*nbins);
  M_YtTAn_ = (struct M_handle ***)wkspace_all0c(sizeof(struct M_handle **)*nbins);
  for (nb=0;nb<nbins;nb++){
    M_WtSZn_[nb] = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*M_St[nb]->rpop_j);
    M_YtTAn_[nb] = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*M_Tt[nb]->rpop_j);
    for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){
      M_WtSZn_[nb][ns_j] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Wt[nb]->nrows,M_Zt[nb]->nrows);
      M_YtTAn_[nb][ns_j] = M_handle_w_make(BITJ,GLOBAL_B_MLT,M_Yt[nb]->nrows,M_At[nb]->nrows);
      /* for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){ } */}
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_Yn_WtSZn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Yn_WtSZn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Yn_YtTAn_ww = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Yn_YtTAn_ww = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_Yn_WtSZn_ww[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_St[nb]->nrows*(unsigned long long int)M_Zt[nb]->nrows);
    lf_Yn_YtTAn_ww[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)M_At[nb]->nrows);
    /* for (nb=0;nb<nbins;nb++){ } */}
  lf_Yn_WtSZn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Yn_WtSZn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  lf_Yn_YtTAn_uu = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins); length_Yn_YtTAn_uu = (int *)wkspace_all0c(sizeof(int)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_Yn_WtSZn_uu[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_St[nb]->nrows*(unsigned long long int)M_Zt[nb]->nrows);
    lf_Yn_YtTAn_uu[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)M_Tt[nb]->nrows*(unsigned long long int)M_At[nb]->nrows);
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
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Yn[nb],Y_ajdk,&(lf_Yn_ajdk[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
  	  wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Wn[nb],Y_ajdk,&(lf_Wn_ajdk[nb]));
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
	  wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_At[nb],M_Tt[nb],M_Yt[nb],A_ajdk,Y_ajdk,lf_AtTn[nb],lf_YtTn[nb],&(lf_AtTYn[nb]));
	  GLOBAL_pthread_toc();
	  GLOBAL_pthread_tic();
	  wrap_AtTYn_vv__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_A,n_spacing_A,n_spacing_A,M_Zt[nb],M_St[nb],M_Wt[nb],A_ajdk,Y_ajdk,lf_ZtSn[nb],lf_WtSn[nb],&(lf_ZtSWn[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% AnZt and YnWt: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
	  for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){
	    ns_a = M_St[nb]->m_a_[ns_j]; ns_b = M_St[nb]->m_b_[ns_j];
	    GLOBAL_pthread_tic(); 
	    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Zt[nb]->mr_j,M_Zt[nb]->mr_b,M_Wt[nb]->mr_j,M_Wt[nb]->mr_b,lf_ZtSWn[nb],L3_lf_get(lf_ZtSWn[nb],lf_ZtSWn[nb]->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_Zt[nb]->nrows),&(M_Wt[nb]->nrows),&(M_ZtSWn_[nb][ns_j]),&(GLOBAL_B_MLT),(addressable_0));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Wt[nb]->mr_j,M_Wt[nb]->mr_b,M_Zt[nb]->mr_j,M_Zt[nb]->mr_b,lf_ZtSWn[nb],L3_lf_get(lf_ZtSWn[nb],lf_ZtSWn[nb]->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_Wt[nb]->nrows),&(M_Zt[nb]->nrows),&(M_WtSZn_[nb][ns_j]),&(GLOBAL_B_MLT),(addressable_1));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_At[nb]->mr_j,M_At[nb]->mr_b,M_Yt[nb]->mr_j,M_Yt[nb]->mr_b,lf_AtTYn[nb],L3_lf_get(lf_AtTYn[nb],lf_AtTYn[nb]->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_At[nb]->nrows),&(M_Yt[nb]->nrows),&(M_AtTYn_[nb][ns_j]),&(GLOBAL_B_MLT),(addressable_0));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    wrap_xcalc(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),M_Yt[nb]->mr_j,M_Yt[nb]->mr_b,M_At[nb]->mr_j,M_At[nb]->mr_b,lf_AtTYn[nb],L3_lf_get(lf_AtTYn[nb],lf_AtTYn[nb]->lf,0,0,0,0,0,0,ns_j,ns_b,ns_a),&(M_Yt[nb]->nrows),&(M_At[nb]->nrows),&(M_YtTAn_[nb][ns_j]),&(GLOBAL_B_MLT),(addressable_1));
	    GLOBAL_pthread_toc();
	    /* for (ns_j=0;ns_j<M_St[nb]->rpop_j;ns_j++){ } */}
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% xcalc: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
	  length_An_ZtSWn_ww[nb] = wrap_An_ZtSWn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_An[nb],M_St[nb],lf_ZtSWn[nb],M_ZtSWn_[nb],M_Wt[nb],A_ajdk,(addressable_0),&(lf_An_ZtSWn_uu[nb]),&(lf_An_ZtSWn_ww[nb]));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  length_Yn_WtSZn_ww[nb] = wrap_An_ZtSWn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_Yn[nb],M_St[nb],lf_ZtSWn[nb],M_WtSZn_[nb],M_Zt[nb],Y_ajdk,(addressable_1),&(lf_Yn_WtSZn_uu[nb]),&(lf_Yn_WtSZn_ww[nb]));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  length_An_AtTYn_ww[nb] = wrap_An_ZtSWn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_An[nb],M_Tt[nb],lf_AtTYn[nb],M_AtTYn_[nb],M_Yt[nb],A_ajdk,(addressable_0),&(lf_An_AtTYn_uu[nb]),&(lf_An_AtTYn_ww[nb]));
	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
	  length_Yn_YtTAn_ww[nb] = wrap_An_ZtSWn_ww__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_spacing_B,n_spacing_B,n_spacing_B,M_Yn[nb],M_Tt[nb],lf_AtTYn[nb],M_YtTAn_[nb],M_At[nb],Y_ajdk,(addressable_1),&(lf_Yn_YtTAn_uu[nb]),&(lf_Yn_YtTAn_ww[nb]));
	  GLOBAL_pthread_toc();
  	  /* for (nb=0;nb<nbins;nb++){ } */}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_ZtSWn Yn_WtSZn: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_B=0;n_spacing_B<=2;n_spacing_B++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }}} */}}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_An_ZtSWn_ww_test]\n");}
}
