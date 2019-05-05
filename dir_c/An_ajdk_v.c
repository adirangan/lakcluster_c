#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void calc_A_ajdk(double *A_p,int A_pcols,double **A_ajdk_p)
{
  /* A_ajdk[(0:A_pcols-1) + AJDK_j_k*A_pcols] stores alpha.^j.*D.^k, where alpha=p-q and D=1/(4*p*q);
     For example: A_ajdk[(0:A_pcols-1) + AJDK_1_0*A_pcols] stores alpha ;
  */
  int verbose=0,nc=0;
  double *A_ajdk=NULL;
  if (*A_ajdk_p==NULL){ (*A_ajdk_p) = (double *) wkspace_all0c(AJDK_TOT*A_pcols*sizeof(double));} 
  A_ajdk = *A_ajdk_p;
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_0_0*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),0)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),0);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_1_0*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),1)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),0);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_2_0*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),2)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),0);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_3_0*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),3)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),0);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_4_0*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),4)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),0);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_0_1*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),0)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),1);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_1_1*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),1)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),1);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_2_1*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),2)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),1);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_3_1*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),3)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),1);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_4_1*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),4)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),1);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_0_2*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),0)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),2);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_1_2*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),1)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),2);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_2_2*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),2)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),2);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_3_2*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),3)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),2);}
  for (nc=0;nc<A_pcols;nc++){ A_ajdk[nc + AJDK_4_2*A_pcols] = pow(A_p[nc]-(1-A_p[nc]),4)*pow(1.0/maximum(0.01,(4.0*A_p[nc]*(1-A_p[nc]))),2);}
  if (verbose){ raprintf(*A_ajdk_p,"double",A_pcols,AJDK_TOT," %% A_ajdk: ");} 
}

void wrap_M_setup_test_excerpt_0(char *error_vs_speed,double mrnd,int nrows,int rpop_b,int ncols,int cpop_b,unsigned char *bn,unsigned char *bt)
{
  int brows = bsize(nrows);
  int bcols = bsize(ncols);
  unsigned char *bn_tmp=NULL,*bn_end=NULL,*bt_tmp=NULL,*bt_end=NULL;
  unsigned char nx=0;
  int nr=0,nc=0;
  double dtmp=0;
  mrnd = maximum(0,minimum(1,mrnd));
  if (strcmp(error_vs_speed,"speed")==0){ 
    bn_end = bn; for (nr=0;nr<rpop_b;nr++){ bn_end += bcols;}
    nx=0; bn_tmp = bn; while (bn_tmp<bn_end){ *(bn_tmp++)=nx++;}
    bt_end = bt; for (nc=0;nc<cpop_b;nc++){ bt_end += brows;}
    nx=0; bt_tmp = bt; while (bt_tmp<bt_end){ *(bt_tmp++)=nx++;}
    /* if (strcmp(error_vs_speed,"speed")==0){ } */}
  if (strcmp(error_vs_speed,"error")==0){ 
    for (nc=0;nc<ncols;nc++){ bt_tmp = &(bt[nc*brows]); for (nr=0;nr<nrows;nr++){ bn_tmp = &(bn[nr*bcols]); 
	dtmp = rand01; if (dtmp>mrnd){ bset__on(bn_tmp,nc); bset__on(bt_tmp,nr);} else /* if (dtmp<=mrnd) */{ bset_off(bn_tmp,nc); bset_off(bt_tmp,nr);}
	/* for (nc=0;nc<ncols;nc++){ for (nr=0;nr<nrows;nr++){ }} */}}
    /* if (strcmp(error_vs_speed,"error")==0){ } */}
}

void wrap_M_setup_test_excerpt_1(char *error_vs_speed,double mrnd,int nrows,int ncols,unsigned char *mr_b,unsigned char *mc_b,unsigned char **bn_ra_p,unsigned char **bt_ra_p)
{
  int verbose=0;
  int brows = bsize(nrows);
  int bcols = bsize(ncols);
  unsigned char *bn_ra=NULL,*bt_ra=NULL;
  unsigned char *bn=NULL,*bt=NULL;
  int nr_a=0,nr_b=0,nr_x=0,nc_a=0,nc_b=0,nc_x=0;
  int rpop_b=0,cpop_b=0;
  unsigned char *wkspace_mark=NULL;
  if (verbose){ printf(" %% [entering wrap_M_setup_test_excerpt_1], %s mrnd %0.2f nrows %d ncols %d\n",error_vs_speed,mrnd,nrows,ncols);}
  rpop_b = popcount_uchar_array(mr_b,brows); cpop_b = popcount_uchar_array(mc_b,bcols);
  *bn_ra_p = wkspace_all0c(rpop_b*bcols); bn_ra = *bn_ra_p;
  *bt_ra_p = wkspace_all0c(cpop_b*brows); bt_ra = *bt_ra_p;
  if (strcmp(error_vs_speed,"error")==0){
    GLOBAL_wkspace_point = 0; wkspace_mark = wkspace_base;    
    bn = wkspace_all0c(nrows*bcols); bt = wkspace_all0c(ncols*brows);
    wrap_M_setup_test_excerpt_0("error",mrnd,nrows,rpop_b,ncols,cpop_b,bn,bt);
    nr_b=0; nr_a=0;
    while (nr_b<rpop_b && nr_a<nrows){
      if (bget__on(mr_b,nr_a)){ for (nc_x=0;nc_x<bcols;nc_x++){ bn_ra[nc_x+nr_b*bcols] = bn[nc_x+nr_a*bcols];} nr_b++; /* if (bget__on(mr_b,nr_b)){ } */}
      nr_a++; /* while (nr_b<rpop_b && nr_a<nrows){ } */}
    if (verbose>1){
      printf(" %% nrows %d; rpop_b: %d\n",nrows,rpop_b);
      bprintf(mr_b ,POPLENGTH,1     ,nrows," %% mr_b : ");
      bprintf(bn   ,POPLENGTH,nrows ,ncols," %% bn   : ");
      bprintf(bn_ra,POPLENGTH,rpop_b,ncols," %% bn_ra: ");
      /* if (verbose){ } */}
    nc_b=0; nc_a=0;
    while (nc_b<cpop_b && nc_a<ncols){
      if (bget__on(mc_b,nc_a)){ for (nr_x=0;nr_x<brows;nr_x++){ bt_ra[nr_x+nc_b*brows] = bt[nr_x+nc_a*brows];} nc_b++; /* if (bget__on(mc_b,nc_b)){ } */}
      nc_a++; /* while (nc_b<cpop_b && nc_a<ncols){ } */}
    if (verbose>1){
      printf(" %% ncols %d; cpop_b: %d\n",ncols,cpop_b);
      bprintf(mr_b ,POPLENGTH,1     ,ncols," %% mc_b : ");
      bprintf(bt   ,POPLENGTH,ncols ,nrows," %% bt   : ");
      bprintf(bt_ra,POPLENGTH,cpop_b,nrows," %% bt_ra: ");
      /* if (verbose){ } */}
    wkspace_reset(wkspace_mark); GLOBAL_wkspace_point = 1;
    /* if (strcmp(error_vs_speed,"error")==0){ } */}
  else if (strcmp(error_vs_speed,"speed")==0){
    wrap_M_setup_test_excerpt_0("speed",mrnd,nrows,rpop_b,ncols,cpop_b,bn_ra,bt_ra);
    /* else if (strcmp(error_vs_speed,"speed")==0){ } */}
  if (verbose){ printf(" %% [finished wrap_M_setup_test_excerpt_1]\n");}
}

void wrap_M_setup_test_excerpt_2(int nrows,int ncols,unsigned char *bn,unsigned char *bt,unsigned char *mr_b,unsigned char *mr_j,unsigned char *mc_b,unsigned char *mc_j,struct M_handle **Mn_p,struct M_handle **Mt_p)
{
  int nr=0,nc=0;
  struct M_handle *Mn=NULL,*Mt=NULL;
  *Mn_p = M_handle_v_make(BITJ,nrows,ncols,NULL,bn,mr_b,mc_b); Mn = (*Mn_p);
  for (nr=0;nr<bsize(nrows);nr++){ Mn->mr_j[nr] = mr_j[nr];} for (nc=0;nc<bsize(ncols);nc++){ Mn->mc_j[nc] = mc_j[nc];}
  *Mt_p = M_handle_v_make(BITJ,ncols,nrows,NULL,bt,mc_b,mr_b); Mt = (*Mt_p);
  for (nr=0;nr<bsize(ncols);nr++){ Mt->mr_j[nr] = mc_j[nr];} for (nc=0;nc<bsize(nrows);nc++){ Mt->mc_j[nc] = mr_j[nc];}
  M_mxget(Mn); M_mxget(Mt);
}

void wrap_M_setup_test(char *error_vs_speed,double mrnd,int nbins,int *A_n_rows,int A_n_cols,int *Z_n_rows,int Y_n_cols,int T_n_cols,struct M_handle ***M_An_p,struct M_handle ***M_At_p,struct M_handle ***M_Zn_p,struct M_handle ***M_Zt_p,struct M_handle ***M_Yn_p,struct M_handle ***M_Yt_p,struct M_handle ***M_Wn_p,struct M_handle ***M_Wt_p,struct M_handle ***M_Tn_p,struct M_handle ***M_Tt_p,struct M_handle ***M_Sn_p,struct M_handle ***M_St_p,double **A_p_p,double **A_ajdk_p,struct L_handle ***lf_An_ajdk_p,struct L_handle ***lf_Zn_ajdk_p,double **Y_p_p,double **Y_ajdk_p,struct L_handle ***lf_Yn_ajdk_p,struct L_handle ***lf_Wn_ajdk_p)
{
  int verbose=1;
  unsigned char **b_An=NULL,**b_At=NULL,**b_An_mr_b=NULL,**b_An_mr_j=NULL,*b_An_mc_b=NULL,*b_An_mc_j=NULL;
  unsigned char **b_Zn=NULL,**b_Zt=NULL,**b_Zn_mr_b=NULL,**b_Zn_mr_j=NULL;
  unsigned char **b_Yn=NULL,**b_Yt=NULL,*b_Yn_mc_b=NULL,*b_Yn_mc_j=NULL;
  unsigned char **b_Tn=NULL,**b_Tt=NULL,*b_Tn_mc_b=NULL,*b_Tn_mc_j=NULL;
  unsigned char **b_Wn=NULL,**b_Wt=NULL,**b_Sn=NULL,**b_St=NULL;
  unsigned char *bb=NULL,*bj=NULL,*bn=NULL,*bt=NULL;
  int bcols_A=0,bcols_Y=0,bcols_T=0;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double dtmp=0;
  int A_pcols=0,Y_pcols=0;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nb=0,nx=0,nr=0,nc=0;
  double mrnT=0;
  mrnd = maximum(0,minimum(1,mrnd)); mrnT = minimum(0.5,mrnd);
  if (verbose){ printf(" %% [entering wrap_M_setup_test], %s mrnd %0.2f nbins %d A_n_cols %d Y_n_cols %d T_n_cols %d\n",error_vs_speed,mrnd,nbins,A_n_cols,Y_n_cols,T_n_cols);}
  if (verbose){ raprintf(A_n_rows,"int",1,nbins," %% A_n_rows: ");} if (verbose){ raprintf(Z_n_rows,"int",1,nbins," %% Z_n_rows: ");}
  bcols_A = bsize(A_n_cols); bcols_Y = bsize(Y_n_cols); bcols_T = bsize(T_n_cols);
  if (verbose){ printf(" %% bcols_A %d bcols_Y %d bcols_T %d\n",bcols_A,bcols_Y,bcols_T);}
  if (verbose){ printf(" %% initializing masks: \n");}
  b_An_mr_b = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins); b_An_mr_j = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    b_An_mr_b[nb] = (unsigned char *)wkspace_all0c(bsize(A_n_rows[nb])); b_An_mr_j[nb] = (unsigned char *)wkspace_all0c(bsize(A_n_rows[nb])); bb = b_An_mr_b[nb]; bj = b_An_mr_j[nb];
    for (nx=0;nx<A_n_rows[nb];nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
    /* for (nb=0;nb<nbins;nb++){ } */}
  b_An_mc_b = (unsigned char *)wkspace_all0c(bcols_A); b_An_mc_j = (unsigned char *)wkspace_all0c(bcols_A);
  bb = b_An_mc_b; bj = b_An_mc_j; for (nx=0;nx<A_n_cols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  b_Zn_mr_b = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins); b_Zn_mr_j = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    b_Zn_mr_b[nb] = (unsigned char *)wkspace_all0c(bsize(Z_n_rows[nb])); b_Zn_mr_j[nb] = (unsigned char *)wkspace_all0c(bsize(Z_n_rows[nb])); bb = b_Zn_mr_b[nb]; bj = b_Zn_mr_j[nb];
    for (nx=0;nx<Z_n_rows[nb];nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
    /* for (nb=0;nb<nbins;nb++){ } */}
  b_Yn_mc_b = (unsigned char *)wkspace_all0c(bcols_Y); b_Yn_mc_j = (unsigned char *)wkspace_all0c(bcols_Y);
  bb = b_Yn_mc_b; bj = b_Yn_mc_j; for (nx=0;nx<Y_n_cols;nx++){ dtmp = rand01; if (dtmp>mrnd*mrnd){ bset__on(bb,nx);} if (dtmp>mrnd){ bset__on(bj,nx);}}
  b_Tn_mc_b = (unsigned char *)wkspace_all0c(bcols_T); b_Tn_mc_j = (unsigned char *)wkspace_all0c(bcols_T);
  bb = b_Tn_mc_b; bj = b_Tn_mc_j; for (nx=0;nx<T_n_cols;nx++){ dtmp = rand01; if (dtmp>mrnT*mrnT){ bset__on(bb,nx);} if (dtmp>mrnT){ bset__on(bj,nx);}} bset__on(bb,0); bset__on(bj,0); bset__on(bb,T_n_cols-1), bset__on(bj,T_n_cols-1);
  if (verbose){ printf(" %% initializing arrays: \n");}
  b_An = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins); b_At = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,A_n_rows[nb],A_n_cols,b_An_mr_b[nb],b_An_mc_b,&(b_An[nb]),&(b_At[nb]));
    /* for (nb=0;nb<nbins;nb++){ } */}
  b_Zn = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins); b_Zt = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    b_Zn[nb] = (unsigned char *)wkspace_all0c(Z_n_rows[nb]*bcols_A); b_Zt[nb] = (unsigned char *)wkspace_all0c(bsize(Z_n_rows[nb])*A_n_cols);
    wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,Z_n_rows[nb],A_n_cols,b_Zn_mr_b[nb],b_An_mc_b,&(b_Zn[nb]),&(b_Zt[nb]));
    /* for (nb=0;nb<nbins;nb++){ } */}
  b_Yn = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins); b_Yt = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    b_Yn[nb] = (unsigned char *)wkspace_all0c(A_n_rows[nb]*bcols_Y); b_Yt[nb] = (unsigned char *)wkspace_all0c(bsize(A_n_rows[nb])*Y_n_cols);
    wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,A_n_rows[nb],Y_n_cols,b_An_mr_b[nb],b_Yn_mc_b,&(b_Yn[nb]),&(b_Yt[nb]));
    /* for (nb=0;nb<nbins;nb++){ } */}
  b_Wn = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins); b_Wt = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    b_Wn[nb] = (unsigned char *)wkspace_all0c(Z_n_rows[nb]*bcols_Y); b_Wt[nb] = (unsigned char *)wkspace_all0c(bsize(Z_n_rows[nb])*Y_n_cols);
    wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,Z_n_rows[nb],Y_n_cols,b_Zn_mr_b[nb],b_Yn_mc_b,&(b_Wn[nb]),&(b_Wt[nb]));
    /* for (nb=0;nb<nbins;nb++){ } */}
  b_Tn = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins); b_Tt = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    b_Tn[nb] = (unsigned char *)wkspace_all0c(A_n_rows[nb]*bcols_T); b_Tt[nb] = (unsigned char *)wkspace_all0c(bsize(A_n_rows[nb])*T_n_cols);
    wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,A_n_rows[nb],T_n_cols,b_An_mr_b[nb],b_Tn_mc_b,&(b_Tn[nb]),&(b_Tt[nb]));
    for (nx=0;nx<popcount_uchar_array(b_An_mr_b[nb],bsize(A_n_rows[nb]));nx++){ bn = &(b_Tn[nb][nx*bcols_T]); bset__on(bn,0);}
    for (nx=0;nx<A_n_rows[nb];nx++){ bt = &(b_Tt[nb][0*bsize(A_n_rows[nb])]); bset__on(bt,nx);}
    /* for (nb=0;nb<nbins;nb++){ } */}
  b_Sn = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins); b_St = (unsigned char **)wkspace_all0c(sizeof(unsigned char *)*nbins);
  for (nb=0;nb<nbins;nb++){ 
    b_Sn[nb] = (unsigned char *)wkspace_all0c(Z_n_rows[nb]*bcols_T); b_St[nb] = (unsigned char *)wkspace_all0c(bsize(Z_n_rows[nb])*T_n_cols);
    wrap_M_setup_test_excerpt_1(error_vs_speed,0.5,Z_n_rows[nb],T_n_cols,b_Zn_mr_b[nb],b_Tn_mc_b,&(b_Sn[nb]),&(b_St[nb]));
    for (nx=0;nx<popcount_uchar_array(b_Zn_mr_b[nb],bsize(Z_n_rows[nb]));nx++){ bn = &(b_Sn[nb][nx*bcols_T]); bset__on(bn,0);}
    for (nx=0;nx<Z_n_rows[nb];nx++){ bt = &(b_St[nb][0*bsize(Z_n_rows[nb])]); bset__on(bt,nx);}
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% initializing M_handle: \n");}
  M_An = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins); M_At = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){ wrap_M_setup_test_excerpt_2(A_n_rows[nb],A_n_cols,b_An[nb],b_At[nb],b_An_mr_b[nb],b_An_mr_j[nb],b_An_mc_b,b_An_mc_j,&(M_An[nb]),&(M_At[nb])); /* for (nb=0;nb<nbins;nb++){ } */}
  M_Zn = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins); M_Zt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){ wrap_M_setup_test_excerpt_2(Z_n_rows[nb],A_n_cols,b_Zn[nb],b_Zt[nb],b_Zn_mr_b[nb],b_Zn_mr_j[nb],b_An_mc_b,b_An_mc_j,&(M_Zn[nb]),&(M_Zt[nb])); /* for (nb=0;nb<nbins;nb++){ } */}
  M_Yn = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins); M_Yt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){ wrap_M_setup_test_excerpt_2(A_n_rows[nb],Y_n_cols,b_Yn[nb],b_Yt[nb],b_An_mr_b[nb],b_An_mr_j[nb],b_Yn_mc_b,b_Yn_mc_j,&(M_Yn[nb]),&(M_Yt[nb])); /* for (nb=0;nb<nbins;nb++){ } */}
  M_Wn = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins); M_Wt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){ wrap_M_setup_test_excerpt_2(Z_n_rows[nb],Y_n_cols,b_Wn[nb],b_Wt[nb],b_Zn_mr_b[nb],b_Zn_mr_j[nb],b_Yn_mc_b,b_Yn_mc_j,&(M_Wn[nb]),&(M_Wt[nb])); /* for (nb=0;nb<nbins;nb++){ } */}
  M_Tn = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins); M_Tt = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){ wrap_M_setup_test_excerpt_2(A_n_rows[nb],T_n_cols,b_Tn[nb],b_Tt[nb],b_An_mr_b[nb],b_An_mr_j[nb],b_Tn_mc_b,b_Tn_mc_j,&(M_Tn[nb]),&(M_Tt[nb])); /* for (nb=0;nb<nbins;nb++){ } */}
  M_Sn = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins); M_St = (struct M_handle **)wkspace_all0c(sizeof(struct M_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){ wrap_M_setup_test_excerpt_2(Z_n_rows[nb],T_n_cols,b_Sn[nb],b_St[nb],b_Zn_mr_b[nb],b_Zn_mr_j[nb],b_Tn_mc_b,b_Tn_mc_j,&(M_Sn[nb]),&(M_St[nb])); /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% initializing A_ajdk and Y_ajdk: \n");}
  A_pcols = psize(A_n_cols);
  A_p = (double *)wkspace_all0c(A_pcols*sizeof(double));
  for (nc=0;nc<A_pcols;nc++){ A_p[nc] = maximum(0.1,minimum(0.9,rand01));} if (GLOBAL_TEST_sparse==0){ printf(" %% turning off D and a\n"); for (nc=0;nc<A_pcols;nc++){ A_p[nc] = 0.5;}}
  A_ajdk = (double *)wkspace_all0c(AJDK_TOT*A_pcols*sizeof(double)); 
  calc_A_ajdk(A_p,A_pcols,&(A_ajdk));
  lf_An_ajdk = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_Zn_ajdk = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_An_ajdk[nb] = L_handle_make((unsigned long long int)A_n_rows[nb]*(unsigned long long int)AJDK_TOT);
    lf_Zn_ajdk[nb] = L_handle_make((unsigned long long int)Z_n_rows[nb]*(unsigned long long int)AJDK_TOT);
    /* for (nb=0;nb<nbins;nb++){ } */}
  Y_pcols = psize(Y_n_cols);
  Y_p = (double *)wkspace_all0c(Y_pcols*sizeof(double));
  for (nc=0;nc<Y_pcols;nc++){ Y_p[nc] = maximum(0.1,minimum(0.9,rand01));} if (GLOBAL_TEST_sparse==0){ printf(" %% turning off D and a\n"); for (nc=0;nc<Y_pcols;nc++){ Y_p[nc] = 0.5;}}
  Y_ajdk = (double *)wkspace_all0c(AJDK_TOT*Y_pcols*sizeof(double));
  calc_A_ajdk(Y_p,Y_pcols,&(Y_ajdk));
  lf_Yn_ajdk = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  lf_Wn_ajdk = (struct L_handle **)wkspace_all0c(sizeof(struct L_handle *)*nbins);
  for (nb=0;nb<nbins;nb++){
    lf_Yn_ajdk[nb] = L_handle_make((unsigned long long int)A_n_rows[nb]*(unsigned long long int)AJDK_TOT);
    lf_Wn_ajdk[nb] = L_handle_make((unsigned long long int)Z_n_rows[nb]*(unsigned long long int)AJDK_TOT);
    /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% setting output: \n");}
  if (M_An_p!=NULL){ *M_An_p = M_An;} if (M_At_p!=NULL){ *M_At_p = M_At;}
  if (M_Zn_p!=NULL){ *M_Zn_p = M_Zn;} if (M_Zt_p!=NULL){ *M_Zt_p = M_Zt;}
  if (M_Yn_p!=NULL){ *M_Yn_p = M_Yn;} if (M_Yt_p!=NULL){ *M_Yt_p = M_Yt;}
  if (M_Wn_p!=NULL){ *M_Wn_p = M_Wn;} if (M_Wt_p!=NULL){ *M_Wt_p = M_Wt;}
  if (M_Tn_p!=NULL){ *M_Tn_p = M_Tn;} if (M_Tt_p!=NULL){ *M_Tt_p = M_Tt;}
  if (M_Sn_p!=NULL){ *M_Sn_p = M_Sn;} if (M_St_p!=NULL){ *M_St_p = M_St;}
  if (A_p_p!=NULL){ *A_p_p = A_p;} if (A_ajdk_p!=NULL){ *A_ajdk_p = A_ajdk;} if (lf_An_ajdk_p!=NULL){ *lf_An_ajdk_p = lf_An_ajdk;} if (lf_Zn_ajdk_p!=NULL){ *lf_Zn_ajdk_p = lf_Zn_ajdk;}
  if (Y_p_p!=NULL){ *Y_p_p = Y_p;} if (Y_ajdk_p!=NULL){ *Y_ajdk_p = Y_ajdk;} if (lf_Yn_ajdk_p!=NULL){ *lf_Yn_ajdk_p = lf_Yn_ajdk;} if (lf_Wn_ajdk_p!=NULL){ *lf_Wn_ajdk_p = lf_Wn_ajdk;}
  if (verbose){ printf(" %% [finished wrap_M_setup_test]\n");}
}

void wrap_M_setup_test_test()
{
  /* test with input file: wrap_M_setup_test.in ; 
     GLOBAL_verbose= 3;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= wrap_M_setup_test;
     GLOBAL_TEST_TYP2= error;
     GLOBAL_NBINS= 2;
     GLOBAL_TEST_mrand= 0.5;
     GLOBAL_TEST_A_n_rows= 19;
     GLOBAL_TEST_A_n_cols= 20;
     GLOBAL_TEST_Z_n_rows= 21;
     GLOBAL_TEST_Y_n_cols= 18;
     GLOBAL_TEST_T_n_cols= 3;
     END= 0;
  */
  int verbose=GLOBAL_verbose;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk=NULL,**lf_Zn_ajdk=NULL,**lf_Yn_ajdk=NULL,**lf_Wn_ajdk=NULL;
  int nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  if (verbose){ printf(" %% [entering wrap_M_setup_test_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk,&lf_Zn_ajdk,&Y_p,&Y_ajdk,&lf_Yn_ajdk,&lf_Wn_ajdk);
  M_handle_printf(M_An[0],verbose," %% M_An[0]: "); M_handle_printf(M_At[0],verbose," %% M_At[0]: ");
  M_handle_printf(M_An[1],verbose," %% M_An[1]: "); M_handle_printf(M_At[1],verbose," %% M_At[1]: ");
  M_handle_printf(M_Tn[0],verbose," %% M_Tn[0]: "); M_handle_printf(M_Tt[0],verbose," %% M_Tt[0]: ");
  M_handle_printf(M_Tn[1],verbose," %% M_Tn[1]: "); M_handle_printf(M_Tt[1],verbose," %% M_Tt[1]: ");
  M_handle_printf(M_Zn[0],verbose," %% M_Zn[0]: "); M_handle_printf(M_Zt[0],verbose," %% M_Zt[0]: ");
  M_handle_printf(M_Zn[1],verbose," %% M_Zn[1]: "); M_handle_printf(M_Zt[1],verbose," %% M_Zt[1]: ");
  M_handle_printf(M_Sn[0],verbose," %% M_Sn[0]: "); M_handle_printf(M_St[0],verbose," %% M_St[0]: ");
  M_handle_printf(M_Sn[1],verbose," %% M_Sn[1]: "); M_handle_printf(M_St[1],verbose," %% M_St[1]: ");
  raprintf(A_p,"double",1,psize(A_n_cols)," %% A_p: "); raprintf(A_ajdk,"double",psize(A_n_cols),AJDK_TOT," %% A_ajdk: ");
  raprintf(Y_p,"double",1,psize(Y_n_cols)," %% Y_p: "); raprintf(Y_ajdk,"double",psize(Y_n_cols),AJDK_TOT," %% Y_ajdk: ");
  if (verbose){ printf(" %% [finished wrap_M_setup_test_test]\n");}
}

void *get_An_ajdk_v(void *vp)
{
  /* This function serves as a module which calculates output_An_ajdk = An*diag(mc_j)*A_ajdk(j,k), assuming A_ajdk is given ;
     The terms alpha=(p-q) and D = 1/(4*p*q) are given ;
     The calculation is performed using vector operations ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_An_ajdk = (struct L_handle *)(vpra[ip++]);
  int type_flag = *(int *)(vpra[ip++]); 
  int output_spacing_r = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols)/* rup(M_An->ncols + M_An->ncols_extend,POPLENGTH)/POPLENGTH */;
  int type_flag_single = 0;
  int mx_j=0,mx_chunk=0,ajdk_x_x=0,ma_a=0,ma_b=0,ma_j=0,nc=0,tab_r=0,tab_r_stride=0,tab_O;
  double dtmp=0,dtmp_base=0;
  __m128i *wA_tag;
  __m128i *mc_tag,*mc_end;
  double *dinp=NULL;
  char tempchar[FNAMESIZE];
  if (type_flag==TYPE_p0){ printf(" %% Warning! TYPE_p0 not implemented within get_An_ajdk_v\n");}
  if      (type_flag==TYPE_pm){ type_flag_single = TYPE_p_;}
  else if (type_flag==TYPE_00){ type_flag_single = TYPE_0_;}
  if (verbose>1){ printf(" %% [entering get_An_ajdk_v] tidx %d\n",tidx);}  
  if (verbose>2){ bprintf((unsigned char *)M_An->wX,POPLENGTH,M_An->rpop_b,minimum(M_An->mc_length*BIT8,M_An->ncols)," %% M_An->wX: ");}
  unsigned int *ma_b_,*ma_a_;
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  if (verbose>1){ sprintf(tempchar," %% ma_b_[%d]: ",tidx); raprintf(ma_b_,"unsigned int",1,M_An->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% ma_a_[%d]: ",tidx); raprintf(ma_a_,"unsigned int",1,M_An->rpop_j,tempchar);}
  switch (output_spacing_r){ case SPACING_j: tab_r_stride = M_An->rpop_j; break; case SPACING_b: tab_r_stride = M_An->rpop_b; break; case SPACING_a: tab_r_stride = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  if (verbose>1){ printf(" %% calculating output_An_ajdk\n");}
  output_An_ajdk->spacing_row = output_spacing_r; output_An_ajdk->row_stride = tab_r_stride; 
  output_An_ajdk->spacing_col = SPACING_a; output_An_ajdk->col_stride = AJDK_TOT; 
  if (verbose>1){ raprintf(A_ajdk,"double",A_pcols,AJDK_TOT," %% A_ajdk: ");}
  if (GLOBAL_omp_type==GLOBAL_omp_off){
    ajdk_x_x=0;
    while (ajdk_x_x<AJDK_TOT){
      dtmp_base=0; for (nc=0;nc<A_pcols;nc++){ dtmp_base += (double)A_ajdk[nc + ajdk_x_x*A_pcols] * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8); /* for (nc=0;nc<A_pcols;nc++){ } */} GLOBAL_ops_count_one(tidx,A_pcols,A_pcols*POPLENGTH);
      ma_j=0;
      while (ma_j<M_An->rpop_j){
	ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
	tab_O = tab_r + ajdk_x_x*tab_r_stride; dtmp=dtmp_base;
	if      (type_flag_single==TYPE_p_){ output_An_ajdk->lf[tab_O] = -dtmp; }
	else if (type_flag_single==TYPE_0_){ output_An_ajdk->lf[tab_O] = 0;}
	wA_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	mc_tag = (__m128i*)&(M_An->mc_j[0]); mc_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	dinp = &(A_ajdk[0+ajdk_x_x*A_pcols]);
	if      (type_flag_single==TYPE_p_){ dtmp = popcount_lf(&wA_tag,&mc_tag,&mc_end,&dinp);}
	else if (type_flag_single==TYPE_0_){ dtmp = popcount_lf(&wA_tag,&mc_tag,&mc_end,&dinp);}
	if      (type_flag_single==TYPE_p_){ output_An_ajdk->lf[tab_O] += (2*dtmp);}
	else if (type_flag_single==TYPE_0_){ output_An_ajdk->lf[tab_O] += (1*dtmp);}
	ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
      GLOBAL_ops_count_one(tidx,M_An->rpop_j,M_An->rpop_j*M_An->mc_length*BIT8);
      ajdk_x_x++; /* while (ajdk_x_x<AJDK_TOT){ } */}
    /* if (GLOBAL_omp_type==GLOBAL_omp_off){ } */}
  else if (GLOBAL_omp_type==GLOBAL_omp__on){
    mx_chunk=1;
#pragma omp parallel shared(mx_chunk) private(mx_j,ajdk_x_x,ma_j,ma_a,ma_b,tab_r,tab_O,dtmp,dtmp_base,nc,wA_tag,mc_tag,mc_end,dinp)
    { /* begin omp parallel */
      mx_j=0; tab_r=0;
#pragma omp for schedule(dynamic,mx_chunk)
      for (mx_j=0;mx_j<AJDK_TOT;mx_j++){
	ajdk_x_x = mx_j;
	dtmp_base=0; for (nc=0;nc<A_pcols;nc++){ dtmp_base += (double)A_ajdk[nc + ajdk_x_x*A_pcols] * (double)popcount_uchar_array((unsigned char *)&(M_An->mc_j[nc*POPLENGTH/BIT8]),1*POPLENGTH/BIT8); /* for (nc=0;nc<A_pcols;nc++){ } */}
	ma_j=0;
	while (ma_j<M_An->rpop_j){
	  ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
	  switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break; default: break; /* switch (output_spacing_r){ } */}
	  tab_O = tab_r + ajdk_x_x*tab_r_stride; dtmp=dtmp_base;
	  if      (type_flag_single==TYPE_p_){ output_An_ajdk->lf[tab_O] = -dtmp; }
	  else if (type_flag_single==TYPE_0_){ output_An_ajdk->lf[tab_O] = 0;}
	  wA_tag = (__m128i*)&(M_An->wX[ma_b*M_An->mc_length]);
	  mc_tag = (__m128i*)&(M_An->mc_j[0]); mc_end = (__m128i*)&(M_An->mc_j[M_An->mc_length]);
	  dinp = &(A_ajdk[0+ajdk_x_x*A_pcols]);
	  if      (type_flag_single==TYPE_p_){ dtmp = popcount_lf(&wA_tag,&mc_tag,&mc_end,&dinp);}
	  else if (type_flag_single==TYPE_0_){ dtmp = popcount_lf(&wA_tag,&mc_tag,&mc_end,&dinp);}
#pragma omp critical
	  {/* begin omp critical */
	    if      (type_flag_single==TYPE_p_){ output_An_ajdk->lf[tab_O] += (2*dtmp);}
	    else if (type_flag_single==TYPE_0_){ output_An_ajdk->lf[tab_O] += (dtmp);}
	    /* end omp critical */}
	  ma_j++; /* while (ma_j<M_An->rpop_j){ } */}
	/* for (mx_j=0;mx_j<AJDK_TOT;mx_j++){ } */}
      /* end omp parallel */}
    GLOBAL_ops_count_one(tidx,A_pcols + M_An->rpop_j,A_pcols*POPLENGTH + M_An->rpop_j*M_An->mc_length*BIT8);
    /* else if (GLOBAL_omp_type==GLOBAL_omp__on){ } */}
  return NULL;
}

int wrap_An_ajdk_v(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,struct M_handle *M_An,double *A_ajdk,struct L_handle **output_An_ajdk_p)
{
  /* This function calls An_ajdk_v: calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 6) */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length=0,length_r=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_An_ajdk_v] *tidx %d\n",*tidx);}
  switch (output_spacing_r){ case SPACING_j: length_r = M_An->rpop_j; break; case SPACING_b: length_r = M_An->rpop_b; break; case SPACING_a: length_r = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  length = length_r*AJDK_TOT; if (*output_An_ajdk_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_ajdk_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  length = length_r*AJDK_TOT; if ((*output_An_ajdk_p)->length<length){ printf(" %% Error! length %llu<%llu in wrap_An_ajdk_v\n",(*output_An_ajdk_p)->length,length);} memset((*output_An_ajdk_p)->lf,0,length*sizeof(double)); 
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = A_ajdk; vpra[ip++] = *output_An_ajdk_p;
  switch (type_flag){ case TYPE_p0: vpra[ip++] = &addressable_type_p0; break; case TYPE_pm: vpra[ip++] = &addressable_type_pm; break; case TYPE_00: vpra[ip++] = &addressable_type_00; break; default: break; /* switch (type_flag){ } */}
  switch (output_spacing_r){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_ajdk_v,vpra)){ printf("Warning! cannot create thread %d in wrap_An_ajdk_v\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_ajdk_v(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_An_ajdk_v] *tidx %d\n",*tidx);}
  return length;
}  


void *get_An_ajdk_u(void *vp)
{
  /* This function serves as a module which calculates output_An_ajdk = An*diag(mc_j)*A_ajdk(j,k), assuming A_ajdk is given ;
     The terms alpha=(p-q) and D = 1/(4*p*q) are given ;
     calculation performed using element-wise operations ;
  */
  int verbose=0;
  int ip=0;
  void **vpra=(void **)vp;
  int tidx = *(int *)(vpra[ip++]);
  struct M_handle *M_An = (struct M_handle *)(vpra[ip++]);
  double *A_ajdk = (double *)(vpra[ip++]);
  struct L_handle *output_An_ajdk = (struct L_handle *)(vpra[ip++]);
  int type_flag = *(int *)(vpra[ip++]); 
  int output_spacing_r = *(int *)(vpra[ip++]);
  int A_pcols = psize(M_An->ncols)/* rup(M_An->ncols + M_An->ncols_extend,POPLENGTH)/POPLENGTH */;
  int type_flag_single = 0;
  int ajdk_x_x=0,ma_a=0,ma_b=0,ma_j=0,nc_a=0,tab_r=0,tab_r_stride=0,tab_O;
  double dtmp=0;
  unsigned char *A_tag=NULL;
  char tempchar[FNAMESIZE];
  if (type_flag==TYPE_p0){ printf(" %% Warning! TYPE_p0 not implemented within get_An_ajdk_u\n");}
  if      (type_flag==TYPE_pm){ type_flag_single = TYPE_p_;}
  else if (type_flag==TYPE_00){ type_flag_single = TYPE_0_;}
  if (verbose>1){ printf(" %% [entering get_An_ajdk_u] tidx %d\n",tidx);}  
  if (verbose>2){ bprintf((unsigned char *)M_An->wX,POPLENGTH,M_An->rpop_b,minimum(M_An->mc_length*BIT8,M_An->ncols)," %% M_An->wX: ");}
  unsigned int *ma_b_,*ma_a_;
  ma_b_ = M_An->m_b_; ma_a_ = M_An->m_a_;
  if (verbose>1){ sprintf(tempchar," %% ma_b_[%d]: ",tidx); raprintf(ma_b_,"unsigned int",1,M_An->rpop_j,tempchar);}
  if (verbose>1){ sprintf(tempchar," %% ma_a_[%d]: ",tidx); raprintf(ma_a_,"unsigned int",1,M_An->rpop_j,tempchar);}
  switch (output_spacing_r){ case SPACING_j: tab_r_stride = M_An->rpop_j; break; case SPACING_b: tab_r_stride = M_An->rpop_b; break; case SPACING_a: tab_r_stride = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  if (verbose>1){ printf(" %% calculating output_An_ajdk\n");}
  output_An_ajdk->spacing_row = output_spacing_r; output_An_ajdk->row_stride = tab_r_stride;
  output_An_ajdk->spacing_col = SPACING_a; output_An_ajdk->col_stride = AJDK_TOT; 
  if (verbose>1){ raprintf(A_ajdk,"double",A_pcols,AJDK_TOT," %% A_ajdk: ");}
  if (verbose>1){ raprintf(output_An_ajdk->lf,"double",tab_r_stride,AJDK_TOT," %% output_An_ajdk->lf: ");}
  memset(output_An_ajdk->lf,0,tab_r_stride*AJDK_TOT*sizeof(double));
  ajdk_x_x=0;
  while (ajdk_x_x<AJDK_TOT){
    for (ma_j=0;ma_j<M_An->rpop_j;ma_j++){
      ma_a = ma_a_[ma_j]; ma_b = ma_b_[ma_j];
      switch (output_spacing_r){ case SPACING_j: tab_r=ma_j; break; case SPACING_b: tab_r=ma_b; break; case SPACING_a: tab_r=ma_a; break;default: break; /* switch (output_spacing_r){ } */}
      tab_O = tab_r + ajdk_x_x*tab_r_stride;
      A_tag = (unsigned char *)(&(M_An->wX[ma_b*M_An->mc_length]));
      if      (type_flag==TYPE_pm){ dtmp=0; for (nc_a=0;nc_a<M_An->ncols;nc_a++){ dtmp += bget____(A_tag,nc_a)*bget__on(M_An->mc_j,nc_a)*A_ajdk[nc_a/POPLENGTH + ajdk_x_x*A_pcols];}}
      else if (type_flag==TYPE_00){ dtmp=0; for (nc_a=0;nc_a<M_An->ncols;nc_a++){ dtmp += bget__on(A_tag,nc_a)*bget__on(M_An->mc_j,nc_a)*A_ajdk[nc_a/POPLENGTH + ajdk_x_x*A_pcols];}}
      output_An_ajdk->lf[tab_O] = dtmp;
      /* for (ma_j=0;ma_j<M_An->rpop_j;ma_j++){ } */}
    ajdk_x_x++; /* while (ajdk_x_x<AJDK_TOT){ } */}
  GLOBAL_ops_count_one(tidx,AJDK_TOT*M_An->rpop_j*M_An->ncols,0);
  if (verbose>1){ raprintf(output_An_ajdk->lf,"double",tab_r_stride,AJDK_TOT," %% output_An_ajdk->lf: ");}
  if (verbose>1){ printf(" %% [finished get_An_ajdk_u] tidx %d\n",tidx);}  
  return NULL;
}

int wrap_An_ajdk_u(int *tidx,void **vpra,pthread_t *thread_in,int type_flag,int output_spacing_r,struct M_handle *M_An,double *A_ajdk,struct L_handle **output_An_ajdk_p)
{
  /* This function uses the M_handle M_An to set up An_ajdk;
     calculation performed in thread *thread_in ; thread number *tidx ;
     variable space in **vpra (should be at least size 6) */
  int verbose=0;
  /* unsigned char *wkspace_mark=NULL; */
  unsigned long long int length=0,length_r=0,ip=0;
  if (verbose){ printf(" %% [entering wrap_An_ajdk_u] *tidx %d\n",*tidx);}
  switch (output_spacing_r){ case SPACING_j: length_r = M_An->rpop_j; break; case SPACING_b: length_r = M_An->rpop_b; break; case SPACING_a: length_r = M_An->nrows; break; default: break; /* switch (output_spacing_r){ } */}
  length = length_r*AJDK_TOT; if (*output_An_ajdk_p==NULL){ if (verbose){ printf(" %% allocating output size %llu*%d\n",length,(int)sizeof(double));} *output_An_ajdk_p = L_handle_make(length);}
  if (verbose>2){ bprintf(M_An->mr_b,M_An->bitj,1,M_An->nrows," %% M_An->mr_b: ");}
  if (verbose>2){ bprintf(M_An->mr_j,M_An->bitj,1,M_An->nrows," %% M_An->mr_j: ");}
  if (verbose>2){ bprintf(M_An->mc_b,M_An->bitj,1,M_An->ncols," %% M_An->mc_b: ");}
  if (verbose>2){ bprintf(M_An->mc_j,M_An->bitj,1,M_An->ncols," %% M_An->mc_j: ");}
  length = length_r*AJDK_TOT; if ((*output_An_ajdk_p)->length<length){ printf(" %% Warning! length %llu<%llu in wrap_An_ajdk_u\n",(*output_An_ajdk_p)->length,length);} memset((*output_An_ajdk_p)->lf,0,length*sizeof(double)); 
  ip=0; vpra[ip++] = tidx; vpra[ip++] = M_An; vpra[ip++] = A_ajdk; vpra[ip++] = *output_An_ajdk_p;
  switch (type_flag){ case TYPE_p0: vpra[ip++] = &addressable_type_p0; break; case TYPE_pm: vpra[ip++] = &addressable_type_pm; break; case TYPE_00: vpra[ip++] = &addressable_type_00; break; default: break; /* switch (type_flag){ } */}
  switch (output_spacing_r){ case SPACING_j: vpra[ip++] = &addressable_spacing_j; break; case SPACING_b: vpra[ip++] = &addressable_spacing_b; break; case SPACING_a: vpra[ip++] = &addressable_spacing_a; break; default: break; /* switch (output_spacing){ } */}
  if (*tidx>0){ if (pthread_create(thread_in,NULL,&get_An_ajdk_u,vpra)){ printf("Warning! cannot create thread %d in wrap_An_ajdk_u\n",*tidx);}}
  else /* if (*tidx<=0) */{ get_An_ajdk_u(vpra);} /* must join threads later */;
  if (verbose){ printf(" %% [finished wrap_An_ajdk_u] *tidx %d\n",*tidx);}
  return length;
}  

void wrap_An_ajdk_v_test()
{
  /* test for errors with input file: An_ajdk_v_error.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= An_ajdk_v;
     GLOBAL_TEST_TYP2= error;
     GLOBAL_NBINS= 4;
     GLOBAL_TEST_mrand= 0.25;
     GLOBAL_TEST_A_n_rows= 19;
     GLOBAL_TEST_A_n_cols= 4090;
     GLOBAL_TEST_Z_n_rows= 21;
     GLOBAL_TEST_Y_n_cols= 18;
     GLOBAL_TEST_T_n_cols= 13;
     GLOBAL_TEST_niter= 5;
     END= 0;
  */
  /* test for speed with input file: An_ajdk_v_speed.in ;
     GLOBAL_verbose= 1;
     GLOBAL_thread_count= 8;
     GLOBAL_omp_type= 0;
     GLOBAL_TEST_TYPE= An_ajdk_v;
     GLOBAL_TEST_TYP2= speed;
     GLOBAL_NBINS= 7;
     GLOBAL_TEST_mrand= 0.95;
     GLOBAL_TEST_A_n_rows= 1930;
     GLOBAL_TEST_A_n_cols= 43900;
     GLOBAL_TEST_Z_n_rows= 2120;
     GLOBAL_TEST_Y_n_cols= 21000;
     GLOBAL_TEST_T_n_cols= 1;
     GLOBAL_TEST_niter= 1;
     END= 0;
  */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0);
  int iteration_max = GLOBAL_TEST_niter;
  int nbins = GLOBAL_NBINS,*A_n_rows=NULL,A_n_cols = GLOBAL_TEST_A_n_cols,*Z_n_rows=NULL,Y_n_cols = GLOBAL_TEST_Y_n_cols,T_n_cols = GLOBAL_TEST_T_n_cols;
  struct M_handle **M_An=NULL,**M_At=NULL,**M_Zn=NULL,**M_Zt=NULL,**M_Yn=NULL,**M_Yt=NULL,**M_Wn=NULL,**M_Wt=NULL,**M_Tn=NULL,**M_Tt=NULL,**M_Sn=NULL,**M_St=NULL;
  double *A_p=NULL,*A_ajdk=NULL,*Y_p=NULL,*Y_ajdk=NULL;
  struct L_handle **lf_An_ajdk_v=NULL,**lf_Zn_ajdk_v=NULL,**lf_Yn_ajdk_v=NULL,**lf_Wn_ajdk_v=NULL;
  int nl=0,nb=0,n_type=0,n_spacing_A=0,n_spacing_B=0;
  int *length_An_ajdk_v=NULL,*length_Zn_ajdk_v=NULL,*length_Yn_ajdk_v=NULL,*length_Wn_ajdk_v=NULL;
  struct L_handle **lf_An_ajdk_u=NULL,**lf_Zn_ajdk_u=NULL,**lf_Yn_ajdk_u=NULL,**lf_Wn_ajdk_u=NULL;
  int *length_An_ajdk_u=NULL,*length_Zn_ajdk_u=NULL,*length_Yn_ajdk_u=NULL,*length_Wn_ajdk_u=NULL;
  if (verbose){ printf(" %% [entering wrap_An_ajdk_v_test]\n");}
  A_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ A_n_rows[nb] = maximum(1,GLOBAL_TEST_A_n_rows + (nb%3) - 1);}
  Z_n_rows = (int *) wkspace_all0c(sizeof(int)*nbins); for (nb=0;nb<nbins;nb++){ Z_n_rows[nb] = maximum(1,GLOBAL_TEST_Z_n_rows + (nb%5) - 2);}
  wrap_M_setup_test(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,nbins,A_n_rows,A_n_cols,Z_n_rows,Y_n_cols,T_n_cols,&M_An,&M_At,&M_Zn,&M_Zt,&M_Yn,&M_Yt,&M_Wn,&M_Wt,&M_Tn,&M_Tt,&M_Sn,&M_St,&A_p,&A_ajdk,&lf_An_ajdk_v,&lf_Zn_ajdk_v,&Y_p,&Y_ajdk,&lf_Yn_ajdk_v,&lf_Wn_ajdk_v);
  length_An_ajdk_v = (int *) wkspace_all0c(sizeof(int)*nbins);
  length_Zn_ajdk_v = (int *) wkspace_all0c(sizeof(int)*nbins);
  length_Yn_ajdk_v = (int *) wkspace_all0c(sizeof(int)*nbins);
  length_Wn_ajdk_v = (int *) wkspace_all0c(sizeof(int)*nbins);
  if (error_check){
    length_An_ajdk_u = (int *) wkspace_all0c(sizeof(int)*nbins);
    length_Zn_ajdk_u = (int *) wkspace_all0c(sizeof(int)*nbins);
    length_Yn_ajdk_u = (int *) wkspace_all0c(sizeof(int)*nbins);
    length_Wn_ajdk_u = (int *) wkspace_all0c(sizeof(int)*nbins);
    lf_An_ajdk_u = (struct L_handle **) wkspace_all0c(sizeof(struct L_handle *)*nbins);
    lf_Zn_ajdk_u = (struct L_handle **) wkspace_all0c(sizeof(struct L_handle *)*nbins);
    lf_Yn_ajdk_u = (struct L_handle **) wkspace_all0c(sizeof(struct L_handle *)*nbins);
    lf_Wn_ajdk_u = (struct L_handle **) wkspace_all0c(sizeof(struct L_handle *)*nbins);
    for (nb=0;nb<nbins;nb++){
      lf_An_ajdk_u[nb] = L_handle_make((unsigned long long int)M_An[nb]->nrows*(unsigned long long int)AJDK_TOT);
      lf_Zn_ajdk_u[nb] = L_handle_make((unsigned long long int)M_Zn[nb]->nrows*(unsigned long long int)AJDK_TOT);
      lf_Yn_ajdk_u[nb] = L_handle_make((unsigned long long int)M_Yn[nb]->nrows*(unsigned long long int)AJDK_TOT);
      lf_Wn_ajdk_u[nb] = L_handle_make((unsigned long long int)M_Wn[nb]->nrows*(unsigned long long int)AJDK_TOT);
      /* for (nb=0;nb<nbins;nb++){ } */}
    /* if (error_check){ } */}
  for (n_type=1;n_type<=2;n_type++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){
  	if (verbose){ printf(" %% %s; %s\n",TYPE_name[n_type],SPACING_name[n_spacing_A]);}
  	GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
  	GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  	for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){
  	  GLOBAL_pthread_tic();
  	  length_An_ajdk_v[nb] = wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_An[nb],A_ajdk,&(lf_An_ajdk_v[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
  	  length_Zn_ajdk_v[nb] = wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Zn[nb],A_ajdk,&(lf_Zn_ajdk_v[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
  	  length_Yn_ajdk_v[nb] = wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Yn[nb],Y_ajdk,&(lf_Yn_ajdk_v[nb]));
  	  GLOBAL_pthread_toc();
  	  GLOBAL_pthread_tic();
  	  length_Wn_ajdk_v[nb] = wrap_An_ajdk_v(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Wn[nb],Y_ajdk,&(lf_Wn_ajdk_v[nb]));
  	  GLOBAL_pthread_toc();
  	  /* for (nl=0;nl<iteration_max;nl++){ for (nb=0;nb<nbins;nb++){ }} */}}
  	GLOBAL_pthread_tuc();
	GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_ajdk_v: ");
	GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	if (error_check){
	  GLOBAL_tic(0); GLOBAL_ops_reset_all(); GLOBAL_ops_f_sum=0; GLOBAL_ops_b_sum=0;
	  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
	  for (nb=0;nb<nbins;nb++){
	    GLOBAL_pthread_tic();
	    length_An_ajdk_u[nb] = wrap_An_ajdk_u(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_An[nb],A_ajdk,&(lf_An_ajdk_u[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_Zn_ajdk_u[nb] = wrap_An_ajdk_u(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Zn[nb],A_ajdk,&(lf_Zn_ajdk_u[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_Yn_ajdk_u[nb] = wrap_An_ajdk_u(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Yn[nb],Y_ajdk,&(lf_Yn_ajdk_u[nb]));
	    GLOBAL_pthread_toc();
	    GLOBAL_pthread_tic();
	    length_Wn_ajdk_u[nb] = wrap_An_ajdk_u(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),n_type,n_spacing_A,M_Wn[nb],Y_ajdk,&(lf_Wn_ajdk_u[nb]));
	    GLOBAL_pthread_toc();
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  GLOBAL_pthread_tuc();
	  GLOBAL_ops_addup_all(); GLOBAL_ops_printf_all(verbose && !error_check," %% An_ajdk_u: ");
	  GLOBAL_ops_toc(-1,0,verbose && !error_check," %% total time: ");
	  for (nb=0;nb<nbins;nb++){
	    printf(" %% lf_An_ajdk_v[%d] error %0.16f\n",nb,dra_diff(lf_An_ajdk_v[nb]->lf,lf_An_ajdk_u[nb]->lf,length_An_ajdk_v[nb],1));
	    printf(" %% lf_Zn_ajdk_v[%d] error %0.16f\n",nb,dra_diff(lf_Zn_ajdk_v[nb]->lf,lf_Zn_ajdk_u[nb]->lf,length_Zn_ajdk_v[nb],1));
	    printf(" %% lf_Yn_ajdk_v[%d] error %0.16f\n",nb,dra_diff(lf_Yn_ajdk_v[nb]->lf,lf_Yn_ajdk_u[nb]->lf,length_Yn_ajdk_v[nb],1));
	    printf(" %% lf_Wn_ajdk_v[%d] error %0.16f\n",nb,dra_diff(lf_Wn_ajdk_v[nb]->lf,lf_Wn_ajdk_u[nb]->lf,length_Wn_ajdk_v[nb],1));
	    /* for (nb=0;nb<nbins;nb++){ } */}
	  /* if (error_check){ } */}
  	/* for (n_type=1;n_type<=2;n_type++){ for (n_spacing_A=0;n_spacing_A<=2;n_spacing_A++){ }} */}}
  wkspace_printf();
  if (verbose){ printf(" %% [finished wrap_An_ajdk_v_test]\n");}
}  
