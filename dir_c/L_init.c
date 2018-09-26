#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

struct L_handle *L_handle_make(unsigned long long int length)
{
  struct L_handle *L=NULL;
  L = (struct L_handle *)wkspace_all0c(sizeof(struct L_handle)*1);
  L->length = length;
  L->lf = (double *)wkspace_all0c(sizeof(double)*length);
  L->spacing_row = SPACING_b; L->row_stride=0;
  L->spacing_col = SPACING_b; L->col_stride=0;
  L->spacing_lyr = SPACING_b; L->lyr_stride=0;
  return L;
}

void L_handle_copy(struct L_handle *L,struct L_handle *L_in)
{
  if (L!=NULL && L_in!=NULL){
    L->length = L_in->length;
    memcpy(L->lf,L_in->lf,sizeof(double)*L->length);
    L->spacing_row=L_in->spacing_row;
    L->row_stride=L_in->row_stride;
    L->spacing_col=L_in->spacing_col;
    L->col_stride=L_in->col_stride;
    L->spacing_lyr=L_in->spacing_lyr;
    L->lyr_stride=L_in->lyr_stride;
    /* if (L!=NULL && L_in!=NULL){ } */}
}

void L_zero(struct L_handle *L){ unsigned long long int nl=0; for (nl=0;nl<L->length;nl++){ L->lf[nl]=0;}}

void L2_transpose(struct L_handle *L1,struct L_handle *L2)
{
  unsigned long long int nr=0,nc=0;
  L1->spacing_row = L2->spacing_col; L1->spacing_col = L2->spacing_row; L1->spacing_lyr = L2->spacing_lyr;
  L1->row_stride = L2->col_stride; L1->col_stride = L2->row_stride; L1->lyr_stride = L2->lyr_stride;
  for (nc=0;nc<L2->col_stride;nc++){ for (nr=0;nr<L2->row_stride;nr++){ L1->lf[nc+nr*L1->row_stride] = L2->lf[nr+nc*L2->row_stride];}}
}

void L2_duplicate(struct L_handle *L1,struct L_handle *L2)
{
  unsigned long long int nr=0,nc=0;
  L1->spacing_row = L2->spacing_row; L1->spacing_col = L2->spacing_col; L1->spacing_lyr = L2->spacing_lyr;
  L1->row_stride = L2->row_stride; L1->col_stride = L2->col_stride; L1->lyr_stride = L2->lyr_stride;
  for (nc=0;nc<L2->col_stride;nc++){ for (nr=0;nr<L2->row_stride;nr++){ L1->lf[nr+nc*L1->row_stride] = L2->lf[nr+nc*L2->row_stride];}}
}

void L3_duplicate(struct L_handle *L1,struct L_handle *L2)
{
  unsigned long long int nr=0,nc=0,nl=0;
  L1->spacing_row = L2->spacing_row; L1->spacing_col = L2->spacing_col; L1->spacing_lyr = L2->spacing_lyr;
  L1->row_stride = L2->row_stride; L1->col_stride = L2->col_stride; L1->lyr_stride = L2->lyr_stride;
  for (nl=0;nl<L2->lyr_stride;nl++){ for (nc=0;nc<L2->col_stride;nc++){ for (nr=0;nr<L2->row_stride;nr++){ L1->lf[nr+nc*L1->row_stride+nl*L1->row_stride*L1->col_stride] = L2->lf[nr+nc*L2->row_stride+nl*L2->row_stride*L2->col_stride];}}}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double *L3_get(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a)
{
  unsigned long long int tab_r=0,tab_c=0,tab_l=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L3_get\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L3_get\n");}
  if (nl_j>nl_b || nl_j>nl_a || nl_b>nl_a){ printf(" %% Warning! improper indices in L3_get\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  switch (L->spacing_lyr){ case SPACING_j: tab_l = nl_j; break; case SPACING_b: tab_l = nl_b; break; case SPACING_a: tab_l = nl_a; break; default: break; /* switch (L->spacing_lyr){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L3_get\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L3_get\n",tab_c,L->col_stride);}
  if (tab_l<0 || tab_l>=L->lyr_stride){ printf(" %% Warning! tab_l %llu/%llu in L3_get\n",tab_l,L->lyr_stride);}
  return &(L->lf[tab_r + tab_c*L->row_stride + tab_l*L->row_stride*L->col_stride]);
}

double *L2_get(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a)
{
  unsigned long long int tab_r=0,tab_c=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L2_get: nr: %llu,%llu,%llu \n",nr_j,nr_b,nr_a);}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L2_get: nc: %llu,%llu,%llu \n",nc_j,nc_b,nc_a);}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L2_get\n",tab_r,L->row_stride); exit(1);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L2_get\n",tab_c,L->col_stride); exit(1);}
  return &(L->lf[tab_r + tab_c*L->row_stride]);
}

double *L1_get(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a)
{
  unsigned long long int tab_r=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L1_get\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L1_get\n",tab_r,L->row_stride);}
  return &(L->lf[tab_r]);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void L3_set(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a,double val)
{
  unsigned long long int tab_r=0,tab_c=0,tab_l=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L3_set\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L3_set\n");}
  if (nl_j>nl_b || nl_j>nl_a || nl_b>nl_a){ printf(" %% Warning! improper indices in L3_set\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  switch (L->spacing_lyr){ case SPACING_j: tab_l = nl_j; break; case SPACING_b: tab_l = nl_b; break; case SPACING_a: tab_l = nl_a; break; default: break; /* switch (L->spacing_lyr){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L3_set\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L3_set\n",tab_c,L->col_stride);}
  if (tab_l<0 || tab_l>=L->lyr_stride){ printf(" %% Warning! tab_l %llu/%llu in L3_set\n",tab_l,L->lyr_stride);}
  L->lf[tab_r + tab_c*L->row_stride + tab_l*L->row_stride*L->col_stride] = val;
}

void L2_set(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,double val)
{
  unsigned long long int tab_r=0,tab_c=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L2_set\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L2_set\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L2_set\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L2_set\n",tab_c,L->col_stride);}
  L->lf[tab_r + tab_c*L->row_stride] = val;
}

void L1_set(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,double val)
{
  unsigned long long int tab_r=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L1_set\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L1_set\n",tab_r,L->row_stride);}
  L->lf[tab_r] = val;
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void L3_plusequals(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a,double val)
{
  unsigned long long int tab_r=0,tab_c=0,tab_l=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L3_plusequals\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L3_plusequals\n");}
  if (nl_j>nl_b || nl_j>nl_a || nl_b>nl_a){ printf(" %% Warning! improper indices in L3_plusequals\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  switch (L->spacing_lyr){ case SPACING_j: tab_l = nl_j; break; case SPACING_b: tab_l = nl_b; break; case SPACING_a: tab_l = nl_a; break; default: break; /* switch (L->spacing_lyr){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L3_plusequals\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L3_plusequals\n",tab_c,L->col_stride);}
  if (tab_l<0 || tab_l>=L->lyr_stride){ printf(" %% Warning! tab_l %llu/%llu in L3_plusequals\n",tab_l,L->lyr_stride);}
  L->lf[tab_r + tab_c*L->row_stride + tab_l*L->row_stride*L->col_stride] += val;
}

void L2_plusequals(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,double val)
{
  unsigned long long int tab_r=0,tab_c=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L2_plusequals\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L2_plusequals\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L2_plusequals\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L2_plusequals\n",tab_c,L->col_stride);}
  L->lf[tab_r + tab_c*L->row_stride] += val;
}

void L1_plusequals(struct L_handle *L,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,double val)
{
  unsigned long long int tab_r=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L1_plusequals\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L1_plusequals\n",tab_r,L->row_stride);}
  L->lf[tab_r] += val;
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double *L3_lf_get(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a)
{
  unsigned long long int tab_r=0,tab_c=0,tab_l=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L3_lf_get\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L3_lf_get\n");}
  if (nl_j>nl_b || nl_j>nl_a || nl_b>nl_a){ printf(" %% Warning! improper indices in L3_lf_get\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  switch (L->spacing_lyr){ case SPACING_j: tab_l = nl_j; break; case SPACING_b: tab_l = nl_b; break; case SPACING_a: tab_l = nl_a; break; default: break; /* switch (L->spacing_lyr){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L3_lf_get\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L3_lf_get\n",tab_c,L->col_stride);}
  if (tab_l<0 || tab_l>=L->lyr_stride){ printf(" %% Warning! tab_l %llu/%llu in L3_lf_get\n",tab_l,L->lyr_stride);}
  return &(lf[tab_r + tab_c*L->row_stride + tab_l*L->row_stride*L->col_stride]);
}

double *L2_lf_get(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a)
{
  unsigned long long int tab_r=0,tab_c=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L2_lf_get: nr: %llu,%llu,%llu \n",nr_j,nr_b,nr_a);}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L2_lf_get: nc: %llu,%llu,%llu \n",nc_j,nc_b,nc_a);}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L2_lf_get\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L2_lf_get\n",tab_c,L->col_stride);}
  return &(lf[tab_r + tab_c*L->row_stride]);
}

double *L1_lf_get(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a)
{
  unsigned long long int tab_r=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L1_lf_get\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L1_lf_get\n",tab_r,L->row_stride);}
  return &(lf[tab_r]);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void L3_lf_set(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a,double val)
{
  unsigned long long int tab_r=0,tab_c=0,tab_l=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L3_lf_set\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L3_lf_set\n");}
  if (nl_j>nl_b || nl_j>nl_a || nl_b>nl_a){ printf(" %% Warning! improper indices in L3_lf_set\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  switch (L->spacing_lyr){ case SPACING_j: tab_l = nl_j; break; case SPACING_b: tab_l = nl_b; break; case SPACING_a: tab_l = nl_a; break; default: break; /* switch (L->spacing_lyr){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L3_lf_set\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L3_lf_set\n",tab_c,L->col_stride);}
  if (tab_l<0 || tab_l>=L->lyr_stride){ printf(" %% Warning! tab_l %llu/%llu in L3_lf_set\n",tab_l,L->lyr_stride);}
  lf[tab_r + tab_c*L->row_stride + tab_l*L->row_stride*L->col_stride] = val;
}

void L2_lf_set(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,double val)
{
  unsigned long long int tab_r=0,tab_c=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L2_lf_set\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L2_lf_set\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L2_lf_set\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L2_lf_set\n",tab_c,L->col_stride);}
  lf[tab_r + tab_c*L->row_stride] = val;
}

void L1_lf_set(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,double val)
{
  unsigned long long int tab_r=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L1_lf_set\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L1_lf_set\n",tab_r,L->row_stride);}
  lf[tab_r] = val;
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void L3_lf_plusequals(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,unsigned long long int nl_j,unsigned long long int nl_b,unsigned long long int nl_a,double val)
{
  unsigned long long int tab_r=0,tab_c=0,tab_l=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L3_lf_plusequals\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L3_lf_plusequals\n");}
  if (nl_j>nl_b || nl_j>nl_a || nl_b>nl_a){ printf(" %% Warning! improper indices in L3_lf_plusequals\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  switch (L->spacing_lyr){ case SPACING_j: tab_l = nl_j; break; case SPACING_b: tab_l = nl_b; break; case SPACING_a: tab_l = nl_a; break; default: break; /* switch (L->spacing_lyr){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L3_lf_plusequals\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L3_lf_plusequals\n",tab_c,L->col_stride);}
  if (tab_l<0 || tab_l>=L->lyr_stride){ printf(" %% Warning! tab_l %llu/%llu in L3_lf_plusequals\n",tab_l,L->lyr_stride);}
  lf[tab_r + tab_c*L->row_stride + tab_l*L->row_stride*L->col_stride] += val;
}

void L2_lf_plusequals(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,unsigned long long int nc_j,unsigned long long int nc_b,unsigned long long int nc_a,double val)
{
  unsigned long long int tab_r=0,tab_c=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L2_lf_plusequals\n");}
  if (nc_j>nc_b || nc_j>nc_a || nc_b>nc_a){ printf(" %% Warning! improper indices in L2_lf_plusequals\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break; /* switch (L->spacing_col){ } */}  
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L2_lf_plusequals\n",tab_r,L->row_stride);}
  if (tab_c<0 || tab_c>=L->col_stride){ printf(" %% Warning! tab_c %llu/%llu in L2_lf_plusequals\n",tab_c,L->col_stride);}
  lf[tab_r + tab_c*L->row_stride] += val;
}

void L1_lf_plusequals(struct L_handle *L,double *lf,unsigned long long int nr_j,unsigned long long int nr_b,unsigned long long int nr_a,double val)
{
  unsigned long long int tab_r=0;
  if (nr_j>nr_b || nr_j>nr_a || nr_b>nr_a){ printf(" %% Warning! improper indices in L1_lf_plusequals\n");}
  switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break; /* switch (L->spacing_row){ } */}
  if (tab_r<0 || tab_r>=L->row_stride){ printf(" %% Warning! tab_r %llu/%llu in L1_lf_plusequals\n",tab_r,L->row_stride);}
  lf[tab_r] += val;
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void L2_clean(struct L_handle *L,unsigned long long int nrows,unsigned char *mr_b,unsigned char *mr_j,unsigned long long int ncols,unsigned char *mc_b,unsigned char *mc_j)
{
  unsigned long long int nr_j=0,nr_b=0,nr_a=0; unsigned long long int nc_j=0,nc_b=0,nc_a=0; unsigned long long int tab_r=0,tab_c=0,tab_x=0;
  nc_j=0;nc_b=0;nc_a=0;
  while (nc_a<ncols){
    if (bget__on(mc_b,nc_a)){
      if (bget__on(mc_j,nc_a)){
	switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break;}
	nr_j=0;nr_b=0;nr_a=0;
	while (nr_a<nrows){
	  if (bget__on(mr_b,nr_a)){
	    if (bget__on(mr_j,nr_a)){
	      switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break;}
	      tab_x = tab_r + tab_c*L->row_stride;
	      /* do nothing */
	      nr_j++; /* if (bget__on(mr_j,nr_a)){ } */}
	    else /* if not in mr_j */{
	      switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break;}
	      tab_x = tab_r + tab_c*L->row_stride;
	      if (L->spacing_row==SPACING_b || L->spacing_row==SPACING_a){ L->lf[tab_x]=0;}
	      /* if not in mr_j */}
	    nr_b++; /* if (bget__on(mr_b,nr_a)){ } */}
	  else /* if not in mr_b */{
	    switch (L->spacing_row){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break;}
	    tab_x = tab_r + tab_c*L->row_stride;
	    if (L->spacing_row==SPACING_a){ L->lf[tab_x]=0;}
	    /* if not in mr_b */}
	  nr_a++; /* while (nr_a<nrows){ } */}
	nc_j++; /* if (bget__on(mc_j,nc_a)){ } */}
      else /* if not in mc_j */{
	switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break;}
	if (L->spacing_col==SPACING_b || L->spacing_col==SPACING_a){ for (nr_a=0;nr_a<L->row_stride;nr_a++){ tab_x = nr_a + tab_c*L->row_stride; L->lf[tab_x]=0;}}
	/* if not in bmc_j */}
      nc_b++; /* if (bget__on(bmc_b,nc_a)){ } */}
    else /* if not in bmc_b */{
      switch (L->spacing_col){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break;}
      if (L->spacing_col==SPACING_a){ for (nr_a=0;nr_a<L->row_stride;nr_a++){ tab_x = nr_a + tab_c*L->row_stride; L->lf[tab_x]=0;}}
      /* if not in bmc_b */}
    nc_a++;/* while (nc_a<ncols){ } */}
}

void LL2_plustimesequals(struct L_handle *LA,struct L_handle *LB,double dtmp,unsigned long long int nrows,unsigned char *mr_b,unsigned char *mr_j,unsigned long long int ncols,unsigned char *mc_b,unsigned char *mc_j)
{
  unsigned long long int nr_a=0,nr_b=0,nr_j=0,nc_a=0,nc_b=0,nc_j=0;
  nc_j=0;nc_b=0;nc_a=0;
  while (nc_a<ncols){
    if (bget__on(mc_b,nc_a)){
      if (bget__on(mc_j,nc_a)){
	nr_j=0;nr_b=0;nr_a=0;
	while (nr_a<nrows){
	  if (bget__on(mr_b,nr_a)){
	    if (bget__on(mr_j,nr_a)){
	      L2_plusequals(LA , nr_j,nr_b,nr_a , nc_j,nc_b,nc_a , dtmp*(*L2_get(LB , nr_j,nr_b,nr_a , nc_j,nc_b,nc_a)));
	      nr_j++; /* if (bget__on(mr_j,nr_a)){ } */}
	    else /* if not in mr_j */{
	      /* do nothing */
	      /* if not in mr_j */}
	    nr_b++; /* if (bget__on(mr_b,nr_a)){ } */}
	  else /* if not in mr_b */{
	    /* do nothing */
	    /* if not in mr_b */}
	  nr_a++; /* while (nr_a<nrows){ } */}
	nc_j++; /* if (bget__on(mc_j,nc_a)){ } */}
      else /* if not in mc_j */{
	/* do nothing */
	/* if not in bmc_j */}
      nc_b++; /* if (bget__on(bmc_b,nc_a)){ } */}
    else /* if not in bmc_b */{
      /* do nothing */
      /* if not in bmc_b */}
    nc_a++;/* while (nc_a<ncols){ } */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

void lfprintf(struct L_handle *L,char *prefix)
{
  unsigned long long int nl=0;
  char tmpchar[FNAMESIZE];
  if (L->lyr_stride>1){
    for (nl=0;nl<L->lyr_stride;nl++){
      sprintf(tmpchar,"%s [%llu]: ",prefix,nl);
      raprintf(&(L->lf[nl*L->row_stride*L->col_stride]),"double",L->row_stride,L->col_stride,tmpchar);
	/* for (nl=0;nl<L->lyr_stride;nl++){ } */}
    /* if (L->lyr_stride>1){ } */}
  else if (L->lyr_stride<=1){
    if (L->row_stride>0 && L->col_stride>1){
      raprintf(L->lf,"double",L->row_stride,L->col_stride,prefix);
      /* if (L->row_stride>0 && L->col_stride>1){ } */}
    else if (L->row_stride>0 && L->col_stride<=1){
      raprintf(L->lf,"double",1,L->row_stride,prefix);
      /* else if (L->row_stride>0 && L->col_stride==1){ } */}
    /* else if (L->lyr_stride<=1){ } */}
}
