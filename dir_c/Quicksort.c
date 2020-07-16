#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

int iPartition(int *i_,int stride,int l,int r) 
{
  int pivot=0,tmpi=0;
  int i=0,j=0;
  pivot = i_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i_[stride*i] <= pivot && i <= r );
    do{ j--;} while( i_[stride*j] > pivot );
    if( i >= j ) break;
    tmpi = i_[stride*i]; i_[stride*i] = i_[stride*j]; i_[stride*j] = tmpi;
  }while(1);
  tmpi = i_[stride*l]; i_[stride*l] = i_[stride*j]; i_[stride*j] = tmpi;
  return j;
}

void iQuickSort(unsigned int nn,int *i_,int stride,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  /* test with: */
  /*
    int i_[] = { 7 , 8 , 6 , 3 , 2 , 9 , 10 , 11 , 6 , 7 , 6 };
    raprintf(i_,"int",1,11," %% i_: ");
    iQuickSort(0,i_,1,0,11-1);
    raprintf(i_,"int",1,11," %% i_: ");
  */
  int j=0; 
  if( l < r ) { 
    if (nn<GLOBAL_recursion_limit){ 
      j = iPartition(i_,stride,l,r); 
      iQuickSort(nn+1,i_,stride,l,j-1); 
      iQuickSort(nn+1,i_,stride,j+1,r);
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* if (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in iQuickSort\n",nn);
      /* recursion_limit breached */}  
    /* if( l < r ) { } */}
}

int dPartition(double *d_,int stride,int l,int r) 
{
  double pivot=0,tmpd=0;
  int i=0,j=0;
  pivot = d_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( d_[stride*i] <= pivot && i <= r );
    do{ j--;} while( d_[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = d_[stride*i]; d_[stride*i] = d_[stride*j]; d_[stride*j] = tmpd;
  }while(1);
  tmpd = d_[stride*l]; d_[stride*l] = d_[stride*j]; d_[stride*j] = tmpd;
  return j;
}

void dQuickSort(unsigned int nn,double *d_,int stride,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  /* test with: */
  /*
    double d_[] = { 7.0 , 8.0 , 6.0 , 3.0 , 2.0 , 9.0 , 10.0 , 11.0 , 6.0 , 7.0 , 6.0 };
    raprintf(d_,"double",1,11," %% d_: ");
    dQuickSort(0,d_,1,0,11-1);
    raprintf(d_,"double",1,11," %% d_: ");
  */
  int j=0; 
  if( l < r ) { 
    if (nn<GLOBAL_recursion_limit){ 
      j = dPartition(d_,stride,l,r); 
      dQuickSort(nn+1,d_,stride,l,j-1); 
      dQuickSort(nn+1,d_,stride,j+1,r);
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* if (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in dQuickSort\n",nn);
      /* recursion_limit breached */}  
    /* if( l < r ) { } */}
}

int dPartition_xij(double *d_,int stride,int *i_,int *i_b_,int *i_j_,int *lnb,int *mr_a,int *mr_b,int l,int r) 
{
  double pivot=0,tmpd=0;
  int i=0,j=0,tmpi=0;
  pivot = d_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( d_[stride*i] <= pivot && i <= r );
    do{ j--;} while( d_[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = d_[stride*i]; d_[stride*i] = d_[stride*j]; d_[stride*j] = tmpd;
    if (i_!=NULL){ tmpi = i_[i]; i_[i] = i_[j]; i_[j] = tmpi;}
    if (i_b_!=NULL){ tmpi = i_b_[i]; i_b_[i] = i_b_[j]; i_b_[j] = tmpi;}
    if (i_j_!=NULL){ tmpi = i_j_[i]; i_j_[i] = i_j_[j]; i_j_[j] = tmpi;}
    if (lnb!=NULL){ tmpi = lnb[i]; lnb[i] = lnb[j]; lnb[j] = tmpi;}
    if (mr_a!=NULL){ tmpi = mr_a[i]; mr_a[i] = mr_a[j]; mr_a[j] = tmpi;}
    if (mr_b!=NULL){ tmpi = mr_b[i]; mr_b[i] = mr_b[j]; mr_b[j] = tmpi;}
  }while(1);
  tmpd = d_[stride*l]; d_[stride*l] = d_[stride*j]; d_[stride*j] = tmpd;
  if (i_!=NULL){ tmpi = i_[l]; i_[l] = i_[j]; i_[j] = tmpi;}
  if (i_b_!=NULL){ tmpi = i_b_[l]; i_b_[l] = i_b_[j]; i_b_[j] = tmpi;}
  if (i_j_!=NULL){ tmpi = i_j_[l]; i_j_[l] = i_j_[j]; i_j_[j] = tmpi;}
  if (lnb!=NULL){ tmpi = lnb[l]; lnb[l] = lnb[j]; lnb[j] = tmpi;}
  if (mr_a!=NULL){ tmpi = mr_a[l]; mr_a[l] = mr_a[j]; mr_a[j] = tmpi;}
  if (mr_b!=NULL){ tmpi = mr_b[l]; mr_b[l] = mr_b[j]; mr_b[j] = tmpi;}
  return j;
}

unsigned int dQuickSort_xij(unsigned int nn,double *d_,int stride,int *i_,int *i_b_,int *i_j_,int *lnb,int *mr_a,int *mr_b,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  /* test with: */
  /*
    double d_[] = { 7.0 , 8.0 , 6.0 , 3.0 , 2.0 , 9.0 , 10.0 , 11.0 , 6.0 , 7.0 , 6.0 };
    int    i_[] = {   0 ,   1 ,   2 ,   3 ,   4 ,   5 ,    6 ,    7 ,   8 ,   9 ,  10 };
    int  i_b_[] = {   0 ,   1 ,   2 ,   3 ,   4 ,   5 ,    6 ,    7 ,   8 ,   9 ,  10 };
    int  i_j_[] = {   0 ,   1 ,   2 ,   3 ,   4 ,   5 ,    6 ,    7 ,   8 ,   9 ,  10 };
    int   lnb[] = {   0 ,   1 ,   0 ,   1 ,   0 ,   1 ,    0 ,    1 ,   0 ,   1 ,   0 };
    raprintf(d_,"double",1,11," %% d_: ");
    raprintf(i_,"int",1,11," %% i_: ");
    raprintf(lnb,"int",1,11," %% lnb: ");
    dQuickSort(0,d_,1,i_,i_b_,i_j_,lnb,NULL,NULL,0,11-1);
    raprintf(d_,"double",1,11," %% d_: ");
    raprintf(i_,"int",1,11," %% i_: ");
    raprintf(lnb,"int",1,11," %% lnb: ");
  */
  int j=0; unsigned int n1=nn,n2=nn;
  if( l < r ) { 
    j = dPartition_xij(d_,stride,i_,i_b_,i_j_,lnb,mr_a,mr_b,l,r); 
    if (nn<GLOBAL_recursion_limit){ 
      n1 = dQuickSort_xij(nn+1,d_,stride,i_,i_b_,i_j_,lnb,mr_a,mr_b,l,j-1); 
      n2 = dQuickSort_xij(nn+1,d_,stride,i_,i_b_,i_j_,lnb,mr_a,mr_b,j+1,r); 
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in dQuickSort_xij.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",nn);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void irandperm(unsigned int nn,int *i_)
{
  int verbose=0;
  double d_[nn];
  int nj=0;
  for (nj=0;nj<nn;nj++){ d_[nj] = rand01;}
  if (verbose){ printf(" %% [entering irandperm] nn %d\n",nn);}
  dQuickSort_xij(0,d_,1,i_,NULL,NULL,NULL,NULL,NULL,0,nn-1);
  if (verbose){ raprintf(i_,"int",1,nn," %% i_: ");}
  if (verbose){ printf(" %% [finished irandperm] nn %d\n",nn);}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int ulliPartition_index(unsigned long long int *ulli_,int stride,int *i_,int l,int r) 
{
  unsigned long long int pivot=0,tmpulli=0;
  int i=0,j=0,tmpi=0;
  pivot = ulli_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( ulli_[stride*i] <= pivot && i <= r );
    do{ j--;} while( ulli_[stride*j] > pivot );
    if( i >= j ) break;
    tmpulli = ulli_[stride*i]; ulli_[stride*i] = ulli_[stride*j]; ulli_[stride*j] = tmpulli;
    if (i_!=NULL){ tmpi = i_[i]; i_[i] = i_[j]; i_[j] = tmpi;}
  }while(1);
  tmpulli = ulli_[stride*l]; ulli_[stride*l] = ulli_[stride*j]; ulli_[stride*j] = tmpulli;
  if (i_!=NULL){ tmpi = i_[l]; i_[l] = i_[j]; i_[j] = tmpi;}
  return j;
}

unsigned int ulliQuickSort_index(unsigned int nn,unsigned long long int *ulli_,int stride,int *i_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  int j=0; unsigned int n1=nn,n2=nn;
  if( l < r ) { 
    j = ulliPartition_index(ulli_,stride,i_,l,r); 
    if (nn<GLOBAL_recursion_limit){ 
      n1 = ulliQuickSort_index(nn+1,ulli_,stride,i_,l,j-1); 
      n2 = ulliQuickSort_index(nn+1,ulli_,stride,i_,j+1,r); 
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in ulliQuickSort_index.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",nn);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void ulliQuickSort_index_driver(int n_ulli,unsigned long long int *ulli_,int stride,unsigned long long int *ulli_workspace_,int *i_)
{
  /* finds index listing i_ so that ulli_[stride*i_[.]] is sorted. */
  /* test with: ;
  int n_ulli = 5;
  int i_[n_ulli];
  unsigned long long int ulli_[2*n_ulli];
  unsigned long long int ulli_workspace_[1*n_ulli];
  ulli_[0] = 15;
  ulli_[2] = 25;
  ulli_[4] = 13;
  ulli_[6] = 18;
  ulli_[8] = 20;
  ulliQuickSort_index_driver(n_ulli,ulli_,2,ulli_workspace_,i_);
  raprintf(ulli_workspace_,"unsigned long long int",1,n_ulli," %% ulli_workspace_: ");
  raprintf(i_,"int",1,n_ulli," %% i_: ");
  exit(0);
  */
  int nulli=0;
  for (nulli=0;nulli<n_ulli;nulli++){ ulli_workspace_[nulli] = ulli_[stride*nulli]; i_[nulli]=nulli;}
  ulliQuickSort_index(0,ulli_workspace_,1,i_,0,n_ulli-1);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int dPartition_index(double *d_,int stride,int *i_,int l,int r) 
{
  double pivot=0,tmpd=0;
  int i=0,j=0,tmpi=0;
  pivot = d_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( d_[stride*i] <= pivot && i <= r );
    do{ j--;} while( d_[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = d_[stride*i]; d_[stride*i] = d_[stride*j]; d_[stride*j] = tmpd;
    if (i_!=NULL){ tmpi = i_[i]; i_[i] = i_[j]; i_[j] = tmpi;}
  }while(1);
  tmpd = d_[stride*l]; d_[stride*l] = d_[stride*j]; d_[stride*j] = tmpd;
  if (i_!=NULL){ tmpi = i_[l]; i_[l] = i_[j]; i_[j] = tmpi;}
  return j;
}

unsigned int dQuickSort_index(unsigned int nn,double *d_,int stride,int *i_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  int j=0; unsigned int n1=nn,n2=nn;
  if( l < r ) { 
    j = dPartition_index(d_,stride,i_,l,r); 
    if (nn<GLOBAL_recursion_limit){ 
      n1 = dQuickSort_index(nn+1,d_,stride,i_,l,j-1); 
      n2 = dQuickSort_index(nn+1,d_,stride,i_,j+1,r); 
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in dQuickSort_index.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",nn);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void dQuickSort_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *i_)
{
  /* finds index listing i_ so that d_[stride*i_[.]] is sorted. */
  /* test with: ;
  int n_d = 5;
  int i_[n_d];
  double d_[2*n_d];
  double d_workspace_[1*n_d];
  d_[0] = 15;
  d_[2] = 25;
  d_[4] = 13;
  d_[6] = 18;
  d_[8] = 20;
  dQuickSort_index_driver(n_d,d_,2,d_workspace_,i_);
  raprintf(d_workspace_,"double",1,n_d," %% d_workspace_: ");
  raprintf(i_,"int",1,n_d," %% i_: ");
  exit(0);
  */
  int nd=0;
  for (nd=0;nd<n_d;nd++){ d_workspace_[nd] = d_[stride*nd]; i_[nd]=nd;}
  dQuickSort_index(0,d_workspace_,1,i_,0,n_d-1);
}

void dQuickSort_index_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *i_orig_from_sort_,int *i_workspace_,int *i_sort_from_orig_)
{
  /* finds index listing i_orig_from_sort_ so that d_[stride*i_orig_from_sort_[.]] is sorted. */
  /* test with: ;
  int n_d = 5;
  int i_orig_from_sort_[n_d];
  int i_sort_from_orig_[n_d];
  double d_[2*n_d];
  double d_workspace_[1*n_d];
  int i_workspace_[1*n_d];
  d_[0] = 15;
  d_[2] = 25;
  d_[4] = 13;
  d_[6] = 18;
  d_[8] = 20;
  dQuickSort_index_index_driver(n_d,d_,2,d_workspace_,i_orig_from_sort_,i_workspace_,i_sort_from_orig_);
  raprintf(d_,"double",1,2*n_d," %% d_: ");
  raprintf(d_workspace_,"double",1,n_d," %% d_workspace_: ");
  raprintf(i_orig_from_sort_,"int",1,n_d," %% i_orig_from_sort_: ");
  raprintf(i_sort_from_orig_,"int",1,n_d," %% i_sort_from_orig_: ");
  exit(0);
  */
  int nd=0;
  for (nd=0;nd<n_d;nd++){ d_workspace_[nd] = d_[stride*nd]; i_orig_from_sort_[nd]=nd;}
  dQuickSort_index(0,d_workspace_,1,i_orig_from_sort_,0,n_d-1);
  for (nd=0;nd<n_d;nd++){ i_workspace_[nd] = i_orig_from_sort_[nd]; i_sort_from_orig_[nd]=nd;}
  lQuickSort_index(0,i_workspace_,1,i_sort_from_orig_,0,n_d-1);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int lPartition_index(int *l_,int stride,int *i_,int l,int r) 
{
  int pivot=0,tmpl=0;
  int i=0,j=0,tmpi=0;
  pivot = l_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( l_[stride*i] <= pivot && i <= r );
    do{ j--;} while( l_[stride*j] > pivot );
    if( i >= j ) break;
    tmpl = l_[stride*i]; l_[stride*i] = l_[stride*j]; l_[stride*j] = tmpl;
    if (i_!=NULL){ tmpi = i_[i]; i_[i] = i_[j]; i_[j] = tmpi;}
  }while(1);
  tmpl = l_[stride*l]; l_[stride*l] = l_[stride*j]; l_[stride*j] = tmpl;
  if (i_!=NULL){ tmpi = i_[l]; i_[l] = i_[j]; i_[j] = tmpi;}
  return j;
}

unsigned int lQuickSort_index(unsigned int nn,int *l_,int stride,int *i_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  int j=0; unsigned int n1=nn,n2=nn;
  if( l < r ) { 
    j = lPartition_index(l_,stride,i_,l,r); 
    if (nn<GLOBAL_recursion_limit){ 
      n1 = lQuickSort_index(nn+1,l_,stride,i_,l,j-1); 
      n2 = lQuickSort_index(nn+1,l_,stride,i_,j+1,r); 
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in lQuickSort_index.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",nn);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void lQuickSort_index_driver(int n_l,int *l_,int stride,int *l_workspace_,int *i_)
{
  /* finds index listing i_ so that l_[stride*i_[.]] is sorted. */
  /* test with: ;
  int n_l = 5;
  int i_[n_l];
  int l_[2*n_l];
  int l_workspace_[1*n_l];
  l_[0] = 15;
  l_[2] = 25;
  l_[4] = 13;
  l_[6] = 18;
  l_[8] = 20;
  lQuickSort_index_driver(n_l,l_,2,l_workspace_,i_);
  raprintf(l_workspace_,"int",1,n_l," %% l_workspace_: ");
  raprintf(i_,"int",1,n_l," %% i_: ");
  exit(0);
  */
  int nl=0;
  for (nl=0;nl<n_l;nl++){ l_workspace_[nl] = l_[stride*nl]; i_[nl]=nl;}
  lQuickSort_index(0,l_workspace_,1,i_,0,n_l-1);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
