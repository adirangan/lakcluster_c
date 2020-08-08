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
    do{ i++;} while( i<=r && i_[stride*i] <= pivot);
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
    do{ i++;} while( i<=r && d_[stride*i] <= pivot);
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
    do{ i++;} while( i<=r && d_[stride*i] <= pivot );
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

int ulliPartition_index(unsigned long long int *ulli_,int stride,int *index_,int l,int r) 
{
  unsigned long long int pivot=0,tmpulli=0;
  int i=0,j=0,tmpi=0;
  pivot = ulli_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i<=r && ulli_[stride*i] <= pivot );
    do{ j--;} while( ulli_[stride*j] > pivot );
    if( i >= j ) break;
    tmpulli = ulli_[stride*i]; ulli_[stride*i] = ulli_[stride*j]; ulli_[stride*j] = tmpulli;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  tmpulli = ulli_[stride*l]; ulli_[stride*l] = ulli_[stride*j]; ulli_[stride*j] = tmpulli;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  return j;
}

unsigned int ulliQuickSort_index(unsigned int nn,unsigned long long int *ulli_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  int j=0; unsigned int n1=nn,n2=nn;
  if( l < r ) { 
    j = ulliPartition_index(ulli_,stride,index_,l,r); 
    if (nn<GLOBAL_recursion_limit){ 
      n1 = ulliQuickSort_index(nn+1,ulli_,stride,index_,l,j-1); 
      n2 = ulliQuickSort_index(nn+1,ulli_,stride,index_,j+1,r); 
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in ulliQuickSort_index.\n",nn);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void ulliQuickSort_index_driver(int n_ulli,unsigned long long int *ulli_,int stride,unsigned long long int *ulli_workspace_,int *index_)
{
  /* finds index listing index_ so that ulli_[stride*index_[.]] is sorted. */
  /* test with: ;
  int n_ulli = 5;
  int index_[n_ulli];
  unsigned long long int ulli_[2*n_ulli];
  unsigned long long int ulli_workspace_[1*n_ulli];
  ulli_[0] = 15;
  ulli_[2] = 25;
  ulli_[4] = 13;
  ulli_[6] = 18;
  ulli_[8] = 20;
  ulliQuickSort_index_driver(n_ulli,ulli_,2,ulli_workspace_,index_);
  raprintf(ulli_workspace_,"unsigned long long int",1,n_ulli," %% ulli_workspace_: ");
  raprintf(index_,"int",1,n_ulli," %% index_: ");
  */
  int nulli=0;
  for (nulli=0;nulli<n_ulli;nulli++){ ulli_workspace_[nulli] = ulli_[stride*nulli]; index_[nulli]=nulli;}
  ulliQuickSort_index(0,ulli_workspace_,1,index_,0,n_ulli-1);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int dPartition_index(double *d_,int stride,int *index_,int l,int r) 
{
  double pivot=0,tmpd=0;
  int i=0,j=0,tmpi=0;
  pivot = d_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i<=r && d_[stride*i] <= pivot );
    do{ j--;} while( d_[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = d_[stride*i]; d_[stride*i] = d_[stride*j]; d_[stride*j] = tmpd;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  tmpd = d_[stride*l]; d_[stride*l] = d_[stride*j]; d_[stride*j] = tmpd;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  return j;
}

unsigned int dQuickSort_index(unsigned int nn,double *d_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  int j=0; unsigned int n1=nn,n2=nn;
  if( l < r ) { 
    j = dPartition_index(d_,stride,index_,l,r); 
    if (nn<GLOBAL_recursion_limit){ 
      n1 = dQuickSort_index(nn+1,d_,stride,index_,l,j-1); 
      n2 = dQuickSort_index(nn+1,d_,stride,index_,j+1,r); 
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in dQuickSort_index.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",nn);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void dQuickSort_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_)
{
  /* finds index listing index_ so that d_[stride*index_[.]] is sorted. */
  /* test with: ;
  int n_d = 5;
  int index_[n_d];
  double d_[2*n_d];
  double d_workspace_[1*n_d];
  d_[0] = 15;
  d_[2] = 25;
  d_[4] = 13;
  d_[6] = 18;
  d_[8] = 20;
  dQuickSort_index_driver(n_d,d_,2,d_workspace_,index_);
  raprintf(d_workspace_,"double",1,n_d," %% d_workspace_: ");
  raprintf(index_,"int",1,n_d," %% index_: ");
  */
  int nd=0;
  for (nd=0;nd<n_d;nd++){ d_workspace_[nd] = d_[stride*nd]; index_[nd]=nd;}
  dQuickSort_index(0,d_workspace_,1,index_,0,n_d-1);
}

void dQuickSort_index_index_driver(int n_d,double *d_,int stride,double *d_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_)
{
  /* finds index listing index_orig_from_sort_ so that d_[stride*index_orig_from_sort_[.]] is sorted. */
  /* test with: ;
  int n_d = 5;
  int index_orig_from_sort_[n_d];
  int index_sort_from_orig_[n_d];
  double d_[2*n_d];
  double d_workspace_[1*n_d];
  int index_workspace_[1*n_d];
  d_[0] = 15;
  d_[2] = 25;
  d_[4] = 13;
  d_[6] = 18;
  d_[8] = 20;
  dQuickSort_index_index_driver(n_d,d_,2,d_workspace_,index_orig_from_sort_,index_sort_from_orig_,index_workspace_);
  raprintf(d_,"double",1,2*n_d," %% d_: ");
  raprintf(d_workspace_,"double",1,n_d," %% d_workspace_: ");
  raprintf(index_orig_from_sort_,"int",1,n_d," %% index_orig_from_sort_: ");
  raprintf(index_sort_from_orig_,"int",1,n_d," %% index_sort_from_orig_: ");
  */
  int nd=0;
  for (nd=0;nd<n_d;nd++){ d_workspace_[nd] = d_[stride*nd]; if (index_orig_from_sort_!=NULL){ index_orig_from_sort_[nd]=nd;}}
  dQuickSort_index(0,d_workspace_,1,index_orig_from_sort_,0,n_d-1);
  if (index_sort_from_orig_!=NULL){
    for (nd=0;nd<n_d;nd++){ index_workspace_[nd] = index_orig_from_sort_[nd]; index_sort_from_orig_[nd]=nd;}
    lQuickSort_index(0,index_workspace_,1,index_sort_from_orig_,0,n_d-1);
    /* if (index_sort_from_orig_!=NULL){ } */}
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int lPartition_index(int *l_,int stride,int *index_,int l,int r) 
{
  int pivot=0,tmpl=0;
  int i=0,j=0,tmpi=0;
  pivot = l_[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( i<=r && l_[stride*i] <= pivot );
    do{ j--;} while( l_[stride*j] > pivot );
    if( i >= j ) break;
    tmpl = l_[stride*i]; l_[stride*i] = l_[stride*j]; l_[stride*j] = tmpl;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  tmpl = l_[stride*l]; l_[stride*l] = l_[stride*j]; l_[stride*j] = tmpl;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  return j;
}

unsigned int lQuickSort_index(unsigned int nn,int *l_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  int j=0; unsigned int n1=nn,n2=nn;
  if( l < r ) { 
    j = lPartition_index(l_,stride,index_,l,r); 
    if (nn<GLOBAL_recursion_limit){ 
      n1 = lQuickSort_index(nn+1,l_,stride,index_,l,j-1); 
      n2 = lQuickSort_index(nn+1,l_,stride,index_,j+1,r); 
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in lQuickSort_index.\n",nn);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}

void lQuickSort_index_driver(int n_l,int *l_,int stride,int *l_workspace_,int *index_)
{
  /* finds index listing index_ so that l_[stride*index_[.]] is sorted. */
  /* test with: ;
  int n_l = 5;
  int index_[n_l];
  int l_[2*n_l];
  int l_workspace_[1*n_l];
  l_[0] = 15;
  l_[2] = 25;
  l_[4] = 13;
  l_[6] = 18;
  l_[8] = 20;
  lQuickSort_index_driver(n_l,l_,2,l_workspace_,index_);
  raprintf(l_workspace_,"int",1,n_l," %% l_workspace_: ");
  raprintf(index_,"int",1,n_l," %% index_: ");
  */
  int nl=0;
  for (nl=0;nl<n_l;nl++){ l_workspace_[nl] = l_[stride*nl]; index_[nl]=nl;}
  lQuickSort_index(0,l_workspace_,1,index_,0,n_l-1);
}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int charpPartition_index(char **charp_,int stride,int *index_,int l,int r) 
{
  int verbose=0;
  char *pivot=NULL,*tmpcharp=NULL;
  int i=0,j=0,tmpi=0;
  if (verbose){ printf(" %% [entering charpPartition_index] stride %d l %d r %d\n",stride,l,r);}
  pivot = charp_[stride*l];
  if (verbose){ printf(" %% pivot %s\n",pivot);}
  i = l; j = r+1;
  if (verbose){ printf(" %% entering loop\n");}
  do{
    if (verbose){ printf(" %% increase i\n");}
    do{ i++; if (verbose){ printf(" %% i %d: stride*i %d --> %s vs %s\n",i,stride*i,charp_[stride*i],pivot);} } while( i<=r && strcmp(charp_[stride*i],pivot)<=0);
    if (verbose){ printf(" %% decrease j\n");}
    do{ j--; if (verbose){ printf(" %% j %d: stride*j %d --> %s vs %s\n",j,stride*j,charp_[stride*j],pivot);} } while( strcmp(charp_[stride*j],pivot)> 0 );
    if( i >= j ) break;
    if (verbose){ printf(" %% i %d j %d swapping\n",i,j);}
    tmpcharp = charp_[stride*i]; charp_[stride*i] = charp_[stride*j]; charp_[stride*j] = tmpcharp;
    if (index_!=NULL){ tmpi = index_[i]; index_[i] = index_[j]; index_[j] = tmpi;}
  }while(1);
  if (verbose){ printf(" %% l %d j %d swapping\n",l,j);}
  tmpcharp = charp_[stride*l]; charp_[stride*l] = charp_[stride*j]; charp_[stride*j] = tmpcharp;
  if (index_!=NULL){ tmpi = index_[l]; index_[l] = index_[j]; index_[j] = tmpi;}
  if (verbose){ printf(" %% [finished charpPartition_index] returning %d\n",j);}
  return j;
}

unsigned int charpQuickSort_index(unsigned int nn,char **charp_,int stride,int *index_,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  int verbose=0;
  int j=0; unsigned int n1=nn,n2=nn;
  if (verbose){ printf(" %% [entering charpQuickSort_index] nn %d stride %d l %d r %d\n",nn,stride,l,r);}
  if( l < r ) {
    if (verbose){ printf(" %% l<r: calling charpPartition_index\n");}
    j = charpPartition_index(charp_,stride,index_,l,r); 
    if (verbose){ printf(" %% done calling charpPartition_index\n");}
    if (nn<GLOBAL_recursion_limit){ 
      if (verbose){ printf(" %% splitting...\n");}
      n1 = charpQuickSort_index(nn+1,charp_,stride,index_,l,j-1); 
      n2 = charpQuickSort_index(nn+1,charp_,stride,index_,j+1,r); 
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in charpQuickSort_index.\n",nn);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  if (verbose){ printf(" %% [finished charpQuickSort_index] returning %d\n",maximum(n1,n2));}
  return maximum(n1,n2);
}

void charpQuickSort_index_driver(int n_charp,char **charp_,int stride,char **charp_workspace_,int *index_)
{
  /* finds index listing index_ so that charp_[stride*index_[.]] is sorted. */
  int ncharp=0;
  for (ncharp=0;ncharp<n_charp;ncharp++){ charp_workspace_[ncharp] = charp_[stride*ncharp]; index_[ncharp]=ncharp;}
  charpQuickSort_index(0,charp_workspace_,1,index_,0,n_charp-1);
}

void charpQuickSort_index_test()
{
  int n_charp = 5,ncharp=0;
  int index_[n_charp];
  char *charp_[2*n_charp];
  char *charp_workspace_[1*n_charp];
  charp_[0] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[0],"dog");
  charp_[1] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[1]," ");
  charp_[2] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[2],"cat");
  charp_[3] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[3]," ");
  charp_[4] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[4],"bird");
  charp_[5] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[5]," ");
  charp_[6] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[6],"pig");
  charp_[7] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[7]," ");
  charp_[8] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[8],"bug");
  charp_[9] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[9]," ");
  charpQuickSort_index_driver(n_charp,charp_,2,charp_workspace_,index_);
  for (ncharp=0;ncharp<2*n_charp;ncharp++){ printf(" %% ncharp %d: %s\n",ncharp,charp_[ncharp]);}
  raprintf(index_,"int",1,n_charp," %% index_: ");
}

void charpQuickSort_index_index_driver(int n_charp,char **charp_,int stride,char **charp_workspace_,int *index_orig_from_sort_,int *index_sort_from_orig_,int *index_workspace_)
{
  /* finds index listing index_orig_from_sort_ so that charp_[stride*index_orig_from_sort_[.]] is sorted. */
  int verbose=0;
  int ncharp=0;
  if (verbose){ printf(" %% [entering charpQuickSort_index_index_driver]\n");}
  for (ncharp=0;ncharp<n_charp;ncharp++){
    if (verbose){ printf(" %% ncharp %d/%d\n",ncharp,n_charp);}
    charp_workspace_[ncharp] = charp_[stride*ncharp];
    if (index_orig_from_sort_!=NULL){ index_orig_from_sort_[ncharp]=ncharp;}
    /* for (ncharp=0;ncharp<n_charp;ncharp++){ } */}
  if (verbose){ printf(" %% calling charpQuickSort_index\n");}
  charpQuickSort_index(0,charp_workspace_,1,index_orig_from_sort_,0,n_charp-1);
  if (verbose){ printf(" %% finished charpQuickSort_index\n");}
  if (index_sort_from_orig_!=NULL){
    if (verbose){ printf(" %% setting index_sort_from_orig\n");}
    for (ncharp=0;ncharp<n_charp;ncharp++){ index_workspace_[ncharp] = index_orig_from_sort_[ncharp]; index_sort_from_orig_[ncharp]=ncharp;}
    lQuickSort_index(0,index_workspace_,1,index_sort_from_orig_,0,n_charp-1);
    /* if (index_sort_from_orig_!=NULL){ } */}
  if (verbose){ printf(" %% [finished charpQuickSort_index_index_driver]\n");}
}

void charpQuickSort_index_index_test()
{
  int n_charp = 5,ncharp=0;
  int index_orig_from_sort_[n_charp];
  int index_sort_from_orig_[n_charp];
  char *charp_[2*n_charp];
  char *charp_workspace_[1*n_charp];
  int *index_workspace_[1*n_charp];
  charp_[0] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[0],"dog");
  charp_[1] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[1]," ");
  charp_[2] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[2],"cat");
  charp_[3] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[3]," ");
  charp_[4] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[4],"bird");
  charp_[5] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[5]," ");
  charp_[6] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[6],"pig");
  charp_[7] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[7]," ");
  charp_[8] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[8],"bug");
  charp_[9] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_[9]," ");
  charpQuickSort_index_index_driver(n_charp,charp_,2,charp_workspace_,index_orig_from_sort_,index_sort_from_orig_,index_workspace_);
  for (ncharp=0;ncharp<2*n_charp;ncharp++){ printf(" %% ncharp %d: %s\n",ncharp,charp_[ncharp]);}
  raprintf(index_orig_from_sort_,"int",1,n_charp," %% index_orig_from_sort_: ");
  raprintf(index_sort_from_orig_,"int",1,n_charp," %% index_sort_from_orig_: ");
}

void charpIntersect_index_index_driver(int n_charp_0,char **charp_0_,int stride_0,int n_charp_1,char **charp_1_,int stride_1,char **charp_workspace_,int *n_intersect_p,int *index_0_orig_from_sort_,int *index_0_orig_from_intersect_,int *index_0_intersect_from_orig_,int *index_1_orig_from_sort_,int *index_1_orig_from_intersect_,int *index_1_intersect_from_orig_)
{
  /* This intersection retains repeated entries with a multiplicity given by the minimum multiplicity across inputs. */
  /* The indexed entries point to the first instances of each element in the intersection */
  /* We expect charp_workspace_ to be of size maximum(n_charp_0,n_charp_1). */
  /* charp_workspace_ ends up holding the intersection. */
  int verbose=1;
  int n_intersect=0,nintersect=0;
  int n0_sort=0,n1_sort=0;
  int n0_orig=0,n1_orig=0;
  char *charp_0=NULL,*charp_1=NULL;
  int tmp_compare=0;
  if (verbose){ printf(" %% [entering charpIntersect_index_index_driver]\n");}
  /* first sort charp_0_ and charp_1_ */
  charpQuickSort_index_index_driver(n_charp_0,charp_0_,stride_0,charp_workspace_,index_0_orig_from_sort_,NULL,NULL);
  if (verbose>1){ for (n0_sort=0;n0_sort<n_charp_0;n0_sort++){ n0_orig = index_0_orig_from_sort_[n0_sort]; printf(" %% (%d,%d) <-- %s\n",n0_sort,n0_orig,charp_0_[stride_0*n0_orig]);} /* if (verbose){ } */}
  charpQuickSort_index_index_driver(n_charp_1,charp_1_,stride_1,charp_workspace_,index_1_orig_from_sort_,NULL,NULL);
  if (verbose>1){ for (n1_sort=0;n1_sort<n_charp_1;n1_sort++){ n1_orig = index_1_orig_from_sort_[n1_sort]; printf(" %% (%d,%d) <-- %s\n",n1_sort,n1_orig,charp_1_[stride_1*n1_orig]);} /* if (verbose){ } */}
  /* now step through both charp_0_ and charp_1_ to count number of identical entries */
  if (index_0_intersect_from_orig_!=NULL){ for (n0_orig=0;n0_orig<n_charp_0;n0_orig++){ index_0_intersect_from_orig_[n0_orig] = -1;}}
  if (index_1_intersect_from_orig_!=NULL){ for (n1_orig=0;n1_orig<n_charp_1;n1_orig++){ index_1_intersect_from_orig_[n1_orig] = -1;}}
  n_intersect = 0;
  n0_sort = 0; if (n0_sort<n_charp_0){ n0_orig = index_0_orig_from_sort_[n0_sort]; charp_0 = charp_0_[stride_0*n0_orig];}
  n1_sort = 0; if (n1_sort<n_charp_1){ n1_orig = index_1_orig_from_sort_[n1_sort]; charp_1 = charp_1_[stride_1*n1_orig];}
  while ((n0_sort<n_charp_0) && (n1_sort<n_charp_1)){
    if (verbose>1){ printf(" %% %d %s %d %s\n",n0_sort,charp_0,n1_sort,charp_1);}
    tmp_compare = strcmp(charp_0,charp_1);
    if (tmp_compare< 0){ /* charp_0 <  charp_1 */
      //charp_workspace_[n_intersect] = charp_0; //<-- use for union. ;
      //index_0_orig_from_intersect_[n_intersect] = n0_orig; //<-- use for union. ;
      n0_sort += 1; if (n0_sort<n_charp_0){ n0_orig = index_0_orig_from_sort_[n0_sort]; charp_0 = charp_0_[stride_0*n0_orig];}
      //n_intersect++; //<-- use for union. ;
      /* if (tmp_compare<=0){ } */}
    if (tmp_compare==0){ /* charp_0 == charp_1 */
      charp_workspace_[n_intersect] = charp_0;
      if (index_0_orig_from_intersect_!=NULL){ index_0_orig_from_intersect_[n_intersect] = n0_orig;}
      if (index_0_intersect_from_orig_!=NULL){ index_0_intersect_from_orig_[n0_orig] = n_intersect;}
      if (index_1_orig_from_intersect_!=NULL){ index_1_orig_from_intersect_[n_intersect] = n1_orig;}
      if (index_1_intersect_from_orig_!=NULL){ index_1_intersect_from_orig_[n1_orig] = n_intersect;}
      n0_sort += 1; if (n0_sort<n_charp_0){ n0_orig = index_0_orig_from_sort_[n0_sort]; charp_0 = charp_0_[stride_0*n0_orig];}
      n1_sort += 1; if (n1_sort<n_charp_1){ n1_orig = index_1_orig_from_sort_[n1_sort]; charp_1 = charp_1_[stride_1*n1_orig];}
      n_intersect++;
      /* if (tmp_compare==0){ } */}
    if (tmp_compare> 0){ /* charp_0 >  charp_1 */
      //charp_workspace_[n_intersect] = charp_1; //<-- use for union. ;
      //index_1_orig_from_intersect_[n_intersect] = n1_orig; //<-- use for union. ;
      n1_sort += 1; if (n1_sort<n_charp_1){ n1_orig = index_1_orig_from_sort_[n1_sort]; charp_1 = charp_1_[stride_1*n1_orig];}
      //n_intersect++; //<-- use for union. ;
      /* if (tmp_compare> 0){ } */}
    /* while ((n0_sort<n_charp_0) && (n1_sort<n_charp_1)){ } */}
  if (n_intersect_p!=NULL){ *n_intersect_p = n_intersect;}
  if (verbose>0){
    if (verbose>0){ printf(" %% Testing intersection... \n");}
    for (nintersect=0;nintersect<n_intersect;nintersect++){
      n0_orig = index_0_orig_from_intersect_[nintersect];
      n1_orig = index_1_orig_from_intersect_[nintersect];
      if ( (strcmp(charp_0_[stride_0*n0_orig],charp_1_[stride_1*n1_orig])!=0) || (strcmp(charp_workspace_[nintersect],charp_0_[stride_0*n0_orig])!=0) || (strcmp(charp_workspace_[nintersect],charp_1_[stride_1*n1_orig])!=0) ){
	printf(" %% Warning, nintersect %d (%s) n0_orig %d (%s) n1_orig %d (%s)\n",nintersect,charp_workspace_[nintersect],n0_orig,charp_0_[stride_0*n0_orig],n1_orig,charp_1_[stride_1*n1_orig]);
	/* if error */}
      if ( (strcmp(charp_workspace_[index_0_intersect_from_orig_[n0_orig]],charp_0_[stride_0*n0_orig])!=0) || (strcmp(charp_workspace_[index_1_intersect_from_orig_[n1_orig]],charp_1_[stride_1*n1_orig])!=0) ){
	printf(" %% Warning, index_0_intersect_from_orig_[n0_orig] %d (%s) n0_orig %d (%s)\n",index_0_intersect_from_orig_[n0_orig],charp_workspace_[index_0_intersect_from_orig_[n0_orig]],n0_orig,charp_0_[stride_0*n0_orig]);
	printf(" %% Warning, index_1_intersect_from_orig_[n1_orig] %d (%s) n1_orig %d (%s)\n",index_1_intersect_from_orig_[n1_orig],charp_workspace_[index_1_intersect_from_orig_[n1_orig]],n1_orig,charp_1_[stride_1*n1_orig]);
	/* if error */}
      /* for (nintersect=0;nintersect<n_intersect;nintersect++){ } */}
    if (verbose>0){ printf(" %% Finished testing %d elements. \n",n_intersect);}
    /* if (verbose>0){ } */}
  wkspace_printf();
  if (verbose){ printf(" %% [finished charpIntersect_index_index_driver]\n");}
}

void charpIntersect_index_index_test()
{
  int n_charp_0 = 5,ncharp_0=0,stride_0=2;
  int *index_0_orig_from_sort_=NULL;
  int *index_0_orig_from_intersect_=NULL;
  int *index_0_intersect_from_orig_=NULL;
  char **charp_0_=NULL;
  int n_charp_1 = 6,ncharp_1=0,stride_1=3;
  int *index_1_orig_from_sort_=NULL;
  int *index_1_orig_from_intersect_=NULL;
  int *index_1_intersect_from_orig_=NULL;
  char **charp_1_=NULL;
  char **charp_workspace_=NULL;
  int nl=0;
  int n_intersect=0,nintersect=0;
  charp_workspace_ = (char **)wkspace_all0c(sizeof(char *)*maximum(n_charp_0,n_charp_1));
  index_0_orig_from_sort_ = (int *) wkspace_all0c(sizeof(int)*n_charp_0);
  index_0_orig_from_intersect_ = (int *) wkspace_all0c(sizeof(int)*n_charp_0);
  index_0_intersect_from_orig_ = (int *) wkspace_all0c(sizeof(int)*n_charp_0);
  charp_0_ = (char **)wkspace_all0c(sizeof(char *)*stride_0*n_charp_0);
  nl=0;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"dog"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"cat"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"bird"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"bug"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"bug"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  index_1_orig_from_sort_ = (int *) wkspace_all0c(sizeof(int)*n_charp_1);
  index_1_orig_from_intersect_ = (int *) wkspace_all0c(sizeof(int)*n_charp_1);
  index_1_intersect_from_orig_ = (int *) wkspace_all0c(sizeof(int)*n_charp_1);
  charp_1_ = (char **)wkspace_all0c(sizeof(char *)*stride_1*n_charp_1);
  nl=0;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"pear"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"cat"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"bug"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"wolf"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"bug"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"dog"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charpIntersect_index_index_driver(n_charp_0,charp_0_,stride_0,n_charp_1,charp_1_,stride_1,charp_workspace_,&n_intersect,index_0_orig_from_sort_,index_0_orig_from_intersect_,index_0_intersect_from_orig_,index_1_orig_from_sort_,index_1_orig_from_intersect_,index_1_intersect_from_orig_);
  for (ncharp_0=0;ncharp_0<stride_0*n_charp_0;ncharp_0++){ printf(" %% ncharp_0 %d: %s\n",ncharp_0,charp_0_[ncharp_0]);}
  for (ncharp_1=0;ncharp_1<stride_1*n_charp_1;ncharp_1++){ printf(" %% ncharp_1 %d: %s\n",ncharp_1,charp_1_[ncharp_1]);}
  for (nintersect=0;nintersect<n_intersect;nintersect++){ printf(" %% nintersect %d: %s\n",nintersect,charp_workspace_[nintersect]);}
  raprintf(index_0_orig_from_intersect_,"int",1,n_intersect," %% index_0_orig_from_intersect_: ");
  raprintf(index_0_intersect_from_orig_,"int",1,n_charp_0," %% index_0_intersect_from_orig_: ");
  raprintf(index_1_orig_from_intersect_,"int",1,n_intersect," %% index_1_orig_from_intersect_: ");
  raprintf(index_1_intersect_from_orig_,"int",1,n_charp_1," %% index_1_intersect_from_orig_: ");
  exit(0);
}

void charpIntersectall_index_index_driver(int n_list,int *n_charp_,char ***charp__,int *stride_,char **charp_workspace_,int *n_intersect_p,int **index_orig_from_sort__,int **index_orig_from_intersect__,int **index_intersect_from_orig__)
{
  /* This intersection retains repeated entries with a multiplicity given by the minimum multiplicity across inputs. */
  /* The indexed entries point to the first instances of each element in the intersection */
  /* We expect charp_workspace_ to be of size maximum(n_charp_). */
  /* charp_workspace_ ends up holding the intersection. */
  int verbose=1;
  int nlist=0;
  int n_intersect=0,nintersect=0;
  int n_orig=0,n_sort=0;
  int n_orig_[n_list];
  int n_sort_[n_list];
  char *charp_[n_list];
  int flag_continue=0;
  int flag_found=0;
  int tmp_compare=0;
  int nlist_min=0;
  if (verbose){ printf(" %% [entering charpIntersectall_index_index_driver]\n");}
  /* first sort each charp__[nlist] */
  for (nlist=0;nlist<n_list;nlist++){
    if (verbose>1){ printf(" %% nlist %d/%d calling QuickSort\n",nlist,n_list);}
    charpQuickSort_index_index_driver(n_charp_[nlist],charp__[nlist],stride_[nlist],charp_workspace_,index_orig_from_sort__[nlist],NULL,NULL);
    if (verbose>1){ for (n_sort=0;n_sort<minimum(10,n_charp_[nlist]);n_sort++){ n_orig = index_orig_from_sort__[nlist][n_sort]; printf(" %% nlist %d: (%d,%d) <-- %s\n",nlist,n_sort,n_orig,charp__[nlist][stride_[nlist]*n_orig]);} /* if (verbose){ } */}
    /* for (nlist=0;nlist<n_list;nlist++){ } */}
  /* now step through all the charp__[nlist] to count number of identical entries */
  for (nlist=0;nlist<n_list;nlist++){
    if (verbose>1){ printf(" %% nlist %d/%d setting index_intersect_from_orig__\n",nlist,n_list);}
    if (index_intersect_from_orig__[nlist]!=NULL){ for (n_orig=0;n_orig<n_charp_[nlist];n_orig++){ index_intersect_from_orig__[nlist][n_orig] = -1;}}
    /* for (nlist=0;nlist<n_list;nlist++){ } */}
  if (verbose>1){ printf(" %% stepping through sorted lists.\n");}
  n_intersect = 0;
  for (nlist=0;nlist<n_list;nlist++){
    n_sort_[nlist] = 0;
    if (n_sort_[nlist]<n_charp_[nlist]){
      n_orig_[nlist] = index_orig_from_sort__[nlist][n_sort_[nlist]]; charp_[nlist] = charp__[nlist][stride_[nlist]*n_orig_[nlist]];
      /* if (n_sort_[nlist]<n_charp_[nlist]){ } */}
    /* for (nlist=0;nlist<n_list;nlist++){ } */}
  flag_continue=1;
  for (nlist=0;nlist<n_list;nlist++){
    flag_continue *= (n_sort_[nlist]<n_charp_[nlist] ? 1 : 0);
    /* for (nlist=0;nlist<n_list;nlist++){ } */}
  while (flag_continue){
    nlist_min=0; flag_found=0;
    for (nlist=0;nlist<n_list;nlist++){
      tmp_compare = strcmp(charp_[nlist_min],charp_[nlist]);
      if (tmp_compare>0){ flag_found=1; nlist_min=nlist;}
      if (tmp_compare<0){ flag_found=1; nlist_min=nlist_min;}
      /* for (nlist=0;nlist<n_list;nlist++){ } */}
    if (flag_found==1){ n_sort_[nlist_min] += 1; if (n_sort_[nlist_min]<n_charp_[nlist_min]){ n_orig_[nlist_min] = index_orig_from_sort__[nlist_min][n_sort_[nlist_min]]; charp_[nlist_min] = charp__[nlist_min][stride_[nlist_min]*n_orig_[nlist_min]];}}
    if (flag_found==0){
      charp_workspace_[n_intersect] = charp_[0];
      for (nlist=0;nlist<n_list;nlist++){
	if (index_orig_from_intersect__[nlist]!=NULL){ index_orig_from_intersect__[nlist][n_intersect] = n_orig_[nlist];}
	if (index_intersect_from_orig__[nlist]!=NULL){ index_intersect_from_orig__[nlist][n_orig_[nlist]] = n_intersect;}
	n_sort_[nlist] += 1; if (n_sort_[nlist]<n_charp_[nlist]){ n_orig_[nlist] = index_orig_from_sort__[nlist][n_sort_[nlist]]; charp_[nlist] = charp__[nlist][stride_[nlist]*n_orig_[nlist]];}
	/* for (nlist=0;nlist<n_list;nlist++){ } */}
      n_intersect++;
      /* if (flag_found==0){ } */}
    flag_continue=1;
    for (nlist=0;nlist<n_list;nlist++){
      flag_continue *= (n_sort_[nlist]<n_charp_[nlist] ? 1 : 0);
      /* for (nlist=0;nlist<n_list;nlist++){ } */}
    /* while (flag_continue){ } */}
  if (n_intersect_p!=NULL){ *n_intersect_p = n_intersect;}
  if (verbose>1){ printf(" %% found n_intersect %d.\n",n_intersect);}
  if (verbose>0){
    printf(" %% Testing intersection... \n");
    for (nintersect=0;nintersect<n_intersect;nintersect++){
      for (nlist=0;nlist<n_list;nlist++){
	n_orig_[nlist] = index_orig_from_intersect__[nlist][nintersect];
	if (strcmp(charp_workspace_[nintersect],charp__[nlist][stride_[nlist]*n_orig_[nlist]])!=0){
	  printf(" %% Warning, nintersect %d (%s) nlist %d n_orig %d (%s)\n",nintersect,charp_workspace_[nintersect],nlist,n_orig_[nlist],charp__[nlist][stride_[nlist]*n_orig_[nlist]]);
	  /* if error */}
	if (strcmp(charp_workspace_[index_intersect_from_orig__[nlist][n_orig_[nlist]]],charp__[nlist][stride_[nlist]*n_orig_[nlist]])!=0){
	  printf(" %% Warning, index_intersect_from_orig__[nlist][n_orig_[nlist]] %d (%s) n_orig_[nlist] %d (%s)\n",index_intersect_from_orig__[nlist][n_orig_[nlist]],charp_workspace_[index_intersect_from_orig__[nlist][n_orig_[nlist]]],n_orig_[nlist],charp__[nlist][stride_[nlist]*n_orig_[nlist]]);
	  /* if error */}
	/* for (nlist=0;nlist<n_list;nlist++){ } */}
      /* for (nintersect=0;nintersect<n_intersect;nintersect++){ } */}
    printf(" %% Finished testing %d elements. \n",n_intersect);
    /* if (verbose>0){ } */}
  wkspace_printf();
  if (verbose){ printf(" %% [finished charpIntersectall_index_index_driver]\n");}
}

void charpIntersectall_index_index_test()
{
  int n_list = 3,nlist=0;
  int n_charp_0 = 5,ncharp_0=0,stride_0=2;
  int index_0_orig_from_sort_[n_charp_0];
  int index_0_orig_from_intersect_[n_charp_0];
  int index_0_intersect_from_orig_[n_charp_0];
  char *charp_0_[2*n_charp_0];
  int n_charp_1 = 6,ncharp_1=0,stride_1=3;
  int index_1_orig_from_sort_[n_charp_1];
  int index_1_orig_from_intersect_[n_charp_1];
  int index_1_intersect_from_orig_[n_charp_1];
  char *charp_1_[3*n_charp_1];
  int n_charp_2 = 4,ncharp_2=0,stride_2=2;
  int index_2_orig_from_sort_[n_charp_2];
  int index_2_orig_from_intersect_[n_charp_2];
  int index_2_intersect_from_orig_[n_charp_2];
  char *charp_2_[stride_2*n_charp_2];
  char *charp_workspace_[1*maximum(maximum(n_charp_0,n_charp_1),n_charp_2)];
  int nl=0;
  int n_intersect=0,nintersect=0;
  int *index_orig_from_sort__[n_list];
  int *index_orig_from_intersect__[n_list];
  int *index_intersect_from_orig__[n_list];
  int n_charp_[n_list];
  int stride_[n_list];
  char **charp__[n_list];
  int n_orig=0,n_sort=0;
  nl=0;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"dog"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"cat"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"bird"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"bug"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl],"bug"); nl++;
  charp_0_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_0_[nl]," "); nl++;
  nl=0;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"pear"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"cat"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"bug"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"wolf"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"bug"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl],"dog"); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  charp_1_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_1_[nl]," "); nl++;
  nl=0;
  charp_2_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_2_[nl],"bug"); nl++;
  charp_2_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_2_[nl]," "); nl++;
  charp_2_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_2_[nl],"bird"); nl++;
  charp_2_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_2_[nl]," "); nl++;
  charp_2_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_2_[nl],"bug"); nl++;
  charp_2_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_2_[nl]," "); nl++;
  charp_2_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_2_[nl],"dog"); nl++;
  charp_2_[nl] = (char *)wkspace_all0c(sizeof(char)*8); sprintf(charp_2_[nl]," "); nl++;
  index_orig_from_sort__[0] = index_0_orig_from_sort_;
  index_orig_from_sort__[1] = index_1_orig_from_sort_;
  index_orig_from_sort__[2] = index_2_orig_from_sort_;
  index_orig_from_intersect__[0] = index_0_orig_from_intersect_;
  index_orig_from_intersect__[1] = index_1_orig_from_intersect_;
  index_orig_from_intersect__[2] = index_2_orig_from_intersect_;
  index_intersect_from_orig__[0] = index_0_intersect_from_orig_;
  index_intersect_from_orig__[1] = index_1_intersect_from_orig_;
  index_intersect_from_orig__[2] = index_2_intersect_from_orig_;
  n_charp_[0] = n_charp_0;
  n_charp_[1] = n_charp_1;
  n_charp_[2] = n_charp_2;
  stride_[0] = stride_0;
  stride_[1] = stride_1;
  stride_[2] = stride_2;
  charp__[0] = charp_0_;
  charp__[1] = charp_1_;
  charp__[2] = charp_2_;
  charpIntersectall_index_index_driver(n_list,n_charp_,charp__,stride_,charp_workspace_,&n_intersect,index_orig_from_sort__,index_orig_from_intersect__,index_intersect_from_orig__);
  for (nlist=0;nlist<n_list;nlist++){
    for (n_orig=0;n_orig<n_charp_[nlist];n_orig++){ printf(" %% nlist %d n_orig %d: %s\n",nlist,n_orig,charp__[nlist][stride_[nlist]*n_orig]);}
    raprintf(index_orig_from_intersect__[nlist],"int",1,n_intersect," %% index_orig_from_intersect__: ");
    raprintf(index_intersect_from_orig__[nlist],"int",1,n_charp_[nlist]," %% index_intersect_from_orig__: ");
    /* for (nlist=0;nlist<n_list;nlist++){ } */}
  for (nintersect=0;nintersect<n_intersect;nintersect++){ printf(" %% nintersect %d: %s\n",nintersect,charp_workspace_[nintersect]);}
  exit(0);
}
