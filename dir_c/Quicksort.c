
int dPartition(double *dra,int stride,int l,int r) 
{
  double pivot=0,tmpd=0;
  int i=0,j=0;
  pivot = dra[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( dra[stride*i] <= pivot && i <= r );
    do{ j--;} while( dra[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = dra[stride*i]; dra[stride*i] = dra[stride*j]; dra[stride*j] = tmpd;
  }while(1);
  tmpd = dra[stride*l]; dra[stride*l] = dra[stride*j]; dra[stride*j] = tmpd;
  return j;
}

void dQuickSort(unsigned int nn,double *dra,int stride,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  /* test with: */
  /*
    double dra[] = { 7.0 , 8.0 , 6.0 , 3.0 , 2.0 , 9.0 , 10.0 , 11.0 , 6.0 , 7.0 , 6.0 };
    raprintf(dra,"double",1,11," %% dra: ");
    dQuickSort(dra,ira,ira_b,ira_j,0,10);
    raprintf(dra,"double",1,11," %% dra: ");
  */
  int j=0; 
  if( l < r ) { 
    if (nn<GLOBAL_recursion_limit){ 
      j = dPartition(dra,stride,l,r); 
      dQuickSort(nn+1,dra,stride,l,j-1); 
      dQuickSort(nn+1,dra,stride,j+1,r);
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* if (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in dQuickSort\n",nn);
      /* recursion_limit breached */}  
    /* if( l < r ) { } */}
}

int dPartition_xij(double *dra,int stride,int *ira,int *ira_b,int *ira_j,int *lnb,int *mr_a,int *mr_b,int l,int r) 
{
  double pivot=0,tmpd=0;
  int i=0,j=0,tmpi=0;
  pivot = dra[stride*l];
  i = l; j = r+1;
  do{
    do{ i++;} while( dra[stride*i] <= pivot && i <= r );
    do{ j--;} while( dra[stride*j] > pivot );
    if( i >= j ) break;
    tmpd = dra[stride*i]; dra[stride*i] = dra[stride*j]; dra[stride*j] = tmpd;
    if (ira!=NULL){ tmpi = ira[i]; ira[i] = ira[j]; ira[j] = tmpi;}
    if (ira_b!=NULL){ tmpi = ira_b[i]; ira_b[i] = ira_b[j]; ira_b[j] = tmpi;}
    if (ira_j!=NULL){ tmpi = ira_j[i]; ira_j[i] = ira_j[j]; ira_j[j] = tmpi;}
    if (lnb!=NULL){ tmpi = lnb[i]; lnb[i] = lnb[j]; lnb[j] = tmpi;}
    if (mr_a!=NULL){ tmpi = mr_a[i]; mr_a[i] = mr_a[j]; mr_a[j] = tmpi;}
    if (mr_b!=NULL){ tmpi = mr_b[i]; mr_b[i] = mr_b[j]; mr_b[j] = tmpi;}
  }while(1);
  tmpd = dra[stride*l]; dra[stride*l] = dra[stride*j]; dra[stride*j] = tmpd;
  if (ira!=NULL){ tmpi = ira[l]; ira[l] = ira[j]; ira[j] = tmpi;}
  if (ira_b!=NULL){ tmpi = ira_b[l]; ira_b[l] = ira_b[j]; ira_b[j] = tmpi;}
  if (ira_j!=NULL){ tmpi = ira_j[l]; ira_j[l] = ira_j[j]; ira_j[j] = tmpi;}
  if (lnb!=NULL){ tmpi = lnb[l]; lnb[l] = lnb[j]; lnb[j] = tmpi;}
  if (mr_a!=NULL){ tmpi = mr_a[l]; mr_a[l] = mr_a[j]; mr_a[j] = tmpi;}
  if (mr_b!=NULL){ tmpi = mr_b[l]; mr_b[l] = mr_b[j]; mr_b[j] = tmpi;}
  return j;
}

unsigned int dQuickSort_xij(unsigned int nn,double *dra,int stride,int *ira,int *ira_b,int *ira_j,int *lnb,int *mr_a,int *mr_b,int l,int r)
{
  /* modified from http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c */
  /* test with: */
  /*
    double dra[] = { 7.0 , 8.0 , 6.0 , 3.0 , 2.0 , 9.0 , 10.0 , 11.0 , 6.0 , 7.0 , 6.0 };
    int    ira[] = {   0 ,   1 ,   2 ,   3 ,   4 ,   5 ,    6 ,    7 ,   8 ,   9 ,  10 };
    int  ira_b[] = {   0 ,   1 ,   2 ,   3 ,   4 ,   5 ,    6 ,    7 ,   8 ,   9 ,  10 };
    int  ira_j[] = {   0 ,   1 ,   2 ,   3 ,   4 ,   5 ,    6 ,    7 ,   8 ,   9 ,  10 };
    int   lnb[] = {   0 ,   1 ,   0 ,   1 ,   0 ,   1 ,    0 ,    1 ,   0 ,   1 ,   0 };
    raprintf(dra,"double",1,11," %% dra: ");
    raprintf(ira,"int",1,11," %% ira: ");
    raprintf(lnb,"int",1,11," %% lnb: ");
    dQuickSort_xij(dra,ira,ira_b,ira_j,0,10);
    raprintf(dra,"double",1,11," %% dra: ");
    raprintf(ira,"int",1,11," %% ira: ");
    raprintf(lnb,"int",1,11," %% lnb: ");
  */
  int j=0; unsigned int n1=nn,n2=nn;
  if( l < r ) { 
    j = dPartition_xij(dra,stride,ira,ira_b,ira_j,lnb,mr_a,mr_b,l,r); 
    if (nn<GLOBAL_recursion_limit){ 
      n1 = dQuickSort_xij(nn+1,dra,stride,ira,ira_b,ira_j,lnb,mr_a,mr_b,l,j-1); 
      n2 = dQuickSort_xij(nn+1,dra,stride,ira,ira_b,ira_j,lnb,mr_a,mr_b,j+1,r); 
      /* if (nn<GLOBAL_recursion_limit){ } */}
    else /* (nn>=GLOBAL_recursion_limit) */{ 
      printf(" %% Warning! recursion limit %d reached in dQuickSort_xij.\n %% This means that many row- and col-scores are identical.\n %% It is likely that your data is mostly zero or constant.\n",nn);
      /* recursion_limit breached */}
    /* if( l < r ) { } */}
  return maximum(n1,n2);
}
