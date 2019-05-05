#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

unsigned long int lrand()
{
  /* generates a unsigned long int (64 bit) */
  int verbose=0;
  int sr=(int)round(log(RAND_MAX)/log(2)),sl=8*sizeof(unsigned long int),sld=8*sizeof(long double);
  int nshift=sl/sr + (sl - sr*(sl/sr) > 0),ns=0;
  unsigned long int ld=0;
  if (verbose){ printf(" %% rand() size %d unsigned long int size %d, factor %d, long double size %d\n",sr,sl,nshift,sld);}
  ld=0;
  for (ns=0;ns<nshift;ns++){ 
    ld <<= sr; ld+=abs(rand()); 
    if (verbose){ printf(" %% ld set to %ld 2^(%0.1f) out of %d\n",ld,(double)logl((long double)ld)/log(2),sl);}
    /* for (ns=0;ns<nshift;ns++){ } */}
  return ld;
}

double randn()
{
  /* box muller polar form */
  double u=0,v=0,s=0;
  while (s==0 || s>1){ u=2*rand01-1;v=2*rand01-1;s=u*u+v*v;}
  return u*sqrt(-2*log(s)/s);
}

unsigned long int RGET(unsigned long int *rseed)
{
  /* basically:
     RNEXT = (RPREV*((unsigned long int)(pow(2,RPOW)+RADD))%((unsigned long int) pow(2,2*RPOW-1)));
     return RNEXT */
  /* test via:
     unsigned long int rs=lrand();
     int nbins=64,htau=1,tab=0,nda=0,ndb=0;
     double *dra=NULL,da=0,db=0;
     struct matrix *h1ra=NULL,*h2ra=NULL;
     h1ra = mmake();minit0(h1ra,1,nbins);
     h2ra = mmake();minit0(h2ra,nbins,nbins);
     dra = (double *) tcalloc(htau+1,sizeof(double));
     for (nl=0;nl<htau+1;nl++){ dra[nl]=R01GET(&rs);} tab=0;
     for (nl=0;nl<1000000;nl++){ 
     dra[tab]=R01GET(&rs); da=dra[tab]; tab+=1; if (tab==htau+1){tab=0;} db=dra[tab]; 
     nda = crop((int)floor(nbins*da),0,nbins-1); ndb = crop((int)floor(nbins*db),0,nbins-1);
     sentry(h1ra,0,nda,gentry(h1ra,0,nda)+1);
     sentry(h2ra,nda,ndb,gentry(h2ra,nda,ndb)+1);} 
     raprintf(h1ra->mtrx,"double",1,nbins," %% h1ra: ");
     raprintf(h2ra->mtrx,"double",nbins,nbins," %% h2ra: ");
     mspy(h2ra,1.0);
     exit(0);
  */
  *rseed = (*rseed*POW2RPOWPLUSRADD%POW22RPOWMINUSONE);
  return *rseed;
}

double R01GET(unsigned long int *rseed)
{
 /* basically:
    RNEXT = (RPREV*((unsigned long int)(pow(2,RPOW)+RADD))%((unsigned long int) pow(2,2*RPOW-1)));
    return = (double)RNEXT/(double)pow(2,2*RPOW-1); */
  int verbose=0;
  *rseed = (*rseed*POW2RPOWPLUSRADD%POW22RPOWMINUSONE);
  if (verbose){ printf(" %% %ld, %ld/%ld, %lf\n",*rseed,*rseed,POW22RPOWMINUSONE,(double)((long double)*rseed/(long double)POW22RPOWMINUSONE));}
  return (double)((long double)*rseed/(long double)POW22RPOWMINUSONE);
}

double RNGET(unsigned long int *rseed)
{
  /* box muller polar form */
  double u=0,v=0,s=0;
  while (s==0 || s>1){ u=2*R01GET(rseed)-1;v=2*R01GET(rseed)-1;s=u*u+v*v;}
  return u*sqrt(-2*log(s)/s);
}

double RISIGET(unsigned long int *rseed,double rate)
{
  double r = R01GET(rseed);
  if (r==0){ r=1;}
  return -log(r)/rate;
}
