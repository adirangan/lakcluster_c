#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void ira_cumulative_sum(int *ira,unsigned long long int length,int **jra_p)
{ 
  /* cumulative sum */
  unsigned long long int i=0;
  int sum=0;
  if (*jra_p==NULL){ *jra_p = (int *) wkspace_all0c(length*sizeof(int));}
  for (i=0;i<length;i++){ sum += ira[i]; (*jra_p)[i]=sum;}
}

void ra_stats(void *vra,char *type,unsigned long long int length,void *max_p,void *min_p,double *mean_p,double *stdev_p)
{
  /* finds the stats of an array */
  unsigned long long int i=0;
  double mean=0,stdev=0;
  double dmax=0,dmin=0;
  long long int llmax=0,llmin=0;
  int imax=0,imin=0;
  if (vra!=NULL){
    if (strcmp(type,"double")==0){
      for (i=0;i<length;i++){ mean += ((double *)vra)[i];}
      mean /= (double)length;
      dmax=((double *)vra)[0]; dmin=((double *)vra)[0];
      for (i=0;i<length;i++){
	if (((double *)vra)[i]>dmax){ dmax=((double *)vra)[i];}
	if (((double *)vra)[i]<dmin){ dmin=((double *)vra)[i];}
	stdev += (((double *)vra)[i]-mean)*(((double *)vra)[i]-mean);}
      if (max_p!=NULL){*(double *)max_p=dmax;}
      if (min_p!=NULL){*(double *)min_p=dmin;}
      stdev /= (double)length;
      stdev = sqrt(stdev);}
    else if (strcmp(type,"long long int")==0){
      for (i=0;i<length;i++){ mean += ((long long int *)vra)[i];}
      mean /= (double)length;
      llmax=((long long int *)vra)[0]; llmin=((long long int *)vra)[0];
      for (i=0;i<length;i++){
	if (((long long int *)vra)[i]>llmax){ llmax=((long long int *)vra)[i];}
	if (((long long int *)vra)[i]<llmin){ llmin=((long long int *)vra)[i];}
	stdev += (((long long int *)vra)[i]-mean)*(((long long int *)vra)[i]-mean);}
      if (max_p!=NULL){*(long long int *)max_p=llmax;}
      if (min_p!=NULL){*(long long int *)min_p=llmin;}
      stdev /= (double)length;
      stdev = sqrt(stdev);}
    else if (strcmp(type,"int")==0){
      for (i=0;i<length;i++){ mean += ((int *)vra)[i];}
      mean /= (double)length;
      imax=((int *)vra)[0]; imin=((int *)vra)[0];
      for (i=0;i<length;i++){
	if (((int *)vra)[i]>imax){ imax=((int *)vra)[i];}
	if (((int *)vra)[i]<imin){ imin=((int *)vra)[i];}
	stdev += (((int *)vra)[i]-mean)*(((int *)vra)[i]-mean);}
      if (max_p!=NULL){*(int *)max_p=imax;}
      if (min_p!=NULL){*(int *)min_p=imin;}
      stdev /= (double)length;
      stdev = sqrt(stdev);}
    if (mean_p!=NULL){*mean_p=mean;}
    if (stdev_p!=NULL){*stdev_p=stdev;}}
}

void dra_plustimesequals_s___m_m(double *dra1,unsigned long long int nrows,unsigned long long int ncols,double *dra2,double multiplier,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* adds multiplier*dra2 to dra1 along masks bmr_b,bmc_b assuming SPACING_a */ 
  unsigned long long int nr_a=0; unsigned long long int nc_a=0; 
  if (dra2!=NULL && dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL){ for (nc_a=0;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ dra1[nr_a+nc_a*nrows] += multiplier*dra2[nr_a+nc_a*nrows];}}}}}
}
void dra_plusdivequals_s___m_m(double *dra1,unsigned long long int nrows,unsigned long long int ncols,double *dra2,double denominator,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* adds dra2/denominator to dra1 along masks bmr_b,bmc_b assuming SPACING_a */ 
  unsigned long long int nr_a=0; unsigned long long int nc_a=0; 
  double d_use = denominator ;
  if (fabs(d_use)<GLOBAL_tolerance){ d_use=1;}
  if (dra2!=NULL && dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL){ for (nc_a=0;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ dra1[nr_a+nc_a*nrows] += dra2[nr_a+nc_a*nrows]/d_use;}}}}}
}
void dra_mds_nrm_s___m_m(double *dra1,unsigned long long int nrows,unsigned long long int ncols,double denominator,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* replaces first column of dra1 with dra1(:,1) - sum(dra1(:,2:end),2)/denominator */
  /* assuming SPACING_a */
  unsigned long long int nr_a=0; unsigned long long int nc_a=0;
  if (dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL && denominator>0){ 
    for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ for (nc_a=1;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ dra1[nr_a+0*nrows] -= dra1[nr_a+nc_a*nrows]/denominator;}}}}
    /* if denominator>0 */}  
}
void dra_mds_pow_s___m_m(double *dra1,unsigned long long int nrows,unsigned long long int ncols,double denominator,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* if dra1(j,1)>0 then replaces dra1(j,1) with sgn(tmp)*sqrt(abs(tmp)), with tmp = dra1(j,1).^2 - sum(dra1(j,2:end).^2,2)/denominator */
  /* assuming SPACING_a */
  unsigned long long int nr_a=0; unsigned long long int nc_a=0;
  if (dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL && denominator>0){ 
    for (nr_a=0;nr_a<nrows;nr_a++){ 
      if (bget__on(bmr_b,nr_a)){ 
	if (dra1[nr_a+0*nrows]>=0){ 
	  dra1[nr_a+0*nrows] *= dra1[nr_a+0*nrows]; 
	  for (nc_a=1;nc_a<ncols;nc_a++){ 
	    if (bget__on(bmc_b,nc_a)){ dra1[nr_a+0*nrows] -= dra1[nr_a+nc_a*nrows]*dra1[nr_a+nc_a*nrows]/denominator;}
	    /* for (nc_a=1;nc_a<ncols;nc_a++){ } */}
	  dra1[nr_a+0*nrows] = dra1[nr_a+0*nrows]>=0 ? sqrt(dra1[nr_a+0*nrows]) : -sqrt(-dra1[nr_a+0*nrows]);
	  /* if (dra1[nr_a+0*nrows]>=0){ } */}
	/* if (bget__on(bmr_b,nr_a)){ } */}
      /* for (nr_a=0;nr_a<nrows;nr_a++){ } */}
    /* if denominator>0 */}
}

void dra_sumx(unsigned long long int nx,double *dra,unsigned long long int stride,unsigned long long int length)
{
  /* sums length entries into entry nx, separated by stride */ unsigned long long int nr=0; for (nr=0;nr<nx;nr++){ dra[nx] += dra[nr*stride];} for (nr=nx+1;nr<length;nr++){ dra[nx] += dra[nr*stride];}
}

void dra_plustimesequals(double *dra1,unsigned long long int length,double *dra2,double multiplier)
{
  /* adds multiplier*dra2 to dra1 */ unsigned long long int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] += multiplier*dra2[nr];}
}
void dra_plusequals(double *dra1,unsigned long long int length,double *dra2)
{
  /* adds dra2 to dra1 */ unsigned long long int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] += dra2[nr];}
}
void dra_subtequals(double *dra1,unsigned long long int length,double *dra2)
{
  /* subtracts dra2 from dra1 */ unsigned long long int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] -= dra2[nr];}
}
void dra_plus(double *dra1,unsigned long long int length,double adder)
{
  /* adds adder to dra1 */ unsigned long long int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] += adder;}
}
void dra_times(double *dra,unsigned long long int length,double multiplier)
{
  /* multiplies dra by multiplier */ unsigned long long int nr=0; for (nr=0;nr<length;nr++){ dra[nr] *= multiplier;}
}
void dra_dup(double *dra1,unsigned long long int length,double *dra2)
{
  /* sets dra1 to dra2 */ unsigned long long int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] = dra2[nr];}
}
void ura_dup(unsigned char *ura1,unsigned long long int length,unsigned char *ura2)
{
  /* sets ura1 to ura2 */ unsigned long long int nr=0; for (nr=0;nr<length;nr++){ ura1[nr] = ura2[nr];}
}

double dra_diff(double *dra,double *drb,unsigned long long int length,unsigned long long int stride)
{
  double output=0,dtmp=0,*dstop=&(dra[length*stride]);
  while(dra<dstop){ dtmp=(*(dra))-(*(drb)); output+=dtmp*dtmp; dra+=stride;drb+=stride;}
  return sqrt(output);
}
