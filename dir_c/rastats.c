void ra_stats(void *vra,char *type,int length,void *max_p,void *min_p,double *mean_p,double *stdev_p)
{
  /* finds the stats of an array */
  int i=0;
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

double dra_sum(double *dra,int length)
{
  /* sums array of double type */
  double *dstop=NULL,dsum=0;
  dstop=&(dra[length]); dsum=0; while (dra<dstop){ dsum+=*(dra++);}
  return dsum;
}

void dra_clean_s_x_m_m(double *dra1,int nrows,int rpop_b,int rpop_j,int ncols,int cpop_b,int cpop_j,unsigned char *bmr_b,unsigned char *bmr_j,unsigned char *bmc_b,unsigned char *bmc_j,int spacing_r,int spacing_c)
{
  /* clears dra1 outside of masks bmr_j,bmc_j assuming spacing_r,spacing_c */ 
  int nr_j=0,nr_b=0,nr_a=0; int nc_j=0,nc_b=0,nc_a=0; int tab_r=0,tab_c=0,tab_r_stride=0,tab_x=0;
  if (dra1!=NULL && bmr_b!=NULL && bmr_j!=NULL && bmc_b!=NULL && bmc_j!=NULL){ 
    switch (spacing_r){ case SPACING_j: tab_r_stride = rpop_j; break; case SPACING_b: tab_r_stride = rpop_b; break; case SPACING_a: tab_r_stride = nrows; break; default: break;}
    /* switch (spacing_c){ case SPACING_j: tab_c_stride = cpop_j; break; case SPACING_b: tab_c_stride = cpop_b; break; case SPACING_a: tab_c_stride = ncols; break; default: break;} */
    nc_j=0;nc_b=0;nc_a=0;
    while (nc_a<ncols){
      if (bget__on(bmc_b,nc_a)){
	if (bget__on(bmc_j,nc_a)){
	  switch (spacing_c){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break;}
	  nr_j=0;nr_b=0;nr_a=0;
	  while (nr_a<nrows){
	    if (bget__on(bmr_b,nr_a)){
	      if (bget__on(bmr_j,nr_a)){
		switch (spacing_r){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break;}
		tab_x = tab_r + tab_c*tab_r_stride;
		/* do nothing */
		nr_j++; /* if (bget__on(bmr_j,nr_a)){ } */}
	      else /* if not in bmr_j */{
		switch (spacing_r){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break;}
		tab_x = tab_r + tab_c*tab_r_stride;
		if (spacing_r==SPACING_b || spacing_r==SPACING_a){ dra1[tab_x]=0;}
		/* if not in bmr_j */}
	      nr_b++; /* if (bget__on(bmr_b,nr_a)){ } */}
	    else /* if not in bmr_b */{
	      switch (spacing_r){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break;}
	      tab_x = tab_r + tab_c*tab_r_stride;
	      if (spacing_r==SPACING_a){ dra1[tab_x]=0;}
	      /* if not in bmr_b */}
	    nr_a++; /* while (nr_a<nrows){ } */}
	  nc_j++; /* if (bget__on(bmc_j,nc_a)){ } */}
	else /* if not in bmc_j */{
	  switch (spacing_c){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break;}
	  if (spacing_c==SPACING_b || spacing_c==SPACING_a){ for (nr_a=0;nr_a<tab_r_stride;nr_a++){ tab_x = nr_a + tab_c*tab_r_stride; dra1[tab_x]=0;}}
	  /* if not in bmc_j */}
	nc_b++; /* if (bget__on(bmc_b,nc_a)){ } */}
      else /* if not in bmc_b */{
	switch (spacing_c){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break;}
	if (spacing_c==SPACING_a){ for (nr_a=0;nr_a<tab_r_stride;nr_a++){ tab_x = nr_a + tab_c*tab_r_stride; dra1[tab_x]=0;}}
	/* if not in bmc_b */}
      nc_a++;/* while (nc_a<ncols){ } */}
    /* if (dra1!=NULL && bmr_b!=NULL && bmr_j!=NULL && bmc_b!=NULL && bmc_j!=NULL){ } */}
}
void dra_plustimesequals_s_x_m_m(double *dra1,int nrows,int rpop_b,int rpop_j,int ncols,int cpop_b,int cpop_j,double *dra2,double multiplier,unsigned char *bmr_b,unsigned char *bmr_j,unsigned char *bmc_b,unsigned char *bmc_j,int spacing_r,int spacing_c)
{
  /* adds multiplier*dra2 to dra1 along masks bmr_j,bmc_j assuming spacing_r,spacing_c */ 
  int nr_j=0,nr_b=0,nr_a=0; int nc_j=0,nc_b=0,nc_a=0; int tab_r=0,tab_c=0,tab_r_stride=0,tab_x=0;
  if (dra2!=NULL && dra1!=NULL && bmr_b!=NULL && bmr_j!=NULL && bmc_b!=NULL && bmc_j!=NULL){ 
    switch (spacing_r){ case SPACING_j: tab_r_stride = rpop_j; break; case SPACING_b: tab_r_stride = rpop_b; break; case SPACING_a: tab_r_stride = nrows; break; default: break;}
    /* switch (spacing_c){ case SPACING_j: tab_c_stride = cpop_j; break; case SPACING_b: tab_c_stride = cpop_b; break; case SPACING_a: tab_c_stride = ncols; break; default: break;} */
    nc_j=0;nc_b=0;nc_a=0;
    while (nc_a<ncols){
      if (bget__on(bmc_b,nc_a)){
	if (bget__on(bmc_j,nc_a)){
	  switch (spacing_c){ case SPACING_j: tab_c = nc_j; break; case SPACING_b: tab_c = nc_b; break; case SPACING_a: tab_c = nc_a; break; default: break;}
	  nr_j=0;nr_b=0;nr_a=0;
	  while (nr_a<nrows){
	    if (bget__on(bmr_b,nr_a)){
	      if (bget__on(bmr_j,nr_a)){
		switch (spacing_r){ case SPACING_j: tab_r = nr_j; break; case SPACING_b: tab_r = nr_b; break; case SPACING_a: tab_r = nr_a; break; default: break;}
		tab_x = tab_r + tab_c*tab_r_stride;
		dra1[tab_x] += multiplier*dra2[tab_x];
		nr_j++; /* if (bget__on(bmr_j,nr_a)){ } */}
	      nr_b++; /* if (bget__on(bmr_b,nr_a)){ } */}
	    nr_a++; /* while (nr_a<nrows){ } */}
	  nc_j++; /* if (bget__on(bmc_j,nc_a)){ } */}
	nc_b++; /* if (bget__on(bmc_b,nc_a)){ } */}
      nc_a++;/* while (nc_a<ncols){ } */}
    /* if (dra2!=NULL && dra1!=NULL && bmr_b!=NULL && bmr_j!=NULL && bmc_b!=NULL && bmc_j!=NULL){ } */}
}
void dra_plustimesequals_s___m_m(double *dra1,int nrows,int ncols,double *dra2,double multiplier,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* adds multiplier*dra2 to dra1 along masks bmr_b,bmc_b assuming SPACING_a */ 
  int nr_a=0; int nc_a=0; 
  if (dra2!=NULL && dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL){ for (nc_a=0;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ dra1[nr_a+nc_a*nrows] += multiplier*dra2[nr_a+nc_a*nrows];}}}}}
}
void dra_plusdivequals_s___m_m(double *dra1,int nrows,int ncols,double *dra2,double denominator,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* adds dra2/denominator to dra1 along masks bmr_b,bmc_b assuming SPACING_a */ 
  int nr_a=0; int nc_a=0; 
  double d_use = denominator ;
  if (fabs(d_use)<GLOBAL_tolerance){ d_use=1;}
  if (dra2!=NULL && dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL){ for (nc_a=0;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ dra1[nr_a+nc_a*nrows] += dra2[nr_a+nc_a*nrows]/d_use;}}}}}
}
void dra_aucflip_s___m_m(double *dra1,int nrows,int ncols,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* replaces dra1 with max(dra1,1-dra1) along masks bmr_b,bmc_b assuming SPACING_a */ 
  int nr_a=0; int nc_a=0; 
  if (dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL){ for (nc_a=0;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ dra1[nr_a+nc_a*nrows] = maximum(dra1[nr_a+nc_a*nrows],1-dra1[nr_a+nc_a*nrows]);}}}}}
}
void dra_min1_s___m_m(double *dra1,int nrows,int ncols,double *dra2,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* replaces dra1 with min(dra1,dra2) along masks bmr_b,bmc_b assuming SPACING_a */ 
  int nr_a=0; int nc_a=0; 
  if (dra2!=NULL && dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL){ for (nc_a=0;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ dra1[nr_a+nc_a*nrows] = minimum(dra1[nr_a+nc_a*nrows],dra2[nr_a+nc_a*nrows]);}}}}}
}
void dra_pow2r_s___m_m(double *dra1,int nrows,int ncols,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* replaces dra1 with pow(max(dra1,0),2) along masks bmr_b,bmc_b assuming SPACING_a */ 
  int nr_a=0; int nc_a=0; 
  if (dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL){ for (nc_a=0;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ dra1[nr_a+nc_a*nrows] = maximum(dra1[nr_a+nc_a*nrows],0)*maximum(dra1[nr_a+nc_a*nrows],0);}}}}}
}
void dra_mdsnrm_s___m_m(double *dra1,int nrows,int ncols,double denominator,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* replaces first column of dra1 with dra1(:,1) - sum(dra1(:,2:end),2)/denominator */
  /* assuming SPACING_a */
  int nr_a=0; int nc_a=0;
  if (dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL && denominator>0){ 
    for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ for (nc_a=1;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ dra1[nr_a+0*nrows] -= dra1[nr_a+nc_a*nrows]/denominator;}}}}
    /* if denominator>0 */}  
}
void dra_mdspow_s___m_m(double *dra1,int nrows,int ncols,double denominator,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* if dra1(j,1)>0 then replaces dra1(j,1) with sgn(tmp)*sqrt(abs(tmp)), with tmp = dra1(j,1).^2 - sum(dra1(j,2:end).^2,2)/denominator */
  /* assuming SPACING_a */
  int nr_a=0; int nc_a=0;
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
void dra_mdsabs_s___m_m(double *dra1,int nrows,int ncols,double denominator,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* replaces dra1(j,1) with dra1(j,1) - sqrt(mean(dra1(j,2:end).^2))/denominator */
  /* assuming SPACING_a */
  int nr_a=0; int nc_a=0; int ntmp=0; double dtmp=0;
  if (dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL && denominator>0){ 
    for (nr_a=0;nr_a<nrows;nr_a++){ 
      if (bget__on(bmr_b,nr_a)){ 
	ntmp=0; dtmp=0;
	for (nc_a=1;nc_a<ncols;nc_a++){ 
	  if (bget__on(bmc_b,nc_a)){ ntmp++; dtmp += dra1[nr_a+nc_a*nrows]*dra1[nr_a+nc_a*nrows];}
	  /* for (nc_a=1;nc_a<ncols;nc_a++){ } */}
	if (ntmp>0){ dtmp /= ntmp;} dtmp = sqrt(dtmp); 
	dra1[nr_a+0*nrows] -= dtmp/denominator;
	/* if (bget__on(bmr_b,nr_a)){ } */}
      /* for (nr_a=0;nr_a<nrows;nr_a++){ } */}
    /* if denominator>0 */}
}
void dra_sqrtr_s___m_m(double *dra1,int nrows,int ncols,unsigned char *bmr_b,unsigned char *bmc_b)
{
  /* replaces dra1 with sqrt(max(dra1,0)) along masks bmr_b,bmc_b assuming SPACING_a */ 
  int nr_a=0; int nc_a=0; 
  if (dra1!=NULL && bmr_b!=NULL && bmc_b!=NULL){ for (nc_a=0;nc_a<ncols;nc_a++){ if (bget__on(bmc_b,nc_a)){ for (nr_a=0;nr_a<nrows;nr_a++){ if (bget__on(bmr_b,nr_a)){ dra1[nr_a+nc_a*nrows] = sqrt(maximum(dra1[nr_a+nc_a*nrows],0));}}}}}
}
void dra_plustimesequals(double *dra1,int length,double *dra2,double multiplier)
{
  /* adds multiplier*dra2 to dra1 */ int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] += multiplier*dra2[nr];}
}
void dra_plusequals(double *dra1,int length,double *dra2)
{
  /* adds dra2 to dra1 */ int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] += dra2[nr];}
}
void dra_subtequals(double *dra1,int length,double *dra2)
{
  /* subtracts dra2 from dra1 */ int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] -= dra2[nr];}
}
void dra_plus(double *dra1,int length,double adder)
{
  /* adds adder to dra1 */ int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] += adder;}
}
void dra_times(double *dra,int length,double multiplier)
{
  /* multiplies dra by multiplier */ int nr=0; for (nr=0;nr<length;nr++){ dra[nr] *= multiplier;}
}
void dra_dup(double *dra1,int length,double *dra2)
{
  /* sets dra1 to dra2 */ int nr=0; for (nr=0;nr<length;nr++){ dra1[nr] = dra2[nr];}
}
void ura_dup(unsigned char *ura1,int length,unsigned char *ura2)
{
  /* sets ura1 to ura2 */ int nr=0; for (nr=0;nr<length;nr++){ ura1[nr] = ura2[nr];}
}

void dra_normalize(int type_flag,double *Tra,double *Sra,int nrows_T,int nrows_S,int ncols,double *p_,double *q_,double *a_,double *b_)
{
  /* type_flag==TYPE_CAT: categorical normalization; each column beta-centered ;
     type_flag==TYPE_CNT: continuous normalization; each column normalized ;
  */
  int nc=0;
  double sum_T=0,sum_S=0,p=0,q=0,alpha=0,beta=0,mean_T=0,stdev_T=0,mean_S=0,stdev_S=0,mean=0,stdev=0;
  switch (type_flag){
  case TYPE_CAT: /* categorical */
    for (nc=0;nc<ncols;nc++){
      /* ra_sum(&(Tra[0+nc*nrows_T]),"double",nrows_T,&sum_T); */
      /* ra_sum(&(Sra[0+nc*nrows_S]),"double",nrows_S,&sum_S); */
      sum_T = dra_sum(&(Tra[0+nc*nrows_T]),nrows_T);
      sum_S = dra_sum(&(Sra[0+nc*nrows_S]),nrows_S);
      p = (sum_T+sum_S)/(nrows_T+nrows_S); q = 1-p; alpha = p-q; beta = 2*p*q;
      dra_times(&(Tra[0+nc*nrows_T]),nrows_T,2);
      dra_plus(&(Tra[0+nc*nrows_T]),nrows_T,-1-alpha);
      dra_times(&(Sra[0+nc*nrows_S]),nrows_S,2);
      dra_plus(&(Sra[0+nc*nrows_S]),nrows_S,-1-alpha);
      if (p_!=NULL){ p_[nc]=p;} if (q_!=NULL){ q_[nc]=q;}
      if (a_!=NULL){ a_[nc]=alpha;} if (b_!=NULL){ b_[nc]= beta;}
      /* for (nc=0;nc<ncols;nc++){ } */}
    break;
  case TYPE_CNT: /* continuous */
    for (nc=0;nc<ncols;nc++){
      ra_stats(&(Tra[0+nc*nrows_T]),"double",nrows_T,NULL,NULL,&mean_T,&stdev_T);
      ra_stats(&(Sra[0+nc*nrows_S]),"double",nrows_S,NULL,NULL,&mean_S,&stdev_S);
      mean = (nrows_T*mean_T + nrows_S*mean_S)/(nrows_T+nrows_S);
      stdev = sqrt( ((stdev_T*stdev_T + mean_T*mean_T)*nrows_T + (stdev_S*stdev_S + mean_S*mean_S)*nrows_S)/(nrows_T+nrows_S) - mean*mean);
      dra_plus(&(Tra[0+nc*nrows_T]),nrows_T,-mean);
      dra_plus(&(Sra[0+nc*nrows_S]),nrows_S,-mean);
      dra_times(&(Tra[0+nc*nrows_T]),nrows_T,1.0/stdev);
      dra_times(&(Sra[0+nc*nrows_S]),nrows_S,1.0/stdev);
      if (p_!=NULL){ p_[nc]=0;} if (q_!=NULL){ q_[nc]=0;}
      if (a_!=NULL){ a_[nc]=0;} if (b_!=NULL){ b_[nc]=1;}
      /* for (nc=0;nc<ncols;nc++){ } */}
    break;
  default: printf(" %% warning! unknown type %d in dra_normalize\n",type_flag); break; /* switch (type_flag){ } */}
}

void dra_to_llra(double *Tra,double *Sra,int nrows_T,int nrows_S,int ncols,long long int *Tll,long long int *Sll,long long int *denominator_p)
{
  double dmin_T=0,dmax_T=0,dmin_S,dmax_S=0,dex=0;
  int ldex=0,nl;
  long long int denominator=0;
  ra_stats(Tra,"double",nrows_T*ncols,&dmax_T,&dmin_T,NULL,NULL);
  ra_stats(Sra,"double",nrows_S*ncols,&dmax_S,&dmin_S,NULL,NULL);
  dex = maximum(maximum(fabs(dmax_T),fabs(dmin_T)),maximum(fabs(dmax_S),fabs(dmin_S)));
  ldex = maximum(16,ceil(log2(dex)));
  denominator = 1; while (0<ldex--){ denominator *= 2;}
  for (nl=0;nl<nrows_T*ncols;nl++){ Tll[nl] = llround(Tra[nl]*denominator);}
  for (nl=0;nl<nrows_S*ncols;nl++){ Sll[nl] = llround(Sra[nl]*denominator);}
  if (denominator_p!=NULL){ *denominator_p=denominator;}  
}

double dra_norm(double *dra,int length,int stride)
{
  double output=0,*dstop=&(dra[length*stride]);
  while (dra<dstop){ output+=(*(dra))*(*(dra)); dra+=stride;}
  return sqrt(output);
}

double dra_diff(double *dra,double *drb,int length,int stride)
{
  double output=0,dtmp=0,*dstop=&(dra[length*stride]);
  while(dra<dstop){ dtmp=(*(dra))-(*(drb)); output+=dtmp*dtmp; dra+=stride;drb+=stride;}
  return sqrt(output);
}
