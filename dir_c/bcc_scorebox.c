#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void bcc_scorebox_mxA(struct bcc_ajdk *D,int rdrop,int cdrop)
{
  /* reading rdrop and cdrop from input: ;
     copies D->A_bmc_j to D->A_umc_j ; 
     define D->A_umc_j_rmv and D->A_umc_j_rtn ; removing first cdrop entries of D->QC_index_local_mc_a ;
     copies E->A_bmr_j to E->A_umr_j ; 
     define E->A_umr_j_rmv and E->A_umr_j_rtn ; removing first rdrop entries of D->QR_index_local_mr_a ; using D->QR_index_local_nb to index nb ;
     call bcc_sumscores_mxB ; copying all instances of umc and umr to associated bmc and bmr respectively ;
   */
  int verbose=0; 
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_single *E=NULL;
  int nb1=0,nr=0,nc=0;
  if (verbose){ printf(" %% [entering bcc_scorebox_mxA]\n");}
  for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = D->A_umc_j[nc]; D->A_umc_j_rmv[nc] = 0;}
  for (nc=0;nc<minimum(cdrop,D->A_cpop_j);nc++){ D->A_umc_j_rmv[D->QC_index_local_mc_a[nc]]=1; D->A_umc_j_rtn[D->QC_index_local_mc_a[nc]]=0;}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = E->A_umr_j[nr]; E->A_umr_j_rmv[nr] = 0;}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nr=0;nr<minimum(rdrop,D->A_rpop_j_total);nr++){ E_[D->QR_index_local_nb[nr]]->A_umr_j_rmv[D->QR_index_local_mr_a[nr]]=1; E_[D->QR_index_local_nb[nr]]->A_umr_j_rtn[D->QR_index_local_mr_a[nr]]=0;}
  if (verbose>1){
    if (verbose){ raprintf(D->A_umc_j    ,"char",1,D->A_ncols," %% D->A_umc_j    : ");}
    if (verbose){ raprintf(D->A_umc_j_rmv,"char",1,D->A_ncols," %% D->A_umc_j_rmv: ");}
    if (verbose){ raprintf(D->A_umc_j_rtn,"char",1,D->A_ncols," %% D->A_umc_j_rtn: ");}
    for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j    : ",nb1); raprintf(E->A_umr_j    ,"char",1,E->A_nrows,D->tmpAnchar);
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j_rmv: ",nb1); raprintf(E->A_umr_j_rmv,"char",1,E->A_nrows,D->tmpAnchar);
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j_rtn: ",nb1); raprintf(E->A_umr_j_rtn,"char",1,E->A_nrows,D->tmpAnchar);
      /* for (nb1=0;nb1<nbins;nb1++){ } */}
    /* if (verbose>1){ } */}
  if (verbose>1){ printf(" %% calling bcc_sumscores_mxB\n");}
  bcc_sumscores_mxB(D); 
  if (verbose>1){ printf(" %% finished bcc_sumscores_mxB\n");}
  if (verbose>1){
    sprintf(D->tmpAnchar," %%%% D->A_bmc_j    : "); bprintf(D->A_bmc_j    ,D->bitj,1,D->A_ncols,D->tmpAnchar);
    sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rmv: "); bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols,D->tmpAnchar);
    sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rtn: "); bprintf(D->A_bmc_j_rtn,D->bitj,1,D->A_ncols,D->tmpAnchar);
    for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j    : ",nb1); bprintf(E->A_bmr_j    ,D->bitj,1,E->A_nrows,D->tmpAnchar);
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j_rmv: ",nb1); bprintf(E->A_bmr_j_rmv,D->bitj,1,E->A_nrows,D->tmpAnchar);
      sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j_rtn: ",nb1); bprintf(E->A_bmr_j_rtn,D->bitj,1,E->A_nrows,D->tmpAnchar);
      /* for (nb1=0;nb1<nbins;nb1++){ } */}
    /* if (verbose>1){ } */}
  if (verbose){ printf(" %% [finished bcc_scorebox_mxA]\n");}
}

void bcc_scorebox_svalue(struct bcc_ajdk *D)
{
  /* Calculates D->QR_svalue and D->QC_svalue from E->QR_AnAtTAnAt_nrm and D->QC_AtTAnAtTAn_nrm, respectively. */
  int verbose=0; 
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb1=0; struct bcc_single *E=NULL;
  int na_a=0,na_b=0,na_j=0,ma_a=0,ma_b=0,ma_j=0;
  int Irow=0,Icol=0;
  if (verbose){ printf(" %% [entering bcc_scorebox_svalue]\n");}
  D->Irem=0; for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; D->Irem += (E->A_rpop_j>0?1:0);} 
  xcc_get_Ireq("empty low",nbins,D->Irem,&(D->Ireq),&(Irow),&(Icol));
  if (verbose){ printf(" %% found D->Irem %d/%d, using D->Ireq %d; Irow %d Icol %d\n",D->Irem,nbins,D->Ireq,Irow,Icol);}
  if (verbose){
    printf(" %% D->A_ncols %d D->A_cpop_b %d D->A_cpop_j %d\n",D->A_ncols,D->A_cpop_b,D->A_cpop_j);
    if (verbose>1){ sprintf(D->tmpAnchar," %%%% D->A_bmc_b: "); bprintf(D->A_bmc_b,D->bitj,1,D->A_ncols,D->tmpAnchar);}
    if (verbose>1){ sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
    /* if (verbose){ } */}
  if (D->A_cbother && D->A_rpop_b_total){
    na_a=0;na_b=0;na_j=0;
    while (na_a<D->A_ncols){
      if (bget__on(D->A_bmc_b,na_a)){
	if (bget__on(D->A_bmc_j,na_a)){
	  D->QC_svalue[na_j] = D->QC_AtTAnAtTAn_nrm[na_a + 0*D->A_ncols + Icol*D->A_ncols*D->T_ncols];
	  na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	na_b++; /* if (bget__on(D->A_bmc_b,na_a)){ } */}
      na_a++;/* while (na_a<D->A_ncols){ } */}
    if (verbose>1){
      raprintf(  D->QC_svalue,"double",1,D->A_cpop_j," %%   D->QC_svalue :");
      /* if (verbose){ } */}
    /* if (D->A_cbother && D->A_rpop_b_total){ } */}
  if (D->A_rpop_b_total && D->A_cbother){
    D->A_nrows_total=0;D->A_rpop_b_total=0;D->A_rpop_j_total=0;
    for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; 
      if (verbose){
	printf(" %% nb1 %d E->A_nrows %d E->A_rpop_b %d E->A_rpop_j %d\n",nb1,E->A_nrows,E->A_rpop_b,E->A_rpop_j);
	if (verbose>1){ sprintf(D->tmpAnchar," %%%% E->A_bmr_b: "); bprintf(E->A_bmr_b,D->bitj,1,E->A_nrows,D->tmpAnchar);}
	if (verbose>1){ sprintf(D->tmpAnchar," %%%% E->A_bmr_j: "); bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
	/* if (verbose){ } */}
      ma_a=0;ma_b=0;ma_j=0;
      while (ma_a<E->A_nrows){
	if (bget__on(E->A_bmr_b,ma_a)){
	  if (bget__on(E->A_bmr_j,ma_a)){
	    D->QR_svalue[D->A_rpop_j_total] = E->QR_AnAtTAnAt_nrm[ma_a + 0*E->A_nrows + Irow*E->A_nrows*D->T_ncols];
	    ma_j++; D->A_rpop_j_total++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	  ma_b++; D->A_rpop_b_total++; /* while (ma_a<E->A_nrows){ } */}
	ma_a++; D->A_nrows_total++; /* while (ma_a<E->A_nrows){ } */}
      /* for (nb1=0;nb1<nbins;nb1++){ } */}
    if (verbose>1){
      raprintf(  D->QR_svalue,"double",1,D->A_rpop_j_total," %% D->QR_svalue :");
      /* if (verbose){ } */}
    /* if (D->A_rpop_b_total && D->A_cbother){ } */}
  if (verbose){ printf(" %% [finished bcc_scorebox_svalue]\n");}
}

void bcc_scorebox_srt(struct bcc_ajdk *D,int nrows,int *mr_index_local_nb,int *mr_index_local_mr,int ncols,int *mc_index_sort)
{
  /* Copies last D->A_rpop_j_total entries from mr_index_sort (i.e., mr_index_local_nb and mr_index_local_mr) ;
     into D->QR_index_local_nb and D->QR_index_local_mr_a, respectively. ;
     Copies last D->A_cpop_j entries from mc_index_sort ;
     into D->QC_index_local_mc_a. ;
  */
  int verbose=0; 
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb1=0; struct bcc_single *E=NULL;
  int na_j=0,nr=0;
  if (verbose){ printf(" %% [entering bcc_scorebox_srt]\n");}
  if (mc_index_sort!=NULL && ncols>=D->A_cpop_j){ 
    for (na_j=0;na_j<D->A_cpop_j;na_j++){ D->QC_index_local_mc_a[na_j] = mc_index_sort[ncols-D->A_cpop_j+na_j];} 
    /* if (mc_index_sort!=NULL){ } */}
  if (mr_index_local_nb!=NULL && mr_index_local_mr!=NULL && nrows>=D->A_rpop_j_total){ 
    for (nr=0;nr<D->A_rpop_j_total;nr++){ 
      D->QR_index_local_mr_a[nr] = mr_index_local_mr[nrows-D->A_rpop_j_total+nr];
      D->QR_index_local_nb[nr] = mr_index_local_nb[nrows-D->A_rpop_j_total+nr];
      /* for (nr=0;nr<D->A_rpop_j_total;nr++){ } */}
    /* if (mr_index_local_nb!=NULL && mr_index_local_mr!=NULL){ } */};
  if (verbose>1){ printf(" %% finished copying mr_index_local_nb,mr_index_local_mr and mc_index_sort\n");}    
  if (verbose>1){
    raprintf(D->QR_index_local_mr_a,   "int",1,D->A_rpop_j_total," %% D->QR_index_local_mr_a pos :");
    raprintf(  D->QR_index_local_nb,   "int",1,D->A_rpop_j_total," %%   D->QR_index_local_nb pos :");
    raprintf(D->QC_index_local_mc_a,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_a pos :");
    /* if (verbose){ } */}
  if (verbose){ printf(" %% [finished bcc_scorebox_srt]\n");}
}

void bcc_out_xdrop_lkp(struct bcc_ajdk *D,int nrows,int *mr_index_sort,int **mr_index_local_nb_,int **mr_index_local_mr_)
{
  /* Converts mr_index_sort to lnb (local nb) and index_local_mr (local mr) indices. ;
     Note that all indices are of 'a' type (i.e., referring to mr_a). */
  int verbose=0;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nl=0,tmpd=0,nb=0,mx_a_tot=0;
  struct bcc_single *E=NULL;
  if (verbose){ printf(" %% [entering bcc_out_xdrop_lkp]\n");}
  mx_a_tot = 0; for (nb=0;nb<nbins;nb++){ E = E_[nb]; mx_a_tot += E->A_nrows; /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% nbins %d, mx_a_tot %d, nrows %d\n",nbins,mx_a_tot,nrows);}
  if (verbose){ printf(" %% allocating mr_index_local_nb_\n");}
  if (*mr_index_local_nb_==NULL){ (*mr_index_local_nb_) = (int *) wkspace_all0c(nrows*sizeof(int));}
  if (verbose){ printf(" %% allocating mr_index_local_mr_\n");}
  if (*mr_index_local_mr_==NULL){ (*mr_index_local_mr_) = (int *) wkspace_all0c(nrows*sizeof(int));}
  if (verbose){ printf(" %% setting mr_index_local_xx_\n");}
  for (nl=0;nl<nrows;nl++){
    if (verbose){ printf(" %% %% nl %d/%d\n",nl,nrows);}
    tmpd = mr_index_sort[nl];
    if (verbose){ printf(" %% %% nl %d/%d: tmpd %d\n",nl,nrows,tmpd);}
    nb=0; while (tmpd>=E_[nb]->A_nrows){ tmpd-=E_[nb]->A_nrows; nb++;}
    if (verbose){ printf(" %% %% nl %d/%d: tmpd %d nb %d\n",nl,nrows,tmpd,nb);}
    (*mr_index_local_nb_)[nl] = nb; (*mr_index_local_mr_)[nl] = tmpd;
    /* for (nl=0;nl<nrows;nl++){ } */}
  if (verbose>1){ raprintf(*mr_index_local_nb_,"int",1,nrows," %% mr_index_local_nb: ");}
  if (verbose>1){ raprintf(*mr_index_local_mr_,"int",1,nrows," %% mr_index_local_mr: ");}
  if (verbose){ printf(" %% [finished bcc_out_xdrop_lkp]\n");}
}

void xcc_out_xdrop_sort(int nlines,int *out_xdrop,int *nrows_,int **mr_index_sort_,int *ncols_,int **mc_index_sort_)
{
  /* converts out_xdrop into sorted lists of indices */
  int verbose=0;
  int nl=0,nrmax=0,ncmax=0,nr=0,nc=0;
  if (verbose){ printf(" %% [entering xcc_out_xdrop_sort]\n");}
  (*nrows_)=0;(*ncols_)=0; nrmax=0;ncmax=0; for (nl=0;nl<nlines;nl++){ if (out_xdrop[0 + nl*2]>=0){ nrmax = maximum(nrmax,out_xdrop[0 + nl*2]); (*nrows_)++;} if (out_xdrop[1 + nl*2]>=0){ ncmax = maximum(ncmax,out_xdrop[1 + nl*2]); (*ncols_)++;}} if ((*nrows_)!=nrmax+1 || (*ncols_)!=ncmax+1){ if (verbose>-1){ printf(" %% Note: (*nrows_) %d (*ncols_) %d nrmax %d ncmax %d in xcc_out_xdrop_sort\n",(*nrows_),(*ncols_),nrmax,ncmax);}}
  if (verbose){ printf(" %% (*nrows_) %d (*ncols_) %d\n",(*nrows_),(*ncols_));}
  if (*mr_index_sort_==NULL){ (*mr_index_sort_) = (int *) wkspace_all0c((*nrows_)*sizeof(int));}
  if (*mc_index_sort_==NULL){ (*mc_index_sort_) = (int *) wkspace_all0c((*ncols_)*sizeof(int));}
  nr=0; nl=0; while (nr<(*nrows_) && nl<nlines){ if (out_xdrop[0 + nl*2]>=0){ (*mr_index_sort_)[nr] = out_xdrop[0 + nl*2]; nr++;} nl++;}
  nc=0; nl=0; while (nc<(*ncols_) && nl<nlines){ if (out_xdrop[1 + nl*2]>=0){ (*mc_index_sort_)[nc] = out_xdrop[1 + nl*2]; nc++;} nl++;}
  if (verbose>1){ raprintf(*mr_index_sort_,"int",1,(*nrows_)," %% mr_index_sort: ");}
  if (verbose>1){ raprintf(*mc_index_sort_,"int",1,(*ncols_)," %% mc_index_sort: ");}
  if (verbose){ printf(" %% [finished xcc_out_xdrop_sort]\n");}
}

void xcc_out_xdrop_load(char *fname,int *nlines_,int **out_xdrop_)
{
  /* recover out_xdrop from disc */
  int verbose=0;
  int nl=0; char c[128];
  FILE *fp=NULL; 
  if (verbose){ printf(" %% [entering xcc_out_xdrop_load]\n");}
  if (out_xdrop_!=NULL){
    if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s when reading from disc in xcc_out_xdrop_load.\n",fname); exit(RET_READ_FAIL);}
    (*nlines_)=0; while(!feof(fp)){ (*nlines_) += (fgetc(fp)=='\n');} if (verbose){ printf(" %% found %d lines in %s\n",(*nlines_),fname);}
    if (*out_xdrop_==NULL){ (*out_xdrop_) = (int *) wkspace_all0c(2*(*nlines_)*sizeof(int));}
    fseeko(fp,(off_t)(0),SEEK_SET);
    for (nl=0;nl<(*nlines_);nl++){ fscanf(fp,"%d",&((*out_xdrop_)[0 + nl*2])); fscanf(fp,"%d",&((*out_xdrop_)[1 + nl*2])); fscanf(fp,"%c",c);}
    fclose(fp);fp=NULL; 
    if (verbose>1){ raprintf(*out_xdrop_,"int",2,(*nlines_)," %% out_xdrop: ");}
    /* if (out_xdrop_!=NULL){ } */}
  if (verbose){ printf(" %% [finished xcc_out_xdrop_load]\n");}
}

void xcc_out_xdrop_ini2(int rpop_b_total,int *global_mr_a,int cpop_b,int *local_mc_a,int *nlines_,int **out_xdrop_)
{
  /* initalize out_xdrop in order */
  int verbose=2;
  int nl=0,nlr=0,nlc=0;
  if (verbose){ printf(" %% [entering xcc_out_xdrop_ini2]\n");}
  if (verbose>1){ printf(" %% rpop_b_total %d cpop_b %d\n",rpop_b_total,cpop_b);}
  if (verbose>2){
    raprintf(global_mr_a,"int",1,rpop_b_total," %% global_mr_a: ");
    raprintf(local_mc_a,"int",1,cpop_b," %% local_mc_a: ");
    /* if (verbose>1){ } */}
  if (out_xdrop_!=NULL){
    (*nlines_) = rpop_b_total + cpop_b; if (verbose){ printf(" %% setting %d lines\n",(*nlines_));}
    if (verbose>1){ printf(" %% allocating out_xdrop_\n");}
    if (*out_xdrop_==NULL){ (*out_xdrop_) = (int *) wkspace_all0c(2*(*nlines_)*sizeof(int));}
    if (verbose>1){ printf(" %% creating out_xdrop_\n");}
    nl=0; 
    if (verbose>1){ printf(" %% planting rows\n");}
    nlr=0;
    while (nl<(*nlines_) && nlr<rpop_b_total){
      (*out_xdrop_)[0 + nl*2] = global_mr_a[nlr];
      (*out_xdrop_)[1 + nl*2] = -1;
      nl++; nlr++; /* while rows */}
    if (verbose>1){ printf(" %% nl %d nlr %d\n",nl,nlr);}
    if (verbose>1){ printf(" %% planting cols\n");}
    nlc=0;
    while (nl<(*nlines_) && nlc<cpop_b){
      (*out_xdrop_)[0 + nl*2] = -1;
      (*out_xdrop_)[1 + nl*2] = local_mc_a[nlc];
      nl++; nlc++; /* while cols */}
    if (verbose>1){ printf(" %% nl %d nlc %d\n",nl,nlc);}
    if (verbose>2){ raprintf(*out_xdrop_,"int",2,(*nlines_)," %% out_xdrop: ");}
    /* if (out_xdrop_!=NULL){ } */}
  if (verbose){ printf(" %% [finished xcc_out_xdrop_ini2]\n");}
}

void xcc_out_xdrop_init(int rpop_b_total,int *global_mr_a,int cpop_b,int *local_mc_a,int *nlines_,int **out_xdrop_)
{
  /* initalize out_xdrop randomly */
  int verbose=0;
  int rij[rpop_b_total]; int cij[cpop_b]; int iij[rpop_b_total+cpop_b];
  int nl=0,nlr=0,nlc=0; 
  if (verbose){ printf(" %% [entering xcc_out_xdrop_init]\n");}
  if (verbose>1){ printf(" %% rpop_b_total %d cpop_b %d\n",rpop_b_total,cpop_b);}
  if (verbose>2){
    raprintf(global_mr_a,"int",1,rpop_b_total," %% global_mr_a: ");
    raprintf(local_mc_a,"int",1,cpop_b," %% local_mc_a: ");
    /* if (verbose>1){ } */}
  if (verbose>1){ printf(" %% setting rij\n");}
  for (nl=0;nl<rpop_b_total;nl++){ rij[nl] = global_mr_a[nl];} 
  if (verbose>1){ printf(" %% permuting rij\n");}
  irandperm(rpop_b_total,rij);
  if (verbose>1){ printf(" %% setting cij\n");}
  for (nl=0;nl<cpop_b;nl++){ cij[nl] = local_mc_a[nl];} 
  if (verbose>1){ printf(" %% permuting cij\n");}
  irandperm(cpop_b,cij);
  if (verbose>1){ printf(" %% setting iij 0\n");}
  for (nl=0;nl<rpop_b_total;nl++){ iij[nl] = 0;}
  if (verbose>1){ printf(" %% setting iij 1\n");}
  for (nl=rpop_b_total;nl<rpop_b_total+cpop_b;nl++){ iij[nl] = 1;}
  if (verbose>1){ printf(" %% permuting iij\n");}
  irandperm(rpop_b_total+cpop_b,iij);
  if (out_xdrop_!=NULL){
    (*nlines_) = rpop_b_total + cpop_b; if (verbose){ printf(" %% setting %d lines\n",(*nlines_));}
    if (verbose>1){ printf(" %% allocating out_xdrop_\n");}
    if (*out_xdrop_==NULL){ (*out_xdrop_) = (int *) wkspace_all0c(2*(*nlines_)*sizeof(int));}
    if (verbose>1){ printf(" %% creating out_xdrop_\n");}
    nl=0;
    while (nl<(*nlines_)){
      if (iij[nl]==0){
	(*out_xdrop_)[0 + nl*2] = rij[nlr];
	(*out_xdrop_)[1 + nl*2] = -1;
	nlr++; /* if (iij[nl]==0){ } */}
      if (iij[nl]==1){
	(*out_xdrop_)[0 + nl*2] = -1;
	(*out_xdrop_)[1 + nl*2] = cij[nlc];
	nlc++; /* if (iij[nl]==1){ } */}
      nl++; /* while (nl<(*nlines_)){ } */}
    if (verbose>1){ raprintf(*out_xdrop_,"int",2,(*nlines_)," %% out_xdrop: ");}
    /* if (out_xdrop_!=NULL){ } */}
  if (verbose){ printf(" %% [finished xcc_out_xdrop_init]\n");}
}

void xcc_scorebox_xdrop_array(int nxols,int x_max,int x_num,int **xval_,int **xdrop_)
{
  /* given maximum index x_max (out of nxols) and number of values x_num, ;
     we set up an array of index-values xval (from 0 to nxols-1) ;
     as well as an array of differences xdrop (number of indices to remove at each step). ;
     For example, if: ;
     nxols = 16, x_max = 7, x_num = 3 ;
     we set: ;
     xval = round([7,3.5,0]) = [7,4,0] ;
     xdrop = [8,3,4] ;
     Note that the first value in xdrop is the number of indices between (nxols-1) and xval(0). ;
     Subsequent values xdrop[j] = xval[j-1]-xval[j]. ;
   */
  int verbose=0;
  int nx=0,xdrop_initial=0;
  double d=0.0;
  if (verbose){ printf(" %% [entering xcc_scorebox_xdrop_array] nxols %d x_max %d x_num %d\n",nxols,x_max,x_num);}
  if (*xval_==NULL){ (*xval_) = (int *) wkspace_all0c(x_num*sizeof(int));}
  if (*xdrop_==NULL){ (*xdrop_) = (int *) wkspace_all0c(x_num*sizeof(int));}
  d = ((double)(x_max - 0))/((double)maximum(1,x_num-1));
  for (nx=0;nx<x_num;nx++){ (*xval_)[nx] = (int)round(x_max - nx*d);}
  (*xdrop_)[0] = maximum(0,nxols - (x_max+1));
  for (nx=1;nx<x_num;nx++){ (*xdrop_)[nx] = (*xval_)[nx-1]-(*xval_)[nx];}
  if (verbose>1){ raprintf((*xval_),"int",1,x_num," %% xval: ");}
  if (verbose>1){ raprintf((*xdrop_),"int",1,x_num," %% xdrop: ");}
  if (verbose){ printf(" %% [finished xcc_scorebox_xdrop_array]\n");}  
}
