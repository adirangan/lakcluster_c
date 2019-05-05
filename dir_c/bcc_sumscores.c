#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */


void xcc_out_trace_get(char *fname,int *nlines_,int **out_trace_ni_,int **out_trace_nr_,int **out_trace_nc_,double **out_trace_QR_,double **out_trace_QC_,int **out_trace_nb_)
{
  /* recover out_trace from disc */
  int verbose=0;
  int nl=0; char c[128];
  FILE *fp=NULL; 
  if (verbose){ printf(" %% [entering xcc_out_trace_get]\n");}
  if ((fp=fopen(fname,"r"))==NULL){ printf(" %% Warning! could not open %s when reading from disc.\n",fname); exit(RET_READ_FAIL);}
  (*nlines_)=0; while(!feof(fp)){ (*nlines_) += (fgetc(fp)=='\n');} if (verbose){ printf(" %% found %d lines in %s\n",(*nlines_),fname);}
  if (*out_trace_ni_==NULL){ (*out_trace_ni_) = (int *) wkspace_all0c(1*(*nlines_)*sizeof(int));}
  if (*out_trace_nr_==NULL){ (*out_trace_nr_) = (int *) wkspace_all0c(1*(*nlines_)*sizeof(int));}
  if (*out_trace_nc_==NULL){ (*out_trace_nc_) = (int *) wkspace_all0c(1*(*nlines_)*sizeof(int));}
  if (*out_trace_QR_==NULL){ (*out_trace_QR_) = (double *) wkspace_all0c(1*(*nlines_)*sizeof(double));}
  if (*out_trace_QC_==NULL){ (*out_trace_QC_) = (double *) wkspace_all0c(1*(*nlines_)*sizeof(double));}
  if (*out_trace_nb_==NULL){ (*out_trace_nb_) = (int *) wkspace_all0c(1*(*nlines_)*sizeof(int));}
  fseeko(fp,(off_t)(0),SEEK_SET);
  for (nl=0;nl<(*nlines_);nl++){ 
    fscanf(fp,"%d",&((*out_trace_ni_)[nl])); fscanf(fp,"%d",&((*out_trace_nr_)[nl])); fscanf(fp,"%d",&((*out_trace_nc_)[nl])); 
    fscanf(fp,"%lf",&((*out_trace_QR_)[nl])); fscanf(fp,"%lf",&((*out_trace_QC_)[nl])); 
    fscanf(fp,"%d",&((*out_trace_nb_)[nl])); 
    fscanf(fp,"%c",c); /* for (nl=0;nl<(*nlines_);nl++){ } */}
  fclose(fp);fp=NULL; 
  if (verbose>1){ raprintf(*out_trace_ni_,"int",1,(*nlines_)," %% out_trace_ni: ");}
  if (verbose>1){ raprintf(*out_trace_nr_,"int",1,(*nlines_)," %% out_trace_nr: ");}
  if (verbose>1){ raprintf(*out_trace_nc_,"int",1,(*nlines_)," %% out_trace_nc: ");}
  if (verbose>1){ raprintf(*out_trace_QR_,"double",1,(*nlines_)," %% out_trace_QR: ");}
  if (verbose>1){ raprintf(*out_trace_QC_,"double",1,(*nlines_)," %% out_trace_QC: ");}
  if (verbose>1){ raprintf(*out_trace_nb_,"int",1,(*nlines_)," %% out_trace_nb: ");}
  if (verbose){ printf(" %% [finished xcc_out_trace_get]\n");}
}

void bcc_sumscores_mxB(struct bcc_ajdk *D)
{
  /* writes D->A_umc_j into D->A_bmc_j  */
  /* writes D->A_umc_j_rmv into D->A_bmc_j_rmv  */
  /* writes D->A_umc_j_rtn into D->A_bmc_j_rtn */
  /* writes E->A_umr_j into E->A_bmr_j  */
  /* writes E->A_umr_j_rmv into E->A_bmr_j_rmv  */
  /* writes E->A_umr_j_rtn into E->A_bmr_j_rtn */
  /* finally, calls bcc_lrup_mxset, 
     which copies all of the 
     A_bmr_j_rtn etc (i.e., retained) and 
     A_bmr_j_rmv etc (i.e., removed) bitmasks 
     into the temporary variables used for lrup 
     (e.g., E->M_an used in low-rank update) 
  */
  int verbose=0;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb1=0,nr=0,nc=0; struct bcc_single *E=NULL;
  if (verbose){ raprintf(D->A_umc_j,"char",1,D->A_ncols," %% D->A_umc_j: ");}
  if (verbose){ raprintf(D->A_umc_j_rmv,"char",1,D->A_ncols," %% D->A_umc_j_rmv: ");}
  if (verbose){ raprintf(D->A_umc_j_rtn,"char",1,D->A_ncols," %% D->A_umc_j_rtn: ");}
  fill_uchar_zero(D->A_bmc_j,bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */); for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j,D->A_umc_j,nc);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j: "); if (verbose){ bprintf(D->A_bmc_j,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  fill_uchar_zero(D->A_bmc_j_rmv,bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */); for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j_rmv,D->A_umc_j_rmv,nc);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rmv: "); if (verbose){ bprintf(D->A_bmc_j_rmv,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  fill_uchar_zero(D->A_bmc_j_rtn,bsize(D->A_ncols)/* rup(D->A_ncols+D->A_ncols_extend,POPLENGTH)/BIT8 */); for (nc=0;nc<D->A_ncols;nc++){ b_copy_u(D->A_bmc_j_rtn,D->A_umc_j_rtn,nc);}
  sprintf(D->tmpAnchar," %%%% D->A_bmc_j_rtn: "); if (verbose){ bprintf(D->A_bmc_j_rtn,D->bitj,1,D->A_ncols,D->tmpAnchar);}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j: ",nb1); if (verbose){ raprintf(E->A_umr_j,"char",1,E->A_nrows,D->tmpAnchar);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j_rmv: ",nb1); if (verbose){ raprintf(E->A_umr_j_rmv,"char",1,E->A_nrows,D->tmpAnchar);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_umr_j_rtn: ",nb1); if (verbose){ raprintf(E->A_umr_j_rtn,"char",1,E->A_nrows,D->tmpAnchar);}
    fill_uchar_zero(E->A_bmr_j,bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */); for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j,E->A_umr_j,nr);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j: ",nb1); if (verbose){ bprintf(E->A_bmr_j,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    fill_uchar_zero(E->A_bmr_j_rmv,bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */); for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j_rmv,E->A_umr_j_rmv,nr);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j_rmv: ",nb1); if (verbose){ bprintf(E->A_bmr_j_rmv,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    fill_uchar_zero(E->A_bmr_j_rtn,bsize(E->A_nrows)/* rup(E->A_nrows+E->A_nrows_extend,POPLENGTH)/BIT8 */); for (nr=0;nr<E->A_nrows;nr++){ b_copy_u(E->A_bmr_j_rtn,E->A_umr_j_rtn,nr);}
    sprintf(D->tmpAnchar," %%%% E_[%.2d]->A_bmr_j_rtn: ",nb1); if (verbose){ bprintf(E->A_bmr_j_rtn,D->bitj,1,E->A_nrows,D->tmpAnchar);}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  D->A_cpop_j_rmv = popcount_uchar_array(D->A_bmc_j_rmv,D->A_mc_length);
  D->A_cpop_j_rtn = popcount_uchar_array(D->A_bmc_j_rtn,D->A_mc_length);
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    E->A_rpop_j_rmv = popcount_uchar_array(E->A_bmr_j_rmv,E->A_mr_length);
    E->A_rpop_j_rtn = popcount_uchar_array(E->A_bmr_j_rtn,E->A_mr_length);
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  bcc_lrup_mxset(D);
  if (verbose){ printf(" %% [finished bcc_sumscores_mxB]\n");}
}

void bcc_sumscores_mxC(struct bcc_ajdk *D)
{
  /* writes data to out_trace ; 
     0: out_iteration ;
     1: A_rpop_j_total ; rows remaining
     2: A_cpop_j ; columns remaining
     3: average of QR_svalue ; average row-score
     4: average of QC_svalue ; average col-score
     5: nb_rem ; covariate categories remaining
  */
  int verbose=0;
  double tmp_mean=0;
  /* 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MA = nrows_rem; NA = ncols_rem; MZ = size(Z_orig,1); NZ = ncols_rem;
    tmp_R = sum(QR); R_d = log(MA) + log(NA) + log(MA) + log(NA);
    tmp_C = sum(QC); C_d = log(NA) + log(MA) + log(NA) + log(MA);
    out_trace(iteration,:) = [iteration , MA , NA , sign(tmp_R)*exp(log(abs(tmp_R))-R_d) , sign(tmp_C)*exp(log(abs(tmp_C))-C_d) sum(cov_bother_flag_AAAA)+sum(cov_bother_flag_AZZA)];
    for nr=1:length(rdrop);
    out_xdrop(out_xdrop_ij,:) = [A_n_rind_vals_lookup(rdrop(end+1-nr))-1 , -1];
    out_xdrop_ij = out_xdrop_ij+1;
    end;%for nr=1:length(rdrop);
    for nc=1:length(cdrop);
    out_xdrop(out_xdrop_ij,:) = [-1 , A_n_cind_vals_lookup(cdrop(end+1-nc))-1];
    out_xdrop_ij = out_xdrop_ij+1;
    end;%for nc=1:length(cdrop);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */
  D->out_trace[0 + D->out_iteration*D->out_trace_length] = D->out_iteration;
  D->out_trace[1 + D->out_iteration*D->out_trace_length] = D->A_rpop_j_total;
  D->out_trace[2 + D->out_iteration*D->out_trace_length] = D->A_cpop_j;
  ra_stats(D->QR_svalue,"double",D->A_rpop_j_total,NULL,NULL,&tmp_mean,NULL); if (!isfinite(tmp_mean)){ tmp_mean = 1.0;}
  D->out_trace[3 + D->out_iteration*D->out_trace_length] = tmp_mean; 
  ra_stats(D->QC_svalue,"double",D->A_cpop_j,NULL,NULL,&tmp_mean,NULL); if (!isfinite(tmp_mean)){ tmp_mean = 1.0;}
  D->out_trace[4 + D->out_iteration*D->out_trace_length] = tmp_mean;
  D->out_trace[5 + D->out_iteration*D->out_trace_length] = D->Irem;
  D->out_iteration++;
  if (verbose){ printf(" %% [finished bcc_sumscores_mxC]\n");}
}

void bcc_sumscores_dmp(struct bcc_ajdk *D)
{
  /* dump temporary variables used within bcc to disc */
  /* Note that out_trace holds:
     0: out_iteration ;
     1: A_rpop_j_total ; rows remaining
     2: A_cpop_j ; columns remaining
     3: average of QR_svalue ; average row-score
     4: average of QC_svalue ; average col-score
     5: nb_rem ; covariate categories remaining
  */
  int verbose=0;
  int nl=0; FILE *fp=NULL; char tempchar[FNAMESIZE];
  FILE *fp_a=NULL; char tempchar_a[FNAMESIZE];
  FILE *fp_b=NULL; char tempchar_b[FNAMESIZE];
  if (verbose){ printf(" %% [entering bcc_sumscores_dmp]\n");}
  sprintf(tempchar,"%s/out_trace.txt",GLOBAL_DIR_NAME); if ((fp=fopen(tempchar,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",tempchar); exit(RET_READ_FAIL);}
  for (nl=0;nl<D->out_iteration;nl++){
    fprintf(fp,"%d %d %d %0.16f %0.16f %d\n",(int)D->out_trace[0 + nl*D->out_trace_length],(int)D->out_trace[1 + nl*D->out_trace_length],(int)D->out_trace[2 + nl*D->out_trace_length],D->out_trace[3 + nl*D->out_trace_length],D->out_trace[4 + nl*D->out_trace_length],(int)D->out_trace[5 + nl*D->out_trace_length]);
    /* for (nl=0;nl<D->out_iteration;nl++){ } */}
  fclose(fp);fp=NULL; 
  sprintf(tempchar_a,"%s/out_xdrop_a.txt",GLOBAL_DIR_NAME); if ((fp_a=fopen(tempchar_a,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",tempchar_a); exit(RET_READ_FAIL);}
  sprintf(tempchar_b,"%s/out_xdrop_b.txt",GLOBAL_DIR_NAME); if ((fp_b=fopen(tempchar_b,"w"))==NULL){ printf(" %% Warning! could not open %s when writing to disc.\n",tempchar_b); exit(RET_READ_FAIL);}
  for (nl=0;nl<D->out_xdrop_ij;nl++){
    fprintf(fp_a,"%d %d\n",(int)D->out_xdrop_a[0 + nl*2],(int)D->out_xdrop_a[1 + nl*2]);
    fprintf(fp_b,"%d %d\n",(int)D->out_xdrop_b[0 + nl*2],(int)D->out_xdrop_b[1 + nl*2]);
    /* for (nl=0;nl<D->out_xdrop_ij;nl++){ } */}
  fclose(fp_a);fp_a=NULL; 
  fclose(fp_b);fp_b=NULL; 
  if (verbose){ printf(" %% [finished bcc_sumscores_dmp]\n");}
}

void bcc_sumscores_mxA(struct bcc_ajdk *D)
{
  /* Obtains rdrop and cdrop from get_xdrop (based on GLOBAL_gamma). ;
     Copies D->A_bmc_j to D->A_umc_j. ; 
     Defines D->A_umc_j_rmv and D->A_umc_j_rtn ; removing first cdrop entries of D->QC_index_local_mc_a. ;
     Copies E->A_bmr_j to E->A_umr_j. ; 
     Defines E->A_umr_j_rmv and E->A_umr_j_rtn ; removing first rdrop entries of D->QR_index_local_mr_a ; using D->QR_index_local_nb to index nb. ;
     Call bcc_sumscores_mxB ; 
      Copies all instances of umc and umr to associated bmc and bmr respectively ;
      Note that bcc_sumscores_mxB also calls 
      bcc_lrup_mxset:
       Copies all of the 
       A_bmr_j_rtn etc (i.e., retained) and 
       A_bmr_j_rmv etc (i.e., removed) bitmasks 
       into the temporary variables used for lrup 
       (e.g., E->M_an used in low-rank update).
     Call bcc_sumscores_mxC:
      Writing data to out_trace ;
     Finally: add to D->out_xdrop_a and D->out_xdrop_b (using D->QR_index_global_mr_x for rows, and D->QC_index_local_mc for cols). ; 
     Note: D->out_xdrop_x lists row-indices first, then col-indices ;
   */
  int verbose=0; 
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_single *E=NULL;
  int rdrop=0;int cdrop=0;
  int nb1=0,nr=0,nc=0;
  if (verbose){ printf(" %% [entering bcc_sumscores_mxA]\n");}
  get_xdrop(D->A_rpop_j_total,D->A_cpop_j,&rdrop,&cdrop);
  for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j[nc] = bget__on(D->A_bmc_j,nc);} for (nc=0;nc<D->A_ncols;nc++){ D->A_umc_j_rtn[nc] = D->A_umc_j[nc]; D->A_umc_j_rmv[nc] = 0;}
  for (nc=0;nc<minimum(cdrop,D->A_cpop_j);nc++){ D->A_umc_j_rmv[D->QC_index_local_mc_a[nc]]=1; D->A_umc_j_rtn[D->QC_index_local_mc_a[nc]]=0;}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];
    for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j[nr] = bget__on(E->A_bmr_j,nr);} for (nr=0;nr<E->A_nrows;nr++){ E->A_umr_j_rtn[nr] = E->A_umr_j[nr]; E->A_umr_j_rmv[nr] = 0;}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  for (nr=0;nr<minimum(rdrop,D->A_rpop_j_total);nr++){ E_[D->QR_index_local_nb[nr]]->A_umr_j_rmv[D->QR_index_local_mr_a[nr]]=1; E_[D->QR_index_local_nb[nr]]->A_umr_j_rtn[D->QR_index_local_mr_a[nr]]=0;}
  bcc_sumscores_mxB(D); bcc_sumscores_mxC(D);
  if (verbose){ raprintf(D->out_trace,"double_trn",D->out_trace_length,D->out_iteration," %% D->out_trace: ");}
  for (nr=0;nr<minimum(rdrop,D->A_rpop_j_total);nr++){ 
    D->out_xdrop_a[0 + D->out_xdrop_ij*2] = D->QR_index_global_mr_a[nr]; D->out_xdrop_a[1 + D->out_xdrop_ij*2] = -1;
    D->out_xdrop_b[0 + D->out_xdrop_ij*2] = D->QR_index_global_mr_b[nr]; D->out_xdrop_b[1 + D->out_xdrop_ij*2] = -1;
    D->out_xdrop_ij++; /* for (nr=0;nr<minimum(rdrop,D->A_rpop_j_total);nr++){ } */}
  for (nc=0;nc<minimum(cdrop,D->A_cpop_j);nc++){ 
    D->out_xdrop_a[1 + D->out_xdrop_ij*2] = D->QC_index_local_mc_a[nc]; D->out_xdrop_a[0 + D->out_xdrop_ij*2] = -1; 
    D->out_xdrop_b[1 + D->out_xdrop_ij*2] = D->QC_index_local_mc_b[nc]; D->out_xdrop_b[0 + D->out_xdrop_ij*2] = -1; 
    D->out_xdrop_ij++; /* for (nc=0;nc<minimum(cdrop,D->A_cpop_j);nc++){ } */}
  if (verbose){ raprintf(D->out_xdrop_a,"int",2,D->out_xdrop_ij," %% D->out_xdrop_a: ");}
  if (verbose){ raprintf(D->out_xdrop_b,"int",2,D->out_xdrop_ij," %% D->out_xdrop_b: ");}
  if (verbose){ printf(" %% [finished bcc_sumscores_mxA]\n");}
}

void xcc_get_Ireq(char *strategy,int nbins,int Irem,int *Ireq_p,int *Irow_p,int *Icol_p)
{
  int Ireq=0,Irow=0,Icol=0;
  if (strstr(strategy,"empty high")){ /* If empty bins are sorted high */
    Ireq = maximum(0,minimum(nbins,minimum(Irem,GLOBAL_Ireq)));
    if (Ireq<=0){ Irow=0; Icol=0;} 
    if (Ireq>0){ Irow = Irem - Ireq; Icol = Irem*Irem - Ireq*Ireq;}
    /* if (strstr(strategy,"empty high")){ } */}
  if (strstr(strategy,"empty low")){ /* If empty bins are sorted low */
    Ireq = maximum(0,minimum(nbins,GLOBAL_Ireq));
    if (Ireq<=0){ Irow=0; Icol=0;} 
    if (Ireq>0){ Irow = nbins - Ireq; Icol = nbins*nbins - Ireq*Ireq;}
    /* if (strstr(strategy,"empty low")){ } */}
  if (Ireq_p!=NULL){ *Ireq_p = Ireq;}
  if (Irow_p!=NULL){ *Irow_p = Irow;}
  if (Icol_p!=NULL){ *Icol_p = Icol;}
}

void bcc_sumscores_xij(struct bcc_ajdk *D)
{
  /* Sets D->QC_index_local_mc_x.
     These arrays store the indices (i.e., index locations) of the masks 
     M_An->mc_x==D->A_bmc_x.
     Then sorts col-scores in increasing order across columns.
     At this point the first few entries of D->QC_index_local_mc_x correspond to the columns that are least likely to participate in a bicluster. ;
     Sets D->QR_index_local_mr_x and D->QR_index_local_nb.
     These arrays store the indices (i.e., index locations) of the masks 
     M_An[nb]->mr==E->A_bmr_x.
     Also copies D->QR_index_global_mr_x from E->QR_index_globalL_mr_x.
     These arrays store the overall index locations of the masks M_An->mr accumulated across bins. ;
     Then sorts row-scores in increasing order across rows.
     At this point the first few entries of D->QR_index_global_mr_x correspond to the rows that are least likely to participate in a bicluster. ;
     Also redefines D->A_rpop_j_total using E->A_bmr_b and E->A_bmr_j. 
  */
  int verbose=0; 
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb1=0; struct bcc_single *E=NULL;
  int na_a=0,na_b=0,na_j=0,ma_a=0,ma_b=0,ma_j=0;
  int Irow=0,Icol=0;
  unsigned int nn=0;
  if (verbose){ printf(" %% [entering bcc_sumscores_xij]\n");}
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
	  D->QC_index_local_mc_a[na_j] = na_a; D->QC_index_local_mc_b[na_j] = na_b; D->QC_index_local_mc_j[na_j] = na_j;
	  na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	na_b++; /* if (bget__on(D->A_bmc_b,na_a)){ } */}
      na_a++;/* while (na_a<D->A_ncols){ } */}
    if (verbose>1){
      raprintf(  D->QC_svalue,"double",1,D->A_cpop_j," %%   D->QC_svalue pre :");
      raprintf(D->QC_index_local_mc_a,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_a pre :");
      raprintf(D->QC_index_local_mc_b,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_b pre :");
      raprintf(D->QC_index_local_mc_j,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_j pre :");
      /* if (verbose){ } */}
    nn = dQuickSort_xij(0,D->QC_svalue,1,D->QC_index_local_mc_a,D->QC_index_local_mc_b,D->QC_index_local_mc_j,NULL,NULL,NULL,0,D->A_cpop_j-1);
    if (verbose>1){ printf(" %% finished sorting QC_svalue, maximum recursion_level %d\n",nn);}
    if (verbose>1){
      raprintf(  D->QC_svalue,"double",1,D->A_cpop_j," %%   D->QC_svalue pos :");
      raprintf(D->QC_index_local_mc_a,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_a pos :");
      raprintf(D->QC_index_local_mc_b,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_b pos :");
      raprintf(D->QC_index_local_mc_j,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_j pos :");
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
	    D->QR_index_local_mr_a[D->A_rpop_j_total] = ma_a; D->QR_index_local_mr_b[D->A_rpop_j_total] = ma_b; D->QR_index_local_mr_j[D->A_rpop_j_total] = ma_j; D->QR_index_local_nb[D->A_rpop_j_total] = nb1; 
	    D->QR_index_global_mr_a[D->A_rpop_j_total] = E->QR_index_global_mr_a[ma_a]; D->QR_index_global_mr_b[D->A_rpop_j_total] = E->QR_index_global_mr_b[ma_a];
	    ma_j++; D->A_rpop_j_total++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	  ma_b++; D->A_rpop_b_total++; /* while (ma_a<E->A_nrows){ } */}
	ma_a++; D->A_nrows_total++; /* while (ma_a<E->A_nrows){ } */}
      /* for (nb1=0;nb1<nbins;nb1++){ } */}
    if (verbose>1){
      raprintf(  D->QR_svalue,"double",1,D->A_rpop_j_total," %%   D->QR_svalue pre :");
      raprintf(D->QR_index_local_mr_a,   "int",1,D->A_rpop_j_total," %% D->QR_index_local_mr_a pre :");
      raprintf(D->QR_index_local_mr_b,   "int",1,D->A_rpop_j_total," %% D->QR_index_local_mr_b pre :");
      raprintf(D->QR_index_local_mr_j,   "int",1,D->A_rpop_j_total," %% D->QR_index_local_mr_j pre :");
      raprintf( D->QR_index_local_nb,   "int",1,D->A_rpop_j_total," %%  D->QR_index_local_nb pre :");
      raprintf( D->QR_index_global_mr_a,   "int",1,D->A_rpop_j_total," %%  D->QR_index_global_mr_a pre :");
      raprintf( D->QR_index_global_mr_b,   "int",1,D->A_rpop_j_total," %%  D->QR_index_global_mr_b pre :");
      /* if (verbose){ } */}
    nn = dQuickSort_xij(0,D->QR_svalue,1,D->QR_index_local_mr_a,D->QR_index_local_mr_b,D->QR_index_local_mr_j,D->QR_index_local_nb,D->QR_index_global_mr_a,D->QR_index_global_mr_b,0,D->A_rpop_j_total-1);
    if (verbose>1){ printf(" %% finished sorting QR_svalue, maximum recursion_level %d\n",nn);}    
    if (verbose>1){
      raprintf(  D->QR_svalue,"double",1,D->A_rpop_j_total," %%   D->QR_svalue pos :");
      raprintf(D->QR_index_local_mr_a,   "int",1,D->A_rpop_j_total," %% D->QR_index_local_mr_a pos :");
      raprintf(D->QR_index_local_mr_b,   "int",1,D->A_rpop_j_total," %% D->QR_index_local_mr_b pos :");
      raprintf(D->QR_index_local_mr_j,   "int",1,D->A_rpop_j_total," %% D->QR_index_local_mr_j pos :");
      raprintf( D->QR_index_local_nb,   "int",1,D->A_rpop_j_total," %%  D->QR_index_local_nb pos :");
      raprintf( D->QR_index_global_mr_a,   "int",1,D->A_rpop_j_total," %%  D->QR_index_global_mr_a pos :");
      raprintf( D->QR_index_global_mr_b,   "int",1,D->A_rpop_j_total," %%  D->QR_index_global_mr_b pos :");
      /* if (verbose){ } */}
    /* if (D->A_rpop_b_total && D->A_cbother){ } */}
  if (verbose){ printf(" %% [finished bcc_sumscores_xij]\n");}
}

void bcc_sumscores_cmb(struct bcc_ajdk *D)
{
  /* Combines scores such as AtTAnAtTAn and AtTYnYtTAn to form final score. ;
   */
  int verbose=0;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb1=0,nb2=0,nbx=0; struct bcc_single *E=NULL; unsigned char bmc1=255;
  if (verbose){ printf(" %% [entering bcc_sumscores_cmb]\n");}
  /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    QC = QC_AAAA_min - QC_AAZZ_min;
    for nb1=0:nbins-1;
    QR{1+nb1} = QR_AAAA{1+nb1} - QR_AZZA{1+nb1};
    end;%for nb1=0:nbins-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */
  if (verbose){ printf(" %% [entering bcc_sumscores_cmb]\n");}
  if (D->A_cbother && D->A_rpop_b_total){
    if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->Z_rpop_b_total){ for (nbx=0;nbx<nbins*nbins;nbx++){ dra_plustimesequals_s___m_m(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,1/*D->T_ncols*/,&(D->QC_AtTYnWtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),+1.0,D->A_bmc_j,&bmc1/*D->T_bmc_j*/);}}
    if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ for (nbx=0;nbx<nbins*nbins;nbx++){ dra_plustimesequals_s___m_m(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,1/*D->T_ncols*/,&(D->QC_AtTAnZtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),-1.0,D->A_bmc_j,&bmc1/*D->T_bmc_j*/);}}
    if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->A_rpop_b_total){ for (nbx=0;nbx<nbins*nbins;nbx++){ dra_plustimesequals_s___m_m(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,1/*D->T_ncols*/,&(D->QC_AtTYnYtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),-1.0,D->A_bmc_j,&bmc1/*D->T_bmc_j*/);}}
    if (verbose){ for (nbx=0;nbx<nbins*nbins;nbx++){ sprintf(D->tmpAnchar," %%%% D->QC_AtTAnAtTAn_nrm[%03d]: ",nbx); raprintf(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,1,D->tmpAnchar);}}
    /* if (D->A_cbother && D->A_rpop_b_total){ } */}
  for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; 
    if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){
      if (E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother){ for (nb2=0;nb2<nbins;nb2++){ dra_plustimesequals_s___m_m(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,1/*D->T_ncols*/,&(E->QR_AnZtSWnYt_nrm[nb2*E->A_nrows*D->T_ncols]),+1.0,E->A_bmr_j,&bmc1/*D->T_bmc_j*/);}}
      if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ for (nb2=0;nb2<nbins;nb2++){ dra_plustimesequals_s___m_m(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,1/*D->T_ncols*/,&(E->QR_AnZtSZnAt_nrm[nb2*E->A_nrows*D->T_ncols]),-1.0,E->A_bmr_j,&bmc1/*D->T_bmc_j*/);}}
      if (E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother){ for (nb2=0;nb2<nbins;nb2++){ dra_plustimesequals_s___m_m(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,1/*D->T_ncols*/,&(E->QR_AnAtTYnYt_nrm[nb2*E->A_nrows*D->T_ncols]),-1.0,E->A_bmr_j,&bmc1/*D->T_bmc_j*/);}}
      if (verbose){ for (nb2=0;nb2<nbins;nb2++){ sprintf(D->tmpAnchar," %%%% E_[%d]->QR_AnAtTAnAt_nrm[%d]: ",nb1,nb2); raprintf(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,1,D->tmpAnchar);}}
      /* if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ } */}
    /* for (nb1=0;nb1<nbins;nb1++){ } */}
  if (verbose){ printf(" %% [finished bcc_sumscores_cmb]\n");}
}

void bcc_sumscores_srt(struct bcc_ajdk *D)
{
  /* If GLOBAL_Ireq>0, then we sort the scores by bin. ;
     Otherwise, if GLOBAL_Ireq<=0, we use the summed scores over all bins (already stored in nb==0).
  */
  int verbose=0; 
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb1=0,nb2=0,nbx=0; struct bcc_single *E=NULL; 
  /* unsigned char bmc1=255; */
  int ma_j=0,ma_b=0,ma_a=0,na_j=0,na_b=0,na_a=0;
  if (verbose){ printf(" %% [entering bcc_sumscores_srt]\n");}
  /* 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     QC_AAAA_nrm = min(QC_AAAA_nrm,2);
     QC_AAZZ_nrm = min(QC_AAZZ_nrm,2);
     for nb1=0:nbins-1;
     QR_AAAA_nrm{1+nb1} = sort(QR_AAAA_nrm{1+nb1},2);
     QR_AZZA_nrm{1+nb1} = sort(QR_AZZA_nrm{1+nb1},2);
     end;%for nb1=0:nbins-1;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */
  if (GLOBAL_Ireq>0){
    for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];     
      if (verbose){ printf(" %% nb1 %d\n",nb1);}
      if (E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother){ 
	ma_a=0; ma_b=0; ma_j=0; 
	while (ma_a<E->A_nrows){
	  if (bget__on(E->A_bmr_b,ma_a)){
	    if (bget__on(E->A_bmr_j,ma_a)){
	      if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnZtSWnYt_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	      if (GLOBAL_Ireq>0){ dQuickSort(0,&(E->QR_AnZtSWnYt_nrm[ma_a+0*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,0,nbins-1);}
	      if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnZtSWnYt_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	      ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	    ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	  ma_a++; /* while (ma_a<E->A_nrows){ } */}
	/* if (E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother){ } */}
      if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ 
	ma_a=0; ma_b=0; ma_j=0; 
	while (ma_a<E->A_nrows){
	  if (bget__on(E->A_bmr_b,ma_a)){
	    if (bget__on(E->A_bmr_j,ma_a)){
	      if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnZtSZnAt_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	      if (GLOBAL_Ireq>0){ dQuickSort(0,&(E->QR_AnZtSZnAt_nrm[ma_a+0*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,0,nbins-1);}
	      if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnZtSZnAt_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	      ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	    ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	  ma_a++; /* while (ma_a<E->A_nrows){ } */}
	/* if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ } */}
      if (E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother){ 
	ma_a=0; ma_b=0; ma_j=0; 
	while (ma_a<E->A_nrows){
	  if (bget__on(E->A_bmr_b,ma_a)){
	    if (bget__on(E->A_bmr_j,ma_a)){
	      if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnAtTYnYt_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	      if (GLOBAL_Ireq>0){ dQuickSort(0,&(E->QR_AnAtTYnYt_nrm[ma_a+0*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,0,nbins-1);}
	      if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnAtTYnYt_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	      ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	    ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	  ma_a++; /* while (ma_a<E->A_nrows){ } */}
	/* if (E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother){ } */}
      if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ 
	ma_a=0; ma_b=0; ma_j=0; 
	while (ma_a<E->A_nrows){
	  if (bget__on(E->A_bmr_b,ma_a)){
	    if (bget__on(E->A_bmr_j,ma_a)){
	      if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnAtTAnAt_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	      if (GLOBAL_Ireq>0){ dQuickSort(0,&(E->QR_AnAtTAnAt_nrm[ma_a+0*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,0,nbins-1);}
	      if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnAtTAnAt_nrm[ma_a+0*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
	      ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
	    ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	  ma_a++; /* while (ma_a<E->A_nrows){ } */}
	/* if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ } */}
      /* for (nb1=0;nb1<nbins;nb1++){ } */}
    if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->Z_rpop_b_total){ 
      na_a=0; na_b=0; na_j=0; 
      while (na_a<D->A_ncols){
	if (bget__on(D->A_bmc_b,na_a)){
	  if (bget__on(D->A_bmc_j,na_a)){
	    if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTYnWtSZn_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	    if (GLOBAL_Ireq>0){ dQuickSort(0,&(D->QC_AtTYnWtSZn_nrm[na_a+0*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,0,nbins*nbins-1);}
	    if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTYnWtSZn_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	    na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	  na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
	na_a++; /* while (na_a<D->A_ncols){ } */}
      /* if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->Z_rpop_b_total){ } */}
    if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ 
      na_a=0; na_b=0; na_j=0; 
      while (na_a<D->A_ncols){
	if (bget__on(D->A_bmc_b,na_a)){
	  if (bget__on(D->A_bmc_j,na_a)){
	    if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTAnZtSZn_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	    if (GLOBAL_Ireq>0){ dQuickSort(0,&(D->QC_AtTAnZtSZn_nrm[na_a+0*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,0,nbins*nbins-1);}
	    if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTAnZtSZn_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	    na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	  na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
	na_a++; /* while (na_a<D->A_ncols){ } */}
      /* if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ } */}
    if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->A_rpop_b_total){ 
      na_a=0; na_b=0; na_j=0; 
      while (na_a<D->A_ncols){
	if (bget__on(D->A_bmc_b,na_a)){
	  if (bget__on(D->A_bmc_j,na_a)){
	    if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTYnYtTAn_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	    if (GLOBAL_Ireq>0){ dQuickSort(0,&(D->QC_AtTYnYtTAn_nrm[na_a+0*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,0,nbins*nbins-1);}
	    if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTYnYtTAn_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	    na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	  na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
	na_a++; /* while (na_a<D->A_ncols){ } */}
      /* if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->A_rpop_b_total){ } */}
    if (D->A_cbother && D->A_rpop_b_total){ 
      na_a=0; na_b=0; na_j=0; 
      while (na_a<D->A_ncols){
	if (bget__on(D->A_bmc_b,na_a)){
	  if (bget__on(D->A_bmc_j,na_a)){
	    if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTAnAtTAn_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	    if (GLOBAL_Ireq>0){ dQuickSort(0,&(D->QC_AtTAnAtTAn_nrm[na_a+0*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,0,nbins*nbins-1);}
	    if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTAnAtTAn_nrm[na_a+0*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
	    na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	  na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
	na_a++; /* while (na_a<D->A_ncols){ } */}
      /* if (D->A_cbother && D->A_rpop_b_total){ } */}
    /* if (GLOBAL_Ireq>0){ } */}
  if (verbose){ printf(" %% [finished bcc_sumscores_srt]\n");}
}

void bcc_sumscores_ifT(struct bcc_ajdk *D)
{
  /* If GLOBAL_Ireq>0, we correct each E_[nb1]->QR_AnAtTAnAt_nrm[0+nb2*A_nrows*T_ncols] for T, across all values of nb1 and nb2. ;
     Otherwise, if GLOBAL_Ireq<=0 also correct each E_[nb1]->QR_AnAtTAnAt_nrm[0+nb2*A_nrows*T_ncols] for T, but only use the nb2==0 term for each nb1.
  */
  int verbose=0; 
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb1=0,nb2=0,nbx=0; struct bcc_single *E=NULL; 
  int n_mds = maximum(0,D->T_cpop_j-1); int n_mds_max = 6; double mds_scale_factor=1.0;
  if (verbose){ printf(" %% [entering bcc_sumscores_ifT], n_mds %d\n",n_mds);}
  /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_mds = 5; mds_scale_factor = [1.00,0.30,0.12,0.10,0.09];
    QC_all = temp_D_At_T_AnAt_T_An{1+0}/((A_nrows_(1+0) - 1)*A_ncols*(A_ncols-1));
    QR_all = temp_AnAt_T_AnAt{1+0}/(A_nrows_(1+0)*(A_nrows_(1+0) - 1)*A_ncols);
    QC = QC_all(:,1) - sum(QC_all(:,2:end),2)/n_mds/mds_scale_factor(n_mds);
    QR = (QR_all(:,1).^2 - sum(QR_all(:,2:end).^2,2)/n_mds/mds_scale_factor(n_mds));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */
  if (n_mds>n_mds_max){ printf(" %% Warning! n_mds>n_mds_max %d\n",n_mds_max); n_mds = n_mds_max;} 
  if (n_mds<=0){
    if (verbose){ printf(" %% n_mds==0, do nothing\n");}
    /* if (n_mds<=0){ } */}
  else /* if (n_mds>0) */{
    mds_scale_factor = GLOBAL_kappa_squared_loop_scale_factor_[n_mds-1];
    if (GLOBAL_kappa_squared>0);{ mds_scale_factor = GLOBAL_kappa_squared;}
    if (verbose){ printf(" %% n_mds==%d, mds_scale_factor=%0.3f\n",n_mds,mds_scale_factor);}
    for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1]; 
      for (nb2=0;nb2<(GLOBAL_Ireq>0 ? nbins : 1);nb2++){ nbx = nb1+nb2*nbins; 
	if (verbose){ printf(" %% nb1 %d nb2 %d nbx %d\n",nb1,nb2,nbx);}
	if (E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother){ 
	  dra_mds_pow_s___m_m(&(E->QR_AnZtSWnYt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,n_mds*mds_scale_factor,E->A_bmr_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(E->QR_AnZtSWnYt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnZtSWnYt_nrm: ");}
	  /* if (E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother){ } */}
	if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ 
	  dra_mds_pow_s___m_m(&(E->QR_AnZtSZnAt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,n_mds*mds_scale_factor,E->A_bmr_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(E->QR_AnZtSZnAt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnZtSZnAt_nrm: ");}
	  /* if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ } */}
	if (E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother){ 
	  dra_mds_pow_s___m_m(&(E->QR_AnAtTYnYt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,n_mds*mds_scale_factor,E->A_bmr_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(E->QR_AnAtTYnYt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnAtTYnYt_nrm: ");}
	  /* if (E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother){ } */}
	if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ 
	  dra_mds_pow_s___m_m(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,n_mds*mds_scale_factor,E->A_bmr_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnAtTAnAt_nrm: ");}
	  /* if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ } */}
	if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->Z_rpop_b_total){ 
	  dra_mds_nrm_s___m_m(&(D->QC_AtTYnWtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,n_mds*mds_scale_factor,D->A_bmc_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(D->QC_AtTYnWtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTYnWtSZn_nrm: ");}
	  /* if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->Z_rpop_b_total){ } */}
	if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ 
	  dra_mds_nrm_s___m_m(&(D->QC_AtTAnZtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,n_mds*mds_scale_factor,D->A_bmc_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(D->QC_AtTAnZtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTAnZtSZn_nrm: ");}
	  /* if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ } */}
	if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->A_rpop_b_total){ 
	  dra_mds_nrm_s___m_m(&(D->QC_AtTYnYtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,n_mds*mds_scale_factor,D->A_bmc_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(D->QC_AtTYnYtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTYnYtTAn_nrm: ");}
	  /* if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->A_rpop_b_total){ } */}
	if (D->A_cbother && D->A_rpop_b_total){ 
	  dra_mds_nrm_s___m_m(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,n_mds*mds_scale_factor,D->A_bmc_j,D->T_bmc_j); 
	  if (verbose){ raprintf(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTAnAtTAn_nrm: ");}
	  /* if (D->A_cbother && D->A_rpop_b_total){ } */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<(GLOBAL_Ireq>0 ? nbins : 1);nb2++){ }} */}}    
    /* if (n_mds>0){ } */}
  if (verbose){ printf(" %% [finished bcc_sumscores_ifT]\n");}
}

void bcc_sumscores_nrm(struct bcc_ajdk *D)
{
  /* If GLOBAL_Ireq>0, we normalize by the number of summands in each bin. ;
     Otherwise, if GLOBAL_Ireq<=0 we normalize by the number of summands across all bins, and then sum the results. ;
  */
  int verbose=0; 
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_double **F_ = D->F_;
  int nb1=0,nb2=0,nbx=0; long long int lld=0; struct bcc_single *E=NULL; struct bcc_double *F=NULL;
  int ns_a=0;
  int ma_a=0,ma_b=0,ma_j=0,na_a=0,na_b=0,na_j=0;
  if (verbose){ printf(" %% [entering bcc_sumscores_nrm]\n");}
  /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nb1=0:nbins-1;for nb2=0:nbins-1;nb_tab = nb1+nb2*nbins;
    QR_AAAA{1+nb_tab} = temp_AnAt_T_AnAt{1+nb_tab}/((A_nrows_(1+nb2) - 1*(nb1==nb2))*A_ncols*(A_ncols-1));
    QR_AZZA{1+nb_tab} = temp_AnZt_S_ZnAt{1+nb_tab}/((Z_nrows_(1+nb2) - 0*(nb1==nb2))*A_ncols*(A_ncols-1));
    QC_AAAA{1+nb_tab} = temp_D_At_T_AnAt_T_An{1+nb_tab}/(A_nrows_(1+nb1)*(A_nrows_(1+nb2) - 1*(nb1==nb2))*A_ncols);
    QC_AAZZ{1+nb_tab} = temp_D_At_T_AnZt_S_Zn{1+nb_tab}/(A_nrows_(1+nb1)*(Z_nrows_(1+nb2) - 0*(nb1==nb2))*A_ncols);
    end;end;%for nb1=0:nbins-1;for nb2=0:nbins-1;nb_tab = nb1+nb2*nbins;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    QC_AAAA_nrm = zeros(A_ncols,nbins*nbins);
    QC_AAZZ_nrm = zeros(A_ncols,nbins*nbins);
    for nb1=0:nbins-1;
    QR_AAAA_nrm{1+nb1} = zeros(A_nrows_(1+nb1),nbins);
    QR_AZZA_nrm{1+nb1} = zeros(A_nrows_(1+nb1),nbins);
    for nb2=0:nbins-1;nb_tab = nb1+nb2*nbins;
    QR_AAAA_nrm{1+nb1}(:,1+nb2) = QR_AAAA{1+nb_tab};
    QR_AZZA_nrm{1+nb1}(:,1+nb2) = QR_AZZA{1+nb_tab};
    QC_AAAA_nrm(:,1+nb_tab) = QC_AAAA{1+nb_tab};
    QC_AAZZ_nrm(:,1+nb_tab) = QC_AAZZ{1+nb_tab};
    end;end;%for nb1=0:nbins-1;for nb2=0:nbins-1;nb_tab = nb1+nb2*nbins;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  */
  if (GLOBAL_Ireq>0){
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins;
	if (verbose){ printf(" %% nb1 %d nb2 %d nbx %d\n",nb1,nb2,nbx);}
	E = E_[nb1]; F = F_[nbx];
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->Y_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_AnZtSWnYt_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double)); 
	  lld = (long long int)(E_[nb2]->M_Zn->rpop_j - 0*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_Zn->cpop_j) * (long long int)(E_[nb2]->M_Wn->cpop_j);
	  dra_plusdivequals_s___m_m(&(E->QR_AnZtSWnYt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_AnZtSWnYt->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_AnZtSWnYt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnZtSWnYt_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->Y_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->A_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_AnZtSZnAt_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb2]->M_Zn->rpop_j - 0*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_Zn->cpop_j) * (long long int)(E_[nb2]->M_Zn->cpop_j - 1);
	  dra_plusdivequals_s___m_m(&(E->QR_AnZtSZnAt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_AnZtSZnAt->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_AnZtSZnAt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnZtSZnAt_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->A_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->Y_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_AnAtTYnYt_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb2]->M_An->rpop_j - 1*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_An->cpop_j) * (long long int)(E_[nb2]->M_Yn->cpop_j);
	  dra_plusdivequals_s___m_m(&(E->QR_AnAtTYnYt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_AnAtTYnYt->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_AnAtTYnYt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnAtTYnYt_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->Y_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->A_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb2]->M_An->rpop_j - 1*(nb1==nb2?1:0)) * (long long int)(E_[nb2]->M_An->cpop_j) * (long long int)(E_[nb2]->M_An->cpop_j - 1);
	  dra_plusdivequals_s___m_m(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_AnAtTAnAt->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnAtTAnAt_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->A_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->Z_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_AtTYnWtSZn_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb1]->M_Yn->rpop_j) * (long long int)(E_[nb1]->M_Yn->cpop_j - 0) * (long long int)(E_[nb2]->M_Wn->rpop_j);
	  dra_plusdivequals_s___m_m(&(D->QC_AtTYnWtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_AtTYnWtSZn->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_AtTYnWtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTYnWtSZn_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->Z_rbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->Z_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_AtTAnZtSZn_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb1]->M_An->rpop_j) * (long long int)(E_[nb1]->M_An->cpop_j - 1) * (long long int)(E_[nb2]->M_Zn->rpop_j);
	  dra_plusdivequals_s___m_m(&(D->QC_AtTAnZtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_AtTAnZtSZn->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_AtTAnZtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTAnZtSZn_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->Z_rbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->A_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_AtTYnYtTAn_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb1]->M_Yn->rpop_j) * (long long int)(E_[nb1]->M_Yn->cpop_j - 0) * (long long int)(E_[nb2]->M_Yn->rpop_j - 1*(nb1==nb2?1:0));
	  dra_plusdivequals_s___m_m(&(D->QC_AtTYnYtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_AtTYnYtTAn->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_AtTYnYtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTYnYtTAn_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->A_rbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->A_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(E_[nb1]->M_An->rpop_j) * (long long int)(E_[nb1]->M_An->cpop_j - 1) * (long long int)(E_[nb2]->M_An->rpop_j - 1*(nb1==nb2?1:0));
	  dra_plusdivequals_s___m_m(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_AtTAnAtTAn->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTAnAtTAn_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->A_rbother){ } */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    /* if (GLOBAL_Ireq>0){ } */}
  if (GLOBAL_Ireq<=0){
    for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ nbx = nb1+nb2*nbins;
	if (verbose){ printf(" %% nb1 %d nb2 %d nbx %d\n",nb1,nb2,nbx);}
	E = E_[nb1]; F = F_[nbx];
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->Y_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_AnZtSWnYt_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double)); 
	  lld = (long long int)(D->Z_rpop_j_total - 0) * (long long int)(E_[nb2]->M_Zn->cpop_j) * (long long int)(E_[nb2]->M_Wn->cpop_j);
	  dra_plusdivequals_s___m_m(&(E->QR_AnZtSWnYt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_AnZtSWnYt->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_AnZtSWnYt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnZtSWnYt_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->Y_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->A_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_AnZtSZnAt_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(D->Z_rpop_j_total - 0) * (long long int)(E_[nb2]->M_Zn->cpop_j) * (long long int)(E_[nb2]->M_Zn->cpop_j - 1);
	  dra_plusdivequals_s___m_m(&(E->QR_AnZtSZnAt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_AnZtSZnAt->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_AnZtSZnAt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnZtSZnAt_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->Z_rbother && D->A_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->Y_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_AnAtTYnYt_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(D->A_rpop_j_total - 1) * (long long int)(E_[nb2]->M_An->cpop_j) * (long long int)(E_[nb2]->M_Yn->cpop_j);
	  dra_plusdivequals_s___m_m(&(E->QR_AnAtTYnYt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_AnAtTYnYt->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_AnAtTYnYt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnAtTYnYt_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->Y_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->A_cbother){ fill_uchar_zero((unsigned char *)(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols])),E->A_nrows*D->T_ncols*sizeof(double));
	  lld = (long long int)(D->A_rpop_j_total - 1) * (long long int)(E_[nb2]->M_An->cpop_j) * (long long int)(E_[nb2]->M_An->cpop_j - 1);
	  dra_plusdivequals_s___m_m(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),E->A_nrows,D->T_ncols,F->QR_AnAtTAnAt->lf,(double)lld,E->A_bmr_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(E->QR_AnAtTAnAt_nrm[nb2*E->A_nrows*D->T_ncols]),"double_trn",E->A_nrows,D->T_ncols," %%%% E->QR_AnAtTAnAt_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && E_[nb2]->A_rbother && D->A_cbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->Z_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_AtTYnWtSZn_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(D->A_rpop_j_total) * (long long int)(E_[nb1]->M_Yn->cpop_j - 0) * (long long int)(D->Z_rpop_j_total);
	  dra_plusdivequals_s___m_m(&(D->QC_AtTYnWtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_AtTYnWtSZn->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_AtTYnWtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTYnWtSZn_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->Z_rbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->Z_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_AtTAnZtSZn_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(D->A_rpop_j_total) * (long long int)(E_[nb1]->M_An->cpop_j - 1) * (long long int)(D->Z_rpop_j_total);
	  dra_plusdivequals_s___m_m(&(D->QC_AtTAnZtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_AtTAnZtSZn->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_AtTAnZtSZn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTAnZtSZn_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->Z_rbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->A_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_AtTYnYtTAn_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(D->A_rpop_j_total) * (long long int)(E_[nb1]->M_Yn->cpop_j - 0) * (long long int)(D->A_rpop_j_total - 1);
	  dra_plusdivequals_s___m_m(&(D->QC_AtTYnYtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_AtTYnYtTAn->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_AtTYnYtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTYnYtTAn_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->Y_cbother && E_[nb2]->A_rbother){ } */}
	if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->A_rbother){ fill_uchar_zero((unsigned char *)(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols])),D->A_ncols*D->T_ncols*sizeof(double));
	  lld = (long long int)(D->A_rpop_j_total) * (long long int)(E_[nb1]->M_An->cpop_j - 1) * (long long int)(D->A_rpop_j_total - 1);
	  dra_plusdivequals_s___m_m(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),D->A_ncols,D->T_ncols,F->QC_AtTAnAtTAn->lf,(double)lld,D->A_bmc_j,D->T_bmc_j);
	  if (verbose){ raprintf(&(D->QC_AtTAnAtTAn_nrm[nbx*D->A_ncols*D->T_ncols]),"double_trn",D->A_ncols,D->T_ncols," %%%% D->QC_AtTAnAtTAn_nrm: ");}
	  /* if (E_[nb1]->A_rbother && D->A_cbother && D->A_cbother && E_[nb2]->A_rbother){ } */}
	/* for (nb1=0;nb1<nbins;nb1++){ for (nb2=0;nb2<nbins;nb2++){ }} */}}
    for (ns_a=0;ns_a<D->T_ncols;ns_a++){
      if (bget__on(D->T_bmc_j,ns_a)){
	for (nb1=0;nb1<nbins;nb1++){ E = E_[nb1];     
	  if (verbose){ printf(" %% nb1 %d\n",nb1);}
	  if (E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother){ 
	    ma_a=0; ma_b=0; ma_j=0; 
	    while (ma_a<E->A_nrows){
	      if (bget__on(E->A_bmr_b,ma_a)){
		if (bget__on(E->A_bmr_j,ma_a)){
		  if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnZtSWnYt_nrm[ma_a+ns_a*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
		  dra_sumx(0,&(E->QR_AnZtSWnYt_nrm[ma_a+ns_a*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,nbins);
		  if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnZtSWnYt_nrm[ma_a+ns_a*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
		  ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
		ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	      ma_a++; /* while (ma_a<E->A_nrows){ } */}
	    /* if (E->A_rbother && D->A_cbother && E->Z_rbother && D->Y_cbother){ } */}
	  if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ 
	    ma_a=0; ma_b=0; ma_j=0; 
	    while (ma_a<E->A_nrows){
	      if (bget__on(E->A_bmr_b,ma_a)){
		if (bget__on(E->A_bmr_j,ma_a)){
		  if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnZtSZnAt_nrm[ma_a+ns_a*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
		  dra_sumx(0,&(E->QR_AnZtSZnAt_nrm[ma_a+ns_a*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,nbins);
		  if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnZtSZnAt_nrm[ma_a+ns_a*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
		  ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
		ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	      ma_a++; /* while (ma_a<E->A_nrows){ } */}
	    /* if (E->A_rbother && D->A_cbother && E->Z_rbother && D->A_cbother){ } */}
	  if (E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother){ 
	    ma_a=0; ma_b=0; ma_j=0; 
	    while (ma_a<E->A_nrows){
	      if (bget__on(E->A_bmr_b,ma_a)){
		if (bget__on(E->A_bmr_j,ma_a)){
		  if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnAtTYnYt_nrm[ma_a+ns_a*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
		  dra_sumx(0,&(E->QR_AnAtTYnYt_nrm[ma_a+ns_a*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,nbins);
		  if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnAtTYnYt_nrm[ma_a+ns_a*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
		  ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
		ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	      ma_a++; /* while (ma_a<E->A_nrows){ } */}
	    /* if (E->A_rbother && D->A_cbother && E->A_rbother && D->Y_cbother){ } */}
	  if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ 
	    ma_a=0; ma_b=0; ma_j=0; 
	    while (ma_a<E->A_nrows){
	      if (bget__on(E->A_bmr_b,ma_a)){
		if (bget__on(E->A_bmr_j,ma_a)){
		  if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pre: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnAtTAnAt_nrm[ma_a+ns_a*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
		  dra_sumx(0,&(E->QR_AnAtTAnAt_nrm[ma_a+ns_a*E->A_nrows+0*E->A_nrows*D->T_ncols]),E->A_nrows*D->T_ncols,nbins);
		  if (verbose){ printf(" %% ma_a %d ma_b %d ma_j %d; pos: ",ma_a,ma_b,ma_j); for (nb2=0;nb2<nbins;nb2++){ printf(" %+07.3f",E->QR_AnAtTAnAt_nrm[ma_a+ns_a*E->A_nrows+nb2*E->A_nrows*D->T_ncols]);} printf("\n");}
		  ma_j++; /* if (bget__on(E->A_bmr_j,ma_a)){ } */}
		ma_b++; /* if (bget__on(E->A_bmr_b,ma_a)){ } */}
	      ma_a++; /* while (ma_a<E->A_nrows){ } */}
	    /* if (E->A_rbother && D->A_cbother && E->A_rbother && D->A_cbother){ } */}
	  /* for (nb1=0;nb1<nbins;nb1++){ } */}
	if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->Z_rpop_b_total){ 
	  na_a=0; na_b=0; na_j=0; 
	  while (na_a<D->A_ncols){
	    if (bget__on(D->A_bmc_b,na_a)){
	      if (bget__on(D->A_bmc_j,na_a)){
		if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTYnWtSZn_nrm[na_a+ns_a*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
		dra_sumx(0,&(D->QC_AtTYnWtSZn_nrm[na_a+ns_a*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,nbins*nbins);
		if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTYnWtSZn_nrm[na_a+ns_a*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
		na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	      na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
	    na_a++; /* while (na_a<D->A_ncols){ } */}
	  /* if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->Z_rpop_b_total){ } */}
	if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ 
	  na_a=0; na_b=0; na_j=0; 
	  while (na_a<D->A_ncols){
	    if (bget__on(D->A_bmc_b,na_a)){
	      if (bget__on(D->A_bmc_j,na_a)){
		if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTAnZtSZn_nrm[na_a+ns_a*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
		dra_sumx(0,&(D->QC_AtTAnZtSZn_nrm[na_a+ns_a*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,nbins*nbins);
		if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTAnZtSZn_nrm[na_a+ns_a*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
		na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	      na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
	    na_a++; /* while (na_a<D->A_ncols){ } */}
	  /* if (D->A_rpop_b_total && D->A_cbother && D->A_cbother && D->Z_rpop_b_total){ } */}
	if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->A_rpop_b_total){ 
	  na_a=0; na_b=0; na_j=0; 
	  while (na_a<D->A_ncols){
	    if (bget__on(D->A_bmc_b,na_a)){
	      if (bget__on(D->A_bmc_j,na_a)){
		if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTYnYtTAn_nrm[na_a+ns_a*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
		dra_sumx(0,&(D->QC_AtTYnYtTAn_nrm[na_a+ns_a*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,nbins*nbins);
		if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTYnYtTAn_nrm[na_a+ns_a*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
		na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	      na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
	    na_a++; /* while (na_a<D->A_ncols){ } */}
	  /* if (D->A_rpop_b_total && D->A_cbother && D->Y_cbother && D->A_rpop_b_total){ } */}
	if (D->A_cbother && D->A_rpop_b_total){ 
	  na_a=0; na_b=0; na_j=0; 
	  while (na_a<D->A_ncols){
	    if (bget__on(D->A_bmc_b,na_a)){
	      if (bget__on(D->A_bmc_j,na_a)){
		if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pre: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTAnAtTAn_nrm[na_a+ns_a*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
		dra_sumx(0,&(D->QC_AtTAnAtTAn_nrm[na_a+ns_a*D->A_ncols+0*D->A_ncols*D->T_ncols]),D->A_ncols*D->T_ncols,nbins*nbins);
		if (verbose){ printf(" %% na_a %d na_b %d na_j %d; pos: ",na_a,na_b,na_j); for (nbx=0;nbx<nbins*nbins;nbx++){ printf(" %+07.3f",D->QC_AtTAnAtTAn_nrm[na_a+ns_a*D->A_ncols+nbx*D->A_ncols*D->T_ncols]);} printf("\n");}
		na_j++; /* if (bget__on(D->A_bmc_j,na_a)){ } */}
	      na_b++; /* if (bget__on(E->A_bmc_b,na_a)){ } */}
	    na_a++; /* while (na_a<D->A_ncols){ } */}
	  /* if (D->A_cbother && D->A_rpop_b_total){ } */}
	/* if (bget__on(D->T_bmc_j,ns_a)){ } */}
      /* for (ns_a=0;ns_a<D->T_ncols;ns_a++){ } */}
    /* if (GLOBAL_Ireq<=0){ } */}
  if (verbose){ printf(" %% [finished bcc_sumscores_nrm]\n");}
}

void bcc_sumscores_test()
{
  /* tests a combination of lrup, flattenloop and sumscores */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0); int iteration_max = GLOBAL_TEST_niter;
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  int nl=0; double ct=0,rt=0,r=0;
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  bcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D);
  GLOBAL_toc(1,1+verbose," %% loading time: ");
  if (verbose>1){ printf(" %% calculating initial loop-subscores.\n");}
  GLOBAL_tic(1);
  bcc_An_ajdk(D);
  bcc_lf_ZtSn(D);
  bcc_lf_ZtSWn(D);
  bcc_lf_YnWt(D);
  bcc_M_ZtSWn(D);
  bcc_M_YnWt(D);
  bcc_lf_AtTYn____WtSZn(D);
  bcc_lf_At_T_YnWt_S_Zn(D);
  if (error_check){ bcc_QC_AtTYnWtSZn_uu(D);}
  bcc_lf_AnZt_S_WnYt(D);
  bcc_lf_An_ZtSWn_Yt(D);
  if (error_check){ bcc_QR_AnZtSWnYt_uu(D);}
  GLOBAL_toc(1,1+verbose," %% initial subscore time: ");
  if (verbose>1){ printf(" %% correcting for collapsed-loops.\n");}
  GLOBAL_tic(1);
  bcc_singlestudy_ww(D);
  bcc_doublestudy_ww(D);
  bcc_flattenloop(D);
  if (error_check){ bcc_QC_AtTYnWtSZn_error(verbose,D);}
  if (error_check){ bcc_QR_AnZtSWnYt_error(verbose,D);}
  GLOBAL_toc(1,1+verbose," %% collapsed-loop correction time: ");
  if (verbose>1){ printf(" %% beginning iteration.\n");}
  nl=0;
  while ((iteration_max<=0 || nl<iteration_max) && (D->A_rpop_j_total>0 || D->A_cpop_j>0)){
    GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %0.5fs(%2.1fh)/%0.5fs(%2.1fh) = %.1f\n",nl,iteration_max,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r);}
    if (verbose>1){ printf(" %% combining subscores to form initial loop-scores.\n");}
    bcc_sumscores_nrm(D);
    bcc_sumscores_ifT(D);
    bcc_sumscores_srt(D);
    bcc_sumscores_cmb(D); 
    if (verbose>1){ printf(" %% finding rows and columns with low scores.\n");}
    bcc_sumscores_xij(D);
    bcc_sumscores_mxA(D);
    bcc_sumscores_dmp(D);
    if (verbose>1){ printf(" %% use low-rank-update to calculate new loop-subscores after removing rows and columns with low scores.\n");}
    if (strstr(D->QR_strategy,"condense")){ bcc_lrup_QR_YnWt_stage_2(D); /* if strategy */}
    else{ /* use AnZt_vv */ bcc_lrup_QR_YnWt_stage_0(D); bcc_lrup_QR_YnWt_stage_1(D); /* if strategy */}
    bcc_lrup_QR_YnWt_stage_3(D);
    bcc_lrup_QC_ZtSWn_stage_0(D);
    bcc_lrup_QC_ZtSWn_stage_1(D);
    bcc_lrup_QC_ZtSWn_stage_2(D);
    bcc_lrup_QC_YnWt_stage_a0(D);
    bcc_lrup_QC_YnWt_stage_a1(D);
    bcc_lrup_QC_YnWt_stage_a2(D);
    bcc_lrup_QC_YnWt_stage_a3(D);
    bcc_lrup_QC_YnWt_stage_b0(D);
    if (strstr(D->QC_strategy,"store one")){ bcc_lrup_QC_YnWt_stage_b3(D); /* if strategy */}
    else{ /* store all */ bcc_lrup_QC_YnWt_stage_b1(D); bcc_lrup_QC_YnWt_stage_b2(D); /* if strategy */}
    bcc_lrup_mxdup(D);
    bcc_M_mxset(D);
    bcc_M_ZtSWn(D);
    bcc_M_YnWt(D);
    bcc_lf_AtTYn____WtSZn(D);
    /* bcc_lf_At_T_YnWt_S_Zn(D); */
    bcc_lrup_QC_YnWt_stage_c(D);
    if (error_check){ bcc_QC_AtTYnWtSZn_uu(D);}
    bcc_lf_AnZt_S_WnYt(D);
    /* fix later */ /* bcc_lf_An_ZtSWn_Yt(D); */
    if (error_check){ bcc_QR_AnZtSWnYt_uu(D);}
    if (verbose>1){ printf(" %% correcting for collapsed-loops once again.\n");}
    bcc_An_ajdk(D);
    bcc_lf_ZtSn(D);
    bcc_singlestudy_ww(D);
    bcc_doublestudy_ww(D);
    bcc_flattenloop(D);
    if (error_check){ bcc_QC_AtTYnWtSZn_error(verbose,D);}
    if (error_check){ bcc_QR_AnZtSWnYt_error(verbose,D);}
    nl++; /* while (nl<iteration_max && D->A_cpop_j>2 && D->A_rpop_j_total>2){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  wkspace_printf();
}

void bcc_lrup_sumscores_test()
{
  /* tests a combination of lrup, flattenloop and sumscores */
  int verbose=GLOBAL_verbose; int error_check = (strcmp(GLOBAL_TEST_TYP2,"error")==0); int iteration_max = GLOBAL_TEST_niter; int xdrop_length = 0;
  int *rdrop=NULL,*cdrop=NULL,*rkeep=NULL,*ckeep=NULL,*rcsum=NULL,*ccsum=NULL;
  double tau_c_est=0,tau_r_est=0;
  struct bcc_ajdk *D=NULL;struct bcc_single **E_=NULL; struct bcc_double **F_=NULL;
  int nl=0; double ct=0,rt=0,r=0,it=0,et=0;
  if (error_check){ printf(" %% skipping inital subscore calculation; not checking for errors\n");}
  GLOBAL_tic(2);
  GLOBAL_tic(1);
  bcc_init(GLOBAL_TEST_TYP2,GLOBAL_TEST_mrand,&D,&E_,&F_,GLOBAL_QR_strategy,GLOBAL_QC_strategy); bcc_init_QX(D);
  GLOBAL_toc(1,1+verbose," %% loading time: ");
  get_xdrop_array(D->A_rpop_j_total,D->A_cpop_j,&xdrop_length,&rdrop,&cdrop,&rkeep,&ckeep);
  ira_cumulative_sum(rkeep,xdrop_length,&rcsum); ira_cumulative_sum(ckeep,xdrop_length,&ccsum);
  sprintf(GLOBAL_skip,"AnZt_vv AnAt_vv AtTYn_vv AtTAn_vv xcalc AtTYn____WtSZn_vv At_T_YnWt_S_Zn_ww AnZt_S_WnYt_vv An_ZtSWn_Yt_ww");
  if (verbose>1){ printf(" %% skipping initial loop-subscores.\n");}
  GLOBAL_tic(1);
  bcc_An_ajdk(D);
  bcc_lf_ZtSn(D);
  bcc_lf_ZtSWn(D);
  bcc_lf_YnWt(D);
  bcc_M_ZtSWn(D);
  bcc_M_YnWt(D);
  bcc_lf_AtTYn____WtSZn(D);
  bcc_lf_At_T_YnWt_S_Zn(D);
  bcc_lf_AnZt_S_WnYt(D);
  bcc_lf_An_ZtSWn_Yt(D);
  GLOBAL_toc(1,1+verbose," %% skipping initial subscore : ");
  sprintf(GLOBAL_skip,"\0");
  if (verbose>1){ printf(" %% correcting for collapsed-loops.\n");}
  GLOBAL_tic(1);
  GLOBAL_tic(4); bcc_singlestudy_ww(D); GLOBAL_toc(4,1+verbose," %% singlestudy_ww: ");
  GLOBAL_tic(4); bcc_doublestudy_ww(D); GLOBAL_toc(4,1+verbose," %% doublestudy_ww: ");
  GLOBAL_tic(4); bcc_flattenloop(D); GLOBAL_toc(4,1+verbose," %% flattenloop: ");
  GLOBAL_toc(1,1+verbose," %% collapsed-loop correction time: ");
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt); it = rt;
  if (verbose>1){ printf(" %% elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; beginning iteration.\n",ct,ct/3600,rt,rt/3600,r);}
  GLOBAL_tic(3);
  nl=0;
  while ((iteration_max<=0 || nl<iteration_max) && (D->A_rpop_j_total>0 || D->A_cpop_j>0)){
    GLOBAL_toc(3,0,""); ct = GLOBAL_elct[3]; rt = GLOBAL_elrt[3]; r=ct/maximum(1,rt); tau_c_est = rt/ccsum[nl]; tau_r_est = rt/rcsum[nl]; 
    //et = maximum(rcsum[xdrop_length-1]*tau_r_est,ccsum[xdrop_length-1]*tau_c_est);
    et = ccsum[xdrop_length-1]*tau_c_est;
    if (verbose>-1){ printf(" %% iteration %.5d/%.5d, D %.5d-x-%.7d, elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; estimated total time %6.1fs(%2.1fh) + %6.1fs(%2.1fh)\n",nl,xdrop_length,D->A_rpop_j_total,D->A_cpop_j,ct,ct/3600,rt,rt/3600,r,it,it/3600,et,et/3600);}
    if (verbose>1){ printf(" %% combining subscores to form initial loop-scores.\n");}
    GLOBAL_tic(4); bcc_sumscores_nrm(D); GLOBAL_toc(4,1+verbose," %% sumscores_nrm: ");
    GLOBAL_tic(4); bcc_sumscores_ifT(D); GLOBAL_toc(4,1+verbose," %% sumscores_ifT: ");
    GLOBAL_tic(4); bcc_sumscores_srt(D); GLOBAL_toc(4,1+verbose," %% sumscores_srt: ");
    GLOBAL_tic(4); bcc_sumscores_cmb(D); GLOBAL_toc(4,1+verbose," %% sumscores_cmb: ");
    if (verbose>1){ printf(" %% finding rows and columns with low scores.\n");}
    GLOBAL_tic(4); bcc_sumscores_xij(D); GLOBAL_toc(4,1+verbose," %% sumscores_xij: ");
    GLOBAL_tic(4); bcc_sumscores_mxA(D); GLOBAL_toc(4,1+verbose," %% sumscores_mxA: ");
    GLOBAL_tic(4); bcc_sumscores_dmp(D); GLOBAL_toc(4,1+verbose," %% sumscores_dmp: ");
    if (verbose>1){ printf(" %% use low-rank-update to calculate new loop-subscores after removing rows and columns with low scores.\n");}
    if (strstr(D->QR_strategy,"condense")){ 
      GLOBAL_tic(4); bcc_lrup_QR_YnWt_stage_2(D); GLOBAL_toc(4,1+verbose," %% lrup_QR_YnWt_stage_2: ");
      /* if strategy */}
    else{ /* use AnZt_vv */ 
      GLOBAL_tic(4); bcc_lrup_QR_YnWt_stage_0(D); GLOBAL_toc(4,1+verbose," %% lrup_QR_YnWt_stage_0: ");
      GLOBAL_tic(4); bcc_lrup_QR_YnWt_stage_1(D); GLOBAL_toc(4,1+verbose," %% lrup_QR_YnWt_stage_1: ");
      /* if strategy */}
    GLOBAL_tic(4); bcc_lrup_QR_YnWt_stage_3(D); GLOBAL_toc(4,1+verbose," %% lrup_QR_YnWt_stage_3: ");
    GLOBAL_tic(4); bcc_lrup_QC_ZtSWn_stage_0(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_ZtSWn_stage_0: ");
    GLOBAL_tic(4); bcc_lrup_QC_ZtSWn_stage_1(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_ZtSWn_stage_1: ");
    GLOBAL_tic(4); bcc_lrup_QC_ZtSWn_stage_2(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_ZtSWn_stage_2: ");
    GLOBAL_tic(4); bcc_lrup_QC_YnWt_stage_a0(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_YnWt_stage_a0: ");
    GLOBAL_tic(4); bcc_lrup_QC_YnWt_stage_a1(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_YnWt_stage_a1: ");
    GLOBAL_tic(4); bcc_lrup_QC_YnWt_stage_a2(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_YnWt_stage_a2: ");
    GLOBAL_tic(4); bcc_lrup_QC_YnWt_stage_a3(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_YnWt_stage_a3: ");
    GLOBAL_tic(4); bcc_lrup_QC_YnWt_stage_b0(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_YnWt_stage_b0: ");
    if (strstr(D->QC_strategy,"store one")){ 
      GLOBAL_tic(4); bcc_lrup_QC_YnWt_stage_b3(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_YnWt_stage_b3: ");
      /* if strategy */}
    else{ /* store all */ 
      GLOBAL_tic(4); bcc_lrup_QC_YnWt_stage_b1(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_YnWt_stage_b1: ");
      GLOBAL_tic(4); bcc_lrup_QC_YnWt_stage_b2(D); GLOBAL_toc(4,1+verbose," %% lrup_QC_YnWt_stage_b2: ");
      /* if strategy */}
    GLOBAL_tic(4); bcc_lrup_mxdup(D); GLOBAL_toc(4,1+verbose," %% lrup_mxdup: ");
    GLOBAL_tic(4); bcc_M_mxset(D); GLOBAL_toc(4,1+verbose," %% M_mxset: ");
    GLOBAL_tic(4); bcc_M_ZtSWn(D); GLOBAL_toc(4,1+verbose," %% M_ZtSWn: ");
    GLOBAL_tic(4); bcc_M_YnWt(D); GLOBAL_toc(4,1+verbose," %% M_YnWt: ");
    GLOBAL_tic(4); bcc_lf_AtTYn____WtSZn(D); GLOBAL_toc(4,1+verbose," %% lf_AtTYn____WtSZn: ");
    /* bcc_lf_At_T_YnWt_S_Zn(D); */
    GLOBAL_tic(4); bcc_lrup_QC_YnWt_stage_c(D); GLOBAL_toc(4,1+verbose," %% QC_YnWt_stage_c: ");
    GLOBAL_tic(4); bcc_lf_AnZt_S_WnYt(D); GLOBAL_toc(4,1+verbose," %% lf_AnZt_S_WnYt: ");
    /* fix later */ /* bcc_lf_An_ZtSWn_Yt(D); */
    if (verbose>1){ printf(" %% correcting for collapsed-loops once again.\n");}
    GLOBAL_tic(4); bcc_An_ajdk(D); GLOBAL_toc(4,1+verbose," %% An_ajdk: ");
    GLOBAL_tic(4); bcc_lf_ZtSn(D); GLOBAL_toc(4,1+verbose," %% lf_ZtSn: ");
    GLOBAL_tic(4); bcc_singlestudy_ww(D); GLOBAL_toc(4,1+verbose," %% singlestudy_ww: ");
    GLOBAL_tic(4); bcc_doublestudy_ww(D); GLOBAL_toc(4,1+verbose," %% doublestudy_ww: ");
    GLOBAL_tic(4); bcc_flattenloop(D); GLOBAL_toc(4,1+verbose," %% flattenloop: ");
    nl++; /* while (nl<iteration_max && D->A_cpop_j>2 && D->A_rpop_j_total>2){ } */}
  GLOBAL_toc(2,0,""); ct = GLOBAL_elct[2]; rt = GLOBAL_elrt[2]; r=ct/maximum(1,rt);
  if (verbose>-1){ printf(" %% total elapsed time ct/rt %6.1fs(%2.1fh)/%6.1fs(%2.1fh) = %2.1f; finished.\n",ct,ct/3600,rt,rt/3600,r);}
  wkspace_printf();
}



