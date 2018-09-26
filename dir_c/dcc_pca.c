#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void dcc_pca_mxA(struct dcc_ajdk *D,int rdrop,int cdrop)
{
  /* reading rdrop and cdrop from input: ;
     copies D->A_bmc_j to D->A_umc_j ; 
     define D->A_umc_j_rmv and D->A_umc_j_rtn ; removing first cdrop entries of D->QC_index_local_mc_a ;
     copies E->A_bmr_j to E->A_umr_j ; 
     define E->A_umr_j_rmv and E->A_umr_j_rtn ; removing first rdrop entries of D->QR_index_local_mr_a ; using D->QR_index_local_nb to index nb ;
     call dcc_sumscores_mxB ; copying all instances of umc and umr to associated bmc and bmr respectively ;
   */
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_single *E=NULL;
  int nb1=0,nr=0,nc=0;
  if (verbose){ printf(" %% [entering dcc_pca_mxA]\n");}
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
  if (verbose>1){ printf(" %% calling dcc_sumscores_mxB\n");}
  dcc_sumscores_mxB(D); 
  if (verbose>1){ printf(" %% finished dcc_sumscores_mxB\n");}
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
  if (verbose){ printf(" %% [finished dcc_pca_mxA]\n");}
}

void dcc_pca_srt(struct dcc_ajdk *D,int nrows,int *mr_index_local_nb,int *mr_index_local_mr,int ncols,int *mc_index_sort)
{
  /* Copies last D->A_rpop_j_total entries from mr_index_sort (i.e., mr_index_local_nb and mr_index_local_mr) ;
     into D->QR_index_local_nb and D->QR_index_local_mr_a, respectively. ;
     Copies last D->A_cpop_j entries from mc_index_sort ;
     into D->QC_index_local_mc_a. ;
  */
  int verbose=0; 
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb1=0; struct dcc_single *E=NULL;
  int na_j=0,nr=0;
  if (verbose){ printf(" %% [entering dcc_pca_srt]\n");}
  if (mc_index_sort!=NULL && ncols>=D->A_cpop_j){ 
    for (na_j=0;na_j<D->A_cpop_j;na_j++){ D->QC_index_local_mc_a[na_j] = mc_index_sort[ncols-D->A_cpop_j+na_j];} 
    /* if (mc_index_sort!=NULL){ } */}
  if (mr_index_local_nb!=NULL && mr_index_local_mr!=NULL && nrows>=D->A_rpop_j_total){ 
    for (nr=0;nr<D->A_rpop_j_total;nr++){ 
      D->QR_index_local_mr_a[nr] = mr_index_local_mr[nrows-D->A_rpop_j_total+nr];
      D->QR_index_local_nb[nr] = mr_index_local_nb[nrows-D->A_rpop_j_total+nr];
      /* for (nr=0;nr<D->A_rpop_j_total;nr++){ } */}
    /* if (mc_index_sort!=NULL){ } */};
  if (verbose>1){ printf(" %% finished copying mr_index_local_nb,mr_index_local_mr and mc_index_sort\n");}    
  if (verbose>1){
    raprintf(D->QR_index_local_mr_a,   "int",1,D->A_rpop_j_total," %% D->QR_index_local_mr_a pos :");
    raprintf(  D->QR_index_local_nb,   "int",1,D->A_rpop_j_total," %%   D->QR_index_local_nb pos :");
    raprintf(D->QC_index_local_mc_a,   "int",1,D->A_cpop_j," %% D->QC_index_local_mc_a pos :");
    /* if (verbose){ } */}
  if (verbose){ printf(" %% [finished dcc_pca_srt]\n");}
}

void dcc_out_xdrop_lkp(struct dcc_ajdk *D,int nrows,int *mr_index_sort,int **mr_index_local_nb_,int **mr_index_local_mr_)
{
  /* Converts mr_index_sort to lnb (local nb) and index_local_mr (local mr) indices. ;
     Note that all indices are of 'a' type (i.e., referring to mr_a). */
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nl=0,tmpd=0,nb=0,mx_a_tot=0;
  struct dcc_single *E=NULL;
  if (verbose){ printf(" %% [entering dcc_out_xdrop_lkp]\n");}
  mx_a_tot = 0; for (nb=0;nb<nbins;nb++){ E = E_[nb]; mx_a_tot += E->A_nrows; /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% nbins %d, mx_a_tot %d, nrows %d\n",nbins,mx_a_tot,nrows);}
  if (*mr_index_local_nb_==NULL){ (*mr_index_local_nb_) = (int *) wkspace_all0c(nrows*sizeof(int));}
  if (*mr_index_local_mr_==NULL){ (*mr_index_local_mr_) = (int *) wkspace_all0c(nrows*sizeof(int));}
  for (nl=0;nl<nrows;nl++){  
    tmpd = mr_index_sort[nl];
    nb=0; while (tmpd>=E_[nb]->A_nrows){ tmpd-=E_[nb]->A_nrows; nb++;}
    (*mr_index_local_nb_)[nl] = nb; (*mr_index_local_mr_)[nl] = tmpd;
    /* for (nl=0;nl<nrows;nl++){ } */}
  if (verbose>1){ raprintf(*mr_index_local_nb_,"int",1,nrows," %% mr_index_local_nb: ");}
  if (verbose>1){ raprintf(*mr_index_local_mr_,"int",1,nrows," %% mr_index_local_mr: ");}
  if (verbose){ printf(" %% [finished dcc_out_xdrop_lkp]\n");}
}

