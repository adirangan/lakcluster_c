#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

struct R_handle *R_handle_make(char *out_xdrop_fname,unsigned long int rseed)
{
  int verbose=0;
  struct R_handle *R=NULL;
  R = (struct R_handle *)wkspace_all0c(sizeof(struct R_handle)*1);
  if (verbose){ printf(" %% [entering R_handle_make]\n");}
  if (out_xdrop_fname!=NULL){ sprintf(R->out_xdrop_name,out_xdrop_fname);} else{ sprintf(R->out_xdrop_name,"\0");}
  if (verbose>1){ printf(" %% R->out_xdrop_name %s\n",R->out_xdrop_name);}
  R->out_xdrop_nlines=0; /* number of lines in out_xdrop */
  R->out_xdrop_nrows=0; /* number of row-indices stored in out_xdrop */
  R->out_xdrop_ncols=0; /* number of col-indices stored in out_xdrop */
  R->out_xdrop=NULL; /* actual array for out_xdrop */
  R->mr_index_sort=NULL; /* row-indices from out_xdrop */
  R->mc_index_sort=NULL; /* col-indices from out_xdrop */
  R->mr_index_local_nb=NULL; /* local nb from row-indices */
  R->mr_index_local_mr=NULL; /* local mr from row-indices */
  R->rseed = rseed; /* seed to use when scrambling */
  if (R->out_xdrop_name!=NULL &&  strcmp(R->out_xdrop_name,"\0")){ if (verbose>1){ printf(" %% calling xcc_out_xdrop_load\n");} xcc_out_xdrop_load(R->out_xdrop_name,&(R->out_xdrop_nlines),&(R->out_xdrop));}
  if (R->out_xdrop_name==NULL || !strcmp(R->out_xdrop_name,"\0")){ printf(" %% Warning! missing out_xdrop_name in R_handle_make\n");}
  if (verbose>1){ printf(" %% calling xcc_out_xdrop_sort\n");}
  xcc_out_xdrop_sort(R->out_xdrop_nlines,R->out_xdrop,&(R->out_xdrop_nrows),&(R->mr_index_sort),&(R->out_xdrop_ncols),&(R->mc_index_sort));
  iQuickSort(0,R->mr_index_sort,1,0,R->out_xdrop_nrows-1); iQuickSort(0,R->mc_index_sort,1,0,R->out_xdrop_ncols-1);
  if (verbose){ printf(" %% [finished R_handle_make]\n");}
  return R;
}

void R_handle_printf(int verbose,struct R_handle *R,char *prefix)
{
  int nl=0,nx=0;
  char tmpchar[FNAMESIZE];
  if (verbose>0){ printf("%s out_xdrop_name: %s\n",prefix,R->out_xdrop_name);}
  if (verbose>0){ printf("%s nlines %d, nrows %d, ncols %d\n",prefix,R->out_xdrop_nlines,R->out_xdrop_nrows,R->out_xdrop_ncols);}
  if (verbose>0){ sprintf(tmpchar,"%s mc_index_sort: ",prefix); raprintf(R->mc_index_sort,"int",1,R->out_xdrop_ncols,prefix);}
  if (verbose>0){ sprintf(tmpchar,"%s mr_index_sort: ",prefix); raprintf(R->mr_index_sort,"int",1,R->out_xdrop_nrows,tmpchar);}
  if (verbose>0){ sprintf(tmpchar,"%s mr_index_local_nb: ",prefix); raprintf(R->mr_index_local_nb,"int",1,R->out_xdrop_nrows,tmpchar);}
  if (verbose>0){ sprintf(tmpchar,"%s mr_index_local_mr: ",prefix); raprintf(R->mr_index_local_mr,"int",1,R->out_xdrop_nrows,tmpchar);}
}

void dcc_out_xdrop_prune(struct dcc_ajdk *D,int *nrows_,int *mr_index_sort,int *ncols_,int *mc_index_sort)
{
  /* removes indices which are not in mr_b and mc_b masks */
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_;
  int nb=0,mx_a_tot=0;
  struct dcc_single *E=NULL;
  int mx_pre=0,nx_pre=0;
  int mx_pos=0,nx_pos=0;
  int nrows=0,ncols=0;
  int tmpd=0,tmpd2=0;
  if (verbose){ printf(" %% [entering dcc_out_xdrop_prune]\n");}
  if (verbose>1){ raprintf(mr_index_sort,"int",1,*nrows_," %% mr_index_sort: ");}
  if (verbose>1){ raprintf(mc_index_sort,"int",1,*ncols_," %% mc_index_sort: ");}
  mx_a_tot = 0; for (nb=0;nb<nbins;nb++){ E = E_[nb]; mx_a_tot += E->A_nrows; /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% nbins %d, mx_a_tot %d, nrows %d\n",nbins,mx_a_tot,nrows);}
  nx_pre=0; ncols=0;
  for (nx_pos=0;nx_pos<*ncols_;nx_pos++){
    tmpd = mc_index_sort[nx_pos];
    if (tmpd< 0 || tmpd>=D->A_ncols){ if (verbose>1){ printf(" %% nx_pos %d tmpd %d out of bounds\n",nx_pos,tmpd);}}
    if (tmpd>=0 && tmpd< D->A_ncols){
      if (!D->A_umc_b[tmpd]){ if (verbose>1){ printf(" %% nx_pos %d tmpd %d off\n",nx_pos,tmpd);}}
      if ( D->A_umc_b[tmpd]){ if (verbose>1){ printf(" %% nx_pos %d tmpd %d  on; retaining\n",nx_pos,tmpd);} mc_index_sort[nx_pre] = tmpd; nx_pre++; ncols++;}
      /* if (tmpd>=0 && tmpd< D->A_ncols){ } */}
    /* for (nx_pos=0;nx_pos<*ncols_;nx_pos++){ } */}
  *ncols_ = ncols;
  mx_pre=0; nrows=0;
  for (mx_pos=0;mx_pos<*nrows_;mx_pos++){
    tmpd = mr_index_sort[mx_pos];
    if (tmpd< 0 || tmpd>=D->A_nrows_total){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d out of bounds\n",mx_pos,tmpd);}}
    if (tmpd>=0 && tmpd< D->A_nrows_total){
      tmpd2=tmpd; nb=0; while (tmpd2>=E_[nb]->A_nrows){ tmpd2-=E_[nb]->A_nrows;nb++;}
      if (verbose>2){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d\n",mx_pos,tmpd,nb,tmpd2);}
      E = E_[nb];
      if (tmpd2< 0 || tmpd2>=E->A_nrows){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d out of bounds\n",mx_pos,tmpd,nb,tmpd2);}}
      if (tmpd2>=0 && tmpd2< E->A_nrows){
	if (!E->A_umr_b[tmpd2]){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d off\n",mx_pos,tmpd,nb,tmpd2);}}
	if ( E->A_umr_b[tmpd2]){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d  on; retaining\n",mx_pos,tmpd,nb,tmpd2);} mr_index_sort[mx_pre] = tmpd; mx_pre++; nrows++;}
	/* if (tmpd2>=0 && tmpd2< E->A_nrows){ } */}
      /* if (tmpd>=0 && tmpd<D->A_nrows_total){ } */}
    /* for (mx_pos=0;mx_pos<*nrows_;mx_pos++){ } */}
  *nrows_ = nrows;
  if (verbose>1){ raprintf(mr_index_sort,"int",1,*nrows_," %% mr_index_sort: ");}
  if (verbose>1){ raprintf(mc_index_sort,"int",1,*ncols_," %% mc_index_sort: ");}
  if (verbose){ printf(" %% [finished dcc_out_xdrop_prune]\n");}
}

void dcc_scramble(struct dcc_ajdk *D,struct R_handle *R)
{
  int verbose=0;
  int nbins = D->nbins; struct dcc_single **E_ = D->E_; struct dcc_single *E = NULL;
  int nb=0,ma_x=0,na_x=0;
  int ma_a_[nbins],ma_b_[nbins],ma_j_[nbins],ma_a=0;
  int na_a=0,na_b=0,na_j=0;
  unsigned char *An_tag=NULL,*At_tag=NULL;
  int tmpd=0;
  if (verbose){ printf(" %% [entering dcc_scramble]\n");}
  if (verbose){ R_handle_printf(verbose,R," %% R: ");}
  if (verbose){ printf(" %% calculating E_[nb]->lf_Zt_rsum.\n"); wkspace_printf();}
  for (nb=0;nb<nbins;nb++){ L_zero(E_[nb]->lf_Zt_rsum);}
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){
    GLOBAL_pthread_tic(); wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),TYPE_00,SPACING_a,E_[nb]->M_Zt,&(E_[nb]->lf_Zt_rsum)); GLOBAL_pthread_toc();
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  if (verbose){ printf(" %% finished calculating E_[nb]->lf_Zt_rsum.\n"); wkspace_printf();}
  for (nb=0;nb<nbins;nb++){ ma_a_[nb]=0; ma_b_[nb]=0; ma_j_[nb]=0;} 
  nb=0; 
  for (ma_x=0;ma_x<R->out_xdrop_nrows;ma_x++){
    while(nb<R->mr_index_local_nb[ma_x]){ nb++;}
    E = E_[nb];
    while(E->M_An->m_a_[ma_j_[nb]] < R->mr_index_local_mr[ma_x]){ ma_j_[nb] += 1;}
    ma_a_[nb] = E->M_An->m_a_[ma_j_[nb]]; ma_b_[nb] = E->M_An->m_b_[ma_j_[nb]];
    ma_a = ma_a_[nb];
    if (verbose>2){ printf(" %% ma_x %d R->mr_index_local_nb[ma_x] %d nb %d  ma_j %d ma_b %d ma_a %d\n",ma_x,R->mr_index_local_nb[ma_x],nb,ma_j_[nb],ma_b_[nb],ma_a_[nb]);}
    An_tag = (unsigned char *)(&(E->M_An->wX[ma_b_[nb]*E->M_An->mc_length]));
    na_j=0;na_b=0;na_a=0;
    for (na_x=0;na_x<R->out_xdrop_ncols;na_x++){
      while(E->M_An->n_a_[na_j] < R->mc_index_sort[na_x]){ na_j++;}
      na_a = E->M_An->n_a_[na_j]; na_b = E->M_An->n_b_[na_j];
      if (verbose>2){ printf(" %% na_j %d na_b %d na_a %d\n",na_j,na_b,na_a);}
      At_tag = (unsigned char *)(&(E->M_At->wX[na_b*E->M_At->mc_length]));
      /* tmpd = (R01GET(&(R->rseed))<D->A_p[na_a/POPLENGTH]) ? 1 : 0; */
      tmpd = (R01GET(&(R->rseed)) < E->lf_Zt_rsum->lf[na_a]/(double)E->M_Zn->rpop_b) ? 1 : 0;
      if (tmpd==1){ bset__on(An_tag,na_a); bset__on(At_tag,ma_a); /* if (tmpd==1){ } */}
      if (tmpd==0){ bset_off(An_tag,na_a); bset_off(At_tag,ma_a); /* if (tmpd==0){ } */}
      /* for (na_x=0;na_x<R->out_xdrop_ncols;na_x++){ } */}
    /* for (ma_x=0;ma_x<R->out_xdrop_nrows;ma_x++){ } */}  
  if (verbose){ printf(" %% [finished dcc_scramble]\n");}
}

void wrap_dcc_scramble(struct dcc_ajdk *D)
{
  int verbose=0;
  struct R_handle *R=NULL;
  int nb=0,nscramble=0; char fname[FNAMESIZE];
  if (verbose){ printf(" %% [entering wrap_dcc_scramble]\n");}
  if (verbose){
    for (nb=0;nb<D->nbins;nb++){ 
      sprintf(fname,"%s/At_%.2d_scramble_pre.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_An->ncols,D->E_[nb]->M_An->rpop_b,D->E_[nb]->M_An->wX);
      sprintf(fname,"%s/An_%.2d_scramble_pre.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_At->ncols,D->E_[nb]->M_At->rpop_b,D->E_[nb]->M_At->wX);
      /* for (nb=0;nb<D->nbins;nb++){ } */}
    /* if (verbose){ } */}
  for (nscramble=0;nscramble<GLOBAL_scramble_num;nscramble++){
    R = R_handle_make(GLOBAL_scramble_out_xdrop_[nscramble],GLOBAL_scramble_rseed_[nscramble]);
    dcc_out_xdrop_prune(D,&(R->out_xdrop_nrows),R->mr_index_sort,&(R->out_xdrop_ncols),R->mc_index_sort);
    dcc_out_xdrop_lkp(D,R->out_xdrop_nrows,R->mr_index_sort,&(R->mr_index_local_nb),&(R->mr_index_local_mr));
    dcc_scramble(D,R);
    /* for (nscramble=0;nscramble<GLOBAL_scramble_num;nscramble++){ } */}
  if (verbose){
    for (nb=0;nb<D->nbins;nb++){ 
      sprintf(fname,"%s/At_%.2d_scramble_pos.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_An->ncols,D->E_[nb]->M_An->rpop_b,D->E_[nb]->M_An->wX);
      sprintf(fname,"%s/An_%.2d_scramble_pos.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_At->ncols,D->E_[nb]->M_At->rpop_b,D->E_[nb]->M_At->wX);
      /* for (nb=0;nb<D->nbins;nb++){ } */}
    /* if (verbose){ } */}
  if (verbose){ printf(" %% [finished wrap_dcc_scramble]\n");}
}

void bcc_scramble(struct bcc_ajdk *D,struct R_handle *R)
{
  int verbose=0;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_; struct bcc_single *E = NULL;
  int nb=0,ma_x=0,na_x=0;
  int ma_a_[nbins],ma_b_[nbins],ma_j_[nbins],ma_a=0;
  int na_a=0,na_b=0,na_j=0;
  unsigned char *An_tag=NULL,*At_tag=NULL;
  int tmpd=0;
  if (verbose){ printf(" %% [entering bcc_scramble]\n");}
  if (verbose){ R_handle_printf(verbose,R," %% R: ");}
  if (verbose){ printf(" %% calculating E_[nb]->lf_Zt_rsum.\n"); wkspace_printf();}
  for (nb=0;nb<nbins;nb++){ L_zero(E_[nb]->lf_Zt_rsum);}
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){
    GLOBAL_pthread_tic(); wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),TYPE_00,SPACING_a,E_[nb]->M_Zt,&(E_[nb]->lf_Zt_rsum)); GLOBAL_pthread_toc();
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  if (verbose){ printf(" %% finished calculating E_[nb]->lf_Zt_rsum.\n"); wkspace_printf();}
  for (nb=0;nb<nbins;nb++){ ma_a_[nb]=0; ma_b_[nb]=0; ma_j_[nb]=0;} 
  nb=0; 
  for (ma_x=0;ma_x<R->out_xdrop_nrows;ma_x++){
    while(nb<R->mr_index_local_nb[ma_x]){ nb++;}
    E = E_[nb];
    while(E->M_An->m_a_[ma_j_[nb]] < R->mr_index_local_mr[ma_x]){ ma_j_[nb] += 1;}
    ma_a_[nb] = E->M_An->m_a_[ma_j_[nb]]; ma_b_[nb] = E->M_An->m_b_[ma_j_[nb]];
    ma_a = ma_a_[nb];
    if (verbose>2){ printf(" %% ma_x %d R->mr_index_local_nb[ma_x] %d nb %d  ma_j %d ma_b %d ma_a %d\n",ma_x,R->mr_index_local_nb[ma_x],nb,ma_j_[nb],ma_b_[nb],ma_a_[nb]);}
    An_tag = (unsigned char *)(&(E->M_An->wX[ma_b_[nb]*E->M_An->mc_length]));
    na_j=0;na_b=0;na_a=0;
    for (na_x=0;na_x<R->out_xdrop_ncols;na_x++){
      while(E->M_An->n_a_[na_j] < R->mc_index_sort[na_x]){ na_j++;}
      na_a = E->M_An->n_a_[na_j]; na_b = E->M_An->n_b_[na_j];
      if (verbose>2){ printf(" %% na_j %d na_b %d na_a %d\n",na_j,na_b,na_a);}
      At_tag = (unsigned char *)(&(E->M_At->wX[na_b*E->M_At->mc_length]));
      /* tmpd = (R01GET(&(R->rseed))<D->A_p[na_a/POPLENGTH]) ? 1 : 0; */
      tmpd = (R01GET(&(R->rseed)) < E->lf_Zt_rsum->lf[na_a]/(double)E->M_Zn->rpop_b) ? 1 : 0;
      if (tmpd==1){ bset__on(An_tag,na_a); bset__on(At_tag,ma_a); /* if (tmpd==1){ } */}
      if (tmpd==0){ bset_off(An_tag,na_a); bset_off(At_tag,ma_a); /* if (tmpd==0){ } */}
      /* for (na_x=0;na_x<R->out_xdrop_ncols;na_x++){ } */}
    /* for (ma_x=0;ma_x<R->out_xdrop_nrows;ma_x++){ } */}  
  if (verbose){ printf(" %% [finished bcc_scramble]\n");}
}

void bcc_out_xdrop_prune(struct bcc_ajdk *D,int *nrows_,int *mr_index_sort,int *ncols_,int *mc_index_sort)
{
  /* removes indices which are not in mr_b and mc_b masks */
  int verbose=0;
  int nbins = D->nbins; struct bcc_single **E_ = D->E_;
  int nb=0,mx_a_tot=0;
  struct bcc_single *E=NULL;
  int mx_pre=0,nx_pre=0;
  int mx_pos=0,nx_pos=0;
  int nrows=0,ncols=0;
  int tmpd=0,tmpd2=0;
  if (verbose){ printf(" %% [entering bcc_out_xdrop_prune]\n");}
  if (verbose>1){ raprintf(mr_index_sort,"int",1,*nrows_," %% mr_index_sort: ");}
  if (verbose>1){ raprintf(mc_index_sort,"int",1,*ncols_," %% mc_index_sort: ");}
  mx_a_tot = 0; for (nb=0;nb<nbins;nb++){ E = E_[nb]; mx_a_tot += E->A_nrows; /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% nbins %d, mx_a_tot %d, nrows %d\n",nbins,mx_a_tot,nrows);}
  nx_pre=0; ncols=0;
  for (nx_pos=0;nx_pos<*ncols_;nx_pos++){
    tmpd = mc_index_sort[nx_pos];
    if (tmpd< 0 || tmpd>=D->A_ncols){ if (verbose>1){ printf(" %% nx_pos %d tmpd %d out of bounds\n",nx_pos,tmpd);}}
    if (tmpd>=0 && tmpd< D->A_ncols){
      if (!D->A_umc_b[tmpd]){ if (verbose>1){ printf(" %% nx_pos %d tmpd %d off\n",nx_pos,tmpd);}}
      if ( D->A_umc_b[tmpd]){ if (verbose>1){ printf(" %% nx_pos %d tmpd %d  on; retaining\n",nx_pos,tmpd);} mc_index_sort[nx_pre] = tmpd; nx_pre++; ncols++;}
      /* if (tmpd>=0 && tmpd< D->A_ncols){ } */}
    /* for (nx_pos=0;nx_pos<*ncols_;nx_pos++){ } */}
  *ncols_ = ncols;
  mx_pre=0; nrows=0;
  for (mx_pos=0;mx_pos<*nrows_;mx_pos++){
    tmpd = mr_index_sort[mx_pos];
    if (tmpd< 0 || tmpd>=D->A_nrows_total){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d out of bounds\n",mx_pos,tmpd);}}
    if (tmpd>=0 && tmpd< D->A_nrows_total){
      tmpd2=tmpd; nb=0; while (tmpd2>=E_[nb]->A_nrows){ tmpd2-=E_[nb]->A_nrows;nb++;}
      if (verbose>2){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d\n",mx_pos,tmpd,nb,tmpd2);}
      E = E_[nb];
      if (tmpd2< 0 || tmpd2>=E->A_nrows){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d out of bounds\n",mx_pos,tmpd,nb,tmpd2);}}
      if (tmpd2>=0 && tmpd2< E->A_nrows){
	if (!E->A_umr_b[tmpd2]){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d off\n",mx_pos,tmpd,nb,tmpd2);}}
	if ( E->A_umr_b[tmpd2]){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d  on; retaining\n",mx_pos,tmpd,nb,tmpd2);} mr_index_sort[mx_pre] = tmpd; mx_pre++; nrows++;}
	/* if (tmpd2>=0 && tmpd2< E->A_nrows){ } */}
      /* if (tmpd>=0 && tmpd<D->A_nrows_total){ } */}
    /* for (mx_pos=0;mx_pos<*nrows_;mx_pos++){ } */}
  *nrows_ = nrows;
  if (verbose>1){ raprintf(mr_index_sort,"int",1,*nrows_," %% mr_index_sort: ");}
  if (verbose>1){ raprintf(mc_index_sort,"int",1,*ncols_," %% mc_index_sort: ");}
  if (verbose){ printf(" %% [finished bcc_out_xdrop_prune]\n");}
}

void wrap_bcc_scramble(struct bcc_ajdk *D)
{
  int verbose=1;
  struct R_handle *R=NULL;
  int nb=0,nscramble=0; char fname[FNAMESIZE];
  if (verbose){ printf(" %% [entering wrap_bcc_scramble]\n");}
  if (verbose){
    for (nb=0;nb<D->nbins;nb++){ 
      sprintf(fname,"%s/At_%.2d_scramble_pre.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_An->ncols,D->E_[nb]->M_An->rpop_b,D->E_[nb]->M_An->wX);
      sprintf(fname,"%s/An_%.2d_scramble_pre.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_At->ncols,D->E_[nb]->M_At->rpop_b,D->E_[nb]->M_At->wX);
      /* for (nb=0;nb<D->nbins;nb++){ } */}
    /* if (verbose){ } */}
  for (nscramble=0;nscramble<GLOBAL_scramble_num;nscramble++){
    R = R_handle_make(GLOBAL_scramble_out_xdrop_[nscramble],GLOBAL_scramble_rseed_[nscramble]);
    bcc_out_xdrop_prune(D,&(R->out_xdrop_nrows),R->mr_index_sort,&(R->out_xdrop_ncols),R->mc_index_sort);
    bcc_out_xdrop_lkp(D,R->out_xdrop_nrows,R->mr_index_sort,&(R->mr_index_local_nb),&(R->mr_index_local_mr));
    bcc_scramble(D,R);
    /* for (nscramble=0;nscramble<GLOBAL_scramble_num;nscramble++){ } */}
  if (verbose){
    for (nb=0;nb<D->nbins;nb++){ 
      sprintf(fname,"%s/At_%.2d_scramble_pos.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_An->ncols,D->E_[nb]->M_An->rpop_b,D->E_[nb]->M_An->wX);
      sprintf(fname,"%s/An_%.2d_scramble_pos.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_At->ncols,D->E_[nb]->M_At->rpop_b,D->E_[nb]->M_At->wX);
      /* for (nb=0;nb<D->nbins;nb++){ } */}
    /* if (verbose){ } */}
  if (verbose){ printf(" %% [finished wrap_bcc_scramble]\n");}
}

void dcg_out_xdrop_prune(struct dcg_ajdk *D,int *nrows_,int *mr_index_sort,int *ncols_,int *mc_index_sort)
{
  /* removes indices which are not in mr_b and mc_b masks */
  int verbose=0;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_;
  int nb=0,mx_a_tot=0;
  struct dcg_single *E=NULL;
  int mx_pre=0,nx_pre=0;
  int mx_pos=0,nx_pos=0;
  int nrows=0,ncols=0;
  int tmpd=0,tmpd2=0;
  if (verbose){ printf(" %% [entering dcg_out_xdrop_prune]\n");}
  if (verbose>1){ raprintf(mr_index_sort,"int",1,*nrows_," %% mr_index_sort: ");}
  if (verbose>1){ raprintf(mc_index_sort,"int",1,*ncols_," %% mc_index_sort: ");}
  mx_a_tot = 0; for (nb=0;nb<nbins;nb++){ E = E_[nb]; mx_a_tot += E->A_nrows; /* for (nb=0;nb<nbins;nb++){ } */}
  if (verbose){ printf(" %% nbins %d, mx_a_tot %d, nrows %d\n",nbins,mx_a_tot,nrows);}
  nx_pre=0; ncols=0;
  for (nx_pos=0;nx_pos<*ncols_;nx_pos++){
    tmpd = mc_index_sort[nx_pos];
    if (tmpd< 0 || tmpd>=D->A_ncols){ if (verbose>1){ printf(" %% nx_pos %d tmpd %d out of bounds\n",nx_pos,tmpd);}}
    if (tmpd>=0 && tmpd< D->A_ncols){
      if (!D->A_umc_b[tmpd]){ if (verbose>1){ printf(" %% nx_pos %d tmpd %d off\n",nx_pos,tmpd);}}
      if ( D->A_umc_b[tmpd]){ if (verbose>1){ printf(" %% nx_pos %d tmpd %d  on; retaining\n",nx_pos,tmpd);} mc_index_sort[nx_pre] = tmpd; nx_pre++; ncols++;}
      /* if (tmpd>=0 && tmpd< D->A_ncols){ } */}
    /* for (nx_pos=0;nx_pos<*ncols_;nx_pos++){ } */}
  *ncols_ = ncols;
  mx_pre=0; nrows=0;
  for (mx_pos=0;mx_pos<*nrows_;mx_pos++){
    tmpd = mr_index_sort[mx_pos];
    if (tmpd< 0 || tmpd>=D->A_nrows_total){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d out of bounds\n",mx_pos,tmpd);}}
    if (tmpd>=0 && tmpd< D->A_nrows_total){
      tmpd2=tmpd; nb=0; while (tmpd2>=E_[nb]->A_nrows){ tmpd2-=E_[nb]->A_nrows;nb++;}
      if (verbose>2){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d\n",mx_pos,tmpd,nb,tmpd2);}
      E = E_[nb];
      if (tmpd2< 0 || tmpd2>=E->A_nrows){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d out of bounds\n",mx_pos,tmpd,nb,tmpd2);}}
      if (tmpd2>=0 && tmpd2< E->A_nrows){
	if (!E->A_umr_b[tmpd2]){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d off\n",mx_pos,tmpd,nb,tmpd2);}}
	if ( E->A_umr_b[tmpd2]){ if (verbose>1){ printf(" %% mx_pos %d tmpd %d nb %d tmpd2 %d  on; retaining\n",mx_pos,tmpd,nb,tmpd2);} mr_index_sort[mx_pre] = tmpd; mx_pre++; nrows++;}
	/* if (tmpd2>=0 && tmpd2< E->A_nrows){ } */}
      /* if (tmpd>=0 && tmpd<D->A_nrows_total){ } */}
    /* for (mx_pos=0;mx_pos<*nrows_;mx_pos++){ } */}
  *nrows_ = nrows;
  if (verbose>1){ raprintf(mr_index_sort,"int",1,*nrows_," %% mr_index_sort: ");}
  if (verbose>1){ raprintf(mc_index_sort,"int",1,*ncols_," %% mc_index_sort: ");}
  if (verbose){ printf(" %% [finished dcg_out_xdrop_prune]\n");}
}

void dcg_scramble(struct dcg_ajdk *D,struct R_handle *R)
{
  int verbose=0;
  int nbins = D->nbins; struct dcg_single **E_ = D->E_; struct dcg_single *E = NULL;
  int nb=0,ma_x=0,na_x=0;
  int ma_a_[nbins],ma_b_[nbins],ma_j_[nbins],ma_a=0;
  int na_a=0,na_b=0,na_j=0;
  unsigned char *An_tag=NULL,*At_tag=NULL;
  int tmpd=0;
  if (verbose){ printf(" %% [entering dcg_scramble]\n");}
  if (verbose){ R_handle_printf(verbose,R," %% R: ");}
  if (verbose){ printf(" %% calculating E_[nb]->lf_Zt_rsum.\n"); wkspace_printf();}
  for (nb=0;nb<nbins;nb++){ L_zero(E_[nb]->lf_Zt_rsum);}
  GLOBAL_nf_cur=0; GLOBAL_nf_opn=0;
  for (nb=0;nb<nbins;nb++){
    GLOBAL_pthread_tic(); wrap_An_v__run(&(GLOBAL_tint[GLOBAL_nf_cur]),GLOBAL_tvp[GLOBAL_nf_cur],&(GLOBAL_threads[GLOBAL_nf_cur]),TYPE_00,SPACING_a,E_[nb]->M_Zt,&(E_[nb]->lf_Zt_rsum)); GLOBAL_pthread_toc();
    /* for (nb=0;nb<nbins;nb++){ } */}
  GLOBAL_pthread_tuc();
  if (verbose){ printf(" %% finished calculating E_[nb]->lf_Zt_rsum.\n"); wkspace_printf();}
  for (nb=0;nb<nbins;nb++){ ma_a_[nb]=0; ma_b_[nb]=0; ma_j_[nb]=0;} 
  nb=0; 
  for (ma_x=0;ma_x<R->out_xdrop_nrows;ma_x++){
    while(nb<R->mr_index_local_nb[ma_x]){ nb++;}
    E = E_[nb];
    while(E->M_An->m_a_[ma_j_[nb]] < R->mr_index_local_mr[ma_x]){ ma_j_[nb] += 1;}
    ma_a_[nb] = E->M_An->m_a_[ma_j_[nb]]; ma_b_[nb] = E->M_An->m_b_[ma_j_[nb]];
    ma_a = ma_a_[nb];
    if (verbose>2){ printf(" %% ma_x %d R->mr_index_local_nb[ma_x] %d nb %d  ma_j %d ma_b %d ma_a %d\n",ma_x,R->mr_index_local_nb[ma_x],nb,ma_j_[nb],ma_b_[nb],ma_a_[nb]);}
    An_tag = (unsigned char *)(&(E->M_An->wX[ma_b_[nb]*E->M_An->mc_length]));
    na_j=0;na_b=0;na_a=0;
    for (na_x=0;na_x<R->out_xdrop_ncols;na_x++){
      while(E->M_An->n_a_[na_j] < R->mc_index_sort[na_x]){ na_j++;}
      na_a = E->M_An->n_a_[na_j]; na_b = E->M_An->n_b_[na_j];
      if (verbose>2){ printf(" %% na_j %d na_b %d na_a %d\n",na_j,na_b,na_a);}
      At_tag = (unsigned char *)(&(E->M_At->wX[na_b*E->M_At->mc_length]));
      /* tmpd = (R01GET(&(R->rseed))<D->A_p[na_a/POPLENGTH]) ? 1 : 0; */
      tmpd = (R01GET(&(R->rseed)) < E->lf_Zt_rsum->lf[na_a]/(double)E->M_Zn->rpop_b) ? 1 : 0;
      if (tmpd==1){ bset__on(An_tag,na_a); bset__on(At_tag,ma_a); /* if (tmpd==1){ } */}
      if (tmpd==0){ bset_off(An_tag,na_a); bset_off(At_tag,ma_a); /* if (tmpd==0){ } */}
      /* for (na_x=0;na_x<R->out_xdrop_ncols;na_x++){ } */}
    /* for (ma_x=0;ma_x<R->out_xdrop_nrows;ma_x++){ } */}  
  if (verbose){ printf(" %% [finished dcg_scramble]\n");}
}

void wrap_dcg_scramble(struct dcg_ajdk *D)
{
  int verbose=0;
  struct R_handle *R=NULL;
  int nb=0,nscramble=0; char fname[FNAMESIZE];
  if (verbose){ printf(" %% [entering wrap_dcg_scramble]\n");}
  if (verbose){
    for (nb=0;nb<D->nbins;nb++){ 
      sprintf(fname,"%s/At_%.2d_scramble_pre.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_An->ncols,D->E_[nb]->M_An->rpop_b,D->E_[nb]->M_An->wX);
      sprintf(fname,"%s/An_%.2d_scramble_pre.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_At->ncols,D->E_[nb]->M_At->rpop_b,D->E_[nb]->M_At->wX);
      /* for (nb=0;nb<D->nbins;nb++){ } */}
    /* if (verbose){ } */}
  for (nscramble=0;nscramble<GLOBAL_scramble_num;nscramble++){
    R = R_handle_make(GLOBAL_scramble_out_xdrop_[nscramble],GLOBAL_scramble_rseed_[nscramble]);
    dcg_out_xdrop_prune(D,&(R->out_xdrop_nrows),R->mr_index_sort,&(R->out_xdrop_ncols),R->mc_index_sort);
    dcg_out_xdrop_lkp(D,R->out_xdrop_nrows,R->mr_index_sort,&(R->mr_index_local_nb),&(R->mr_index_local_mr));
    dcg_scramble(D,R);
    /* for (nscramble=0;nscramble<GLOBAL_scramble_num;nscramble++){ } */}
  if (verbose){
    for (nb=0;nb<D->nbins;nb++){ 
      sprintf(fname,"%s/At_%.2d_scramble_pos.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_An->ncols,D->E_[nb]->M_An->rpop_b,D->E_[nb]->M_An->wX);
      sprintf(fname,"%s/An_%.2d_scramble_pos.b16",GLOBAL_DIR_NAME,nb);
      binary_write(fname,BITJ,D->E_[nb]->M_At->ncols,D->E_[nb]->M_At->rpop_b,D->E_[nb]->M_At->wX);
      /* for (nb=0;nb<D->nbins;nb++){ } */}
    /* if (verbose){ } */}
  if (verbose){ printf(" %% [finished wrap_dcg_scramble]\n");}
}
