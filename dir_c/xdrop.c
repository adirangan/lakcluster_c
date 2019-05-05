#ifndef _MONOLITH
#include "lakcluster_header.h"
#endif /* _MONOLITH */

void get_xdrop(double lrij,double lcij,int *rdrop_p,int *cdrop_p)
{
  int verbose=0; double dbl1=1.000000-0.000001;
  double gamma = GLOBAL_gamma,gamma_tmp_row=0,ammag_tmp_row=0,gamma_tmp_col=0,ammag_tmp_col=0; int rdrop=0;int cdrop=0;
  if (verbose){ printf(" %% [entering get_xdrop] gamma %0.3f, lrij %d lcij %d\n",gamma,(int)lrij,(int)lcij);}
  /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if gamma=0, then we set up gamma_tmp_row to remove a single row ;
    if (gamma>0); gamma_tmp_row = max(gamma,(1-1e-6)/length(tmp_rij)); ammag_tmp_row = min(1-gamma,(length(tmp_rij)-1)/length(tmp_rij));
    elseif gamma<=0; gamma_tmp_row = (1-1e-6)/length(tmp_rij); ammag_tmp_row = (length(tmp_rij)-1)/length(tmp_rij);
    end;%if gamma==0;
    % setting up ammag_tmp_col to remove as many cols as necessary so that log(ncols_pos)/log(nrows_pos) = log(ncols_pre)/log(nrows_pre) ;
    % i.e., log(ammag_tmp_col*ncols_pre)/log(ammag_tmp_row*nrows_pre) = log(ncols_pre)/log(nrows_pre) ;
    % i.e., log(ammag_tmp_col*ncols_pre) = (log(ammag_tmp_row) + log(nrows_pre))*log(ncols_pre)/log(nrows_pre) ;
    % i.e., log(ammag_tmp_col) = log(ammag_tmp_row)*log(ncols_pre)/log(nrows_pre) ;
    % i.e., ammag_tmp_col = exp(log(ammag_tmp_row)*log(ncols_pre)/log(nrows_pre));
    ammag_tmp_col = exp(log(ammag_tmp_row)*log(length(tmp_cij))/log(length(tmp_rij)));
    gamma_tmp_col = 1-ammag_tmp_col-1e-6;
    rdrop = r_rem{iteration}(tmp_rij(1:ceil(gamma_tmp_row*end))); cdrop = c_rem{iteration}(tmp_cij(1:ceil(gamma_tmp_col*end)));
    rkeep = r_rem{iteration}(tmp_rij(ceil(gamma_tmp_row*end)+1:end)); ckeep = c_rem{iteration}(tmp_cij(ceil(gamma_tmp_col*end)+1:end));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Note that we can perform similar row and col removal using sqrt or linear scaling, e.g., :
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M=5e3;N=300e3;g=0.1;nr=M;nc=N;ii=1;nr_(ii)=nr;nc_(ii)=nc;
    while(nr>1); gr=g;ar=min(1-g,(nr-1)/nr);ac=exp(log(ar)*log(nc)/log(nr));gc=1-1e-6-ac;
    rdrop=ceil(gr*nr);cdrop=ceil(gc*nc);nr=nr-rdrop;nc=nc-cdrop;ii=ii+1;nr_(ii)=nr;nc_(ii)=nc;end; % log-scaling ;
    figure;subplot(2,3,1);plot(nr_,nc_,'o-');subplot(2,3,2);plot(log(nr_),log(nc_),'o-'); subplot(2,3,3);plot(sqrt(nr_),sqrt(nc_),'o-'); 
    subplot(2,3,4);plot(1:ii,nr_./nc_,'o-');subplot(2,3,5);plot(1:ii,log(nr_)./log(nc_),'o-');subplot(2,3,6);plot(1:ii,sqrt(nr_)./sqrt(nc_),'o-'); 
    M=5e3;N=300e3;g=0.1;nr=M;nc=N;ii=1;nr_(ii)=nr;nc_(ii)=nc;
    while(nr>1); gr=g;gc=g;
    rdrop=ceil(gr*nr);cdrop=ceil(gc*nc);nr=nr-rdrop;nc=nc-cdrop;ii=ii+1;nr_(ii)=nr;nc_(ii)=nc;end; % linear scaling ;
    figure;subplot(2,3,1);plot(nr_,nc_,'o-');subplot(2,3,2);plot(log(nr_),log(nc_),'o-'); subplot(2,3,3);plot(sqrt(nr_),sqrt(nc_),'o-'); 
    subplot(2,3,4);plot(1:ii,nr_./nc_,'o-');subplot(2,3,5);plot(1:ii,log(nr_)./log(nc_),'o-');subplot(2,3,6);plot(1:ii,sqrt(nr_)./sqrt(nc_),'o-'); 
    M=5e3;N=300e3;g=0.1;nr=M;nc=N;ii=1;nr_(ii)=nr;nc_(ii)=nc;
    while(nr>1); gr=g;ar=min(1-g,(nr-1)/nr);ac=(1/nc)*(sqrt(ar*nr)*sqrt(nc)/sqrt(nr)).^2;gc=1-1e-6-ac;
    rdrop=ceil(gr*nr);cdrop=ceil(gc*nc);nr=nr-rdrop;nc=nc-cdrop;ii=ii+1;nr_(ii)=nr;nc_(ii)=nc;end; % sqrt-scaling ;
    figure;subplot(2,3,1);plot(nr_,nc_,'o-');subplot(2,3,2);plot(log(nr_),log(nc_),'o-'); subplot(2,3,3);plot(sqrt(nr_),sqrt(nc_),'o-'); 
    subplot(2,3,4);plot(1:ii,nr_./nc_,'o-');subplot(2,3,5);plot(1:ii,log(nr_)./log(nc_),'o-');subplot(2,3,6);plot(1:ii,sqrt(nr_)./sqrt(nc_),'o-'); 
  */
  if (lrij<=2){ /* drop everything */ rdrop = lrij; cdrop = lcij;}
  else /* if lrij>2 */{
    if (gamma>0){ gamma_tmp_row = maximum(gamma,(dbl1)/maximum(1,lrij)); ammag_tmp_row = minimum(1-gamma,(lrij-1)/maximum(1,lrij));}
    else /* if gamma<=0 */{ gamma_tmp_row = (dbl1)/maximum(1,lrij); ammag_tmp_row = (lrij-1)/maximum(1,lrij);}
    switch (GLOBAL_gamma_type){
    case GLOBAL_gamma__log: ammag_tmp_col = exp(log(ammag_tmp_row)*log(lcij)/maximum(1,log(maximum(1,lrij)))); break;
    case GLOBAL_gamma_sqrt: ammag_tmp_col = sqrt(ammag_tmp_row*lrij)*sqrt(lcij)/maximum(1,sqrt(maximum(1,lrij)))/((double)maximum(1,lcij)); break;
    case GLOBAL_gamma__lin: ammag_tmp_col = ammag_tmp_row; break;
    case GLOBAL_gamma_rows: ammag_tmp_col = 0; break;
    default: printf(" %% Warning! GLOBAL_gamma_type %d\n",GLOBAL_gamma_type); break; /* switch (GLOBAL_gamma_type){ } */}
    gamma_tmp_col = dbl1-ammag_tmp_col; rdrop = ceil(gamma_tmp_row*lrij); cdrop = ceil(gamma_tmp_col*lcij); if (GLOBAL_gamma_type==GLOBAL_gamma_rows){ cdrop=0;}
    /* if lrij>2 */}
  rdrop = minimum(rdrop,lrij); cdrop = minimum(cdrop,lcij);
  if (verbose>0){ printf(" %% gamma_tmp_row %0.3f ammag_tmp_row %0.3f rdrop %d/%d gamma_tmp_col %0.3f ammag_tmp_col %0.3f cdrop %d/%d\n",gamma_tmp_row,ammag_tmp_row,rdrop,(int)lrij,gamma_tmp_col,ammag_tmp_col,cdrop,(int)lcij);}
  if (rdrop_p!=NULL){ *rdrop_p=rdrop;}
  if (cdrop_p!=NULL){ *cdrop_p=cdrop;}
  if (verbose){ printf(" %% [finished get_xdrop] rdrop %d cdrop %d\n",rdrop,cdrop);}
}

int get_xdrop_length(double lrij,double lcij)
{
  int verbose=0;
  double lrij_tmp=lrij,lcij_tmp=lcij;
  int rdrop=0,cdrop=0;
  int continue_flag=1;
  int length=0;
  if (verbose){ printf(" %% [entering get_xdrop_length] lrij %0.2f lcij %0.2f\n",lrij,lcij);}
  while (continue_flag){
    length++;
    get_xdrop(lrij_tmp,lcij_tmp,&rdrop,&cdrop);
    lrij_tmp-=rdrop; lcij_tmp-=cdrop;
    if (verbose>1){ printf(" %% length %d, rdrop %d cdrop %d lrij_tmp %0.2f lcij_tmp %0.2f\n",length,rdrop,cdrop,lrij_tmp,lcij_tmp);}
    continue_flag = (lrij_tmp>0 && lcij_tmp>0);
    /* while (continue_flag){ } */}
  if (verbose){ printf(" %% [finished get_xdrop_length]\n");}
  return length;
}

void get_xdrop_array(double lrij,double lcij,int *length_p,int **rdrop_p_,int **cdrop_p_,int **rkeep_p_,int **ckeep_p_)
{
  /* Stores rdrop and cdrop values across all iterations.
     Also stores number remaining (rkeep and ckeep).
  */
  int verbose=0;
  int length=0,nl=0,rdrop=0,cdrop=0;
  double lrij_tmp=0,lcij_tmp=0;
  if (verbose){ printf(" %% [entering get_xdrop_array] lrij %d lcij %d\n",(int)lrij,(int)lcij);}
  length = get_xdrop_length(lrij,lcij);
  if (length_p!=NULL){ *length_p = length;}
  if (*rdrop_p_==NULL){ (*rdrop_p_) = (int *) wkspace_all0c(length*sizeof(int));}
  if (*cdrop_p_==NULL){ (*cdrop_p_) = (int *) wkspace_all0c(length*sizeof(int));}
  if (*rkeep_p_==NULL){ (*rkeep_p_) = (int *) wkspace_all0c(length*sizeof(int));}
  if (*ckeep_p_==NULL){ (*ckeep_p_) = (int *) wkspace_all0c(length*sizeof(int));}
  lrij_tmp=lrij; lcij_tmp=lcij;
  for (nl=0;nl<length;nl++){
    get_xdrop(lrij_tmp,lcij_tmp,&rdrop,&cdrop);
    (*rdrop_p_)[nl] = rdrop; (*cdrop_p_)[nl] = cdrop;
    lrij_tmp-=rdrop; lcij_tmp-=cdrop;
    (*rkeep_p_)[nl] = round(lrij_tmp); (*ckeep_p_)[nl] = round(lcij_tmp);
    /* for (nl=0;nl<length;nl++){ } */}  
  if (verbose>1){ raprintf((*rdrop_p_),"int",1,length," %% rdrop_: ");}
  if (verbose>1){ raprintf((*rkeep_p_),"int",1,length," %% rkeep_: ");}
  if (verbose>1){ raprintf((*cdrop_p_),"int",1,length," %% cdrop_: ");}
  if (verbose>1){ raprintf((*ckeep_p_),"int",1,length," %% ckeep_: ");}
  if (verbose){ printf(" %% [finished get_xdrop_array]\n");}  
}

void get_xdrop_array_sub(double lrij,double lcij,int iteration_num,int iteration_max,int iteration_min,int *length_p,int **rdrop_p_,int **cdrop_p_,int **rkeep_p_,int **ckeep_p_)
{
  /* Stores rdrop and cdrop values across a subset of iterations.
     Also stores number remaining (rkeep and ckeep).
     The subset of iterations to be stored starts at iteration_min, and runs up to iteration_max, both values in the range [0,..,length-1].
     The total number of iterations to be stored is iteration_num.
     The jth iteration stored is round(iteration_min + (iteration_max-iteration_min)*j/(iteration_num-1)).
  */
  int verbose=0;
  int iteration_[iteration_num],ni=0;
  int length=0,nl=0,rdrop=0,cdrop=0;
  double lrij_tmp=0,lcij_tmp=0;
  if (verbose){ printf(" %% [entering get_xdrop_array_sub] gamma %0.16f lrij %d lcij %d iteration_num %d iteration_max %d iteration_min %d\n",GLOBAL_gamma,(int)lrij,(int)lcij,iteration_num,iteration_max,iteration_min);}
  length = get_xdrop_length(lrij,lcij); if (length_p!=NULL){ *length_p = length;}
  if (length<iteration_num){ printf(" %% Warning! length %d < iteration_num %d in get_xdrop_array_sub\n",length,iteration_num);}
  if (iteration_num==1){ 
    iteration_[0] = minimum(iteration_min,iteration_max);
    /* if (iteration_num==1){ } */}
  if (iteration_num>1){
    for (nl=0;nl<iteration_num;nl++){
      iteration_[nl] = round(iteration_min + (iteration_max-iteration_min)*(double)nl/(double)(iteration_num-1));
      /* for (nl=0;nl<iteration_num;nl++){ } */}
    /* if (iteration_num>1){ } */}
  if (*rdrop_p_==NULL){ (*rdrop_p_) = (int *) wkspace_all0c(iteration_num*sizeof(int));}
  if (*cdrop_p_==NULL){ (*cdrop_p_) = (int *) wkspace_all0c(iteration_num*sizeof(int));}
  if (*rkeep_p_==NULL){ (*rkeep_p_) = (int *) wkspace_all0c(iteration_num*sizeof(int));}
  if (*ckeep_p_==NULL){ (*ckeep_p_) = (int *) wkspace_all0c(iteration_num*sizeof(int));}
  lrij_tmp=lrij; lcij_tmp=lcij; (*rdrop_p_)[0] = 0; (*cdrop_p_)[0] = 0;
  for (nl=0;nl<iteration_min;nl++){
    get_xdrop(lrij_tmp,lcij_tmp,&rdrop,&cdrop);
    lrij_tmp-=rdrop; lcij_tmp-=cdrop;
    (*rdrop_p_)[0] += rdrop; (*cdrop_p_)[0] += cdrop;
    /* for (nl=0;nl<iteration_min;nl++){ } */}
  (*rkeep_p_)[0] = round(lrij_tmp); (*ckeep_p_)[0] = round(lcij_tmp);
  for (ni=1;ni<iteration_num;ni++){
    (*rdrop_p_)[ni] = 0; (*cdrop_p_)[ni] = 0;
    for (nl=iteration_[ni-1];nl<iteration_[ni];nl++){
      get_xdrop(lrij_tmp,lcij_tmp,&rdrop,&cdrop);
      lrij_tmp-=rdrop; lcij_tmp-=cdrop;
      (*rdrop_p_)[ni] += rdrop; (*cdrop_p_)[ni] += cdrop;
      /* for (nl=iteration_[ni-1];nl<iteration_[ni];nl++){ } */}
    (*rkeep_p_)[ni] = round(lrij_tmp); (*ckeep_p_)[ni] = round(lcij_tmp);
    /* for (ni=1;ni<iteration_num;ni++){ } */}
  if (verbose>1){ raprintf((*rdrop_p_),"int",1,iteration_num," %% rdrop_: ");}
  if (verbose>1){ raprintf((*rkeep_p_),"int",1,iteration_num," %% rkeep_: ");}
  if (verbose>1){ raprintf((*cdrop_p_),"int",1,iteration_num," %% cdrop_: ");}
  if (verbose>1){ raprintf((*ckeep_p_),"int",1,iteration_num," %% ckeep_: ");}
  if (verbose){ printf(" %% [finished get_xdrop_array_sub]\n");}  
}
