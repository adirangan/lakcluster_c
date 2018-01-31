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

int get_xdrop_total(double lrij,double lcij)
{
  double lrij_tmp=lrij,lcij_tmp=lcij;
  int rdrop=0,cdrop=0;
  int continue_flag=1;
  int length=0;
  while (continue_flag){
    length++;
    get_xdrop(lrij_tmp,lcij_tmp,&rdrop,&cdrop);
    lrij_tmp-=rdrop; lcij_tmp-=cdrop;
    continue_flag = (lrij_tmp>0 && lcij_tmp>0);
    /* while (continue_flag){ } */}
  return length;
}
