if ~exist('index_sort_','var');
[~,index_sort_] = sort(mean(B_,1),'descend'); index_sort_ = index_sort_-1;
end;%if ~exist('index_sort_','var');

ng0 = 10; ng1 = 11; 

n_h = 32;
%%%%;
b0_ = B_(:,1+index_sort_(1+ng0)); %<-- one gene across all samples. ;
b1_ = B_(:,1+index_sort_(1+ng1)); %<-- one gene across all samples. ;
blim0_ = prctile(b0_,[0,99]);
blim1_ = prctile(b1_,[0,99]);
hb01__ = hist2d_0(b0_,b1_,n_h,n_h,blim0_,blim1_);
%%%%;
c0_ = log(1+b0_);
c1_ = log(1+b1_);
clim0_ = prctile(c0_,[0,99]);
clim1_ = prctile(c1_,[0,99]);
hc01__ = hist2d_0(c0_,c1_,n_h,n_h,clim0_,clim1_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
figure(1);clf;figbig;figbeach();
p_row = 2; p_col = 3; ns=0;
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
plot(b0_,b1_,'.');
xlim(blim0_); ylim(blim1_);
title('unnormalized scatterplot');
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
imagesc(hb01__);colorbar;
set(gca,'Ydir','normal');
title('unnormalized histogram');
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
imagesc(log(1+hb01__));colorbar;
set(gca,'Ydir','normal');
title('unnormalized log-histogram');
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
plot(c0_,c1_,'.');
xlim(clim0_); ylim(clim1_);
title('log-normalized scatterplot');
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
imagesc(hc01__);colorbar;
set(gca,'Ydir','normal');
title('log-normalized histogram');
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
imagesc(log(1+hc01__));colorbar;
set(gca,'Ydir','normal');
title('log-normalized log-histogram');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
[n_smp,n_gene] = size(B_);
rng(1);
figure(1);clf;figbig;figbeach();
p_row = 6; p_col = 8; ns=0; 
n_h = 32;
n_iteration = floor(p_row*p_col/2);
for niteration=0:n_iteration-1;
ng0 = max(0,min(n_gene-1,floor(n_gene/100*rand())));
ng1 = max(0,min(n_gene-1,floor(n_gene/100*rand())));
%%%%;
b0_ = B_(:,1+index_sort_(1+ng0)); %<-- one gene across all samples. ;
b1_ = B_(:,1+index_sort_(1+ng1)); %<-- one gene across all samples. ;
blim0_ = prctile(b0_,[0,99]);
blim1_ = prctile(b1_,[0,99]);
hb01__ = hist2d_0(b0_,b1_,n_h,n_h,blim0_,blim1_);
%%%%;
c0_ = log(1+b0_);
c1_ = log(1+b1_);
clim0_ = prctile(c0_,[0,99]);
clim1_ = prctile(c1_,[0,99]);
hc01__ = hist2d_0(c0_,c1_,n_h,n_h,clim0_,clim1_);
%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
imagesc(hc01__);colorbar;
set(gca,'Ydir','normal');
title(sprintf('h (%d,%d)',ng0,ng1));
%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
imagesc(log(1+hc01__));colorbar;
set(gca,'Ydir','normal');
title(sprintf('lh (%d,%d)',ng0,ng1));
%%%%;
drawnow;
%%%%;
end;%for niteration=0:n_iteration-1;
