function [kld___,n_rank_,kld__,DVS__,SAD2__] = SampleSVD_5(D_,Sample_Label_,alpha_);
% cluster-supervised pca. ;
% Clusters are each distinct gaussians with different means and (slightly) different standard-deviations. ;
% We also modulate the ratio 'alpha' used to define Alpha_. ;
% For each alpha we calculate lp_iso ;

setup;
verbose=0;
[M,N] = size(D_);
u_Sample_Label_ = unique(Sample_Label_); n_Sample_Label = length(u_Sample_Label_);
n_Sample_Label_ = zeros(n_Sample_Label,1);
Sample_Label_ij_ = cell(n_Sample_Label);
for nSample_Label=1:n_Sample_Label;
Sample_Label_ij_{nSample_Label} = find(Sample_Label_==u_Sample_Label_(nSample_Label));
n_Sample_Label_(nSample_Label) = length(Sample_Label_ij_{nSample_Label});
end;%for nSample_Label=1:n_Sample_Label;
n_cluster = n_Sample_Label;
%%%%%%%%;

%%%%%%%%;
kld__ = kld_iso__1(D_,Sample_Label_); kld_50_ = median(kld__,1);
n_alpha = numel(alpha_);
n_rank_ = zeros(n_alpha,1);
kld___ = zeros(n_cluster*(n_cluster-1),N,n_alpha);
SAD2__ = cell(n_alpha,1);
DVS__ = cell(n_alpha,1);
for nalpha = 1:n_alpha;
alpha = alpha_(nalpha);
disp(sprintf(' %% nalpha %d/%d alpha %0.6f',nalpha,n_alpha,alpha));
%%%%%%%%;
Alpha_ = zeros(M,M);
for ncluster1=1:n_cluster; for ncluster2=1:n_cluster;
tmp = (1/alpha - (ncluster1==ncluster2)) / (n_Sample_Label_(ncluster1) * n_Sample_Label_(ncluster2)) ;
Alpha_(Sample_Label_ij_{ncluster1},Sample_Label_ij_{ncluster2}) = tmp;
end;end;%for ncluster1=1:n_cluster; for ncluster2=1:n_cluster;
A2_ = Alpha_ - diag(sum(Alpha_,2)); E_A2_ = eig(A2_);
n_rank_A2 = length(find(abs(E_A2_/E_A2_(1))>1e-6));
disp(sprintf(' %% E_A2_ [%+0.6f,%+0.6f]; pos %d neg %d rank %d',min(E_A2_),max(E_A2_),length(find(E_A2_>0)),length(find(E_A2_<0)),n_rank_A2));
A2fun = @(x,str) transpose(D_)*(A2_*(D_*x));
n_rank = min(M,n_rank_A2); n_rank_(nalpha) = n_rank;
[tmp_VAD_,tmp_SAD2_] = eigs(A2fun,N,n_rank,'largestabs','IsFunctionSymmetric',1,'Display',0);
tmp_SAD2_ = diag(tmp_SAD2_);
tmp_SAD_ = sqrt(abs(tmp_SAD2_)); %<-- pretend we took eigenvalues of positive definite matrix [transpose(D_)*A2_*D_*transpose(D_)*A2_*D_]. ;
tmp_DVS_ = D_*tmp_VAD_(:,1:n_rank)*diag(1./tmp_SAD_(1:n_rank));
%%%%%%%%;
SAD2__{nalpha} = tmp_SAD2_;
DVS__{nalpha} = tmp_DVS_;
kld___(:,1:n_rank,nalpha) = kld_iso__1(tmp_DVS_,Sample_Label_);
%%%%%%%%;
flag_plot=verbose;
if flag_plot;
figure(1);subplot(3,5,nalpha)
scatter(tmp_DVS_(:,1),tmp_DVS_(:,2),15,Sample_Label_,'filled'); colormap(colormap_beach());
axisnotick;
title(sprintf('alpha %0.2f -kld50 %+0.2f',alpha,-median(kld___(:,2,nalpha))));
figbig;drawnow();
end;%if flag_plot;
%%%%%%%%;
end;%for nalpha = 1:n_alpha;
n_rank = max(n_rank_);
kld_50__ = reshape(median(kld___(:,1:n_rank,:),1),n_rank,n_alpha);
%%%%%%%%;

%%%%%%%%;
plot_flag=verbose;
if plot_flag;
figure(2);
c_ = colormap_beach(); n_c = size(c_,1);
hold on;
plot(log2(1:n_rank),-kld_50_(1:n_rank)./[1:n_rank],'k.-');
for nalpha=1:n_alpha;
nc = max(1,min(n_c,floor(n_c*nalpha/n_alpha)));
plot(log2(1:n_rank),-kld_50__(1:n_rank,nalpha)./transpose([1:n_rank]),'o-','Color',c_(nc,:));
end;%for nalpha=1:n_alpha;
xlabel('log2(rank)'); ylabel('kld/dim');
legend(num2str(transpose([0,alpha_])),'Location','NorthEast');
title('kld/dim versus rank');
figbig;
end;%if plot_flag;
%%%%%%%%;




