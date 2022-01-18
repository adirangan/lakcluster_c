% Here we test the kmeans clustering algorithm. ; 
% Surprisingly, this algorithm seems very good! ;

M = 1781; N = 2e4; n_rank = 24;
snr = 1.0;
n_cluster = 4;
X_ = linspace(0.3,0.7,n_cluster);
MX_ = ceil(M.^X_); NX_ = ceil(N.^X_); mu_ = sqrt(N)*snr./sqrt(MX_.*NX_);
snr_ = sqrt(MX_.*NX_).*mu_/sqrt(N);

A_n_ = randn(M,N); 
for nj=1:2; 
pf_{nj} = randperm(size(A_n_,nj)); 
[~,pi_{nj}] = sort(pf_{nj}); [~,pi_{nj}] = sort(pi_{nj}); 
end;% for nj=1:2;
B_n_ = cell(n_cluster,1);
B_label_ = zeros(M,1); 
MX_sum = 0; NX_sum = 0;
for ncluster=1:n_cluster;
MX = MX_(ncluster); NX = NX_(ncluster); mu = mu_(ncluster);
B_n_{ncluster}=randn(MX,NX)+mu*ones(MX,1)*(2*(rand(1,NX)>0.5) - 1); 
A_n_(pf_{1}(MX_sum + (1:MX)),pf_{2}(NX_sum + (1:NX))) = B_n_{ncluster};
B_label_(pi_{1}(MX_sum + (1:MX))) = ncluster;
MX_sum = MX_sum + MX; NX_sum = NX_sum + NX;
end;%for ncluster=1:n_cluster;
[label_A_,n_label_A_] = label_num_to_enum_0(B_label_);

[tmp_U_,tmp_S_,tmp_V_] = svds(A_n_,n_rank);

n_k = max(8,2*n_cluster);
lP0__ = zeros(n_rank,n_k); label_B___ = cell(n_rank,n_k); cap___ = cell(n_rank,n_k);
for nrank=1:n_rank;
A_n_sub_ = tmp_U_(:,1:nrank);
for nk=1:n_k;
label_B___{nrank,nk} = kmeans(A_n_sub_,nk);
[lP0__(nrank,nk),cap___{nrank,nk}] = label_to_label_enrichment_lP0(label_A_,label_B___{nrank,nk});
end;%for nk=1:n_k;
end;%for nrank=1:n_rank;

figure(1);clf;
colormap(colormap_beach());
imagesc(lP0__);
xlabel('nk'); ylabel('nrank');
title('lP0');
colorbar;

figure(2);clf;
n_row = 6; n_col = 12;
nrank_ = round(linspace(1,n_rank,n_row));
nk_ = round(linspace(2,n_k,n_col));
for nrow=1:n_row; for ncol=1:n_col;
subplot(n_row,n_col,ncol + (nrow-1)*n_col);
nrank = nrank_(nrow); nk = nk_(ncol);
imagesc(cap___{nrank,nk}./repmat(n_label_A_,1,nk),[0,1]);
%xlabel('discovered'); ylabel('planted');
set(gca,'XTick',[],'YTick',[]);
%title(sprintf('nrank %d nk %d intersection matrix',nrank,nk));
title(sprintf('%d,%d',nrank,nk));
%colorbar;
end;end;%for nrow=1:n_row; for ncol=1:n_col;
figbig;
