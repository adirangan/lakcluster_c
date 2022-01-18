% Here we test the kmeans clustering algorithm. ; 
% Surprisingly, this algorithm seems very good! ;

M = 1781; N = 2e4; n_rank = 19; %n_rank = ceil(sqrt(N));
snr = 0.4;
n_cluster = 8; rng(1);
[A_n_,label_A_,n_label_A_,pf_,pi_,snr_,MX_,NX_,mu_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);

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

flag_plot=0;
if flag_plot;
figure(1);clf;
colormap(colormap_beach());
imagesc(-lP0__);colorbar;
xlabel('nk'); ylabel('nrank');
title('lP0');
end;%if flag_plot;

flag_plot=0;
if flag_plot;
figure(2);clf;
colormap(colormap_beach());
n_row = 6; n_col = 12;
nk_ = round(linspace(2,n_k,n_row));
nrank_ = round(linspace(1,n_rank,n_col));
for nrow=1:n_row; for ncol=1:n_col;
subplot(n_row,n_col,ncol + (nrow-1)*n_col);
nk = nk_(nrow); nrank = nrank_(ncol); 
imagesc(cap___{nrank,nk}./repmat(n_label_A_,1,nk),[0,1]);
%xlabel('discovered'); ylabel('planted');
set(gca,'XTick',[],'YTick',[]);
%title(sprintf('nk %d nrank %d intersection matrix',nk,nrank));
title(sprintf('k%d r%d',nk,nrank));
%colorbar;
end;end;%for nrow=1:n_row; for ncol=1:n_col;
figbig;
end;%if flag_plot;

flag_plot=1;
if flag_plot;
figure(3);clf;
colormap(colormap_beach());
n_row = 6; n_col = 12;
nk_ = round(linspace(2,n_k,n_row));
nrank_ = round(linspace(1,n_rank,n_col));
for nrow=1:n_row; for ncol=1:n_col;
subplot(n_row,n_col,ncol + (nrow-1)*n_col);
nk = nk_(nrow); nrank = nrank_(ncol); 
tmp_cap_ = cap___{nrank,nk};
tmp_nlA_ = n_label_A_;
tmp_nlB_ = sum(tmp_cap_,1);
tmp_MX_ = tmp_cap_;
tmp_NX_ = repmat(transpose([N,NX_]),1,nk);
tmp_mu_ = repmat(transpose([0,mu_]),1,nk);
%tmp_snr_ = sqrt(tmp_MX_.*tmp_NX_).*tmp_mu_/sqrt(N);
%imagesc(tmp_snr_,snr + 0.15*[-1,+1]);
tmp_snr_ = sqrt(tmp_MX_.*tmp_NX_).*tmp_mu_./sqrt(repmat(tmp_nlB_,1+n_cluster,1));
imagesc(tmp_snr_,snr*sqrt(N)/sqrt(M) + 0.5*[-1,+1]);
%xlabel('discovered'); ylabel('planted');
set(gca,'XTick',[],'YTick',[]);
%title(sprintf('nk %d nrank %d intersection matrix',nk,nrank));
title(sprintf('k%d r%d',nk,nrank));
%colorbar;
end;end;%for nrow=1:n_row; for ncol=1:n_col;
figbig;
end;%if flag_plot;

