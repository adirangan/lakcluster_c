% testing hdbscan. ;

%M = 1781; N = 2e4; n_cluster = 8; snr = 1.0;
M = 178; N = 2e3; n_cluster = 1; snr = 2.0;
rng(1);
[A_n_,label_A_,n_label_A_,pf_,pi_,snr_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);
[A_n_sub10_, ~, ~]=run_umap(A_n_,'verbose','none','n_neighbors',30,'min_dist',0.051,'n_components',10);
hdbscan_sub10 = HDBSCAN(A_n_sub10_); 
hdbscan_sub10.minpts = 10; hdbscan_sub10.fit_model(); hdbscan_sub10.get_best_clusters(); hdbscan_sub10.get_membership();
[A_n_sub02_, ~, ~]=run_umap(A_n_,'verbose','none','n_neighbors',30,'min_dist',0.051,'n_components',02);
hdbscan_sub02 = HDBSCAN(A_n_sub02_); 
hdbscan_sub02.minpts = 10; hdbscan_sub02.fit_model(); hdbscan_sub02.get_best_clusters(); hdbscan_sub02.get_membership();

subplot(1,3,1);
scatter(A_n_sub02_(:,1),A_n_sub02_(:,2),25,label_A_); title('true labels');
subplot(1,3,2);
scatter(A_n_sub02_(:,1),A_n_sub02_(:,2),25,hdbscan_sub10.labels); title('sub10 labels');
subplot(1,3,3);
scatter(A_n_sub02_(:,1),A_n_sub02_(:,2),25,hdbscan_sub02.labels); title('sub02 labels');
