function ...
test_loader_helper_covariate_mahalanobis_cluster_2( ...
 verbose ...
,flag_force_create ...
,flag_force_replot ...
,dir_trunk ...
,dir_mat ...
,dir_jpg ...
,dir_cluster ...
,H_infix ...
,n_HCOV ...
,H_ ...
,H_VariableName_ ...
,n_u ...
,X_infix ...
,rank_estimate_X ...
,X_ ...
,n_Sample_Label ...
,str_Sample_Label_ ...
);

[tmp_UH_,tmp_SH_,tmp_VH_] = svds(H_,n_HCOV);
%[tmp_UX_,tmp_SX_,tmp_VX_] = svds(X_,n_u); %<-- use full rank of X_. ;
[tmp_UX_,tmp_SX_,tmp_VX_] = svds(X_,rank_estimate_X); %<-- use full rank of X_. ;
UX_H_ = inv(tmp_SX_)*transpose(tmp_UX_)*mean_center_0(H_); %<-- explicit mean-centering not required for dexcluster, but perhaps for others. ;
[tmp_U_UX_H_,tmp_S_UX_H_,tmp_V_UX_H_] = svds(UX_H_,n_HCOV);
%%%%%%%%;
% Identify mahalanobis rank of H_. ;
%%%%%%%%;
tmp_fname_mat = sprintf('%s/S_call_U%s_%s_mahalanobis.mat',dir_mat,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(tmp_fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
[ ...
 rank_estimate_sample ...
,s_B_nlp_ ...
,svd_sample__ ...
,eig_tw__ ...
,eig_B_ ...
,kta_opt__ ...
,h_x__ ...
,h_eig__ ...
,h_opt__ ...
] = ...
rank_estimate_sample_0( [ ones(n_u,1) , UX_H_ ] );
save(tmp_fname_mat ...
,'rank_estimate_sample' ...
,'s_B_nlp_' ...
,'svd_sample__' ...
,'eig_tw__' ...
,'eig_B_' ...
,'kta_opt__' ...
,'h_x__' ...
,'h_eig__' ...
,'h_opt__' ...
);
end;%if (~exist(tmp_fname_mat,'file'));
%%%%%%%%;
% spectral-->isosplit5;
%%%%%%%%;
fname_fig = sprintf('%s/U%s_%s_spectral_scatter_FIGA',dir_jpg,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
n_rank_HC = 6;
for nrank_HC=1:n_rank_HC;
opts_isosplit5 = struct('K_init',200,'isocut_threshold',1.0);
label_H__{nrank_HC} = transpose(isosplit5(transpose(tmp_V_UX_H_(:,1:nrank_HC)*tmp_S_UX_H_(1:nrank_HC,1:nrank_HC)),opts_isosplit5));
end;%for nrank_HC=1:n_rank_HC;
for nrank_HC=1:n_rank_HC;
subplot(2,3,nrank_HC); scatter(tmp_V_UX_H_(:,1),tmp_V_UX_H_(:,2),15,label_H__{nrank_HC}); colormap(colormap_beach());
end;%for nrank_HC=1:n_rank_HC;
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% dexnb;
%%%%%%%%;
gamma = 0.00; n_shuffle = 64; p_set = 0.05; n_member_lob = 2;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('U%s_%s_%s',X_infix,H_infix,str_xfix); 
dir_out = []; E_array_base_ = transpose(UX_H_); E_array_r_ij_ = []; E_array_c_ij_ = [];
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_4( ...
 dir_cluster ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,[] ...
,flag_force_create ...
);
%%%%%%%%;
fname_fig = sprintf('%s/U%s_%s_dexnb_scatter_FIGA',dir_jpg,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
p_set_ = [0.05]*0.5.^[0:5];
for np=1:length(p_set_);
tmp_p_set = p_set_(np);
tmp_ZRimax_label_ = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_4( ...
 dir_cluster ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,tmp_p_set ...
,n_member_lob ...
,[] ...
);
subplot(2,3,np);
scatter(tmp_V_UX_H_(:,1),tmp_V_UX_H_(:,2),15,label_str_to_enum_0(tmp_ZRimax_label_),'filled');
colormap(colormap_beach());
title(sprintf('p %0.6f #=%d',tmp_p_set,length(unique(tmp_ZRimax_label_))));
axisnotick;
end;%for np=1:length(p_set_);
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Print output: ;
%%%%%%%%;
p_set = 0.05; 
[ ...
 ZRimax_output_label_ ...
,ZRimax_lpFmax_label_ ...
,ZRimax_lpnext_label_ ...
] = ...
test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_4( ...
 dir_cluster ...
,dir_out ...
,prefix_base ...
,E_array_base_ ...
,E_array_r_ij_ ...
,E_array_c_ij_ ...
,gamma ...
,n_shuffle ...
,p_set ...
,n_member_lob ...
,[] ...
);
ZRimax_enum_ = label_str_to_enum_0(ZRimax_output_label_);
n_QCluster = length(unique(ZRimax_enum_));
fname_tsv = sprintf('%s/n_QCluster_U%s_%s_mahalanobis.tsv',dir_mat,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(fname_tsv,'file'));
disp(sprintf(' %% %s not found, creating',fname_tsv));
save(fname_tsv ...
,'n_QCluster' ...
,'-ascii','-tabs' ...
);
end;%if (~exist(fname_tsv,'file'));
%%%%%%%%;
% plot QC-clusters. ;
%%%%%%%%;
fname_fig = sprintf('%s/U%s_%s_dexnb_scatter_FIGA0',dir_jpg,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
figure(1);clf;
label_plot_recursive_0(ZRimax_output_label_,ZRimax_lpFmax_label_,H_VariableName_);
title(sprintf('U%s_%s_dexnb',X_infix,H_infix),'Interpreter','none');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% plot QC-clusters in PCA-space;
%%%%%%%%;
fname_fig = sprintf('%s/U%s_%s_dexnb_scatter_FIGA1',dir_jpg,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
figure(1);clf;
tmp_Legend_ = num2str(transpose(1:n_QCluster));
tmp_symbol_ = cell(n_QCluster,1);
for nQCluster=0:n_QCluster-1;
if (mod(nQCluster,5)==0); tmp_symbol = 'o'; end;
if (mod(nQCluster,5)==1); tmp_symbol = '^'; end;
if (mod(nQCluster,5)==2); tmp_symbol = 's'; end;
if (mod(nQCluster,5)==3); tmp_symbol = 'p'; end;
if (mod(nQCluster,5)==4); tmp_symbol = 'h'; end;
tmp_symbol_{1+nQCluster} = tmp_symbol;
end;%for nQCluster=0:n_QCluster-1;
subplot(1,3,1);
colormap(colormap_beach());
hold on;
for nQCluster=0:n_QCluster-1;
tmp_ij_ = find(ZRimax_enum_==1+nQCluster);
tmp_symbol = tmp_symbol_{1+nQCluster};
scatter(tmp_SH_(1,1)*tmp_VH_(tmp_ij_,1),tmp_SH_(2,2)*tmp_VH_(tmp_ij_,2),65,ZRimax_enum_(tmp_ij_),tmp_symbol,'filled','MarkerEdgeColor','k');
end;%for nQCluster=0:n_QCluster-1;
title(sprintf('QC clusters (%s orig)',H_infix));
axis equal; axisnotick;
legend(tmp_Legend_,'Location','SouthWest');
subplot(1,3,2);
colormap(colormap_beach());
hold on;
for nQCluster=0:n_QCluster-1;
tmp_ij_ = find(ZRimax_enum_==1+nQCluster);
tmp_symbol = tmp_symbol_{1+nQCluster};
scatter(tmp_V_UX_H_(tmp_ij_,1),tmp_V_UX_H_(tmp_ij_,2),65,ZRimax_enum_(tmp_ij_),tmp_symbol,'filled','MarkerEdgeColor','k');
end;%for nQCluster=0:n_QCluster-1;
title(sprintf('QC clusters (U%s_%s)',X_infix,H_infix),'Interpreter','none');
axis equal; axisnotick;
legend(tmp_Legend_,'Location','SouthEast');
subplot(1,3,3);
[tmp_AUL_n_] = label_pca_plot_0(transpose(UX_H_),label_str_to_enum_0(ZRimax_output_label_));
colormap(colormap_beach());
hold on;
for nQCluster=0:n_QCluster-1;
tmp_ij_ = find(ZRimax_enum_==1+nQCluster);
tmp_symbol = tmp_symbol_{1+nQCluster};
scatter(tmp_AUL_n_(tmp_ij_,1),tmp_AUL_n_(tmp_ij_,2),65,ZRimax_enum_(tmp_ij_),tmp_symbol,'filled','MarkerEdgeColor','k');
end;%for nQCluster=0:n_QCluster-1;
title(sprintf('QC clusters (labeled pca) (U%s_%s)',X_infix,H_infix),'Interpreter','none');
axis equal; axisnotick;
legend(tmp_Legend_,'Location','SouthEast');
set(gcf,'Position',1+[0,0,512*3,512*2]);
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
if (verbose>1);
disp(sprintf(' n_QCluster: %d',n_QCluster));
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
for nQCluster=0:n_QCluster-1;
disp(sprintf(' QCluster #%d:',1+nQCluster));
disp(sprintf('\t%s\n',H_VariableName_{find(ZRimax_enum_==1+nQCluster)}));
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
end;%for nQCluster=0:n_QCluster-1;
end;%if (verbose>1);
%%%%%%%%;
% Identify mahalanobis rank of each QC-cluster. ;
%%%%%%%%;
QCluster_rank_estimate_mahalanobis_ = zeros(n_QCluster,1);
for nQCluster=0:n_QCluster-1;
tmp_h_index_ = efind(ZRimax_enum_==1+nQCluster);
tmp_H_ = H_(:,1+tmp_h_index_);
tmp_UXH_ = inv(tmp_SX_)*transpose(tmp_UX_)*mean_center_0(tmp_H_);
tmp_n_HCOV = numel(tmp_h_index_);
%%%%%%%%;
tmp_fname_mat = sprintf('%s/S_c%d_U%s_%s_mahalanobis.mat',dir_mat,1+nQCluster,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(tmp_fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
[ ...
 rank_estimate_sample ...
,s_B_nlp_ ...
,svd_sample__ ...
,eig_tw__ ...
,eig_B_ ...
,kta_opt__ ...
,h_x__ ...
,h_eig__ ...
,h_opt__ ...
] = ...
rank_estimate_sample_0( [ ones(n_u,1) ,  tmp_UXH_ ] );
save(tmp_fname_mat ...
,'nQCluster' ...
,'rank_estimate_sample' ...
,'s_B_nlp_' ...
,'svd_sample__' ...
,'eig_tw__' ...
,'eig_B_' ...
,'kta_opt__' ...
,'h_x__' ...
,'h_eig__' ...
,'h_opt__' ...
);
end;%if (~exist(tmp_fname_mat,'file'));
tmp_ = load(tmp_fname_mat);
QCluster_rank_estimate_mahalanobis_(1+nQCluster) = tmp_.rank_estimate_sample;
clear tmp_;
%%%%%%%%;
end;%for nQCluster=0:n_QCluster-1;
%%%%%%%%;
% Identify original rank of each QC-cluster. ;
%%%%%%%%;
QCluster_rank_estimate_original_ = zeros(n_QCluster,1);
for nQCluster=0:n_QCluster-1;
tmp_h_index_ = efind(ZRimax_enum_==1+nQCluster);
tmp_H_ = H_(:,1+tmp_h_index_);
tmp_UXH_ = inv(tmp_SX_)*transpose(tmp_UX_)*mean_center_0(tmp_H_);
tmp_n_HCOV = numel(tmp_h_index_);
%%%%%%%%;
tmp_fname_mat = sprintf('%s/S_c%d_U%s_%s_original.mat',dir_mat,1+nQCluster,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(tmp_fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
[ ...
 rank_estimate_sample ...
,s_B_nlp_ ...
,svd_sample__ ...
,eig_tw__ ...
,eig_B_ ...
,kta_opt__ ...
,h_x__ ...
,h_eig__ ...
,h_opt__ ...
] = ...
rank_estimate_sample_0( [ ones(n_u,1) , mean_center_0(tmp_H_) ] );
save(tmp_fname_mat ...
,'nQCluster' ...
,'rank_estimate_sample' ...
,'s_B_nlp_' ...
,'svd_sample__' ...
,'eig_tw__' ...
,'eig_B_' ...
,'kta_opt__' ...
,'h_x__' ...
,'h_eig__' ...
,'h_opt__' ...
);
end;%if (~exist(tmp_fname_mat,'file'));
tmp_ = load(tmp_fname_mat);
QCluster_rank_estimate_original_(1+nQCluster) = tmp_.rank_estimate_sample;
clear tmp_;
%%%%%%%%;
end;%for nQCluster=0:n_QCluster-1;
%%%%%%%%;
% identify and store QC-clusters. ;
%%%%%%%%;
for nQCluster=0:n_QCluster-1;
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
tmp_infix = sprintf('c%d_U%s_%s__',1+nQCluster,X_infix,H_infix);
tmp_h_index_ = efind(ZRimax_enum_==1+nQCluster);
tmp_n_HCOV = numel(tmp_h_index_);
disp(sprintf(' %% nQCluster %d/%d --> %s(%d)',nQCluster,n_QCluster,tmp_infix,tmp_n_HCOV));
disp(sprintf(' %s_%s: QCluster #%d:',X_infix,H_infix,1+nQCluster));
disp(sprintf('\t[rank: mahalanobis %d original %d]',QCluster_rank_estimate_mahalanobis_(1+nQCluster),QCluster_rank_estimate_original_(1+nQCluster)));
disp(sprintf('\t%s\n',H_VariableName_{1+tmp_h_index_}));
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
tmp_H_ = H_(:,1+tmp_h_index_);
fname_tsv = sprintf('%s/%s.tsv',dir_mat,tmp_infix);
if (flag_force_create | flag_force_replot | ~exist(fname_tsv,'file')); 
disp(sprintf(' %% writing %s',fname_tsv));
save(fname_tsv,'tmp_H_','-ascii','-tabs');
end;%if (~exist(fname_tsv,'file'));
clear tmp_H_;
end;%for nQCluster=0:n_QCluster-1;
%%%%%%%%;
% tsne00-->isosplit5;
%%%%%%%%;
fname_fig = sprintf('%s/U%s_%s_tsne00_scatter_FIGA',dir_jpg,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
figure(1);clf;
tmp_Q_ = eye(n_HCOV); %[tmp_Q_,~] = qr(randn(n_HCOV));
UX_H_sub_ = fast_tsne_dr_0(transpose(UX_H_*tmp_Q_),struct('rand_seed',1,'no_dims',2,'theta',0.0));
opts_isosplit5 = struct('K_init',n_HCOV,'isocut_threshold',1.0);
label_UX_H_ = transpose(isosplit5(transpose(UX_H_sub_),opts_isosplit5));
subplot(1,1,1); scatter(UX_H_sub_(:,1),UX_H_sub_(:,2),15,label_UX_H_); colormap(colormap_beach());
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% tsne50-->isosplit5;
%%%%%%%%;
fname_fig = sprintf('%s/U%s_%s_tsne50_scatter_FIGA',dir_jpg,X_infix,H_infix);
if (flag_force_create | flag_force_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
figure(1);clf;
tmp_Q_ = eye(n_HCOV); %[tmp_Q_,~] = qr(randn(n_HCOV));
UX_H_sub_ = fast_tsne_dr_0(transpose(UX_H_*tmp_Q_),struct('rand_seed',1,'no_dims',2,'theta',0.5));
opts_isosplit5 = struct('K_init',n_HCOV,'isocut_threshold',1.0);
label_UX_H_ = transpose(isosplit5(transpose(UX_H_sub_),opts_isosplit5));
subplot(1,1,1); scatter(UX_H_sub_(:,1),UX_H_sub_(:,2),15,label_UX_H_); colormap(colormap_beach());
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;