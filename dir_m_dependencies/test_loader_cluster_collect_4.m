function test_loader_cluster_collect_4();
% collect results from test_loader_cluster_wrap_0.m ;

setup;
%%%%%%%%;
dir_trunk = '/data/rangan/dir_bcc/dir_jamison';
str_C_tsv_ = sprintf('%s/dir_mat/C_rank_.tsv',dir_trunk);
C_rank_ = textread(str_C_tsv_); C_rank_ = C_rank_(:,1:end-1);
str_C_nsv_ = sprintf('%s/dir_mat/C_VariableName_.nsv',dir_trunk);
fp = fopen(str_C_nsv_,'r'); C_VariableName_ = textscan(fp,'%s\n'); fclose(fp); C_VariableName_ = C_VariableName_{1};
fp_label_A_ = fopen(sprintf('%s/dir_mat/str_CLabel_sub_.nsv',dir_trunk),'r');
str_label_A_ = textscan(fp_label_A_,'%s'); fclose(fp_label_A_);
n_u = numel(str_label_A_{1});
prefix_normalization = 'li16f';
str_E_tsv_ = sprintf('%s/dir_mat/E_%s_.tsv',dir_trunk,prefix_normalization);
fp = fopen(str_E_tsv_,'r'); tmp_ = fgetl(fp); fclose(fp);
n_E_GENE = numel(str2num(tmp_)); clear tmp_;
str_I_tsv_ = sprintf('%s/dir_mat/I_%s_.tsv',dir_trunk,prefix_normalization);
fp = fopen(str_I_tsv_,'r'); tmp_ = fgetl(fp); fclose(fp);
n_I_GENE = numel(str2num(tmp_)); clear tmp_;
label_A_ = label_str_to_enum_0(str_label_A_{1});
%%%%%%%%;

%{
gamma = 0.01; n_shuffle = 64; p_set = 0.05; n_member_lob = 3;
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
p_set_d = floor(p_set*1000); str_p = sprintf('_p%.3d',p_set_d);
str_nml = sprintf('_nml%d',n_member_lob);
str_xfix = sprintf('dexcluster_nonbinary_trace_ZRimax%s%s%s',str_g,str_p,str_nml);
prefix_base = sprintf('tmp_%s',str_xfix);
dir_base = sprintf('%s/dir_loader_cluster/dir_E_li16f_cluster',dir_trunk);
dir_out = []; E_array_base_ = sparse([],[],[],n_u,n_E_GENE,0); E_array_r_ij_ = []; E_array_c_ij_ = [];
ZRimax_label_ = test_loader_dexcluster_nonbinary_trace_ZRimax_recursive_2(dir_base,dir_out,prefix_base,E_array_base_,E_array_r_ij_,E_array_c_ij_,gamma,n_shuffle,p_set,n_member_lob);
% Note that setting p_set=0.05 reduces the number of groups to 26 or so. ;
 %}

for str_X_ = {'E'};%for str_X_ = {'E','I','EI'};
str_X = str_X_{1};
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '));
disp(sprintf(' %% %s: ',str_X));
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '));
%%%%%%%%;
str_Z_ = {'dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3','spectral_isosplit5','tsne50_isosplit5','tsne00_isosplit5','umap00_isosplit5','umap00_hdbscan'};
n_index_ = [2,6,2,2,1,1];
n_Z = numel(str_Z_);
for nZ=1:n_Z;
str_Z = str_Z_{nZ};
n_index = n_index_(nZ);
%%%%%%%%;
%for str_Y_ = {'li16f','li16f_cPhr_r3','li16f_cPhr_r6','li16f_cPhr_r9'};
%for str_Y_ = {'li16f','li16f_c127_r3','li16f_c127_r6','li16f_c127_r9'};
for str_Y_ = {'li16f'};
str_Y = str_Y_{1};
%%%%%%%%;
str_mat = sprintf('/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_%s_%s_cluster/%s.mat',str_X,str_Y,str_Z);
tmp_mat_ = load(str_mat);
for nindex=1:n_index;
tmp_n_cluster = length(unique(tmp_mat_.label_B__{nindex}));
tmp_lpv = tmp_mat_.lpv_(nindex);
disp(sprintf(' %% %s %s %s: index %d: num clusters %d; p_value exp(%0.6f)',str_X,str_Y,str_Z,nindex,tmp_n_cluster,tmp_lpv));
end;%for nindex=1:n_index;
%%%%%%%%;
end;%for str_Y_ = {'li16f','li16f_cPhr_r3','li16f_cPhr_r6','li16f_cPhr_r9'};
%%%%%%%%;
end;%for nZ=1:n_Z;
%%%%%%%%;
end;%for str_X_ = {'E','I','EI'};

 %%%%%%%%%%%%%%%%%%%%%% 
 % Output: 
 %%%%%%%%%%%%%%%%%%%%%% 
 % E: 
 %%%%%%%%%%%%%%%%%%%%%% 
 % E li16f dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3: index 1: num clusters 23; p_value exp(-616.143198)
 % E li16f dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3: index 2: num clusters 23; p_value exp(-773.968332)
 % E li16f spectral_isosplit5: index 1: num clusters 1; p_value exp(0.000000)
 % E li16f spectral_isosplit5: index 2: num clusters 1; p_value exp(0.000000)
 % E li16f spectral_isosplit5: index 3: num clusters 1; p_value exp(0.000000)
 % E li16f spectral_isosplit5: index 4: num clusters 2; p_value exp(-857.726838)
 % E li16f spectral_isosplit5: index 5: num clusters 3; p_value exp(-1129.727346)
 % E li16f spectral_isosplit5: index 6: num clusters 3; p_value exp(-1141.103937)
 % E li16f tsne50_isosplit5: index 1: num clusters 1; p_value exp(0.000000)
 % E li16f tsne50_isosplit5: index 2: num clusters 1; p_value exp(0.000000)
 % E li16f tsne00_isosplit5: index 1: num clusters 1; p_value exp(0.000000)
 % E li16f tsne00_isosplit5: index 2: num clusters 1; p_value exp(0.000000)
 % E li16f umap00_isosplit5: index 1: num clusters 1; p_value exp(0.000000)
 % E li16f umap00_hdbscan: index 1: num clusters 2; p_value exp(-0.670799)
 %%%%%%%%%%%%%%%%%%%%%% 
 % I: 
 %%%%%%%%%%%%%%%%%%%%%% 
 % I li16f dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3: index 1: num clusters 2; p_value exp(-1.601470)
 % I li16f dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3: index 2: num clusters 2; p_value exp(-3.753418)
 % I li16f spectral_isosplit5: index 1: num clusters 1; p_value exp(0.000000)
 % I li16f spectral_isosplit5: index 2: num clusters 1; p_value exp(0.000000)
 % I li16f spectral_isosplit5: index 3: num clusters 2; p_value exp(-69.214223)
 % I li16f spectral_isosplit5: index 4: num clusters 3; p_value exp(-69.688773)
 % I li16f spectral_isosplit5: index 5: num clusters 3; p_value exp(-72.654802)
 % I li16f spectral_isosplit5: index 6: num clusters 4; p_value exp(-83.036552)
 % I li16f tsne50_isosplit5: index 1: num clusters 3; p_value exp(-82.682014)
 % I li16f tsne50_isosplit5: index 2: num clusters 5; p_value exp(-139.068369)
 % I li16f tsne00_isosplit5: index 1: num clusters 5; p_value exp(-107.146088)
 % I li16f tsne00_isosplit5: index 2: num clusters 5; p_value exp(-136.540016)
 % I li16f umap00_isosplit5: index 1: num clusters 2; p_value exp(-15.369407)
 % I li16f umap00_hdbscan: index 1: num clusters 2; p_value exp(-14.084283)
 %%%%%%%%%%%%%%%%%%%%%% 
 % EI: 
 %%%%%%%%%%%%%%%%%%%%%% 
 % EI li16f dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3: index 1: num clusters 8; p_value exp(-407.512286)
 % EI li16f dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3: index 2: num clusters 8; p_value exp(-427.963913)
 % EI li16f spectral_isosplit5: index 1: num clusters 1; p_value exp(0.000000)
 % EI li16f spectral_isosplit5: index 2: num clusters 1; p_value exp(0.000000)
 % EI li16f spectral_isosplit5: index 3: num clusters 1; p_value exp(0.000000)
 % EI li16f spectral_isosplit5: index 4: num clusters 1; p_value exp(0.000000)
 % EI li16f spectral_isosplit5: index 5: num clusters 2; p_value exp(-72.840032)
 % EI li16f spectral_isosplit5: index 6: num clusters 3; p_value exp(-418.916639)
 % EI li16f tsne50_isosplit5: index 1: num clusters 1; p_value exp(0.000000)
 % EI li16f tsne50_isosplit5: index 2: num clusters 1; p_value exp(0.000000)
 % EI li16f tsne00_isosplit5: index 1: num clusters 4; p_value exp(-258.918404)
 % EI li16f tsne00_isosplit5: index 2: num clusters 1; p_value exp(0.000000)
 % EI li16f umap00_isosplit5: index 1: num clusters 1; p_value exp(0.000000)
 % EI li16f umap00_hdbscan: index 1: num clusters 2; p_value exp(-35.839565)

str_X = 'I';
str_Y = 'li16f';
nZ = 3; str_Z = str_Z_{nZ}; n_index = n_index_(nZ);
str_mat = sprintf('/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_%s_%s_cluster/%s.mat',str_X,str_Y,str_Z);
tmp_mat_ = load(str_mat);
for nindex=1:n_index;
tmp_n_cluster = length(unique(tmp_mat_.label_B__{nindex}));
tmp_lpv = tmp_mat_.lpv_(nindex);
disp(sprintf(' %% %s %s %s: index %d: num clusters %d; p_value exp(%0.6f)',str_X,str_Y,str_Z,nindex,tmp_n_cluster,tmp_lpv));
end;%for nindex=1:n_index;
%%%%%%%%;
nindex = n_index;
num_CLabel_ = tmp_mat_.label_B__{nindex}; 
for nu=1:n_u; str_CLabel_{nu} = num2str(num_CLabel_(nu)); end; 
[AB_] = test_loader_cluster_C_AB_0(1024,str_CLabel_,C_rank_);

%%%%%%%%;
tmp2_ = AB_.AB_C_rawZ_; 
[tmp_U_,tmp_S_,tmp_V_] = svds(tmp2_,3); [~,tmp_U_ij_] = sort(tmp_U_(:,1)); [~,tmp_V_ij_] = sort(tmp_V_(:,1)); 
figure(1); clf;
subplot(1,1,1);
imagesc(tmp2_(tmp_U_ij_,tmp_V_ij_)); colormap(colormap_beach()); colorbar; 
set(gca,'XTick',1:AB_.n_CCOV,'XTickLabel',C_VariableName_(tmp_V_ij_)); xtickangle(90);
set(gca,'YTick',[]); ylabel('pair index');
set(gca,'FontSize',8);
title(sprintf('%s %s %s: AB_C_rawZ_',str_X,str_Y,str_Z),'Interpreter','none');
figbig;
figure(2); clf;
subplot(1,1,1);
scatter3(tmp_V_(:,1),tmp_V_(:,2),tmp_V_(:,3),15,'x'); axis vis3d;
xlabel('pca1'); ylabel('pca2'); zlabel('pca3');
title(sprintf('%s %s %s: AB_C_rawZ_',str_X,str_Y,str_Z),'Interpreter','none');
figbig;
%%%%%%%%;
tmp3_ = corr(tmp2_); 
[tmp_U_,tmp_S_,tmp_V_] = svds(tmp3_,3); [~,tmp_U_ij_] = sort(tmp_U_(:,1)); [~,tmp_V_ij_] = sort(tmp_V_(:,1)); 
figure(3); clf;
subplot(1,1,1);
 imagesc(tmp3_(tmp_V_ij_,tmp_V_ij_),[-1,+1]); colormap(colormap_beach()); colorbar;
set(gca,'XTick',1:AB_.n_CCOV,'XTickLabel',C_VariableName_(tmp_V_ij_)); xtickangle(90);
set(gca,'YTick',1:AB_.n_CCOV,'YTickLabel',C_VariableName_(tmp_U_ij_));
set(gca,'FontSize',8);
title(sprintf('%s %s %s: corr(AB_C_rawZ_)',str_X,str_Y,str_Z),'Interpreter','none');
figbig;
figure(4); clf;
subplot(1,1,1);
scatter3(tmp_V_(:,1),tmp_V_(:,2),tmp_V_(:,3),15,'x'); axis vis3d;
xlabel('pca1'); ylabel('pca2'); zlabel('pca3');
title(sprintf('%s %s %s: corr(AB_C_rawZ_)',str_X,str_Y,str_Z),'Interpreter','none');
figbig;






