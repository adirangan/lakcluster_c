%function test_loader_cluster_collect_6();
% collect results from test_loader_cluster_wrap_1.m ;

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

%%%%%%%%;
u_label_A_ = unique(str_label_A_{1});
n_label_A = length(u_label_A_);
label_A_each__ = zeros(n_u,n_label_A);
for nlabel_A=1:n_label_A;
label_A_each__(:,nlabel_A) = zeros(n_u,1);
tmp_ij_ = find(strcmp(str_label_A_{1},u_label_A_{nlabel_A}));
label_A_each__(tmp_ij_,nlabel_A) = 1;
end;%for nlabel_A=1:n_label_A;
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

str_X_ = {'E'}; n_X = numel(str_X_);
str_Z_ = {'dexcluster_nonbinary_trace_ZRimax_g010_p050_nml3','spectral_isosplit5','tsne50_isosplit5','tsne00_isosplit5','umap00_isosplit5','umap00_hdbscan'}; n_Z = numel(str_Z_);
n_Z_index_ = [2,6,2,2,1,1]; assert(length(n_Z_index_)>=n_Z); n_Z_index_max = max(n_Z_index_);
str_Y_ = {'li16f','li16f_c127_r9','li16f_cRRR_r9','li16f_c012_r9'}; n_Y = numel(str_Y_);
nlpv_XYZ_ = -ones(1+n_label_A,n_X,n_Y,n_Z,n_Z_index_max);
n_cluster_XYZ_ = -ones(n_X,n_Y,n_Z,n_Z_index_max);
flag_rerun_each = 0;
for nX=1:n_X;
str_X = str_X_{nX};
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '));
disp(sprintf(' %% %s: ',str_X));
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '));
%%%%%%%%;
for nZ=1:n_Z;
str_Z = str_Z_{nZ};
n_Z_index = n_Z_index_(nZ);
%%%%%%%%;
for nY=1:n_Y;
str_Y = str_Y_{nY};
%%%%%%%%;
str_mat = sprintf('/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_%s_%s_cluster/%s.mat',str_X,str_Y,str_Z);
if ( exist(str_mat,'file'));
disp(sprintf(' %% %s found, not creating',str_mat));
tmp_mat_ = load(str_mat);
for nZ_index=1:n_Z_index;
tmp_n_cluster = length(unique(tmp_mat_.label_B__{nZ_index}));
tmp_lpv = tmp_mat_.lpv_(nZ_index);
disp(sprintf(' %% %s %s %s: index %d: num clusters %d; p_value exp(%0.6f)',str_X,str_Y,str_Z,nZ_index,tmp_n_cluster,tmp_lpv));
nlpv_XYZ_(1,nX,nY,nZ,nZ_index) = -tmp_lpv;
n_cluster_XYZ_(nX,nY,nZ,nZ_index) = tmp_n_cluster;
end;%for nZ_index=1:n_Z_index;
%%%%%%%%;
str_mat_each = sprintf('/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_%s_%s_cluster/%s_each.mat',str_X,str_Y,str_Z);
if (~flag_rerun_each &&  exist(str_mat_each,'file')); 
%disp(sprintf(' %% %s found, not creating',str_mat_each)); 
tmp_mat_each_ = load(str_mat_each); 
lpv_each_ = tmp_mat_each_.lpv_each_;
lp0_each_ = tmp_mat_each_.lp0_each_;
flag_method_each_ = tmp_mat_each_.flag_method_each_;
end;%if ( exist(str_mat_each,'file')); 
if ( flag_rerun_each || ~exist(str_mat_each,'file')); 
%disp(sprintf(' %% %s not found, creating',str_mat_each)); 
lpv_each_ = zeros(n_label_A,n_Z_index);
lp0_each_ = zeros(n_label_A,n_Z_index);
flag_method_each_ = zeros(n_label_A,n_Z_index);
for nZ_index=1:n_Z_index;
for nlabel_A=1:n_label_A;
disp(sprintf(' %% %s %s %s: each for nZ_index %d nlabel_A %d (str %s)',str_X,str_Y,str_Z,nZ_index,nlabel_A,u_label_A_{nlabel_A}));
[lpv_each_(nlabel_A,nZ_index),lP0_each_(nlabel_A,nZ_index),flag_method_each_(nlabel_A,nZ_index)] = label_to_label_enrichment_2(label_A_each__(:,nlabel_A),tmp_mat_.label_B__{nZ_index});
disp(sprintf(' %% lpv exp(%0.6f) lp0 exp(%0.6f) flag_method %d',lpv_each_(nlabel_A,nZ_index),lP0_each_(nlabel_A,nZ_index),flag_method_each_(nlabel_A,nZ_index)));
end;%for nlabel_A=1:n_label_A;
end;%for nZ_index=1:n_Z_index;
save(str_mat_each,'lpv_each_','lp0_each_','flag_method_each_');
end;%if (~exist(str_mat_each,'file')); 
for nZ_index=1:n_Z_index;
tmp_n_cluster = length(unique(tmp_mat_.label_B__{nZ_index}));
tmp_lpv = tmp_mat_.lpv_(nZ_index);
disp(sprintf(' %% %s %s %s: index %d: num clusters %d; p_value exp(%0.6f)',str_X,str_Y,str_Z,nZ_index,tmp_n_cluster,tmp_lpv));
tmp_str1_ = ''; tmp_str2_ = '';
for nlabel_A=1:n_label_A;
tmp_str1_ = sprintf('%s %.3d',tmp_str1_,nlabel_A);
tmp_str2_ = sprintf('%s %.3d',tmp_str2_,round(-lpv_each_(nlabel_A,nZ_index)));
nlpv_XYZ_(1+nlabel_A,nX,nY,nZ,nZ_index) = -lpv_each_(nlabel_A,nZ_index);
end;%for nlabel_A=1:n_label_A;
disp(sprintf(' %% %s',tmp_str1_));
disp(sprintf(' %% %s',tmp_str2_));
end;%for nZ_index=1:n_Z_index;
end;%if ( exist(str_mat,'file'));
if (~exist(str_mat,'file'));
disp(sprintf(' %% %s not found, skipping',str_mat));
end;%if (~exist(str_mat,'file'));
%%%%%%%%;
end;%for nY=1:n_Y;
%%%%%%%%;
end;%for nZ=1:n_Z;
%%%%%%%%;
end;%for nX=1:n_X;

%%%%%%%%;
% Now create table: ;
%%%%%%%%;
for nX=1:n_X;
str_X = str_X_{nX};
for nZ=1:n_Z;
str_Z = str_Z_{nZ};
n_Z_index = n_Z_index_(nZ);
for nZ_index=1:n_Z_index;
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
disp(sprintf(' %% %s %s nZ_index %d',str_X,str_Z,nZ_index));
disp(sprintf(' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
for nY=1:n_Y;
str_Y = str_Y_{nY}; if (numel(str_Y)<7); str_Y = sprintf('%s\t',str_Y); end;
str_ncl_ = '      '; if (n_cluster_XYZ_(nX,nY,nZ,nZ_index)>-1); str_ncl_ = sprintf('%6.0f',1.0*n_cluster_XYZ_(nX,nY,nZ,nZ_index)); end;
str_lpv_ = '      '; if (nlpv_XYZ_(1,nX,nY,nZ,nZ_index)>-1); str_lpv_ = sprintf('%6.0f',max(0,nlpv_XYZ_(1,nX,nY,nZ,nZ_index))); end;
str_out_ = sprintf(' %% %s:\t #C:%s\t ALL: %s\t EACH: ',str_Y,str_ncl_,str_lpv_);
for nlabel_A=1:n_label_A;
str_lpv_ = '       '; if (nlpv_XYZ_(1+nlabel_A,nX,nY,nZ,nZ_index)>-1); str_lpv_ = sprintf('%6.0f',max(0,nlpv_XYZ_(1+nlabel_A,nX,nY,nZ,nZ_index))); end;
str_out_ = sprintf('%s %s',str_out_,str_lpv_);
end;%for nlabel_A=1:n_label_A;
disp(sprintf('%s',str_out_));
end;%for nY=1:n_Y;
end;%for nZ_index=1:n_Z_index;
end;%for nZ=1:n_Z;
end;%for nX=1:n_X;






