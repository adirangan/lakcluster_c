function...
[] = ...
test_loader_cluster_postprocess_1( ...
 dir_base ...
,str_label_A_ ...
,A_n_ ...
,str_xfix ...
,date_diff_threshold ...
,flag_force_create_mat ...
,flag_force_create_tmp ...
);

pre_fname_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
if ( exist(pre_fname_mat,'file'));
disp(sprintf(' %% %s found, postprocessing',pre_fname_mat));

try_fname_pre = sprintf('%s/%s_markergene_B___',dir_base,str_xfix);
[try_flag_skip,try_fname_mat] = open_fname_tmp(try_fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~try_flag_skip;
try;

%%%%%%%%;
% Collect output. ;
%%%%%%%%;
verbose=1;
n_u = numel(str_label_A_{1});
label_A_ = label_str_to_enum_0(str_label_A_{1});
u_label_A_ = unique(str_label_A_{1});
n_label_A = length(u_label_A_);
label_A_each__ = zeros(n_u,n_label_A);
for nlabel_A=0:n_label_A-1;
label_A_each__(:,1+nlabel_A) = zeros(n_u,1);
tmp_index_ = efind(strcmp(str_label_A_{1},u_label_A_{1+nlabel_A}));
label_A_each__(1+tmp_index_,1+nlabel_A) = 1;
end;%for nlabel_A=0:n_label_A-1;
u_label_A_enum_ = zeros(n_label_A,1); 
u_label_A_enum_(1:n_label_A-1) = cellfun(@str2num,u_label_A_(1:n_label_A-1));
u_label_A_enum_(end) = n_label_A;
[~,index_label_A_] = sort(u_label_A_enum_,'ascend'); index_label_A_ = index_label_A_ - 1;
%%%%%%%%;
X_n_ = A_n_; n_X_GENE = size(X_n_,2);
assert(size(X_n_,1)==n_u);
%%%%;
tmp_fname_mat = sprintf('%s/markergene_A__.mat',dir_base);
if (~exist(tmp_fname_mat,'file') |  flag_force_create_mat);
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
%%%%%%%%;
markergene_auc_A__ = zeros(n_X_GENE,n_label_A);
markergene_auz_A__ = zeros(n_X_GENE,n_label_A); %<-- bates-equivalent z-score. ;
markergene_nlp_A__ = zeros(n_X_GENE,n_label_A); %<-- bates-equivalent z-score. ;
for nlabel_A=0:n_label_A-1;
tmp_index_0on_ = efind(label_A_each__(:,1+nlabel_A)==1);
tmp_index_off_ = efind(label_A_each__(:,1+nlabel_A)==0);
tmp_N = min(numel(tmp_index_0on_),numel(tmp_index_off_));
if (verbose); disp(sprintf(' %% nlabel_A %d/%d --> N %d/%d',nlabel_A,n_label_A,tmp_N,n_u)); end;
for nX_GENE=0:n_X_GENE-1;
tmp_auc = auc_0(X_n_(1+tmp_index_off_,1+nX_GENE),X_n_(1+tmp_index_0on_,1+nX_GENE));
tmp_auz = (tmp_auc - 0.5)*sqrt(12*tmp_N);
tmp_nlp = -logp_auc_0(tmp_auc,tmp_N);
markergene_auc_A__(1+nX_GENE,1+nlabel_A) = tmp_auc;
markergene_auz_A__(1+nX_GENE,1+nlabel_A) = tmp_auz;
markergene_nlp_A__(1+nX_GENE,1+nlabel_A) = tmp_nlp;
end;%for nX_GENE=0:n_X_GENE-1;
end;%for nlabel_A=0:n_label_A-1;
%%%%%%%%;
save(tmp_fname_mat ...
,'str_label_A_','n_u','label_A_' ...
,'u_label_A_','n_label_A','label_A_each__' ...
,'u_label_A_enum_','index_label_A_' ...
,'n_X_GENE' ...
,'markergene_auc_A__' ...
,'markergene_auz_A__' ...
,'markergene_nlp_A__' ...
);
end;%if (~exist(tmp_fname_mat,'file') |  flag_force_create_mat);
load(tmp_fname_mat);
%%%%%%%%;
tmp_fname_mat = sprintf('%s/%s.mat',dir_base,str_xfix);
tmp_ = load(tmp_fname_mat);
n_Z_index = numel(tmp_.label_B__);
label_B__ = tmp_.label_B__;
%%%%%%%%;
tmp_fname_mat = sprintf('%s/%s_markergene_B___.mat',dir_base,str_xfix);
if (~exist(tmp_fname_mat,'file') |  flag_force_create_mat);
disp(sprintf(' %% %s not found, creating',tmp_fname_mat));
%%%%%%%%;
lp0_A_each_vs_B__ = zeros(n_label_A,n_Z_index);
lpv_A_each_vs_B__ = zeros(n_label_A,n_Z_index);
for nZ_index=0:n_Z_index-1;
for nlabel_A=0:n_label_A-1;
lp0_A_each_vs_B__(1+nlabel_A,1+nZ_index) = label_to_label_enrichment_lP0(label_A_each__(:,1+nlabel_A),label_B__{1+nZ_index});
lpv_A_each_vs_B__(1+nlabel_A,1+nZ_index) = label_to_label_enrichment_2(label_A_each__(:,1+nlabel_A),label_B__{1+nZ_index});
end;%for nlabel_A=0:n_label_A-1;
end;%for nZ_index=0:n_Z_index-1;
%%%%%%%%;
u_label_B__ = cell(n_Z_index,1);
n_label_B_ = zeros(n_Z_index,1);
label_B_each___ = cell(n_Z_index,1);
nlpv_each___ = cell(n_Z_index,1);
for nZ_index=0:n_Z_index-1;
%%%%%%%%;
u_label_B__{1+nZ_index} = unique(label_B__{1+nZ_index});
n_label_B_(1+nZ_index) = length(u_label_B__{1+nZ_index});
label_B_each___{1+nZ_index} = zeros(n_u,n_label_B_(1+nZ_index));
for nlabel_B=0:n_label_B_(1+nZ_index)-1;
label_B_each___{1+nZ_index}(:,1+nlabel_B) = zeros(n_u,1);
tmp_index_ = efind(label_B__{1+nZ_index}==u_label_B__{1+nZ_index}(1+nlabel_B));
label_B_each___{1+nZ_index}(1+tmp_index_,1+nlabel_B) = 1;
end;%for nlabel_B=0:n_label_B_(1+nZ_index)-1;
%%%%%%%%;
nlpv_each__ = zeros(n_label_A,n_label_B_(1+nZ_index));
for nlabel_A=0:n_label_A-1;
for nlabel_B=0:n_label_B_(1+nZ_index)-1;
nlpv_each__(1+nlabel_A,1+nlabel_B) = -label_pair_enrichment(label_A_each__(:,1+nlabel_A),label_B_each___{1+nZ_index}(:,1+nlabel_B));
end;%for nlabel_B=0:n_label_B_(1+nZ_index)-1;
end;%for nlabel_A=0:n_label_A-1;
nlpv_each___{1+nZ_index} = nlpv_each__;
%%%%%%%%;
markergene_auc_B___{1+nZ_index} = zeros(n_X_GENE,n_label_B_(1+nZ_index));
markergene_auz_B___{1+nZ_index} = zeros(n_X_GENE,n_label_B_(1+nZ_index)); %<-- bates-equivalent z-score. ;
markergene_nlp_B___{1+nZ_index} = zeros(n_X_GENE,n_label_B_(1+nZ_index)); %<-- bates-equivalent z-score. ;
for nlabel_B=0:n_label_B_(1+nZ_index)-1;
tmp_index_0on_ = efind(label_B_each___{1+nZ_index}(:,1+nlabel_B)==1);
tmp_index_off_ = efind(label_B_each___{1+nZ_index}(:,1+nlabel_B)==0);
tmp_N = min(numel(tmp_index_0on_),numel(tmp_index_off_));
if (verbose); disp(sprintf(' %% nlabel_B %d/%d --> N %d/%d',nlabel_B,n_label_B_(1+nZ_index),tmp_N,n_u)); end;
for nX_GENE=0:n_X_GENE-1;
tmp_auc = auc_0(X_n_(1+tmp_index_off_,1+nX_GENE),X_n_(1+tmp_index_0on_,1+nX_GENE));
tmp_auz = (tmp_auc - 0.5)*sqrt(12*tmp_N);
tmp_nlp = -logp_auc_0(tmp_auc,tmp_N);
markergene_auc_B___{1+nZ_index}(1+nX_GENE,1+nlabel_B) = tmp_auc;
markergene_auz_B___{1+nZ_index}(1+nX_GENE,1+nlabel_B) = tmp_auz;
markergene_nlp_B___{1+nZ_index}(1+nX_GENE,1+nlabel_B) = tmp_nlp;
end;%for nX_GENE=0:n_X_GENE-1;
end;%for nlabel_B=0:n_label_B_(1+nZ_index)-1;
%%%%%%%%;
end;%for nZ_index=0:n_Z_index-1;
%%%%%%%%;
save(tmp_fname_mat ...
,'str_label_A_','n_u','label_A_' ...
,'u_label_A_','n_label_A','label_A_each__' ...
,'u_label_A_enum_','index_label_A_' ...
,'n_X_GENE' ...
,'n_Z_index','label_B__' ...
,'u_label_B__','n_label_B_','label_B_each___','nlpv_each___','nlpv_each__' ...
,'markergene_auc_B___' ...
,'markergene_auz_B___' ...
,'markergene_nlp_B___' ...
,'lp0_A_each_vs_B__' ...
,'lpv_A_each_vs_B__' ...
);
end;%if (~exist(tmp_fname_mat,'file') |  flag_force_create_mat);
load(tmp_fname_mat);
%%%%%%%%;

for nZ_index=0:n_Z_index-1;
%%%%%%%%;
tmp_dir_jpg = sprintf('%s/dir_jpg',dir_base);
if (~exist(tmp_dir_jpg,'dir')); sprintf(' %% mkdir %s',tmp_dir_jpg); mkdir(tmp_dir_jpg); end;
fname_fig = sprintf('%s/dir_jpg/%s_z%d_lpv_FIGA',dir_base,str_xfix,1+nZ_index);
if (~exist(sprintf('%s.jpg',fname_fig),'file') |  flag_force_create_mat);
disp(sprintf('%% %s not found, creating',fname_fig));
%%%%%%%%;
[tmp_lP0,tmp_cap_] = label_to_label_enrichment_lP0(label_A_,label_B__{1+nZ_index});
%%%%%%%%;
tmp_index_label_B_ = 0:n_label_B_(1+nZ_index)-1;
index_label_B_ = zeros(n_label_B_(1+nZ_index),1);
for nlabel_A=0:min(n_label_A,n_label_B_(1+nZ_index))-1;
[~,tmp_index] = max(nlpv_each___{1+nZ_index}(1+index_label_A_(1+nlabel_A),1+tmp_index_label_B_)); tmp_index = tmp_index - 1;
index_label_B_(1+nlabel_A) = tmp_index_label_B_(1+tmp_index);
tmp_index_label_B_ = setdiff(tmp_index_label_B_,index_label_B_(1+nlabel_A));
end;%for nlabel_A=0:min(n_label_A,n_label_B_(1+nZ_index))-1;
%%%%%%%%;
figure(1);clf;
figbeach;
subplot(1,2,1);
imagesc(tmp_cap_(1+index_label_A_,1+index_label_B_),[0,100]); colorbar;
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gca,'XTick',1:n_label_B_(1+nZ_index),'XTickLabel',u_label_B__{1+nZ_index}(1+index_label_B_)); xtickangle(90);
axis image; title('cap');
subplot(1,2,2);
imagesc(nlpv_each___{1+nZ_index}(1+index_label_A_,1+index_label_B_),[0,9]); colorbar;
set(gca,'YTick',1:n_label_A,'YTickLabel',u_label_A_(1+index_label_A_)); %xtickangle(90);
set(gca,'XTick',1:n_label_B_(1+nZ_index),'XTickLabel',u_label_B__{1+nZ_index}(1+index_label_B_)); xtickangle(90);
colorbar;
axis image; title('nlpv');
sgtitle(sprintf('%s_z%d',str_xfix,1+nZ_index),'Interpreter','none');
figbig;
disp(sprintf(' %% writing %s.jpg',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file') |  flag_force_create_mat);
%%%%%%%%;
end;%for nZ_index=0:n_Z_index-1;

catch; disp(sprintf(' %% WARNING: error generating %s',try_fname_mat)); end;%try;
close_fname_tmp(try_fname_pre);
end;%if ~try_flag_skip;
clear try_fname_pre try_flag_skip try_fname_mat;

end;%if ( exist(pre_fname_mat,'file'));
