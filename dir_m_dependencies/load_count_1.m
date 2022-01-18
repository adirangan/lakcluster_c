clear;

platform = 'rusty';%platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

verbose=1;
flag_recalc = 0;
flag_replot = 0;
flag_center = 0;
tolerance_master = 1e-2;
nf=0;

dir_code = sprintf('/%s/rangan/dir_bcc/dir_halfloop_dev',string_root);
str_code = sprintf('%s/halfloop_dev',dir_code);
dir_weilin = sprintf('/%s/rangan/dir_bcc/dir_weilin',string_root);
dir_mat = sprintf('%s/dir_mat',dir_weilin); if (~exist(dir_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_mat)); mkdir(dir_mat); end;
dir_jpg = sprintf('%s/dir_jpg',dir_weilin); if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
fname_count_sparse = sprintf('%s/count_sparse_.mat',dir_mat);
if (~exist('B_','var')); load(fname_count_sparse); end;
fname_count_logn_sparse = sprintf('%s/count_logn_sparse_.mat',dir_mat);
if (~exist('B_logn_','var')); load(fname_count_logn_sparse); end;

n_smp = size(B_,1);
n_umi = size(B_,2);
B_sum_1_ = sum(B_,1);
B_sum_2_ = sum(B_,2);

fname_count_pasca = sprintf('%s/count_pasca_.csv',dir_weilin);
fp = fopen(fname_count_pasca);
label_pasca_ = textscan(fp,'%s','Delimiter','\n','headerlines',1);
fclose(fp);
label_pasca_ = label_pasca_{1};
label_pasca_ = label_pasca_(1:n_smp);
[ ...
 label_pasca_enum_ ...
,n_u_label_pasca ...
,u_label_pasca_ ...
,index_nu_pasca_from_nall_ ...
,n_u_label_pasca_ ...
,index_nall_from_nu_pasca__ ...
] = ...
label_str_to_enum_1( ...
 label_pasca_ ...
);

fname_count_stage = sprintf('%s/count_stage_.csv',dir_weilin);
fp = fopen(fname_count_stage);
label_stage_ = textscan(fp,'%s','Delimiter','\n','headerlines',1);
fclose(fp);
label_stage_ = label_stage_{1};
label_stage_ = label_stage_(1:n_smp);
[ ...
 label_stage_enum_ ...
,n_u_label_stage ...
,u_label_stage_ ...
,index_nu_stage_from_nall_ ...
,n_u_label_stage_ ...
,index_nall_from_nu_stage__ ...
] = ...
label_str_to_enum_1( ...
 label_stage_ ...
);

fname_count_Isogenic_type = sprintf('%s/count_Isogenic_type_.csv',dir_weilin);
fp = fopen(fname_count_Isogenic_type);
label_Isogenic_type_ = textscan(fp,'%s','Delimiter','\n','headerlines',1);
fclose(fp);
label_Isogenic_type_ = label_Isogenic_type_{1};
label_Isogenic_type_ = label_Isogenic_type_(1:n_smp);
[ ...
 label_Isogenic_type_enum_ ...
,n_u_label_Isogenic_type ...
,u_label_Isogenic_type_ ...
,index_nu_Isogenic_type_from_nall_ ...
,n_u_label_Isogenic_type_ ...
,index_nall_from_nu_Isogenic_type__ ...
] = ...
label_str_to_enum_1( ...
 label_Isogenic_type_ ...
);

label_pasca_Isog_stage_ = cell(n_smp,1);
for nsmp=0:n_smp-1;
label_pasca_Isog_stage_{1+nsmp} = strcat(label_pasca_{1+nsmp},'_',label_Isogenic_type_{1+nsmp},'_',label_stage_{1+nsmp});
end;%for nsmp=0:n_smp-1;
[ ...
 label_pasca_Isog_stage_enum_ ...
,n_u_label_pasca_Isog_stage ...
,u_label_pasca_Isog_stage_ ...
,index_nu_pasca_Isog_stage_from_nall_ ...
,n_u_label_pasca_Isog_stage_ ...
,index_nall_from_nu_pasca_Isog_stage__ ...
] = ...
label_str_to_enum_1( ...
 label_pasca_Isog_stage_ ...
);

index_stage_3mo_ = efind(label_stage_enum_==1);
fname_fig = sprintf('%s/pca_3mo_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
[UB__,SB__,VB__] = svds(B_logn_,3); SB_ = diag(SB__);
figure(1+nf);nf=nf+1;figbig;figbeach();
subplot(2,2,1);
imagesc(log(1+hist2d_0(UB__(1+index_stage_3mo_,1),UB__(1+index_stage_3mo_,2),32,32)));
title('pca');
subplot(2,2,2);
scatter3(UB__(1+index_stage_3mo_,1),UB__(1+index_stage_3mo_,2),UB__(1+index_stage_3mo_,3),12,label_pasca_enum_(1+index_stage_3mo_),'filled');fig80s;
title('pasca');
subplot(2,2,3);
scatter3(UB__(1+index_stage_3mo_,1),UB__(1+index_stage_3mo_,2),UB__(1+index_stage_3mo_,3),12,label_Isogenic_type_enum_(1+index_stage_3mo_),'filled');fig80s;
title('Isog');
subplot(2,2,4);
scatter3(UB__(1+index_stage_3mo_,1),UB__(1+index_stage_3mo_,2),UB__(1+index_stage_3mo_,3),12,label_pasca_Isog_stage_enum_(1+index_stage_3mo_),'filled');fig80s;
title('pasca+Isog');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

n_svd = 36;
n_smp_3mo = numel(index_stage_3mo_);
B_3mo_ = B_(1+index_stage_3mo_,:);
B_3mo_logn_ = B_logn_(1+index_stage_3mo_,:);

fname_mat = sprintf('%s/SB_3mo_logn__.mat',dir_mat);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
[UB_3mo_logn__,SB_3mo_logn__,VB_3mo_logn__] = svds(B_3mo_logn_,n_svd); SB_3mo_logn_ = diag(SB_3mo_logn__);
save(fname_mat ...
     ,'n_svd' ...
     ,'UB_3mo_logn__' ...
     ,'SB_3mo_logn_' ...
     ,'VB_3mo_logn__' ...
     );
end;%if (~exist(fname_mat,'file'));
load(fname_mat);

fname_mat = sprintf('%s/SB_3mo_logn_shuffle__.mat',dir_mat);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
n_shuffle = 32;
SB_3mo_logn_shuffle__ = zeros(n_svd,n_shuffle);
for nshuffle=0:n_shuffle-1;
rng(nshuffle*1024);
tmp_B_3mo_logn_ = B_3mo_logn_;
for numi=0:n_umi-1;
if (mod(numi,1024)==0); disp(sprintf(' %% nshuffle %d/%d, numi %d/%d',nshuffle,n_shuffle,numi,n_umi)); end;
tmp_ij_ = randperm(n_smp_3mo);
tmp_B_3mo_logn_(:,1+numi) = tmp_B_3mo_logn_(tmp_ij_,1+numi);
end;%for numi=0:n_umi-1;
[tmp_UB_3mo_logn__,tmp_SB_3mo_logn__,tmp_VB_3mo_logn__] = svds(tmp_B_3mo_logn_,n_svd); tmp_SB_3mo_logn_ = diag(tmp_SB_3mo_logn__);
SB_3mo_logn_shuffle__(:,1+nshuffle) = tmp_SB_3mo_logn_;
end;%for nshuffle=0:n_shuffle-1;
save(fname_mat ...
     ,'n_svd' ...
     ,'n_shuffle' ...
     ,'SB_3mo_logn_shuffle__' ...
     );
end;%if (~exist(fname_mat,'file'));
load(fname_mat);

%%%%%%%%;
% Apparent rank is 21. ;
%%%%%%%%;
fname_fig = sprintf('%s/SB_3mo_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;
semilogy(1:n_svd,SB_3mo_logn_,'kx',1:n_svd,SB_3mo_logn_shuffle__,'ro');
title('SB_3mo_logn_','Interpreter','none');
xlabel('rank'); ylabel('sigma');
xlim([0.5,n_svd+0.5]);
grid on;
set(gcf,'Position',1+[0,0,512,1024]);
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));

rank_estimate_B_3mo_logn = 21;

fname_fig = sprintf('%s/nnz_B_3mo_hist_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figsml;
n_N = 64; tmp_h_ = hist(B_3mo_(:),[0:n_N]);
loglog(0:n_N,tmp_h_,'o','MarkerFaceColor','k');
xlabel('count'); ylabel('histogram');
grid on;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% apparent cutoff is at 1, not 2. ;
%%%%%%%%

fname_mat = sprintf('%s/B_3mo_logn_sg__.mat',dir_mat);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
nnz_3mo_gene_ = sum(B_3mo_~=0,1);
index_gene_retain_ = efind(nnz_3mo_gene_>=1); %<-- at least 1 nonzero entries, since 1 does not seem to deviate from algebraic decay. ;
B_3mo_logn_sg__ = B_3mo_logn_(:,1+index_gene_retain_);
[UB_3mo_logn__,SB_3mo_logn__,VB_3mo_logn__] = svds(B_3mo_logn_sg__,2); SB_3mo_logn_ = diag(SB_3mo_logn__);
save(fname_mat ...
     ,'nnz_3mo_gene_' ...
     ,'index_gene_retain_' ...
     ,'B_3mo_logn_sg__' ...
     ,'n_svd' ...
     ,'UB_3mo_logn__' ...
     ,'SB_3mo_logn_' ...
     ,'VB_3mo_logn__' ...
     );
end;%if (~exist(fname_mat,'file'));
load(fname_mat);

label_B_3mo_ = label_pasca_Isog_stage_(1+index_stage_3mo_);
[ ...
 label_B_3mo_enum_ ...
,n_u_label_B_3mo ...
,u_label_B_3mo_ ...
,index_nu_B_3mo_from_nall_ ...
,n_u_label_B_3mo_ ...
,index_nall_from_nu_B_3mo__ ...
] = ...
label_str_to_enum_1( ...
 label_B_3mo_ ...
);
%%%%%%%%;
[~,index_sort_u_label_B_3mo_] = sort(n_u_label_B_3mo_,'descend'); index_sort_u_label_B_3mo_ = index_sort_u_label_B_3mo_ - 1;
index_nu_srt_from_nu_ori_ = index_sort_u_label_B_3mo_;
n_all = sum(n_u_label_B_3mo_);
for nu_label_B_3mo=0:n_u_label_B_3mo-1;
tmp_n_u = n_u_label_B_3mo_(1+index_sort_u_label_B_3mo_(1+nu_label_B_3mo));
u_label_B_3mo = u_label_B_3mo_{1+index_sort_u_label_B_3mo_(1+nu_label_B_3mo)};
disp(sprintf(' %% nu %.2d/%.2d: %.4d/%.4d = %0.2f = %s',nu_label_B_3mo,n_u_label_B_3mo,tmp_n_u,n_all,tmp_n_u/n_all,u_label_B_3mo));
end;%for nu_label_B_3mo=0:n_u_label_B_3mo-1;
%%%%%%%%;
fname_fig = sprintf('%s/B_3mo_logn_sg_pca_FIGA',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+512,1024]);colormap('lines');
markersize_final = 12;
markersize_start = 4;
str_symbol_ = {'o','^','s','p','h'};
c_lines__ = colormap('lines'); n_c_lines = size(c_lines__,1);
hold on;
for nl=0:numel(index_sort_u_label_B_3mo_)-1;
nc_lines = nl;
str_symbol = str_symbol_{1+mod(nl,5)};
markersize_use = round(markersize_start + (markersize_final-markersize_start)*nl/(numel(index_sort_u_label_B_3mo_)-1));
tmp_nu_label_B_3mo = index_sort_u_label_B_3mo_(1+nl);
tmp_index_ = index_nall_from_nu_B_3mo__{1+tmp_nu_label_B_3mo};
plot(UB_3mo_logn__(1+tmp_index_,1),UB_3mo_logn__(1+tmp_index_,2),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_lines__(1+nc_lines,:),'MarkerEdgeColor','k');
end;%for nl=0:numel(index_sort_u_label_B_3mo_)-1;
axis equal; axisnotick;
legend(u_label_B_3mo_(1+index_sort_u_label_B_3mo_));
title('pca');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% cluster full data-set. ;
%%%%%%%%;
flag_calculate = strcmp(platform,'access1');
if flag_calculate;
test_dexcluster_driver_5( ...
 str_code ...
,B_3mo_logn_sg__ ...
,index_nu_B_3mo_from_nall_ ...
,dir_mat ...
,rank_estimate_B_3mo_logn ...
);
end;%if flag_calculate;

disp('returning'); return;

%%%%%%%%;
% create label which only refers to pasca, ignoring Isogenic_type. ;
%%%%%%%%;
label_B_3mo_pasca_ = label_pasca_(1+index_stage_3mo_);
[ ...
 label_B_3mo_pasca_enum_ ...
,n_u_label_B_3mo_pasca ...
,u_label_B_3mo_pasca_ ...
,index_nu_B_3mo_pasca_from_nall_ ...
,n_u_label_B_3mo_pasca_ ...
,index_nall_from_nu_B_3mo_pasca__ ...
] = ...
label_str_to_enum_1( ...
 label_B_3mo_pasca_ ...
);
%%%%%%%%;
[~,index_sort_u_label_B_3mo_pasca_] = sort(n_u_label_B_3mo_pasca_,'descend'); index_sort_u_label_B_3mo_pasca_ = index_sort_u_label_B_3mo_pasca_ - 1;
index_nu_srt_from_nu_ori_ = index_sort_u_label_B_3mo_pasca_;
n_all = sum(n_u_label_B_3mo_pasca_);
for nu_label_B_3mo_pasca=0:n_u_label_B_3mo_pasca-1;
tmp_n_u = n_u_label_B_3mo_pasca_(1+index_sort_u_label_B_3mo_pasca_(1+nu_label_B_3mo_pasca));
u_label_B_3mo_pasca = u_label_B_3mo_pasca_{1+index_sort_u_label_B_3mo_pasca_(1+nu_label_B_3mo_pasca)};
disp(sprintf(' %% nu %.2d/%.2d: %.4d/%.4d = %0.2f = %s',nu_label_B_3mo_pasca,n_u_label_B_3mo_pasca,tmp_n_u,n_all,tmp_n_u/n_all,u_label_B_3mo_pasca));
end;%for nu_label_B_3mo_pasca=0:n_u_label_B_3mo_pasca-1;
%%%%%%%%;
% create label which only refers to Isogenic_type, ignoring pasca. ;
%%%%%%%%;
label_B_3mo_Isog_ = label_Isogenic_type_(1+index_stage_3mo_);
[ ...
 label_B_3mo_Isog_enum_ ...
,n_u_label_B_3mo_Isog ...
,u_label_B_3mo_Isog_ ...
,index_nu_B_3mo_Isog_from_nall_ ...
,n_u_label_B_3mo_Isog_ ...
,index_nall_from_nu_B_3mo_Isog__ ...
] = ...
label_str_to_enum_1( ...
 label_B_3mo_Isog_ ...
);
%%%%%%%%;
[~,index_sort_u_label_B_3mo_Isog_] = sort(n_u_label_B_3mo_Isog_,'descend'); index_sort_u_label_B_3mo_Isog_ = index_sort_u_label_B_3mo_Isog_ - 1;
index_nu_srt_from_nu_ori_ = index_sort_u_label_B_3mo_Isog_;
n_all = sum(n_u_label_B_3mo_Isog_);
for nu_label_B_3mo_Isog=0:n_u_label_B_3mo_Isog-1;
tmp_n_u = n_u_label_B_3mo_Isog_(1+index_sort_u_label_B_3mo_Isog_(1+nu_label_B_3mo_Isog));
u_label_B_3mo_Isog = u_label_B_3mo_Isog_{1+index_sort_u_label_B_3mo_Isog_(1+nu_label_B_3mo_Isog)};
disp(sprintf(' %% nu %.2d/%.2d: %.4d/%.4d = %0.2f = %s',nu_label_B_3mo_Isog,n_u_label_B_3mo_Isog,tmp_n_u,n_all,tmp_n_u/n_all,u_label_B_3mo_Isog));
end;%for nu_label_B_3mo_Isog=0:n_u_label_B_3mo_Isog-1;
%%%%%%%%;

%%%%%%%%;
% collect statistics from full data-set. ;
%%%%%%%%;
[ ...
 flag_found_ ...
,lpv_quad_ ...
,n_cluster_found_ ...
,label_B__ ...
,prefix_method_ ...
,legend_method_ ...
,nrank_method_ ...
,symbol_method_ ...
] = ...
test_dexcluster_driver_collect_5( ...
 str_code ...
,index_nu_B_3mo_from_nall_ ...
,dir_mat ...
,rank_estimate_B_3mo_logn ...
);
%%%%%%%%;
[ ...
 ~ ...
,lpv_quad_pasca_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
test_dexcluster_driver_collect_5( ...
 str_code ...
,index_nu_B_3mo_pasca_from_nall_ ...
,dir_mat ...
,rank_estimate_B_3mo_logn ...
);
%%%%%%%%;
[ ...
 ~ ...
,lpv_quad_Isog_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
test_dexcluster_driver_collect_5( ...
 str_code ...
,index_nu_B_3mo_Isog_from_nall_ ...
,dir_mat ...
,rank_estimate_B_3mo_logn ...
);
%%%%%%%%;

if ~exist('USB_3mo_logn__','var'); [tmp_U_,tmp_S_,tmp_V_] = svds(B_3mo_logn_sg__,rank_estimate_B_3mo_logn); USB_3mo_logn__ = tmp_U_*tmp_S_; clear tmp_U_ tmp_S_ tmp_V_; end;

%%%%%%%%;
% loading hnbr results from full data-set. ;
%%%%%%%%;
infix_method_0 = 'hnbrtZRgumb_r21_r1'; infix_method_1 = sprintf('%s',infix_method_0); E_array_base_ = USB_3mo_logn__;
%infix_method_0 = 'hnbtZRgumb_r1'; infix_method_1 = sprintf('%s_g010',infix_method_0);
dir_mat_method = sprintf('%s/dir_tmp_%s/dir_tmp_%s',dir_mat,infix_method_0,infix_method_1);
fname_output_label = sprintf('%s/output_label__.txt',dir_mat_method);
fp = fopen(fname_output_label); output_label__ = textscan(fp,'%s','Delimiter','\n'); output_label__ = output_label__{1}; fclose(fp);
fname_nlpbra_label = sprintf('%s/nlpbra_label__.txt',dir_mat_method);
fp = fopen(fname_nlpbra_label); nlpbra_label__ = textscan(fp,'%s','Delimiter','\n'); nlpbra_label__ = nlpbra_label__{1}; fclose(fp);
fname_nlpnex_label = sprintf('%s/nlpnex_label__.txt',dir_mat_method);
fp = fopen(fname_nlpnex_label); nlpnex_label__ = textscan(fp,'%s','Delimiter','\n'); nlpnex_label__ = nlpnex_label__{1}; fclose(fp);
[lpv,lP0,fla,cap_,cup_] = label_to_label_enrichment_quad_4(index_nu_B_3mo_from_nall_,label_str_to_enum_0(output_label__));
%%%%%%%%;
label_tree ...
= ...
label_to_tree_1( ...
 0 ...
,[] ...
,[] ...
,[] ...
,output_label__ ...
,nlpbra_label__ ...
,nlpnex_label__ ...
);
%%%%%%%%%%%%%%%%;
nlpnex_threshold_ = 3:0.5:27; n_nlpnex_threshold = numel(nlpnex_threshold_);
ouput_label_pid___ = cell(n_nlpnex_threshold,1);
lpv_method_0_ = zeros(n_nlpnex_threshold,1);
lpv_method_0_pasca_ = zeros(n_nlpnex_threshold,1);
lpv_method_0_Isog_ = zeros(n_nlpnex_threshold,1);
n_cluster_method_0_ = zeros(n_nlpnex_threshold,1);
lpv_method_1_ = zeros(n_nlpnex_threshold,1);
lpv_method_1_pasca_ = zeros(n_nlpnex_threshold,1);
lpv_method_1_Isog_ = zeros(n_nlpnex_threshold,1);
n_cluster_method_1_ = zeros(n_nlpnex_threshold,1);
for nnlpnex_threshold=0:n_nlpnex_threshold-1;
nlpnex_threshold = nlpnex_threshold_(1+nnlpnex_threshold);
if (verbose); disp(sprintf(' %% nnlpnex_threshold %d/%d nlpnex_threshold %0.2f',nnlpnex_threshold,n_nlpnex_threshold,nlpnex_threshold)); end;
tmp_output_label__ = label_tree_output_label_from_nlpnex_0(label_tree,nlpnex_threshold,[]);
output_label_pid___{1+nnlpnex_threshold} = tmp_output_label__;
label_B_0_ = label_str_to_enum_1(tmp_output_label__);
tmp_lpv = label_to_label_enrichment_quad_4(index_nu_B_3mo_from_nall_,label_B_0_);
tmp_lpv_pasca = label_to_label_enrichment_quad_4(index_nu_B_3mo_pasca_from_nall_,label_B_0_);
tmp_lpv_Isog = label_to_label_enrichment_quad_4(index_nu_B_3mo_Isog_from_nall_,label_B_0_);
lpv_method_0_(1+nnlpnex_threshold) = tmp_lpv;
lpv_method_0_pasca_(1+nnlpnex_threshold) = tmp_lpv_pasca;
lpv_method_0_Isog_(1+nnlpnex_threshold) = tmp_lpv_Isog;
n_cluster_method_0_(1+nnlpnex_threshold) = numel(unique(label_B_0_));
label_B_1_ = label_rearrange_0(label_B_0_,E_array_base_);
tmp_lpv = label_to_label_enrichment_quad_4(index_nu_B_3mo_from_nall_,label_B_1_);
tmp_lpv_pasca = label_to_label_enrichment_quad_4(index_nu_B_3mo_pasca_from_nall_,label_B_1_);
tmp_lpv_Isog = label_to_label_enrichment_quad_4(index_nu_B_3mo_Isog_from_nall_,label_B_1_);
lpv_method_1_(1+nnlpnex_threshold) = tmp_lpv;
lpv_method_1_pasca_(1+nnlpnex_threshold) = tmp_lpv_pasca;
lpv_method_1_Isog_(1+nnlpnex_threshold) = tmp_lpv_Isog;
n_cluster_method_1_(1+nnlpnex_threshold) = numel(unique(label_B_1_));
end;%for nnlpnex_threshold=0:n_nlpnex_threshold-1;
%%%%%%%%%%%%%%%%;

%%%%%%%%;
fname_fig = sprintf('%s/lpv_nlpnex_threshold_%s_FIGA',dir_jpg,infix_method_1);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figmed;
index_louvain00pr = efind(strcmp(prefix_method_,'louvain00pr_default'));
index_umap00pr = efind(strcmp(prefix_method_,'umap00pr_default'));
index_hnbrtZRgumb = efind(strcmp(prefix_method_,'hnbrtZRgumb'));
%%%%%%%%;
subplot(1,3,1); hold on;
%plot(nlpnex_threshold_,-lpv_method_0_,'mo-');
plot(nlpnex_threshold_,-lpv_method_1_,'ms-');
plot(nlpnex_threshold_,-lpv_quad_(1+index_louvain00pr)*ones(n_nlpnex_threshold,1),'rx-');
plot(nlpnex_threshold_,-lpv_quad_(1+index_umap00pr)*ones(n_nlpnex_threshold,1),'k^-');
plot(nlpnex_threshold_,-lpv_quad_(1+index_hnbrtZRgumb)*ones(n_nlpnex_threshold,1),'m-');
xlabel('nlp'); ylabel('entropy');
title('pasca_Isog','Interpreter','none');
%%%%%%%%;
subplot(1,3,2); hold on;
%plot(nlpnex_threshold_,-lpv_method_0_pasca_,'mo-');
plot(nlpnex_threshold_,-lpv_method_1_pasca_,'ms-');
plot(nlpnex_threshold_,-lpv_quad_pasca_(1+index_louvain00pr)*ones(n_nlpnex_threshold,1),'rx-');
plot(nlpnex_threshold_,-lpv_quad_pasca_(1+index_umap00pr)*ones(n_nlpnex_threshold,1),'k^-');
plot(nlpnex_threshold_,-lpv_quad_pasca_(1+index_hnbrtZRgumb)*ones(n_nlpnex_threshold,1),'m-');
xlabel('nlp'); ylabel('entropy');
title('pasca','Interpreter','none');
%%%%%%%%;
subplot(1,3,3); hold on;
%plot(nlpnex_threshold_,-lpv_method_0_Isog_,'mo-');
plot(nlpnex_threshold_,-lpv_method_1_Isog_,'ms-');
plot(nlpnex_threshold_,-lpv_quad_Isog_(1+index_louvain00pr)*ones(n_nlpnex_threshold,1),'rx-');
plot(nlpnex_threshold_,-lpv_quad_Isog_(1+index_umap00pr)*ones(n_nlpnex_threshold,1),'k^-');
plot(nlpnex_threshold_,-lpv_quad_Isog_(1+index_hnbrtZRgumb)*ones(n_nlpnex_threshold,1),'m-');
xlabel('nlp'); ylabel('entropy');
title('Isog','Interpreter','none');
%%%%%%%%;
sgtitle(sprintf('%s',dir_mat_method),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
fname_fig = sprintf('%s/label_plot_recursive_nlpscale_%s_FIGA',dir_jpg,infix_method_1);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(1);figbig;;
c_nlpvt__ = colormap_nlpvt(); n_c_nlpvt = size(c_nlpvt__,1); nlp_lim_ = [0,27];
subplot(1,9,1);
cla;
tmp_x__ = zeros(4,n_nlpnex_threshold);
tmp_y__ = zeros(4,n_nlpnex_threshold);
tmp_c___ = zeros(1,n_nlpnex_threshold,3);
for nnlpnex_threshold=0:n_nlpnex_threshold-1;
nlpnex_threshold = nlpnex_threshold_(1+nnlpnex_threshold);
nlpv = -lpv_method_(1+nnlpnex_threshold);
tmp_x__(:,1+nnlpnex_threshold) = [0;nlpv;nlpv;0];
tmp_y__(:,1+nnlpnex_threshold) = nlpnex_threshold + 0.5*[-1;-1;+1;+1];
nc_nlpvt = max(0,min(n_c_nlpvt-1,floor(n_c_nlpvt*(nlpnex_threshold-min(nlp_lim_))/diff(nlp_lim_))));
tmp_c___(1,1+nnlpnex_threshold,:) = c_nlpvt__(1+nc_nlpvt,:);
end;%for nnlpnex_threshold=0:n_nlpnex_threshold-1;
hold on;
patch(tmp_x__,tmp_y__,tmp_c___);
hold off;
ylim([min(nlp_lim_),max(nlp_lim_)+1]);
set(gca,'Ytick',0:max(nlp_lim_)+1);
xlim([0,max(-lpv_method_)*1.25]);
xlabel('entropy');
ylabel('nlp');
grid on;
%%%%%%%%;
subplot(1,9,[2:9]);
label_plot_recursive_nlpscale_2(output_label__,nlpbra_label__,nlpnex_label__,[],[],[]);
sgtitle(sprintf('%s',dir_mat_method),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;


