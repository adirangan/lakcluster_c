%%%%%%%%;
% Now try and plot the base auc (i.e., all patients) for the various prs-scores. ;
%%%%%%%%;
clear;
dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PGC_20190328');
dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c_dev');
%%%%%%%%;
lisa_arm1_r0 = lisa_struct_make_ver0('','',4,'dex',0.004,34,2,0,5,2,0,0); % arguments: mr_string,mc_string,cl_num,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
lisa_arm1_r0 = lisa_struct_prefix_ver0(lisa_arm1_r0,dir_code,dir_trunk); 
lisa_arm1_r0.nshuffle = 0;  lisa_arm1_r0 = lisa_struct_names_ver0(lisa_arm1_r0); 
lisa_arm1_r0 = lisa_struct_xdrop_ver0(lisa_arm1_r0); lisa_arm1_r0 = lisa_struct_mdsfam_ver0(lisa_arm1_r0); 
%disp(sprintf(' %% loading bim for %s',lisa_arm1_r0.string_name_s0));
lisa_arm1_r0 = lisa_struct_bim_ver0(lisa_arm1_r0); %<-- this is large and takes a while to load. ;
lisa_arm1_r0 = lisa_struct_mx_ver0(lisa_arm1_r0); lisa_arm1_r0 = lisa_struct_studyindex_ver0(lisa_arm1_r0); 
lisa_arm1_r0 = lisa_struct_trace_ver0(lisa_arm1_r0);
%%%%%%%%;
cl_num_arm1 = 4; cl_num_arm2 = 1; 
%%%%%%%%;
prs_prefix = 'ori'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 175; ni_r1_arm2 = 163;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver3_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni000_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'keep_bc1_ni175'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 175; ni_r1_arm2 = 163;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver3_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni175_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'keep_bc1_ni375'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 175; ni_r1_arm2 = 163;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver3_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni375_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'keep_bc2_ni163'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 175; ni_r1_arm2 = 163;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver3_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni163_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'tree_plt4_bicl1_ni188'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 238; ni_r1_arm2 = 263;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver5_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni188_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'tree_plt4_bicl1_ni200'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 238; ni_r1_arm2 = 263;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver5_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni200_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'tree_plt4_bicl1_ni213'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 238; ni_r1_arm2 = 263;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver5_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni213_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'tree_plt4_bicl1_ni225'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 238; ni_r1_arm2 = 263;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver5_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni225_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'tree_plt4_bicl1_ni238'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 238; ni_r1_arm2 = 263;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver5_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni238_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'tree_plt4_bicl1_ni250'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 238; ni_r1_arm2 = 263;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver5_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni250_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'tree_plt4_bicl1_ni263'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 238; ni_r1_arm2 = 263;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver5_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni263_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;
prs_prefix = 'tree_plt4_bicl1_ni275'; cl_num_arm1_string = ''; cl_num_arm2_string = ''; ni_r0_arm2 = 238; ni_r1_arm2 = 263;
lisa_arm1_r0.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1_r0.dir_out_s0,lisa_arm1_r0.cl_num,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string);
if (~exist(lisa_arm1_r0.dir_out_s0_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_prs_scatterplot_ver3',lisa_arm1_r0.dir_out_s0_pca)); end;
dir_out_s0_pca_r0_arm1 = sprintf('%s/dir_prs_mat',lisa_arm1_r0.dir_out_s0_pca);
fname_mat = sprintf('%s/prs_scatterplot_ver5_%s_trn%d%s_tst%d%s_ni%.3d_ni%.3d.mat',dir_out_s0_pca_r0_arm1,prs_prefix,cl_num_arm1,cl_num_arm1_string,cl_num_arm2,cl_num_arm2_string,ni_r0_arm2,ni_r1_arm2);
tmp_auc19_prs_t_ = load(fname_mat,'auc19_prs_t_');
auc19_prs_t_ni275_ = squeeze(tmp_auc19_prs_t_.auc19_prs_t_(1,1,:));
%%%%%%%%;

n_xxx = 10;
auc19_prs_t__ = cell(n_xxx,1);
str_prs_t_ = cell(n_xxx,1);
n_xxx = 1; 
auc19_prs_t__{n_xxx} = auc19_prs_t_ni163_; str_prs_t_{n_xxx} = '163'; n_xxx = n_xxx + 1;
auc19_prs_t__{n_xxx} = auc19_prs_t_ni175_; str_prs_t_{n_xxx} = '175'; n_xxx = n_xxx + 1;
auc19_prs_t__{n_xxx} = auc19_prs_t_ni188_; str_prs_t_{n_xxx} = '188'; n_xxx = n_xxx + 1;
auc19_prs_t__{n_xxx} = auc19_prs_t_ni200_; str_prs_t_{n_xxx} = '200'; n_xxx = n_xxx + 1;
auc19_prs_t__{n_xxx} = auc19_prs_t_ni213_; str_prs_t_{n_xxx} = '213'; n_xxx = n_xxx + 1;
auc19_prs_t__{n_xxx} = auc19_prs_t_ni225_; str_prs_t_{n_xxx} = '225'; n_xxx = n_xxx + 1;
auc19_prs_t__{n_xxx} = auc19_prs_t_ni238_; str_prs_t_{n_xxx} = '238'; n_xxx = n_xxx + 1;
auc19_prs_t__{n_xxx} = auc19_prs_t_ni250_; str_prs_t_{n_xxx} = '250'; n_xxx = n_xxx + 1;
auc19_prs_t__{n_xxx} = auc19_prs_t_ni263_; str_prs_t_{n_xxx} = '263'; n_xxx = n_xxx + 1;
auc19_prs_t__{n_xxx} = auc19_prs_t_ni275_; str_prs_t_{n_xxx} = '275'; n_xxx = n_xxx + 1;
n_xxx = n_xxx-1;

%%%%%%%%;
% These prs thresholds apply to original prs score only, as well as last 16 indices for other prs scores. ;
%%%%%%%%;
str_prs_ = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
ij_ori_ = 1:16;
ij_xxx_ = 5:20;
c_ = colormap('cool'); n_c = size(c_,1); clim_ = [163,275];
for nxxx=1:n_xxx;
subplot(2,5,nxxx);
hold on;
nc = max(1,min(n_c,floor(n_c*nxxx/n_xxx))); 
plot(1:16,auc19_prs_t_ni000_(ij_ori_),'.-','Color',[0,0,0],'LineWidth',2,'MarkerSize',25);
plot(1:16,auc19_prs_t__{nxxx}(ij_xxx_),'.-','Color',c_(nc,:),'LineWidth',2,'MarkerSize',25);
plot(1:16,0.5*ones(1,16),'k:');
hold off;
ylim([0.48,0.58]); ylabel('cauc');
xlim([1,16]);
set(gca,'XTick',1:16,'XTickLabel',str_prs_cutoff_); xtickangle(90); xlabel('prs threshold');
legend({'000',str_prs_t_{nxxx}},'Location','NorthWest');
title(sprintf('prs ni%s',str_prs_t_{nxxx}));
end;%for nxxx=1:n_xxx;
set(gcf,'Position',1+[0,0,512*4,512*2]);
fname_base = sprintf('%s/dir_misc/prs_comparison_trn%d_tst%d_C',dir_trunk,cl_num_arm1,cl_num_arm2);
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));


%{
%%
%
%%
%%%%%%%%;
% These prs thresholds apply to original prs score only, as well as last 16 indices for other prs scores. ;
%%%%%%%%;
str_prs_ = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'}; 
str_prs_cutoff_ = {'1e-7.5','1e-7','1e-6.5','1e-6','1e-5.5','1e-5','1e-4.5','1e-4','1e-3.5','1e-3','1e-2.5','1e-2','1e-1.5','1e-1','1e-0.5','1e0'}; 
ij_ori_ = 1:16;
ij_xxx_ = 5:20;
c_ = colormap(colormap_beach()); n_c = size(c_,1); clim_ = [163,275];
subplot(1,1,1);
hold on;
%ni =   0; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni000_(ij_ori_),'o-','Color',c_(nc,:));
ni =   0; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni000_(ij_ori_),'.-','Color',[0,0,0]);
ni = 163; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni163_(ij_xxx_),'.-','Color',c_(nc,:));
ni = 175; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni175_(ij_xxx_),'.-','Color',c_(nc,:));
ni = 188; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni188_(ij_xxx_),'.-','Color',c_(nc,:));
ni = 200; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni200_(ij_xxx_),'.-','Color',c_(nc,:));
ni = 213; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni213_(ij_xxx_),'.-','Color',c_(nc,:));
ni = 225; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni225_(ij_xxx_),'.-','Color',c_(nc,:));
ni = 238; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni238_(ij_xxx_),'.-','Color',c_(nc,:));
ni = 250; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni250_(ij_xxx_),'.-','Color',c_(nc,:));
ni = 263; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni263_(ij_xxx_),'.-','Color',c_(nc,:));
ni = 275; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni275_(ij_xxx_),'.-','Color',c_(nc,:));
%ni = 375; nc = max(1,min(n_c,floor(n_c*(ni-min(clim_))/diff(clim_)))); plot(1:16,auc19_prs_t_ni375_(ij_xxx_),'o-','Color',c_(nc,:));
plot(1:16,0.5*ones(1,16),'k:');
hold off;
ylim([0.48,0.58]); ylabel('cauc');
xlim([1,16]);
set(gca,'XTick',1:16,'XTickLabel',str_prs_cutoff_); xtickangle(90); xlabel('prs threshold');
%legend({'000','163','175','188','200','213','225','238','250','263','275','375'},'Location','NorthWest');
legend({'000','163','175','188','200','213','225','238','250','263','275'},'Location','NorthWest');
title(sprintf('prs cauc trn%d tst%d',cl_num_arm1,cl_num_arm2));
 %}


%{
clf;
%%%%%%%%;
%%%%%%%%;
%Sn3 	10^(-9.5) 	3.16E-10 	0 	0
%Sn2 	10^(-9) 	1E-09 		0 	1
%Sn1 	10^(-8.5) 	3.16E-09 	0 	1
%S0 	10^(-8) 	1E-08 		0 	1
%%%%%%%%;
%S1 	10^(-7.5) 	3.16E-08 	0 	1
%S2 	10^(-7) 	1E-07 		0 	1
%S3 	10^(-6.5) 	3.16E-07 	2 	3
%S4 	10^(-6) 	0.000001 	3 	4
%S5 	10^(-5.5) 	3.16E-06 	13 	8
%S6 	10^(-5) 	0.00001 	23 	15
%S7 	10^(-4.5) 	3.16E-05 	64 	38
%S8 	10^(-4) 	0.0001 		175 	107
%S9 	10^(-3.5) 	0.000316 	460 	314
%S10 	10^(-3) 	0.001 		1224 	810
%S11 	10^(-2.5) 	0.00316 	3295 	2151
%S12 	10^(-2) 	0.01 		8707 	5298
%S13 	10^(-1.5) 	0.0316 		21723 	12659
%S14 	10^(-1) 	0.1 		52636 	29157
%S15 	10^(-0.5) 	0.316 		119024 	61337
%S16 	10^(-0) 	1 		227175 	102286
%%%%%%%%;
subplot(1,1,1);
n_snp_A_ = [0,0,2,3,13,23,64,175,460,1224,3295,8707,21723,52636,119024,227175];
n_snp_B_ = [1,1,3,4, 8,15,38,107,314, 810,2151,5298,12659,29157, 61337,102286];
hold on;
plot(log10(n_snp_A_),auc19_prs_t_A_(ij_ori_),'ko-');
plot(log10(n_snp_B_),auc19_prs_t_B_(ij_xxx_),'ro-');
plot([0,log10(max(n_snp_A_))],0.5*ones(1,2),'k:');
hold off;
ylim([0.48,0.58]); ylabel('cauc');
xlim([0,log10(max(n_snp_A_))]);
xlabel('log10(number of snps)');
legend('ori','bc1 175','Location','NorthWest');
title(sprintf('prs cauc trn%d tst%d',cl_num_arm1,cl_num_arm2));
%%%%%%%%;
set(gcf,'Position',1+[0,0,512*1,512]);
fname_base = sprintf('%s/dir_misc/prs_comparison_trn%d_tst%d_C',dir_trunk,cl_num_arm1,cl_num_arm2);
print('-djpeg',sprintf('%s.jpg',fname_base));
print('-depsc',sprintf('%s.eps',fname_base));
%%%%%%%%;
 %}
