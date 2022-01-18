function dolphin_dXX_9(dXX,n_shuffle,age_percentile_range_,sex_constrain,isU_constrain,nvar_recalc_);
% Note that age-range is 0.64 --> 54.35. ;
% try: ;
%{

[~,M_10_p] = min(abs(sort(age_ori_M_) - 10)); M_10_p = (M_10_p-0)/numel(age_ori_M_);
[~,M_30_p] = min(abs(sort(age_ori_M_) - 30)); M_30_p = (M_30_p-0)/numel(age_ori_M_);
[~,F_10_p] = min(abs(sort(age_ori_F_) - 10)); F_10_p = (F_10_p-0)/numel(age_ori_F_);
[~,F_30_p] = min(abs(sort(age_ori_F_) - 30)); F_30_p = (F_30_p-0)/numel(age_ori_F_);
disp(sprintf(' %% age_ori_M_ prctile %0.2f --> %0.2f',100*M_10_p,prctile(age_ori_M_,100*M_10_p)));
disp(sprintf(' %% age_ori_M_ prctile %0.2f --> %0.2f',100*M_30_p,prctile(age_ori_M_,100*M_30_p)));
disp(sprintf(' %% age_ori_F_ prctile %0.2f --> %0.2f',100*F_10_p,prctile(age_ori_F_,100*F_10_p)));
disp(sprintf(' %% age_ori_F_ prctile %0.2f --> %0.2f',100*F_30_p,prctile(age_ori_F_,100*F_30_p)));

%}

if (nargin<1);
nvar_recalc_ = [21,43]; %<-- Creatinine and GFR. ;
for nU=2:-1:0;
if nU==0; isU_constrain = ''; end;
if nU==1; isU_constrain = 'lo'; end;
if nU==2; isU_constrain = 'hi'; end;
%%%%;
dXX = 0; n_shuffle = 256; sex_constrain = '';
dolphin_dXX_9(dXX,n_shuffle,[ 0, 20],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[20, 70],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[70,100],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[      ],sex_constrain,isU_constrain,nvar_recalc_);
%%%%;
dXX = 0; n_shuffle = 256; sex_constrain = '';
dolphin_dXX_9(dXX,n_shuffle,[ 0, 20],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[20, 70],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[70,100],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[      ],sex_constrain,isU_constrain,nvar_recalc_);
%%%%;
dXX = 0; n_shuffle = 256; sex_constrain = 'M';
dolphin_dXX_9(dXX,n_shuffle,[ 0, 21],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[21, 73],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[73,100],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[      ],sex_constrain,isU_constrain,nvar_recalc_);
%%%%;
dXX = 0; n_shuffle = 256; sex_constrain = 'M';
dolphin_dXX_9(dXX,n_shuffle,[ 0, 21],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[21, 73],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[73,100],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[      ],sex_constrain,isU_constrain,nvar_recalc_);
%%%%;
dXX = 0; n_shuffle = 256; sex_constrain = 'F';
dolphin_dXX_9(dXX,n_shuffle,[ 0, 14],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[14, 68],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[68,100],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[      ],sex_constrain,isU_constrain,nvar_recalc_);
%%%%;
dXX = 0; n_shuffle = 256; sex_constrain = 'F';
dolphin_dXX_9(dXX,n_shuffle,[ 0, 14],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[14, 68],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[68,100],sex_constrain,isU_constrain,nvar_recalc_);
dolphin_dXX_9(dXX,n_shuffle,[      ],sex_constrain,isU_constrain,nvar_recalc_);
%%%%;
end;%for nU=0:2;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); dXX = []; end; na=na+1;
if (nargin<1+na); n_shuffle = []; end; na=na+1;
if (nargin<1+na); age_percentile_range_ = []; end; na=na+1;
if (nargin<1+na); sex_constrain = []; end; na=na+1;
if (nargin<1+na); isU_constrain = []; end; na=na+1;
if (nargin<1+na); nvar_recalc_ = []; end; na=na+1;

% test with: ;
% dXX = 0; n_shuffle = 256; age_percentile_range_ = []; sex_constrain = ''; isU_constrain = 'hi';

if (isempty(dXX)); dXX = 0; end;
if (isempty(n_shuffle)); n_shuffle = 32; end;
if (isempty(sex_constrain)); sex_constrain = ''; end;
if (isempty(isU_constrain)); isU_constrain = ''; end;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_dolphin',string_root);
dir_jpg = sprintf('%s/dir_jpg_ind_7',dir_trunk);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
dir_mat = sprintf('%s/dir_mat_ind_7',dir_trunk);
if (~exist(dir_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_mat)); mkdir(dir_mat); end;
dir_shuffle_mat = sprintf('%s/dir_shuffle%d_mat_ind_7',dir_trunk,n_shuffle);
if (~exist(dir_shuffle_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_shuffle_mat)); mkdir(dir_shuffle_mat); end;

tolerance_master = 1e-6;
flag_replot=0;
date_diff_threshold = 0.5;
flag_force_create_mat = 0;
flag_force_create_tmp = 0;

%%%%;
str_sex_constrain = '';
if ~isempty(sex_constrain);
str_sex_constrain = sprintf('_%s',sex_constrain);
end;%if ~isempty(sex_constrain);
%%%%;
str_isU_constrain = '';
if ~isempty(isU_constrain);
str_isU_constrain = sprintf('_U%s',isU_constrain);
end;%if ~isempty(isU_constrain);
%%%%;
str_age_constrain = '';
if (~isempty(age_percentile_range_));
str_age_constrain = sprintf('_age_%.2d%.2d',min(age_percentile_range_),min(99,max(age_percentile_range_)));
end;%if (~isempty(age_percentile_range_));
%%%%;
fname_infix = ...
sprintf( ...
 'aid%.2d%s%s%s' ...
,dXX ...
,str_age_constrain ...
,str_sex_constrain ...
,str_isU_constrain ...
);
disp(sprintf(' %% fname_infix %s',fname_infix));

dolphin_dXX_load_7;
string_dat_name_ = string_dat_ori_name_;
n_var = n_var_ori;
isM_ = isM_ori_;
isF_ = isF_ori_;
isU_bbar_lo_ = isU_bbar_lo_ori_;
isU_bbar_hi_ = isU_bbar_hi_ori_;
aid_ = aid_ori_;
age_ = age_ori_;
dat__ = res_ind_nrm__;
tmp_index_use_ = efind( 1 ...
&isfinite(age_) ...
& (sum(isfinite(dat__),2)> 0) ...
);
isM_ = isM_(1+tmp_index_use_);
isF_ = isF_(1+tmp_index_use_);
isU_bbar_lo_ = isU_bbar_lo_(1+tmp_index_use_);
isU_bbar_hi_ = isU_bbar_hi_(1+tmp_index_use_);
aid_ = aid_(1+tmp_index_use_);
age_ = age_(1+tmp_index_use_);
dat__ = dat__(1+tmp_index_use_,:);

%%%%%%%%;
% select sex. ;
%%%%%%%%;
if ~isempty(sex_constrain);
if strcmp(sex_constrain,'M'); tmp_index_ = efind(isM_); end;%if strcmp(sex_constrain,'M');
if strcmp(sex_constrain,'F'); tmp_index_ = efind(isF_); end;%if strcmp(sex_constrain,'F');
isM_ = isM_(1+tmp_index_);
isF_ = isF_(1+tmp_index_);
isU_bbar_lo_ = isU_bbar_lo_(1+tmp_index_);
isU_bbar_hi_ = isU_bbar_hi_(1+tmp_index_);
aid_ = aid_(1+tmp_index_);
age_ = age_(1+tmp_index_);
dat__ = dat__(1+tmp_index_,:);
end;%if ~isempty(sex_constrain);

%%%%%%%%;
% select isU. ;
%%%%%%%%;
if ~isempty(isU_constrain);
if strcmp(isU_constrain,'lo'); tmp_index_ = efind(isU_bbar_lo_); end;%if strcmp(isU_constrain,'lo');
if strcmp(isU_constrain,'hi'); tmp_index_ = efind(isU_bbar_hi_); end;%if strcmp(isU_constrain,'hi');
isM_ = isM_(1+tmp_index_);
isF_ = isF_(1+tmp_index_);
isU_bbar_lo_ = isU_bbar_lo_(1+tmp_index_);
isU_bbar_hi_ = isU_bbar_hi_(1+tmp_index_);
aid_ = aid_(1+tmp_index_);
age_ = age_(1+tmp_index_);
dat__ = dat__(1+tmp_index_,:);
end;%if ~isempty(isU_constrain);

%%%%%%%%;
% select age range. ;
%%%%%%%%;
if (~isempty(age_percentile_range_));
tmp_index_ = efind(age_>=prctile(age_,min(age_percentile_range_)) & age_<=prctile(age_,max(age_percentile_range_)));
disp(sprintf(' %% age_percentile_range_ [%0.2f %0.2f] --> [%0.2f %0.2f] <-- %d found',age_percentile_range_,prctile(age_,age_percentile_range_),numel(tmp_index_)));
isM_ = isM_(1+tmp_index_);
isF_ = isF_(1+tmp_index_);
isU_bbar_lo_ = isU_bbar_lo_(1+tmp_index_);
isU_bbar_hi_ = isU_bbar_hi_(1+tmp_index_);
aid_ = aid_(1+tmp_index_);
age_ = age_(1+tmp_index_);
dat__ = dat__(1+tmp_index_,:);
end;%if (~isempty(age_percentile_range_));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% processing full data-matrix. ;
%%%%%%%%;
fname_pre = sprintf('%s/dolphin_%s_s0000',dir_mat,fname_infix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~flag_skip;
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
%%%%%%%%;
[ ...
 parameter ...
,a_ ...
,A__ ...
,BB_inv__ ...
,CC_inv__ ...
,L ...
,niteration ...
,a__ ...
,A___ ...
,BB_inv___ ...
,CC_inv___ ...
,L_ ...
,index_min ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ ...
,age_ ...
,dat__ ...
);
save(fname_mat ...
,'parameter' ...
,'a_' ...
,'A__' ...
,'BB_inv__' ...
,'CC_inv__' ...
,'L' ...
,'niteration' ...
,'a__' ...
,'A___' ...
,'BB_inv___' ...
,'CC_inv___' ...
,'L_' ...
,'index_min' ...
);
%%%%%%%%;
close_fname_tmp(fname_pre);
end;%if ~flag_skip;

fname_mat = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, loading',fname_mat));
load(fname_mat);
index_mss_ = efind(~isfinite(dat__));
dat_imp__ = dolphin_impute_aA_0(aid_,age_,dat__,a_,A__,index_mss_);
%%%%%%%%;
% plot correlations. ;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_correlation_0',dir_jpg,fname_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1); clf; figbeach;
imagesc(corr(dat_imp__),[-1,+1]);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',12);
colorbar();
axis image;
figbig;
sgtitle('correlation matrix (after imputation)','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% plot SDE-regression results. ;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_aABC_4',dir_jpg,fname_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
subplot(2,2,1);
bar(a_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
title('drift a_','Interpreter','none');
subplot(2,2,2);
colormap(colormap_beach());
imagesc(A__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('interaction A__','Interpreter','none');
subplot(2,2,3);
colormap(colormap_beach());
imagesc(BB_inv__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('SDE noise BB_inv__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(CC_inv__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('observation noise CC_inv__','Interpreter','none');
figbig;
sgtitle('coefficient matrices','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% Now repeat for several permuted trials. ;
%%%%%%%%;
fname_pre = sprintf('%s/dolphin_%s_sxxxx',dir_shuffle_mat,fname_infix);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~flag_skip;
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
%%%%%%%%;
a_prm__ = zeros(n_var,1+n_shuffle);
A_prm___ = zeros(n_var,n_var,1+n_shuffle);
BB_inv_prm___ = zeros(n_var,n_var,1+n_shuffle);
CC_inv_prm___ = zeros(n_var,n_var,1+n_shuffle);
L_prm_ = zeros(1+n_shuffle,1);
for nshuffle=0:n_shuffle;
if (nshuffle==0); dat_prm__ = dat__; end;
if (nshuffle> 0); dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle); end;
if (mod(nshuffle,32)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
%%%%%%%%;
[ ...
 parameter ...
,a_prm_ ...
,A_prm__ ...
,BB_inv_prm__ ...
,CC_inv_prm__ ...
,L_prm ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_prm__ ...
);
a_prm__(:,1+nshuffle) = a_prm_;
A_prm___(:,:,1+nshuffle) = A_prm__;
BB_inv_prm___(:,:,1+nshuffle) = BB_inv_prm__;
CC_inv_prm___(:,:,1+nshuffle) = CC_inv_prm__;
L_prm_(1+nshuffle) = L_prm;
end;%for nshuffle=0:n_shuffle;
save(fname_mat ...
,'parameter' ...
,'a_prm__' ...
,'A_prm___' ...
,'BB_inv_prm___' ...
,'CC_inv_prm___' ...
,'L_prm_' ...
);
close_fname_tmp(fname_pre);
end;% if ~flag_skip;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% processing variable pairs. ;
%%%%%%%%;
fname_pre = sprintf('%s/dolphin_%s_cmb_s0000',dir_mat,fname_infix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ( isempty(nvar_recalc_));
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~flag_skip;
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
%%%%%%%%;
a_cmb_ = zeros(n_var,1);
A_cmb__ = zeros(n_var,n_var);
BB_inv_cmb__ = zeros(n_var,n_var);
CC_inv_cmb__ = zeros(n_var,n_var);
L_cmb = 0;
a_each__ = zeros(n_var,n_var);
A_each__ = zeros(n_var,n_var);
BB_inv_each__ = zeros(n_var,n_var);
CC_inv_each__ = zeros(n_var,n_var);
L_each__ = zeros(n_var,n_var);
for nvar0=0:n_var-1;
disp(sprintf(' %% nvar0 %d/%d',nvar0,n_var));
for nvar1=nvar0+1:n_var-1;
dat_pair__ = dat__(:,1+[nvar0,nvar1]);
[ ...
 parameter ...
,a_pair_ ...
,A_pair__ ...
,BB_inv_pair__ ...
,CC_inv_pair__ ...
,L_pair ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_pair__ ...
);
a_cmb_(1+[nvar0;nvar1]) = a_cmb_(1+[nvar0;nvar1]) + a_pair_(:);
A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + A_pair__;
BB_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = BB_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + BB_inv_pair__;
CC_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = CC_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + CC_inv_pair__;
L_cmb = L_cmb + L_pair;
a_each__(1+nvar0,1+nvar1) = a_pair_(1+0); a_each__(1+nvar1,1+nvar0) = a_pair_(1+1);
A_each__(1+nvar0,1+nvar1) = A_pair__(1+0,1+0); A_each__(1+nvar1,1+nvar0) = A_pair__(1+1,1+1);
BB_inv_each__(1+nvar0,1+nvar1) = BB_inv_pair__(1+0,1+0); BB_inv_each__(1+nvar1,1+nvar0) = BB_inv_pair__(1+1,1+1);
CC_inv_each__(1+nvar0,1+nvar1) = CC_inv_pair__(1+0,1+0); CC_inv_each__(1+nvar1,1+nvar0) = CC_inv_pair__(1+1,1+1);
L_each__(1+nvar0,1+nvar1) = L_pair; L_each__(1+nvar1,1+nvar0) = L_pair;
end;%for nvar1=nvar0:n_var-1;
end;%for nvar0=0:n_var-1;
save(fname_mat ...
,'parameter' ...
,'a_cmb_' ...
,'A_cmb__' ...
,'BB_inv_cmb__' ...
,'CC_inv_cmb__' ...
,'L_cmb' ...
,'a_each__' ...
,'A_each__' ...
,'BB_inv_each__' ...
,'CC_inv_each__' ...
,'L_each__' ...
);
%%%%%%%%;
close_fname_tmp(fname_pre);
end;% if ~flag_skip;
end;%if ( isempty(nvar_recalc_));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (~isempty(nvar_recalc_));
fname_recalc_pre = sprintf('%s/dolphin_recalc_%s_cmb_s0000',dir_mat,fname_infix);
[flag_recalc_skip,fname_recalc_mat] = open_fname_tmp(fname_recalc_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~flag_recalc_skip;
nvar_recalc_ = sort(unique(nvar_recalc_));
n_nvar_recalc = numel(nvar_recalc_);
fname_mat = sprintf('%s.mat',fname_pre);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, updating',fname_mat));
load(fname_mat);
parameter.tolerance_master = tolerance_master;
%%%%;
for nnvar_recalc0=0:n_nvar_recalc-1;
for nnvar_recalc1=0:n_var-1;
nvar0 = min(nvar_recalc_(1+nnvar_recalc0),nnvar_recalc1);
nvar1 = max(nvar_recalc_(1+nnvar_recalc0),nnvar_recalc1);
if (nvar0<nvar1);
disp(sprintf(' %% nvar0 %d/%d nvar1 %d/%d',nvar0,n_var,nvar1,n_var));
dat_pair__ = dat__(:,1+[nvar0,nvar1]);
[ ...
 parameter ...
,a_pair_ ...
,A_pair__ ...
,BB_inv_pair__ ...
,CC_inv_pair__ ...
,L_pair ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_pair__ ...
);
a_cmb_(1+[nvar0;nvar1]) = a_pair_(:);
A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = A_pair__;
BB_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = BB_inv_pair__;
CC_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = CC_inv_pair__;
L_cmb = L_pair;
a_each__(1+nvar0,1+nvar1) = a_pair_(1+0); a_each__(1+nvar1,1+nvar0) = a_pair_(1+1);
A_each__(1+nvar0,1+nvar1) = A_pair__(1+0,1+0); A_each__(1+nvar1,1+nvar0) = A_pair__(1+1,1+1);
BB_inv_each__(1+nvar0,1+nvar1) = BB_inv_pair__(1+0,1+0); BB_inv_each__(1+nvar1,1+nvar0) = BB_inv_pair__(1+1,1+1);
CC_inv_each__(1+nvar0,1+nvar1) = CC_inv_pair__(1+0,1+0); CC_inv_each__(1+nvar1,1+nvar0) = CC_inv_pair__(1+1,1+1);
L_each__(1+nvar0,1+nvar1) = L_pair; L_each__(1+nvar1,1+nvar0) = L_pair;
end;%if (nvar0<nvar1);
end;%for nnvar_recalc1=0:n_var-1;
end;%for nnvar_recalc0=0:n_nvar_recalc-1;
%%%%;
save(fname_mat ...
,'parameter' ...
,'a_cmb_' ...
,'A_cmb__' ...
,'BB_inv_cmb__' ...
,'CC_inv_cmb__' ...
,'L_cmb' ...
,'a_each__' ...
,'A_each__' ...
,'BB_inv_each__' ...
,'CC_inv_each__' ...
,'L_each__' ...
);
%%%%%%%%;
end;%if ( exist(fname_mat,'file'));
save(fname_recalc_mat,'fname_recalc_mat');
close_fname_tmp(fname_recalc_pre);
end;%if ~flag_recalc_skip;
end;%if (~isempty(nvar_recalc_));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

fname_mat = sprintf('%s/dolphin_%s_cmb_s0000.mat',dir_mat,fname_infix);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, loading',fname_mat));
load(fname_mat);
%%%%%%%%;
% plot SDE-regression results. ;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_cmb_aABC_4',dir_jpg,fname_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
subplot(2,2,1);
bar(a_cmb_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
title('drift a_cmb_','Interpreter','none');
subplot(2,2,2);
colormap(colormap_beach());
imagesc(A_cmb__ - diag(diag(A_cmb__))); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('interaction A_cmb__','Interpreter','none');
subplot(2,2,3);
colormap(colormap_beach());
imagesc(BB_inv_cmb__ - diag(diag(BB_inv_cmb__))); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('SDE noise BB_inv_cmb__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(CC_inv_cmb__ - diag(diag(CC_inv_cmb__))); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('observation noise CC_inv_cmb__','Interpreter','none');
figbig;
sgtitle('coefficient matrices','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% Now repeat for several permuted trials. ;
%%%%%%%%;
fname_pre = sprintf('%s/dolphin_%s_cmb_sxxxx',dir_shuffle_mat,fname_infix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ( isempty(nvar_recalc_));
[flag_skip,fname_mat] = open_fname_tmp(fname_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~flag_skip;
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
%%%%%%%%;
a_cmb_prm__ = zeros(n_var,1+n_shuffle);
A_cmb_prm___ = zeros(n_var,n_var,1+n_shuffle);
BB_inv_cmb_prm___ = zeros(n_var,n_var,1+n_shuffle);
CC_inv_cmb_prm___ = zeros(n_var,n_var,1+n_shuffle);
L_cmb_prm_ = zeros(1+n_shuffle,1);
a_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
A_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
BB_inv_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
CC_inv_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
L_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
for nshuffle=0:n_shuffle;
if (nshuffle==0); dat_prm__ = dat__; end;
if (nshuffle> 0); dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle); end;
if (mod(nshuffle,32)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
%%%%%%%%;
a_cmb_ = zeros(n_var,1);
A_cmb__ = zeros(n_var,n_var);
BB_inv_cmb__ = zeros(n_var,n_var);
CC_inv_cmb__ = zeros(n_var,n_var);
L_cmb = 0;
a_each__ = zeros(n_var,n_var);
A_each__ = zeros(n_var,n_var);
BB_inv_each__ = zeros(n_var,n_var);
CC_inv_each__ = zeros(n_var,n_var);
L_each__ = zeros(n_var,n_var);
for nvar0=0:n_var-1;
disp(sprintf(' %% nvar0 %d/%d',nvar0,n_var));
for nvar1=nvar0+1:n_var-1;
dat_pair__ = dat_prm__(:,1+[nvar0,nvar1]);
[ ...
 parameter ...
,a_pair_ ...
,A_pair__ ...
,BB_inv_pair__ ...
,CC_inv_pair__ ...
,L_pair ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_pair__ ...
);
a_cmb_(1+[nvar0;nvar1]) = a_cmb_(1+[nvar0;nvar1]) + a_pair_(:);
A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + A_pair__;
BB_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = BB_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + BB_inv_pair__;
CC_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = CC_inv_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + CC_inv_pair__;
L_cmb = L_cmb + L_pair;
a_each__(1+nvar0,1+nvar1) = a_pair_(1+0); a_each__(1+nvar1,1+nvar0) = a_pair_(1+1);
A_each__(1+nvar0,1+nvar1) = A_pair__(1+0,1+0); A_each__(1+nvar1,1+nvar0) = A_pair__(1+1,1+1);
BB_inv_each__(1+nvar0,1+nvar1) = BB_inv_pair__(1+0,1+0); BB_inv_each__(1+nvar1,1+nvar0) = BB_inv_pair__(1+1,1+1);
CC_inv_each__(1+nvar0,1+nvar1) = CC_inv_pair__(1+0,1+0); CC_inv_each__(1+nvar1,1+nvar0) = CC_inv_pair__(1+1,1+1);
L_each__(1+nvar0,1+nvar1) = L_pair; L_each__(1+nvar1,1+nvar0) = L_pair;
end;%for nvar1=nvar0:n_var-1;
end;%for nvar0=0:n_var-1;
a_cmb_prm__(:,1+nshuffle) = a_cmb_;
A_cmb_prm___(:,:,1+nshuffle) = A_cmb__;
BB_inv_cmb_prm___(:,:,1+nshuffle) = BB_inv_cmb__;
CC_inv_cmb_prm___(:,:,1+nshuffle) = CC_inv_cmb__;
L_cmb_prm_(1+nshuffle) = L_cmb;
a_each_prm___(:,:,1+nshuffle) = a_each__;
A_each_prm___(:,:,1+nshuffle) = A_each__;
BB_inv_each_prm___(:,:,1+nshuffle) = BB_inv_each__;
CC_inv_each_prm___(:,:,1+nshuffle) = CC_inv_each__;
L_each_prm___(:,:,1+nshuffle) = L_each__;
end;%for nshuffle=0:n_shuffle;
save(fname_mat ...
,'parameter' ...
,'a_cmb_prm__' ...
,'A_cmb_prm___' ...
,'BB_inv_cmb_prm___' ...
,'CC_inv_cmb_prm___' ...
,'L_cmb_prm_' ...
,'a_each_prm___' ...
,'A_each_prm___' ...
,'BB_inv_each_prm___' ...
,'CC_inv_each_prm___' ...
,'L_each_prm___' ...
);
close_fname_tmp(fname_pre);
end;% if ~flag_skip;
end;%if ( isempty(nvar_recalc_));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (~isempty(nvar_recalc_));
fname_recalc_pre = sprintf('%s/dolphin_recalc_%s_cmb_sxxxx',dir_mat,fname_infix);
[flag_recalc_skip,fname_recalc_mat] = open_fname_tmp(fname_recalc_pre,date_diff_threshold,flag_force_create_mat,flag_force_create_tmp);
if ~flag_recalc_skip;
nvar_recalc_ = sort(unique(nvar_recalc_));
n_nvar_recalc = numel(nvar_recalc_);
fname_mat = sprintf('%s.mat',fname_pre);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, updating',fname_mat));
load(fname_mat);
parameter.tolerance_master = tolerance_master;
%%%%%%%%;
for nshuffle=0:n_shuffle;
if (nshuffle==0); dat_prm__ = dat__; end;
if (nshuffle> 0); dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle); end;
if (mod(nshuffle,32)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
%%%%%%%%;
for nnvar_recalc0=0:n_nvar_recalc-1;
for nnvar_recalc1=0:n_var-1;
nvar0 = min(nvar_recalc_(1+nnvar_recalc0),nnvar_recalc1);
nvar1 = max(nvar_recalc_(1+nnvar_recalc0),nnvar_recalc1);
if (nvar0<nvar1);
disp(sprintf(' %% nvar0 %d/%d nvar1 %d/%d',nvar0,n_var,nvar1,n_var));
dat_pair__ = dat_prm__(:,1+[nvar0,nvar1]);
[ ...
 parameter ...
,a_pair_ ...
,A_pair__ ...
,BB_inv_pair__ ...
,CC_inv_pair__ ...
,L_pair ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ ...
,age_ ...
,dat_pair__ ...
);
A_cmb_prm___(1+[nvar0;nvar1],1+[nvar0;nvar1],1+nshuffle) = A_pair__;
BB_inv_cmb_prm___(1+[nvar0;nvar1],1+[nvar0;nvar1],1+nshuffle) = BB_inv_pair__;
CC_inv_cmb_prm___(1+[nvar0;nvar1],1+[nvar0;nvar1],1+nshuffle) = CC_inv_pair__;
a_each_prm___(1+nvar0,1+nvar1,1+nshuffle) = a_pair_(1+0);
a_each_prm___(1+nvar1,1+nvar0,1+nshuffle) = a_pair_(1+1);
A_each_prm___(1+nvar0,1+nvar1,1+nshuffle) = A_pair__(1+0,1+0);
A_each_prm___(1+nvar1,1+nvar0,1+nshuffle) = A_pair__(1+1,1+1);
BB_inv_each_prm___(1+nvar0,1+nvar1,1+nshuffle) = BB_inv_pair__(1+0,1+0);
BB_inv_each_prm___(1+nvar1,1+nvar0,1+nshuffle) = BB_inv_pair__(1+1,1+1);
CC_inv_each_prm___(1+nvar0,1+nvar1,1+nshuffle) = CC_inv_pair__(1+0,1+0);
CC_inv_each_prm___(1+nvar1,1+nvar0,1+nshuffle) = CC_inv_pair__(1+1,1+1);
L_each_prm__(1+nvar0,1+nvar1,1+nshuffle) = L_pair;
L_each_prm__(1+nvar1,1+nvar0,1+nshuffle) = L_pair;
end;%if (nvar0<nvar1);
end;%for nnvar_recalc1=0:n_var-1;
end;%for nnvar_recalc0=0:n_nvar_recalc-1;
end;%for nshuffle=0:n_shuffle;
save(fname_mat ...
,'parameter' ...
,'a_cmb_prm__' ...
,'A_cmb_prm___' ...
,'BB_inv_cmb_prm___' ...
,'CC_inv_cmb_prm___' ...
,'L_cmb_prm_' ...
,'a_each_prm___' ...
,'A_each_prm___' ...
,'BB_inv_each_prm___' ...
,'CC_inv_each_prm___' ...
,'L_each_prm___' ...
);
%%%%%%%%;
end;%if ( exist(fname_mat,'file'));
save(fname_recalc_mat,'fname_recalc_mat');
close_fname_tmp(fname_recalc_pre);
end;%if ~flag_recalc_skip;
end;%if (~isempty(nvar_recalc_));

disp('returning'); return;

