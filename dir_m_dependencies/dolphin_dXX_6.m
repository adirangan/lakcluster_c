function dolphin_dXX_6(dXX,n_shuffle,age_percentile_range_);
%{
note that age-range is 0.64 --> 54.35. ;
try: ;
%}

if (nargin<1);
dXX = 0; n_shuffle = 256;
dolphin_dXX_6(dXX,n_shuffle,[ 0, 50]);
dolphin_dXX_6(dXX,n_shuffle,[50,100]);
dolphin_dXX_6(dXX,n_shuffle);
%for dXX=0:145;
%n_shuffle = 64;
%try;dolphin_dXX_6(dXX,n_shuffle);catch;disp(sprintf(' %% error processing dXX %.3d',dXX)); end;%try;
%end;%for dXX=0:145;
disp('returning'); return;
end;%if (nargin<1);

if (nargin<2); n_shuffle = []; end;
if (nargin<3); age_percentile_range_ = []; end;
if (isempty(n_shuffle)); n_shuffle = 32; end;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_dolphin',string_root);
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);
%dir_bootstrap_mat = sprintf('%s/dir_bootstrap_mat',dir_trunk);
dir_shuffle_mat = sprintf('%s/dir_shuffle%d_mat',dir_trunk,n_shuffle);
if (~exist(dir_shuffle_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_shuffle_mat)); mkdir(dir_shuffle_mat); end;

flag_replot=0;
flag_recalc=0;

flag_drift0=1;
flag_center=1;
flag_constrain_CC=1;

if ( isempty(age_percentile_range_));
fname_infix = sprintf('aid%.2dd%dc%dC%d',dXX,flag_drift0,flag_center,flag_constrain_CC);
end;%if ( isempty(age_percentile_range_));
if (~isempty(age_percentile_range_));
fname_infix = sprintf('aid%.2d_age%.2d%.2d_d%dc%dC%d',dXX,min(age_percentile_range_),min(99,max(age_percentile_range_)),flag_drift0,flag_center,flag_constrain_CC);
end;%if (~isempty(age_percentile_range_));
disp(sprintf(' %% fname_infix %s',fname_infix));

dolphin_dXX_load_3;

%%%%%%%%;
% select age range. ;
%%%%%%%%;
if (~isempty(age_percentile_range_));
tmp_index_ = efind(age_>=prctile(age_,min(age_percentile_range_)) & age_<=prctile(age_,max(age_percentile_range_)));
disp(sprintf(' %% age_percentile_range_ [%0.2f %0.2f] --> [%0.2f %0.2f] <-- %d found',age_percentile_range_,prctile(age_,age_percentile_range_),numel(tmp_index_)));
aid_ = aid_(1+tmp_index_);
age_ = age_(1+tmp_index_);
dat__ = dat__(1+tmp_index_,:);
end;%if (~isempty(age_percentile_range_));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% processing full data-matrix. ;
%%%%%%%%;
fname_mat = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix);
fname_tmp = sprintf('%s/dolphin_%s_s0000.tmp',dir_mat,fname_infix);
if (flag_recalc | (~exist(fname_mat,'file') & ~exist(fname_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_tmp,'fname_mat');
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.n_step = 1;
parameter.n_iteration = 4;
parameter.flag_constrain_CC = flag_constrain_CC;
parameter.flag_crude = 0;
%%%%%%%%;
[ ...
 a_ ...
,A__ ...
,BB__ ...
,CC__ ...
,L ...
,niteration ...
,a__ ...
,A___ ...
,BB___ ...
,CC___ ...
,L_ ...
,index_min ...
] = ...
dolphin_estimate_aABC_2( ...
 aid_ ...
,age_ ...
,dat__ ...
,parameter ...
);
save(fname_mat ...
,'parameter' ...
,'a_' ...
,'A__' ...
,'BB__' ...
,'CC__' ...
,'L' ...
,'niteration' ...
,'a__' ...
,'A___' ...
,'BB___' ...
,'CC___' ...
,'L_' ...
,'index_min' ...
);
%%%%%%%%;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
end;%if ( exist(fname_mat,'file'));

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
fname_fig = sprintf('%s/dolphin_%s_aABC_2',dir_jpg,fname_infix);
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
imagesc(BB__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('SDE noise BB__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(CC__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('observation noise CC__','Interpreter','none');
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
fname_mat = sprintf('%s/dolphin_%s_sxxxx.mat',dir_shuffle_mat,fname_infix);
fname_tmp = sprintf('%s/dolphin_%s_sxxxx.tmp',dir_shuffle_mat,fname_infix);
if (flag_recalc | (~exist(fname_mat,'file') & ~exist(fname_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_tmp,'fname_mat');
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.n_step = 1;
parameter.n_iteration = 4;
parameter.flag_constrain_CC = flag_constrain_CC;
parameter.flag_crude = 0;
%%%%%%%%;
a_prm__ = zeros(n_var,1+n_shuffle);
A_prm___ = zeros(n_var,n_var,1+n_shuffle);
BB_prm___ = zeros(n_var,n_var,1+n_shuffle);
CC_prm___ = zeros(n_var,n_var,1+n_shuffle);
L_prm_ = zeros(1+n_shuffle,1);
for nshuffle=0:n_shuffle;
if (nshuffle==0); dat_prm__ = dat__; end;
if (nshuffle> 0); dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle); end;
if (mod(nshuffle,32)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
%%%%%%%%;
[ ...
 a_prm_ ...
,A_prm__ ...
,BB_prm__ ...
,CC_prm__ ...
,L_prm ...
] = ...
dolphin_estimate_aABC_2( ...
 aid_ ...
,age_ ...
,dat_prm__ ...
,parameter ...
);
a_prm__(:,1+nshuffle) = a_prm_;
A_prm___(:,:,1+nshuffle) = A_prm__;
BB_prm___(:,:,1+nshuffle) = BB_prm__;
CC_prm___(:,:,1+nshuffle) = CC_prm__;
L_prm_(1+nshuffle) = L_prm;
end;%for nshuffle=0:n_shuffle;
save(fname_mat ...
,'parameter' ...
,'a_prm__' ...
,'A_prm___' ...
,'BB_prm___' ...
,'CC_prm___' ...
,'L_prm_' ...
);
delete(fname_tmp);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_mat));
end;%if ( exist(fname_mat,'file'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% processing variable pairs. ;
%%%%%%%%;
fname_mat = sprintf('%s/dolphin_%s_cmb_s0000.mat',dir_mat,fname_infix);
fname_tmp = sprintf('%s/dolphin_%s_cmb_s0000.tmp',dir_mat,fname_infix);
if (flag_recalc | (~exist(fname_mat,'file') & ~exist(fname_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_tmp,'fname_mat');
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.n_step = 1;
parameter.n_iteration = 4;
parameter.flag_constrain_CC = flag_constrain_CC;
parameter.flag_crude = 0;
%%%%%%%%;
a_cmb_ = zeros(n_var,1);
A_cmb__ = zeros(n_var,n_var);
BB_cmb__ = zeros(n_var,n_var);
CC_cmb__ = zeros(n_var,n_var);
L_cmb = 0;
a_each__ = zeros(n_var,n_var);
A_each__ = zeros(n_var,n_var);
BB_each__ = zeros(n_var,n_var);
CC_each__ = zeros(n_var,n_var);
L_each__ = zeros(n_var,n_var);
for nvar0=0:n_var-1;
disp(sprintf(' %% nvar0 %d/%d',nvar0,n_var));
for nvar1=nvar0+1:n_var-1;
dat_pair__ = dat__(:,1+[nvar0,nvar1]);
[ ...
 a_pair_ ...
,A_pair__ ...
,BB_pair__ ...
,CC_pair__ ...
,L_pair ...
] = ...
dolphin_estimate_aABC_2( ...
 aid_ ...
,age_ ...
,dat_pair__ ...
,parameter ...
);
a_cmb_(1+[nvar0;nvar1]) = a_cmb_(1+[nvar0;nvar1]) + a_pair_(:);
A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + A_pair__;
BB_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = BB_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + BB_pair__;
CC_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = CC_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + CC_pair__;
L_cmb = L_cmb + L_pair;
a_each__(1+nvar0,1+nvar1) = a_pair_(1+0); a_each__(1+nvar1,1+nvar0) = a_pair_(1+1);
A_each__(1+nvar0,1+nvar1) = A_pair__(1+0,1+0); A_each__(1+nvar1,1+nvar0) = A_pair__(1+1,1+1);
BB_each__(1+nvar0,1+nvar1) = BB_pair__(1+0,1+0); BB_each__(1+nvar1,1+nvar0) = BB_pair__(1+1,1+1);
CC_each__(1+nvar0,1+nvar1) = CC_pair__(1+0,1+0); CC_each__(1+nvar1,1+nvar0) = CC_pair__(1+1,1+1);
L_each__(1+nvar0,1+nvar1) = L_pair; L_each__(1+nvar1,1+nvar0) = L_pair;
end;%for nvar1=nvar0:n_var-1;
end;%for nvar0=0:n_var-1;
save(fname_mat ...
,'parameter' ...
,'a_cmb_' ...
,'A_cmb__' ...
,'BB_cmb__' ...
,'CC_cmb__' ...
,'L_cmb' ...
,'a_each__' ...
,'A_each__' ...
,'BB_each__' ...
,'CC_each__' ...
,'L_each__' ...
);
%%%%%%%%;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
end;%if ( exist(fname_mat,'file'));

fname_mat = sprintf('%s/dolphin_%s_cmb_s0000.mat',dir_mat,fname_infix);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, loading',fname_mat));
load(fname_mat);
%%%%%%%%;
% plot SDE-regression results. ;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_cmb_aABC_2',dir_jpg,fname_infix);
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
imagesc(BB_cmb__ - diag(diag(BB_cmb__))); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('SDE noise BB_cmb__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(CC_cmb__ - diag(diag(CC_cmb__))); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('observation noise CC_cmb__','Interpreter','none');
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
fname_mat = sprintf('%s/dolphin_%s_cmb_sxxxx.mat',dir_shuffle_mat,fname_infix);
fname_tmp = sprintf('%s/dolphin_%s_cmb_sxxxx.tmp',dir_shuffle_mat,fname_infix);
if (flag_recalc | (~exist(fname_mat,'file') & ~exist(fname_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_tmp,'fname_mat');
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.n_step = 1;
parameter.n_iteration = 4;
parameter.flag_constrain_CC = flag_constrain_CC;
parameter.flag_crude = 0;
%%%%%%%%;
a_cmb_prm__ = zeros(n_var,1+n_shuffle);
A_cmb_prm___ = zeros(n_var,n_var,1+n_shuffle);
BB_cmb_prm___ = zeros(n_var,n_var,1+n_shuffle);
CC_cmb_prm___ = zeros(n_var,n_var,1+n_shuffle);
L_cmb_prm_ = zeros(1+n_shuffle,1);
a_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
A_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
BB_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
CC_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
L_each_prm___ = zeros(n_var,n_var,1+n_shuffle);
for nshuffle=0:n_shuffle;
if (nshuffle==0); dat_prm__ = dat__; end;
if (nshuffle> 0); dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle); end;
if (mod(nshuffle,32)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
%%%%%%%%;
a_cmb_ = zeros(n_var,1);
A_cmb__ = zeros(n_var,n_var);
BB_cmb__ = zeros(n_var,n_var);
CC_cmb__ = zeros(n_var,n_var);
L_cmb = 0;
a_each__ = zeros(n_var,n_var);
A_each__ = zeros(n_var,n_var);
BB_each__ = zeros(n_var,n_var);
CC_each__ = zeros(n_var,n_var);
L_each__ = zeros(n_var,n_var);
for nvar0=0:n_var-1;
disp(sprintf(' %% nvar0 %d/%d',nvar0,n_var));
for nvar1=nvar0+1:n_var-1;
dat_pair__ = dat_prm__(:,1+[nvar0,nvar1]);
[ ...
 a_pair_ ...
,A_pair__ ...
,BB_pair__ ...
,CC_pair__ ...
,L_pair ...
] = ...
dolphin_estimate_aABC_2( ...
 aid_ ...
,age_ ...
,dat_pair__ ...
,parameter ...
);
a_cmb_(1+[nvar0;nvar1]) = a_cmb_(1+[nvar0;nvar1]) + a_pair_(:);
A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = A_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + A_pair__;
BB_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = BB_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + BB_pair__;
CC_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) = CC_cmb__(1+[nvar0;nvar1],1+[nvar0;nvar1]) + CC_pair__;
L_cmb = L_cmb + L_pair;
a_each__(1+nvar0,1+nvar1) = a_pair_(1+0); a_each__(1+nvar1,1+nvar0) = a_pair_(1+1);
A_each__(1+nvar0,1+nvar1) = A_pair__(1+0,1+0); A_each__(1+nvar1,1+nvar0) = A_pair__(1+1,1+1);
BB_each__(1+nvar0,1+nvar1) = BB_pair__(1+0,1+0); BB_each__(1+nvar1,1+nvar0) = BB_pair__(1+1,1+1);
CC_each__(1+nvar0,1+nvar1) = CC_pair__(1+0,1+0); CC_each__(1+nvar1,1+nvar0) = CC_pair__(1+1,1+1);
L_each__(1+nvar0,1+nvar1) = L_pair; L_each__(1+nvar1,1+nvar0) = L_pair;
end;%for nvar1=nvar0:n_var-1;
end;%for nvar0=0:n_var-1;
a_cmb_prm__(:,1+nshuffle) = a_cmb_;
A_cmb_prm___(:,:,1+nshuffle) = A_cmb__;
BB_cmb_prm___(:,:,1+nshuffle) = BB_cmb__;
CC_cmb_prm___(:,:,1+nshuffle) = CC_cmb__;
L_cmb_prm_(1+nshuffle) = L_cmb;
a_each_prm___(:,:,1+nshuffle) = a_each__;
A_each_prm___(:,:,1+nshuffle) = A_each__;
BB_each_prm___(:,:,1+nshuffle) = BB_each__;
CC_each_prm___(:,:,1+nshuffle) = CC_each__;
L_each_prm___(:,:,1+nshuffle) = L_each__;
end;%for nshuffle=0:n_shuffle;
save(fname_mat ...
,'parameter' ...
,'a_cmb_prm__' ...
,'A_cmb_prm___' ...
,'BB_cmb_prm___' ...
,'CC_cmb_prm___' ...
,'L_cmb_prm_' ...
,'a_each_prm___' ...
,'A_each_prm___' ...
,'BB_each_prm___' ...
,'CC_each_prm___' ...
,'L_each_prm___' ...
);
delete(fname_tmp);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_mat));
end;%if ( exist(fname_mat,'file'));

disp('returning'); return;

