function dolphin_dXX_bootstrap_3(dXX,n_shuffle,age_percentile_range_);
%{
note that age-range is 0.64 --> 54.35. ;
try: ;

dXX = 0; n_shuffle = 64;
dolphin_dXX_bootstrap_3(dXX,n_shuffle,[ 0, 50]);
dolphin_dXX_bootstrap_3(dXX,n_shuffle,[25, 75]);
dolphin_dXX_bootstrap_3(dXX,n_shuffle,[50,100]);
dolphin_dXX_bootstrap_3(dXX,n_shuffle);

%}

if (nargin<2); n_shuffle = []; end;
if (nargin<3); age_percentile_range_ = []; end;
if (isempty(n_shuffle)); n_shuffle = 32; end;

setup_OptiPlex; dir_trunk = '/home/rangan/dir_bcc/dir_dolphin';
%setup_access1; dir_trunk = '/data/rangan/dir_bcc/dir_dolphin';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);
dir_bootstrap_mat = sprintf('%s/dir_bootstrap_mat',dir_trunk);
dir_shuffle_mat = sprintf('%s/dir_shuffle_mat',dir_trunk);

flag_replot=1;
flag_recalc=0;

flag_drift0=1;
flag_a_is_0=0;
flag_center=1;
%dXX = 0;

if ( isempty(age_percentile_range_));
fname_infix = sprintf('aid%.2dd%dc%d',dXX,flag_drift0,flag_center);
end;%if ( isempty(age_percentile_range_));
if (~isempty(age_percentile_range_));
fname_infix = sprintf('aid%.2d_age%.2d%.2d_d%dc%d',dXX,min(age_percentile_range_),min(99,max(age_percentile_range_)),flag_drift0,flag_center);
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

index_mss_ = efind(~isfinite(dat__));

%%%%%%%%;
% processing full data-matrix. ;
%%%%%%%%;
fname_mat = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix);
fname_tmp = sprintf('%s/dolphin_%s_s0000.tmp',dir_mat,fname_infix);
if (flag_recalc | (~exist(fname_mat,'file') & ~exist(fname_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_tmp,'fname_mat');
%%%%%%%%;
BB__ = []; CC__ = [];
%%%%%%%%;
n_step = 1; n_iteration = 4;
%%%%%%%%;
if flag_a_is_0==0;
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_0(aid_,age_,dat__,n_step,n_iteration);
save(fname_mat ...
,'n_step' ...
,'n_iteration' ...
,'a_' ...
,'A__' ...
,'BB__' ...
,'CC__' ...
,'L_' ...
,'niteration' ...
,'a__' ...
,'A___' ...
,'BB___' ...
,'CC___' ...
);
end;%if flag_a_is_0==0;
if flag_a_is_0==1;
[A__,BB__,CC__,L_,niteration,A___,BB___,CC___] = dolphin_estimate_ABC_0(aid_,age_,dat__,n_step,n_iteration);
save(fname_mat ...
,'n_step' ...
,'n_iteration' ...
,'A__' ...
,'BB__' ...
,'CC__' ...
,'L_' ...
,'niteration' ...
,'A___' ...
,'BB___' ...
,'CC___' ...
);
end;%if flag_a_is_0==1;
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
fname_fig = sprintf('%s/dolphin_%s_aABC_0',dir_jpg,fname_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
if flag_a_is_0==0;
subplot(2,2,1);
bar(a_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
title('drift a_','Interpreter','none');
end;%if flag_a_is_0==0;
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
% Now repeat for several bootstrapped trials. ;
% default bootstrap_fraction = 0.5;
%%%%%%%%;
n_step = 1; n_iteration = 4;
for nshuffle=1:n_shuffle;
[aid_0_,age_0_,dat_0__,aid_1_,age_1_,dat_1__,n_smp_0,n_smp_1,prm_aid__,index_aid_0__,index_aid_1__] = dolphin_bootstrap_0(aid_,age_,dat__,nshuffle);
%%%%%%%%;
fname_0_mat = sprintf('%s/dolphin_%s_b0_%.4d.mat',dir_bootstrap_mat,fname_infix,nshuffle);
fname_0_tmp = sprintf('%s/dolphin_%s_b0_%.4d.tmp',dir_bootstrap_mat,fname_infix,nshuffle);
if (flag_recalc | (~exist(fname_0_mat,'file') & ~exist(fname_0_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_0_mat));
save(fname_0_tmp,'fname_0_mat');
if flag_a_is_0==0;
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_0(aid_0_,age_0_,dat_0__,n_step,n_iteration);
save(fname_0_mat ...
,'nshuffle' ... 
,'n_step' ... 
,'n_iteration' ... 
,'a_' ... 
,'A__' ... 
,'BB__' ... 
,'CC__' ... 
,'L_' ...
,'niteration' ...
,'a__' ... 
,'A___' ... 
,'BB___' ... 
,'CC___' ... 
     );
end;%if flag_a_is_0==0;
if flag_a_is_0==1;
[A__,BB__,CC__,L_,niteration,A___,BB___,CC___] = dolphin_estimate_ABC_0(aid_0_,age_0_,dat_0__,n_step,n_iteration);
save(fname_0_mat ...
,'nshuffle' ... 
,'n_step' ... 
,'n_iteration' ... 
,'A__' ... 
,'BB__' ... 
,'CC__' ... 
,'L_' ...
,'niteration' ...
,'A___' ... 
,'BB___' ... 
,'CC___' ... 
     );
end;%if flag_a_is_0==1;
delete(fname_0_tmp);
end;%if (~exist(fname_0_mat,'file'));
if ( exist(fname_0_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_0_mat));
end;%if ( exist(fname_0_mat,'file'));
%%%%%%%%;
fname_1_mat = sprintf('%s/dolphin_%s_b1_%.4d.mat',dir_bootstrap_mat,fname_infix,nshuffle);
fname_1_tmp = sprintf('%s/dolphin_%s_b1_%.4d.tmp',dir_bootstrap_mat,fname_infix,nshuffle);
if (flag_recalc | (~exist(fname_1_mat,'file') & ~exist(fname_1_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_1_mat));
save(fname_1_tmp,'fname_1_mat');
if flag_a_is_0==0;
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_0(aid_1_,age_1_,dat_1__,n_step,n_iteration);
save(fname_1_mat ...
,'nshuffle' ... 
,'n_step' ... 
,'n_iteration' ... 
,'a_' ... 
,'A__' ... 
,'BB__' ... 
,'CC__' ... 
,'L_' ...
,'niteration' ...
,'a__' ... 
,'A___' ... 
,'BB___' ... 
,'CC___' ... 
     );
end;%if flag_a_is_0==0;
if flag_a_is_0==1;
[A__,BB__,CC__,L_,niteration,A___,BB___,CC___] = dolphin_estimate_ABC_0(aid_1_,age_1_,dat_1__,n_step,n_iteration);
save(fname_1_mat ...
,'nshuffle' ... 
,'n_step' ... 
,'n_iteration' ... 
,'A__' ... 
,'BB__' ... 
,'CC__' ... 
,'L_' ...
,'niteration' ...
,'A___' ... 
,'BB___' ... 
,'CC___' ... 
     );
end;%if flag_a_is_0==1;
delete(fname_1_tmp);
end;%if (~exist(fname_1_mat,'file'));
if ( exist(fname_1_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_1_mat));
end;%if ( exist(fname_1_mat,'file'));
%%%%%%%%;
end;%for nshuffle=1:n_shuffle;

%%%%%%%%;
% Now look at bootstrapped iterates. ;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_bootstrap_A_FIGA',dir_jpg,fname_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;figbeach();
A_lim_ = [-1,+1];
np=0;
%subplot(4,8,1+np);np=np+1;
%fname_mat = sprintf('%s/dolphin_%s_s0000.mat',dir_bootstrap_mat,fname_infix);
%tmp_ = load(fname_mat); A_ori__ = tmp_.A__; imagesc(A_ori__,A_lim_); axis image; axisnotick; title('ori');
for nshuffle=1:min(16,n_shuffle);
subplot(4,8,1+np);np=np+1;
fname_mat = sprintf('%s/dolphin_%s_b0_%.4d.mat',dir_bootstrap_mat,fname_infix,nshuffle);
if ( exist(fname_mat,'file'));
tmp_ = load(fname_mat); A_b0__ = tmp_.A__; imagesc(A_b0__,A_lim_); axis image; axisnotick; title(sprintf('b%d A',nshuffle));
end;%if ( exist(fname_mat,'file'));
subplot(4,8,1+np);np=np+1;
fname_mat = sprintf('%s/dolphin_%s_b1_%.4d.mat',dir_bootstrap_mat,fname_infix,nshuffle);
if ( exist(fname_mat,'file'));
tmp_ = load(fname_mat); A_b1__ = tmp_.A__; imagesc(A_b1__,A_lim_); axis image; axisnotick; title(sprintf('b%d B',nshuffle));
end;%if ( exist(fname_mat,'file'));
end;%for nshuffle=1:14;
figbig;
sgtitle(sprintf('bootstrap iterates (%s)',fname_infix),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now estimate nlp for each interaction-term using bootstrapped iterates. ;
%%%%%%%%;
fname_mat = sprintf('%s/dolphin_%s_s0000.mat',dir_bootstrap_mat,fname_infix);
if ( exist(fname_mat,'file'));
%%%%%%%%;
tmp_ = load(fname_mat); A_ori__ = tmp_.A__; 
A_avg__ = zeros(size(A_ori__));
A_std__ = zeros(size(A_ori__));
n_found = 0;
for nshuffle=1:n_shuffle;
fname_mat_0 = sprintf('%s/dolphin_%s_b0_%.4d.mat',dir_bootstrap_mat,fname_infix,nshuffle);
if ( exist(fname_mat_0,'file'));
tmp_ = load(fname_mat_0); A_b0__ = tmp_.A__; A_avg__ = A_avg__ + A_b0__.^1; A_std__ = A_std__ + A_b0__.^2; n_found = n_found+1;
end;%if ( exist(fname_mat_0,'file'));
fname_mat_1 = sprintf('%s/dolphin_%s_b1_%.4d.mat',dir_bootstrap_mat,fname_infix,nshuffle);
if ( exist(fname_mat_1,'file'));
tmp_ = load(fname_mat_1); A_b1__ = tmp_.A__; A_avg__ = A_avg__ + A_b1__.^1; A_std__ = A_std__ + A_b1__.^2; n_found = n_found+1;
end;%if ( exist(fname_mat_1,'file'));
end;%for nshuffle=1:n_shuffle;
A_avg__ = A_avg__/n_found; A_std__ = sqrt(A_std__/n_found - A_avg__.^2);
A_nlp__ = -z_to_lp(A_avg__./A_std__);
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_bootstrap_A_FIGB',dir_jpg,fname_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); figure(1);
figure(1);clf;figbig;figbeach();fontsize_use = 8;
subplot_(1) = subplot(2,2,1); imagesc(A_ori__,A_lim_); axis image; axisnotick; title('ori'); colorbar;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',fontsize_use);
subplot_(2) = subplot(2,2,2); imagesc(A_avg__,A_lim_); axis image; axisnotick; title('avg'); colorbar;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',fontsize_use);
subplot_(3) = subplot(2,2,3); imagesc(A_std__,A_lim_); axis image; axisnotick; title('std'); colorbar;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',fontsize_use);
subplot_(4) = subplot(2,2,4); imagesc(A_nlp__,[0,27]); axis image; axisnotick; title('nlp (local)'); colorbar;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',fontsize_use);
sgtitle('original vs bootstrap');
colormap(subplot_(1),colormap_beach);
colormap(subplot_(2),colormap_beach);
colormap(subplot_(3),colormap_beach);
colormap(subplot_(4),colormap_nlpvt(64));
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

disp('returning'); return;

%%%%%%%%;
% Now repeat for several permuted trials. ;
%%%%%%%%;
n_step = 1; n_iteration = 4;
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_%s_s%.4d.mat',dir_shuffle_mat,fname_infix,nshuffle);
fname_tmp = sprintf('%s/dolphin_%s_s%.4d.tmp',dir_shuffle_mat,fname_infix,nshuffle);
if (flag_recalc | (~exist(fname_mat,'file') & ~exist(fname_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_tmp,'fname_mat');
dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle);
if flag_a_is_0==0;
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_0(aid_,age_,dat_prm__,n_step,n_iteration);
save(fname_mat ...
,'nshuffle' ... 
,'n_step' ... 
,'n_iteration' ... 
,'a_' ... 
,'A__' ... 
,'BB__' ... 
,'CC__' ... 
,'L_' ...
,'niteration' ...
,'a__' ... 
,'A___' ... 
,'BB___' ... 
,'CC___' ... 
     );
end;%if flag_a_is_0==0;
if flag_a_is_0==1;
[A__,BB__,CC__,L_,niteration,A___,BB___,CC___] = dolphin_estimate_ABC_0(aid_,age_,dat_prm__,n_step,n_iteration);
save(fname_mat ...
,'nshuffle' ... 
,'n_step' ... 
,'n_iteration' ... 
,'A__' ... 
,'BB__' ... 
,'CC__' ... 
,'L_' ...
,'niteration' ...
,'A___' ... 
,'BB___' ... 
,'CC___' ... 
     );
end;%if flag_a_is_0==1;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_mat));
end;%if ( exist(fname_mat,'file'));
end;%for nshuffle=1:n_shuffle;

flag_all = 1;
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_%s_s%.4d.mat',dir_shuffle_mat,fname_infix,nshuffle);
flag_all = flag_all * exist(fname_mat,'file');
end;%for nshuffle=1:n_shuffle;
if flag_all;
%%%%%%%%;
% collect null distribution. ;
%%%%%%%%;
if flag_a_is_0==0; a_null_avg_ = zeros(n_var,1); end;
A_null_avg__ = zeros(n_var,n_var);
BB_null_avg__ = zeros(n_var,n_var);
CC_null_avg__ = zeros(n_var,n_var);
if flag_a_is_0==0; a_null_std_ = zeros(n_var,1); end;
A_null_std__ = zeros(n_var,n_var);
BB_null_std__ = zeros(n_var,n_var);
CC_null_std__ = zeros(n_var,n_var);
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_%s_s%.4d.mat',dir_shuffle_mat,fname_infix,nshuffle);
tmp_ = load(fname_mat);
[~,tmp_index] = min(tmp_.L_(1:1+tmp_.niteration)); tmp_index = tmp_index-1;
if flag_a_is_0==0; a_null_avg_ = a_null_avg_ + tmp_.a__{1+tmp_index}; end;
A_null_avg__ = A_null_avg__ + tmp_.A___{1+tmp_index};
BB_null_avg__ = BB_null_avg__ + tmp_.BB___{1+tmp_index};
CC_null_avg__ = CC_null_avg__ + tmp_.CC___{1+tmp_index};
if flag_a_is_0==0; a_null_std_ = a_null_std_ + tmp_.a__{1+tmp_index}.^2; end;
A_null_std__ = A_null_std__ + tmp_.A___{1+tmp_index}.^2;
BB_null_std__ = BB_null_std__ + tmp_.BB___{1+tmp_index}.^2;
CC_null_std__ = CC_null_std__ + tmp_.CC___{1+tmp_index}.^2;
clear tmp_;
end;%for nshuffle=1:n_shuffle;
if flag_a_is_0==0; a_null_avg_ = a_null_avg_ / n_shuffle; end;
A_null_avg__ = A_null_avg__ / n_shuffle;
BB_null_avg__ = BB_null_avg__ / n_shuffle;
CC_null_avg__ = CC_null_avg__ / n_shuffle;
if flag_a_is_0==0; a_null_std_ = sqrt(a_null_std_/n_shuffle - a_null_avg_.^2); end;
A_null_std__ = sqrt(A_null_std__/n_shuffle - A_null_avg__.^2);
BB_null_std__ = sqrt(BB_null_std__/n_shuffle - BB_null_avg__.^2);
CC_null_std__ = sqrt(CC_null_std__/n_shuffle - CC_null_avg__.^2);
fname_mat = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix);
tmp_ = load(fname_mat);
[~,tmp_index] = min(tmp_.L_(1:1+2)); tmp_index = tmp_index-1;
if flag_a_is_0==0; a_ = tmp_.a__{1+tmp_index}; end;
A__ = tmp_.A___{1+tmp_index};
BB__ = tmp_.BB___{1+tmp_index};
CC__ = tmp_.CC___{1+tmp_index};
clear tmp_;
if flag_a_is_0==0; a_Z_ = real((a_ - a_null_avg_)./a_null_std_); end;
A_Z__ = real((A__ - A_null_avg__)./A_null_std__);
BB_Z__ = real((BB__ - BB_null_avg__)./BB_null_std__);
CC_Z__ = real((CC__ - CC_null_avg__)./CC_null_std__);
if flag_a_is_0==0; a_nlp_ = -z_to_lp(a_Z_); end;
A_nlp__ = -z_to_lp(A_Z__);
BB_nlp__ = -z_to_lp(BB_Z__);
CC_nlp__ = -z_to_lp(CC_Z__);
%%%%%%%%;   
fname_fig = sprintf('%s/dolphin_%s_aABC_nlp_0',dir_jpg,fname_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
if flag_a_is_0==0;
subplot(2,2,1);
bar(a_nlp_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
title('-log p for drift a_','Interpreter','none');
end;%if flag_a_is_0==0;
subplot(2,2,2);
colormap(colormap_beach());
imagesc(A_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('-log p for interaction A__','Interpreter','none');
subplot(2,2,3);
colormap(colormap_beach());
imagesc(BB_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('-log p for SDE noise BB__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(CC_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('-log p for observation noise CC__','Interpreter','none');
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
end;%if flag_all;

%%%%%%%%;
% Now step through every pair of variables and check for pairwise causality. ;
%%%%%%%%;
n_step = 1; n_iteration = 4;
fname_mat = sprintf('%s/dolphin_%s_pair_cmb_s0000.mat',dir_mat,fname_infix);
fname_tmp = sprintf('%s/dolphin_%s_pair_cmb_s0000.tmp',dir_mat,fname_infix);
if (flag_recalc | (~exist(fname_mat,'file') & ~exist(fname_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_tmp,'fname_mat');
%%%%%%%%;
if flag_a_is_0==0; pair_a_cmb_ = zeros(n_var,1); end;
pair_A_cmb__ = zeros(n_var,n_var);
pair_BB_cmb__ = zeros(n_var,n_var);
pair_CC_cmb__ = zeros(n_var,n_var);
n_pair = n_var*(n_var-1)/2;
for nvar_0=0:n_var-1;
for nvar_1=nvar_0+1:n_var-1;
disp(sprintf(' %% nvar_0 %d nvar_1 %d',nvar_0,nvar_1));
pair_dat__ = dat__(:,1+[nvar_0,nvar_1]);
if flag_a_is_0==0;
[pair_a_,pair_A__,pair_BB__,pair_CC__,pair_L_,pair_niteration,pair_a__,pair_A___,pair_BB___,pair_CC___] = dolphin_estimate_aABC_0(aid_,age_,pair_dat__,n_step,n_iteration);
end;%if flag_a_is_0==0;
if flag_a_is_0==1;
[pair_A__,pair_BB__,pair_CC__,pair_L_,pair_niteration,pair_A___,pair_BB___,pair_CC___] = dolphin_estimate_ABC_0(aid_,age_,pair_dat__,n_step,n_iteration);end;%if flag_a_is_0==1;
[~,tmp_index] = min(pair_L_(1:1+pair_niteration)); tmp_index = tmp_index-1;
if flag_a_is_0==0; pair_a_cmb_(1+nvar_0) = pair_a_cmb_(1+nvar_0) + pair_a__{1+tmp_index}(1+0)/n_pair; end;
if flag_a_is_0==0; pair_a_cmb_(1+nvar_1) = pair_a_cmb_(1+nvar_1) + pair_a__{1+tmp_index}(1+1)/n_pair; end;
pair_A_cmb__(1+nvar_0,1+nvar_1) = pair_A___{1+tmp_index}(1+0,1+1);
pair_A_cmb__(1+nvar_1,1+nvar_0) = pair_A___{1+tmp_index}(1+1,1+0);
pair_A_cmb__(1+nvar_0,1+nvar_0) = pair_A_cmb__(1+nvar_0,1+nvar_0) + pair_A___{1+tmp_index}(1+0,1+0)/n_pair;
pair_A_cmb__(1+nvar_1,1+nvar_1) = pair_A_cmb__(1+nvar_1,1+nvar_1) + pair_A___{1+tmp_index}(1+1,1+1)/n_pair;
pair_BB_cmb__(1+nvar_0,1+nvar_1) = pair_BB___{1+tmp_index}(1+0,1+1);
pair_BB_cmb__(1+nvar_1,1+nvar_0) = pair_BB___{1+tmp_index}(1+1,1+0);
pair_BB_cmb__(1+nvar_0,1+nvar_0) = pair_BB_cmb__(1+nvar_0,1+nvar_0) + pair_BB___{1+tmp_index}(1+0,1+0)/n_pair;
pair_BB_cmb__(1+nvar_1,1+nvar_1) = pair_BB_cmb__(1+nvar_1,1+nvar_1) + pair_BB___{1+tmp_index}(1+1,1+1)/n_pair;
pair_CC_cmb__(1+nvar_0,1+nvar_1) = pair_CC___{1+tmp_index}(1+0,1+1);
pair_CC_cmb__(1+nvar_1,1+nvar_0) = pair_CC___{1+tmp_index}(1+1,1+0);
pair_CC_cmb__(1+nvar_0,1+nvar_0) = pair_CC_cmb__(1+nvar_0,1+nvar_0) + pair_CC___{1+tmp_index}(1+0,1+0)/n_pair;
pair_CC_cmb__(1+nvar_1,1+nvar_1) = pair_CC_cmb__(1+nvar_1,1+nvar_1) + pair_CC___{1+tmp_index}(1+1,1+1)/n_pair;
end;%for nvar_1=nvar_0+1:n_var-1;
end;%for nvar_0=0:n_var-1;
%%%%%%%%;
if flag_a_is_0==0;
save(fname_mat ...
,'n_step' ... 
,'n_iteration' ... 
,'pair_a_cmb_' ... 
,'pair_A_cmb__' ... 
,'pair_BB_cmb__' ... 
,'pair_CC_cmb__' ... 
     );
end;%if flag_a_is_0==0;
if flag_a_is_0==1;
save(fname_mat ...
,'n_step' ... 
,'n_iteration' ... 
,'pair_A_cmb__' ... 
,'pair_BB_cmb__' ... 
,'pair_CC_cmb__' ... 
     );
end;%if flag_a_is_0==1;
%%%%%%%%;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_mat));
end;%if ( exist(fname_mat,'file'));

fname_mat = sprintf('%s/dolphin_%s_pair_cmb_s0000.mat',dir_mat,fname_infix);
if ( exist(fname_mat,'file'));
load(fname_mat);
fname_fig = sprintf('%s/dolphin_%s_pair_aABC_cmb_0',dir_jpg,fname_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
if flag_a_is_0==0;
subplot(2,2,1);
bar(pair_a_cmb_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
title('drift pair_a_cmb_','Interpreter','none');
end;%if flag_a_is_0==0;
subplot(2,2,2);
colormap(colormap_beach());
imagesc(pair_A_cmb__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('interaction pair_A_cmb__','Interpreter','none');
subplot(2,2,3);
colormap(colormap_beach());
imagesc(pair_BB_cmb__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('SDE noise pair_BB_cmb__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(pair_CC_cmb__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('observation noise pair_CC_cmb__','Interpreter','none');
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
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% Now repeat for several permuted trials. ;
%%%%%%%%;
n_step = 1; n_iteration = 4;
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_%s_pair_cmb_s%.4d.mat',dir_shuffle_mat,fname_infix,nshuffle);
fname_tmp = sprintf('%s/dolphin_%s_pair_cmb_s%.4d.tmp',dir_shuffle_mat,fname_infix,nshuffle);
if (flag_recalc | (~exist(fname_mat,'file') & ~exist(fname_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_tmp,'fname_mat');
dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle);
%%%%%%%%;
if flag_a_is_0==0; pair_a_cmb_ = zeros(n_var,1); end;
pair_A_cmb__ = zeros(n_var,n_var);
pair_BB_cmb__ = zeros(n_var,n_var);
pair_CC_cmb__ = zeros(n_var,n_var);
n_pair = n_var*(n_var-1)/2;
for nvar_0=0:n_var-1;
for nvar_1=nvar_0+1:n_var-1;
disp(sprintf(' %% nvar_0 %d nvar_1 %d',nvar_0,nvar_1));
pair_dat_prm__ = dat_prm__(:,1+[nvar_0,nvar_1]);
if flag_a_is_0==0;
[pair_a_,pair_A__,pair_BB__,pair_CC__,pair_L_,pair_niteration,pair_a__,pair_A___,pair_BB___,pair_CC___] = dolphin_estimate_aABC_0(aid_,age_,pair_dat_prm__,n_step,n_iteration);
end;%if flag_a_is_0==0;
if flag_a_is_0==1;
[pair_A__,pair_BB__,pair_CC__,pair_L_,pair_niteration,pair_A___,pair_BB___,pair_CC___] = dolphin_estimate_ABC_0(aid_,age_,pair_dat_prm__,n_step,n_iteration);
end;%if flag_a_is_0==1;
[~,tmp_index] = min(pair_L_(1:1+pair_niteration)); tmp_index = tmp_index-1;
if flag_a_is_0==0; pair_a_cmb_(1+nvar_0) = pair_a_cmb_(1+nvar_0) + pair_a__{1+tmp_index}(1+0)/n_pair; end;
if flag_a_is_0==0; pair_a_cmb_(1+nvar_1) = pair_a_cmb_(1+nvar_1) + pair_a__{1+tmp_index}(1+1)/n_pair; end;
pair_A_cmb__(1+nvar_0,1+nvar_1) = pair_A___{1+tmp_index}(1+0,1+1);
pair_A_cmb__(1+nvar_1,1+nvar_0) = pair_A___{1+tmp_index}(1+1,1+0);
pair_A_cmb__(1+nvar_0,1+nvar_0) = pair_A_cmb__(1+nvar_0,1+nvar_0) + pair_A___{1+tmp_index}(1+0,1+0)/n_pair;
pair_A_cmb__(1+nvar_1,1+nvar_1) = pair_A_cmb__(1+nvar_1,1+nvar_1) + pair_A___{1+tmp_index}(1+1,1+1)/n_pair;
pair_BB_cmb__(1+nvar_0,1+nvar_1) = pair_BB___{1+tmp_index}(1+0,1+1);
pair_BB_cmb__(1+nvar_1,1+nvar_0) = pair_BB___{1+tmp_index}(1+1,1+0);
pair_BB_cmb__(1+nvar_0,1+nvar_0) = pair_BB_cmb__(1+nvar_0,1+nvar_0) + pair_BB___{1+tmp_index}(1+0,1+0)/n_pair;
pair_BB_cmb__(1+nvar_1,1+nvar_1) = pair_BB_cmb__(1+nvar_1,1+nvar_1) + pair_BB___{1+tmp_index}(1+1,1+1)/n_pair;
pair_CC_cmb__(1+nvar_0,1+nvar_1) = pair_CC___{1+tmp_index}(1+0,1+1);
pair_CC_cmb__(1+nvar_1,1+nvar_0) = pair_CC___{1+tmp_index}(1+1,1+0);
pair_CC_cmb__(1+nvar_0,1+nvar_0) = pair_CC_cmb__(1+nvar_0,1+nvar_0) + pair_CC___{1+tmp_index}(1+0,1+0)/n_pair;
pair_CC_cmb__(1+nvar_1,1+nvar_1) = pair_CC_cmb__(1+nvar_1,1+nvar_1) + pair_CC___{1+tmp_index}(1+1,1+1)/n_pair;
end;%for nvar_1=nvar_0+1:n_var-1;
end;%for nvar_0=0:n_var-1;
%%%%%%%%;
if flag_a_is_0==0;
save(fname_mat ...
,'n_step' ... 
,'n_iteration' ...
,'nshuffle' ...
,'pair_a_cmb_' ... 
,'pair_A_cmb__' ... 
,'pair_BB_cmb__' ... 
,'pair_CC_cmb__' ... 
     );
end;%if flag_a_is_0==0;
if flag_a_is_0==1;
save(fname_mat ...
,'n_step' ... 
,'n_iteration' ...
,'nshuffle' ...
,'pair_A_cmb__' ... 
,'pair_BB_cmb__' ... 
,'pair_CC_cmb__' ... 
     );
end;%if flag_a_is_0==1;
%%%%%%%%;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_mat));
end;%if ( exist(fname_mat,'file'));
end;%for nshuffle=1:n_shuffle;

flag_all = 1;
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_%s_pair_cmb_s%.4d.mat',dir_shuffle_mat,fname_infix,nshuffle);
flag_all = flag_all * exist(fname_mat,'file');
end;%for nshuffle=1:n_shuffle;
if flag_all;
%%%%%%%%;
% collect null distribution. ;
%%%%%%%%;
if flag_a_is_0==0; pair_a_cmb_null_avg_ = zeros(n_var,1); end;
pair_A_cmb_null_avg__ = zeros(n_var,n_var);
pair_BB_cmb_null_avg__ = zeros(n_var,n_var);
pair_CC_cmb_null_avg__ = zeros(n_var,n_var);
if flag_a_is_0==0; pair_a_cmb_null_std_ = zeros(n_var,1); end;
pair_A_cmb_null_std__ = zeros(n_var,n_var);
pair_BB_cmb_null_std__ = zeros(n_var,n_var);
pair_CC_cmb_null_std__ = zeros(n_var,n_var);
n_shuffle_use = 32;
for nshuffle=1:n_shuffle_use;
fname_mat = sprintf('%s/dolphin_%s_pair_cmb_s%.4d.mat',dir_shuffle_mat,fname_infix,nshuffle);
tmp_ = load(fname_mat);
disp(sprintf(' %% norms: pair_A__ %0.2f pair_BB__ %0.2f pair_CC__ %0.2f',fnorm(tmp_.pair_A_cmb__),fnorm(tmp_.pair_BB_cmb__),fnorm(tmp_.pair_CC_cmb__)));
if flag_a_is_0==0; pair_a_cmb_null_avg_ = pair_a_cmb_null_avg_ + tmp_.pair_a_cmb_; end;
pair_A_cmb_null_avg__ = pair_A_cmb_null_avg__ + tmp_.pair_A_cmb__;
pair_BB_cmb_null_avg__ = pair_BB_cmb_null_avg__ + tmp_.pair_BB_cmb__;
pair_CC_cmb_null_avg__ = pair_CC_cmb_null_avg__ + tmp_.pair_CC_cmb__;
if flag_a_is_0==0; pair_a_cmb_null_std_ = pair_a_cmb_null_std_ + tmp_.pair_a_cmb_.^2; end;
pair_A_cmb_null_std__ = pair_A_cmb_null_std__ + tmp_.pair_A_cmb__.^2;
pair_BB_cmb_null_std__ = pair_BB_cmb_null_std__ + tmp_.pair_BB_cmb__.^2;
pair_CC_cmb_null_std__ = pair_CC_cmb_null_std__ + tmp_.pair_CC_cmb__.^2;
clear tmp_;
end;%for nshuffle=1:n_shuffle_use;
if flag_a_is_0==0; pair_a_cmb_null_avg_ = pair_a_cmb_null_avg_ / n_shuffle_use; end;
pair_A_cmb_null_avg__ = pair_A_cmb_null_avg__ / n_shuffle_use;
pair_BB_cmb_null_avg__ = pair_BB_cmb_null_avg__ / n_shuffle_use;
pair_CC_cmb_null_avg__ = pair_CC_cmb_null_avg__ / n_shuffle_use;
if flag_a_is_0==0; pair_a_cmb_null_std_ = sqrt(pair_a_cmb_null_std_/n_shuffle_use - pair_a_cmb_null_avg_.^2); end;
pair_A_cmb_null_std__ = sqrt(pair_A_cmb_null_std__/n_shuffle_use - pair_A_cmb_null_avg__.^2);
pair_BB_cmb_null_std__ = sqrt(pair_BB_cmb_null_std__/n_shuffle_use - pair_BB_cmb_null_avg__.^2);
pair_CC_cmb_null_std__ = sqrt(pair_CC_cmb_null_std__/n_shuffle_use - pair_CC_cmb_null_avg__.^2);
fname_mat = sprintf('%s/dolphin_%s_pair_cmb_s0000.mat',dir_mat,fname_infix);
tmp_ = load(fname_mat);
if flag_a_is_0==0; pair_a_cmb_ = tmp_.pair_a_cmb_; end;
pair_A_cmb__ = tmp_.pair_A_cmb__;
pair_BB_cmb__ = tmp_.pair_BB_cmb__;
pair_CC_cmb__ = tmp_.pair_CC_cmb__;
clear tmp_;
if flag_a_is_0==0; pair_a_cmb_Z_ = real((pair_a_cmb_ - pair_a_cmb_null_avg_)./pair_a_cmb_null_std_); end;
pair_A_cmb_Z__ = real((pair_A_cmb__ - pair_A_cmb_null_avg__)./pair_A_cmb_null_std__);
pair_BB_cmb_Z__ = real((pair_BB_cmb__ - pair_BB_cmb_null_avg__)./pair_BB_cmb_null_std__);
pair_CC_cmb_Z__ = real((pair_CC_cmb__ - pair_CC_cmb_null_avg__)./pair_CC_cmb_null_std__);
if flag_a_is_0==0; pair_a_cmb_nlp_ = -z_to_lp(pair_a_cmb_Z_); end;
pair_A_cmb_nlp__ = -z_to_lp(pair_A_cmb_Z__);
pair_BB_cmb_nlp__ = -z_to_lp(pair_BB_cmb_Z__);
pair_CC_cmb_nlp__ = -z_to_lp(pair_CC_cmb_Z__);
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_%s_pair_aABC_cmb_nlp_0',dir_jpg,fname_infix);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
if flag_a_is_0==0;
subplot(2,2,1);
bar(pair_a_cmb_nlp_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
title('-log p for drift pair_a_cmb_','Interpreter','none');
end;%if flag_a_is_0==0;
subplot(2,2,2);
colormap(colormap_beach());
imagesc(pair_A_cmb_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('-log p for interaction pair_A_cmb__','Interpreter','none');
subplot(2,2,3);
colormap(colormap_beach());
imagesc(pair_BB_cmb_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('-log p for SDE noise pair_BB_cmb__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(pair_CC_cmb_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
%colorbar;
title('-log p for observation noise pair_CC_cmb__','Interpreter','none');
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
end;%if flag_all;

%{
  %%%%%%%%;
  % step through each variable and examine the entries. ;
  %%%%%%%%;
  n_step = 1;
  dt_all_ = [];
  for naid=0:n_aid-1;
  tmp_index_aid_ = index_aid__{1+naid};
  tmp_age_ = age_(1+tmp_index_aid_);
  [tmp_age_sort_,tmp_index_] = sort(tmp_age_,'ascend'); tmp_index_ = tmp_index_-1;
  for nstep=0;n_step-1;
  dt_ = tmp_age_sort_(1+(1+nstep):end)-tmp_age_sort_(1:end-(1+nstep));
  dt_all_ = [dt_all_;dt_];
  end;%for nstep=0;n_step-1;
  end;%for naid=0:n_aid-1;
  %%%%%%%%;
  for nvar=0:n_var-1;
  c_ = {'r','k'};
  subplot(2,1,1); cla;
  hold on;
  age_base=0;
  for naid=0:n_aid-1;
  tmp_index_aid_ = index_aid__{1+naid};
  tmp_age_ = age_(1+tmp_index_aid_);
  [tmp_age_,tmp_index_] = sort(tmp_age_,'ascend'); tmp_index_ = tmp_index_-1;
  tmp_var_ = dat_imp__(1+tmp_index_aid_,1+nvar); tmp_var_ = tmp_var_(1+tmp_index_);
  plot(age_base+tmp_age_,tmp_var_,'.','Color',c_{1+mod(naid,2)});
  age_base = age_base + max(tmp_age_) + 1;
  end;%for naid=0:n_aid-1;
  hold off;
  subplot(2,1,2); cla;
  hold on;
  dt_all_ = [];
  dv_all_ = [];
  for naid=0:n_aid-1;
  tmp_index_aid_ = index_aid__{1+naid};
  tmp_age_ = age_(1+tmp_index_aid_);
  [tmp_age_,tmp_index_] = sort(tmp_age_,'ascend'); tmp_index_ = tmp_index_-1;
  tmp_var_ = dat_imp__(1+tmp_index_aid_,1+nvar); tmp_var_ = tmp_var_(1+tmp_index_);
  dt_ = diff(tmp_age_);
  dv_ = diff(tmp_var_);
  dt_all_ = [dt_all_;dt_];
  dv_all_ = [dv_all_;dv_];
  end;%for naid=0:n_aid-1;
  plot(dt_all_,dv_all_.^2,'k.');
  hold off;
  sgtitle(sprintf('%d %s',nvar,string_dat_name_{1+nvar}));
  drawnow();
  pause();
  end;%for nvar=0:n_var-1;
  %}

%{

  %%%%%%%%;
  % plot single variables and their drift. ;
  %%%%%%%%;
  [a_,b_] = dolphin_estimate_ab_0(aid_,age_,dat__);
  for nvar=0:n_var-1;
  c_ = {'r','k'};
  subplot(10,4,1+nvar); cla;
  hold on;
  age_base=0;
  for naid=0:n_aid-1;
  tmp_index_aid_ = index_aid__{1+naid};
  tmp_age_ = age_(1+tmp_index_aid_);
  [tmp_age_,tmp_index_] = sort(tmp_age_,'ascend'); tmp_index_ = tmp_index_-1;
  tmp_var_ = dat__(1+tmp_index_aid_,1+nvar); tmp_var_ = tmp_var_(1+tmp_index_);
  plot(0*age_base+tmp_age_,tmp_var_,'.','Color',c_{1+mod(naid,2)});
  age_base = age_base + max(tmp_age_) + 1;
  end;%for naid=0:n_aid-1;
  tmp_t_ = 1:50;
  plot(tmp_t_,b_(1+nvar) + tmp_t_*a_(1+nvar),'g-','LineWidth',3);
  hold off;
  title(sprintf('%d %s',nvar,string_dat_name_{1+nvar}));
  figbig;
  drawnow();
  end;%for nvar=0:n_var-1;


  %}
