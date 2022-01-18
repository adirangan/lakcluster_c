function dolphin_dXX_5(dXX,n_shuffle,age_percentile_range_);
%{
note that age-range is 0.64 --> 54.35. ;
try: ;

dXX = 0; n_shuffle = 256;
dolphin_dXX_5(dXX,n_shuffle,[ 0, 50]);
dolphin_dXX_5(dXX,n_shuffle,[50,100]);
dolphin_dXX_5(dXX,n_shuffle);

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
flag_constrain_CC=1;
%dXX = 0;

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
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_1(aid_,age_,dat__,n_step,n_iteration,flag_constrain_CC);
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
fname_fig = sprintf('%s/dolphin_%s_aABC_1',dir_jpg,fname_infix);
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

%{
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
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_1(aid_0_,age_0_,dat_0__,n_step,n_iteration,flag_constrain_CC);
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
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_1(aid_1_,age_1_,dat_1__,n_step,n_iteration,flag_constrain_CC);
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
%}

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
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_1(aid_,age_,dat_prm__,n_step,n_iteration,flag_constrain_CC);
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

disp('returning'); return;

