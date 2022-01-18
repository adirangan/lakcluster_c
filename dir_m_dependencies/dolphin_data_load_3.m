clear;
setup_OptiPlex;
dir_trunk = '/home/rangan/dir_bcc/dir_dolphin';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);

flag_replot=1;

fname_data = sprintf('%s/Dolphin_Data_201918755_s2.xlsx',dir_trunk);
T_ = readtable(fname_data);
aid_ori_ = T_{:,1};
age_ori_ = T_{:,4};
nvar_base = 7;
dat_ori__ = T_{:,(1+nvar_base):end};
n_smp_ori = size(dat_ori__,1);
n_var_ori = size(dat_ori__,2);
string_dat_ori_name_ = cell(n_var_ori,1);
for nvar=0:n_var_ori-1;
string_dat_ori_name_{1+nvar} = T_.Properties.VariableNames{1+nvar_base+nvar};
end;%for nvar=0:n_var_ori-1;

dat_ret__ = dat_ori__;
mss_var_threshold = 500;
index_var_rmv_ = efind(sum(~isfinite(dat_ret__),1)>mss_var_threshold);
index_var_rmv_ = union(index_var_rmv_,efind(strcmp(string_dat_ori_name_,'RBCDist')));
index_var_rmv_ = union(index_var_rmv_,efind(strcmp(string_dat_ori_name_,'Albumin')));
index_var_ret_ = setdiff(0:n_var_ori-1,index_var_rmv_);
dat_ret__ = dat_ret__(:,1+index_var_ret_);
n_var_ret = size(dat_ret__,2);
string_dat_ret_name_ = cell(n_var_ret,1);
for nvar=0:n_var_ret-1;
string_dat_ret_name_{1+nvar} = string_dat_ori_name_{1+index_var_ret_(1+nvar)};
end;%for nvar=0:n_var_ret-1;

mss_smp_threshold = 7;
index_smp_ret_ = setdiff(0:n_smp_ori-1,efind(sum(~isfinite(dat_ret__),2)>mss_smp_threshold));
dat_ret__ = dat_ret__(1+index_smp_ret_,:);
n_smp_ret = size(dat_ret__,1);

fname_fig = sprintf('%s/dolphin_ori_vs_ret',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
subplot(2,2,1); stairs(cumsum(sort(sum(~isfinite(dat_ori__),1)))/numel(dat_ori__)); title('cumsum(sum(~isfinite(dat_ori__),1))','Interpreter','none');
subplot(2,2,2); stairs(cumsum(sort(sum(~isfinite(dat_ori__),2)))/numel(dat_ori__)); title('cumsum(sum(~isfinite(dat_ori__),2))','Interpreter','none');
subplot(2,2,3); stairs(cumsum(sort(sum(~isfinite(dat_ret__),1)))/numel(dat_ret__)); title('cumsum(sum(~isfinite(dat_ret__),1))','Interpreter','none');
subplot(2,2,4); stairs(cumsum(sort(sum(~isfinite(dat_ret__),2)))/numel(dat_ret__)); title('cumsum(sum(~isfinite(dat_ret__),2))','Interpreter','none');
figbig;
sgtitle('original vs retained')
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

aid_ret_ = aid_ori_(1+index_smp_ret_);
age_ret_ = age_ori_(1+index_smp_ret_);

aid_ = aid_ret_;
age_ = age_ret_;
dat_nrm__ = dat_ret__;
[n_smp,n_var] = size(dat_nrm__);
string_dat_name_ = string_dat_ret_name_;
prctile_threshold = 0.01;
for nvar=0:n_var-1;
tmp_index_ = efind(isfinite(dat_nrm__(:,1+nvar)));
tmp_lo = prctile(dat_nrm__(1+tmp_index_,1+nvar),100*prctile_threshold);
tmp_hi = prctile(dat_nrm__(1+tmp_index_,1+nvar),100*(1-prctile_threshold));
dat_nrm__(1+tmp_index_,1+nvar) = max(tmp_lo,min(tmp_hi,dat_nrm__(1+tmp_index_,1+nvar)));
end;%for nvar=0:n_var-1;

u_aid_ = unique(aid_);
n_aid = numel(u_aid_);
n_aid_ = zeros(n_aid,1);
index_aid__ = cell(n_aid,1);
for naid=0:n_aid-1;
tmp_index_aid_ = efind(aid_==u_aid_(1+naid));
index_aid__{1+naid} = tmp_index_aid_;
n_aid_(1+naid) = numel(index_aid__{1+naid});
tmp_dat_nrm__ = normalize(dat_nrm__(1+tmp_index_aid_,:),'zscore');
dat_nrm__(1+tmp_index_aid_,:) = tmp_dat_nrm__;
end;%for naid=0:n_aid-1;

dat__ = dat_nrm__;

fname_fig = sprintf('%s/dolphin_dt_scatter_A',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
c_ = colormap_beach(); n_c = size(c_,1);
n_plot_row = 5;
n_plot_col = 10;
nplot=0;
for nvar=0:24;%for nvar=0:n_var-1;
dt_all_ = [];
dv_all_ = [];
for naid=0:n_aid-1;
nc = max(0,min(n_c-1,floor(n_c*naid/n_aid)));
tmp_index_finite_ = efind(isfinite(dat__(1+index_aid__{1+naid},1+nvar)));
tmp_index_aid_ = index_aid__{1+naid}(1+tmp_index_finite_);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_,tmp_index_] = sort(tmp_age_,'ascend'); tmp_index_ = tmp_index_-1;
tmp_var_ = dat__(1+tmp_index_aid_,1+nvar); tmp_var_ = tmp_var_(1+tmp_index_);
dt_ = diff(tmp_age_);
dv_ = diff(tmp_var_);
dt_all_ = [dt_all_;dt_];
dv_all_ = [dv_all_;dv_];
subplot(n_plot_row,n_plot_col,1+nplot+0);
hold on;
plot(sqrt(dt_),dv_,'.','Color',c_(1+nc,:));
hold off;
end;%for naid=0:n_aid-1;
xlim_ = sqrt([0,5]);
ylim_ = prctile(dv_all_,[1,99]);
%subplot(n_plot_row,n_plot_col,1+nplot+1);
%plot(tmp_age_,tmp_var_,'.-','Color',c_(1+nc,:));
subplot(n_plot_row,n_plot_col,1+nplot+0);
xlim(xlim_); ylim(ylim_); title(sprintf('%s',string_dat_name_{1+nvar}));
subplot(n_plot_row,n_plot_col,1+nplot+1);
n_h = 32;
tmp_h__ = hist2d_0(sqrt(dt_all_),dv_all_,n_h,n_h,xlim_,ylim_);
tmp_e_ = linspace(ylim_(1),ylim_(2),n_h);
tmp_b_ = 1:n_h;
tmp_avg_ = zeros(n_h,1);
tmp_std_ = zeros(n_h,1);
for nh=0:n_h-1;
tmp_h__(:,1+nh) = tmp_h__(:,1+nh)/max(1,sum(tmp_h__(:,1+nh)));
tmp_avg_(1+nh) = dot(tmp_b_,tmp_h__(:,1+nh));
tmp_std_(1+nh) = sqrt(dot(tmp_b_.*tmp_b_,tmp_h__(:,1+nh)) - tmp_avg_(1+nh).^2);
end;%for nh=0:31;
hold on;
imagesc(tmp_h__,[0.02,0.20]);
colormap(colormap_beach());
plot(tmp_b_,tmp_avg_,'k.');
plot(tmp_b_,tmp_avg_+1.0*tmp_std_,'k.');
plot(tmp_b_,tmp_avg_-1.0*tmp_std_,'k.');
hold off;
xlim([1,n_h]); ylim([1,n_h]);
axisnotick;
title(sprintf('%s',string_dat_name_{1+nvar}));
nplot=nplot+2;
end;%for nvar=0:n_var-1;
figbig;
sgtitle('sqrt(dt) horizontal, dv vertical')
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/dolphin_dt_scatter_B',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
c_ = colormap_beach(); n_c = size(c_,1);
n_plot_row = 5;
n_plot_col = 10;
nplot=0;
for nvar=25:n_var-1;%for nvar=0:n_var-1;
dt_all_ = [];
dv_all_ = [];
for naid=0:n_aid-1;
nc = max(0,min(n_c-1,floor(n_c*naid/n_aid)));
tmp_index_finite_ = efind(isfinite(dat__(1+index_aid__{1+naid},1+nvar)));
tmp_index_aid_ = index_aid__{1+naid}(1+tmp_index_finite_);
tmp_age_ = age_(1+tmp_index_aid_);
[tmp_age_,tmp_index_] = sort(tmp_age_,'ascend'); tmp_index_ = tmp_index_-1;
tmp_var_ = dat__(1+tmp_index_aid_,1+nvar); tmp_var_ = tmp_var_(1+tmp_index_);
dt_ = diff(tmp_age_);
dv_ = diff(tmp_var_);
dt_all_ = [dt_all_;dt_];
dv_all_ = [dv_all_;dv_];
subplot(n_plot_row,n_plot_col,1+nplot+0);
hold on;
plot(sqrt(dt_),dv_,'.','Color',c_(1+nc,:));
hold off;
end;%for naid=0:n_aid-1;
xlim_ = sqrt([0,5]);
ylim_ = prctile(dv_all_,[1,99]);
%subplot(n_plot_row,n_plot_col,1+nplot+1);
%plot(tmp_age_,tmp_var_,'.-','Color',c_(1+nc,:));
subplot(n_plot_row,n_plot_col,1+nplot+0);
xlim(xlim_); ylim(ylim_); title(sprintf('%s',string_dat_name_{1+nvar}));
subplot(n_plot_row,n_plot_col,1+nplot+1);
n_h = 32;
tmp_h__ = hist2d_0(sqrt(dt_all_),dv_all_,n_h,n_h,xlim_,ylim_);
tmp_e_ = linspace(ylim_(1),ylim_(2),n_h);
tmp_b_ = 1:n_h;
tmp_avg_ = zeros(n_h,1);
tmp_std_ = zeros(n_h,1);
for nh=0:n_h-1;
tmp_h__(:,1+nh) = tmp_h__(:,1+nh)/max(1,sum(tmp_h__(:,1+nh)));
tmp_avg_(1+nh) = dot(tmp_b_,tmp_h__(:,1+nh));
tmp_std_(1+nh) = sqrt(dot(tmp_b_.*tmp_b_,tmp_h__(:,1+nh)) - tmp_avg_(1+nh).^2);
end;%for nh=0:31;
hold on;
imagesc(tmp_h__,[0.02,0.20]);
colormap(colormap_beach());
plot(tmp_b_,tmp_avg_,'k.');
plot(tmp_b_,tmp_avg_+1.0*tmp_std_,'k.');
plot(tmp_b_,tmp_avg_-1.0*tmp_std_,'k.');
hold off;
xlim([1,n_h]); ylim([1,n_h]);
axisnotick;
title(sprintf('%s',string_dat_name_{1+nvar}));
nplot=nplot+2;
end;%for nvar=0:n_var-1;
figbig;
sgtitle('sqrt(dt) horizontal, dv vertical')
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

fname_mat = sprintf('%s/dolphin_data_s0000.mat',dir_mat);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
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
a_ = zeros(n_var,1);
A__ = zeros(n_var,n_var);
%%%%%%%%;
index_mss_ = efind(~isfinite(dat__));
BB__ = []; CC__ = [];
%%%%%%%%;
% Now impute data using A. ;
%%%%%%%%;
dat_imp__ = dolphin_impute_aA_0(aid_,age_,dat__,a_,A__,index_mss_);
%%%%%%%%;
% Now estimate initial B,C. ;
%%%%%%%%;
[BB__,CC__,l2_R__,sum_1,sum_dt,sum_dtdt,sum_DDj__,sum_DDjdt__] = dolphin_estimate_BC_from_aA_0(aid_,age_,dat_imp__,a_,A__,n_step);
[BB_crude__,CC_crude__] = dolphin_estimate_BC_from_aA_crude_0(aid_,age_,dat_imp__,a_,A__,n_step);
disp(sprintf(' %% BB__ vs BB_crude__: %0.16f',fnorm(BB__ - BB_crude__)/fnorm(BB__)));
disp(sprintf(' %% CC__ vs CC_crude__: %0.16f',fnorm(CC__ - CC_crude__)/fnorm(CC__)));
BB__ = BB_crude__; CC__ = CC_crude__; %<-- use crude method. ;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_BCD_0',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
colormap(colormap_beach());
subplot(2,3,1+0);
imagesc(BB__); axis image; title('BB__','Interpreter','none');
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_);xtickangle(90);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
subplot(2,3,2+0);
imagesc(CC__); axis image; title('CC__','Interpreter','none');
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_);xtickangle(90);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
subplot(2,3,3+0);
imagesc(log10(l2_R__)); axis image; title('log10(l2_R__)','Interpreter','none');
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_);xtickangle(90);
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
subplot(2,3,1+3);
plot(diag(BB__),'ko'); title('BB__','Interpreter','none');
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_);xtickangle(90);
subplot(2,3,2+3);
plot(diag(CC__),'ko'); title('CC__','Interpreter','none');
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_);xtickangle(90);
subplot(2,3,3+3);
plot(diag(log10(l2_R__)),'ko'); title('log10(l2_R__)','Interpreter','none');
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_);xtickangle(90);
figbig;
sgtitle('initial BB__, CC__ and log10(l2_R__)','Interpreter','none')
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Now estimate a_ and A__ from BB__ and CC__. ;
%%%%%%%%;
dt_lim_ = [0;max(dt_all_)];
[a_,A__,L] = dolphin_estimate_aA_from_BC_0(aid_,age_,dat_imp__,BB__,CC__,dt_lim_,n_step);
disp(sprintf(' %% initial: negative-log-likelihood %0.16f',L));
%%%%%%%%;
% Now iterate a few times. ;
%%%%%%%%;
n_iteration = 4;
L_ = zeros(n_iteration+1,1);
a__ = cell(1+n_iteration,1);
A___ = cell(1+n_iteration,1);
BB___ = cell(1+n_iteration,1);
CC___ = cell(1+n_iteration,1);
niteration=0;
L_old = L;
L_(1+niteration) = L_old;
a__{1+niteration} = a_;
A___{1+niteration} = A__;
BB___{1+niteration} = BB__;
CC___{1+niteration} = CC__;
flag_continue=1;
while (flag_continue);
% Re-impute data using A. ;
dat_imp__ = dolphin_impute_aA_0(aid_,age_,dat__,a_,A__,index_mss_);
% Re-estimate B,C. ;
[BB__,CC__,l2_R__,sum_1,sum_dt,sum_dtdt,sum_DDj__,sum_DDjdt__] = dolphin_estimate_BC_from_aA_0(aid_,age_,dat_imp__,a_,A__,n_step);
[BB_crude__,CC_crude__] = dolphin_estimate_BC_from_aA_crude_0(aid_,age_,dat_imp__,a_,A__,n_step);
disp(sprintf(' %% BB__ vs BB_crude__: %0.16f',fnorm(BB__ - BB_crude__)/fnorm(BB__)));
disp(sprintf(' %% CC__ vs CC_crude__: %0.16f',fnorm(CC__ - CC_crude__)/fnorm(CC__)));
BB__ = BB_crude__; CC__ = CC_crude__; %<-- use crude method. ;
% Re-estimate a_ and A__ from BB__ and CC__. ;
[a_,A__,L] = dolphin_estimate_aA_from_BC_0(aid_,age_,dat_imp__,BB__,CC__,dt_lim_,n_step);
L_new = L; L_(1+niteration+1)=L;
a__{1+niteration+1} = a_;
A___{1+niteration+1} = A__;
BB___{1+niteration+1} = BB__;
CC___{1+niteration+1} = CC__;
disp(sprintf(' %% iteration %d: negative-log-likelihood %0.16f',niteration,L));
flag_continue=0;
niteration=niteration+1;
if (niteration<n_iteration & fnorm(L_old-L_new)/fnorm(L_old)>1e-3); flag_continue=1; end;
L_old = L_new;
end;%while (flag_continue);
%%%%%%%%;
save(fname_mat,'n_step','n_iteration','a_','A__','BB__','CC__','L','a__','A___','BB___','CC___','L_');
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_aABC_0',dir_jpg);
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
colorbar;
title('interaction A__','Interpreter','none');
subplot(2,2,3);
colormap(colormap_beach());
imagesc(BB__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
title('SDE noise BB__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(CC__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
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
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% Now repeat for several permuted trials. ;
%%%%%%%%;
n_shuffle = 32; n_step = 1; n_iteration = 4;
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_data_s%.4d.mat',dir_mat,nshuffle);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle);
[a_null_,A_null__,BB_null__,CC_null__,L_null_,niteration_null,a_null__,A_null___,BB_null___,CC_null___] = dolphin_estimate_aABC_0(aid_,age_,dat_prm__,n_step,n_iteration);
save(fname_mat ...
,'n_step' ... 
,'n_iteration' ... 
,'nshuffle' ... 
,'a_null_' ... 
,'A_null__' ... 
,'BB_null__' ... 
,'CC_null__' ... 
,'L_null_' ...
,'niteration_null' ...
,'a_null__' ... 
,'A_null___' ... 
,'BB_null___' ... 
,'CC_null___' ... 
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_mat));
end;%if ( exist(fname_mat,'file'));
end;%for nshuffle=1:n_shuffle;

%%%%%%%%;
% collect null distribution. ;
%%%%%%%%;
a_null_avg_ = zeros(n_var,1);
A_null_avg__ = zeros(n_var,n_var);
BB_null_avg__ = zeros(n_var,n_var);
CC_null_avg__ = zeros(n_var,n_var);
a_null_std_ = zeros(n_var,1);
A_null_std__ = zeros(n_var,n_var);
BB_null_std__ = zeros(n_var,n_var);
CC_null_std__ = zeros(n_var,n_var);
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_data_s%.4d.mat',dir_mat,nshuffle);
tmp_ = load(fname_mat);
[~,tmp_index] = min(tmp_.L_null_(1:1+tmp_.niteration_null)); tmp_index = tmp_index-1;
a_null_avg_ = a_null_avg_ + tmp_.a_null__{1+tmp_index};
A_null_avg__ = A_null_avg__ + tmp_.A_null___{1+tmp_index};
BB_null_avg__ = BB_null_avg__ + tmp_.BB_null___{1+tmp_index};
CC_null_avg__ = CC_null_avg__ + tmp_.CC_null___{1+tmp_index};
a_null_std_ = a_null_std_ + tmp_.a_null__{1+tmp_index}.^2;
A_null_std__ = A_null_std__ + tmp_.A_null___{1+tmp_index}.^2;
BB_null_std__ = BB_null_std__ + tmp_.BB_null___{1+tmp_index}.^2;
CC_null_std__ = CC_null_std__ + tmp_.CC_null___{1+tmp_index}.^2;
clear tmp_;
end;%for nshuffle=1:n_shuffle;
a_null_avg_ = a_null_avg_ / n_shuffle;
A_null_avg__ = A_null_avg__ / n_shuffle;
BB_null_avg__ = BB_null_avg__ / n_shuffle;
CC_null_avg__ = CC_null_avg__ / n_shuffle;
a_null_std_ = sqrt(a_null_std_/n_shuffle - a_null_avg_.^2);
A_null_std___ = sqrt(A_null_std__/n_shuffle - A_null_avg__.^2);
BB_null_std___ = sqrt(BB_null_std__/n_shuffle - BB_null_avg__.^2);
CC_null_std___ = sqrt(CC_null_std__/n_shuffle - CC_null_avg__.^2);
fname_mat = sprintf('%s/dolphin_data_s0000.mat',dir_mat);
tmp_ = load(fname_mat);
[~,tmp_index] = min(tmp_.L_(1:1+3)); tmp_index = tmp_index-1;
a_ = tmp_.a__{1+tmp_index};
A__ = tmp_.A___{1+tmp_index};
BB__ = tmp_.BB___{1+tmp_index};
CC__ = tmp_.CC___{1+tmp_index};
clear tmp_;
a_Z_ = real((a_ - a_null_avg_)./a_null_std_);
A_Z__ = real((A__ - A_null_avg__)./A_null_std__);
BB_Z__ = real((BB__ - BB_null_avg__)./BB_null_std__);
CC_Z__ = real((CC__ - CC_null_avg__)./CC_null_std__);
a_nlp_ = -z_to_lp(a_Z_);
A_nlp__ = -z_to_lp(A_Z__);
BB_nlp__ = -z_to_lp(BB_Z__);
CC_nlp__ = -z_to_lp(CC_Z__);
   
fname_fig = sprintf('%s/dolphin_aABC_nlp_0',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
subplot(2,2,1);
bar(a_nlp_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
title('-log p for drift a_','Interpreter','none');
subplot(2,2,2);
colormap(colormap_beach());
imagesc(A_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
title('-log p for interaction A__','Interpreter','none');
subplot(2,2,3);
colormap(colormap_beach());
imagesc(BB_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
title('-log p for SDE noise BB__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(CC_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
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
% Now step through every pair of variables and check for pairwise causality. ;
%%%%%%%%;
n_step = 1; n_iteration = 4;
fname_mat = sprintf('%s/dolphin_data_pair_cmb_s0000.mat',dir_mat);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
pair_a_cmb_ = zeros(n_var,1);
pair_A_cmb__ = zeros(n_var,n_var);
pair_BB_cmb__ = zeros(n_var,n_var);
pair_CC_cmb__ = zeros(n_var,n_var);
n_pair = n_var*(n_var-1)/2;
for nvar_0=0:n_var-1;
for nvar_1=nvar_0+1:n_var-1;
disp(sprintf(' %% nvar_0 %d nvar_1 %d',nvar_0,nvar_1));
pair_dat__ = dat__(:,1+[nvar_0,nvar_1]);
[pair_a_,pair_A__,pair_BB__,pair_CC__,pair_L_,pair_niteration,pair_a__,pair_A___,pair_BB___,pair_CC___] = dolphin_estimate_aABC_0(aid_,age_,pair_dat__,n_step,n_iteration);
[~,tmp_index] = min(pair_L_(1:1+pair_niteration)); tmp_index = tmp_index-1;
pair_a_cmb_(1+nvar_0) = pair_a_cmb_(1+nvar_0) + pair_a__{1+tmp_index}(1+0)/n_pair;
pair_a_cmb_(1+nvar_1) = pair_a_cmb_(1+nvar_1) + pair_a__{1+tmp_index}(1+1)/n_pair;
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
save(fname_mat ...
,'n_step' ... 
,'n_iteration' ... 
,'pair_a_cmb_' ...
,'pair_A_cmb__' ... 
,'pair_BB_cmb__' ... 
,'pair_CC_cmb__' ... 
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_mat));
end;%if ( exist(fname_mat,'file'));

fname_mat = sprintf('%s/dolphin_data_pair_cmb_s0000.mat',dir_mat);
load(fname_mat);
fname_fig = sprintf('%s/dolphin_pair_aABC_cmb_0',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
subplot(2,2,1);
bar(pair_a_cmb_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
title('drift pair_a_cmb_','Interpreter','none');
subplot(2,2,2);
colormap(colormap_beach());
imagesc(pair_A_cmb__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
title('interaction pair_A_cmb__','Interpreter','none');
subplot(2,2,3);
colormap(colormap_beach());
imagesc(pair_BB_cmb__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
title('SDE noise pair_BB_cmb__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(pair_CC_cmb__); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
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

%%%%%%%%;
% Now repeat for several permuted trials. ;
%%%%%%%%;
n_shuffle = 32; n_step = 1; n_iteration = 4;
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_data_pair_cmb_s%.4d.mat',dir_mat,nshuffle);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle);
%%%%%%%%;
pair_a_cmb_ = zeros(n_var,1);
pair_A_cmb__ = zeros(n_var,n_var);
pair_BB_cmb__ = zeros(n_var,n_var);
pair_CC_cmb__ = zeros(n_var,n_var);
n_pair = n_var*(n_var-1)/2;
for nvar_0=0:n_var-1;
for nvar_1=nvar_0+1:n_var-1;
disp(sprintf(' %% nvar_0 %d nvar_1 %d',nvar_0,nvar_1));
pair_dat_prm__ = dat_prm__(:,1+[nvar_0,nvar_1]);
[pair_a_,pair_A__,pair_BB__,pair_CC__,pair_L_,pair_niteration,pair_a__,pair_A___,pair_BB___,pair_CC___] = dolphin_estimate_aABC_0(aid_,age_,pair_dat_prm__,n_step,n_iteration);
[~,tmp_index] = min(pair_L_(1:1+pair_niteration)); tmp_index = tmp_index-1;
pair_a_cmb_(1+nvar_0) = pair_a_cmb_(1+nvar_0) + pair_a__{1+tmp_index}(1+0)/n_pair;
pair_a_cmb_(1+nvar_1) = pair_a_cmb_(1+nvar_1) + pair_a__{1+tmp_index}(1+1)/n_pair;
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
save(fname_mat ...
,'n_step' ... 
,'n_iteration' ...
,'nshuffle' ...
,'pair_a_cmb_' ...
,'pair_A_cmb__' ... 
,'pair_BB_cmb__' ... 
,'pair_CC_cmb__' ... 
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not loading',fname_mat));
end;%if ( exist(fname_mat,'file'));
end;%for nshuffle=1:n_shuffle;

%%%%%%%%;
% collect null distribution. ;
%%%%%%%%;
pair_a_cmb_null_avg_ = zeros(n_var,1);
pair_A_cmb_null_avg__ = zeros(n_var,n_var);
pair_BB_cmb_null_avg__ = zeros(n_var,n_var);
pair_CC_cmb_null_avg__ = zeros(n_var,n_var);
pair_a_cmb_null_std_ = zeros(n_var,1);
pair_A_cmb_null_std__ = zeros(n_var,n_var);
pair_BB_cmb_null_std__ = zeros(n_var,n_var);
pair_CC_cmb_null_std__ = zeros(n_var,n_var);
n_shuffle_use = 32;
for nshuffle=1:n_shuffle_use;
fname_mat = sprintf('%s/dolphin_data_pair_cmb_s%.4d.mat',dir_mat,nshuffle);
tmp_ = load(fname_mat);
disp(sprintf(' %% norms: pair_a_ %0.2f pair_A__ %0.2f pair_BB__ %0.2f pair_CC__ %0.2f',fnorm(tmp_.pair_a_cmb_),fnorm(tmp_.pair_A_cmb__),fnorm(tmp_.pair_BB_cmb__),fnorm(tmp_.pair_CC_cmb__)));
pair_a_cmb_null_avg_ = pair_a_cmb_null_avg_ + tmp_.pair_a_cmb_;
pair_A_cmb_null_avg__ = pair_A_cmb_null_avg__ + tmp_.pair_A_cmb__;
pair_BB_cmb_null_avg__ = pair_BB_cmb_null_avg__ + tmp_.pair_BB_cmb__;
pair_CC_cmb_null_avg__ = pair_CC_cmb_null_avg__ + tmp_.pair_CC_cmb__;
pair_a_cmb_null_std_ = pair_a_cmb_null_std_ + tmp_.pair_a_cmb_.^2;
pair_A_cmb_null_std__ = pair_A_cmb_null_std__ + tmp_.pair_A_cmb__.^2;
pair_BB_cmb_null_std__ = pair_BB_cmb_null_std__ + tmp_.pair_BB_cmb__.^2;
pair_CC_cmb_null_std__ = pair_CC_cmb_null_std__ + tmp_.pair_CC_cmb__.^2;
clear tmp_;
end;%for nshuffle=1:n_shuffle_use;
pair_a_cmb_null_avg_ = pair_a_cmb_null_avg_ / n_shuffle_use;
pair_A_cmb_null_avg__ = pair_A_cmb_null_avg__ / n_shuffle_use;
pair_BB_cmb_null_avg__ = pair_BB_cmb_null_avg__ / n_shuffle_use;
pair_CC_cmb_null_avg__ = pair_CC_cmb_null_avg__ / n_shuffle_use;
pair_a_cmb_null_std_ = sqrt(pair_a_cmb_null_std_/n_shuffle_use - pair_a_cmb_null_avg_.^2);
pair_A_cmb_null_std___ = sqrt(pair_A_cmb_null_std__/n_shuffle_use - pair_A_cmb_null_avg__.^2);
pair_BB_cmb_null_std___ = sqrt(pair_BB_cmb_null_std__/n_shuffle_use - pair_BB_cmb_null_avg__.^2);
pair_CC_cmb_null_std___ = sqrt(pair_CC_cmb_null_std__/n_shuffle_use - pair_CC_cmb_null_avg__.^2);
fname_mat = sprintf('%s/dolphin_data_pair_cmb_s0000.mat',dir_mat);
tmp_ = load(fname_mat);
pair_a_cmb_ = tmp_.pair_a_cmb_;
pair_A_cmb__ = tmp_.pair_A_cmb__;
pair_BB_cmb__ = tmp_.pair_BB_cmb__;
pair_CC_cmb__ = tmp_.pair_CC_cmb__;
clear tmp_;
pair_a_cmb_Z_ = real((pair_a_cmb_ - pair_a_cmb_null_avg_)./pair_a_cmb_null_std_);
pair_A_cmb_Z__ = real((pair_A_cmb__ - pair_A_cmb_null_avg__)./pair_A_cmb_null_std__);
pair_BB_cmb_Z__ = real((pair_BB_cmb__ - pair_BB_cmb_null_avg__)./pair_BB_cmb_null_std__);
pair_CC_cmb_Z__ = real((pair_CC_cmb__ - pair_CC_cmb_null_avg__)./pair_CC_cmb_null_std__);
pair_a_cmb_nlp_ = -z_to_lp(pair_a_cmb_Z_);
pair_A_cmb_nlp__ = -z_to_lp(pair_A_cmb_Z__);
pair_BB_cmb_nlp__ = -z_to_lp(pair_BB_cmb_Z__);
pair_CC_cmb_nlp__ = -z_to_lp(pair_CC_cmb_Z__);

fname_fig = sprintf('%s/dolphin_pair_aABC_cmb_nlp_0',dir_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);
clf;
subplot(2,2,1);
bar(pair_a_cmb_nlp_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
title('-log p for drift pair_a_cmb_','Interpreter','none');
subplot(2,2,2);
colormap(colormap_beach());
imagesc(pair_A_cmb_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
title('-log p for interaction pair_A_cmb__','Interpreter','none');
subplot(2,2,3);
colormap(colormap_beach());
imagesc(pair_BB_cmb_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
title('-log p for SDE noise pair_BB_cmb__','Interpreter','none');
subplot(2,2,4);
colormap(colormap_beach());
imagesc(pair_CC_cmb_nlp__,[0,6]); axis image;
set(gca,'YTick',1:n_var,'YTickLabel',string_dat_name_);
set(gca,'XTick',1:n_var,'XTickLabel',string_dat_name_); xtickangle(90);
set(gca,'TickLength',[0,0]);
set(gca,'fontsize',6);
colorbar;
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

