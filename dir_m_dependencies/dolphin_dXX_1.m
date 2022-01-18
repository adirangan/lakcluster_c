clear;
%setup_OptiPlex; dir_trunk = '/home/rangan/dir_bcc/dir_dolphin';
setup_access1; dir_trunk = '/data/rangan/dir_bcc/dir_dolphin';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);

flag_replot=1;
flag_recalc=0;

flag_drift0=1;
flag_a_is_0=0;
flag_center=1;
dXX = 0;

fname_infix = sprintf('aid%.2dd%dc%d',dXX,flag_drift0,flag_center);

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

%{

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
  % RBCDist and Albumin look strange. ;
  %%%%%%%%;
  figure(1);clf;
  c_ = {'r','k'};
  nvar = efind(strcmp(string_dat_ori_name_,'RBCDist'));
  subplot(2,1,1); cla;
  hold on;
  age_base=0;
  for naid=0:n_aid-1;
  tmp_index_aid_ = index_aid__{1+naid};
  tmp_age_ = age_(1+tmp_index_aid_);
  [tmp_age_,tmp_index_] = sort(tmp_age_,'ascend'); tmp_index_ = tmp_index_-1;
  tmp_var_ = dat_ori__(1+tmp_index_aid_,1+nvar); tmp_var_ = tmp_var_(1+tmp_index_);
  plot(age_base+tmp_age_,tmp_var_,'.','Color',c_{1+mod(naid,2)});
  age_base = age_base + max(tmp_age_) + 1;
  end;%for naid=0:n_aid-1;
  hold off;
  title(string_dat_ori_name_{1+nvar});
  %%%%%%%%;
  nvar = efind(strcmp(string_dat_ori_name_,'Albumin'));
  subplot(2,1,2); cla;
  hold on;
  age_base=0;
  for naid=0:n_aid-1;
  tmp_index_aid_ = index_aid__{1+naid};
  tmp_age_ = age_(1+tmp_index_aid_);
  [tmp_age_,tmp_index_] = sort(tmp_age_,'ascend'); tmp_index_ = tmp_index_-1;
  tmp_var_ = dat_ori__(1+tmp_index_aid_,1+nvar); tmp_var_ = tmp_var_(1+tmp_index_);
  plot(age_base+tmp_age_,tmp_var_,'.','Color',c_{1+mod(naid,2)});
  age_base = age_base + max(tmp_age_) + 1;
  end;%for naid=0:n_aid-1;
  hold off;
  title(string_dat_ori_name_{1+nvar});
  %%%%%%%%;
  figbig;
  print('-depsc',sprintf('%s/dolphin_ori_RBCDist_Albumin_0.eps',dir_jpg));
  print('-djpeg',sprintf('%s/dolphin_ori_RBCDist_Albumin_0.jpg',dir_jpg));

  %%%%%%%%;
  % looks like we should log-normalize SED60 and GGT. ;
  %%%%%%%%;
  figure(1);clf;
  nvar = efind(strcmp(string_dat_ori_name_,'NRBC'));
  subplot(2,3,1); hist(dat_ori__(:,1+nvar),32); title(string_dat_ori_name_{1+nvar});
  subplot(2,3,4); hist(log(1+dat_ori__(:,1+nvar)),32); title(sprintf('log(%s)',string_dat_ori_name_{1+nvar}));
  nvar = efind(strcmp(string_dat_ori_name_,'SED60'));
  subplot(2,3,2); hist(dat_ori__(:,1+nvar),32); title(string_dat_ori_name_{1+nvar});
  subplot(2,3,5); hist(log(1+dat_ori__(:,1+nvar)),32); title(sprintf('log(%s)',string_dat_ori_name_{1+nvar}));
  nvar = efind(strcmp(string_dat_ori_name_,'GGT'));
  subplot(2,3,3); hist(dat_ori__(:,1+nvar),32); title(string_dat_ori_name_{1+nvar});
  subplot(2,3,6); hist(log(1+dat_ori__(:,1+nvar)),32); title(sprintf('log(%s)',string_dat_ori_name_{1+nvar}));
  print('-depsc',sprintf('%s/dolphin_ori_hist_NRBC_SED60_GGT_0.eps',dir_jpg));
  print('-djpeg',sprintf('%s/dolphin_ori_hist_NRBC_SED60_GGT_0.jpg',dir_jpg));
  %}

nvar = efind(strcmp(string_dat_ori_name_,'SED60'));
dat_ori__(:,1+nvar) = log(1+dat_ori__(:,1+nvar));
string_dat_ori_name_{1+nvar} = sprintf('log(1+%s)',string_dat_ori_name_{1+nvar});
nvar = efind(strcmp(string_dat_ori_name_,'GGT'));
dat_ori__(:,1+nvar) = log(1+dat_ori__(:,1+nvar));
string_dat_ori_name_{1+nvar} = sprintf('log(1+%s)',string_dat_ori_name_{1+nvar});

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

[a_drift_,b_drift_] = dolphin_estimate_ab_0(aid_,age_,dat_nrm__);

n_step = 1; dt_all_ = dolphin_dt_all_0(aid_,age_,n_step);

u_aid_ = unique(aid_);
n_aid = numel(u_aid_);
n_aid_ = zeros(n_aid,1);
index_aid__ = cell(n_aid,1);
for naid=0:n_aid-1;
tmp_index_aid_ = efind(aid_==u_aid_(1+naid));
index_aid__{1+naid} = tmp_index_aid_;
n_aid_(1+naid) = numel(index_aid__{1+naid});
tmp_dat_nrm__ = dat_nrm__(1+tmp_index_aid_,:);
tmp_age_ = age_(1+tmp_index_aid_);
if (flag_drift0==1); tmp_dat_nrm__ = tmp_dat_nrm__ - (ones(numel(tmp_age_),1)*reshape(b_drift_,[1,n_var]) + tmp_age_*reshape(a_drift_,[1,n_var])); end;
if (flag_drift0==0); tmp_dat_nrm__ = tmp_dat_nrm__; end;
if (flag_drift0==0 & flag_center==1); tmp_dat_nrm__ = normalize(tmp_dat_nrm__,'zscore'); end;
if (flag_drift0==0 & flag_center==0); tmp_dat_nrm__ = tmp_dat_nrm__; end;
dat_nrm__(1+tmp_index_aid_,:) = tmp_dat_nrm__;
end;%for naid=0:n_aid-1;
if (flag_drift0==1 & flag_center==1); dat_nrm__ = normalize(dat_nrm__,'zscore'); end;
if (flag_drift0==1 & flag_center==0); dat_nrm__ = dat_nrm__; end;

n_step = 1; relative_variation = dolphin_relative_variation_0(aid_,age_,dat_nrm__,n_step);

%%%%%%%%;
% select dolphin. ;
%%%%%%%%;
if dXX>0; 
naid = efind(u_aid_==dXX);
tmp_index_aid_ = efind(aid_==u_aid_(1+naid));
aid_ = aid_(1+tmp_index_aid_);
age_ = age_(1+tmp_index_aid_);
dat_nrm__ = dat_nrm__(1+tmp_index_aid_,:);
end;%if dXX>0; 

dat__ = dat_nrm__;

fname_mat = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix);
fname_tmp = sprintf('%s/dolphin_%s_s0000.tmp',dir_mat,fname_infix);
if (flag_recalc | (~exist(fname_mat,'file') & ~exist(fname_tmp,'file')));
disp(sprintf(' %% %s not found, creating',fname_mat));
save(fname_tmp,'fname_mat');
%%%%%%%%;
index_mss_ = efind(~isfinite(dat__));
BB__ = []; CC__ = [];
%%%%%%%%;
n_step = 1; n_iteration = 4;
%%%%%%%%;
if flag_a_is_0==0;
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_0(aid_,age_,dat__,n_step,n_iteration);
save(fname_mat,'n_step','n_iteration','a_','A__','BB__','CC__','L_','niteration','a__','A___','BB___','CC___');
end;%if flag_a_is_0==0;
if flag_a_is_0==1;
[A__,BB__,CC__,L_,niteration,A___,BB___,CC___] = dolphin_estimate_ABC_0(aid_,age_,dat__,n_step,n_iteration);
save(fname_mat,'n_step','n_iteration','A__','BB__','CC__','L_','niteration','A___','BB___','CC___');
end;%if flag_a_is_0==1;
%%%%%%%%;
delete(fname_tmp);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
end;%if ( exist(fname_mat,'file'));

fname_mat = sprintf('%s/dolphin_%s_s0000.mat',dir_mat,fname_infix);
if ( exist(fname_mat,'file'));
load(fname_mat);
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
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% Now repeat for several permuted trials. ;
%%%%%%%%;
n_shuffle = 32; n_step = 1; n_iteration = 4;
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_%s_s%.4d.mat',dir_mat,fname_infix,nshuffle);
fname_tmp = sprintf('%s/dolphin_%s_s%.4d.tmp',dir_mat,fname_infix,nshuffle);
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
fname_mat = sprintf('%s/dolphin_%s_s%.4d.mat',dir_mat,fname_infix,nshuffle);
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
fname_mat = sprintf('%s/dolphin_%s_s%.4d.mat',dir_mat,fname_infix,nshuffle);
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
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% Now repeat for several permuted trials. ;
%%%%%%%%;
n_shuffle = 32; n_step = 1; n_iteration = 4;
for nshuffle=1:n_shuffle;
fname_mat = sprintf('%s/dolphin_%s_pair_cmb_s%.4d.mat',dir_mat,fname_infix,nshuffle);
fname_tmp = sprintf('%s/dolphin_%s_pair_cmb_s%.4d.tmp',dir_mat,fname_infix,nshuffle);
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
fname_mat = sprintf('%s/dolphin_%s_s%.4d.mat',dir_mat,fname_infix,nshuffle);
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
fname_mat = sprintf('%s/dolphin_%s_pair_cmb_s%.4d.mat',dir_mat,fname_infix,nshuffle);
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
