%%%%%%%%;
% step through each of the 44 variables. ;
%%%%%%%%;
nf=0;
%dXX = 0;
%age_percentile_range_ = [];
%sex_constrain = [];
%%%%%%%%;
platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (~exist('dir_trunk','var')); 
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_dolphin',string_root);
end;%if (~exist('dir_trunk','var')); 
dir_jpg_base = sprintf('%s/dir_jpg',dir_trunk);
%dir_mat = sprintf('%s/dir_mat',dir_trunk);
%dir_shuffle_mat = sprintf('%s/dir_shuffle_mat',dir_trunk);
%%%%%%%%;
if (~exist('flag_replot','var')); flag_replot=0; end;
%%%%%%%%;
fname_data = sprintf('%s/Dolphin_Data_201918755_s2.xlsx',dir_trunk);
T_ = readtable(fname_data);
%%%%%%%%;
aid_ori_ = T_{:,1};
u_aid_ori_ = unique(aid_ori_); n_u_aid_ori = numel(u_aid_ori_);
n_u_aid_ori_ = zeros(n_u_aid_ori,1);
index_aid_ori__ = cell(n_u_aid_ori,1);
for nu_aid_ori=0:n_u_aid_ori-1;
index_aid_ori_ = efind(aid_ori_==u_aid_ori_(1+nu_aid_ori));
n_u_aid_ori_(1+nu_aid_ori) = numel(index_aid_ori_);
index_aid_ori__{1+nu_aid_ori} = index_aid_ori_;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
%%%%%%%%;
spe_ori_ = T_{:,3}; %<-- retain TT species. ;
u_spe_ori_ = unique(spe_ori_); n_u_spe_ori = numel(u_spe_ori_);
n_u_spe_ori_ = zeros(n_u_spe_ori,1);
index_spe_ori__ = cell(n_u_spe_ori,1);
for nu_spe_ori=0:n_u_spe_ori-1;
index_spe_ori_ = efind(strcmp(spe_ori_,u_spe_ori_{1+nu_spe_ori}));
n_u_spe_ori_(1+nu_spe_ori) = numel(index_spe_ori_);
index_spe_ori__{1+nu_spe_ori} = index_spe_ori_;
end;%for nu_spe_ori=0:n_u_spe_ori-1;
isTT_ori_ = cellfun(@(x)strcmp(x,'TT'),spe_ori_);
index_smp_TT_ = efind(isTT_ori_);
%%%%%%%%;
lab_ori_ = T_{:,7}; %<-- may need to ignore lab_code, as some dolphins have multiple lab_codes. ;
lab_let_ori_ = cast(lab_ori_ - min(lab_ori_) + 'A','char');
u_lab_ori_ = unique(lab_ori_); n_u_lab_ori = numel(u_lab_ori_);
n_u_lab_ori_ = zeros(n_u_lab_ori,1);
index_lab_ori__ = cell(n_u_lab_ori,1);
for nu_lab_ori=0:n_u_lab_ori-1;
index_lab_ori_ = efind(lab_ori_==u_lab_ori_(1+nu_lab_ori));
n_u_lab_ori_(1+nu_lab_ori) = numel(index_lab_ori_);
index_lab_ori__{1+nu_lab_ori} = index_lab_ori_;
end;%for nu_lab_ori=0:n_u_lab_ori-1;
u_lab_color__ = [ ...
  0.5 0.1 0.8 ...
; 0.8 0.1 0.5 ...
; 1.0 0.9 0.9 ...
; 0.0 0.0 0.0 ...
; 0.0 0.0 0.0 ...
; 0.9 0.9 1.0 ...
; 0.0 0.0 0.0 ...
; 0.0 0.0 0.0 ...
; 0.0 0.0 0.0 ...
; 0.0 0.0 0.0 ...
];
lab_ori_color__ = zeros(size(lab_ori_,1),3);
for nu_lab_ori=0:n_u_lab_ori-1;
index_lab_ori_ = index_lab_ori__{1+nu_lab_ori};
n_u_lab_ori = n_u_lab_ori_(1+nu_lab_ori);
lab_ori_color__(1+index_lab_ori_,:) = repmat(u_lab_color__(1+nu_lab_ori,:),[n_u_lab_ori,1]);
end;%for nu_lab_ori=0:n_u_lab_ori-1;
%%%%%%%%;
sex_ori_ = T_{:,2};
isM_ori_ = cellfun(@(x)strcmp(x,'M'),sex_ori_);
isF_ori_ = cellfun(@(x)strcmp(x,'F'),sex_ori_);
age_ori_ = T_{:,4};
age_ori_M_ = age_ori_(1+efind(isM_ori_));
age_ori_F_ = age_ori_(1+efind(isF_ori_));
disp(sprintf(' %% A: %d prctile [ 5 15 20 30 65 70 85 95 ] %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f ',numel(age_ori_),prctile(age_ori_,[ 5 15 20 30 65 70 85 95 ])));
disp(sprintf(' %% M: %d prctile [ 5 15 20 30 65 70 85 95 ] %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f ',sum(isM_ori_),prctile(age_ori_M_,[ 5 15 20 30 65 70 85 95 ])));
disp(sprintf(' %% F: %d prctile [ 5 15 20 30 65 70 85 95 ] %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f ',sum(isF_ori_),prctile(age_ori_F_,[ 5 15 20 30 65 70 85 95 ])));
h_age_ = 1:65;
h_age_M_ = hist(age_ori_M_,h_age_);
h_age_F_ = hist(age_ori_F_,h_age_);
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_age_hist_MF',dir_jpg_base);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1); clf; figmed; figbeach;
subplot(1,1,1);
bar(h_age_,[h_age_F_ ; h_age_M_],'grouped','BarWidth',1.5);
xlabel('age');
ylabel('number');
legend({'F','M'});
title('Histogram of ages for Female (blue) and Male (red) dolphins');
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
nvar_base = 7;
dat_ori__ = T_{:,(1+nvar_base):end};
n_smp_ori = size(dat_ori__,1);
n_var_ori = size(dat_ori__,2);
string_dat_ori_name_ = cell(n_var_ori,1);
for nvar=0:n_var_ori-1;
string_dat_ori_name_{1+nvar} = T_.Properties.VariableNames{1+nvar_base+nvar};
end;%for nvar=0:n_var_ori-1;
%%%%%%%%;

%%%%%%%%;
% Set defaults for variables. ;
%%%%%%%%;
flag_zero_as_missing_ = ones(n_var_ori,1);
flag_zero_as_missing_(1+efind(strcmp(string_dat_ori_name_,'NRBC'))) = 0;
flag_zero_as_missing_(1+efind(strcmp(string_dat_ori_name_,'Monocytes'))) = 0;
flag_zero_as_missing_(1+efind(strcmp(string_dat_ori_name_,'ACMonocytes'))) = 0;
flag_zero_as_missing_(1+efind(strcmp(string_dat_ori_name_,'UricAcid'))) = 0;
flag_zero_as_missing_(1+efind(strcmp(string_dat_ori_name_,'Bilirubin'))) = 0;
flag_zero_as_missing_(1+efind(strcmp(string_dat_ori_name_,'SED60'))) = 0;
flag_zero_as_missing_(1+efind(strcmp(string_dat_ori_name_,'GFR'))) = 0;
%%%%;
flag_pseudocount_ = zeros(n_var_ori,1);
val_pseudocount_ = zeros(n_var_ori,1);
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'NRBC'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'NRBC'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Monocytes'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Monocytes'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'ACMonocytes'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'ACMonocytes'))) = 0.01;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'EOS'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'EOS'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Creatinine'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Creatinine'))) = 0.1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'UricAcid'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'UricAcid'))) = 0.1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Potassium'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Potassium'))) = 0.1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'CO2'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'CO2'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Protein'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Protein'))) = 0.1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Albumin'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Albumin'))) = 0.1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Calcium'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Calcium'))) = 0.1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'InorgPhos'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'InorgPhos'))) = 0.1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'AlkPhos'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'AlkPhos'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'ALT'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'ALT'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'GGT'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'GGT'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Bilirubin'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Bilirubin'))) = 0.1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Triglyceride'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Triglyceride'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Iron'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Iron'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'CPK'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'CPK'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'SED60'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'SED60'))) = 1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Mg'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'Mg'))) = 0.1;
flag_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'GFR'))) = 1;
val_pseudocount_(1+efind(strcmp(string_dat_ori_name_,'GFR'))) = 4e-4;
%%%%;
flag_inv_transform_ = zeros(n_var_ori,1);
flag_inv_transform_(1+efind(strcmp(string_dat_ori_name_,'GFR'))) = 1;
%%%%;
flag_upb_threshold_ = zeros(n_var_ori,1);
val_upb_threshold_ = zeros(n_var_ori,1);
flag_upb_threshold_(1+efind(strcmp(string_dat_ori_name_,'RBCDist'))) = 1;
val_upb_threshold_(1+efind(strcmp(string_dat_ori_name_,'RBCDist'))) = 40;
%%%%;
flag_lob_threshold_ = zeros(n_var_ori,1);
val_lob_threshold_ = zeros(n_var_ori,1);
flag_lob_threshold_(1+efind(strcmp(string_dat_ori_name_,'Albumin'))) = 1;
val_lob_threshold_(1+efind(strcmp(string_dat_ori_name_,'Albumin'))) = 2;
flag_lob_threshold_(1+efind(strcmp(string_dat_ori_name_,'GFR'))) = 1;
val_lob_threshold_(1+efind(strcmp(string_dat_ori_name_,'GFR'))) = 100;
%%%%;

%%%%%%%%;
% setting Creatinine to be a function of GFR. ;
% Afterwards we no longer need bounds on GFR. ;
%%%%%%%%;
nvar_ori_0 = efind(strcmp(string_dat_ori_name_,'GFR'));
nvar_ori_1 = efind(strcmp(string_dat_ori_name_,'Creatinine'));
for nsmp_ori=0:n_smp_ori-1;
x = dat_ori__(1+nsmp_ori,1+nvar_ori_0);
y = dat_ori__(1+nsmp_ori,1+nvar_ori_1);
if ( isfinite(y)); x = 300/max(1e-12,y.^(2/sqrt(3))); end;
if (~isfinite(y) & isfinite(x) & (x>0)); y = round((300/max(1e-12,x)).^(sqrt(3)/2),2); end;
dat_ori__(1+nsmp_ori,1+nvar_ori_0) = x;
dat_ori__(1+nsmp_ori,1+nvar_ori_1) = y;
tmp_x_(1+nsmp_ori) = x;
tmp_y_(1+nsmp_ori) = y;
end;%for nsmp_ori=0:n_smp_ori-1;
%%%%%%%%;
%flag_lob_threshold_(1+efind(strcmp(string_dat_ori_name_,'GFR'))) = 0;
%val_lob_threshold_(1+efind(strcmp(string_dat_ori_name_,'GFR'))) = 0;
%%%%%%%%;
flag_plot=0;
if flag_plot;
nvar_ori_0 = efind(strcmp(string_dat_ori_name_,'GFR')); nvar_ori_1 = efind(strcmp(string_dat_ori_name_,'Creatinine')); tmp_index_use_ = efind(isTT_ori_); x_ = dat_ori__(1+tmp_index_use_,1+nvar_ori_0); y_ = dat_ori__(1+tmp_index_use_,1+nvar_ori_1); figure;clf;figbig; subplot(1,2,1); loglog(x_,y_,'ro',x_,round((300./x_).^(sqrt(3)/2),2),'gx'); xlabel(string_dat_ori_name_{1+nvar_ori_0}); ylabel(string_dat_ori_name_{1+nvar_ori_1}); subplot(1,2,2); plot(res_nrm__(:,1+nvar_ori_0),res_nrm__(:,1+nvar_ori_1),'ko'); xlabel(string_dat_ori_name_{1+nvar_ori_0}); ylabel(string_dat_ori_name_{1+nvar_ori_1});
end;%if flag_plot;

%%%%%%%%;
% Determine skewness and select log-transform. ;
% Limiting to TT species. ;
%%%%%%%%;
flag_log_transform_ = zeros(n_var_ori,1);
age_lob =  0; age_upb = 55;
for nvar_ori=0:n_var_ori-1;
string_dat_ori_name = string_dat_ori_name_{1+nvar_ori};
tmp_aid_ori_ = aid_ori_(1+index_smp_TT_);
tmp_lab_ori_ = lab_ori_(1+index_smp_TT_);
tmp_lab_ori_color__ = lab_ori_color__(1+index_smp_TT_,:);
tmp_sex_ori_ = sex_ori_(1+index_smp_TT_);
tmp_age_ori_ = age_ori_(1+index_smp_TT_);
tmp_dat_ori_ = dat_ori__(1+index_smp_TT_,1+nvar_ori);
flag_zero_as_missing = flag_zero_as_missing_(1+nvar_ori);
if flag_zero_as_missing; tmp_dat_ori_(1+efind(tmp_dat_ori_==0)) = NaN; end;%if flag_zero_as_missing;
flag_upb_threshold = flag_upb_threshold_(1+nvar_ori);
if flag_upb_threshold;
val_upb_threshold = val_upb_threshold_(1+nvar_ori);
tmp_dat_ori_(1+efind(tmp_dat_ori_> val_upb_threshold)) = NaN;
end;%if flag_upb_threshold;
flag_lob_threshold = flag_lob_threshold_(1+nvar_ori);
if flag_lob_threshold;
val_lob_threshold = val_lob_threshold_(1+nvar_ori);
tmp_dat_ori_(1+efind(tmp_dat_ori_< val_lob_threshold)) = NaN;
end;%if flag_lob_threshold;
flag_inv_transform = flag_inv_transform_(1+nvar_ori);
if flag_inv_transform; tmp_dat_ori_ = 1./tmp_dat_ori_; end;
tmp_index_use_ = efind( isfinite(tmp_age_ori_) & isfinite(tmp_dat_ori_) & (tmp_age_ori_>=age_lob) & (tmp_age_ori_<=age_upb) );
p_use = numel(tmp_index_use_)/numel(index_smp_TT_);
tmp_aid_use_ = tmp_aid_ori_(1+tmp_index_use_);
tmp_lab_use_ = tmp_lab_ori_(1+tmp_index_use_);
tmp_lab_color_use__ = tmp_lab_ori_color__(1+tmp_index_use_,:);
tmp_sex_use_ = tmp_sex_ori_(1+tmp_index_use_);
tmp_age_use_ = tmp_age_ori_(1+tmp_index_use_);
tmp_dat_use_ = tmp_dat_ori_(1+tmp_index_use_);
tmp_dat_p01 = prctile(tmp_dat_use_,01);
tmp_dat_p99 = prctile(tmp_dat_use_,99);
ddat = 0;
flag_pseudocount = flag_pseudocount_(1+nvar_ori);
if flag_pseudocount; ddat = val_pseudocount_(1+nvar_ori); end;%if flag_pseudocount;
tmp_log_use_ = log(ddat + tmp_dat_use_);
tmp_log_p01 = prctile(tmp_log_use_,01);
tmp_log_p99 = prctile(tmp_log_use_,99);
tmp_index_dat_use_mid_ = efind( (tmp_dat_use_>=tmp_dat_p01) & (tmp_dat_use_<=tmp_dat_p99) );
skewness_dat  = skewness(normalize(tmp_dat_use_(1+tmp_index_dat_use_mid_)));
skewness_log  = skewness(normalize(tmp_log_use_(1+tmp_index_dat_use_mid_)));
if (abs(skewness_dat)<=abs(skewness_log)); flag_log_transform = 0; end;
if (abs(skewness_dat)> abs(skewness_log)); flag_log_transform = 1; end;
disp(sprintf('nvar_ori %d: %s --> %0.2f %0.2f --> %d',nvar_ori,string_dat_ori_name,skewness_dat,skewness_log,flag_log_transform));
flag_log_transform_(1+nvar_ori) = flag_log_transform;
%%%%%%%%%%%%%%%%;
fname_fig = sprintf('%s/fig_log_transform_nvar%.2d_%s_FIGA',dir_jpg_base,nvar_ori,string_dat_ori_name);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
sgtitle(sprintf('%s --> %0.2f %0.2f --> %d',string_dat_ori_name,skewness_dat,skewness_log,flag_log_transform));
ylim_ = [tmp_dat_p01,tmp_dat_p99]; ylim_ = mean(ylim_) + 0.5*1.5*diff(ylim_)*[-1,+1];
subplot(2,3,[1,2]); plot(tmp_age_use_,tmp_dat_use_,'k.');
xlim([0,55]); ylim(ylim_); xlabel('age'); ylabel(string_dat_ori_name); grid on;
subplot(2,3,3); hist(tmp_dat_use_,linspace(tmp_dat_p01,tmp_dat_p99,128)); title('hist');
xlim(ylim_);
ylim_ = [tmp_log_p01,tmp_log_p99]; ylim_ = mean(ylim_) + 0.5*1.5*diff(ylim_)*[-1,+1];
subplot(2,3,[4,5]); plot(tmp_age_use_,tmp_log_use_,'k.');
xlim([0,55]); ylim(ylim_); xlabel('age'); ylabel(string_dat_ori_name); grid on;
subplot(2,3,6); hist(tmp_log_use_,linspace(tmp_log_p01,tmp_log_p99,128)); title('hist log');
xlim(ylim_);
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%%%%%%%%%;
string_dat_ori_name_title = sprintf('%s',string_dat_ori_name);
if flag_log_transform;
string_dat_ori_name_title = sprintf('log(%s)',string_dat_ori_name);
tmp_dat_use_ = tmp_log_use_;
tmp_dat_p01 = tmp_log_p01;
tmp_dat_p99 = tmp_log_p99;
end;%if flag_log_transform;
ylim_ = [tmp_dat_p01,tmp_dat_p99]; ylim_ = mean(ylim_) + 0.5*1.5*diff(ylim_)*[-1,+1];
tmp_index_dat_use_mid_ = efind( (tmp_dat_use_>=tmp_dat_p01) & (tmp_dat_use_<=tmp_dat_p99) );
fname_fig = sprintf('%s/fig_chebfit_nvar%.2d_%s_FIGA',dir_jpg_base,nvar_ori,string_dat_ori_name);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
n_degree = 1;
[c_obj,x_value_use_,y_value_use_] = chebfit_safe_1(tmp_age_use_(1+tmp_index_dat_use_mid_),tmp_dat_use_(1+tmp_index_dat_use_mid_),n_degree);
figure(1+nf);nf=nf+1;clf;figbig;
sgtitle(sprintf('%s p %0.2f',string_dat_ori_name_title,p_use));
p_row = 3; p_col = 4; ns=0;
markersize_use = 4;
fontsize_use = 8;
%%%%;
subplot(p_row,p_col,1+ns+[0:2]);ns=ns+3;
hold on;
tmp_index_ = efind(strcmp(tmp_sex_use_,'F'));
plot(tmp_age_use_(1+tmp_index_),tmp_dat_use_(1+tmp_index_),'ko','MarkerSize',markersize_use,'MarkerFaceColor',[0,0.5,1.0]);
tmp_index_ = efind(strcmp(tmp_sex_use_,'M'));
plot(tmp_age_use_(1+tmp_index_),tmp_dat_use_(1+tmp_index_),'ko','MarkerSize',markersize_use,'MarkerFaceColor',[1.0,0.5,0]);
plot(sort(tmp_age_use_),c_obj(sort(tmp_age_use_)),'r-','LineWidth',3);
hold off;
xlim([age_lob,age_upb]); xlabel('age');
ylim(ylim_);
ylabel(string_dat_ori_name_title,'Interpreter','none');
title('sex');
grid on;
%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
hist(tmp_dat_use_,linspace(tmp_dat_p01,tmp_dat_p99,128));
xlim(ylim_);
title('histogram');
%%%%;
subplot(p_row,p_col,1+1*p_col+[0:p_col-1]);
hold on;
na=0; na_ = [0]; tmp_c_ = [1,0,0];
str_symbol_ = {'o','^','s','p','h'}; nsymbol=0;
for naid=0:n_u_aid_ori-1;
tmp_index_ = efind(tmp_aid_use_==u_aid_ori_(1+naid));
if (numel(tmp_index_)>0);
tmp_age_min = min(tmp_age_use_(1+tmp_index_));
tmp_age_max = max(tmp_age_use_(1+tmp_index_));
str_symbol = str_symbol_{1+nsymbol};
plot( ...
 1+na+(tmp_age_use_(1+tmp_index_)-tmp_age_min) ...
,tmp_dat_use_(1+tmp_index_) ...
,str_symbol ...
,'MarkerEdgeColor','k' ...
,'MarkerSize',markersize_use ...
,'MarkerFaceColor',tmp_c_ ...
);
na = na + (tmp_age_max - tmp_age_min) + 1;
na_ = [na_,na];
nsymbol = nsymbol+1; if (nsymbol>4); nsymbol=0; end;
tmp_c_ = 1-tmp_c_;
end;%if (numel(tmp_index_)>0);
end;%for naid=0:n_u_aid_ori-1;
hold off;
xlim([0,1+na]); xlabel('age'); set(gca,'XTick',na_);
ylim(ylim_);
ylabel(string_dat_ori_name_title,'Interpreter','none');
title('aid');
grid on;
%%%%;
subplot(p_row,p_col,1+2*p_col+[0:p_col-1]);
hold on;
na=0; na_ = [0]; tmp_c_ = [1,0,0];
str_symbol_ = {'o','^','s','p','h'}; nsymbol=0;
for naid=0:n_u_aid_ori-1;
tmp_index_ = efind(tmp_aid_use_==u_aid_ori_(1+naid));
if (numel(tmp_index_)>0);
tmp_age_min = min(tmp_age_use_(1+tmp_index_));
tmp_age_max = max(tmp_age_use_(1+tmp_index_));
str_symbol = str_symbol_{1+nsymbol};
scatter( ...
 1+na+(tmp_age_use_(1+tmp_index_)-tmp_age_min) ...
,tmp_dat_use_(1+tmp_index_) ...
,markersize_use ...
,tmp_lab_color_use__(1+tmp_index_,:) ...
,'filled' ...
);
na = na + (tmp_age_max - tmp_age_min) + 1;
na_ = [na_,na];
nsymbol = nsymbol+1; if (nsymbol>4); nsymbol=0; end;
tmp_c_ = 1-tmp_c_;
end;%if (numel(tmp_index_)>0);
end;%for naid=0:n_u_aid_ori-1;
hold off;
xlim([0,1+na]); xlabel('age'); set(gca,'XTick',na_);
ylim(ylim_);
ylabel(string_dat_ori_name_title,'Interpreter','none');
title('lab');
grid on;
%%%%;
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for nvar_ori=0:n_var_ori-1;
%%%%%%%%;

%%%%%%%%;
% Determine slope and calculate thresholded normalized residual. ;
% Limiting to TT species. ;
% res_ is segregated by sex. ;
%%%%%%%%;
index_smp_TT_isM_ = intersect(index_smp_TT_,efind(isM_ori_));
index_smp_TT_isF_ = intersect(index_smp_TT_,efind(isF_ori_));
n_degree = 1;
polyfit_M__ = zeros(n_var_ori,1+n_degree);
polyfit_F__ = zeros(n_var_ori,1+n_degree);
che_lob =  5; %<-- do not include very young dolphins when determining baseline. ;
che_upb = 55;
prctile_upb = 99;
prctile_lob =  1;
dat_nrm__ = NaN*ones(size(dat_ori__));
res_nrm__ = NaN*ones(size(dat_ori__));
for nvar_ori=0:n_var_ori-1;
string_dat_ori_name = string_dat_ori_name_{1+nvar_ori};
for nn=0:1;
if nn==0; flag_isM=1; flag_isF=0; index_smp_TT_isX_ = index_smp_TT_isM_; end;
if nn==1; flag_isM=0; flag_isF=1; index_smp_TT_isX_ = index_smp_TT_isF_; end;
tmp_aid_ori_ = aid_ori_(1+index_smp_TT_isX_);
tmp_sex_ori_ = sex_ori_(1+index_smp_TT_isX_);
tmp_age_ori_ = age_ori_(1+index_smp_TT_isX_);
tmp_dat_ori_ = dat_ori__(1+index_smp_TT_isX_,1+nvar_ori);
flag_zero_as_missing = flag_zero_as_missing_(1+nvar_ori);
if flag_zero_as_missing; tmp_dat_ori_(1+efind(tmp_dat_ori_==0)) = NaN; end;%if flag_zero_as_missing;
flag_upb_threshold = flag_upb_threshold_(1+nvar_ori);
if flag_upb_threshold;
val_upb_threshold = val_upb_threshold_(1+nvar_ori);
tmp_dat_ori_(1+efind(tmp_dat_ori_> val_upb_threshold)) = NaN;
end;%if flag_upb_threshold;
flag_lob_threshold = flag_lob_threshold_(1+nvar_ori);
if flag_lob_threshold;
val_lob_threshold = val_lob_threshold_(1+nvar_ori);
tmp_dat_ori_(1+efind(tmp_dat_ori_< val_lob_threshold)) = NaN;
end;%if flag_lob_threshold;
flag_inv_transform = flag_inv_transform_(1+nvar_ori);
if flag_inv_transform; tmp_dat_ori_ = 1./tmp_dat_ori_; end;
flag_log_transform = flag_log_transform_(1+nvar_ori);
if flag_log_transform;
ddat = 0;
flag_pseudocount = flag_pseudocount_(1+nvar_ori);
if flag_pseudocount; ddat = val_pseudocount_(1+nvar_ori); end;%if flag_pseudocount;
tmp_dat_ori_ = log(ddat + tmp_dat_ori_);
end;%if flag_log_transform;
tmp_dat_lob = prctile(tmp_dat_ori_,prctile_lob);
tmp_dat_upb = prctile(tmp_dat_ori_,prctile_upb);
tmp_index_use_ = efind( 1 ...
&isfinite(tmp_age_ori_) ...
&isfinite(tmp_dat_ori_) ...
&(tmp_age_ori_>=che_lob) ...
&(tmp_age_ori_<=che_upb) ...
&(tmp_dat_ori_>=tmp_dat_lob) ...
&(tmp_dat_ori_<=tmp_dat_upb) ...
);
tmp_age_use_ = tmp_age_ori_(1+tmp_index_use_);
tmp_dat_use_ = tmp_dat_ori_(1+tmp_index_use_);
polyfit_X_ = polyfit_safe_0(tmp_age_use_,tmp_dat_use_,n_degree);
if flag_isM; polyfit_M__(1+nvar_ori,:) = polyfit_X_; end;
if flag_isF; polyfit_F__(1+nvar_ori,:) = polyfit_X_; end;
tmp_index_use_ = efind( 1 ...
&isfinite(tmp_age_ori_) ...
&isfinite(tmp_dat_ori_) ...
);
tmp_age_use_ = tmp_age_ori_(1+tmp_index_use_);
tmp_dat_use_ = tmp_dat_ori_(1+tmp_index_use_);
tmp_res_use_ = tmp_dat_use_ - polyval(polyfit_X_,tmp_age_use_);
tmp_res_cut_ = max(prctile(tmp_res_use_,prctile_lob),min(prctile(tmp_res_use_,prctile_upb),tmp_res_use_));
tmp_dat_cut_ = max(prctile(tmp_dat_use_,prctile_lob),min(prctile(tmp_dat_use_,prctile_upb),tmp_dat_use_));
tmp_res_nrm_ = normalize(tmp_res_cut_);
tmp_dat_nrm_ = normalize(tmp_dat_cut_);
res_nrm__(1+index_smp_TT_isX_(1+tmp_index_use_),1+nvar_ori) = tmp_res_nrm_;
dat_nrm__(1+index_smp_TT_isX_(1+tmp_index_use_),1+nvar_ori) = tmp_dat_nrm_;
end;%for nn=0:1;
end;%for nvar_ori=0:n_var_ori-1;
%%%%%%%%;

%%%%%%%%;
fname_fig = sprintf('%s/fig_res_plot_FIGA',dir_jpg_base);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 5; p_col = 9; ns=0;
for nvar_ori=0:n_var_ori-1;
string_dat_ori_name = string_dat_ori_name_{1+nvar_ori};
string_dat_ori_name_title = sprintf('%s',string_dat_ori_name);
flag_log_transform = flag_log_transform_(1+nvar_ori);
if flag_log_transform;
string_dat_ori_name_title = sprintf('log(%s)',string_dat_ori_name);
end;%if flag_log_transform;
subplot(p_row,p_col,1+ns);ns=ns+1;
hold on;
plot(age_ori_(1+index_smp_TT_isM_),res_nrm__(1+index_smp_TT_isM_,1+nvar_ori),'r.');
plot(age_ori_(1+index_smp_TT_isF_),res_nrm__(1+index_smp_TT_isF_,1+nvar_ori),'g.');
hold off;
xlim([0,55]);
title(string_dat_ori_name_title);
end;%for nvar_ori=0:n_var_ori-1;
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
 
%%%%%%%%;
fname_fig = sprintf('%s/fig_res_hist_FIGA',dir_jpg_base);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 5; p_col = 9; ns=0;
for nvar_ori=0:n_var_ori-1;
string_dat_ori_name = string_dat_ori_name_{1+nvar_ori};
string_dat_ori_name_title = sprintf('%s',string_dat_ori_name);
flag_log_transform = flag_log_transform_(1+nvar_ori);
if flag_log_transform;
string_dat_ori_name_title = sprintf('log(%s)',string_dat_ori_name);
end;%if flag_log_transform;
h_x_ = (-3:0.1:+3);
hM_ = hist(res_nrm__(1+index_smp_TT_isM_,1+nvar_ori),h_x_);
hF_ = hist(res_nrm__(1+index_smp_TT_isF_,1+nvar_ori),h_x_);
subplot(p_row,p_col,1+ns);ns=ns+1;
hold on;
bar(h_x_,+hF_);
bar(h_x_,-hM_);
hold off;
title(string_dat_ori_name_title);
end;%for nvar_ori=0:n_var_ori-1;
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% res_ind_ is further segregated by individual (for determining baseline). ;
%%%%%%%%;
res_ind_nrm__ = res_nrm__;
polyfit_res_ind_pnv___ = zeros(1+n_degree,n_u_aid_ori,n_var_ori);
for nu_aid_ori=0:n_u_aid_ori-1;
u_aid_ori = u_aid_ori_(1+nu_aid_ori);
index_aid_ori_ = index_aid_ori__{1+nu_aid_ori};
tmp_age_ori_ = age_ori_(1+index_aid_ori_);
for nvar_ori=0:n_var_ori-1;
tmp_res_nrm_ = res_nrm__(1+index_aid_ori_,1+nvar_ori);
tmp_index_use_ = efind( 1 ...
&isfinite(tmp_age_ori_) ...
&(tmp_age_ori_>=che_lob) ...
&(tmp_age_ori_<=che_upb) ...
&isfinite(tmp_res_nrm_) ...
);
if numel(tmp_index_use_)> (1+n_degree);
tmp_age_use_ = tmp_age_ori_(1+tmp_index_use_);
tmp_res_use_ = tmp_res_nrm_(1+tmp_index_use_);
polyfit_res_ind_ = polyfit_safe_0(tmp_age_use_,tmp_res_use_,n_degree);
polyfit_res_ind_pnv___(:,1+nu_aid_ori,1+nvar_ori) = polyfit_res_ind_;
tmp_index_use_ = efind( 1 ...
&isfinite(tmp_age_ori_) ...
&isfinite(tmp_res_nrm_) ...
);
tmp_res_ind_use_ = tmp_res_nrm_(1+tmp_index_use_) - polyval(polyfit_res_ind_,tmp_age_ori_(1+tmp_index_use_));
tmp_res_ind_nrm_ = normalize(tmp_res_ind_use_);
res_ind_nrm__(1+index_aid_ori_(1+tmp_index_use_),1+nvar_ori) = tmp_res_ind_nrm_;
end;%if numel(tmp_index_use_)> (1+n_degree);
end;%for nvar_ori=0:n_var_ori-1;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
%%%%%%%%;

%%%%%%%%;
fname_fig = sprintf('%s/fig_res_ind_plot_FIGA',dir_jpg_base);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 5; p_col = 9; ns=0;
for nvar_ori=0:n_var_ori-1;
string_dat_ori_name = string_dat_ori_name_{1+nvar_ori};
string_dat_ori_name_title = sprintf('%s',string_dat_ori_name);
flag_log_transform = flag_log_transform_(1+nvar_ori);
if flag_log_transform;
string_dat_ori_name_title = sprintf('log(%s)',string_dat_ori_name);
end;%if flag_log_transform;
subplot(p_row,p_col,1+ns);ns=ns+1;
hold on;
plot(age_ori_(1+index_smp_TT_isM_),res_ind_nrm__(1+index_smp_TT_isM_,1+nvar_ori),'r.');
plot(age_ori_(1+index_smp_TT_isF_),res_ind_nrm__(1+index_smp_TT_isF_,1+nvar_ori),'g.');
hold off;
xlim([0,55]);
title(string_dat_ori_name_title);
end;%for nvar_ori=0:n_var_ori-1;
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
 
%%%%%%%%;
fname_fig = sprintf('%s/fig_res_ind_hist_FIGA',dir_jpg_base);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 5; p_col = 9; ns=0;
for nvar_ori=0:n_var_ori-1;
string_dat_ori_name = string_dat_ori_name_{1+nvar_ori};
string_dat_ori_name_title = sprintf('%s',string_dat_ori_name);
flag_log_transform = flag_log_transform_(1+nvar_ori);
if flag_log_transform;
string_dat_ori_name_title = sprintf('log(%s)',string_dat_ori_name);
end;%if flag_log_transform;
h_x_ = (-3:0.1:+3);
hM_ = hist(res_ind_nrm__(1+index_smp_TT_isM_,1+nvar_ori),h_x_);
hF_ = hist(res_ind_nrm__(1+index_smp_TT_isF_,1+nvar_ori),h_x_);
subplot(p_row,p_col,1+ns);ns=ns+1;
hold on;
bar(h_x_,+hF_);
bar(h_x_,-hM_);
hold off;
title(string_dat_ori_name_title);
end;%for nvar_ori=0:n_var_ori-1;
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%%%%%;
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now examine increments. ;
%%%%%%%%;
for nvar_ori=0:n_var_ori-1;
string_dat_ori_name = string_dat_ori_name_{1+nvar_ori};
string_dat_ori_name_title = sprintf('%s',string_dat_ori_name);
flag_log_transform = flag_log_transform_(1+nvar_ori);
if flag_log_transform;
string_dat_ori_name_title = sprintf('log(%s)',string_dat_ori_name);
end;%if flag_log_transform;
%%%%%%%%;
n_h_x = 32;
n_h_y = 12;
nstep_ = 1:3; n_nstep = numel(nstep_);
dt_all__ = cell(n_nstep,1); dv_all__ = cell(n_nstep,1); tab_ = zeros(n_nstep,1);
for nnstep=0:n_nstep-1;
dt_all__{1+nnstep} = zeros(n_smp_ori,1); dv_all__{1+nnstep} = zeros(n_smp_ori,1);
end;%for nnstep=0:n_nstep-1;
for nu_aid_ori=0:n_u_aid_ori-1;
index_aid_ori_ = index_aid_ori__{1+nu_aid_ori};
tmp_age_ori_ = age_ori_(1+index_aid_ori_);
tmp_res_ori_ = res_nrm__(1+index_aid_ori_,1+nvar_ori);
tmp_index_finite_ = efind( isfinite(tmp_age_ori_) & isfinite(tmp_res_ori_) );
tmp_age_use_ = tmp_age_ori_(1+tmp_index_finite_);
tmp_res_use_ = tmp_res_ori_(1+tmp_index_finite_);
[tmp_age_srt_,tmp_index_srt_] = sort(tmp_age_use_,'ascend'); tmp_index_srt_ = tmp_index_srt_-1;
tmp_res_srt_ = tmp_res_use_(1+tmp_index_srt_);
for nnstep=0:n_nstep-1;
nstep = nstep_(1+nnstep);
tab = tab_(1+nnstep);
dt_ = tmp_age_srt_(1+nstep:end) - tmp_age_srt_(1:end-nstep);
dv_ = tmp_res_srt_(1+nstep:end) - tmp_res_srt_(1:end-nstep);
tab_add = numel(dt_);
dt_all__{1+nnstep}(1+tab+[0:tab_add-1]) = dt_;
dv_all__{1+nnstep}(1+tab+[0:tab_add-1]) = dv_;
tab_(1+nnstep) = tab_(1+nnstep) + tab_add;
end;%for nnstep=0:n_nstep-1;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
for nnstep=0:n_nstep-1;
dt_all__{1+nnstep} = dt_all__{1+nnstep}(1:tab_(1+nnstep)); dv_all__{1+nnstep} = dv_all__{1+nnstep}(1:tab_(1+nnstep));
end;%for nnstep=0:n_nstep-1;
tab_sum = sum(tab_);
tab_csum_ = cumsum([0;tab_]);
dt_all_all_ = zeros(tab_sum,1); dv_all_all_ = zeros(tab_sum,1);
for nnstep=0:n_nstep-1;
dt_all_all_(1+tab_csum_(1+nnstep)+[0:tab_(1+nnstep)-1]) = dt_all__{1+nnstep};
dv_all_all_(1+tab_csum_(1+nnstep)+[0:tab_(1+nnstep)-1]) = dv_all__{1+nnstep};
end;%for nnstep=0:n_nstep-1;
xlim_ = [0,sqrt(max(dt_all_all_))];
ylim_ = prctile(dv_all_all_,[1,99]);
tmp_h___ = zeros(n_h_y,n_h_x,n_nstep+1);
tmp_avg__ = zeros(n_h_x,n_nstep+1);
tmp_std__ = zeros(n_h_x,n_nstep+1);
for nnstep=0:n_nstep+1-1;
if (nnstep< n_nstep); tmp_dt_all_ = dt_all__{1+nnstep}; tmp_dv_all_ = dv_all__{1+nnstep}; end;
if (nnstep==n_nstep); tmp_dt_all_ = dt_all_all_; tmp_dv_all_ = dv_all_all_; end;
tmp_h__ = hist2d_0(sqrt(tmp_dt_all_),tmp_dv_all_,n_h_x,n_h_y,xlim_,ylim_);
tmp_e_ = linspace(ylim_(1),ylim_(2),n_h_x);
tmp_b_ = 1:n_h_y;
tmp_avg_ = zeros(n_h_x,1);
tmp_std_ = zeros(n_h_x,1);
for nh=0:n_h_x-1;
tmp_h__(:,1+nh) = tmp_h__(:,1+nh)/max(1,sum(tmp_h__(:,1+nh)));
tmp_avg_(1+nh) = dot(tmp_b_,tmp_h__(:,1+nh));
tmp_std_(1+nh) = sqrt(dot(tmp_b_.*tmp_b_,tmp_h__(:,1+nh)) - tmp_avg_(1+nh).^2);
end;%for nh=0:n_h-1;
tmp_h___(:,:,1+nnstep) = tmp_h__;
tmp_avg__(:,1+nnstep) = tmp_avg_;
tmp_std__(:,1+nnstep) = tmp_std_;
end;%for nnstep=0:n_nstep+1-1;
%%%%%%%%;
flag_plot=0;
if flag_plot;
fname_fig = sprintf('%s/dolphin_drift0_dt_scatter_nvar%.2d_%s_FIGD',dir_jpg_base,nvar_ori,string_dat_ori_name_title);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1); clf; figbig; fig80s();
markersize_use = 32;
n_plot_row = 2;
n_plot_col = n_nstep+1;
for nnstep=0:n_nstep+1-1;
if (nnstep< n_nstep); tmp_dt_all_ = dt_all__{1+nnstep}; tmp_dv_all_ = dv_all__{1+nnstep}; end;
if (nnstep==n_nstep); tmp_dt_all_ = dt_all_all_; tmp_dv_all_ = dv_all_all_; end;
subplot(n_plot_row,n_plot_col,1+nnstep+0*n_plot_col);
plot(sqrt(tmp_dt_all_),tmp_dv_all_,'k.');
xlim(xlim_);ylim(ylim_); xlabel('sqrt(dt)'); ylabel('dv');
if (nnstep< n_nstep); title(sprintf('%d-step',nstep_(1+nnstep))); end;
if (nnstep==n_nstep); title(sprintf('sum-step')); end;
subplot(n_plot_row,n_plot_col,1+nnstep+1*n_plot_col);
hold on;
imagesc(tmp_h___(:,:,1+nnstep),[0.02,0.20]);
tmp_b_ = 1:n_h_x; tmp_avg_ = tmp_avg__(:,1+nnstep); tmp_std_ = tmp_std__(:,1+nnstep);
plot(tmp_b_,tmp_avg_,'w.','MarkerSize',markersize_use);
plot(tmp_b_,tmp_avg_+1.0*tmp_std_,'y.','MarkerSize',markersize_use);
plot(tmp_b_,tmp_avg_-1.0*tmp_std_,'y.','MarkerSize',markersize_use);
hold off;
xlim([1,n_h_x]); ylim([1,n_h_y]);
axisnotick;
xlabel('sqrt(dt)'); ylabel('dv');
end;%for nnstep=0:n_nstep+1-1;
sgtitle(sprintf('%s',string_dat_ori_name_title));
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%for nvar_ori=0:n_var_ori-1;

%%%%%%%%;
% Now plot a few increments. ;
%%%%%%%%;
for nvar_ori=0:n_var_ori-1;%for nvar_ori=efind(strcmp(string_dat_ori_name_,'MCH'));
string_dat_ori_name = string_dat_ori_name_{1+nvar_ori};
string_dat_ori_name_title = sprintf('%s',string_dat_ori_name);
flag_log_transform = flag_log_transform_(1+nvar_ori);
if flag_log_transform;
string_dat_ori_name_title = sprintf('log(%s)',string_dat_ori_name);
end;%if flag_log_transform;
%%%%%%%%;
n_h_x = 12;
n_h_y = 18;
dt_all_ = zeros(n_smp_ori,1); dv_all_ = zeros(n_smp_ori,1); tab=0;
for nu_aid_ori=0:n_u_aid_ori-1;
index_aid_ori_ = index_aid_ori__{1+nu_aid_ori};
tmp_age_ori_ = age_ori_(1+index_aid_ori_);
tmp_res_ori_ = res_nrm__(1+index_aid_ori_,1+nvar_ori);
tmp_index_finite_ = efind( isfinite(tmp_age_ori_) & isfinite(tmp_res_ori_) );
tmp_age_use_ = tmp_age_ori_(1+tmp_index_finite_);
tmp_res_use_ = tmp_res_ori_(1+tmp_index_finite_);
[tmp_age_srt_,tmp_index_srt_] = sort(tmp_age_use_,'ascend'); tmp_index_srt_ = tmp_index_srt_-1;
tmp_res_srt_ = tmp_res_use_(1+tmp_index_srt_);
nstep=1;
dt_ = tmp_age_srt_(1+nstep:end) - tmp_age_srt_(1:end-nstep);
dv_ = tmp_res_srt_(1+nstep:end) - tmp_res_srt_(1:end-nstep);
tab_add = numel(dt_);
dt_all_(1+tab+[0:tab_add-1]) = dt_;
dv_all_(1+tab+[0:tab_add-1]) = dv_;
tab = tab + tab_add;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
dt_all_ = dt_all_(1:tab);
dv_all_ = dv_all_(1:tab);
%xlim_ = [0,sqrt(max(dt_all_))];
xlim_ = [0,sqrt(prctile(dt_all_,95))];
ylim_ = prctile(dv_all_,[ 1,99]);
tmp_h__ = zeros(n_h_y,n_h_x);
tmp_h__ = hist2d_0(sqrt(dt_all_),dv_all_,n_h_x,n_h_y,xlim_,ylim_);
tmp_g__ = tmp_h__;
tmp_e_ = linspace(ylim_(1),ylim_(2),n_h_x);
tmp_b_ = 1:n_h_y;
tmp_avg_ = zeros(n_h_x,1);
tmp_avg_dv_ = zeros(n_h_x,1);
tmp_std_ = zeros(n_h_x,1);
tmp_std_dv_ = zeros(n_h_x,1);
for nh=0:n_h_x-1;
tmp_h__(:,1+nh) = tmp_h__(:,1+nh)/max(1,sum(tmp_h__(:,1+nh)));
tmp_g__(:,1+nh) = tmp_g__(:,1+nh)/max(1,max(tmp_g__(:,1+nh)));
tmp_avg_(1+nh) = dot(tmp_b_,tmp_h__(:,1+nh));
tmp_avg_dv_(1+nh) = (tmp_avg_(1+nh) - 1 - (n_h_y-1)/2) / (n_h_y-1) * diff(ylim_) + mean(ylim_);
tmp_std_(1+nh) = sqrt(dot(tmp_b_.*tmp_b_,tmp_h__(:,1+nh)) - tmp_avg_(1+nh).^2);
tmp_std_dv_(1+nh) = tmp_std_(1+nh) / (n_h_y-1) * diff(ylim_) + mean(ylim_);
end;%for nh=0:n_h-1;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_drift0_dt_scatter_nvar%.2d_%s_FIGE',dir_jpg_base,nvar_ori,string_dat_ori_name_title);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1); clf; figsml; fig80s();
c_80s__ = colormap_80s;
sgtitle(string_dat_ori_name_title);
markersize_use = 8; linewidth_use = 3; markeredgecolor_use = 0.65*[1,1,1];
n_plot_block = 5;
n_plot_row = 1; n_plot_col = 2*n_plot_block+1; ns=0;
subplot(n_plot_row,n_plot_col,1+ns+[0:n_plot_block-1]);ns=ns+n_plot_block;
hold on;
plot(sqrt(dt_all_),dv_all_,'k.');
tmp_bv_ = linspace(min(xlim_),max(xlim_),1+n_h_x); tmp_bv_ = 0.5*(tmp_bv_(1:end-1)+tmp_bv_(2:end));
plot(tmp_bv_,tmp_avg_dv_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','w','LineWidth',linewidth_use,'MarkerEdgeColor',markeredgecolor_use);
plot(tmp_bv_,tmp_avg_dv_+1.0*tmp_std_dv_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','y','LineWidth',linewidth_use,'MarkerEdgeColor',markeredgecolor_use);
plot(tmp_bv_,tmp_avg_dv_-1.0*tmp_std_dv_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','y','LineWidth',linewidth_use,'MarkerEdgeColor',markeredgecolor_use);
hold off;
xlim(xlim_);ylim(ylim_); xlabel('sqrt(dt)'); ylabel('dv');
subplot(n_plot_row,n_plot_col,1+ns+[0:n_plot_block-1]);ns=ns+n_plot_block;
hold on;
%imagesc(tmp_h__,[0.02,0.20]);
%imagesc(tmp_g__,[0.00,exp(-1/2)]); %<-- one standard-deviation away from mean. ;
imagesc(tmp_g__,[0.00,1.00]); %<-- full range. ;
tmp_b_ = 1:n_h_x;
plot(tmp_b_,tmp_avg_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','w','LineWidth',linewidth_use,'MarkerEdgeColor',markeredgecolor_use);
plot(tmp_b_,tmp_avg_+1.0*tmp_std_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','y','LineWidth',linewidth_use,'MarkerEdgeColor',markeredgecolor_use);
plot(tmp_b_,tmp_avg_-1.0*tmp_std_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','y','LineWidth',linewidth_use,'MarkerEdgeColor',markeredgecolor_use);
hold off;
xlim([1,n_h_x]); ylim([1,n_h_y]);
axisnotick;
xlabel('sqrt(dt)'); ylabel('dv');
subplot(n_plot_row,n_plot_col,n_plot_col);
imagesc(transpose(1:64));
set(gca,'Ydir','normal');
set(gca,'YTick',[1,64],'YTickLabel',[0,1]);
set(gca,'TickLength',[0,0]);
set(gca,'XTick',[]);
sgtitle(sprintf('%s',string_dat_ori_name_title));
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%for nvar_ori=0:n_var_ori-1;

flag_check=0;
if flag_check;
%%%%%%%%;
% Note: at this point the residual-matrix is quite noisy: ;
% There is not an obvious age-related signal in the dominant 2 principal-components. ;
% and the total rank is close to 40. ;
%%%%%%%%;
res_fin__ = isfinite(res_nrm__);
tmp_index_fin_ = efind(sum(res_fin__,2));
res_nrm_sub__ = res_nrm__(1+tmp_index_fin_,:);
res_nrm_sub__ = transpose(res_nrm_sub__);
ij_missed_ = find(~isfinite(res_nrm_sub__));
res_nrm_sub_imp__ = res_nrm_sub__; res_nrm_sub_imp__(ij_missed_)=0;
%rank_estimate_res = rank_estimate_sample_1(transpose(res_nrm_sub_imp__),0.05);
rank_estimate_res = 40;
[res_nrm_sub_imp__,res_nrm_sub_fit__] = svd_impute_fit_3(res_nrm_sub__,ij_missed_,rank_estimate_res);
age_sub_ = age_ori_(1+tmp_index_fin_);
[tmp_U__,tmp_S__,tmp_V__] = svds(res_nrm_sub_fit__,2);
%%%%%%%%;
figure(1+nf);nf=nf+1; figbig;
markersize_use = 6;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
age_sub_min = min(age_sub_);
age_sub_max = max(age_sub_);
hold on;
for nl=0:numel(age_sub_)-1;
age_sub = age_sub_(1+nl);
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*(age_sub - age_sub_min)/(age_sub_max - age_sub_min))));
plot(tmp_V__(1+nl,1),tmp_V__(1+nl,2),'ko','MarkerSize',markersize_use,'MarkerFaceColor',c_80s__(1+nc_80s,:));
end;%for nl=0:numel(age_sub_)-1;
hold off;
%%%%%%%%;
end;%if flag_check;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now examine slopes of individual dolphins for the blood-based aging-rate biomarkers ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
bbar_upb_vs_lob_ = ...
[ ...
 0 ...
,0 ...
,0 ...
,0 ...
,1 ...
,1 ...
];
bbar_lob_ = ...
[ ...
  12 ...
,163 ...
, 49 ...
,  6 ...
,1.1 ...
,6.4 ...
];
bbar_upb_ = ...
[ ...
  16 ...
,470 ...
,130 ...
, 23 ...
,1.8 ...
,7.8 ...
];
index_nvar_bbar_ = ...
[ ...
 efind(strcmp(string_dat_ori_name_,'HGB')) ...
,efind(strcmp(string_dat_ori_name_,'AlkPhos')) ...
,efind(strcmp(string_dat_ori_name_,'Platelets')) ...
,efind(strcmp(string_dat_ori_name_,'Lymphs')) ...
,efind(strcmp(string_dat_ori_name_,'Creatinine')) ...
,efind(strcmp(string_dat_ori_name_,'Protein')) ...
];
n_var_bbar = numel(index_nvar_bbar_);
polyfit_dat_ind_pnv___ = NaN*ones(1+n_degree,n_u_aid_ori,n_var_bbar);
for nu_aid_ori=0:n_u_aid_ori-1;
u_aid_ori = u_aid_ori_(1+nu_aid_ori);
index_aid_ori_ = index_aid_ori__{1+nu_aid_ori};
tmp_age_ori_ = age_ori_(1+index_aid_ori_);
for nvar_bbar=0:n_var_bbar-1;
nvar_ori = index_nvar_bbar_(1+nvar_bbar);
tmp_dat_nrm_ = dat_nrm__(1+index_aid_ori_,1+nvar_ori);
tmp_index_use_ = efind( 1 ...
&isfinite(tmp_age_ori_) ...
&(tmp_age_ori_>=che_lob) ...
&(tmp_age_ori_<=che_upb) ...
&isfinite(tmp_dat_nrm_) ...
);
if numel(tmp_index_use_)> (1+n_degree);
tmp_age_use_ = tmp_age_ori_(1+tmp_index_use_);
tmp_dat_use_ = tmp_dat_nrm_(1+tmp_index_use_);
polyfit_dat_ind_ = polyfit_safe_0(tmp_age_use_,tmp_dat_use_,n_degree);
polyfit_dat_ind_pnv___(:,1+nu_aid_ori,1+nvar_bbar) = polyfit_dat_ind_;
end;%if numel(tmp_index_use_)> (1+n_degree);
end;%for nvar_bbar=0:n_var_bbar-1;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
%%%%%%%%;
% restrict attention only to dolphins with at least 20 data-points, ;
% and find principal-components of slope. ;
%%%%%%%%;
polyfit_dat_ind_slope_nv__ = reshape(polyfit_dat_ind_pnv___(1,:,:),[n_u_aid_ori,n_var_bbar]);
tmp_index_use_ = efind( (n_u_aid_ori_>=20) & isfinite(sum(polyfit_dat_ind_slope_nv__,2)) );
polyfit_dat_ind_slope_use_nv__ = polyfit_dat_ind_slope_nv__(1+tmp_index_use_,:);
[U_bbar__,S_bbar__,V_bbar__] = svds(polyfit_dat_ind_slope_use_nv__,6);
U_bbar_all__ = polyfit_dat_ind_slope_nv__*V_bbar__*inv(S_bbar__);
%plot(U_bbar__(:,1),U_bbar__(:,2),'rx',U_bbar_all__(:,1),U_bbar_all__(:,2),'bo');
%%%%%%%%;
% pinpoint dolphin in the lower 10th percentile, ;
% or upper 90th percentile for the blood-based-aging-rate biomarkers. ;
%%%%%%%%;
age_cut = 10;
is_bbar_lob_ab__ = zeros(n_u_aid_ori,n_var_bbar);
is_bbar_upb_ab__ = zeros(n_u_aid_ori,n_var_bbar);
is_bbar_old_ab__ = zeros(n_u_aid_ori,n_var_bbar);
at_bbar_lob_ab__ = NaN*ones(n_u_aid_ori,n_var_bbar);
at_bbar_upb_ab__ = NaN*ones(n_u_aid_ori,n_var_bbar);
at_bbar_old_ab__ = NaN*ones(n_u_aid_ori,n_var_bbar);
age_lob_a_ = NaN*ones(n_u_aid_ori,1);
age_upb_a_ = NaN*ones(n_u_aid_ori,1);
for nu_aid_ori=0:n_u_aid_ori-1;
u_aid_ori = u_aid_ori_(1+nu_aid_ori);
index_aid_ori_ = index_aid_ori__{1+nu_aid_ori};
tmp_age_ori_ = age_ori_(1+index_aid_ori_);
age_lob_a_(1+nu_aid_ori) = min(tmp_age_ori_);
age_upb_a_(1+nu_aid_ori) = max(tmp_age_ori_);
for nvar_bbar=0:n_var_bbar-1;
nvar_ori = index_nvar_bbar_(1+nvar_bbar);
tmp_dat_ori_ = dat_ori__(1+index_aid_ori_,1+nvar_ori);
tmp_index_use_ = efind( isfinite(tmp_age_ori_) & isfinite(tmp_dat_ori_) );
tmp_age_use_ = tmp_age_ori_(1+tmp_index_use_);
tmp_dat_use_ = tmp_dat_ori_(1+tmp_index_use_);
tmp_index_lob_ = efind(tmp_dat_use_<=bbar_lob_(1+nvar_bbar));
tmp_is_lob = 0;
tmp_at_lob = NaN;
if ~isempty(tmp_index_lob_);
if max(tmp_age_use_(1+tmp_index_lob_))>age_cut;
tmp_is_lob = 1;
%tmp_at_lob = max(age_cut,min(tmp_age_use_(1+tmp_index_lob_)));
tmp_at_lob = max(age_lob,min(tmp_age_use_(1+tmp_index_lob_)));
end;%if max(tmp_age_use_(1+tmp_index_lob_))>age_cut;
end;%if ~isempty(tmp_index_lob_);
tmp_index_upb_ = efind(tmp_dat_use_>=bbar_upb_(1+nvar_bbar));
tmp_is_upb = 0;
if ~isempty(tmp_index_upb_);
if max(tmp_age_use_(1+tmp_index_upb_))>age_cut;
tmp_is_upb = 1;
%tmp_at_upb = max(age_cut,min(tmp_age_use_(1+tmp_index_upb_)));
tmp_at_upb = max(age_lob,min(tmp_age_use_(1+tmp_index_upb_)));
end;%if max(tmp_age_use_(1+tmp_index_upb_))>age_cut;
end;%if ~isempty(tmp_index_upb_);
is_bbar_lob_ab__(1+nu_aid_ori,1+nvar_bbar) = tmp_is_lob;
is_bbar_upb_ab__(1+nu_aid_ori,1+nvar_bbar) = tmp_is_upb;
at_bbar_lob_ab__(1+nu_aid_ori,1+nvar_bbar) = tmp_at_lob;
at_bbar_upb_ab__(1+nu_aid_ori,1+nvar_bbar) = tmp_at_upb;
if (bbar_upb_vs_lob_(1+nvar_bbar)==0);
is_bbar_old_ab__(1+nu_aid_ori,1+nvar_bbar) = tmp_is_lob;
at_bbar_old_ab__(1+nu_aid_ori,1+nvar_bbar) = tmp_at_lob;
end;%if (bbar_upb_vs_lob_(1+nvar_bbar)==0);
if (bbar_upb_vs_lob_(1+nvar_bbar)==1);
is_bbar_old_ab__(1+nu_aid_ori,1+nvar_bbar) = tmp_is_upb;
at_bbar_old_ab__(1+nu_aid_ori,1+nvar_bbar) = tmp_at_upb;
end;%if (bbar_upb_vs_lob_(1+nvar_bbar)==1);
end;%for nvar_bbar=0:n_var_bbar-1;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
%%%%%%%%;
% Now associate U_bbar_all__ with individual measurements. ;
%%%%%%%%;
U_bbar_ori__ = NaN*ones(n_u_aid_ori,n_var_bbar);
tmp_a_ = []; na=0;
tmp_U_ = []; na=0;
for nu_aid_ori=0:n_u_aid_ori-1;
u_aid_ori = u_aid_ori_(1+nu_aid_ori);
index_aid_ori_ = index_aid_ori__{1+nu_aid_ori};
tmp_n = n_u_aid_ori_(1+nu_aid_ori);
U_bbar_all = U_bbar_all__(1+nu_aid_ori,1);
if ~isfinite(U_bbar_all);
U_bbar_ori__(1+index_aid_ori_,:) = NaN;
end;%if ~isfinite(U_bbar_all);
if  isfinite(U_bbar_all);
U_bbar_ori__(1+index_aid_ori_,:) = repmat(U_bbar_all__(1+nu_aid_ori,:),[tmp_n,1]);
tmp_a_ = [tmp_a_;u_aid_ori*ones(tmp_n,1)];
tmp_U_ = [tmp_U_;U_bbar_all*ones(tmp_n,1)];
na=na+tmp_n;
end;%if  isfinite(U_bbar_all);
end;%for nu_aid_ori=0:n_u_aid_ori-1;
U_p50 = median(tmp_U_);

%%%%%%%%;
% Now make a scatterplot of the individual dolphins projected on bbar pcs. ;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_bbar_old_FIGG',dir_jpg_base);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(1);clf;figbig;
p_row=3;p_col=2;ns=0;
markersize_use = 16;
linewidth_use = 2;
ulim_ = prctile(U_bbar_all__(1+efind(n_u_aid_ori_>50),1),[  0,100]);
ulim_ = mean(ulim_) + 1.25*0.5*diff(ulim_)*[-1,+1];
for nvar_bbar=0:n_var_bbar-1;
nvar_ori = index_nvar_bbar_(1+nvar_bbar);
subplot(p_row,p_col,1+nvar_bbar);
hold on;
plot(U_p50*[1;1],[age_lob;age_upb],'-','Color',0.65*[1,1,1],'LineWidth',3);
plot(ulim_,age_cut*[1,1],':','Color',0.65*[1,1,1],'LineWidth',1);
for nu_aid_ori=0:n_u_aid_ori-1;
is_bbar_old = is_bbar_old_ab__(1+nu_aid_ori,1+nvar_bbar);
U_bbar_all = U_bbar_all__(1+nu_aid_ori,1);
tmp_age_lob = age_lob_a_(1+nu_aid_ori);
tmp_age_upb = age_upb_a_(1+nu_aid_ori);
if ~is_bbar_old;
plot(U_bbar_all*[1;1],[tmp_age_lob;tmp_age_upb],'g-','LineWidth',linewidth_use);
end;%if ~is_bbar_old;
if  is_bbar_old;
tmp_age_mid = at_bbar_old_ab__(1+nu_aid_ori,1+nvar_bbar);
plot(U_bbar_all*[1;1],[tmp_age_lob;tmp_age_mid],'g-','LineWidth',linewidth_use);
plot(U_bbar_all*[1;1],[tmp_age_mid;tmp_age_upb],'r-','LineWidth',linewidth_use);
plot(U_bbar_all*[1],[tmp_age_mid],'k+','LineWidth',linewidth_use);
end;%if  is_bbar_old;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
hold off;
xlim(ulim_); xlabel('principal-component');
%ylim([age_cut,age_upb]); ylabel('age');
ylim([age_lob,age_upb]); ylabel('age');
title(string_dat_ori_name_{1+nvar_ori});
end;%for nvar_bbar=0:n_var_bbar-1;
%%%%%%%%;
%sgtitle(sprintf('%s',fname_fig),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Check histograms of U_bbar_ori__ by age and sex. ;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_bbar_hist_FIGI',dir_jpg_base);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(1);clf;figbig;
ulim_ = prctile(U_bbar_all__(1+efind(n_u_aid_ori_>50),1),[  0,100]);
p_row=3;p_col=2;ns=0;
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
tmp_index_ = efind( isfinite(U_bbar_ori__(:,1)) & isM_ori_ );
tmp_U_bbar_ori_ = U_bbar_ori__(1+tmp_index_,1);
hist(tmp_U_bbar_ori_,linspace(min(ulim_),max(ulim_),128));
xlim(ulim_);
tmp_f_lo = numel(efind(tmp_U_bbar_ori_<=U_p50))/numel(tmp_U_bbar_ori_);
tmp_f_hi = numel(efind(tmp_U_bbar_ori_> U_p50))/numel(tmp_U_bbar_ori_);
title(sprintf('U isM: %0.2f %0.2f',tmp_f_lo,tmp_f_hi));
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
tmp_index_ = efind( isfinite(U_bbar_ori__(:,1)) & isF_ori_ );
tmp_U_bbar_ori_ = U_bbar_ori__(1+tmp_index_,1);
hist(tmp_U_bbar_ori_,linspace(min(ulim_),max(ulim_),128));
xlim(ulim_);
tmp_f_lo = numel(efind(tmp_U_bbar_ori_<=U_p50))/numel(tmp_U_bbar_ori_);
tmp_f_hi = numel(efind(tmp_U_bbar_ori_> U_p50))/numel(tmp_U_bbar_ori_);
title(sprintf('U isF: %0.2f %0.2f',tmp_f_lo,tmp_f_hi));
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
tmp_index_ = efind( isfinite(U_bbar_ori__(:,1)) & (age_ori_>= 0) & (age_ori_< 10) );
tmp_U_bbar_ori_ = U_bbar_ori__(1+tmp_index_,1);
hist(tmp_U_bbar_ori_,linspace(min(ulim_),max(ulim_),128));
xlim(ulim_);
tmp_f_lo = numel(efind(tmp_U_bbar_ori_<=U_p50))/numel(tmp_U_bbar_ori_);
tmp_f_hi = numel(efind(tmp_U_bbar_ori_> U_p50))/numel(tmp_U_bbar_ori_);
title(sprintf('U age  0-10: %0.2f %0.2f',tmp_f_lo,tmp_f_hi));
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
tmp_index_ = efind( isfinite(U_bbar_ori__(:,1)) & (age_ori_>=10) & (age_ori_< 30) );
tmp_U_bbar_ori_ = U_bbar_ori__(1+tmp_index_,1);
hist(tmp_U_bbar_ori_,linspace(min(ulim_),max(ulim_),128));
xlim(ulim_);
tmp_f_lo = numel(efind(tmp_U_bbar_ori_<=U_p50))/numel(tmp_U_bbar_ori_);
tmp_f_hi = numel(efind(tmp_U_bbar_ori_> U_p50))/numel(tmp_U_bbar_ori_);
title(sprintf('U age 10-30: %0.2f %0.2f',tmp_f_lo,tmp_f_hi));
%%%%%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
tmp_index_ = efind( isfinite(U_bbar_ori__(:,1)) & (age_ori_>=30) & (age_ori_< 55) );
tmp_U_bbar_ori_ = U_bbar_ori__(1+tmp_index_,1);
hist(tmp_U_bbar_ori_,linspace(min(ulim_),max(ulim_),128));
xlim(ulim_);
tmp_f_lo = numel(efind(tmp_U_bbar_ori_<=U_p50))/numel(tmp_U_bbar_ori_);
tmp_f_hi = numel(efind(tmp_U_bbar_ori_> U_p50))/numel(tmp_U_bbar_ori_);
title(sprintf('U age 30-55: %0.2f %0.2f',tmp_f_lo,tmp_f_hi));
%%%%%%%%;
sgtitle(sprintf('%s',fname_fig),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

isU_bbar_hi_ori_ = zeros(n_smp_ori,1);
tmp_index_ = efind( isfinite(U_bbar_ori__(:,1)) & (U_bbar_ori__(:,1)> U_p50) );
isU_bbar_hi_ori_(1+tmp_index_) = 1;
isU_bbar_lo_ori_ = zeros(n_smp_ori,1);
tmp_index_ = efind( isfinite(U_bbar_ori__(:,1)) & (U_bbar_ori__(:,1)< U_p50) );
isU_bbar_lo_ori_(1+tmp_index_) = 1;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now represent a 2-variable case-study using res_ind. ;
% Here we plot A__. ;
%%%%%%%%;
nvar_pair_ = ...
{ ...
  { 'SED60' , 'CPK' } ...
 ,{ 'SED60' , 'Albumin' } ...
 ,{ 'CPK' , 'Creatinine' } ...
 ,{ 'GFR' , 'Creatinine' } ...
 ,{ 'Iron' , 'MCH' } ...
 ,{ 'Bilirubin' , 'MCHC' } ...
 ,{ 'Bilirubin' , 'HCT' } ...
 ,{ 'CPK' , 'Creatinine' } ...
 ,{ 'HGB' , 'MCH' } ...
 ,{ 'RBCDist' , 'MCH' } ...
 ,{ 'MPV' , 'UricAcid' } ...
 ,{ 'RBC' , 'Platelets' } ...
};
n_pair = numel(nvar_pair_);
%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npair=0:n_pair-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
str_pair_0 = nvar_pair_{1+npair}{1+0};
str_pair_1 = nvar_pair_{1+npair}{1+1};
nvar_ori_0 = efind(strcmp(string_dat_ori_name_,str_pair_0));
nvar_ori_1 = efind(strcmp(string_dat_ori_name_,str_pair_1));
res_ind_nrm_sub__ = res_ind_nrm__(:,1+[nvar_ori_0,nvar_ori_1]);
n_var_sub = 2;
%%%%%%%%;
n_shuffle = 0;
a_ind_sub_prm__ = zeros(n_var_sub,1+n_shuffle);
A_ind_sub_prm___ = zeros(n_var_sub,n_var_sub,1+n_shuffle);
BB_inv_ind_sub_prm___ = zeros(n_var_sub,n_var_sub,1+n_shuffle);
CC_inv_ind_sub_prm___ = zeros(n_var_sub,n_var_sub,1+n_shuffle);
for nshuffle=0:n_shuffle;
if (mod(nshuffle,16)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
res_ind_nrm_sub_prm__ = res_ind_nrm_sub__;
if (nshuffle>0);
[res_ind_nrm_sub_prm__] = dolphin_permute_0(aid_ori_,age_ori_,res_ind_nrm_sub__,nshuffle);
end;%if (nshuffle>0);
parameter = struct('type','parameter');
[ ...
 parameter ...
,a_ind_sub_prm_ ...
,A_ind_sub_prm__ ...
,BB_inv_ind_sub_prm__ ...
,CC_inv_ind_sub_prm__ ...
,L_ind_sub_prm ...
,niteration_ind_sub_prm ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ori_ ...
,age_ori_ ...
,res_ind_nrm_sub_prm__ ...
);
a_ind_sub_prm__(:,1+nshuffle) = a_ind_sub_prm_;
A_ind_sub_prm___(:,:,1+nshuffle) = A_ind_sub_prm__;
BB_inv_ind_sub_prm___(:,:,1+nshuffle) = BB_inv_ind_sub_prm__;
CC_inv_ind_sub_prm___(:,:,1+nshuffle) = CC_inv_ind_sub_prm__;
end;%for nshuffle=0:n_shuffle;
%%%%%%%%;
a_ind_sub_avg_ = mean(a_ind_sub_prm__(:,2:end),2);
a_ind_sub_std_ = std(a_ind_sub_prm__(:,2:end),1,2);
A_ind_sub_avg__ = mean(A_ind_sub_prm___(:,:,2:end),3);
A_ind_sub_std__ = std(A_ind_sub_prm___(:,:,2:end),1,3);
BB_inv_ind_sub_avg__ = mean(BB_inv_ind_sub_prm___(:,:,2:end),3);
BB_inv_ind_sub_std__ = std(BB_inv_ind_sub_prm___(:,:,2:end),1,3);
CC_inv_ind_sub_avg__ = mean(CC_inv_ind_sub_prm___(:,:,2:end),3);
CC_inv_ind_sub_std__ = std(CC_inv_ind_sub_prm___(:,:,2:end),1,3);
a_ind_sub_z__ = (a_ind_sub_prm__ - repmat(a_ind_sub_avg_,[1,1+n_shuffle]))./repmat(a_ind_sub_std_,[1,1+n_shuffle]);
a_ind_sub_nlp__ = -z_to_p_twosided_0(a_ind_sub_z__);
a_ind_sub_snlp__ = sign(a_ind_sub_prm__ - repmat(a_ind_sub_avg_,[1,1+n_shuffle])).*a_ind_sub_nlp__;
A_ind_sub_z___ = (A_ind_sub_prm___ - repmat(A_ind_sub_avg__,[1,1,1+n_shuffle]))./repmat(A_ind_sub_std__,[1,1,1+n_shuffle]);
A_ind_sub_nlp___ = -z_to_p_twosided_0(A_ind_sub_z___);
A_ind_sub_snlp___ = sign(A_ind_sub_prm___ - repmat(A_ind_sub_avg__,[1,1,1+n_shuffle])).*A_ind_sub_nlp___;
BB_inv_ind_sub_z___ = (BB_inv_ind_sub_prm___ - repmat(BB_inv_ind_sub_avg__,[1,1,1+n_shuffle]))./repmat(BB_inv_ind_sub_std__,[1,1,1+n_shuffle]);
BB_inv_ind_sub_nlp___ = -z_to_p_twosided_0(BB_inv_ind_sub_z___);
BB_inv_ind_sub_snlp___ = sign(BB_inv_ind_sub_prm___ - repmat(BB_inv_ind_sub_avg__,[1,1,1+n_shuffle])).*BB_inv_ind_sub_nlp___;
CC_inv_ind_sub_z___ = (CC_inv_ind_sub_prm___ - repmat(CC_inv_ind_sub_avg__,[1,1,1+n_shuffle]))./repmat(CC_inv_ind_sub_std__,[1,1,1+n_shuffle]);
CC_inv_ind_sub_nlp___ = -z_to_p_twosided_0(CC_inv_ind_sub_z___);
CC_inv_ind_sub_snlp___ = sign(CC_inv_ind_sub_prm___ - repmat(CC_inv_ind_sub_avg__,[1,1,1+n_shuffle])).*CC_inv_ind_sub_nlp___;
%%%%%%%%;
%%%%;
%%%%%%%%;
a_ind_sub_ = a_ind_sub_prm__(:,1);
A_ind_sub__ = A_ind_sub_prm___(:,:,1);
BB_inv_ind_sub__ = BB_inv_ind_sub_prm___(:,:,1);
CC_inv_ind_sub__ = CC_inv_ind_sub_prm___(:,:,1);
index_mss_ = efind(~isfinite(res_ind_nrm_sub__));
res_ind_nrm_sub_imp__ = dolphin_impute_aA_0(aid_ori_,age_ori_,res_ind_nrm_sub__,a_ind_sub_,A_ind_sub__,index_mss_);
%%%%%%%%;
[ ...
 parameter ...
,res_ind_nrm_sub_fil__ ...
,P_pos_nnv___ ...
] = ...
dolphin_estimate_kalman_filter_0( ...
 parameter ...
,aid_ori_ ...
,age_ori_ ...
,res_ind_nrm_sub_imp__ ...
,a_ind_sub_ ...
,A_ind_sub__ ...
,1*BB_inv_ind_sub__ ...
,1*CC_inv_ind_sub__ ...
);
[ ...
 parameter ...
,res_ind_nrm_sub_fil__ ...
,P_pos_nnv___ ...
] = ...
dolphin_estimate_kalman_filter_0( ...
 parameter ...
,aid_ori_ ...
,age_ori_ ...
,res_ind_nrm_sub_imp__ ...
,a_ind_sub_ ...
,A_ind_sub__ ...
,1*BB_inv_ind_sub__ ...
,1*CC_inv_ind_sub__ ...
,P_pos_nnv___ ...
);
%%%%%%%%;
%%%%;
%%%%%%%%;
parameter_gen = struct('type','parameter');
parameter_gen.tolerance_master = 1e-6;
parameter_gen.dt_avg = 1/64;
parameter_gen.T_max = age_upb;
parameter_gen.rseed = 0;
[ ...
 parameter_gen ...
,age_ind_gen_ ...
,X_ind_gen__ ...
,Y_ind_gen__ ...
] = ...
SDE_generate_data_0( ...
 parameter_gen ...
,a_ind_sub_ ...
,A_ind_sub__ ...
,[] ...
,BB_inv_ind_sub__ ...
,[] ...
,CC_inv_ind_sub__ ...
,zeros(2,1) ...
);
%%%%%%%%;
tmp_res_ind_gen__ = transpose(X_ind_gen__);
glim_ = prctile(X_ind_gen__,[ 1,99],'all');
glim_ = 0*mean(glim_) + 1.25*0.5*diff(glim_)*[-1,+1];
%%%%%%%%;
n_index_use_ = zeros(n_u_aid_ori,1);
for nu_aid_ori=0:n_u_aid_ori-1;
index_aid_ori_ = intersect( index_aid_ori__{1+nu_aid_ori} , index_smp_TT_ );
index_use_ = efind( isfinite(age_ori_(1+index_aid_ori_)) & isfinite(sum(res_ind_nrm_sub_fil__(1+index_aid_ori_,:),2)) );
n_index_use = numel(index_use_);
n_index_use_(1+nu_aid_ori) = n_index_use;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
[~,tmp_index_srt_] = sort(n_index_use_,'descend'); tmp_index_srt_ = tmp_index_srt_-1;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_res_ind_nrm_sub_fil_nvar%.2d_%s_nvar%.2d_%s_FIGF',dir_jpg_base,nvar_ori_0,string_dat_ori_name_{1+nvar_ori_0},nvar_ori_1,string_dat_ori_name_{1+nvar_ori_1});
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;
p_row = 2; p_col = 6; p_pers = 1; ns=0;
p_big_ = [ 0 , 1 , p_col , p_col+1 ];
p_sml_ = setdiff([0:p_row*p_col-1],p_big_);
linewidth_use = 3;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
%%%%%%%%;
parameter_gen = struct('type','parameter');
parameter_gen.x0_lim_ = glim_;
parameter_gen.x1_lim_ = glim_;
subplot(p_row,p_col,1 + p_big_);cla;
vector_field_0(parameter_gen,A_ind_sub__); axisnotick; axis image;
hold on;
s=surfline_0(tmp_res_ind_gen__(:,1+0),tmp_res_ind_gen__(:,1+1),age_ind_gen_);
set(s,'LineWidth',linewidth_use);
hold off;
set(gca,'Clim',[age_lob,age_upb]);
xlabel(string_dat_ori_name_{1+nvar_ori_0});
ylabel(string_dat_ori_name_{1+nvar_ori_1});
title('SDE sample');
%%%%%%%%;
for nl=0:n_u_aid_ori-1;
nu_aid_ori = tmp_index_srt_(1+nl);
if ns<numel(p_sml_);
subplot(p_row,p_col,1+p_sml_(1+ns));
%%%%%%%%;
index_aid_ori_ = intersect( index_aid_ori__{1+nu_aid_ori} , index_smp_TT_ );
index_use_ = efind( isfinite(age_ori_(1+index_aid_ori_)) & isfinite(sum(res_ind_nrm_sub_fil__(1+index_aid_ori_,:),2)) );
[tmp_age_srt_,tmp_ij] = sort(age_ori_(1+index_aid_ori_(1+index_use_)),'ascend');
tmp_res_ind_srt__ = res_ind_nrm_sub_fil__(1+index_aid_ori_(1+index_use_(tmp_ij)),:);
tmp_age_srt_ = tmp_age_srt_(2:end); tmp_res_ind_srt__ = tmp_res_ind_srt__(2:end,:);
%%%%;
xlim_ = prctile(tmp_res_ind_srt__(:,1+0),[5,95],'all');
xlim_ = mean(xlim_) + 1.25*0.5*diff(xlim_)*[-1,+1];
ylim_ = prctile(tmp_res_ind_srt__(:,1+1),[5,95],'all');
ylim_ = mean(ylim_) + 1.25*0.5*diff(ylim_)*[-1,+1];
hold on;
s=surfline_0(tmp_res_ind_srt__(:,1+0),tmp_res_ind_srt__(:,1+1),tmp_age_srt_);
set(s,'LineWidth',linewidth_use);
hold off;
set(gca,'Clim',[age_lob,age_upb]);
%%%%%%%%;
xlim(xlim_);xlabel(string_dat_ori_name_{1+nvar_ori_0});
ylim(ylim_);ylabel(string_dat_ori_name_{1+nvar_ori_1});
axis square; axis image;
title(sprintf('naid %d',u_aid_ori_(1+nu_aid_ori)));
grid on;
end;%if ns<numel(p_sml_)-1;
if (mod(nl,p_pers)==0); ns=ns+1; end;
end;%for nl=0:n_u_aid_ori-1;
fig80s;
sgtitle(sprintf('%s vs %s',string_dat_ori_name_{1+nvar_ori_0},string_dat_ori_name_{1+nvar_ori_1}));
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for npair=0:n_pair-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now represent a 2-variable case-study using res_ind. ;
% Here we plot B__. ;
%%%%%%%%;
nvar_pair_ = ...
{ ...
  { 'MCH' , 'Protein' } ...
 ,{ 'MCH' , 'Protein' } ...
};
n_pair = numel(nvar_pair_);
%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npair=0:n_pair-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
str_pair_0 = nvar_pair_{1+npair}{1+0};
str_pair_1 = nvar_pair_{1+npair}{1+1};
nvar_ori_0 = efind(strcmp(string_dat_ori_name_,str_pair_0));
nvar_ori_1 = efind(strcmp(string_dat_ori_name_,str_pair_1));
res_ind_nrm_sub__ = res_ind_nrm__(:,1+[nvar_ori_0,nvar_ori_1]);
n_var_sub = 2;
%%%%%%%%;
n_shuffle = 0;
a_ind_sub_prm__ = zeros(n_var_sub,1+n_shuffle);
A_ind_sub_prm___ = zeros(n_var_sub,n_var_sub,1+n_shuffle);
BB_inv_ind_sub_prm___ = zeros(n_var_sub,n_var_sub,1+n_shuffle);
CC_inv_ind_sub_prm___ = zeros(n_var_sub,n_var_sub,1+n_shuffle);
for nshuffle=0:n_shuffle;
if (mod(nshuffle,16)==0); disp(sprintf(' %% nshuffle %d/%d',nshuffle,n_shuffle)); end;
res_ind_nrm_sub_prm__ = res_ind_nrm_sub__;
if (nshuffle>0);
[res_ind_nrm_sub_prm__] = dolphin_permute_0(aid_ori_,age_ori_,res_ind_nrm_sub__,nshuffle);
end;%if (nshuffle>0);
parameter = struct('type','parameter');
[ ...
 parameter ...
,a_ind_sub_prm_ ...
,A_ind_sub_prm__ ...
,BB_inv_ind_sub_prm__ ...
,CC_inv_ind_sub_prm__ ...
,L_ind_sub_prm ...
,niteration_ind_sub_prm ...
] = ...
dolphin_estimate_aABC_4( ...
 parameter ...
,aid_ori_ ...
,age_ori_ ...
,res_ind_nrm_sub_prm__ ...
);
a_ind_sub_prm__(:,1+nshuffle) = a_ind_sub_prm_;
A_ind_sub_prm___(:,:,1+nshuffle) = A_ind_sub_prm__;
BB_inv_ind_sub_prm___(:,:,1+nshuffle) = BB_inv_ind_sub_prm__;
CC_inv_ind_sub_prm___(:,:,1+nshuffle) = CC_inv_ind_sub_prm__;
end;%for nshuffle=0:n_shuffle;
%%%%%%%%;
a_ind_sub_avg_ = mean(a_ind_sub_prm__(:,2:end),2);
a_ind_sub_std_ = std(a_ind_sub_prm__(:,2:end),1,2);
A_ind_sub_avg__ = mean(A_ind_sub_prm___(:,:,2:end),3);
A_ind_sub_std__ = std(A_ind_sub_prm___(:,:,2:end),1,3);
BB_inv_ind_sub_avg__ = mean(BB_inv_ind_sub_prm___(:,:,2:end),3);
BB_inv_ind_sub_std__ = std(BB_inv_ind_sub_prm___(:,:,2:end),1,3);
CC_inv_ind_sub_avg__ = mean(CC_inv_ind_sub_prm___(:,:,2:end),3);
CC_inv_ind_sub_std__ = std(CC_inv_ind_sub_prm___(:,:,2:end),1,3);
a_ind_sub_z__ = (a_ind_sub_prm__ - repmat(a_ind_sub_avg_,[1,1+n_shuffle]))./repmat(a_ind_sub_std_,[1,1+n_shuffle]);
a_ind_sub_nlp__ = -z_to_p_twosided_0(a_ind_sub_z__);
a_ind_sub_snlp__ = sign(a_ind_sub_prm__ - repmat(a_ind_sub_avg_,[1,1+n_shuffle])).*a_ind_sub_nlp__;
A_ind_sub_z___ = (A_ind_sub_prm___ - repmat(A_ind_sub_avg__,[1,1,1+n_shuffle]))./repmat(A_ind_sub_std__,[1,1,1+n_shuffle]);
A_ind_sub_nlp___ = -z_to_p_twosided_0(A_ind_sub_z___);
A_ind_sub_snlp___ = sign(A_ind_sub_prm___ - repmat(A_ind_sub_avg__,[1,1,1+n_shuffle])).*A_ind_sub_nlp___;
BB_inv_ind_sub_z___ = (BB_inv_ind_sub_prm___ - repmat(BB_inv_ind_sub_avg__,[1,1,1+n_shuffle]))./repmat(BB_inv_ind_sub_std__,[1,1,1+n_shuffle]);
BB_inv_ind_sub_nlp___ = -z_to_p_twosided_0(BB_inv_ind_sub_z___);
BB_inv_ind_sub_snlp___ = sign(BB_inv_ind_sub_prm___ - repmat(BB_inv_ind_sub_avg__,[1,1,1+n_shuffle])).*BB_inv_ind_sub_nlp___;
CC_inv_ind_sub_z___ = (CC_inv_ind_sub_prm___ - repmat(CC_inv_ind_sub_avg__,[1,1,1+n_shuffle]))./repmat(CC_inv_ind_sub_std__,[1,1,1+n_shuffle]);
CC_inv_ind_sub_nlp___ = -z_to_p_twosided_0(CC_inv_ind_sub_z___);
CC_inv_ind_sub_snlp___ = sign(CC_inv_ind_sub_prm___ - repmat(CC_inv_ind_sub_avg__,[1,1,1+n_shuffle])).*CC_inv_ind_sub_nlp___;
%%%%%%%%;
%%%%;
%%%%%%%%;
a_ind_sub_ = a_ind_sub_prm__(:,1);
A_ind_sub__ = A_ind_sub_prm___(:,:,1);
BB_inv_ind_sub__ = BB_inv_ind_sub_prm___(:,:,1);
CC_inv_ind_sub__ = CC_inv_ind_sub_prm___(:,:,1);
index_mss_ = efind(~isfinite(res_ind_nrm_sub__));
res_ind_nrm_sub_imp__ = dolphin_impute_aA_0(aid_ori_,age_ori_,res_ind_nrm_sub__,a_ind_sub_,A_ind_sub__,index_mss_);
%%%%%%%%;
[ ...
 parameter ...
,res_ind_nrm_sub_fil__ ...
,P_pos_nnv___ ...
] = ...
dolphin_estimate_kalman_filter_0( ...
 parameter ...
,aid_ori_ ...
,age_ori_ ...
,res_ind_nrm_sub_imp__ ...
,a_ind_sub_ ...
,A_ind_sub__ ...
,1*BB_inv_ind_sub__ ...
,1*CC_inv_ind_sub__ ...
);
[ ...
 parameter ...
,res_ind_nrm_sub_fil__ ...
,P_pos_nnv___ ...
] = ...
dolphin_estimate_kalman_filter_0( ...
 parameter ...
,aid_ori_ ...
,age_ori_ ...
,res_ind_nrm_sub_imp__ ...
,a_ind_sub_ ...
,A_ind_sub__ ...
,1*BB_inv_ind_sub__ ...
,1*CC_inv_ind_sub__ ...
,P_pos_nnv___ ...
);
%%%%%%%%;
%%%%;
%%%%%%%%;
parameter_gen = struct('type','parameter');
parameter_gen.tolerance_master = 1e-6;
parameter_gen.dt_avg = 1/64;
parameter_gen.T_max = age_upb;
parameter_gen.rseed = 0;
[ ...
 parameter_gen ...
,age_ind_gen_ ...
,X_ind_gen__ ...
,Y_ind_gen__ ...
] = ...
SDE_generate_data_0( ...
 parameter_gen ...
,a_ind_sub_ ...
,A_ind_sub__ ...
,[] ...
,BB_inv_ind_sub__ ...
,[] ...
,CC_inv_ind_sub__ ...
,zeros(2,1) ...
);
%%%%%%%%;
tmp_res_ind_gen__ = transpose(X_ind_gen__);
glim_ = prctile(X_ind_gen__,[ 1,99],'all');
glim_ = 0*mean(glim_) + 1.25*0.5*diff(glim_)*[-1,+1];
%%%%%%%%;
n_index_use_ = zeros(n_u_aid_ori,1);
for nu_aid_ori=0:n_u_aid_ori-1;
index_aid_ori_ = intersect( index_aid_ori__{1+nu_aid_ori} , index_smp_TT_ );
index_use_ = efind( isfinite(age_ori_(1+index_aid_ori_)) & isfinite(sum(res_ind_nrm_sub_fil__(1+index_aid_ori_,:),2)) );
n_index_use = numel(index_use_);
n_index_use_(1+nu_aid_ori) = n_index_use;
end;%for nu_aid_ori=0:n_u_aid_ori-1;
[~,tmp_index_srt_] = sort(n_index_use_,'descend'); tmp_index_srt_ = tmp_index_srt_-1;
%%%%%%%%;
fname_fig = sprintf('%s/dolphin_res_ind_nrm_sub_fil_nvar%.2d_%s_nvar%.2d_%s_FIGF',dir_jpg_base,nvar_ori_0,string_dat_ori_name_{1+nvar_ori_0},nvar_ori_1,string_dat_ori_name_{1+nvar_ori_1});
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;
p_row = 2; p_col = 6; p_pers = 1; ns=0;
p_big_ = [ 0 , 1 , p_col , p_col+1 ];
p_sml_ = setdiff([0:p_row*p_col-1],p_big_);
linewidth_use = 3;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
%%%%%%%%;
parameter_gen = struct('type','parameter');
parameter_gen.x0_lim_ = glim_;
parameter_gen.x1_lim_ = glim_;
subplot(p_row,p_col,1 + p_big_);cla;
vector_field_0(parameter_gen,A_ind_sub__); axisnotick; axis image;
hold on;
s=surfline_0(tmp_res_ind_gen__(:,1+0),tmp_res_ind_gen__(:,1+1),age_ind_gen_);
set(s,'LineWidth',linewidth_use);
hold off;
set(gca,'Clim',[age_lob,age_upb]);
xlabel(string_dat_ori_name_{1+nvar_ori_0});
ylabel(string_dat_ori_name_{1+nvar_ori_1});
title('SDE sample');
%%%%%%%%;
for nl=0:n_u_aid_ori-1;
nu_aid_ori = tmp_index_srt_(1+nl);
if ns<numel(p_sml_);
subplot(p_row,p_col,1+p_sml_(1+ns));
%%%%%%%%;
index_aid_ori_ = intersect( index_aid_ori__{1+nu_aid_ori} , index_smp_TT_ );
index_use_ = efind( isfinite(age_ori_(1+index_aid_ori_)) & isfinite(sum(res_ind_nrm_sub_fil__(1+index_aid_ori_,:),2)) );
[tmp_age_srt_,tmp_ij] = sort(age_ori_(1+index_aid_ori_(1+index_use_)),'ascend');
tmp_res_ind_srt__ = res_ind_nrm_sub_fil__(1+index_aid_ori_(1+index_use_(tmp_ij)),:);
tmp_age_srt_ = tmp_age_srt_(2:end); tmp_res_ind_srt__ = tmp_res_ind_srt__(2:end,:);
%%%%;
xlim_ = prctile(tmp_res_ind_srt__(:,1+0),[5,95],'all');
xlim_ = mean(xlim_) + 1.25*0.5*diff(xlim_)*[-1,+1];
ylim_ = prctile(tmp_res_ind_srt__(:,1+1),[5,95],'all');
ylim_ = mean(ylim_) + 1.25*0.5*diff(ylim_)*[-1,+1];
hold on;
s=surfline_0(tmp_res_ind_srt__(:,1+0),tmp_res_ind_srt__(:,1+1),tmp_age_srt_);
set(s,'LineWidth',linewidth_use);
hold off;
set(gca,'Clim',[age_lob,age_upb]);
%%%%%%%%;
xlim(xlim_);xlabel(string_dat_ori_name_{1+nvar_ori_0});
ylim(ylim_);ylabel(string_dat_ori_name_{1+nvar_ori_1});
axis square; axis image;
title(sprintf('naid %d',u_aid_ori_(1+nu_aid_ori)));
grid on;
end;%if ns<numel(p_sml_)-1;
if (mod(nl,p_pers)==0); ns=ns+1; end;
end;%for nl=0:n_u_aid_ori-1;
fig80s;
sgtitle(sprintf('%s vs %s',string_dat_ori_name_{1+nvar_ori_0},string_dat_ori_name_{1+nvar_ori_1}));
disp(sprintf(' %% writing %s',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for npair=0:n_pair-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;















