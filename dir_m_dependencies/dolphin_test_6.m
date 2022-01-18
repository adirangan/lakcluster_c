function [tru_,est_,null_] = dolphin_test_6(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
% Generates data with drift and tests_dolphin_estimate_aABC_2. ;
% Assumes exponential distribution of dt with mean dt_avg. ;
% Simulates SDE recovery, alongside irg granger_test. ;
% Verdict seems clear: flag_crude==0 is better than flag_crude==1. ;
%{
%try: ;

dXX = 0;
n_var = 2; %<-- pair. ;
dt_avg = 0.257; %<-- matches d00 data. ;
relative_variation = 0.8973; %<-- matches d00 data. ;
for n_T_0in = [128,256,512,1024]; %<-- matches d01 data. ;%n_T_0in = 5664; %<-- matches d00 data. ;
snr_ = 2.^[-3:+0.5:+3]; n_snr = numel(snr_);
rseed_ = 0:16; n_rseed = numel(rseed_);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
for nsnr=0:n_snr-1;
snr = snr_(1+nsnr);
[tru_,est_,null_] = dolphin_test_6(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
end;%for nsnr=0:n_snr-1;
end;%for nrseed=0:n_rseed-1;
end;%for n_T_0in = [128,256,512,1024];

%}

date_diff_threshold = 0.5;
flag_verbose = 0;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); string_root = 'data'; end;
if (strcmp(platform,'eval1')); string_root = 'home'; end;
if (strcmp(platform,'OptiPlex')); string_root = 'home'; end;
dir_trunk = sprintf('/%s/rangan/dir_bcc/dir_dolphin',string_root);
%setup_OptiPlex; dir_trunk = '/home/rangan/dir_bcc/dir_dolphin';
%setup_access1; dir_trunk = '/data/rangan/dir_bcc/dir_dolphin';

dir_R = sprintf('%s/dir_R',dir_trunk);
dir_irg_csv = sprintf('%s/dir_irg_csv',dir_R);
if (~exist(dir_irg_csv,'dir')); disp(sprintf(' %% mkdir %s',dir_irg_csv)); mkdir(dir_irg_csv); end;
dir_irg_mat = sprintf('%s/dir_irg_mat',dir_R);
if (~exist(dir_irg_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_irg_mat)); mkdir(dir_irg_mat); end;
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

flag_replot=1;
flag_a_is_0=0;

string_infix = dolphin_test_6_infix(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
fname_dat_mat = sprintf('%s/%s_dat__.mat',dir_irg_mat,string_infix);
fname_age_csv = sprintf('%s/%s_age_.csv',dir_irg_csv,string_infix);
fname_dat_csv = sprintf('%s/%s_dat__.csv',dir_irg_csv,string_infix);
fname_c0_out_mat = sprintf('%s/%s_c0_out.mat',dir_irg_mat,string_infix);
fname_c1_out_mat = sprintf('%s/%s_c1_out.mat',dir_irg_mat,string_infix);
fname_irg_R = sprintf('%s/%s_irg.R',dir_R,string_infix);
fname_irg_out = sprintf('%s/%s_irg_out__.csv',dir_irg_csv,string_infix);

flag_skip = 0;
if ( exist(fname_irg_out,'file')); flag_skip = 1; end;
if (~exist(fname_irg_out,'file'));
flag_skip=0;
fname_irg_tmp = sprintf('%s/%s_irg_out__.tmp',dir_irg_csv,string_infix);
if ( exist(fname_irg_tmp,'file'));
tmp_date_diff = datenum(clock) - datenum(dir(fname_irg_tmp).date);
if (tmp_date_diff< date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = recent, skipping',fname_irg_tmp,tmp_date_diff));
flag_skip=1;
end;%if (tmp_date_diff< date_diff_threshold);
if (tmp_date_diff>=date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = stale, deleting',fname_irg_tmp,tmp_date_diff));
delete(fname_irg_tmp);
flag_skip=0;
end;%if (tmp_date_diff>=date_diff_threshold);
end;%if ( exist(fname_irg_tmp,'file'));
end;%if (~exist(fname_irg_out,'file'));

if (~flag_skip);
disp(sprintf(' %% %s not found, creating',fname_irg_out));
save(fname_irg_tmp,'fname_irg_tmp');

%%%%%%%%;
% Generate data. ;
%%%%%%%%;
rng(rseed);
%%%%%%%%;
if ( exist(fname_dat_mat,'file')); disp(sprintf(' %% %s found, not creating',fname_dat_mat)); end;
if (~exist(fname_dat_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_dat_mat));
%%%%%%%%;
a_tru_ = 1.0*randn(n_var,1);
A_tru__ = randn(n_var,n_var);
[Psi__,Lambda__] = eig(A_tru__); Psi_inv__ = inv(Psi__);
Lambda_ = diag(Lambda__);
for nvar=0:n_var-1; Lambda_(1+nvar) = min(-0.1,real(Lambda_(1+nvar))) + i*imag(Lambda_(1+nvar)); end; %<-- rectify spectrum of A. ;
Lambda__ = diag(Lambda_);
A_tru__ = real(Psi__*Lambda__*Psi_inv__);
C_tru__  = 1.0*randn(n_var,n_var); CC_tru__ = transpose(C_tru__)*C_tru__;
B_tru__  = 1.0*randn(n_var,n_var); BB_tru__ = transpose(B_tru__)*B_tru__;
A_inv_a_ = inv(A_tru__)*a_tru_;
tmp_n_T = 1024*1;
R_diff = @(log_amplitude) evolve_SDE(rseed,tmp_n_T,n_var,dt_avg,A_inv_a_,Psi__,Lambda_,Psi_inv__,B_tru__,C_tru__,snr,log_amplitude) - relative_variation ;
log_amplitude = fzero(R_diff,1.0);
R_est = evolve_SDE(rseed,tmp_n_T,n_var,dt_avg,A_inv_a_,Psi__,Lambda_,Psi_inv__,B_tru__,C_tru__,snr,log_amplitude);
if (flag_verbose); disp(sprintf(' %% estimated log_amplitude %0.4f --> amplitude %0.4f --> R_est %0.4f',log_amplitude,exp(log_amplitude),R_est)); end;
%%%%%%%%;
[R_avg,dt_,age_,Y__] = evolve_SDE(rseed,n_T_0in,n_var,dt_avg,A_inv_a_,Psi__,Lambda_,Psi_inv__,B_tru__,C_tru__,snr,log_amplitude);
amplitude = exp(log_amplitude);
%%%%%%%%;
aid_ = ones(n_T_0in,1);
dat__ = transpose(Y__);
if (flag_verbose); plot(Y__(1,:),Y__(2,:),'.-'); end;
[n_smp,n_var] = size(dat__);
string_dat_name_ = {'x';'y'};
dt_all_ = diff(age_);
n_step = 1;
dt_lim_ = [0;max(age_(1+n_step:end)-age_(1:end-n_step))];
%%%%%%%%;
save(fname_dat_mat ...
     , 'rseed', 'n_var', 'dt_avg', 'n_T_0in', 'relative_variation', 'snr' ...
     , 'a_tru_', 'A_tru__', 'C_tru__', 'CC_tru__', 'B_tru__', 'BB_tru__' ...
     , 'A_inv_a_', 'R_diff', 'log_amplitude', 'R_est', 'R_avg', 'dt_', 'age_', 'aid_', 'dat__' ...
     );
writematrix(age_,fname_age_csv);
writematrix(dat__,fname_dat_csv);
%%%%%%%%;
end;%if (~exist(fname_dat_mat,'file'));
load(fname_dat_mat);

%%%%%%%%;
% Now estimate drift. ;
%%%%%%%%;
dat_nrm__ = dat__;
[a_drift_,b_drift_] = dolphin_estimate_ab_0(aid_,age_,dat_nrm__);
a_est_ = a_drift_;
b_est_ = b_drift_;
amplitude = exp(log_amplitude);
corr_a = corr(a_est_(:),amplitude*a_tru_(:));
dat_nrm__ = dat_nrm__ - (ones(numel(age_),1)*reshape(b_drift_,[1,n_var]) + age_*reshape(a_drift_,[1,n_var]));

if ( exist(fname_c0_out_mat,'file')); disp(sprintf(' %% %s found, not creating',fname_c0_out_mat)); end;
if (~exist(fname_c0_out_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_c0_out_mat));
%%%%%%%%;
% then estimate SDE for parameter_c0;
%%%%%%%%;
parameter_c0 = struct('type','parameter');
parameter_c0.flag_verbose = 0;
parameter_c0.n_step = 1;
parameter_c0.n_iteration = 4;
parameter_c0.flag_constrain_CC = 0;
parameter_c0.n_shuffle = 32;
parameter_c0.flag_crude = 0;
%%%%%%%%;
parameter = parameter_c0;
dolphin_test_6_helper;
save(fname_c0_out_mat ...
     , 'tru_', 'est_', 'null_' ...
     );
%%%%%%%%;
end;%if (~exist(fname_c0_out_mat,'file'));
load(fname_c0_out_mat);

if ( exist(fname_c1_out_mat,'file')); disp(sprintf(' %% %s found, not creating',fname_c1_out_mat)); end;
if (~exist(fname_c1_out_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_c1_out_mat));
%%%%%%%%;
% then estimate SDE for parameter_c1;
%%%%%%%%;
parameter_c1 = struct('type','parameter');
parameter_c1.flag_verbose = 0;
parameter_c1.n_step = 1;
parameter_c1.n_iteration = 4;
parameter_c1.flag_constrain_CC = 0;
parameter_c1.n_shuffle = 32;
parameter_c1.flag_crude = 1;
%%%%%%%%;
parameter = parameter_c1;
dolphin_test_6_helper;
save(fname_c1_out_mat ...
     , 'tru_', 'est_', 'null_' ...
     );
%%%%%%%%;
end;%if (~exist(fname_c1_out_mat,'file'));
load(fname_c1_out_mat);

%{

%%%%%%%%;
% Generate irg R file. ;
%%%%%%%%;
if ( exist(fname_irg_R,'file')); disp(sprintf(' %% %s found, not creating',fname_irg_R)); end;
if (~exist(fname_irg_R,'file')); disp(sprintf(' %% %s not found, creating',fname_irg_R));
%%%%%%%%;
fp = fopen(fname_irg_R,'w');
fprintf(fp,'# install.packages("remotes");\n');
fprintf(fp,'# remotes::install_github("SMAC-Group/irg");\n');
fprintf(fp,'setwd("%s");\n',dir_R);
fprintf(fp,'library(irg);\n');
fprintf(fp,'library(mgcv);\n');
fprintf(fp,'################;\n');
fprintf(fp,'# Average signals at matching time points and then standardize;\n');
fprintf(fp,'avg <- function(signals) {\n');
fprintf(fp,'  signal <- rep(NA, length(signals)/3)\n');
fprintf(fp,'  for(i in 1:(length(signals)/3)) {\n');
fprintf(fp,'    signal[i] <- mean(signals[(3*(i-1) + 1) : (3*i)])\n');
fprintf(fp,'  }\n');
fprintf(fp,'  signal_std <- (signal - mean(signal))/sd(signal)\n');
fprintf(fp,'  return(signal_std)\n');
fprintf(fp,'}\n');
fprintf(fp,'# Pre-processing function uses GAM to detrend the signals and then uses previous function to aggregate and normalize\n');
fprintf(fp,'pre_process_single <- function(dat_0_,age_) {\n');
fprintf(fp,'  mod_0 <- gam(dat_0_ ~ s(age_))\n');
fprintf(fp,'  resid_0 <- avg(resid(mod_0))\n');
fprintf(fp,'  return(resid_0)\n');
fprintf(fp,'}\n');
fprintf(fp,'################;\n');
fprintf(fp,'age_ <- read.table("%s",sep=",");\n',fname_age_csv);
fprintf(fp,'dat__ <- read.table("%s",sep=",");\n',fname_dat_csv);
fprintf(fp,'age_nrm_ <- rep(NA,length(age_[,1])/3);\n\n');
fprintf(fp,'for(i in 1:(length(age_[,1])/3)) {\n');
fprintf(fp,'  age_nrm_[i] <- mean(age_[(3*(i-1) + 1) : (3*i),1])\n');
fprintf(fp,'}#for(i in 1:(length(age_[,1])/3));\n');
fprintf(fp,'################;\n');
fprintf(fp,'n_smp <- dim(dat__)[1];\n');
fprintf(fp,'n_var <- dim(dat__)[2];\n');
fprintf(fp,'################;\n');
fprintf(fp,'for (nvar0 in 0:(n_var-1)){');
fprintf(fp,'if (nvar0<n_var-1){');
fprintf(fp,'for (nvar1 in (nvar0+1):(n_var-1)){');
fprintf(fp,'################;\n');
fprintf(fp,'string_out <- c(nvar0,nvar1);\n');
fprintf(fp,'write(string_out,"");\n');
fprintf(fp,'dat_0_ <- dat__[,1+nvar0];\n');
fprintf(fp,'dat_1_ <- dat__[,1+nvar1];\n');
fprintf(fp,'################;\n');
fprintf(fp,'dat_nrm_0_ <- pre_process_single(dat_0_,age_[,1]);\n');
fprintf(fp,'dat_nrm_1_ <- pre_process_single(dat_1_,age_[,1]);\n');
fprintf(fp,'################;\n');
fprintf(fp,'dat_nrm__ <- data.frame(dat_nrm_0_,dat_nrm_1_,age_nrm_);\n');
fprintf(fp,'result_nrm_ <- granger_test(dat_nrm__[,1],dat_nrm__[,2],dat_nrm__[,3],theta=NULL,alternative="twodir",H=100,seed=0,showprogress=TRUE);\n');
fprintf(fp,'output_nrm_ <- c(nvar0,nvar1,result_nrm_$pvalue,result_nrm_$parameters);\n');
fprintf(fp,'write(output_nrm_,"%s",append=TRUE);\n',fname_irg_out);
fprintf(fp,'################;\n');
fprintf(fp,'}#for (nvar1 in (nvar0+1):(n_var-1));\n');
fprintf(fp,'}#if (nvar0<n_var-1);\n');
fprintf(fp,'}#for (nvar0 in 0:(n_var-1));\n');
fclose(fp);
%%%%%%%%;
end;%if (~exist(fname_irg_R,'file'));

%%%%%%%%;
% Generate irg R file. ;
%%%%%%%%;
if ( exist(fname_irg_out,'file')); disp(sprintf(' %% %s found, not creating',fname_irg_out)); end;
if (~exist(fname_irg_out,'file')); disp(sprintf(' %% %s not found, creating',fname_irg_out));
string_command = sprintf(' Rscript %s;',fname_irg_R);
system(string_command);
end;%if (~exist(fname_irg_out,'file'));

%}

delete(fname_irg_tmp);
end;%if (~flag_skip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function string_infix = dolphin_test_6_infix(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
str_dt = sprintf('dt%.4d',floor(1000*dt_avg));
str_n_T = sprintf('T%.4d',floor(n_T_0in));
str_rv = sprintf('rv%.4d',floor(1000*relative_variation));
str_snr = sprintf('snr%.4d',floor(1000*snr));
string_infix = sprintf('dt6_rng%d_n%d_%s_%s_%s_%s',rseed,n_var,str_dt,str_n_T,str_rv,str_snr);

function [R_avg,dt_,age_,Y__] = evolve_SDE(rseed,n_T,n_var,dt_avg,A_inv_a_,Psi__,Lambda_,Psi_inv__,B__,C__,snr,log_amplitude);
flag_verbose=0;
if flag_verbose; disp(sprintf(' %% [entering evolve_SDE] log_amplitude %0.2f --> amplitude %0.2f',log_amplitude,exp(log_amplitude))); end;
rng(rseed);
Y__ = zeros(n_var,n_T);
Y__(:,1+0) = randn(n_var,1);
age_ = zeros(n_T,1);
dt_ = zeros(n_T-1,1);
R_ = zeros(n_T-1,1);
for nT=1:n_T-1;
dt=-dt_avg * log(1-rand());
dt_(1+nT-1) = dt;
age_(1+nT) = age_(1+nT-1) + dt;
Y_pre_ = Y__(:,1+nT-1);
tmp_exp_Adt__ = Psi__ * diag(exp(+exp(log_amplitude)*Lambda_*dt)) * Psi_inv__;
Y_pos_ = real( -A_inv_a_ + tmp_exp_Adt__ * (Y_pre_ + A_inv_a_ + SDE_sample_int_expAsBdW(dt,Psi__,exp(log_amplitude)*Lambda_,exp(log_amplitude)/snr*Psi_inv__*B__,1e-2)) );
Y__(:,1+nT+0) = Y_pos_ + exp(log_amplitude)/snr*C__*randn(n_var,1);
R_(1+nT-1) = fnorm(Y__(:,1+nT+0) - Y__(:,1+nT-1))/fnorm(Y__(:,1+nT-1));
end;%for nT=1:n_T-1;
R_avg = mean(R_);
if flag_verbose; disp(sprintf(' %% [finished evolve_SDE] R_avg %0.2f',R_avg)); end;
