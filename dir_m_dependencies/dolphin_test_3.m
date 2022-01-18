function [corr_a,corr_A,corr_BB,corr_CC] = dolphin_test_3(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
% Assumes exponential distribution of dt with mean dt_avg. ;
%{
%try: ;

clear;
setup_OptiPlex;
dir_trunk = '/home/rangan/dir_bcc/dir_dolphin';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);

dXX = 0;
n_var = 2; %<-- pair. ;
dt_avg = 0.257; %<-- matches d00 data. ;
n_T_0in = 5664; %<-- matches d00 data. ;
relative_variation = 0.8973; %<-- matches d00 data. ;
snr_ = 2.^[-3:1:+3]; n_snr = numel(snr_);
rseed_ = 0:8; n_rseed = numel(rseed_);
fname_mat = sprintf('%s/dolphin_test_aid%.2dn%d.mat',dir_mat,dXX,n_var);
if (~exist(fname_mat,'file'));
disp(sprintf(' %s not found, creating',fname_mat));
corr_a__ = zeros(n_snr,n_rseed);
corr_A__ = zeros(n_snr,n_rseed);
corr_BB__ = zeros(n_snr,n_rseed);
corr_CC__ = zeros(n_snr,n_rseed);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
for nsnr=0:n_snr-1;
snr = snr_(1+nsnr);
[corr_a,corr_A,corr_BB,corr_CC] = dolphin_test_3(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
corr_a__(1+nsnr,1+nrseed) = corr_a;
corr_A__(1+nsnr,1+nrseed) = corr_A;
corr_BB__(1+nsnr,1+nrseed) = corr_BB;
corr_CC__(1+nsnr,1+nrseed) = corr_CC;
disp(sprintf(' %% nsnr %d/%d nrseed %d/%d --> corr_a %0.2f corr_A %0.2f corr_BB %0.2f corr_CC %0.2f',nsnr,n_snr,nrseed,n_rseed,corr_a,corr_A,corr_BB,corr_CC));
end;%for nsnr=0:n_snr-1;
end;%for nrseed=0:n_rseed-1;
save(fname_mat ...
,'n_var','dt_avg','n_T_0in','relative_variation','snr_','n_snr','rseed_','n_rseed' ...
,'corr_a__' ...
,'corr_A__' ...
,'corr_BB__' ...
,'corr_CC__' ...
);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

dXX = 32;
n_var = 2; %<-- pair. ;
dt_avg = 0.1809; %<-- matches d32 data. ;
n_T_0in = 117; %<-- matches d32 data. ;
relative_variation = 0.9365; %<-- matches d32 data. ;
snr_ = 2.^[-3:1:+3]; n_snr = numel(snr_);
rseed_ = 0:8; n_rseed = numel(rseed_);
fname_mat = sprintf('%s/dolphin_test_aid%.2dn%d.mat',dir_mat,dXX,n_var);
if (~exist(fname_mat,'file'));
disp(sprintf(' %s not found, creating',fname_mat));
corr_a__ = zeros(n_snr,n_rseed);
corr_A__ = zeros(n_snr,n_rseed);
corr_BB__ = zeros(n_snr,n_rseed);
corr_CC__ = zeros(n_snr,n_rseed);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
for nsnr=0:n_snr-1;
snr = snr_(1+nsnr);
[corr_a,corr_A,corr_BB,corr_CC] = dolphin_test_3(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
corr_a__(1+nsnr,1+nrseed) = corr_a;
corr_A__(1+nsnr,1+nrseed) = corr_A;
corr_BB__(1+nsnr,1+nrseed) = corr_BB;
corr_CC__(1+nsnr,1+nrseed) = corr_CC;
disp(sprintf(' %% nsnr %d/%d nrseed %d/%d --> corr_a %0.2f corr_A %0.2f corr_BB %0.2f corr_CC %0.2f',nsnr,n_snr,nrseed,n_rseed,corr_a,corr_A,corr_BB,corr_CC));
end;%for nsnr=0:n_snr-1;
end;%for nrseed=0:n_rseed-1;
save(fname_mat ...
,'n_var','dt_avg','n_T_0in','relative_variation','snr_','n_snr','rseed_','n_rseed' ...
,'corr_a__' ...
,'corr_A__' ...
,'corr_BB__' ...
,'corr_CC__' ...
);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

dXX = 0;
for n_var = [2,4,8,16,32]; %n_var = 8; %<-- more than 2 variables. ;
dt_avg = 0.257; %<-- matches d00 data. ;
n_T_0in = 5664; %<-- matches d00 data. ;
relative_variation = 0.8973; %<-- matches d00 data. ;
snr_ = 2.^[-3:1:+3]; n_snr = numel(snr_);
rseed_ = 0:8; n_rseed = numel(rseed_);
fname_mat = sprintf('%s/dolphin_test_aid%.2dn%d.mat',dir_mat,dXX,n_var);
if (~exist(fname_mat,'file'));
disp(sprintf(' %s not found, creating',fname_mat));
corr_a__ = zeros(n_snr,n_rseed);
corr_A__ = zeros(n_snr,n_rseed);
corr_BB__ = zeros(n_snr,n_rseed);
corr_CC__ = zeros(n_snr,n_rseed);
for nrseed=0:n_rseed-1;
rseed = rseed_(1+nrseed);
for nsnr=0:n_snr-1;
snr = snr_(1+nsnr);
[corr_a,corr_A,corr_BB,corr_CC] = dolphin_test_3(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
corr_a__(1+nsnr,1+nrseed) = corr_a;
corr_A__(1+nsnr,1+nrseed) = corr_A;
corr_BB__(1+nsnr,1+nrseed) = corr_BB;
corr_CC__(1+nsnr,1+nrseed) = corr_CC;
disp(sprintf(' %% nsnr %d/%d nrseed %d/%d --> corr_a %0.2f corr_A %0.2f corr_BB %0.2f corr_CC %0.2f',nsnr,n_snr,nrseed,n_rseed,corr_a,corr_A,corr_BB,corr_CC));
end;%for nsnr=0:n_snr-1;
end;%for nrseed=0:n_rseed-1;
save(fname_mat ...
,'n_var','dt_avg','n_T_0in','relative_variation','snr_','n_snr','rseed_','n_rseed' ...
,'corr_a__' ...
,'corr_A__' ...
,'corr_BB__' ...
,'corr_CC__' ...
);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
end;%for n_var = [2,4,8,16,32]; %n_var = 8; %<-- more than 2 variables. ;

%}
flag_verbose = 0;

setup_OptiPlex;
dir_trunk = '/home/rangan/dir_bcc/dir_dolphin';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);

flag_replot=1;
flag_crude=1;

%%%%%%%%;
% Generate data. ;
%%%%%%%%;
rng(rseed);

a_tru_ = 1.0*randn(n_var,1);
A_tru__ = randn(n_var,n_var);
[Psi__,Lambda__] = eig(A_tru__); Psi_inv__ = inv(Psi__);
Lambda_ = diag(Lambda__);
for nvar=0:n_var-1; Lambda_(1+nvar) = min(-0.1,real(Lambda_(1+nvar))) + i*imag(Lambda_(1+nvar)); end; %<-- rectify spectrum of A. ;
Lambda__ = diag(Lambda_);
A_tru__ = real(Psi__*Lambda__*Psi_inv__);
C_tru__  = 1.0*randn(n_var,n_var); CC_tru__ = C_tru__*transpose(C_tru__);
B_tru__  = 1.0*randn(n_var,n_var); BB_tru__ = B_tru__*transpose(B_tru__);
A_inv_a_ = inv(A_tru__)*a_tru_;
tmp_n_T = 1024*1;
R_diff = @(log_amplitude) evolve_SDE(rseed,tmp_n_T,n_var,dt_avg,A_inv_a_,Psi__,Lambda_,Psi_inv__,B_tru__,C_tru__,snr,log_amplitude) - relative_variation ;
log_amplitude = fzero(R_diff,1.0);
R_est = evolve_SDE(rseed,tmp_n_T,n_var,dt_avg,A_inv_a_,Psi__,Lambda_,Psi_inv__,B_tru__,C_tru__,snr,log_amplitude);
if (flag_verbose); disp(sprintf(' %% estimated log_amplitude %0.4f --> amplitude %0.4f --> R_est %0.4f',log_amplitude,exp(log_amplitude),R_est)); end;

[R_avg,dt_,age_,Y__] = evolve_SDE(rseed,n_T_0in,n_var,dt_avg,A_inv_a_,Psi__,Lambda_,Psi_inv__,B_tru__,C_tru__,snr,log_amplitude);
amplitude = exp(log_amplitude);

aid_ = ones(n_T_0in,1);
dat__ = transpose(Y__);
if (flag_verbose); plot(Y__(1,:),Y__(2,:),'.-'); end;
[n_smp,n_var] = size(dat__);
string_dat_name_ = {'x';'y'};
dt_all_ = diff(age_);
n_step = 1;
dt_lim_ = [0;max(age_(1+n_step:end)-age_(1:end-n_step))];

a_est_ = zeros(n_var,1); A_est__ = zeros(n_var,n_var);
%%%%%%%%;
% Now estimate initial B,C. ;
%%%%%%%%;
[BB_est__,CC_est__,l2_R__,sum_1,sum_dt,sum_dtdt,sum_DDj__,sum_DDjdt__] = dolphin_estimate_BC_from_aA_0(aid_,age_,dat__,a_est_,A_est__,n_step);
[BB_crude__,CC_crude__] = dolphin_estimate_BC_from_aA_crude_0(aid_,age_,dat__,a_est_,A_est__,n_step);
if (flag_verbose); disp(sprintf(' %% BB_est__ vs BB_crude__: %0.16f',fnorm(BB_est__ - BB_crude__)/fnorm(BB_est__))); end;
if (flag_verbose); disp(sprintf(' %% CC_est__ vs CC_crude__: %0.16f',fnorm(CC_est__ - CC_crude__)/fnorm(CC_est__))); end;
if (flag_crude); BB_est__ = BB_crude__; CC_est__ = CC_crude__; end; %<-- use crude method. ;
%%%%%%%%;
% Now estimate a_est_ and A_est__. ;
%%%%%%%%;
[a_est_,A_est__,L] = dolphin_estimate_aA_from_BC_0(aid_,age_,dat__,BB_est__,CC_est__,dt_lim_,n_step);
if (flag_verbose); disp(sprintf(' %% initial: negative-log-likelihood %0.16f',L)); end;
%%%%%%%%;
% Now iterate a few times. ;
%%%%%%%%;
n_iteration = 4;
L_ = zeros(n_iteration+1,1);
a_est__ = cell(1+n_iteration,1);
A_est___ = cell(1+n_iteration,1);
BB_est___ = cell(1+n_iteration,1);
CC_est___ = cell(1+n_iteration,1);
niteration=0;
L_old = L;
L_(1+niteration) = L_old;
a_est__{1+niteration} = a_est_;
A_est___{1+niteration} = A_est__;
BB_est___{1+niteration} = BB_est__;
CC_est___{1+niteration} = CC_est__;
flag_continue=1;
while (flag_continue);
% Re-estimate B,C. ;
[BB_est__,CC_est__,l2_R__,sum_1,sum_dt,sum_dtdt,sum_DDj__,sum_DDjdt__] = dolphin_estimate_BC_from_aA_0(aid_,age_,dat__,a_est_,A_est__,n_step);
[BB_crude__,CC_crude__] = dolphin_estimate_BC_from_aA_crude_0(aid_,age_,dat__,a_est_,A_est__,n_step);
if (flag_verbose); disp(sprintf(' %% BB_est__ vs BB_crude__: %0.16f',fnorm(BB_est__ - BB_crude__)/fnorm(BB_est__))); end;
if (flag_verbose); disp(sprintf(' %% CC_est__ vs CC_crude__: %0.16f',fnorm(CC_est__ - CC_crude__)/fnorm(CC_est__))); end;
if (flag_crude); BB_est__ = BB_crude__; CC_est__ = CC_crude__; end; %<-- use crude method. ;
% Re-estimate a_est_ and A_est__ from BB_est__ and CC_est__. ;
[a_est_,A_est__,L] = dolphin_estimate_aA_from_BC_0(aid_,age_,dat__,BB_est__,CC_est__,dt_lim_,n_step);
L_new = L; L_(1+niteration+1)=L;
a_est__{1+niteration+1} = a_est_;
A_est___{1+niteration+1} = A_est__;
BB_est___{1+niteration+1} = BB_est__;
CC_est___{1+niteration+1} = CC_est__;
if (flag_verbose); disp(sprintf(' %% iteration %d: negative-log-likelihood %0.16f',niteration,L)); end;
f_a_ = fnorm(a_est_ - amplitude*a_tru_)/fnorm(amplitude*a_tru_);
f_A__ = fnorm(A_est__ - amplitude*A_tru__)/fnorm(amplitude*A_tru__);
f_BB__ = fnorm(BB_est__ - amplitude/snr*BB_tru__)/fnorm(amplitude/snr*BB_tru__);
f_CC__ = fnorm(CC_est__ - amplitude/snr*CC_tru__)/fnorm(amplitude/snr*CC_tru__);
if (flag_verbose); disp(sprintf(' %% iteration %d: f_a_ %0.2f f_A__ %0.2f f_BB__ %0.2f f_CC__ %0.2f',niteration,f_a_,f_A__,f_BB__,f_CC__)); end;
flag_continue=0;
niteration=niteration+1;
if (niteration<n_iteration & fnorm(L_old-L_new)/fnorm(L_old)>1e-3); flag_continue=1; end;
L_old = L_new;
end;%while (flag_continue);
%%%%%%%%;

[~,tmp_index_min] = min(L_(1:1+niteration)); tmp_index_min = tmp_index_min - 1;
if (flag_verbose); disp(sprintf(' %% tmp_index_min %d/%d --> L %0.2f',tmp_index_min,niteration,L_(1+tmp_index_min))); end;
a_est_ = a_est__{1+tmp_index_min};
A_est__ = A_est___{1+tmp_index_min};
BB_est__ = BB_est___{1+tmp_index_min};
CC_est__ = CC_est___{1+tmp_index_min};
if (flag_verbose>2);
  a_est_,;amplitude*a_tru_,;
  A_est__,;amplitude*A_tru__,;
  BB_est__,;amplitude/snr*BB_tru__,;
  CC_est__,;amplitude/snr*CC_tru__,;
end;%if (flag_verbose);
corr_a = corr(a_est_(:),amplitude*a_tru_(:));
corr_A = corr(A_est__(:),amplitude*A_tru__(:));
corr_BB = corr(BB_est__(:),amplitude/snr*BB_tru__(:));
corr_CC = corr(CC_est__(:),amplitude/snr*CC_tru__(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

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
