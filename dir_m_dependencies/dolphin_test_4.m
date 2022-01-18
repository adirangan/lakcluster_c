function [corr_a,corr_A,corr_BB,corr_CC,tru_,est_,null_] = dolphin_test_4(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
% Assumes exponential distribution of dt with mean dt_avg. ;
%{
%try: ;

clear;
setup_OptiPlex;
dir_trunk = '/home/rangan/dir_bcc/dir_dolphin';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);
dir_datmat = sprintf('%s/dir_datmat',dir_trunk);

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
[corr_a,corr_A,corr_BB,corr_CC] = dolphin_test_4(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
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
[corr_a,corr_A,corr_BB,corr_CC] = dolphin_test_4(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
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
[corr_a,corr_A,corr_BB,corr_CC] = dolphin_test_4(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
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

fname_fig = sprintf('%s/dolphin_test_4_FIGA',dir_jpg);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;
c_ = colormap_beach(); n_c = size(c_,1);
for n_var = [2,4,8,16,32];
nc = max(0,min(n_c-1,floor(n_c*log2(n_var)/5)));
tmp_ = load(sprintf('../dir_mat/dolphin_test_aid00n%d.mat',n_var));
subplot(2,2,1); hold on; plot(1:tmp_.n_snr,median(tmp_.corr_a__,2),'.-','Color',c_(1+nc,:),'LineWidth',2,'MarkerSize',25); hold off;
subplot(2,2,2); hold on; plot(1:tmp_.n_snr,median(tmp_.corr_A__,2),'.-','Color',c_(1+nc,:),'LineWidth',2,'MarkerSize',25); hold off;
subplot(2,2,3); hold on; plot(1:tmp_.n_snr,median(tmp_.corr_BB__,2),'.-','Color',c_(1+nc,:),'LineWidth',2,'MarkerSize',25); hold off;
subplot(2,2,4); hold on; plot(1:tmp_.n_snr,median(tmp_.corr_CC__,2),'.-','Color',c_(1+nc,:),'LineWidth',2,'MarkerSize',25); hold off;
end;%for n_var = [2,4,8,16,32];
subplot(2,2,1); xlabel('snr'); set(gca,'XTick',1:tmp_.n_snr,'XTickLabel',tmp_.snr_); xtickangle(90); ylabel('correlation'); xlim([0,tmp_.n_snr+1]); ylim([-0.25,+1.25]); grid on; title('a'); legend({'n2','n4','n8','n16','n32'},'Location','SouthEast');
subplot(2,2,2); xlabel('snr'); set(gca,'XTick',1:tmp_.n_snr,'XTickLabel',tmp_.snr_); xtickangle(90); ylabel('correlation'); xlim([0,tmp_.n_snr+1]); ylim([-0.25,+1.25]); grid on; title('A'); legend({'n2','n4','n8','n16','n32'},'Location','SouthEast');
subplot(2,2,3); xlabel('snr'); set(gca,'XTick',1:tmp_.n_snr,'XTickLabel',tmp_.snr_); xtickangle(90); ylabel('correlation'); xlim([0,tmp_.n_snr+1]); ylim([-0.25,+1.25]); grid on; title('BB'); legend({'n2','n4','n8','n16','n32'},'Location','SouthEast');
subplot(2,2,4); xlabel('snr'); set(gca,'XTick',1:tmp_.n_snr,'XTickLabel',tmp_.snr_); xtickangle(90); ylabel('correlation'); xlim([0,tmp_.n_snr+1]); ylim([-0.25,+1.25]); grid on; title('CC'); legend({'n2','n4','n8','n16','n32'},'Location','SouthEast');
sgtitle('Correlation of recovered SDE parameters');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
clear tmp_;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% generate some data. ;
%%%%%%%%;
clear;
setup_access1;
dir_trunk = '/data/rangan/dir_bcc/dir_dolphin';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);
dir_mat = sprintf('%s/dir_mat',dir_trunk);
dir_datmat = sprintf('%s/dir_datmat',dir_trunk);
flag_verbose=0;
dXX = 0;
for n_var = [2,4,8,16,32]; %n_var = 8; %<-- more than 2 variables. ;
dt_avg = 0.257; %<-- matches d00 data. ;
n_T_0in = 5664; %<-- matches d00 data. ;
relative_variation = 0.8973; %<-- matches d00 data. ;
snr_ = 2.^[-3:1:+3]; n_snr = numel(snr_);
rseed=0;
for nsnr=0:n_snr-1;
snr = snr_(1+nsnr);
fname_base = sprintf('%s/dolphin_test_gen_aid%.2dn%dls%drs%d',dir_datmat,dXX,n_var,3+log2(snr),rseed);
fname_age_csv = sprintf('%s_age_.csv',fname_base);
fname_dat__csv = sprintf('%s_dat__.csv',fname_base);
fname_a_csv = sprintf('%s_a_.csv',fname_base);
fname_A__csv = sprintf('%s_A__.csv',fname_base);
fname_BB__csv = sprintf('%s_BB__.csv',fname_base);
fname_CC__csv = sprintf('%s_CC__.csv',fname_base);
fname_mat = sprintf('%s.mat',fname_base);
if (~exist(fname_mat,'file') | ~exist(fname_dat__csv,'file'));
disp(sprintf(' %s not found, creating',fname_mat));
%%%%%%%%;
rng(rseed);
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
writematrix(age_,fname_age_csv);
writematrix(dat__,fname_dat__csv);
writematrix(a_tru_,fname_a_csv);
writematrix(A_tru__,fname_A__csv); 
writematrix(BB_tru__,fname_BB__csv); 
writematrix(CC_tru__,fname_CC__csv); 
save(fname_mat ...
,'n_smp','n_var','age_','dat__' ...
,'a_tru_','A_tru__','BB_tru__','CC_tru__' ...
);
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %s found, not creating',fname_mat));
end;%if ( exist(fname_mat,'file'));
end;%for nsnr=0:n_snr-1;
end;%for n_var = [2,4,8,16,32]; %n_var = 8; %<-- more than 2 variables. ;


%}
flag_verbose = 0;

setup_OptiPlex;
dir_trunk = '/home/rangan/dir_bcc/dir_dolphin';
dir_jpg = sprintf('%s/dir_jpg',dir_trunk);

flag_replot=1;
flag_crude=1;
flag_a_is_0=0;

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
C_tru__  = 1.0*randn(n_var,n_var); CC_tru__ = transpose(C_tru__)*C_tru__;
B_tru__  = 1.0*randn(n_var,n_var); BB_tru__ = transpose(B_tru__)*B_tru__;
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

n_iteration = 32;
[a_est_,A_est__,BB_est__,CC_est__,L_,niteration,a_est__,A_est___,BB_est___,CC_est___] = dolphin_estimate_aABC_0(aid_,age_,dat__,n_step,n_iteration);

[~,tmp_index_min] = min(L_(1:1+niteration)); tmp_index_min = tmp_index_min - 1;
if (flag_verbose); disp(sprintf(' %% tmp_index_min %d/%d --> L %0.2f',tmp_index_min,niteration,L_(1+tmp_index_min))); end;
a_est_ = a_est__{1+tmp_index_min};
A_est__ = A_est___{1+tmp_index_min};
BB_est__ = BB_est___{1+tmp_index_min};
CC_est__ = CC_est___{1+tmp_index_min};
if (flag_verbose>2);
  a_est_,;amplitude*a_tru_,;
  A_est__,;amplitude*A_tru__,;
  BB_est__,;(amplitude/snr)^2*BB_tru__,;
  CC_est__,;(amplitude/snr)^2*CC_tru__,;
end;%if (flag_verbose);
corr_a = corr(a_est_(:),amplitude*a_tru_(:));
corr_A = corr(A_est__(:),amplitude*A_tru__(:));
corr_BB = corr(BB_est__(:),(amplitude/snr)^2*BB_tru__(:));
corr_CC = corr(CC_est__(:),(amplitude/snr)^2*CC_tru__(:));

n_shuffle = 32;
a_null__ = cell(n_shuffle,1);
A_null___ = cell(n_shuffle,1);
BB_null___ = cell(n_shuffle,1);
CC_null___ = cell(n_shuffle,1);
n_step = 1; n_iteration = 4;
for nshuffle=1:n_shuffle;
dat_prm__ = dolphin_permute_0(aid_,age_,dat__,nshuffle);
[a_,A__,BB__,CC__,L_,niteration,a__,A___,BB___,CC___] = dolphin_estimate_aABC_0(aid_,age_,dat_prm__,n_step,n_iteration);
[~,tmp_index_min] = min(L_(1:1+niteration)); tmp_index_min = tmp_index_min-1;
a_null__{1+nshuffle-1} = a__{1+tmp_index_min};
A_null___{1+nshuffle-1} = A___{1+tmp_index_min};
BB_null___{1+nshuffle-1} = BB___{1+tmp_index_min};
CC_null___{1+nshuffle-1} = CC___{1+tmp_index_min};
L_null_{1+nshuffle-1} = L_(1+tmp_index_min);
end;%for nshuffle=1:n_shuffle;

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
%%%%%%%%;
for nshuffle=0:n_shuffle-1;
if flag_a_is_0==0; a_null_avg_ = a_null_avg_ + a_null__{1+nshuffle}; end;
A_null_avg__ = A_null_avg__ + A_null___{1+nshuffle};
BB_null_avg__ = BB_null_avg__ + BB_null___{1+nshuffle};
CC_null_avg__ = CC_null_avg__ + CC_null___{1+nshuffle};
if flag_a_is_0==0; a_null_std_ = a_null_std_ + a_null__{1+nshuffle}.^2; end;
A_null_std__ = A_null_std__ + A_null___{1+nshuffle}.^2;
BB_null_std__ = BB_null_std__ + BB_null___{1+nshuffle}.^2;
CC_null_std__ = CC_null_std__ + CC_null___{1+nshuffle}.^2;
end;%for nshuffle=0:n_shuffle-1;
%%%%%%%%;
if flag_a_is_0==0; a_null_avg_ = a_null_avg_ / n_shuffle; end;
A_null_avg__ = A_null_avg__ / n_shuffle;
BB_null_avg__ = BB_null_avg__ / n_shuffle;
CC_null_avg__ = CC_null_avg__ / n_shuffle;
if flag_a_is_0==0; a_null_std_ = sqrt(a_null_std_/n_shuffle - a_null_avg_.^2); end;
A_null_std__ = sqrt(A_null_std__/n_shuffle - A_null_avg__.^2);
BB_null_std__ = sqrt(BB_null_std__/n_shuffle - BB_null_avg__.^2);
CC_null_std__ = sqrt(CC_null_std__/n_shuffle - CC_null_avg__.^2);
%%%%%%%%;
if flag_a_is_0==0; a_Z_ = real((a_est_ - a_null_avg_)./a_null_std_); end;
A_Z__ = real((A_est__ - A_null_avg__)./A_null_std__);
BB_Z__ = real((BB_est__ - BB_null_avg__)./BB_null_std__);
CC_Z__ = real((CC_est__ - CC_null_avg__)./CC_null_std__);
if flag_a_is_0==0; a_nlp_ = -z_to_lp(a_Z_); end;
A_nlp__ = -z_to_lp(A_Z__);
BB_nlp__ = -z_to_lp(BB_Z__);
CC_nlp__ = -z_to_lp(CC_Z__);

tru_ = struct('type','tru');
tru_.a_tru_ = a_tru_;
tru_.A_tru__ = A_tru__;
tru_.BB_tru__ = BB_tru__;
tru_.CC_tru__ = CC_tru__;

est_ = struct('type','est');
est_.a_est_ = a_est_;
est_.A_est__ = A_est__;
est_.BB_est__ = BB_est__;
est_.CC_est__ = CC_est__;

null_ = struct('type','null');
null_.a_null__ = a_null__;
null_.A_null___ = A_null___;
null_.BB_null___ = BB_null___;
null_.CC_null___ = CC_null___;
null_.a_null_avg_ = a_null_avg_;
null_.A_null_avg__ = A_null_avg__;
null_.BB_null_avg__ = BB_null_avg__;
null_.CC_null_avg__ = CC_null_avg__;
null_.a_null_std_ = a_null_std_;
null_.A_null_std__ = A_null_std__;
null_.BB_null_std__ = BB_null_std__;
null_.CC_null_std__ = CC_null_std__;
null_.a_Z_ = a_Z_;
null_.A_Z__ = A_Z__;
null_.BB_Z__ = BB_Z__;
null_.CC_Z__ = CC_Z__;
null_.a_nlp_ = a_nlp_;
null_.A_nlp__ = A_nlp__;
null_.BB_nlp__ = BB_nlp__;
null_.CC_nlp__ = CC_nlp__;



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
