% test out a loop-counting message-passing algorithm. ;
clear;
rng(1);
M = 2048*3; N = M+1; snr = 0.95; n_cluster = 1;
[A_n_,label_A_,n_label_A_,pf_,pi_,snr_,MX_,NX_,mu_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);
D_n_base_ = A_n_(pi_{1},pi_{2});
D_n_ = D_n_base_;
B_n_X_ = zeros(M,N); B_n_X_(1:MX_,1:NX_) = 1;
B_n_ij_ = find(B_n_X_(:)); n_B_n = numel(B_n_ij_);
D_n_ij_ = setdiff(1:M*N,B_n_ij_); n_D_n = numel(D_n_ij_);

n_tau = 16; tau_ = 2.^[1:n_tau];
n_h = 128;
n_iteration = 16;
A_ = zeros(n_tau,n_iteration);
for niteration=1:n_iteration;
%%%%%%%%;
tmp_Lr_ = repmat(sum(D_n_.^2,2),1,N); tmp_Lc_ = repmat(sum(D_n_.^2,1),M,1);
tmp_DDDD_ = D_n_.*(D_n_*transpose(D_n_)*D_n_) - D_n_.^2 .* (tmp_Lr_+tmp_Lc_) + D_n_.^4;
E_n_ = tmp_DDDD_; E_n_ = (E_n_-mean(E_n_))/std(E_n_(:));
for ntau=1:n_tau;
tau = tau_(ntau);
tmp_F_n_ = exp( (E_n_ - mean(E_n_(:)))/tau ); tmp_F_n_ = (tmp_F_n_-mean(tmp_F_n_))/std(tmp_F_n_(:));
A_(ntau,niteration) = auc_0(tmp_F_n_(D_n_ij_),tmp_F_n_(B_n_ij_));
clear tmp_F_n_;
end;%for ntau=1:n_tau;
%%%%%%%%;
sg_D_n_ = std(D_n_(:)); ylim_D_n_ = 0.0 + 7.5*sg_D_n_*[-1,+1];
h_bins_D_n_ = linspace(ylim_D_n_(1),ylim_D_n_(2),n_h);
h_D_n_sub_ = hist(D_n_(B_n_ij_),h_bins_D_n_); h_D_n_all_ = hist(D_n_(D_n_ij_),h_bins_D_n_);
%%%%%%%%;
sg_E_n_ = std(E_n_(:)); ylim_E_n_ = 0.0 + 7.5*sg_E_n_*[-1,+1];
h_bins_E_n_ = linspace(ylim_E_n_(1),ylim_E_n_(2),n_h);
h_E_n_sub_ = hist(E_n_(B_n_ij_),h_bins_E_n_); h_E_n_all_ = hist(E_n_(D_n_ij_),h_bins_E_n_);
%%%%%%%%;
flag_plot=1;
if flag_plot;
subplot(1,3,1); cla; hold on;
stairs(h_bins_D_n_,h_D_n_all_/sum(h_D_n_all_),'b-');
stairs(h_bins_D_n_,h_D_n_sub_/sum(h_D_n_sub_),'r-');
subplot(1,3,2); cla; hold on;
stairs(h_bins_E_n_,h_E_n_all_/sum(h_E_n_all_),'b-');
stairs(h_bins_E_n_,h_E_n_sub_/sum(h_E_n_sub_),'r-');
subplot(1,3,3); cla; hold on;
plot(log2(tau_),A_(:,niteration)); xlabel('log2(tau)'); ylabel('auc');
figbig;drawnow;pause();
end;%if flag_plot;
%%%%%%%%;
D_n_ = E_n_;
%%%%%%%%;
end;%for niteration=1:n_iteration;

