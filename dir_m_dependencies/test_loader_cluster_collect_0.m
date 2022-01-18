% collect results from test_loader_cluster_wrap_0.m ;

str_li30p = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/dexcluster_nonbinary_trace_ZRmax_g010_p150_nml10.mat';
str_li30p_cPhr_r1 = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cPhr_r1_cluster/dexcluster_nonbinary_trace_ZRmax_g010_p150_nml10.mat';

tmp_li30p_ = load(str_li30p);
tmp_li30p_cPhr_r1_ = load(str_li30p_cPhr_r1);


tmp_dir = '/data/rangan/dir_bcc/dir_jamison/dir_loader_cluster/dir_E_li30p_cluster/dir_tmp_dexcluster_nonbinary_trace_ZRmax_g010_p150_nml10/dir_tmp_dexcluster_nonbinary_trace_ZRmax_g010_p150_nml10_dexnb_trace_ZRmax_g010/dir_trace';
n_shuffle = 64;
for nshuffle=0:n_shuffle;
tmp_trace_{1+nshuffle} = textread(sprintf('%s/out_trace_s%.4d.txt',tmp_dir,nshuffle));
tmp_trace_F_{1+nshuffle} = textread(sprintf('%s/out_trace_F_s%.4d.txt',tmp_dir,nshuffle));
end;%for nshuffle=0:n_shuffle;

n_iteration = size(tmp_trace_{1},1);
n_iter_ = tmp_trace_{1}(:,1);
r_rem_ = tmp_trace_{1}(:,2);
c_rem_ = tmp_trace_{1}(:,3);
%%%%%%%%;
tmp_QR_ = zeros(n_iteration,1+n_shuffle);
tmp_QC_ = zeros(n_iteration,1+n_shuffle);
for nshuffle=0:n_shuffle;
tmp_QR_(:,1+nshuffle) = tmp_trace_{1+nshuffle}(:,4);
tmp_QC_(:,1+nshuffle) = tmp_trace_{1+nshuffle}(:,5);
end;%for nshuffle=0:n_shuffle;
%%%%%%%%;
tmp_QR_avg_ = mean(tmp_QR_(:,2:end),2);
tmp_QR_std_ = std(tmp_QR_(:,2:end),1,2);
tmp_ZR_ = zeros(n_iteration,1+n_shuffle);
for nshuffle=0:n_shuffle;
tmp_ZR_(:,1+nshuffle) = (tmp_QR_(:,1+nshuffle) - tmp_QR_avg_)./max(1e-12,tmp_QR_std_);
end;%for nshuffle=0:n_shuffle;
%%%%%%%%;
tmp_QC_avg_ = mean(tmp_QC_(:,2:end),2);
tmp_QC_std_ = std(tmp_QC_(:,2:end),1,2);
tmp_ZC_ = zeros(n_iteration,1+n_shuffle);
for nshuffle=0:n_shuffle;
tmp_ZC_(:,1+nshuffle) = (tmp_QC_(:,1+nshuffle) - tmp_QC_avg_)./max(1e-12,tmp_QC_std_);
end;%for nshuffle=0:n_shuffle;
%%%%%%%%;
n_iteration_F = size(tmp_trace_F_{1},1);
n_iter_F_ = tmp_trace_F_{1}(:,1);
r_rem_F_ = tmp_trace_F_{1}(:,2);
c_rem_F_ = tmp_trace_F_{1}(:,3);
%%%%%%%%;
tmp_QR_F_ = zeros(n_iteration_F,1+n_shuffle);
tmp_QC_F_ = zeros(n_iteration_F,1+n_shuffle);
for nshuffle=0:n_shuffle;
tmp_QR_F_(:,1+nshuffle) = tmp_trace_F_{1+nshuffle}(:,4);
tmp_QC_F_(:,1+nshuffle) = tmp_trace_F_{1+nshuffle}(:,5);
end;%for nshuffle=0:n_shuffle;
%%%%%%%%;
tmp_QR_F_avg_ = mean(tmp_QR_F_(:,2:end),2);
tmp_QR_F_std_ = std(tmp_QR_F_(:,2:end),1,2);
tmp_ZR_F_ = zeros(n_iteration_F,1+n_shuffle);
for nshuffle=0:n_shuffle;
tmp_ZR_F_(:,1+nshuffle) = (tmp_QR_F_(:,1+nshuffle) - tmp_QR_F_avg_)./max(1e-12,tmp_QR_F_std_);
end;%for nshuffle=0:n_shuffle;
%%%%%%%%;
tmp_QC_F_avg_ = mean(tmp_QC_F_(:,2:end),2);
tmp_QC_F_std_ = std(tmp_QC_F_(:,2:end),1,2);
tmp_ZC_F_ = zeros(n_iteration_F,1+n_shuffle);
for nshuffle=0:n_shuffle;
tmp_ZC_F_(:,1+nshuffle) = (tmp_QC_F_(:,1+nshuffle) - tmp_QC_F_avg_)./max(1e-12,tmp_QC_F_std_);
end;%for nshuffle=0:n_shuffle;
%%%%%%%%;

subplot(2,2,1);
hold on;
plot(n_iter_,tmp_ZR_(:,2:end),'b-','LineWidth',1);
plot(n_iter_,tmp_ZR_(:,1),'r-','LineWidth',2);
subplot(2,2,2);
hold on;
plot(n_iter_,tmp_ZC_(:,2:end),'b-','LineWidth',1);
plot(n_iter_,tmp_ZC_(:,1),'r-','LineWidth',2);
subplot(2,2,3);
hold on;
plot(n_iter_F_,tmp_ZR_F_(:,2:end),'b-','LineWidth',1);
plot(n_iter_F_,tmp_ZR_F_(:,1),'r-','LineWidth',2);
subplot(2,2,4);
hold on;
plot(n_iter_F_,tmp_ZC_F_(:,2:end),'b-','LineWidth',1);
plot(n_iter_F_,tmp_ZC_F_(:,1),'r-','LineWidth',2);

tmp_lp_R_F_ = z_to_p_0(tmp_ZR_F_);
tmp_lp_r_rem_F_ = lnchoosek(max(r_rem_F_),r_rem_F_);
% test out probability of particular cutpoint: ;
plot(n_iter_F_,tmp_lp_r_rem_F_(:) - tmp_lp_R_F_(:,1));


