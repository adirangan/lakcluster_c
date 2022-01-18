function test_dexcluster_multi_csv_2(M,N,snr,n_cluster,gamma,n_rank,n_iteration,n_shuffle,flag_rerun);
% saves data used for test_dexcluster_multi_2. ;
%{

  M = 178; N = 2e3; n_cluster = 1; gamma = 0.00; n_rank = 6; n_iteration = 1; n_shuffle = 256; flag_rerun = 0;
  n_step = 21;
  snr_ = linspace(0.5,1.5,n_step);
  for nstep=1:n_step;
  snr = snr_(nstep);
  test_dexcluster_multi_csv_2(M,N,snr,n_cluster,gamma,n_rank,n_iteration,n_shuffle,flag_rerun);
  end;%for nstep=1:n_step;

  %}

setup();
n_vrank = 2;
dir_trunk = sprintf('/data/rangan/dir_bcc/dir_jamison');
dir_pvclust = sprintf('/data/rangan/dir_bcc/dir_pvclust');
dir_base = sprintf('%s/dir_dexcluster_multi_csv_2',dir_pvclust);
if (~exist(dir_base,'dir')); disp(sprintf(' %% mkdir %s',dir_base)); mkdir(dir_base); end;
%%%%%%%%;
for niteration=1:n_iteration;
disp(sprintf(' %% niteration %d/%d',niteration,n_iteration));
%%%%%%%%;
rng(niteration);
[A_n_,label_A_,n_label_A_,pf_,pi_,snr_] = random_matrix_planted_cluster_0(M,N,snr,n_cluster);
str_xfix = test_dexcluster_multi_xfix_2('base',M,N,snr,n_cluster,gamma,niteration);
str_base = sprintf('%s/%s.csv',dir_base,str_xfix);
if (~exist(str_base,'file')); 
disp(sprintf(' %% %s not found, creating',str_base)); 
writematrix(A_n_(pi_{1},pi_{2}),str_base);
end;%if (~exist(str_base,'file')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for niteration=1:n_iteration;

