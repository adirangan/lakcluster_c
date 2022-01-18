clear;
setup();

rng(1);
N=196; X = ceil(N^0.75);
A = randn(N,N+1); B=randn(X,X+1)+2.5*ones(X,1)*(2*(rand(1,X+1)>0.5) - 1); 
%A = randn(N,N+1); B=randn(X,X+1)+1.5*ones(X,1)*ones(1,X+1);
for nj=1:3; 
%pf_{nj} = randperm(size(A,nj)); 
pf_{nj} = 1:size(A,nj);
[~,pi_{nj}] = sort(pf_{nj}); [~,pi_{nj}] = sort(pi_{nj}); 
end;% for nj=1:3;
A(pf_{1}(1:X),pf_{2}(1:X+1)) = B;
%A = -A;

dir_code = '/data/rangan/dir_bcc/dir_lakcluster_c_dev';
dir_trunk = pwd;
prefix = 'test';
M_n_ = {A};
rev_flag = 0;
A_n_rind_ = {[1:size(A,1)]};
A_n_cind = [1:size(A,2)];
Z_n_rind_ = {[]};
T_n_ = {ones(size(A,1),1)};
T_n_cind = [1];
GLOBAL_TEST_sparse = 0;
GLOBAL_kappa_squared = 0;
gamma = 0.0;
B_MLT = 34;
Ireq = 0;
shuffle_num = 0;
verbose_flag = 0;
force_create_flag = 1;
slurm_walltime = 0;
slurm_nnodes = 1;
slurm_tpn = 1;
slurm_memdecl = 0;

[tmp_U_,tmp_S_,tmp_V_] = svds(A,1); [~,tmp_U_ij_] = sort(abs(tmp_U_),'descend'); [~,tmp_V_ij_] = sort(abs(tmp_V_),'descend');
rkeep_s_ = tmp_U_ij_; ckeep_s_ = tmp_V_ij_;
[~,tmp_Bn_] = intersect(rkeep_s_,pi_{1}(1:size(B,1)),'stable'); [~,tmp_Bc_] = intersect(rkeep_s_,pi_{1}(size(B,1)+1:end),'stable'); Auc_s_(1) = auc_0(tmp_Bn_,tmp_Bc_);
[~,tmp_Bn_] = intersect(ckeep_s_,pi_{2}(1:size(B,2)),'stable'); [~,tmp_Bc_] = intersect(ckeep_s_,pi_{2}(size(B,2)+1:end),'stable'); Auc_s_(2) = auc_0(tmp_Bn_,tmp_Bc_);
Auc_s_avg = mean(Auc_s_); disp(sprintf(' %% spectral Auc [%0.2f %0.2f] --> %0.2f',Auc_s_,Auc_s_avg));

dexcluster_uADZSZDA_dr_2(dir_code,dir_trunk,prefix,M_n_,rev_flag,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind,GLOBAL_TEST_sparse,GLOBAL_kappa_squared,gamma,B_MLT,Ireq,shuffle_num,verbose_flag,force_create_flag,slurm_walltime,slurm_nnodes,slurm_tpn,slurm_memdecl);
%dexcluster_uADZSZDA_dr_0(dir_code,dir_trunk,prefix,M_n_,rev_flag,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind,gamma,B_MLT,Ireq,shuffle_num,0,verbose_flag,force_create_flag);
out_xdrop_a_ = textread('/data/rangan/dir_bcc/dir_jamison/dir_test/dir_test_dex_uADZSZDA_D_g000_B34_I0_s0000/out_xdrop_a.txt');
rdrop_a_ = out_xdrop_a_(:,1); rdrop_a_ = 1+rdrop_a_(find(rdrop_a_>-1)); rkeep_a_ = rdrop_a_(end:-1:1); 
cdrop_a_ = out_xdrop_a_(:,2); cdrop_a_ = 1+cdrop_a_(find(cdrop_a_>-1)); ckeep_a_ = cdrop_a_(end:-1:1);
[~,tmp_Bn_] = intersect(rkeep_a_,pi_{1}(1:size(B,1)),'stable'); [~,tmp_Bc_] = intersect(rkeep_a_,pi_{1}(size(B,1)+1:end),'stable'); Auc_a_(1) = auc_0(tmp_Bn_,tmp_Bc_);
[~,tmp_Bn_] = intersect(ckeep_a_,pi_{2}(1:size(B,2)),'stable'); [~,tmp_Bc_] = intersect(ckeep_a_,pi_{2}(size(B,2)+1:end),'stable'); Auc_a_(2) = auc_0(tmp_Bn_,tmp_Bc_);
Auc_a_avg = mean(Auc_a_); disp(sprintf(' %% dexcluster: Auc [%0.2f %0.2f] --> %0.2f',Auc_a_,Auc_a_avg));

A_n_ = tutorial_binary_uncompress('./dir_test/test_M_0_n.b16'); A_t_ = transpose(A_n_); e_n_= ones(size(A_n_,1),1);
out_xdrop_m_ = test_dcc_AAAA(A_n_);
rdrop_m_ = out_xdrop_m_(:,1); rdrop_m_ = 1+rdrop_m_(find(rdrop_m_>-1)); 
cdrop_m_ = out_xdrop_m_(:,2); cdrop_m_ = 1+cdrop_m_(find(cdrop_m_>-1)); 
rkeep_m_ = rdrop_m_(end:-1:1); ckeep_m_ = cdrop_m_(end:-1:1);
[~,tmp_Bn_] = intersect(rkeep_m_,pi_{1}(1:size(B,1)),'stable'); [~,tmp_Bc_] = intersect(rkeep_m_,pi_{1}(size(B,1)+1:end),'stable'); Auc_m_(1) = auc_0(tmp_Bn_,tmp_Bc_);
[~,tmp_Bn_] = intersect(ckeep_m_,pi_{2}(1:size(B,2)),'stable'); [~,tmp_Bc_] = intersect(ckeep_m_,pi_{2}(size(B,2)+1:end),'stable'); Auc_m_(2) = auc_0(tmp_Bn_,tmp_Bc_);
Auc_m_avg = mean(Auc_m_); disp(sprintf(' %% test_dcc Auc [%0.2f %0.2f] --> %0.2f',Auc_m_,Auc_m_avg));

[xdrop_condensed_,Z_avg_,Z_avg_condensed_,xdrop_,xrmv] = tensor_auc_ver3(A);
rdrop_t_ = xdrop_condensed_{1}; cdrop_t_ = xdrop_condensed_{2};
rkeep_t_ = rdrop_t_(end:-1:1); ckeep_t_ = cdrop_t_(end:-1:1);
[~,tmp_Bn_] = intersect(rkeep_t_,pi_{1}(1:size(B,1)),'stable'); [~,tmp_Bc_] = intersect(rkeep_t_,pi_{1}(size(B,1)+1:end),'stable'); Auc_t_(1) = auc_0(tmp_Bn_,tmp_Bc_);
[~,tmp_Bn_] = intersect(ckeep_t_,pi_{2}(1:size(B,2)),'stable'); [~,tmp_Bc_] = intersect(ckeep_t_,pi_{2}(size(B,2)+1:end),'stable'); Auc_t_(2) = auc_0(tmp_Bn_,tmp_Bc_);
Auc_t_avg = mean(Auc_t_); disp(sprintf(' %% tensor_auc Auc [%0.2f %0.2f] --> %0.2f',Auc_t_,Auc_t_avg));

colormap(colormap_beach());
subplot(2,3,1); imagesc(A(pi_{1},pi_{2}),[-1,1]); set(gca,'XTick',[],'Ytick',[]); title('A orig');
subplot(2,3,2); imagesc(A_n_(pi_{1},pi_{2}),[-1,1]); set(gca,'XTick',[],'Ytick',[]); title('A binarized');
subplot(2,3,3); imagesc(A(rkeep_a_,ckeep_a_),[-1,1]); set(gca,'XTick',[],'Ytick',[]); title(sprintf('dexcluster A%0.2f',Auc_a_avg));
subplot(2,3,4); imagesc(A_n_(rkeep_m_,ckeep_m_),[-1,1]); set(gca,'XTick',[],'Ytick',[]); title(sprintf('matlab dcc A%0.2f',Auc_m_avg));
subplot(2,3,5); imagesc(A(rkeep_t_,ckeep_t_),[-1,1]); set(gca,'XTick',[],'Ytick',[]); title(sprintf('tensor auc A%0.2f',Auc_t_avg));
subplot(2,3,6); imagesc(A(rkeep_s_,ckeep_s_),[-1,1]); set(gca,'XTick',[],'Ytick',[]); title(sprintf('spectral A%0.2f',Auc_s_avg));
set(gcf,'Position',1+[0,0,1024*2,1024]);
