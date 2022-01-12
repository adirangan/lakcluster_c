function ....
[ ...
 mr_A_prm_ ...
,mr_Z_prm_ ...
] = ...
xxxcluster_uADZSZDA_shuffle_2( ...
 shuffle_num ...
,M_n_ ...
,A_n_rij_ ...
,A_n_cij ...
,Z_n_rij_ ...
,T_n_ ...
,T_n_cij ...
);
% creates shuffled row masks respecting categorical-covariates and continuous-covariate sectors ;
% test with: ;
%{
  xxxcluster_uADZSZDA_shuffle_2();
  %}

if (nargin<1);
clear M_n_ A_n_rij_ Z_n_rij_ T_n_  ;
nbins = 3;
M_n_ = cell(3,1);
A_n_rij_ = cell(3,1);
Z_n_rij_ = cell(3,1);
T_n_ = cell(3,1);
A_n_cols = 2;
A_n_cij = [1];
T_n_cols = 20; T_n_cij = [1,3,7];
%T_n_cols = 1; T_n_cij = [1];
n_sector = 2^(length(T_n_cij)-1);
for (nb1=0:nbins-1);
M_n_{1+nb1} = zeros(128*(1+nb1),A_n_cols);
T_n_{1+nb1} = 2*(randn(128*(1+nb1),T_n_cols)>0)-1;
T_n_{1+nb1}(:,1) = 1;
prm = randperm(size(M_n_{1+nb1},1));
lA = round(size(M_n_{1+nb1},1)/3);
lZ = round(size(M_n_{1+nb1},1)/2);
A_n_rij_{1+nb1} = prm(1:lA);
Z_n_rij_{1+nb1} = prm(lA + (1:lZ));
M_n_{1+nb1}(A_n_rij_{1+nb1},1) = 2*nb1 + 0;
T_tmp = T_n_{1+nb1}(A_n_rij_{1+nb1},T_n_cij(2:end)); nsec = 2^(length(T_n_cij(2:end)));
A_sector_{1+nb1} = (T_tmp>0)*transpose(2.^[0:length(T_n_cij(2:end))-1]);
M_n_{1+nb1}(A_n_rij_{1+nb1},2) = A_sector_{1+nb1};
M_n_{1+nb1}(Z_n_rij_{1+nb1},1) = 2*nb1 + 1;
T_tmp = T_n_{1+nb1}(Z_n_rij_{1+nb1},T_n_cij(2:end)); nsec = 2^(length(T_n_cij(2:end)));
Z_sector_{1+nb1} = (T_tmp>0)*transpose(2.^[0:length(T_n_cij(2:end))-1]);
M_n_{1+nb1}(Z_n_rij_{1+nb1},2) = Z_sector_{1+nb1};
end;%for (nb1=0:nbins-1);
shuffle_num = 1;
[mr_A_prm_,mr_Z_prm_] = xxxcluster_uADZSZDA_shuffle_2(shuffle_num,M_n_,A_n_rij_,A_n_cij,Z_n_rij_,T_n_,T_n_cij);
prows = 4;
for nb1=0:nbins-1;
lM = size(M_n_{1+nb1},1);
lA = length(A_n_rij_{1+nb1});
lZ = length(Z_n_rij_{1+nb1});
A_p_rij_{1+nb1} = find(mr_A_prm_{1+nb1}); Z_p_rij_{1+nb1} = find(mr_Z_prm_{1+nb1});
%b_ = 2*nb1 + (0:1);
b_ = 0:2*nbins-1;
s_ = 0:n_sector-1; if (n_sector==1); s_ = [-1:+1]; end;
h_tmp = [hist(M_n_{1+nb1}(A_n_rij_{1+nb1},1),b_) ; hist(M_n_{1+nb1}(A_p_rij_{1+nb1},1),b_)];
subplot(nbins,prows,1 + nb1*prows); bar(b_,transpose(h_tmp)); title(sprintf('A-label 1+nb1 %d',nb1));
h_tmp = [hist(M_n_{1+nb1}(A_n_rij_{1+nb1},2),s_) ; hist(M_n_{1+nb1}(A_p_rij_{1+nb1},2),s_)];
subplot(nbins,prows,2 + nb1*prows); bar(s_,transpose(h_tmp)); title(sprintf('A-sector 1+nb1 %d',nb1));
h_tmp = [hist(M_n_{1+nb1}(Z_n_rij_{1+nb1},1),b_) ; hist(M_n_{1+nb1}(Z_p_rij_{1+nb1},1),b_)];
subplot(nbins,prows,3 + nb1*prows); bar(b_,transpose(h_tmp)); title(sprintf('Z-label 1+nb1 %d',nb1));
h_tmp = [hist(M_n_{1+nb1}(Z_n_rij_{1+nb1},2),s_) ; hist(M_n_{1+nb1}(Z_p_rij_{1+nb1},2),s_)];
subplot(nbins,prows,4 + nb1*prows); bar(s_,transpose(h_tmp)); title(sprintf('Z-sector 1+nb1 %d',nb1));
end;%for nb1=0:nbins-1;
disp(sprintf('returning')); return;
end;%if (nargin<1);

nbins = length(M_n_);
[A_n_cols,Y_n_cols,T_n_cols,~,~] = xxxcluster_uADZSZDA_check_1(shuffle_num,M_n_,A_n_rij_,A_n_cij,Z_n_rij_,T_n_,T_n_cij);

for nb1=0:nbins-1;
mr_M = zeros(size(M_n_{1+nb1},1),1);
mr_A_ori_{1+nb1} = mr_M; mr_A_ori_{1+nb1}(A_n_rij_{1+nb1})=1;
mr_Z_ori_{1+nb1} = mr_M; mr_Z_ori_{1+nb1}(Z_n_rij_{1+nb1})=1;
mr_A_prm_{1+nb1} = mr_A_ori_{1+nb1};
mr_Z_prm_{1+nb1} = mr_Z_ori_{1+nb1};
end;%for nb1=0:nbins-1;

if (shuffle_num>0); % performing covariate-respecting shuffle ;
rng(shuffle_num); 
for nb1=0:nbins-1;
mr_M_tmp = zeros(size(M_n_{1+nb1},1),1); mr_M_tmp(A_n_rij_{1+nb1})=1; mr_M_tmp(Z_n_rij_{1+nb1})=1;
disp(sprintf(' %% nb1 %d, sum(mr_M_tmp) %d',1+nb1,sum(mr_M_tmp)));
mr_M__on = find(mr_M_tmp);
T_tmp = T_n_{1+nb1}(mr_M__on,T_n_cij(2:end)); n_sector = 2^(length(T_n_cij)-1);
M_sector_{1+nb1} = (T_tmp>0)*transpose(2.^[0:length(T_n_cij)-2]);
for nsector=0:n_sector-1; 
sec_rij_{1+nsector} = find(M_sector_{1+nb1}==nsector);
mr_X_tmp_{1+nsector} = zeros(size(M_n_{1+nb1},1),1); mr_X_tmp_{1+nsector}(mr_M__on(sec_rij_{1+nsector}))=1;
end;%for nsector=0:n_sector-1; 
for nsector=0:n_sector-1;
lA_(1+nsector) = length(intersect(A_n_rij_{1+nb1},find(mr_X_tmp_{1+nsector})));
lZ_(1+nsector) = length(intersect(Z_n_rij_{1+nb1},find(mr_X_tmp_{1+nsector})));
prm_{1+nsector} = randperm(length(sec_rij_{1+nsector}));
prm_A_{1+nsector} = sec_rij_{1+nsector}(prm_{1+nsector}(1:lA_(1+nsector))); 
prm_Z_{1+nsector} = sec_rij_{1+nsector}(prm_{1+nsector}(lA_(1+nsector) + (1:lZ_(1+nsector))));
end;%for nsector=0:n_sector-1;
mr_A_prm_{1+nb1} = zeros(size(M_n_{1+nb1},1),1); mr_Z_prm_{1+nb1} = zeros(size(M_n_{1+nb1},1),1);
for nsector=0:n_sector-1;
mr_A_prm_{1+nb1}(mr_M__on(prm_A_{1+nsector}))=1; mr_Z_prm_{1+nb1}(mr_M__on(prm_Z_{1+nsector}))=1;
end;%for nsector=0:n_sector-1;
disp_flag=1;
if disp_flag;
tmp_A_prm = find(mr_A_prm_{1+nb1}); tmp_Z_prm = find(mr_Z_prm_{1+nb1});
for nsector=0:n_sector-1; 
lA_prm_(1+nsector) = length(intersect(tmp_A_prm,find(mr_X_tmp_{1+nsector}))); 
lZ_prm_(1+nsector) = length(intersect(tmp_Z_prm,find(mr_X_tmp_{1+nsector}))); 
lX_prm_(1+nsector) = length(intersect(tmp_A_prm,tmp_Z_prm));
end;%for nsector=0:n_sector-1;
disp(sprintf(' %% shuffle: nb1 %d',nb1));
disp([lA_ ; lA_prm_ ; lZ_ ; lZ_prm_ ; lX_prm_]); 
disp(sprintf(' %% mr_A_prm %d/%d; mr_Z_prm %d/%d',sum(mr_A_prm_{1+nb1}),length(A_n_rij_{1+nb1}),sum(mr_Z_prm_{1+nb1}),length(Z_n_rij_{1+nb1})));
end;% if disp_flag;
end;%for nb1=0:nbins-1;
end;% if (shuffle_num>0); % performing covariate-respecting shuffle ;
