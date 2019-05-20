function dcgxpander_test_uAZSZA_ver0(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num)
% trying to test dcgxpander with single bins and mrnd=0 ; 
% we only embed one bicluster: ;
% X0 is balanced across covariate categories and T, and is exclusive to the cases ;
% try: ;
%{
  
  dir_trunk=pwd;N=1024*1;X_factor=0.575;X_esm=0.0;gamma=0.25;B_MLT=32;rng_num=1;
  dcgxpander_test_uAZSZA_ver0(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

  %}


if (nargin<8); pt_num = -1; end;
if (nargin<7); rng_num = 1; end;
if (nargin<4); B_MLT = 32; end;
if (nargin<5); gamma = 0.025; end;
if (nargin<4); X_esm = 0.1; end;
if (nargin<3); X_factor = 0.65; end;
if (nargin<2); N = 1024; end;
if (nargin<1); dir_trunk = pwd; end;

nbins = 1;
nrows_A_ = ceil(N/nbins*ones(1,nbins)) ;
nrows_Z_ = 0;%ceil(2*N/nbins*ones(1,nbins)) ;
ncols_A = N+1 ;
nmds = 0;
ncols_T = 1+nmds; 
if nmds>0;
for nb1=0:nbins-1; 
MDT_{1+nb1} = randn(nrows_A_(1+nb1),nmds); 
MDT_{1+nb1}(1:floor(nrows_A_(1+nb1)/2),1) = -abs(MDT_{1+nb1}(1:floor(nrows_A_(1+nb1)/2),1));
MDT_{1+nb1}(floor(nrows_A_(1+nb1)/2):end,1) = +abs(MDT_{1+nb1}(floor(nrows_A_(1+nb1)/2):end,1));
MDS_{1+nb1} = randn(nrows_Z_(1+nb1),nmds); 
MDS_{1+nb1}(1:floor(nrows_Z_(1+nb1)/2),1) = -abs(MDS_{1+nb1}(1:floor(nrows_Z_(1+nb1)/2),1));
MDS_{1+nb1}(floor(nrows_Z_(1+nb1)/2):end,1) = +abs(MDS_{1+nb1}(floor(nrows_Z_(1+nb1)/2):end,1));
end;%for nb1=0:nbins-1;
end;%if nmds>0;
if nmds==0;
for nb1=0:nbins-1; 
MDT_{1+nb1} = zeros(nrows_A_(1+nb1),nmds); 
MDS_{1+nb1} = zeros(nrows_Z_(1+nb1),nmds); 
end;%for nb1=0:nbins-1; 
end;%if nmds==0;
mrnd = 0.0;
verbose_flag=0;
%Ireq = 0; % 112716 ;
%Ireq = 1; % 112916 ;
Ireq = 0; % 010817 ;

nrows_X = ceil(N.^(X_factor)); nrows_XX = ceil(1.5*nrows_X);
ncols_X = ceil((N+1).^(X_factor)); ncols_XX = ceil(1.5*ncols_X);
X_eps = X_esm/sqrt(nrows_X); % spectral error for X;
X_pp = tutorial_g1em(X_esm); disp(sprintf(' %% X_esm %0.4f, log10(X_esm) %0.4f; X_pp %0.4f, log10(X_pp) %0.4f',X_esm,log10(X_esm),X_pp,log10(X_pp)));

test_string = 'dcgxpander_test_uAZSZA';
dir_base = sprintf('%s/dir_%s',dir_trunk,test_string); if ~exist(dir_base,'dir'); mkdir(dir_base); end;
[posfix_use,X_factor_d,X_esm_d,gamma_d] = xfix_gen(N,X_factor,X_esm,gamma,B_MLT,rng_num);
in__name = sprintf('%s_%s__in',test_string,posfix_use);
out_name = sprintf('%s_%s_out',test_string,posfix_use);
dir__in = sprintf('%s/dir_%s',dir_base,in__name); if ~exist(dir__in,'dir'); mkdir(dir__in); end;
dir_out = sprintf('%s/dir_%s',dir_base,out_name); if ~exist(dir_out,'dir'); mkdir(dir_out); end; 
%path(path,dir__in); path(path,dir_out);
rng(rng_num);

clear An At Zn Zt Tn Tt Sn St mr_A_ mr_Z_ mc_A mc_T bitj;
bitj = 16;

mc_A = (rand(1,ncols_A)>mrnd);
mc_T = (rand(1,ncols_T)>mrnd); mc_T(1)=1; % ensure first column of [T;S] is all ones ;
for nb1=0:nbins-1; 
mr_A_{1+nb1} = (rand(nrows_A_(1+nb1),1)>mrnd);
mr_Z_{1+nb1} = (rand(nrows_Z_(1+nb1),1)>mrnd);
end;%for nb1=0:nbins-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate random matrices ;
for nb1=0:nbins-1; 
An{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_A)<0.5)-1;
Zn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_A)<0.5)-1;
Tn{1+nb1} = 2*([ones(nrows_A_(1+nb1),1) , MDT_{1+nb1}]>0)-1;
Sn{1+nb1} = 2*([ones(nrows_Z_(1+nb1),1) , MDS_{1+nb1}]>0)-1;
end;%for nb1=0:nbins-1; 
% generate clusters ;
% XW is a simple nonspecific cluster ;
XW = 2*(rand(nrows_X,ncols_X)<X_pp)-1;
for nb1=0:nbins-1; XW_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,XW_pr_i_{1+nb1}] = sort(XW_pr_{1+nb1}); [tmp,XW_pr_i_{1+nb1}] = sort(XW_pr_i_{1+nb1}); end;
XW_pc = randperm(ncols_A);
tmp_c = [XW_pc]; [tmp,XW_pc_i] = sort(tmp_c); [tmp,XW_pc_i] = sort(XW_pc_i);

%u_mds = randn(nmds,1); u_mds = u_mds/norm(u_mds);
u_mds = 1;
for nb1=0:nbins-1;
T_u_{1+nb1} = MDT_{1+nb1}*u_mds;
S_u_{1+nb1} = MDS_{1+nb1}*u_mds;
T_u_pos_{1+nb1} = find(T_u_{1+nb1}>0);
T_u_neg_{1+nb1} = find(T_u_{1+nb1}<=0);
end;%for nb1=0:nbins-1;
% generate clusters ;

% X0 is balanced across covariate categories and T, and is exclusive to the cases ;
X0 = 2*(rand(nrows_X,ncols_X)<X_pp)-1; X0(:,floor(end/2):end) = -X0(:,floor(end/2):end);
for nb1=0:nbins-1; X0_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X0_pr_i_{1+nb1}] = sort(X0_pr_{1+nb1}); [tmp,X0_pr_i_{1+nb1}] = sort(X0_pr_i_{1+nb1}); end;
X0_pc = randperm(ncols_A); 
% X4 is balanced across covariate categories and T, but extends to include cases and controls ;
X4 = 2*(rand(3*nrows_XX,ncols_XX)<X_pp)-1; X4(:,floor(end/2):end) = -X4(:,floor(end/2):end);
for nb1=0:nbins-1; 
X4_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X4_pr_i_{1+nb1}] = sort(X4_pr_{1+nb1}); [tmp,X4_pr_i_{1+nb1}] = sort(X4_pr_i_{1+nb1}); 
X4_pz_{1+nb1} = randperm(nrows_Z_(1+nb1)); [tmp,X4_pz_i_{1+nb1}] = sort(X4_pz_{1+nb1}); [tmp,X4_pz_i_{1+nb1}] = sort(X4_pz_i_{1+nb1}); 
end;%for nb1=0:nbins-1; 
X4_pc = randperm(ncols_A); 
% X1 is balanced across continuous-covariates and is exclusive to the cases, but is restricted to category-2;
X1 = 2*(rand(nrows_XX,ncols_XX)<X_pp)-1; X1(:,floor(end/2):end) = -X1(:,floor(end/2):end);
for nb1=0:nbins-1; X1_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X1_pr_i_{1+nb1}] = sort(X1_pr_{1+nb1}); [tmp,X1_pr_i_{1+nb1}] = sort(X1_pr_i_{1+nb1}); end;
X1_pc = randperm(ncols_A); 
% X2 is balanced across continuous-covariates and is exclusive to the cases, but is restricted to category-0;
X2 = 2*(rand(nrows_XX,ncols_XX)<X_pp)-1; X2(:,floor(end/2):end) = -X2(:,floor(end/2):end);
for nb1=0:nbins-1; X2_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X2_pr_i_{1+nb1}] = sort(X2_pr_{1+nb1}); [tmp,X2_pr_i_{1+nb1}] = sort(X2_pr_i_{1+nb1}); end;
X2_pc = randperm(ncols_A); 
% X3 is balanced across continuous-covariates and is exclusive to the cases, but is restricted to category-1;
X3 = 2*(rand(nrows_XX,ncols_XX)<X_pp)-1; X3(:,floor(end/2):end) = -X3(:,floor(end/2):end);
for nb1=0:nbins-1; X3_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X3_pr_i_{1+nb1}] = sort(X3_pr_{1+nb1}); [tmp,X3_pr_i_{1+nb1}] = sort(X3_pr_i_{1+nb1}); end;
X3_pc = randperm(ncols_A); 
% X5 is balanced across covariate categories and is exclusive to the cases, but is restricted to positive T_u ;
X5 = 2*(rand(nrows_XX,ncols_XX)<X_pp)-1;%tutorial_binary_makelr2(nrows_XX,ncols_XX,X_pp);
X5(:,floor(end/2):end) = -X5(:,floor(end/2):end);
X5_pc = randperm(ncols_A); 
for nb1=0:nbins-1; 
X5_pr_{1+nb1} = [ T_u_pos_{1+nb1}(randperm(length(T_u_pos_{1+nb1}))) ; T_u_neg_{1+nb1}(randperm(length(T_u_neg_{1+nb1}))) ];
[tmp,X5_pr_i_{1+nb1}] = sort(X5_pr_{1+nb1}); [tmp,X5_pr_i_{1+nb1}] = sort(X5_pr_i_{1+nb1});
end;%for nb1=0:nbins-1; 
% X6 is balanced across covariate categories and is exclusive to the cases, but is restricted to negative T_u ;
X6 = 2*(rand(nrows_XX,ncols_XX)<X_pp)-1;%tutorial_binary_makelr2(nrows_XX,ncols_XX,X_pp);
X6(:,floor(end/2):end) = -X6(:,floor(end/2):end);
X6_pc = randperm(ncols_A); 
for nb1=0:nbins-1; 
X6_pr_{1+nb1} = [ T_u_neg_{1+nb1}(randperm(length(T_u_neg_{1+nb1}))) ; T_u_pos_{1+nb1}(randperm(length(T_u_pos_{1+nb1}))) ];
[tmp,X6_pr_i_{1+nb1}] = sort(X6_pr_{1+nb1}); [tmp,X6_pr_i_{1+nb1}] = sort(X6_pr_i_{1+nb1});
end;%for nb1=0:nbins-1; 

[tmp,X0_pc_i] = sort(X0_pc); [tmp,X0_pc_i] = sort(X0_pc_i);
[tmp,X1_pc_i] = sort(X1_pc); [tmp,X1_pc_i] = sort(X1_pc_i);
[tmp,X2_pc_i] = sort(X2_pc); [tmp,X2_pc_i] = sort(X2_pc_i);
[tmp,X3_pc_i] = sort(X3_pc); [tmp,X3_pc_i] = sort(X3_pc_i);
[tmp,X4_pc_i] = sort(X4_pc); [tmp,X4_pc_i] = sort(X4_pc_i);
[tmp,X5_pc_i] = sort(X5_pc); [tmp,X5_pc_i] = sort(X5_pc_i);
[tmp,X6_pc_i] = sort(X6_pc); [tmp,X6_pc_i] = sort(X6_pc_i);

%insert clusters;
% do not insert XW ;
% for (nb1=0:nbins-1); An{1+nb1}(XW_pr_{1+nb1}(0*floor(nrows_X/nbins) + (1:floor(nrows_X/nbins))),XW_pc(1:ncols_X)) = XW(nb1*floor(nrows_X/nbins) + (1:floor(nrows_X/nbins)),1:ncols_X); end;

% insert X0;
for nb1=0:nbins-1; An{1+nb1}(X0_pr_{1+nb1}(0*floor(nrows_X/nbins) + (1:floor(nrows_X/nbins))),X0_pc(1:ncols_X)) = X0(nb1*floor(nrows_X/nbins) + (1:floor(nrows_X/nbins)),1:ncols_X); end;
%{
% insert X4;
for nb1=0:nbins-1; 
An{1+nb1}(X4_pr_{1+nb1}(0*floor(nrows_XX/nbins) + (1:floor(nrows_XX/nbins))),X4_pc(1:ncols_XX)) = X4(nb1*floor(3*nrows_XX/nbins) + (1:floor(nrows_XX/nbins)),1:ncols_XX); 
Zn{1+nb1}(X4_pz_{1+nb1}(0*floor(nrows_XX/nbins) + (1:floor(2*nrows_XX/nbins))),X4_pc(1:ncols_XX)) = X4(nb1*floor(3*nrows_XX/nbins) + floor(nrows_XX/nbins) + (1:floor(2*nrows_XX/nbins)),1:ncols_XX); 
end;%for nb1=0:nbins-1; 
% insert X1;
nb1=2-1; An{1+nb1}(X1_pr_{1+nb1}(0*floor(nrows_XX/1) + (1:floor(nrows_XX/2/1))),X1_pc(1:ncols_XX)) = X1(0*floor(nrows_XX/2/1) + (1:floor(nrows_XX/2/1)),1:ncols_XX);
nb1=3-1; An{1+nb1}(X1_pr_{1+nb1}(0*floor(nrows_XX/1) + (1:floor(nrows_XX/2/1))),X1_pc(1:ncols_XX)) = X1(1*floor(nrows_XX/2/1) + (1:floor(nrows_XX/2/1)),1:ncols_XX);
% insert X2;
nb1=1-1; An{1+nb1}(X2_pr_{1+nb1}(0*floor(nrows_XX/1) + (1:floor(nrows_XX/2/1))),X2_pc(1:ncols_XX)) = X2(0*floor(nrows_XX/2/1) + (1:floor(nrows_XX/2/1)),1:ncols_XX);
nb1=3-1; An{1+nb1}(X2_pr_{1+nb1}(0*floor(nrows_XX/1) + (1:floor(nrows_XX/2/1))),X2_pc(1:ncols_XX)) = X2(1*floor(nrows_XX/2/1) + (1:floor(nrows_XX/2/1)),1:ncols_XX);
% insert X3;
nb1=1-1; An{1+nb1}(X3_pr_{1+nb1}(0*floor(nrows_XX/1) + (1:floor(nrows_XX/2/1))),X3_pc(1:ncols_XX)) = X3(0*floor(nrows_XX/2/1) + (1:floor(nrows_XX/2/1)),1:ncols_XX);
nb1=2-1; An{1+nb1}(X3_pr_{1+nb1}(0*floor(nrows_XX/1) + (1:floor(nrows_XX/2/1))),X3_pc(1:ncols_XX)) = X3(1*floor(nrows_XX/2/1) + (1:floor(nrows_XX/2/1)),1:ncols_XX);
% insert X5 ;
for nb1=0:nbins-1; An{1+nb1}(X5_pr_{1+nb1}(0*floor(nrows_XX/nbins) + (1:floor(nrows_XX/nbins))),X5_pc(1:ncols_XX)) = X5(nb1*floor(nrows_XX/nbins) + (1:floor(nrows_XX/nbins)),1:ncols_XX); end;
% insert X6 ;
for nb1=0:nbins-1; An{1+nb1}(X6_pr_{1+nb1}(0*floor(nrows_XX/nbins) + (1:floor(nrows_XX/nbins))),X6_pc(1:ncols_XX)) = X6(nb1*floor(nrows_XX/nbins) + (1:floor(nrows_XX/nbins)),1:ncols_XX); end;
 %}

% generate transposes ;
for nb1=0:nbins-1;
At{1+nb1} = transpose(An{1+nb1}); Zt{1+nb1} = transpose(Zn{1+nb1});
Tt{1+nb1} = transpose(Tn{1+nb1}); St{1+nb1} = transpose(Sn{1+nb1});
end;%for nb1=0:nbins-1;

disp_flag=0;
if disp_flag;
figure; clf;
for nb1=0:nbins-1;
% X4;
subplot(3,nbins,1+nb1+0*nbins);
tmp_c = X4_pc_i;
tmp_r = X4_pr_i_{1+nb1};
imagesc([ ...
	  Tn{1+nb1}(X4_pr_i_{1+nb1},:) , zeros(nrows_A_(1+nb1),8) , repmat(T_u_{1+nb1}(X4_pr_i_{1+nb1}),1,8) , zeros(nrows_A_(1+nb1),8) , An{1+nb1}(X4_pr_i_{1+nb1},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X4, A, nb1 %d',nb1));
subplot(3,nbins,[1+nb1+1*nbins , 1+nb1+2*nbins]);
tmp_c = X4_pc_i;
tmp_r = X4_pz_i_{1+nb1};
imagesc([ ...
	  Sn{1+nb1}(X4_pz_i_{1+nb1},:) , zeros(nrows_Z_(1+nb1),8) , repmat(S_u_{1+nb1}(X4_pz_i_{1+nb1}),1,8) , zeros(nrows_Z_(1+nb1),8) , Zn{1+nb1}(X4_pz_i_{1+nb1},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X4, Z, nb1 %d',nb1));
end;%for nb1=0:nbins-1;
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
figure; clf;
for nb1=0:nbins-1;
% X1;
subplot(3,nbins,1+nb1+0*nbins);
tmp_c = X1_pc_i;
tmp_r = X1_pr_i_{1+nb1};
imagesc([ ...
	  Tn{1+nb1}(X1_pr_i_{1+nb1},:) , zeros(nrows_A_(1+nb1),8) , repmat(T_u_{1+nb1}(X1_pr_i_{1+nb1}),1,8) , zeros(nrows_A_(1+nb1),8) , An{1+nb1}(X1_pr_i_{1+nb1},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X1, nb1 %d',nb1));
% X2;
subplot(3,nbins,1+nb1+1*nbins);
tmp_c = X2_pc_i;
tmp_r = X2_pr_i_{1+nb1};
imagesc([ ...
	  Tn{1+nb1}(X2_pr_i_{1+nb1},:) , zeros(nrows_A_(1+nb1),8) , repmat(T_u_{1+nb1}(X2_pr_i_{1+nb1}),1,8) , zeros(nrows_A_(1+nb1),8) , An{1+nb1}(X2_pr_i_{1+nb1},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X2, nb1 %d',nb1));
% X3;
subplot(3,nbins,1+nb1+2*nbins);
tmp_c = X3_pc_i;
tmp_r = X3_pr_i_{1+nb1};
imagesc([ ...
	  Tn{1+nb1}(X3_pr_i_{1+nb1},:) , zeros(nrows_A_(1+nb1),8) , repmat(T_u_{1+nb1}(X3_pr_i_{1+nb1}),1,8) , zeros(nrows_A_(1+nb1),8) , An{1+nb1}(X3_pr_i_{1+nb1},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X3, nb1 %d',nb1));
end;%for nb1=0:nbins-1;
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
figure; clf;
for nb1=0:nbins-1;
% X5;
subplot(3,nbins,1+nb1+0*nbins);
tmp_c = X5_pc_i;
tmp_r = X5_pr_i_{1+nb1};
imagesc([ ...
	  Tn{1+nb1}(X5_pr_i_{1+nb1},:) , zeros(nrows_A_(1+nb1),8) , repmat(T_u_{1+nb1}(X5_pr_i_{1+nb1}),1,8) , zeros(nrows_A_(1+nb1),8) , An{1+nb1}(X5_pr_i_{1+nb1},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X5, nb1 %d',nb1));
% X6;
subplot(3,nbins,1+nb1+1*nbins);
tmp_c = X6_pc_i;
tmp_r = X6_pr_i_{1+nb1};
imagesc([ ...
	  Tn{1+nb1}(X6_pr_i_{1+nb1},:) , zeros(nrows_A_(1+nb1),8) , repmat(T_u_{1+nb1}(X6_pr_i_{1+nb1}),1,8) , zeros(nrows_A_(1+nb1),8) , An{1+nb1}(X6_pr_i_{1+nb1},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X6, nb1 %d',nb1));
% X0;
subplot(3,nbins,1+nb1+2*nbins);
tmp_c = X0_pc_i;
tmp_r = X0_pr_i_{1+nb1};
imagesc([ ...
	  Tn{1+nb1}(X0_pr_i_{1+nb1},:) , zeros(nrows_A_(1+nb1),8) , repmat(T_u_{1+nb1}(X0_pr_i_{1+nb1}),1,8) , zeros(nrows_A_(1+nb1),8) , An{1+nb1}(X0_pr_i_{1+nb1},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X0, nb1 %d',nb1));
end;%for nb1=0:nbins-1;
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
figure; clf;
for nb1=0:nbins-1;
% X0;
subplot(1,nbins,1+nb1);
tmp_c = X0_pc_i;
tmp_r = X0_pr_i_{1+nb1};
imagesc([ ...
	  An{1+nb1}(tmp_r,tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X0, nb1 %d',nb1));
end;%for nb1=0:nbins-1;
return;
end;%if disp_flag;

for nb1=0:nbins-1; 
At{1+nb1} = transpose(An{1+nb1}); Zt{1+nb1} = transpose(Zn{1+nb1});
Tt{1+nb1} = transpose(Tn{1+nb1}); St{1+nb1} = transpose(Sn{1+nb1});
end;%for nb1=0:nbins-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear mc_K mr_J_;
mc_K = mc_A; mc_K(X0_pc(1:floor(ncols_X/2)))=0;
for nb1=0:nbins-1;
mr_J_{1+nb1} = mr_A_{1+nb1}; mr_J_{1+nb1}(X0_pr_i_{1+nb1}(1:floor(nrows_X/nbins/2)))=0;
end;%for nb1=0:nbins-1;
tmpKtchar = sprintf('%s/mc_K.b16',dir__in);tutorial_binary_compress(bitj,mc_K(:)>0,tmpKtchar); %disp(sprintf('mc_K: %s',num2str(mc_K)));
for nb1=0:nbins-1;
tmpJtchar = sprintf('%s/mr_J_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_J_{1+nb1}(:)>0,tmpJtchar); %disp(sprintf('mr_J_[%d]: %s',nb1,num2str(mr_J_{1+nb1})));
end;%for nb1=0:nbins-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpAtchar = sprintf('%s/mc_A.b16',dir__in);tutorial_binary_compress(bitj,mc_A(:)>0,tmpAtchar); %disp(sprintf('mc_A: %s',num2str(mc_A)));
tmpTtchar = sprintf('%s/mc_T.b16',dir__in);tutorial_binary_compress(bitj,mc_T(:)>0,tmpTtchar); %disp(sprintf('mc_T: %s',num2str(mc_T)));
for nb1=0:nbins-1;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,An{1+nb1}>0,tmpAnchar); 
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,At{1+nb1}>0,tmpAtchar); 
tmpZnchar = sprintf('%s/Z_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Zn{1+nb1}>0,tmpZnchar); 
tmpZtchar = sprintf('%s/Z_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Zt{1+nb1}>0,tmpZtchar); 
tmpTnchar = sprintf('%s/T_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Tn{1+nb1}>0,tmpTnchar); 
tmpTtchar = sprintf('%s/T_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Tt{1+nb1}>0,tmpTtchar); 
tmpSnchar = sprintf('%s/S_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Sn{1+nb1}>0,tmpSnchar); 
tmpStchar = sprintf('%s/S_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,St{1+nb1}>0,tmpStchar); 
tmpAtchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_A_{1+nb1}(:)>0,tmpAtchar); %disp(sprintf('mr_A_[%d]: %s',nb1,num2str(mr_A_{1+nb1})));
tmpZtchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_Z_{1+nb1}(:)>0,tmpZtchar); %disp(sprintf('mr_Z_[%d]: %s',nb1,num2str(mr_Z_{1+nb1})));
end;%for nb1=0:nbins-1;

fname__in = sprintf('%s/%s.in',dir__in,in__name);
fp = fopen(fname__in,'w');
fprintf(fp,'GLOBAL_verbose= 0;\n');
fprintf(fp,'GLOBAL_thread_count= 8;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','dcgxpander_driver');
fprintf(fp,'GLOBAL_NBINS= %d;\n',nbins);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',B_MLT);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
%fprintf(fp,'GLOBAL_TEST_sparse= 0;\n');
fprintf(fp,'GLOBAL_Ireq= 0;\n');
fprintf(fp,'GLOBAL_DIR_XPRE= %s;\n',dir__in);
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_A_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',ncols_A);
fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'mr_A_%.2d.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cind= mc_A.b16;\n',dir__in);
fprintf(fp,'GLOBAL_Z_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'Z_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'Z_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_Z_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'mr_Z_%.2d.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'T_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'T_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_T_n_cind= mc_T.b16;\n',dir__in);
fprintf(fp,'GLOBAL_S_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'S_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_S_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'S_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_J_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'mr_J_%.2d.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_K_n_cind= mc_K.b16;\n',dir__in);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by dcgxpander_test_uAZSZA_ver0.m on %s;\n',date);
fclose(fp);

call_flag=1;%call_flag = input(' call? 1=yes (default), 0=no:'); if isempty(call_flag); call_flag=1; end;
if call_flag;
disp(sprintf('%s/../lakcluster_ver18 < %s',dir_trunk,fname__in));
system(sprintf('%s/../lakcluster_ver18 < %s',dir_trunk,fname__in));
end;%if call_flag;

%disp(sprintf('%s/dims_X',dir_out));
run(sprintf('%s/timing',dir_out)); %run(sprintf('%s/dims_X',dir_out));
tmpAtchar = sprintf('%s/mc_A.b16',dir__in);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAtchar); mc_A_tmp = tutorial_binary_uncompress(tmpAtchar,1:nrows,1:ncols)>0; 
for nb1=0:nbins-1; 
tmpAnchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAnchar); mr_A_tmp_{1+nb1} = tutorial_binary_uncompress(tmpAnchar,1:nrows,1:ncols)>0; 
tmpZnchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpZnchar); mr_Z_tmp_{1+nb1} = tutorial_binary_uncompress(tmpZnchar,1:nrows,1:ncols)>0; 
end;%for nb1=0:nbins-1; 

out_xgrow_a = textread(sprintf('%s/out_xgrow_a.txt',dir_out)); out_xgrow_b = textread(sprintf('%s/out_xgrow_b.txt',dir_out));
r_ij = out_xgrow_a(find(out_xgrow_a(:,1)>-1),1); c_ij = out_xgrow_a(find(out_xgrow_a(:,2)>-1),2);

clear r_on_ z_on_ c_on y_on ;
for nb1=0:nbins-1; r_on_{1+nb1} = find(mr_A_tmp_{1+nb1}); z_on_{1+nb1} = find(mr_Z_tmp_{1+nb1}); end;%for nb1=0:nbins-1;
c_on = find(mc_A_tmp); 
r_on_x = r_on_{1}; for nb1=1:nbins-1; r_on_x = unionall({r_on_x,r_on_{1+nb1}+sum(nrows_A_(1:1+nb1-1))}); end;
z_on_x = z_on_{1}; for nb1=1:nbins-1; z_on_x = unionall({z_on_x,z_on_{1+nb1}+sum(nrows_Z_(1:1+nb1-1))}); end;

XW_pc_h = intersect(c_on,c_ij_i(XW_pc(1:ncols_X))); AW_pc_h = setdiff(c_on,XW_pc_h); 
auc_c_W = get_auc(AW_pc_h,XW_pc_h);

X0_pc_h = intersect(c_on,c_ij_i(X0_pc(1:ncols_X))); A0_pc_h = setdiff(c_on,X0_pc_h); 
auc_c_0 = get_auc(A0_pc_h,X0_pc_h);
%{
X1_pc_h = intersect(c_on,c_ij_i(X1_pc(1:ncols_XX))); A1_pc_h = setdiff(c_on,X1_pc_h); 
auc_c_1 = get_auc(A1_pc_h,X1_pc_h);
X2_pc_h = intersect(c_on,c_ij_i(X2_pc(1:ncols_XX))); A2_pc_h = setdiff(c_on,X2_pc_h); 
auc_c_2 = get_auc(A2_pc_h,X2_pc_h);
X3_pc_h = intersect(c_on,c_ij_i(X3_pc(1:ncols_XX))); A3_pc_h = setdiff(c_on,X3_pc_h); 
auc_c_3 = get_auc(A3_pc_h,X3_pc_h);
X4_pc_h = intersect(c_on,c_ij_i(X4_pc(1:ncols_XX))); A4_pc_h = setdiff(c_on,X4_pc_h); 
auc_c_4 = get_auc(A4_pc_h,X4_pc_h);
X5_pc_h = intersect(c_on,c_ij_i(X5_pc(1:ncols_XX))); A5_pc_h = setdiff(c_on,X5_pc_h); 
auc_c_5 = get_auc(A5_pc_h,X5_pc_h);
X6_pc_h = intersect(c_on,c_ij_i(X6_pc(1:ncols_XX))); A6_pc_h = setdiff(c_on,X6_pc_h); 
auc_c_6 = get_auc(A6_pc_h,X6_pc_h);
 %}

nrows_A_sum = cumsum([0,nrows_A_]);
for nb1=0:nbins-1;
XW_pr_g_{1+nb1} = XW_pr_{1+nb1} + nrows_A_sum(1+nb1);
X0_pr_g_{1+nb1} = X0_pr_{1+nb1} + nrows_A_sum(1+nb1);
%X1_pr_g_{1+nb1} = X1_pr_{1+nb1} + nrows_A_sum(1+nb1);
%X2_pr_g_{1+nb1} = X2_pr_{1+nb1} + nrows_A_sum(1+nb1);
%X3_pr_g_{1+nb1} = X3_pr_{1+nb1} + nrows_A_sum(1+nb1);
%X4_pr_g_{1+nb1} = X4_pr_{1+nb1} + nrows_A_sum(1+nb1);
%X5_pr_g_{1+nb1} = X5_pr_{1+nb1} + nrows_A_sum(1+nb1);
%X6_pr_g_{1+nb1} = X6_pr_{1+nb1} + nrows_A_sum(1+nb1);
end;%for nb1=0:nbins-1;

%{
% XW;
for nb1=0:nbins-1;
nb1=1-1; XW_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(XW_pr_g_{1+nb1}(1:floor(nrows_X/3)))); AW_pr_h_{1+nb1} = setdiff(r_on_x,XW_pr_h_{1+nb1});
auc_r_W(1+nb1) = get_auc(AW_pr_h_{1+nb1},XW_pr_h_{1+nb1});
end;%for nb1=0:nbins-1;
auc_r_W_x = mean(auc_r_W);
 %}

%{
% X0;
nb1=1-1; X0_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X0_pr_g_{1+nb1}(1:floor(nrows_X/3)))); A0_pr_h_{1+nb1} = setdiff(r_on_x,X0_pr_h_{1+nb1});
auc_r_0(1+nb1) = get_auc(A0_pr_h_{1+nb1},X0_pr_h_{1+nb1});
nb1=2-1; X0_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X0_pr_g_{1+nb1}(1:floor(nrows_X/3)))); A0_pr_h_{1+nb1} = setdiff(r_on_x,X0_pr_h_{1+nb1});
auc_r_0(1+nb1) = get_auc(A0_pr_h_{1+nb1},X0_pr_h_{1+nb1});
nb1=3-1; X0_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X0_pr_g_{1+nb1}(1:floor(nrows_X/3)))); A0_pr_h_{1+nb1} = setdiff(r_on_x,X0_pr_h_{1+nb1});
auc_r_0(1+nb1) = get_auc(A0_pr_h_{1+nb1},X0_pr_h_{1+nb1});
auc_r_0_x = squeeze((auc_r_0(1) + auc_r_0(2) + auc_r_0(3))/3);
% X1;
nb1=2-1; X1_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X1_pr_g_{1+nb1}(1:floor(nrows_XX/2)))); A1_pr_h_{1+nb1} = setdiff(r_on_x,X1_pr_h_{1+nb1});
auc_r_1(1+nb1) = get_auc(A1_pr_h_{1+nb1},X1_pr_h_{1+nb1});
nb1=3-1; X1_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X1_pr_g_{1+nb1}(1:floor(nrows_XX/2)))); A1_pr_h_{1+nb1} = setdiff(r_on_x,X1_pr_h_{1+nb1});
auc_r_1(1+nb1) = get_auc(A1_pr_h_{1+nb1},X1_pr_h_{1+nb1});
auc_r_1_x = squeeze((auc_r_1(2) + auc_r_1(3))/2);
% X2;
nb1=1-1; X2_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X2_pr_g_{1+nb1}(1:floor(nrows_XX/2)))); A2_pr_h_{1+nb1} = setdiff(r_on_x,X2_pr_h_{1+nb1});
auc_r_2(1+nb1) = get_auc(A2_pr_h_{1+nb1},X2_pr_h_{1+nb1});
nb1=3-1; X2_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X2_pr_g_{1+nb1}(1:floor(nrows_XX/2)))); A2_pr_h_{1+nb1} = setdiff(r_on_x,X2_pr_h_{1+nb1});
auc_r_2(1+nb1) = get_auc(A2_pr_h_{1+nb1},X2_pr_h_{1+nb1});
auc_r_2_x = squeeze((auc_r_2(1) + auc_r_2(3))/2);
% X3;
nb1=1-1; X3_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X3_pr_g_{1+nb1}(1:floor(nrows_XX/2)))); A3_pr_h_{1+nb1} = setdiff(r_on_x,X3_pr_h_{1+nb1});
auc_r_3(1+nb1) = get_auc(A3_pr_h_{1+nb1},X3_pr_h_{1+nb1});
nb1=2-1; X3_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X3_pr_g_{1+nb1}(1:floor(nrows_XX/2)))); A3_pr_h_{1+nb1} = setdiff(r_on_x,X3_pr_h_{1+nb1});
auc_r_3(1+nb1) = get_auc(A3_pr_h_{1+nb1},X3_pr_h_{1+nb1});
auc_r_3_x = squeeze((auc_r_3(1) + auc_r_3(2))/2);
% X4;
nb1=1-1; X4_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X4_pr_g_{1+nb1}(1:floor(nrows_XX/3)))); A4_pr_h_{1+nb1} = setdiff(r_on_x,X4_pr_h_{1+nb1});
auc_r_4(1+nb1) = get_auc(A4_pr_h_{1+nb1},X4_pr_h_{1+nb1});
nb1=2-1; X4_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X4_pr_g_{1+nb1}(1:floor(nrows_XX/3)))); A4_pr_h_{1+nb1} = setdiff(r_on_x,X4_pr_h_{1+nb1});
auc_r_4(1+nb1) = get_auc(A4_pr_h_{1+nb1},X4_pr_h_{1+nb1});
nb1=3-1; X4_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X4_pr_g_{1+nb1}(1:floor(nrows_XX/3)))); A4_pr_h_{1+nb1} = setdiff(r_on_x,X4_pr_h_{1+nb1});
auc_r_4(1+nb1) = get_auc(A4_pr_h_{1+nb1},X4_pr_h_{1+nb1});
auc_r_4_x = squeeze((auc_r_4(1) + auc_r_4(2) + auc_r_4(3))/3);
% X5;
nb1=1-1; X5_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X5_pr_g_{1+nb1}(1:floor(nrows_XX/3)))); A5_pr_h_{1+nb1} = setdiff(r_on_x,X5_pr_h_{1+nb1});
auc_r_5(1+nb1) = get_auc(A5_pr_h_{1+nb1},X5_pr_h_{1+nb1});
nb1=2-1; X5_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X5_pr_g_{1+nb1}(1:floor(nrows_XX/3)))); A5_pr_h_{1+nb1} = setdiff(r_on_x,X5_pr_h_{1+nb1});
auc_r_5(1+nb1) = get_auc(A5_pr_h_{1+nb1},X5_pr_h_{1+nb1});
nb1=3-1; X5_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X5_pr_g_{1+nb1}(1:floor(nrows_XX/3)))); A5_pr_h_{1+nb1} = setdiff(r_on_x,X5_pr_h_{1+nb1});
auc_r_5(1+nb1) = get_auc(A5_pr_h_{1+nb1},X5_pr_h_{1+nb1});
auc_r_5_x = squeeze((auc_r_5(1) + auc_r_5(2) + auc_r_5(3))/3);
% X6;
nb1=1-1; X6_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X6_pr_g_{1+nb1}(1:floor(nrows_XX/3)))); A6_pr_h_{1+nb1} = setdiff(r_on_x,X6_pr_h_{1+nb1});
auc_r_6(1+nb1) = get_auc(A6_pr_h_{1+nb1},X6_pr_h_{1+nb1});
nb1=2-1; X6_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X6_pr_g_{1+nb1}(1:floor(nrows_XX/3)))); A6_pr_h_{1+nb1} = setdiff(r_on_x,X6_pr_h_{1+nb1});
auc_r_6(1+nb1) = get_auc(A6_pr_h_{1+nb1},X6_pr_h_{1+nb1});
nb1=3-1; X6_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X6_pr_g_{1+nb1}(1:floor(nrows_XX/3)))); A6_pr_h_{1+nb1} = setdiff(r_on_x,X6_pr_h_{1+nb1});
auc_r_6(1+nb1) = get_auc(A6_pr_h_{1+nb1},X6_pr_h_{1+nb1});
auc_r_6_x = squeeze((auc_r_6(1) + auc_r_6(2) + auc_r_6(3))/3);

auc_x_0 = 0.5*(auc_r_0_x + auc_c_0) ;
auc_x_1 = 0.5*(auc_r_1_x + auc_c_1) ;
auc_x_2 = 0.5*(auc_r_2_x + auc_c_2) ;
auc_x_3 = 0.5*(auc_r_3_x + auc_c_3) ;
auc_x_4 = 0.5*(auc_r_4_x + auc_c_4) ;
auc_x_5 = 0.5*(auc_r_5_x + auc_c_5) ;
auc_x_6 = 0.5*(auc_r_6_x + auc_c_6) ;

xx_ra = [ ...
          auc_x_0   , auc_x_1   , auc_x_2   , auc_x_3   , auc_x_4   , auc_x_5   , auc_x_6   ; ...
          auc_r_0_x , auc_r_1_x , auc_r_2_x , auc_r_3_x , auc_r_4_x , auc_r_5_x , auc_r_6_x ; ...
          auc_c_0   , auc_c_1   , auc_c_2   , auc_c_3   , auc_c_4   , auc_c_5   , auc_c_6   ; ...
	  ];
disp(num2str(xx_ra));
 %}

% X0;
nb1=1-1; X0_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X0_pr_g_{1+nb1}(1:floor(nrows_X/nbins)))); A0_pr_h_{1+nb1} = setdiff(r_on_x,X0_pr_h_{1+nb1});
auc_r_0(1+nb1) = get_auc(A0_pr_h_{1+nb1},X0_pr_h_{1+nb1});
auc_r_0_x = auc_r_0(1);

auc_x_0 = 0.5*(auc_r_0_x + auc_c_0) ;
xx_ra = [ auc_x_0 ; auc_r_0_x ; auc_c_0 ];
disp(num2str(xx_ra));

xx_ra = xx_ra(:);

if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

clean_flag=1;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;
