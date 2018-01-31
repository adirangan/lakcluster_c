function lakcluster_test_uADZSZDA_Ireq2_ver7(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num,verbose_flag)
% creating a larger example involving multiple categorical and continuous covariates, as well as cases and controls (but not involving Y). ;
% uses lakcluster_ver16 ;
% try: ;
%{
  
  N=1024;X_factor=0.55;X_esm=0.1;gamma=0.01;B_MLT=32;rng_num=1;lakcluster_test_uADZSZDA_Ireq2_ver7(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

  %}

if (nargin<9); verbose_flag=0; end;
if (nargin<8); pt_num = -1; end;
if (nargin<7); rng_num = 1; end;
if (nargin<4); B_MLT = 32; end;
if (nargin<5); gamma = 0.025; end;
if (nargin<4); X_esm = 0.1; end;
if (nargin<3); X_factor = 0.65; end;
if (nargin<2); N = 1024; end;
if (nargin<1); dir_trunk = pwd; end;

nbins = 3;
nrows_A_ = ceil(N/nbins*ones(1,nbins));
nrows_Z_ = ceil(2*N/nbins*ones(1,nbins));
ncols_A3 = ceil(N/1920)*1920;
ncols_A = 3*ncols_A3;
ncols_Y = 0;
nmds = 2; 
ncols_T = nmds+1; 
for nb1=0:nbins-1; MDT_{1+nb1} = randn(nrows_A_(1+nb1),nmds); MDS_{1+nb1} = randn(nrows_Z_(1+nb1),nmds); end;%for nb1=0:nbins-1;
mrnd = 0;
QR_strategy = 'YnWt';
QC_strategy = 'YnWt store one';

nrows_X = ceil(N.^(X_factor)); nrows_XX = ceil(1.5*nrows_X);
ncols_X = ceil((ncols_A).^(X_factor)); ncols_XX = ceil(1.5*ncols_X);
X_eps = X_esm/sqrt(nrows_X); % spectral error for X;
X_pp = tutorial_g1em(X_esm);
Ireq = 2;

test_string = 'lakcluster_test_uADZSZDA_Ireq2';
dir_base = sprintf('%s/dir_%s',dir_trunk,test_string); if ~exist(dir_base,'dir'); mkdir(dir_base); end;
[posfix_use,X_factor_d,X_esm_d,gamma_d] = xfix_gen(N,X_factor,X_esm,gamma,B_MLT,rng_num);
in__name = sprintf('%s_%s__in',test_string,posfix_use);
out_name = sprintf('%s_%s_out',test_string,posfix_use);
dir__in = sprintf('%s/dir_%s',dir_base,in__name); if ~exist(dir__in,'dir'); mkdir(dir__in); end;
dir_out = sprintf('%s/dir_%s',dir_base,out_name); if ~exist(dir_out,'dir'); mkdir(dir_out); end; 
%path(path,dir__in); path(path,dir_out);
rng(rng_num);

clear An At Zn Zt Yn Yt Wn Wt Tn Tt Sn St mr_A mr_Z mc_A mc_Y mc_T bitj;
bitj = 16;

mc_A = (rand(1,ncols_A)>mrnd);
mc_Y = (rand(1,ncols_Y)>mrnd);
mc_T = (rand(1,ncols_T)>mrnd); mc_T(1)=1; % ensure first column of [T;S] is all ones ;
for nb1=0:nbins-1; 
mr_A{1+nb1} = (rand(nrows_A_(1+nb1),1)>mrnd);
mr_Z{1+nb1} = (rand(nrows_Z_(1+nb1),1)>mrnd);
end;%for nb1=0:nbins-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate random matrices ;
sx_ra(1) = 0.125; sx_start(1) = 0*ncols_A3;
sx_ra(2) = 0.250; sx_start(2) = 1*ncols_A3;
sx_ra(3) = 0.500; sx_start(3) = 2*ncols_A3;
sx_mid = 0.250;
for nb1=0:nbins-1; 
An{1+nb1} = [2*(rand(nrows_A_(1+nb1),ncols_A3)<sx_ra(1))-1 , 2*(rand(nrows_A_(1+nb1),ncols_A3)<sx_ra(2))-1 , 2*(rand(nrows_A_(1+nb1),ncols_A3)<sx_ra(3))-1];
Zn{1+nb1} = [2*(rand(nrows_Z_(1+nb1),ncols_A3)<sx_ra(1))-1 , 2*(rand(nrows_Z_(1+nb1),ncols_A3)<sx_ra(2))-1 , 2*(rand(nrows_Z_(1+nb1),ncols_A3)<sx_ra(3))-1];
Tn{1+nb1} = [ones(nrows_A_(1+nb1),1) , 2*(MDT_{1+nb1}>0)-1];
Sn{1+nb1} = [ones(nrows_Z_(1+nb1),1) , 2*(MDS_{1+nb1}>0)-1];
end;%for nb1=0:nbins-1; 
% generate clusters ;
% XW is a simple nonspecific cluster ;
XW = tutorial_binary_makelr3(nrows_X,ncols_X,X_pp,0.5);
for nb1=0:nbins-1; XW_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,XW_pr_i_{1+nb1}] = sort(XW_pr_{1+nb1}); [tmp,XW_pr_i_{1+nb1}] = sort(XW_pr_i_{1+nb1}); end;
XW_pc = randperm(ncols_A);
tmp_c = [XW_pc]; [tmp,XW_pc_i] = sort(tmp_c); [tmp,XW_pc_i] = sort(XW_pc_i);

u_mds = randn(nmds,1); u_mds = u_mds/norm(u_mds);
for nb1=0:nbins-1;
T_u_{1+nb1} = MDT_{1+nb1}*u_mds; T_u_pos_{1+nb1} = find(T_u_{1+nb1}>0); T_u_neg_{1+nb1} = find(T_u_{1+nb1}<=0);
S_u_{1+nb1} = MDS_{1+nb1}*u_mds; S_u_pos_{1+nb1} = find(S_u_{1+nb1}>0); S_u_neg_{1+nb1} = find(S_u_{1+nb1}<=0);
end;%for nb1=0:nbins-1;
% generate clusters ;
% X0 is balanced across covariate categories and T, and is exclusive to the cases ;
X0 = tutorial_binary_makelr3(nrows_X,ncols_X,X_pp,sx_mid);
for nb1=0:nbins-1; X0_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X0_pr_i_{1+nb1}] = sort(X0_pr_{1+nb1}); [tmp,X0_pr_i_{1+nb1}] = sort(X0_pr_i_{1+nb1}); end;
X0_pc = randperm(ncols_A); [tmp,X0_pc_i] = sort(X0_pc); [tmp,X0_pc_i] = sort(X0_pc_i);
% X1 is balanced across T and is exclusive to the cases, but is restricted to covariate categories 2,3 ;
X1 = tutorial_binary_makelr3(nrows_XX,ncols_XX,X_pp,sx_mid);
for nb1=0:nbins-1; X1_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X1_pr_i_{1+nb1}] = sort(X1_pr_{1+nb1}); [tmp,X1_pr_i_{1+nb1}] = sort(X1_pr_i_{1+nb1}); end;
X1_pc = randperm(ncols_A); [tmp,X1_pc_i] = sort(X1_pc); [tmp,X1_pc_i] = sort(X1_pc_i);
% X2 is balanced across T and is exclusive to the cases, but is restricted to covariate categories 1,3 ;
X2 = tutorial_binary_makelr3(nrows_XX,ncols_XX,X_pp,sx_mid);
for nb1=0:nbins-1; X2_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X2_pr_i_{1+nb1}] = sort(X2_pr_{1+nb1}); [tmp,X2_pr_i_{1+nb1}] = sort(X2_pr_i_{1+nb1}); end;
X2_pc = randperm(ncols_A); [tmp,X2_pc_i] = sort(X2_pc); [tmp,X2_pc_i] = sort(X2_pc_i);
% X3 is balanced across T and is exclusive to the cases, but is restricted to covariate categories 1,2 ;
X3 = tutorial_binary_makelr3(nrows_XX,ncols_XX,X_pp,sx_mid);
for nb1=0:nbins-1; X3_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X3_pr_i_{1+nb1}] = sort(X3_pr_{1+nb1}); [tmp,X3_pr_i_{1+nb1}] = sort(X3_pr_i_{1+nb1}); end;
X3_pc = randperm(ncols_A); [tmp,X3_pc_i] = sort(X3_pc); [tmp,X3_pc_i] = sort(X3_pc_i);
% X4 is balanced across covariate categories and T, but spans both the cases and controls ;
X4 = tutorial_binary_makelr3(3*nrows_XX,ncols_XX,X_pp,sx_mid);
for nb1=0:nbins-1; X4_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X4_pr_i_{1+nb1}] = sort(X4_pr_{1+nb1}); [tmp,X4_pr_i_{1+nb1}] = sort(X4_pr_i_{1+nb1}); end;
for nb1=0:nbins-1; X4_pz_{1+nb1} = randperm(nrows_Z_(1+nb1)); [tmp,X4_pz_i_{1+nb1}] = sort(X4_pz_{1+nb1}); [tmp,X4_pz_i_{1+nb1}] = sort(X4_pz_i_{1+nb1}); end;
X4_pc = randperm(ncols_A); [tmp,X4_pc_i] = sort(X4_pc); [tmp,X4_pc_i] = sort(X4_pc_i);
% X5 is balanced across covariate categories and is exclusive to the cases, but is restricted to positive T_u ;
X5 = tutorial_binary_makelr3(nrows_XX,ncols_XX,X_pp,sx_mid);
X5_pc = randperm(ncols_A); [tmp,X5_pc_i] = sort(X5_pc); [tmp,X5_pc_i] = sort(X5_pc_i);
for nb1=0:nbins-1; 
X5_pr_{1+nb1} = [ T_u_pos_{1+nb1}(randperm(length(T_u_pos_{1+nb1}))) ; T_u_neg_{1+nb1}(randperm(length(T_u_neg_{1+nb1}))) ];
[tmp,X5_pr_i_{1+nb1}] = sort(X5_pr_{1+nb1}); [tmp,X5_pr_i_{1+nb1}] = sort(X5_pr_i_{1+nb1});
end;%for nb1=0:nbins-1; 
% X6 is balanced across covariate categories and is exclusive to the cases, but is restricted to negative T_u ;
X6 = tutorial_binary_makelr3(nrows_XX,ncols_XX,X_pp,sx_mid);
X6_pc = randperm(ncols_A); [tmp,X6_pc_i] = sort(X6_pc); [tmp,X6_pc_i] = sort(X6_pc_i);
for nb1=0:nbins-1; 
X6_pr_{1+nb1} = [ T_u_neg_{1+nb1}(randperm(length(T_u_neg_{1+nb1}))) ; T_u_pos_{1+nb1}(randperm(length(T_u_pos_{1+nb1}))) ];
[tmp,X6_pr_i_{1+nb1}] = sort(X6_pr_{1+nb1}); [tmp,X6_pr_i_{1+nb1}] = sort(X6_pr_i_{1+nb1});
end;%for nb1=0:nbins-1; 

%insert clusters;
% insert XW ;
% for (nb1=0:nbins-1); An{1+nb1}(XW_pr_{1+nb1}(0*floor(nrows_X/nbins) + (1:floor(nrows_X/nbins))),XW_pc(1:ncols_X)) = XW(nb1*floor(nrows_X/nbins) + (1:floor(nrows_X/nbins)),1:ncols_X); end;

% insert X0;
nb1=1-1; An{1+nb1}(X0_pr_{1+nb1}(0*floor(nrows_X/3) + (1:floor(nrows_X/3))),X0_pc(1:ncols_X)) = X0(0*floor(nrows_X/3) + (1:floor(nrows_X/3)),1:ncols_X);
nb1=2-1; An{1+nb1}(X0_pr_{1+nb1}(0*floor(nrows_X/3) + (1:floor(nrows_X/3))),X0_pc(1:ncols_X)) = X0(1*floor(nrows_X/3) + (1:floor(nrows_X/3)),1:ncols_X);
nb1=3-1; An{1+nb1}(X0_pr_{1+nb1}(0*floor(nrows_X/3) + (1:floor(nrows_X/3))),X0_pc(1:ncols_X)) = X0(2*floor(nrows_X/3) + (1:floor(nrows_X/3)),1:ncols_X);
% insert X1 ;
nb1=2-1; An{1+nb1}(X1_pr_{1+nb1}(0*floor(nrows_XX/2) + (1:floor(nrows_XX/2))),X1_pc(1:ncols_XX)) = X1(0*floor(nrows_XX/2) + (1:floor(nrows_XX/2)),1:ncols_XX);
nb1=3-1; An{1+nb1}(X1_pr_{1+nb1}(0*floor(nrows_XX/2) + (1:floor(nrows_XX/2))),X1_pc(1:ncols_XX)) = X1(1*floor(nrows_XX/2) + (1:floor(nrows_XX/2)),1:ncols_XX);
% insert X2 ;
nb1=1-1; An{1+nb1}(X2_pr_{1+nb1}(0*floor(nrows_XX/2) + (1:floor(nrows_XX/2))),X2_pc(1:ncols_XX)) = X2(0*floor(nrows_XX/2) + (1:floor(nrows_XX/2)),1:ncols_XX);
nb1=3-1; An{1+nb1}(X2_pr_{1+nb1}(0*floor(nrows_XX/2) + (1:floor(nrows_XX/2))),X2_pc(1:ncols_XX)) = X2(1*floor(nrows_XX/2) + (1:floor(nrows_XX/2)),1:ncols_XX);
% insert X3 ;
nb1=1-1; An{1+nb1}(X3_pr_{1+nb1}(0*floor(nrows_XX/2) + (1:floor(nrows_XX/2))),X3_pc(1:ncols_XX)) = X3(0*floor(nrows_XX/2) + (1:floor(nrows_XX/2)),1:ncols_XX);
nb1=2-1; An{1+nb1}(X3_pr_{1+nb1}(0*floor(nrows_XX/2) + (1:floor(nrows_XX/2))),X3_pc(1:ncols_XX)) = X3(1*floor(nrows_XX/2) + (1:floor(nrows_XX/2)),1:ncols_XX);
% insert X4 ;
nb1=1-1; An{1+nb1}(X4_pr_{1+nb1}(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3))),X4_pc(1:ncols_XX)) = X4(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3)),1:ncols_XX);
nb1=2-1; An{1+nb1}(X4_pr_{1+nb1}(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3))),X4_pc(1:ncols_XX)) = X4(1*floor(nrows_XX/3) + (1:floor(nrows_XX/3)),1:ncols_XX);
nb1=3-1; An{1+nb1}(X4_pr_{1+nb1}(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3))),X4_pc(1:ncols_XX)) = X4(2*floor(nrows_XX/3) + (1:floor(nrows_XX/3)),1:ncols_XX);
nb1=1-1; Zn{1+nb1}(X4_pz_{1+nb1}(0*floor(nrows_XX*2/3) + (1:floor(nrows_XX*2/3))),X4_pc(1:ncols_XX)) = X4(3*floor(nrows_XX/3) + (1:floor(nrows_XX*2/3)),1:ncols_XX);
nb1=2-1; Zn{1+nb1}(X4_pz_{1+nb1}(0*floor(nrows_XX*2/3) + (1:floor(nrows_XX*2/3))),X4_pc(1:ncols_XX)) = X4(5*floor(nrows_XX/3) + (1:floor(nrows_XX*2/3)),1:ncols_XX);
nb1=3-1; Zn{1+nb1}(X4_pz_{1+nb1}(0*floor(nrows_XX*2/3) + (1:floor(nrows_XX*2/3))),X4_pc(1:ncols_XX)) = X4(7*floor(nrows_XX/3) + (1:floor(nrows_XX*2/3)),1:ncols_XX);
% insert X5 ;
nb1=1-1; An{1+nb1}(X5_pr_{1+nb1}(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3))),X5_pc(1:ncols_XX)) = X5(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3)),1:ncols_XX);
nb1=2-1; An{1+nb1}(X5_pr_{1+nb1}(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3))),X5_pc(1:ncols_XX)) = X5(1*floor(nrows_XX/3) + (1:floor(nrows_XX/3)),1:ncols_XX);
nb1=3-1; An{1+nb1}(X5_pr_{1+nb1}(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3))),X5_pc(1:ncols_XX)) = X5(2*floor(nrows_XX/3) + (1:floor(nrows_XX/3)),1:ncols_XX);
% insert X6 ;
nb1=1-1; An{1+nb1}(X6_pr_{1+nb1}(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3))),X6_pc(1:ncols_XX)) = X6(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3)),1:ncols_XX);
nb1=2-1; An{1+nb1}(X6_pr_{1+nb1}(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3))),X6_pc(1:ncols_XX)) = X6(1*floor(nrows_XX/3) + (1:floor(nrows_XX/3)),1:ncols_XX);
nb1=3-1; An{1+nb1}(X6_pr_{1+nb1}(0*floor(nrows_XX/3) + (1:floor(nrows_XX/3))),X6_pc(1:ncols_XX)) = X6(2*floor(nrows_XX/3) + (1:floor(nrows_XX/3)),1:ncols_XX);

% generate transposes ;
for nb1=0:nbins-1;
At{1+nb1} = transpose(An{1+nb1}); Zt{1+nb1} = transpose(Zn{1+nb1});
Tt{1+nb1} = transpose(Tn{1+nb1}); St{1+nb1} = transpose(Sn{1+nb1});
end;%for nb1=0:nbins-1;

disp_flag=0;
if disp_flag;
% X1, X2, X3;
tmp_c = X1_pc_i;
figure; 
imagesc([ ...
	  Tn{1+0}(X1_pr_i_{1+0},:) , zeros(nrows_A_(1+0),8) , repmat(T_u_{1+0}(X1_pr_i_{1+0}),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(X1_pr_i_{1+0},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+1}(X1_pr_i_{1+1},:) , zeros(nrows_A_(1+1),8) , repmat(T_u_{1+1}(X1_pr_i_{1+1}),1,8) , zeros(nrows_A_(1+1),8) , An{1+1}(X1_pr_i_{1+1},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+2}(X1_pr_i_{1+2},:) , zeros(nrows_A_(1+2),8) , repmat(T_u_{1+2}(X1_pr_i_{1+2}),1,8) , zeros(nrows_A_(1+2),8) , An{1+2}(X1_pr_i_{1+2},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X1'));
tmp_c = X2_pc_i;
figure; 
imagesc([ ...
	  Tn{1+0}(X2_pr_i_{1+0},:) , zeros(nrows_A_(1+0),8) , repmat(T_u_{1+0}(X2_pr_i_{1+0}),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(X2_pr_i_{1+0},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+1}(X2_pr_i_{1+1},:) , zeros(nrows_A_(1+1),8) , repmat(T_u_{1+1}(X2_pr_i_{1+1}),1,8) , zeros(nrows_A_(1+1),8) , An{1+1}(X2_pr_i_{1+1},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+2}(X2_pr_i_{1+2},:) , zeros(nrows_A_(1+2),8) , repmat(T_u_{1+2}(X2_pr_i_{1+2}),1,8) , zeros(nrows_A_(1+2),8) , An{1+2}(X2_pr_i_{1+2},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X2'));
tmp_c = X3_pc_i;
figure; 
imagesc([ ...
	  Tn{1+0}(X3_pr_i_{1+0},:) , zeros(nrows_A_(1+0),8) , repmat(T_u_{1+0}(X3_pr_i_{1+0}),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(X3_pr_i_{1+0},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+1}(X3_pr_i_{1+1},:) , zeros(nrows_A_(1+1),8) , repmat(T_u_{1+1}(X3_pr_i_{1+1}),1,8) , zeros(nrows_A_(1+1),8) , An{1+1}(X3_pr_i_{1+1},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+2}(X3_pr_i_{1+2},:) , zeros(nrows_A_(1+2),8) , repmat(T_u_{1+2}(X3_pr_i_{1+2}),1,8) , zeros(nrows_A_(1+2),8) , An{1+2}(X3_pr_i_{1+2},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X3'));
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
% X4;
tmp_c = X4_pc_i;
figure; 
imagesc([ ...
	  Tn{1+0}(X4_pr_i_{1+0},:) , zeros(nrows_A_(1+0),8) , repmat(T_u_{1+0}(X4_pr_i_{1+0}),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(X4_pr_i_{1+0},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+1}(X4_pr_i_{1+1},:) , zeros(nrows_A_(1+1),8) , repmat(T_u_{1+1}(X4_pr_i_{1+1}),1,8) , zeros(nrows_A_(1+1),8) , An{1+1}(X4_pr_i_{1+1},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+2}(X4_pr_i_{1+2},:) , zeros(nrows_A_(1+2),8) , repmat(T_u_{1+2}(X4_pr_i_{1+2}),1,8) , zeros(nrows_A_(1+2),8) , An{1+2}(X4_pr_i_{1+2},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Sn{1+0}(X4_pz_i_{1+0},:) , zeros(nrows_Z_(1+0),8) , repmat(S_u_{1+0}(X4_pz_i_{1+0}),1,8) , zeros(nrows_Z_(1+0),8) , Zn{1+0}(X4_pz_i_{1+0},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Sn{1+1}(X4_pz_i_{1+1},:) , zeros(nrows_Z_(1+1),8) , repmat(S_u_{1+1}(X4_pz_i_{1+1}),1,8) , zeros(nrows_Z_(1+1),8) , Zn{1+1}(X4_pz_i_{1+1},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Sn{1+2}(X4_pz_i_{1+2},:) , zeros(nrows_Z_(1+2),8) , repmat(S_u_{1+2}(X4_pz_i_{1+2}),1,8) , zeros(nrows_Z_(1+2),8) , Zn{1+2}(X4_pz_i_{1+2},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X4'));
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
% X5 and X6;
tmp_c = X5_pc_i;
figure; 
imagesc([ ...
	  Tn{1+0}(X5_pr_i_{1+0},:) , zeros(nrows_A_(1+0),8) , repmat(T_u_{1+0}(X5_pr_i_{1+0}),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(X5_pr_i_{1+0},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+1}(X5_pr_i_{1+1},:) , zeros(nrows_A_(1+1),8) , repmat(T_u_{1+1}(X5_pr_i_{1+1}),1,8) , zeros(nrows_A_(1+1),8) , An{1+1}(X5_pr_i_{1+1},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+2}(X5_pr_i_{1+2},:) , zeros(nrows_A_(1+2),8) , repmat(T_u_{1+2}(X5_pr_i_{1+2}),1,8) , zeros(nrows_A_(1+2),8) , An{1+2}(X5_pr_i_{1+2},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X5'));
tmp_c = X6_pc_i;
figure; 
imagesc([ ...
	  Tn{1+0}(X6_pr_i_{1+0},:) , zeros(nrows_A_(1+0),8) , repmat(T_u_{1+0}(X6_pr_i_{1+0}),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(X6_pr_i_{1+0},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+1}(X6_pr_i_{1+1},:) , zeros(nrows_A_(1+1),8) , repmat(T_u_{1+1}(X6_pr_i_{1+1}),1,8) , zeros(nrows_A_(1+1),8) , An{1+1}(X6_pr_i_{1+1},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+2}(X6_pr_i_{1+2},:) , zeros(nrows_A_(1+2),8) , repmat(T_u_{1+2}(X6_pr_i_{1+2}),1,8) , zeros(nrows_A_(1+2),8) , An{1+2}(X6_pr_i_{1+2},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X6'));
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
% X0;
tmp_c = X0_pc_i;
figure; 
imagesc([ ...
	  Tn{1+0}(X0_pr_i_{1+0},:) , zeros(nrows_A_(1+0),8) , repmat(T_u_{1+0}(X0_pr_i_{1+0}),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(X0_pr_i_{1+0},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+1}(X0_pr_i_{1+1},:) , zeros(nrows_A_(1+1),8) , repmat(T_u_{1+1}(X0_pr_i_{1+1}),1,8) , zeros(nrows_A_(1+1),8) , An{1+1}(X0_pr_i_{1+1},tmp_c) ; ...
	    zeros(8,ncols_T) , zeros(8,8) , zeros(8,8) , zeros(8,8) , zeros(8,ncols_A) ; ...
	  Tn{1+2}(X0_pr_i_{1+2},:) , zeros(nrows_A_(1+2),8) , repmat(T_u_{1+2}(X0_pr_i_{1+2}),1,8) , zeros(nrows_A_(1+2),8) , An{1+2}(X0_pr_i_{1+2},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X0'));
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
% XX;
figure; 
imagesc([ An{1}(XW_pr_i_{1},XW_pc_i) ; ...
	  ] , [0,1]);
title(sprintf('XW'));
return;
end;%if disp_flag;

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
tmpAtchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_A{1+nb1}(:)>0,tmpAtchar); %disp(sprintf('mr_A_[%d]: %s',nb1,num2str(mr_A{1+nb1})));
tmpZtchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_Z{1+nb1}(:)>0,tmpZtchar); %disp(sprintf('mr_Z_[%d]: %s',nb1,num2str(mr_Z{1+nb1})));
end;%for nb1=0:nbins-1;

fname__in = sprintf('%s/%s.in',dir__in,in__name);
fp = fopen(fname__in,'w');
fprintf(fp,'GLOBAL_verbose= 0;\n');
fprintf(fp,'GLOBAL_thread_count= 8;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','lakcluster_driver');
fprintf(fp,'GLOBAL_QR_strategy= %s;\n',QR_strategy);
fprintf(fp,'GLOBAL_QC_strategy= %s;\n',QC_strategy);
fprintf(fp,'GLOBAL_NBINS= %d;\n',nbins);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',B_MLT);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
%fprintf(fp,'GLOBAL_TEST_sparse= 0;\n');
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
fprintf(fp,'GLOBAL_DIR_XPRE= %s;\n',dir__in);
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_A_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',ncols_A);
fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'mr_A_%.2d.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cind= mc_A.b16;\n');
fprintf(fp,'GLOBAL_Z_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'Z_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'Z_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_Z_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'mr_Z_%.2d.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',ncols_Y);
fprintf(fp,'GLOBAL_T_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'T_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'T_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_T_n_cind= mc_T.b16;\n');
fprintf(fp,'GLOBAL_S_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'S_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_S_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'S_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lakcluster_test_uADZSZDA_Ireq2_ver7.m on %s;\n',date);
fclose(fp);

call_flag=1;%call_flag = input(' call? 1=yes (default), 0=no:'); if isempty(call_flag); call_flag=1; end;
if call_flag;
disp(sprintf('%s/../lakcluster_ver18 < %s/%s.in',dir_trunk,dir__in,in__name));
system(sprintf('%s/../lakcluster_ver18 < %s/%s.in',dir_trunk,dir__in,in__name));
end;%if call_flag;

run(sprintf('%s/timing',dir_out)); %run(sprintf('%s/dims_X',dir_out));
tmpAtchar = sprintf('%s/mc_A.b16',dir__in);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAtchar); mc_A_tmp = tutorial_binary_uncompress(tmpAtchar,1:nrows,1:ncols)>0; 
for nb1=0:nbins-1; 
tmpAnchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAnchar); mr_A_tmp{1+nb1} = tutorial_binary_uncompress(tmpAnchar,1:nrows,1:ncols)>0; 
tmpZnchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpZnchar); mr_Z_tmp{1+nb1} = tutorial_binary_uncompress(tmpZnchar,1:nrows,1:ncols)>0; 
end;%for nb1=0:nbins-1; 

out_xdrop_a = textread(sprintf('%s/out_xdrop_a.txt',dir_out)); out_xdrop_b = textread(sprintf('%s/out_xdrop_b.txt',dir_out));
r_ij = out_xdrop_a(find(out_xdrop_a(:,1)>-1),1); c_ij = out_xdrop_a(find(out_xdrop_a(:,2)>-1),2);
[r_ij_ord,r_ij_i] = sort(1+r_ij); [c_ij_ord,c_ij_i] = sort(1+c_ij);
clear r_on_ z_on_ c_on y_on ;
for nb1=0:nbins-1; r_on_{1+nb1} = find(mr_A_tmp{1+nb1}); z_on_{1+nb1} = find(mr_Z_tmp{1+nb1}); end;%for nb1=0:nbins-1;
c_on = find(mc_A_tmp); 
r_on_x = r_on_{1}; for nb1=1:nbins-1; r_on_x = unionall({r_on_x,r_on_{1+nb1}+sum(nrows_A_(1:1+nb1-1))}); end;
z_on_x = z_on_{1}; for nb1=1:nbins-1; z_on_x = unionall({z_on_x,z_on_{1+nb1}+sum(nrows_Z_(1:1+nb1-1))}); end;

XW_pc_h = intersect(c_on,c_ij_i(XW_pc(1:ncols_X))); AW_pc_h = setdiff(c_on,XW_pc_h); 
auc_c_W = get_auc(AW_pc_h,XW_pc_h);

X0_pc_h = intersect(c_on,c_ij_i(X0_pc(1:ncols_X))); A0_pc_h = setdiff(c_on,X0_pc_h); 
auc_c_0 = get_auc(A0_pc_h,X0_pc_h);
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

nrows_A_sum = cumsum([0,nrows_A_]);
for nb1=0:nbins-1;
XW_pr_g_{1+nb1} = XW_pr_{1+nb1} + nrows_A_sum(1+nb1);
X0_pr_g_{1+nb1} = X0_pr_{1+nb1} + nrows_A_sum(1+nb1);
X1_pr_g_{1+nb1} = X1_pr_{1+nb1} + nrows_A_sum(1+nb1);
X2_pr_g_{1+nb1} = X2_pr_{1+nb1} + nrows_A_sum(1+nb1);
X3_pr_g_{1+nb1} = X3_pr_{1+nb1} + nrows_A_sum(1+nb1);
X4_pr_g_{1+nb1} = X4_pr_{1+nb1} + nrows_A_sum(1+nb1);
X5_pr_g_{1+nb1} = X5_pr_{1+nb1} + nrows_A_sum(1+nb1);
X6_pr_g_{1+nb1} = X6_pr_{1+nb1} + nrows_A_sum(1+nb1);
end;%for nb1=0:nbins-1;

% XW;
for nb1=0:nbins-1;
nb1=1-1; XW_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(XW_pr_g_{1+nb1}(1:floor(nrows_X/3)))); AW_pr_h_{1+nb1} = setdiff(r_on_x,XW_pr_h_{1+nb1});
auc_r_W(1+nb1) = get_auc(AW_pr_h_{1+nb1},XW_pr_h_{1+nb1});
end;%for nb1=0:nbins-1;
auc_r_W_x = mean(auc_r_W);

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
%auc_x_X = (auc_x_1 + auc_x_2 + auc_x_3)/3;
%auc_x_Y = (auc_x_5 + auc_x_6)/2;

xx_ra = [ ...
          auc_x_0   , auc_x_1   , auc_x_2   , auc_x_3   , auc_x_4   , auc_x_5   , auc_x_6   ; ...
          auc_r_0_x , auc_r_1_x , auc_r_2_x , auc_r_3_x , auc_r_4_x , auc_r_5_x , auc_r_6_x ; ...
          auc_c_0   , auc_c_1   , auc_c_2   , auc_c_3   , auc_c_4   , auc_c_5   , auc_c_6   ; ...
	  ];
disp(num2str(xx_ra));
xx_ra = xx_ra(:);

if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

clean_flag=0;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;

