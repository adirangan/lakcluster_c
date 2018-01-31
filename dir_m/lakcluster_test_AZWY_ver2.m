function lakcluster_test_AZWY_ver2(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num)
% starting with basic cluster X1n embedded in An{nbins}, and another disjoint basic cluster X2n embedded across An{nbins} and Zn{nbins} : ;
% As a nod to the real world, we make nrows_Z_(nbins) twice as large as nrows_A_(nbins) (i.e., more controls than cases). ;
% try: ;
%{
  
  N=512;X_factor=0.65;X_esm=0.1;gamma=0.05;B_MLT=0;rng_num=1;lakcluster_test_AZWY_ver2(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

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

nbins = 1;
nrows_A_ = ceil(N/nbins*ones(1,nbins));
nrows_Z_ = ceil(2*N/nbins*ones(1,nbins));
ncols_A = N+1;
ncols_Y = N+2;
nmds = 0; 
ncols_T = 1;
mrnd = 0;
QR_strategy = 'YnWt';
QC_strategy = 'ZtSWn';

nrows_X = ceil(N.^(X_factor)); nrows_XX = ceil(1.5*nrows_X);
ncols_X = ceil((ncols_A).^(X_factor)); ncols_XX = ceil(1.5*ncols_X);
X_eps = X_esm/sqrt(nrows_X); % spectral error for X;
X_pp = tutorial_g1em(X_esm);

test_string = 'lakcluster_test_AZWY';
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
for nb1=0:nbins-1; 
An{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_A)<0.5)-1;
Yn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_Y)<0.5)-1;
Zn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_A)<0.5)-1;
Wn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_Y)<0.5)-1;
end;%for nb1=0:nbins-1; 

% generate clusters ;
% X0 is exclusive to A ;
X0 = tutorial_binary_makelr2(nrows_X,ncols_X,X_pp);
for nb1=0:nbins-1; X0_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X0_pr_i_{1+nb1}] = sort(X0_pr_{1+nb1}); [tmp,X0_pr_i_{1+nb1}] = sort(X0_pr_i_{1+nb1}); end;
X0_pc = randperm(ncols_A); [tmp,X0_pc_i] = sort(X0_pc); [tmp,X0_pc_i] = sort(X0_pc_i);
% X1 spans A and Z ; 
X1 = tutorial_binary_makelr2(3*nrows_XX,ncols_XX,X_pp);
for nb1=0:nbins-1; X1_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X1_pr_i_{1+nb1}] = sort(X1_pr_{1+nb1}); [tmp,X1_pr_i_{1+nb1}] = sort(X1_pr_i_{1+nb1}); end;
for nb1=0:nbins-1; X1_pz_{1+nb1} = randperm(nrows_Z_(1+nb1)); [tmp,X1_pz_i_{1+nb1}] = sort(X1_pz_{1+nb1}); [tmp,X1_pz_i_{1+nb1}] = sort(X1_pz_i_{1+nb1}); end;
X1_pc = randperm(ncols_A); [tmp,X1_pc_i] = sort(X1_pc); [tmp,X1_pc_i] = sort(X1_pc_i);
% X2 spans A and Y ; 
X2 = tutorial_binary_makelr2(nrows_XX,2*ncols_XX,X_pp);
for nb1=0:nbins-1; X2_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X2_pr_i_{1+nb1}] = sort(X2_pr_{1+nb1}); [tmp,X2_pr_i_{1+nb1}] = sort(X2_pr_i_{1+nb1}); end;
X2_pc = randperm(ncols_A); [tmp,X2_pc_i] = sort(X2_pc); [tmp,X2_pc_i] = sort(X2_pc_i);
X2_py = randperm(ncols_Y); [tmp,X2_py_i] = sort(X2_py); [tmp,X2_py_i] = sort(X2_py_i);
% X3 spans A, Z, W and Y ; 
X3 = tutorial_binary_makelr2(3*nrows_XX,2*ncols_XX,X_pp);
for nb1=0:nbins-1; X3_pr_{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,X3_pr_i_{1+nb1}] = sort(X3_pr_{1+nb1}); [tmp,X3_pr_i_{1+nb1}] = sort(X3_pr_i_{1+nb1}); end;
for nb1=0:nbins-1; X3_pz_{1+nb1} = randperm(nrows_Z_(1+nb1)); [tmp,X3_pz_i_{1+nb1}] = sort(X3_pz_{1+nb1}); [tmp,X3_pz_i_{1+nb1}] = sort(X3_pz_i_{1+nb1}); end;
X3_pc = randperm(ncols_A); [tmp,X3_pc_i] = sort(X3_pc); [tmp,X3_pc_i] = sort(X3_pc_i);
X3_py = randperm(ncols_Y); [tmp,X3_py_i] = sort(X3_py); [tmp,X3_py_i] = sort(X3_py_i);

%insert clusters;
% insert X0;
nb1=1-1; An{1+nb1}(X0_pr_{1+nb1}(0*floor(nrows_X) + (1:floor(nrows_X))),X0_pc(1:ncols_X)) = X0(0*floor(nrows_X) + (1:floor(nrows_X)),1:ncols_X);
% insert X1 ;
nb1=1-1; An{1+nb1}(X1_pr_{1+nb1}(0*floor(nrows_XX) + (1:floor(nrows_XX))),X1_pc(1:ncols_XX)) = X1(0*floor(nrows_XX) + (1:floor(nrows_XX)),1:ncols_XX);
nb1=1-1; Zn{1+nb1}(X1_pz_{1+nb1}(0*floor(nrows_XX*2) + (1:floor(nrows_XX*2))),X1_pc(1:ncols_XX)) = X1(1*floor(nrows_XX) + (1:floor(nrows_XX*2)),1:ncols_XX);
% insert X2 ;
nb1=1-1; An{1+nb1}(X2_pr_{1+nb1}(0*floor(nrows_XX) + (1:floor(nrows_XX))),X2_pc(1:ncols_XX)) = X2(0*floor(nrows_XX) + (1:floor(nrows_XX)),0*ncols_XX + (1:ncols_XX));
nb1=1-1; Yn{1+nb1}(X2_pr_{1+nb1}(0*floor(nrows_XX) + (1:floor(nrows_XX))),X2_py(1:ncols_XX)) = X2(0*floor(nrows_XX) + (1:floor(nrows_XX)),1*ncols_XX + (1:ncols_XX));
% insert X3 ;
nb1=1-1; An{1+nb1}(X3_pr_{1+nb1}(0*floor(nrows_XX) + (1:floor(nrows_XX))),X3_pc(1:ncols_XX)) = X3(0*floor(nrows_XX) + (1:floor(nrows_XX)),0*ncols_XX + (1:ncols_XX));
nb1=1-1; Yn{1+nb1}(X3_pr_{1+nb1}(0*floor(nrows_XX) + (1:floor(nrows_XX))),X3_py(1:ncols_XX)) = X3(0*floor(nrows_XX) + (1:floor(nrows_XX)),1*ncols_XX + (1:ncols_XX));
nb1=1-1; Zn{1+nb1}(X3_pz_{1+nb1}(0*floor(nrows_XX) + (1:floor(nrows_XX*2))),X3_pc(1:ncols_XX)) = X3(1*floor(nrows_XX) + (1:floor(nrows_XX*2)),0*ncols_XX + (1:ncols_XX));
nb1=1-1; Wn{1+nb1}(X3_pz_{1+nb1}(0*floor(nrows_XX) + (1:floor(nrows_XX*2))),X3_py(1:ncols_XX)) = X3(1*floor(nrows_XX) + (1:floor(nrows_XX*2)),1*ncols_XX + (1:ncols_XX));

% generate transposes ;
for nb1=0:nbins-1;
At{1+nb1} = transpose(An{1+nb1}); Zt{1+nb1} = transpose(Zn{1+nb1});
Yt{1+nb1} = transpose(Yn{1+nb1}); Wt{1+nb1} = transpose(Wn{1+nb1});
end;%for nb1=0:nbins-1;

disp_flag=0;
if disp_flag;
% X3 ;
figure; pgap = 8;
imagesc([ ...
	An{1+0}(X3_pr_i_{1+0},X3_pc_i) , zeros(nrows_A_(1+0),pgap) , Yn{1+0}(X3_pr_i_{1+0},X3_py_i) ; ...
        zeros(pgap,ncols_A) , zeros(pgap,pgap) , zeros(pgap,ncols_Y) ; ...
	Zn{1+0}(X3_pz_i_{1+0},X3_pc_i) , zeros(nrows_Z_(1+0),pgap) , Wn{1+0}(X3_pz_i_{1+0},X3_py_i) ; ...
	  ] , [0,1]);
title(sprintf('X3'));
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
% X0 ;
figure; pgap = 8;
imagesc([ ...
	An{1+0}(X0_pr_i_{1+0},X0_pc_i)  ; ...
	  ] , [0,1]);
title(sprintf('X0'));
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
% X1 ;
figure; pgap = 8;
imagesc([ ...
	An{1+0}(X1_pr_i_{1+0},X1_pc_i)  ; ...
        zeros(pgap,ncols_A)  ; ...
	Zn{1+0}(X1_pz_i_{1+0},X1_pc_i)  ; ...
	  ] , [0,1]);
title(sprintf('X1'));
return;
end;%if disp_flag;

disp_flag=0;
if disp_flag;
% X2 ;
figure; pgap = 8;
imagesc([ ...
	An{1+0}(X2_pr_i_{1+0},X2_pc_i) , zeros(nrows_A_(1+0),pgap) , Yn{1+0}(X2_pr_i_{1+0},X2_py_i) ; ...
	  ] , [0,1]);
title(sprintf('X2'));
return;
end;%if disp_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpAtchar = sprintf('%s/mc_A.b16',dir__in);tutorial_binary_compress(bitj,mc_A(:)>0,tmpAtchar); %disp(sprintf('mc_A: %s',num2str(mc_A)));
tmpYtchar = sprintf('%s/mc_Y.b16',dir__in);tutorial_binary_compress(bitj,mc_Y(:)>0,tmpYtchar); %disp(sprintf('mc_Y: %s',num2str(mc_Y)));
for nb1=0:nbins-1;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,An{1+nb1}>0,tmpAnchar); 
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,At{1+nb1}>0,tmpAtchar); 
tmpZnchar = sprintf('%s/Z_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Zn{1+nb1}>0,tmpZnchar); 
tmpZtchar = sprintf('%s/Z_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Zt{1+nb1}>0,tmpZtchar); 
tmpYnchar = sprintf('%s/Y_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Yn{1+nb1}>0,tmpYnchar); 
tmpYtchar = sprintf('%s/Y_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Yt{1+nb1}>0,tmpYtchar); 
tmpWnchar = sprintf('%s/W_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Wn{1+nb1}>0,tmpWnchar); 
tmpWtchar = sprintf('%s/W_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Wt{1+nb1}>0,tmpWtchar); 
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
fprintf(fp,'GLOBAL_Y_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'Y_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Y_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'Y_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',ncols_Y);
fprintf(fp,'GLOBAL_Y_n_cind= mc_Y.b16;\n');
fprintf(fp,'GLOBAL_W_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'W_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_W_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'W_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lakcluster_test_AZWY_ver2.m on %s;\n',date);
fclose(fp);

call_flag=1;%call_flag = input(' call? 1=yes (default), 0=no:'); if isempty(call_flag); call_flag=1; end;
if call_flag;
disp(sprintf('%s/../lakcluster_ver18 < %s/%s.in',dir_trunk,dir__in,in__name));
system(sprintf('%s/../lakcluster_ver18 < %s/%s.in',dir_trunk,dir__in,in__name));
end;%if call_flag;

run(sprintf('%s/timing',dir_out)); 
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

X0_pc_h = intersect(c_on,c_ij_i(X0_pc(1:ncols_X))); A0_pc_h = setdiff(c_on,X0_pc_h); 
auc_c_0 = get_auc(A0_pc_h,X0_pc_h);
X1_pc_h = intersect(c_on,c_ij_i(X1_pc(1:ncols_XX))); A1_pc_h = setdiff(c_on,X1_pc_h); 
auc_c_1 = get_auc(A1_pc_h,X1_pc_h);
X2_pc_h = intersect(c_on,c_ij_i(X2_pc(1:ncols_XX))); A2_pc_h = setdiff(c_on,X2_pc_h); 
auc_c_2 = get_auc(A2_pc_h,X2_pc_h);
X3_pc_h = intersect(c_on,c_ij_i(X3_pc(1:ncols_XX))); A3_pc_h = setdiff(c_on,X3_pc_h); 
auc_c_3 = get_auc(A3_pc_h,X3_pc_h);

nrows_A_sum = cumsum([0,nrows_A_]);
for nb1=0:nbins-1;
X0_pr_g_{1+nb1} = X0_pr_{1+nb1} + nrows_A_sum(1+nb1);
X1_pr_g_{1+nb1} = X1_pr_{1+nb1} + nrows_A_sum(1+nb1);
X2_pr_g_{1+nb1} = X2_pr_{1+nb1} + nrows_A_sum(1+nb1);
X3_pr_g_{1+nb1} = X3_pr_{1+nb1} + nrows_A_sum(1+nb1);
end;%for nb1=0:nbins-1;

% X0;
nb1=1-1; X0_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X0_pr_g_{1+nb1}(1:floor(nrows_X)))); A0_pr_h_{1+nb1} = setdiff(r_on_x,X0_pr_h_{1+nb1});
auc_r_0_x(1+nb1) = get_auc(A0_pr_h_{1+nb1},X0_pr_h_{1+nb1});
% X1;
nb1=1-1; X1_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X1_pr_g_{1+nb1}(1:floor(nrows_XX)))); A1_pr_h_{1+nb1} = setdiff(r_on_x,X1_pr_h_{1+nb1});
auc_r_1_x(1+nb1) = get_auc(A1_pr_h_{1+nb1},X1_pr_h_{1+nb1});
% X2;
nb1=1-1; X2_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X2_pr_g_{1+nb1}(1:floor(nrows_XX)))); A2_pr_h_{1+nb1} = setdiff(r_on_x,X2_pr_h_{1+nb1});
auc_r_2_x(1+nb1) = get_auc(A2_pr_h_{1+nb1},X2_pr_h_{1+nb1});
% X3;
nb1=1-1; X3_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X3_pr_g_{1+nb1}(1:floor(nrows_XX)))); A3_pr_h_{1+nb1} = setdiff(r_on_x,X3_pr_h_{1+nb1});
auc_r_3_x(1+nb1) = get_auc(A3_pr_h_{1+nb1},X3_pr_h_{1+nb1});

auc_x_0 = 0.5*(auc_r_0_x + auc_c_0) ;
auc_x_1 = 0.5*(auc_r_1_x + auc_c_1) ;
auc_x_2 = 0.5*(auc_r_2_x + auc_c_2) ;
auc_x_3 = 0.5*(auc_r_3_x + auc_c_3) ;

xx_ra = [ ...
          auc_x_0   , auc_x_1   , auc_x_2   , auc_x_3   ; ...
          auc_r_0_x , auc_r_1_x , auc_r_2_x , auc_r_3_x ; ...
          auc_c_0   , auc_c_1   , auc_c_2   , auc_c_3   ; ...
	  ];
disp(num2str(xx_ra));
xx_ra = xx_ra(:);

disp_flag=1;
if disp_flag;
subplot(1,2,1);imagesc(An{1+0}(X0_pr_{1+0},X0_pc),[-1,1]);
subplot(1,2,2);imagesc(An{1+0}(1+r_ij(end:-1:1),1+c_ij(end:-1:1)),[-1,1]);
end;%if disp_flag;

if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

clean_flag=1;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;
