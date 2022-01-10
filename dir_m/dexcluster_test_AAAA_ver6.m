function dexcluster_test_AAAA_ver6(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num)
% trying to test dexcluster with single bins and mrnd=0 ; 
% this time we embed only a single balanced cluster. no controls! ;
% try: ;
%{
  
  N=4024*1;X_factor=0.55;X_esm=0.1;gamma=0.05;B_MLT= 32;rng_num=1;dexcluster_test_AAAA_ver6(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

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
nrows_Z_ = 0 ; 
ncols_A = 8*N+1 ;
nmds = 0; % 010917;
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
Ireq = 0; % 010917 ;

nrows_X = ceil(N.^(X_factor)); nrows_XX = ceil(1.5*nrows_X);
ncols_X = ceil((N+1).^(X_factor)); ncols_XX = ceil(1.5*ncols_X);
X_eps = X_esm/sqrt(nrows_X); % spectral error for X;
X_pp = tutorial_g1em(X_esm); disp(sprintf(' %% X_esm %0.4f, log10(X_esm) %0.4f; X_pp %0.4f, log10(X_pp) %0.4f',X_esm,log10(X_esm),X_pp,log10(X_pp)));

test_string = 'dexcluster_test_AAAA';
dir_base = sprintf('%s/dir_%s',dir_trunk,test_string); if ~exist(dir_base,'dir'); mkdir(dir_base); end;
[posfix_use,X_factor_d,X_esm_d,gamma_d] = xfix_gen(N,X_factor,X_esm,gamma,B_MLT,rng_num);
in__name = sprintf('%s_%s__in',test_string,posfix_use);
out_name = sprintf('%s_%s_out',test_string,posfix_use);
dir__in = sprintf('%s/dir_%s',dir_base,in__name); if ~exist(dir__in,'dir'); mkdir(dir__in); end;
dir_out = sprintf('%s/dir_%s',dir_base,out_name); if ~exist(dir_out,'dir'); mkdir(dir_out); end; 
%path(path,dir__in); path(path,dir_out);
rng(rng_num);

clear An At Zn Zt Tn Tt Sn St mr_A mr_Z mc_A mc_T bitj;
bitj = 16;

mc_A = (rand(1,ncols_A)>mrnd);
mc_T = (rand(1,ncols_T)>mrnd); mc_T(1)=1; % ensure first column of [T;S] is all ones ;
for nb1=0:nbins-1; 
mr_A{1+nb1} = (rand(nrows_A_(1+nb1),1)>mrnd);
mr_Z{1+nb1} = (rand(nrows_Z_(1+nb1),1)>mrnd);
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

[tmp,X0_pc_i] = sort(X0_pc); [tmp,X0_pc_i] = sort(X0_pc_i);

%insert clusters;

% insert X0;
for nb1=0:nbins-1; An{1+nb1}(X0_pr_{1+nb1}(0*floor(nrows_X/nbins) + (1:floor(nrows_X/nbins))),X0_pc(1:ncols_X)) = X0(nb1*floor(nrows_X/nbins) + (1:floor(nrows_X/nbins)),1:ncols_X); end;

% generate transposes ;
for nb1=0:nbins-1;
At{1+nb1} = transpose(An{1+nb1}); Zt{1+nb1} = transpose(Zn{1+nb1});
Tt{1+nb1} = transpose(Tn{1+nb1}); St{1+nb1} = transpose(Sn{1+nb1});
end;%for nb1=0:nbins-1;


disp_flag=0;
if disp_flag;
% X0;
nb1=0;
tmp_c = X0_pc_i;
tmp_r = X0_pr_i_{1+nb1};
imagesc([ ...
	  An{1+nb1}(X0_pr_i_{1+nb1},tmp_c) ; ...
	  ] , [0,1]);
title(sprintf('X0, nb1 %d',nb1));
return;
end;%if disp_flag;

for nb1=0:nbins-1; 
At{1+nb1} = transpose(An{1+nb1}); Zt{1+nb1} = transpose(Zn{1+nb1});
Tt{1+nb1} = transpose(Tn{1+nb1}); St{1+nb1} = transpose(Sn{1+nb1});
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
tmpAtchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_A{1+nb1}(:)>0,tmpAtchar); %disp(sprintf('mr_A_[%d]: %s',nb1,num2str(mr_A{1+nb1})));
tmpZtchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_Z{1+nb1}(:)>0,tmpZtchar); %disp(sprintf('mr_Z_[%d]: %s',nb1,num2str(mr_Z{1+nb1})));
end;%for nb1=0:nbins-1;

fname__in = sprintf('%s/%s.in',dir__in,in__name);
fp = fopen(fname__in,'w');
fprintf(fp,'GLOBAL_verbose= 0;\n');
fprintf(fp,'GLOBAL_thread_count= 8;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','dexcluster_driver');
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
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_Z_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by dexcluster_test_AAAA_ver6.m on %s;\n',date);
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

nrows_A_sum = cumsum([0,nrows_A_]);
for nb1=0:nbins-1;
X0_pr_g_{1+nb1} = X0_pr_{1+nb1} + nrows_A_sum(1+nb1);
end;%for nb1=0:nbins-1;

% X0;
nb1=1-1; X0_pr_h_{1+nb1} = intersect(r_on_x,r_ij_i(X0_pr_g_{1+nb1}(1:floor(nrows_X/nbins)))); A0_pr_h_{1+nb1} = setdiff(r_on_x,X0_pr_h_{1+nb1});
auc_r_0(1+nb1) = get_auc(A0_pr_h_{1+nb1},X0_pr_h_{1+nb1});
auc_r_0_x = auc_r_0(1);

auc_x_0 = 0.5*(auc_r_0_x + auc_c_0) ;
xx_ra = [ auc_x_0 ; auc_r_0_x ; auc_c_0 ];
disp(num2str(xx_ra));

xx_ra = xx_ra(:);

plot_flag=1;
if plot_flag;
figure;
nb1 = 1;
subplot(1,2,1); imagesc(An{nb1}(X0_pr_i_{nb1},X0_pc_i),[-1,+1]);
subplot(1,2,2); imagesc(An{nb1}(1+r_ij(end:-1:1),1+c_ij(end:-1:1)),[-1,+1]);
end;%if plot_flag;

if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

clean_flag=1;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;
