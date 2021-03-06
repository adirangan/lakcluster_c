function lakcluster_test_uAZZA_ver4(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num)
% This case considers four clusters embedded within a data-set possessing cases and controls, ;
% as well as two categorical covariate categories. ;
% The first cluster straddles both covariate categories within the cases, ;
% whereas the next two clusters are each restricted to a single category within the cases. ;
% The fourth cluster is even larger and straddles both covariates, but also straddles cases and controls. ;
% Here we make nrows_Z_(nbins) twice as large as nrows_A_(nbins) (i.e., twice as many controls as cases). ;
% uses lakcluster_ver16;
% try: ;
%{
  
  N=1024;X_factor=0.55;X_esm=0.1;gamma=0.05;B_MLT=32;rng_num=1;lakcluster_test_uAZZA_ver4(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

  nbins=2; % for uAZZA;
  N_ra = 2.^[6:12]; X_factor_ra = [0.45:0.025:0.65]; X_esm_ra = 10.^[-3:0.5:1]; gamma_ra = [0,0.5,0.15,0.5]; B_ra = [0,4,8,16]; rng_num_ra = 1:4;
  for nN=1:length(N_ra); N = N_ra(nN);
  for nf=1:length(X_factor_ra); X_factor = X_factor_ra(nf);
  for ne=1:length(X_esm_ra); X_esm = X_esm_ra(ne);
  for ng=1:length(gamma_ra); gamma = gamma_ra(ng);
  for nB=1:length(B_ra); B_MLT = B_ra(nB);
  for nr=1:length(rng_num_ra); rng_num = rng_num_ra(nr);
  lakcluster_test_uAZZA_ver4(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);
  xfix_del(pwd,'lakcluster_test_uAZZA',N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
  end;end;end;end;end;end;%for N,f,esm,B,gamma,rng;

  %}


if (nargin<8); pt_num = -1; end;
if (nargin<7); rng_num = 1; end;
if (nargin<4); B_MLT = 32; end;
if (nargin<5); gamma = 0.025; end;
if (nargin<4); X_esm = 0.1; end;
if (nargin<3); X_factor = 0.65; end;
if (nargin<2); N = 1024; end;
if (nargin<1); dir_trunk = pwd; end;

nbins = 2;
nrows_A_ = ceil([N/2 , N/2+1]);
nrows_Z_ = ceil([N/1 , N/1+1]);
ncols_A = N+1;
ncols_Y = 0;
ncols_T = 1;
mrnd = 0;
QR_strategy = 'YnWt';
QC_strategy = 'YnWt';
verbose_flag=0;

nrows_X = ceil(N.^(X_factor));
ncols_X = ceil((N+1).^(X_factor));
X_eps = X_esm/sqrt(nrows_X); % spectral error for X;
X_pp = tutorial_g1em(X_esm);

test_string = 'lakcluster_test_uAZZA';
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

X_1 = tutorial_binary_makelr2(1*nrows_X,1*ncols_X,X_pp);
nrows_X_11 = ceil(nrows_X/2); 
nrows_X_12 = nrows_X - nrows_X_11;
X_11 = X_1(1:nrows_X_11,:);
X_12 = X_1(nrows_X_11 + (1:nrows_X_12),:);
X_2 = tutorial_binary_makelr2(1*nrows_X,1*ncols_X,X_pp);
X_3 = tutorial_binary_makelr2(1*nrows_X,1*ncols_X,X_pp);
X_4 = tutorial_binary_makelr2(3*nrows_X,2*ncols_X,X_pp);
nrows_X_41 = ceil(nrows_X/2); 
nrows_X_42 = ceil(nrows_X/2);
nrows_X_43 = nrows_X;
nrows_X_44 = 3*nrows_X - nrows_X_41 - nrows_X_42 - nrows_X_43;
X_41 = X_4(1:nrows_X_41,:);
X_42 = X_4(nrows_X_41 + (1:nrows_X_42),:);
X_43 = X_4(nrows_X_41 + nrows_X_42 + (1:nrows_X_43),:);
X_44 = X_4(nrows_X_41 + nrows_X_42 + nrows_X_43 + (1:nrows_X_44),:);

for nb1=0:nbins-1; 
An{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_A)>0.5)-1; 
Zn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_A)>0.5)-1; 
Yn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_Y)>0.5)-1; 
Wn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_Y)>0.5)-1; 
Tn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_T)>0.5)-1; Tn{1+nb1}(:,1) = 1; 
Sn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_T)>0.5)-1; Sn{1+nb1}(:,1) = 1; 
end;%for nb1=0:nbins-1;

%insert clusters;
pr11_f = randperm(nrows_A_(1+0)); [tmp,pr11_i] = sort(pr11_f); [tmp,pr11_i] = sort(pr11_i);
pr12_f = randperm(nrows_A_(1+1)); [tmp,pr12_i] = sort(pr12_f); [tmp,pr12_i] = sort(pr12_i);
pc1x_f = randperm(ncols_A);       [tmp,pc1x_i] = sort(pc1x_f); [tmp,pc1x_i] = sort(pc1x_i);
pr2__f = randperm(nrows_A_(1+0)); [tmp,pr2__i] = sort(pr2__f); [tmp,pr2__i] = sort(pr2__i);
pc2__f = randperm(ncols_A);       [tmp,pc2__i] = sort(pc2__f); [tmp,pc2__i] = sort(pc2__i);
pr3__f = randperm(nrows_A_(1+1)); [tmp,pr3__i] = sort(pr3__f); [tmp,pr3__i] = sort(pr3__i);
pc3__f = randperm(ncols_A);       [tmp,pc3__i] = sort(pc3__f); [tmp,pc3__i] = sort(pc3__i);
pr41_f = randperm(nrows_A_(1+0)); [tmp,pr41_i] = sort(pr41_f); [tmp,pr41_i] = sort(pr41_i);
pr42_f = randperm(nrows_A_(1+1)); [tmp,pr42_i] = sort(pr42_f); [tmp,pr42_i] = sort(pr42_i);
pz41_f = randperm(nrows_Z_(1+0)); [tmp,pz41_i] = sort(pz41_f); [tmp,pz41_i] = sort(pz41_i);
pz42_f = randperm(nrows_Z_(1+1)); [tmp,pz42_i] = sort(pz42_f); [tmp,pz42_i] = sort(pz42_i);
pc4x_f = randperm(ncols_A);       [tmp,pc4x_i] = sort(pc4x_f); [tmp,pc4x_i] = sort(pc4x_i);
An{1+0}(pr41_f(1:nrows_X_41),pc4x_f(1:2*ncols_X)) = X_41;
An{1+1}(pr42_f(1:nrows_X_42),pc4x_f(1:2*ncols_X)) = X_42;
Zn{1+0}(pz41_f(1:nrows_X_43),pc4x_f(1:2*ncols_X)) = X_43;
Zn{1+1}(pz42_f(1:nrows_X_44),pc4x_f(1:2*ncols_X)) = X_44;
An{1+1}(pr3__f(1:nrows_X),pc3__f(1:ncols_X)) = X_3;
An{1+0}(pr2__f(1:nrows_X),pc2__f(1:ncols_X)) = X_2;
An{1+0}(pr11_f(1:nrows_X_11),pc1x_f(1:ncols_X)) = X_11;
An{1+1}(pr12_f(1:nrows_X_12),pc1x_f(1:ncols_X)) = X_12;

disp_flag=0;
if disp_flag;
figure;imagesc([ An{1+0}(pr11_i,pc1x_i) ; zeros(8,ncols_A) ; An{1+1}(pr12_i,pc1x_i) ; zeros(8,ncols_A) ; Zn{1+0}(:,pc1x_i) ; zeros(8,ncols_A) ; Zn{1+1}(:,pc1x_i) ]); title('X_1');
figure;imagesc([ An{1+0}(pr2__i,pc2__i) ; zeros(8,ncols_A) ; An{1+1}(:,pc2__i) ; zeros(8,ncols_A) ; Zn{1+0}(:,pc2__i) ; zeros(8,ncols_A) ; Zn{1+1}(:,pc2__i) ]); title('X_2');
figure;imagesc([ An{1+0}(:,pc3__i) ; zeros(8,ncols_A) ; An{1+1}(pr3__i,pc3__i) ; zeros(8,ncols_A) ; Zn{1+0}(:,pc3__i) ; zeros(8,ncols_A) ; Zn{1+1}(:,pc3__i) ]); title('X_3');
figure;imagesc([ An{1+0}(pr41_i,pc4x_i) ; zeros(8,ncols_A) ; An{1+1}(pr42_i,pc4x_i) ; zeros(8,ncols_A) ; Zn{1+0}(pz41_i,pc4x_i) ; zeros(8,ncols_A) ; Zn{1+1}(pz42_i,pc4x_i) ]); title('X_4');
return;
end;%if disp_flag;

for nb1=0:nbins-1; 
At{1+nb1} = transpose(An{1+nb1}); Zt{1+nb1} = transpose(Zn{1+nb1});
Yt{1+nb1} = transpose(Yn{1+nb1}); Wt{1+nb1} = transpose(Wn{1+nb1});
Tt{1+nb1} = transpose(Tn{1+nb1}); St{1+nb1} = transpose(Sn{1+nb1});
end;%for nb1=0:nbins-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpAtchar = sprintf('%s/mc_A.b16',dir__in);tutorial_binary_compress(bitj,mc_A(:)>0,tmpAtchar); %disp(sprintf('mc_A: %s',num2str(mc_A)));
for nb1=0:nbins-1;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,An{1+nb1}>0,tmpAnchar); 
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,At{1+nb1}>0,tmpAtchar); 
tmpZnchar = sprintf('%s/Z_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Zn{1+nb1}>0,tmpZnchar); 
tmpZtchar = sprintf('%s/Z_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Zt{1+nb1}>0,tmpZtchar); 
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
fprintf(fp,'GLOBAL_Ireq= 0;\n');
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
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lakcluster_test_uAZZA_ver4.m on %s;\n',date);
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
r_on_x = unionall({r_on_{1},r_on_{2}+nrows_A_(1)});
z_on_x = unionall({z_on_{1},z_on_{2}+nrows_Z_(1)});

X_pc1x = intersect(c_on,c_ij_i(pc1x_f(1:ncols_X))); A_pc1x = setdiff(c_on,X_pc1x);
X_pc2_ = intersect(c_on,c_ij_i(pc2__f(1:ncols_X))); A_pc2_ = setdiff(c_on,X_pc2_);
X_pc3_ = intersect(c_on,c_ij_i(pc3__f(1:ncols_X))); A_pc3_ = setdiff(c_on,X_pc3_);
X_pc4x = intersect(c_on,c_ij_i(pc4x_f(1:2*ncols_X))); A_pc4x = setdiff(c_on,X_pc4x);
pr11_g = pr11_f;
pr12_g = pr12_f + nrows_A_(1);
pr2__g = pr2__f;
pr3__g = pr3__f + nrows_A_(1);
pr41_g = pr41_f;
pr42_g = pr42_f + nrows_A_(1);

X_pr11 = intersect(r_on_x,r_ij_i(pr11_g(1:nrows_X_11))); A_pr11 = setdiff(r_on_x,X_pr11);
X_pr12 = intersect(r_on_x,r_ij_i(pr12_g(1:nrows_X_12))); A_pr12 = setdiff(r_on_x,X_pr12);
X_pr2_ = intersect(r_on_x,r_ij_i(pr2__g(1:nrows_X))); A_pr2_ = setdiff(r_on_x,X_pr2_);
X_pr3_ = intersect(r_on_x,r_ij_i(pr3__g(1:nrows_X))); A_pr3_ = setdiff(r_on_x,X_pr3_);
X_pr41 = intersect(r_on_x,r_ij_i(pr41_g(1:nrows_X_41))); A_pr41 = setdiff(r_on_x,X_pr41);
X_pr42 = intersect(r_on_x,r_ij_i(pr42_g(1:nrows_X_42))); A_pr42 = setdiff(r_on_x,X_pr42);

auc_c_1x = get_auc((A_pc1x),(X_pc1x));
auc_r_11 = get_auc((A_pr11),(X_pr11)); 
auc_r_12 = get_auc((A_pr12),(X_pr12)); 
auc_c_2 = get_auc((A_pc2_),(X_pc2_));
auc_r_2 = get_auc((A_pr2_),(X_pr2_));
auc_c_3 = get_auc((A_pc3_),(X_pc3_));
auc_r_3 = get_auc((A_pr3_),(X_pr3_));
auc_c_4x = get_auc((A_pc4x),(X_pc4x));
auc_r_41 = get_auc((A_pr41),(X_pr41)); 
auc_r_42 = get_auc((A_pr42),(X_pr42)); 
xx_ra = [auc_c_1x,auc_r_11,auc_r_12,auc_c_2,auc_r_2,auc_c_3,auc_r_3,auc_c_4x,auc_r_41,auc_r_42];
disp(num2str(xx_ra));

if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

clean_flag=1;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;

