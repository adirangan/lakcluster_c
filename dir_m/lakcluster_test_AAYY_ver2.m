function lakcluster_test_AAYY_ver2(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num)
% starting with basic cluster X1n embedded in An{nbins}, and another disjoint basic cluster X2n embedded across An{nbins} and Yn{nbins} : ;
% As a nod to the real world, we make ncols_Y==5 (i.e., only a few MDS components for each patient). ;
% try: ;
%{
  
  N=1024;X_factor=0.55;X_esm=0.1;gamma=0.05;B_MLT=32;rng_num=1;lakcluster_test_AAYY_ver2(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

  nbins=1; % for AAYY;
  N_ra = 2.^[6:12]; X_factor_ra = [0.45:0.025:0.65]; X_esm_ra = 10.^[-3:0.5:1]; gamma_ra = [0,0.5,0.15,0.5]; B_ra = [0,4,8,16]; rng_num_ra = 1:4;
  for nN=1:length(N_ra); N = N_ra(nN);
  for nf=1:length(X_factor_ra); X_factor = X_factor_ra(nf);
  for ne=1:length(X_esm_ra); X_esm = X_esm_ra(ne);
  for ng=1:length(gamma_ra); gamma = gamma_ra(ng);
  for nB=1:length(B_ra); B_MLT = B_ra(nB);
  for nr=1:length(rng_num_ra); rng_num = rng_num_ra(nr);
  lakcluster_test_AAYY_ver2(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);
  xfix_del(pwd,'lakcluster_test_AAYY',N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
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

nbins = 1;
nrows_A_ = [N];
nrows_Z_ = [0];
ncols_A = N+1;
ncols_Y = N;
ncols_T = 1;
mrnd = 0;
QR_strategy = 'YnWt';
QC_strategy = 'ZtSWn';
verbose_flag=0;

nbx = nbins;% bin to implant X;
nrows_X = ceil(nrows_A_(nbx).^(X_factor));
ncols_X = ceil(ncols_A.^(X_factor));
X_eps = X_esm/sqrt(nrows_X); % spectral error for X;
X_pp = tutorial_g1em(X_esm);

test_string = 'lakcluster_test_AAYY';
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

for nb1=0:nbins-1; 
An{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_A)>0.5)-1; 
Zn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_A)>0.5)-1; 
Yn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_Y)>0.5)-1; 
Wn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_Y)>0.5)-1; 
Tn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_T)>0.5)-1; Tn{1+nb1}(:,1) = 1; 
Sn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_T)>0.5)-1; Sn{1+nb1}(:,1) = 1; 
if (1+nb1)==nbx; % insert clusters ;
X__ = tutorial_binary_makelr2(nrows_X,ncols_X,X_pp);
XAY = tutorial_binary_makelr2(nrows_X,ncols_X + ncols_Y,X_pp);
XA_ = XAY(:,1:ncols_X);
X_Y = XAY(:,ncols_X+1:ncols_X+ncols_Y);
pr2_f{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,pr2_i{1+nb1}] = sort(pr2_f{1+nb1}); [tmp,pr2_i{1+nb1}] = sort(pr2_i{1+nb1});
pc2_f{1+nb1} = randperm(ncols_A);         [tmp,pc2_i{1+nb1}] = sort(pc2_f{1+nb1}); [tmp,pc2_i{1+nb1}] = sort(pc2_i{1+nb1});
An{1+nb1}(pr2_f{1+nb1}(1:1*nrows_X),pc2_f{1+nb1}(1:ncols_X)) = XA_;
pc3_f{1+nb1} = randperm(ncols_Y);         [tmp,pc3_i{1+nb1}] = sort(pc3_f{1+nb1}); [tmp,pc3_i{1+nb1}] = sort(pc3_i{1+nb1});
Yn{1+nb1}(pr2_f{1+nb1}(1:1*nrows_X),pc3_f{1+nb1}(1:min(ncols_X,ncols_Y))) = X_Y(:,1:min(ncols_X,ncols_Y));
pr1_f{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,pr1_i{1+nb1}] = sort(pr1_f{1+nb1}); [tmp,pr1_i{1+nb1}] = sort(pr1_i{1+nb1});
pc1_f{1+nb1} = randperm(ncols_A);         [tmp,pc1_i{1+nb1}] = sort(pc1_f{1+nb1}); [tmp,pc1_i{1+nb1}] = sort(pc1_i{1+nb1});
An{1+nb1}(pr1_f{1+nb1}(1:1*nrows_X),pc1_f{1+nb1}(1:ncols_X)) = X__;
end;%if (1+nb1)==nbins; 
At{1+nb1} = transpose(An{1+nb1}); Zt{1+nb1} = transpose(Zn{1+nb1});
Yt{1+nb1} = transpose(Yn{1+nb1}); Wt{1+nb1} = transpose(Wn{1+nb1});
Tt{1+nb1} = transpose(Tn{1+nb1}); St{1+nb1} = transpose(Sn{1+nb1});
end;%for nb1=0:nbins-1;

disp_flag=0;
if disp_flag;
figure;imagesc([ An{1+0}(pr1_i{1+0},pc1_i{1+0}) , zeros(nrows_A_(1+0),8) , Yn{1+0}(pr1_i{1+0},         :) ],[-1,1]); title('X00');
figure;imagesc([ An{1+0}(pr2_i{1+0},pc2_i{1+0}) , zeros(nrows_A_(1+0),8) , Yn{1+0}(pr2_i{1+0},pc3_i{1+0}) ],[-1,1]); title('XAY');
return;
end;%if disp_flag;

tmpAtchar = sprintf('%s/mc_A.b16',dir__in);tutorial_binary_compress(bitj,mc_A(:)>0,tmpAtchar); %disp(sprintf('mc_A: %s',num2str(mc_A)));
tmpYtchar = sprintf('%s/mc_Y.b16',dir__in);tutorial_binary_compress(bitj,mc_Y(:)>0,tmpYtchar); %disp(sprintf('mc_Y: %s',num2str(mc_Y)));
tmpTtchar = sprintf('%s/mc_T.b16',dir__in);tutorial_binary_compress(bitj,mc_T(:)>0,tmpTtchar); %disp(sprintf('mc_T: %s',num2str(mc_T)));
for nb1=0:nbins-1;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,An{1+nb1}>0,tmpAnchar); 
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,At{1+nb1}>0,tmpAtchar); 
tmpZnchar = sprintf('%s/Z_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Zn{1+nb1}>0,tmpZnchar); 
tmpZtchar = sprintf('%s/Z_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Zt{1+nb1}>0,tmpZtchar); 
tmpYnchar = sprintf('%s/Y_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Yn{1+nb1}>0,tmpYnchar); 
tmpYtchar = sprintf('%s/Y_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Yt{1+nb1}>0,tmpYtchar); 
tmpWnchar = sprintf('%s/W_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Wn{1+nb1}>0,tmpWnchar); 
tmpWtchar = sprintf('%s/W_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Wt{1+nb1}>0,tmpWtchar); 
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
fprintf(fp,'GLOBAL_DIR_XPRE= %s;\n',dir__in);
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'A_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_A_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',ncols_A);
fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'mr_A_%.2d.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cind= mc_A.b16;\n');
fprintf(fp,'GLOBAL_Y_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'Y_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Y_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'Y_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',ncols_Y);
fprintf(fp,'GLOBAL_Y_n_cind= mc_Y.b16;\n');
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lakcluster_test_AAYY_ver2.m on %s;\n',date);
fclose(fp);

call_flag=1;%call_flag = input(' call? 1=yes (default), 0=no:'); if isempty(call_flag); call_flag=1; end;
if call_flag;
disp(sprintf('%s/../lakcluster_ver18 < %s/%s.in',dir_trunk,dir__in,in__name));
system(sprintf('%s/../lakcluster_ver18 < %s/%s.in',dir_trunk,dir__in,in__name));
end;%if call_flag;

run(sprintf('%s/timing',dir_out));
tmpAtchar = sprintf('%s/mc_A.b16',dir__in);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAtchar); mc_A_tmp = tutorial_binary_uncompress(tmpAtchar,1:nrows,1:ncols)>0; 
tmpYtchar = sprintf('%s/mc_Y.b16',dir__in);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpYtchar); mc_Y_tmp = tutorial_binary_uncompress(tmpYtchar,1:nrows,1:ncols)>0; 
for nb1=0:nbins-1; 
tmpAnchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAnchar); mr_A_tmp{1+nb1} = tutorial_binary_uncompress(tmpAnchar,1:nrows,1:ncols)>0; 
tmpZnchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpZnchar); mr_Z_tmp{1+nb1} = tutorial_binary_uncompress(tmpZnchar,1:nrows,1:ncols)>0; 
end;%for nb1=0:nbins-1; 

out_xdrop_a = textread(sprintf('%s/out_xdrop_a.txt',dir_out)); out_xdrop_b = textread(sprintf('%s/out_xdrop_b.txt',dir_out));
r_ij = out_xdrop_a(find(out_xdrop_a(:,1)>-1),1); c_ij = out_xdrop_a(find(out_xdrop_a(:,2)>-1),2);
[r_ij_ord,r_ij_i] = sort(1+r_ij); [c_ij_ord,c_ij_i] = sort(1+c_ij);
nb1 = nbins-1; nb2 = nbins-1; nb_tab = nb1+nb2*nbins;
r_on = find(mr_A_tmp{1+nb1}); c_on = find(mc_A_tmp);
z_on = find(mr_Z_tmp{1+nb1}); y_on = find(mc_Y_tmp);

A_pr1_i = intersect(r_on,r_ij_i(pr1_f{1+nb1}(nrows_X+1:nrows_A_(1+nb1)))); lAr1{1+nb1} = length(A_pr1_i);
X_pr1_i = intersect(r_on,r_ij_i(pr1_f{1+nb1}(1:nrows_X))); lXr1{1+nb1} = length(X_pr1_i); 
A_pc1_i = intersect(c_on,c_ij_i(pc1_f{1+nb1}(ncols_X+1:ncols_A))); lAc1 = length(A_pc1_i);
X_pc1_i = intersect(c_on,c_ij_i(pc1_f{1+nb1}(1:ncols_X))); lXc1 = length(X_pc1_i);
A_pr2_i = intersect(r_on,r_ij_i(pr2_f{1+nb1}(nrows_X+1:nrows_A_(1+nb1)))); lAr2{1+nb1} = length(A_pr2_i);
X_pr2_i = intersect(r_on,r_ij_i(pr2_f{1+nb1}(1:nrows_X))); lXr2{1+nb1} = length(X_pr2_i); 
A_pc2_i = intersect(c_on,c_ij_i(pc2_f{1+nb1}(ncols_X+1:ncols_A))); lAc2 = length(A_pc2_i);
X_pc2_i = intersect(c_on,c_ij_i(pc2_f{1+nb1}(1:ncols_X))); lXc2 = length(X_pc2_i);

auc_r1 = get_auc(A_pr1_i,X_pr1_i); auc_c1 = get_auc(A_pc1_i,X_pc1_i);
auc_r2 = get_auc(A_pr2_i,X_pr2_i); auc_c2 = get_auc(A_pc2_i,X_pc2_i);
disp(sprintf(' %% rows: (X1/A1) = %d/%d = %0.0f%% --> auc_r1 %0.2f',lXr1{1+nb1},lAr1{1+nb1},100*lXr1{1+nb1}/lAr1{1+nb1},auc_r1));
disp(sprintf(' %% cols: (X1/A1) = %d/%d = %0.0f%% --> auc_c1 %0.2f',lXc1,lAc1,100*lXc1/lAc1,auc_c1));
disp(sprintf(' %% rows: (X2/A2) = %d/%d = %0.0f%% --> auc_r2 %0.2f',lXr2{1+nb1},lAr2{1+nb1},100*lXr2{1+nb1}/lAr2{1+nb1},auc_r2));
disp(sprintf(' %% cols: (X2/A2) = %d/%d = %0.0f%% --> auc_c2 %0.2f',lXc2,lAc2,100*lXc2/lAc2,auc_c2));

xx_ra = [auc_r1,auc_c1,0.5*(auc_r1+auc_c1),auc_r2,auc_c2,0.5*(auc_r2+auc_c2)];
if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

clean_flag=1;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;

