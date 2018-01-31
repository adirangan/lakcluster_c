function lakcluster_test_ADAADA_ver2(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num)
% starting with 9 basic clusters embedded in 3-level sparse An{nbins};
% try: ;
%{
  
  N=1024;X_factor=0.5;X_esm=0.1;gamma=0.05;B_MLT=0;rng_num=1;lakcluster_test_ADAADA_ver2(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

  nbins=1; % for ADAADA;
  N_ra = 2.^[6:12]; X_factor_ra = [0.45:0.025:0.65]; X_esm_ra = 10.^[-3:0.5:1]; gamma_ra = [0,0.5,0.15,0.5]; B_ra = [0,4,8,16]; rng_num_ra = 1:4;
  for nN=1:length(N_ra); N = N_ra(nN);
  for nf=1:length(X_factor_ra); X_factor = X_factor_ra(nf);
  for ne=1:length(X_esm_ra); X_esm = X_esm_ra(ne);
  for ng=1:length(gamma_ra); gamma = gamma_ra(ng);
  for nB=1:length(B_ra); B_MLT = B_ra(nB);
  for nr=1:length(rng_num_ra); rng_num = rng_num_ra(nr);
  lakcluster_test_ADAADA_ver2(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);
  xfix_del(pwd,'lakcluster_test_ADAADA',N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
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
ncols_A = 3*ceil(N/1920)*1920;
ncols_Y = 0;
ncols_T = 1;
mrnd = 0.0;
QR_strategy = 'YnWt';
QC_strategy = 'YnWt';
verbose_flag=0;

nbx = nbins;% bin to implant X;
nrows_X = ceil(nrows_A_(nbx).^(X_factor));
ncols_X = ceil(ncols_A.^(X_factor));
X_eps = X_esm/sqrt(nrows_X); % spectral error for X;
X_pp = tutorial_g1em(X_esm);

test_string = 'lakcluster_test_ADAADA';
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

sx_ra(1) = 0.125;
sx_ra(2) = 0.500;
sx_ra(3) = 0.875;
for nb1=0:nbins-1; 
An{1+nb1} = [2*(rand(nrows_A_(1+nb1),ncols_A/3)<sx_ra(1))-1 , 2*(rand(nrows_A_(1+nb1),ncols_A/3)<sx_ra(2))-1 , 2*(rand(nrows_A_(1+nb1),ncols_A/3)<sx_ra(3))-1];
Zn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_A)>0.5)-1; 
Yn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_Y)>0.5)-1; 
Wn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_Y)>0.5)-1; 
Tn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_T)>0.5)-1; Tn{1+nb1}(:,1) = 1; 
Sn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_T)>0.5)-1; Sn{1+nb1}(:,1) = 1; 
if (1+nb1)==nbx; % insert clusters ;
for np1=0:2;for np2=0:2;
X_{1+np1,1+np2} = tutorial_binary_makelr3(nrows_X,ncols_X,X_pp,sx_ra(1+np1));
Xrij{1+np1,1+np2} = min(nrows_A_(1+nb1)-nrows_X,np1*nrows_X + np2*3*nrows_X) + (1:nrows_X);
Xcij{1+np1,1+np2} = min(ncols_A-ncols_X,np1*ncols_X + np2*ncols_A/3) + (1:ncols_X);
An{1+nb1}(Xrij{1+np1,1+np2},Xcij{1+np1,1+np2}) = X_{1+np1,1+np2};
end;end;%for np1=0:2;for np2=0:2;
end;%if (1+nb1)==nbins; 
At{1+nb1} = transpose(An{1+nb1}); Zt{1+nb1} = transpose(Zn{1+nb1});
Yt{1+nb1} = transpose(Yn{1+nb1}); Wt{1+nb1} = transpose(Wn{1+nb1});
Tt{1+nb1} = transpose(Tn{1+nb1}); St{1+nb1} = transpose(Sn{1+nb1});
end;%for nb1=0:nbins-1;

disp_flag=0;
if disp_flag;
imagesc(An{1});
return;
end;%if disp_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

check_flag=0;%check_flag=input(' check input? 1=yes, 0=no (default)'); if isempty(check_flag); check_flag=0; end;
if (check_flag);
tmpAtchar = sprintf('%s/mc_A.b16',dir__in);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAtchar); mc_A_tmp = tutorial_binary_uncompress(tmpAtchar,1:nrows,1:ncols)>0; disp(sprintf('mc_A: %s',num2str(mc_A_tmp)));
tmpYtchar = sprintf('%s/mc_Y.b16',dir__in);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpYtchar); mc_Y_tmp = tutorial_binary_uncompress(tmpYtchar,1:nrows,1:ncols)>0; disp(sprintf('mc_Y: %s',num2str(mc_Y_tmp)));
tmpTtchar = sprintf('%s/mc_T.b16',dir__in);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpTtchar); mc_T_tmp = tutorial_binary_uncompress(tmpTtchar,1:nrows,1:ncols)>0; disp(sprintf('mc_T: %s',num2str(mc_T_tmp)));
for nb1=0:nbins-1;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAnchar); An_tmp{1+nb1} = tutorial_binary_uncompress(tmpAnchar,1:nrows,1:ncols); nrows_A_tmp_(1+nb1) = nrows; ncols_A_tmp = ncols;
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAtchar); At_tmp{1+nb1} = tutorial_binary_uncompress(tmpAtchar,1:nrows,1:ncols);
tmpZnchar = sprintf('%s/Z_%.2d_n.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpZnchar); Zn_tmp{1+nb1} = tutorial_binary_uncompress(tmpZnchar,1:nrows,1:ncols); nrows_Z_tmp_(1+nb1) = nrows;
tmpZtchar = sprintf('%s/Z_%.2d_t.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpZtchar); Zt_tmp{1+nb1} = tutorial_binary_uncompress(tmpZtchar,1:nrows,1:ncols);
tmpYnchar = sprintf('%s/Y_%.2d_n.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpYnchar); Yn_tmp{1+nb1} = tutorial_binary_uncompress(tmpYnchar,1:nrows,1:ncols); ncols_Y_tmp = ncols;
tmpYtchar = sprintf('%s/Y_%.2d_t.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpYtchar); Yt_tmp{1+nb1} = tutorial_binary_uncompress(tmpYtchar,1:nrows,1:ncols);
tmpWnchar = sprintf('%s/W_%.2d_n.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpWnchar); Wn_tmp{1+nb1} = tutorial_binary_uncompress(tmpWnchar,1:nrows,1:ncols);
tmpWtchar = sprintf('%s/W_%.2d_t.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpWtchar); Wt_tmp{1+nb1} = tutorial_binary_uncompress(tmpWtchar,1:nrows,1:ncols);
tmpTnchar = sprintf('%s/T_%.2d_n.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpTnchar); Tn_tmp{1+nb1} = tutorial_binary_uncompress(tmpTnchar,1:nrows,1:ncols); ncols_T_tmp = ncols;
tmpTtchar = sprintf('%s/T_%.2d_t.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpTtchar); Tt_tmp{1+nb1} = tutorial_binary_uncompress(tmpTtchar,1:nrows,1:ncols);
tmpSnchar = sprintf('%s/S_%.2d_n.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpSnchar); Sn_tmp{1+nb1} = tutorial_binary_uncompress(tmpSnchar,1:nrows,1:ncols);
tmpStchar = sprintf('%s/S_%.2d_t.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpStchar); St_tmp{1+nb1} = tutorial_binary_uncompress(tmpStchar,1:nrows,1:ncols);
tmpAnchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAnchar); mr_A_tmp{1+nb1} = tutorial_binary_uncompress(tmpAnchar,1:nrows,1:ncols)>0; disp(sprintf('mr_A[%.2d]: %s',nb1,num2str(mr_A_tmp{1+nb1})));
tmpZnchar = sprintf('%s/mr_Z_%.2d.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpZnchar); mr_Z_tmp{1+nb1} = tutorial_binary_uncompress(tmpZnchar,1:nrows,1:ncols)>0; disp(sprintf('mr_Z[%.2d]: %s',nb1,num2str(mr_Z_tmp{1+nb1})));
end;%for nb1=0:nbins-1;
disp(sprintf(' %% mc_A: %f',norm(mc_A(:)-mc_A_tmp(:))));
disp(sprintf(' %% mc_Y: %f',norm(mc_Y(:)-mc_Y_tmp(:))));
disp(sprintf(' %% mc_T: %f',norm(mc_T(:)-mc_T_tmp(:))));
for nb1=0:nbins-1;
disp(sprintf(' %% An(%d): %f',nb1,norm(An{1+nb1}-An_tmp{1+nb1})));
disp(sprintf(' %% At(%d): %f',nb1,norm(At{1+nb1}-At_tmp{1+nb1})));
disp(sprintf(' %% Zn(%d): %f',nb1,norm(Zn{1+nb1}-Zn_tmp{1+nb1})));
disp(sprintf(' %% Zt(%d): %f',nb1,norm(Zt{1+nb1}-Zt_tmp{1+nb1})));
disp(sprintf(' %% Yn(%d): %f',nb1,norm(Yn{1+nb1}-Yn_tmp{1+nb1})));
disp(sprintf(' %% Yt(%d): %f',nb1,norm(Yt{1+nb1}-Yt_tmp{1+nb1})));
disp(sprintf(' %% Wn(%d): %f',nb1,norm(Wn{1+nb1}-Wn_tmp{1+nb1})));
disp(sprintf(' %% Wt(%d): %f',nb1,norm(Wt{1+nb1}-Wt_tmp{1+nb1})));
disp(sprintf(' %% Tn(%d): %f',nb1,norm(Tn{1+nb1}-Tn_tmp{1+nb1})));
disp(sprintf(' %% Tt(%d): %f',nb1,norm(Tt{1+nb1}-Tt_tmp{1+nb1})));
disp(sprintf(' %% Sn(%d): %f',nb1,norm(Sn{1+nb1}-Sn_tmp{1+nb1})));
disp(sprintf(' %% St(%d): %f',nb1,norm(St{1+nb1}-St_tmp{1+nb1})));
disp(sprintf(' %% mr_A(%d): %f',nb1,norm(mr_A{1+nb1}(:)-mr_A_tmp{1+nb1}(:))));
disp(sprintf(' %% mr_Z(%d): %f',nb1,norm(mr_Z{1+nb1}(:)-mr_Z_tmp{1+nb1}(:))));
end;%for nb1=0:nbins-1;
end;%if (check_flag);

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
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s/A_%.2d_n.b16',dir__in,nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s/A_%.2d_t.b16',dir__in,nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_A_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',ncols_A);
fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s/mr_A_%.2d.b16',dir__in,nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cind= %s/mc_A.b16;\n',dir__in);
fprintf(fp,'GLOBAL_T_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s/T_%.2d_n.b16',dir__in,nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s/T_%.2d_t.b16',dir__in,nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_T_n_cind= %s/mc_T.b16;\n',dir__in);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lakcluster_test_ADAADA_ver2.m on %s;\n',date);
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

for np1=0:2;for np2=0:2;
X_pr_{1+np1,1+np2} = intersect(r_on,r_ij_i(Xrij{1+np1,1+np2})); lX_pr_(1+np1,1+np2) = length(X_pr_{1+np1,1+np2});
X_pc_{1+np1,1+np2} = intersect(c_on,c_ij_i(Xcij{1+np1,1+np2})); lX_pc_(1+np1,1+np2) = length(X_pc_{1+np1,1+np2});
A_pr_{1+np1,1+np2} = setdiff(r_on,X_pr_{1+np1,1+np2}); lA_pr_{1+np1,1+np2} = length(A_pr_{1+np1,1+np2});
A_pc_{1+np1,1+np2} = setdiff(c_on,X_pc_{1+np1,1+np2}); lA_pc_{1+np1,1+np2} = length(A_pc_{1+np1,1+np2});
end;end;%for np1=0:2;for np2=0:2;
A_pr_x = setdiff(r_on,unionall({X_pr_{1,1},X_pr_{1,2},X_pr_{1,3},X_pr_{2,1},X_pr_{2,2},X_pr_{2,3},X_pr_{3,1},X_pr_{3,2},X_pr_{3,3}})); 
A_pc_x = setdiff(r_on,unionall({X_pc_{1,1},X_pc_{1,2},X_pc_{1,3},X_pc_{2,1},X_pc_{2,2},X_pc_{2,3},X_pc_{3,1},X_pc_{3,2},X_pc_{3,3}})); 

% compares each X_pr to the other indices, even those in other biclusters ;
%{
ARout = [    get_auc((A_pr_{1,1}),(X_pr_{1,1})) , get_auc((A_pr_{1,2}),(X_pr_{1,2})) , get_auc((A_pr_{1,3}),(X_pr_{1,3})) ;...
	     get_auc((A_pr_{2,1}),(X_pr_{2,1})) , get_auc((A_pr_{2,2}),(X_pr_{2,2})) , get_auc((A_pr_{2,3}),(X_pr_{2,3})) ;...
	     get_auc((A_pr_{3,1}),(X_pr_{3,1})) , get_auc((A_pr_{3,2}),(X_pr_{3,2})) , get_auc((A_pr_{3,3}),(X_pr_{3,3})) ;...
	 ];
ACout = [    get_auc((A_pc_{1,1}),(X_pc_{1,1})) , get_auc((A_pc_{1,2}),(X_pc_{1,2})) , get_auc((A_pc_{1,3}),(X_pc_{1,3})) ;...
	     get_auc((A_pc_{2,1}),(X_pc_{2,1})) , get_auc((A_pc_{2,2}),(X_pc_{2,2})) , get_auc((A_pc_{2,3}),(X_pc_{2,3})) ;...
	     get_auc((A_pc_{3,1}),(X_pc_{3,1})) , get_auc((A_pc_{3,2}),(X_pc_{3,2})) , get_auc((A_pc_{3,3}),(X_pc_{3,3})) ;...
	 ];
 %}

% compares each X_pr only to the other indices which are not in other biclusters ;
ARout = [    get_auc((A_pr_x),(X_pr_{1,1})) , get_auc((A_pr_x),(X_pr_{1,2})) , get_auc((A_pr_x),(X_pr_{1,3})) ;...
	     get_auc((A_pr_x),(X_pr_{2,1})) , get_auc((A_pr_x),(X_pr_{2,2})) , get_auc((A_pr_x),(X_pr_{2,3})) ;...
	     get_auc((A_pr_x),(X_pr_{3,1})) , get_auc((A_pr_x),(X_pr_{3,2})) , get_auc((A_pr_x),(X_pr_{3,3})) ;...
	 ];
ACout = [    get_auc((A_pc_x),(X_pc_{1,1})) , get_auc((A_pc_x),(X_pc_{1,2})) , get_auc((A_pc_x),(X_pc_{1,3})) ;...
	     get_auc((A_pc_x),(X_pc_{2,1})) , get_auc((A_pc_x),(X_pc_{2,2})) , get_auc((A_pc_x),(X_pc_{2,3})) ;...
	     get_auc((A_pc_x),(X_pc_{3,1})) , get_auc((A_pc_x),(X_pc_{3,2})) , get_auc((A_pc_x),(X_pc_{3,3})) ;...
	 ];

% compares each X_pr to each other X_pr in the same sparsity class;
XRout = [    get_auc((X_pr_{1,1}),(X_pr_{2,1})) , get_auc((X_pr_{1,1}),(X_pr_{3,1})) , get_auc((X_pr_{2,1}),(X_pr_{3,1})) ;...
             get_auc((X_pr_{1,2}),(X_pr_{2,2})) , get_auc((X_pr_{1,2}),(X_pr_{3,2})) , get_auc((X_pr_{2,2}),(X_pr_{3,2})) ;...
             get_auc((X_pr_{1,3}),(X_pr_{2,3})) , get_auc((X_pr_{1,3}),(X_pr_{3,3})) , get_auc((X_pr_{2,3}),(X_pr_{3,3})) ;...
	 ];
XCout = [    get_auc((X_pc_{1,1}),(X_pc_{2,1})) , get_auc((X_pc_{1,1}),(X_pc_{3,1})) , get_auc((X_pc_{2,1}),(X_pc_{3,1})) ;...
             get_auc((X_pc_{1,2}),(X_pc_{2,2})) , get_auc((X_pc_{1,2}),(X_pc_{3,2})) , get_auc((X_pc_{2,2}),(X_pc_{3,2})) ;...
             get_auc((X_pc_{1,3}),(X_pc_{2,3})) , get_auc((X_pc_{1,3}),(X_pc_{3,3})) , get_auc((X_pc_{2,3}),(X_pc_{3,3})) ;...
	 ];

disp(sprintf(' %% ARout: ')); disp(num2str(ARout));
disp(sprintf(' %% ACout: ')); disp(num2str(ACout));
disp(sprintf(' %% XRout: ')); disp(num2str(XRout));
disp(sprintf(' %% XCout: ')); disp(num2str(XCout));

xx_ra = transpose([ARout(:);ACout(:);XRout(:);XCout(:)]);
if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

clean_flag=1;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;

