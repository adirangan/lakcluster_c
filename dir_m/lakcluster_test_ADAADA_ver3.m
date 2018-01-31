function lakcluster_test_ADAADA_ver3(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num) ;
% Two biclusters, both embedded in an A with sparsity sx_ra(1) = 0.0625;
% X_1 = tutorial_binary_makelr3(1*nrows_X,1*ncols_X,X_pp,2*sx_ra(1)); base size, double sparsity ;
% X_2 = tutorial_binary_makelr3(2*nrows_X,2*ncols_X,X_pp,1*sx_ra(1)); double size, base sparsity ;
% try: ;
%{
  
  N=128;X_factor=0.65;X_esm=0.01;gamma=0.05;B_MLT=32;rng_num=1;lakcluster_test_ADAADA_ver3(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

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

%sx_ra(1) = 0.125;
sx_ra(1) = 0.0625;
X_1 = tutorial_binary_makelr3(1*nrows_X,1*ncols_X,X_pp,2*sx_ra(1));
X_2 = tutorial_binary_makelr3(2*nrows_X,2*ncols_X,X_pp,1*sx_ra(1));

for nb1=0:nbins-1; 
An{1+nb1} = [2*(rand(nrows_A_(1+nb1),ncols_A)<sx_ra(1))-1 ];
if (1+nb1)==nbx; % insert clusters ;
pr1__f = randperm(nrows_A_(1)); [tmp,pr1__i] = sort(pr1__f); [tmp,pr1__i] = sort(pr1__i);
pc1__f = randperm(ncols_A);       [tmp,pc1__i] = sort(pc1__f); [tmp,pc1__i] = sort(pc1__i);
pr2__f = randperm(nrows_A_(1)); [tmp,pr2__i] = sort(pr2__f); [tmp,pr2__i] = sort(pr2__i);
pc2__f = randperm(ncols_A);       [tmp,pc2__i] = sort(pc2__f); [tmp,pc2__i] = sort(pc2__i);
An{1+0}(pr2__f(1:2*nrows_X),pc2__f(1:2*ncols_X)) = X_2;
An{1+0}(pr1__f(1:1*nrows_X),pc1__f(1:1*ncols_X)) = X_1;
[tmp_r,tmp_r2,tmp_r1] = intersect(pr2__f(1:2*nrows_X),pr1__f(1:1*nrows_X));
[tmp_c,tmp_c2,tmp_c1] = intersect(pc2__f(1:2*ncols_X),pc1__f(1:1*ncols_X));
for nr=1:length(tmp_r);
for nc=1:length(tmp_c);
if (rand()>0.5); An{1+0}(tmp_r(nr),tmp_c(nc)) = X_1(tmp_r1(nr),tmp_c1(nc));
 else; An{1+0}(tmp_r(nr),tmp_c(nc)) = X_2(tmp_r2(nr),tmp_c2(nc));
end;%if ;
end;%for nc=1:length(tmp_c);
end;%for nr=1:length(tmp_r);
end;%if (1+nb1)==nbins; 
At{1+nb1} = transpose(An{1+nb1}); 
end;%for nb1=0:nbins-1;

disp_flag=0;
if disp_flag;
figure;imagesc([ An{1+0}(pr1__i,pc1__i) ]); title('X_1');
figure;imagesc([ An{1+0}(pr2__i,pc2__i) ]); title('X_2');
return;
end;%if disp_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpAtchar = sprintf('%s/mc_A.b16',dir__in);tutorial_binary_compress(bitj,mc_A(:)>0,tmpAtchar); %disp(sprintf('mc_A: %s',num2str(mc_A)));
for nb1=0:nbins-1;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,An{1+nb1}>0,tmpAnchar); 
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,At{1+nb1}>0,tmpAtchar); 
tmpAtchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);tutorial_binary_compress(bitj,mr_A{1+nb1}(:)>0,tmpAtchar); %disp(sprintf('mr_A_[%d]: %s',nb1,num2str(mr_A{1+nb1})));
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
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',nrows_Z_(1+nb1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',ncols_Y);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lakcluster_test_ADAADA_ver3.m on %s;\n',date);
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
end;%for nb1=0:nbins-1; 

out_xdrop_a = textread(sprintf('%s/out_xdrop_a.txt',dir_out)); out_xdrop_b = textread(sprintf('%s/out_xdrop_b.txt',dir_out));
r_ij = out_xdrop_a(find(out_xdrop_a(:,1)>-1),1); c_ij = out_xdrop_a(find(out_xdrop_a(:,2)>-1),2);
[r_ij_ord,r_ij_i] = sort(1+r_ij); [c_ij_ord,c_ij_i] = sort(1+c_ij);
clear r_on_ c_on ;
for nb1=0:nbins-1; r_on_{1+nb1} = find(mr_A_tmp{1+nb1}); end;
c_on = find(mc_A_tmp); 
r_on_x = r_on_{1};

X_pc1_ = intersect(c_on,c_ij_i(pc1__f(1:1*ncols_X))); A_pc1_ = setdiff(c_on,X_pc1_);
X_pc2_ = intersect(c_on,c_ij_i(pc2__f(1:2*ncols_X))); A_pc2_ = setdiff(c_on,X_pc2_);
pr1__g = pr1__f;
pr2__g = pr2__f;

X_pr1_ = intersect(r_on_x,r_ij_i(pr1__g(1:1*nrows_X))); A_pr1_ = setdiff(r_on_x,X_pr1_);
X_pr2_ = intersect(r_on_x,r_ij_i(pr2__g(1:2*nrows_X))); A_pr2_ = setdiff(r_on_x,X_pr2_);

auc_c_1 = get_auc((A_pc1_),(X_pc1_)); auc_r_1 = get_auc((A_pr1_),(X_pr1_)); auc_x_1 = 0.5*(auc_c_1 + auc_r_1);
auc_c_2 = get_auc((A_pc2_),(X_pc2_)); auc_r_2 = get_auc((A_pr2_),(X_pr2_)); auc_x_2 = 0.5*(auc_c_2 + auc_r_2);
xx_ra = [auc_c_1,auc_r_1,auc_x_1,auc_c_2,auc_r_2,auc_x_2];
disp(num2str(xx_ra));

if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

clean_flag=0;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;

