function lakcluster_test_AAAA_ver5(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num)
% starting with basic cluster Xn embedded in An{nbins}: ; testing lakcluster_ver18 ;
% try: ;
%{
  
  N=1024;X_factor=0.55;X_esm=0.1;gamma=0.05;B_MLT=32;rng_num=1;lakcluster_test_AAAA_ver5(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

  nbins=1; % for AAAA;
  N_ra = 2.^[6:12]; X_factor_ra = [0.45:0.025:0.65]; X_esm_ra = 10.^[-3:0.5:1]; gamma_ra = [0,0.5,0.15,0.5]; B_ra = [0,4,8,16]; rng_num_ra = 1:4;
  for nN=1:length(N_ra); N = N_ra(nN);
  for nf=1:length(X_factor_ra); X_factor = X_factor_ra(nf);
  for ne=1:length(X_esm_ra); X_esm = X_esm_ra(ne);
  for ng=1:length(gamma_ra); gamma = gamma_ra(ng);
  for nB=1:length(B_ra); B_MLT = B_ra(nB);
  for nr=1:length(rng_num_ra); rng_num = rng_num_ra(nr);
  lakcluster_test_AAAA_ver5(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);
  xfix_del(pwd,'lakcluster_test_AAAA',N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
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
ncols_A = 2*N+1;
ncols_Y = 0;
ncols_T = 1;
mrnd = 0;
QR_strategy = 'YnWt';
QC_strategy = 'YnWt';
verbose_flag=0;

nbc = nbins;% bin to implant X;
nrows_X = ceil(nrows_A_(nbc).^(X_factor));
ncols_X = ceil(ncols_A.^(X_factor));
X_eps = X_esm/sqrt(nrows_X); % spectral error for X;
X_pp = tutorial_g1em(X_esm);

test_string = 'lakcluster_test_AAAA';
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
if (1+nb1)==nbc; % insert cluster ;
X = tutorial_binary_makelr2(nrows_X,ncols_X,X_pp);
pr_f{1+nb1} = randperm(nrows_A_(1+nb1)); [tmp,pr_i{1+nb1}] = sort(pr_f{1+nb1}); [tmp,pr_i{1+nb1}] = sort(pr_i{1+nb1});
pc_f{1+nb1} = randperm(ncols_A);         [tmp,pc_i{1+nb1}] = sort(pc_f{1+nb1}); [tmp,pc_i{1+nb1}] = sort(pc_i{1+nb1});
An{1+nb1}(pr_f{1+nb1}(1:nrows_X),pc_f{1+nb1}(1:ncols_X)) = X;
end;%if (1+nb1)==nbins; 
At{1+nb1} = transpose(An{1+nb1});
Zn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_A)>0.5)-1; Zt{1+nb1} = transpose(Zn{1+nb1});
Yn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_Y)>0.5)-1; Yt{1+nb1} = transpose(Yn{1+nb1});
Wn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_Y)>0.5)-1; Wt{1+nb1} = transpose(Wn{1+nb1});
Tn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_T)>0.5)-1; Tn{1+nb1}(:,1) = 1; Tt{1+nb1} = transpose(Tn{1+nb1});
Sn{1+nb1} = 2*(rand(nrows_Z_(1+nb1),ncols_T)>0.5)-1; Sn{1+nb1}(:,1) = 1; St{1+nb1} = transpose(Sn{1+nb1});
end;%for nb1=0:nbins-1;

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
%fprintf(fp,'GLOBAL_TEST_sparse= 0;\n');
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
fprintf(fp,'%% generated by lakcluster_test_AAAA_ver5.m on %s;\n',date);
fclose(fp);

call_flag=1;%call_flag = input(' call? 1=yes (default), 0=no:'); if isempty(call_flag); call_flag=1; end;
if call_flag;
disp(sprintf('%s/../lakcluster_ver18 < %s/%s.in',dir_trunk,dir__in,in__name));
system(sprintf('%s/../lakcluster_ver18 < %s/%s.in',dir_trunk,dir__in,in__name));
end;%if call_flag;

run(sprintf('%s/timing',dir_out));
tmpAtchar = sprintf('%s/mc_A.b16',dir__in);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAtchar); mc_A_tmp = tutorial_binary_uncompress(tmpAtchar,1:nrows,1:ncols)>0; 
for nb1=0:nbins-1; tmpAnchar = sprintf('%s/mr_A_%.2d.b16',dir__in,nb1);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpAnchar); mr_A_tmp{1+nb1} = tutorial_binary_uncompress(tmpAnchar,1:nrows,1:ncols)>0; end;

out_xdrop_a = textread(sprintf('%s/out_xdrop_a.txt',dir_out)); out_xdrop_b = textread(sprintf('%s/out_xdrop_b.txt',dir_out));
r_ij = out_xdrop_a(find(out_xdrop_a(:,1)>-1),1); c_ij = out_xdrop_a(find(out_xdrop_a(:,2)>-1),2);
[r_ij_ord,r_ij_i] = sort(1+r_ij); [c_ij_ord,c_ij_i] = sort(1+c_ij);
nb1 = nbins-1; nb2 = nbins-1; nb_tab = nb1+nb2*nbins;
r_on = find(mr_A_tmp{1+nb1}); c_on = find(mc_A_tmp);
A_pr_i = intersect(r_on,r_ij_i(pr_f{1+nb1}(nrows_X+1:nrows_A_(1+nb1)))); lAr{1+nb1} = length(A_pr_i);
X_pr_i = intersect(r_on,r_ij_i(pr_f{1+nb1}(1:nrows_X))); lXr{1+nb1} = length(X_pr_i); 
A_pc_i = intersect(c_on,c_ij_i(pc_f{1+nb1}(ncols_X+1:ncols_A))); lAc = length(A_pc_i);
X_pc_i = intersect(c_on,c_ij_i(pc_f{1+nb1}(1:ncols_X))); lXc = length(X_pc_i);

auc_r = get_auc(A_pr_i,X_pr_i);
auc_c = get_auc(A_pc_i,X_pc_i);
disp(sprintf(' %% rows: (X/A) = %d/%d = %0.0f%% --> auc_r %0.2f',lXr{1+nb1},lAr{1+nb1},100*lXr{1+nb1}/lAr{1+nb1},auc_r));
disp(sprintf(' %% cols: (X/A) = %d/%d = %0.0f%% --> auc_c %0.2f',lXc,lAc,100*lXc/lAc,auc_c));

xx_ra = [auc_r,auc_c,0.5*(auc_r+auc_c)];
if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

plot_flag=0;
if plot_flag;
figure;
nb1 = 1;
subplot(1,2,1); imagesc(An{nb1}(pr_i{nb1},pc_i{nb1}),[-1,+1]);
subplot(1,2,2); imagesc(An{nb1}(1+r_ij(end:-1:1),1+c_ij(end:-1:1)),[-1,+1]);
end;%if plot_flag;

clean_flag=1;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;

