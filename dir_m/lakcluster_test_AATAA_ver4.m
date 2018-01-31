function lakcluster_test_AATAA_ver4(dir_trunk,N,X_factor,X_esm,gamma,B_MLT,rng_num,pt_num)
% Embedding a single mds-balanced cluster in An{nbins}, balanced w.r.t. 2 mds-components. ;
% Adding two other mds-imbalanced clusters (of double and quadruple the size, respectively) to An{nbins}.;
% using lakcluster_ver18 ;
% try: ;
%{
  
  N=1024;X_factor=0.55;X_esm=0.1;gamma=0.05;B_MLT=0;rng_num=1;lakcluster_test_AATAA_ver4(pwd,N,X_factor,X_esm,gamma,B_MLT,rng_num);

  %}


if (nargin<8); pt_num = -1; end;
if (nargin<7); rng_num = 1; end;
if (nargin<4); B_MLT = 32; end;
if (nargin<5); gamma = 0.025; end;
if (nargin<4); X_esm = 0.1; end;
if (nargin<3); X_factor = 0.65; end;
if (nargin<2); N = 1024; end;
if (nargin<1); dir_trunk = pwd; end;

rng(rng_num);
nbins = 1;
nrows_A_ = [N];
nrows_Z_ = [0];
ncols_A = N+1;
ncols_Y = 0;
nmds = 2; 
ncols_T = nmds+1; MDS = randn(N,nmds);
mrnd = 0;
QR_strategy = 'YnWt';
QC_strategy = 'ZtSWn';
verbose_flag=0;

nbx = nbins;% bin to implant X;
nrows_X = ceil(nrows_A_(nbx).^(X_factor));
ncols_X = ceil(ncols_A.^(X_factor));
X_eps = X_esm/sqrt(nrows_X); % spectral error for X;
X_pp = tutorial_g1em(X_esm);

test_string = 'lakcluster_test_AATAA';
dir_base = sprintf('%s/dir_%s',dir_trunk,test_string); if ~exist(dir_base,'dir'); mkdir(dir_base); end;
[posfix_use,X_factor_d,X_esm_d,gamma_d] = xfix_gen(N,X_factor,X_esm,gamma,B_MLT,rng_num);
in__name = sprintf('%s_%s__in',test_string,posfix_use);
out_name = sprintf('%s_%s_out',test_string,posfix_use);
dir__in = sprintf('%s/dir_%s',dir_base,in__name); if ~exist(dir__in,'dir'); mkdir(dir__in); end;
dir_out = sprintf('%s/dir_%s',dir_base,out_name); if ~exist(dir_out,'dir'); mkdir(dir_out); end; 
%path(path,dir__in); path(path,dir_out);

clear An At Zn Zt Yn Yt Wn Wt Tn Tt Sn St mr_A mr_Z mc_A mc_Y mc_T bitj;
bitj = 16;

mc_A = (rand(1,ncols_A)>mrnd);
mc_Y = (rand(1,ncols_Y)>mrnd);
mc_T = (rand(1,ncols_T)>mrnd); mc_T(1)=1; % ensure first column of [T;S] is all ones ;
mc_T = ones(1,ncols_T);
for nb1=0:nbins-1; 
mr_A{1+nb1} = (rand(nrows_A_(1+nb1),1)>mrnd);
end;%for nb1=0:nbins-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%u1 = 1; u2 = 2; u3 = 4; 
%u1 = 1; u2 = 0; u3 = 1; % 082016
%u1 = 1; u2 = 0; u3 = 2; % 111716
u1 = 1; u2 = 2; u3 = 4; % 112116
X_1 = tutorial_binary_makelr2(u1*nrows_X,u1*ncols_X,X_pp);
X_2 = tutorial_binary_makelr2(u2*nrows_X,u2*ncols_X,X_pp);
X_3 = tutorial_binary_makelr2(u3*nrows_X,u3*ncols_X,X_pp);

for nb1=0:nbins-1; 
An{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_A)>0.5)-1; 
Tn{1+nb1} = 2*(rand(nrows_A_(1+nb1),ncols_T)>0.5)-1; Tn{1+nb1}(:,1) = 1; 
if nb1==0; Tn{1+0} = [ones(nrows_A_(1+0),1) , 2*(MDS>0)-1]; end;%nb1==0;
end;%for nb1=0:nbins-1;

u_mds = randn(nmds,1); u_mds = u_mds/norm(u_mds);
T_u = MDS*u_mds;
T_u_pos = find(T_u>0);
T_u_neg = find(T_u<=0);

%insert clusters;
pc1_f = randperm(ncols_A);       [tmp,pc1_i] = sort(pc1_f); [tmp,pc1_i] = sort(pc1_i);
pr1_f = randperm(nrows_A_(1+0)); [tmp,pr1_i] = sort(pr1_f); [tmp,pr1_i] = sort(pr1_i);
pc2_f = randperm(ncols_A);       [tmp,pc2_i] = sort(pc2_f); [tmp,pc2_i] = sort(pc2_i);
pr2_f = [T_u_pos(randperm(length(T_u_pos)));T_u_neg(randperm(length(T_u_neg)))]; [tmp,pr2_i] = sort(pr2_f); [tmp,pr2_i] = sort(pr2_i);
pc3_f = randperm(ncols_A);       [tmp,pc3_i] = sort(pc3_f); [tmp,pc3_i] = sort(pc3_i);
pr3_f = [T_u_neg(randperm(length(T_u_neg)));T_u_pos(randperm(length(T_u_pos)))]; [tmp,pr3_i] = sort(pr3_f); [tmp,pr3_i] = sort(pr3_i);

An{1+0}(pr1_f(1:u1*nrows_X),pc1_f(1:u1*ncols_X)) = X_1;
An{1+0}(pr2_f(1:u2*nrows_X),pc2_f(1:u2*ncols_X)) = X_2;
An{1+0}(pr3_f(1:u3*nrows_X),pc3_f(1:u3*ncols_X)) = X_3;

disp_flag=0;
if disp_flag;
figure;imagesc([ Tn{1+0}(pr1_i,:) , zeros(nrows_A_(1+0),8) , repmat(T_u(pr1_i),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(pr1_i,pc1_i) ],[-1,1]); title('X_1');
figure;imagesc([ Tn{1+0}(pr2_i,:) , zeros(nrows_A_(1+0),8) , repmat(T_u(pr2_i),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(pr2_i,pc2_i) ],[-1,1]); title('X_2');
figure;imagesc([ Tn{1+0}(pr3_i,:) , zeros(nrows_A_(1+0),8) , repmat(T_u(pr3_i),1,8) , zeros(nrows_A_(1+0),8) , An{1+0}(pr3_i,pc3_i) ],[-1,1]); title('X_3');
return;
end;%if disp_flag;

for nb1=0:nbins-1; 
At{1+nb1} = transpose(An{1+nb1}); 
Tt{1+nb1} = transpose(Tn{1+nb1}); 
end;%for nb1=0:nbins-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% expanding T ; used for testing algorithm ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncols_T_expand = 1 + 2*nmds; prm_T = randperm(2*nmds); mds_ij = find(prm_T<=nmds); %disp(prm_T); disp(mds_ij); % T_expand used to spread out columns of T ;
T_ij = [1 , 1 + (mds_ij)];
for nb1=0:nbins-1;
Tn_expand{1+nb1} = randn(size(Tn{1+nb1},1),ncols_T_expand);
Tn_expand{1+nb1}(:,T_ij) = Tn{1+nb1};
Tt_expand{1+nb1} = transpose(Tn_expand{1+nb1});
end;%for nb1=0:nbins-1;
mc_T_expand = zeros(1,ncols_T_expand); mc_T_expand(T_ij) = mc_T;
tmpAtchar = sprintf('%s/mc_A.b16',dir__in);tutorial_binary_compress(bitj,mc_A(:)>0,tmpAtchar); %disp(sprintf('mc_A: %s',num2str(mc_A)));
tmpTtchar = sprintf('%s/mc_T.b16',dir__in);tutorial_binary_compress(bitj,mc_T(:)>0,tmpTtchar); %disp(sprintf('mc_T: %s',num2str(mc_T)));
tmpTtchar = sprintf('%s/mc_T_expand.b16',dir__in);tutorial_binary_compress(bitj,mc_T_expand(:)>0,tmpTtchar); %disp(sprintf('mc_T_expand: %s',num2str(mc_T_expand)));
%disp(sprintf(' %% error %0.16f',norm(Tn_expand{1}(:,T_ij) - Tn{1}(:,:))));
for nb1=0:nbins-1;
tmpAnchar = sprintf('%s/A_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,An{1+nb1}>0,tmpAnchar); 
tmpAtchar = sprintf('%s/A_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,At{1+nb1}>0,tmpAtchar); 
tmpTnchar = sprintf('%s/T_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Tn{1+nb1}>0,tmpTnchar); 
tmpTtchar = sprintf('%s/T_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Tt{1+nb1}>0,tmpTtchar); 
tmpTnchar = sprintf('%s/T_expand_%.2d_n.b16',dir__in,nb1);tutorial_binary_compress(bitj,Tn_expand{1+nb1}>0,tmpTnchar); 
tmpTtchar = sprintf('%s/T_expand_%.2d_t.b16',dir__in,nb1);tutorial_binary_compress(bitj,Tt_expand{1+nb1}>0,tmpTtchar); 
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
use_T_expand = 0;%use_T_expand = input(' use T_expand [0,1]?');
if use_T_expand;
fprintf(fp,'GLOBAL_T_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'T_expand_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'T_expand_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T_expand);
fprintf(fp,'GLOBAL_T_n_cind= mc_T_expand.b16;\n');
 else;
fprintf(fp,'GLOBAL_T_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'T_%.2d_n.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'T_%.2d_t.b16',nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T);
fprintf(fp,'GLOBAL_T_n_cind= mc_T.b16;\n');
end;%if use_T_expand; 
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lakcluster_test_AATAA_ver4.m on %s;\n',date);
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
nb1 = nbins-1; nb2 = nbins-1; nb_tab = nb1+nb2*nbins;
r_on = find(mr_A_tmp{1+nb1}); c_on = find(mc_A_tmp);

X_pr1_i = intersect(r_on,r_ij_i(pr1_f(1:u1*nrows_X))); lXr1 = length(X_pr1_i); 
X_pc1_i = intersect(c_on,c_ij_i(pc1_f(1:u1*ncols_X))); lXc1 = length(X_pc1_i);
X_pr2_i = intersect(r_on,r_ij_i(pr2_f(1:u2*nrows_X))); lXr2 = length(X_pr2_i); 
X_pc2_i = intersect(c_on,c_ij_i(pc2_f(1:u2*ncols_X))); lXc2 = length(X_pc2_i);
X_pr3_i = intersect(r_on,r_ij_i(pr3_f(1:u3*nrows_X))); lXr3 = length(X_pr3_i); 
X_pc3_i = intersect(c_on,c_ij_i(pc3_f(1:u3*ncols_X))); lXc3 = length(X_pc3_i);

A_pr1_i = setdiff(r_on,X_pr1_i); lAr1 = length(A_pr1_i);
A_pr2_i = setdiff(r_on,X_pr2_i); lAr2 = length(A_pr2_i);
A_pr3_i = setdiff(r_on,X_pr3_i); lAr3 = length(A_pr3_i);
A_pc1_i = setdiff(c_on,X_pc1_i); lAc1 = length(A_pc1_i);
A_pc2_i = setdiff(c_on,X_pc2_i); lAc2 = length(A_pc2_i);
A_pc3_i = setdiff(c_on,X_pc3_i); lAc3 = length(A_pc3_i);

auc_c_x1 = get_auc(A_pc1_i,X_pc1_i); auc_r_x1 = get_auc(A_pr1_i,X_pr1_i);
auc_c_x2 = get_auc(A_pc2_i,X_pc2_i); auc_r_x2 = get_auc(A_pr2_i,X_pr2_i);
auc_c_x3 = get_auc(A_pc3_i,X_pc3_i); auc_r_x3 = get_auc(A_pr3_i,X_pr3_i);

disp(sprintf(' %% rows: (X1/A1) = %d/%d = %0.0f%% --> auc_r_x1 %0.2f',lXr1,lAr1,100*lXr1/lAr1,auc_r_x1));
disp(sprintf(' %% cols: (X1/A1) = %d/%d = %0.0f%% --> auc_c_x1 %0.2f',lXc1,lAc1,100*lXc1/lAc1,auc_c_x1));
disp(sprintf(' %% rows: (X2/A2) = %d/%d = %0.0f%% --> auc_r_x2 %0.2f',lXr2,lAr2,100*lXr2/lAr2,auc_r_x2));
disp(sprintf(' %% cols: (X2/A2) = %d/%d = %0.0f%% --> auc_c_x2 %0.2f',lXc2,lAc2,100*lXc2/lAc2,auc_c_x2));
disp(sprintf(' %% rows: (X3/A3) = %d/%d = %0.0f%% --> auc_r_x3 %0.2f',lXr3,lAr3,100*lXr3/lAr3,auc_r_x3));
disp(sprintf(' %% cols: (X3/A3) = %d/%d = %0.0f%% --> auc_c_x3 %0.2f',lXc3,lAc3,100*lXc3/lAc3,auc_c_x3));

xx_ra = [auc_r_x1,auc_c_x1,0.5*(auc_r_x1+auc_c_x1), auc_r_x2,auc_c_x2,0.5*(auc_r_x2+auc_c_x2), auc_r_x3,auc_c_x3,0.5*(auc_r_x3+auc_c_x3)];
if pt_num>0; fname_xx = sprintf('%s/xx_out.txt.pt%.3d',dir_base,pt_num); else; fname_xx = sprintf('%s/xx_out.txt',dir_base); end;
fid = fopen(fname_xx,'a'); fprintf(fid,'%.5d,%.2d,%.2d,%.2d,%.2d,%.3d,%0.2f,%0.2f',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num,elrt,elct); fprintf(fid,',%0.16f',xx_ra); fprintf(fid,'\n'); fclose(fid);

clean_flag=1;
if clean_flag;
xfix_del(dir_trunk,test_string,N,X_factor,X_esm,gamma,B_MLT,rng_num,nbins);
end;%if clean_flag;

