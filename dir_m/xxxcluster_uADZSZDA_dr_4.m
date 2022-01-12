function ...
[ ...
 parameter ...
] = ...
xxxcluster_uADZSZDA_dr_4( ...
 parameter ...
,str_prefix ...
,M_n_ ...
,A_n_rij_ ...
,A_n_cij ...
,Z_n_rij_ ...
,T_n_ ...
,T_n_cij ...
);
% removing references to Z and S if empty. ;
% test with: ;
%{
  xxxcluster_uADZSZDA_dr_4([]);
  %}

if (nargin<1);
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); str_prefix=[]; end; na=na+1;
if (nargin<1+na); M_n_=[]; end; na=na+1;
if (nargin<1+na); A_n_rij_=[]; end; na=na+1;
if (nargin<1+na); A_n_cij=[]; end; na=na+1;
if (nargin<1+na); Z_n_rij_=[]; end; na=na+1;
if (nargin<1+na); T_n_=[]; end; na=na+1;
if (nargin<1+na); T_n_cij=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if isempty(str_prefix); str_prefix = 'test'; end;
if isempty(M_n_); M_n_ = {randn(1024,1024)}; end;
if isempty(A_n_rij_); A_n_rij_ = {[1+0*512+[0:511]]}; end;
if isempty(A_n_cij); A_n_cij = [1:512]; end;
if isempty(Z_n_rij_); Z_n_rij_ = {[1+1*512+[0:511]]}; end;
if isempty(T_n_); T_n_ = {ones(1024,1)}; end;
if isempty(T_n_cij); T_n_cij = 1; end;

if ~isfield(parameter,'tolerance_master'); parameter.tolerance_master = 10^(-9.5); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'dir_trunk'); parameter.dir_trunk = pwd; end;
if ~isfield(parameter,'dir_code'); parameter.dir_code = sprintf('%s/..',pwd); end;
if ~isfield(parameter,'str_code'); parameter.str_code = 'lakcluster_ver18'; end;
if ~isfield(parameter,'str_lak_vs_dex'); parameter.str_lak_vs_dex = 'dex'; end;
if ~isfield(parameter,'flag_reverse'); parameter.flag_reverse = 0; end;
if ~isfield(parameter,'GLOBAL_TEST_sparse'); parameter.GLOBAL_TEST_sparse = 0; end;
if ~isfield(parameter,'GLOBAL_kappa_squared'); parameter.GLOBAL_kappa_squared = 0; end;
if ~isfield(parameter,'QR_strategy'); parameter.QR_strategy = 'YnWt condense'; end;
if ~isfield(parameter,'QC_strategy'); parameter.QC_strategy = 'YnWt store one'; end;
if ~isfield(parameter,'gamma'); parameter.gamma = 0.025; end;
if ~isfield(parameter,'B_MLT'); parameter.B_MLT = abs(floor(log2(parameter.tolerance_master))); end;
if ~isfield(parameter,'Ireq'); parameter.Ireq = 0; end;
if ~isfield(parameter,'shuffle_num'); parameter.shuffle_num = 0; end;
if ~isfield(parameter,'pt_num'); parameter.pt_num = -1; end;
if ~isfield(parameter,'shuffle_num'); parameter.shuffle_num = 0; end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
if ~isfield(parameter,'flag_force_create'); parameter.flag_force_create = 0; end;
if ~isfield(parameter,'slurm_walltime'); parameter.slurm_walltime = 0; end;
if ~isfield(parameter,'slurm_nnodes'); parameter.slurm_nnodes = 1; end;
if ~isfield(parameter,'slurm_tpn'); parameter.slurm_tpn = 1; end;
if ~isfield(parameter,'slurm_memdecl'); parameter.slurm_memdecl = 1; end;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;
dir_trunk = parameter.dir_trunk;
dir_code = parameter.dir_code;
str_code = parameter.str_code;
str_lak_vs_dex = parameter.str_lak_vs_dex;
flag_reverse = parameter.flag_reverse;
GLOBAL_TEST_sparse = parameter.GLOBAL_TEST_sparse;
GLOBAL_kappa_squared = parameter.GLOBAL_kappa_squared;
QR_strategy = parameter.QR_strategy;
QC_strategy = parameter.QC_strategy;
gamma = parameter.gamma;
B_MLT = parameter.B_MLT;
Ireq = parameter.Ireq;
shuffle_num = parameter.shuffle_num;
pt_num = parameter.pt_num;
shuffle_num = parameter.shuffle_num;
flag_verbose = parameter.flag_verbose;
flag_force_create = parameter.flag_force_create;
slurm_walltime = parameter.slurm_walltime;
slurm_nnodes = parameter.slurm_nnodes;
slurm_tpn = parameter.slurm_tpn;
slurm_memdecl = parameter.slurm_memdecl;

nbins=length(M_n_);
[ ...
A_n_cols ...
,Y_n_cols ...
,T_n_cols ...
,~ ...
,~ ...
] = ...
xxxcluster_uADZSZDA_check_1( ...
 shuffle_num ...
,M_n_ ...
,A_n_rij_ ...
,A_n_cij ...
,Z_n_rij_ ...
,T_n_ ...
,T_n_cij ...
);

bitj = 16;

Z_bother = 0;
for nb1=0:nbins-1;
if (length(Z_n_rij_{1+nb1})>0); Z_bother = 1; end;
end;%for nb1=0:nbins-1;
if ~Z_bother; flag_reverse = 0; end;

str_out_suffix = sprintf('%s_%s',str_prefix,xxxcluster_uADZSZDA_xfix_gen_ver1(str_lak_vs_dex,flag_reverse,A_n_rij_,Z_n_rij_,T_n_cij,GLOBAL_TEST_sparse,gamma,B_MLT,Ireq,shuffle_num));
dir_0in = sprintf('%s/dir_%s',dir_trunk,str_prefix);
dir_out = sprintf('%s/dir_%s',dir_0in,str_out_suffix); 
disp(sprintf(' str_out_suffix: %s',str_out_suffix));
disp(sprintf(' dir_0in: %s',dir_0in));
disp(sprintf(' dir_out: %s',dir_out));
if ~exist(dir_out,'dir'); disp(sprintf(' %% creating %s',dir_out)); mkdir(dir_out); 
 else disp(sprintf(' %% directory %s already exists, not creating.',dir_out)); end;
flag_trace_found = 0; 
tmpchar_trace = sprintf('%s/out_trace.txt',dir_out);
if exist(tmpchar_trace,'file');
tmp_trace = textread(tmpchar_trace);
if length(tmp_trace)> 6; disp(sprintf(' %% found %s of length %d, not rerunning.',tmpchar_trace,length(tmp_trace))); flag_trace_found = 1; end;
if (flag_force_create & length(tmp_trace)> 6); disp(sprintf(' %% found %s of length %d, actually, rerunning anyway.',tmpchar_trace,length(tmp_trace))); flag_trace_found = 0; end;
if length(tmp_trace)<=6; disp(sprintf(' %% found %s of length %d,     rerunning.',tmpchar_trace,length(tmp_trace))); flag_trace_found = 0; end;
end;%if exist(tmpchar_trace,'file');

if ~flag_trace_found;

dir_0in_plus_prefix = sprintf('%s/%s',dir_0in,str_prefix); 
dir_out_plus_prefix = sprintf('%s/%s',dir_out,str_out_suffix);

for nb1=0:nbins-1;
mr_M = zeros(size(M_n_{1+nb1},1),1);
mr_A_ori_{1+nb1} = mr_M; mr_A_ori_{1+nb1}(A_n_rij_{1+nb1})=1;
tmpchar = sprintf('%s_mr_A_%d.b16',dir_out_plus_prefix,0+nb1); 
if (flag_force_create | ~exist(tmpchar,'file'));
disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,length(mr_A_ori_{1+nb1}),1)); 
binary_compress(bitj,mr_A_ori_{1+nb1}(:)>0,tmpchar);
end;% if exist;
if Z_bother;
mr_Z_ori_{1+nb1} = mr_M; mr_Z_ori_{1+nb1}(Z_n_rij_{1+nb1})=1;
tmpchar = sprintf('%s_mr_Z_%d.b16',dir_out_plus_prefix,0+nb1); 
if (flag_force_create | ~exist(tmpchar,'file'));
disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,length(mr_Z_ori_{1+nb1}),1)); 
binary_compress(bitj,mr_Z_ori_{1+nb1}(:)>0,tmpchar);
end;% if exist;
end;%if Z_bother;
end;%for nb1=0:nbins-1;

if (shuffle_num>0); % performing covariate-respecting shuffle ;
if Z_bother;
[mr_A_prm_,mr_Z_prm_] = xxxcluster_uADZSZDA_shuffle_2(shuffle_num,M_n_,A_n_rij_,A_n_cij,Z_n_rij_,T_n_,T_n_cij);
for nb1=0:nbins-1;
tmpchar = sprintf('%s_mr_A_%d.b16',dir_out_plus_prefix,0+nb1); binary_compress(bitj,mr_A_prm_{1+nb1}(:)>0,tmpchar);
if Z_bother;
tmpchar = sprintf('%s_mr_Z_%d.b16',dir_out_plus_prefix,0+nb1); binary_compress(bitj,mr_Z_prm_{1+nb1}(:)>0,tmpchar);
end;%if Z_bother;
end;%for nb1=0:nbins-1;
end;%if Z_bother;
if ~Z_bother;
rng(shuffle_num);
for nb1=0:nbins-1;
[tmp_Q_,~] = qr(randn(size(M_n_{1+nb1},1)));
L_n_{1+nb1} = tmp_Q_*M_n_{1+nb1};
end;%for nb1=0:nbins-1;
end;%if ~Z_bother;
end;% if (shuffle_num>0); % performing covariate-respecting shuffle ;

% write M_n_ ;
for nb1=0:nbins-1;
tmpchar = sprintf('%s_M_%d_n.b16',dir_0in_plus_prefix,0+nb1); 
if (flag_force_create | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of length %d-x-%d',tmpchar,size(M_n_{1+nb1}))); binary_compress(bitj,M_n_{1+nb1}>0,tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
if ~Z_bother & shuffle_num>0;
tmpchar = sprintf('%s_L_%d_n.b16',dir_out_plus_prefix,0+nb1); 
if (flag_force_create | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of length %d-x-%d',tmpchar,size(L_n_{1+nb1}))); binary_compress(bitj,L_n_{1+nb1}>0,tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
end;%if ~Z_bother & shuffle_num>0;
tmpchar = sprintf('%s_M_%d_t.b16',dir_0in_plus_prefix,0+nb1); 
if (flag_force_create | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,size(transpose(M_n_{1+nb1})))); binary_compress(bitj,transpose(M_n_{1+nb1}>0),tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
if ~Z_bother & shuffle_num>0;
tmpchar = sprintf('%s_L_%d_t.b16',dir_out_plus_prefix,0+nb1); 
if (flag_force_create | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of length %d-x-%d',tmpchar,size(transpose(L_n_{1+nb1})))); binary_compress(bitj,transpose(L_n_{1+nb1}>0),tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
end;%if ~Z_bother & shuffle_num>0;
% write T_n_ ;
tmpchar = sprintf('%s_T_%d_n.b16',dir_0in_plus_prefix,0+nb1); 
if (flag_force_create | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,size(T_n_{1+nb1}))); binary_compress(bitj,T_n_{1+nb1}>0,tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
tmpchar = sprintf('%s_T_%d_t.b16',dir_0in_plus_prefix,0+nb1); 
if (flag_force_create | ~exist(tmpchar,'file')); disp(sprintf(' %% creating %s of size %d-x-%d',tmpchar,size(transpose(T_n_{1+nb1})))); binary_compress(bitj,transpose(T_n_{1+nb1}>0),tmpchar); else; disp(sprintf(' %% file %s already exists, not creating.',tmpchar)); end;
end;%for nb1=0:nbins-1;

% write A_n_cij and T_n_cij ;
mc_A = zeros(A_n_cols,1); mc_A(A_n_cij)=1;
tmpchar = sprintf('%s_mc_A.b16',dir_out_plus_prefix); binary_compress(bitj,mc_A(:)>0,tmpchar);
mc_T = zeros(T_n_cols,1); mc_T(T_n_cij)=1;
tmpchar = sprintf('%s_mc_T.b16',dir_out_plus_prefix); binary_compress(bitj,mc_T(:)>0,tmpchar);

fname_0in = sprintf('%s.in',dir_out_plus_prefix);
fp = fopen(fname_0in,'w');
fprintf(fp,'GLOBAL_verbose= %d;\n',flag_verbose);
fprintf(fp,'GLOBAL_thread_count= 8;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s%s;\n',str_lak_vs_dex,'cluster_driver');
fprintf(fp,'GLOBAL_QR_strategy= %s;\n',QR_strategy);
fprintf(fp,'GLOBAL_QC_strategy= %s;\n',QC_strategy);
fprintf(fp,'GLOBAL_NBINS= %d;\n',nbins);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',B_MLT);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_TEST_sparse= %d;\n',GLOBAL_TEST_sparse);
fprintf(fp,'GLOBAL_kappa_squared= %0.16f;\n',GLOBAL_kappa_squared);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
if Z_bother | shuffle_num==0;
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',dir_0in_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',dir_0in_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
end;%if Z_bother | shuffle_num==0;
if ~Z_bother & shuffle_num>0;
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_L_%d_n.b16',dir_out_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_L_%d_t.b16',dir_out_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
end;%if ~Z_bother & shuffle_num>0;
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',size(M_n_{1+nb1},1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',A_n_cols);
if (flag_reverse==1); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_Z_%d.b16',dir_out_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (flag_reverse==0); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_A_%d.b16',dir_out_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A.b16;\n',dir_out_plus_prefix);
if Z_bother;
fprintf(fp,'GLOBAL_Z_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',dir_0in_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',dir_0in_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',size(M_n_{1+nb1},1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
if (flag_reverse==1); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_A_%d.b16',dir_out_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (flag_reverse==0); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_Z_%d.b16',dir_out_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
end;%if Z_bother;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_cols); 
fprintf(fp,'GLOBAL_T_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_n.b16',dir_0in_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_t.b16',dir_0in_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T.b16;\n',dir_out_plus_prefix);
if Z_bother;
fprintf(fp,'GLOBAL_S_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_n.b16',dir_0in_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_S_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_t.b16',dir_0in_plus_prefix,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
end;%if Z_bother;
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by xxxcluster_uADZSZDA_dr_4.m on %s;\n',date);
fclose(fp);

flag_call=1*(slurm_walltime<=0);
if flag_call;
disp(sprintf('%s/%s < %s',dir_code,str_code,fname_0in));
system(sprintf('%s/%s < %s',dir_code,str_code,fname_0in));
end;%if flag_call;

flag_slurm = 1*(slurm_walltime>0);
if flag_slurm;
slurm_fname = sprintf('%s.slurm',dir_out_plus_prefix);
slurm_fp = fopen(slurm_fname,'w');
fprintf(slurm_fp,'#!/bin/sh \n');
fprintf(slurm_fp,'#\n');
fprintf(slurm_fp,'#SBATCH --verbose\n');
fprintf(slurm_fp,'#SBATCH --job-name=%s\n',fname_0in);
fprintf(slurm_fp,'#SBATCH --output=%s_output.log\n',dir_out_plus_prefix);
fprintf(slurm_fp,'#SBATCH --error=%s_error.log\n',dir_out_plus_prefix);
slurm_walltime_h = floor(slurm_walltime);
slurm_walltime_m = min(59,ceil(60*(slurm_walltime - slurm_walltime_h)));
sprintf(' %% slurm_walltime=%d:%.2d:59',slurm_walltime_h,slurm_walltime_m);
fprintf(slurm_fp,'#SBATCH --time=%d:%.2d:59\n',slurm_walltime_h,slurm_walltime_m);
fprintf(slurm_fp,'#SBATCH --nodes=%d --ntasks-per-node=%d\n',slurm_nnodes,slurm_tpn);
fprintf(slurm_fp,'#SBATCH --mem=%dGB\n',slurm_memdecl);
fprintf(slurm_fp,'\n');
fprintf(slurm_fp,'/bin/hostname\n');
fprintf(slurm_fp,'/bin/pwd\n');
fprintf(slurm_fp,'module load matlab/2016b\n');
fprintf(slurm_fp,'%s/%s < %s\n',dir_code,str_code,fname_0in);
fclose(slurm_fp);
type(slurm_fname);
system(sprintf('sbatch %s;\n',slurm_fname));
end;%if flag_slurm;

end;%if ~flag_trace_found;
