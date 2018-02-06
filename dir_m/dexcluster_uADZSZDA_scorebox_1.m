function dexcluster_uADZSZDA_scorebox_1(dir_code,dir_trunk,prefix,M_n_,rev_flag,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind,GLOBAL_TEST_sparse,gamma,B_MLT,Ireq,shuffle_num,verbose_flag,force_create_flag,scorebox_out_xdrop,scorebox_row_max,scorebox_row_num,scorebox_col_max,scorebox_col_num,slurm_walltime,slurm_nnodes,slurm_tpn,slurm_memdecl)
% using dexcluster_ver18 with slurm ; 
% test with: ;
%{
  dexcluster_uADZSZDA_scorebox_1();
  %}

if (nargin<1);
disp('returning'); return;
end;%if (nargin<1);

na=1;
if (nargin<na); dir_code = pwd; end; na=na+1;
if (nargin<na); dir_trunk = pwd; end; na=na+1;
if (nargin<na); prefix = 'test'; end; na=na+1;
if (nargin<na); M_n_{1} = randn(1024); end; na=na+1;
if (nargin<na); rev_flag = 0; end; na=na+1;
if (nargin<na); A_n_rind_ = {1:512}; end; na=na+1;
if (nargin<na); A_n_cind = 1:512; end; na=na+1;
if (nargin<na); Z_n_rind_ = {1:512}; end; na=na+1;
if (nargin<na); T_n_{1} = ones(1024,1); end; na=na+1;
if (nargin<na); T_n_cind = 1; end; na=na+1;
if (nargin<na); GLOBAL_TEST_sparse = 1; end; na=na+1;
if (nargin<na); gamma = 0; end; na=na+1;
if (nargin<na); B_MLT = 32; end; na=na+1;
if (nargin<na); Ireq = 1; end; na=na+1;
if (nargin<na); shuffle_num = 0; end; na=na+1;
if (nargin<na); verbose_flag = 0; end; na=na+1;
if (nargin<na); force_create_flag = 1; end; na=na+1;
if (nargin<na); slurm_walltime = 0; end; na=na+1;
if (nargin<na); slurm_nnodes = 1; end; na=na+1;
if (nargin<na); slurm_tpn = 1; end; na=na+1;
if (nargin<na); slurm_memdecl = 1; end; na=na+1;

nbins=length(M_n_);
[A_n_cols,Y_n_cols,T_n_cols,~,~] = lakcluster_uADZSZDA_check_0(shuffle_num,M_n_,A_n_rind_,A_n_cind,Z_n_rind_,T_n_,T_n_cind);

QR_strategy = 'YnWt';
QC_strategy = 'YnWt store one';
%QC_strategy = 'ZtSWn';
bitj = 16;

test_string = sprintf('%s_%s',prefix,dexcluster_uADZSZDA_xfix_gen_ver0(rev_flag,gamma,B_MLT,Ireq,shuffle_num));
dir__in = sprintf('%s/dir_%s',dir_trunk,prefix);
dir_out = sprintf('%s/dir_%s',dir__in,test_string); 
disp(sprintf(' test_string: %s',test_string));
disp(sprintf(' dir__in: %s',dir__in));
disp(sprintf(' dir_out: %s',dir_out));
if ~exist(dir_out,'dir'); disp(sprintf(' %% creating %s',dir_out)); mkdir(dir_out); 
 else disp(sprintf(' %% directory %s already exists, not creating.',dir_out)); end;
found_scorebox_flag = 0; 

scorebox_prefix = sprintf('%s/scorebox_r%dn%d_c%dn%d',dir_out,scorebox_row_max,scorebox_row_num,scorebox_col_max,scorebox_col_num);
scorebox_name = sprintf('%s_A_rpop_j_total.mda',scorebox_prefix); %S_A_rpop_j_total = mda_read_i4(scorebox_name);
if (exist(scorebox_name,'file')); disp(sprintf(' %% Found %s',scorebox_name)); found_scorebox_flag = found_scorebox_flag+1; end;
scorebox_name = sprintf('%s_A_cpop_j.mda',scorebox_prefix); %S_A_cpop_j = mda_read_i4(scorebox_name);
if (exist(scorebox_name,'file')); disp(sprintf(' %% Found %s',scorebox_name)); found_scorebox_flag = found_scorebox_flag+1; end;
scorebox_name = sprintf('%s_Irem.mda',scorebox_prefix); %S_Irem = mda_read_i4(scorebox_name);
if (exist(scorebox_name,'file')); disp(sprintf(' %% Found %s',scorebox_name)); found_scorebox_flag = found_scorebox_flag+1; end;
scorebox_name = sprintf('%s_QR_avg.mda',scorebox_prefix); %S_QR_avg = mda_read_r8(scorebox_name);
if (exist(scorebox_name,'file')); disp(sprintf(' %% Found %s',scorebox_name)); found_scorebox_flag = found_scorebox_flag+1; end;
scorebox_name = sprintf('%s_QC_avg.mda',scorebox_prefix); %S_QC_avg = mda_read_r8(scorebox_name);
if (exist(scorebox_name,'file')); disp(sprintf(' %% Found %s',scorebox_name)); found_scorebox_flag = found_scorebox_flag+1; end;
disp(sprintf(' %% Found %d/%d',found_scorebox_flag,5));
found_scorebox_flag=(found_scorebox_flag==5);
if found_scorebox_flag; disp(sprintf(' %% found scorebox, not rerunning.')); end;
if (force_create_flag & found_scorebox_flag); disp(sprintf(' %% found scorebox, rerunning anyway.')); found_scorebox_flag=0; end;
if (~found_scorebox_flag); disp(sprintf(' %% rerunning.')); found_scorebox_flag = 0; end;

if ~found_scorebox_flag;

d_inpre = sprintf('%s/%s',dir__in,prefix); 
d_oupre = sprintf('%s/%s',dir_out,test_string);

fname__in = sprintf('%s_scorebox.in',d_oupre);
fp = fopen(fname__in,'w');
fprintf(fp,'GLOBAL_verbose= %d;\n',verbose_flag);
fprintf(fp,'GLOBAL_thread_count= 8;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','dexcluster_scorebox');
fprintf(fp,'GLOBAL_NBINS= %d;\n',nbins);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',B_MLT);
fprintf(fp,'GLOBAL_TEST_sparse= %d;\n',GLOBAL_TEST_sparse);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',size(M_n_{1+nb1},1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',A_n_cols);
if (rev_flag==1); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_Z_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_A_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A.b16;\n',d_oupre);
fprintf(fp,'GLOBAL_Z_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',size(M_n_{1+nb1},1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
if (rev_flag==1); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_A_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_Z_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_cols); 
fprintf(fp,'GLOBAL_T_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T.b16;\n',d_oupre);
fprintf(fp,'GLOBAL_S_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_S_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_T_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_scorebox_out_xdrop= %s;\n',scorebox_out_xdrop);
fprintf(fp,'GLOBAL_scorebox_row_max= %d;\n',scorebox_row_max);
fprintf(fp,'GLOBAL_scorebox_row_num= %d;\n',scorebox_row_num);
fprintf(fp,'GLOBAL_scorebox_col_max= %d;\n',scorebox_col_max);
fprintf(fp,'GLOBAL_scorebox_col_num= %d;\n',scorebox_col_num);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by dexcluster_uADZSZDA_scorebox_1.m on %s;\n',date);
fclose(fp);

call_flag=1*(slurm_walltime<=0);
if call_flag;
disp(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
system(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
end;%if call_flag;

slurm_flag = 1*(slurm_walltime>0);
if slurm_flag;
slurm_fname = sprintf('%s.slurm',d_oupre);
slurm_fp = fopen(slurm_fname,'w');
fprintf(slurm_fp,'#!/bin/sh \n');
fprintf(slurm_fp,'#\n');
fprintf(slurm_fp,'#SBATCH --verbose\n');
fprintf(slurm_fp,'#SBATCH --job-name=%s\n',fname__in);
fprintf(slurm_fp,'#SBATCH --output=%s_output.log\n',d_oupre);
fprintf(slurm_fp,'#SBATCH --error=%s_error.log\n',d_oupre);
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
fprintf(slurm_fp,'%s/lakcluster_ver18 < %s\n',dir_code,fname__in);
fclose(slurm_fp);
type(slurm_fname);
system(sprintf('sbatch %s;\n',slurm_fname));
end;%if slurm_flag;

end;%if ~found_scorebox_flag;
