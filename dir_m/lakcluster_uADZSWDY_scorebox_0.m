function lakcluster_uADZSWDY_scorebox_0(dir_code,dir_trunk,prefix,M_n_,rev_flag,A_n_rind_,A_n_cind,Z_n_rind_,Y_n_cind,T_n_,T_n_cind,gamma,B_MLT,Ireq,shuffle_num,walltime,verbose_flag,force_create_flag,scorebox_out_xdrop,scorebox_row_max,scorebox_row_num,scorebox_col_max,scorebox_col_num)
% using lakcluster_ver18 ; 
% test with: ;
%{
  lakcluster_uADZSWDY_scorebox_0();
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
if (nargin<na); A_n_cind = 1:256; end; na=na+1;
if (nargin<na); Z_n_rind_ = {1:512}; end; na=na+1;
if (nargin<na); Y_n_rind_ = {257:512}; end; na=na+1;
if (nargin<na); T_n_{1} = ones(1024,1); end; na=na+1;
if (nargin<na); T_n_cind = 1; end; na=na+1;
if (nargin<na); gamma = 0; end; na=na+1;
if (nargin<na); B_MLT = 32; end; na=na+1;
if (nargin<na); Ireq = 1; end; na=na+1;
if (nargin<na); shuffle_num = 0; end; na=na+1;
if (nargin<na); walltime = 0; end; na=na+1;
if (nargin<na); verbose_flag = 0; end; na=na+1;
if (nargin<na); force_create_flag = 1; end; na=na+1;

nbins=length(M_n_);
[M_n_cols,A_n_cols,Y_n_cols,T_n_cols,~,~] = lakcluster_uADZSWDY_check_0(shuffle_num,M_n_,A_n_rind_,A_n_cind,Z_n_rind_,Y_n_cind,T_n_,T_n_cind);

QR_strategy = 'YnWt';
QC_strategy = 'YnWt store one';
%QC_strategy = 'ZtSWn';
bitj = 16;

test_string = sprintf('%s_%s',prefix,lakcluster_uADZSWDY_xfix_gen_ver0(rev_flag,gamma,B_MLT,Ireq,shuffle_num));
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
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','lakcluster_scorebox');
fprintf(fp,'GLOBAL_QR_strategy= %s;\n',QR_strategy);
fprintf(fp,'GLOBAL_QC_strategy= %s;\n',QC_strategy);
fprintf(fp,'GLOBAL_NBINS= %d;\n',nbins);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',B_MLT);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
fprintf(fp,'GLOBAL_A_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',size(M_n_{1+nb1},1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',M_n_cols);
if (rev_flag==1); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_Z_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_A_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A.b16;\n',d_oupre);
fprintf(fp,'GLOBAL_Z_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb1=0:nbins-1; fprintf(fp,'%d',size(M_n_{1+nb1},1)); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
if (rev_flag==1); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_A_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb1=0:nbins-1; fprintf(fp,'%s_mr_Z_%d.b16',d_oupre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_Y_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Y_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',M_n_cols);
fprintf(fp,'GLOBAL_Y_n_cind= %s_mc_Y.b16;\n',d_oupre);
fprintf(fp,'GLOBAL_W_n_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_n.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_W_t_name_= '); for nb1=0:nbins-1; fprintf(fp,'%s_M_%d_t.b16',d_inpre,0+nb1); if nb1<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
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
fprintf(fp,'%% generated by lakcluster_uADZSWDY_scorebox_0.m on %s;\n',date);
fclose(fp);

call_flag=1*(walltime<=0);
if call_flag;
disp(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
system(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
end;%if call_flag;

pbs_flag=1*(walltime>0);
if pbs_flag;
fname_pbs = sprintf('%s.pbs',d_oupre);
if ~found_scorebox_flag;
fp = fopen(fname_pbs,'w');
fprintf(fp,'#PBS -S /bin/bash\n');
fprintf(fp,'#PBS -lnodes=1:ppn=15\n');
%fprintf(fp,'#PBS -lnodes=1:cpu2+cpu3\n');
fprintf(fp,'#PBS -lwalltime=%d:00:00\n',walltime * (~found_scorebox_flag));
fprintf(fp,'#PBS -e %s_pbs_error.log\n',d_oupre);
fprintf(fp,'#PBS -o %s_pbs_output.log\n',d_oupre);
if ~found_scorebox_flag;
fprintf(fp,'%s/lakcluster_ver18 < %s\n',dir_code,fname__in);
end;%if ~found_scorebox_flag;
fclose(fp);
type(fname_pbs);
system(sprintf('qsub %s;\n',fname_pbs));
end;%if ~found_scorebox_flag;
end;%if pbs_flag;

end;%if ~found_scorebox_flag;
