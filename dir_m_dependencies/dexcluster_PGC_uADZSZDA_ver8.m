function dexcluster_PGC_uADZSZDA_ver8(dir_code,dir_trunk,prefix,p_threshold,rev_flag,nbins,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle,pbs_walltime,row_factor,col_factor,verbose_flag)
% set up to allow pbs_walltime as input ; using lakcluster_ver18 ; no replication ; reverse-search ;
% try: ;
%{

  !scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_lakcluster_c/dir_m/?*.m ~/dir_lakcluster_c/dir_m/ ;
  !cd ~/dir_lakcluster_c; make lakcluster_ver18; cd ~/dir_PGC_20180304/dir_m ;
  !scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_PGC_20180304/dir_m/dexcluster_PGC_uADZSZDA_ver8.m ~/dir_PGC_20180304/dir_m/ ;
  !scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_PGC_20180304/dir_m/xxxcluster_PGC_uADZSZDA_xfix_ver8.m ~/dir_PGC_20180304/dir_m/ ;
  !scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_PGC_20180304/dir_m/xxxcluster_PGC_getbim_ver2.m ~/dir_PGC_20180304/dir_m/ ;
  !scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_PGC_20180304/dir_m/xxxcluster_uADZSZDA_check_ver1.m ~/dir_PGC_20180304/dir_m/ ;
  !scp -p rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_PGC_20180304/dir_m/xxxcluster_uADZSZDA_shuffle_ver2.m ~/dir_PGC_20180304/dir_m/ ;
 
  % interactive submit with something like: ;
  % ls ~/dir_PGC_20180304/dir_PGC_cl3_maf01_analyze/dir_*?/?*.pbs | xargs -p -n 1 qsub ;
  % submit with something like: ;
  % ls ~/dir_PGC_20180304/dir_PGC_cl3_maf01_analyze/dir_*?/?*.pbs | xargs -n 1 qsub ;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % trunk ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  dir_trunk = '/home/arangan/dir_PGC_20180304';
  dir_code = sprintf('/home/arangan/dir_lakcluster_c');
  %dir_trunk = '/data/rangan/dir_bcc/dir_PGC_20180304';
  %dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % test ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  p_threshold_ = [0.01 , 0.05 , 0.10 , 0.25]; gamma=0.004;B_MLT=34;
  pbs_walltime=0;row_factor=1.0;col_factor=1.0; 
  verbose_flag=0;
  cl_num=3;%cl_num = input('cluster number? [1,2,3]');
  if cl_num==1; prefix = 'PGC_cl1_maf01'; nbins = 9; Ireq = floor(nbins/2); end;
  if cl_num==2; prefix = 'PGC_cl2_maf01'; nbins = 7; Ireq = floor(nbins/2); end;
  if cl_num==3; prefix = 'PGC_cl3_maf01'; nbins = 10; Ireq = floor(nbins/2); end;
  for p_threshold = p_threshold_;
  for rev_flag=[0,1]; 
  for mds_tab=[1]; 
  n_mds = 20; if mds_tab==1; mds_used_ = [1:2]; mds_repl = 1; else mds_used_ = []; mds_repl = 1; end; 
  nshuffle=0; dexcluster_PGC_uADZSZDA_ver8(dir_code,dir_trunk,prefix,p_threshold,rev_flag,nbins,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle,pbs_walltime,row_factor,col_factor,verbose_flag);
  for nshuffle=1:8; dexcluster_PGC_uADZSZDA_ver8(dir_code,dir_trunk,prefix,p_threshold,rev_flag,nbins,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle,pbs_walltime,row_factor,col_factor,verbose_flag); end;
  end;%for mds_tab=[1]; 
  end;%for rev_flag=[0,1]; 
  end;%for p_threshold = p_threshold_;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % run ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  p_threshold_ = [0.01 , 0.05 , 0.10 , 0.25]; gamma=0.004;B_MLT=34;
  %pbs_walltime_ = [1.50 , 1.25 , 1.00 , 0.75];
  pbs_walltime_ = [4.50 , 4.25 , 4.00 , 3.75];
  row_factor=1.0;col_factor=1.0; 
  verbose_flag=0;
  cl_num=3;%cl_num = input('cluster number? [1,2,3]');
  if cl_num==1; prefix = 'PGC_cl1_maf01'; nbins = 9; Ireq = floor(nbins/2); end;
  if cl_num==2; prefix = 'PGC_cl2_maf01'; nbins = 7; Ireq = floor(nbins/2); end;
  if cl_num==3; prefix = 'PGC_cl3_maf01'; nbins = 10; Ireq = floor(nbins/2); end;
  for n_p=1:length(p_threshold_);
  p_threshold = p_threshold_(n_p); pbs_walltime = pbs_walltime_(n_p);
  for rev_flag = [0:1]; for mds_tab = [1];
  n_mds = 20; if mds_tab==1; mds_used_ = [1:2]; mds_repl = 1; else mds_used_ = []; mds_repl = 1; end; 
  nshuffle=0; dexcluster_PGC_uADZSZDA_ver8(dir_code,dir_trunk,prefix,p_threshold,rev_flag,nbins,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle,pbs_walltime,row_factor,col_factor,verbose_flag);
  for nshuffle=1:8; dexcluster_PGC_uADZSZDA_ver8(dir_code,dir_trunk,prefix,p_threshold,rev_flag,nbins,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle,pbs_walltime,row_factor,col_factor,verbose_flag); end;
  end;end;%for mds_tab = [1];%for rev_flag = [0:1]; 
  end;%  for n_p = 1:length(p_threshold_);

  %}

na=1;
if (nargin<na); dir_code = pwd; end; na=na+1;
if (nargin<na); dir_trunk = pwd; end; na=na+1;
if (nargin<na); prefix = 'clx'; end; na=na+1;
if (nargin<na); p_threshold = 0; end; na=na+1;
if (nargin<na); rev_flag = 0; end; na=na+1;
if (nargin<na); nbins = 1; end; na=na+1;
if (nargin<na); n_mds = 20; end; na=na+1;
if (nargin<na); mds_used_ = [1:2]; end; na=na+1;
if (nargin<na); mds_repl = 1; end; na=na+1;
if (nargin<na); gamma = 0.002; end; na=na+1;
if (nargin<na); B_MLT = 8; end; na=na+1;
if (nargin<na); Ireq = 0; end; na=na+1;
if (nargin<na); nshuffle = 0; end; na=na+1;
if (nargin<na); pbs_walltime = 13; end; na=na+1;
if (nargin<na); row_factor = 1.0; end; na=na+1;
if (nargin<na); col_factor = 1.0; end; na=na+1;
if (nargin<na); verbose_flag = 0; end; na=na+1;

bitj = 16;

name_string = sprintf('%s_%s',prefix,xxxcluster_PGC_uADZSZDA_xfix_ver8('dex',p_threshold,rev_flag,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle));
name_suffix = sprintf('%s','analyze');
disp(sprintf(' name_string: %s',name_string));
dir__in = sprintf('%s/dir_%s',dir_trunk,prefix);
dir_tmp = sprintf('%s_%s',dir__in,name_suffix); if ~exist(dir_tmp,'dir'); mkdir(dir_tmp); end;
dir_out = sprintf('%s_%s/dir_%s',dir__in,name_suffix,name_string); if ~exist(dir_out,'dir'); mkdir(dir_out); end;
found_trace_flag = 0; 
tmpchar = sprintf('%s/out_trace.txt',dir_out);
if exist(tmpchar,'file');
tmp = textread(tmpchar);
if length(tmp)> 970; disp(sprintf(' %% found %s of length %d, not rerunning.',tmpchar,length(tmp))); found_trace_flag = 1; end;
%if length(tmp)> 970; disp(sprintf(' %% found %s of length %d, actually, rerunning anyway.',tmpchar,length(tmp))); found_trace_flag = 0; end;
if length(tmp)<=970; disp(sprintf(' %% found %s of length %d,     rerunning.',tmpchar,length(tmp))); found_trace_flag = 0; end;
end;%if exist(tmpchar,'file');

if ~found_trace_flag

d_inpre = sprintf('%s/%s',dir__in,prefix); 
d_oupre = sprintf('%s/%s',dir_out,name_string);

Y_n_cols=0;

% reading original row-masks for A and Z ;
mr_A_ori_ = cell(nbins,1);
mr_Z_ori_ = cell(nbins,1);
A_n_rind_ = cell(nbins,1);
Z_n_rind_ = cell(nbins,1);
for nb=0:nbins-1;
tmpchar = sprintf('%s_mr_A_%0.2d.b16',d_inpre,1+nb);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); 
mr_A_ori_{1+nb} = (tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols)>0);
A_n_rind_{1+nb} = find(mr_A_ori_{1+nb});
tmpchar = sprintf('%s_mr_Z_%0.2d.b16',d_inpre,1+nb);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); 
mr_Z_ori_{1+nb} = (tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols)>0);
Z_n_rind_{1+nb} = find(mr_Z_ori_{1+nb});
end;%for nb=0:nbins-1;

% reading original col-masks for A ;
tmpchar = sprintf('%s_mc_A.b16',d_inpre);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); 
mc_A = (tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols)>0);
A_n_cind = find(mc_A);

M_n_rows_ = zeros(nbins,1);
A_n_rows_ = zeros(nbins,1);
Z_n_rows_ = zeros(nbins,1);
for nb=0:nbins-1;
M_n_rows_(1+nb) = length(mr_A_ori_{1+nb});
A_n_rows_(1+nb) = sum(mr_A_ori_{1+nb});
Z_n_rows_(1+nb) = sum(mr_Z_ori_{1+nb});
end;%for nb=0:nbins-1;
M_n_cols = length(mc_A);
disp(sprintf(' %% nbins %d; total vs cases vs controls',nbins));
disp(num2str([M_n_rows_ , A_n_rows_ , Z_n_rows_]));

% compressing T_n_ ;
T_n_ = cell(nbins,1);
T_n_crop_ = cell(nbins,1);
T_t_crop_ = cell(nbins,1);
T_n_crop_cols = 1+length(mds_used_);
for nb=0:nbins-1;
tmpchar = sprintf('%s_T_%0.2d_n.b16',d_inpre,1+nb); [bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); disp(sprintf(' %% reading %s = (%d,%d)',tmpchar,nrows,ncols));
T_n_{1+nb} = tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols); 
T_n_crop_{1+nb} = T_n_{1+nb}(:,[1,(1+(mds_used_))]); 
T_t_crop_{1+nb} = transpose(T_n_crop_{1+nb});
tmpchar = sprintf('%s_T_crop_%0.2d_n.b16',d_oupre,1+nb); disp(sprintf(' %% writing %s = (%d,%d)',tmpchar,size(T_n_crop_{1+nb})));
tutorial_binary_compress(bitj,T_n_crop_{1+nb}>0,tmpchar);
tmpchar = sprintf('%s_T_crop_%0.2d_t.b16',d_oupre,1+nb); disp(sprintf(' %% writing %s = (%d,%d)',tmpchar,size(T_t_crop_{1+nb})));
tutorial_binary_compress(bitj,T_t_crop_{1+nb}>0,tmpchar);
end;%for nb=0:nbins-1;
mc_T_crop = ones(T_n_crop_cols,1);
tmpchar = sprintf('%s_mc_T_crop.b16',d_oupre); disp(sprintf(' %% writing %s = (%d,%d)',tmpchar,size(mc_T_crop))); tutorial_binary_compress(bitj,mc_T_crop(:)>0,tmpchar);
disp(sprintf('mc_T_crop:'));disp(num2str(transpose(mc_T_crop)));
T_n_crop_cind = 1:T_n_crop_cols;

% checking consistency ;
flag_error = xxxcluster_uADZSZDA_check_ver1(M_n_rows_,M_n_cols,A_n_rind_,A_n_cind,Z_n_rind_,T_n_crop_cols,T_n_crop_,T_n_crop_cind); 
if (flag_error); disp(sprintf(' %% Warning! incorrect dimensions in dexcluster_PGC_uADZSZDA_ver8')); return; end;

% performing covariate-respecting shuffle ;
[mr_A_prm_,mr_Z_prm_] = xxxcluster_uADZSZDA_shuffle_ver2(nshuffle,M_n_rows_,M_n_cols,A_n_rind_,A_n_cind,Z_n_rind_,T_n_crop_cols,T_n_crop_,T_n_crop_cind);

% writing row-masks ;
A_n_rows_used=0;
Z_n_rows_used=0;
mr_A_use_ = mr_A_ori_; mr_Z_use_ = mr_Z_ori_;
if (nshuffle>0); mr_A_use_ = mr_A_prm_; mr_Z_use_ = mr_Z_prm_; end;%if (nshuffle>0); 
for nb=0:nbins-1;
if (row_factor<1); mr_A_use_{1+nb} = mr_A_use_{1+nb}.*(rand(size(mr_A_use_{1+nb}))<row_factor); mr_Z_use_{1+nb} = mr_Z_use_{1+nb}.*(rand(size(mr_Z_use_{1+nb}))<row_factor); end;%if (row_factor<1);
disp(sprintf('nb %.2d : mr_A_ori_ npats %.5d ncase %.4d, mr_A_use_ npats %.5d ncase %.4d, overlap %.4d',nb,length(mr_A_ori_{1+nb}),sum(mr_A_ori_{1+nb}),length(mr_A_use_{1+nb}),sum(mr_A_use_{1+nb}),sum(mr_A_ori_{1+nb}.*mr_A_use_{1+nb})));
disp(sprintf('nb %.2d : mr_Z_ori_ npats %.5d nctrl %.4d, mr_Z_use_ npats %.5d nctrl %.4d, overlap %.4d',nb,length(mr_Z_ori_{1+nb}),sum(mr_Z_ori_{1+nb}),length(mr_Z_use_{1+nb}),sum(mr_Z_use_{1+nb}),sum(mr_Z_ori_{1+nb}.*mr_Z_use_{1+nb})));
tmpchar = sprintf('%s_mr_A_%0.2d.b16',d_oupre,1+nb);tutorial_binary_compress(bitj,mr_A_use_{1+nb}(:)>0,tmpchar);
tmpchar = sprintf('%s_mr_Z_%0.2d.b16',d_oupre,1+nb);tutorial_binary_compress(bitj,mr_Z_use_{1+nb}(:)>0,tmpchar);
A_n_rows_used = A_n_rows_used + sum(mr_A_use_{1+nb}(:)>0);
Z_n_rows_used = Z_n_rows_used + sum(mr_Z_use_{1+nb}(:)>0);
end;%for nb=0:nbins-1;

% writing col-mask ;
mc_A_use = ones(M_n_cols,1);
fname_bim = sprintf('%s/%s_bim.ext',dir__in,prefix);
mc_A_bim = xxxcluster_PGC_getbim_ver2(fname_bim,p_threshold,p_threshold^2);
mc_A_use = mc_A_use .* mc_A_bim;
disp(sprintf(' %% p_threshold %0.2f, retaining %d mc values, but setting %d mc values to 0',p_threshold,sum(mc_A_use),sum(~mc_A_use)));
if (col_factor<1);
mc_A_use = mc_A_use.*(rand(size(mc_A_use))<col_factor); 
disp(sprintf(' %% col_factor %0.2f, retaining %d mc values, but setting %d mc values to 0',col_factor,sum(mc_A_use),sum(~mc_A_use)));
end;%if (col_factor<1);
fname_mc_A = sprintf('%s_mc_A.b16',d_oupre);
disp(sprintf(' %% creating %s',fname_mc_A));
tmpchar_out = fname_mc_A;tutorial_binary_compress(bitj,mc_A_use(:)>0,tmpchar_out); 
disp(sprintf('mc_A nsnps %.9d/%.9d',sum(mc_A_use(:)>0),length(mc_A_use)));
A_n_cols_used = sum(mc_A_use(:)>0);

disp(sprintf(' %% A_n_rows_used %d Z_n_rows_used %d A_n_cols_used %d',A_n_rows_used,Z_n_rows_used,A_n_cols_used));

fname__in = sprintf('%s.in',d_oupre);
fp = fopen(fname__in,'w');
fprintf(fp,'GLOBAL_verbose= %d;\n',verbose_flag);
fprintf(fp,'GLOBAL_thread_count= 15;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','dexcluster_driver');
fprintf(fp,'GLOBAL_NBINS= %d;\n',nbins);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',B_MLT);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
fprintf(fp,'GLOBAL_A_n_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_A_%0.2d_n.b16',d_inpre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_A_%0.2d_t.b16',d_inpre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_rows_= '); for nb=0:nbins-1; fprintf(fp,'%d',M_n_rows_(1+nb)); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',M_n_cols);
if (rev_flag==1); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb=0:nbins-1; fprintf(fp,'%s_mr_Z_%0.2d.b16',d_oupre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb=0:nbins-1; fprintf(fp,'%s_mr_A_%0.2d.b16',d_oupre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A.b16;\n',d_oupre);
fprintf(fp,'GLOBAL_Z_n_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_A_%0.2d_n.b16',d_inpre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_t_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_A_%0.2d_t.b16',d_inpre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb=0:nbins-1; fprintf(fp,'%d',M_n_rows_(1+nb)); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
if (rev_flag==1); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb=0:nbins-1; fprintf(fp,'%s_mr_A_%0.2d.b16',d_oupre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb=0:nbins-1; fprintf(fp,'%s_mr_Z_%0.2d.b16',d_oupre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_crop_cols); 
fprintf(fp,'GLOBAL_T_n_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_T_crop_%0.2d_n.b16',d_oupre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_T_crop_%0.2d_t.b16',d_oupre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T_crop.b16;\n',d_oupre);
fprintf(fp,'GLOBAL_S_n_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_T_crop_%0.2d_n.b16',d_oupre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_S_t_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_T_crop_%0.2d_t.b16',d_oupre,1+nb); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by dexcluster_PGC_uADZSZDA_ver8.m on %s;\n',date);
fclose(fp);
type(fname__in);

call_flag=1*(pbs_walltime<=0);%call_flag = input(' call? 1=yes (default), 0=no:'); if isempty(call_flag); call_flag=1; end;
if call_flag;
disp(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
system(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
end;%if call_flag;

pbs_flag=1*(pbs_walltime>0);
if pbs_flag;
fname_pbs = sprintf('%s.pbs',d_oupre);
fp = fopen(fname_pbs,'w');
fprintf(fp,'#PBS -S /bin/bash\n');
fprintf(fp,'#PBS -lnodes=1:ppn=15\n');
%fprintf(fp,'#PBS -lnodes=1:cpu2+cpu3\n');
pbs_walltime_h = floor(pbs_walltime); pbs_walltime_m = min(59,ceil(60*(pbs_walltime - pbs_walltime_h))); 
sprintf(' %% pbs_walltime=%d:%.2d:59',pbs_walltime_h,pbs_walltime_m);
fprintf(fp,'#PBS -lwalltime=%d:%.2d:59\n',pbs_walltime_h,pbs_walltime_m);
fprintf(fp,'#PBS -e %s_pbs_error.log\n',d_oupre);
fprintf(fp,'#PBS -o %s_pbs_output.log\n',d_oupre);
fprintf(fp,'%s/lakcluster_ver18 < %s\n',dir_code,fname__in);
fclose(fp);
type(fname_pbs);
system(sprintf('qsub %s;\n',fname_pbs));
end;%if pbs_flag;

end;%if ~found_trace_flag;

