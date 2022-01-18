function dexcluster_PGC_ADZSZDA_ver7(dir_code,dir_trunk,prefix,p_threshold,rev_flag,n_mds,mds_used,gamma,B_MLT,Ireq,nshuffle,pbs_walltime,row_factor,col_factor,verbose_flag)
% set up to allow pbs_walltime as input ; using lakcluster_ver18 ; no replication ; reverse-search ;
% try: ;
%{

  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_lakcluster_c ~/ ;
  !cd ~/dir_lakcluster_c; make lakcluster_ver18; cd ~/dir_PGC_20180304/dir_m ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_lakcluster_c/dir_m/dexcluster_PGC_ADZSZDA_ver7.m ~/dir_lakcluster_c/dir_m/ ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_lakcluster_c/dir_m/dexcluster_PGC_xfix_gen_ver7.m ~/dir_lakcluster_c/dir_m/ ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_lakcluster_c/dir_m/lakcluster_PGC_getbim_ver1.m ~/dir_lakcluster_c/dir_m/ ;
  
  % interactive submit with something like: ;
  % ls ~/dir_PGC_20180304/dir_PGC_cl3_maf01_analyze/dir_*?/?*.pbs | xargs -p -n 1 qsub ;
  % submit with something like: ;
  % ls ~/dir_PGC_20180304/dir_PGC_cl3_maf01_analyze/dir_*?/?*.pbs | xargs -n 1 qsub ;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % test ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  dir_code = sprintf('/home/arangan/dir_lakcluster_c');dir_trunk = sprintf('/home/arangan/dir_PGC_20180304'); 
  p_threshold_ = [0.01 , 0.05 , 0.10 , 0.25]; gamma=0.004;B_MLT=34;
  pbs_walltime=0;row_factor=1.0;col_factor=1.0; 
  verbose_flag=0;
  cl_num=3;%cl_num = input('cluster number? [1,2,3]');
  if cl_num==1; prefix = 'PGC_cl1_maf01'; nbins = 11; Ireq = floor(nbins/2); end;
  if cl_num==2; prefix = 'PGC_cl2_maf01'; nbins = 7; Ireq = floor(nbins/2); end;
  if cl_num==3; prefix = 'PGC_cl3_maf01'; nbins = 10; Ireq = floor(nbins/2); end;
  Ireq = 0;
  for p_threshold = p_threshold_;
  rev_flag=1; mds_tab=1; n_mds = 20; if mds_tab==1; mds_used = [1:2]; else mds_used = []; end; 
  nshuffle=0; dexcluster_PGC_ADZSZDA_ver7(dir_code,dir_trunk,prefix,p_threshold,rev_flag,n_mds,mds_used,gamma,B_MLT,0,nshuffle,pbs_walltime,row_factor,col_factor,verbose_flag);
  end;%for p_threshold = p_threshold_;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % run ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  dir_code = sprintf('/home/arangan/dir_lakcluster_c');dir_trunk = sprintf('/home/arangan/dir_PGC_20180304'); 
  p_threshold_ = [0.01 , 0.05 , 0.10 , 0.25]; gamma=0.004;B_MLT=34;
  pbs_walltime_ = [1.50 , 1.25 , 1.00 , 0.75];
  row_factor=1.0;col_factor=1.0; 
  verbose_flag=0;
  cl_num=3;%cl_num = input('cluster number? [1,2,3]');
  if cl_num==1; prefix = 'PGC_cl1_maf01'; nbins = 11; Ireq = floor(nbins/2); end;
  if cl_num==2; prefix = 'PGC_cl2_maf01'; nbins = 7; Ireq = floor(nbins/2); end;
  if cl_num==3; prefix = 'PGC_cl3_maf01'; nbins = 10; Ireq = floor(nbins/2); end;
  Ireq = 0;
  for n_p=1:length(p_threshold_);
  p_threshold = p_threshold_(n_p); pbs_walltime = pbs_walltime_(n_p);
  for rev_flag = [0:1]; for mds_tab = [1];
  n_mds = 20; if mds_tab==1; mds_used = [1:2]; else mds_used = []; end; 
  nshuffle=0; dexcluster_PGC_ADZSZDA_ver7(dir_code,dir_trunk,prefix,p_threshold,rev_flag,n_mds,mds_used,gamma,B_MLT,0,nshuffle,pbs_walltime,row_factor,col_factor,verbose_flag);
  for nshuffle=1:8; dexcluster_PGC_ADZSZDA_ver7(dir_code,dir_trunk,prefix,p_threshold,rev_flag,n_mds,mds_used,gamma,B_MLT,0,nshuffle,pbs_walltime,row_factor,col_factor,verbose_flag); end;
  end;end;%for mds_tab = 0;%for rev_flag = [0:1]; 
  end;%  for n_p = 1:length(p_threshold_);

  %}

na=1;
if (nargin<na); dir_trunk = pwd; end; na=na+1;
if (nargin<na); prefix = 'clx'; end; na=na+1;
if (nargin<na); p_threshold = 0; end; na=na+1;
if (nargin<na); rev_flag = 0; end; na=na+1;
if (nargin<na); n_mds = 20; end; na=na+1;
if (nargin<na); mds_used = [1:2]; end; na=na+1;
if (nargin<na); gamma = 0.002; end; na=na+1;
if (nargin<na); B_MLT = 8; end; na=na+1;
if (nargin<na); Ireq = 0; end; na=na+1;
if (nargin<na); nshuffle = 0; end; na=na+1;
if (nargin<na); pbs_walltime = 13; end; na=na+1;
if (nargin<na); row_factor = 1.0; end; na=na+1;
if (nargin<na); col_factor = 1.0; end; na=na+1;
if (nargin<na); verbose_flag = 0; end; na=na+1;

nbins=1;
mrnd = 0;
bitj = 16;

name_string = sprintf('%s_%s',prefix,dexcluster_PGC_xfix_gen_ver7(p_threshold,rev_flag,n_mds,mds_used,gamma,B_MLT,Ireq,nshuffle));
name_suffix = sprintf('%s','analyze');
disp(sprintf(' name_string: %s',name_string));
dir__in = sprintf('%s/dir_%s',dir_trunk,prefix);
dir_tmp = sprintf('%s_%s',dir__in,name_suffix); if ~exist(dir_tmp,'dir'); mkdir(dir_tmp); end;
dir_out = sprintf('%s_%s/dir_%s',dir__in,name_suffix,name_string); if ~exist(dir_out,'dir'); mkdir(dir_out); end;
found_trace_flag = 0; 
tmpchar = sprintf('%s/out_trace.txt',dir_out);
if exist(tmpchar,'file');
tmp = textread(tmpchar);
%if length(tmp)> 970; disp(sprintf(' %% found %s of length %d, not rerunning.',tmpchar,length(tmp))); found_trace_flag = 1; end;
if length(tmp)> 970; disp(sprintf(' %% found %s of length %d, actually, rerunning anyway.',tmpchar,length(tmp))); found_trace_flag = 0; end;
if length(tmp)<=970; disp(sprintf(' %% found %s of length %d,     rerunning.',tmpchar,length(tmp))); found_trace_flag = 0; end;
end;%if exist(tmpchar,'file');

if ~found_trace_flag

d_inpre = sprintf('%s/%s',dir__in,prefix); 
d_oupre = sprintf('%s/%s',dir_out,name_string);

ncols_Y=0; %ncols_T=21;

nb=1; % only 1 bin ;
% A and Z;
tmpchar = sprintf('%s_mr_A_full.b16',d_inpre);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); 
mr_A = (tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols)>0);
tmpchar = sprintf('%s_mr_Z_full.b16',d_inpre);[bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); 
mr_Z = (tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols)>0);
tmp_A_rij = find(mr_A); tmp_Z_rij = find(mr_Z);

if (nshuffle>0); % performing mds-sector-respecting shuffle ;
rng(nshuffle); 
mr_A_prm = zeros(size(mr_A)); mr_Z_prm = zeros(size(mr_Z));
nb=1; tmpchar = sprintf('%s_T_full_n.b16',d_inpre); [bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); tmp_T = tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols); 
disp(sprintf(' %% reading %s = (%d,%d)',tmpchar,nrows,ncols));
tmp_T = tmp_T(:,1 + (mds_used))>0; nsec = 2^size(tmp_T,2);
tmp_T = tmp_T * transpose(2.^[0:size(tmp_T,2)-1]);
for ns = 1:nsec; sec_rij_{ns} = find(tmp_T==ns-1); end;
for ns = 1:nsec; 
lA_(ns) = length(intersect(tmp_A_rij,sec_rij_{ns})); lZ_(ns) = length(intersect(tmp_Z_rij,sec_rij_{ns}));
prm_{ns} = randperm(length(sec_rij_{ns}));
prm_A_{ns} = sec_rij_{ns}(prm_{ns}(1:lA_(ns))); prm_Z_{ns} = sec_rij_{ns}(prm_{ns}(lA_(ns) + (1:lZ_(ns))));
end;%for ns=1:nsec;
mr_A_prm = zeros(size(mr_A)); mr_Z_prm = zeros(size(mr_Z));
for ns=1:nsec;
mr_A_prm(prm_A_{ns}) = 1; mr_Z_prm(prm_Z_{ns}) = 1;
end;%for ns=1:nsec;
disp_flag=0;
if disp_flag;
tmp_A_prm = find(mr_A_prm); tmp_Z_prm = find(mr_Z_prm);
for ns=1:nsec; lA_prm_(ns) = length(intersect(tmp_A_prm,sec_rij_{ns})); lZ_prm_(ns) = length(intersect(tmp_Z_prm,sec_rij_{ns})); end;
for ns=1:nsec; lA_cap_(ns) = length(intersect(intersect(tmp_A_rij,sec_rij_{ns}),intersect(tmp_A_prm,sec_rij_{ns}))); lZ_cap_(ns) = length(intersect(intersect(tmp_Z_rij,sec_rij_{ns}),intersect(tmp_Z_prm,sec_rij_{ns}))); end;
disp([lA_ ; lA_prm_ ; lA_cap_]); disp([lZ_ ; lZ_prm_ ; lZ_cap_]); 
disp(sprintf('mr_A %d --> mr_A_prm %d; mr_Z %d --> mr_Z_prm %d',sum(mr_A),sum(mr_A_prm),sum(mr_Z),sum(mr_Z_prm)));
end;% if disp_flag;
end;%if (nshuffle>0); 

mr_A_tmp = mr_A; nrows_A_(nb) = length(mr_A_tmp);
mr_Z_tmp = mr_Z; nrows_Z_(nb) = length(mr_Z_tmp);
if (nshuffle>0); 
disp(sprintf(' %% nshuffle %d... shuffling: ',nshuffle));
mr_A_tmp = mr_A_prm;
mr_Z_tmp = mr_Z_prm;
end;%if (nshuffle>0); 
if (row_factor<1); mr_A_tmp = mr_A_tmp.*(rand(size(mr_A_tmp))<row_factor); end;%if (row_factor<1);
disp(sprintf('nb %.2d : mr_A npats %.5d ncase %.4d, mr_A_tmp npats %.5d ncase %.4d, overlap %.4d',nb,length(mr_A),sum(mr_A),length(mr_A_tmp),sum(mr_A_tmp),sum(mr_A.*mr_A_tmp)));
disp(sprintf('nb %.2d : mr_Z npats %.5d nctrl %.4d, mr_Z_tmp npats %.5d nctrl %.4d, overlap %.4d',nb,length(mr_Z),sum(mr_Z),length(mr_Z_tmp),sum(mr_Z_tmp),sum(mr_Z.*mr_Z_tmp)));
if ~found_trace_flag;
tmpchar = sprintf('%s_mr_A_full.b16',d_oupre);tutorial_binary_compress(bitj,mr_A_tmp(:)>0,tmpchar);
tmpchar = sprintf('%s_mr_Z_full.b16',d_oupre);tutorial_binary_compress(bitj,mr_Z_tmp(:)>0,tmpchar);
end;%if ~found_trace_flag;
rows_used = sum(mr_A_tmp(:)>0);

% A ;
tmpchar__in = sprintf('%s_mc_A.b16',d_inpre); [bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar__in); ncols_A = max(nrows,ncols);
mc_A_tmp = ones(ncols_A,1); if (col_factor<1); mc_A_tmp = mc_A_tmp.*(rand(size(mc_A_tmp))<col_factor); end;%if (col_factor<1);
fname_bim = sprintf('%s_bim.ext',d_inpre);
fname_mc_A_pre = sprintf('%s_bim_p%.2d_mc_A.b16',d_inpre,floor(100*p_threshold));
fname_mc_A_pos = sprintf('%s_mc_A.b16',d_oupre);
if exist(fname_mc_A_pre,'file');
disp(sprintf(' %% copying %s to %s',fname_mc_A_pre,fname_mc_A_pos));
copyfile(fname_mc_A_pre,fname_mc_A_pos);
[bitj,nr_tmp,nc_tmp] = tutorial_binary_getsize(fname_mc_A_pre); [mc_A_tmp] = tutorial_binary_uncompress(fname_mc_A_pre,1:nr_tmp,1:nc_tmp);
disp(sprintf('mc_A nsnps %.9d/%.9d',sum(mc_A_tmp(:)>0),length(mc_A_tmp)));
cols_used = sum(mc_A_tmp(:)>0);
 else;
disp(sprintf(' %% creating %s',fname_mc_A_pre));
mc_A_bim = lakcluster_PGC_getbim_ver1(fname_bim,p_threshold,p_threshold^2);
disp(sprintf(' %% p_threshold %0.2f, retaining %d mc values, but setting %d mc values to 0',p_threshold,sum(mc_A_bim),sum(~mc_A_bim)));
mc_A_tmp = mc_A_bim;
if ~found_trace_flag; 
tmpchar_out = fname_mc_A_pre;tutorial_binary_compress(bitj,mc_A_tmp(:)>0,tmpchar_out); 
tmpchar_out = fname_mc_A_pos;tutorial_binary_compress(bitj,mc_A_tmp(:)>0,tmpchar_out); 
end;%if ~found_trace_flag; 
disp(sprintf('mc_A nsnps %.9d/%.9d',sum(mc_A_tmp(:)>0),length(mc_A_tmp)));
cols_used = sum(mc_A_tmp(:)>0);
end;%if exist(fname_mc_A_pre,'file');

% T ;
tmpchar = sprintf('%s_T_full_n.b16',d_inpre); [bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); disp(sprintf(' %% reading %s = (%d,%d)',tmpchar,nrows,ncols));
tmp_T = tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols); 
Tn_crop = tmp_T(:,[1,(1+(mds_used))]); Tt_crop = transpose(Tn_crop); ncols_T_crop = size(Tn_crop,2);
tmpchar = sprintf('%s_T_crop_n.b16',d_oupre); disp(sprintf(' %% writing %s = (%d,%d)',tmpchar,size(Tn_crop)));
tutorial_binary_compress(bitj,Tn_crop>0,tmpchar);
tmpchar = sprintf('%s_T_crop_t.b16',d_oupre); disp(sprintf(' %% writing %s = (%d,%d)',tmpchar,size(Tn_crop,2),size(Tn_crop,1)));
tutorial_binary_compress(bitj,Tt_crop>0,tmpchar);
tmpchar = sprintf('%s_mc_T.b16',d_inpre); [bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); ncols_T = max(nrows,ncols);
%mc_T_tmp = zeros(ncols_T,1); mc_T_tmp(1)=1; mc_T_tmp(1 + (mds_used)) = ones(size(mds_used));
mc_T_crop = ones(ncols_T_crop,1);
if ~found_trace_flag; tmpchar = sprintf('%s_mc_T_crop.b16',d_oupre);tutorial_binary_compress(bitj,mc_T_crop(:)>0,tmpchar); end;
disp(sprintf('mc_T_crop:'));disp(num2str(transpose(mc_T_crop)));

disp(sprintf(' %% rows_used %d cols_used %d',rows_used,cols_used));

fname__in = sprintf('%s.in',d_oupre);
if ~found_trace_flag;
fp = fopen(fname__in,'w');
fprintf(fp,'GLOBAL_verbose= %d;\n',verbose_flag);
fprintf(fp,'GLOBAL_thread_count= 15;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','dexcluster_driver');
fprintf(fp,'GLOBAL_NBINS= %d;\n',nbins);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',B_MLT);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
fprintf(fp,'GLOBAL_A_n_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_A_full_n.b16',d_inpre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_A_t_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_A_full_t.b16',d_inpre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
if (rev_flag==1); fprintf(fp,'GLOBAL_A_n_rows_= '); for nb=0:nbins-1; fprintf(fp,'%d',nrows_Z_(1+nb)); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_A_n_rows_= '); for nb=0:nbins-1; fprintf(fp,'%d',nrows_A_(1+nb)); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',ncols_A);
if (rev_flag==1); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb=0:nbins-1; fprintf(fp,'%s_mr_Z_full.b16',d_oupre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_A_n_rind_= '); for nb=0:nbins-1; fprintf(fp,'%s_mr_A_full.b16',d_oupre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A.b16;\n',d_oupre);
fprintf(fp,'GLOBAL_Z_n_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_A_full_n.b16',d_inpre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_Z_t_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_A_full_t.b16',d_inpre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
if (rev_flag==1); fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb=0:nbins-1; fprintf(fp,'%d',nrows_A_(1+nb)); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_Z_n_rows_= '); for nb=0:nbins-1; fprintf(fp,'%d',nrows_Z_(1+nb)); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==1); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb=0:nbins-1; fprintf(fp,'%s_mr_A_full.b16',d_oupre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
if (rev_flag==0); fprintf(fp,'GLOBAL_Z_n_rind_= '); for nb=0:nbins-1; fprintf(fp,'%s_mr_Z_full.b16',d_oupre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n'); end;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',ncols_Y);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',ncols_T_crop); 
fprintf(fp,'GLOBAL_T_n_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_T_crop_n.b16',d_oupre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_t_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_T_crop_t.b16',d_oupre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T_crop.b16;\n',d_oupre);
fprintf(fp,'GLOBAL_S_n_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_T_crop_n.b16',d_oupre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_S_t_name_= '); for nb=0:nbins-1; fprintf(fp,'%s_T_crop_t.b16',d_oupre); if nb<nbins-1; fprintf(fp,','); else fprintf(fp,';'); end; end; fprintf(fp,'\n');
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by dexcluster_PGC_ADZSZDA_ver7.m on %s;\n',date);
fclose(fp);
end;%if ~found_trace_flag;

call_flag=1*(pbs_walltime<=0);%call_flag = input(' call? 1=yes (default), 0=no:'); if isempty(call_flag); call_flag=1; end;
if call_flag;
if ~found_trace_flag;
disp(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
system(sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in));
end;%if ~found_trace_flag;
end;%if call_flag;

pbs_flag=1*(pbs_walltime>0);
if pbs_flag;
fname_pbs = sprintf('%s.pbs',d_oupre);
if ~found_trace_flag;
fp = fopen(fname_pbs,'w');
fprintf(fp,'#PBS -S /bin/bash\n');
fprintf(fp,'#PBS -lnodes=1:ppn=15\n');
%fprintf(fp,'#PBS -lnodes=1:cpu2+cpu3\n');
pbs_walltime_h = floor(pbs_walltime); pbs_walltime_m = min(59,ceil(60*(pbs_walltime - pbs_walltime_h))); 
sprintf(' %% pbs_walltime=%d:%.2d:59',pbs_walltime_h,pbs_walltime_m);
fprintf(fp,'#PBS -lwalltime=%d:%.2d:59\n',pbs_walltime_h,pbs_walltime_m);
fprintf(fp,'#PBS -e %s_pbs_error.log\n',d_oupre);
fprintf(fp,'#PBS -o %s_pbs_output.log\n',d_oupre);
if ~found_trace_flag;
fprintf(fp,'%s/lakcluster_ver18 < %s\n',dir_code,fname__in);
end;%if ~found_trace_flag;
fclose(fp);
type(fname_pbs);
system(sprintf('qsub %s;\n',fname_pbs));
end;%if ~found_trace_flag;
end;%if pbs_flag;

end;%if ~found_trace_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output_string,gamma_d] = dexcluster_PGC_xfix_gen_ver7(p_threshold,rev_flag,n_mds,mds_used,gamma,B_MLT,Ireq,nshuffle);
p_str = sprintf('p%.2d',floor(100*p_threshold));
if (rev_flag==1); rev_str = 'X'; else rev_str = 'D'; end;
mds_tmp = zeros(1,n_mds); mds_tmp(mds_used)=1; mds_code = dot(2.^[0:n_mds-1],mds_tmp);
gamma_d = floor(gamma*1000);
output_string = sprintf('dex_ADZSZDA_%s_%s_m%d_g%.3d_B%.2d_n%d_s%.4d',p_str,rev_str,mds_code,gamma_d,B_MLT,Ireq,nshuffle);
