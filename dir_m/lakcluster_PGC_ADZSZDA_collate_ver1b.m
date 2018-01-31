function lakcluster_PGC_ADZSZDA_collate_ver1b(dir_trunk,prefix,infix,rev_flag,nmds,mds_used,gamma_,B_MLT_,Ireq_,shuffle_num_)
% runs through the directories indexed by PGC_xfix_gen_ver1b and collects output_trace ;
% try with: 
%{
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/?*.[m] ./ ;
  !scp -rp rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/lakcluster_PGC_ADZSZDA_collate_ver1b.m ./ ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % collating rerun of PGC_cl3 ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  dir_trunk = sprintf('%s',pwd); cl_num=3; 
  if cl_num==1; prefix = 'PGC_cl1'; nbins = 11; Ireq = floor(nbins/2); end;
  if cl_num==2; prefix = 'PGC_cl2'; nbins = 7; Ireq = floor(nbins/2); end;
  if cl_num==3; prefix = 'PGC_cl3'; nbins = 10; Ireq = floor(nbins/2); end;
  %infix = 'lak_ADZSZDA';
  infix = 'dex_ADZSZDA';
  %infix = '';
  for rev_flag = [1]; for mds_tab = [1];
  %try;
  nmds = 20; if mds_tab==1; mds_used = [1:2]; else mds_used = []; end; 
  gamma=0.05;B_MLT=34;Ireq=0;shuffle_num=0;shuffle_num_=[0:64];
  lakcluster_PGC_ADZSZDA_collate_ver1b(dir_trunk,prefix,infix,rev_flag,nmds,mds_used,[gamma],[B_MLT],[Ireq],shuffle_num_);
  %catch; disp(' error'); end;%try;
  end;end;%for mds_tab = 0:1;%for rev_flag = 0:1; 
  !scp -p PGC_cl3_dex_ADZSZDA_X_m3_g050_B34_n0_sXXXX.mat rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_code_081717/dir_m/dir_PGC_cl3/ ;

 %}

na=1;
if (nargin<na); dir_trunk = pwd; end; na=na+1;
if (nargin<na); prefix = 'PGC_cl3_maf01'; end; na=na+1;
if (nargin<na); prefix = 'lak_ADZSZDA'; end; na=na+1;
if (nargin<na); rev_flag = 0; end; na=na+1;
if (nargin<na); nmds = 20; end; na=na+1;
if (nargin<na); mds_used = [1:2]; end; na=na+1;
if (nargin<na); gamma = 0.1; end; na=na+1;
if (nargin<na); B_MLT = 8; end; na=na+1;
if (nargin<na); Ireq = 0; end; na=na+1;
if (nargin<na); shuffle_num = 0; end; na=na+1;
if (nargin<na); walltime = 13; end; na=na+1;
if (nargin<na); row_factor = 2.0; end; na=na+1;
if (nargin<na); col_factor = 2.0; end; na=na+1;

dir_base = sprintf('%s/dir_%s/dir_%s',dir_trunk,prefix,prefix);
clear xa_trn__ xb_trn__ tr_trn__ t_trn_elct t_trn_elrt ;
n_found_trn=0; n_tot = length(gamma_)*length(B_MLT_)*length(Ireq_)*length(shuffle_num_);
for n_gamma=1:length(gamma_);
for n_B_MLT=1:length(B_MLT_);
for n_Ireq=1:length(Ireq_);
for n_shuffle_num=1:length(shuffle_num_);
% get xdrop_a and xdrop_b ;
if (n_shuffle_num==1);
gamma = gamma_(n_gamma); B_MLT = B_MLT_(n_B_MLT); Ireq = Ireq_(n_Ireq); shuffle_num = shuffle_num_(n_shuffle_num);
xfix = PGC_xfix_gen_ver1b(infix,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,shuffle_num);
tmpchar_trn = sprintf('%s_%s/out_xdrop_a.txt',dir_base,xfix);
if exist(tmpchar_trn,'file'); 
tmp = textread(tmpchar_trn);
disp(sprintf('found %s of size %d-x-%d',tmpchar_trn,size(tmp)));
xa_trn__{n_gamma}{n_B_MLT}{n_Ireq}{n_shuffle_num} = tmp;
end;% if exists;
if ~exist(tmpchar_trn,'file'); disp(sprintf('could not locate %s',tmpchar_trn)); end;
tmpchar_trn = sprintf('%s_%s/out_xdrop_b.txt',dir_base,xfix);
if exist(tmpchar_trn,'file'); 
tmp = textread(tmpchar_trn);
disp(sprintf('found %s of size %d-x-%d',tmpchar_trn,size(tmp)));
xb_trn__{n_gamma}{n_B_MLT}{n_Ireq}{n_shuffle_num} = tmp;
end;% if exists;
if ~exist(tmpchar_trn,'file'); disp(sprintf('could not locate %s',tmpchar_trn)); end;
end;%if (n_shuffle_num==1);
% get trace ;
gamma = gamma_(n_gamma); B_MLT = B_MLT_(n_B_MLT); Ireq = Ireq_(n_Ireq); shuffle_num = shuffle_num_(n_shuffle_num);
xfix = PGC_xfix_gen_ver1b(infix,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,shuffle_num);
tmpchar_trn = sprintf('%s_%s/out_trace.txt',dir_base,xfix);
if exist(tmpchar_trn,'file'); 
tmp = textread(tmpchar_trn);
disp(sprintf('found %s of size %d-x-%d',tmpchar_trn,size(tmp))); n_found_trn = n_found_trn+1;
tr_trn__{n_gamma}{n_B_MLT}{n_Ireq}{n_shuffle_num} = tmp;
end;% if exists;
if ~exist(tmpchar_trn,'file'); disp(sprintf('could not locate %s',tmpchar_trn)); end;
% get timing ;
gamma = gamma_(n_gamma); B_MLT = B_MLT_(n_B_MLT); Ireq = Ireq_(n_Ireq); shuffle_num = shuffle_num_(n_shuffle_num);
xfix = PGC_xfix_gen_ver1b(infix,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,shuffle_num);
tmpchar_trn = sprintf('%s_%s/timing.m',dir_base,xfix);
if exist(tmpchar_trn,'file');
tmp_fp = fopen(tmpchar_trn); [tmp_a] = textscan(tmp_fp,'elct = %f; elrt = %f;'); fclose(tmp_fp);
t_trn_elct(n_gamma,n_B_MLT,n_Ireq,n_shuffle_num) = tmp_a{1};
t_trn_elrt(n_gamma,n_B_MLT,n_Ireq,n_shuffle_num) = tmp_a{2};
end;% if exists;
if ~exist(tmpchar_trn,'file'); disp(sprintf('could not locate %s',tmpchar_trn)); end;
end;end;end;end;% input arrays ;
disp(sprintf('found trn %d/%d ',n_found_trn,n_tot));
savename = sprintf('%s_%s.mat',prefix,PGC_savename_gen_ver1b(infix,rev_flag,nmds,mds_used,gamma_,B_MLT_,Ireq_,shuffle_num_));
disp(sprintf('saving to %s',savename));
save(savename,'xa_trn__','xb_trn__','tr_trn__','t_trn_elct','t_trn_elrt','dir_trunk','prefix','gamma_','B_MLT_','Ireq_','shuffle_num_','dir_base');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output_string,gamma_d] = PGC_xfix_gen_ver1b(infix,rev_flag,nmds,mds_used,gamma,B_MLT,Ireq,shuffle_num);
if (rev_flag==1); rev_str = 'X'; else rev_str = 'D'; end;
mds_tmp = zeros(1,nmds); mds_tmp(mds_used)=1; mds_code = dot(2.^[0:nmds-1],mds_tmp);
%gamma_d = floor(gamma*1000);
gamma_d = floor(gamma*100);
if (length(infix)>0);
output_string = sprintf('%s_%s_m%d_g%.2d_B%.2d_n%d_s%.4d',infix,rev_str,mds_code,gamma_d,B_MLT,Ireq,shuffle_num);
end;%if (length(infix)>0);
if (length(infix)==0);
output_string = sprintf('%s_m%d_g%.2d_B%.2d_n%d_s%.4d',rev_str,mds_code,gamma_d,B_MLT,Ireq,shuffle_num);
end;%if (length(infix)==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output_string,gamma_d] = PGC_savename_gen_ver1b(infix,rev_flag,nmds,mds_used,gamma_,B_MLT_,Ireq_,shuffle_num_);
if (rev_flag==1); rev_str = 'X'; else rev_str = 'D'; end;
mds_tmp = zeros(1,nmds); mds_tmp(mds_used)=1; mds_code = dot(2.^[0:nmds-1],mds_tmp);
gamma_d = floor(gamma_*1000);
if (length(infix)>0);
output_string = sprintf('%s_%s_m%d',infix,rev_str,mds_code);
end;%if (length(infix)>0);
if (length(infix)==0);
output_string = sprintf('%s_m%d',rev_str,mds_code);
end;%if (length(infix)==0);
if length(gamma_)>1; output_string = sprintf('%s_gXXX',output_string); else output_string = sprintf('%s_g%.3d',output_string,floor(gamma_(1)*1000)); end;
if length(B_MLT_)>1; output_string = sprintf('%s_BXX',output_string); else output_string = sprintf('%s_B%.2d',output_string,B_MLT_(1)); end;
if length(Ireq_)>1; output_string = sprintf('%s_nX',output_string); else output_string = sprintf('%s_n%.1d',output_string,Ireq_(1)); end;
if length(shuffle_num_)>1; output_string = sprintf('%s_sXXXX',output_string); else output_string = sprintf('%s_s%.4d',output_string,shuffle_num_(1)); end;
