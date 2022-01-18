function lisa_xdropplot_initialize_ver0(str_input,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_shuffle) ;
%{

  str_input = 'pcaplot'; gamma = 0.004; B_MLT = 34; n_mds = 20; n_shuffle = 64;
  flag_dex_vs_lak = 'dex'; flag_reverse = 1; 
  cl_num = 4; 
  for cl_num_arm2=1:4;
  if (cl_num_arm2~=cl_num);
  for n_maf = 1:5;
  n_cov = 1; 
  lisa_xdropplot_initialize_ver0(str_input,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_shuffle) ;
  end;%for n_maf = 1:5;
  end;%if (cl_num_arm2~=cl_num);
  end;%for cl_num_arm2=1:4;

  str_input = 'pcaplot'; gamma = 0.004; B_MLT = 34; n_mds = 20; n_shuffle = 64;
  flag_dex_vs_lak = 'dex'; flag_reverse = 0; 
  cl_num = 1; cl_num_arm2 = 2; n_maf = 2; n_cov = 1; 
  lisa_xdropplot_initialize_ver0(str_input,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_shuffle) ;
    
  %}
ni=1;
if (nargin<ni); str_input = 'figplot'; end; ni=ni+1;
if (nargin<ni); cl_num = 1; end; ni=ni+1;
if (nargin<ni); cl_num_arm2 = 2; end; ni=ni+1;
if (nargin<ni); flag_dex_vs_lak = 'dex'; end; ni=ni+1;
if (nargin<ni); gamma = 0.004; end; ni=ni+1;
if (nargin<ni); B_MLT = 34; end; ni=ni+1;
if (nargin<ni); n_mds = 20; end; ni=ni+1;
if (nargin<ni); flag_reverse = 0; end; ni=ni+1;
if (nargin<ni); n_maf = 5; end; ni=ni+1;
if (nargin<ni); n_cov = 1; end; ni=ni+1;
if (nargin<ni); n_shuffle = 64; end; ni=ni+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% choose directory names. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c'); 
dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c_dev'); 
dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PGC_20180304');
path(path,sprintf('%s/dir_m',dir_code)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% select parameters used during analysis. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%cl_num = 4; 
[study_trunk_,study_name_,n_study] = lisa_define_study_ver0(cl_num);
Ireq_half = floor(n_study/2); Icat_full = n_study;
string_prefix = sprintf('PGC_cl%d_maf01',cl_num);
%flag_dex_vs_lak = 'dex'; 
%gamma=0.004;B_MLT=34;n_mds=20;
%flag_reverse = 0;
%n_maf = 5;
%n_cov = 1;
%%%%%%%%%%%%%%%% ;
if n_maf==1; maf_lo_threshold = 0.01; maf_hi_threshold = 0.50; end; % all (since we thresholded maf at 0.01 to start with) ;
if n_maf==2; maf_lo_threshold = 0.01; maf_hi_threshold = 0.05; end; % extreme snps (differential expression only, no loop-counting later on) ;
if n_maf==3; maf_lo_threshold = 0.01; maf_hi_threshold = 0.10; end; % extreme snps (differential expression only, no loop-counting later on) ;
if n_maf==4; maf_lo_threshold = 0.10; maf_hi_threshold = 0.50; end; % more balanced snps (differential expression to start, as well as loop-counting later on) ;
if n_maf==5; maf_lo_threshold = 0.25; maf_hi_threshold = 0.50; end; % more balanced snps (differential expression to start, as well as loop-counting later on) ;
%%%%%%%%%%%%%%%% ;
if n_cov==1; Ireq=0; Icat=1; mds_used_=[1:2]; mds_repl=2; end; % only correct for mds-components [1:2], but replicated twice, so that each 'sector' is 45 degrees ;
if n_cov==2; Ireq=Ireq_half; Icat=Icat_full; mds_used_=[]; mds_repl=0; end; % correct for study, requiring at least Ireq_half, but not for mds-components. ;
if n_cov==3; Ireq=Ireq_half; Icat=Icat_full; mds_used_=[1:2]; mds_repl=1; end; % correct for study, requiring at least Ireq_half, as well as correct for mds-components [1:2], but replicated only once, so that each 'sector' is 90 degrees ;
if n_cov==4; Ireq=Ireq_half; Icat=Icat_full; mds_used_=[1:2]; mds_repl=2; end; % correct for study, requiring at least Ireq_half, as well as correct for mds-components [1:2], replicated twice, so that each 'sector' is 45 degrees ;
%%%%%%%%%%%%%%%% ;
%n_shuffle = 64;
dir__in = sprintf('%s/dir_%s',dir_trunk,string_prefix);
dir_out = sprintf('%s/dir_%s_analyze',dir_trunk,string_prefix);
disp(sprintf(' %% dir__in: %s',dir__in));
disp(sprintf(' %% dir_out: %s',dir_out));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% generate figplot. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if strfind(str_input,'figplot');
lisa_xdropplot_loadtrace_ver0;
lisa_xdropplot_figplot_ver0;
end;%if strfind(str_input,'figplot');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% perform pca. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if strfind(str_input,'pcaplot');
lisa_xdropplot_loadtrace_ver0;
lisa_xdropplot_pcaplot_ver0;
end;%if strfind(str_input,'pcaplot');
