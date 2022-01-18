function [lisa] = lisa_struct_fig_wrap_ver0(specification) ;
%{
 %try: ;

  dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PGC_20190328');
  dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c_dev');
  flag_dex_vs_lak = 'dex'; %<-- differentially expressed clustering. ;
  cl_num = 4; %<-- train on platform 4. ;
  flag_reverse = 0; %<-- forward bicluster (i.e., case-specific). ;
  if (strcmp(flag_dex_vs_lak,'dex')); gamma = [0.004]; end; %<-- gamma is the fraction eliminated per iteration. 000 implies a single patient eliminated per iteration. ;
  if (strcmp(flag_dex_vs_lak,'lak')); gamma = [0.001]; end; %<-- gamma is the fraction eliminated per iteration. 000 implies a single patient eliminated per iteration. ;
  B_MLT = 34; n_mds = 20; mr_string = '';mc_string = ''; %<-- accurate to 2^(-34), 20 total mds components (but only 2 used). ; No special mc_string. ;
  n_maf = 5; n_cov = 2; %<-- minor-allele-frequency cutoff 25-50, 2 covariates (mds-components) used, repeated twice. ;
  n_scramble = 0; n_shuffle = 128; %<-- no previous bicluster extracted/scrambled first, no random shuffling. ;
  flag_rerun=0; %<-- regenerate mat-file.; 
  %mr_string = 'BX1'; mc_string = 'BX1';
  for flag_reverse = 0;%for flag_reverse = 0:1;
  for n_scramble = 0;%for n_scramble = 0:1;
  for cl_num = 4;%for cl_num = 1:4;
  specification = struct();
  specification.dir_trunk = dir_trunk;
  specification.dir_code = dir_code;
  specification.mr_string = mr_string;
  specification.mc_string = mc_string;
  specification.cl_num = cl_num;
  specification.flag_dex_vs_lak = flag_dex_vs_lak;
  specification.gamma = gamma;
  specification.B_MLT = B_MLT;
  specification.n_mds = n_mds;
  specification.flag_reverse = flag_reverse;
  specification.n_maf = n_maf;
  specification.n_cov = n_cov;
  specification.n_scramble = n_scramble;
  specification.n_shuffle = n_shuffle;
  specification.flag_rerun = flag_rerun;
  lisa_struct_fig_wrap_ver0(specification) ; close all;
  end;%for cl_num = 1:4;
  end;%for n_scramble = 0:1;
  end;%for flag_reverse = 0:1;

 %}

setup;

dir_trunk = specification.dir_trunk;
dir_code = specification.dir_code;
mr_string = specification.mr_string;
mc_string = specification.mc_string;
cl_num = specification.cl_num;
flag_dex_vs_lak = specification.flag_dex_vs_lak;
gamma = specification.gamma;
B_MLT = specification.B_MLT;
n_mds = specification.n_mds;
flag_reverse = specification.flag_reverse;
n_maf = specification.n_maf;
n_cov = specification.n_cov;
n_scramble = specification.n_scramble;
n_shuffle = specification.n_shuffle;
flag_rerun = specification.flag_rerun;

%%%%%%%%;
lisa = lisa_struct_make_ver0(mr_string,mc_string,cl_num,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
lisa = lisa_struct_prefix_ver0(lisa,dir_code,dir_trunk); 
lisa.nshuffle = 0;  lisa = lisa_struct_names_ver0(lisa); 
lisa = lisa_struct_xdrop_ver0(lisa); lisa = lisa_struct_mdsfam_ver0(lisa); 
lisa = lisa_struct_bim_ver0(lisa); %<-- this is large and takes a while to load. ;
lisa = lisa_struct_mx_ver0(lisa); lisa = lisa_struct_studyindex_ver0(lisa); 
lisa = lisa_struct_trace_ver0(lisa);
lisa = lisa_struct_BD_ver0(lisa);
%%%%%%%%;
lisa = lisa_struct_fig_ver0(lisa);

