function lisa_xdropplot_initialize_ver2(str_input,mc_string,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
%{

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % check for lak / dex ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  str_input = 'checktrace'; B_MLT = 34; n_mds = 20; n_scramble = 0; n_shuffle = 64;
  for flag_dex_vs_lak = {'lak'}; %for flag_dex_vs_lak = {'dex','lak'};
  %if (strcmp(flag_dex_vs_lak{1},'dex')); gamma = 0.004; end; if (strcmp(flag_dex_vs_lak{1},'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
  if (strcmp(flag_dex_vs_lak{1},'dex')); gamma = 0.004; end; if (strcmp(flag_dex_vs_lak{1},'lak')); gamma = 0.001; mc_string = ''; end;
  for flag_reverse = [0:1]; %for flag_reverse = [0:1]; 
  for cl_num = 1:4; %for cl_num = 1:4; 
  cl_num_arm2 = []; 
  for n_maf = [4,5]; %for n_maf = 1:5; 
  for n_cov = 2; %for n_cov = 0:2; 
  lisa_xdropplot_initialize_ver2(str_input,mc_string,cl_num,cl_num_arm2,flag_dex_vs_lak{1},gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
  end;%for n_cov = 0:5; 
  end;%for n_maf = 1:5; 
  end;%for cl_num = 1:4; cl_num_arm2 = []; 
  end;%for flag_reverse = [0:1]; 
  end;%for flag_dex_vs_lak = {'dex','lak'}; 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % generate figures for dex ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  str_input = 'figplot'; B_MLT = 34; n_mds = 20; n_scramble = 1; n_shuffle = 64;
  for flag_dex_vs_lak = {'dex'}; 
  if (strcmp(flag_dex_vs_lak{1},'dex')); gamma = 0.004; end; if (strcmp(flag_dex_vs_lak{1},'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
  for flag_reverse = 0; 
  for cl_num = 4; 
  cl_num_arm2 = []; 
  for n_maf = 1:5; 
  for n_cov = 2; 
  lisa_xdropplot_initialize_ver2(str_input,mc_string,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
  end;%for n_cov = 0:2; 
  end;%for n_maf = 1:5; 
  end;%for cl_num = 1:4; 
  end;%for flag_reverse = 0:1;
  end;%for flag_dex_vs_lak = {'dex','lak'}; 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % generate pca projection for dex ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  str_input = 'pcaplot'; gamma = 0.004; B_MLT = 34; n_mds = 20; n_scramble = 1; n_shuffle = 64;
  flag_dex_vs_lak = 'dex'; 
  if (strcmp(flag_dex_vs_lak,'dex')); gamma = 0.004; end; if (strcmp(flag_dex_vs_lak,'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
  flag_reverse = 0; 
  cl_num = 4; cl_num_arm2 = 1; n_maf = 5; n_cov = 2;
  lisa_xdropplot_initialize_ver2(str_input,mc_string,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
  % generate pca projection for all cl_num_arm2 (dex) ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
  str_input = 'pcaplot'; B_MLT = 34; n_mds = 20; n_shuffle = 64;
  flag_dex_vs_lak = 'dex'; 
  if (strcmp(flag_dex_vs_lak,'dex')); gamma = 0.004; mc_string = ''; end; if (strcmp(flag_dex_vs_lak,'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
  for flag_reverse = 0;%for flag_reverse = 0:1;
  for n_scramble = 1;%for n_scramble = 0:1; 
  cl_num = 4; 
  for cl_num_arm2=1:4;
  if (cl_num_arm2==cl_num);%if (cl_num_arm2~=cl_num);
  for n_maf = [5];
  for n_cov = 2; 
  lisa_xdropplot_initialize_ver2(str_input,mc_string,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
  end;%for n_cov = 2; 
  end;%for n_maf = 1:5;
  end;%if (cl_num_arm2~=cl_num);
  end;%for cl_num_arm2=1:4;
  end;%for n_scramble = 0:1; 
  end;%for flag_reverse = 0:1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % generate figures for lak ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  str_input = 'figplot'; B_MLT = 34; n_mds = 20; n_scramble = 0; n_shuffle = 64;
  flag_dex_vs_lak = 'lak'; 
  if (strcmp(flag_dex_vs_lak,'dex')); gamma = 0.004; end; if (strcmp(flag_dex_vs_lak,'lak')); gamma = 0.001; mc_string = 'm20_p85'; mc_string = ''; end;
  for flag_reverse = 1; %for flag_reverse = 0:1; 
  cl_num = 3; cl_num_arm2 = []; n_maf = 4; 
  for n_cov = 2;%for n_cov = 0:2; 
  lisa_xdropplot_initialize_ver2(str_input,mc_string,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
  end;%for n_cov = 0:2; 
  end;%for flag_reverse = 0:1; 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  % generate pca projection for lak ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
  str_input = 'pcaplot'; B_MLT = 34; n_mds = 20; n_scramble = 1; n_shuffle = 0;
  flag_dex_vs_lak = 'lak'; 
  if (strcmp(flag_dex_vs_lak,'dex')); gamma = 0.004; end; if (strcmp(flag_dex_vs_lak,'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
  flag_reverse = 0; cl_num = 4; cl_num_arm2 = 1; n_maf = 5; n_cov = 2; 
  lisa_xdropplot_initialize_ver2(str_input,mc_string,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
  % generate pca projection for all cl_num_arm2 (lak) ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
  str_input = 'pcaplot'; B_MLT = 34; n_mds = 20; n_scramble = 0; n_shuffle = 64;
  flag_dex_vs_lak = 'lak'; 
  if (strcmp(flag_dex_vs_lak,'dex')); gamma = 0.004; end; if (strcmp(flag_dex_vs_lak,'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
  for flag_reverse = 0:1;%for flag_reverse = 0:1;
  cl_num = 2; 
  for cl_num_arm2=1:4;
  if (cl_num_arm2~=cl_num);
  for n_maf = 4;%for n_maf = 4:5; 
  for n_cov = 2;%for n_cov = 0:2; 
  lisa_xdropplot_initialize_ver2(str_input,mc_string,cl_num,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
  end;%for n_cov = 0:2; 
  end;%for n_maf = 1:5;
  end;%if (cl_num_arm2~=cl_num);
  end;%for cl_num_arm2=1:4;
  end;%for flag_reverse = 0:1;

  %%%%%%%%%%%%%%%%;
  % clean ;
  %%%%%%%%%%%%%%%%;
  str_input = 'clean'; gamma = 0.004; B_MLT = 34; n_mds = 20; n_scramble = 0; n_shuffle = 64;
  for flag_dex_vs_lak = {'dex','lak'}; 
  if (strcmp(flag_dex_vs_lak{1},'dex')); gamma = 0.004; end; if (strcmp(flag_dex_vs_lak{1},'lak')); gamma = 0.001; mc_string = 'm20_p85'; end;
  for flag_reverse = 0:1;
  for cl_num = 1:4;
  for n_maf = 1:5;
  for n_cov = 1:5; 
  lisa_xdropplot_initialize_ver2(str_input,mc_string,cl_num,[],flag_dex_vs_lak{1},gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
  end;%for n_cov = 1:5; 
  end;%for n_maf = 1:5;
  end;%for cl_num = 1:4;
  end;%for flag_reverse = 0:1;
  end;%for flag_dex_vs_lak = {'dex','lak'}; 
    
  %}

ni=1;
if (nargin<ni); str_input = 'figplot'; end; ni=ni+1;
if (nargin<ni); mc_string = ''; end; ni=ni+1;
if (nargin<ni); cl_num = 1; end; ni=ni+1;
if (nargin<ni); cl_num_arm2 = 2; end; ni=ni+1;
if (nargin<ni); flag_dex_vs_lak = 'dex'; end; ni=ni+1;
if (nargin<ni); gamma = 0.004; end; ni=ni+1;
if (nargin<ni); B_MLT = 34; end; ni=ni+1;
if (nargin<ni); n_mds = 20; end; ni=ni+1;
if (nargin<ni); flag_reverse = 0; end; ni=ni+1;
if (nargin<ni); n_maf = 5; end; ni=ni+1;
if (nargin<ni); n_cov = 2; end; ni=ni+1;
if (nargin<ni); n_scramble = 0; end; ni=ni+1;
if (nargin<ni); n_shuffle = 64; end; ni=ni+1;

verbose = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% choose directory names and set prefix used during analysis. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
lisa_setprefix_ver2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% check for completeness. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if strfind(str_input,'checktrace');
lisa_xdropplot_checktrace_ver2;
end;%if strfind(str_input,'check');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% generate figplot. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if strfind(str_input,'figplot');
lisa_xdropplot_loadtrace_ver2;
lisa_xdropplot_figplot_ver0;
end;%if strfind(str_input,'figplot');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% perform pca. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if strfind(str_input,'pcaplot');
lisa_xdropplot_loadtrace_ver2;
lisa_xdropplot_pcaplot_ver3;
end;%if strfind(str_input,'pcaplot');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% clean. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if strfind(str_input,'clean');
lisa_xdropplot_clean_ver1;
end;%if strfind(str_input,'clean');