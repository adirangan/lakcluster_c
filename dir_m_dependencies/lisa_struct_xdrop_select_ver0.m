function lisa = lisa_struct_xdrop_select_ver0(specification);
%{
  dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PGC_20190328');
  dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c_dev');
  flag_dex_vs_lak = 'dex'; %<-- differentially expressed clustering. ;
  cl_num = 4; %<-- train on platform 4. ;
  flag_reverse = 0; %<-- forward bicluster (i.e., case-specific). ;
  gamma = [0.004]; %<-- gamma is the fraction eliminated per iteration. 000 implies a single patient eliminated per iteration. ;
  B_MLT = 34; n_mds = 20; mr_string = '';mc_string = ''; %<-- accurate to 2^(-34), 20 total mds components (but only 2 used). ; No special mc_string. ;
  n_maf = 5; n_cov = 2; %<-- minor-allele-frequency cutoff 25-50, 2 covariates (mds-components) used, repeated twice. ;
  n_scramble = 0; n_shuffle = 0; %<-- no previous bicluster extracted/scrambled first, no random shuffling. ;
  flag_rerun=0; %<-- regenerate file.; 
  niteration = 175; %<-- iteration to select. ;
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
  specification.niteration = niteration;
  %{
  cl_num = 1; flag_reverse = 0; niteration = 475; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 2; flag_reverse = 0; niteration = 325; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 3; flag_reverse = 0; niteration = 450; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 4; flag_reverse = 0; niteration = 325; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 1; flag_reverse = 1; niteration = 525; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 2; flag_reverse = 1; niteration = 325; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 3; flag_reverse = 1; niteration = 525; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  %cl_num = 4; flag_reverse = 1; niteration =  25; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  %}
  %{
  cl_num = 1; flag_reverse = 0; niteration = 300; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 2; flag_reverse = 0; niteration = 200; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 3; flag_reverse = 0; niteration = 325; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 4; flag_reverse = 0; niteration = 250; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 1; flag_reverse = 1; niteration = 300; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 2; flag_reverse = 1; niteration = 250; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  cl_num = 3; flag_reverse = 1; niteration = 425; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  %cl_num = 4; flag_reverse = 1; niteration =  25; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  %}
  %{
    cl_num = 4; flag_reverse = 0; n_maf = 1; niteration =  50; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.n_maf = n_maf; specification.niteration = niteration; 
    lisa_struct_xdrop_select_ver0(specification) ;
    cl_num = 4; flag_reverse = 0; n_maf = 4; niteration =  200; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.n_maf = n_maf; specification.niteration = niteration; 
    lisa_struct_xdrop_select_ver0(specification) ;
    cl_num = 4; flag_reverse = 0; n_maf = 6; niteration =  150; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.n_maf = n_maf; specification.niteration = niteration; 
    lisa_struct_xdrop_select_ver0(specification) ;
    %}
  %{
    cl_num = 4; flag_reverse = 0; n_maf = 1; niteration =  125; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.n_maf = n_maf; specification.niteration = niteration; 
    mr_string = 'BDX'; mc_string = 'BDX'; specification.mr_string = mr_string; specification.mc_string = mc_string;
    lisa_struct_xdrop_select_ver0(specification) ;
    cl_num = 4; flag_reverse = 0; n_maf = 6; niteration =  138; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.n_maf = n_maf; specification.niteration = niteration; 
    mr_string = 'BDX'; mc_string = 'BDX'; specification.mr_string = mr_string; specification.mc_string = mc_string;
    lisa_struct_xdrop_select_ver0(specification) ;
    cl_num = 4; flag_reverse = 0; n_maf = 4; niteration =  163; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.n_maf = n_maf; specification.niteration = niteration; 
    mr_string = 'BDX'; mc_string = 'BDX'; specification.mr_string = mr_string; specification.mc_string = mc_string;
    lisa_struct_xdrop_select_ver0(specification) ;
    cl_num = 4; flag_reverse = 0; n_maf = 5; niteration =  150; specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.n_maf = n_maf; specification.niteration = niteration; 
    mr_string = 'BDX'; mc_string = 'BDX'; specification.mr_string = mr_string; specification.mc_string = mc_string;
    lisa_struct_xdrop_select_ver0(specification) ;    
    %}
  %{
    cl_num = 4; flag_reverse = 0; n_maf = 5; n_scramble = 1; niteration =  263; 
    specification.cl_num = cl_num; specification.flag_reverse = flag_reverse; specification.n_maf = n_maf; specification.n_scramble = n_scramble; specification.niteration = niteration; 
    lisa_struct_xdrop_select_ver0(specification) ;
    
    %}
  for niteration=round(150:25/2:375);
  specification.niteration = niteration; lisa_struct_xdrop_select_ver0(specification) ;
  end;%for niteration=round(150:25/2:375);
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
niteration = specification.niteration;

%%%%%%%%;
lisa = lisa_struct_make_ver0(mr_string,mc_string,cl_num,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
lisa = lisa_struct_prefix_ver0(lisa,dir_code,dir_trunk); 
lisa.nshuffle = 0;  lisa = lisa_struct_names_ver0(lisa); 
lisa.dir_out_s0_select = sprintf('%s/dir_select',lisa.dir_out_s0);
if (~exist(lisa.dir_out_s0_select,'dir')); disp(sprintf(' %% creating %s',lisa.dir_out_s0_select)); mkdir(lisa.dir_out_s0_select); end;
fname_xdrop = sprintf('%s/out_xdrop_ni%d.txt',lisa.dir_out_s0_select,niteration);
fname_rdrop = sprintf('%s/out_rdrop_ni%d.txt',lisa.dir_out_s0_select,niteration);
fname_cdrop = sprintf('%s/out_cdrop_ni%d.txt',lisa.dir_out_s0_select,niteration);
flag_exist = (exist(fname_xdrop,'file') & exist(fname_rdrop,'file') & exist(fname_cdrop,'file'));
if (flag_exist & ~flag_rerun); disp(sprintf(' %% found %s',fname_xdrop)); end;
if (~flag_exist | (flag_exist & flag_rerun));
disp(sprintf(' %% could not find %s, creating',fname_xdrop));
%%%%%%%%;
lisa = lisa_struct_xdrop_ver0(lisa); lisa = lisa_struct_mdsfam_ver0(lisa); 
lisa = lisa_struct_bim_ver0(lisa); %<-- this is large and takes a while to load. ;
lisa = lisa_struct_mx_ver0(lisa); lisa = lisa_struct_studyindex_ver0(lisa); 
lisa = lisa_struct_trace_ver0(lisa);
%%%%%%%%;
r_rem = lisa.r_rem_(niteration); c_rem = lisa.c_rem_(niteration);
rdrop_a_rem_ = lisa.rdrop_a_(end-r_rem+1:end); cdrop_a_rem_ = lisa.cdrop_a_(end-c_rem+1:end);
%%%%%%%%;
scramble_out_xdrop = zeros(r_rem + c_rem,2);
scramble_out_xdrop(1:r_rem,1) = rdrop_a_rem_(:);
scramble_out_xdrop(r_rem + (1:c_rem),2) = cdrop_a_rem_(:);
%fname_xdrop = sprintf('%s/out_xdrop_ni%d.txt',lisa.dir_out_s0_select,niteration);
fp = fopen(fname_xdrop,'w'); fprintf(fp,'%d %d\n',transpose(scramble_out_xdrop-1)); fclose(fp);
disp(sprintf(' %% printed to:\n%s',fname_xdrop));
%%%%%%%%;
%fname_rdrop = sprintf('%s/out_rdrop_ni%d.txt',lisa.dir_out_s0_select,niteration);
fid=fopen(fname_rdrop,'w');
for nr=1:r_rem;
ma = rdrop_a_rem_(nr);
fprintf(fid,'%s %s %s %s %d %d %s %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f %0.16f\n',lisa.fam_{1}{ma},lisa.fam_{2}{ma},lisa.fam_{3}{ma},lisa.fam_{4}{ma},lisa.fam_{5}(ma),lisa.fam_{6}(ma),lisa.fam_{7}{ma},lisa.mds_sort_(ma,:));
end;%for nr=1:r_rem;
fclose(fid);
disp(sprintf('%s',fname_rdrop));
%%%%%%%%;
%fname_cdrop = sprintf('%s/out_cdrop_ni%d.txt',lisa.dir_out_s0_select,niteration);
fid=fopen(fname_cdrop,'w');
for nc=1:c_rem;
na = cdrop_a_rem_(nc);
fprintf(fid,'%s\t%c\t%c\t%s\t%f\t%f\t%f\n',lisa.bim_{1}{na},lisa.bim_{2}(na),lisa.bim_{3}(na),lisa.bim_{4}{na},lisa.bim_{5}(na),lisa.bim_{6}(na),lisa.bim_{7}(na));
end;%for nc=1:c_rem;
fclose(fid);
disp(sprintf('%s',fname_cdrop));
%%%%%%%%;
end;%if (~flag_exist | (flag_exist & flag_rerun));


