%%%%%%%%%%%%%%%%;
% Multiple runs. ;
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% declare trunk ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
flag_platform = 2;  flag_scramble=1;
if flag_platform==1; string_platform = '/home/arangan'; end;
if flag_platform==2; string_platform = '/data/rangan/dir_bcc'; end;
if flag_platform==3; string_platform = '/scratch/avr209'; end;
dir_trunk = sprintf('%s/dir_PGC_20190328',string_platform);
dir_code = sprintf('%s/dir_lakcluster_c_dev',string_platform);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
flag_rerun = 0;
B_MLT=34;n_mds=20; n_shuffle = 0;
slurm_nnodes = 1; slurm_tpn = 15; slurm_memdecl = 64; % Note, if 64+GB requested we should run on snappy. ;
row_factor=1.0;col_factor=1.0; 
flag_verbose=0;
for flag_dex_vs_lak = {'dex'}; %for flag_dex_vs_lak = {'dex','lak'};
%if (strcmp(flag_dex_vs_lak{1},'dex')); gamma=0.004; mr_string = 'BDX'; mc_string = 'BDX'; end;
if (strcmp(flag_dex_vs_lak{1},'dex')); gamma=0.004; mr_string = ''; mc_string = ''; end;
if (strcmp(flag_dex_vs_lak{1},'lak')); gamma=0.001; mr_string = ''; mc_string = ''; end;
for flag_reverse = [0];%for flag_reverse=[0,1]; 
for cl_num = 4; %for cl_num=1:4;
if (strcmp(flag_dex_vs_lak{1},'dex') & cl_num> 1); pbs_walltime = 13.50*strcmp(string_platform,'/home/arangan'); slurm_walltime = 13.50*strcmp(string_platform,'/scratch/avr209'); end;
if (strcmp(flag_dex_vs_lak{1},'lak') & cl_num> 1); pbs_walltime = 64.50*strcmp(string_platform,'/home/arangan'); slurm_walltime = 64.50*strcmp(string_platform,'/scratch/avr209'); end;
if (strcmp(flag_dex_vs_lak{1},'dex') & cl_num==1); pbs_walltime = 13.50*strcmp(string_platform,'/home/arangan'); slurm_walltime = 13.50*strcmp(string_platform,'/scratch/avr209'); end;
if (strcmp(flag_dex_vs_lak{1},'lak') & cl_num==1); pbs_walltime = 180.50*strcmp(string_platform,'/home/arangan'); slurm_walltime = 180.50*strcmp(string_platform,'/scratch/avr209'); end;
[study_trunk_,study_name_,n_study] = lisa_define_study_ver0(cl_num); Icat_full = n_study; Ireq_half = floor(Icat_full/2);
string_prefix = sprintf('PGC_cl%d_maf01',cl_num);
for n_maf = 5;%for n_maf = [1,6,5,4];
[maf_lo_threshold,maf_hi_threshold] = lisa_set_maf_xx_threshold(n_maf);
for n_cov=2;%for n_cov=[0,1,2];
[Ireq,Icat,mds_used_,mds_repl] = lisa_set_mds_used(n_cov);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% original run (no shuffle) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_platform==1 | flag_platform==2 | flag_platform==3);
nshuffle=0; 
n_scramble=0;
xxxcluster_PGC_uADZSZDA_ver15(dir_code,dir_trunk,string_prefix,maf_lo_threshold,maf_hi_threshold,mr_string,mc_string,flag_dex_vs_lak{1},flag_reverse,Icat,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,{},[],nshuffle,flag_rerun,pbs_walltime,slurm_walltime,slurm_nnodes,slurm_tpn,slurm_memdecl,row_factor,col_factor,flag_verbose);
%%%%%%%%;
if flag_scramble;
n_scramble=1; scramble_rseed_ = [1]; flag_run = 0;
if (cl_num==1 & flag_reverse==1); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl1_maf01_analyze/dir_PGC_cl1_maf01_dex_p25_X_m2r2_g004/dir_select/out_xdrop_BX1.txt',dir_trunk)}; flag_run=1;
%%%%%%%%;
elseif (cl_num==4 & flag_reverse==0 & n_maf==1 & strcmp(mr_string,'') & strcmp(mc_string,'')); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p01_D_m2r2_g004/dir_select/out_xdrop_ni50.txt',dir_trunk)}; flag_run=1;
elseif (cl_num==4 & flag_reverse==0 & n_maf==6 & strcmp(mr_string,'') & strcmp(mc_string,'')); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p05_D_m2r2_g004/dir_select/out_xdrop_ni150.txt',dir_trunk)}; flag_run=1;
elseif (cl_num==4 & flag_reverse==0 & n_maf==4 & strcmp(mr_string,'') & strcmp(mc_string,'')); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p10_D_m2r2_g004/dir_select/out_xdrop_ni200.txt',dir_trunk)}; flag_run=1;
%elseif (cl_num==4 & flag_reverse==0 & n_maf==5 & strcmp(mr_string,'') & strcmp(mc_string,'')); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p25_D_m2r2_g004/dir_select/out_xdrop_ni175.txt',dir_trunk)}; flag_run=1;
elseif (cl_num==4 & flag_reverse==0 & n_maf==5 & strcmp(mr_string,'') & strcmp(mc_string,'')); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p25_D_m2r2_g004/dir_select/out_xdrop_ni238.txt',dir_trunk)}; flag_run=1;
%%%%%%%%;
elseif (cl_num==4 & flag_reverse==0 & n_maf==1 & strcmp(mr_string,'BDX') & strcmp(mc_string,'BDX')); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p01_BDX_BDX_D_m2r2_g004/dir_select/out_xdrop_ni125.txt',dir_trunk)}; flag_run=1;
elseif (cl_num==4 & flag_reverse==0 & n_maf==6 & strcmp(mr_string,'BDX') & strcmp(mc_string,'BDX')); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p05_BDX_BDX_D_m2r2_g004/dir_select/out_xdrop_ni138.txt',dir_trunk)}; flag_run=1;
elseif (cl_num==4 & flag_reverse==0 & n_maf==4 & strcmp(mr_string,'BDX') & strcmp(mc_string,'BDX')); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p10_BDX_BDX_D_m2r2_g004/dir_select/out_xdrop_ni163.txt',dir_trunk)}; flag_run=1;
elseif (cl_num==4 & flag_reverse==0 & n_maf==5 & strcmp(mr_string,'BDX') & strcmp(mc_string,'BDX')); scramble_out_xdrop_ = {sprintf('%s/dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p25_BDX_BDX_D_m2r2_g004/dir_select/out_xdrop_ni150.txt',dir_trunk)}; flag_run=1;
%%%%%%%%;
 else;
disp(sprintf(' %% Warning! scramble not implemented in xxxcluster_PGC_uADZSZDA_ver15_dr.m'));
end;%if implemented ;
if flag_run; xxxcluster_PGC_uADZSZDA_ver15(dir_code,dir_trunk,string_prefix,maf_lo_threshold,maf_hi_threshold,mr_string,mc_string,flag_dex_vs_lak{1},flag_reverse,Icat,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,scramble_out_xdrop_,scramble_rseed_,nshuffle,flag_rerun,pbs_walltime,slurm_walltime,slurm_nnodes,slurm_tpn,slurm_memdecl,row_factor,col_factor,flag_verbose); end;
end;%if flag_scramble;
%%%%%%%%;
end;%if (flag_platform==1 through 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% shuffled runs. ;
% Note that there is no need to shuffle the scrambled runs. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_platform==2 | flag_platform==3);
for nshuffle=1:n_shuffle;
n_scramble=0;
xxxcluster_PGC_uADZSZDA_ver15(dir_code,dir_trunk,string_prefix,maf_lo_threshold,maf_hi_threshold,mr_string,mc_string,flag_dex_vs_lak{1},flag_reverse,Icat,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,{},[],nshuffle,flag_rerun,pbs_walltime,slurm_walltime,slurm_nnodes,slurm_tpn,slurm_memdecl,row_factor,col_factor,flag_verbose); 
end;%for nshuffle=1:n_shuffle; 
end;%if (flag_platform==2 | flag_platform==3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for n_cov=[1];
end;%for n_maf=[1:5];
end;%for cl_num=1:4;
end;%for flag_reverse=[0,1]; 
end;%for flag_dex_vs_lak = {'dex','lak'};
