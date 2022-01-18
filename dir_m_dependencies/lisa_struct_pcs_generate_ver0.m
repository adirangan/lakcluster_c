function [lisa_arm1] = lisa_struct_pcs_generate_ver0(specification) ;
%{

  %%%%%%%%;
  % Basic run. ;
  %%%%%%%%;
  dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PGC_20190328');
  dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c_dev');
  flag_dex_vs_lak = 'dex'; %<-- differentially expressed clustering. ;
  cl_num_arm1 = 4; %<-- train on platform 4. ;
  cl_num_arm1_mc_A = []; cl_num_arm1_mr_A_full_ = []; cl_num_arm1_mr_Z_full_ = []; cl_num_arm1_string = []; 
  flag_reverse = 0; %<-- forward bicluster (i.e., case-specific). ;
  gamma = [0.004]; %<-- gamma is the fraction eliminated per iteration. 000 implies a single patient eliminated per iteration. ;
  B_MLT = 34; n_mds = 20; %<-- accurate to 2^(-34), 20 total mds components (but only 2 used). ; 
  mr_string_arm1 = '';mc_string_arm1 = ''; %No special mc_string. ;
  for n_maf = [5]; %for n_maf = [1,6,4,5]; %<-- minor-allele-frequency cutoff 25-50;
  n_cov = 2; %<-- 2 covariates (mds-components) used, repeated twice. ;
  n_scramble = 0; %<-- no previous bicluster extracted/scrambled first ;
  n_shuffle = 128; %<-- number of random shuffles. ;
  flag_rerun=0; %<-- regenerate files.; 
  pca_b_mlt = 44; pca_tolerance = 1e-2; pca_rank = 2;
  n_iteration_stride = 25/2;
  for flag_reverse = [0];%for flag_reverse = 0:1;
  specification = struct();
  specification.dir_trunk = dir_trunk;
  specification.dir_code = dir_code;
  specification.mr_string_arm1 = mr_string_arm1;
  specification.mc_string_arm1 = mc_string_arm1;
  specification.cl_num_arm1 = cl_num_arm1;
  specification.cl_num_arm1_mc_A = cl_num_arm1_mc_A;
  specification.cl_num_arm1_mr_A_full_ = cl_num_arm1_mr_A_full_;
  specification.cl_num_arm1_mr_Z_full_ = cl_num_arm1_mr_Z_full_;
  specification.cl_num_arm1_string = cl_num_arm1_string;
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
  specification.pca_b_mlt = pca_b_mlt;
  specification.pca_tolerance = pca_tolerance;
  specification.pca_rank = pca_rank;
  specification.n_iteration_stride = n_iteration_stride;
  lisa_struct_pcs_generate_ver0(specification) ;
  end;%for flag_reverse = 0:1;
  end;%for n_maf = [1,4,5,6]; %n_maf = 5; %<-- minor-allele-frequency cutoff 25-50;

 %}

dir_trunk = specification.dir_trunk;
dir_code = specification.dir_code;
mr_string_arm1 = specification.mr_string_arm1;
mc_string_arm1 = specification.mc_string_arm1;
cl_num_arm1 = specification.cl_num_arm1;
cl_num_arm1_mc_A = specification.cl_num_arm1_mc_A;
cl_num_arm1_mr_A_full_ = specification.cl_num_arm1_mr_A_full_;
cl_num_arm1_mr_Z_full_ = specification.cl_num_arm1_mr_Z_full_;
cl_num_arm1_string = specification.cl_num_arm1_string;
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
pca_b_mlt = specification.pca_b_mlt;
pca_tolerance = specification.pca_tolerance;
pca_rank = specification.pca_rank;
n_iteration_stride = specification.n_iteration_stride;

%%%%%%%%;
lisa_arm1 = lisa_struct_make_ver0(mr_string_arm1,mc_string_arm1,cl_num_arm1,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
lisa_arm1 = lisa_struct_prefix_ver0(lisa_arm1,dir_code,dir_trunk); 
lisa_arm1.nshuffle = 0;  lisa_arm1 = lisa_struct_names_ver0(lisa_arm1); 
lisa_arm1 = lisa_struct_xdrop_ver0(lisa_arm1); lisa_arm1 = lisa_struct_mdsfam_ver0(lisa_arm1); 
%lisa_arm1 = lisa_struct_bim_ver0(lisa_arm1); %<-- this is large and takes a while to load. ;
lisa_arm1 = lisa_struct_mx_ver0(lisa_arm1); lisa_arm1 = lisa_struct_studyindex_ver0(lisa_arm1); 
lisa_arm1 = lisa_struct_trace_ver0(lisa_arm1);
%%%%%%%%;

%%%%%%%%;
lisa_arm1.dir_out_s0_pcs = sprintf('%s/dir_pcs%s',lisa_arm1.dir_out_s0,cl_num_arm1_string);
if (~exist(lisa_arm1.dir_out_s0_pcs,'dir')); disp(sprintf(' %% mkdir %s',lisa_arm1.dir_out_s0_pcs)); mkdir(lisa_arm1.dir_out_s0_pcs); end;
%%%%%%%%;
lisa_arm1.dir_out_s0_pcs_jpg = sprintf('%s/dir_pcs_jpg',lisa_arm1.dir_out_s0_pcs);
if (~exist(lisa_arm1.dir_out_s0_pcs_jpg,'dir')); disp(sprintf(' %% creating %s',lisa_arm1.dir_out_s0_pcs_jpg)); mkdir(lisa_arm1.dir_out_s0_pcs_jpg); end;
%%%%%%%%;
lisa_arm1.dir_out_s0_pcs_b16 = sprintf('%s/dir_pcs_b16',lisa_arm1.dir_out_s0_pcs);
if (~exist(lisa_arm1.dir_out_s0_pcs_b16,'dir')); disp(sprintf(' %% creating %s',lisa_arm1.dir_out_s0_pcs_b16)); mkdir(lisa_arm1.dir_out_s0_pcs_b16); end;
%%%%%%%%;
lisa_arm1.dir_out_s0_pcs_inp = sprintf('%s/dir_pcs_inp',lisa_arm1.dir_out_s0_pcs);
if (~exist(lisa_arm1.dir_out_s0_pcs_inp,'dir')); disp(sprintf(' %% creating %s',lisa_arm1.dir_out_s0_pcs_inp)); mkdir(lisa_arm1.dir_out_s0_pcs_inp); end;
%%%%%%%%;
lisa_arm1.dir_out_s0_pcs_mda = sprintf('%s/dir_pcs_mda',lisa_arm1.dir_out_s0_pcs);
if (~exist(lisa_arm1.dir_out_s0_pcs_mda,'dir')); disp(sprintf(' %% creating %s',lisa_arm1.dir_out_s0_pcs_mda)); mkdir(lisa_arm1.dir_out_s0_pcs_mda); end;
%%%%%%%%;
lisa_arm1.dir_out_s0_pcs_mat = sprintf('%s/dir_pcs_mat',lisa_arm1.dir_out_s0_pcs);
if (~exist(lisa_arm1.dir_out_s0_pcs_mat,'dir')); disp(sprintf(' %% creating %s',lisa_arm1.dir_out_s0_pcs_mat)); mkdir(lisa_arm1.dir_out_s0_pcs_mat); end;
%%%%%%%%;
lisa_arm1.dir_out_s0_pcs_shuffle = sprintf('%s/dir_pcs_mat',lisa_arm1.dir_out_s0_pcs);
if (~exist(lisa_arm1.dir_out_s0_pcs_mat,'dir')); disp(sprintf(' %% creating %s',lisa_arm1.dir_out_s0_pcs_mat)); mkdir(lisa_arm1.dir_out_s0_pcs_mat); end;
%%%%%%%%;

n_iteration = lisa_arm1.n_iteration;
niteration_ = round(0:n_iteration_stride:n_iteration); 
niteration_(1)=1; niteration_(end) = min(niteration_(end),n_iteration-1);
niteration_ = niteration_(length(niteration_):-1:1);

flag_verbose=0;
flag_found_all = 1;
for nshuffle=0:n_shuffle;
pcs_infix = sprintf('s%.4d',nshuffle);
pcs_fname_V_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_V_.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_AnV_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_AnV_.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_ZnV_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_ZnV_.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_rkeep_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_rkeep.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_ckeep_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_ckeep.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_auc19_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_auc19.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
if (...
    exist(pcs_fname_AnV_arm1_,'file') & ...
    exist(pcs_fname_ZnV_arm1_,'file') & ...
    exist(pcs_fname_ckeep_arm1_,'file') & ...
    exist(pcs_fname_auc19_arm1_,'file') & ...
    1);
if (flag_verbose>1); disp(sprintf(' %% nshuffle %d; found:\n%% %s\n%% %s\n%% %s\n%% %s\n%% %s',nshuffle,pcs_fname_V_arm1_,pcs_fname_AnV_arm1_,pcs_fname_ZnV_arm1_,pcs_fname_ckeep_arm1_,pcs_fname_auc19_arm1_)); end;% if (flag_verbose>1);
 else;
flag_found_all = 0;
disp(sprintf(' %% nshuffle %d; missing files',nshuffle));
end;% if exist;
end;%for nshuffle=0:n_shuffle;

if (flag_found_all==1); disp(sprintf(' %% pcs %s found all files',lisa_arm1.string_name_s0)); end;% if (flag_found_all==1);

if (flag_rerun==1 | flag_found_all==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Step through random shuffles. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
for nshuffle=0:n_shuffle;

flag_found_local = 0;
pcs_infix = sprintf('s%.4d',nshuffle);
pcs_fname_V_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_V_.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_AnV_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_AnV_.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_ZnV_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_ZnV_.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_rkeep_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_rkeep.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_ckeep_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_ckeep.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
pcs_fname_auc19_arm1_ = sprintf('%s/pcs_%s_k%d_B%d_auc19.mda',lisa_arm1.dir_out_s0_pcs_mda,pcs_infix,pca_rank,pca_b_mlt);
if (...
    exist(pcs_fname_AnV_arm1_,'file') & ...
    exist(pcs_fname_ZnV_arm1_,'file') & ...
    exist(pcs_fname_ckeep_arm1_,'file') & ...
    exist(pcs_fname_auc19_arm1_,'file') & ...
    1);
flag_found_local=1;
if (flag_verbose>1); disp(sprintf(' %% nshuffle %d; found:\n%% %s\n%% %s\n%% %s\n%% %s\n%% %s',nshuffle,pcs_fname_V_arm1_,pcs_fname_AnV_arm1_,pcs_fname_ZnV_arm1_,pcs_fname_ckeep_arm1_,pcs_fname_auc19_arm1_)); end;% if (flag_verbose>1);
 else;
flag_found_local = 0;
disp(sprintf(' %% nshuffle %d; missing files',nshuffle));
end;% if exist;

if (flag_rerun==0 & flag_found_local==1); disp(sprintf(' %% nshuffle %d: files found, not rerunning',nshuffle)); end;
if (flag_rerun==1 | flag_found_local==0); disp(sprintf(' %% nshuffle %d: files not found, rerunning',nshuffle));

lisa_arm1 = lisa_struct_xdrop_shuffle_ver0(lisa_arm1,nshuffle);

tmp_mc_A_arm1_ = zeros(size(lisa_arm1.mc_A_));
tmp_mc_A_arm1_(lisa_arm1.cdrop_a_shuffle_) = 1;
bitj=16;
tmpchar = sprintf('%s/%s_mc_A_arm1.b16',lisa_arm1.dir_out_s0_pcs_b16,lisa_arm1.string_name_s0);
if (~exist(tmpchar,'file')); disp(sprintf(' %% writing %s',tmpchar)); tutorial_binary_compress(bitj,tmp_mc_A_arm1_(:)>0,tmpchar); end;

tmp_mr_X_full_ = lisa_arm1.mr_A_full_ + lisa_arm1.mr_Z_full_;
tmp_mr_D_full_ = zeros(size(tmp_mr_X_full_)); 
tmp_mr_D_full_(lisa_arm1.rdrop_a_shuffle_) = 1;
tmp_mr_X_full_ = tmp_mr_X_full_ - tmp_mr_D_full_;
bitj=16;
if (flag_reverse==0); 
tmpchar = sprintf('%s/%s_mr_A_arm1_full_%s.b16',lisa_arm1.dir_out_s0_pcs_b16,lisa_arm1.string_name_s0,pcs_infix);
if (~exist(tmpchar,'file')); disp(sprintf(' %% writing %s',tmpchar)); tutorial_binary_compress(bitj,tmp_mr_D_full_(:)>0,tmpchar); end;
tmpchar = sprintf('%s/%s_mr_Z_arm1_full_%s.b16',lisa_arm1.dir_out_s0_pcs_b16,lisa_arm1.string_name_s0,pcs_infix);
if (~exist(tmpchar,'file')); disp(sprintf(' %% writing %s',tmpchar)); tutorial_binary_compress(bitj,tmp_mr_X_full_(:)>0,tmpchar); end;
end;%if (flag_reverse==0); 
if (flag_reverse==1); 
tmpchar = sprintf('%s/%s_mr_Z_arm1_full_%s.b16',lisa_arm1.dir_out_s0_pcs_b16,lisa_arm1.string_name_s0,pcs_infix);
if (~exist(tmpchar,'file')); disp(sprintf(' %% writing %s',tmpchar)); tutorial_binary_compress(bitj,tmp_mr_D_full_(:)>0,tmpchar); end;
tmpchar = sprintf('%s/%s_mr_A_arm1_full_%s.b16',lisa_arm1.dir_out_s0_pcs_b16,lisa_arm1.string_name_s0,pcs_infix);
if (~exist(tmpchar,'file')); disp(sprintf(' %% writing %s',tmpchar)); tutorial_binary_compress(bitj,tmp_mr_X_full_(:)>0,tmpchar); end;
end;%if (flag_reverse==1); 

if (flag_rerun==0 &  exist(pcs_fname_AnV_arm1_,'file') &  exist(pcs_fname_ZnV_arm1_,'file'));
disp(sprintf(' %% loading %s and %s',pcs_fname_AnV_arm1_,pcs_fname_ZnV_arm1_));
end;%if ( exist(pcs_fname_AnV_arm1_,'file') &  exist(pcs_fname_ZnV_arm1_,'file'));
if (flag_rerun==1 | ~exist(pcs_fname_AnV_arm1_,'file') | ~exist(pcs_fname_ZnV_arm1_,'file'));
disp(sprintf(' %% calculating %s, %s and %s',pcs_fname_V_arm1_,pcs_fname_AnV_arm1_,pcs_fname_ZnV_arm1_));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% setting up input file for pca_driver. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
d_inpre = sprintf('%s/%s',lisa_arm1.dir__in,lisa_arm1.string_prefix); 
d_oupre = sprintf('%s/%s',lisa_arm1.dir_out_s0,lisa_arm1.string_name_s0);
d_oupre_pcs = sprintf('%s/%s',lisa_arm1.dir_out_s0_pcs_b16,lisa_arm1.string_name_s0);
%pca_b_mlt = 44;
%pca_tolerance = 1e-2;
%pca_rank = 2;
flag_T = 0;
if (lisa_arm1.mds_repl< 1 | length(lisa_arm1.mds_used_)~=2); flag_T = 0; T_n_crop_cols = 1+length(lisa_arm1.mds_used_); end;
if (lisa_arm1.mds_repl>=1 & length(lisa_arm1.mds_used_)==2); flag_T = 1; end;
lisa_arm1.mds_str = sprintf('m%dr%d',length(lisa_arm1.mds_used_),lisa_arm1.mds_repl);
mds_kappa_squared = 1.0;
if (length(lisa_arm1.mds_used_)>0); mds_kappa_squared = textread(sprintf('%s_T_%s_kappa.txt',d_inpre,lisa_arm1.mds_str)); end;
T_n_cols = 1+length(lisa_arm1.mds_used_)*lisa_arm1.mds_repl;
Y_n_cols=0;
M_n_rows_ = zeros(lisa_arm1.Icat,1); M_n_cols = length(tmp_mc_A_arm1_);
for nb=0:lisa_arm1.Icat-1; M_n_rows_(1+nb) = length(lisa_arm1.mr_A__{1+nb});end;%for nb=0:lisa_arm1.Icat-1;
fname__in = sprintf('%s/%s_%s.in',lisa_arm1.dir_out_s0_pcs_inp,lisa_arm1.string_name_s0,pcs_infix);
fp = fopen(fname__in,'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fprintf(fp,'GLOBAL_verbose= %d;\n',flag_verbose);
fprintf(fp,'GLOBAL_thread_count= 15;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','pca_driver');
fprintf(fp,'GLOBAL_pca_infix= pcs_%s;\n',pcs_infix);
if (nshuffle==0); fprintf(fp,'GLOBAL_pca_out_xdrop= %s/%s;\n',lisa_arm1.dir_out_s0,'out_xdrop_a.txt'); end;
if (nshuffle> 0); fprintf(fp,'GLOBAL_pca_out_xdrop= %s/out_xdrop_a_s%.4d.txt;\n',lisa_arm1.dir_out_trace,nshuffle); end;
fprintf(fp,'GLOBAL_pca_iteration_num= %d;\n',length(niteration_));
fprintf(fp,'GLOBAL_pca_iteration_max= %d;\n',niteration_(1));
fprintf(fp,'GLOBAL_pca_iteration_min= %d;\n',niteration_(end));
fprintf(fp,'GLOBAL_pca_rank= %d;\n',pca_rank);
fprintf(fp,'GLOBAL_pca_tolerance= %0.16f;\n',pca_tolerance);
fprintf(fp,'GLOBAL_TEST_niter= 1024;\n');
fprintf(fp,'GLOBAL_NBINS= %d;\n',lisa_arm1.Icat);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',pca_b_mlt);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',lisa_arm1.Ireq);
if (flag_T==1); fprintf(fp,'GLOBAL_kappa_squared= %0.16f;\n',mds_kappa_squared); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%if (nshuffle==0); A_name_string = '_A_'; end;%if (nshuffle==0);
%if (nshuffle>0); A_name_string = '_A_shuffle_'; end;%if (nshuffle>0);
A_name_string = '_A_';
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_name_= ',lisa_arm1.Icat,d_inpre,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_t_name_= ',lisa_arm1.Icat,d_inpre,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_A_n_rows_= ',lisa_arm1.Icat,M_n_rows_);
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',M_n_cols);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',lisa_arm1.Icat,d_oupre_pcs,'_mr_Z_arm1_',sprintf('_%s.b16',pcs_infix)); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',lisa_arm1.Icat,d_oupre_pcs,'_mr_A_arm1_',sprintf('_%s.b16',pcs_infix)); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A_arm1.b16;\n',d_oupre_pcs);
fprintf(fp,'GLOBAL_A_p_name= %s/A_p.mda;\n',lisa_arm1.dir_out_s0);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_name_= ',lisa_arm1.Icat,d_inpre,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_t_name_= ',lisa_arm1.Icat,d_inpre,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_Z_n_rows_= ',lisa_arm1.Icat,M_n_rows_);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',lisa_arm1.Icat,d_oupre,'_mr_A_','.b16'); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',lisa_arm1.Icat,d_oupre,'_mr_Z_','.b16'); end;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
if (flag_T==0);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_crop_cols); 
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_T_n_name_= ',lisa_arm1.Icat,d_oupre,'_T_crop_','_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_T_t_name_= ',lisa_arm1.Icat,d_oupre,'_T_crop_','_t.b16');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T_crop.b16;\n',d_oupre);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_S_n_name_= ',lisa_arm1.Icat,d_oupre,'_T_crop_','_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_S_t_name_= ',lisa_arm1.Icat,d_oupre,'_T_crop_','_t.b16');
end;%if (flag_T==0);
if (flag_T==1);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_cols); 
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_T_n_name_= ',lisa_arm1.Icat,d_inpre,sprintf('_T_%s_',lisa_arm1.mds_str),'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_T_t_name_= ',lisa_arm1.Icat,d_inpre,sprintf('_T_%s_',lisa_arm1.mds_str),'_t.b16');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T_%s.b16;\n',d_oupre,lisa_arm1.mds_str);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_S_n_name_= ',lisa_arm1.Icat,d_inpre,sprintf('_T_%s_',lisa_arm1.mds_str),'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_S_t_name_= ',lisa_arm1.Icat,d_inpre,sprintf('_T_%s_',lisa_arm1.mds_str),'_t.b16');
end;%if (flag_T==1);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',lisa_arm1.dir_out_s0_pcs_mda);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lisa_struct_pcs_generate_ver0 on %s;\n',date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fclose(fp);
type(fname__in);
disp(sprintf(' %% fname__in:\n%s',fname__in));
command_str = sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in);
system(command_str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% setting up input file for pca_proj arm1. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
d_inpre = sprintf('%s/%s',lisa_arm1.dir__in,lisa_arm1.string_prefix); 
d_oupre = sprintf('%s/%s',lisa_arm1.dir_out_s0,lisa_arm1.string_name_s0);
d_oupre_pcs = sprintf('%s/%s',lisa_arm1.dir_out_s0_pcs_b16,lisa_arm1.string_name_s0);
%pca_b_mlt = 44;
%pca_tolerance = 1e-2;
%pca_rank = 2;
lisa_arm1.mds_str = sprintf('m%dr%d',length(lisa_arm1.mds_used_),lisa_arm1.mds_repl);
Y_n_cols=0;
M_n_rows_ = zeros(lisa_arm1.Icat,1); M_n_cols = length(tmp_mc_A_arm1_);
for nb=0:lisa_arm1.Icat-1; M_n_rows_(1+nb) = length(lisa_arm1.mr_A__{1+nb});end;%for nb=0:lisa_arm1.Icat-1;
fname__in = sprintf('%s/%s_%s.in',lisa_arm1.dir_out_s0_pcs_inp,lisa_arm1.string_name_s0,pcs_infix);
fp = fopen(fname__in,'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fprintf(fp,'GLOBAL_verbose= %d;\n',flag_verbose);
fprintf(fp,'GLOBAL_thread_count= 15;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','pca_proj_driver');
fprintf(fp,'GLOBAL_pca_infix= pcs_%s;\n',pcs_infix);
if (nshuffle==0); fprintf(fp,'GLOBAL_pca_out_xdrop= %s/%s;\n',lisa_arm1.dir_out_s0,'out_xdrop_a.txt'); end;
if (nshuffle> 0); fprintf(fp,'GLOBAL_pca_out_xdrop= %s/out_xdrop_a_s%.4d.txt;\n',lisa_arm1.dir_out_trace,nshuffle); end;
fprintf(fp,'GLOBAL_pca_V_= %s;\n',pcs_fname_V_arm1_);
fprintf(fp,'GLOBAL_pca_iteration_num= %d;\n',length(niteration_));
fprintf(fp,'GLOBAL_pca_iteration_max= %d;\n',niteration_(1));
fprintf(fp,'GLOBAL_pca_iteration_min= %d;\n',niteration_(end));
fprintf(fp,'GLOBAL_pca_rank= %d;\n',pca_rank);
fprintf(fp,'GLOBAL_pca_tolerance= %0.16f;\n',pca_tolerance);
fprintf(fp,'GLOBAL_TEST_niter= 1024;\n');
fprintf(fp,'GLOBAL_NBINS= %d;\n',lisa_arm1.Icat);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',pca_b_mlt);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',lisa_arm1.Ireq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%if (nshuffle==0); A_name_string = '_A_'; end;%if (nshuffle==0);
%if (nshuffle>0); A_name_string = '_A_shuffle_'; end;%if (nshuffle>0);
A_name_string = '_A_';
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_name_= ',lisa_arm1.Icat,d_inpre,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_t_name_= ',lisa_arm1.Icat,d_inpre,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_A_n_rows_= ',lisa_arm1.Icat,M_n_rows_);
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',M_n_cols);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',lisa_arm1.Icat,d_oupre,'_mr_Z_','.b16'); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',lisa_arm1.Icat,d_oupre,'_mr_A_','.b16'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A_arm1.b16;\n',d_oupre_pcs);
fprintf(fp,'GLOBAL_A_p_name= %s/A_p.mda;\n',lisa_arm1.dir_out_s0);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_name_= ',lisa_arm1.Icat,d_inpre,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_t_name_= ',lisa_arm1.Icat,d_inpre,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_Z_n_rows_= ',lisa_arm1.Icat,M_n_rows_);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',lisa_arm1.Icat,d_oupre,'_mr_A_','.b16'); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',lisa_arm1.Icat,d_oupre,'_mr_Z_','.b16'); end;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',lisa_arm1.dir_out_s0_pcs_mda);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by lisa_struct_pcs_generate_ver0 on %s;\n',date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fclose(fp);
type(fname__in);
disp(sprintf(' %% fname__in:\n%s',fname__in));
command_str = sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in);
system(command_str);
end;%if (~exist(pcs_fname_AnV_arm1_,'file') | ~exist(pcs_fname_ZnV_arm1_,'file'));

if ( exist(pcs_fname_V_arm1_,'file') &  exist(pcs_fname_AnV_arm1_,'file') &  exist(pcs_fname_ZnV_arm1_,'file')); 
disp(sprintf(' %% found:\n %% %s\n %% %s\n %% %s\n %% deleting:\n %% %s',pcs_fname_V_arm1_,pcs_fname_AnV_arm1_,pcs_fname_ZnV_arm1_,pcs_fname_V_arm1_)); 
command_str = sprintf(' rm -rf %s',pcs_fname_V_arm1_);
disp(command_str);
system(command_str);
end;%if ( exist(pcs_fname_V_arm1_,'file'));

pcs_DnV_arm1_ = mda_read_r8(pcs_fname_AnV_arm1_);
pcs_XnV_arm1_ = mda_read_r8(pcs_fname_ZnV_arm1_);
pcs_rkeep_arm1_ = mda_read_i4(pcs_fname_rkeep_arm1_);
pcs_ckeep_arm1_ = mda_read_i4(pcs_fname_ckeep_arm1_);
%%%%%%%%;
flag_plot=0;
if flag_plot;
for ni=1:length(pcs_ckeep_arm1_);
ckeep = pcs_ckeep_arm1_(ni);
niteration = niteration_(ni);
clf;
hold on;
tmp_B_ij_ = lisa_arm1.rdrop_a_shuffle_(end+1-lisa_arm1.r_rem_(niteration):end);
tmp_D_ij_ = lisa_arm1.rdrop_a_shuffle_(1:end-lisa_arm1.r_rem_(niteration));
plot(pcs_DnV_arm1_(tmp_B_ij_,1,ni),pcs_DnV_arm1_(tmp_B_ij_,2,ni),'r.');
plot(pcs_DnV_arm1_(tmp_D_ij_,1,ni),pcs_DnV_arm1_(tmp_D_ij_,2,ni),'g.');
plot(pcs_XnV_arm1_(find(tmp_mr_X_full_),1,ni),pcs_XnV_arm1_(find(tmp_mr_X_full_),2,ni),'b.');
hold off;
title(sprintf('ni %d niteration %d <-- ckeep %d',ni,niteration,ckeep));
drawnow();
pause();
end;%for ni=1:length(pcs_ckeep_arm1_);
end;%if flag_plot;
%%%%%%%%;
pcs_auc19_arm1_ = zeros(length(niteration_),1);
tmp_D_ij_ = find(tmp_mr_D_full_);
tmp_X_ij_ = find(tmp_mr_X_full_);
tmp_mds_X_ = lisa_arm1.mds_sort_(tmp_X_ij_,[1:6,19]); 
tmp_mds_D_ = lisa_arm1.mds_sort_(tmp_D_ij_,[1:6,19]);
for ni=1:length(niteration_);
niteration = niteration_(ni);
tmp_auc19_arm1 = cauc_0(pcs_XnV_arm1_(tmp_X_ij_,1,ni),pcs_DnV_arm1_(tmp_D_ij_,1,ni),tmp_mds_X_,tmp_mds_D_);
pcs_auc19_arm1_(ni) = tmp_auc19_arm1;
end;%for ni=1:length(niteration_);
%%%%%%%%;
flag_plot=0;
if flag_plot;
plot(niteration_,pcs_auc19_arm1_,'k.-','LineWidth',1);
xlim([1,n_iteration]); ylim([0.45,1.0]);
end;%if flag_plot;
%%%%%%%%;
disp(sprintf(' %% writing %s',pcs_fname_auc19_arm1_));
mda_write_d3_r8(pcs_auc19_arm1_,pcs_fname_auc19_arm1_);

end;%if (flag_rerun==1 | flag_found_local==0); disp(sprintf(' %% nshuffle %s, files not found, rerunning',nshuffle));

end;%for nshuffle=0:n_shuffle;

end;%if (flag_rerun==1 | flag_found_all==0);