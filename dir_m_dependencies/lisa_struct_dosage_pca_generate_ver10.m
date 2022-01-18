function lisa_struct_dosage_pca_generate_ver10(specification) ;
% This function loads dosage data (line-by-line) and calculates a running-dot-product between the dosage data and a pca-vector V calculated using the genotyped data from arm1. ;
% This is updated to consider either n_scramble==0 (i.e. r0) or n_scramble==1 (i.e., r1). ;
% In addition, the genotyped data from arm2 alone is extracted from the dosage data (i.e., geno) and compared against the pca-vector V calculated using only the genotyped data from intersect(arm1,arm2). ;
% try: ;
%{

  %%%%%%%%;
  % Basic run (not on lisa). ;
  %%%%%%%%;
  dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PGC_20190328');
  dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c_dev');
  flag_dex_vs_lak = 'dex'; %<-- differentially expressed clustering. ;
  cl_num_arm1 = 4; cl_num_arm2 = 1; %<-- train on platform 4, replicate on platform 1. ;
  cl_num_arm1_mc_A = []; cl_num_arm1_mr_A_full_ = []; cl_num_arm1_mr_Z_full_ = []; cl_num_arm1_string = []; 
  cl_num_arm2_mc_A = []; cl_num_arm2_mr_A_full_ = []; cl_num_arm2_mr_Z_full_ = []; cl_num_arm2_string = []; 
  flag_reverse = 0; %<-- forward bicluster (i.e., case-specific). ;
  gamma = [0.004]; %<-- gamma is the fraction eliminated per iteration. 000 implies a single patient eliminated per iteration. ;
  B_MLT = 34; n_mds = 20; mr_string_arm1 = '';mc_string_arm1 = ''; mr_string_arm2 = '';mc_string_arm2 = ''; %<-- accurate to 2^(-34), 20 total mds components (but only 2 used). ; No special mc_string. ;
  n_maf = 5; n_cov = 2; %<-- minor-allele-frequency cutoff 25-50, 2 covariates (mds-components) used, repeated twice. ;
  n_scramble = 1; n_shuffle = 0; %<-- no previous bicluster extracted/scrambled first, no random shuffling. ;
  flag_rerun=0; %<-- regenerate files.; 
  pca_b_mlt = 44; pca_tolerance = 1e-2; pca_rank = 2;
  n_iteration_stride = 25/2;
  flag_stage = 0;
  flag_lisa = 0;
  if n_scramble==0; niteration = 175; end;
  if n_scramble==1; niteration = 163; end;
  for flag_reverse = 0;%for flag_reverse = 0:1;
  for cl_num_arm2 = [2,3];%for cl_num_arm2 = 1:4;
  specification = struct();
  specification.dir_trunk = dir_trunk;
  specification.dir_code = dir_code;
  specification.mr_string_arm1 = mr_string_arm1;
  specification.mc_string_arm1 = mc_string_arm1;
  specification.mr_string_arm2 = mr_string_arm2;
  specification.mc_string_arm2 = mc_string_arm2;
  specification.cl_num_arm1 = cl_num_arm1;
  specification.cl_num_arm1_mc_A = cl_num_arm1_mc_A;
  specification.cl_num_arm1_mr_A_full_ = cl_num_arm1_mr_A_full_;
  specification.cl_num_arm1_mr_Z_full_ = cl_num_arm1_mr_Z_full_;
  specification.cl_num_arm1_string = cl_num_arm1_string;
  specification.cl_num_arm2 = cl_num_arm2;
  specification.cl_num_arm2_mc_A = cl_num_arm2_mc_A;
  specification.cl_num_arm2_mr_A_full_ = cl_num_arm2_mr_A_full_;
  specification.cl_num_arm2_mr_Z_full_ = cl_num_arm2_mr_Z_full_;
  specification.cl_num_arm2_string = cl_num_arm2_string;
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
  specification.flag_stage = flag_stage;
  specification.flag_lisa = flag_lisa;
  specification.niteration = niteration;
  lisa_struct_dosage_pca_generate_ver10(specification) ;
  end;%for cl_num_arm2 = 1:4;
  end;%for flag_reverse = 0:1;

  %%%%%%%%;
  % Basic run (on lisa). ;
  %%%%%%%%;
  dir_trunk = sprintf('/home/grouther/dir_PGC_20190328');
  dir_code = sprintf('/home/grouther/dir_lakcluster_c_dev');
  flag_dex_vs_lak = 'dex'; %<-- differentially expressed clustering. ;
  cl_num_arm1 = 4; cl_num_arm2 = 1; %<-- train on platform 4, replicate on platform 1. ;
  cl_num_arm1_mc_A = []; cl_num_arm1_mr_A_full_ = []; cl_num_arm1_mr_Z_full_ = []; cl_num_arm1_string = []; 
  cl_num_arm2_mc_A = []; cl_num_arm2_mr_A_full_ = []; cl_num_arm2_mr_Z_full_ = []; cl_num_arm2_string = []; 
  flag_reverse = 0; %<-- forward bicluster (i.e., case-specific). ;
  gamma = [0.004]; %<-- gamma is the fraction eliminated per iteration. 000 implies a single patient eliminated per iteration. ;
  B_MLT = 34; n_mds = 20; mr_string_arm1 = '';mc_string_arm1 = ''; mr_string_arm2 = '';mc_string_arm2 = ''; %<-- accurate to 2^(-34), 20 total mds components (but only 2 used). ; No special mc_string. ;
  n_maf = 5; n_cov = 2; %<-- minor-allele-frequency cutoff 25-50, 2 covariates (mds-components) used, repeated twice. ;
  n_scramble = 1; n_shuffle = 0; %<-- no previous bicluster extracted/scrambled first, no random shuffling. ;
  flag_rerun=0; %<-- regenerate files.; 
  pca_b_mlt = 44; pca_tolerance = 1e-2; pca_rank = 2;
  n_iteration_stride = 25/2;
  flag_stage = 1;
  flag_lisa = 1;
  niteration = 163;
  for flag_reverse = 0;%for flag_reverse = 0:1;
  for cl_num_arm2 = 1;%for cl_num_arm2 = 1:4;
  specification = struct();
  specification.dir_trunk = dir_trunk;
  specification.dir_code = dir_code;
  specification.mr_string_arm1 = mr_string_arm1;
  specification.mc_string_arm1 = mc_string_arm1;
  specification.mr_string_arm2 = mr_string_arm2;
  specification.mc_string_arm2 = mc_string_arm2;
  specification.cl_num_arm1 = cl_num_arm1;
  specification.cl_num_arm1_mc_A = cl_num_arm1_mc_A;
  specification.cl_num_arm1_mr_A_full_ = cl_num_arm1_mr_A_full_;
  specification.cl_num_arm1_mr_Z_full_ = cl_num_arm1_mr_Z_full_;
  specification.cl_num_arm1_string = cl_num_arm1_string;
  specification.cl_num_arm2 = cl_num_arm2;
  specification.cl_num_arm2_mc_A = cl_num_arm2_mc_A;
  specification.cl_num_arm2_mr_A_full_ = cl_num_arm2_mr_A_full_;
  specification.cl_num_arm2_mr_Z_full_ = cl_num_arm2_mr_Z_full_;
  specification.cl_num_arm2_string = cl_num_arm2_string;
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
  specification.flag_stage = flag_stage;
  specification.flag_lisa = flag_lisa;
  specification.niteration = niteration;
  lisa_struct_dosage_pca_generate_ver10(specification) ;
  end;%for cl_num_arm2 = 1:4;
  end;%for flag_reverse = 0:1;

  %}

% Note the following asymmetry regarding sparse matrices: ;
% summing over columns is *much* faster than summing over rows. ;
% (i.e., about 25 times faster for the example below). ;
%{
tmp_N = 1024*8;
tmp_ = sparse(randperm(tmp_N),randperm(tmp_N),1,tmp_N,tmp_N);
disp('sum(*,2)'); tic;for nl=1:1024;sum(tmp_(:,max(1,min(tmp_N,floor(tmp_N*rand())))),2); end;toc;
disp('sum(*,1)'); tic;for nl=1:1024;sum(tmp_(max(1,min(tmp_N,floor(tmp_N*rand()))),:),1); end;toc;
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
mr_string_arm2 = specification.mr_string_arm2;
mc_string_arm2 = specification.mc_string_arm2;
cl_num_arm2 = specification.cl_num_arm2;
cl_num_arm2_mc_A = specification.cl_num_arm2_mc_A;
cl_num_arm2_mr_A_full_ = specification.cl_num_arm2_mr_A_full_;
cl_num_arm2_mr_Z_full_ = specification.cl_num_arm2_mr_Z_full_;
cl_num_arm2_string = specification.cl_num_arm2_string;
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
flag_stage = specification.flag_stage;
flag_lisa = specification.flag_lisa;
flag_lisa_use = specification.flag_lisa;
niteration = specification.niteration;
n_figure = 1;

%%%%%%%%;
lisa_arm1 = lisa_struct_make_ver0(mr_string_arm1,mc_string_arm1,cl_num_arm1,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,n_scramble,n_shuffle) ;
lisa_arm1 = lisa_struct_prefix_ver0(lisa_arm1,dir_code,dir_trunk); 
lisa_arm1.nshuffle = 0;  lisa_arm1 = lisa_struct_names_ver0(lisa_arm1); 
lisa_arm1 = lisa_struct_xdrop_ver0(lisa_arm1); lisa_arm1 = lisa_struct_mdsfam_ver0(lisa_arm1); 
if (flag_stage==0); disp(sprintf(' %% loading bim')); lisa_arm1 = lisa_struct_bim_ver0(lisa_arm1); end; %<-- this is large and takes a while to load. ;
lisa_arm1 = lisa_struct_mx_ver0(lisa_arm1); lisa_arm1 = lisa_struct_studyindex_ver0(lisa_arm1); 
lisa_arm1 = lisa_struct_trace_ver0(lisa_arm1);
%%%%%%%%;
lisa_arm2 = lisa_struct_make_ver0(mr_string_arm2,mc_string_arm2,cl_num_arm2,flag_dex_vs_lak,gamma,B_MLT,n_mds,flag_reverse,n_maf,n_cov,0,n_shuffle) ;
lisa_arm2 = lisa_struct_prefix_ver0(lisa_arm2,dir_code,dir_trunk); 
lisa_arm2.nshuffle = 0;  lisa_arm2 = lisa_struct_names_ver0(lisa_arm2); 
lisa_arm2 = lisa_struct_xdrop_ver0(lisa_arm2); lisa_arm2 = lisa_struct_mdsfam_ver0(lisa_arm2); 
if (flag_stage==0); disp(sprintf(' %% loading bim')); lisa_arm2 = lisa_struct_bim_ver0(lisa_arm2); end; %<-- this is large and takes a while to load. ;
lisa_arm2 = lisa_struct_mx_ver0(lisa_arm2); lisa_arm2 = lisa_struct_studyindex_ver0(lisa_arm2); 
lisa_arm2 = lisa_struct_trace_ver0(lisa_arm2);
%%%%%%%%;

lisa_arm1.dir_out_s0_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1.dir_out_s0,lisa_arm1.cl_num,cl_num_arm1_string,lisa_arm2.cl_num,cl_num_arm2_string);
if (~exist(lisa_arm1.dir_out_s0_pca,'dir')); disp(sprintf(' %% mkdir %s',lisa_arm1.dir_out_s0_pca)); mkdir(lisa_arm1.dir_out_s0_pca); end;
lisa_arm1.dir_out_s0_pca_jpg = sprintf('%s/dir_jpg',lisa_arm1.dir_out_s0_pca);
if (~exist(lisa_arm1.dir_out_s0_pca_jpg,'dir')); disp(sprintf(' %% creating %s',lisa_arm1.dir_out_s0_pca_jpg)); mkdir(lisa_arm1.dir_out_s0_pca_jpg); end;
n_iteration = lisa_arm1.n_iteration;
niteration_ = round(0:n_iteration_stride:n_iteration); 
niteration_(1)=1; niteration_(end) = min(niteration_(end),n_iteration-1);
niteration_ = niteration_(length(niteration_):-1:1);

flag_verbose=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_stage==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
%extract Ap ;
%%%%%%%%%%%%%%%%;
fname_tmp = sprintf('%s/A_p.mda',lisa_arm1.dir_out_s0); fcheck(fname_tmp);
lisa_arm1.A_p_ = mda_read_r8(fname_tmp);
fname_tmp = sprintf('%s/AZ_rsum.mda',lisa_arm1.dir_out_s0); fcheck(fname_tmp);
lisa_arm1.AZ_rsum_ = mda_read_r8(fname_tmp);
flag_check=0;
if flag_check;
disp(sprintf(' %% comparing atfr_, AZ_rsum_, A_p_ for arm1'));
tmp_ = 1:length(lisa_arm1.bim_{1}); plot(tmp_,lisa_arm1.bim_{6},'k-',tmp_,lisa_arm1.AZ_rsum_/lisa_arm1.n_patient,'ro',tmp_,lisa_arm1.A_p_(1+floor((tmp_-1)/1920)),'gx'); %<-- These should all overlap. ;
end;%if flag_check;
fname_tmp = sprintf('%s/A_p.mda',lisa_arm2.dir_out_s0); fcheck(fname_tmp);
lisa_arm2.A_p_ = mda_read_r8(fname_tmp);
disp(sprintf(' %% Warning! using unstable A_p'));
lisa_arm2.A_p_pca_ = mda_read_r8(sprintf('%s/A_p.mda',lisa_arm1.dir_out_s0_pca));
fname_tmp = sprintf('%s/AZ_rsum.mda',lisa_arm2.dir_out_s0); fcheck(fname_tmp);
lisa_arm2.AZ_rsum_ = mda_read_r8(fname_tmp);
disp(sprintf(' %% Warning! using unstable AZ_rsum'));
lisa_arm2.AZ_rsum_pca_ = mda_read_r8(sprintf('%s/AZ_rsum.mda',lisa_arm1.dir_out_s0_pca));
flag_check=0;
if flag_check;
disp(sprintf(' %% comparing atfr_, AZ_rsum_, A_p_ for arm2'));
subplot(1,2,1);tmp_ = 1:length(lisa_arm2.bim_{1}); plot(tmp_,lisa_arm2.bim_{6},'k-',tmp_,lisa_arm2.AZ_rsum_/lisa_arm2.n_patient,'ro',tmp_,lisa_arm2.A_p_(1+floor((tmp_-1)/1920)),'gx'); %<-- These should all overlap. ;
subplot(1,2,2);tmp_ = 1:length(lisa_arm2.bim_{1}); plot(tmp_,lisa_arm2.bim_{6},'k-',tmp_,lisa_arm2.AZ_rsum_pca_/lisa_arm2.n_patient,'ro',tmp_,lisa_arm2.A_p_pca_(1+floor((tmp_-1)/1920)),'gx'); %<-- These should all overlap. ;
end;%if flag_check;

%%%%%%%%%%%%%%%%;
%pick niteration. ;
%%%%%%%%%%%%%%%%;
disp(sprintf(' %% niteration %d',niteration));
%%%%%%%%;

tmp_dir_pca = sprintf('%s/dir_pca_trn%d%s_tst%d%s',lisa_arm1.dir_out_s0,lisa_arm1.cl_num,cl_num_arm1_string,lisa_arm1.cl_num,cl_num_arm1_string);
if (~exist(tmp_dir_pca,'dir')); disp(sprintf(' %% Warning! %s does not exist in lisa_struct_dosage_pca_generate_ver10.m ',tmp_dir_pca)); end;
lisa_arm1.pca_infix = sprintf('pca_ni%d_tst%d',niteration,cl_num_arm1);
V_trnx_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_V_.mda',tmp_dir_pca,lisa_arm1.pca_infix,pca_rank,pca_b_mlt)); %<-- This is the V-vector associated with (all) the snps in arm1, ordered to apply to arm1. ;
%%%%%%%%;
pca_infix_arm1 = sprintf('pca_ni%d_tst%d',niteration,cl_num_arm2);
V_tstx_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_V_.mda',lisa_arm1.dir_out_s0_pca,pca_infix_arm1,pca_rank,pca_b_mlt)); %<-- This is the V-vector associated with the snps in both arm1 and arm2, ordered to apply to arm1. ;
V_tstx_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_V_arm2_.mda',lisa_arm1.dir_out_s0_pca,pca_infix_arm1,pca_rank,pca_b_mlt)); %<-- This is the V-vector associated with the snps in both arm1 and arm2, ordered to apply to arm2. ;
% Note that V_tstx_arm1_ and V_tstx_arm2_ are simply permutations of one another. tmp_ij1_ = find(V_tstx_arm1_(:,1)~=0); tmp_ij2_ = find(V_tstx_arm2_(:,1)~=0); plot(sort(V_tstx_arm1_(tmp_ij1_,1))-sort(V_tstx_arm2_(tmp_ij2_,1))); clear tmp_ij1_ tmp_ij2_ ; %<-- should be all zeros. ;
% While V_trnx_arm1_ is not a permutation of V_tstx_arm1_: The former is calculated using the entire bicluster, whereas the latter is calculated using only the snps in both arm1 and arm2. Nevertheless: these two vectors are quite correlated. ; plot(V_trnx_arm1_(:,1),V_tstx_arm1_(:,1),'.') ; %<-- should be highly correlated. ;
%%%%%%%%;
pca_proj_infix_arm1 = sprintf('pca_proj_arm1_ni%d_tst%d',niteration,cl_num_arm2); %<-- Note that here we use the projection of arm1, with snps limited by overlap with arm two ;
DnV_tstx_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_AnV_.mda',lisa_arm1.dir_out_s0_pca,pca_proj_infix_arm1,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca from arm one, since this is where the bicluster is. ;
XnV_tstx_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_ZnV_.mda',lisa_arm1.dir_out_s0_pca,pca_proj_infix_arm1,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca from arm one, since this is where the bicluster is. ;
%if (strcmp(flag_dex_vs_lak,'lak')); DnV_tstx_arm1_ = abs(DnV_tstx_arm1_); XnV_tstx_arm1_ = abs(XnV_tstx_arm1_); end;%if (strcmp(flag_dex_vs_lak,'lak'));
mr_D_rmv_arm1_ = lisa_arm1.mr_D_*0; mr_D_rmv_arm1_(lisa_arm1.rdrop_a_(1:end-lisa_arm1.trace_{1}(niteration,2)))=1;
mr_D_ret_arm1_ = lisa_arm1.mr_D_*0; mr_D_ret_arm1_(lisa_arm1.rdrop_a_(end-lisa_arm1.trace_{1}(niteration,2)+1:end))=1;
%%%%%%%%;
pca_proj_infix_arm2 = sprintf('pca_proj_arm2_ni%d_tst%d',niteration,cl_num_arm2); %<-- Note that here we use the projection of arm2, with snps limited by overlap with arm two ;
DnV_tstx_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_AnV_.mda',lisa_arm1.dir_out_s0_pca,pca_proj_infix_arm2,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca from arm one, since this is where the bicluster is. ;
XnV_tstx_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_ZnV_.mda',lisa_arm1.dir_out_s0_pca,pca_proj_infix_arm2,pca_rank,pca_b_mlt)); %<-- Note that we use dir_out_s0000_pca from arm one, since this is where the bicluster is. ;
%if (strcmp(flag_dex_vs_lak,'lak')); DnV_tstx_arm2_ = abs(DnV_tstx_arm2_); XnV_tstx_arm2_ = abs(XnV_tstx_arm2_); end;%if (strcmp(flag_dex_vs_lak,'lak')); %<-- this is a crude way of collapsing multidimensional information. ;

%%%%%%%%;
% establish snp-indices for each allelic combination. ;
%%%%%%%%;
%%%%%%%%;
% First arm1: ;
%%%%%%%%;
bo1_ = lisa_arm1.bim_{1}; allele1_arm1_ = lisa_arm1.bim_{2}; allele2_arm1_ = lisa_arm1.bim_{3}; alleletype_arm1_ = lisa_arm1.bim_{4}; atfr_arm1_ = lisa_arm1.bim_{6};
[bu1_,bo1_to_bu1_,bu1_to_bo1_] = unique(bo1_);
bo1_by_bu1_xref_ = sparse(1:length(bo1_),bu1_to_bo1_,1,length(bo1_),length(bu1_));
%%%%%%%%;
flag_check=0;
if flag_check;
disp(sprintf(' %% checking bo1_by_bu1_xref_'));
nbu1=max(1,min(length(bu1_),floor(length(bu1_)*rand())));
disp(sprintf(' %% nbu1 %.9d: bu1_{nbu1} %s',nbu1,bu1_{nbu1}));
tmp_bo1_ = find(sum(bo1_by_bu1_xref_(:,nbu1),2));
for nl01=1:length(tmp_bo1_);
nbo1 = tmp_bo1_(nl01);
disp(sprintf(' %% nbo1 %.9d: bo1_{nbo1} %s (%s)',nbo1,bo1_{nbo1},alleletype_arm1_{nbo1}));
assert(strcmp(bo1_{nbo1},bu1_{nbu1}));
end;%for nl01=1:length(tmp_bo1_);
end;%if flag_check;
%%%%%%%%;
tmp_ = 1:length(bo1_); POPLENGTH = 1920; tmp_A_p_arm1_ = lisa_arm1.A_p_(1+floor((tmp_-1)/POPLENGTH));
flag_check = 0; [atfr_arm1_and_,A_p_arm1_and_,V_trnx_arm1_and_,V_tstx_arm1_and_] = lisa_dosage_pca_helper_2('and',alleletype_arm1_,bo1_,bu1_,atfr_arm1_,tmp_A_p_arm1_,flag_check,bo1_by_bu1_xref_,V_trnx_arm1_,V_tstx_arm1_);
flag_check = 0; [atfr_arm1_xor_,A_p_arm1_xor_,V_trnx_arm1_xor_,V_tstx_arm1_xor_] = lisa_dosage_pca_helper_2('xor',alleletype_arm1_,bo1_,bu1_,atfr_arm1_,tmp_A_p_arm1_,flag_check,bo1_by_bu1_xref_,V_trnx_arm1_,V_tstx_arm1_);
flag_check = 0; [atfr_arm1_nor_,A_p_arm1_nor_,V_trnx_arm1_nor_,V_tstx_arm1_nor_] = lisa_dosage_pca_helper_2('nor',alleletype_arm1_,bo1_,bu1_,atfr_arm1_,tmp_A_p_arm1_,flag_check,bo1_by_bu1_xref_,V_trnx_arm1_,V_tstx_arm1_);
%%%%%%%%;
% Plotting allelic-types against one another. ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
tmp_ij_ = find( ~ ( (V_trnx_arm1_and_(:,1)==0) & (V_trnx_arm1_xor_(:,1)==0) & (V_trnx_arm1_nor_(:,1)==0) ) );
plot3(V_trnx_arm1_and_(tmp_ij_,1),V_trnx_arm1_xor_(tmp_ij_,1),V_trnx_arm1_nor_(tmp_ij_,1),'.');
xlabel('and'); ylabel('xor'); zlabel('nor'); axis vis3d;
end;% if flag_plot;
flag_plot=0;
if flag_plot;
tmp_ij_ = find( ~ ( (V_trnx_arm1_and_(:,1)==0) & (V_trnx_arm1_xor_(:,1)==0) & (V_trnx_arm1_nor_(:,1)==0) ) );
subplot(3,3,1); plot(V_trnx_arm1_and_(tmp_ij_,1),V_trnx_arm1_xor_(tmp_ij_,1),'b.'); xlabel('and'); ylabel('xor'); title('and-vs-xor 1');
subplot(3,3,2); plot(V_trnx_arm1_and_(tmp_ij_,1),V_trnx_arm1_nor_(tmp_ij_,1),'b.'); xlabel('and'); ylabel('nor'); title('and-vs-nor 1');
subplot(3,3,3); plot(V_trnx_arm1_xor_(tmp_ij_,1),V_trnx_arm1_nor_(tmp_ij_,1),'b.'); xlabel('xor'); ylabel('nor'); title('xor-vs-nor 1');
subplot(3,3,4); plot(V_trnx_arm1_and_(tmp_ij_,2),V_trnx_arm1_xor_(tmp_ij_,2),'r.'); xlabel('and'); ylabel('xor'); title('and-vs-xor 2');
subplot(3,3,5); plot(V_trnx_arm1_and_(tmp_ij_,2),V_trnx_arm1_nor_(tmp_ij_,2),'r.'); xlabel('and'); ylabel('nor'); title('and-vs-nor 2');
subplot(3,3,6); plot(V_trnx_arm1_xor_(tmp_ij_,2),V_trnx_arm1_nor_(tmp_ij_,2),'r.'); xlabel('xor'); ylabel('nor'); title('xor-vs-nor 2');
subplot(3,3,7); plot(V_trnx_arm1_and_(tmp_ij_,1),V_trnx_arm1_and_(tmp_ij_,2),'g.'); xlabel('and 1'); ylabel('and 2'); title('and 1-vs-2');
subplot(3,3,8); plot(V_trnx_arm1_xor_(tmp_ij_,1),V_trnx_arm1_xor_(tmp_ij_,2),'g.'); xlabel('xor 1'); ylabel('xor 2'); title('xor 1-vs-2');
subplot(3,3,9); plot(V_trnx_arm1_nor_(tmp_ij_,1),V_trnx_arm1_nor_(tmp_ij_,2),'g.'); xlabel('nor 1'); ylabel('nor 2'); title('nor 1-vs-2');
end;% if flag_plot;
flag_plot=0;
if flag_plot;
tmp_ij_ = find( ~ ( (V_trnx_arm1_and_(:,1)==0) & (V_trnx_arm1_xor_(:,1)==0) & (V_trnx_arm1_nor_(:,1)==0) ) );
subplot(3,1,1); plot(atfr_arm1_and_(tmp_ij_,1),atfr_arm1_xor_(tmp_ij_,1),'b.'); xlabel('and'); ylabel('xor'); title('and-vs-xor');
subplot(3,1,2); plot(atfr_arm1_and_(tmp_ij_,1),atfr_arm1_nor_(tmp_ij_,1),'b.'); xlabel('and'); ylabel('nor'); title('and-vs-nor');
subplot(3,1,3); plot(atfr_arm1_xor_(tmp_ij_,1),atfr_arm1_nor_(tmp_ij_,1),'b.'); xlabel('xor'); ylabel('nor'); title('xor-vs-nor');
end;% if flag_plot;
%%%%%%%%;
% Now arm2: ;
%%%%%%%%;
bo2_ = lisa_arm2.bim_{1}; allele1_arm2_ = lisa_arm2.bim_{2}; allele2_arm2_ = lisa_arm2.bim_{3}; alleletype_arm2_ = lisa_arm2.bim_{4}; atfr_arm2_ = lisa_arm2.bim_{6};
[bu2_,bo2_to_bu2_,bu2_to_bo2_] = unique(bo2_);
bo2_by_bu2_xref_ = sparse(1:length(bo2_),bu2_to_bo2_,1,length(bo2_),length(bu2_));
%%%%%%%%;
flag_check=0;
if flag_check;
disp(sprintf(' %% checking bo2_by_bu2_xref_'));
nbu2=max(1,min(length(bu2_),floor(length(bu2_)*rand())));
disp(sprintf(' %% nbu2 %.9d: bu2_{nbu2} %s',nbu2,bu2_{nbu2}));
tmp_bo2_ = find(sum(bo2_by_bu2_xref_(:,nbu2),2));
for nlo2=1:length(tmp_bo2_);
nbo2 = tmp_bo2_(nlo2);
disp(sprintf(' %% nbo2 %.9d: bo2_{nbo2} %s (%s)',nbo2,bo2_{nbo2},alleletype_arm2_{nbo2}));
assert(strcmp(bo2_{nbo2},bu2_{nbu2}));
end;%for nlo2=1:length(tmp_bo2_);
end;%if flag_check;
%%%%%%%%;
tmp_ = 1:length(bo2_); POPLENGTH = 1920; tmp_A_p_arm2_ = lisa_arm2.A_p_(1+floor((tmp_-1)/POPLENGTH));
flag_check = 0; [atfr_arm2_and_,A_p_arm2_and_] = lisa_dosage_pca_helper_2('and',alleletype_arm2_,bo2_,bu2_,atfr_arm2_,tmp_A_p_arm2_,flag_check,bo2_by_bu2_xref_);
flag_check = 0; [atfr_arm2_xor_,A_p_arm2_xor_] = lisa_dosage_pca_helper_2('xor',alleletype_arm2_,bo2_,bu2_,atfr_arm2_,tmp_A_p_arm2_,flag_check,bo2_by_bu2_xref_);
flag_check = 0; [atfr_arm2_nor_,A_p_arm2_nor_] = lisa_dosage_pca_helper_2('nor',alleletype_arm2_,bo2_,bu2_,atfr_arm2_,tmp_A_p_arm2_,flag_check,bo2_by_bu2_xref_);
%%%%%%%%;
% Plotting allelic-types against one another. ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
subplot(3,1,1); plot(atfr_arm2_and_(:,1),atfr_arm2_xor_(:,1),'b.'); xlabel('and'); ylabel('xor'); title('and-vs-xor');
subplot(3,1,2); plot(atfr_arm2_and_(:,1),atfr_arm2_nor_(:,1),'b.'); xlabel('and'); ylabel('nor'); title('and-vs-nor');
subplot(3,1,3); plot(atfr_arm2_xor_(:,1),atfr_arm2_nor_(:,1),'b.'); xlabel('xor'); ylabel('nor'); title('xor-vs-nor');
end;% if flag_plot;
%%%%%%%%;
% Now intersect arm1 and arm2;
%%%%%%%%;
[bu3_,bu3_to_bu1_,bu3_to_bu2_] = intersect(bu1_,bu2_,'stable');
bu1_by_bu3_xref_ = sparse(bu3_to_bu1_,1:length(bu3_),1,length(bu1_),length(bu3_));
bu3_by_bu1_xref_ = sparse(1:length(bu3_),bu3_to_bu1_,1,length(bu3_),length(bu1_));
bu2_by_bu3_xref_ = sparse(bu3_to_bu2_,1:length(bu3_),1,length(bu2_),length(bu3_));
bu3_by_bu2_xref_ = sparse(1:length(bu3_),bu3_to_bu2_,1,length(bu3_),length(bu2_));
flag_check=0;
if flag_check;
disp(sprintf(' %% checking bu3_ vs bu1_ and bu2_'));
%%%%%%%%;
nbu3=max(1,min(length(bu3_),floor(length(bu3_)*rand())));
disp(sprintf(' %% nbu3 %.9d: bu3_{nbu3} %s',nbu3,bu3_{nbu3}));
%%%%%%%%;
tmp_bu1_ = find(sum(bu1_by_bu3_xref_(:,nbu3),2));
for nlu1=1:length(tmp_bu1_);
nbu1 = tmp_bu1_(nlu1);
disp(sprintf(' %% nbu1 %.9d: bu1_{nbu1} %s',nbu1,bu1_{nbu1}));
assert(strcmp(bu1_{nbu1},bu3_{nbu3}));
tmp_bo1_ = find(sum(bo1_by_bu1_xref_(:,nbu1),2));
for nlo1=1:length(tmp_bo1_);
nbo1 = tmp_bo1_(nlo1);
disp(sprintf(' %% nbo1 %.9d: bo1_{nbo1} %s (%s)',nbo1,bo1_{nbo1},alleletype_arm1_{nbo1}));
assert(strcmp(bo1_{nbo1},bu3_{nbu3}));
end;%for nlo1=1:length(tmp_bo1_);
end;%for nlu1=1:length(tmp_bu1_);
%%%%%%%%;
tmp_bu2_ = find(sum(bu2_by_bu3_xref_(:,nbu3),2));
for nlu2=1:length(tmp_bu2_);
nbu2 = tmp_bu2_(nlu2);
disp(sprintf(' %% nbu2 %.9d: bu2_{nbu2} %s',nbu2,bu2_{nbu2}));
assert(strcmp(bu2_{nbu2},bu3_{nbu3}));
tmp_bo2_ = find(sum(bo2_by_bu2_xref_(:,nbu2),2));
for nlo2=1:length(tmp_bo2_);
nbo2 = tmp_bo2_(nlo2);
disp(sprintf(' %% nbo2 %.9d: bo2_{nbo2} %s (%s)',nbo2,bo2_{nbo2},alleletype_arm2_{nbo2}));
assert(strcmp(bo2_{nbo2},bu3_{nbu3}));
end;%for nlo2=1:length(tmp_bo2_);
end;%for nlu2=1:length(tmp_bu2_);
%%%%%%%%;
end;%if flag_check;

flag_check=0;
if flag_check;
disp(sprintf(' %% checking V_tstx_arm1_ against V_tstx_arm2_'));
tmp_tstx_arm1_and_ = zeros(length(bu3_),2); tmp_tstx_arm1_xor_ = zeros(length(bu3_),2); tmp_tstx_arm1_nor_ = zeros(length(bu3_),2);
tmp_tstx_arm2_and_ = zeros(length(bu3_),2); tmp_tstx_arm2_xor_ = zeros(length(bu3_),2); tmp_tstx_arm2_nor_ = zeros(length(bu3_),2);
for nbu3=1:length(bu3_);
%%%%%%%%;
nbu1 = find(sum(bu1_by_bu3_xref_(:,nbu3),2)); assert(length(nbu1)==1);
tmp_bo1_ = find(sum(bo1_by_bu1_xref_(:,nbu1),2)); assert(length(tmp_bo1_)>=1);
for nlo1=1:length(tmp_bo1_);
nbo1 = tmp_bo1_(nlo1);
if strcmp(alleletype_arm1_{nbo1},'and'); tmp_tstx_arm1_and_(nbu3,:)=V_tstx_arm1_(nbo1,:); end;
if strcmp(alleletype_arm1_{nbo1},'xor'); tmp_tstx_arm1_xor_(nbu3,:)=V_tstx_arm1_(nbo1,:); end;
if strcmp(alleletype_arm1_{nbo1},'nor'); tmp_tstx_arm1_nor_(nbu3,:)=V_tstx_arm1_(nbo1,:); end;
end;%for nlo1=1:length(tmp_bo1_);
%%%%%%%%;
nbu2 = find(sum(bu2_by_bu3_xref_(:,nbu3),2)); assert(length(nbu2)==1);
tmp_bo2_ = find(sum(bo2_by_bu2_xref_(:,nbu2),2)); assert(length(tmp_bo2_)>=1);
for nlo2=1:length(tmp_bo2_);
nbo2 = tmp_bo2_(nlo2);
if strcmp(alleletype_arm2_{nbo2},'and'); tmp_tstx_arm2_and_(nbu3,:)=V_tstx_arm2_(nbo2,:); end;
if strcmp(alleletype_arm2_{nbo2},'xor'); tmp_tstx_arm2_xor_(nbu3,:)=V_tstx_arm2_(nbo2,:); end;
if strcmp(alleletype_arm2_{nbo2},'nor'); tmp_tstx_arm2_nor_(nbu3,:)=V_tstx_arm2_(nbo2,:); end;
end;%for nlo2=1:length(tmp_bo2_);
end;%for nbu3=1:length(bu3_);
disp(sprintf(' %% V_tstx_arm1_ - V_tstx_arm2_ (and) error %0.16f',norm(tmp_tstx_arm1_and_ - tmp_tstx_arm2_and_,'fro')));
disp(sprintf(' %% V_tstx_arm1_ - V_tstx_arm2_ (xor) error %0.16f',norm(tmp_tstx_arm1_xor_ - tmp_tstx_arm2_xor_,'fro')));
disp(sprintf(' %% V_tstx_arm1_ - V_tstx_arm2_ (nor) error %0.16f',norm(tmp_tstx_arm1_nor_ - tmp_tstx_arm2_nor_,'fro')));
figure(n_figure); n_figure = n_figure + 1;
subplot(2,3,1); plot(tmp_tstx_arm1_and_(:,1),tmp_tstx_arm2_and_(:,1),'.'); title('and 1');
subplot(2,3,2); plot(tmp_tstx_arm1_xor_(:,1),tmp_tstx_arm2_xor_(:,1),'.'); title('xor 1');
subplot(2,3,3); plot(tmp_tstx_arm1_nor_(:,1),tmp_tstx_arm2_nor_(:,1),'.'); title('nor 1');
subplot(2,3,4); plot(tmp_tstx_arm1_and_(:,2),tmp_tstx_arm2_and_(:,2),'.'); title('and 2');
subplot(2,3,5); plot(tmp_tstx_arm1_xor_(:,2),tmp_tstx_arm2_xor_(:,2),'.'); title('xor 2');
subplot(2,3,6); plot(tmp_tstx_arm1_nor_(:,2),tmp_tstx_arm2_nor_(:,2),'.'); title('nor 2');
end;%if flag_check;

tmp_dir = sprintf('%s/dir_dosage',lisa_arm1.dir_out_s0_pca);
if (~exist(tmp_dir,'dir')); disp(sprintf(' %% mkdir %s',tmp_dir)); mkdir(tmp_dir); end;
fname_stage_0 = sprintf('%s/%s_ni%d_dosage_stage_0.mat',tmp_dir,lisa_arm1.string_name_s0,niteration); 
fname_stage_1 = sprintf('%s/%s_ni%d_dosage_stage_1.mat',tmp_dir,lisa_arm1.string_name_s0,niteration); 
lisa_arm2_n_patient = lisa_arm2.n_patient;
lisa_arm1_bim_ = lisa_arm1.bim_;
lisa_arm1_A_p_ = lisa_arm1.A_p_;
lisa_arm2_A_p_ = lisa_arm2.A_p_;
mr_D_arm2_ = lisa_arm2.mr_D_;
mr_X_arm2_ = lisa_arm2.mr_X_;
study_name_arm2_ = lisa_arm2.study_name_;
fam_name_arm1_ = lisa_arm1.fam_name_;
fam_name_arm2_ = lisa_arm2.fam_name_;
tmp = 0;
save(fname_stage_0,...
'A_p_arm1_and_',...
'A_p_arm1_nor_',...
'A_p_arm1_xor_',...
'A_p_arm2_and_',...
'A_p_arm2_nor_',...
'A_p_arm2_xor_',...
'B_MLT',...
'DnV_tstx_arm1_',...
'DnV_tstx_arm2_',...
'POPLENGTH',...
'V_trnx_arm1_',...
'V_trnx_arm1_and_',...
'V_trnx_arm1_nor_',...
'V_trnx_arm1_xor_',...
'V_tstx_arm1_',...
'V_tstx_arm1_and_',...
'V_tstx_arm1_nor_',...
'V_tstx_arm1_xor_',...
'V_tstx_arm2_',...
'XnV_tstx_arm1_',...
'XnV_tstx_arm2_',...
'allele1_arm1_',...
'allele1_arm2_',...
'allele2_arm1_',...
'allele2_arm2_',...
'alleletype_arm1_',...
'alleletype_arm2_',...
'atfr_arm1_',...
'atfr_arm1_and_',...
'atfr_arm1_nor_',...
'atfr_arm1_xor_',...
'atfr_arm2_',...
'atfr_arm2_and_',...
'atfr_arm2_nor_',...
'atfr_arm2_xor_',...
'bo1_',...
'bo1_by_bu1_xref_',...
'bo1_to_bu1_',...
'bo2_',...
'bo2_by_bu2_xref_',...
'bo2_to_bu2_',...
'bu1_',...
'bu1_by_bu3_xref_',...
'bu1_to_bo1_',...
'bu2_',...
'bu2_by_bu3_xref_',...
'bu2_to_bo2_',...
'bu3_',...
'bu3_by_bu1_xref_',...
'bu3_by_bu2_xref_',...
'bu3_to_bu1_',...
'bu3_to_bu2_',...
'cl_num_arm1',...
'cl_num_arm2',...
'fam_name_arm1_',...
'fam_name_arm2_',...
'flag_check',...
'flag_dex_vs_lak',...
'flag_lisa',...
'flag_plot',...
'flag_rerun',...
'flag_reverse',...
'flag_stage',...
'flag_verbose',...
'fname_stage_0',...
'fname_stage_1',...
'gamma',...
'lisa_arm1_A_p_',...
'lisa_arm1_bim_',...
'lisa_arm2_A_p_',...
'lisa_arm2_n_patient',...
'mr_D_arm2_',...
'mr_D_ret_arm1_',...
'mr_D_rmv_arm1_',...
'mr_X_arm2_',...
'n_cov',...
'n_figure',...
'n_iteration',...
'n_iteration_stride',...
'n_maf',...
'n_mds',...
'n_scramble',...
'n_shuffle',...
'niteration',...
'niteration_',...
'pca_b_mlt',...
'pca_infix_arm1',...
'pca_proj_infix_arm1',...
'pca_proj_infix_arm2',...
'pca_rank',...
'pca_tolerance',...
'study_name_arm2_',...
'tmp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_stage==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_stage==1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

verbose=1; 
tmp_dir = sprintf('%s/dir_dosage',lisa_arm1.dir_out_s0_pca);
if (~exist(tmp_dir,'dir')); disp(sprintf(' %% mkdir %s',tmp_dir)); mkdir(tmp_dir); end;
fname_stage_0 = sprintf('%s/%s_ni%d_dosage_stage_0.mat',tmp_dir,lisa_arm1.string_name_s0,niteration); 
fname_stage_0_use = fname_stage_0;
fname_stage_1 = sprintf('%s/%s_ni%d_dosage_stage_1.mat',tmp_dir,lisa_arm1.string_name_s0,niteration); 
fname_stage_1_use = fname_stage_1;
load(fname_stage_0);
flag_lisa = flag_lisa_use; %<-- not necessarily the variable which was loaded. ;
fname_stage_0 = fname_stage_0_use;
fname_stage_1 = fname_stage_1_use;
flag_stage = 1;
if (flag_lisa==0); dir_dosage_data = sprintf('%s/dir_icuk_bicluster_dosage',dir_trunk); end;
if (flag_lisa==1); dir_dosage_data = sprintf('/home/grouther/dir_icuk_bicluster_dosage'); end;

dosage_imputed_DnV_arm2_ = zeros(lisa_arm2_n_patient,2);
dosage_imputed_XnV_arm2_ = zeros(lisa_arm2_n_patient,2);
dosage_imputed_DnV_arm2__ = cell(length(study_name_arm2_),1);
dosage_imputed_XnV_arm2__ = cell(length(study_name_arm2_),1);
dosage_rounded_DnV_arm2_ = zeros(lisa_arm2_n_patient,2);
dosage_rounded_XnV_arm2_ = zeros(lisa_arm2_n_patient,2);
dosage_rounded_DnV_arm2__ = cell(length(study_name_arm2_),1);
dosage_rounded_XnV_arm2__ = cell(length(study_name_arm2_),1);
dosage_genotyp_DnV_arm2_ = zeros(lisa_arm2_n_patient,2);
dosage_genotyp_XnV_arm2_ = zeros(lisa_arm2_n_patient,2);
dosage_genotyp_DnV_arm2__ = cell(length(study_name_arm2_),1);
dosage_genotyp_XnV_arm2__ = cell(length(study_name_arm2_),1);
mismatch_d1__ = cell(length(study_name_arm2_),1);

for nstudy_arm2=1:length(study_name_arm2_);
%%%%%%%%%%%%%%%%;
study_name_arm2 = study_name_arm2_{nstudy_arm2};
tmp_ij = strfind(study_name_arm2,'bip_');
study_name_arm2_prefix = study_name_arm2(tmp_ij+(4:7));
if flag_lisa==1;
fname_dosage    = sprintf('%s/%s_eur_sr-qc.hg19.ch.fl.out.dosage_%s_out_cdrop_bicl%d_ni%d.out.dosage',dir_dosage_data,study_name_arm2_prefix,lisa_arm1.string_name_s0,1+n_scramble,niteration);
tmp_command = sprintf('gunzip -f -k %s.gz',fname_dosage); disp(tmp_command); system(tmp_command);
end;%if flag_lisa==1;
if flag_lisa==0;
fname_dosage = sprintf('%s/%s_test1.cropped','/data/rangan/dir_bcc/dir_PGC_20190328/dir_icuk_bicluster_test',study_name_arm2_prefix);
end;%if flag_lisa==0;
fcheck(fname_dosage);
n_dosage_snp = wc_0(fname_dosage)-1;

fid = fopen(fname_dosage);
tmp_line_ = fgetl(fid);
tmp_dosage_patient_ = strsplit(tmp_line_);
assert(strcmp(tmp_dosage_patient_{1},'SNP'));
assert(strcmp(tmp_dosage_patient_{2},'A1'));
assert(strcmp(tmp_dosage_patient_{3},'A2'));
n_dosage_patient = floor((length(tmp_dosage_patient_)-3)/2);
dosage_fam_name_ = cell(n_dosage_patient,1);
for ndosage_patient=1:n_dosage_patient;
dosage_fam_name_{ndosage_patient} = sprintf('%s&%s',strtrim(tmp_dosage_patient_{2 + 2*ndosage_patient}),strtrim(tmp_dosage_patient_{3 + 2*ndosage_patient}));
end;%for ndosage_patient=1:n_dosage_patient;
dosage_bim_snp_ = cell(n_dosage_snp,1);
allele1_dosage_ = zeros(n_dosage_snp,1);
allele2_dosage_ = zeros(n_dosage_snp,1);
for ndosage_snp=1:n_dosage_snp;
if (mod(ndosage_snp,1000)==0); disp(sprintf(' %% ndosage_snp %d/%d',ndosage_snp,n_dosage_snp)); end;
tmp_line_ = fgetl(fid);
tmp_dosage_snp_ = textscan(tmp_line_,'%s',3);%tmp_dosage_snp_ = strsplit(tmp_line_);
tmp_dosage_snp_ = tmp_dosage_snp_{1};
flag_check=0;
if flag_check;
tmp_ij_ = find(strcmp(lisa_arm1_bim_{1},tmp_dosage_snp_{1}));
for nl=1:length(tmp_ij_);
disp(sprintf(' %% found snp %s (%s,%s) <-- bim %s (%c,%c): %s',tmp_dosage_snp_{1},tmp_dosage_snp_{2},tmp_dosage_snp_{3},lisa_arm1_bim_{1}{tmp_ij_(nl)},lisa_arm1_bim_{2}(tmp_ij_(nl)),lisa_arm1_bim_{3}(tmp_ij_(nl)),alleletype_arm1_{tmp_ij_(nl)}));
end;%for nl=1:length(tmp_ij_);
end;%if flag_check;
dosage_bim_snp_{ndosage_snp} = tmp_dosage_snp_{1};
allele1_dosage_(ndosage_snp) = tmp_dosage_snp_{2};
allele2_dosage_(ndosage_snp) = tmp_dosage_snp_{3};
end;%for ndosage_snp=1:n_dosage_snp;
fclose(fid);

%%%%%%%%;
dud_ = dosage_bim_snp_ ; %<-- we assume this is already unique. ;
%%%%%%%%;
[eu1_,eu1_to_dud_,eu1_to_bu1_] = intersect(dud_,bu1_,'stable');
dud_by_eu1_xref_ = sparse(eu1_to_dud_,1:length(eu1_),1,length(dud_),length(eu1_));
eu1_by_dud_xref_ = sparse(1:length(eu1_),eu1_to_dud_,1,length(eu1_),length(dud_));
bu1_by_eu1_xref_ = sparse(eu1_to_bu1_,1:length(eu1_),1,length(bu1_),length(eu1_));
%%%%%%%%;
flag_check=0;
if flag_check;
neu1=max(1,min(length(eu1_),floor(length(eu1_)*rand())));
disp(sprintf(' %% neu1 %.9d: eu1_{neu1} %s',neu1,eu1_{neu1}));
tmp_dud_ = find(sum(dud_by_eu1_xref_(:,neu1),2));
for nl1d=1:length(tmp_dud_);
ndud = tmp_dud_(nl1d);
disp(sprintf(' %% ndud %.9d: dud_{ndud} %s',ndud,dud_{ndud}));
assert(strcmp(dud_{ndud},eu1_{neu1}));
end;%for nl1d=1:length(tmp_dud_);
tmp_bu1_ = find(sum(bu1_by_eu1_xref_(:,neu1),2));
for nlu1=1:length(tmp_bu1_);
nbu1 = tmp_bu1_(nlu1);
disp(sprintf(' %% nbu1 %.9d: bu1_{nbu1} %s and (%0.6f %0.6f) , xor (%0.6f %0.6f) , nor (%0.6f %0.6f)',nbu1,bu1_{nbu1},V_trnx_arm1_and_(nbu1,:),V_trnx_arm1_xor_(nbu1,:),V_trnx_arm1_nor_(nbu1,:)));
assert(strcmp(bu1_{nbu1},eu1_{neu1}));
end;%for nlu1=1:length(tmp_bu1_);
tmp_bo1_ = find(sum(bo1_by_bu1_xref_(:,nbu1),2));
for nlo1=1:length(tmp_bo1_);
nbo1 = tmp_bo1_(nlo1);
disp(sprintf(' %% nbo1 %.9d: bo1_{nbo1} %s (%s) (%0.6f %0.6f)',nbo1,bo1_{nbo1},alleletype_arm1_{nbo1},V_trnx_arm1_(nbo1,:)));
assert(strcmp(bo1_{nbo1},eu1_{neu1}));
end;%for nlo1=1:length(tmp_bo1_);
end;%if flag_check;
%%%%%%%%;
[eu2_,eu2_to_dud_,eu2_to_bu2_] = intersect(dud_,bu2_,'stable');
dud_by_eu2_xref_ = sparse(eu2_to_dud_,1:length(eu2_),1,length(dud_),length(eu2_));
eu2_by_dud_xref_ = sparse(1:length(eu2_),eu2_to_dud_,1,length(eu2_),length(dud_));
bu2_by_eu2_xref_ = sparse(eu2_to_bu2_,1:length(eu2_),1,length(bu2_),length(eu2_));
%%%%%%%%;
flag_check=0;
if flag_check;
neu2=max(1,min(length(eu2_),floor(length(eu2_)*rand())));
disp(sprintf(' %% neu2 %.9d: eu2_{neu2} %s',neu2,eu2_{neu2}));
tmp_dud_ = find(sum(dud_by_eu2_xref_(:,neu2),2));
for nl1d=1:length(tmp_dud_);
ndud = tmp_dud_(nl1d);
disp(sprintf(' %% ndud %.9d: dud_{ndud} %s',ndud,dud_{ndud}));
assert(strcmp(dud_{ndud},eu2_{neu2}));
end;%for nl1d=1:length(tmp_dud_);
tmp_bu2_ = find(sum(bu2_by_eu2_xref_(:,neu2),2));
for nl12=1:length(tmp_bu2_);
nbu2 = tmp_bu2_(nl12);
disp(sprintf(' %% nbu2 %.9d: bu2_{nbu2} %s',nbu2,bu2_{nbu2}));
assert(strcmp(bu2_{nbu2},eu2_{neu2}));
end;%for nl12=1:length(tmp_bu2_);
tmp_bo2_ = find(sum(bo2_by_bu2_xref_(:,nbu2),2));
for nlo2=1:length(tmp_bo2_);
nbo2 = tmp_bo2_(nlo2);
disp(sprintf(' %% nbo2 %.9d: bo2_{nbo2} %s (%s)',nbo2,bo2_{nbo2},alleletype_arm2_{nbo2}));
assert(strcmp(bo2_{nbo2},eu2_{neu2}));
end;%for nlo2=1:length(tmp_bo2_);
end;%if flag_check;

% Ordering V_trnx_arm1_xxx_ and V_tstx_arm1_xxx_ to align with dosage file. ;
clear V_trnx_perm_and_ V_trnx_perm_xor_ V_trnx_perm_nor_;
clear V_tstx_perm_and_ V_tstx_perm_xor_ V_tstx_perm_nor_;
clear atfr_arm1_perm_and_ atfr_arm1_perm_xor_ atfr_arm1_perm_nor_;
clear atfr_arm2_perm_and_ atfr_arm2_perm_xor_ atfr_arm2_perm_nor_;
clear A_p_arm1_perm_and_ A_p_arm1_perm_xor_ A_p_arm1_perm_nor_;
clear A_p_arm2_perm_and_ A_p_arm2_perm_xor_ A_p_arm2_perm_nor_;
%%%%%%%%;
V_trnx_perm_and_ = zeros(length(dud_),2);
V_trnx_perm_xor_ = zeros(length(dud_),2);
V_trnx_perm_nor_ = zeros(length(dud_),2);
V_tstx_perm_and_ = zeros(length(dud_),2);
V_tstx_perm_xor_ = zeros(length(dud_),2);
V_tstx_perm_nor_ = zeros(length(dud_),2);
atfr_arm1_perm_and_ = zeros(length(dud_),1);
atfr_arm1_perm_xor_ = zeros(length(dud_),1);
atfr_arm1_perm_nor_ = zeros(length(dud_),1);
A_p_arm1_perm_and_ = zeros(length(dud_),1);
A_p_arm1_perm_xor_ = zeros(length(dud_),1);
A_p_arm1_perm_nor_ = zeros(length(dud_),1);
for neu1=1:length(eu1_);
if (mod(neu1,1000)==0); disp(sprintf(' %% neu1 %d/%d',neu1,length(eu1_))); end;
tmp_dud_ = find(sum(dud_by_eu1_xref_(:,neu1),2));
if (length(tmp_dud_)>1); disp(sprintf(' %% Warning! tmp_dud_ length %d in lisa_dosage_pca_ver9',length(tmp_dud_))); end;
ndud = tmp_dud_(1);
tmp_bu1_ = find(sum(bu1_by_eu1_xref_(:,neu1),2));
if (length(tmp_bu1_)>1); disp(sprintf(' %% Warning! tmp_bu1_ length %d in lisa_dosage_pca_ver9',length(tmp_bu1_))); end;
nbu1 = tmp_bu1_(1);
V_trnx_perm_and_(ndud,:) = V_trnx_arm1_and_(nbu1,:);
V_trnx_perm_xor_(ndud,:) = V_trnx_arm1_xor_(nbu1,:);
V_trnx_perm_nor_(ndud,:) = V_trnx_arm1_nor_(nbu1,:);
V_tstx_perm_and_(ndud,:) = V_tstx_arm1_and_(nbu1,:);
V_tstx_perm_xor_(ndud,:) = V_tstx_arm1_xor_(nbu1,:);
V_tstx_perm_nor_(ndud,:) = V_tstx_arm1_nor_(nbu1,:);
atfr_arm1_perm_and_(ndud) = atfr_arm1_and_(nbu1);
atfr_arm1_perm_xor_(ndud) = atfr_arm1_xor_(nbu1);
atfr_arm1_perm_nor_(ndud) = atfr_arm1_nor_(nbu1);
A_p_arm1_perm_and_(ndud) = A_p_arm1_and_(nbu1);
A_p_arm1_perm_xor_(ndud) = A_p_arm1_xor_(nbu1);
A_p_arm1_perm_nor_(ndud) = A_p_arm1_nor_(nbu1);
end;%for neu1=1:length(eu1_);
atfr_arm1_perm_and_alpha_ = atfr_arm1_perm_and_ - (1 - atfr_arm1_perm_and_); atfr_arm1_perm_and_delta_ = 1./( 4 .* atfr_arm1_perm_and_ .* (1 - atfr_arm1_perm_and_) );
atfr_arm1_perm_xor_alpha_ = atfr_arm1_perm_xor_ - (1 - atfr_arm1_perm_xor_); atfr_arm1_perm_xor_delta_ = 1./( 4 .* atfr_arm1_perm_xor_ .* (1 - atfr_arm1_perm_xor_) );
atfr_arm1_perm_nor_alpha_ = atfr_arm1_perm_nor_ - (1 - atfr_arm1_perm_nor_); atfr_arm1_perm_nor_delta_ = 1./( 4 .* atfr_arm1_perm_nor_ .* (1 - atfr_arm1_perm_nor_) );
A_p_arm1_perm_and_alpha_ = A_p_arm1_perm_and_ - (1 - A_p_arm1_perm_and_); A_p_arm1_perm_and_delta_ = 1./( 4 .* A_p_arm1_perm_and_ .* (1 - A_p_arm1_perm_and_) );
A_p_arm1_perm_xor_alpha_ = A_p_arm1_perm_xor_ - (1 - A_p_arm1_perm_xor_); A_p_arm1_perm_xor_delta_ = 1./( 4 .* A_p_arm1_perm_xor_ .* (1 - A_p_arm1_perm_xor_) );
A_p_arm1_perm_nor_alpha_ = A_p_arm1_perm_nor_ - (1 - A_p_arm1_perm_nor_); A_p_arm1_perm_nor_delta_ = 1./( 4 .* A_p_arm1_perm_nor_ .* (1 - A_p_arm1_perm_nor_) );
%%%%%%%%;
atfr_arm2_perm_and_ = zeros(length(dud_),1);
atfr_arm2_perm_xor_ = zeros(length(dud_),1);
atfr_arm2_perm_nor_ = zeros(length(dud_),1);
A_p_arm2_perm_and_ = zeros(length(dud_),1);
A_p_arm2_perm_xor_ = zeros(length(dud_),1);
A_p_arm2_perm_nor_ = zeros(length(dud_),1);
for neu2=1:length(eu2_);
if (mod(neu2,1000)==0); disp(sprintf(' %% neu2 %d/%d',neu2,length(eu2_))); end;
tmp_dud_ = find(sum(dud_by_eu2_xref_(:,neu2),2));
if (length(tmp_dud_)>1); disp(sprintf(' %% Warning! tmp_dud_ length %d in lisa_dosage_pca_ver9',length(tmp_dud_))); end;
ndud = tmp_dud_(1);
tmp_bu2_ = find(sum(bu2_by_eu2_xref_(:,neu2),2));
if (length(tmp_bu2_)>1); disp(sprintf(' %% Warning! tmp_bu2_ length %d in lisa_dosage_pca_ver9',length(tmp_bu2_))); end;
nbu2 = tmp_bu2_(1);
atfr_arm2_perm_and_(ndud) = atfr_arm2_and_(nbu2);
atfr_arm2_perm_xor_(ndud) = atfr_arm2_xor_(nbu2);
atfr_arm2_perm_nor_(ndud) = atfr_arm2_nor_(nbu2);
A_p_arm2_perm_and_(ndud) = A_p_arm2_and_(nbu2);
A_p_arm2_perm_xor_(ndud) = A_p_arm2_xor_(nbu2);
A_p_arm2_perm_nor_(ndud) = A_p_arm2_nor_(nbu2);
end;%for neu2=1:length(eu2_);
atfr_arm2_perm_and_alpha_ = atfr_arm2_perm_and_ - (1 - atfr_arm2_perm_and_); atfr_arm2_perm_and_delta_ = 1./( 4 .* atfr_arm2_perm_and_ .* (1 - atfr_arm2_perm_and_) );
atfr_arm2_perm_xor_alpha_ = atfr_arm2_perm_xor_ - (1 - atfr_arm2_perm_xor_); atfr_arm2_perm_xor_delta_ = 1./( 4 .* atfr_arm2_perm_xor_ .* (1 - atfr_arm2_perm_xor_) );
atfr_arm2_perm_nor_alpha_ = atfr_arm2_perm_nor_ - (1 - atfr_arm2_perm_nor_); atfr_arm2_perm_nor_delta_ = 1./( 4 .* atfr_arm2_perm_nor_ .* (1 - atfr_arm2_perm_nor_) );
A_p_arm2_perm_and_alpha_ = A_p_arm2_perm_and_ - (1 - A_p_arm2_perm_and_); A_p_arm2_perm_and_delta_ = 1./( 4 .* A_p_arm2_perm_and_ .* (1 - A_p_arm2_perm_and_) );
A_p_arm2_perm_xor_alpha_ = A_p_arm2_perm_xor_ - (1 - A_p_arm2_perm_xor_); A_p_arm2_perm_xor_delta_ = 1./( 4 .* A_p_arm2_perm_xor_ .* (1 - A_p_arm2_perm_xor_) );
A_p_arm2_perm_nor_alpha_ = A_p_arm2_perm_nor_ - (1 - A_p_arm2_perm_nor_); A_p_arm2_perm_nor_delta_ = 1./( 4 .* A_p_arm2_perm_nor_ .* (1 - A_p_arm2_perm_nor_) );
flag_check=0;
if flag_check;
subplot(2,3,1); plot(atfr_arm1_perm_and_,atfr_arm2_perm_and_,'.'); title('atfr and');
subplot(2,3,2); plot(atfr_arm1_perm_xor_,atfr_arm2_perm_xor_,'.'); title('atfr xor');
subplot(2,3,3); plot(atfr_arm1_perm_nor_,atfr_arm2_perm_nor_,'.'); title('atfr nor');
subplot(2,3,4); plot(A_p_arm1_perm_and_,A_p_arm2_perm_and_,'.'); title('A_p and');
subplot(2,3,5); plot(A_p_arm1_perm_xor_,A_p_arm2_perm_xor_,'.'); title('A_p xor');
subplot(2,3,6); plot(A_p_arm1_perm_nor_,A_p_arm2_perm_nor_,'.'); title('A_p nor');
end;%if flag_check;

fu1_ = fam_name_arm1_; %<-- We assume this is unique. ;
fu2_ = fam_name_arm2_; %<-- We assume this is unique. ;
gud_ = dosage_fam_name_; %<-- We assume this is unique. ;
%%%%%%%%;
[hu1_,hu1_to_gud_,hu1_to_fu1_] = intersect(gud_,fu1_,'stable');
gud_by_hu1_xref_ = sparse(hu1_to_gud_,1:length(hu1_),1,length(gud_),length(hu1_));
fu1_by_hu1_xref_ = sparse(hu1_to_fu1_,1:length(hu1_),1,length(fu1_),length(hu1_));
gud_to_fu1_ = zeros(length(gud_),1);
for nhu1=1:length(hu1_);
tmp_fu1_ = find(sum(fu1_by_hu1_xref_(:,nhu1),2));
if (length(tmp_fu1_)>1); disp(sprintf(' %% Warning! tmp_fu1_ length %d in lisa_dosage_pca_ver9',length(tmp_fu1_))); end;
nfu1 = tmp_fu1_(1);
tmp_gud_ = find(sum(gud_by_hu1_xref_(:,nhu1),2));
if (length(tmp_gud_)>1); disp(sprintf(' %% Warning! tmp_gud_ length %d in lisa_dosage_pca_ver9',length(tmp_gud_))); end;
ngud = tmp_gud_(1);
assert(strcmp(fu1_{nfu1},gud_{ngud}));
gud_to_fu1_(ngud) = nfu1;
end;%for nhu1=1:length(hu1_);
%%%%%%%%;
[hu2_,hu2_to_gud_,hu2_to_fu2_] = intersect(gud_,fu2_,'stable');
gud_by_hu2_xref_ = sparse(hu2_to_gud_,1:length(hu2_),1,length(gud_),length(hu2_));
fu2_by_hu2_xref_ = sparse(hu2_to_fu2_,1:length(hu2_),1,length(fu2_),length(hu2_));
gud_to_fu2_ = zeros(length(gud_),1);
for nhu2=1:length(hu2_);
tmp_fu2_ = find(sum(fu2_by_hu2_xref_(:,nhu2),2));
if (length(tmp_fu2_)>1); disp(sprintf(' %% Warning! tmp_fu2_ length %d in lisa_dosage_pca_ver9',length(tmp_fu2_))); end;
nfu2 = tmp_fu2_(1);
tmp_gud_ = find(sum(gud_by_hu2_xref_(:,nhu2),2));
if (length(tmp_gud_)>1); disp(sprintf(' %% Warning! tmp_gud_ length %d in lisa_dosage_pca_ver9',length(tmp_gud_))); end;
ngud = tmp_gud_(1);
assert(strcmp(fu2_{nfu2},gud_{ngud}));
gud_to_fu2_(ngud) = nfu2;
end;%for nhu2=1:length(hu2_);

%%%%%%%%;
char2num_ = zeros(84,1); 
char2num_(65) = 1; %<-- A ;
char2num_(67) = 2; %<-- C ;
char2num_(71) = 3; %<-- G ;
char2num_(84) = 4; %<-- T ;
%%%%%%%%;
allele_dosage__ = zeros(4,4);
for nl=1:n_dosage_snp;
allele_dosage__(char2num_(allele1_dosage_(nl)),char2num_(allele2_dosage_(nl))) = allele_dosage__(char2num_(allele1_dosage_(nl)),char2num_(allele2_dosage_(nl))) + 1;
end;%for nl=1:n_dosage_snp;
%%%%%%%%;
allele_arm1__ = zeros(4,4);
for nl=1:length(allele1_arm1_);
allele_arm1__(char2num_(allele1_arm1_(nl)),char2num_(allele2_arm1_(nl))) = allele_arm1__(char2num_(allele1_arm1_(nl)),char2num_(allele2_arm1_(nl))) +1;
end;%for nl=1:length(allele1_arm1_);
%%%%%%%%;
allele_arm2__ = zeros(4,4);
for nl=1:length(allele1_arm2_);
allele_arm2__(char2num_(allele1_arm2_(nl)),char2num_(allele2_arm2_(nl))) = allele_arm2__(char2num_(allele1_arm2_(nl)),char2num_(allele2_arm2_(nl))) +1;
end;%for nl=1:length(allele1_arm2_);
%%%%%%%%;
disp('allele_dosage__'); disp(allele_dosage__);
disp('allele_arm1__'); disp(allele_arm1__);
disp('allele_arm2__'); disp(allele_arm2__);

%%%%%%%%%%%%%%%%;
% We define different kinds of allelic mismatch errors: ;
%%%%%%%%%%%%%%%%;
error_code_null = 0;      % <-- No error; alleles match and minor-allele-frequency is commensurate. ;
error_code_synonym = 1;   % <-- Synonym; allele names must be interchanged with synonyms (i.e,. A<-->T and C<-->G). ;
error_code_dominance = 2; % <-- Dominance; the alleles must be switched (i.e., XY --> YX) so that the first allele listed corresponds to the minor (i.e., less frequent) allele. ;
error_code_both = 3;      % <-- Both; the allele names must be interchanged *and* the order must be switched. ;
error_code_ambiguous = 8; % <-- Ambiguity; I cannot determine the correct alignment. ;
error_code_missing = 9;   % <-- Missing; this snp is not in the overlap between the two sets. ;
%%%%%%%%%%%%%%%%;
flag_check=0;
if flag_check;
%%%%%%%%;
mismatch_d1_ = 9*ones(n_dosage_snp,1);
for ndosage_snp=1:n_dosage_snp;
if (mod(ndosage_snp,1000)==0); disp(sprintf(' %% ndosage_snp %d/%d',ndosage_snp,n_dosage_snp)); end;
ndud = ndosage_snp;
neu1 = find(sum(eu1_by_dud_xref_(:,ndud),2)); assert(length(neu1)<=1);
if (length(neu1)==1);
nbu1 = find(sum(bu1_by_eu1_xref_(:,neu1),2)); assert(length(nbu1)==1);
tmp_bo1_ = find(sum(bo1_by_bu1_xref_(:,nbu1),2)); assert(length(tmp_bo1_)>=1);
for nlo1=1:length(tmp_bo1_);
nbo1 = tmp_bo1_(nlo1);
assert(strcmp(bo1_{nbo1},dud_{ndud}));
assert(strcmp(bo1_{nbo1},eu1_{neu1}));
assert(strcmp(bo1_{nbo1},bu1_{nbu1}));
end;%for nlo1=1:length(tmp_bo1_);
nbo1 = tmp_bo1_(1);
mismatch_d1_(ndud) = lisa_allele_mismatch_ver0(allele1_dosage_(ndud),allele2_dosage_(ndud),allele1_arm1_(nbo1),allele2_arm1_(nbo1));
end;%if (length(neu1)==1);
end;%for ndosage_snp=1:n_dosage_snp;
%%%%%%%%;
mismatch_d2_ = 9*ones(n_dosage_snp,1);
for ndosage_snp=1:n_dosage_snp;
if (mod(ndosage_snp,1000)==0); disp(sprintf(' %% ndosage_snp %d/%d',ndosage_snp,n_dosage_snp)); end;
ndud = ndosage_snp;
neu2 = find(sum(eu2_by_dud_xref_(:,ndud),2)); assert(length(neu2)<=1);
if (length(neu2)==1);
nbu2 = find(sum(bu2_by_eu2_xref_(:,neu2),2)); assert(length(nbu2)==1);
tmp_bo2_ = find(sum(bo2_by_bu2_xref_(:,nbu2),2)); assert(length(tmp_bo2_)>=1);
for nlo2=1:length(tmp_bo2_);
nbo2 = tmp_bo2_(nlo2);
assert(strcmp(bo2_{nbo2},dud_{ndud}));
assert(strcmp(bo2_{nbo2},eu2_{neu2}));
assert(strcmp(bo2_{nbo2},bu2_{nbu2}));
end;%for nlo2=1:length(tmp_bo2_);
nbo2 = tmp_bo2_(1);
mismatch_d2_(ndud) = lisa_allele_mismatch_ver0(allele1_dosage_(ndud),allele2_dosage_(ndud),allele1_arm2_(nbo2),allele2_arm2_(nbo2));
end;%if (length(neu2)==1);
end;%for ndosage_snp=1:n_dosage_snp;
%%%%%%%%;
end%if flag_check;

%%%%%%%%%%%%%%%%;
% Now step through the dosage file line by line. ;
% updating the inner-products as we go. ;
%%%%%%%%%%%%%%%%;
ent_cutoff = 0.03; 
flag_V_trnx = 1; %<-- use V from training set (all snps from arm1). ;
flag_V_tstx = 0; %<-- use V from testing set (only snps in intersection of arm1 and arm2). ;
flag_normalize_atfr = 0; %<-- use atfr to normalize. ;
flag_normalize_A_p = 1; %<-- use A_p to normalize. ;
assert(~(flag_normalize_atfr & flag_normalize_A_p));
flag_normalize_arm1 = 1; %<-- use arm1 to normalize. ;
flag_normalize_arm2 = 0; %<-- use arm2 to normalize. ;
assert(~(flag_normalize_arm1 & flag_normalize_arm2));
flag_match_arm1 = 1; %<-- swap dosage data to match arm1. ;
flag_match_arm2 = 0; %<-- swap dosage data to match arm2. ;
assert(~(flag_match_arm1 & flag_match_arm2));
tmp_dosage_imputed_AnV_ = zeros(n_dosage_patient,2);
tmp_dosage_rounded_AnV_ = zeros(n_dosage_patient,2);
tmp_dosage_genotyp_AnV_ = zeros(n_dosage_patient,2);
mismatch_d1_imputed__{nstudy_arm2} = 9*ones(n_dosage_snp,1);
mismatch_d1_rounded__{nstudy_arm2} = 9*ones(n_dosage_snp,1);
mismatch_d1_genotyp__{nstudy_arm2} = 9*ones(n_dosage_snp,1);
mismatch_d2_imputed__{nstudy_arm2} = 9*ones(n_dosage_snp,1);
mismatch_d2_rounded__{nstudy_arm2} = 9*ones(n_dosage_snp,1);
mismatch_d2_genotyp__{nstudy_arm2} = 9*ones(n_dosage_snp,1);
flag_bo1_all_sum = 0; flag_bo2_and_sum = 0; flag_bo2_xor_sum = 0; flag_bo2_nor_sum = 0; n_z_a = [0,0];
%%%%%%%%;
fid = fopen(fname_dosage);
tmp_line_ = fgetl(fid); %<-- skip header line. ;
for ndosage_snp=1:n_dosage_snp;
if (mod(ndosage_snp,1000)==0); disp(sprintf(' %% ndosage_snp %d/%d',ndosage_snp,n_dosage_snp)); end;
tmp_line_ = fgetl(fid);
tmp_dosage_snp_ = textscan(tmp_line_,'%s',3);
tmp_dosage_al1_ = tmp_dosage_snp_{1}{2};
tmp_dosage_al2_ = tmp_dosage_snp_{1}{3};
ndud = ndosage_snp;
assert(strcmp(tmp_dosage_snp_{1}{1},dud_{ndud}));
assert(tmp_dosage_al1_==allele1_dosage_(ndud));
assert(tmp_dosage_al2_==allele2_dosage_(ndud));
%%%%%%%%;
% check arm1 for snp;
%%%%%%%%;
neu1 = find(sum(eu1_by_dud_xref_(:,ndud),2)); assert(length(neu1)<=1);
if (length(neu1)==1);
flag_bo1_all_sum = flag_bo1_all_sum + 1;
nbu1 = find(sum(bu1_by_eu1_xref_(:,neu1),2)); assert(length(nbu1)==1);
tmp_bo1_ = find(sum(bo1_by_bu1_xref_(:,nbu1),2)); assert(length(tmp_bo1_)>=1);
nbo1 = tmp_bo1_(1);
assert(strcmp(bo1_{nbo1},dud_{ndud}));
assert(strcmp(bo1_{nbo1},eu1_{neu1}));
assert(strcmp(bo1_{nbo1},bu1_{nbu1}));
tmp_length = sum(cellfun('length',tmp_dosage_snp_{1})) + length(tmp_dosage_snp_{1});
tmp_dosage_imputed_val_ = sscanf(tmp_line_(tmp_length:end),'%f');
assert(length(tmp_dosage_imputed_val_)/2==n_dosage_patient);
tmp_dosage_rounded_val_ = max(0,min(2,round(tmp_dosage_imputed_val_)));
tmp_dosage_imputed_val_nor_ = tmp_dosage_imputed_val_(1:2:end-1);
tmp_dosage_imputed_val_xor_ = tmp_dosage_imputed_val_(2:2:end-0);
tmp_dosage_imputed_val_and_ = ones(n_dosage_patient,1) - tmp_dosage_imputed_val_nor_ - tmp_dosage_imputed_val_xor_;
tmp_dosage_rounded_val_nor_ = tmp_dosage_rounded_val_(1:2:end-1);
tmp_dosage_rounded_val_xor_ = tmp_dosage_rounded_val_(2:2:end-0);
tmp_dosage_rounded_val_and_ = ones(n_dosage_patient,1) - tmp_dosage_rounded_val_nor_ - tmp_dosage_rounded_val_xor_;
%%%%%%%%;
% check arm2 for snp;
%%%%%%%%;
tmp_error_code_d2 = error_code_missing;
flag_bo2_and = 0; flag_bo2_xor = 0; flag_bo2_nor = 0;
nbu3 = find(sum(bu3_by_bu1_xref_(:,nbu1),2)); assert(length(nbu3)<=1);
if (~isempty(nbu3));
nbu2 = find(sum(bu2_by_bu3_xref_(:,nbu3),2)); assert(length(nbu2)==1);
tmp_bo2_ = find(sum(bo2_by_bu2_xref_(:,nbu2),2)); assert(length(tmp_bo2_)>=1);
%%%%%%%%;
% check for mismatch in d2. ;
%%%%%%%%;
nbo2 = tmp_bo2_(1);
assert(strcmp(bo2_{nbo2},dud_{ndud}));
assert(strcmp(bo2_{nbo2},bu1_{nbu1}));
tmp_error_code_d2 = lisa_allele_mismatch_ver0(...
  allele1_dosage_(ndud),allele2_dosage_(ndud),allele1_arm2_(nbo2),allele2_arm2_(nbo2)...
  ,sum(abs(tmp_dosage_imputed_val_and_)),sum(abs(tmp_dosage_imputed_val_xor_)),sum(abs(tmp_dosage_imputed_val_nor_))...
  ,atfr_arm2_perm_and_(ndud),atfr_arm2_perm_xor_(ndud),atfr_arm2_perm_nor_(ndud)...
  ,ent_cutoff...
  );
mismatch_d2_genotyp__{nstudy_arm2}(ndud) = tmp_error_code_d2;
mismatch_d2_imputed__{nstudy_arm2}(ndud) = tmp_error_code_d2;
%%%%%%%%;
for nlo2=1:length(tmp_bo2_);
nbo2 = tmp_bo2_(nlo2);
if strcmp(alleletype_arm2_{nbo2},'and'); flag_bo2_and = 1; flag_bo2_and_sum = flag_bo2_and_sum + 1; end;
if strcmp(alleletype_arm2_{nbo2},'xor'); flag_bo2_xor = 1; flag_bo2_xor_sum = flag_bo2_xor_sum + 1; end;
if strcmp(alleletype_arm2_{nbo2},'nor'); flag_bo2_nor = 1; flag_bo2_nor_sum = flag_bo2_nor_sum + 1; end;
end;%for nlo2=1:length(tmp_bo2_);
end;%if (~isempty(nbu3));
tmp_dosage_genotyp_val_nor_ = flag_bo2_nor * tmp_dosage_rounded_val_nor_;
tmp_dosage_genotyp_val_xor_ = flag_bo2_xor * tmp_dosage_rounded_val_xor_;
tmp_dosage_genotyp_val_and_ = flag_bo2_and * tmp_dosage_rounded_val_and_;
%%%%%%%%;
% check for mismatch in d1. ;
%%%%%%%%;
tmp_error_code_d1 = lisa_allele_mismatch_ver0(...
  allele1_dosage_(ndud),allele2_dosage_(ndud),allele1_arm1_(nbo1),allele2_arm1_(nbo1)...
  ,sum(abs(tmp_dosage_imputed_val_and_)),sum(abs(tmp_dosage_imputed_val_xor_)),sum(abs(tmp_dosage_imputed_val_nor_))...
  ,atfr_arm1_perm_and_(ndud),atfr_arm1_perm_xor_(ndud),atfr_arm1_perm_nor_(ndud)...
  ,ent_cutoff...
  );
mismatch_d1_genotyp__{nstudy_arm2}(ndud) = tmp_error_code_d1;
if (flag_bo2_and || flag_bo2_xor || flag_bo2_nor); mismatch_d1_imputed__{nstudy_arm2}(ndud) = tmp_error_code_d1; end;
%%%%%%%%;
% match to either arm1 or arm2 if necessary. ;
%%%%%%%%;
if ( flag_match_arm1 & ( (tmp_error_code_d1==error_code_dominance) | (tmp_error_code_d1==error_code_both) ) )...
 | ( flag_match_arm2 & ( (tmp_error_code_d2==error_code_dominance) | (tmp_error_code_d2==error_code_both) ) );
tmp_ = tmp_dosage_imputed_val_nor_ ; tmp_dosage_imputed_val_nor_ = tmp_dosage_imputed_val_and_; tmp_dosage_imputed_val_and_ = tmp_ ; %<-- switching dominant allele ;
tmp_ = tmp_dosage_rounded_val_nor_ ; tmp_dosage_rounded_val_nor_ = tmp_dosage_rounded_val_and_; tmp_dosage_rounded_val_and_ = tmp_ ; %<-- switching dominant allele ;
tmp_ = tmp_dosage_genotyp_val_nor_ ; tmp_dosage_genotyp_val_nor_ = tmp_dosage_genotyp_val_and_; tmp_dosage_genotyp_val_and_ = tmp_ ; %<-- switching dominant allele ;
end;%if match. ;
%%%%%%%%;
% set to zero if necessary. ;
%%%%%%%%;
if ( flag_match_arm1 & ( 0*(tmp_error_code_d1==error_code_ambiguous) | 1*(tmp_error_code_d1==error_code_missing) ) ) ...
 | ( flag_match_arm2 & ( 0*(tmp_error_code_d2==error_code_ambiguous) | 1*(tmp_error_code_d2==error_code_missing) ) ) ;
tmp_dosage_imputed_val_and_normalized_ = zeros(length(gud_),1);
tmp_dosage_imputed_val_xor_normalized_ = zeros(length(gud_),1);
tmp_dosage_imputed_val_nor_normalized_ = zeros(length(gud_),1);
tmp_dosage_rounded_val_and_normalized_ = zeros(length(gud_),1);
tmp_dosage_rounded_val_xor_normalized_ = zeros(length(gud_),1);
tmp_dosage_rounded_val_nor_normalized_ = zeros(length(gud_),1);
tmp_dosage_genotyp_val_and_normalized_ = zeros(length(gud_),1);
tmp_dosage_genotyp_val_xor_normalized_ = zeros(length(gud_),1);
tmp_dosage_genotyp_val_nor_normalized_ = zeros(length(gud_),1);
%%%%%%%%;
% otherwise if nonzero: ;
%%%%%%%%;
else;
if flag_normalize_atfr & flag_normalize_arm1
assert(isfinite(atfr_arm1_perm_and_delta_(ndud)) & isfinite(atfr_arm1_perm_and_alpha_(ndud)));
tmp_dosage_imputed_val_and_normalized_ = (2*tmp_dosage_imputed_val_and_ - 1 - atfr_arm1_perm_and_alpha_(ndud)) * sqrt(atfr_arm1_perm_and_delta_(ndud));
tmp_dosage_imputed_val_xor_normalized_ = (2*tmp_dosage_imputed_val_xor_ - 1 - atfr_arm1_perm_xor_alpha_(ndud)) * sqrt(atfr_arm1_perm_xor_delta_(ndud));
tmp_dosage_imputed_val_nor_normalized_ = (2*tmp_dosage_imputed_val_nor_ - 1 - atfr_arm1_perm_nor_alpha_(ndud)) * sqrt(atfr_arm1_perm_nor_delta_(ndud));
tmp_dosage_rounded_val_and_normalized_ = (2*tmp_dosage_rounded_val_and_ - 1 - atfr_arm1_perm_and_alpha_(ndud)) * sqrt(atfr_arm1_perm_and_delta_(ndud));
tmp_dosage_rounded_val_xor_normalized_ = (2*tmp_dosage_rounded_val_xor_ - 1 - atfr_arm1_perm_xor_alpha_(ndud)) * sqrt(atfr_arm1_perm_xor_delta_(ndud));
tmp_dosage_rounded_val_nor_normalized_ = (2*tmp_dosage_rounded_val_nor_ - 1 - atfr_arm1_perm_nor_alpha_(ndud)) * sqrt(atfr_arm1_perm_nor_delta_(ndud));
tmp_dosage_genotyp_val_and_normalized_ = (2*tmp_dosage_genotyp_val_and_ - 1 - atfr_arm1_perm_and_alpha_(ndud)) * sqrt(atfr_arm1_perm_and_delta_(ndud));
tmp_dosage_genotyp_val_xor_normalized_ = (2*tmp_dosage_genotyp_val_xor_ - 1 - atfr_arm1_perm_xor_alpha_(ndud)) * sqrt(atfr_arm1_perm_xor_delta_(ndud));
tmp_dosage_genotyp_val_nor_normalized_ = (2*tmp_dosage_genotyp_val_nor_ - 1 - atfr_arm1_perm_nor_alpha_(ndud)) * sqrt(atfr_arm1_perm_nor_delta_(ndud));
end;%if flag_normalize_atfr & flag_normalize_arm1
if flag_normalize_A_p & flag_normalize_arm1
assert(isfinite(A_p_arm1_perm_and_delta_(ndud)) & isfinite(A_p_arm1_perm_and_alpha_(ndud)));
tmp_dosage_imputed_val_and_normalized_ = (2*tmp_dosage_imputed_val_and_ - 1 - A_p_arm1_perm_and_alpha_(ndud)) * sqrt(A_p_arm1_perm_and_delta_(ndud));
tmp_dosage_imputed_val_xor_normalized_ = (2*tmp_dosage_imputed_val_xor_ - 1 - A_p_arm1_perm_xor_alpha_(ndud)) * sqrt(A_p_arm1_perm_xor_delta_(ndud));
tmp_dosage_imputed_val_nor_normalized_ = (2*tmp_dosage_imputed_val_nor_ - 1 - A_p_arm1_perm_nor_alpha_(ndud)) * sqrt(A_p_arm1_perm_nor_delta_(ndud));
tmp_dosage_rounded_val_and_normalized_ = (2*tmp_dosage_rounded_val_and_ - 1 - A_p_arm1_perm_and_alpha_(ndud)) * sqrt(A_p_arm1_perm_and_delta_(ndud));
tmp_dosage_rounded_val_xor_normalized_ = (2*tmp_dosage_rounded_val_xor_ - 1 - A_p_arm1_perm_xor_alpha_(ndud)) * sqrt(A_p_arm1_perm_xor_delta_(ndud));
tmp_dosage_rounded_val_nor_normalized_ = (2*tmp_dosage_rounded_val_nor_ - 1 - A_p_arm1_perm_nor_alpha_(ndud)) * sqrt(A_p_arm1_perm_nor_delta_(ndud));
tmp_dosage_genotyp_val_and_normalized_ = (2*tmp_dosage_genotyp_val_and_ - 1 - A_p_arm1_perm_and_alpha_(ndud)) * sqrt(A_p_arm1_perm_and_delta_(ndud));
tmp_dosage_genotyp_val_xor_normalized_ = (2*tmp_dosage_genotyp_val_xor_ - 1 - A_p_arm1_perm_xor_alpha_(ndud)) * sqrt(A_p_arm1_perm_xor_delta_(ndud));
tmp_dosage_genotyp_val_nor_normalized_ = (2*tmp_dosage_genotyp_val_nor_ - 1 - A_p_arm1_perm_nor_alpha_(ndud)) * sqrt(A_p_arm1_perm_nor_delta_(ndud));
end;%if flag_normalize_A_p & flag_normalize_arm1
%%%%%%%%;
if flag_normalize_atfr & flag_normalize_arm2;
assert(isfinite(atfr_arm2_perm_and_delta_(ndud)) & isfinite(atfr_arm2_perm_and_alpha_(ndud)));
tmp_dosage_imputed_val_and_normalized_ = (2*tmp_dosage_imputed_val_and_ - 1 - atfr_arm2_perm_and_alpha_(ndud)) * sqrt(atfr_arm2_perm_and_delta_(ndud));
tmp_dosage_imputed_val_xor_normalized_ = (2*tmp_dosage_imputed_val_xor_ - 1 - atfr_arm2_perm_xor_alpha_(ndud)) * sqrt(atfr_arm2_perm_xor_delta_(ndud));
tmp_dosage_imputed_val_nor_normalized_ = (2*tmp_dosage_imputed_val_nor_ - 1 - atfr_arm2_perm_nor_alpha_(ndud)) * sqrt(atfr_arm2_perm_nor_delta_(ndud));
tmp_dosage_rounded_val_and_normalized_ = (2*tmp_dosage_rounded_val_and_ - 1 - atfr_arm2_perm_and_alpha_(ndud)) * sqrt(atfr_arm2_perm_and_delta_(ndud));
tmp_dosage_rounded_val_xor_normalized_ = (2*tmp_dosage_rounded_val_xor_ - 1 - atfr_arm2_perm_xor_alpha_(ndud)) * sqrt(atfr_arm2_perm_xor_delta_(ndud));
tmp_dosage_rounded_val_nor_normalized_ = (2*tmp_dosage_rounded_val_nor_ - 1 - atfr_arm2_perm_nor_alpha_(ndud)) * sqrt(atfr_arm2_perm_nor_delta_(ndud));
tmp_dosage_genotyp_val_and_normalized_ = (2*tmp_dosage_genotyp_val_and_ - 1 - atfr_arm2_perm_and_alpha_(ndud)) * sqrt(atfr_arm2_perm_and_delta_(ndud));
tmp_dosage_genotyp_val_xor_normalized_ = (2*tmp_dosage_genotyp_val_xor_ - 1 - atfr_arm2_perm_xor_alpha_(ndud)) * sqrt(atfr_arm2_perm_xor_delta_(ndud));
tmp_dosage_genotyp_val_nor_normalized_ = (2*tmp_dosage_genotyp_val_nor_ - 1 - atfr_arm2_perm_nor_alpha_(ndud)) * sqrt(atfr_arm2_perm_nor_delta_(ndud));
end;%if flag_normalize_atfr & flag_normalize_arm2
if flag_normalize_A_p & flag_normalize_arm2;
assert(isfinite(A_p_arm2_perm_and_delta_(ndud)) & isfinite(A_p_arm2_perm_and_alpha_(ndud)));
tmp_dosage_imputed_val_and_normalized_ = (2*tmp_dosage_imputed_val_and_ - 1 - A_p_arm2_perm_and_alpha_(ndud)) * sqrt(A_p_arm2_perm_and_delta_(ndud));
tmp_dosage_imputed_val_xor_normalized_ = (2*tmp_dosage_imputed_val_xor_ - 1 - A_p_arm2_perm_xor_alpha_(ndud)) * sqrt(A_p_arm2_perm_xor_delta_(ndud));
tmp_dosage_imputed_val_nor_normalized_ = (2*tmp_dosage_imputed_val_nor_ - 1 - A_p_arm2_perm_nor_alpha_(ndud)) * sqrt(A_p_arm2_perm_nor_delta_(ndud));
tmp_dosage_rounded_val_and_normalized_ = (2*tmp_dosage_rounded_val_and_ - 1 - A_p_arm2_perm_and_alpha_(ndud)) * sqrt(A_p_arm2_perm_and_delta_(ndud));
tmp_dosage_rounded_val_xor_normalized_ = (2*tmp_dosage_rounded_val_xor_ - 1 - A_p_arm2_perm_xor_alpha_(ndud)) * sqrt(A_p_arm2_perm_xor_delta_(ndud));
tmp_dosage_rounded_val_nor_normalized_ = (2*tmp_dosage_rounded_val_nor_ - 1 - A_p_arm2_perm_nor_alpha_(ndud)) * sqrt(A_p_arm2_perm_nor_delta_(ndud));
tmp_dosage_genotyp_val_and_normalized_ = (2*tmp_dosage_genotyp_val_and_ - 1 - A_p_arm2_perm_and_alpha_(ndud)) * sqrt(A_p_arm2_perm_and_delta_(ndud));
tmp_dosage_genotyp_val_xor_normalized_ = (2*tmp_dosage_genotyp_val_xor_ - 1 - A_p_arm2_perm_xor_alpha_(ndud)) * sqrt(A_p_arm2_perm_xor_delta_(ndud));
tmp_dosage_genotyp_val_nor_normalized_ = (2*tmp_dosage_genotyp_val_nor_ - 1 - A_p_arm2_perm_nor_alpha_(ndud)) * sqrt(A_p_arm2_perm_nor_delta_(ndud));
end;%if flag_normalize_A_p & flag_normalize_arm2
end;%if set to zero. ;
%%%%%%%%%%%%%%%%;
% check against b16 file. ;
%%%%%%%%%%%%%%%%;
flag_check=0;
if flag_check;
%%%%%%%%;
tmp_bo2_ = find(strcmp(bo2_,dud_{ndud}));
if (length(tmp_bo2_)>0);
%lisa_dosage_pca_helper_0;
fname_A_n_arm2 = sprintf('%s/%s_A_full_n.b16',lisa_arm2.dir__in,lisa_arm2.string_prefix); fcheck(fname_A_n_arm2);
[tmp_nbins,tmp_nrows,tmp_ncols] = tutorial_binary_getsize(fname_A_n_arm2);
%%%%%%%%;
tmp_bo2_and = find(strcmp(bo2_,dud_{ndud}) & strcmp(alleletype_arm2_,'and')); if (length(tmp_bo2_and)>1); disp(sprintf(' %% Warning! tmp_bo2_and length %d in lisa_dosage_pca_ver9.m',length(tmp_bo2_and))); end;
if (length(tmp_bo2_and)>0);
tmp_arm2_and_ = tutorial_binary_uncompress(fname_A_n_arm2,1:tmp_nrows,min(tmp_ncols,tmp_bo2_and));
tmp_arm2_perm_and_ = zeros(length(gud_),1);
for nhu2=1:length(hu2_);
tmp_fu2_ = find(sum(fu2_by_hu2_xref_(:,nhu2),2)); if (length(tmp_fu2_)>1); disp(sprintf(' %% Warning! tmp_fu2_ length %d in lisa_dosage_pca_ver9',length(tmp_fu2_))); end; nfu2 = tmp_fu2_(1);
tmp_gud_ = find(sum(gud_by_hu2_xref_(:,nhu2),2)); if (length(tmp_gud_)>1); disp(sprintf(' %% Warning! tmp_gud_ length %d in lisa_dosage_pca_ver9',length(tmp_gud_))); end; ngud = tmp_gud_(1);
assert(strcmp(fu2_{nfu2},gud_{ngud}));
tmp_arm2_perm_and_(ngud) = tmp_arm2_and_(nfu2);
end;%for nhu2=1:length(hu2_);
tmp_str_and_and = sprintf('arm2_and vs dosage_imputed_val_and %.10d',sum(abs((tmp_arm2_perm_and_>0) - tmp_dosage_imputed_val_and_)));
tmp_str_and_xor = sprintf('arm2_and vs dosage_imputed_val_xor %.10d',sum(abs((tmp_arm2_perm_and_>0) - tmp_dosage_imputed_val_xor_)));
tmp_str_and_nor = sprintf('arm2_and vs dosage_imputed_val_nor %.10d',sum(abs((tmp_arm2_perm_and_>0) - tmp_dosage_imputed_val_nor_)));
disp(sprintf(' %% %% ndud %.10d tmp_bo2_and %.10d %s (%s) , %s , %s , %s',ndud,tmp_bo2_and,dud_{ndud},alleletype_arm2_{tmp_bo2_and},tmp_str_and_and,tmp_str_and_xor,tmp_str_and_nor));
end;%if (length(tmp_bo2_and)>0);
%%%%%%%%;
tmp_bo2_xor = find(strcmp(bo2_,dud_{ndud}) & strcmp(alleletype_arm2_,'xor')); if (length(tmp_bo2_xor)>1); disp(sprintf(' %% Warning! tmp_bo2_xor length %d in lisa_dosage_pca_ver9.m',length(tmp_bo2_xor))); end;
if (length(tmp_bo2_xor)>0);
tmp_arm2_xor_ = tutorial_binary_uncompress(fname_A_n_arm2,1:tmp_nrows,min(tmp_ncols,tmp_bo2_xor));
tmp_arm2_perm_xor_ = zeros(length(gud_),1);
for nhu2=1:length(hu2_);
tmp_fu2_ = find(sum(fu2_by_hu2_xref_(:,nhu2),2)); if (length(tmp_fu2_)>1); disp(sprintf(' %% Warning! tmp_fu2_ length %d in lisa_dosage_pca_ver9',length(tmp_fu2_))); end; nfu2 = tmp_fu2_(1);
tmp_gud_ = find(sum(gud_by_hu2_xref_(:,nhu2),2)); if (length(tmp_gud_)>1); disp(sprintf(' %% Warning! tmp_gud_ length %d in lisa_dosage_pca_ver9',length(tmp_gud_))); end; ngud = tmp_gud_(1);
assert(strcmp(fu2_{nfu2},gud_{ngud}));
tmp_arm2_perm_xor_(ngud) = tmp_arm2_xor_(nfu2);
end;%for nhu2=1:length(hu2_);
tmp_str_xor_and = sprintf('arm2_xor vs dosage_imputed_val_and %.10d',sum(abs((tmp_arm2_perm_xor_>0) - tmp_dosage_imputed_val_and_)));
tmp_str_xor_xor = sprintf('arm2_xor vs dosage_imputed_val_xor %.10d',sum(abs((tmp_arm2_perm_xor_>0) - tmp_dosage_imputed_val_xor_)));
tmp_str_xor_nor = sprintf('arm2_xor vs dosage_imputed_val_nor %.10d',sum(abs((tmp_arm2_perm_xor_>0) - tmp_dosage_imputed_val_nor_)));
disp(sprintf(' %% %% ndud %.10d tmp_bo2_xor %.10d %s (%s) , %s , %s , %s',ndud,tmp_bo2_xor,dud_{ndud},alleletype_arm2_{tmp_bo2_xor},tmp_str_xor_and,tmp_str_xor_xor,tmp_str_xor_nor));
end;%if (length(tmp_bo2_xor)>0);
%%%%%%%%;
tmp_bo2_nor = find(strcmp(bo2_,dud_{ndud}) & strcmp(alleletype_arm2_,'nor')); if (length(tmp_bo2_nor)>1); disp(sprintf(' %% Warning! tmp_bo2_nor length %d in lisa_dosage_pca_ver9.m',length(tmp_bo2_nor))); end;
if (length(tmp_bo2_nor)>0);
tmp_arm2_nor_ = tutorial_binary_uncompress(fname_A_n_arm2,1:tmp_nrows,min(tmp_ncols,tmp_bo2_nor));
tmp_arm2_perm_nor_ = zeros(length(gud_),1);
for nhu2=1:length(hu2_);
tmp_fu2_ = find(sum(fu2_by_hu2_xref_(:,nhu2),2)); if (length(tmp_fu2_)>1); disp(sprintf(' %% Warning! tmp_fu2_ length %d in lisa_dosage_pca_ver9',length(tmp_fu2_))); end; nfu2 = tmp_fu2_(1);
tmp_gud_ = find(sum(gud_by_hu2_xref_(:,nhu2),2)); if (length(tmp_gud_)>1); disp(sprintf(' %% Warning! tmp_gud_ length %d in lisa_dosage_pca_ver9',length(tmp_gud_))); end; ngud = tmp_gud_(1);
assert(strcmp(fu2_{nfu2},gud_{ngud}));
tmp_arm2_perm_nor_(ngud) = tmp_arm2_nor_(nfu2);
end;%for nhu2=1:length(hu2_);
tmp_str_nor_and = sprintf('arm2_nor vs dosage_imputed_val_and %.10d',sum(abs((tmp_arm2_perm_nor_>0) - tmp_dosage_imputed_val_and_)));
tmp_str_nor_xor = sprintf('arm2_nor vs dosage_imputed_val_xor %.10d',sum(abs((tmp_arm2_perm_nor_>0) - tmp_dosage_imputed_val_xor_)));
tmp_str_nor_nor = sprintf('arm2_nor vs dosage_imputed_val_nor %.10d',sum(abs((tmp_arm2_perm_nor_>0) - tmp_dosage_imputed_val_nor_)));
disp(sprintf(' %% %% ndud %.10d tmp_bo2_nor %.10d %s (%s) , %s , %s , %s',ndud,tmp_bo2_nor,dud_{ndud},alleletype_arm2_{tmp_bo2_nor},tmp_str_nor_and,tmp_str_nor_xor,tmp_str_nor_nor));
end;%if (length(tmp_bo2_nor)>0);
%%%%%%%%;
end;%if (length(tmp_bo2_)>0);
end;%if flag_check;
%%%%%%%%;
if flag_V_trnx;
tmp_dosage_imputed_AnV_ = tmp_dosage_imputed_AnV_ ...
   + tmp_dosage_imputed_val_and_normalized_*V_trnx_perm_and_(ndud,:) ...
   + tmp_dosage_imputed_val_xor_normalized_*V_trnx_perm_xor_(ndud,:) ...
   + tmp_dosage_imputed_val_nor_normalized_*V_trnx_perm_nor_(ndud,:) ...
  ;
tmp_dosage_rounded_AnV_ = tmp_dosage_rounded_AnV_ ...
   + tmp_dosage_rounded_val_and_normalized_*V_trnx_perm_and_(ndud,:) ...
   + tmp_dosage_rounded_val_xor_normalized_*V_trnx_perm_xor_(ndud,:) ...
   + tmp_dosage_rounded_val_nor_normalized_*V_trnx_perm_nor_(ndud,:) ...
  ;
tmp_dosage_genotyp_AnV_ = tmp_dosage_genotyp_AnV_ ...
   + tmp_dosage_genotyp_val_and_normalized_*V_trnx_perm_and_(ndud,:) ...
   + tmp_dosage_genotyp_val_xor_normalized_*V_trnx_perm_xor_(ndud,:) ...
   + tmp_dosage_genotyp_val_nor_normalized_*V_trnx_perm_nor_(ndud,:) ...
  ;
n_z_a = n_z_a + (V_trnx_perm_and_(ndud,:)~=0) + (V_trnx_perm_xor_(ndud,:)~=0) + (V_trnx_perm_nor_(ndud,:)~=0);
end;%if flag_V_trnx;
if flag_V_tstx;
tmp_dosage_imputed_AnV_ = tmp_dosage_imputed_AnV_ ...
   + tmp_dosage_imputed_val_and_normalized_*V_tstx_perm_and_(ndud,:) ...
   + tmp_dosage_imputed_val_xor_normalized_*V_tstx_perm_xor_(ndud,:) ...
   + tmp_dosage_imputed_val_nor_normalized_*V_tstx_perm_nor_(ndud,:) ...
  ;
tmp_dosage_rounded_AnV_ = tmp_dosage_rounded_AnV_ ...
   + tmp_dosage_rounded_val_and_normalized_*V_tstx_perm_and_(ndud,:) ...
   + tmp_dosage_rounded_val_xor_normalized_*V_tstx_perm_xor_(ndud,:) ...
   + tmp_dosage_rounded_val_nor_normalized_*V_tstx_perm_nor_(ndud,:) ...
  ;
tmp_dosage_genotyp_AnV_ = tmp_dosage_genotyp_AnV_ ...
   + tmp_dosage_genotyp_val_and_normalized_*V_tstx_perm_and_(ndud,:) ...
   + tmp_dosage_genotyp_val_xor_normalized_*V_tstx_perm_xor_(ndud,:) ...
   + tmp_dosage_genotyp_val_nor_normalized_*V_tstx_perm_nor_(ndud,:) ...
  ;
n_z_a = n_z_a + (V_tstx_perm_and_(ndud,:)~=0) + (V_tstx_perm_xor_(ndud,:)~=0) + (V_tstx_perm_nor_(ndud,:)~=0);
end;%if flag_V_tstx;
end;%if (length(neu1)==1);
end;%for ndosage_snp=1:n_dosage_snp;
fclose(fid);
%%%%%%%%;
disp(sprintf(' %% flag_bo1_all_sum = %d; flag_bo2_and_sum = %d; flag_bo2_xor_sum = %d; flag_bo2_nor_sum = %d',flag_bo1_all_sum,flag_bo2_and_sum,flag_bo2_xor_sum,flag_bo2_nor_sum));
tmp_h_mismatch_d1_imputed_ = hist(mismatch_d1_imputed__{nstudy_arm2},[0,1,2,3,8,9]);
tmp_h_mismatch_d1_genotyp_ = hist(mismatch_d1_genotyp__{nstudy_arm2},[0,1,2,3,8,9]);
disp(sprintf(' %% tmp_h_d1_mismatched_imputed_:       null    synonym  dominance       both  ambiguous   missing '));
disp(sprintf(' %% tmp_h_d1_mismatched_imputed_: %0.10d %0.10d %0.10d %0.10d %0.10d %0.10d',tmp_h_mismatch_d1_imputed_));
disp(sprintf(' %% tmp_h_d1_mismatched_genotyp_: %0.10d %0.10d %0.10d %0.10d %0.10d %0.10d',tmp_h_mismatch_d1_genotyp_));
tmp_h_mismatch_d2_imputed_ = hist(mismatch_d2_imputed__{nstudy_arm2},[0,1,2,3,8,9]);
tmp_h_mismatch_d2_genotyp_ = hist(mismatch_d2_genotyp__{nstudy_arm2},[0,1,2,3,8,9]);
disp(sprintf(' %% tmp_h_d1_mismatched_imputed_:       null    synonym  dominance       both  ambiguous   missing '));
disp(sprintf(' %% tmp_h_d2_mismatched_imputed_: %0.10d %0.10d %0.10d %0.10d %0.10d %0.10d',tmp_h_mismatch_d2_imputed_));
disp(sprintf(' %% tmp_h_d2_mismatched_genotyp_: %0.10d %0.10d %0.10d %0.10d %0.10d %0.10d',tmp_h_mismatch_d2_genotyp_));
tmp_mismatch_d1_imputed_swap_ = ( (mismatch_d1_imputed__{nstudy_arm2}==error_code_dominance) | (mismatch_d1_imputed__{nstudy_arm2}==error_code_both) );
tmp_mismatch_d1_genotyp_swap_ = ( (mismatch_d1_genotyp__{nstudy_arm2}==error_code_dominance) | (mismatch_d1_genotyp__{nstudy_arm2}==error_code_both) );
tmp_mismatch_d2_imputed_swap_ = ( (mismatch_d2_imputed__{nstudy_arm2}==error_code_dominance) | (mismatch_d2_imputed__{nstudy_arm2}==error_code_both) );
tmp_mismatch_d2_genotyp_swap_ = ( (mismatch_d2_genotyp__{nstudy_arm2}==error_code_dominance) | (mismatch_d2_genotyp__{nstudy_arm2}==error_code_both) );
tmp_mismatch_d1_imputed_keep_ = ( (mismatch_d1_imputed__{nstudy_arm2}==error_code_null) | (mismatch_d1_imputed__{nstudy_arm2}==error_code_synonym) );
tmp_mismatch_d1_genotyp_keep_ = ( (mismatch_d1_genotyp__{nstudy_arm2}==error_code_null) | (mismatch_d1_genotyp__{nstudy_arm2}==error_code_synonym) );
tmp_mismatch_d2_imputed_keep_ = ( (mismatch_d2_imputed__{nstudy_arm2}==error_code_null) | (mismatch_d2_imputed__{nstudy_arm2}==error_code_synonym) );
tmp_mismatch_d2_genotyp_keep_ = ( (mismatch_d2_genotyp__{nstudy_arm2}==error_code_null) | (mismatch_d2_genotyp__{nstudy_arm2}==error_code_synonym) );
tmp_mismatch_dx_imputed_keep_ = tmp_mismatch_d1_imputed_keep_ & tmp_mismatch_d2_imputed_keep_ ;
tmp_mismatch_dx_genotyp_keep_ = tmp_mismatch_d1_genotyp_keep_ & tmp_mismatch_d2_genotyp_keep_ ;
tmp_mismatch_dx_imputed_swap_ = tmp_mismatch_d1_imputed_swap_ & tmp_mismatch_d2_imputed_swap_ ;
tmp_mismatch_dx_genotyp_swap_ = tmp_mismatch_d1_genotyp_swap_ & tmp_mismatch_d2_genotyp_swap_ ;
tmp_mismatch_dx_imputed_miss_ = tmp_mismatch_d1_imputed_swap_ & tmp_mismatch_d2_imputed_keep_ ;
tmp_mismatch_dx_genotyp_miss_ = tmp_mismatch_d1_genotyp_swap_ & tmp_mismatch_d2_genotyp_keep_ ;
assert(sum( tmp_mismatch_dx_imputed_keep_ & tmp_mismatch_dx_imputed_swap_ )==0);
assert(sum( tmp_mismatch_dx_imputed_keep_ & tmp_mismatch_dx_imputed_miss_ )==0);
assert(sum( tmp_mismatch_dx_imputed_swap_ & tmp_mismatch_dx_imputed_miss_ )==0);
assert(sum( tmp_mismatch_dx_genotyp_keep_ & tmp_mismatch_dx_genotyp_swap_ )==0);
assert(sum( tmp_mismatch_dx_genotyp_keep_ & tmp_mismatch_dx_genotyp_miss_ )==0);
assert(sum( tmp_mismatch_dx_genotyp_swap_ & tmp_mismatch_dx_genotyp_miss_ )==0);
disp(sprintf(' %% imputed: keep %d swap %d miss %d; genotyp: keep %d swap %d miss %d'...
,sum(tmp_mismatch_dx_imputed_keep_),sum(tmp_mismatch_dx_imputed_swap_),sum(tmp_mismatch_dx_imputed_miss_)...
,sum(tmp_mismatch_dx_genotyp_keep_),sum(tmp_mismatch_dx_genotyp_swap_),sum(tmp_mismatch_dx_genotyp_miss_)...
	     ));
%%%%%%%%;

tmp_imputed_AnV_arm2_ = zeros(length(fu2_),2);
tmp_rounded_AnV_arm2_ = zeros(length(fu2_),2);
tmp_genotyp_AnV_arm2_ = zeros(length(fu2_),2);
for nhu2=1:length(hu2_);
tmp_fu2_ = find(sum(fu2_by_hu2_xref_(:,nhu2),2));
if (length(tmp_fu2_)>1); disp(sprintf(' %% Warning! tmp_fu2_ length %d in lisa_dosage_pca_ver9',length(tmp_fu2_))); end;
nfu2 = tmp_fu2_(1);
tmp_gud_ = find(sum(gud_by_hu2_xref_(:,nhu2),2));
if (length(tmp_gud_)>1); disp(sprintf(' %% Warning! tmp_gud_ length %d in lisa_dosage_pca_ver9',length(tmp_gud_))); end;
ngud = tmp_gud_(1);
assert(strcmp(fu2_{nfu2},gud_{ngud}));
tmp_imputed_AnV_arm2_(nfu2,:) = tmp_dosage_imputed_AnV_(ngud,:);
tmp_rounded_AnV_arm2_(nfu2,:) = tmp_dosage_rounded_AnV_(ngud,:);
tmp_genotyp_AnV_arm2_(nfu2,:) = tmp_dosage_genotyp_AnV_(ngud,:);
end;%for nhu2=1:length(hu2_);
tmp_imputed_DnV_arm2_ = zeros(length(fu2_),2); tmp_Dij_ = find(mr_D_arm2_); tmp_imputed_DnV_arm2_(tmp_Dij_,:) = tmp_imputed_AnV_arm2_(tmp_Dij_,:);
tmp_imputed_XnV_arm2_ = zeros(length(fu2_),2); tmp_Xij_ = find(mr_X_arm2_); tmp_imputed_XnV_arm2_(tmp_Xij_,:) = tmp_imputed_AnV_arm2_(tmp_Xij_,:);
tmp_rounded_DnV_arm2_ = zeros(length(fu2_),2); tmp_Dij_ = find(mr_D_arm2_); tmp_rounded_DnV_arm2_(tmp_Dij_,:) = tmp_rounded_AnV_arm2_(tmp_Dij_,:);
tmp_rounded_XnV_arm2_ = zeros(length(fu2_),2); tmp_Xij_ = find(mr_X_arm2_); tmp_rounded_XnV_arm2_(tmp_Xij_,:) = tmp_rounded_AnV_arm2_(tmp_Xij_,:);
tmp_genotyp_DnV_arm2_ = zeros(length(fu2_),2); tmp_Dij_ = find(mr_D_arm2_); tmp_genotyp_DnV_arm2_(tmp_Dij_,:) = tmp_genotyp_AnV_arm2_(tmp_Dij_,:);
tmp_genotyp_XnV_arm2_ = zeros(length(fu2_),2); tmp_Xij_ = find(mr_X_arm2_); tmp_genotyp_XnV_arm2_(tmp_Xij_,:) = tmp_genotyp_AnV_arm2_(tmp_Xij_,:);

flag_plot= (0 & ~flag_lisa);
if flag_plot;
figure(n_figure); n_figure = n_figure+1; clf;
subplot(2,2,1); hold on;
tmp_Dij_ = intersect(gud_to_fu2_,find(mr_D_arm2_));
plot(DnV_tstx_arm2_(tmp_Dij_,1),tmp_imputed_DnV_arm2_(tmp_Dij_,1),'bo'); ylim([-5,5]);
plot(DnV_tstx_arm2_(tmp_Dij_,1),tmp_rounded_DnV_arm2_(tmp_Dij_,1),'g+'); ylim([-5,5]);
plot(DnV_tstx_arm2_(tmp_Dij_,1),tmp_genotyp_DnV_arm2_(tmp_Dij_,1),'rx'); ylim([-5,5]);
title('D 1');
subplot(2,2,2); hold on;
tmp_Xij_ = intersect(gud_to_fu2_,find(mr_X_arm2_));
plot(XnV_tstx_arm2_(tmp_Xij_,1),tmp_imputed_XnV_arm2_(tmp_Xij_,1),'bo'); ylim([-5,5]);
plot(XnV_tstx_arm2_(tmp_Xij_,1),tmp_rounded_XnV_arm2_(tmp_Xij_,1),'g+'); ylim([-5,5]);
plot(XnV_tstx_arm2_(tmp_Xij_,1),tmp_genotyp_XnV_arm2_(tmp_Xij_,1),'rx'); ylim([-5,5]);
title('X 1');
subplot(2,2,3); hold on;
tmp_Dij_ = intersect(gud_to_fu2_,find(mr_D_arm2_));
plot(DnV_tstx_arm2_(tmp_Dij_,2),tmp_imputed_DnV_arm2_(tmp_Dij_,2),'bo'); ylim([-5,5]);
plot(DnV_tstx_arm2_(tmp_Dij_,2),tmp_rounded_DnV_arm2_(tmp_Dij_,2),'g+'); ylim([-5,5]);
plot(DnV_tstx_arm2_(tmp_Dij_,2),tmp_genotyp_DnV_arm2_(tmp_Dij_,2),'rx'); ylim([-5,5]);
title('D 2');
subplot(2,2,4); hold on;
tmp_Xij_ = intersect(gud_to_fu2_,find(mr_X_arm2_));
plot(XnV_tstx_arm2_(tmp_Xij_,2),tmp_imputed_XnV_arm2_(tmp_Xij_,2),'bo'); ylim([-5,5]);
plot(XnV_tstx_arm2_(tmp_Xij_,2),tmp_rounded_XnV_arm2_(tmp_Xij_,2),'g+'); ylim([-5,5]);
plot(XnV_tstx_arm2_(tmp_Xij_,2),tmp_genotyp_XnV_arm2_(tmp_Xij_,2),'rx'); ylim([-5,5]);
title('X 2'); hold on;
end;%if flag_plot;

dosage_imputed_DnV_arm2__{nstudy_arm2} = tmp_imputed_DnV_arm2_ ;
dosage_imputed_XnV_arm2__{nstudy_arm2} = tmp_imputed_XnV_arm2_ ;
dosage_imputed_DnV_arm2_ = dosage_imputed_DnV_arm2_ + tmp_imputed_DnV_arm2_ ;
dosage_imputed_XnV_arm2_ = dosage_imputed_XnV_arm2_ + tmp_imputed_XnV_arm2_ ;
dosage_rounded_DnV_arm2__{nstudy_arm2} = tmp_rounded_DnV_arm2_ ;
dosage_rounded_XnV_arm2__{nstudy_arm2} = tmp_rounded_XnV_arm2_ ;
dosage_rounded_DnV_arm2_ = dosage_rounded_DnV_arm2_ + tmp_rounded_DnV_arm2_ ;
dosage_rounded_XnV_arm2_ = dosage_rounded_XnV_arm2_ + tmp_rounded_XnV_arm2_ ;
dosage_genotyp_DnV_arm2__{nstudy_arm2} = tmp_genotyp_DnV_arm2_ ;
dosage_genotyp_XnV_arm2__{nstudy_arm2} = tmp_genotyp_XnV_arm2_ ;
dosage_genotyp_DnV_arm2_ = dosage_genotyp_DnV_arm2_ + tmp_genotyp_DnV_arm2_ ;
dosage_genotyp_XnV_arm2_ = dosage_genotyp_XnV_arm2_ + tmp_genotyp_XnV_arm2_ ;

if (flag_lisa==1); tmp_command = sprintf('rm -rf %s',fname_dosage); disp(tmp_command); system(tmp_command); end;

save(fname_stage_1,...
      'dosage_imputed_DnV_arm2__','dosage_imputed_XnV_arm2__','dosage_imputed_DnV_arm2_','dosage_imputed_XnV_arm2_'...
     ,'dosage_rounded_DnV_arm2__','dosage_rounded_XnV_arm2__','dosage_rounded_DnV_arm2_','dosage_rounded_XnV_arm2_'...
     ,'dosage_genotyp_DnV_arm2__','dosage_genotyp_XnV_arm2__','dosage_genotyp_DnV_arm2_','dosage_genotyp_XnV_arm2_'...
     ,'mismatch_d1_imputed__','mismatch_d1_genotyp__'...
     ,'mismatch_d2_imputed__','mismatch_d2_genotyp__'...
     );

%%%%%%%%%%%%%%%%;
end;%for nstudy_arm2=1:length(study_name_arm2_);

save(fname_stage_1,...
      'dosage_imputed_DnV_arm2__','dosage_imputed_XnV_arm2__','dosage_imputed_DnV_arm2_','dosage_imputed_XnV_arm2_'...
     ,'dosage_rounded_DnV_arm2__','dosage_rounded_XnV_arm2__','dosage_rounded_DnV_arm2_','dosage_rounded_XnV_arm2_'...
     ,'dosage_genotyp_DnV_arm2__','dosage_genotyp_XnV_arm2__','dosage_genotyp_DnV_arm2_','dosage_genotyp_XnV_arm2_'...
     ,'mismatch_d1_imputed__','mismatch_d1_genotyp__'...
     ,'mismatch_d2_imputed__','mismatch_d2_genotyp__'...
     );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_stage==1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
