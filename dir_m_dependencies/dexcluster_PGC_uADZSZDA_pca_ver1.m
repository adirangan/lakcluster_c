flag_setup=0;

if flag_setup;
% Now set up row- and col-masks which correspond to a particular iteration from PGC_cl4_maf01_dex_p25_D_m2r2_g004. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% First load PGC_cl4_maf01_dex_p25_D_m2r2_g004. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
dir_code = sprintf('/data/rangan/dir_bcc/dir_lakcluster_c_dev'); 
dir_trunk = sprintf('/data/rangan/dir_bcc/dir_PGC_20180304');
path(path,sprintf('%s/dir_m',dir_code)); 
%%%%%%%%%%%%%%%% ;
cl_num = 4; 
[study_trunk_,study_name_,n_study] = lisa_define_study_ver0(cl_num);
Ireq_half = floor(n_study/2); Icat_full = n_study;
string_prefix = sprintf('PGC_cl%d_maf01',cl_num);
flag_dex_vs_lak = 'dex'; 
gamma=0.004;B_MLT=34;n_mds=20;
flag_reverse = 0; 
n_maf = 5;
n_cov = 1;
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
n_shuffle = 64;
dir__in = sprintf('%s/dir_%s',dir_trunk,string_prefix);
dir_out = sprintf('%s/dir_%s_analyze',dir_trunk,string_prefix);
disp(sprintf(' %% dir_out: %s',dir_out));
%%%%%%%%%%%%%%%% ;
lisa_xdropplot_loadtrace_ver0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Now intersect columns in arm1 with columns from arm2 (replication arm). ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
n_snp = length(bim_{1});
bim_condensed_ = cell(n_snp,1);
for nsnp=1:n_snp;
bim_condensed_{nsnp} = sprintf('%d_%s_%s',bim_{1}(nsnp),bim_{2}{nsnp},bim_{7}{nsnp});
end;%for nsnp=1:n_snp;
%%%%%%%%%%%%%%%% ;
cl_num_arm2 = 1;
string_prefix_arm2 = sprintf('PGC_cl%d_maf01',cl_num_arm2);
dir__in_arm2 = sprintf('%s/dir_%s',dir_trunk,string_prefix_arm2);
dir_out_arm2 = sprintf('%s/dir_%s_analyze',dir_trunk,string_prefix_arm2);
disp(sprintf(' %% cl_num_arm2 %d --> %s',cl_num_arm2,string_prefix_arm2));
[study_trunk_arm2_,study_name_arm2_,n_study_arm2] = lisa_define_study_ver0(cl_num_arm2);
fcheck(sprintf('%s/%s_bim.ext',dir__in_arm2,string_prefix_arm2));
fid = fopen(sprintf('%s/%s_bim.ext',dir__in_arm2,string_prefix_arm2)); 
bim_arm2_ = textscan(fid,'%d\t%s\t%d\t%d\t%c\t%c\t%s\t%f\t%f\t%f\t%f\n','headerlines',0); fclose(fid);
n_snp_arm2 = length(bim_arm2_{1});
bim_arm2_condensed_ = cell(n_snp_arm2,1);
for nsnp=1:n_snp_arm2;
bim_arm2_condensed_{nsnp} = sprintf('%d_%s_%s',bim_arm2_{1}(nsnp),bim_arm2_{2}{nsnp},bim_arm2_{7}{nsnp});
end;%for nsnp=1:n_snp_arm2;
%%%%%%%%%%%%%%%% ;
[~,i_arm1,i_arm2] = intersect(bim_condensed_,bim_arm2_condensed_,'stable');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;%if flag_setup;

niteration_ = 0:50:750; niteration(1)=1; niteration_(end) = length(r_rem_)-1;
niteration_ = niteration_(length(niteration_):-1:1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Now set niteration and generate row- and col-masks. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
for ni=1:length(niteration_);
niteration = niteration_(ni);

disp(sprintf(' %% ni %d/%d, niteration %d',ni,length(niteration_),niteration));
%niteration = 400; % drawn from ../dir_PGC_cl4_maf01_analyze/dir_PGC_cl4_maf01_dex_p25_D_m2r2_g004/PGC_cl4_maf01_dex_p25_D_m2r2_g004_trace_A.jpg. ;
mr_A_tmp = mr_A_full*0; mr_A_tmp(rdrop_a_(length(rdrop_a_)-r_rem_(niteration)+1:length(rdrop_a_)))=1;
tmpchar = sprintf('%s/%s_mc_A.b16',dir_out_s0000,string_name_s0000); 
[bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); mc_A = tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols)>0;
mc_A_tmp = mc_A*0; mc_A_tmp(cdrop_a_(length(cdrop_a_)-c_rem_(niteration)+1:length(cdrop_a_)))=1;
%%%%%%%%%%%%%%%% ;
% mc_A_tm2 contains columns for initial arm. ;
%%%%%%%%%%%%%%%% ;
[~,ii] = intersect(i_arm1,find(mc_A_tmp),'stable');
for nc=1:length(ii); assert(strcmp(bim_arm2_condensed_{i_arm2(ii(nc))},bim_condensed_{i_arm1(ii(nc))})); end;%for nc=1:length(ii);
mc_A_tm2 = mc_A*0; mc_A_tm2(i_arm1(ii))=1;
%%%%%%%%%%%%%%%% ;
% mc_A_tm2_arm2 contains columns for replication arm. ;
%%%%%%%%%%%%%%%% ;
tmpchar = sprintf('%s/PGC_cl%d_maf01_mc_A.b16',dir__in_arm2,cl_num_arm2);
[bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); mc_A_arm2 = tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols)>0;
mc_A_tm2_arm2 = mc_A_arm2*0; mc_A_tm2_arm2(i_arm2(ii))=1;
j_ = find(mc_A_tm2); j_arm2_ = find(mc_A_tm2_arm2); assert(length(j_)==length(j_arm2_));
assert(length(intersect(bim_arm2_condensed_(j_arm2_),bim_condensed_(j_)))==length(j_));
%%%%%%%%%%%%%%%% ;
pca_infix = sprintf('pca_ni%d_clb%d',niteration,cl_num_arm2);
%%%%%%%%%%%%%%%% ;
% Assuming Icat==1 (i.e., full) for single row mask. ;
tmpchar = sprintf('%s/%s_mr_A_full_%s_arm1.b16',dir_out_s0000,string_name_s0000,pca_infix); %delete(tmpchar);
tutorial_binary_compress(bitj,mr_A_tmp(:)>0,tmpchar);
tmpchar = sprintf('%s/%s_mc_A_%s_arm1.b16',dir_out_s0000,string_name_s0000,pca_infix); %delete(tmpchar);
tutorial_binary_compress(bitj,mc_A_tm2(:)>0,tmpchar);
tmpchar = sprintf('%s/%s_mc_A_%s_arm2.b16',dir_out_s0000,string_name_s0000,pca_infix); %delete(tmpchar);
tutorial_binary_compress(bitj,mc_A_tm2_arm2(:)>0,tmpchar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% setting up input file for pca_driver. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
nshuffle=0;
flag_verbose=0;
d_inpre = sprintf('%s/%s',dir__in,string_prefix); 
d_oupre = sprintf('%s/%s',dir_out_s0000,string_name_s0000);
pca_b_mlt = 44;
pca_tolerance = 1e-2;
pca_rank = 2;
flag_T = 1;
mds_str = sprintf('m%dr%d',length(mds_used_),mds_repl);
mds_kappa_squared = textread(sprintf('%s_T_%s_kappa.txt',d_inpre,mds_str));
T_n_cols = 1+length(mds_used_)*mds_repl;
Y_n_cols=0;
M_n_rows_ = zeros(Icat,1); M_n_cols = length(mc_A);
for nb=0:Icat-1; M_n_rows_(1+nb) = length(mr_A_{1+nb});end;%for nb=0:Icat-1;
pca_infix = sprintf('pca_ni%d_clb%d',niteration,cl_num_arm2);
fname__in = sprintf('%s/%s_%s.in',dir_out_s0000,string_name_s0000,pca_infix);
fp = fopen(fname__in,'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fprintf(fp,'GLOBAL_verbose= %d;\n',flag_verbose);
fprintf(fp,'GLOBAL_thread_count= 15;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','pca_driver');
fprintf(fp,'GLOBAL_pca_infix= %s;\n',pca_infix);
fprintf(fp,'GLOBAL_pca_iteration_num= 1;\n');
fprintf(fp,'GLOBAL_pca_iteration_max= 0;\n');
fprintf(fp,'GLOBAL_pca_iteration_min= 0;\n');
fprintf(fp,'GLOBAL_pca_rank= %d;\n',pca_rank);
fprintf(fp,'GLOBAL_pca_tolerance= %0.16f;\n',pca_tolerance);
fprintf(fp,'GLOBAL_TEST_niter= 1024;\n');
fprintf(fp,'GLOBAL_NBINS= %d;\n',Icat);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',pca_b_mlt);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
if (flag_T==1); fprintf(fp,'GLOBAL_kappa_squared= %0.16f;\n',mds_kappa_squared); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if (nshuffle==0); A_name_string = '_A_'; end;%if (nshuffle==0);
if (nshuffle>0); A_name_string = '_A_shuffle_'; end;%if (nshuffle>0);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_name_= ',Icat,d_inpre,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_t_name_= ',Icat,d_inpre,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_A_n_rows_= ',Icat,M_n_rows_);
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',M_n_cols);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',Icat,d_oupre,'_mr_Z_',sprintf('_%s_arm1.b16',pca_infix)); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',Icat,d_oupre,'_mr_A_',sprintf('_%s_arm1.b16',pca_infix)); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A_%s_arm1.b16;\n',d_oupre,pca_infix);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_name_= ',Icat,d_inpre,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_t_name_= ',Icat,d_inpre,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_Z_n_rows_= ',Icat,M_n_rows_);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',Icat,d_oupre,'_mr_A_','.b16'); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',Icat,d_oupre,'_mr_Z_','.b16'); end;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
if (flag_T==0);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_crop_cols); 
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_T_n_name_= ',Icat,d_oupre,'_T_crop_','_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_T_t_name_= ',Icat,d_oupre,'_T_crop_','_t.b16');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T_crop.b16;\n',d_oupre);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_S_n_name_= ',Icat,d_oupre,'_T_crop_','_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_S_t_name_= ',Icat,d_oupre,'_T_crop_','_t.b16');
end;%if (flag_T==0);
if (flag_T==1);
fprintf(fp,'GLOBAL_T_n_cols= %d;\n',T_n_cols); 
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_T_n_name_= ',Icat,d_inpre,sprintf('_T_%s_',mds_str),'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_T_t_name_= ',Icat,d_inpre,sprintf('_T_%s_',mds_str),'_t.b16');
fprintf(fp,'GLOBAL_T_n_cind= %s_mc_T_%s.b16;\n',d_oupre,mds_str);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_S_n_name_= ',Icat,d_inpre,sprintf('_T_%s_',mds_str),'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_S_t_name_= ',Icat,d_inpre,sprintf('_T_%s_',mds_str),'_t.b16');
end;%if (flag_T==1);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out_s0000);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by pca_PGC_uADZSZDA_ver0.m on %s;\n',date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fclose(fp);
type(fname__in);
disp(sprintf(' %% fname__in:\n%s',fname__in));
command_str = sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in);
system(command_str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% Now load V_ and generate pca for arm2. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
pca_fname_V_ = sprintf('%s/%s_k%d_B%d_V_.mda',dir_out_s0000,pca_infix,pca_rank,pca_b_mlt);
V_ = mda_read_r8(pca_fname_V_);
for (nr=1:pca_rank); assert(sum(V_(:,nr)~=0)<=length(mc_A_tm2)); assert(sum(V_(:,nr))==sum(V_(:,nr).*mc_A_tm2)); end;
V_arm2_ = zeros(length(mc_A_arm2),pca_rank);
for nc=1:length(ii); V_arm2_(i_arm2(ii(nc)),:) = V_(i_arm1(ii(nc)),:); end;%for nc=1:length(ii);
for (nr=1:pca_rank); assert(sum(V_arm2_(:,nr)~=0)<=length(mc_A_tm2_arm2)); assert(sum(V_arm2_(:,nr))==sum(V_arm2_(:,nr).*mc_A_tm2_arm2)); end;
pca_fname_V_arm2_ = sprintf('%s/%s_k%d_B%d_V_arm2_.mda',dir_out_s0000,pca_infix,pca_rank,pca_b_mlt);
mda_write_d3_r8(V_arm2_,pca_fname_V_arm2_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% setting up input file for pca_proj arm1. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
nshuffle=0;
flag_verbose=0;
d_inpre = sprintf('%s/%s',dir__in,string_prefix); 
d_oupre = sprintf('%s/%s',dir_out_s0000,string_name_s0000);
pca_b_mlt = 44;
pca_tolerance = 1e-2;
pca_rank = 2;
mds_str = sprintf('m%dr%d',length(mds_used_),mds_repl);
Y_n_cols=0;
M_n_rows_ = zeros(Icat,1); M_n_cols = length(mc_A);
for nb=0:Icat-1; M_n_rows_(1+nb) = length(mr_A_{1+nb});end;%for nb=0:Icat-1;
pca_proj_arm1_infix = sprintf('pca_proj_arm1_ni%d_clb%d',niteration,cl_num_arm2);
fname__in = sprintf('%s/%s_%s.in',dir_out_s0000,string_name_s0000,pca_proj_arm1_infix);
fp = fopen(fname__in,'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fprintf(fp,'GLOBAL_verbose= %d;\n',flag_verbose);
fprintf(fp,'GLOBAL_thread_count= 15;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','pca_proj_driver');
fprintf(fp,'GLOBAL_pca_infix= %s;\n',pca_proj_arm1_infix);
fprintf(fp,'GLOBAL_pca_V_= %s;\n',pca_fname_V_);
fprintf(fp,'GLOBAL_pca_iteration_num= 1;\n');
fprintf(fp,'GLOBAL_pca_iteration_max= 0;\n');
fprintf(fp,'GLOBAL_pca_iteration_min= 0;\n');
fprintf(fp,'GLOBAL_pca_rank= %d;\n',pca_rank);
fprintf(fp,'GLOBAL_pca_tolerance= %0.16f;\n',pca_tolerance);
fprintf(fp,'GLOBAL_TEST_niter= 1024;\n');
fprintf(fp,'GLOBAL_NBINS= %d;\n',Icat);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',pca_b_mlt);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if (nshuffle==0); A_name_string = '_A_'; end;%if (nshuffle==0);
if (nshuffle>0); A_name_string = '_A_shuffle_'; end;%if (nshuffle>0);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_name_= ',Icat,d_inpre,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_t_name_= ',Icat,d_inpre,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_A_n_rows_= ',Icat,M_n_rows_);
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',M_n_cols);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',Icat,d_oupre,'_mr_Z_','.b16'); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',Icat,d_oupre,'_mr_A_','.b16'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A_%s_arm1.b16;\n',d_oupre,pca_infix);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_name_= ',Icat,d_inpre,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_t_name_= ',Icat,d_inpre,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_Z_n_rows_= ',Icat,M_n_rows_);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',Icat,d_oupre,'_mr_A_','.b16'); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',Icat,d_oupre,'_mr_Z_','.b16'); end;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out_s0000);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by pca_PGC_uADZSZDA_ver0.m on %s;\n',date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fclose(fp);
type(fname__in);
disp(sprintf(' %% fname__in:\n%s',fname__in));
command_str = sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in);
system(command_str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% setting up input file for pca_proj arm2. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
nshuffle=0;
flag_verbose=0;
d_inpre_arm2 = sprintf('%s/%s',dir__in_arm2,string_prefix_arm2); 
d_oupre = sprintf('%s/%s',dir_out_s0000,string_name_s0000);
pca_b_mlt = 44;
pca_tolerance = 1e-2;
pca_rank = 2;
mds_str = sprintf('m%dr%d',length(mds_used_),mds_repl);
Y_n_cols=0;
tmpchar = sprintf('%s/PGC_cl%d_maf01_mr_A_full.b16',dir__in_arm2,cl_num_arm2);
[bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); mr_A_full_arm2 = tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols)>0;
tmpchar = sprintf('%s/PGC_cl%d_maf01_mr_Z_full.b16',dir__in_arm2,cl_num_arm2);
[bitj,nrows,ncols] = tutorial_binary_getsize(tmpchar); mr_Z_full_arm2 = tutorial_binary_uncompress(tmpchar,1:nrows,1:ncols)>0;
Icat = 1; M_n_rows_ = length(mr_A_full_arm2);
M_n_cols = length(mc_A_arm2);
pca_proj_arm2_infix = sprintf('pca_proj_arm2_ni%d_clb%d',niteration,cl_num_arm2);
fname__in = sprintf('%s/%s_%s.in',dir_out_s0000,string_name_s0000,pca_proj_arm2_infix);
fp = fopen(fname__in,'w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fprintf(fp,'GLOBAL_verbose= %d;\n',flag_verbose);
fprintf(fp,'GLOBAL_thread_count= 15;\n');
fprintf(fp,'GLOBAL_omp_type= 1;\n');
fprintf(fp,'GLOBAL_TEST_TYPE= %s;\n','pca_proj_driver');
fprintf(fp,'GLOBAL_pca_infix= %s;\n',pca_proj_arm2_infix);
fprintf(fp,'GLOBAL_pca_V_= %s;\n',pca_fname_V_arm2_);
fprintf(fp,'GLOBAL_pca_iteration_num= 1;\n');
fprintf(fp,'GLOBAL_pca_iteration_max= 0;\n');
fprintf(fp,'GLOBAL_pca_iteration_min= 0;\n');
fprintf(fp,'GLOBAL_pca_rank= %d;\n',pca_rank);
fprintf(fp,'GLOBAL_pca_tolerance= %0.16f;\n',pca_tolerance);
fprintf(fp,'GLOBAL_TEST_niter= 1024;\n');
fprintf(fp,'GLOBAL_NBINS= %d;\n',Icat);
fprintf(fp,'GLOBAL_B_MLT= %d;\n',pca_b_mlt);
fprintf(fp,'GLOBAL_gamma= %0.4f;\n',gamma);
fprintf(fp,'GLOBAL_Ireq= %d;\n',Ireq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
if (nshuffle==0); A_name_string = '_A_'; end;%if (nshuffle==0);
if (nshuffle>0); A_name_string = '_A_shuffle_'; end;%if (nshuffle>0);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_name_= ',Icat,d_inpre_arm2,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_t_name_= ',Icat,d_inpre_arm2,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_A_n_rows_= ',Icat,M_n_rows_);
fprintf(fp,'GLOBAL_A_n_cols= %d;\n',M_n_cols);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',Icat,d_inpre_arm2,'_mr_Z_','.b16'); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_A_n_rind_= ',Icat,d_inpre_arm2,'_mr_A_','.b16'); end;
fprintf(fp,'GLOBAL_A_n_cind= %s_mc_A_%s_arm2.b16;\n',d_oupre,pca_infix);
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_name_= ',Icat,d_inpre_arm2,A_name_string,'_n.b16');
dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_t_name_= ',Icat,d_inpre_arm2,A_name_string,'_t.b16');
dexcluster_PGC_uADZSZDA_excerpt_1(fp,'GLOBAL_Z_n_rows_= ',Icat,M_n_rows_);
if (flag_reverse==1); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',Icat,d_inpre_arm2,'_mr_A_','.b16'); end;
if (flag_reverse==0); dexcluster_PGC_uADZSZDA_excerpt_0(fp,'GLOBAL_Z_n_rind_= ',Icat,d_inpre_arm2,'_mr_Z_','.b16'); end;
fprintf(fp,'GLOBAL_Y_n_cols= %d;\n',Y_n_cols);
fprintf(fp,'GLOBAL_DIR_NAME= %s;\n',dir_out_s0000);
fprintf(fp,'END= 0;\n');
fprintf(fp,'%% generated by pca_PGC_uADZSZDA_ver0.m on %s;\n',date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fclose(fp);
type(fname__in);
disp(sprintf(' %% fname__in:\n%s',fname__in));
command_str = sprintf('%s/lakcluster_ver18 < %s',dir_code,fname__in);
system(command_str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% plot results. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
pca_proj_arm1_infix = sprintf('pca_proj_arm1_ni%d_clb%d',niteration,cl_num_arm2);
pca_proj_arm2_infix = sprintf('pca_proj_arm2_ni%d_clb%d',niteration,cl_num_arm2);
mr_A_tmp = mr_A_full*0; mr_A_tmp(rdrop_a_(length(rdrop_a_)-r_rem_(niteration)+1:length(rdrop_a_)))=1;
figure(ni);clf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
AnV_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_AnV_.mda',dir_out_s0000,pca_proj_arm1_infix,pca_rank,pca_b_mlt));
ZnV_arm1_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_ZnV_.mda',dir_out_s0000,pca_proj_arm1_infix,pca_rank,pca_b_mlt));
subplot(2,2,1);hold on; 
plot(ZnV_arm1_(find(mr_Z_full),1),ZnV_arm1_(find(mr_Z_full),2),'k.','MarkerSize',15);
plot(AnV_arm1_(find(mr_A_full-mr_A_tmp),1),AnV_arm1_(find(mr_A_full-mr_A_tmp),2),'g.','MarkerSize',15);
plot(AnV_arm1_(find(mr_A_tmp),1),AnV_arm1_(find(mr_A_tmp),2),'r.','MarkerSize',25);
hold off;
xlabel('PC1');ylabel('PC2'); title(sprintf('ni %d arm1',niteration));
subplot(2,2,2); hold on;
hlim = mean(union(ZnV_arm1_(:,1),AnV_arm1_(:,1))) + [-1,1]*3.5*std(union(ZnV_arm1_(:,1),AnV_arm1_(:,1)));
hbins=linspace(min(hlim),max(hlim),64);
h_Z_full = hist(ZnV_arm1_(find(mr_Z_full),1),hbins); h_Z_full = h_Z_full/sum(mr_Z_full);
h_A_full = hist(AnV_arm1_(find(mr_A_full),1),hbins); h_A_full = h_A_full/sum(mr_A_full);
h_A_exc = hist(AnV_arm1_(find(mr_A_full-mr_A_tmp),1),hbins); h_A_exc = h_A_exc/sum(mr_A_full);
h_A_tmp = hist(AnV_arm1_(find(mr_A_tmp),1),hbins); h_A_tmp = h_A_tmp/sum(mr_A_full);
stairs(hbins,h_Z_full,'k-','LineWidth',2);
stairs(hbins,h_A_full,'m-','LineWidth',2);
stairs(hbins,h_A_exc,'g-','LineWidth',2);
stairs(hbins,h_A_tmp,'r-','LineWidth',2);
hold off;
xlim(hlim);xlabel('PC1');ylabel('hist'); title(sprintf('ni %d arm1',niteration));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
AnV_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_AnV_.mda',dir_out_s0000,pca_proj_arm2_infix,pca_rank,pca_b_mlt));
ZnV_arm2_ = mda_read_r8(sprintf('%s/%s_k%d_B%d_ZnV_.mda',dir_out_s0000,pca_proj_arm2_infix,pca_rank,pca_b_mlt));
subplot(2,2,3);hold on; 
plot(ZnV_arm2_(find(mr_Z_full_arm2),1),ZnV_arm2_(find(mr_Z_full_arm2),2),'k.','MarkerSize',15);
plot(AnV_arm2_(find(mr_A_full_arm2),1),AnV_arm2_(find(mr_A_full_arm2),2),'m.','MarkerSize',15);
hold off;
xlabel('PC1');ylabel('PC2'); title(sprintf('ni %d arm2',niteration));
subplot(2,2,4); hold on;
hlim = mean(union(ZnV_arm2_(:,1),AnV_arm2_(:,1))) + [-1,1]*3.5*std(union(ZnV_arm2_(:,1),AnV_arm2_(:,1)));
hbins=linspace(min(hlim),max(hlim),64);
h_Z_full = hist(ZnV_arm2_(find(mr_Z_full_arm2),1),hbins); h_Z_full = h_Z_full/sum(mr_Z_full_arm2);
h_A_full = hist(AnV_arm2_(find(mr_A_full_arm2),1),hbins); h_A_full = h_A_full/sum(mr_A_full_arm2);
stairs(hbins,h_Z_full,'k-','LineWidth',2);
stairs(hbins,h_A_full,'m-','LineWidth',2);
hold off;
xlim(hlim);xlabel('PC1');ylabel('hist'); title(sprintf('ni %d arm2',niteration));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
drawnow();
fname_fig = sprintf('%s/%s_pca_proj_ni%d_clb%d.jpg',dir_out_s0000,string_name_s0000,niteration,cl_num_arm2); print('-djpeg',fname_fig);

end;%for ni=1:length(niteration_);