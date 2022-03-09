setup_0;

flag_verbose=1;
flag_disp=flag_verbose; nf=0;

if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Testing a version of xxxcluster_fromdisk_uADZSZDA_ver16_dr_0 ;')); end;
if (flag_verbose); disp(sprintf(' %% with a modified maf_lo_threshold. ;')); end;
if (flag_verbose); disp(sprintf(' %% First running test_stripped_xxxcluster_fromdisk_uADZSZDA_ver16. ;')); end;
dir_trunk = sprintf('%s/dir_test_xxxcluster_fromdisk_ver16',pwd);
test_stripped_xxxcluster_fromdisk_ver16(struct('flag_verbose',0));
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now running xxxcluster_fromdisk_uADZSZDA_ver16_dr_0 on same files. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note that, this time, we return the ''parameter'' structure. ;')); end;
if (flag_verbose); disp(sprintf(' %% We will later alter this structure to run a new analysis ;')); end;
if (flag_verbose); disp(sprintf(' %% on a subset of the snps (selected based on minor-allele-frequency). ;')); end;
parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.dir_trunk = dir_trunk;
parameter.study_name_of_branch_s_ = {'dir_study00'};
parameter.study_name_without_extension_s_ = {'study00'};
parameter.mds_name_plus_extension_s_ = {'study00_mds_tsv.txt'};
parameter.ent_cutoff = 0.0045022835561375; %<-- taken from test_stripped_xxxcluster_fromdisk_uADZSZDA_ver16;
parameter = xxxcluster_fromdisk_uADZSZDA_ver16_dr_0(parameter);
if (flag_verbose); disp(sprintf(' %% ;')); end;

if (flag_verbose); disp(sprintf(' %% Now we examine the bim.ext file. ;')); end;
fname_bim_ext = sprintf('%s/%s_bim.ext',parameter.dir_0in,parameter.str_prefix);
[ ...
 n_snp ...
,bim_khr_ ...
,bim_vid_ ...
,bim_gdi_ ...
,bim_pdi_ ...
,bim_al1_ ...
,bim_al2_ ...
,bim_alt_ ...
,bim_ent_ ...
,bim_frq_ ...
,bim_mss_ ...
,bim_maf_ ...
,bim_name_ ...
,bim_ ...
] = ...
load_bimext_ver1( ...
fname_bim_ext ...
);
bim_maf_p35 = prctile(bim_maf_,35);
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% For this synthetic data-set the minor-allele-frequencies (maf) ;')); end;
if (flag_verbose); disp(sprintf(' %% range from: %0.6f to %0.6f. ;',min(bim_maf_),max(bim_maf_))); end;
if (flag_verbose); disp(sprintf(' %% The 35-th percentile for the maf is: %0.6f ;',bim_maf_p35)); end;
if (flag_verbose); disp(sprintf(' %% As a test, we will set the maf_lo_threshold to equal this value. ;')); end;
if (flag_verbose); disp(sprintf(' %% Consequently, the number of used columns will shrink. ;')); end;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% First we copy the parameter structure, ;')); end;
if (flag_verbose); disp(sprintf(' %% and specify the new maf_lo_threshold ;')); end;
if (flag_verbose); disp(sprintf(' %% (which was previously 0.10 by default). ;')); end;
parameter_maflo = parameter;
parameter_maflo.maf_lo_threshold = bim_maf_p35; %<-- In a real data-set you might replace bim_maf_p35 with, say, 0.05. ;
for nshuffle=0:parameter_maflo.n_shuffle-1+1;
parameter_maflo.nshuffle = nshuffle;
parameter_maflo = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_maflo);
end;%for nshuffle=0:parameter_maflo.n_shuffle-1+1;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now the results are stored in a new output directory, named: ;')); end;
if (flag_verbose); disp(sprintf(' %% %s ;',parameter_maflo.dir_out_trace)); end;
if (flag_verbose); disp(sprintf(' %% with the new infix: ;')); end;
if (flag_verbose); disp(sprintf(' %% %s ;',parameter_maflo.str_name_s0000)); end;
if (flag_verbose); disp(sprintf(' %% Note that this is different from the previous infix: ;')); end;
if (flag_verbose); disp(sprintf(' %% %s ;',parameter.str_name_s0000)); end;
if (flag_verbose); disp(sprintf(' %% The ''pXX'' sequence indicates the maf_lo_threshold. ;')); end;
if (flag_verbose); disp(sprintf(' %% If you were to set a maf_hi_threshold as well (e.g., 0.49) ;')); end;
if (flag_verbose); disp(sprintf(' %% then the sequence would read ''p46q49''. ;')); end;
if (flag_verbose); disp(sprintf(' %% These sequences are generated via: ;')); end;
if (flag_verbose); disp(sprintf(' %% xxxcluster_fromdisk_uADZSZDA_xfix_gen_ver16.m ;')); end;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now the actual maf thresholding is done by: ;')); end;
if (flag_verbose); disp(sprintf(' %% mc_from_bim_ext_ver5.m ;')); end;
if (flag_verbose); disp(sprintf(' %% which actually thresholds based on *both* the maf, ;')); end;
if (flag_verbose); disp(sprintf(' %% as well as the allele-frequency. ;')); end;
if (flag_verbose); disp(sprintf(' %% Specifically, allele-combinations are used only if: ;')); end;
if (flag_verbose); disp(sprintf(' %% Their corresponding maf is above the cutoff, ;')); end;
if (flag_verbose); disp(sprintf(' %% and their corresponding frequency is above the cutoff^2. ;')); end;
if (flag_verbose); disp(sprintf(' %% This may not be exactly what you want! ;')); end;

flag_continue=input(' %% continue? 1 or 0[default]'); if isempty(flag_continue); flag_continue=0; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_continue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now, if desired, we can load the traces from each run. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note that we use parameter_maflo instead of parameter.')); end;
trace__ = load_trace__from_dir_ver0(parameter_maflo.dir_out_trace);
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now we display the original trace (red) against the label-shuffled traces (black). ;')); end;
[tmp_nlpR,ij_nlpR] = max(trace__.nlpR_s0000_);
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
linewidth_sml = 0.5;
linewidth_big = 3;
markersize_big = 16;
%%%%;
subplot(1,2,1);
hold on;
plot(trace__.niter_s0000_,trace__.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(trace__.niter_s0000_,trace__.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(trace__.niter_s0000_(ij_nlpR),trace__.ZR_s0000_(ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([min(trace__.niter_s0000_),max(trace__.niter_s0000_)]); xlabel('iteration');
ylabel('negative-log-p');
%%%%;
subplot(1,2,2);
hold on;
r_rem_max = max(trace__.r_rem_s0000_);
plot(r_rem_max-trace__.r_rem_s0000_,trace__.ZR_is__,'k-','LineWidth',linewidth_sml);
plot(r_rem_max-trace__.r_rem_s0000_,trace__.ZR_s0000_,'r-','LineWidth',linewidth_big);
plot(r_rem_max-trace__.r_rem_s0000_(ij_nlpR),trace__.ZR_s0000_(ij_nlpR),'ko','MarkerFaceColor','r','MarkerSize',markersize_big);
hold off;
xlim([0,r_rem_max]);
xlabel('patients excluded');
ylabel('negative-log-p');
%%%%;
end;%if flag_disp;

if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now we calculate principal-components. ;')); end;
if (flag_verbose); disp(sprintf(' %% Once again, we must use parameter_maflo. ;')); end;
mx__ = load_mx__from_parameter_ver0(parameter_maflo);
mr_dvx_ = zeros(size(mx__.mr_A_default_full_));
mr_dvx_(find(mx__.mr_A_full_))=2;
mr_dvx_(find(mx__.mr_Z_full_))=1;
index_use_ = efind(mr_dvx_);
parameter_maflo.str_out_xdrop_a_s0000 = sprintf('%s/out_xdrop_a.txt',parameter_maflo.dir_out_s0000);
[ ...
 xdrop_ ...
] = ...
load_out_xdrop_from_str_ver0( ...
parameter_maflo.str_out_xdrop_a_s0000 ...
);
index_rdrop_ = xdrop_.index_rdrop_;
index_cdrop_ = xdrop_.index_cdrop_;
if flag_disp;
figure(1+nf);nf=nf+1;clf;
figmed;fig80s;
end;%if flag_disp;
p_row = 1; p_col = 2; np=0;
markersize_sml = 4;
markersize_med = 8;
markersize_big = 16;
xdrop_ij_rkeep_ = xdrop_.ij_rkeep_;
xdrop_ij_ckeep_ = xdrop_.ij_ckeep_;
%%%%;
[A_p_c_,A_p_0_,AZ_rsum_] = load_A_p_c_from_dir_0(parameter_maflo.dir_out_s0000);
alpha_c_ = A_p_c_ - (1-A_p_c_);
D_c_ = sqrt(1./max(0.01,4.0*A_p_c_.*(1-A_p_c_)));
mx__ = load_mx__from_parameter_ver0(parameter_maflo);
%%%%%%%%;
ni=0;
[ ...
 AZnV_ni0_driver_pd__ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_ni_ver16( ...
 parameter_maflo ...
,ni ...
,xdrop_ ...
,trace__ ...
,mx__ ...
);
if flag_disp;
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_ni0_driver_pd__(1+index_use_,1),AZnV_ni0_driver_pd__(1+index_use_,2),16,mr_dvx_(1+index_use_),'filled','MarkerEdgeColor','k');
axisnotick; title('full data-set (driver)'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;
%%%%%%%%;
ni=ij_nlpR-1;
[ ...
 AZnV_nix_driver_pd__ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_ni_ver16( ...
 parameter_maflo ...
,ni ...
,xdrop_ ...
,trace__ ...
,mx__ ...
);
if flag_disp;
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_nix_driver_pd__(1+index_use_,1),AZnV_nix_driver_pd__(1+index_use_,2),16,mr_dvx_(1+index_use_),'filled','MarkerEdgeColor','k');
axisnotick; title('bicluster-informed (driver)'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_continue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

