setup_0;

flag_verbose=1;
flag_disp=flag_verbose; nf=0;

if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Testing subsampled version of xxxcluster_fromdisk_uADZSZDA_ver16_dr_0. ;')); end;
if (flag_verbose); disp(sprintf(' %% First running test_stripped_xxxcluster_fromdisk_uADZSZDA_ver16. ;')); end;
dir_trunk = sprintf('%s/dir_test_xxxcluster_fromdisk_ver16',pwd);
test_stripped_xxxcluster_fromdisk_ver16(struct('flag_verbose',0));
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now running xxxcluster_fromdisk_uADZSZDA_ver16_dr_0 on same files. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note that, this time, we return the ''parameter'' structure. ;')); end;
if (flag_verbose); disp(sprintf(' %% We will later alter this structure to run a new analysis ;')); end;
if (flag_verbose); disp(sprintf(' %% on a subset of the patients. ;')); end;
parameter = struct('type','parameter');
parameter.flag_verbose = 0;
parameter.dir_trunk = dir_trunk;
parameter.study_name_of_branch_s_ = {'dir_study00'};
parameter.study_name_without_extension_s_ = {'study00'};
parameter.mds_name_plus_extension_s_ = {'study00_mds_tsv.txt'};
parameter.ent_cutoff = 0.0045022835561375; %<-- test_stripped_xxxcluster_fromdisk_uADZSZDA_ver16;
parameter = xxxcluster_fromdisk_uADZSZDA_ver16_dr_0(parameter);
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now loading the row- and col-masks used in the original run. ;')); end;
if (flag_verbose); disp(sprintf(' %% The row-masks refer to patients: ;')); end;
if (flag_verbose); disp(sprintf(' %% A refers to case-patients, Z to ctrls. ;')); end;
if (flag_verbose); disp(sprintf(' %% The ''full'' infix refers to all the studies combined. ;')); end;
if (flag_verbose); disp(sprintf(' %% Because we have only one study, we can load the full row masks. ;')); end;
if (flag_verbose); disp(sprintf(' %% The col-mask refers to snps (or rather, allele-combinations). ;')); end;
if (flag_verbose); disp(sprintf(' %% The cols are the same across studies, so there is only one col-mask. ;')); end;
fname_mr_A_full = sprintf('%s/%s_mr_A_full.b16',parameter.dir_0in,parameter.str_prefix);
fname_mr_Z_full = sprintf('%s/%s_mr_Z_full.b16',parameter.dir_0in,parameter.str_prefix);
fname_mc_A = sprintf('%s/%s_mc_A.b16',parameter.dir_0in,parameter.str_prefix);
mr_A_ori_ = binary_uncompress(fname_mr_A_full)>0; %<-- thresholded to avoid negative values. ;
mr_Z_ori_ = binary_uncompress(fname_mr_Z_full)>0; %<-- thresholded to avoid negative values. ;
mc_A_ori_ = binary_uncompress(fname_mc_A)>0; %<-- thresholded to avoid negative values. ;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now we can alter these row- and col-masks. ;')); end;
n_patient = numel(mr_A_ori_);
mr_A_alt_ = mr_A_ori_; mr_A_alt_(end-10:2:end)=0; %<-- remove a few cases. ;
mr_Z_alt_ = mr_Z_ori_; mr_Z_alt_(1:2:10)=0; %<-- remove a few ctrls. ;
n_snp = numel(mc_A_ori_);
mc_A_alt_ = mc_A_ori_; mc_A_alt_(1:2:20)=0; %<-- remove a few cols. ;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now we specify a row- and col-mask infix for naming. ;')); end;
if (flag_verbose); disp(sprintf(' %% Then we save the altered masks as files in the input directory. ;')); end;
str_mr_0in = 'rsub';
str_mc_0in = 'csub';
fname_mr_A_alt_full = sprintf('%s/%s_mr_A_%s_full.b16',parameter.dir_0in,parameter.str_prefix,str_mr_0in);
fname_mr_Z_alt_full = sprintf('%s/%s_mr_Z_%s_full.b16',parameter.dir_0in,parameter.str_prefix,str_mr_0in);
fname_mc_A_alt = sprintf('%s/%s_mc_A_%s.b16',parameter.dir_0in,parameter.str_prefix,str_mc_0in);
bitj=16;
binary_compress(bitj,mr_A_alt_>0,fname_mr_A_alt_full);
binary_compress(bitj,mr_Z_alt_>0,fname_mr_Z_alt_full);
binary_compress(bitj,mc_A_alt_>0,fname_mc_A_alt);
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% We do the same for each study. ;')); end;
fname_mr_A_alt_01 = sprintf('%s/%s_mr_A_%s_01.b16',parameter.dir_0in,parameter.str_prefix,str_mr_0in);
fname_mr_Z_alt_01 = sprintf('%s/%s_mr_Z_%s_01.b16',parameter.dir_0in,parameter.str_prefix,str_mr_0in);
bitj=16;
binary_compress(bitj,mr_A_alt_>0,fname_mr_A_alt_01);
binary_compress(bitj,mr_Z_alt_>0,fname_mr_Z_alt_01);
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now we copy the parameter structure, ;')); end;
if (flag_verbose); disp(sprintf(' %% and add in the altered mask infixes ;')); end;
if (flag_verbose); disp(sprintf(' %% (which were previously empty by default). ;')); end;
if (flag_verbose); disp(sprintf(' %% The nshuffle=0 run does not permute the case-ctrl labels. ;')); end;
if (flag_verbose); disp(sprintf(' %% The nshuffle>0 runs operate on label-shuffled data. ;')); end;
if (flag_verbose); disp(sprintf(' %% This label-shuffling permutes only the cases and ctrls ;')); end;
if (flag_verbose); disp(sprintf(' %% included in the altered row-masks. ;')); end;
if (flag_verbose); disp(sprintf(' %% Of course these permutations also respect study-boundaries and mds-sectors ;')); end;
if (flag_verbose); disp(sprintf(' %% if multiple studies and/or mds-components are included. ;')); end;
parameter_alt = parameter;
parameter_alt.str_mr_0in = str_mr_0in;
parameter_alt.str_mc_0in = str_mc_0in;
for nshuffle=0:parameter_alt.n_shuffle-1+1;
parameter_alt.nshuffle = nshuffle;
parameter_alt = xxxcluster_fromdisk_uADZSZDA_ver16(parameter_alt);
end;%for nshuffle=0:parameter_alt.n_shuffle-1+1;
if (flag_verbose); disp(sprintf(' %% ;')); end;

if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now we load the traces from each run. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note that we use parameter_alt instead of parameter.')); end;
trace__ = load_trace__from_dir_ver0(parameter_alt.dir_out_trace);
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
if (flag_verbose); disp(sprintf(' %% Once again, we must use parameter_alt. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note also that we subselect the relevant rows: ')); end;
if (flag_verbose); disp(sprintf(' %% ''index_rsub_'' is limited to the alternate mask. ;')); end;
mr_dvx_ = zeros(n_patient,1);
mr_dvx_(find(mr_A_alt_))=2;
mr_dvx_(find(mr_Z_alt_))=1;
index_rsub_ = efind(mr_dvx_);
parameter_alt.str_out_xdrop_a_s0000 = sprintf('%s/out_xdrop_a.txt',parameter_alt.dir_out_s0000);
[ ...
 xdrop_ ...
] = ...
load_out_xdrop_from_str_ver0( ...
parameter_alt.str_out_xdrop_a_s0000 ...
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
[A_p_c_,A_p_0_,AZ_rsum_] = load_A_p_c_from_dir_0(parameter_alt.dir_out_s0000);
alpha_c_ = A_p_c_ - (1-A_p_c_);
D_c_ = sqrt(1./max(0.01,4.0*A_p_c_.*(1-A_p_c_)));
mx__ = load_mx__from_parameter_ver0(parameter_alt);
%%%%%%%%;
ni=0;
[ ...
 AZnV_ni0_driver_pd__ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_ni_ver16( ...
 parameter_alt ...
,ni ...
,xdrop_ ...
,trace__ ...
,mx__ ...
);
if flag_disp;
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_ni0_driver_pd__(1+index_rsub_,1),AZnV_ni0_driver_pd__(1+index_rsub_,2),16,mr_dvx_(1+index_rsub_),'filled','MarkerEdgeColor','k');
axisnotick; title('full data-set (driver)'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;
%%%%%%%%;
ni=ij_nlpR-1;
[ ...
 AZnV_nix_driver_pd__ ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_D_from_ni_ver16( ...
 parameter_alt ...
,ni ...
,xdrop_ ...
,trace__ ...
,mx__ ...
);
if flag_disp;
subplot(p_row,p_col,1+np);np=np+1;
scatter(AZnV_nix_driver_pd__(1+index_rsub_,1),AZnV_nix_driver_pd__(1+index_rsub_,2),16,mr_dvx_(1+index_rsub_),'filled','MarkerEdgeColor','k');
axisnotick; title('bicluster-informed (driver)'); xlabel('PC1'); ylabel('PC2');
end;%if flag_disp;

