clear;
fname = '/data/rangan/dir_bcc/dir_seek_072916/data/human/go_bp_iea.txt';
fcheck(fname);
fid=fopen(fname);
n_line = 0; flag_continue=1;
flag_continue=1;
while (flag_continue); tmp_ = fgetl(fid); if tmp_~=-1; n_line = n_line+1; else flag_continue=0; end; end;
fclose(fid);
disp(sprintf(' %% reading %s: found %d lines',fname,n_line));
n_pathway = floor(n_line/2);
pathway_name_ = cell(n_pathway,1);
pathway_EZid__ = cell(n_pathway,1);
pathway_size_ = zeros(n_pathway,1);
fid=fopen(fname);
for npathway=1:n_pathway;
tmp_ = fgetl(fid);
pathway_name_{npathway} = tmp_;
tmp_ = fgetl(fid);
pathway_EZid__{npathway} = str2num(tmp_);
pathway_size_(npathway) = length(pathway_EZid__{npathway});
end;%for npathway=1:n_pathway;
fclose(fid);
%%%%%%%%;
flag_cap=0;
if flag_cap;
pathway_cap_ = zeros(n_pathway,n_pathway);
for np1=1:n_pathway; 
pathway_cap_(np1,np1) = 1.0;
if (mod(np1,100)==0); disp(sprintf(' %% np1 %d/%d',np1,n_pathway)); end;
for np2=np1+1:n_pathway;
tmp = length(intersect(pathway_EZid__{np1},pathway_EZid__{np2}))/max(1,min(pathway_size_(np1),pathway_size_(np2)));
pathway_cap_(np1,np2) = tmp; pathway_cap_(np2,np1) = tmp;
end;%for np2=1:n_pathway;
end;%for np1=1:n_pathway; 
end;%if flag_cap;
%%%%%%%%;
% Now load annot file created with: ;
% lisa_annot_to_EZ_ver0('/data/rangan/dir_bcc/dir_PGC_20180304/dir_misc/bip_all_bim_uniq_38_annotat_35_10.genes.annot');
% This file contains: ; 
% EZ__ a cell-array associating snps with genes. ;
% EZ__{nsnp} gives the list of (string) EZids associated with the snp with (string) rsid listed in rs_sort_{nsnp}. ;
%%%%%%%%;
fname_mat = '/data/rangan/dir_bcc/dir_PGC_20180304/dir_misc/bip_all_bim_uniq_38_annotat_35_10.genes.annot.mat'; fcheck(fname_mat); load(fname_mat);
%%%%%%%%;
% Now loading xdrop file as well as bim file. ;
% Note that bim fields are: ;
%  1: chromosome number. ;
%  2: rsid. ;
%  3: genetic distance. ;
%  4: physical distance. ;
%  5: allele 1. ;
%  6: allele 2. ;
%  7: allele type. ;
%  8: total relative entropy (across nonexcluded patients). ;
%  9: frequency of allele type (sparsity of column). ;
% 10: frequency of missing data for snp. ;
% 11: minor-allele-frequency for that snp. ;
%%%%%%%%;
verbose=1; cl_num = 4; n_maf = 5; n_cov = 2; flag_dex_vs_lak = 'dex'; if (strcmp(flag_dex_vs_lak,'lak')); gamma = 0.001; else; gamma = 0.004; end; mc_string = ''; B_MLT=34;n_mds=20;
flag_reverse=0; n_shuffle=0; n_scramble=0;
lisa_setprefix_ver2;
lisa_setnames_ver2;
lisa_xdropextract_ver2; %<-- entries in cdrop_a run from 1 to length(mc_A)=A_ncols;
lisa_loadmdsfambim_ver2;
lisa_mxcheck_ver2;
lisa_studyindex_ver2;
lisa_traceextract_ver2; %<-- trace_{1} contains the trace associated with xdrop_a. ;
%%%%%%%%;
% Now extracting unique snp-rsids. ;
% Note that, while most snps in cdrop_a are listed three times, ;
% a few are listed only 2 times. ;
% This is because we exclude allelic-combinations with a frequency ;
% lower than maf_cutoff^2 (usually 0.25^2 = 0.0625). ;
%%%%%%%%;
flag_test=0;
if flag_test;
rs_drop_a_ = bim_{2}(cdrop_a_);
[~,tmp0_,tmp1_] = unique(rs_drop_a_);
tmp2_ = sort(tmp1_);
tmp3_ = diff(tmp2_);
tmp4_ = diff(find(diff(tmp3_)>0));
tmp_example = 1+min(find(tmp4_==2));
disp(sprintf(' %% rs_drop_a_(find(tmp1_==tmp_example)):')); disp(rs_drop_a_(find(tmp1_==tmp_example)));
disp(sprintf(' %% bim_{2}(cdrop_a_(find(tmp1_==tmp_example))): chromosome number '));
disp(bim_{2}(cdrop_a_(find(tmp1_==tmp_example))));
disp(sprintf(' %% bim_{4}(cdrop_a_(find(tmp1_==tmp_example))): physical distance (base pair) '));
disp(bim_{4}(cdrop_a_(find(tmp1_==tmp_example))));
disp(sprintf(' %% bim_{7}(cdrop_a_(find(tmp1_==tmp_example))): allele type '));
disp(bim_{7}(cdrop_a_(find(tmp1_==tmp_example))));
disp(sprintf(' %% bim_{9}(cdrop_a_(find(tmp1_==tmp_example))): allele frequency '));
disp(bim_{9}(cdrop_a_(find(tmp1_==tmp_example))));
disp(sprintf(' %% bim_{11}(cdrop_a_(find(tmp1_==tmp_example))): snp frequency '));
disp(bim_{11}(cdrop_a_(find(tmp1_==tmp_example))));
disp(sprintf(' %% bim_{11}(cdrop_a_(find(tmp1_==tmp_example)))^2: min(snp_frequency^2,(1-snp_frequency)^2) '));
disp(min(bim_{11}(cdrop_a_(find(tmp1_==tmp_example))).^2,(1-bim_{11}(cdrop_a_(find(tmp1_==tmp_example)))).^2));
tmp_ij_ = find(strcmp(bim_{2},rs_drop_a_(min(find(tmp1_==tmp_example)))));
disp(sprintf(' %% bim_{2}(tmp_ij_): chromosome number '));
disp(bim_{2}(tmp_ij_));
disp(sprintf(' %% bim_{4}(tmp_ij_): physical distance (base pair) '));
disp(bim_{4}(tmp_ij_));
disp(sprintf(' %% bim_{7}(tmp_ij_): allele type '));
disp(bim_{7}(tmp_ij_));
disp(sprintf(' %% bim_{9}(tmp_ij_): allele frequency '));
disp(bim_{9}(tmp_ij_));
disp(sprintf(' %% bim_{11}(tmp_ij_): snp frequency '));
disp(bim_{11}(tmp_ij_));
disp(sprintf(' %% bim_{11}(tmp_ij_)^2: min(snp_frequency^2,(1-snp_frequency)^2) '));
disp(min(bim_{11}(tmp_ij_).^2,(1-bim_{11}(tmp_ij_)).^2));
end;%flag_test=1;
%%%%%%%%;
% snp_drop_rs_ contains the names of the snps listed in cdrop_a_. ;
% cap_snp_uni_ is the minimal set of snp-names which is in the intersection of both rs_sort_ and snp_drop_rs_. ;
% snp_drop_xref_ is a sparse matrix that can be used to cross-reference the entries of cdrop_a_ with the unique snp-names stored in snp_drop_rs_uni_. ;
% i.e., : given a list of unique indices snp_drop_rs_uni_(tmp_ij1_), ;
% we have tmp_ij2_ = find(sum(snp_drop_xref_(:,tmp_ij1_),2)>0) returns the corresponding indices from within cdrop_a_. ;
% That is to say, unique(snp_drop_rs_(tmp_ij2_)) should match unique(snp_drop_rs_uni_(tmp_ij1_)). ;
%%%%%%%%;
snp_drop_rs_ = bim_{2}(cdrop_a_); %<-- rs index for each snp in bim_. ;
snp_drop_at_ = bim_{7}(cdrop_a_); %<-- allele type (and, xor, nor) for each snp in bim_. ;
[snp_drop_rs_uni_,snp_drop_rs_get_,snp_drop_rs_put_] = unique(snp_drop_rs_);
snp_drop_xref_ = sparse(1:length(snp_drop_rs_put_),snp_drop_rs_put_,1,length(snp_drop_rs_),length(snp_drop_rs_uni_));
flag_test=0;
if flag_test;
tmp_ij1_ = randperm(length(snp_drop_rs_uni_)); tmp_ij1_ = tmp_ij1_(1:floor(end/3));
tmp_ij2_ = find(sum(snp_drop_xref_(:,tmp_ij1_),2)>0);
assert(length(intersect(unique(snp_drop_rs_(tmp_ij2_)),unique(snp_drop_rs_uni_(tmp_ij1_))))==length(tmp_ij1_));
end;%if flag_test;
[cap_snp_uni_,cap_snp_drop_get_,cap_snp_sort_get_] = intersect(snp_drop_rs_uni_,rs_sort_,'stable');
cap_snp_sort_get_xref_ = sparse(1:length(cap_snp_sort_get_),cap_snp_sort_get_,1,length(cap_snp_sort_get_),length(rs_sort_)); %<-- this can be used to cross-reference indices referencing snps in rs_sort_ with those referencing entries in cap_snp_sort_get_. ;
flag_test=1;
if flag_test;
n_found = 0;
for nsnp1=1:length(rs_sort_);
nsnp2 = find(cap_snp_sort_get_xref_(:,nsnp1));
if (~isempty(nsnp2));
assert(cap_snp_sort_get_(nsnp2)==nsnp1);
n_found = n_found+1;
end;%if (~isempty(nsnp2));
end;%for nsnp=1:length(rs_sort_);
disp(sprintf(' %% testing cap_snp_sort_get_xref_; n_found %d/%d',n_found,length(cap_snp_uni_)));
end;%if flag_test;
%%%%%%%%;
% EZ3__ is a sparse matrix with rows indexing rs_sort_, columns indexing the unique EZids, ;
% and sparse entries indicating which snps are associated with which genes. ;
%%%%%%%%;
EZ2__ = cell(n_rs_sort,1);
for nEZ=1:n_rs_sort;
tmp_ = zeros(ln_(nEZ),1);
for nln=1:ln_(nEZ);
tmp_(nln) = str2num(EZ__{nEZ}{nln});
EZ2__{nEZ} = tmp_;
end;%for nln=1:ln_(nEZ);
end;%for nEZ=1:n_rs_sort;
%%%%%%%%;
EZ_EZmax = 0; for nEZ=1:n_rs_sort; EZ_EZmax = max(EZ_EZmax,max(EZ2__{nEZ})); end;
n_ln = sum(ln_);
EZ_snp_sort_ = zeros(n_ln,1);
EZ_EZid_ = zeros(n_ln,1);
n_tmp = 0;
for nEZ=1:n_rs_sort;
nln = ln_(nEZ);
tmp_ij_ = n_tmp + (1:nln);
EZ_snp_sort_(tmp_ij_) = nEZ;
EZ_EZid_(tmp_ij_) = EZ2__{nEZ};
n_tmp = n_tmp + nln;
end;%for nEZ=1:n_rs_sort;
%%%%%%%%;
[EZ_EZid_uni_,EZ_EZid_get_,EZ_EZid_put_] = unique(EZ_EZid_);
EZ3__ = sparse(EZ_snp_sort_,EZ_EZid_put_,1); %<-- Now EZ3__ can be used to cross-reference snps and genes. ;
% i.e., if a list of EZids from EZ_EZid_ is given by EZ_EZid_uni_(tmp2_), then tmp3_ = find(sum(EZ3__(:,tmp2_),2)) can be used to find the snps associated with those EZids. ;
% More specifically, rs_sort_{tmp3_} lists the (string) rsids such that each of the EZ2__{tmp3_(nl)} lists a set of EZids such that at least one is within the original list from EZ_EZid_uni_(tmp2_). ;
flag_test=0;
if flag_test;
tmp2_ = randperm(length(EZ_EZid_uni_)); tmp2_ = tmp2_(1:floor(end/30));
tmp3_ = find(sum(EZ3__(:,tmp2_),2)); 
for nl=1:length(tmp3_);
assert(length(intersect(EZ2__{tmp3_(nl)},EZ_EZid_uni_(tmp2_)))>0);
end;%for nl=1:length(tmp3_);
end;%if flag_test;
%%%%%%%%;
% now check;
%%%%%%%%;
h_s2_ = hist(sum(EZ3__,2),1:44);
h_s1_ = hist(sum(EZ3__,1),1:1024); %<-- Seems like, on average, each gene is referenced by 35-40 snps. ;
%%%%%%%%;
% Now determine which snps are associated with each pathway. ;
%%%%%%%%;
pathway_EZid_ = zeros(sum(pathway_size_),1);
tmp=0;
for npathway=1:n_pathway;
pathway_EZid_(tmp + (1:pathway_size_(npathway))) = pathway_EZid__{npathway};
tmp = tmp+pathway_size_(npathway);
end;%for npathway=1:n_pathway;
pathway_EZmax = max(pathway_EZid_);
cup_EZmax = max(EZ_EZmax,pathway_EZmax);
[pathway_EZid_uni_,pathway_EZid_get_,pathway_EZid_put_] = unique(pathway_EZid_);
[cap_EZid_uni_,cap_EZid_EZ_get_,cap_EZid_pathway_get_] = intersect(EZ_EZid_uni_,pathway_EZid_uni_,'stable');
cap_EZid_EZ_put_ = sparse(cap_EZid_uni_,1,cap_EZid_EZ_get_,cup_EZmax,1); %<-- this can be used to cross-reference EZids. i.e., EZ_EZid_uni_(cap_EZid_EZ_put_(cap_EZid_uni_(nEZ))) == cap_EZid_uni_(nEZ);
%plot(cap_EZid_uni_ - EZ_EZid_uni_(cap_EZid_EZ_put_(cap_EZid_uni_)),'.'); %<-- should be all zeros. ;
cap_EZid_pathway_put_ = sparse(cap_EZid_uni_,1,cap_EZid_pathway_get_,cup_EZmax,1); %<-- this can be used to cross-refence EZids. i.e,. pathway_EZid_uni_(cap_EZid_pathway_put_(cap_EZid_uni_(nEZ))) == cap_EZid_uni_(nEZ);
%plot(cap_EZid_uni_ - pathway_EZid_uni_(cap_EZid_pathway_put_(cap_EZid_uni_)),'.'); %<-- should be all zeros. ;
%%%%%%%%;

%%%%%%%%%%%%%%%%;
% link pathways with snps in cdrop_a_ by taking the union overall snp-sets for each gene.; 
%%%%%%%%%%%%%%%%;
pathway_snp_cdrop_a_get__ = cell(n_pathway,1); %<-- this will hold the snp-indices (as entries of cdrop_a_) associated with the snps in each pathway. ; i.e., bim_{2}(cdrop_a_(pathway_snp_cdrop_a_get__{npathway})) will hold the (string) rsids associated with pathway pathway_name_{npathway}. ;
for npathway=1:n_pathway;
tmp_pathway_name = pathway_name_{npathway};
tmp_pathway_EZid_ = pathway_EZid__{npathway};
tmp_pathway_size = pathway_size_(npathway);
disp(sprintf(' %% npathway %d/%d name %s size %d',npathway,n_pathway,tmp_pathway_name,tmp_pathway_size));
tmp0_ = find(cap_EZid_pathway_put_(tmp_pathway_EZid_)>0);
tmp0_ = pathway_EZid_uni_(cap_EZid_pathway_put_(tmp_pathway_EZid_(tmp0_))); %<-- now tmp0_ stores only the EZid within pathway_EZid_uni_. ;
tmp1_ = find(cap_EZid_EZ_put_(tmp0_)>0); %<-- now tmp1_ indicates which of the the EZid within pathway_EZid_uni_ are also within EZ_EZid_uni_. ;
tmp2_ = cap_EZid_EZ_put_(tmp0_(tmp1_)); %<-- now tmp2_ allows us to reference the EZid within pathway_EZid_uni_ that are also within EZ_EZid_uni. ; That is, EZ_EZid_uni_(tmp2_) will list the actual EZids in this intersection. ;
tmp3_ = EZ_EZid_uni_(tmp2_); %<-- now tmp3_ stores only the EZid within EZ_EZid_uni_ and pathway_EZid_uni_. ;
tmp4_ = find(sum(EZ3__(:,tmp2_),2)); %<-- now tmp4_ indexes the snp ids from rs_sort_ associated with the genes in the pathway. ;
flag_check=0;
if flag_check;
disp(sprintf(' %% checking tmp4_ length %d',length(tmp4_)));
for nl=1:length(tmp4_);
assert(length(intersect(EZ2__{tmp4_(nl)},tmp3_))>0);
end;%for nl=1:length(tmp4_);
end;%if flag_check;
%{
[~,tmp5_,tmp6_] = intersect(cap_snp_sort_get_,tmp4_,'stable'); %<-- tmp4_(tmp6_) contains the rs_sort_ indices associated with snps in cap_snp_sort_uni_. ;
tmp7b_ = cap_snp_drop_get_(tmp5_); %<-- now snp_drop_rs_uni_(tmp7b_) should equal rs_sort_(tmp4_(tmp6_));
flag_check=0;
if flag_check;
disp(sprintf(' %% checking tmp6_ length %d',length(tmp6_)));
for nl=1:length(tmp6_);
assert(strcmp(snp_drop_rs_uni_(tmp7b_(nl)),rs_sort_(tmp4_(tmp6_(nl)))));
end;%for nl=1:length(tmp6_);
end;%if flag_check;
 %}
tmp7_ = cap_snp_drop_get_(find(sum(cap_snp_sort_get_xref_(:,tmp4_),2)));%<-- now snp_drop_rs_uni_(tmp7_) should equal rs_sort_(tmp4_(tmp6_)), where [~,~,tmp6_] = intersect(cap_snp_sort_get_,tmp4_,'stable');
flag_check=0;
if flag_check;
[~,~,tmp6_] = intersect(cap_snp_sort_get_,tmp4_,'stable');
assert(sum(strcmp(snp_drop_rs_uni_(tmp7_),rs_sort_(tmp4_(tmp6_))))==length(tmp7_));
end;%if flag_check;
%plot (tmp7_-tmp7b_,'.'); %<-- this should be all zeros. ;
flag_check=0;
if flag_check;
disp(sprintf(' %% found tmp7_ of length %d',length(tmp7_)));
end;%if flag_check;
tmp8_ = find(sum(snp_drop_xref_(:,tmp7_),2)>0); %<-- now snp_drop_rs_(tmp8_) should contain only names from snp_drop_rs_uni_(tmp7_);
flag_check=0;
if flag_check;
disp(sprintf(' %% checking tmp8_ length %d',length(tmp8_)));
assert(length(intersect(snp_drop_rs_(tmp8_),snp_drop_rs_uni_(tmp7_)))==length(tmp7_));
assert(length(union(snp_drop_rs_(tmp8_),snp_drop_rs_uni_(tmp7_)))==length(tmp7_));
end;%if flag_check;
tmp9_ = cdrop_a_(tmp8_); %<-- these are the snp column-indices listed in cdrop_a which correspond to the pathway. ;
pathway_snp_cdrop_a_get__{npathway} = tmp8_; %<-- we store only the actual rank-order within cdrop_a for each of the snps in the pathway. ;
end;%for npathway=1:n_pathway;
%%%%%%%%;
% a-posteriori check. ;
%%%%%%%%;
flag_check=0;
if flag_check;
npathway = 1+floor(n_pathway*rand());
disp(sprintf(' %% npathway %d name %s size %d',npathway,pathway_name_{npathway},pathway_size_(npathway)));
tmp_rs_ = unique(bim_{2}(cdrop_a_(pathway_snp_cdrop_a_get__{npathway})));
disp(sprintf(' %% snps %d (%d unique): ',length(pathway_snp_cdrop_a_get__{npathway}),length(tmp_rs_)));
for nl=1:length(tmp_rs_);
tmp_ij = find(strcmp(rs_sort_,tmp_rs_{nl}));
%disp(sprintf(' %% nl %d/%d snp %s, EZid intersection %d',nl,length(tmp_rs_),tmp_rs_{nl},length(intersect(EZ2__{tmp_ij},pathway_EZid__{npathway}))));
assert(length(intersect(EZ2__{tmp_ij},pathway_EZid__{npathway}))>0);
end;%for nl=1:length(tmp_rs_);
end;%if flag_check;
%%%%%%%%;
% Create mask (of indices along cdrop_a) ;
%%%%%%%%;
pathway_snp_cdrop_a_mc_ = zeros(size(cdrop_a_));
pathway_snp_mc_ = zeros(n_snp,1);
for npathway=1:n_pathway;
pathway_snp_size_(npathway) = length(pathway_snp_cdrop_a_get__{npathway});
pathway_snp_cdrop_a_mc_(pathway_snp_cdrop_a_get__{npathway}) = pathway_snp_cdrop_a_mc_(pathway_snp_cdrop_a_get__{npathway}) + 1;
pathway_snp_mc_(cdrop_a_(pathway_snp_cdrop_a_get__{npathway})) = 1;
end;%for npathway=1:n_pathway;

%%%%%%%%;
% EZ4__ is a sparse matrix with rows indexing cdrop_a_, columns indexing the unique EZids, ;
% and sparse entries indicating which entries of cdrop_a_ are associated with which genes. ;
% EZ4_XX__ is a sparse matrix with rows indexing cdrop_a_, columns indexing the unique EZids, ;
% and sparse entries indicating which entry of cdrop_a_ is associated with the XX-percentile for each gene. ;
% Note that the actual rsid taken from snp_drop_rs_{find(EZ4_XX__(:,ng_EZ))} will not necessarily correspond to the gene indexed by EZ_EZid_uni_(ng_EZ). ;
% This is because the indices of find(EZ4_XX__(:,ng_EZ)) are designed to be percentiles referencing the distribution of indices within cdrop_a_ associated with EZ_EZid_uni_(ng_EZ). ;
% Thus, find(EZ4_01__(:,ng_EZ)) and find(EZ4_99__(:,ng_EZ)) will correspond to the first and last entries of cdrop_a_ associated with EZ_EZid_uni_(ng_EZ), ;
% but the intermediate find(EZ4_25__(:,ng_EZ)) through find(EZ4_75__(:,ng_EZ)) may not. ;
%%%%%%%%;
n_total = 0; n_total_XX = 0; clear tmp10_ tmp11_ tmp12_ tmp13_ tmp10_01_ tmp10_25_ tmp10_50_ tmp10_75_ tmp10_99_ tmp13_XX_;
for ng0=1:length(cap_EZid_uni_);
ng_pathway = cap_EZid_pathway_get_(ng0); ng_EZ = cap_EZid_EZ_get_(ng0);
tmp_EZ_EZid_ = EZ_EZid_uni_(ng_EZ); %<-- this is merely a single EZid, which should already be within EZ_EZid_uni_. ;
tmp_pathway_EZid_ = pathway_EZid_uni_(ng_pathway); %<-- this is merely a single EZid, which should already be within pathway_EZid_uni_. ;
assert(tmp_EZ_EZid_==tmp_pathway_EZid_);
tmp0_ = find(cap_EZid_pathway_put_(tmp_pathway_EZid_)>0); %<-- not necessary here. ;
tmp0_ = cap_EZid_pathway_put_(tmp_pathway_EZid_(tmp0_)); %<-- this should be nothing more than ng_pathway. ;
tmp0_ = pathway_EZid_uni_(tmp0_); %<-- now tmp0_ stores only the EZid within pathway_EZid_uni_. ;
tmp1_ = find(cap_EZid_EZ_put_(tmp0_)>0); %<-- now tmp1_ indicates which of the the EZid within pathway_EZid_uni_ are also within EZ_EZid_uni_. ;
tmp2_ = cap_EZid_EZ_put_(tmp0_(tmp1_)); %<-- now tmp2_ allows us to reference the EZid within pathway_EZid_uni_ that are also within EZ_EZid_uni. ; That is, EZ_EZid_uni_(tmp2_) will list the actual EZids in this intersection. ;
tmp3_ = EZ_EZid_uni_(tmp2_); %<-- now tmp3_ stores only the EZid within EZ_EZid_uni_ and pathway_EZid_uni_. ;
tmp4_ = find(sum(EZ3__(:,tmp2_),2)); %<-- now tmp4_ indexes the snp ids from rs_sort_ associated with the genes in the pathway. ;
flag_check=0;
if flag_check;
disp(sprintf(' %% checking tmp4_ length %d',length(tmp4_)));
for nl=1:length(tmp4_);
assert(length(intersect(EZ2__{tmp4_(nl)},tmp3_))>0);
end;%for nl=1:length(tmp4_);
end;%if flag_check;
%{
[~,tmp5_,tmp6_] = intersect(cap_snp_sort_get_,tmp4_,'stable'); %<-- tmp4_(tmp6_) contains the rs_sort_ indices associated with snps in cap_snp_sort_uni_. ;
tmp7b_ = cap_snp_drop_get_(tmp5_); %<-- now snp_drop_rs_uni_(tmp7b_) should equal rs_sort_(tmp4_(tmp6_));
flag_check=0;
if flag_check;
disp(sprintf(' %% checking tmp6_ length %d',length(tmp6_)));
for nl=1:length(tmp6_);
assert(strcmp(snp_drop_rs_uni_(tmp7b_(nl)),rs_sort_(tmp4_(tmp6_(nl)))));
end;%for nl=1:length(tmp6_);
end;%if flag_check;
 %}
tmp7_ = cap_snp_drop_get_(find(sum(cap_snp_sort_get_xref_(:,tmp4_),2)));%<-- now snp_drop_rs_uni_(tmp7_) should equal rs_sort_(tmp4_(tmp6_)), where [~,~,tmp6_] = intersect(cap_snp_sort_get_,tmp4_,'stable');
flag_check=0;
if flag_check;
[~,~,tmp6_] = intersect(cap_snp_sort_get_,tmp4_,'stable');
assert(sum(strcmp(snp_drop_rs_uni_(tmp7_),rs_sort_(tmp4_(tmp6_))))==length(tmp7_));
end;%if flag_check;
%plot (tmp7_-tmp7b_,'.'); %<-- this should be all zeros. ;
tmp8_ = find(sum(snp_drop_xref_(:,tmp7_),2)>0); %<-- now snp_drop_rs_(tmp8_) should contain only names from snp_drop_rs_uni_(tmp7_);
flag_check=0;
if flag_check;
disp(sprintf(' %% checking tmp8_ length %d',length(tmp8_)));
assert(length(intersect(snp_drop_rs_(tmp8_),snp_drop_rs_uni_(tmp7_)))==length(tmp7_));
assert(length(union(snp_drop_rs_(tmp8_),snp_drop_rs_uni_(tmp7_)))==length(tmp7_));
end;%if flag_check;
tmp9_ = cdrop_a_(tmp8_); %<-- these are the snp column-indices listed in cdrop_a which correspond to the genes in the pathway. ;
tmp10_(n_total + (1:length(tmp8_))) = tmp8_; %<-- we store only the actual rank-order within cdrop_a for each of the snps associated with the genes in the pathway. ;
tmp11_(n_total + (1:length(tmp8_))) = ng_pathway; %<-- indexing pathway_EZid_uni_. ;
tmp12_(n_total + (1:length(tmp8_))) = tmp_pathway_EZid_; %<-- actual EZid . ;
tmp13_(n_total + (1:length(tmp8_))) = ng_EZ; %<-- index such that EZ_EZid_uni_(ng_EZ)==tmp_pathway_EZid_.;
n_total = n_total + length(tmp8_);
if length(tmp8_)>0;
tmp_c01_ = max(1,min(length(cdrop_a_),round(prctile(tmp8_,  0))));
tmp_c25_ = max(1,min(length(cdrop_a_),round(prctile(tmp8_, 25))));
tmp_c50_ = max(1,min(length(cdrop_a_),round(prctile(tmp8_, 50))));
tmp_c75_ = max(1,min(length(cdrop_a_),round(prctile(tmp8_, 75))));
tmp_c99_ = max(1,min(length(cdrop_a_),round(prctile(tmp8_,100))));
tmp10_01_(n_total_XX + 1) = tmp_c01_;
tmp10_25_(n_total_XX + 1) = tmp_c25_;
tmp10_50_(n_total_XX + 1) = tmp_c50_;
tmp10_75_(n_total_XX + 1) = tmp_c75_;
tmp10_99_(n_total_XX + 1) = tmp_c99_;
tmp13_XX_(n_total_XX + 1) = ng_EZ;
n_total_XX = n_total_XX + 1;
end;%if length(tmp8_)>0;
end;%for ng0=1:length(cap_EZid_uni_);
%%%%%%%%;
EZ4__ = sparse(tmp10_,tmp13_,1,length(cdrop_a_),length(EZ_EZid_uni_)); %<-- Now EZ4__ can be used to cross-reference entries of cdrop_a_ and genes. ;
EZ4_01__ = sparse(tmp10_01_,tmp13_XX_,1,length(cdrop_a_),length(EZ_EZid_uni_));
EZ4_25__ = sparse(tmp10_25_,tmp13_XX_,1,length(cdrop_a_),length(EZ_EZid_uni_));
EZ4_50__ = sparse(tmp10_50_,tmp13_XX_,1,length(cdrop_a_),length(EZ_EZid_uni_));
EZ4_75__ = sparse(tmp10_75_,tmp13_XX_,1,length(cdrop_a_),length(EZ_EZid_uni_));
EZ4_99__ = sparse(tmp10_99_,tmp13_XX_,1,length(cdrop_a_),length(EZ_EZid_uni_));
% i.e., if a list of EZids from EZ_EZid_ is given by EZ_EZid_uni_(tmp2_), then tmp3_ = find(sum(EZ4__(:,tmp2_),2)) can be used to find the entries of cdrop_a_ associated with those EZids. ;
% More specifically, bim_{2}(cdrop_a_(tmp3_)) == snp_drop_rs_(tmp3_) lists the (string) rsids such that, for tmp4_ = find(strcmp(rs_sort_,snp_drop_rs_(tmp3_))), each of the EZ2__{tmp4_(nl)} lists a set of EZids such that at least one is within the original list from EZ_EZid_uni_(tmp2_). ;
flag_test=1;
if flag_test;
tmp2_ = randperm(length(EZ_EZid_uni_)); tmp2_ = tmp2_(1:20);
disp(sprintf(' %% selecting %d genes from within EZ_EZid_uni_. ',length(tmp2_)));
for nl2=1:length(tmp2_);
tmp3_ = find(sum(EZ4__(:,tmp2_(nl2)),2));
if (length(tmp3_)>0);
assert(find(sum(EZ4_99__(:,tmp2_(nl2)),2))==max(tmp3_));
assert(find(sum(EZ4_01__(:,tmp2_(nl2)),2))==min(tmp3_));
end;%if (length(tmp3_)>0);
if (length(tmp3_)==0);
assert(isempty(find(sum(EZ4_99__(:,tmp2_(nl2)),2))));
assert(isempty(find(sum(EZ4_01__(:,tmp2_(nl2)),2))));
end;%if (length(tmp3_)==0);
end;%for nl2=1:length(tmp2_);
tmp3_ = find(sum(EZ4__(:,tmp2_),2)); 
disp(sprintf(' %% found %d corresponding entries of cdrop_a_. ',length(tmp3_)));
for nl1=1:length(tmp3_);
tmp4_ = find(strcmp(rs_sort_,snp_drop_rs_(tmp3_(nl1))));
disp(sprintf(' %% %% entry %d/%d, found %d snps',nl1,length(tmp3_),length(tmp4_)));
for nl2=1:length(tmp4_);
assert(length(intersect(EZ2__{tmp4_(nl2)},EZ_EZid_uni_(tmp2_)))>0);
end;%for nl2=1:length(tmp4_);
end;%for nl1=1:length(tmp3_);
end;%if flag_test;

%%%%%%%%%%%%%%%%;
% link pathways with snps in cdrop_a_ by taking only one of the snps from each snp-set for each gene.; 
%%%%%%%%%%%%%%%%;
%pathway_EZid_cdrop_a_XX__ = cell(n_pathway,1); %<-- this will hold the snp-indices (as entries of cdrop_a_) associated with EZids in each pathway. ; i.e., bim_{2}(cdrop_a_(pathway_EZid_cdrop_a_XX__{npathway})) will hold the (string) rsids associated with first XX-percentile occurrences of the genes in pathway pathway_name_{npathway}. ;
pathway_EZid_01_cdrop_a_get__ = cell(n_pathway,1); %<-- actually the 00th percentile (i.e., first occurrence). ;
pathway_EZid_25_cdrop_a_get__ = cell(n_pathway,1);
pathway_EZid_50_cdrop_a_get__ = cell(n_pathway,1);
pathway_EZid_75_cdrop_a_get__ = cell(n_pathway,1);
pathway_EZid_99_cdrop_a_get__ = cell(n_pathway,1); %<-- actually the 100th percentile (i.e., last occurrence). ;
for npathway=1:n_pathway;
tmp_pathway_name = pathway_name_{npathway};
tmp_pathway_EZid_ = pathway_EZid__{npathway};
tmp_pathway_size = pathway_size_(npathway);
disp(sprintf(' %% npathway %d/%d name %s size %d',npathway,n_pathway,tmp_pathway_name,tmp_pathway_size));
tmp0_ = find(cap_EZid_pathway_put_(tmp_pathway_EZid_)>0);
tmp0_ = pathway_EZid_uni_(cap_EZid_pathway_put_(tmp_pathway_EZid_(tmp0_))); %<-- now tmp0_ stores only the EZid within pathway_EZid_uni_. ;
tmp1_ = find(cap_EZid_EZ_put_(tmp0_)>0); %<-- now tmp1_ indicates which of the the EZid within pathway_EZid_uni_ are also within EZ_EZid_uni_. ;
tmp2_ = cap_EZid_EZ_put_(tmp0_(tmp1_)); %<-- now tmp2_ allows us to reference the EZid within pathway_EZid_uni_ that are also within EZ_EZid_uni. ; That is, EZ_EZid_uni_(tmp2_) will list the actual EZids in this intersection. ;
tmp3_ = EZ_EZid_uni_(tmp2_); %<-- now tmp3_ stores only the EZid within EZ_EZid_uni_ and pathway_EZid_uni_. ;
tmp4_c01_ = find(sum(EZ4_01__(:,tmp2_),2)); %<-- now tmp4_ indexes the first entry from cdrop_a_ associated with each of the genes in tmp3_(nl2). ;
tmp4_c25_ = find(sum(EZ4_25__(:,tmp2_),2));
tmp4_c50_ = find(sum(EZ4_50__(:,tmp2_),2));
tmp4_c75_ = find(sum(EZ4_75__(:,tmp2_),2));
tmp4_c99_ = find(sum(EZ4_99__(:,tmp2_),2)); %<-- now tmp4_ indexes the last entry from cdrop_a_ associated with each of the genes in tmp3_(nl2). ;
pathway_EZid_01_cdrop_a_get__{npathway} = tmp4_c01_;
pathway_EZid_25_cdrop_a_get__{npathway} = tmp4_c25_;
pathway_EZid_50_cdrop_a_get__{npathway} = tmp4_c50_;
pathway_EZid_75_cdrop_a_get__{npathway} = tmp4_c75_;
pathway_EZid_99_cdrop_a_get__{npathway} = tmp4_c99_;
end;%for npathway=1:n_pathway;
%%%%%%%%;
% a-posteriori check. ;
%%%%%%%%;
flag_check=0;
if flag_check;
npathway = 1+floor(n_pathway*rand());
disp(sprintf(' %% npathway %d name %s size %d',npathway,pathway_name_{npathway},pathway_size_(npathway)));
tmp_rs_01_ = unique(bim_{2}(cdrop_a_(pathway_EZid_01_cdrop_a_get__{npathway})));
tmp_rs_25_ = unique(bim_{2}(cdrop_a_(pathway_EZid_25_cdrop_a_get__{npathway})));
tmp_rs_50_ = unique(bim_{2}(cdrop_a_(pathway_EZid_50_cdrop_a_get__{npathway})));
tmp_rs_75_ = unique(bim_{2}(cdrop_a_(pathway_EZid_75_cdrop_a_get__{npathway})));
tmp_rs_99_ = unique(bim_{2}(cdrop_a_(pathway_EZid_99_cdrop_a_get__{npathway})));
disp(sprintf(' %% snps %d (%d unique): ',length(pathway_EZid_01_cdrop_a_get__{npathway}),length(tmp_rs_01_)));
for nl=1:length(tmp_rs_01_);
disp(sprintf(' %% checking tmp_rs_(%d)',nl));
tmp_ij_01_ = find(strcmp(rs_sort_,tmp_rs_01_{nl})); assert(length(intersect(EZ2__{tmp_ij_01_},pathway_EZid__{npathway}))>0);
%tmp_ij_25_ = find(strcmp(rs_sort_,tmp_rs_25_{nl})); assert(length(intersect(EZ2__{tmp_ij_25_},pathway_EZid__{npathway}))>0);
%tmp_ij_50_ = find(strcmp(rs_sort_,tmp_rs_50_{nl})); assert(length(intersect(EZ2__{tmp_ij_50_},pathway_EZid__{npathway}))>0);
%tmp_ij_75_ = find(strcmp(rs_sort_,tmp_rs_75_{nl})); assert(length(intersect(EZ2__{tmp_ij_75_},pathway_EZid__{npathway}))>0);
tmp_ij_99_ = find(strcmp(rs_sort_,tmp_rs_99_{nl})); assert(length(intersect(EZ2__{tmp_ij_99_},pathway_EZid__{npathway}))>0);
end;%for nl=1:length(tmp_rs_01_);
end;%if flag_check;
%%%%%%%%;
% Create masks (of indices along cdrop_a) ;
%%%%%%%%;
pathway_EZid_01_cdrop_a_mc_ = zeros(size(cdrop_a_));
pathway_EZid_25_cdrop_a_mc_ = zeros(size(cdrop_a_));
pathway_EZid_50_cdrop_a_mc_ = zeros(size(cdrop_a_));
pathway_EZid_75_cdrop_a_mc_ = zeros(size(cdrop_a_));
pathway_EZid_99_cdrop_a_mc_ = zeros(size(cdrop_a_));
pathway_EZid_01_mc_ = zeros(n_snp,1);
pathway_EZid_25_mc_ = zeros(n_snp,1);
pathway_EZid_50_mc_ = zeros(n_snp,1);
pathway_EZid_75_mc_ = zeros(n_snp,1);
pathway_EZid_99_mc_ = zeros(n_snp,1);
for npathway=1:n_pathway;
pathway_EZid_01_size_(npathway) = length(pathway_EZid_01_cdrop_a_get__{npathway});
pathway_EZid_25_size_(npathway) = length(pathway_EZid_25_cdrop_a_get__{npathway});
pathway_EZid_50_size_(npathway) = length(pathway_EZid_50_cdrop_a_get__{npathway});
pathway_EZid_75_size_(npathway) = length(pathway_EZid_75_cdrop_a_get__{npathway});
pathway_EZid_99_size_(npathway) = length(pathway_EZid_99_cdrop_a_get__{npathway});
pathway_EZid_01_cdrop_a_mc_(pathway_EZid_01_cdrop_a_get__{npathway}) = pathway_EZid_01_cdrop_a_mc_(pathway_EZid_01_cdrop_a_get__{npathway}) + 1;
pathway_EZid_25_cdrop_a_mc_(pathway_EZid_25_cdrop_a_get__{npathway}) = pathway_EZid_25_cdrop_a_mc_(pathway_EZid_25_cdrop_a_get__{npathway}) + 1;
pathway_EZid_50_cdrop_a_mc_(pathway_EZid_50_cdrop_a_get__{npathway}) = pathway_EZid_50_cdrop_a_mc_(pathway_EZid_50_cdrop_a_get__{npathway}) + 1;
pathway_EZid_75_cdrop_a_mc_(pathway_EZid_75_cdrop_a_get__{npathway}) = pathway_EZid_75_cdrop_a_mc_(pathway_EZid_75_cdrop_a_get__{npathway}) + 1;
pathway_EZid_99_cdrop_a_mc_(pathway_EZid_99_cdrop_a_get__{npathway}) = pathway_EZid_99_cdrop_a_mc_(pathway_EZid_99_cdrop_a_get__{npathway}) + 1;
pathway_EZid_01_mc_(cdrop_a_(pathway_EZid_01_cdrop_a_get__{npathway})) = 1;
pathway_EZid_25_mc_(cdrop_a_(pathway_EZid_25_cdrop_a_get__{npathway})) = 1;
pathway_EZid_50_mc_(cdrop_a_(pathway_EZid_50_cdrop_a_get__{npathway})) = 1;
pathway_EZid_75_mc_(cdrop_a_(pathway_EZid_75_cdrop_a_get__{npathway})) = 1;
pathway_EZid_99_mc_(cdrop_a_(pathway_EZid_99_cdrop_a_get__{npathway})) = 1;
end;%for npathway=1:n_pathway;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now calculate enrichment for all pathways. ;
% This abandons chebfun and just uses a gaussian approximation. ;
%%%%%%%%;
flag_plot=0;
log_pval_99__ = cell(n_pathway,1);
k_99__ = cell(n_pathway,1); K_99_ = zeros(n_pathway,1); M_sub_99_ = zeros(n_pathway,1); M_all_99_ = zeros(n_pathway,1);
m_sub_99__ = cell(n_pathway,1);
m_all_99__ = cell(n_pathway,1);
if flag_plot; hold on;
stairs(niter_,log(0.05)*ones(size(niter_)),'k:');
stairs(niter_,log(0.001)*ones(size(niter_)),'k:');
hold off;
end;%if flag_plot; hold on;
tmp_N_all = length(cdrop_a_); tmp_N_sub = sum((pathway_EZid_99_cdrop_a_mc_>0)); 
[pathway_EZid_99_size_uni_,pathway_EZid_99_size_get_,pathway_EZid_99_size_put_] = unique(pathway_EZid_99_size_);
for nl=length(pathway_EZid_99_size_uni_):-1:1;%for nl=1:length(pathway_EZid_99_size_uni_);
pathway_EZid_99_size = pathway_EZid_99_size_uni_(nl);
tmp_pij_ = find(pathway_EZid_99_size_==pathway_EZid_99_size);
disp(sprintf(' %% nl %d/%d: pathway_EZid_99_size %d, found %d',nl,length(pathway_EZid_99_size_uni_),pathway_EZid_99_size,length(tmp_pij_)));
tmp_K = pathway_EZid_99_size; %tmp_K = length(pathway_EZid_99_cdrop_a_get__{npathway}); 
tmp_M_all = tmp_N_all-tmp_K; tmp_M_sub = tmp_N_sub-tmp_K;
tmp_n_sub_ = cumsum(pathway_EZid_99_cdrop_a_mc_>0,'reverse');
tmp_n_all_ = flip(transpose(1:length(cdrop_a_)));
range_KM_sub_ = [ [tmp_K/2] + tmp_K/2*0.95*[-1,+1] , [tmp_M_sub/2] + tmp_M_sub/2*0.95*[-1,+1] ]; %range_KM_sub_ = [0,tmp_K,0,tmp_M_sub];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nl2=1:length(tmp_pij_);
npathway = tmp_pij_(nl2);
%%%%%%%%%%%%%%%%;
disp(sprintf(' %% %% %% npathway %d name %s EZid size %d snp_cdrop size %d',npathway,pathway_name_{npathway},length(pathway_EZid__{npathway}),length(pathway_EZid_99_cdrop_a_get__{npathway})));
tmp_k_mc_ = zeros(length(cdrop_a_),1); tmp_k_mc_(pathway_EZid_99_cdrop_a_get__{npathway}) = 1;
tmp_k_ = cumsum(tmp_k_mc_,'reverse'); assert(min(tmp_k_)>=0 & max(tmp_k_)<=tmp_K);
tmp_m_sub_ = tmp_n_sub_ - tmp_k_; assert(min(tmp_m_sub_)>=0 & max(tmp_m_sub_)<=tmp_M_sub);
tmp_m_all_ = tmp_n_all_ - tmp_k_; assert(min(tmp_m_all_)>=0 & max(tmp_m_all_)<=tmp_M_all);
tmp_pk_ = tmp_k_/tmp_K; tmp_qk_ = 1 - tmp_pk_;
tmp_pm_sub_ = tmp_m_sub_/tmp_M_sub; tmp_qm_sub_ = 1 - tmp_pm_sub_; 
tmp_pm_all_ = tmp_m_all_/tmp_M_all; tmp_qm_all_ = 1 - tmp_pm_all_; 
%tmp_sg_sub_ = sqrt(1./(tmp_M_sub .* tmp_pm_sub_ .* tmp_qm_sub_) + 1./(tmp_K .* tmp_pk_ .* tmp_qk_));
%tmp_sg_all_ = sqrt(1./(tmp_M_all .* tmp_pm_all_ .* tmp_qm_all_) + 1./(tmp_K .* tmp_pk_ .* tmp_qk_));
%%%%%%%%;
log_pval_99_ = log(0.5)*ones(length(c_rem_),1);
c_ij_ = length(cdrop_a_) - c_rem_ + 1;
tmp_ij_ = find(tmp_pk_(c_ij_)>tmp_pm_sub_(c_ij_));
if (length(tmp_ij_)>0); 
tmp_k_get_ = tmp_k_(c_ij_(tmp_ij_)); tmp_m_sub_get_ = tmp_m_sub_(c_ij_(tmp_ij_));
sh_ = snchoosek(tmp_M_sub,tmp_m_sub_get_) + snchoosek(tmp_K,tmp_k_get_) - snchoosek(tmp_N_sub,tmp_m_sub_get_+tmp_k_get_) ;
dtsh_ = log(tmp_K-tmp_k_get_) - log(tmp_k_get_) - log(tmp_M_sub-tmp_m_sub_get_) + log(tmp_m_sub_get_);
dttsh_ = - ( tmp_M_sub./(tmp_m_sub_get_ .* (tmp_M_sub - tmp_m_sub_get_)) + tmp_K./(tmp_k_get_ .* (tmp_K - tmp_k_get_)) );
tmp_x_ = sqrt(-dttsh_/2) .* ( 0  + dtsh_ ./ (2*dttsh_/2) ) ;
tmp_z_ = erfcln(tmp_x_);
log_eshc_ = sh_ + (-dtsh_ .* dtsh_ ./ (2*dttsh_)) - 0.5*log(-dttsh_/2) + log(sqrt(pi)/2) + tmp_z_;
log_pval_99_(tmp_ij_) = log_eshc_;
end; %if (length(tmp_ij_)>0); 
log_pval_99__{npathway} = log_pval_99_;
K_99_(npathway) = tmp_K; M_sub_99_(npathway) = tmp_M_sub; M_all_99_(npathway) = tmp_M_all;
k_99__{npathway} = tmp_k_(c_ij_); 
m_sub_99__{npathway} = tmp_m_sub_(c_ij_);
m_all_99__{npathway} = tmp_m_all_(c_ij_);
if flag_plot;
hold on; stairs(niter_,log_pval_99_,'k-','LineWidth',0.5); hold off; 
ylim([-15,0]); xlim([1,length(c_rem_)]);
tic_ij_ = max(1,0:25:length(c_rem_));
%set(gca,'XTick',niter_(tic_ij_),'XTickLabel',c_rem_(tic_ij_)); xtickangle(90); xlabel('columns remaining');
set(gca,'XTick',niter_(tic_ij_),'XTickLabel',1+niter_(tic_ij_)); xtickangle(90); xlabel('iteration');
ylabel('log(pvalue)');
drawnow();
end;%if flag_plot;
%%%%%%%%%%%%%%%%;
end;%for nl2=1:length(tmp_pij_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nl=1:length(pathway_EZid_99_size_uni_);
if flag_plot;
stairs(niter_,log(0.05)*ones(size(niter_)),'k:');
stairs(niter_,log(0.001)*ones(size(niter_)),'k:');
hold off;
ylim([-15,0]); xlim([1,length(c_rem_)]);
tic_ij_ = max(1,0:25:length(c_rem_));
%set(gca,'XTick',niter_(tic_ij_),'XTickLabel',c_rem_(tic_ij_)); xtickangle(90); xlabel('columns remaining');
set(gca,'XTick',niter_(tic_ij_),'XTickLabel',1+niter_(tic_ij_)); xtickangle(90); xlabel('iteration');
ylabel('log(pvalue)');
end;% if flag_plot;

%%%%%%%%;
% plot aggregate results. ;
%%%%%%%%;
colormap(colormap_beach());
imagesc( log(1+hist2d_0(repmat(niter_,1,n_pathway),cell2mat(transpose(log_pval_99__)),length(niter_),1501,[0,niter_(end)],[-15,0])) , [0,4] );
colorbar;
tic_ij_ = max(1,0:25:length(c_rem_));
%set(gca,'XTick',niter_(tic_ij_),'XTickLabel',c_rem_(tic_ij_)); xtickangle(90); xlabel('columns remaining');
set(gca,'XTick',niter_(tic_ij_),'XTickLabel',1+niter_(tic_ij_)); xtickangle(90); xlabel('iteration');
tmp_l_ = linspace(-15,0,1501); [~,~,tic_ij_] = intersect(-15:0,tmp_l_,'stable');
set(gca,'YTick',tic_ij_,'YTickLabel',-15:0); ylabel('log(pvalue)');

% list top enrichments. ;
ni = 276;
tmp_pval_99_ = cell2mat(transpose(log_pval_99__));
[~,tmp_ij_] = sort(tmp_pval_99_(ni,:));
for np=1:45;
npathway = tmp_ij_(np);
disp(sprintf(' %% ni %.4d: np %.4d/%.4d: npathway %.4d: (%.4d/%.4d/%.4d) log10_pval %+8.4f \t %s ',ni,np,45,npathway,k_99__{npathway}(ni),m_sub_99__{npathway}(ni),K_99_(npathway),log_pval_99__{npathway}(ni)/log(10),pathway_name_{npathway}));
end;%for np=1:45;

% list enrichment for a particular pathway. ;
ni = 276; npathway = find(pathway_name_{
tmp_pval_99_ = cell2mat(transpose(log_pval_99__));
[~,tmp_ij_] = sort(tmp_pval_99_(ni,:));
for np=1:45;
npathway = tmp_ij_(np);
disp(sprintf(' %% ni %.4d: np %.4d/%.4d: npathway %.4d: (%.4d/%.4d/%.4d) log10_pval %+8.4f \t %s ',ni,np,45,npathway,k_99__{npathway}(ni),m_sub_99__{npathway}(ni),K_99_(npathway),log_pval_99__{npathway}(ni)/log(10),pathway_name_{npathway}));
end;%for np=1:45;
