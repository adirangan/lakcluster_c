function mc = mc_from_bim_ext_ver5(fname,p_lo_threshold,p_hi_threshold);
% creates a column-mask mc which is 0 for all alleles such that: ;
% either ;
% a) the minor-allele-frequency of that snp is below p_lo_threshold, ;
% or ;
% b) the minor-allele-frequency of that snp is above p_hi_threshold, ;
% or ;
% c) the frequency of allele type is below p_lo_threshold.^2 ;
% or ;
% d) the frequency of allele type is above (1-p_lo_threshold).^2 ;

if nargin<3;
p_hi_threshold = 1.0;
end;%if nargin<3;

verbose=0;
if verbose;
disp(sprintf(' %% Bim: ; '));
disp(sprintf(' %% This has the same structure as a typical bim file, ; '));
disp(sprintf(' %%\t Chromosome code (either an integer, or ''X''/''Y''/''XY''/''MT''; ''0'' indicates unknown) or name ;'));
disp(sprintf(' %%\t Variant identifier ;'));
disp(sprintf(' %%\t Position in morgans or centimorgans (safe to use dummy value of ''0'') ;'));
disp(sprintf(' %%\t Base-pair coordinate (1-based; limited to 2^31-2) ;'));
disp(sprintf(' %%\t Allele 1 (corresponding to clear bits in .bed; usually minor) ;'));
disp(sprintf(' %%\t Allele 2 (corresponding to set bits in .bed; usually major) ;'));
disp(sprintf(' %% with five extra fields. ; '));
disp(sprintf(' %% field 6+1: allele type. ; '));
disp(sprintf(' %% field 6+2: total relative entropy across (nonexcluded) patients. ; '));
disp(sprintf(' %% field 6+3: frequency of allele type (i.e., sparsity of column). ; '));
disp(sprintf(' %% field 6+4: frequency of missing data for that snp. ; '));
disp(sprintf(' %% field 6+5: minor-allele-frequency of that snp. ; '));
disp(sprintf(' %% Moreover, the allele-ordering should be sorted in terms of column-sparsity. ; '));

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
fname ...
);

ij_tmp_ = find(bim_maf_>0.5); 
bim_maf_(ij_tmp_) = 1 - bim_maf_(ij_tmp_);
mc = ones(n_snp,1);
p_min = min(p_lo_threshold,p_hi_threshold);
p_max = max(p_lo_threshold,p_hi_threshold);
ij_tmp_ = find( (bim_maf_ < p_min) | (bim_maf_ > p_max) | (bim_frq_ < p_min.^2) | (bim_frq_ > (1-p_min).^2) );
mc(ij_tmp_)=0;
