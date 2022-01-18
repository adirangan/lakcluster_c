function mc = xxxcluster_PGC_getbim_ver2(fname,p_threshold,f_threshold);
% creates a column-mask mc which is 0 for all alleles such that: ;
% either ;
% a) the minor-allele-frequency of that snp is below p_threshold, ;
% or ;
% b) the column-frequency of that allele is below f_threshold. ;
% ;
% Note that we typically expect f_threshold to be close to p_threshold^2. ;

if nargin<3;
f_threshold = min(p_threshold,1-p_threshold)^2;
end;%if nargin<3;

verbose=0;
if verbose;
disp(sprintf('%% Bim: ; '));
disp(sprintf('%% This has the same structure as a typical bim file, ; '));
disp(sprintf('%% with five extra fields. ; '));
disp(sprintf('%% field 6+1: allele type. ; '));
disp(sprintf('%% field 6+2: total relative entropy across (nonexcluded) patients. ; '));
disp(sprintf('%% field 6+3: frequency of allele type (i.e., sparsity of column). ; '));
disp(sprintf('%% field 6+4: frequency of missing data for that snp. ; '));
disp(sprintf('%% field 6+5: minor-allele-frequency of that snp. ; '));
disp(sprintf('%% Moreover, the allele-ordering should be sorted in terms of column-sparsity. ; '));
end;%if verbose;

fp = fopen(fname,'r');
bim_ = textscan(fp,'%d\t%s\t%d\t%d\t%c\t%c\t%s\t%f\t%f\t%f\t%f\n');
fclose(fp);
marker_id_ = bim_{2};
bit_type_ =  bim_{7};
I_opt_tot_ = bim_{8};
fr2_xxx_tot_ = bim_{9};
frq_mss_tot_ = bim_{10};
x_opt_tot_ = bim_{11};
mc = ones(length(marker_id_),1);
tmp = find(x_opt_tot_ < p_threshold | x_opt_tot_ > 1-p_threshold | fr2_xxx_tot_ < f_threshold | fr2_xxx_tot_ > 1-f_threshold);
mc(tmp)=0;
