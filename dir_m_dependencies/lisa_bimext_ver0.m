function [n_snp,bim_chr_,bim__id_,bim_gdi_,bim_pdi_,bim_al1_,bim_al2_,bim_alt_,bim_ent_,bim_frq_,bim_mss_,bim_maf_,bim_name_,bim_] = lisa_bimext_ver0(fname);
% extract fields from bim.ext file. ;

fcheck(fname);
fid = fopen(fname); 
bim_ = textscan(fid,'%d\t%s\t%d\t%d\t%c\t%c\t%s\t%f\t%f\t%f\t%f\n','headerlines',0);
fclose(fid);

n_snp = length(bim_{1});
assert(n_snp==wc_0(fname));

bim_chr_ = bim_{1}; %<-- chromosome;
bim__id_ = bim_{2}; %<-- rsid;
bim_gdi_ = bim_{3}; %<-- genetic distance;
bim_pdi_ = bim_{4}; %<-- physical distance;
bim_al1_ = bim_{5}; %<-- allele 1;
bim_al2_ = bim_{6}; %<-- allele 2;
bim_alt_ = bim_{7}; %<-- allele type;
bim_ent_ = bim_{8}; %<-- relative entropy;
bim_frq_ = bim_{9}; %<-- sparsity of allele type (i.e., frequency of column);
bim_mss_ = bim_{10}; %<-- fraction of missing data;
bim_maf_ = bim_{11}; %<-- frequency of snp (i.e., minor allele frequency);

if nargout>13;
bim_name_ = cell(n_snp,1); 
for nsnp=1:n_snp; 
bim_name_{nsnp} = sprintf('%s%s%s%s%s%s%s',bim_{2}{nsnp},'&',bim_{5}(nsnp),'&',bim_{6}(nsnp),'&',bim_{7}{nsnp}); 
end;%for nsnp=1:n_snp;
end;%if nargout>13;
