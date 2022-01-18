function [n_snp,bim__id_,bim_al1_,bim_al2_,bim_alt_,bim_ent_,bim_frq_,bim_maf_,bim_name_,bim_] = lisa_bimex2_ver0(fname);
% extract fields from bim.ex2 file. ;

fcheck(fname);
fid = fopen(fname); 
bim_ = textscan(fid,'%s\t%c\t%c\t%s\t%f\t%f\t%f\n','headerlines',0);
fclose(fid);

n_snp = length(bim_{1});
assert(n_snp==wc_0(fname));

bim__id_ = bim_{1}; %<-- rsid; (2) ;
bim_al1_ = bim_{2}; %<-- allele 1; (5) ;
bim_al2_ = bim_{3}; %<-- allele 2; (6) ;
bim_alt_ = bim_{4}; %<-- allele type; (7) ;
bim_ent_ = bim_{5}; %<-- relative entropy; (8) ;
bim_frq_ = bim_{6}; %<-- sparsity of allele type (i.e., frequency of column); (9) ;
bim_maf_ = bim_{7}; %<-- frequency of snp (i.e., minor allele frequency); (11) ;

%%%%%%%%;
% original bimext ;
%%%%%%%%;
%bim_chr_ = bim_{1}; %<-- chromosome; (na) ;
%bim__id_ = bim_{2}; %<-- rsid; (1) ;
%bim_gdi_ = bim_{3}; %<-- genetic distance; (na) ;
%bim_pdi_ = bim_{4}; %<-- physical distance; (na) ;
%bim_al1_ = bim_{5}; %<-- allele 1; (2) ;
%bim_al2_ = bim_{6}; %<-- allele 2; (3) ; 
%bim_alt_ = bim_{7}; %<-- allele type; (4) ;
%bim_ent_ = bim_{8}; %<-- relative entropy; (5) ;
%bim_frq_ = bim_{9}; %<-- sparsity of allele type (i.e., frequency of column); (6) ;
%bim_mss_ = bim_{10}; %<-- fraction of missing data; (na) ;
%bim_maf_ = bim_{11}; %<-- frequency of snp (i.e., minor allele frequency); (7) ;

if nargout>8;
bim_name_ = cell(n_snp,1); 
for nsnp=1:n_snp; 
bim_name_{nsnp} = sprintf('%s%s%s%s%s%s%s',bim_{1}{nsnp},'&',bim_{2}(nsnp),'&',bim_{3}(nsnp),'&',bim_{4}{nsnp}); 
end;%for nsnp=1:n_snp;
end;%if nargout>8;
