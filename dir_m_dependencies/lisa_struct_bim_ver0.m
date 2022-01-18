function l = lisa_struct_bim_ver0(l);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% load bim file. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
l.fname_bim = sprintf('%s/%s_bim.ex2',l.dir__in,l.string_prefix); fcheck(l.fname_bim);
[l.n_snp,l.bim__id_,l.bim_al1_,l.bim_al2_,l.bim_alt_,l.bim_ent_,l.bim_frq_,l.bim_maf_,l.bim_name_,l.bim_] = lisa_bimex2_ver0(l.fname_bim);
