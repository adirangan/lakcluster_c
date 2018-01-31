function mc = lakcluster_PGC_getbim_ver0(fname,p_threshold);
fp = fopen(fname,'r');
% taken from bed_to_bim_ver2.m line 434 ; 
% fields are: ; 
%{
  Chromosome
  Marker ID
  Genetic distance
  Physical position
  Allele 1
  Allele 2
  bit type (and vs or) 
  ii_tot (entropy)
  fr2_xxx_tot (sparsity of column)
  frq_mss_tot (missingness of snp)
    %}
bim = textscan(fp,'%d\t%s\t%d\t%d\t%c\t%c\t%s\t%f\t%f\t%f\n');
fclose(fp);
marker_id_ = bim{2};
bit_type_ =  bim{7};
ii_tot_ = bim{8};
fr2_tot_ = bim{9};
mss_tot_ = bim{10};
[C_first_,IA_first_,IC_first_] = unique(marker_id_,'first');
[C_last_,IA_last_,IC_last_] = unique(marker_id_,'last');
fr_and = zeros(length(C_first_),1); fr__or = zeros(length(C_first_),1);
ii_and = zeros(length(C_first_),1); ii__or = zeros(length(C_first_),1);
ms_and = zeros(length(C_first_),1); ms__or = zeros(length(C_first_),1);
AA = zeros(length(C_first_),1);
Aa = zeros(length(C_first_),1);
pp = zeros(length(C_first_),1);
disp_flag=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nc=1:length(C_first_);
C_first = C_first_{nc};
marker_id = marker_id_{IA_first_(nc)};
bit_type = bit_type_{IA_first_(nc)};
ii_tot = ii_tot_(IA_first_(nc));
fr2_tot = fr2_tot_(IA_first_(nc));
mss_tot = mss_tot_(IA_first_(nc));
if strcmp(bit_type,'and'); 
ii_and(nc) = ii_tot; fr_and(nc) = fr2_tot; ms_and(nc) = mss_tot; 
end;%if strcmp(bit_type,'and'); 
if strcmp(bit_type,'or'); 
ii__or(nc) = ii_tot; fr__or(nc) = fr2_tot; ms__or(nc) = mss_tot; 
end;%if strcmp(bit_type,'or'); 
if disp_flag; disp(sprintf(' %% %s-->%s: %s %f',C_first,marker_id,bit_type,fr2_tot)); end;
end;%for nc=1:length(C_first_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nc=1:length(C_last_);
C_last = C_last_{nc};
marker_id = marker_id_{IA_last_(nc)};
bit_type = bit_type_{IA_last_(nc)};
ii_tot = ii_tot_(IA_last_(nc));
fr2_tot = fr2_tot_(IA_last_(nc));
mss_tot = mss_tot_(IA_last_(nc));
if strcmp(bit_type,'and'); 
ii_and(nc) = ii_tot; fr_and(nc) = fr2_tot; ms_and(nc) = mss_tot; 
end;%if strcmp(bit_type,'and'); 
if strcmp(bit_type,'or'); 
ii__or(nc) = ii_tot; fr__or(nc) = fr2_tot; ms__or(nc) = mss_tot; 
end;%if strcmp(bit_type,'or'); 
if disp_flag; disp(sprintf(' %% %s-->%s: %s %f',C_last,marker_id,bit_type,fr2_tot)); end;
end;%for nc=1:length(C_last_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA = fr_and ; Aa = fr__or - fr_and ; pp = AA + 0.5*Aa ; 
mc = ones(length(marker_id_),1);
tmp = find(pp>1-p_threshold);
mc(IA_first_(tmp))=0; mc(IA_last_(tmp))=0;
