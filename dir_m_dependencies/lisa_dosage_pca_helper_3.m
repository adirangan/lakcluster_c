function [A_arm1_xxx_] = lisa_dosage_pca_helper_3(string_alleletype,alleletype_arm1_,bo1_,bu1_,A_arm1_);

co1_ = find(strcmp(alleletype_arm1_,string_alleletype));
[cu1_,co1_to_cu1_,cu1_to_co1_] = unique(bo1_(co1_));
co1_by_cu1_xref_ = sparse(1:length(co1_),cu1_to_co1_,1,length(co1_),length(cu1_));
[du1_,du1_to_bu1_,du1_to_cu1_] = intersect(bu1_,cu1_,'stable');
bu1_by_du1_xref_ = sparse(du1_to_bu1_,1:length(du1_),1,length(bu1_),length(du1_));
cu1_by_du1_xref_ = sparse(du1_to_cu1_,1:length(du1_),1,length(cu1_),length(du1_));
A_arm1_xxx_ = zeros(length(bu1_),size(A_arm1_,2));
for ndu1=1:length(du1_);
nbu1 = find(sum(bu1_by_du1_xref_(:,ndu1),2));
ncu1 = find(sum(cu1_by_du1_xref_(:,ndu1),2));
nco1 = find(sum(co1_by_cu1_xref_(:,ncu1),2));
nbo1 = co1_(nco1);
A_arm1_xxx_(nbu1,:) = A_arm1_(nbo1,:);
end;%for ndu1=1:length(du1_);
clear co1_ cu1_ co1_to_cu1_ cu1_to_co1_ du1_ du1_to_bu1_ du1_to_cu1_ bu1_by_du1_xref_ cu1_by_du1_xref_ ;
