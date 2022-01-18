function ontology = ontology_struct_popuqu_ver0(ontology,annot);

eu_ = annot.eu_;
eu_max = annot.eu_max;

%%%%%%%%;
% Now determine which snps are associated with each pathway. ;
% po__ is a cell array: po__{npathway} contains the EZids associated with pathway_name_{npathway}. ;
% po_ is an array of (nonunique) EZids created from po__. ;
% pu_ is an array of unique EZids in po_;
% qu_ is the intersection of eu_ (from the annotation) and pu_ (from the ontology). ;
%%%%%%%%;
po__ = ontology.pathway_EZid__;
po_ = zeros(sum(ontology.pathway_size_),1);
tmp=0;
for npathway=1:ontology.n_pathway;
po_(tmp + (1:ontology.pathway_size_(npathway))) = po__{npathway};
tmp = tmp+ontology.pathway_size_(npathway);
end;%for npathway=1:ontology.n_pathway;
ontology.pathway_EZmax = max(po_);
ontology.cup_EZmax = max(eu_max,ontology.pathway_EZmax);
[pu_,po_to_pu_,pu_to_po_] = unique(po_);
[qu_,eu_to_qu_,pu_to_qu_] = intersect(eu_,pu_,'stable');
eu_to_qu_xref_ = sparse(qu_,1,eu_to_qu_,ontology.cup_EZmax,1); %<-- this can be used to cross-reference EZids. i.e., eu_(eu_to_qu_xref_(qu_(nEZ))) == qu_(nEZ);
%plot(qu_ - eu_(eu_to_qu_xref_(qu_)),'.'); %<-- should be all zeros. ;
pu_to_qu_xref_ = sparse(qu_,1,pu_to_qu_,ontology.cup_EZmax,1); %<-- this can be used to cross-refence EZids. i.e,. pu_(pu_to_qu_xref_(qu_(nEZ))) == qu_(nEZ);
%plot(qu_ - pu_(pu_to_qu_xref_(qu_)),'.'); %<-- should be all zeros. ;
%%%%%%%%;

ontology.po__ = po__;
ontology.po_ = po_;
ontology.pu_ = pu_;
ontology.po_to_pu_ = po_to_pu_;
ontology.pu_to_po_ = pu_to_po_;
ontology.qu_ = qu_;
ontology.eu_to_qu_ = eu_to_qu_;
ontology.pu_to_qu_ = pu_to_qu_;
ontology.eu_to_qu_xref_ = eu_to_qu_xref_;
ontology.pu_to_qu_xref_ = pu_to_qu_xref_;
