function output_label__ = label_tree_output_label_from_psplit_0(label_tree,psplit_threshold,output_label__);
psplit_tolerance = 1e-15;
if isempty(label_tree.parent); output_label__ = cell(numel(label_tree.item_index_),1); end;
flag_isleaf = isempty(label_tree.childA) & isempty(label_tree.childB);
flag_branch = ~flag_isleaf & (label_tree.psplit < psplit_threshold+psplit_tolerance);
if ( flag_isleaf | ~flag_branch );
n_item = numel(label_tree.item_index_);
for nitem=0:n_item-1;
item_index = label_tree.item_index_(1+nitem);
output_label__{1+item_index} = label_tree.label_string;
end;%for nitem=0:n_item-1;
end;%if ( flag_isleaf | ~flag_branch );
if ( ~flag_isleaf & flag_branch );
output_label__ = label_tree_output_label_from_psplit_0(label_tree.childA,psplit_threshold,output_label__);
output_label__ = label_tree_output_label_from_psplit_0(label_tree.childB,psplit_threshold,output_label__);
end;%if ( ~flag_isleaf & flag_branch );