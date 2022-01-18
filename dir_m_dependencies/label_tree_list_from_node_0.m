function list_node_ = label_tree_list_from_node_0(label_tree,flag_include_leaf);
flag_isleaf = isempty(label_tree.childA) & isempty(label_tree.childB);
%%%%%%%%;
if ( (~flag_isleaf) | (flag_include_leaf) );
list_node_ = cell(1,1);
list_node_{1} = label_tree;
end;%if ( (~flag_isleaf) | (flag_include_leaf) );
if ( (flag_isleaf) & (~flag_include_leaf) );
list_node_ = cell(0,1);
end;%if ( (flag_isleaf) & (~flag_include_leaf) );
%%%%%%%%;
if ( ~flag_isleaf )
list_node_A_ = label_tree_list_from_node_0(label_tree.childA,flag_include_leaf);
list_node_B_ = label_tree_list_from_node_0(label_tree.childB,flag_include_leaf);
list_node_ = [ list_node_ ; list_node_A_ ; list_node_B_ ];
end;%if ( ~flag_isleaf );
