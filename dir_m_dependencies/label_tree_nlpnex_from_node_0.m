function [nlpnex_,psplit_] = label_tree_nlpnex_from_node_0(label_tree);
nlpnex_ = label_tree.nlpnex; psplit_ = label_tree.psplit;
if ( ~isempty(label_tree.childA) & ~isempty(label_tree.childB) );
[childA_nlpnex_,childA_psplit_] = label_tree_nlpnex_from_node_0(label_tree.childA);
[childB_nlpnex_,childB_psplit_] = label_tree_nlpnex_from_node_0(label_tree.childB);
nlpnex_ = [ nlpnex_ ; childA_nlpnex_ ; childB_nlpnex_ ];
psplit_ = [ psplit_ ; childA_psplit_ ; childB_psplit_ ];
end;%if ( ~isempty(label_tree.childA) & ~isempty(label_tree.childB) );
