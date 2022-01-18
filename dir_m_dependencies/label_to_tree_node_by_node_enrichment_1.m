function ...
[ ...
 nlpnex_threshold_ ...
,psplit_threshold_ ...
,output_label_pid___ ...
,lpv__ ...
,n_cluster_ ...
,lpv_rearrange_2b__ ...
,n_cluster_rearrange_2b_ ...
] = ...
label_to_tree_node_by_node_enrichment_1( ...
 label_tree ...
 ,label_A__ ...
 ,A_n_ ...
);

verbose=1;
if (verbose); disp(sprintf(' %% [entering label_to_tree_node_by_node_enrichment_1]')); end;

na=0;
if (nargin<1+na); label_tree = []; end; na=na+1;
if (nargin<1+na); label_A__ = []; end; na=na+1;
if (nargin<1+na); A_n_ = []; end; na=na+1;

list_node_ = label_tree_list_from_node_0(label_tree,0);
n_node = numel(list_node_);
nlpnex_ = zeros(n_node,1);
psplit_ = zeros(n_node,1);
for nnode=0:n_node-1;
nlpnex_(1+nnode) = list_node_{1+nnode}.nlpnex;
psplit_(1+nnode) = list_node_{1+nnode}.psplit;
end;%for nnode=0:n_node-1;
[~,index_sort_] = sort(psplit_,'ascend'); index_sort_ = index_sort_-1;
nlpnex_threshold_ = nlpnex_(1+index_sort_);
psplit_threshold_ = psplit_(1+index_sort_);

n_label_A = numel(label_A__);

ouput_label_pid___ = cell(n_node,1);
lpv__ = zeros(n_node,n_label_A);
n_cluster_ = zeros(n_node,1);
lpv_rearrange_2b__ = zeros(n_node,n_label_A);
n_cluster_rearrange_2b_ = zeros(n_node,1);
time_output_label = 0;
time_enrichment = 0;
time_rearrange = 0;
%%%%%%%%%%%%%%%%;
for nnode=0:n_node-1;
psplit = psplit_(1+index_sort_(1+nnode));
if (verbose); disp(sprintf(' %% nnode %d/%d, psplit %0.16f --> nlp %0.2f',nnode,n_node,psplit,-log(psplit))); end;
tmp_t = tic();
tmp_output_label__ = label_tree_output_label_from_psplit_0(label_tree,psplit,[]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% output_label time %0.2fs',tmp_t)); end;
time_output_label = time_output_label + tmp_t;
output_label_pid___{1+nnode} = tmp_output_label__;
label_B_ = label_str_to_enum_1(tmp_output_label__);
n_cluster_(1+nnode) = numel(unique(label_B_));
tmp_t = tic();
for nlabel_A=0:n_label_A-1;
tmp_lpv = label_to_label_enrichment_quad_4(label_A__{1+nlabel_A},label_B_);
lpv__(1+nnode,1+nlabel_A) = tmp_lpv;
end;%for nlabel_A=0:n_label_A-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% enrichment time %0.2fs',tmp_t)); end;
time_enrichment = time_enrichment + tmp_t;
%%%%;
if ~isempty(A_n_);
%%%%;
if (nargout>=4);
tmp_t = tic();
label_B_rearrange_2b_ = label_rearrange_2(label_B_,A_n_,1);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% label_rearrange_2b time %0.2fs',tmp_t)); end;
time_rearrange = time_rearrange + tmp_t;
n_cluster_rearrange_2b_(1+nnode) = numel(unique(label_B_rearrange_2b_));
tmp_t = tic();
for nlabel_A=0:n_label_A-1;
tmp_lpv = label_to_label_enrichment_quad_4(label_A__{1+nlabel_A},label_B_rearrange_2b_);
lpv_rearrange_2b__(1+nnode,1+nlabel_A) = tmp_lpv;
end;%for nlabel_A=0:n_label_A-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% enrichment time %0.2fs',tmp_t)); end;
time_enrichment = time_enrichment + tmp_t;
end;%if (nargout>=4);
%%%%;
end;%if ~isempty(A_n_);
%%%%;
end;%for nnode=0:n_node-1;
%%%%%%%%%%%%%%%%;

if (verbose); disp(sprintf(' %% output_label time %0.2fs',time_output_label)); end;
if (verbose); disp(sprintf(' %% enrichment time %0.2fs',time_enrichment)); end;
if (verbose); disp(sprintf(' %% rearrange time %0.2fs',time_rearrange)); end;

if (verbose); disp(sprintf(' %% [finished label_to_tree_node_by_node_enrichment_1]')); end;
