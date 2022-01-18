function ...
[ ...
 output_label_pid___ ...
,lpv__ ...
,n_cluster_ ...
,lpv_rearrange_0__ ...
,n_cluster_rearrange_0_ ...
,lpv_rearrange_1__ ...
,n_cluster_rearrange_1_ ...
,lpv_rearrange_2a__ ...
,n_cluster_rearrange_2a_ ...
,lpv_rearrange_2b__ ...
,n_cluster_rearrange_2b_ ...
] = ...
label_to_tree_enrichment_0( ...
 label_tree ...
 ,label_A__ ...
 ,nlpnex_threshold_ ...
 ,A_n_ ...
);

verbose=0;
if (verbose); disp(sprintf(' %% [entering label_to_tree_enrichment_0]')); end;

na=0;
if (nargin<1+na); label_tree = []; end; na=na+1;
if (nargin<1+na); label_A__ = []; end; na=na+1;
if (nargin<1+na); nlpnex_threshold_ = []; end; na=na+1;
if (nargin<1+na); A_n_ = []; end; na=na+1;

if isempty(nlpnex_threshold_); nlpnex_threshold_ = 3:0.5:27; end;
n_nlpnex_threshold = numel(nlpnex_threshold_);
n_label_A = numel(label_A__);

ouput_label_pid___ = cell(n_nlpnex_threshold,1);
lpv__ = zeros(n_nlpnex_threshold,n_label_A);
n_cluster_ = zeros(n_nlpnex_threshold,1);
lpv_rearrange_0__ = zeros(n_nlpnex_threshold,n_label_A);
n_cluster_rearrange_0_ = zeros(n_nlpnex_threshold,1);
lpv_rearrange_1__ = zeros(n_nlpnex_threshold,n_label_A);
n_cluster_rearrange_1_ = zeros(n_nlpnex_threshold,1);
%%%%%%%%%%%%%%%%;
for nnlpnex_threshold=0:n_nlpnex_threshold-1;
nlpnex_threshold = nlpnex_threshold_(1+nnlpnex_threshold);
if (verbose); disp(sprintf(' %% nnlpnex_threshold %d/%d nlpnex_threshold %0.2f',nnlpnex_threshold,n_nlpnex_threshold,nlpnex_threshold)); end;
tmp_output_label__ = label_tree_output_label_from_nlpnex_0(label_tree,nlpnex_threshold,[]);
output_label_pid___{1+nnlpnex_threshold} = tmp_output_label__;
label_B_ = label_str_to_enum_1(tmp_output_label__);
n_cluster_(1+nnlpnex_threshold) = numel(unique(label_B_));
for nlabel_A=0:n_label_A-1;
tmp_lpv = label_to_label_enrichment_quad_4(label_A__{1+nlabel_A},label_B_);
lpv__(1+nnlpnex_threshold,1+nlabel_A) = tmp_lpv;
end;%for nlabel_A=0:n_label_A-1;
%%%%;
if ~isempty(A_n_);
%%%%;
if (nargout>=4);
tmp_t = tic();
label_B_rearrange_0_ = label_rearrange_0(label_B_,A_n_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% label_rearrange_0 time %0.2fs',tmp_t)); end;
n_cluster_rearrange_0_(1+nnlpnex_threshold) = numel(unique(label_B_rearrange_0_));
for nlabel_A=0:n_label_A-1;
tmp_lpv = label_to_label_enrichment_quad_4(label_A__{1+nlabel_A},label_B_rearrange_0_);
lpv_rearrange_0__(1+nnlpnex_threshold,1+nlabel_A) = tmp_lpv;
end;%for nlabel_A=0:n_label_A-1;
end;%if (nargout>=4);
%%;
if (nargout>=6);
tmp_t = tic();
label_B_rearrange_1_ = label_rearrange_1(label_B_,A_n_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% label_rearrange_1 time %0.2fs',tmp_t)); end;
n_cluster_rearrange_1_(1+nnlpnex_threshold) = numel(unique(label_B_rearrange_1_));
for nlabel_A=0:n_label_A-1;
tmp_lpv = label_to_label_enrichment_quad_4(label_A__{1+nlabel_A},label_B_rearrange_1_);
lpv_rearrange_1__(1+nnlpnex_threshold,1+nlabel_A) = tmp_lpv;
end;%for nlabel_A=0:n_label_A-1;
end;%if (nargout>=6);
%%;
if (nargout>=8);
tmp_t = tic();
label_B_rearrange_2a_ = label_rearrange_2(label_B_,A_n_,0);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% label_rearrange_2a time %0.2fs',tmp_t)); end;
n_cluster_rearrange_2a_(1+nnlpnex_threshold) = numel(unique(label_B_rearrange_2a_));
for nlabel_A=0:n_label_A-1;
tmp_lpv = label_to_label_enrichment_quad_4(label_A__{1+nlabel_A},label_B_rearrange_2a_);
lpv_rearrange_2a__(1+nnlpnex_threshold,1+nlabel_A) = tmp_lpv;
end;%for nlabel_A=0:n_label_A-1;
end;%if (nargout>=8);
%%;
if (nargout>=10);
tmp_t = tic();
label_B_rearrange_2b_ = label_rearrange_2(label_B_,A_n_,1);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% label_rearrange_2b time %0.2fs',tmp_t)); end;
n_cluster_rearrange_2b_(1+nnlpnex_threshold) = numel(unique(label_B_rearrange_2b_));
for nlabel_A=0:n_label_A-1;
tmp_lpv = label_to_label_enrichment_quad_4(label_A__{1+nlabel_A},label_B_rearrange_2b_);
lpv_rearrange_2b__(1+nnlpnex_threshold,1+nlabel_A) = tmp_lpv;
end;%for nlabel_A=0:n_label_A-1;
end;%if (nargout>=10);
%%%%;
end;%if ~isempty(A_n_);
%%%%;
end;%for nnlpnex_threshold=0:n_nlpnex_threshold-1;
%%%%%%%%%%%%%%%%;

if (verbose); disp(sprintf(' %% [finished label_to_tree_enrichment_0]')); end;
