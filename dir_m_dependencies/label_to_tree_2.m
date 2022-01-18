function ...
label_tree ...
= ...
label_to_tree_2( ...
 depth ...
,p_nex ...
,label_string ...
,item_index_ ...
,output_label_ ...
,nlpbra_num_ ...
,nlpnex_num_ ...
);

%%%%%%%%;
% Assumes that leaves have nlpbra and nlpnex. ;
%%%%%%%%;

if isempty(depth); depth = 0; end;
if isempty(p_nex); p_nex = 0; end;
if isempty(label_string); label_string = ''; end;

if isempty(item_index_); item_index_ = transpose(0:numel(output_label_)-1); end;
if (~isempty(nlpbra_num_));
if (~isnumeric(nlpbra_num_{1})); nlpbra_num_ = cellfun(@str2num,nlpbra_num_,'UniformOutput',0); end;
end;%if (~isempty(nlpbra_num_));
if (~isempty(nlpnex_num_));
if (~isnumeric(nlpnex_num_{1})); nlpnex_num_ = cellfun(@str2num,nlpnex_num_,'UniformOutput',0); end;
end;%if (~isempty(nlpnex_num_));

label_tree = struct('type','label_tree');
label_tree.psplit = [];
label_tree.depth = depth;
label_tree.label_string = sprintf('%s0',label_string);
label_tree.parent = [];
label_tree.item_index_ = item_index_;
label_tree.nlpbra = [];
label_tree.nlpnex = [];
label_tree.childA = [];
label_tree.childB = [];

n_item = numel(item_index_);
assert(numel(output_label_)==n_item);
assert(numel(nlpbra_num_)==n_item);
assert(numel(nlpnex_num_)==n_item);
output_label_start_ = cell(n_item,1);
nlpbra_num_start_ = zeros(n_item,1);
nlpnex_num_start_ = zeros(n_item,1);
for nitem=0:n_item-1;
output_label_start_{1+nitem} = output_label_{1+nitem}(1+depth);
%%%%%%%%;
if (~strcmp(output_label_start_{1+nitem},'0')); %<-- not leaf. ;
nlpbra_num_start_(1+nitem) = nlpbra_num_{1+nitem}(1+depth);
nlpnex_num_start_(1+nitem) = nlpnex_num_{1+nitem}(1+depth);
end;%if (~strcmp(output_label_start_{1+nitem},'0'));
%%%%%%%%;
% Here we assume that leaves list nlpbra and nlpnex. ;
%%%%%%%%;
if ( strcmp(output_label_start_{1+nitem},'0')); %<-- yes leaf. ;
nlpbra_num_start_(1+nitem) = nlpbra_num_{1+nitem}(1+depth);
nlpnex_num_start_(1+nitem) = nlpnex_num_{1+nitem}(1+depth);
end;%if ( strcmp(output_label_start_{1+nitem},'0'));
%%%%%%%%;
end;%for nitem=0:n_item-1;
tmp_index_A_ = efind(strcmp(output_label_start_,'A'));
tmp_index_B_ = efind(strcmp(output_label_start_,'B'));
tmp_index_0_ = efind(strcmp(output_label_start_,'0'));
n_item_A = numel(tmp_index_A_);
n_item_B = numel(tmp_index_B_);
n_item_0 = numel(tmp_index_0_);
if (n_item_A> 0 & n_item_B> 0); assert(n_item_0==0); end;
if (n_item_0> 0); %<-- yes leaf. ;
assert( (n_item_A==0) & (n_item_B==0) );
assert(numel(strcmp(output_label_,label_tree.label_string))==n_item);
%%%%%%%%;
% Once again we assume that leaves list nlpbra and nlpnex. ;
%%%%%%%%;
label_tree.nlpbra = nlpbra_num_start_(1);
label_tree.nlpnex = nlpnex_num_start_(1);
p_bra = exp(-label_tree.nlpbra);
p_nex = 1 - (1-p_nex)*(1-p_bra);
label_tree.psplit = p_nex;
%%%%%%%%;
end;%if (n_item_0> 0);
if ( (n_item_A> 0) & (n_item_B> 0) ); %<-- not leaf. ;
assert(numel(unique(nlpbra_num_start_))==1);
assert(numel(unique(nlpnex_num_start_))==1);
label_tree.nlpbra = nlpbra_num_start_(1);
label_tree.nlpnex = nlpnex_num_start_(1);
p_bra = exp(-label_tree.nlpbra);
p_nex = 1 - (1-p_nex)*(1-p_bra);
label_tree.psplit = p_nex;
label_tree.childA ...
= ...
label_to_tree_2( ...
 depth+1 ...
,p_nex ...
,sprintf('%sA',label_string) ...
,item_index_(1+tmp_index_A_) ...
,output_label_(1+tmp_index_A_) ...
,nlpbra_num_(1+tmp_index_A_) ...
,nlpnex_num_(1+tmp_index_A_) ...
);
label_tree.childA.parent = label_tree;
label_tree.childB ...
= ...
label_to_tree_2( ...
 depth+1 ...
,p_nex ...
,sprintf('%sB',label_string) ...
,item_index_(1+tmp_index_B_) ...
,output_label_(1+tmp_index_B_) ...
,nlpbra_num_(1+tmp_index_B_) ...
,nlpnex_num_(1+tmp_index_B_) ...
);
label_tree.childB.parent = label_tree;
end;%if ( (n_item_A> 0) & (n_item_B> 0) );

