function label_num_ = label_str_to_num_0(label_str_);
% converts string labels to numeric labels, ;
% assuming that the strings are themselves number. ;
% e.g., '32', '-1', etc.
n_u = numel(label_str_);
u_label_ = unique(label_str_);
n_label = length(u_label_);
n_label_ = zeros(n_label,1);
for nlabel = 1:n_label;
n_label_(nlabel) = length(find(strcmp(label_str_,u_label_(nlabel))));
end;%for nlabel = 1:n_label;
label_num_ = zeros(n_u,1);
for nlabel=1:n_label;
tmp_ij_ = find(strcmp(label_str_,u_label_(nlabel)));
if (~isempty(str2num(u_label_{nlabel})));
tmp_num = str2num(u_label_{nlabel});
if (~isempty(tmp_num)); label_num_(tmp_ij_) = tmp_num; end;
if ( isempty(tmp_num)); label_num_(tmp_ij_) = 0; end;
end;%if (~isempty(str2num(u_label_{nlabel})));
end;%for nlabel=1:n_label;

