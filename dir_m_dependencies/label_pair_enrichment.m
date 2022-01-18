function [lpv,lP_0,cap_,cup_] = label_pair_enrichment(label_A_,label_B_);
% Calculates enrichment p-value for label_A_ and label_B_. ;
% Assumes numeric labels, although labels do not need to be sequential. ;
% Also assumes that each label is binary (i.e., two choices). ;
% lpv = p-value of intersections. ;
% lP_0 = probability of achieving intersections. ;
% cap_ = set of intersections. ;
% cup_ = set of unions. ;
%%%%%%%%;

if (nargin<1);
disp(sprintf(' %% [testing label_pair_enrichment]'));
n_X = 8;
label_A_ = [3;3;3;5;5;5;5;5];
label_B_ = [2;2;2;2;2;2;4;4];
disp(sprintf(' %% label_A_: ')); disp(num2str(transpose(label_A_)));
disp(sprintf(' %% label_B_: ')); disp(num2str(transpose(label_B_)));
[lpv,lP_0,cap_,cup_] = label_pair_enrichment(label_A_,label_B_);
disp(sprintf(' %% lpv %0.2f --> log(%0.2f)',lpv,exp(lpv)));
disp('returning'); return;
end;%if (nargin<1);

%%%%%%%%;
verbose=0;
if (verbose); disp(sprintf(' %% [entering label_pair_enrichment]')); end;
lP_0 = 0;
%%%%%%%%;
n_A = numel(label_A_);
u_label_A_ = unique(label_A_);
n_label_A = length(u_label_A_);
n_label_A_ = zeros(n_label_A,1);
for nlabel_A = 1:n_label_A;
ij_A_{nlabel_A} = find(label_A_ == u_label_A_(nlabel_A));
n_label_A_(nlabel_A) = numel(ij_A_{nlabel_A});
end;%for nlabel_A=1:n_label_A;
%%%%%%%%;
n_B = numel(label_B_);
u_label_B_ = unique(label_B_);
n_label_B = length(u_label_B_);
n_label_B_ = zeros(n_label_B,1);
for nlabel_B = 1:n_label_B;
ij_B_{nlabel_B} = find(label_B_ == u_label_B_(nlabel_B));
n_label_B_(nlabel_B) = numel(ij_B_{nlabel_B});
end;%for nlabel_B=1:n_label_B;
%%%%%%%%;
cap_ = zeros(n_label_A,n_label_B);
cup_ = zeros(n_label_A,n_label_B);
for nlabel_A=1:n_label_A;
for nlabel_B=1:n_label_B;
cap_(nlabel_A,nlabel_B) = numel(intersect(ij_A_{nlabel_A},ij_B_{nlabel_B}));
cup_(nlabel_A,nlabel_B) = numel(union(ij_A_{nlabel_A},ij_B_{nlabel_B}));
end;%for nlabel_B=1:n_label_B;
end;%for nlabel_A=1:n_label_A;
if (verbose); disp(sprintf(' %% cap_: ')); disp(cap_); end;
%%%%%%%%;
if ( (n_label_A<=1) | (n_label_B<=1) );
lpv = 0;
lP_0 = 0;
end;%if ( (n_label_A<=1) | (n_label_B<=1) );
%%%%%%%%;
if ( (n_label_A>=2) & (n_label_B>=2) );
%%%%%%%%;
if (n_label_A~=2); disp(sprintf(' %% Warning, n_label_A %d in label_pair_enrichment',n_label_A)); end;
assert(n_label_A==2);
if (n_label_B~=2); disp(sprintf(' %% Warning, n_label_B %d in label_pair_enrichment',n_label_B)); end;
assert(n_label_B==2);
%%%%%%%%;
if (n_A~=n_B); disp(sprintf(' %% Warning, n_A %d n_B %d in label_pair_enrichment',n_A,n_B)); end;
assert(n_A==n_B);
n_X = n_A;
%%%%%%%%;
% Given a list of specific sets from A, ;
% as well as the cardinality of the sets in B, ;
% the number of ways of finding the precise intersections listed in cap_ is: ;
% lN = sum(gammaln(1+n_label_A_)) - sum(gammaln(1+cap_(:))). ;
% Similarly, the probability of observing cap_ is: ;
% lP = sum(gammaln(1+n_label_A_)) - sum(gammaln(1+cap_(:))) + sum(gammaln(1+n_label_B_)) - gammaln(1+n_X). ;
%%%%%%%%;
num_A_ = zeros(n_A,1);
for nlabel_A=1:n_label_A;
num_A_(ij_A_{nlabel_A}) = nlabel_A;
end;%for nlabel_A=1:n_label_A;
num_B_ = zeros(n_B,1);
for nlabel_B=1:n_label_B;
num_B_(ij_B_{nlabel_B}) = nlabel_B;
end;%for nlabel_B=1:n_label_B;
lP_base = sum(gammaln(1+n_label_A_)) + sum(gammaln(1+n_label_B_)) - gammaln(1+n_X);
cap_0_ = sparse(num_A_,num_B_,1,n_label_A,n_label_B);
tmp_0 = sum(gammaln(1+cap_0_),'all');
lN_0 = +tmp_0;
lP_0 = lP_base - tmp_0;
if (verbose); disp(sprintf(' %% cap_0_: ')); disp(cap_0_); end;
%%%%%%%%;
% Now the p-value is the sum of all probabilities of intersections which are rarer than the precise intersection in cap_. ;
%%%%%%%%;
lp_ = -Inf*ones(1+n_A,1);
for na=0:n_A;
C11 = na;
C12 = n_label_A_(1) - C11;
C21 = n_label_B_(1) - C11;
C22 = n_X - C11 - C12 - C21 ;
tmp_cap_ = [C11 , C12 ; C21 , C22];
if (verbose); disp(sprintf(' %% na %d/%d',na,n_A)); disp(tmp_cap_); end;
if ( isempty(find(tmp_cap_< 0)));
lp_(1+na) = lP_base - sum(gammaln(1+tmp_cap_),'all');
end;%if ( isempty(find(tmp_cap_< 0)));
end;%for na=0:n_A;
p_ = exp(lp_);
lpv = log(sum(p_(find(lp_<=lP_0))));
if (verbose); disp(sprintf(' %% p_: ')); disp(p_); end;
if (verbose); disp(sprintf(' %% sum p_ %0.2f',sum(p_))); end;
if (verbose); disp(sprintf(' %% lpv %0.2f --> log(%0.2f)',lpv,exp(lpv))); end;
%%%%%%%%;
end;%if ( (n_label_A>=2) & (n_label_B>=2) );
%%%%%%%%;
if (verbose); disp(sprintf(' %% [finished label_pair_enrichment]')); end;
