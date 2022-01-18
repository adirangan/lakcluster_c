function [lp,cap_,cup_] = label_to_label_enrichment_0(label_A_,label_B_);
% Calculates enrichment p-value. ;
% Assumes numeric labels, although labels do not need to be sequential. ;
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
if (n_A~=n_B); disp(sprintf(' %% Warning! n_A %d n_B %d in label_to_label_enrichment',n_A,n_B)); end;
assert(n_A==n_B);
n_X = n_A;
%%%%%%%%;
cap_ = zeros(n_label_A,n_label_B);
cup_ = zeros(n_label_A,n_label_B);
for nlabel_A=1:n_label_A;
for nlabel_B=1:n_label_B;
cap_(nlabel_A,nlabel_B) = numel(intersect(ij_A_{nlabel_A},ij_B_{nlabel_B}));
cup_(nlabel_A,nlabel_B) = numel(union(ij_A_{nlabel_A},ij_B_{nlabel_B}));
end;%for nlabel_B=1:n_label_B;
end;%for nlabel_A=1:n_label_A;

%{

  n_X = 30;
  n_label_A = 4;
  n_label_A_ = [4;5;10;11];
  n_label_B = 2;
  n_label_B_ = [9;21];
  tmp_lN = sum(gammaln(1+n_label_A_));
  tmp_lP = sum(gammaln(1+n_label_A_)) + sum(gammaln(1+n_label_B_)) - gammaln(1+n_X);
  cap_0_ = [4,0;5,0;0,10;0,11];
  tmp_0 = sum(gammaln(1+cap_0_),'all');
  lN_0 = tmp_lN - tmp_0;
  lP_0 = tmp_lP - tmp_0;

  cap_0_ = [1,3;2,3;2,8;4,7];
  tmp_0 = sum(gammaln(1+cap_0_),'all');
  lN_0 = tmp_lN - tmp_0;
  lP_0 = tmp_lP - tmp_0;

  n_label_B = 1;
  n_label_B_ = [30];
  tmp_lN = sum(gammaln(1+n_label_A_));
  tmp_lP = sum(gammaln(1+n_label_A_)) + sum(gammaln(1+n_label_B_)) - gammaln(1+n_X);
  cap_0_ = [4;5;10;11];
  tmp_0 = sum(gammaln(1+cap_0_),'all');
  lN_0 = tmp_lN - tmp_0;
  lP_0 = tmp_lP - tmp_0;
  
  %}

num_A_ = zeros(n_A,1);
for nlabel_A=1:n_label_A;
num_A_(ij_A_{nlabel_A}) = nlabel_A;
end;%for nlabel_A=1:n_label_A;
num_B_ = zeros(n_B,1);
for nlabel_B=1:n_label_B;
num_B_(ij_B_{nlabel_B}) = nlabel_B;
end;%for nlabel_B=1:n_label_B;
tmp_lP = sum(gammaln(1+n_label_A_)) + sum(gammaln(1+n_label_B_)) - gammaln(1+n_X);
cap_0_ = sparse(num_A_,num_B_,1,n_label_A,n_label_B);
tmp_0 = sum(gammaln(1+cap_0_),'all');
lN_0 = +tmp_0;
lP_0 = tmp_lP - tmp_0;
cap_A_ = sparse(num_A_,num_A_,1,n_label_A,n_label_A);
tmp_A = sum(gammaln(1+cap_A_),'all');
lN_A = +tmp_A;
lP_A = 2*sum(gammaln(1+n_label_A_)) - gammaln(1+n_X) - tmp_A;
%%%%%%%%;
% Given a list of specific sets from A, ;
% as well as the cardinality of the sets in B, ;
% the number of ways of finding the precise intersections listed in cap_ is: ;
% lN = sum(gammaln(1+n_label_A_)) - sum(gammaln(1+cap_(:))). ;
% Similarly, the probability of observing cap_ is: ;
% lP = sum(gammaln(1+n_label_A_)) - sum(gammaln(1+cap_(:))) + sum(gammaln(1+n_label_B_)) - gammaln(1+n_X). ;
%%%%%%%%;
n_iteration=1024*8;
lN_ = zeros(n_iteration,1);
lP_ = zeros(n_iteration,1);
cap__ = zeros(n_label_A,n_label_B,n_iteration);
for niteration=1:n_iteration;
tmp_cap_ = sparse(num_A_,num_B_(randperm(n_B)),1,n_label_A,n_label_B);
tmp = sum(gammaln(1+tmp_cap_),'all');
lN_(niteration) = +tmp;
lP_(niteration) = tmp_lP - tmp;
cap__(:,:,niteration) = tmp_cap_;
end;%for niteration=1:n_iteration;
clear tmp_cap_;
% compare mean(cap__,3) with (n_label_A_*transpose(n_label_B_)/n_X) ;
% try: ;
%{
  cap_avg__ = mean(cap__,3);
  cap_bar__ = (n_label_A_*transpose(n_label_B_)/n_X);
  disp(sprintf(' %% relative deviation: %0.16f',fnorm(cap_avg__-cap_bar__)/fnorm(cap_bar__)));
  %}
[~,u_cap_ij_] = unique(transpose(reshape(cap__,n_label_A*n_label_B,n_iteration)),'row');
clear cap__;
n_u = numel(u_cap_ij_);
lN_u_ = lN_(u_cap_ij_);
lP_u_ = lP_(u_cap_ij_);

%%%%%%%%;
% Now model the probability of observing lN_u (w.r.t. sampling label_B_) ;
% as a gaussian in lN_u. ;
%%%%%%%%;
mu = mean(lN_u_); sg = std(lN_u_,1);
flag_plot=0;
if flag_plot;
n_h = 33;
hlN_ = mu + sg*4.5*linspace(-1,1,n_h);
h_ = hist(lN_u_,hlN_); h_ = h_/sum(h_)/mean(diff(hlN_));
g_ = 1/sqrt(2*pi) / sg * exp(-(hlN_-mu).^2/(2*sg^2));
subplot(1,2,1); plot(hlN_,g_,'r-',hlN_,h_,'k.'); xlabel('lN_u'); ylabel('p'); title('p');
subplot(1,2,2); plot(hlN_,log(g_),'r-',hlN_,log(h_),'k.'); ylim([min(log(g_)),max(log(g_))]); xlabel('lN_u'); ylabel('log(p)'); title('log(p)');
end;%if flag_plot;
z = (lN_0-mu)/sg;
%f = sqrt(pi)/sqrt(2)*sg*erfc(z/sqrt(2)) / exp(-z^2/2) ; %<-- integral( @(x) exp(-(x-mu).^2/(2*sg^2)) , mu + z*sg , Inf ) ;
lf = log(sqrt(pi)/sqrt(2)*sg) + erfcln(z/sqrt(2)) + z^2/2 ;
lp = z_to_p_0(z);
lpf = lP_0 + lf;




