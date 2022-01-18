function ...
[ ...
 parameter ...
,g_opt_ ...
,nlp_aA_opt_ ...
,nlp_aA_emp_ ...
,p_aA_opt_ ...
,p_aA_emp_ ...
] = ...
gumbel_fit_0( ...
 parameter ...
,A_ ...
,a_ ...
);
% Fits gumbel distribution to the data in A_. ;
% If a_ is provided, estimates one-sided p-value of a_ relative to A_. ;
% post-pends the p-value of the elements of A_. ;
if (nargin<1);
A_ = 0.5 + (0:7);
a_ = [6,7];
[g_opt_,nlp_opt_,nlp_emp_,p_opt_,p_emp_] = gumbel_fit(A_,a_);
disp(sprintf(' %% g_opt_ %0.16f %0.16f',g_opt_));
disp(sprintf(' %% nlp_opt_ %0.16f %0.16f',nlp_opt_));
disp(sprintf(' %% nlp_emp_ %0.16f %0.16f',nlp_emp_));
disp(sprintf(' %% p_opt_ %0.16f %0.16f',p_opt_));
disp(sprintf(' %% p_emp_ %0.16f %0.16f',p_emp_));
disp('returning');return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); A_=[]; end; na=na+1;
if (nargin<1+na); a_=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;

verbose=0;
%%%%%%%%;
g_ori_ = [ 0 ; 1 ]; %<-- mu, sqrt(beta) for gumbel. ;
A_ = A_(:);
n_u = numel(A_);
B_ = (A_ - mean(A_))/std(A_,1); %<-- center and normalize A_. ;
tmp_gumbel_nll = @(g_) sum(gumbel_nll(B_(:),g_));
g_opt_ = g_ori_;
g_opt_ = fminsearch(tmp_gumbel_nll,g_ori_,optimset('MaxIter',1e4,'MaxFunEvals',1e4));
if (verbose);
figure(1);clf;
n_h = 32;
h_x_lim_ = mean(B_) + std(B_,1)*3.5*[-1,+1];
h_x_ = linspace(min(h_x_lim_),max(h_x_lim_),n_h);
h_B_ = hist(B_,h_x_);
h_B_ = h_B_/sum(h_B_)/mean(diff(h_x_));
x_ = linspace(min(h_x_lim_),max(h_x_lim_),1024);
rho_ = gumbel_pdf(x_,g_opt_);
hold on;
stairs(h_x_,h_B_,'k-','LineWidth',2);
plot(x_,rho_,'r-','LineWidth',2);
hold off;
xlabel('x');
ylabel('dist');
title('gumbel_fit','Interpreter','none');
end;%if (verbose);
%%%%%%%%;
nlp_a_opt_ = []; nlp_a_emp_ = [];
p_a_opt_ = []; p_a_emp_ = [];
%%%%%%%%;
if ~isempty(a_);
n_a = numel(a_); n_A = numel(A_);
nlp_a_opt_ = zeros(n_a,1); nlp_a_emp_ = zeros(n_a,1);
p_a_opt_ = zeros(n_a,1); p_a_emp_ = zeros(n_a,1);
for na=0:n_a-1;
a = a_(1+na);
p_a_emp = ( numel(find(A_>a)) + 0.5*numel(find(A_==a)) )/n_u;
nlp_a_emp = -log(p_a_emp);
b = (a - mean(A_))/std(A_,1);
p_a_opt = 1 - gumbel_cdf(b,g_opt_);
nlp_a_opt = -log(p_a_opt);
p_a_emp_(1+na) = p_a_emp;
nlp_a_emp_(1+na) = nlp_a_emp;
p_a_opt_(1+na) = p_a_opt;
nlp_a_opt_(1+na) = nlp_a_opt;
end;%for na=0:n_a-1;
end;%if ~isempty(a_);
%%%%%%%%;
p_A_emp_ = zeros(n_u,1);
[~,ij_uns_from_srt_] = sort(A_(:),'ascend');
p_A_emp_(ij_uns_from_srt_) = (flip(0:n_u-1) + 0.5)/n_u;
nlp_A_emp_ = -log(p_A_emp_);
p_A_opt_ = 1 - gumbel_cdf(B_,g_opt_);
nlp_A_opt_ = -log(p_A_opt_);
%%%%%%%%;
nlp_aA_opt_ = [nlp_a_opt_ ; nlp_A_opt_];
p_aA_opt_ = [p_a_opt_ ; p_A_opt_];
nlp_aA_emp_ = [nlp_a_emp_ ; nlp_A_emp_];
p_aA_emp_ = [p_a_emp_ ; p_A_emp_];
