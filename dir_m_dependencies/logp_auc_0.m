function output = logp_auc_0(A,N);
% returns log of p-value for auc A given N=minimum(ND,NX);
% asymptotic formula assuming both ND,NX are large and ;
% that the probability of identical v-values is small. ;

if nargin<1;
%%%%%%%%;
ND = 32; NX = 1024;
n_iteration = 1024*16;
%%%%%%%%;
a_ = zeros(n_iteration,1);
for niteration=1:n_iteration;
if (mod(niteration,100)==0); disp(sprintf(' %% ni %d/%d',niteration,n_iteration)); end;
a_(niteration) = get_auc(rand(NX,1),rand(ND,1));
end;%for niteration=1:n_iteration;
%%%%%%%%;
n_p = 64; p_ = linspace(0,100,n_p);
p_a_ = zeros(n_p,1);
for np=1:n_p;
p=p_(np); p_a_(np) = prctile(a_,p);
end;%for np=1:n_p;
%%%%%%%%;
n_a = 64; A_ = linspace(0.4,0.6,n_a);
P_A_ = zeros(n_a,1);
for na=1:n_a;
A=A_(na); P_A_(na) = exp(logp_auc_0(A,min(ND,NX)));
end;%for na=1:n_a;
%%%%%%%%;
plot(A_,log(P_A_),'kx-',p_a_,log(1-p_/100),'r.-')
xlabel('auc'); ylabel('probability');legend('analytical','empirical');
%%%%%%%%;
disp('returning'); return;
end;%if nargin<1;

A2 = A-0.5;%A2 = abs(A-0.5); 
A3 = A2*sqrt(6*N);
%logp = log(0.5*erfc(A3));
if (A3<5);
logp = log(0.5*erfc(A3));
else;
logp = log(0.5*g2(A3));
end;%if;
output = logp;

function output = g1(x_);
output = exp(-x_.^2)./x_./sqrt(pi);

function output = g2(x_);
output = exp(-x_.^2)./x_./sqrt(pi).*(1-1./(2.*x_.^2));

