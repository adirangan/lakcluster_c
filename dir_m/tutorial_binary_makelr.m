function [output,p] = tutorial_binary_makelr(nrows,ncols,eps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [output,p] = tutorial_binary_makelr(nrows,ncols,eps);
%
% This function generates a binary (+1/-1) matrix A of size nrows-x-ncols which has ;
% numerical rank 1 with error eps. 
% The inputs are nrows and ncols (integers) dictating the size of A, ;
%     the error eps (a double betwen 0 and 1, default 0), ;
%     representing the ratio between the second and first singular values of A. ;
% The output is the matrix A, as well as the fraction p of elements of A that ;
%     are perturbed from their rank-1 initialization. ;
%
% test by running with no arguments:
% i.e., >> tutorial_binary_makelr();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<2;
disp(sprintf(' testing tutorial_binary_makelr: '));
disp(sprintf(' generating test data A of size 256x512'));
nrows=256;ncols=256;eps=0.2;
A = tutorial_binary_makelr(nrows,ncols,eps);
imagesc(A,[-1,1]); title('heatmap of A');
[u,s,v] = svds(A,2); esm_tmp = s(2,2)/s(1,1) * sqrt(sqrt(nrows*ncols));
disp(sprintf(' requisite esm %0.2f empirical esm %0.2f',eps*sqrt(sqrt(nrows*ncols)),esm_tmp));
return;
end;%if nargin<2;

if nargin<3; eps = 0.00; end;

if eps<=0; p=0;
elseif eps>=1; p=0.5;
elseif eps>0 & eps<1;
sMN = sqrt(nrows*ncols);
%eps_diff = @(p) p_to_eps(p,sMN) - eps;
%p = fzero(eps_diff,0.25);
c = eps^2*sMN/4; d = sqrt(1-c/(1+c));
p = (1-d)/2;
p=max(0,min(0.5,p));
end;
%disp(sprintf(' %% eps %0.2f, p %0.2f',eps,p));

u = 2*(randn(nrows,1)>0)-1; v = 2*(randn(ncols,1)>0)-1;
A = u*transpose(v); B = reshape(A,nrows*ncols,1);
C = zeros(nrows*ncols,1); prm = randperm(nrows*ncols,ceil(p*(nrows*ncols))); C(prm) = -2*B(prm); C = reshape(C,nrows,ncols); % norm 2*p*sqrt(nrows*ncols);
output = reshape(A+C,nrows,ncols); %norm (1-2p)*sqrt(nrows*ncols); s2 close to 4*sqrt(sqrt(nrows*ncols)*p*q);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eps = p_to_eps(p,sMN);
% this is the function: ;
% eps = 4sqrt(sMN*p*q) / (sMN*(p-q)), ;
% where sMN = sqrt(M*N), and q = 1-p. ;
sMN = max(1,sMN); q=1-p;
if p<=0; eps = 0;
elseif p>=0.5; eps=1;
elseif p>0 & p<0.5; 
eps = 4*sqrt(sMN.*p.*q) ./ (sMN.*(q-p)) ;
end;

%{

  clear;
  M = 564; N = 328; 
  p_ra = (2.^[linspace(-9,-1,25)]);
  q_ra = 1-p_ra;
  max_iteration = 8;
  for np=1:length(p_ra);
  p = p_ra(np); q = 1-p;
  disp(sprintf(' %% p %0.2f, q %0.2f',p,q));
  for ni=1:max_iteration;
  %u = 2*(randn(M,1)>0)-1; v = 2*(randn(N,1)>0)-1; % could be all ones;
  %A = u*transpose(v); B = reshape(A,M*N,1);
  A = ones(M,N); B = ones(M*N,1);
  C = zeros(M*N,1); prm = randperm(M*N,ceil(p*(M*N))); C(prm) = -2*B(prm); C = reshape(C,M,N); % norm 2*p*sqrt(M*N);
  D = reshape(A+C,M,N); %norm (1-2p)*sqrt(M*N);
  sA = svds(A,1); % should be simply sqrt(M*N) ;
  sC = svds(C,1); sD = svds(D,2);
  sA_ra(np,ni) = sA;
  sC_ra(np,ni) = sC;
  sD1_ra(np,ni) = sD(1); % close to (1-2*p)*sqrt(M*N) = (q-p)*sqrt(M*N);
  sD2_ra(np,ni) = sD(2); % close to 4*sqrt(sqrt(M*N)*p*q);
  %[U,S,V] = svds(D,2);
  %Un = sign(U(:,1))/sqrt(M); %close to U(:,1);
  %Vn = sign(V(:,1))/sqrt(N); %close to V(:,1);
  %Un = ones(M,1)/sqrt(M); %close to U(:,1);
  %Vn = ones(N,1)/sqrt(N); %close to V(:,1);
  %Dtmp = (D-Un(:,1)*transpose(Vn(:,1))*(1-2*p)*sqrt(M*N)); % norm = sD(2);
  %Dtmp = Dtmp/2; % 11 - C - (1-2p)11 = 2p11-C;  % norm = sD(2)/2;
  %DDtmp = Dtmp*transpose(Dtmp);
  %DDr = DDtmp - diag(diag(DDtmp));
  end;%for ni=1:max_iteration;
  end;%for np=1:length(p_ra);
  mA = exp(mean(log(sA_ra),2));
  mC = exp(mean(log(sC_ra),2));
  mD1 = exp(mean(log(sD1_ra),2));
  mD2 = exp(mean(log(sD2_ra),2));
  subplot(1,2,1);loglog((p_ra),(sD2_ra),'b.-',(p_ra),4*sqrt((sqrt(M*N)*p_ra.*q_ra)),'ro-'); xlabel('log2(p)'); ylabel('log2(4sqrt(sqrt(M*N)pq))');
  subplot(1,2,2);plot(p_ra,mD2./mD1,'b.-',p_ra,4*sqrt(p_ra.*q_ra)./(q_ra-p_ra)/sqrt(sqrt(M*N)),'ro-',p_ra,4*sqrt(p_ra)/sqrt(sqrt(M*N)),'go-');

  %}
