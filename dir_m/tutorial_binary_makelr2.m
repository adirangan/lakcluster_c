function [output] = tutorial_binary_makelr2(nrows,ncols,pp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [output,p] = tutorial_binary_makelr2(nrows,ncols,pp);
%
% This function generates a binary (+1/-1) matrix A of size nrows-x-ncols which has ;
% the requisite probability pp of each loop being rank-1. ;
% The inputs are nrows and ncols (integers) dictating the size of A, ;
%     the probability pp (a double between 0.5 and 1) is the probability that a random ;
%     2x2 submatrix of A will be rank-1 (rather than rank-2). ;
% The output is the matrix A. ;
%
% test by running with no arguments:
% i.e., >> tutorial_binary_makelr2();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2;
disp(sprintf(' testing tutorial_binary_makelr2: '));
nrows=512;ncols=512;esm = 0.75; pp = tutorial_g1em(esm,1);
disp(sprintf(' generating test data A of size %dx%d, esm %0.2f, pp %0.4f',nrows,ncols,esm,pp));
A = tutorial_binary_makelr2(nrows,ncols,pp);
imagesc(A,[-1,1]); title('heatmap of A');
ZC = diag(A'*A*A'*A) - nrows*(ncols + nrows - 1);
T = (nrows*(nrows-1)*ncols*(ncols-1)); f = (sum(ZC) + T)/2/T;
disp(sprintf(' requisite pp %0.4f, empirical pp %0.4f',pp,f));
B = 2*(tutorial_makelr(nrows,ncols,1,esm/sqrt(sqrt(nrows*ncols)),0)>0)-1;
ZC = diag(B'*B*B'*B) - nrows*(ncols + nrows - 1);
T = (nrows*(nrows-1)*ncols*(ncols-1)); f = (sum(ZC) + T)/2/T;
disp(sprintf(' empirical pp of matrix B generated with eps*sqrt(m) = %0.2f: %0.4f',esm,f));
return;
end;%if nargin<2;

if nargin<3; pp = 1.00; end;

pp = min(1.0,max(0.5,pp));
p = pp_to_p(pp);

u = 2*(randn(nrows,1)>0)-1; v = 2*(randn(ncols,1)>0)-1;
A = u*transpose(v); B = reshape(A,nrows*ncols,1);
C = zeros(nrows*ncols,1); prm = randperm(nrows*ncols,ceil(p*(nrows*ncols))); C(prm) = -2*B(prm); C = reshape(C,nrows,ncols); % norm 2*p*sqrt(nrows*ncols);
output = reshape(A+C,nrows,ncols); %norm (1-2p)*sqrt(nrows*ncols); s2 close to 4*sqrt(sqrt(nrows*ncols)*p*q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = pp_to_p(pp);
% pp = probability of a loop being rank-1. ;
% qq = probability of a loop being rank-2. ;
% p = probability of an (entrywise) random replacement. ;
% q = probability of no random replacement. ;
% we see that qq = 4pq^3 + 4qp^3 = 4pq*(q^2+p^2) = 4pq*(1-2pq) ;
pp = min(1.0,max(0.5,pp));
qq = 1-pp;
pq = 0.25 * (1 - 0.25 * sqrt(16 - 32*qq));
p = 0.5 * (1 - sqrt(1 - 4*pq));
