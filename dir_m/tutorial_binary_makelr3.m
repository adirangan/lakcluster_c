function [output] = tutorial_binary_makelr3(nrows,ncols,pp,sc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [output,p] = tutorial_binary_makelr3(nrows,ncols,pp);
%
% This function generates a sparse binary (+1/-1) matrix A of size nrows-x-ncols which has ;
% the requisite probability pp of each loop being rank-1. ;
% The inputs are nrows and ncols (integers) dictating the size of A, ;
%     the probability pp (a double between 0.5 and 1) is the probability that a random ;
%     2x2 submatrix of A will be rank-1 (rather than rank-2), ;
%     and the sparsity-coefficient sc (a probability between 0 and 1). ;
% The output is the matrix A. ;
%
% test by running with no arguments:
% i.e., >> tutorial_binary_makelr3();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<2;
disp(sprintf(' testing tutorial_binary_makelr3: '));
disp(sprintf(' generating test data A of size 256x512'));
nrows=256;ncols=512;pp = 0.65;sc = 0.15;
A = tutorial_binary_makelr3(nrows,ncols,pp,sc);
imagesc(A,[-1,1]); title('heatmap of A');
ZC = diag(A'*A*A'*A) - nrows*(ncols + nrows - 1);
T = (nrows*(nrows-1)*ncols*(ncols-1));
f = (sum(ZC) + T)/2/T;
disp(sprintf(' requisite pp %0.2f, empirical pp %0.2f',pp,f));
sc_tmp = mean(A(:)>0); 
disp(sprintf(' requisite sc %0.2f, empirical sc %0.2f',sc,sc_tmp));
return;
end;%if nargin<2;

if nargin<3; pp = 1.00; end;
if nargin<4; sc = 0.5; end;

sc_flip = min(sc,1-sc); % flipped to have more -1s if sc>0.5;
pp_loop = min(1.0,max(0.5,pp));
p_indiv = pp_to_p(pp_loop);
p_indiv = min(p_indiv,sc_flip);
q_indiv = 1-p_indiv;
sc_flop = (sc_flip-p_indiv)/(q_indiv-p_indiv); % more extreme sparsity before perturbation ;
sc_p = 0.5 - 0.5*sqrt(1 - 2*min(sc_flop,1-sc_flop)); % component sparsity, flipped to have more -1s if sc_flop>0.5. ;

u = 2*(rand(nrows,1)<sc_p)-1; v = 2*(rand(ncols,1)<sc_p)-1;
A = -u*transpose(v); B = reshape(A,nrows*ncols,1); B_up = find(B==+1); B_dn = find(B==-1);
prm = randperm(nrows*ncols,ceil(p_indiv*nrows*ncols));
C = zeros(nrows*ncols,1); C(prm) = -2*B(prm);
C = reshape(C,nrows,ncols); % norm 2*p*sqrt(nrows*ncols);
output = reshape(A+C,nrows,ncols); %norm (1-2p)*sqrt(nrows*ncols); s2 close to 4*sqrt(sqrt(nrows*ncols)*p*q);
if (sc>0.5); output = -output; end; % flipped back to have more +1s if sc>0.5 ;

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
