function output = tutorial_makelr(nrows,ncols,k,eps,nrm_flag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [output] = tutorial_makelr(nrows,ncols,k,eps,nrm_flag);
%
% This function generates a matrix A of size nrows-x-ncols which has ;
% numerical rank k with error eps. A is randomly generated to resemble a ;
% multivariate gaussian distribution. Specifically, the first k singular ;
% values of A are 1, whereas singular values k+1 and beyond equal eps, ;
% with the orientation of the principal components chosen uniformly at random. ;
% The inputs are nrows and ncols (integers) dictating the size of A, ;
%     the numerical rank k (an integer, default 1), ;
%     the error eps (a double betwen 0 and 1, default 0), ;
%     a normalization_flag nrm_flag (either 0 or 1, default 0). ;
% If nrm_flag==1, then then the values of A are normalized to match ;
%     variance with those in a uniformly gaussian random matrix. ;
% The output is A;
%
% test by running with no arguments:
% i.e., >> tutorial_makelr();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% makes an eps-error k-rank block;
% tries to balance the variance ;
% test with: ;
%{

  N=512;nbins=128;dim=3;eps=0.1;
  nrows = N; ncols = 0.25*N;
  sigma = sqrt(dim*(1-eps.^2)/(nrows*ncols) + eps.^2/max(nrows,ncols));
  hbins = linspace(-6*sigma,+6*sigma,nbins);
  A = tutorial_makelr(nrows,ncols,dim,eps);
  h = hist(A(:),hbins); h = h/sum(h)/mean(diff(hbins));
  he = 1/sqrt(2*pi*sigma.^2)*exp(-hbins.^2/2/sigma.^2);
  figure;cla;plot(hbins,h,'ro-',hbins,he,'k-'); xlim([min(hbins),max(hbins)]);
  disp(sprintf('std: %f vs %f',sigma,std(A(:))));
  N = 128;
  nrows = N; ncols = 0.5*N;
  A1 = tutorial_makelr(nrows,ncols,1,0.1,1);
  A2 = tutorial_makelr(nrows,ncols,3,0.1,1);
  A3 = tutorial_makelr(nrows,ncols,1,1.0,1);
  A4 = tutorial_makelr(nrows,ncols,3,1.0,1);
  A = [A1 , A2 , A3 , A4];
  figure;
  subplot(2,1,1);cla;imagesc(A);
  subplot(2,2,3);cla;plot(mean(A));title('mean');
  subplot(2,2,4);cla;plot(var(A));title('var');

 %}

if nargin<2;
disp(sprintf(' testing tutorial_makelr: '));
disp(sprintf(' generating test data A of size 512x512'));
N=512;nbins=128;dim=3;eps=0.1;
nrows = N; ncols = 0.25*N;
sigma = sqrt(dim*(1-eps.^2)/(nrows*ncols) + eps.^2/max(nrows,ncols));
hbins = linspace(-6*sigma,+6*sigma,nbins);
A = tutorial_makelr(nrows,ncols,dim,eps);
h = hist(A(:),hbins); h = h/sum(h)/mean(diff(hbins));
he = 1/sqrt(2*pi*sigma.^2)*exp(-hbins.^2/2/sigma.^2);
figure;cla;
subplot(2,2,1);
disp(sprintf(' plotting distribution of values of A.'));
plot(hbins,h,'ro-',hbins,he,'k-'); xlim([min(hbins),max(hbins)]);
title('distribution of A-values');
disp(sprintf(' sigma: %f vs distribution std %f',sigma,std(A(:))));
disp(sprintf(' generating test data A of size 256x256'));
N = 128;
nrows = N; ncols = 0.5*N;
A1 = tutorial_makelr(nrows,ncols,1,0.1,1);
A2 = tutorial_makelr(nrows,ncols,3,0.1,1);
A3 = tutorial_makelr(nrows,ncols,1,1.0,1);
A4 = tutorial_makelr(nrows,ncols,3,1.0,1);
A = [A1 , A2 , A3 , A4];
disp(sprintf(' plotting mean and variance for each column of A.'));
subplot(2,2,2);cla;imagesc(A);title('heatmap of A');
subplot(2,2,3);cla;plot(mean(A));title('mean');xlim([1,size(A,2)]);
subplot(2,2,4);cla;plot(var(A));title('var');xlim([1,size(A,2)]);
disp(sprintf(' printing figure (see test.jpg)'));
print('-djpeg','./test.jpg');
return;
end;%if nargin<2;

if nargin<5; nrm_flag=0; end;
if nargin<4; eps=0; end;
if nargin<3; k=1; end;

output = randn(nrows,ncols);
[U,S,V] = svds(output,max(1024,4*k)); 
%S=diag(S); S(k+1:end) = eps*S(k+1:end); S=diag(S); 
S=diag(S); S(1:k)=1; S(k+1:end) = eps; S=diag(S); 
sigma=1;
if nrm_flag==1;
sigma = sqrt(k*(1-eps.^2)/(nrows*ncols) + eps.^2/max(nrows,ncols));
end;%if nrm_flag;
output = U*S*transpose(V)/sigma;
