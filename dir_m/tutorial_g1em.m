function pp = tutorial_g1em(eps,m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% pp = tutorial_g1em(eps,m);
%
% This function calculates the probability pp that a 2-x-2 binarized loop within ;
% a large rank-1 matrix is rank-1, given that the large matrix is size m-x-m, ;
% with spectral error eps (i.e., sigma_2/sigma_1). ;
% The value of pp actually only depends on eps*sqrt(m). ;
% 
% The inputs are eps (a double) describing the spectral error, and m (an integer). ;
% If only one input is given, it is treated as eps*sqrt(m) (i.e., m is set to 1). ;
% The output is pp. ;
%
% test by running with no arguments:
% i.e., >> tutorial_g1em();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1;
disp(sprintf(' testing tutorial_g1em: '));
esm_ra = 10.^(-3:0.125:1);
for ne=1:length(esm_ra);
pp_ra(ne) = tutorial_g1em(esm_ra(ne));
end;%for ne=1:length(esm_ra);
plot(log10(esm_ra),pp_ra,'k.-'); xlabel('log10(eps*sqrt(m))'); ylabel('g_{1,eps,m}')
return;
end;%if nargin<1;

if nargin<2; m=1; end;
esm = eps*sqrt(m);
N_w = 128;
w_ra = linspace(0,1,N_w); w_ra = 2*pi*(w_ra(1:end-1) + 1/(2*N_w));
dw = mean(diff(w_ra));
pp_ra = zeros(length(w_ra),1);
for nw=1:length(w_ra);
w = w_ra(nw);
pval_r = @(r) 1/pi * max( atan2( - (esm./r) .* sqrt(1 + (esm./r).^2) , 0.5 * sin(2*w) ) , atan2( (esm./r) .* sqrt(1 + (esm./r).^2) , - 0.5 * sin(2*w) ) );
ival_r = @(r) 1/2/pi * (1 - 2*pval_r(r) + 2*pval_r(r).^2) .* r .* exp(-r.^2 / 2);
pp_ra(nw) = integral(ival_r,0.00001,10);
%pp_ra(nw) = integral(ival_r,0.00000001,10);
end;%for nw=1:length(w_ra);
pp = dw*sum(pp_ra);



