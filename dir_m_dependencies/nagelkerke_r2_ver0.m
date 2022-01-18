function nagelkerke_r2 = nagelkerke_r2_ver0(X_full_,Y_full_,beta_full_);
%%%%%%%%;
% function output = nagelkerke_r2_ver0(X_full_,Y_full_,beta_full_);
% ;
% Calculates the nagelkerke pseudo-r^2 for data X_full_ Y_full_, assuming ;
% That Y_full_ is a categorical outcome with values in [1,2]. ;
% ;
% Inputs: 
% X_full_: N-by-M array to feed into mnrfit. ;
% Y_full_: N-by-1 array to feed into mnrfit. ;
% beta_full_ : M-by-1 outcome from [beta_full_,dev_full,stats_full] = mnrfit(X_full_,Y_full_) ;
% beta_null_ : 1-by-1 outcome from [beta_null_,dev_null,stats_null] = mnrfit(X_null_,Y_null_) ;
%              where X_null_==0 and Y_null_=Y_full_ correspond to the null hypothesis. ;
%              This is calculated to be beta_null_ = log(ND/NX). ;
% ;
% Outputs:
% nagelkerke_r2:  ( 1 - ( L_null / L_full )^(2/N) ) / Z, ;
%                 where the numerator is the Cox and Snell pseudo-R^2, ;
%                 L_full/L_null is the likelihood associated with the full/null model, ;
%                 and the normalizing factor Z in the denominator is: ;
%             Z:  ( 1 - ( L_null )^(2/N) ). ;
%%%%%%%%;

if nargin<1;
X_ = [1;2;3;4;5;6;7;8];
Y_ = [2;1;2;1;2;1;2;2];
disp(sprintf(' %% testing:\n %% r2 = %0.4f',nagelkerke_r2_ver0(X_,Y_)));
disp('returning');return;
end;%if nargin<1;

[N,M] = size(X_full_);
ND = length(find(Y_full_==1)); NX = length(find(Y_full_==2));
beta_null_ = log(max(1,ND)/max(1,NX));
Y_null_ = Y_full_;
X_null_ = zeros(size(Y_null_));

if nargin<3; warning off; beta_full_ = mnrfit(X_full_,Y_full_); warning on; end;

p_full_ = [ones(N,1) , X_full_] * beta_full_ ;
j_full_ = find(Y_full_==1);
l_full = sum(p_full_(j_full_)) - sum(log(1 + exp(p_full_))); %<-- log L_full. ;

p_null_ = [ones(N,1)] * beta_null_ ;
j_null_ = find(Y_null_==1);
l_null = sum(p_null_(j_null_)) - sum(log(1 + exp(p_null_))); %<-- log L_null. ;

Z = 1 - exp(2*l_null/N);
nagelkerke_r2 = (1 - exp(2*l_null/N - 2*l_full/N)) / Z ;


