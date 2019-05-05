function output = popcount_xor_lf(mem1_,mem2_,mask_,mask_end,din_,popcount_,POPLENGTH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%
% function output = popcount_xor_lf(mem1_,mem2_,mask_,mask_end,din_,popcount_,POPLENGTH);
%
% This function calculates: ;
% \sum__{j=0}^{mask_end-1} din_(1+j/POPLENGTH)*popcount(bitand(bitxor(mem1_(1+j),mem2_(1+j)),mask_(1+j))) ;
% which accumulates the xor of mem1_ and mem2_ from mask_(1+0) to mask_(1+mask_end-1), ;
% multiplied by din_. ;
% ;
% We assume that mem1_ and mem2_ and mask_ are of type uint8. ;
% ;
% We assume that mask indices and mask_end are given ;
% ;
% We also assume that each entry of din_ is applied to bit-chunks of size 'POPLENGTH'. ;
% The value of din_ will be replaced with all ones (i.e., '1') if no value is specified. ;
% The value of POPLENGTH = 1920 will be used if no value is specified. ;
% ;
% We also assume that the popcount function is applied through the lookup table popcount_(uint8_input). ;
% The array popcount_ will be generated if no value is specified. ;
% ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

if (nargin<1);

disp(sprintf(' '));
disp('testing popcount_xor_lf: ');
BIT8 = 8; bitj = 16; POPLENGTH = 1920; POPSUB = POPLENGTH/BIT8;
An_nrows = round(2*POPLENGTH + 13); An_ncols = 12;
An_brows = rup(An_nrows,bitj)/BIT8;
An_ = randn(An_nrows,An_ncols)>0;
disp(sprintf(' converting An_ to binary uAt_: '));
uAt_ = tutorial_integer_to_binary(bitj,An_nrows,An_ncols,An_);
disp(sprintf(' uncompressing uAt_ --> iAt_: '));
iAt_ = tutorial_binary_to_integer(bitj,An_nrows,An_ncols,uAt_);
disp(sprintf('error |An_-iAt_| = %f',norm(An_-iAt_)));
An_mr_ = 1.0*(randn(An_nrows,1)>0);
disp(sprintf(' converting An_mr_ to binary uAt_mr_: '));
uAt_mr_ = tutorial_integer_to_binary(bitj,An_nrows,1,An_mr_);
An_prows = rup(An_nrows,POPLENGTH)/POPLENGTH;
din_= rand(An_prows,1);
X1_ = zeros(An_ncols,1); X2_ = zeros(An_ncols,1);
popcount_ = make_popcount_();
for nc1=0:An_ncols-1; for nc2=0:An_ncols-1;
mem1_ = uAt_(:,1+nc1); mem2_ = uAt_(:,1+nc2); mask_ = uAt_mr_; mask_end = An_brows;
X1_(1+nc1,1+nc2) = popcount_xor_lf(mem1_,mem2_,mask_,mask_end,din_,popcount_,POPLENGTH);
tmp_X2 = 0;
for nr1=0:An_nrows-1;
tmp_B21 = iAt_(1+nr1,1+nc1);
tmp_B22 = iAt_(1+nr1,1+nc2);
tmp__or = (tmp_B21 | tmp_B22);
tmp_and = (tmp_B21 & tmp_B22);
tmp_nnd = ~tmp_and;
tmp_xor = tmp__or & tmp_nnd;
tmp_X2 = tmp_X2 + din_(1 + floor(nr1/POPLENGTH)) * (tmp_xor & An_mr_(1+nr1)) ;
end;%for nr1=0:An_nrows-1;
X2_(1+nc1,1+nc2) = tmp_X2;
end;end;%for nc1=0:An_ncols-1; for nc2=0:An_ncols-1;
disp(sprintf('error |X1_-X2_| = %f',norm(X1_-X2_)));

return;
end;%if (nargin<1);

BIT8 = 8;
brows = length(mem1_);
assert(mask_end<=brows);

if (nargin<7);
POPLENGTH = 1920;
end;%if (nargin<7);

if (nargin<6);
popcount_ = make_popcount_();
end;%if (nargin<6);

POPSUB = POPLENGTH/BIT8;
if (nargin<5);
din_ = ones(1,rup(brows,POPSUB)/POPSUB);
end;%if (nargin<5);

din__ = reshape(repmat(reshape(din_,1,length(din_)),POPSUB,1),rup(brows,POPSUB),1);
output = sum(din__(1:brows).*double(intlut(bitand(bitxor(mem1_(1:mask_end),mem2_(1:mask_end)),mask_(1:mask_end)),popcount_)));
