function output_ = AnZt_vv(uAn_,uZn_,A_ajdk_,An_ajdk_,Zn_ajdk_,umc_,popcount_,POPLENGTH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
%
% function output_ = AnZt_vv(uAn_,uZn_,A_ajdk_,An_ajdk_,Zn_ajdk_,umc_,popcount_,POPLENGTH);
%
% This function calculates: ;
% output_(1+nr1,1+nr2) = \sum__{nc=0}^{mask_end-1} (uAn_(1+nr1,1+nc) - A_ajdk(1+nc + AJDK_1_0*A_pcols))*A_AJDK(1+nc + AJDK_1_1*A_pcols)*(uZn_(1+nr1,1+nc) - A_ajdk(1+nc + AJDK_1_0*A_pcols))*umc_(1+nc) ;
% which calculates the (sparsity-corrected) inner-product An*Zt, limited to the mask umc_. ;
% ;
% We assume that uAn_, uZn_ and umc_ are of type uint8. ;
% ;
% We assume that mask indices and mask_end are given ;
% ;
% We also assume that each entry of A_ajdk_ is applied to bit-chunks of size 'POPLENGTH'. ;
% The value of POPLENGTH = 1920 will be used if no value is specified. ;
% ;
% We also assume that the popcount function is applied through the lookup table popcount_(uint8_input). ;
% The array popcount_ will be generated if no value is specified. ;
% ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

if (nargin<1);

disp(sprintf(' '));
disp(' testing AnZt_vv: ');
AJDK_AXDX_define;
popcount_ = make_popcount_();
BIT8 = 8; bitj = 16; POPLENGTH = 1920; POPSUB = POPLENGTH/BIT8;
At_nrows = round(2*POPLENGTH + 13);
At_ncols = 12; Zt_ncols = 17;
At_ = randn(At_nrows,At_ncols)>0;
uAn_ = tutorial_integer_to_binary(bitj,At_nrows,At_ncols,At_);
iAn_ = tutorial_binary_to_integer(bitj,At_nrows,At_ncols,uAn_);
Zt_ = randn(At_nrows,Zt_ncols)>0;
uZn_ = tutorial_integer_to_binary(bitj,At_nrows,Zt_ncols,Zt_);
iZn_ = tutorial_binary_to_integer(bitj,At_nrows,Zt_ncols,uZn_);
disp(sprintf('error |At_-iAn_| = %f',norm(At_-iAn_)));
disp(sprintf('error |Zt_-iZn_| = %f',norm(Zt_-iZn_)));
At_mr_ = 1.0*(randn(At_nrows,1)>0);
uAn_mr_ = tutorial_integer_to_binary(bitj,At_nrows,1,At_mr_);
At_prows = rup(At_nrows,POPLENGTH)/POPLENGTH;
A_p_ = 0.1 + 0.8*rand(At_prows,1);
A_ajdk_ = calc_A_ajdk(At_prows,A_p_);
An_ajdk_ = An_ajdk_v(uAn_,A_ajdk_,uAn_mr_,POPLENGTH);
Zn_ajdk_ = An_ajdk_v(uZn_,A_ajdk_,uAn_mr_,POPLENGTH);
X1_ = AnZt_vv(uAn_,uZn_,A_ajdk_,An_ajdk_,Zn_ajdk_,uAn_mr_,popcount_,POPLENGTH);
X2_ = zeros(At_ncols,Zt_ncols);
for nc1=0:At_ncols-1; for nc2=0:Zt_ncols-1;
tmp_X2 = 0;
for nr1=0:At_nrows-1;
tmp_At = 2*iAn_(1+nr1,1+nc1)-1;
tmp_Zt = 2*iZn_(1+nr1,1+nc2)-1;
tmp_a = A_ajdk_(1 + floor(nr1/POPLENGTH) + AJDK_1_0*At_prows);
tmp_D = A_ajdk_(1 + floor(nr1/POPLENGTH) + AJDK_0_1*At_prows);
tmp_X2 = tmp_X2 +  (tmp_At - tmp_a) * tmp_D * (tmp_Zt - tmp_a) * At_mr_(1+nr1) ;
end;%for nr1=0:At_nrows-1;
X2_(1+nc1,1+nc2) = tmp_X2;
end;end;%for nc1=0:At_ncols-1; for nc2=0:Zt_ncols-1;
disp(sprintf('error |X1_-X2_| = %f',norm(X1_-X2_)));

return;
end;%if (nargin<1);

AJDK_AXDX_define;
BIT8 = 8;
[At_brows,At_ncols] = size(uAn_);
[Zt_brows,Zt_ncols] = size(uZn_);
assert(At_brows==Zt_brows);
An_bcols = At_brows;

if (nargin<8);
POPLENGTH = 1920;
end;%if (nargin<8);

if (nargin<7);
popcount_ = make_popcount_();
end;%if (nargin<7);

POPSUB = POPLENGTH/BIT8;
An_pcols = rup(At_brows,POPSUB)/POPSUB;

if (~isempty(A_ajdk_));
dtmp_a0d1 = 0;
for nc=0:An_bcols-1;
dtmp_a0d1 = dtmp_a0d1 + A_ajdk_(1+floor(nc/POPSUB) + AJDK_0_1*An_pcols)*double(intlut(umc_(1+nc),popcount_));
end;%for nc=0:An_bcols-1;
dtmp_a2d1 = 0;
for nc=0:An_bcols-1;
dtmp_a2d1 = dtmp_a2d1 + A_ajdk_(1+floor(nc/POPSUB) + AJDK_2_1*An_pcols)*double(intlut(umc_(1+nc),popcount_));
end;%for nc=0:An_bcols-1;
end;%if (~isempty(A_ajdk_));

if (~isempty(A_ajdk_) & ~isempty(An_ajdk_) & ~isempty(Zn_ajdk_));
for nca=0:At_ncols-1;
for ncz=0:Zt_ncols-1;
output_(1+nca + ncz*At_ncols) = dtmp_a0d1;
dtmp = popcount_xor_lf(uAn_(:,1+nca),uZn_(:,1+ncz),umc_,At_brows,A_ajdk_(1 + (0:An_pcols-1) + AJDK_0_1*An_pcols),popcount_,POPLENGTH);
output_(1+nca + ncz*At_ncols) = output_(1+nca + ncz*At_ncols) - 2*dtmp;
output_(1+nca + ncz*At_ncols) = output_(1+nca + ncz*At_ncols) - ( An_ajdk_(1+nca + AJDK_1_1*At_ncols) + Zn_ajdk_(1+ncz + AJDK_1_1*Zt_ncols) - dtmp_a2d1 );
end;%for ncz=0:Zt_ncols-1;
end;%for nca=0:At_ncols-1;
end;%if (~isempty(A_ajdk_) & ~isempty(An_ajdk_) & ~isempty(Zn_ajdk_));

if (~isempty(A_ajdk_) & isempty(An_ajdk_) & ~isempty(Zn_ajdk_));
for nca=0:At_ncols-1;
for ncz=0:Zt_ncols-1;
output_(1+nca + ncz*At_ncols) = dtmp_a0d1;
dtmp = popcount_xor_lf(uAn_(:,1+nca),uZn_(:,1+ncz),umc_,At_brows,A_ajdk_(1 + (0:An_pcols-1) + AJDK_0_1*An_pcols),popcount_,POPLENGTH);
output_(1+nca + ncz*At_ncols) = output_(1+nca + ncz*At_ncols) - 2*dtmp;
output_(1+nca + ncz*At_ncols) = output_(1+nca + ncz*At_ncols) - ( Zn_ajdk_(1+ncz + AJDK_1_1*Zt_ncols) - dtmp_a2d1 );
end;%for ncz=0:Zt_ncols-1;
end;%for nca=0:At_ncols-1;
end;%if (~isempty(A_ajdk_) & isempty(An_ajdk_) & ~isempty(Zn_ajdk_));

if (~isempty(A_ajdk_) & ~isempty(An_ajdk_) & isempty(Zn_ajdk_));
for nca=0:At_ncols-1;
for ncz=0:Zt_ncols-1;
output_(1+nca + ncz*At_ncols) = dtmp_a0d1;
dtmp = popcount_xor_lf(uAn_(:,1+nca),uZn_(:,1+ncz),umc_,At_brows,A_ajdk_(1 + (0:An_pcols-1) + AJDK_0_1*An_pcols),popcount_,POPLENGTH);
output_(1+nca + ncz*At_ncols) = output_(1+nca + ncz*At_ncols) - 2*dtmp;
output_(1+nca + ncz*At_ncols) = output_(1+nca + ncz*At_ncols) - ( An_ajdk_(1+nca + AJDK_1_1*At_ncols) - dtmp_a2d1 );
end;%for ncz=0:Zt_ncols-1;
end;%for nca=0:At_ncols-1;
end;%if (~isempty(A_ajdk_) & ~isempty(An_ajdk_) & isempty(Zn_ajdk_));

if (~isempty(A_ajdk_) & isempty(An_ajdk_) & isempty(Zn_ajdk_));
for nca=0:At_ncols-1;
for ncz=0:Zt_ncols-1;
output_(1+nca + ncz*At_ncols) = dtmp_a0d1;
dtmp = popcount_xor_lf(uAn_(:,1+nca),uZn_(:,1+ncz),umc_,At_brows,A_ajdk_(1 + (0:An_pcols-1) + AJDK_0_1*An_pcols),popcount_,POPLENGTH);
output_(1+nca + ncz*At_ncols) = output_(1+nca + ncz*At_ncols) - 2*dtmp;
output_(1+nca + ncz*At_ncols) = output_(1+nca + ncz*At_ncols) - ( - dtmp_a2d1 );
end;%for ncz=0:Zt_ncols-1;
end;%for nca=0:At_ncols-1;
end;%if (~isempty(A_ajdk_) & isempty(An_ajdk_) & isempty(Zn_ajdk_));

if (isempty(A_ajdk_));
tmp_ajdk_ = ones(An_pcols,1);
dtmp_a0d1 = 0;
for nc=0:An_bcols-1;
dtmp_a0d1 = dtmp_a0d1 + double(intlut(umc_(1+nc),popcount_));
end;%for nc=0:An_bcols-1;
for nca=0:At_ncols-1;
for ncz=0:Zt_ncols-1;
output_(1+nca + ncz*At_ncols) = dtmp_a0d1;
dtmp = popcount_xor_lf(uAn_(:,1+nca),uZn_(:,1+ncz),umc_,At_brows,tmp_ajdk_,popcount_,POPLENGTH);
output_(1+nca + ncz*At_ncols) = output_(1+nca + ncz*At_ncols) - 2*dtmp;
end;%for ncz=0:Zt_ncols-1;
end;%for nca=0:At_ncols-1;
end;%if (isempty(A_ajdk_));

output_ = reshape(output_,At_ncols,Zt_ncols);



