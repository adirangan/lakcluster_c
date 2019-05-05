function uAt_ = tutorial_integer_to_binary(bitj,An_nrows,An_ncols,iAt_);

if (nargin<1);
bitj = 16;
An_nrows = round(1920*2.5); 
An_ncols = 13; 
An_ = randn(An_nrows,An_ncols)>0; 
uAt_ = tutorial_integer_to_binary(bitj,An_nrows,An_ncols,An_);
iAt_ = tutorial_binary_to_integer(bitj,An_nrows,An_ncols,uAt_);
disp(sprintf('error |An_ - iAt_| = %f',norm(An_ - iAt_)));
return;
end;% if (nargin<1);

BIT8=8;
At_brows = rup(An_nrows,bitj)/BIT8;
An_nrows_extend = mod(bitj - mod(An_nrows,bitj),bitj);
input = iAt_;
input(An_nrows + (1:An_nrows_extend),[]) = zeros(An_nrows_extend,0);
input = input(:);
lout = length(input);
uAt_ = zeros(lout/BIT8,1,'uint8');
input = reshape(input,BIT8,lout/BIT8);
br = zeros(1,BIT8,'uint8'); br = 2.^(BIT8-1:-1:0);
uAt_ = uint8(br*input); uAt_=reshape(uAt_,At_brows,An_ncols);

