function iAt_ = tutorial_binary_to_integer(bitj,An_nrows,An_ncols,uAt_);

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
An_brows = rup(An_nrows,bitj)/BIT8;
iAt_ = zeros(BIT8*An_brows,An_ncols);
for nb=0:An_brows-1;
iAt_( 1 + nb*BIT8 + (0:7) , : ) = transpose(dec2bin(uAt_(1+nb,:),8)-'0');
end;%for nb=0:brows-1;
iAt_ = iAt_(1:An_nrows,:);
