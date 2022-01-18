function [output] = mda_uncompress(fname);
% open mda file: ;
%{
  A = randn(7,3,5);
  mda_compress(A,'tmp.mda');
  B = mda_uncompress('tmp.mda');
  disp(sprintf(' %% error %f',norm(A(:)-B(:))));
  %}
fp = fopen(fname);
ndims = fread(fp,1,'int');
d_ = zeros(1,ndims); length=1;
for nd=1:ndims;
d_(nd) = fread(fp,1,'int');
length = length*d_(nd);
end;%for nd=1:ndims;
output = reshape(fread(fp,length,'double'),d_);
fclose(fp);
