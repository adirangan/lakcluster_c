function mda_compress(ra,fname);
% create mda file: ;
%{
  A = randn(7,3,5);
  mda_compress(A,'tmp.mda');
  B = mda_uncompress('tmp.mda');
  disp(sprintf(' %% error %f',norm(A(:)-B(:))));
  %}
s = size(ra); ndims = length(s);
fp = fopen(fname,'w');
fwrite(fp,ndims,'int');
for nd=1:ndims; fwrite(fp,s(nd),'int'); end;
fwrite(fp,ra(:),'double');
fclose(fp);
