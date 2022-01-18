function mda_write_d3_r8(ra,fname);
% create mda file: ;
%{
  A = randn(7,3);
  mda_write_d3_r8(A,'tmp.mda');
  B = mda_uncompress('tmp.mda');
  disp(sprintf(' %% error %f',norm(A(:)-B(:))));
  A = randn(7,3,2);
  mda_write_d3_r8(A,'tmp.mda');
  B = mda_uncompress('tmp.mda');
  disp(sprintf(' %% error %f',norm(A(:)-B(:))));
  A = randn(7,3,2,4);
  mda_write_d3_r8(A,'tmp.mda');
  B = mda_uncompress('tmp.mda');
  disp(sprintf(' %% error %f',norm(A(:)-B(:))));
  %}
s = size(ra); ndims = length(s);
fp = fopen(fname,'w');
fwrite(fp,max(3,ndims),'int');
for nd=1:ndims; fwrite(fp,s(nd),'int'); end;
for nd=ndims+1:3; fwrite(fp,1,'int'); end;
fwrite(fp,ra(:),'double');
fclose(fp);
