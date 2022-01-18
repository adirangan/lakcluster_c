function output = mda_read_i4(fname);
% read mda file written by mda_write_i4 in mda_io.c ;
fp = fopen(fname);
if (fp==-1); disp(sprintf(' %% Warning! cannot open %s in mda_read_r8',fname)); end;
n_d = fread(fp,1,'int32');
d_ = zeros(1,n_d);
dd = 1;
for nd=1:n_d;
d_(nd) = fread(fp,1,'int32');
dd = dd*d_(nd);
end;%for nd=1:n_d;
output = fread(fp,dd,'int32');
fp = fopen(fname);
output = reshape(output,d_);