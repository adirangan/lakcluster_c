function c_MDA_write_r4(fname,A);
fp = fopen(fname,'w');
n_d = ndims(A);
d_ = size(A);
fwrite(fp,n_d,'int32');
fwrite(fp,d_,'int32');
fwrite(fp,cast(A,'single'),'float32');
fclose(fp);
