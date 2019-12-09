function output = binary_uncompress(filename_to_read,row_ind,col_ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function output = binary_uncompress(filename_to_read,row_ind,col_ind);
%
% This function extracts a binary (-1/+1) submatrix (output) from the stored file ;
% (named filename_to_read), assuming that this file was written with ;
% binary_compress. ;
% The inputs are filename_to_read (a string) and row_ind and col_ind ;
% (each integer arrays) which indicate the row- and column-indices to read. ;
% Note that we expect row_ind to range from 1 to nrows, and for col_ind to range ;
% from 1 to ncols, where [nrows,ncols] is the size of the stored array. ;
% The output is the submatrix requested. (of size length(row_ind)-by-length(col_ind)). ;
%
% test by running with no arguments:
% i.e., >> binary_uncompress();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1;
disp(sprintf(' '));
disp(' testing binary_uncompress: ');
A = randn(3,12)>0;
binary_compress(16,A,'binary_compress_test.b16');
B = binary_uncompress('binary_compress_test.b16',[1,3],[2,4,6,8,10]);
disp(sprintf(' creating array named A: '));
disp(num2str(A([1,3],[2,4,6,8,10])));
disp(sprintf(' stored as binary_compress_test.b16. '));
disp(sprintf(' reading from this file a recovered array named B: '));
disp(num2str(B>0));
disp(sprintf('error |A-B| = %f',norm(A([1,3],[2,4,6,8,10])-(B>0))));
return;
end;%if nargin<1;

if nargin<3;
[~,nrows,ncols] = binary_getsize(filename_to_read);
output = binary_uncompress(filename_to_read,1:nrows,1:ncols);
return;
end;% if nargin<3;

verbose=0;
BIT8=8;

fcheck(filename_to_read);
fid = fopen(filename_to_read,'r');
bitj = fread(fid,1,'int');
nrows = fread(fid,1,'int');
ncols = fread(fid,1,'int');
if nrows*ncols<=0;
output = zeros(nrows,ncols);
else;

nrows_extend = mod(bitj - mod(nrows,bitj),bitj);
mr_length = (nrows + nrows_extend)/BIT8;
[col_ind_s,col_ind_j] = sort(col_ind-1);
[col_ind_r,col_ind_i] = sort(col_ind_j);
b = zeros(1,mr_length*length(col_ind_s),'uint8');
nc=0;ncc=0; while (col_ind_s(1+ncc)<0); ncc = ncc+1; end;%while;
while (nc<ncols & ncc<length(col_ind_s));
while (nc<ncols & ncc<length(col_ind_s) & nc<col_ind_s(1+ncc)); if (verbose>1); disp(sprintf(' %% jumping column %d',nc)); end; fseek(fid,mr_length,'cof'); nc = nc+1; end;%while;
if (nc<ncols & ncc<length(col_ind_s) & nc==col_ind_s(1+ncc)); if (verbose); disp(sprintf(' %% reading column %d(%d)',nc,ncc)); end; btmp = uint8(fread(fid,mr_length,'uint8')); b(1 + (ncc*mr_length : (ncc+1)*mr_length-1)) = btmp; nc = nc+1; ncc = ncc+1; end;%if;
end;%while (nc<ncols);

br = cast(zeros(1,BIT8),'uint8'); br = 2.^(BIT8-1:-1:0);
br2 = cast(zeros(1,BIT8),'uint8'); br2 = 2.^(BIT8:-1:1);
output = zeros(mr_length*BIT8*length(col_ind_s),1);
for nl=0:mr_length*length(col_ind_s)-1; 
if (mod(nl,100)==0 & verbose); disp(sprintf('\r %% nl %d',nl)); end;%if (verbose);
output(1 + (0:BIT8-1) + nl*BIT8) = transpose(mod(double(b(1+nl)),br2)>=br); 
end;%for nl=0:mr_length*length(col_ind_s)-1; 
output = reshape(output,nrows+nrows_extend,length(col_ind_s));
output = output(row_ind,col_ind_i);
if (verbose); disp(sprintf('recovered %s of size %d,%d',filename_to_read,length(row_ind),length(col_ind))); end;
output = 2*(output)-1;

end;%if nrows*ncols<=0;

fclose(fid);
