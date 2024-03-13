function output = char_uncompress(filename_to_read,ij_row_,ij_col_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function output = char_uncompress(filename_to_read,ij_row_,ij_col_);
%
% This function extracts a char (0,...,255) submatrix (output) ;
% from the stored file (named filename_to_read), ;
% assuming that this file was written with ;
% char_compress. ;
% The inputs are filename_to_read (a string) and ij_row_ and ij_col_ ;
% (each integer arrays) which indicate the row- and column-indices to read. ;
% Note that we expect ij_row_ to range from 1 to nrows, and for ij_col_ to range ;
% from 1 to ncols, where [nrows,ncols] is the size of the stored array. ;
% The output is the submatrix requested. (of size length(ij_row_)-by-length(ij_col_)). ;
%
% test by running with no arguments:
% i.e., >> char_uncompress();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1;
disp(sprintf(' %% '));
disp(sprintf(' %% testing char_uncompress: '));
A = cast(floor(256*rand(3,12)),'uint8');
char_compress(A,'char_compress_test.char');
B = char_uncompress('char_compress_test.char',[1,3],[2,4,6,8,10]);
disp(sprintf(' %% creating array named A: '));
disp(num2str(A([1,3],[2,4,6,8,10])));
disp(sprintf(' %% stored as char_compress_test.char. '));
disp(sprintf(' %% reading from this file a recovered array named B: '));
disp(num2str(B));
disp(sprintf(' %% error |A-B| = %f',norm(cast(A([1,3],[2,4,6,8,10])-B,'double'),'fro')));
return;
end;%if nargin<1;

if nargin<3;
[nrows,ncols] = char_getsize(filename_to_read);
output = char_uncompress(filename_to_read,1:nrows,1:ncols);
return;
end;% if nargin<3;

verbose=0;

fcheck(filename_to_read);
fid = fopen(filename_to_read,'r');
nrows = fread(fid,1,'int');
ncols = fread(fid,1,'int');
if nrows*ncols<=0;
output = zeros(nrows,ncols);
else;

nrows_extend = 0; %nrows_extend = mod(bitj - mod(nrows,bitj),bitj);
mr_length = (nrows + 0)/1; %mr_length = (nrows + nrows_extend)/BIT8;
[ij_col_s_,ij_col_j_] = sort(ij_col_-1);
[ij_col_r_,ij_col_i_] = sort(ij_col_j_);
b = zeros(numel(ij_row_),numel(ij_col_s_),'uint8');
nc=0;ncc=0; while (ij_col_s_(1+ncc)<0); ncc = ncc+1; end;%while;
while (nc<ncols & ncc<length(ij_col_s_));
while (nc<ncols & ncc<length(ij_col_s_) & nc<ij_col_s_(1+ncc));
if (verbose>1); disp(sprintf(' %% jumping column %d',nc)); end;
fseek(fid,mr_length,'cof');
nc = nc+1;
end;%while jump columns. ;
if (nc<ncols & ncc<length(ij_col_s_) & nc==ij_col_s_(1+ncc));
if (verbose); disp(sprintf(' %% reading column %d(%d)',nc,ncc)); end;
btmp = uint8(fread(fid,mr_length,'uint8'));
b(:,1+ncc) = btmp(ij_row_);
clear btmp;
nc = nc+1; ncc = ncc+1;
end;%if;
end;%while (nc<ncols);
output = b;
if (verbose); disp(sprintf('recovered %s of size %d,%d',filename_to_read,length(ij_row_),length(ij_col_))); end;
end;%if nrows*ncols<=0;

fclose(fid);
