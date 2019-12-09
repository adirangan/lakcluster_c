function b = binary_compress(bitj,input,filename_to_save);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function b = binary_compress(bitj,input,filename_to_save);
%
% This function compresses a binary (0/1) array (input) into a filename ;
% (named filename_to_save), ensuring that the number of rows stored is ;
% a multiple of bitj. ;
% The inputs are bitj (an integer), input (an array of 0s and 1s), ;
% and filename_to_save (a string). ;
% The output is b (an array) of uint8s that encode the input. ;
%
% test by running with no arguments:
% i.e., >> binary_compress();
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3;
disp(sprintf(' '));
disp(' testing binary_compress: ');
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
end;%if nargin<3;

%fcheck(filename_to_save);
fid = fopen(filename_to_save,'w');
fwrite(fid,bitj,'int');
[nrows,ncols] = size(input);
fwrite(fid,nrows,'int');
fwrite(fid,ncols,'int');
nrows_extend = mod(bitj - mod(nrows,bitj),bitj);
input(nrows + (1:nrows_extend),[]) = zeros(nrows_extend,0);
input = input(:);
lout = length(input);
bit8=8;
b = zeros(lout/bit8,1,'uint8');
input = reshape(input,bit8,lout/bit8);
br = zeros(1,bit8,'uint8'); br = 2.^(bit8-1:-1:0);
b = br*input; b=b(:);
fwrite(fid,b,'uint8');
fclose(fid);
