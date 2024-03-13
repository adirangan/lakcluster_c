function c = char_compress(input,filename_to_save);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function c = char_compress(input,filename_to_save);
%
% This function compresses a uint8 (0,...,255) array (input) ;
% into a filename (named filename_to_save). ;
% The output is c (an array) of uint8s that encode the input. ;
%
% test by running with no arguments:
% i.e., >> char_compress();
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2;
disp(sprintf(' %% '));
disp(sprintf(' %% testing char_compress: '));
A_r8 = floor(1024*rand(3,12))-512;
A_uint8 = cast(A_r8,'uint8');
char_compress(A_r8,'char_compress_test.char');
B = char_uncompress('char_compress_test.char',[1,3],[2,4,6,8,10]);
disp(sprintf(' %% creating array named A_r8: '));
disp(num2str(A_r8([1,3],[2,4,6,8,10])));
disp(sprintf(' %% casting as array named A_uint8: '));
disp(num2str(A_uint8([1,3],[2,4,6,8,10])));
disp(sprintf(' %% stored as char_compress_test.char. '));
disp(sprintf(' %% reading from this file a recovered array named B: '));
disp(num2str(B));
disp(sprintf(' %% error |A_uint8-B| = %f',norm(cast(A_uint8([1,3],[2,4,6,8,10])-B,'double'),'fro')));
return;
end;%if nargin<2;

%fcheck(filename_to_save);
fid = fopen(filename_to_save,'w');
[nrows,ncols] = size(input);
fwrite(fid,nrows,'int');
fwrite(fid,ncols,'int');
fwrite(fid,cast(input,'uint8'),'uint8');
fclose(fid);
