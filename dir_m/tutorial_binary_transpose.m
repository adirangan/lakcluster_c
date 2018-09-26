function tutorial_binary_transpose(filename_to_read,filename_to_write);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function tutorial_binary_transpose(filename_to_read,filename_to_write);
%
% This function transposes one binary (-1/+1) submatrix from the input file ;
% (named filename_to_read) to the output file (named filename_to_write). ;
%
% test by running with no arguments:
% i.e., >> tutorial_binary_transpose();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2;
disp(sprintf(' '));
disp(' testing tutorial_binary_transpose: ');
A_n = randn(3,12)>0;
tutorial_binary_compress(16,A_n,'tutorial_binary_transpose_test_n.b16');
A_t = transpose(A_n);
tutorial_binary_compress(16,A_t,'tutorial_binary_transpose_test_t.b16');
tutorial_binary_transpose('tutorial_binary_transpose_test_t.b16','tutorial_binary_transpose_test_x.b16')
[bitj,nrows,ncols] = tutorial_binary_getsize('tutorial_binary_transpose_test_x.b16');
A_x = tutorial_binary_uncompress('tutorial_binary_transpose_test_x.b16',1:nrows,1:ncols)>0;
disp(sprintf(' %% error: %0.16f',norm(A_n-A_x)));
return;
end;%if nargin<2;

[bitj,nrows,ncols] = tutorial_binary_getsize(filename_to_read);
A_n = tutorial_binary_uncompress(filename_to_read,1:nrows,1:ncols)>0;
A_t = transpose(A_n);
tutorial_binary_compress(bitj,A_t,filename_to_write);
