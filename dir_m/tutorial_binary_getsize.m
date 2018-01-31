function [j,nrows,ncols] = tutorial_binary_getsize(filename_to_read);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [j,nrows,ncols] = tutorial_binary_getsize(bitj,input,filename_to_save);
%
% This function get the size of a binary array (stored in filename_to read) ;
% The inputs are filename_to_read (a string). ;
% The output is j (i.e., bitj=16), nrows and ncols (all integers). ;
%
% test by running with no arguments:
% i.e., >> tutorial_binary_getsize();
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1;
disp(sprintf(' '));
disp(' testing tutorial_binary_getsize: ');
A = randn(3,12)>0;
tutorial_binary_compress(16,A,'./tutorial_binary_getsize_test.b16');
disp(sprintf(' creating %dx%d array named A: ',3,12));
disp(sprintf(' stored as ./tutorial_binary_getsize_test.b16. '));
[bitj,nr,nc] = tutorial_binary_getsize('./tutorial_binary_getsize_test.b16');
disp(sprintf(' reading from this file sizes %dx%d.',nr,nc));
return;
end;%if nargin<1;

fid = fopen(filename_to_read,'r');
j = fread(fid,1,'int');
nrows = fread(fid,1,'int');
ncols = fread(fid,1,'int');
fclose(fid);
