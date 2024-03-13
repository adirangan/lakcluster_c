function [nrows,ncols] = char_getsize(filename_to_read);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [nrows,ncols] = char_getsize(filename_to_save);
%
% This function gets the size of a char array (stored in filename_to read) ;
% The inputs are filename_to_read (a string). ;
% The output is nrows and ncols (both integers). ;
%
% test by running with no arguments:
% i.e., >> char_getsize();
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1;
disp(sprintf(' %% '));
disp(sprintf(' %% testing char_getsize: '));
A = cast(floor(256*rand(3,12)),'uint8');
char_compress(A,'char_getsize_test.char');
disp(sprintf(' %% creating %dx%d array named A: ',3,12));
disp(sprintf(' %% stored as ./char_getsize_test.char. '));
[nr,nc] = char_getsize('./char_getsize_test.char');
disp(sprintf(' %% reading from this file sizes %dx%d.',nr,nc));
return;
end;%if nargin<1;

fcheck(filename_to_read);
fid = fopen(filename_to_read,'r');
nrows = fread(fid,1,'int');
ncols = fread(fid,1,'int');
fclose(fid);
