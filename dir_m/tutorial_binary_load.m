function [output,nrows,brows,ncols] = tutorial_binary_load(filename_to_read);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function output = tutorial_binary_load(filename_to_read);
%
% This function extracts a binary (-1/+1) submatrix (output) from the stored file ;
% (named filename_to_read), assuming that this file was written with ;
% tutorial_binary_compress. ;
% The inputs are filename_to_read (a string);
%
% test by running with no arguments:
% i.e., >> tutorial_binary_load();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1;
disp(sprintf(' '));
disp(' testing tutorial_binary_load: ');
A = randn(11,12)>0;
tutorial_binary_compress(16,A,'tutorial_binary_compress_test.b16');
B = tutorial_binary_load('tutorial_binary_compress_test.b16');
disp(sprintf(' creating array named A: '));
disp(num2str(A));
disp(sprintf(' stored as tutorial_binary_compress_test.b16. '));
disp(sprintf(' reading from this file a recovered array named B: '));
disp(num2str(B));
B2 = zeros(16,12);
B2(1:8,:) = transpose(dec2bin(B(1,:),8)-'0');
B2(9:16,:) = transpose(dec2bin(B(2,:),8)-'0');
B2 = B2(1:11,:);
disp(sprintf(' uncompressing B --> B2: '));
disp(num2str(B2));
disp(sprintf('error |A-B2| = %f',norm(A-B2)));
return;
end;%if nargin<3;

verbose=1;
BIT8=8;

if (~exist(filename_to_read,'file')); disp(sprintf(' %% Warning! %s does not exist',filename_to_read)); end;
fid = fopen(filename_to_read,'r');
bitj = fread(fid,1,'int');
nrows = fread(fid,1,'int');
ncols = fread(fid,1,'int');
if nrows*ncols<=0; % empty output ;
if (verbose); disp(sprintf(' %% %s empty',filename_to_read)); end;
output = zeros(nrows,ncols);
 else; % read something ;
nrows_extend = mod(bitj - mod(nrows,bitj),bitj);
brows = (nrows + nrows_extend)/BIT8;
output = zeros(1,brows*ncols,'uint8');
output = reshape(uint8(fread(fid,brows*ncols,'uint8')),brows,ncols);
if (verbose); disp(sprintf(' %% %s: read %d(%d)-x-%d',filename_to_read,nrows,brows,ncols)); end;
end; %if nrows*ncols<=0;
fclose(fid);
