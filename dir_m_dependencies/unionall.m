function output = unionall(input);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function output = unionall(input);
% 
% This simple function takes a cell-array of input and outputs ;
% the union of each cell. ;
% Input = cell array. ;
% Output = union of each cell. ;
%
% test by running with no arguments:
% i.e., >> unionall();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1;
disp(sprintf(' '));
disp(sprintf(' testing unionall: '));
disp(sprintf(' First we generate a simple cell array: '));
disp(sprintf(' C = {[1 2],[1 3],[3 6 8],[8 1]};'));
C = {[1 2],[1 3],[3 6 8],[8 1]};
disp(sprintf(' and then take its union: '));
disp(sprintf(' D = unionall(C); '));
disp(sprintf(' D = '));
D = unionall(C);
disp(num2str(transpose(D)));
return;
end;%if nargin<1;

l0=length(input); output = [];
for nl=1:l0; output = union(output,input{nl}); end;
