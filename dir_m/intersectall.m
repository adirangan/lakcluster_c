function output = intersectall(input);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function output = intersectall(input);
% 
% This simple function takes a cell-array of input and outputs ;
% the intersection of all the cells. ;
% Input = cell array. ;
% Output = intersection of all the cells. ;
%
% test by running with no arguments:
% i.e., >> intersectall();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1;
disp(sprintf(' '));
disp(sprintf(' testing intersectall: '));
disp(sprintf(' First we generate a simple cell array: '));
disp(sprintf(' C = {[1 3 2],[1 2 3],[3 6 1 8],[8 1 3]};'));
C = {[1 3 2],[1 2 3],[3 6 1 8],[8 1 3]};
disp(sprintf(' and then take its intersection: '));
disp(sprintf(' D = intersectall(C); '));
disp(sprintf(' D = '));
D = intersectall(C);
disp(num2str((D)));
return;
end;%if nargin<1;

l0=length(input); output = [];
for nl=1:l0; if (nl==1); output = input{nl}; else output = intersect(output,input{nl}); end; end;
