function [ memory_in_GB ] = memory_GB(varname)
% MEMORY_GB uses the WHOS command, ;
% evaluating the BASE workspace and summing up the bytes. ;
% Output displayed in GB. ;

if nargin<1; varname = []; end;

if ~isempty(varname);
memory_in_GB = 0;
mem_elements = evalin('base','whos');
i_found=0;
for i=1:size(mem_elements,1);
if strcmp(mem_elements(i).name,varname); i_found=i; end;
end;%for i=1:size(mem_elements,1);
if i_found>0;
memory_in_GB = mem_elements(i_found).bytes;
memory_in_GB = memory_in_GB/1024/1024/1024;
end;%if i_found>0;
end;%if ~isempty(varname);

if  isempty(varname);
mem_elements = evalin('base','whos');
if (size(mem_elements,1)> 0);
for i = 1:size(mem_elements,1);
memory_array(i) = mem_elements(i).bytes;
end;%for i = 1:size(mem_elements,1);
memory_in_GB = sum(memory_array);
memory_in_GB = memory_in_GB/1024/1024/1024;
 else;
memory_in_GB = 0;
end;%if (size(mem_elements,1)> 0);
end;%if  isempty(varname);
