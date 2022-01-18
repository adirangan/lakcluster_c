function [ memory_in_GB ] = monitor_memory_whos( )
%MONITOR_MEMORY_WHOS uses the WHOS command and evaluates inside the BASE
%workspace and sums up the bytes.  The output is displayed in GB.

mem_elements = evalin('base','whos');
if size(mem_elements,1) > 0

    for i = 1:size(mem_elements,1)
        memory_array(i) = mem_elements(i).bytes;
    end

    memory_in_GB = sum(memory_array);
    memory_in_GB = memory_in_GB/1024/1024/1024;
else
    memory_in_GB = 0;
end
