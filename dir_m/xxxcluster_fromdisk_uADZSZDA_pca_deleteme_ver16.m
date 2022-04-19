function ...
[ ...
 parameter ...
] = ...
xxxcluster_fromdisk_uADZSZDA_pca_deleteme_ver16( ...
 parameter ...
);

str_thisfunction = 'xxxcluster_fromdisk_uADZSZDA_pca_deleteme_ver16';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose = 0; end;
flag_verbose = parameter.flag_verbose;

%%%%;
if ~isfield(parameter,'dir_0in'); disp(sprintf(' %% Warning, parameter.dir_0in undefined in %s',str_thisfunction)); end;
dir_0in = parameter.dir_0in;
str_suffix = sprintf('%s','analyze');
dir_tmp = sprintf('%s_%s',dir_0in,str_suffix);
if ~exist(dir_tmp,'dir'); disp(sprintf(' %% mkdir %s',dir_tmp)); mkdir(dir_tmp); end;
%%%%;
if ~isfield(parameter,'str_name_s0000');
parameter.str_name_s0000 = 'default';
end;%if ~isfield(parameter,'str_name_s0000');
str_name_s0000 = parameter.str_name_s0000;
%%%%;
if ~isfield(parameter,'dir_out_s0000');
str_suffix = sprintf('%s','analyze');
dir_out_s0000 = sprintf('%s_%s/dir_%s',dir_0in,str_suffix,str_name_s0000);
if ~exist(dir_out_s0000,'dir'); disp(sprintf(' %% mkdir %s',dir_out_s0000)); mkdir(dir_out_s0000); end;
parameter.dir_out_s0000 = dir_out_s0000;
end;%if ~isfield(parameter,'dir_out_s0000');
%%%%;
str_name_s0000 = parameter.str_name_s0000;
dir_out_s0000 = parameter.dir_out_s0000;
if (flag_verbose); disp(sprintf(' %% str_name_s0000: %s',str_name_s0000)); end;
if (flag_verbose); disp(sprintf(' %% dir_out_s0000: %s',dir_out_s0000)); end;

str_command = sprintf(' find %s -iname "*deleteme*" ',parameter.dir_out_s0000);
disp(str_command);
system(str_command);
%%%%%%%%
str_command = sprintf(' find %s -iname "*deleteme*" -exec rm {} + ',parameter.dir_out_s0000);
disp(str_command);
system(str_command);
