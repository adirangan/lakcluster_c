function ...
[ ...
 parameter ...
, ...
] = ...
function_template( ...
 parameter ...
, ...
);

str_thisfunction = 'function_template';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));

disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

