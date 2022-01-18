function ...
[ ...
 flag_skip ...
,fname_ext ...
] = ...
open_fname_tmp( ...
 fname_pre ...
,date_diff_threshold ...
,flag_force_create_ext ...
,flag_force_create_tmp ...
,str_ext ...
);
na=0;
if (nargin<1+na); fname_pre=''; end; na=na+1;
if (nargin<1+na); date_diff_threshold=[]; end; na=na+1;
if (nargin<1+na); flag_force_create_ext=[]; end; na=na+1;
if (nargin<1+na); flag_force_create_tmp=[]; end; na=na+1;
if (nargin<1+na); str_ext=[]; end; na=na+1;
%%%%;
if ( isempty(date_diff_threshold)); date_diff_threshold = 1.5; end;
if ( isempty(flag_force_create_ext)); flag_force_create_ext = 0; end;
if ( isempty(flag_force_create_tmp)); flag_force_create_tmp = 0; end;
if ( isempty(str_ext)); str_ext = 'mat'; end;
%%%%;
fname_ext = sprintf('%s.%s',fname_pre,str_ext);
fname_tmp = sprintf('%s.tmp',fname_pre);
if ( exist(fname_ext,'file') & ~flag_force_create_ext);
  flag_skip=1;
end;%if ( exist(fname_ext,'file') & ~flag_force_create_ext);
if (~exist(fname_ext,'file') |  flag_force_create_ext);
flag_skip=0;
if ( exist(fname_tmp,'file'));
tmp_date_diff = datenum(clock) - datenum(dir(fname_tmp).date);
if (tmp_date_diff< date_diff_threshold & ~flag_force_create_tmp);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = recent, skipping',fname_tmp,tmp_date_diff));
flag_skip=1;
end;%if (tmp_date_diff< date_diff_threshold & ~flag_force_create_tmp);
if (tmp_date_diff>=date_diff_threshold |  flag_force_create_tmp);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = stale, deleting',fname_tmp,tmp_date_diff));
delete(fname_tmp);
flag_skip=0;
end;%if (tmp_date_diff>=date_diff_threshold |  flag_force_create_tmp);
end;%if ( exist(fname_tmp,'file'));
if (~flag_skip);
disp(sprintf(' %% %s not found, creating',fname_pre));
save(fname_tmp,'fname_ext');
end;%if (~flag_skip);
end;%if (~exist(fname_ext,'file') |  flag_force_create_ext);
