function close_fname_tmp(fname_pre);
na=0;
if (nargin<1+na); fname_pre=''; end; na=na+1;
%%%%;
fname_tmp = sprintf('%s.tmp',fname_pre);
delete(fname_tmp);
