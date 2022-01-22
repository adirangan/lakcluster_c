function flag = fcheck(fname,verbose);
if (nargin<2); verbose=1; end;
flag=1;
if (~exist(fname,'file')); flag=0; if (verbose); disp(sprintf(' %% Warning! could not find file %s',fname)); end; end;
