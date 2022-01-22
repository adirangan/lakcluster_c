function n = wc_0(fname);
n = 0; 
if (~fcheck(fname)); return; end;
fid = fopen(fname,'r');
tline = fgetl(fid);
while ischar(tline); tline = fgetl(fid); n = n+1; end;
fclose(fid);
