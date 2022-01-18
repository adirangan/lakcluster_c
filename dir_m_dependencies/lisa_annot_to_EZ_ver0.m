function lisa_annot_to_EZ_ver0(fname_annot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% sparse matrix storing rs for each EZ. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fcheck(fname_annot);
n_lines = wc_0(fname_annot); n_header = 2; n_lines = n_lines - n_header;
fid = fopen(fname_annot,'r'); 
for nl=1:n_header; tline = fgetl(fid); end;%for nl=1:n_header;
EZ_ = cell(n_lines,1);
bp_ = cell(n_lines,1);
rs__ = cell(n_lines,1);
for nl=1:n_lines;
if (mod(nl,500)==0); disp(sprintf(' %% nl %d/%d',nl,n_lines)); end;
tline = fgetl(fid);
tmp_ = strsplit(tline);
tmp_len = length(tmp_);
EZ_{nl} = tmp_{1};
bp_{nl} = tmp_{2};
rs__{nl} = tmp_(3:end);
end;%for nl=1:n_lines;
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% count number of uniqe rs. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_rs_all = 0;
for nl=1:n_lines;
n_rs_all = n_rs_all + length(rs__{nl});
end;%for nl=1:n_lines;
rs_all_ = cell(n_rs_all,1);
nrs_all = 0;
for nl=1:n_lines;
tmp_len = length(rs__{nl});
rs_all_(nrs_all + (1:tmp_len)) = rs__{nl};
nrs_all = nrs_all + tmp_len;
end;%for nl=1:n_lines;
[rs_sort_,ij_all_,ij_sort_] = unique(rs_all_);
n_rs_sort = length(rs_sort_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% transpose sparse matrix to store EZ for each rs ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
ln_ = zeros(n_rs_sort,1);
EZ__ = cell(n_rs_sort,1);
nrs_all=0;
for nl=1:n_lines;
if (mod(nl,500)==0); disp(sprintf(' %% nl %d/%d',nl,n_lines)); end;
tmp_len = length(rs__{nl});
for nll=1:tmp_len;
nrs_all = nrs_all + 1;
assert(strcmp(rs__{nl}{nll},rs_all_{nrs_all}));
assert(strcmp(rs_all_{nrs_all},rs_sort_{ij_sort_(nrs_all)}));
ln_(ij_sort_(nrs_all)) = ln_(ij_sort_(nrs_all)) + 1;
EZ__{ij_sort_(nrs_all)}(ln_(ij_sort_(nrs_all))) = {0};
EZ__{ij_sort_(nrs_all)}{ln_(ij_sort_(nrs_all))} = EZ_{nl};
end;%for nll=1:tmp_len;
end;%for nl=1:n_lines;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% check sparse transpose;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
flag_check=0;
if flag_check;
disp(sprintf(' %% checking'));
for nrs_sort=1:n_rs_sort;
if (mod(nrs_sort,500)==0); disp(sprintf(' %% nrs_sort %d/%d',nrs_sort,n_rs_sort)); end;
tmp_rs = rs_sort_{nrs_sort};
tmp_len = length(EZ__{nrs_sort});
for nll=1:tmp_len;
tmp_ij = find(strcmp(EZ_,EZ__{nrs_sort}{nll}));
tmp_rs_ = rs__{tmp_ij};
assert(length(find(strcmp(tmp_rs_,tmp_rs)))>0);
end;%for nll=1:tmp_len;
end;%for nrs_sort=1:n_rs_sort;
end;%if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% save as mat file. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_mat = sprintf('%s.mat',fname_annot);
save(...
     fname_mat...
     ,'EZ_'...
     ,'EZ__'...
     ,'bp_'...
     ,'fname_annot'...
     ,'rs_sort_'...
     ,'ij_all_'...
     ,'ij_sort_'...
     ,'ln_'...
     ,'n_header'...
     ,'n_lines'...
     ,'n_rs_sort'...
     );
