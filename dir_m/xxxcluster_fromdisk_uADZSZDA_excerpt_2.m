function [vname_,str_out__] = xxxcluster_fromdisk_uADZSZDA_excerpt_2(n_bin,str_0in_);
% breaks str_0in_ into a variable name vname_, ;
% as well as a cell-array of n_bin strings. ;

tmp_ij = min(strfind(str_0in_,'='));
vname_ = str_0in_(1:tmp_ij-1);
tmp_ij_ = [tmp_ij+1,strfind(str_0in_,','),strfind(str_0in_,';')];
assert(numel(tmp_ij_)==n_bin+1);
str_out__ = cell(n_bin,1);
for nbin=0:n_bin-1;
str_out__{1+nbin} = str_0in_(tmp_ij_(1+nbin+0)+1:tmp_ij_(1+nbin+1)-1);
end;%for nbin=0:n_bin-1;
