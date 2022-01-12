function xxxcluster_fromdisk_uADZSZDA_excerpt_1(fp,str_0,n_bin,M_);
fprintf(fp,str_0);
for nb=0:n_bin-1; 
fprintf(fp,'%d',M_(1+nb)); 
if nb<n_bin-1; fprintf(fp,','); else fprintf(fp,';'); end; 
end; %for nb=0:n_bin-1; 
fprintf(fp,'\n');
