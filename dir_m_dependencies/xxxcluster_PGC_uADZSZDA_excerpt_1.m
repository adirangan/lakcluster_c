function xxxcluster_PGC_uADZSZDA_excerpt_1(fp,str_0,Icat,M_);
fprintf(fp,str_0);
for nb=0:Icat-1; 
fprintf(fp,'%d',M_(1+nb)); 
if nb<Icat-1; fprintf(fp,','); else fprintf(fp,';'); end; 
end; %for nb=0:Icat-1; 
fprintf(fp,'\n');
