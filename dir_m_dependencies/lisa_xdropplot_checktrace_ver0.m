%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run lisa_xdropplot_initialize first. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

flag_found=1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% set names. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
nshuffle=0;
if (strcmp(flag_dex_vs_lak,'dex')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver10('dex',maf_lo_threshold,maf_hi_threshold,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle)); end;
if (strcmp(flag_dex_vs_lak,'lak')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver10('lak',maf_lo_threshold,maf_hi_threshold,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle)); end;
string_name_s0000 = string_name;
dir_out_s0000 = sprintf('%s/dir_%s',dir_out,string_name_s0000);
dir_out_trace = sprintf('%s/dir_trace',dir_out_s0000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% check for rdrop_a and cdrop_a. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fname_tmp = sprintf('%s/out_xdrop_a.txt',dir_out_s0000); flag_found = flag_found*fcheck(fname_tmp,0);
fname_tmp = sprintf('%s/out_xdrop_b.txt',dir_out_s0000); flag_found = flag_found*fcheck(fname_tmp,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% check for mds and fam and bim file, as well as BP pheno file. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
fcheck(sprintf('%s/mds.mat',dir_trunk),1);
fcheck(sprintf('%s/%s_bim.ext',dir__in,string_prefix),1);
fcheck(sprintf('%s/%s_fam.ext',dir__in,string_prefix),1);
fcheck(sprintf('%s/PGC_BIP32b-BD1.20160612.pheno',dir_trunk),1);
fcheck(sprintf('%s/PGC_BIP32b-BD2.20160612.pheno',dir_trunk),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% check for traces (original stored in s0000) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
trace_found_ = zeros(n_shuffle+1,1);
trace_length_ = zeros(n_shuffle+1,1);
for (nshuffle=0:n_shuffle);
if (strcmp(flag_dex_vs_lak,'dex')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver10('dex',maf_lo_threshold,maf_hi_threshold,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle)); end;
if (strcmp(flag_dex_vs_lak,'lak')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver10('lak',maf_lo_threshold,maf_hi_threshold,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,nshuffle)); end;
trace_found_(1+nshuffle)=0;
trace_length_(1+nshuffle)=0;
if (trace_found_(1+nshuffle)==0);
fname_tmp = sprintf('%s/dir_%s/out_trace.txt',dir_out,string_name);
if (exist(fname_tmp,'file')); 
if (verbose>1); disp(sprintf(' %% found fname_tmp: %s',fname_tmp)); end;
trace_found_(1+nshuffle)=1; end;
end;%if (trace_found_(1+nshuffle)==0);
if (trace_found_(1+nshuffle)==0);
fname_tmp = sprintf('%s/out_trace_s%0.4d.txt',dir_out_trace,nshuffle);
if (exist(fname_tmp,'file')); 
if (verbose>1); disp(sprintf(' %% found fname_tmp: %s',fname_tmp)); end;
trace_found_(1+nshuffle)=1; end;
end;%if (trace_found_(1+nshuffle)==0);
if (trace_found_(1+nshuffle)==1);
tmp_ = textread(fname_tmp);
trace_length_(1+nshuffle) = size(tmp_,1);
end;%if (trace_found_(1+nshuffle)==1);
end;%for (nshuffle=0:n_shuffle);
if (verbose>2); disp(sprintf(' %% found %d/%d traces',sum(trace_found_),1+n_shuffle)); end;
%%%%%%%%%%%%%%%%;
flag_found = flag_found*trace_found_(1);
fraction_found = sum(trace_found_(2:end))/(n_shuffle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% check for pca_proj ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
dname = sprintf('%s/dir_pca',dir_out_s0000);
flag_pca = exist(dname,'dir'); n_pca = 0;
if flag_pca;
tmp_ = dir([dname '/*.*']); tmp = size(tmp_,1);
n_pca = tmp;
end;%if flag_pca;

if (verbose); disp(sprintf(' %% %s: flag_found %d (length %d) fraction_found %0.2f (average length %f) flag_pca %d (number %d)',string_name_s0000,flag_found,trace_length_(1),fraction_found,mean(trace_length_(setdiff(find(trace_found_),[1]))),flag_pca,n_pca)); end;
