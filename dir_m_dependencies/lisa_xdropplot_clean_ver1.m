%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run lisa_xdropplot_initialize first. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% extract rdrop_a and cdrop_a. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
nshuffle=0;
if (strcmp(flag_dex_vs_lak,'dex')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver13('dex',mc_string,maf_lo_threshold,maf_hi_threshold,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,nshuffle)); end;
if (strcmp(flag_dex_vs_lak,'lak')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver13('lak',mc_string,maf_lo_threshold,maf_hi_threshold,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,nshuffle)); end;
string_name_s0000 = string_name;
dir_out_s0000 = sprintf('%s/dir_%s',dir_out,string_name_s0000);
dir_out_trace = sprintf('%s/dir_trace',dir_out_s0000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% extract traces (original stored in s0000) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
for (nshuffle=1:n_shuffle);
if (strcmp(flag_dex_vs_lak,'dex')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver13('dex',mc_string,maf_lo_threshold,maf_hi_threshold,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,nshuffle)); end;
if (strcmp(flag_dex_vs_lak,'lak')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver13('lak',mc_string,maf_lo_threshold,maf_hi_threshold,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,nshuffle)); end;
dirname = sprintf('%s/dir_%s',dir_out,string_name);
if (exist(dirname,'dir')); disp(sprintf(' %% removing %s',dirname)); cmd_str = sprintf('rm -rf %s',dirname); disp(cmd_str); system(cmd_str); end;
end;%for (nshuffle=1:n_shuffle);
