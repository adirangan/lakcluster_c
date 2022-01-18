%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% run lisa_setprefix_ver2 first. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
% set names. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
nshuffle=0;
if (strcmp(flag_dex_vs_lak,'dex')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver13('dex',maf_lo_threshold,maf_hi_threshold,mc_string,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,nshuffle)); end;
if (strcmp(flag_dex_vs_lak,'lak')); string_name = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver13('lak',maf_lo_threshold,maf_hi_threshold,mc_string,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,nshuffle)); end;
string_name_s0000 = string_name;
dir_out_s0000 = sprintf('%s/dir_%s',dir_out,string_name_s0000);
dir_out_trace = sprintf('%s/dir_trace',dir_out_s0000);
nshuffle=0;
if (strcmp(flag_dex_vs_lak,'dex')); string_name_r0 = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver13('dex',maf_lo_threshold,maf_hi_threshold,mc_string,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,0,nshuffle)); end;
if (strcmp(flag_dex_vs_lak,'lak')); string_name_r0 = sprintf('%s_%s',string_prefix,xxxcluster_PGC_uADZSZDA_xfix_ver13('lak',maf_lo_threshold,maf_hi_threshold,mc_string,flag_reverse,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,0,nshuffle)); end;
string_name_s0000_r0 = string_name_r0;
dir_out_s0000_r0 = sprintf('%s/dir_%s',dir_out,string_name_s0000_r0);
dir_out_trace_r0 = sprintf('%s/dir_trace',dir_out_s0000_r0);
