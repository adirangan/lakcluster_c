function l = lisa_struct_xfix_ver0(l)

if max(l.maf_lo_threshold,l.maf_hi_threshold)>=0.5; 
maf_str = sprintf('_p%.2d',floor(100*min(l.maf_lo_threshold,l.maf_hi_threshold)));
else; 
maf_str = sprintf('_p%.2dq%.2d',floor(100*min(l.maf_lo_threshold,l.maf_hi_threshold)),floor(100*max(l.maf_lo_threshold,l.maf_hi_threshold))); 
end;%if max(l.maf_lo_threshold,l.maf_hi_threshold)>=0.5; 
if (length(l.mr_string)>0);
mr_str = sprintf('_%s',l.mr_string);
else;
mr_str = '';
end;%if (length(l.mr_string)>0);
if (length(l.mc_string)>0);
mc_str = sprintf('_%s',l.mc_string);
else;
mc_str = '';
end;%if (length(l.mc_string)>0);
if (l.flag_reverse==1); rev_str = '_X'; else rev_str = '_D'; end;
mds_tmp_ = zeros(1,l.n_mds); 
mds_tmp_(l.mds_used_)=1; 
mds_code = length(l.mds_used_); 
%mds_code = dot(2.^[0:n_mds-1],mds_tmp_);
if (l.mds_repl<1 | length(l.mds_used_)<2); mds_str = sprintf('_m%d',mds_code); else; mds_str = sprintf('_m%dr%d',mds_code,l.mds_repl); end;
gamma_d = floor(l.gamma*1000); gamma_str = sprintf('_g%.3d',gamma_d);
if (strcmp(l.flag_dex_vs_lak,'dex')); B_str = ''; else; B_str = sprintf('_B%.2d',l.B_MLT); end;
if (l.Ireq<=0); Ireq_str = ''; else; Ireq_str = sprintf('_I%d',l.Ireq); end;
if (l.n_scramble<=0); scramble_str = ''; else; scramble_str = sprintf('_r%d',l.n_scramble); end;
if (l.nshuffle<=0); shuffle_str = ''; else; shuffle_str = sprintf('_s%.4d',l.nshuffle); end;
l.output_string = sprintf('%s%s%s%s%s%s%s%s%s%s%s',l.flag_dex_vs_lak,maf_str,mr_str,mc_str,rev_str,mds_str,gamma_str,B_str,Ireq_str,scramble_str,shuffle_str);
scramble_str = '';
l.output_string_r0 = sprintf('%s%s%s%s%s%s%s%s%s%s%s',l.flag_dex_vs_lak,maf_str,mr_str,mc_str,rev_str,mds_str,gamma_str,B_str,Ireq_str,scramble_str,shuffle_str);
if (l.n_scramble<=0); scramble_str = ''; else; scramble_str = sprintf('_r%d',l.n_scramble); end;
shuffle_str = '';
l.output_string_s0 = sprintf('%s%s%s%s%s%s%s%s%s%s%s',l.flag_dex_vs_lak,maf_str,mr_str,mc_str,rev_str,mds_str,gamma_str,B_str,Ireq_str,scramble_str,shuffle_str);
scramble_str = '';
l.output_string_s0_r0 = sprintf('%s%s%s%s%s%s%s%s%s%s%s',l.flag_dex_vs_lak,maf_str,mr_str,mc_str,rev_str,mds_str,gamma_str,B_str,Ireq_str,scramble_str,shuffle_str);

