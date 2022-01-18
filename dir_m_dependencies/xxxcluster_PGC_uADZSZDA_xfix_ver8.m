function [output_string,gamma_d] = xxxcluster_PGC_uADZSZDA_xfix_ver8(prefix,p_threshold,rev_flag,n_mds,mds_used,mds_repl,gamma,B_MLT,Ireq,nshuffle);
p_str = sprintf('p%.2d',floor(100*p_threshold));
if (rev_flag==1); rev_str = 'X'; else rev_str = 'D'; end;
mds_tmp = zeros(1,n_mds); mds_tmp(mds_used)=1; mds_code = dot(2.^[0:n_mds-1],mds_tmp);
gamma_d = floor(gamma*1000);
output_string = sprintf('%s_uADZSZDA_%s_%s_m%dr%d_g%.3d_B%.2d_I%d_s%.4d',prefix,p_str,rev_str,mds_code,mds_repl,gamma_d,B_MLT,Ireq,nshuffle);
