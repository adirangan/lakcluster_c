function [output_string,gamma_d] = xxxcluster_PGC_uADZSZDA_xfix_ver12(prefix,maf_lo_threshold,maf_hi_threshold,rev_flag,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,nshuffle);

if nargin<12;
n_mds = 2; mds_used_ = [1:2];
for maf_lo_threshold = [0.01,0.25];
for maf_hi_threshold = [0.15,0.50];
for rev_flag = [0:1];
for mds_repl = [0,3];
for gamma = [0.004 , 0.05];
for B_MLT = [24,34];
for Ireq = [0,5];
for n_scramble = [0,1];
for nshuffle = [0,8];
disp(xxxcluster_PGC_uADZSZDA_xfix_ver12('dex',maf_lo_threshold,maf_hi_threshold,rev_flag,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,nshuffle));
disp(xxxcluster_PGC_uADZSZDA_xfix_ver12('lak',maf_lo_threshold,maf_hi_threshold,rev_flag,n_mds,mds_used_,mds_repl,gamma,B_MLT,Ireq,n_scramble,nshuffle));
end;end;end;end;end;end;end;end;end;
disp('returning');return;
end;%if nargin<12;

if max(maf_lo_threshold,maf_hi_threshold)>=0.5; 
maf_str = sprintf('_p%.2d',floor(100*min(maf_lo_threshold,maf_hi_threshold)));
else; 
maf_str = sprintf('_p%.2dq%.2d',floor(100*min(maf_lo_threshold,maf_hi_threshold)),floor(100*max(maf_lo_threshold,maf_hi_threshold))); 
end;%if max(maf_lo_threshold,maf_hi_threshold)>=0.5; 
if (rev_flag==1); rev_str = '_X'; else rev_str = '_D'; end;
mds_tmp_ = zeros(1,n_mds); 
mds_tmp_(mds_used_)=1; 
mds_code = length(mds_used_); 
%mds_code = dot(2.^[0:n_mds-1],mds_tmp_);
if (mds_repl<1 | length(mds_used_)<2); mds_str = sprintf('_m%d',mds_code); else; mds_str = sprintf('_m%dr%d',mds_code,mds_repl); end;
gamma_d = floor(gamma*1000); gamma_str = sprintf('_g%.3d',gamma_d);
if (strcmp(prefix,'dex')); B_str = ''; else; B_str = sprintf('_B%.2d',B_MLT); end;
if (Ireq<=0); Ireq_str = ''; else; Ireq_str = sprintf('_I%d',Ireq); end;
if (n_scramble<=0); scramble_str = ''; else; scramble_str = sprintf('_r%d',n_scramble); end;
if (nshuffle<=0); shuffle_str = ''; else; shuffle_str = sprintf('_s%.4d',nshuffle); end;
output_string = sprintf('%s%s%s%s%s%s%s%s%s',prefix,maf_str,rev_str,mds_str,gamma_str,B_str,Ireq_str,scramble_str,shuffle_str);

