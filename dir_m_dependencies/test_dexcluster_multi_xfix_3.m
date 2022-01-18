function str_xfix = test_dexcluster_multi_xfix_3(prefix,nrank,M,N,snr,n_cluster,gamma,niteration);
str_nrank = sprintf('_r%d',nrank);
str_M = sprintf('_M%d',M);
str_N = sprintf('_N%d',N);
str_snr = sprintf('_snr%.3d',round(100*snr));
str_n_cluster = sprintf('_n%d',n_cluster);
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
str_ni = sprintf('_ni%d',niteration);
str_xfix =  sprintf('%s%s%s%s%s%s%s',prefix,str_nrank,str_M,str_N,str_snr,str_n_cluster,str_g,str_ni);
