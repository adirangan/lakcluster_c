function str_xfix = test_dexcluster_multi_xfix_4(prefix,M,N,snr,n_cluster,n_rank,niteration);
str_M = sprintf('_M%d',M);
str_N = sprintf('_N%d',N);
str_snr = sprintf('_snr%.3d',round(100*snr));
str_n_cluster = sprintf('_n%d',n_cluster);
str_n_rank = ''; if ~isempty(n_rank); str_n_rank = sprintf('_r%d',n_rank); end;
str_ni = sprintf('_ni%d',niteration);
str_xfix =  sprintf('%s%s%s%s%s%s%s',prefix,str_M,str_N,str_snr,str_n_cluster,str_n_rank,str_ni);
