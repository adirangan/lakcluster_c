function str_xfix = test_dexcluster_xfix_0(prefix,M,N,gamma,X,mu);
str_M = sprintf('_M%d',M);
str_N = sprintf('_N%d',N);
gamma_d = floor(gamma*1000); str_g = sprintf('_g%.3d',gamma_d);
X_d = floor(X*100); str_X = sprintf('_X%.3d',X_d);
mu_d = floor(mu*100); str_mu = sprintf('_mu%.3d',mu_d);
str_xfix =  sprintf('%s%s%s%s%s%s',prefix,str_M,str_N,str_g,str_X,str_mu);
