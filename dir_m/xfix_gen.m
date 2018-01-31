function [output_string,X_factor_d,X_esm_d,gamma_d] = xfix_gen(N,X_factor,X_esm,gamma,B_MLT,rng_num);

X_factor_d = floor(X_factor*100);
X_esm_d = 10 + floor(-10*log10(X_esm));
gamma_d = floor(gamma*100);
output_string = sprintf('N%.5d_Xf%.2d_Xe%.2d_g%.2d_B%.2d_r%.3d',N,X_factor_d,X_esm_d,gamma_d,B_MLT,rng_num);
