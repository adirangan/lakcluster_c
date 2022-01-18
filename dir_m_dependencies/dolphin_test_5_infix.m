function string_infix = dolphin_test_5_infix(rseed,n_var,dt_avg,n_T_0in,relative_variation,snr);
str_dt = sprintf('dt%.4d',floor(1000*dt_avg));
str_n_T = sprintf('T%.4d',floor(n_T_0in));
str_rv = sprintf('rv%.4d',floor(1000*relative_variation));
str_snr = sprintf('snr%.4d',floor(1000*snr));
string_infix = sprintf('dt5_rng%d_n%d_%s_%s_%s_%s',rseed,n_var,str_dt,str_n_T,str_rv,str_snr);
