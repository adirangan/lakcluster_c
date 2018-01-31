function [output_string,gamma_d] = lakcluster_uADZSWDY_xfix_gen_ver0(rev_flag,gamma,B_MLT,Ireq,shuffle_num);
if (rev_flag==1); rev_str = 'X'; else rev_str = 'D'; end;
gamma_d = floor(gamma*1000);
output_string = sprintf('lak_uADZSWDY_%s_g%.3d_B%.2d_I%d_s%.4d',rev_str,gamma_d,B_MLT,Ireq,shuffle_num);
