function p_abs = abs_p(tmp_upb,p_set,gamma,n_iteration,r_rem_,c_rem_,ZR_s0000_,n_shuffle,tr_top_,tr_avg_);
p_s0000 = tmp_p(gamma,n_iteration,r_rem_,c_rem_,min(tmp_upb,ZR_s0000_),n_shuffle,tr_top_,tr_avg_);
p_abs = abs(p_s0000 - p_set);
