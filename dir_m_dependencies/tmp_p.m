function p_s0000 = tmp_p(gamma,n_iteration,r_rem_,c_rem_,ZR_s0000_,n_shuffle,tr_top_,tr_avg_);

iteration_alo = max(1,floor(n_iteration*gamma));
iteration_ahi = min(n_iteration-1,floor(n_iteration*(1-gamma)));
iteration_avg_ = 1:iteration_ahi; 
dr_ = diff(r_rem_); dc_ = diff(c_rem_); 
dr2_ = transpose(dr_(iteration_avg_))/sum(dr_(iteration_avg_)); 
dc2_ = transpose(dc_(iteration_avg_))/sum(dc_(iteration_avg_));
iteration_top_ = iteration_alo:iteration_ahi;
tr_avg_s0000 = dr2_*ZR_s0000_(iteration_avg_);
tr_top_s0000 = max(ZR_s0000_(iteration_top_));
tr_top_s0000_p = (length(find(tr_top_>tr_top_s0000)) + 0.5*length(find(tr_top_==tr_top_s0000)))/n_shuffle;
tr_avg_s0000_p = (length(find(tr_avg_>tr_avg_s0000)) + 0.5*length(find(tr_avg_==tr_avg_s0000)))/n_shuffle;
tmp_tau = min([tr_top_s0000_p,tr_avg_s0000_p]); tmp_tau = max(0.5/n_shuffle,tmp_tau);
ls_rm = find(tr_top_(:)>=prctile(tr_top_(:),100*(1-tmp_tau)));
ls_ra = find(tr_avg_(:)>=prctile(tr_avg_(:),100*(1-tmp_tau)));
p_s0000 = length(unionall({ls_rm,ls_ra}))/n_shuffle;