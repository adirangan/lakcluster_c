function D_ = svd_impute_fit_helper_0(f_scale_,C_,tmp_C_top_,tmp_C_bot_,ij_missed_);
D_ = C_;
D_(ij_missed_) = f_scale_(1)*tmp_C_top_(ij_missed_) + f_scale_(2)*tmp_C_bot_(ij_missed_);
