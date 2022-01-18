function z_ = a_to_z_to_p_0(x_,a_,s_);
% calculates z-score, then p-value. ;
% x_ = data.; 
% a_ = average.;
% s_ = std. ;
tmp_ij_ = find(s_>0);
z_ = zeros(size(x_));
z_(tmp_ij_) = (x_(tmp_ij_)-a_(tmp_ij_))./s_(tmp_ij_);
z_(tmp_ij_) = erfc( abs(x_(tmp_ij_) - a_(tmp_ij_))./s_(tmp_ij_) / sqrt(2) );
