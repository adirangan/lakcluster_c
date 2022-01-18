function lp_ = z_to_p_0(z_);
% converts z-scores to p-values. ;
% one sided test (i.e., p-value of obtaining a larger z-score). ;
% returns log of p-value. ;
lp_ = log(0.5) + erfcln(z_/sqrt(2));
