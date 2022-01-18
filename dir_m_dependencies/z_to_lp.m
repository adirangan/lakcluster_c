function lp_ = z_to_lp(z_);
% converts z-score to log(p-value). ;
% one sided test (i.e., p-value of obtaining a larger z-score). ;
% returns log of p-value. ;
lp_ = log(0.5) + erfcln(z_/sqrt(2));
