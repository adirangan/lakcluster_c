function output = rup(n,stride) ;
n_extend = mod(stride - mod(n,stride),stride);
output = n+n_extend;
