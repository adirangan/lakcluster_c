function l = slash_diag_0(N);
nx_ = 0:N-1;
ny_ = 0:N-1;
x_ = [nx_+0 , nx_+0 ; nx_+1 , nx_+1];
y_ = [ny_+1 , ny_+0 ; ny_+0 , ny_+1];
l=line(0.5+x_,0.5+y_,1e-3*ones(2,2*N));
