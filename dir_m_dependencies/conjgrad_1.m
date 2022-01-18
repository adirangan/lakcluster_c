function conjgrad_1(A, b, x);

if (nargin<1);
N = 16;
A = randn(N); A = transpose(A)*A;
b = randn(N,1); x = zeros(N,1);
conjgrad_1(A,b,x);
disp('returning'); return;
end;%if (nargin<1);

verbose=1;
eps_tolerance = 1e-9;
N = length(b);
assert(size(A,1)==N);
assert(size(A,2)==N);
assert(length(x)==N);
b = reshape(b,N,1);
x = reshape(x,N,1);

R_ = zeros(N,N);
P_ = zeros(N,N);
X_ = zeros(N,N);
alpha_ = zeros(N,1);

X_(:,1) = x;
R_(:,1) = b - A*X_(:,1);
p = R_(:,1);
p = p/sqrt(transpose(p)*A*p);
P_(:,1) = p;

for ni = 1:2*N-1;
% solve along P_(:,ni) ;
% A*(X_(:,ni) + alpha*P_(:,ni)) = b ; 
% If we want to minimize: \| A*(X_(:,ni) + alpha*P_(:,ni)) - b \|^2, ;
% Then alpha = (transpose(P_(:,ni))*transpose(A)*R_(:,ni)) / norm(A*P_(:,ni)).^2 ;
% But If we want to minimze: 1/2 (X_(:,ni) + alpha*P_(:,ni)) * A * transpose(X_(:,ni) + alpha*P_(:,ni)) - transpose(X_(:,ni) + alpha*P_(:,ni))*b ; 
% Then alpha = transpose(P_(:,ni))*R_(:,ni) / transpose(P_(:,ni))*A*P_(:,ni) ;
alpha_(ni) = ( transpose(P_(:,ni))*R_(:,ni) ) / ( transpose(P_(:,ni))*A*P_(:,ni) ) ;
X_(:,ni+1) = X_(:,ni) + alpha_(ni)*P_(:,ni) ;
R_(:,ni+1) = b - A * X_(:,ni+1) ;
p = R_(:,ni+1);
for nl = 1:ni;
p = p - ( transpose(P_(:,nl))*A*p ) * P_(:,nl); 
end;%for nl = 1:ni;
p = p/sqrt(transpose(p)*A*p);
P_(:,ni+1) = p;
eps_error = norm(A*X_(:,ni+1) - b);
if (verbose);
disp(sprintf(' %% iteration %d/%d; eps_error %0.16f',ni,N,eps_error));
end;%if (verbose);
end; %for ni = 1:N-1;

