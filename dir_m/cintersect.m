function output = cintersect(list_A,list_B);
% Given two permutations of 1:N, ;
% we return the integer-valued cumulative-intersection. ;
% That is to say, output(k) = length(intersect(list_A(1:k),list_B(1:k)));
%
% test with: ;
%{
  % use either: ;
  N = 6; list_A = [5 3 6 4 2 1]; list_B = [6 3 5 1 2 4];
  % or use: ;
  N = 128; list_A = randperm(N); list_B = randperm(N);
  [tmp1,I] = sort(list_A); [tmp2,R] = sort(I); T = I(list_B);
  U = reshape(1:N,size(T));
  S = sparse(max(T,U),U,1);
  C = full(cumsum(sum(S,2)));
  D = zeros(N,1);
  for n=1:N;
  D(n,1) = length(intersect(list_A(1:n),list_B(1:n)));
  end;%for n=1:N;
  disp(sprintf(' %% error: %f',norm(C-D)));
  %}
% Note that we can preprocess the input lists so that they no longer need to be permutations of 1:N ;
%{
  N = 6; list_A = [15 3 16 4 12 1]; list_B = [16 3 15 1 12 4];
  [tmp1,I] = sort(list_A); [tmp2,R] = sort(I); 
  [tmp3,H] = sort(list_B); [tmp4,Q] = sort(H);
  T = I(Q);
  U = reshape(1:N,size(T));
  S = sparse(max(T,U),U,1);
  C = full(cumsum(sum(S,2)));
  D = zeros(N,1);
  for n=1:N;
  D(n,1) = length(intersect(list_A(1:n),list_B(1:n)));
  end;%for n=1:N;
  disp(sprintf(' %% error: %f',norm(C-D)));
  %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(list_A);
[tmp1,I] = sort(list_A); [tmp2,R] = sort(I);
[tmp1,I] = sort(list_A); [tmp2,R] = sort(I); 
[tmp3,H] = sort(list_B); [tmp4,Q] = sort(H);
T = I(Q);
U = reshape(1:N,size(T));
S = sparse(max(T,U),U,1);
C = cumsum(sum(S,2));
output = full(C);
  
  
