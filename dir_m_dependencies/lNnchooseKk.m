function output = lNnchooseKk(N,K,n,k);

output = lnchoosek(K,k) + lnchoosek(N-K,n-k) - lnchoosek(N,n);
