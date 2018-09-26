function output = make_checkered_array(nrows,ncols,blocksize);
output = ones(nrows,ncols);
rij = find(mod(0:nrows-1,2*blocksize)<blocksize);
output(rij,:) = 1-output(rij,:);
cij = find(mod(0:ncols-1,2*blocksize)<blocksize);
output(:,cij) = 1-output(:,cij);
