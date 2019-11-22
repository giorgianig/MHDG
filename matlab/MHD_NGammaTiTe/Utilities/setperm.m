function res = setperm(n,m)

% check
if mod(n,m),error('Error: n must be a multiple of m '), end
res = col(fliplr(reshape(1:n,[m n/m])));