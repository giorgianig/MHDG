function res = expandMatrixA_mult(A,n,mult)
% expand matrix A and M
%  [ A 0 0 0
%    0 A 0 0
%    0 0 A 0
%    0 0 0 A ]
% dimension n
res = zeros([size(A) n n]);
aux = repmat(A, [1 1 n]); 
aux(:,:,1) = aux(:,:,1)*mult(1);
aux(:,:,2) = aux(:,:,2)*mult(2);
aux(:,:,3) = aux(:,:,3)*mult(3);
aux(:,:,4) = aux(:,:,4)*mult(4);
res(:,:,1:n+1:n^2) = aux;
res = permute(res, [3 1 4 2]);
res = reshape(res, n*size(A));