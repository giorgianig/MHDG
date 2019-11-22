function MatIJK = convertMatrixCSRtoIJK(MatCSR)

% Conversion
I = zeros(MatCSR.nnz,1);
J = I;
K = I;

for r = 1:numel(MatCSR.rowptr)-1
    
    for i = MatCSR.rowptr(r):(MatCSR.rowptr(r+1)-1)
    
    I(i) = MatCSR.loc2glob(r);
    J(i) = MatCSR.cols(i);
    K(i) = MatCSR.vals(i);
    end  
end

MatIJK = sparse(I,J,K); 