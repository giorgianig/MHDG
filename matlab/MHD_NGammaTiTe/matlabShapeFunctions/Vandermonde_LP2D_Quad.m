function V = Vandermonde_LP2D_Quad(nDeg,coord1d)

N = size(coord1d,1);%number of nodes/polynomials
if ( N^2 ~= (nDeg+1)^2 )
    error('The number of polynomials does not coincide with the number of nodes')
end

V  = zeros(N^2,N^2);

for j = 1:N
    y = coord1d(j);
    py = orthopoly1D(y,nDeg);
    for i = 1:N
        ind = (j-1)*N+i;
        x = coord1d(i);
        px = orthopoly1D(x,nDeg);
        V(ind,:) = col(px*py');
    end
end

perm = getPermutationsQuads(nDeg);

V(perm,:) = V;