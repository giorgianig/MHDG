function shapeFunctions=computeShapeFunctionsAtPointsQuads(nDeg,coordRef,points)
%    
    coord1d = col(points');
    
    %Vandermonde matrix
    nOfNodes1d = size(coordRef,1);
    V = Vandermonde_LP(nDeg,coordRef);
    [L,U,P] = lu(V');

    % Compute shape functions at interpolation points
    sf1d = zeros(numel(coord1d),nOfNodes1d,2);
    for ipoint = 1:numel(coord1d)
        [p,p_xi] = orthopoly1D_deriv(coord1d(ipoint),nDeg);
        N = U\(L\(P*[p,p_xi]));
        sf1d(ipoint,:,1) = N(:,1);
        sf1d(ipoint,:,2) = N(:,2);
    end
    % permutations
    shapeFunctions = zeros(size(points,1),size(coordRef,1)^2,3);
    perm = getPermutationsQuads(nDeg);
    for i=1:size(points,1)
        ind = (i-1)*2+1;
        aux = col(sf1d(ind,:,1)'*sf1d(ind+1,:,1));
        auxx = col(sf1d(ind,:,2)'*sf1d(ind+1,:,1));
         auxy = col(sf1d(ind,:,1)'*sf1d(ind+1,:,2));
        shapeFunctions(i,perm,1) = aux';
        shapeFunctions(i,perm,2) = auxx';
        shapeFunctions(i,perm,3) = auxy';
    end
    
%     shapeFunctions = createShapeFunctions2dTensor(sf1d,zeros(size(coord1d)),coord1d,perm);
%     shapeFunctions = shapeFunctions(:,:,1)';
