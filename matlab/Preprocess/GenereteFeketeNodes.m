function Xmod = GenereteFeketeNodes(X,T,refEl,nDegRef)

% Generate equal spaced points of the reference elements with same order of
% ez4u mesh

np = size(T,2);
refEl = createReferenceElement(1,np);
nDegRef = refEl.degree;
nodes = [0 0;1 0;0 1];
h = 1/nDegRef;
for i=1:nDegRef-1
    nodes = [nodes; [i 0]*h];
end
for i=1:nDegRef-1
    j = nDegRef - i;
    nodes = [nodes; [j i]*h];
end
for i=1:nDegRef-1
    j = nDegRef - i;
    nodes = [nodes; [0 j]*h];
end
for j = 1:nDegRef-1
    i = (1:nDegRef-j -1)';
    aux = j*ones(size(i));
    nodes = [nodes; [i aux]*h];
end
nodes = 2*nodes - 1;
npoints = size(nodes,1);

% Vandermonde matrix
coordRef = refEl.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = refEl.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');

% Compute shape functions at interpolation points
shapeFunctions = zeros(npoints,nOfNodes);
for ipoint = 1:npoints
    p = orthopoly2D(nodes(ipoint,:),nDeg);
    shapeFunctions(ipoint,:) = (invV*p)';
end

%Generate Fekete nodes from equispacied nodes
nEl = size(T,1);

Xmod = zeros(size(X));
for ielem = 1:nEl
    Te = T(ielem,:);
    Xmod(Te,:) = inv(shapeFunctions)*X(Te,:);
end
end