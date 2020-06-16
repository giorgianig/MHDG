function shock_st = shockCaptPrep(refEl)

nodes = refEl.NodesCoord;
npoints = size(nodes,1);
refEl_lin = createReferenceElement(1,3);

% Vandermonde matrix
coordRef = refEl_lin.NodesCoord;
nOfNodes = size(coordRef,1);
V = Vandermonde_LP(1,coordRef);
invV = inv(V');

% Compute shape functions at interpolation points
shapeFunctions = zeros(npoints,nOfNodes);
for ipoint = 1:npoints
    p = orthopoly2D(nodes(ipoint,:),1);
    shapeFunctions(ipoint,:) = (invV*p)';
end

shock_st.N = shapeFunctions;
