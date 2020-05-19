function u_proj = projectFieldDifferentP(u,p1,p2,Ne)

Nv1 = 0.5*(p1+1)*(p1+2);
Nv2 = 0.5*(p2+1)*(p2+2);
refEl_1 = createReferenceElement(1,Nv1);
refEl_2 = createReferenceElement(1,Nv2);

% Vandermonde matrix
V = Vandermonde_LP(p1,refEl_1.NodesCoord);
invV = inv(V');

% Compute shape functions at interpolation points
shapeFunctions = zeros(Nv2,Nv1);
for ipoint = 1:Nv2
    pp = orthopoly2D(refEl_2.NodesCoord(ipoint,:),p1);
    shapeFunctions(ipoint,:) = (invV*pp)';
end

u = reshape(transpose(reshape(u,2,Ne*Nv1)),Nv1,2*Ne);
u_proj = reshape(transpose(reshape(shapeFunctions*u,Nv2*Ne,2)),2*Ne*Nv2,1);
