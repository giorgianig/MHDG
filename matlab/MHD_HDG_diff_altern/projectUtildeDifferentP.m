function u_proj = projectUtildeDifferentP(u_tilde,p1,p2,Nf)


Nv1 = 0.5*(p1+1)*(p1+2);
Nv2 = 0.5*(p2+1)*(p2+2);
nv1 = p1+1;
nv2 = p2+1;
refEl_1 = createReferenceElement(1,Nv1);
refEl_2 = createReferenceElement(1,Nv2);

% Vandermonde matrix
V = Vandermonde_LP(p1,refEl_1.NodesCoord1d);
invV = inv(V');

% Compute shape functions at interpolation points
shapeFunctions = zeros(nv2,nv1);
for ipoint = 1:nv2
    pp = orthopoly1D(refEl_2.NodesCoord1d(ipoint),p1);
    shapeFunctions(ipoint,:) = (invV*pp)';
end

u_tilde = reshape(transpose(reshape(u_tilde,2,Nf*nv1)),nv1,2*Nf);
u_proj = reshape(transpose(reshape(shapeFunctions*u_tilde,nv2*Nf,2)),2*Nf*nv2,1);



