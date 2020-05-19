function [u_star u_int] = hdg_PostprocessSolution(L,u,K,Bt,int_N,refElv_star,refElv,Ne)

Nv = size(refElv.NodesCoord,1);
coordRef_star = refElv_star.NodesCoord;
Nv_star = size(coordRef_star,1);

% reshape
u = reshape(transpose(reshape(u,2,Ne*Nv)),Nv,2*Ne);  % Nv x 2·Ne
L = reshape(transpose(reshape(L,4,Ne*Nv)),Nv,4*Ne);  % Nv x 4·Ne

% Vandermonde matrix
coordRef = refElv.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = refElv.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');

% Compute shape functions at interpolation points
shapeFunctions = zeros(Nv_star,nOfNodes);
for ipoint = 1:Nv_star
    p = orthopoly2D(coordRef_star(ipoint,:),nDeg);
    shapeFunctions(ipoint,:) = (invV*p)';
end

% interpolate solution at finer mesh
u_int = reshape(transpose(reshape(shapeFunctions*u,Nv_star*Ne,2)),2*Ne*Nv_star,1);
L_int = reshape(transpose(reshape(shapeFunctions*L,Nv_star*Ne,4)),4*Ne*Nv_star,1);

% u star initialization
u_star = zeros(2*Nv_star*Ne,1);

for ielem = 1:Ne

    % index
    ind_u_star = (ielem-1)*2*Nv_star + (1:2*Nv_star);
    ind_L_star = (ielem-1)*4*Nv_star + (1:4*Nv_star);
    
    % elemental matrices
    Ke = K(:,:,ielem);
    Bte = Bt(:,:,ielem);
    int_Ne = int_N(:,:,ielem);
   
    % multiplication
    BtL = Bte*L_int(ind_L_star);
    int_ue_star = int_Ne;
    int_ue = int_Ne'*u_int(ind_u_star);
    
    % Lagrange multipliers
    A = [Ke int_ue_star; int_ue_star' zeros(2)];
    f = [BtL;int_ue];
    
    % elemental solution
    sol = A\f;
    
    % postprocessed solution
    u_star(ind_u_star) = sol(1:end-2);
    
end