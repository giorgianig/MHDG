function sol = initializeSolutionWithL2Projection(X,T,refEl)

sol = zeros(numel(T),2);
Ne = size(T,2);
Nel = size(T,1);
Np = size(X,1);
M = spalloc(Np,Np,3*Np);
fn = zeros(Np,1);
fu = zeros(Np,1);


% generate mass matrix and rhs
for iel = 1:Nel

    Te = T(iel,:);
    Xe = X(Te,:);
    
    % elemental matrices
    [Me,fne,fue] = elementalMatrix(Xe,refEl);
    
    % Assembly
    M(Te,Te) = M(Te,Te) + Me;
    fn(Te) = fn(Te) + fne;    
    fu(Te) = fu(Te) + fue;    
end

ncont = M\fn;
ucont = M\fu;

% generate solution in physical variables
for iel = 1:Nel
    Te = T(iel,:);
    ind = (iel-1)*Ne + (1:Ne);
    sol(ind,1) = ncont(Te);
    sol(ind,2) = ucont(Te);
end

sol = reshape(sol', 2*numel(T),1);




function [Me,fn,fu] = elementalMatrix(Xe,refEl)

np = size(refEl.NodesCoord,1);
fn = zeros(np,1);
fu = fn;
Me = zeros(np,np);

% Information of the reference element
IPw = refEl.IPweights;
N = refEl.N;
Nxi = refEl.Nxi;
Neta = refEl.Neta;
ngauss = length(IPw);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% VOLUME COMPUTATIONS: LOOP IN GAUSS POINTS
for g = 1:ngauss
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];

    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end
    x = N_g*xe;
    y = N_g*ye;

    %x and y derivatives
    invJ = inv(J); 
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    
    %Integration weight
    dvolu=IPw(g)*det(J);
    sol = analyticalSolution([x y]);
    eps = 1e-5;
    %Contribution of the current integration point to the elemental matrix
    Me = Me + (N_g'*N_g+eps*(Nx_g'*Nx_g + Ny_g'*Ny_g) )*dvolu;
    fn = fn + N_g'*sol(1)*dvolu;
    fu = fu + N_g'*sol(2)*dvolu;
end
