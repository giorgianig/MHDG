function sol = initializeSolutionToAnalyticalSol_3D(X,T,refEl)

global theta ntor refElTor neq

N2D = size(T,1);
Np2d    = size(refEl.NodesCoord,1);
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Np      = Np1dTor*Np2d;
Ne      = N2D*ntor;
tdiv    = linspace(0,theta,ntor+1);


sol = zeros(Ne*Np,neq);
% generate solution in conservative variables
for itor = 1:ntor
    te = tdiv(itor)+0.5*(refElTor.NodesCoord1d+1)*(tdiv(itor+1)-tdiv(itor));    
    for iel = 1:N2D
        
        iElem = (itor-1)*N2D+iel;
        ind = (iElem-1)*Np + (1:Np);
        Te = T(iel,:);
        Xe = X(Te,:);
        sole = analyticalSolution3d(Xe,te);
        sol(ind,:) = sole;
    end
end
sol = reshape(sol', neq*Ne*Np,1);
