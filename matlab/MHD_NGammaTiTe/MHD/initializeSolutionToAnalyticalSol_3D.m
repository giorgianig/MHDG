function [sol,solx,soly,solt]= initializeSolutionToAnalyticalSol_3D(X,T,refEl)

global theta ntor refElTor neq

N2D = size(T,1);
Np2d    = size(refEl.NodesCoord,1);
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Np      = Np1dTor*Np2d;
Ne      = N2D*ntor;
tdiv    = linspace(0,theta,ntor+1);


sol = zeros(Ne*Np,neq);
solx = sol;
soly = sol;
solt = sol;
% generate solution in conservative variables
for itor = 1:ntor
    te = tdiv(itor)+0.5*(refElTor.NodesCoord1d+1)*(tdiv(itor+1)-tdiv(itor));    
    for iel = 1:N2D
        
        iElem = (itor-1)*N2D+iel;
        ind = (iElem-1)*Np + (1:Np);
        Te = T(iel,:);
        Xe = X(Te,:);
        [sole, solex, soley, solet]= analyticalSolution_3D(Xe,te);
        sol(ind,:) = sole;
        solx(ind,:) = solex;
        soly(ind,:) = soley;
        solt(ind,:) = solet;
    end
end
sol = reshape(sol', neq*Ne*Np,1);
solx = reshape(solx', neq*Ne*Np,1);
soly = reshape(soly', neq*Ne*Np,1);
solt = reshape(solt', neq*Ne*Np,1);
