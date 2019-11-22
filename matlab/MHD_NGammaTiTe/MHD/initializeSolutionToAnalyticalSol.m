function [sol,solx,soly] = initializeSolutionToAnalyticalSol(X,T)

global neq
sol = zeros(numel(T),neq);
solx = sol;
soly = sol;
Ne = size(T,2);

% generate solution in physical variables
for iel = 1:size(T,1)
    
    ind = (iel-1)*Ne + (1:Ne);
    Te = T(iel,:);
    Xe = X(Te,:);
    [sole,solex,soley] = analyticalSolution(Xe);
    sol(ind,:) = sole;   
    solx(ind,:) = solex;
    soly(ind,:) = soley;
end

sol = reshape(sol', neq*numel(T),1);
solx = reshape(solx', neq*numel(T),1);
soly = reshape(soly', neq*numel(T),1);