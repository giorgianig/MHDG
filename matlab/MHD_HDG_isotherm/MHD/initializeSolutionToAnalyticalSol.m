function sol = initializeSolutionToAnalyticalSol(X,T)

sol = zeros(numel(T),2);
Ne = size(T,2);

% generate solution in physical variables
for iel = 1:size(T,1)
    
    ind = (iel-1)*Ne + (1:Ne);
    Te = T(iel,:);
    Xe = X(Te,:);
    sole = analyticalSolution(Xe);
    sol(ind,:) = sole;   
end

sol = reshape(sol', 2*numel(T),1);
