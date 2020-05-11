function sol = initializeSolution(X,T,F_dir,refEl)

% mesh
neq = 2;
Ne = size(T,1);
Nv = size(T,2);

% sol
sol = zeros(Ne*Nv,neq);
sol(:,1) = 1;
sol(:,2) = 0.1;
%
elements = (1:Ne);
Dir_elem = elements(any(F_dir,2));

% set BC
for iElem = Dir_elem
    Dir_nodes = unique(refEl.faceNodes(F_dir(iElem,:),:));
    ind = (iElem-1)*Nv + Dir_nodes;
    Xf = X(T(iElem,Dir_nodes),:);
    u_ex = analyticalSolution(Xf);
    sol(ind,:) = u_ex;
end
sol = reshape(sol', neq*Ne*Nv,1);