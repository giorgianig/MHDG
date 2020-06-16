function sol = createContinuousSolution_scal(X,T,u)

nOfNodes = size(X,1);
nOfElement = size(T,1);
nOfElementNodes = size(T,2);

% nodal connectivity
N = zeros(nOfNodes,10);
nn = ones(nOfNodes,1);
for ielem = 1:nOfElement
    Te = T(ielem,:);
    nn_Te = nn(Te);
    for kk = 1:nOfElementNodes
        N(Te(kk),nn_Te(kk)) = ielem;
    end
    nn(Te) = nn(Te) +1;
end
N(:,max(nn):end) = [];

% suavized (continuos) solution
sol = zeros(nOfNodes,1);

for iNode = 1:nOfNodes
    elements = N(iNode,(N(iNode,:)~=0));
    n_elements = length(elements);
    aux_sol = zeros(n_elements,1);
    for iElem = 1:n_elements
        elem = elements(iElem);
        Te = T(elem,:);
        local_node = find(Te==iNode,1);
        ind = nOfElementNodes*(elem-1) + local_node;
        aux_sol(iElem) = u(ind);
    end
    sol(iNode) = mean(aux_sol);
end

