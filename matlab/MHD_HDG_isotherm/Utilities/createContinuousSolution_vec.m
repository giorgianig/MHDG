function sol = createContinuousSolution_vec(X,T,u)

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
sol = zeros(2*nOfNodes,1);

for iNode = 1:nOfNodes
    elements = N(iNode,(N(iNode,:)~=0));
    n_elements = length(elements);
    aux_sol_x = zeros(n_elements,1);
    aux_sol_y = zeros(n_elements,1);
    for iElem = 1:n_elements
        elem = elements(iElem);
        Te = T(elem,:);
        local_node = find(Te==iNode,1);
        ind_x = 2*nOfElementNodes*(elem-1) + 2*local_node-1;
        ind_y = 2*nOfElementNodes*(elem-1) + 2*local_node;
        aux_sol_x(iElem) = u(ind_x);
        aux_sol_y(iElem) = u(ind_y);
    end
    ind_node_x = iNode*2-1;
    ind_node_y = iNode*2;
    sol(ind_node_x) = mean(aux_sol_x);
    sol(ind_node_y) = mean(aux_sol_y);
end



