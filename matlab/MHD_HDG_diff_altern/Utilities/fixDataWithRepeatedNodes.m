
% tolerance
tol = 1e9;

allNodes = 1:size(X,1);  % all nodes with the repeated ones
X_round = floor(tol*X);
[aux,j] =unique(X_round,'rows');
rep = setdiff(allNodes,j);  % index of the repeated nodes (the first ones in the X matrix)

% take out boundary nodes from the repeated nodes
% boundaryNodes = [];
% for iboundary = 1:numel(boundaryNames)
%     name = boundaryNames{iboundary}; 
%     eval(['Tb = Tb_' name(4:end) ';'])
%     boundaryNodes = [boundaryNodes;unique(Tb)];
% end
% rep = setdiff(rep,unique(boundaryNodes));

X_round_aux = X_round;
X_round_aux(rep,:) = NaN;
[aux,loc]=ismember(X_round(rep,:),X_round_aux,'rows');
rep_pair = [rep',loc];  % repeated pair index [first index; second index]
if isempty(rep_pair)
    disp('Mesh is OK!')
    return
end
% fix X
disp(['Found ' num2str(numel(rep_pair)) ' repeated nodes'])
notRepeatedNodes = setdiff(allNodes,rep_pair(:,2));
X_fix = X(notRepeatedNodes,:);

% fix T
T_fix = T;
for iNode = 1:numel(rep)
    node = rep_pair(iNode,1);
    node_rep = rep_pair(iNode,2);
    T_fix(T_fix==node_rep) = node;
    T_fix(T_fix>node_rep) = T_fix(T_fix>node_rep) - 1;
    rep_pair(rep_pair>node_rep) = rep_pair(rep_pair>node_rep) - 1;
end

% fix boundaries
aux_change = T_fix~=T;
change_i = T(aux_change(:));
change_j = T_fix(aux_change(:));
change = [change_i,change_j];
[aux,ind] = unique(change(:,1));
change = change(ind,:);
for iboundary = 1:numel(boundaryNames)
    name = boundaryNames{iboundary}; 
    eval(['Tb = Tb_' name(4:end) ';'])
    Tb_fix = Tb;
    [ind, ind2change] = ismember(Tb,change(:,1));
    Tb_fix(ind(:)) = change(ind2change(ind),2);
    auxTb = reshape(Tb_fix,size(Tb));
    Tb_fix = auxTb;
    eval([name '= Tb_fix;']);
end
clear Tb_fix

X = X_fix;
T = T_fix;

