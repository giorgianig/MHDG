function [X_fix,T_fix,rep_pair] = fixMeshWithRepeatedNodes(X,T)
% fix mesh with repeated nodes
% X(:,1) = X(:,1)-X(1,1);
% X(:,2) = X(:,2)-X(1,2);
tol = 1e8;
allNodes = 1:size(X,1);  % all nodes with the repeated ones
X_round = floor(tol*X);
[B,j] =unique(X_round,'rows');
rep = setdiff(allNodes,j);  % index of the repeated nodes (the first ones in the X matrix)
X_round_aux = X_round;
X_round_aux(rep,:) = NaN;
[tf,loc]=ismember(X_round(rep,:),X_round_aux,'rows');
rep_pair = [rep',loc];  % repeated pair index [first index; second index]

if isempty(rep_pair)
    disp('Mesh is OK!')
    return
end
disp(['Found ' num2str(numel(rep_pair)) ' repeated nodes'])
% fix X
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


% X(:,1) = X(:,1)-X(1,1);
% X(:,2) = X(:,2)-X(1,2);