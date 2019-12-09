function [X,T] = sawMeshes(X1,T1,X2,T2,sawCoord,dim)
% X1,T1: mesh 1
% X2,T2: mesh 2
% sawCoord: coordinate at which we make the joining
% dim: along which direction 1-x, 2-y


% number of nodes of the two meshes
nNodes1 = size(X1,1);
nNodes2 = size(X2,1);

% nodes in the meshes
allNodes1 = 1:nNodes1;
allNodes2 = 1:nNodes2;

% Check the break line
tol = 1e-9;
breakNodes1 = allNodes1(abs(X1(:,dim) - sawCoord) < tol);
breakNodes2 = allNodes2(abs(X2(:,dim) - sawCoord) < tol);



if length(breakNodes1) ~= length(breakNodes2)
    error('Number of break nodes in mesh 1 and mesh2 does not coincide')
end
[x1,pos1] = sort(X1(breakNodes1,1));
y1 = X1(breakNodes1(pos1),2);
[x2,pos2] = sort(X2(breakNodes2,1));
y2 = X2(breakNodes2(pos2),2);
xcondition = (x1 < x2 + tol) & (x1 > x2 - tol);
ycondition = (y1 < y2 + tol) & (y1 > y2 - tol);
breakCondition = xcondition & ycondition;
if ~all(breakCondition)
    disp(['The ' breakLine ' nodes...'])
    disp(pos1(~breakCondition))
    disp('do not coincide!!')
    error('Collapse not allowed')
end

% Break nodes position in connectivity matrix
breakNodesValence = 5;
sizeT2 = size(T2);
N2 = getNodalConnectivity(T2);
breakNodes1 = breakNodes1(pos1);
breakNodes2 = breakNodes2(pos2);
nOfBreakNodes = length(breakNodes1);
breakNum = ones(nOfBreakNodes,1);
breakNodes2Pos = zeros(nOfBreakNodes,breakNodesValence);
for inode = 1:nOfBreakNodes
    ibreak2 = breakNodes2(inode);
    elem2 = N2(ibreak2,logical(N2(ibreak2,:)));
    for ielem = 1:length(elem2)
        breakElem = elem2(ielem);
        breakPos = find(T2(breakElem,:) == ibreak2);
        T2(breakElem,breakPos) = nNodes2 + 1; %reset break nodes to maximum value
        breakNodes2Pos(inode,breakNum(inode)) = sub2ind(sizeT2,breakElem,breakPos);
        breakNum(inode) = breakNum(inode) + 1;
    end
end
breakNodes2Pos(:,max(breakNum):end) = [];

% New nodal position
X2(breakNodes2,:) = [];

% New connectivity
T2 = T2 + nNodes1; %with double break nodes
breakNodes2aux = sort(breakNodes2) + nNodes1;
for inode = 1:nOfBreakNodes
    ibreak2aux = breakNodes2aux(inode);
    nodesCondition = T2 > ibreak2aux;
    T2(nodesCondition) = T2(nodesCondition) - 1;
    breakNodes2aux(inode+1:nOfBreakNodes) = breakNodes2aux(inode+1:nOfBreakNodes) - 1;
end
minNum = min(min(T2));
T2 = T2 - minNum + nNodes1 + 1; %without double break nodes

% Add break nodes of mesh1
for inode = 1:nOfBreakNodes
    ibreak1 = breakNodes1(inode);
    ibreak2Pos = breakNodes2Pos(inode,logical(breakNodes2Pos(inode,:)));
    T2(ibreak2Pos) = ibreak1;
end

% New mesh collapsed
T = [T1 ; T2];
X = [X1 ; X2];

%---------------------------------------------------------
    
function N = getNodalConnectivity(T)

nNodes = max(max(T));
nNodesElem = size(T,2);
N = zeros(nNodes,10);
nn = ones(nNodes,1);
for ielem = 1:size(T,1)
    Te = T(ielem,:);
    nn_Te = nn(Te);
    for kk = 1:nNodesElem
        N(Te(kk),nn_Te(kk)) = ielem;
    end
    nn(Te) = nn(Te) + 1;
end
N(:,max(nn):end) = [];

    
    
    
    