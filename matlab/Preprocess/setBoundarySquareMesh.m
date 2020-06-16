function [Tb_UP, Tb_DOWN, Tb_LEFT, Tb_RIGHT, elementFaceInfo] =...
    setBoundarySquareMesh(X,T,faceNodes)

[~, extFaces] = GetFaces(T,faceNodes);

% initialize
nf = size(extFaces,1); 
Tb_UP = zeros(nf,2);
UP = zeros(nf,2);
Tb_DOWN = zeros(nf,2);
DOWN = zeros(nf,2);
Tb_LEFT = zeros(nf,2);
LEFT = zeros(nf,2);
Tb_RIGHT = zeros(nf,2);
RIGHT = zeros(nf,2);

% initilization
ind_UP = 0;
ind_DOWN = 0;
ind_LEFT = 0;
ind_RIGHT = 0;

tol = 1e-9;
for iFace = 1:nf
    iElem = extFaces(iFace,1);
    iface = extFaces(iFace,2);
    nodes = T(iElem,faceNodes(iface,:));
    coord = X(nodes,:);
    xm = mean(coord(:,1));
    ym = mean(coord(:,2));
    
    if ym==1
        % UP
        ind_UP = ind_UP+1;
        Tb_UP(ind_UP,:) = nodes;
        UP(ind_UP,:) = extFaces(iFace,:);
    elseif ym==0
        % DOWN
        ind_DOWN = ind_DOWN+1;
        Tb_DOWN(ind_DOWN,:) = nodes;
        DOWN(ind_DOWN,:) = extFaces(iFace,:);
    elseif xm==0
        % LEFT
        ind_LEFT = ind_LEFT+1;
        Tb_LEFT(ind_LEFT,:) = nodes;
        LEFT(ind_LEFT,:) = extFaces(iFace,:);
    elseif xm==1
        % RIGHT
        ind_RIGHT = ind_RIGHT+1;
        Tb_RIGHT(ind_RIGHT,:) = nodes;
        RIGHT(ind_RIGHT,:) = extFaces(iFace,:);
    else
        error('Error')
    end 
end
Tb_UP(ind_UP+1:end,:) = [];
Tb_DOWN(ind_DOWN+1:end,:) = [];
Tb_LEFT(ind_LEFT+1:end,:) = [];
Tb_RIGHT(ind_RIGHT+1:end,:) = [];
UP(ind_UP+1:end,:) = [];
DOWN(ind_DOWN+1:end,:) = [];
LEFT(ind_LEFT+1:end,:) = [];
RIGHT(ind_RIGHT+1:end,:) = [];

elementFaceInfo.UP = UP;
elementFaceInfo.DOWN = DOWN;
elementFaceInfo.LEFT = LEFT;
elementFaceInfo.RIGHT = RIGHT;

