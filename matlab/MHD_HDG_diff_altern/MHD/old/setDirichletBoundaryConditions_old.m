function [BC init] = setDirichletBoundaryConditions(X,T,boundaryNames,F,infoFaces,refElv)

nv = size(refElv.NodesCoord1d,1);
faceNodes = refElv.faceNodes;
Nv = size(refElv.NodesCoord,1);

% initialization
init = struct('ind',[],'val',[]);
BC = struct('ind',[],'val',[]);
P_cond = true; % condition on the pressure

% loop in boundaries
for iboundary = 1:numel(boundaryNames)
    iname = boundaryNames{iboundary}(4:end);
    infoFace = infoFaces.(['exteriorFaces_' iname]);
    ind = []; val = []; ind4init = [];
    if strcmp(iname,'OUTLET')
        P_cond = setConditionOutlet;
    else
        if strcmp(iname,'INLET')
            [ind,val,ind4init] = setConditionVx(F,infoFace,nv,Nv,faceNodes);
        elseif strcmp(iname,'DOWN')
%             [ind,val,ind4init] = setConditionVx(F,infoFace,nv,Nv,faceNodes);
P_cond = setConditionOutlet;
        elseif strcmp(iname,'LEFT')
            [ind,val,ind4init] = setConditionVx(F,infoFace,nv,Nv,faceNodes);
        elseif strcmp(iname,'RIGHT')
%             [ind,val,ind4init] = setConditionWallLeaky(X,T,F,infoFace,nv,Nv,faceNodes);
        elseif strcmp(iname,'UP')
%             [ind,val,ind4init] = setConditionVx(F,infoFace,nv,Nv,faceNodes);
        elseif strcmp(iname,'WALL')
            [ind,val,ind4init] = setConditionWall(F,infoFace,nv,Nv,faceNodes);
        end
        BC.ind = [BC.ind; ind];
        BC.val = [BC.val; val];
        BC.P_cond = P_cond;
        init.ind = [init.ind; ind4init];
        init.val = [init.val; val];
    end

end

function [ind,val,ind4init] = setConditionVx(F,infoFace,nv,Nv,faceNodes)

% number of faces in the subset
nf = size(infoFace,1);

% initialize
ind = zeros(nf*2*nv,1);
val = ind;
ind4init = ind;

% set the condition
for iface = 1:nf
    iElem = infoFace(iface,1);
    iFace = infoFace(iface,2);
    aux_ind = (iface-1)*2*nv + (1:2*nv);
    Fe = F(iElem,iFace);
    indx = (iface-1)*2*nv + (1:2:2*nv-1);
    ind(aux_ind) = (Fe-1)*2*nv + (1:2*nv);
    val(indx) = 1;

    % index for initialization
    nodes = reshape(bsxfun(@plus,2*faceNodes(iFace,:)-1,(0:1)'),2*nv,1);
    ind4init(aux_ind) = (iElem-1)*2*Nv + nodes;
end

function [ind,val,ind4init] = setConditionVy(F,infoFace,nv,Nv,faceNodes)

% number of faces in the subset
nf = size(infoFace,1);

% initialize
ind = zeros(nf*2*nv,1);
val = ind;
ind4init = ind;

% set the condition
for iface = 1:nf
    iElem = infoFace(iface,1);
    iFace = infoFace(iface,2);
    aux_ind = (iface-1)*2*nv + (1:2*nv);
    Fe = F(iElem,iFace);
    indy = (iface-1)*2*nv + (2:2:2*nv);
    ind(aux_ind) = (Fe-1)*2*nv + (1:2*nv);
    val(indy) = 1;

    % index for initialization
    nodes = reshape(bsxfun(@plus,2*faceNodes(iFace,:)-1,(0:1)'),2*nv,1);
    ind4init(aux_ind) = (iElem-1)*2*Nv + nodes;
end


function [ind,val,ind4init] = setConditionWall(F,infoFace,nv,Nv,faceNodes)

% number of faces in the subset
nf = size(infoFace,1);

% initialize
ind = zeros(nf*2*nv,1);
val = ind;
ind4init = ind;

% set the condition
for iface = 1:nf
    iElem = infoFace(iface,1);
    iFace = infoFace(iface,2);
    aux_ind = (iface-1)*2*nv + (1:2*nv);
    Fe = F(iElem,iFace);
    ind(aux_ind) = (Fe-1)*2*nv + (1:2*nv);

    % index for initialization
    nodes = reshape(bsxfun(@plus,2*faceNodes(iFace,:)-1,(0:1)'),2*nv,1);
    ind4init(aux_ind) = (iElem-1)*2*Nv + nodes;
end

function [ind,val,ind4init] = setConditionCirc(X,T,F,infoFace,nv,Nv,faceNodes)

% number of faces in the subset
nf = size(infoFace,1);

% initialize
ind = zeros(nf*2*nv,1);
val = ind;
ind4init = ind;

% set the condition
for iface = 1:nf
    iElem = infoFace(iface,1);
    iFace = infoFace(iface,2);
    Tf = T(iElem,faceNodes(iFace,:));
    Xf = X(Tf,:);
    xy_s = Xf-0.5;
    aux_ind = (iface-1)*2*nv + (1:2*nv);
    Fe = F(iElem,iFace);
    ind(aux_ind) = (Fe-1)*2*nv + (1:2*nv);
    indx = (iface-1)*2*nv + (1:2:2*nv-1);
    indy = (iface-1)*2*nv + (2:2:2*nv);
    val(indx) = xy_s(:,2);
    val(indy) = -xy_s(:,1);
    % index for initialization
    nodes = reshape(bsxfun(@plus,2*faceNodes(iFace,:)-1,(0:1)'),2*nv,1);
    ind4init(aux_ind) = (iElem-1)*2*Nv + nodes;
end

function [ind,val,ind4init] = setConditionWallLeaky(X,T,F,infoFace,nv,Nv,faceNodes)

% number of faces in the subset
nf = size(infoFace,1);

% initialize
ind = zeros(nf*2*nv,1);
val = ind;
ind4init = ind;

% set the condition
for iface = 1:nf
    iElem = infoFace(iface,1);
    iFace = infoFace(iface,2);
    Tf = T(iElem,faceNodes(iFace,:));
    yf = X(Tf,2);
    aux_ind = (iface-1)*2*nv + (1:2*nv);
    Fe = F(iElem,iFace);
    ind(aux_ind) = (Fe-1)*2*nv + (1:2*nv);
    if any(yf == 1)
        indx = (iface-1)*2*nv + (1:2:2*nv-1);
        yf_min = min(yf);
        val(indx) = (yf-yf_min)/(1-yf_min);
    end
    % index for initialization
    nodes = reshape(bsxfun(@plus,2*faceNodes(iFace,:)-1,(0:1)'),2*nv,1);
    ind4init(aux_ind) = (iElem-1)*2*Nv + nodes;
end

function P_cond = setConditionOutlet

P_cond = false;

