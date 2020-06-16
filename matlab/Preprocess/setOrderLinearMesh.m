function [X, T, Tb_UP, Tb_DOWN, Tb_LEFT, Tb_RIGHT, ...
            elemInfo, elementFaceInfo] = setOrderLinearMesh...
    (X, T, Tb_UP, Tb_DOWN, Tb_LEFT, Tb_RIGHT,...
    elemInfo, elementFaceInfo, p_ref,elemType)


if elemType==1 
    nOfNodesRef = 0.5*(p_ref+1)*(p_ref+2);
    nv = 3;
elseif elemType==0
    nOfNodesRef = (p_ref+1)^2;
    nv = 4;
end
referenceElement = createReferenceElement(elemType,nOfNodesRef);
faceNodes = referenceElement.faceNodes(:,2:end-1); %without vertices
coordRef = referenceElement.NodesCoord;
nOfElementNodes = size(coordRef,1);
nOfFaceNodes = size(faceNodes,2);

[intFaces,elemIntfaceInfo] = GetFaces_mod(T,referenceElement.faceNodes);
nOfInteriorFaces = size(intFaces,1);

nOfElements = size(T,1);
nOfNodes = size(X,1);
Xp = zeros(nOfElementNodes*nOfElements,2);
Tp = zeros(nOfElements,nOfElementNodes);
Xp(1:nOfNodes,:) = X;
Tp(:,1:nv) = T;

%Mesh
elemPos = [1 3];
facePos = [2 4];
auxiliarCoordLogical = true(nOfElementNodes,1);
auxiliarCoordLogical(1:nv) = false;
intFacesMeshed = false(nOfInteriorFaces,1);
ini = nOfNodes + 1;
for elem = 1:nOfElements
    
    coordLogical = auxiliarCoordLogical;
    nOfAlreadyFacesMeshed = 0;
    
    faceInfo = elemIntfaceInfo(elem,:);
    intElemFaces = faceInfo(logical(faceInfo));
    areAlreadyMeshed = intFacesMeshed(intElemFaces,1);
    if any(areAlreadyMeshed)
        facesAlreadyMeshed = intElemFaces(areAlreadyMeshed);
        nOfAlreadyFacesMeshed = length(facesAlreadyMeshed);
        for iface = 1:nOfAlreadyFacesMeshed
            faceMeshed = facesAlreadyMeshed(iface);
            elements = intFaces(faceMeshed,elemPos);
            elementCondition = elements ~= elem;
            elementAlreadyMeshed = elements(elementCondition);
            elementFaceAlreadyMeshed = intFaces(faceMeshed,facePos(elementCondition));
            elementFace = intFaces(faceMeshed,facePos(~elementCondition));
            nodesFaceAlreadyMeshed = faceNodes(elementFaceAlreadyMeshed,:);
            nodesFace = faceNodes(elementFace,:);
            
            Tp(elem,nodesFace) = fliplr(Tp(elementAlreadyMeshed,nodesFaceAlreadyMeshed));
            
            coordLogical(nodesFace) = false;
        end
    end
    
    Te = T(elem,:);
    coordRefMod = coordRef(coordLogical,:);
    if elemType==1
        nodesElemMod = linearMapping(X(Te,:),coordRefMod);
    elseif elemType==0
        nodesElemMod = linearMappingQua(X(Te,:),coordRefMod);
    end
    
    
    ind = ini:ini-1 ...
                 + nOfElementNodes ...
                 - nOfAlreadyFacesMeshed*nOfFaceNodes ...
                 - nv; %without vertices
    
    Xp(ind,:) = nodesElemMod;
    Tp(elem,coordLogical) = ind;
    
    intFacesMeshed(intElemFaces) = true;
    if ~isempty(ind)
        ini = ind(end) + 1;
    end

end

Xp(ini:end,:) = [];

%Boundary
faceNodes = referenceElement.faceNodes;
nOfFaceNodes = size(faceNodes,2);
nameBoundaries = fieldnames(elementFaceInfo);
nOfBoundaries = length(nameBoundaries);
for iboundary = 1:nOfBoundaries
    iname = nameBoundaries{iboundary};
    infoBoundary = elementFaceInfo.(iname);
    nOfBoundaryElements = size(infoBoundary,1);
    Tb = zeros(nOfBoundaryElements,nOfFaceNodes);
    for boundaryElem = 1:nOfBoundaryElements
        iface = infoBoundary(boundaryElem,2);
        ielem = infoBoundary(boundaryElem,1);
        ifaceNodes = faceNodes(iface,:);
        Tb(boundaryElem,:) = Tp(ielem,ifaceNodes);
    end
    
    evalc(['Tb_' iname '= Tb']);
end

%Save
X = Xp;
T = Tp;
elemInfo.nOfNodes = nOfElementNodes;
elemInfo.faceNodes1d = referenceElement.faceNodes1d;
elemInfo.faceNodes = faceNodes;
    
    
    
    
    
    
    
    
    
    
    
    



        