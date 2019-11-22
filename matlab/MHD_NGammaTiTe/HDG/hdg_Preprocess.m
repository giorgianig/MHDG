function [F,F_dir,infoFaces,flipFace] = hdg_preprocess(T,elementFaceInfo,boundaryNames,refEl,Dirichlet_boundaries)


% Creating infoFaces structure
if refEl.elemType==1
    intFaces = GetFaces(T(:,1:3),refEl.faceNodes);
    nf = 3;
elseif refEl.elemType==0
    intFaces = GetFaces(T(:,1:4),refEl.faceNodes);
    nf = 4;
else
    error('Wrong element type')
end
infoFaces.interiorFaces = intFaces;
nOfBoundaries = numel(boundaryNames);
for iboundary = 1:nOfBoundaries
    iname = boundaryNames{iboundary}(4:end);
    infoFaces.(['exteriorFaces_' iname]) = elementFaceInfo.(iname);
end

% Creating faces connectivity
nOfElements = size(T,1);
nOfInteriorFaces = size(intFaces,1);
F = zeros(nOfElements,nf);
F_dir = false(size(F));
flipFace = false(nOfElements,nf);

%% interior faces
for iFace = 1:nOfInteriorFaces
    infoFace = intFaces(iFace,:);
    if F(infoFace(1),infoFace(2))~=0, stop, end
    F(infoFace(1),infoFace(2)) = iFace;
    if F(infoFace(3),infoFace(4))~=0, stop, end
    F(infoFace(3),infoFace(4)) = iFace;
    if infoFace(1)<infoFace(3)
        flipFace(infoFace(3),infoFace(4)) = true;
    else
        flipFace(infoFace(1),infoFace(2)) = true;
    end
end

%% exterior faces
% put Diriclet faces on the bottom
[Dir_sort, renum_bound] = sort(Dirichlet_boundaries);
ind = nOfInteriorFaces;
boundaryNames_sort = boundaryNames(renum_bound);
for iboundary = 1:nOfBoundaries
    iname = boundaryNames_sort{iboundary}(4:end);
    isdir = Dir_sort(iboundary);
    extFaces = infoFaces.(['exteriorFaces_' iname]);
    nOfFaces = size(extFaces,1);
    for iFace = 1:nOfFaces
        infoFace = extFaces(iFace,:);
        if F(infoFace(1),infoFace(2))~=0, stop, end
        F(infoFace(1),infoFace(2)) = iFace + ind;
        F_dir(infoFace(1),infoFace(2)) = isdir;
    end
    ind = ind + nOfFaces;
end

% check 
if ~isempty(setdiff([min(min(F)):max(max(F))],unique(F)))
    error('Something wrong in the face connectivity')
end
