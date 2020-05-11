function plotErrorMap_new(X,T,refEl,error)

% Mesh info
elements = 1:size(T,1);
Ne = numel(elements);

% Creating solution mesh conectivity (tri)
coordRef = refEl.NodesCoord;
elemTriRef = delaunayn(coordRef);
nOfElemTriRef = size(elemTriRef,1);
tri = zeros(nOfElemTriRef*Ne,3);
error_ext = zeros(nOfElemTriRef*Ne,1);
indexElem = 0;
for ielem = 1:Ne
    Te = T(ielem,:);
    for ielemRef = 1:nOfElemTriRef
        indexElemRef = indexElem + ielemRef;
        error_ext(indexElemRef) = error(ielem);
        tri(indexElemRef,:) = Te(elemTriRef(ielemRef,:));
    end
    indexElem = indexElem + nOfElemTriRef;
end
patch('Faces',tri,'Vertices',X,'FaceVertexCData',error_ext,...
    'FaceColor','flat','EdgeAlpha',0);
axis equal


