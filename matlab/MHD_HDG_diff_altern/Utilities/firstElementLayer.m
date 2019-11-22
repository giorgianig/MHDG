function elements = firstElementLayer(F,exteriorFaces,interiorFaces)

elements = exteriorFaces(:,1);
Faces = unique(F(elements,:));
nIntFaces = size(interiorFaces,1);
Faces = sort(Faces(Faces<nIntFaces));
elements = [elements;interiorFaces(Faces,1);interiorFaces(Faces,3)];
elements = unique(elements);
