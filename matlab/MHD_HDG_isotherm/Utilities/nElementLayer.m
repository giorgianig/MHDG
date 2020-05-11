function elements = nElementLayer(F,exteriorFaces,interiorFaces,n)
if n>1
    n =n*2-1;
end

elements = exteriorFaces(:,1);
for ilev = 1:n

    Faces = unique(F(elements,:));
    nIntFaces = size(interiorFaces,1);
    Faces = sort(Faces(Faces<nIntFaces));
    elements = [elements;interiorFaces(Faces,1);interiorFaces(Faces,3)];
    elements = unique(elements);
end