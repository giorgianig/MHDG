aux = zeros(size(infoFaces.exteriorFaces_IN,1),1);

for i=1:size(infoFaces.exteriorFaces_IN,1)
    aux(i) = F(infoFaces.exteriorFaces_IN(i,1),infoFaces.exteriorFaces_IN(i,2));
end

disp('Boundary IN')
aux


aux = zeros(size(infoFaces.exteriorFaces_OUT,1),1);

for i=1:size(infoFaces.exteriorFaces_OUT,1)
    aux(i) = F(infoFaces.exteriorFaces_OUT(i,1),infoFaces.exteriorFaces_OUT(i,2));
end

disp('Boundary OUT')
aux