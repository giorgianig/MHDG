function [Fmod F_per aux_flipFace] = setPeriodicBoundary(X,T,F,F_per,infoFaces,flag1,flag2)

Fmod = F;
aux_flipFace = false(size(F));

% check direction
elem1 = infoFaces.(['exteriorFaces_' flag1])(1,1);
face1 = infoFaces.(['exteriorFaces_' flag1])(1,2);
coord1 = X(T(elem1,findFaceVertices(face1)),:);
norm2bound = find(abs(coord1(1,:) - coord1(2,:)) < 1e-10); 
% 1 - vertical boundary; 2 - horizontal boundary

bound1 = infoFaces.(['exteriorFaces_' flag1]);
bound2 = infoFaces.(['exteriorFaces_' flag2]);

% number of faces in one periodic boundary
nf = size(bound1,1);
for iface1 = 1:nf
    
    el1 = bound1(iface1,1); 
    face1 = bound1(iface1,2);
    F_per(el1,face1) = true;
    coord1 = sort(X(T(el1,findFaceVertices(face1)),3-norm2bound));
    check = true;
    for iface2 = 1:nf 
        
            el2 = bound2(iface2,1); 
            face2 = bound2(iface2,2);
            coord2 = sort(X(T(el2,findFaceVertices(face2)),3-norm2bound));
            if norm(coord1-coord2) < 1e-10
                Fmod(el2,face2) = F(el1,face1);
                F_per(el2,face2) = true;
                if el2>el1
                    aux_flipFace(el1,face1) = true;
                end
                check = false;
                break
            end
    end
    if check
        error('Face not found!')
    end
end



function nodes = findFaceVertices(iface)

switch iface
    case 1
        nodes = [1 2]; 
    case 2
        nodes = [2 3]; 
    case 3
        nodes = [3 1]; 
end