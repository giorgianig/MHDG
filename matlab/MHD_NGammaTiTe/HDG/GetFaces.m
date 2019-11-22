function [intFaces,extFaces] = GetFaces(T,faceNodes)
%
% [intFaces,extFaces] = GetFaces(T)
% (only for triangles and tetrahedrons)
%
% For every face i:
% intFaces(i,:)=[element1 nface1 element2 nface2 node1] for interior faces
% extFaces(i,:)=[element1 nface1] for exterior faces
%
% element1, element2:   number of the elements
% nface1, nface2:       number of face in each element
% node1:  number of node in the 2nd element that matches with the 1st node
%         in the 1st element (with the numbering of the face)
%
% Input:
% T: connectivity of the mesh
%
% Output:
% intFaces,extFaces: interior and exterior faces
%

[nElem,nen] = size(T);
nfaceel     = size(faceNodes,1);


markE    = zeros(nElem,nfaceel);
intFaces = [];
extFaces = [];

%Definition of the faces in the reference element
switch nen
    case 3 %triangle
        Efaces = [1 2; 2 3; 3 1];
    case 4 % quads and tetrahedra
        if size(faceNodes,1) == 4 % quads
            Efaces = [1 2; 2 3; 3 4; 4 1];
        elseif size(faceNodes,1) == 3 % tetra
            Efaces = computeNodesFacesTetra(1);
        end
    case 8 %hexahedra
        Efaces = computeNodesFacesHexa(1);
end

for iElem=1:nElem
    for iFace=1:nfaceel
        if(markE(iElem,iFace)==0)
            markE(iElem,iFace)=1;
            nodesf = T(iElem,Efaces(iFace,:));
            jelem = FindElem(T,iElem,nodesf);
            if(jelem~=0)
                [jface,node1]=FindFace(nodesf,T(jelem,:),Efaces);
                intFaces=[intFaces; iElem,iFace,jelem,jface,node1];
                markE(jelem,jface)=1;
            else
                extFaces=[extFaces; iElem,iFace];
            end
        end
    end
end

%Auxiliar functions
function jelem = FindElem(T,iElem,nodesf)

nen = length(nodesf);

[elems,aux] = find(T==nodesf(1));
elems=elems(elems~=iElem);
T=T(elems,:);
for i=2:nen
    if(~isempty(elems))
        [aux,aux2] = find(T==nodesf(i));
        elems = elems(aux);
        T=T(aux,:);
    end
end

if(isempty(elems))
    jelem=0;
else
    jelem=elems(1);
end

function [jface,node1]=FindFace(nodesf,nodesE,Efaces)

nFaces = size(Efaces,1);
for j=1:nFaces
    nodesj = nodesE(Efaces(j,:));
    if(isempty(setdiff(nodesj,nodesf)))
        jface = j;
        node1 = find(nodesj==(nodesf(1)));
        break;
    end
end


