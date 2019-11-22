function [Xl,Tl,nodalFieldsl]=HDGtoFineLinearApproximation(X,T,nodalFields,nOfSegmentsInEachDirection,referenceElement)
%[Xl,Tl,nodalFieldsl]=HDGtoFineLinearApproximation(X,T,nodalFields,nOfSegmentsInEachDirection,referenceElement)

tol = 1.e-10;

[nOfElements,nOfElementNodes] = size(T);
nodesCoord = referenceElement.NodesCoord;

p = nOfSegmentsInEachDirection;
alpha = linspace(-1,1,p+1);

%Triangularization each element
aux = linspace(-1,1,nOfSegmentsInEachDirection+1);
[x,y]=meshgrid(aux,aux);
x=reshape(x,size(x,1)^2,1); y=reshape(y,size(y,1)^2,1);
i=find(x+y<tol);
Xlref = [x(i),y(i)];
shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,Xlref);
Nstandard = shapeFunctions(:,:,1)';
Tstandard = delaunay(Xlref(:,1),Xlref(:,2));

%Interpolation of coordinates and nodal fields in each element
Xl = []; nodalFieldsl = []; LSl=[]; Tl = [];
for e=1:nOfElements
    Xe = X(T(e,:),:);
    ind = [1:nOfElementNodes]+nOfElementNodes*(e-1);
    Tl = [Tl; Tstandard+size(Xl,1) ];
    Xl = [Xl;Nstandard*Xe]; nodalFieldsl = [nodalFieldsl;Nstandard*nodalFields(ind,:)];
end






