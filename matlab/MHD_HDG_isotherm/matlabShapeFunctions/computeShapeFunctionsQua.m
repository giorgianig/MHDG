function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsQua(nDeg,coord,nOfGaussPoints)
%
% [shapeFunctions,gaussWeights]=computeShapeFunctionsReferenceElement(nDeg,
% coord)
%
% Function to compute the shape functions (& derivatives) at Gauss points
% for a quadrilateral reference element
%
% Input:
% nDeg:  degree of interpolation
% coord: nodal coordinates at the reference element
% nOfGaussPoints: nº of gauss points of the 1D quadrature
%
% Output:
% shapeFunctions: shape functions evaluated at the gauss points
%                 size is nOfNodes X nOfGauss X (nsd + 1)
%                 nsd+1 because shape function (1)
%                 and derivatives (nsd) are stored
% gaussWeights:   weights
%

[z,w] = gaussLegendre(nOfGaussPoints,-1,1);
nIP = length(w); %number of integration points in each direction
nOfGauss = nIP^2;

nOfNodes = size(coord,1);
shapeFunctions = zeros(nOfNodes,nOfGauss,3);
gaussWeights = zeros(nOfGauss,1);
gaussPoints = zeros(nOfGauss,2);

%position and weights of gauss points
iGauss = 1;
for i = 1:nIP
    for j = 1:nIP
        gaussPoints(iGauss,:) = [z(i),z(j)];
        gaussWeights(iGauss) = w(i)*w(j);
        iGauss = iGauss + 1;
    end
end

%Shape functions and derivatives
xi = gaussPoints(:,1); eta = gaussPoints(:,2);
switch nDeg
    case 1 %Q1
        N    = [(1-xi).*(1-eta)/4, (1+xi).*(1-eta)/4, ...
            (1+xi).*(1+eta)/4, (1-xi).*(1+eta)/4];
        Nxi  = [(eta-1)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4];
        Neta = [(xi-1)/4, -(1+xi)/4,   (1+xi)/4,  (1-xi)/4 ];
    case 2 %Q2
        N    = [xi.*(xi-1).*eta.*(eta-1)/4, xi.*(xi+1).*eta.*(eta-1)/4, ...
            xi.*(xi+1).*eta.*(eta+1)/4, xi.*(xi-1).*eta.*(eta+1)/4, ...
            (1-xi.^2).*eta.*(eta-1)/2,  xi.*(xi+1).*(1-eta.^2)/2,   ...
            (1-xi.^2).*eta.*(eta+1)/2,  xi.*(xi-1).*(1-eta.^2)/2,   ...
            (1-xi.^2).*(1-eta.^2)];
        Nxi  = [(xi-1/2).*eta.*(eta-1)/2,   (xi+1/2).*eta.*(eta-1)/2, ...
            (xi+1/2).*eta.*(eta+1)/2,   (xi-1/2).*eta.*(eta+1)/2, ...
            -xi.*eta.*(eta-1),          (xi+1/2).*(1-eta.^2),   ...
            -xi.*eta.*(eta+1),          (xi-1/2).*(1-eta.^2),   ...
            -2*xi.*(1-eta.^2)];
        Neta = [xi.*(xi-1).*(eta-1/2)/2,    xi.*(xi+1).*(eta-1/2)/2, ...
            xi.*(xi+1).*(eta+1/2)/2,    xi.*(xi-1).*(eta+1/2)/2, ...
            (1-xi.^2).*(eta-1/2),       xi.*(xi+1).*(-eta),   ...
            (1-xi.^2).*(eta+1/2),       xi.*(xi-1).*(-eta),   ...
            (1-xi.^2).*(-2*eta)];
    otherwise
        error('Shape functions not implemented yet')
end
shapeFunctions(:,:,1) = N';
shapeFunctions(:,:,2) = Nxi';
shapeFunctions(:,:,3) = Neta';
