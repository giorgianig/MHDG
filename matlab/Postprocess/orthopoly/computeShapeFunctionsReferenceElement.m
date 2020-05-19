function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement(nDeg,coord,nOfGaussPoints,varargin)
%
% [shapeFunctions,gaussWeights]=computeShapeFunctionsReferenceElement(nDeg,
% coord)
%
% Function to compute the shape functions (& derivatives) at Gauss points
%
% Input:
% nDeg:  degree of interpolation
% coord: nodal coordinates at the reference element
% nOfGaussPoints: n� of gauss points of the 1D quadrature
% elementType (optional): 0 for quadrilateral, 1 for triangle. If it isn't
%                         given only triangle or 1D elements are considered
%
% Output:
% shapeFunctions: shape functions evaluated at the gauss points
%                 size is nOfNodes X nOfGauss X (nsd + 1)
%                 nsd+1 because shape function (1)
%                 and derivatives (nsd) are stored
% gaussWeights:   weights
%

if ~isempty(varargin)
    elementType = varargin{:};
    if elementType == 0
        [shapeFunctions,gaussWeights,gaussPoints]=...
            computeShapeFunctionsQua(nDeg,coord,nOfGaussPoints);
        return
    end
end

nsd = size(coord,2);
if nsd==1
    [shapeFunctions,gaussWeights,gaussPoints]=...
        computeShapeFunctionsReferenceElement1D(nDeg,coord,nOfGaussPoints);
elseif nsd==2
    [shapeFunctions,gaussWeights,gaussPoints]=...
        computeShapeFunctionsReferenceElement2D(nDeg,coord,nOfGaussPoints);
elseif nsd==3
    [shapeFunctions,gaussWeights,gaussPoints]=...
        computeShapeFunctionsReferenceElement3D(nDeg,coord,nOfGaussPoints);
else
    error('wrong nsd in computeShapeFunctionsReferenceElement')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement1D(nDeg,coord,nOfGaussPoints)

%number of nodes/polynomials
nOfNodes = nDeg+1;
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsReferenceElement1D')
end

[z,w] = gaussLegendre(nOfGaussPoints,-1,1);
nOfGauss = length(w);

%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

shapeFunctions = zeros(nOfNodes,nOfGauss,2);
gaussWeights = w';
gaussPoints = zeros(nOfGauss,1);

%Integration over [-1,1]
for i = 1:nOfGauss
    x = z(i);
    [p,p_xi] = orthopoly1D_deriv(x,nDeg);
    N = U\(L\(P*[p,p_xi]));
    shapeFunctions(:,i,1) = N(:,1)';
    shapeFunctions(:,i,2) = N(:,2)';
    % only for PFEM
    gaussPoints(i) = x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement2D(nDeg,coord,nOfGaussPoints)

%number of nodes/polynomials
nOfNodes = (nDeg+1)*(nDeg+2)/2;
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsReferenceElement2D')
end
if nDeg <12 && nOfGaussPoints-2<11
    if isequal(nOfGaussPoints,nDeg+2)
        % NEW CALL FOR CUBATURES
        switch nDeg
            case 1
                OrderCubature = 5; % 7 pdG
            case 2
                OrderCubature = 10;
            case 3
                OrderCubature = 10;
            case 4
                OrderCubature = 15;
            case 5
                OrderCubature = 15;
            case 6
                OrderCubature = 15;
            case 7
                OrderCubature = 15;
            case {8,9,10,11}
                OrderCubature = 25;
        end
    else
        switch nOfGaussPoints-2
            case 1
                OrderCubature = 5; % 7 pdG
            case 2
                OrderCubature = 10;
            case 3
                OrderCubature = 10;
            case 4
                OrderCubature = 15;
            case 5
                OrderCubature = 15;
            case 6
                OrderCubature = 15;
            case 7
                OrderCubature = 15;
            case {8,9,10,11}
                OrderCubature = 25;
        end
    end
    [z,w] = GaussLegendreCubature2D(OrderCubature);
    w = 2*w; z = 2*z -1; %mapping onto the normal reference triangle
    nIP = length(w); %number of integration points
    nOfGauss = nIP; % for the cubature
    %Vandermonde matrix
    V = Vandermonde_LP(nDeg,coord);
    [L,U,P] = lu(V');

    shapeFunctions = zeros(nOfNodes,nOfGauss,3);
    gaussWeights = zeros(nOfGauss,1);
    gaussPoints = zeros(nOfGauss,2);

    %Integration over [-1,1]^2 using the cubature
    for i = 1:nIP
        x = z(i,:); % (r,s) coordinates
        [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(x,nDeg);
        N = U\(L\(P*[p,p_xi,p_eta]));
        shapeFunctions(:,i,1) = N(:,1)';
        shapeFunctions(:,i,2) = N(:,2)';
        shapeFunctions(:,i,3) = N(:,3)';
        gaussWeights(i) = w(i);
        gaussPoints(i,:) = x;
    end
else
    %     nOfGaussPoints = nDeg + 2;
    %OLD CALL
    [z,w] = gaussLegendre(nOfGaussPoints,-1,1);
    nIP = length(w); %number of integration points in each direction
    nOfGauss = nIP^2;

    %Vandermonde matrix
    V = Vandermonde_LP(nDeg,coord);
    [L,U,P] = lu(V');

    shapeFunctions = zeros(nOfNodes,nOfGauss,3);
    gaussWeights = zeros(nOfGauss,1);
    gaussPoints = zeros(nOfGauss,2);

    iGauss = 1;
    %Integration over [-1,1]^2
    for i = 1:nIP
        for j = 1:nIP
            x = [z(i),z(j)]; % (r,s) coordinates
            [p,p_xi,p_eta] = orthopoly2D_deriv_rst(x,nDeg);
            N = U\(L\(P*[p,p_xi,p_eta]));
            shapeFunctions(:,iGauss,1) = N(:,1)';
            shapeFunctions(:,iGauss,2) = N(:,2)';
            shapeFunctions(:,iGauss,3) = N(:,3)';
            gaussWeights(iGauss) =(w(i)*w(j))*(1-x(2))/2;
            % only for PFEM
            r = x(1); s = x(2);
            xi = (1+r)*(1-s)/2-1;
            gaussPoints(iGauss,:) = [xi, s];
            iGauss = iGauss + 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement3D(nDeg,coord,nOfGaussPoints)

%number of nodes/polynomials
nOfNodes = (nDeg+1)*(nDeg+2)*(nDeg+3)/6;
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsReferenceElement3D')
end

[z,w] = gaussLegendre(nOfGaussPoints,-1,1);
nIP = length(w); %number of integration points in each direction
nOfGauss = nIP^3;

%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

shapeFunctions = zeros(nOfNodes,nOfGauss,4);
gaussWeights = zeros(nOfGauss,1);
gaussPoints = zeros(nOfGauss,3);

iGauss = 1;
%Integration over [-1,1]^3
for i = 1:nIP
    for j = 1:nIP
        for k = 1:nIP
            x = [z(i),z(j),z(k)]; % (r,s,t) coordinates
            [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_rst(x,nDeg);
            N = U\(L\(P*[p,p_xi,p_eta,p_zeta]));
            shapeFunctions(:,iGauss,1) = N(:,1)';
            shapeFunctions(:,iGauss,2) = N(:,2)';
            shapeFunctions(:,iGauss,3) = N(:,3)';
            shapeFunctions(:,iGauss,4) = N(:,4)';
            gaussWeights(iGauss) =(w(i)*w(j)*w(k))*((1-x(2))/2)*((1-x(3))/2)^2;
            % only for PFEM
            r = x(1); s = x(2); t = x(3);
            eta = (1/2)*(s-s*t-1-t);
            xi = -(1/2)*(r+1)*(eta+t)-1;
            gaussPoints(iGauss,:) = [xi, eta, t];
            iGauss = iGauss + 1;
        end
    end
end
