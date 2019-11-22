function shapeFunctions=computeShapeFunctionsAtPoints(nDeg,coord,points)
%
% [shapeFunctions,gaussWeights]=computeShapeFunctionsReferenceElement(nDeg,
% coord)
%
% Function to compute the shape functions (& derivatives) at given points
%
% Input:
% nDeg:  degree of interpolation
% coord: nodal coordinates at the reference element
% points: points were shape functions are to be evaluated
%
% Output:
% shapeFunctions: shape functions evaluated at the gauss points
%                 size is nOfNodes X nOfGauss X (nsd + 1)
%                 nsd+1 because shape function (1)
%                 and derivatives (nsd) are stored
%

nsd = size(coord,2);
if nsd==1
    shapeFunctions=computeShapeFunctions1D(nDeg,coord,points);
elseif nsd==2
    shapeFunctions=computeShapeFunctions2D(nDeg,coord,points);
else
    error('wrong nsd in computeShapeFunctionsReferenceElement')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shapeFunctions=computeShapeFunctions1D(nDeg,coord,points)

%number of nodes/polynomials
nOfNodes = nDeg+1;
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsReferenceElement1D')
end

z = points;
nOfGauss = length(z);

%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

shapeFunctions = zeros(nOfNodes,nOfGauss,2);

%Integration over [-1,1]
for i = 1:nOfGauss
    x = z(i);
    [p,p_xi] = orthopoly1D_deriv(x,nDeg);
    N = U\(L\(P*[p,p_xi]));
    shapeFunctions(:,i,1) = N(:,1)';
    shapeFunctions(:,i,2) = N(:,2)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shapeFunctions=computeShapeFunctions2D(nDeg,coord,points)

%number of nodes/polynomials
nOfNodes = (nDeg+1)*(nDeg+2)/2;
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsReferenceElement2D')
end

z = points;
nIP = size(z,1); %number of integration points
nOfGauss = nIP; % for the cubature
%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

shapeFunctions = zeros(nOfNodes,nOfGauss,3);

%Integration over [-1,1]^2 using the cubature
for i = 1:nIP
    x = z(i,:); 
    [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(x,nDeg);
    N = U\(L\(P*[p,p_xi,p_eta]));
    shapeFunctions(:,i,1) = N(:,1)';
    shapeFunctions(:,i,2) = N(:,2)';
    shapeFunctions(:,i,3) = N(:,3)';
end


