function X = linearMapping(vertCoord,xiVector)
%
% X = isoparTransf(vertCoord,xiVector)
%
% Linear mapping between local and cartesian ccordinates
%
% Input:
% vertCoord: vertexs of the element
% xiVector:        point in local coordinates
%
% Output:
% X:         point in cartesian coordinates

nOfSpatialDimensions = size(xiVector,2);

switch nOfSpatialDimensions
    case 1
        N = linearShapeFunctions1D(xiVector);
    case 2
        N = linearShapeFunctions2D(xiVector);
    case 3
        N = linearShapeFunctions3D(xiVector);
end

%         
% eval(['N = linearShapeFunctions',int2str(nOfSpatialDimensions),'D(xiVector);']);
X = N*vertCoord;



function N = linearShapeFunctions1D(xiVector)
% Reference interval is [-1,1]
N = 0.5*[1-xiVector, 1+xiVector];

function N = linearShapeFunctions2D(xiVector)
% Reference triangle is [-1,-1; 1,-1; -1,1]
xi = xiVector(:,1);
eta = xiVector(:,2);
N = 0.5*[-xi-eta, 1+xi, 1+eta];

function N = linearShapeFunctions3D(xiVector)
% Reference tetrahedron is [-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
xi = xiVector(:,1);
eta = xiVector(:,2);
zeta = xiVector(:,3);
N = 0.5*[-1-xi-eta-zeta, 1+xi, 1+eta, 1+zeta];
