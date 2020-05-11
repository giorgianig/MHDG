function points_rst = mapXiEtaZeta2rst(points)
%
% points_rst = mapXiEtaZeta2rst(points)
%

nOfSpatialDimensions = size(points, 2);

switch nOfSpatialDimensions
    case 2
        points_rst = mapXiEtaZeta2rst_2D(points);
    case 3
        points_rst = mapXiEtaZeta2rst_3D(points);
    otherwise
        error('mapXiEtaZeta2rst: wrong nOfSpatialDimensions')
end

function points_rst = mapXiEtaZeta2rst_2D(points)

nOfPoints = size(points, 1);
points_rst      = zeros(nOfPoints,2);
points_rst(:,1) = 2*(1+points(:,1))./(1-points(:,2)) - 1;
points_rst(:,2) = points(:,2);

function points_rst = mapXiEtaZeta2rst_3D(points)

nOfPoints = size(points, 1);
points_rst      = zeros(nOfPoints,3);
points_rst(:,1) = -2*(1+points(:,1))./(points(:,2)+points(:,3)) - 1;
points_rst(:,2) =  2*(1+points(:,2))./(1-points(:,3)) - 1;
points_rst(:,3) = points(:,3);