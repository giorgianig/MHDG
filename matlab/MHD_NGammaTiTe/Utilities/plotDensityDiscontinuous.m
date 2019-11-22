function plotDensityDiscontinuous(X,T,sol,referenceElement,nDegRef)

% Check input
if nargin == 4
    nDegRef = 20;
end

sol = transpose(reshape(sol,2,numel(sol)/2));
u = sol(:,1);

% Plotting element (equal spaced points)
nodes = [];
h = 1/nDegRef;
for j = 0:nDegRef
    i = (0:nDegRef-j)';
    aux = j*ones(size(i));
    nodes = [nodes; [i aux]*h];
end
nodes = 2*nodes - 1;
npoints = size(nodes,1);

% Delaunay triangulation of plotting element
elemTriRef = delaunayn(nodes);

% Vandermonde matrix
coordRef = referenceElement.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = referenceElement.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');

nOfElementNodes = size(coordRef,1);
nOfElements = size(T,1);

% Compute shape functions at interpolation points
shapeFunctions = zeros(npoints,nOfNodes);
for ipoint = 1:npoints
    p = orthopoly2D(nodes(ipoint,:),nDeg);
    shapeFunctions(ipoint,:) = (invV*p)';
end

% Loop in elements
patchHandle = zeros(1,nOfElements);
for ielem = 1:nOfElements

    % Interpolate solution and position at interpolation points
    Te = T(ielem,:);
    Xplot = shapeFunctions*X(Te,:);
    ind = ( (ielem-1)*nOfElementNodes+1 ) : (ielem*nOfElementNodes);
    uplot = shapeFunctions*u(ind);

    % Plot interpolated solution in the element
    hold on
    patchHandle(ielem) = patch('Faces',elemTriRef,'Vertices',Xplot,'FaceVertexCData',uplot,...
        'FaceColor','interp','EdgeAlpha',0);
    hold off
end
axis equal
colorbar

