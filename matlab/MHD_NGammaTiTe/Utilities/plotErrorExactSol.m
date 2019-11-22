function plotErrorExactSol(X,T,u,referenceElement,nDegRef)

% reshape u
u = transpose(reshape(u,2,numel(u)/2));

% Check input
if nargin == 4
    nDegRef = 20;
end

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
    ind = (ielem-1)*nOfElementNodes+1:ielem*nOfElementNodes;
    u_ex = analyticalSolution(Xplot);
    uplot(:,1) = shapeFunctions*u(ind,1);
    uplot(:,2) = shapeFunctions*u(ind,2);
    nodalerr = (u_ex-uplot);
    errPlot = sqrt(nodalerr(:,1).^2 + nodalerr(:,2).^2);
    
    % Plot interpolated solution in the element
    hold on
    patchHandle(ielem) = patch('Faces',elemTriRef,'Vertices',Xplot,'FaceVertexCData',errPlot,...
        'FaceColor','interp','EdgeAlpha',0);
    hold off
end
axis equal
colorbar

