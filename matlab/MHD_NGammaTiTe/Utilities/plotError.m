function plotError(X,T,u,refEl,nDegRef)

global Mesh

X = Mesh.lscale*X;
u = transpose(reshape(u,[2 numel(u)/2]));
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
coordRef = refEl.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = refEl.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');

% Compute shape functions at interpolation points
shapeFunctions = zeros(npoints,nOfNodes);
for ipoint = 1:npoints
    p = orthopoly2D(nodes(ipoint,:),nDeg);
    shapeFunctions(ipoint,:) = (invV*p)';
end

%Delaunay's solution mesh   
nOfElemTriRef = size(elemTriRef,1);
nEl = size(T,1);
tri = zeros(nOfElemTriRef*nEl,3);
indexElem = 0;
uplot_dens = zeros(nEl*npoints,1);
uplot_mome = zeros(nEl*npoints,1);

Xplot = zeros(nEl*npoints,2);
for ielem = 1:nEl
    Te = T(ielem,:);
    if size(u,1)==size(X,1)
        % Continuous solution
        ueplot_dens = shapeFunctions*u(Te,1);
        ueplot_mome = shapeFunctions*u(Te,2);
    elseif size(u,1) == numel(T)
        % Discontinuous solution
        ind = ( (ielem-1)*nOfNodes+1 ) : (ielem*nOfNodes);
        ueplot_dens = shapeFunctions*u(ind,1);
        ueplot_mome = shapeFunctions*u(ind,2);
    else
        error('Dimension of u is wrong')
    end
    Xeplot = shapeFunctions*X(Te,:);
    for ielemRef = 1:nOfElemTriRef
        indexElemRef = indexElem + ielemRef;
        tri(indexElemRef,:) = elemTriRef(ielemRef,:)+(ielem-1)*npoints;
        Xplot(tri(indexElemRef,:),:) = Xeplot(elemTriRef(ielemRef,:),:);
        uan = reshape(analyticalSolution(Xeplot(elemTriRef(ielemRef,:),:)),...
            size(Xeplot(elemTriRef(ielemRef,:),:)));
        uplot_dens(tri(indexElemRef,:)) = ueplot_dens(elemTriRef(ielemRef,:))-uan(:,1);
        uplot_mome(tri(indexElemRef,:)) = ueplot_mome(elemTriRef(ielemRef,:))-uan(:,2);
    end
    indexElem = indexElem + nOfElemTriRef;
end

%Plot
figure
patch('Faces',tri,'Vertices',Xplot,'FaceVertexCData',uplot_dens,...
    'FaceColor','interp','EdgeAlpha',0);
axis equal
colormap('jet')
colorbar('location','East');
title('Error in the density')
figure
patch('Faces',tri,'Vertices',Xplot,'FaceVertexCData',uplot_mome,...
    'FaceColor','interp','EdgeAlpha',0);
axis equal
colormap('jet')
colorbar('location','East');
title('Error in the momentum')
