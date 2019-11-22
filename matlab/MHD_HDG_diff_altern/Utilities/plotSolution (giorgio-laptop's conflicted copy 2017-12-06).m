function plotSolution(X,T,u,refEl,nDegRef,cont,logar,plotSuper)

global Mesh

X = Mesh.lscale*X;

if size(u,2) >1
    error('Only 1d solutions')
end
% Check input
if nargin == 4
    nDegRef = 20;
    cont = 0;
    logar=0;
    plotSuper = 0;
elseif nargin == 5
    cont = 0;
    logar=0;
    plotSuper = 0;
elseif nargin == 6
    logar=0;
    plotSuper = 0;
elseif nargin == 7
    plotSuper = 0;
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
uplot = zeros(nEl*npoints,1);
Xplot = zeros(nEl*npoints,2);
for ielem = 1:nEl
    Te = T(ielem,:);
    if numel(u)==size(X,1)
        % Continuous solution
        ueplot = shapeFunctions*u(Te);
    elseif numel(u) == numel(T)
        % Discontinuous solution
        ind = ( (ielem-1)*nOfNodes+1 ) : (ielem*nOfNodes);
        ueplot = shapeFunctions*u(ind);
    else
        error('Dimension of u is wrong')
    end
    Xeplot = shapeFunctions*X(Te,:);
    for ielemRef = 1:nOfElemTriRef
        indexElemRef = indexElem + ielemRef;
        tri(indexElemRef,:) = elemTriRef(ielemRef,:)+(ielem-1)*npoints;
        Xplot(tri(indexElemRef,:),:) = Xeplot(elemTriRef(ielemRef,:),:);
        uplot(tri(indexElemRef,:)) = ueplot(elemTriRef(ielemRef,:));
    end
    indexElem = indexElem + nOfElemTriRef;
end

% if plotSuper
% uplot(abs(uplot)<1) = NaN;
% end


% uplot(abs(uplot)<1e-4) = NaN;

if cont
    tricontour(Xplot,tri,uplot,30);
    hold on
    %     for iname = 1:numel(Mesh.boundaryNames)
    %         boundary.(Mesh.boundaryNames{iname}) = Mesh.(Mesh.boundaryNames{iname});
    %     end
    %     plotBoundary(X,boundary)
    axis equal
    colormap('jet')
    colorbar('location','East');
else
    
    
    %Plot
    if logar
        patch('Faces',tri,'Vertices',Xplot,'FaceVertexCData',log10(uplot),...
            'FaceColor','interp','EdgeAlpha',0);
    else
        patch('Faces',tri,'Vertices',Xplot,'FaceVertexCData',uplot,...
            'FaceColor','interp','EdgeAlpha',0);
    end
    axis equal
    colormap('jet')
    colorbar('location','East','fontname','times new roman',...
        'axislocation','out');
    assi = axis;
end
if plotSuper
    C1= tricontour(Xplot,tri,uplot,[1 1]);
    hold on
    C2 = tricontour(Xplot,tri,uplot,[-1 -1]);
    plot(C1(1,:),C1(2,:),'k.')
    plot(C2(1,:),C2(2,:),'k.')
    load /home/giorgio/Dropbox/Matlab/Meshes/Meshes_2D/WEST/WEST_wall.mat
    plot(Rwall,Zwall,'k-'); axis equal
    axis(assi)
    %     plotMesh(X,T)
end
