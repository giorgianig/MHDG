function plotConstantElementalValue(X,T,val,refEl)

global Mesh

if ~isempty(Mesh)

    X = Mesh.lscale*X;
end

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
uplot = zeros(nEl*npoints,1);
Xplot = zeros(nEl*npoints,2);
for ielem = 1:nEl
    Te = T(ielem,:);
    if val(ielem)==0
    ueplot = NaN*ones(npoints,1);
    else
    ueplot = val(ielem)*ones(npoints,1);
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

%Plot
patch('Faces',tri,'Vertices',Xplot,'FaceVertexCData',uplot,...
    'FaceColor',[0.8 0.8 0.8],'EdgeAlpha',0);
hold on
patch('Faces',tri,'Vertices',Xplot,'FaceVertexCData',uplot,...
    'FaceColor','interp','EdgeAlpha',0);
axis equal
colormap('jet')
colorbar('location','East');
hold on,plotMesh(X,T)