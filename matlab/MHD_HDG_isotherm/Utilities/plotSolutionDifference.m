function plotSolutionDifference(X,T,u,refEl,component,udiff,meshdiff,nDegRef)

global Mesh

X = Mesh.lscale*X;

% reshape solution
u = transpose(reshape(u,2,numel(T)));
if component==1
    u = u(:,1);
elseif component==2
    u = u(:,2)./u(:,1);
end

if size(u,2) >1
    error('Only 1d solutions')
end
% Check input
if nargin == 6
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


% Mesh for the other solution
Sstart  = load(udiff);
MStart = load(meshdiff);
refElstart = createReferenceElement(1,size(MStart.T,2),[]);
u2diff = transpose(reshape(Sstart.u0,2,numel(MStart.T)));
if component==1
    u2diff = u2diff(:,1);
elseif component==2
    u2diff = u2diff(:,2)./u2diff(:,1);
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

disp('Projecting...')
u_proj = evalDGapproximationAtPoints(Xplot,u2diff,MStart.X,MStart.T,refElstart);

%Plot
patch('Faces',tri,'Vertices',Xplot,'FaceVertexCData',(uplot-u_proj)./u_proj,...
    'FaceColor','interp','EdgeAlpha',0);
% patch('Faces',tri,'Vertices',Xplot,'FaceVertexCData',u_proj,...
%     'FaceColor','interp','EdgeAlpha',0);
axis equal
colormap('jet')
colorbar('location','East','fontname','times new roman');

