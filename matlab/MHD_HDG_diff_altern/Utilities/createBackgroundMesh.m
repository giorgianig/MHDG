%% Create Background mesh
clear
close all
global neq 



% Load files to determine the mesh size from a previous solution 
load ../../Meshes/Meshes_2D/Circle_lim_h01_refCor_P10.mat         % the Mesh: must match the solution file
load Saves/u_Circle_lim_h01_refCor_P10.mat_Diffusion_0.00625.mat  % the solution file: must match the mesh
p = 10;                                                           % degree of the reference element used in the solution
nDegRef = 1;                                                     % degree used to create the linear mesh
neq = 2;

%% ********* Start work *********
%%_______________________________

%% Number of elements and element nodes in the original mesh
Ne = size(T,1);
np = 0.5*(p+1)*(p+2);

%% Initialize mesh size (nodal value in the original mesh)
msize_nod = zeros(numel(T),1);

%% Reference element
refEl = createReferenceElement(1,np);

%% Original mesh element size
meshSizeBase = findElementSize(X,T);

%% Nodal connectivity
disp('Generating nodal connectivity')
nNodes = max(max(T));
N = zeros(nNodes,10);
nn = ones(nNodes,1);
for ielem = 1:size(T,1)
    Te = T(ielem,:);
    nn_Te = nn(Te);
    for kk = 1:np
        N(Te(kk),nn_Te(kk)) = ielem;
    end
    nn(Te) = nn(Te) +1;
end
N(:,max(nn):end) = [];

%% Linear reference element for interpolation purpose
nodes = refEl.NodesCoord;
npoints = size(nodes,1);
refEl_lin = createReferenceElement(1,3);

% Vandermonde matrix
coordRef = refEl_lin.NodesCoord;
nOfNodes = size(coordRef,1);
V = Vandermonde_LP(1,coordRef);
invV = inv(V');

% Compute shape functions at interpolation points
shapeFunctions = zeros(npoints,nOfNodes);
for ipoint = 1:npoints
    p = orthopoly2D(nodes(ipoint,:),1);
    shapeFunctions(ipoint,:) = (invV*p)';
end

N_linear = shapeFunctions;




%% Find element size based on the solution
disp('Find element size based on the solution')


% Here I set the element size based on the shock capturing coefficient
msize_elem = meshSizeBase./10.^(findCoeffShockCaptur(u0,X,T,refEl,-7.5));
% msize_elem(msize_elem>meshSizeBase)=meshSizeBase; % we can eliminate this? 

figure
plotConstantElementalValue(X,T,msize_elem,refEl)

 
for ielem = 1:Ne
    
    Te = T(ielem,:);
    Xe = X(Te,:);
    indnodes = (ielem-1)*np + (1:np);
        
    els = N(Te(1:3),:);
    aux_eps_nodal_n = zeros(3,1);
    for i = 1:size(els,1)
        ind = els(i,els(i,:)~=0);
        aux_eps_nodal_n(i) = min(msize_elem(ind));
    end
    msize_nod(indnodes) = N_linear*aux_eps_nodal_n;
end



%% Plotting element (equal spaced points)
disp('Generating plotting element')
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



%% Generate linear mesh and values of mesh size
disp('Generating linear mesh')
nOfElemTriRef = size(elemTriRef,1);
tri = zeros(nOfElemTriRef*Ne,3);
indexElem = 0;
msize = zeros(Ne*npoints,1);
Xplot = zeros(Ne*npoints,2);
for ielem = 1:Ne
    Te = T(ielem,:);
    ind = ( (ielem-1)*nOfNodes+1 ) : (ielem*nOfNodes);
    ueplot = shapeFunctions*msize_nod(ind);
    Xeplot = shapeFunctions*X(Te,:);
    for ielemRef = 1:nOfElemTriRef
        indexElemRef = indexElem + ielemRef;
        tri(indexElemRef,:) = elemTriRef(ielemRef,:)+(ielem-1)*npoints;
        Xplot(tri(indexElemRef,:),:) = Xeplot(elemTriRef(ielemRef,:),:);
        msize(tri(indexElemRef,:)) = ueplot(elemTriRef(ielemRef,:));
    end
    indexElem = indexElem + nOfElemTriRef;
end

figure,
patch('Faces',tri,'Vertices',Xplot,'FaceVertexCData',msize,...
    'FaceColor','interp','EdgeAlpha',0);
axis equal
colormap('jet')
colorbar('location','East');
hold on
plotMesh(X,T)

%% Open file and write 2 lines
disp('Writing background mesh')
fid = fopen('BackgroundMesh.sub','w');
fprintf(fid,'%d\n',size(Xplot,1));
fprintf(fid,'%d\n',size(tri,1));


%% Write nodes
for ip = 1:size(Xplot,1)

    fprintf(fid,'%16.15e %16.15e %16.15e \n',[Xplot(ip,:) msize(ip)]);

end

%% Write elements
for ielem = 1:size(tri,1)
    fprintf(fid,'%d %d %d \n',tri(ielem,:)  );
end

%% Close file
fclose(fid);

disp('Done!')