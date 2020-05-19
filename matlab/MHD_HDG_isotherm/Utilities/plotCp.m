function plotCp(press,X,T,infoFaces,refEl,nint,color)

exteriorFaces = [infoFaces.exteriorFaces_WALL_UP;infoFaces.exteriorFaces_WALL_DOWN];
np = size(refEl.NodesCoord,1);
elements = exteriorFaces(:,1);
faces = exteriorFaces(:,2);
faceNodes = refEl.faceNodes;
np1d = size(faceNodes,2);

% index for the pressure
ind = bsxfun(@plus,np*(elements'-1),faceNodes(faces,:)');

% index for the coordinates
T_t = T';
ind_x = T_t(ind);

% Vandermonde matrix
nDeg = refEl.degree;
V = Vandermonde_LP(nDeg,refEl.NodesCoord1d);
[L,U,P] = lu(V');

% Shape functions
z = -1 : 2/(nint-1) : 1;
shapeFunctions = zeros(nint,np1d,2);
for i = 1:nint
    x = z(i);
    [p,p_xi] = orthopoly1D_deriv(x,nDeg);
    N = U\(L\(P*[p,p_xi]));
    shapeFunctions(i,:,1) = N(:,1);
    shapeFunctions(i,:,2) = N(:,2);
end

% interpolate
x = X(:,1);
p_int = shapeFunctions(:,:,1)*press(ind);
x_int = shapeFunctions(:,:,1)*x(ind_x);

plot(x_int,-2*p_int,[color '--'],'LineWidth',2)
xlim([0 1])
grid on