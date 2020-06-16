function [uint,xyint,indpt] = analyseBoundaryWest(u,X,T,exteriorFaces,refEl,B,nint,proj)

np = size(refEl.NodesCoord,1);
np1d = size(refEl.NodesCoord1d,1);
elements = exteriorFaces(:,1);
faces = exteriorFaces(:,2);
faceNodes = refEl.faceNodes;
nDeg = refEl.degree;
Nf = size(exteriorFaces,1);

% Notable points
pt(1,:) = [2.4410,  0.7870];
pt(2,:) = [1.9360,  0.5840];
pt(3,:) = [1.9510, -0.5920];
pt(4,:) = [2.3580, -0.7550];
pt(5,:) = [2.7977, -0.5128];
pt(6,:) = [2.8180,  0.4868];

% Magnetic field
b = zeros(size(B,1),2);
Br = B(:,1);
Bz = B(:,2);
Bt = B(:,3);
Bmod   = sqrt(Br.^2+Bz.^2+Bt.^2);
b(:,1) = Br./Bmod;
b(:,2) = Bz./Bmod;

% index for the solution
ind = bsxfun(@plus,np*(elements'-1),faceNodes(faces,:)');

% index for the nodes
T_t = T';
indc = T_t(ind);

% Start from point 1
aux_x = abs(X(indc(1,:),1)-pt(1,1))==min(abs(X(indc(1,:),1)-pt(1,1)));
aux_y = abs(X(indc(1,:),2)-pt(1,2))==min(abs(X(indc(1,:),2)-pt(1,2)));
ind1 = find(all([aux_x,aux_y],2),1);
i = ind1;

% Ordinate indeces for speed
indfix = zeros(size(ind));
indfixc = indfix;
% [~,i] = min(ind(:)); 
% i = (i-1)/size(ind,1)+1;
%check
% if min(min(ind))~=ind(1,i),error('something wrong'),end
indfixc(:,1) = indc(:,i);
indfix(:,1) = ind(:,i);
for i = 2:Nf
    aux = indfixc(end,i-1)==indc(1,:);
    indfixc(:,i) = indc(:,aux);
    indfix(:,i) = ind(:,aux);
end

ind = indfix;
indc = indfixc;


% find position of notable points
indpt = zeros(size(pt,1),1);
for i = 1:size(pt,1)
    aux_x = abs(X(indc(1,:),1)-pt(i,1))==min(abs(X(indc(1,:),1)-pt(i,1)));
    aux_y = abs(X(indc(1,:),2)-pt(i,2))==min(abs(X(indc(1,:),2)-pt(i,2)));
    indpt(i) = find(all([aux_x,aux_y],2),1);
end

% Plot only on a part of the boundary
ind = ind(:,indpt(1):indpt(4));
indc= indc(:,indpt(1):indpt(4));
indpt = indpt(1:4);
Nf = size(ind,2);

% Vandermonde matrix
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

x = X(:,1); y = X(:,2);
b1 = b(:,1); b2 = b(:,2);
u_int = shapeFunctions(:,:,1)*u(ind);
x_int = shapeFunctions(:,:,1)*x(indc);
y_int = shapeFunctions(:,:,1)*y(indc);
x_der = shapeFunctions(:,:,2)*x(indc);
y_der = shapeFunctions(:,:,2)*y(indc);
b1e    = shapeFunctions(:,:,1)*b1(indc);
b2e    = shapeFunctions(:,:,1)*b2(indc);

% Loop in faces

% ind_vect =1:nint/2:nint;
% point_x = zeros(nint,Nf);
% point_y = point_x;
uint = zeros(nint,Nf);
for iFace = 1:Nf

    % Integration weight
    xyDerNorm = sqrt(x_der(:,iFace).^2 + y_der(:,iFace).^2);

    % Tangent vector
    t = [ x_der(:,iFace)./xyDerNorm y_der(:,iFace)./xyDerNorm];
    n = [t(:,2) -t(:,1)]; 
    
    bn = b1e(:,iFace).*n(:,1)+b2e(:,iFace).*n(:,2);

    % p face
    pf = u_int(:,iFace);
    
    if proj
        pf = pf.*bn;
    end
    
    uint(:,iFace) =  pf;
%     point_x(:,iFace) = x_int(:,iFace) + vect_x;
%     point_y(:,iFace) = y_int(:,iFace) + vect_y;
end

uint = uint(:);
xyint = [x_int(:),y_int(:)];
indpt = (indpt-1)*nint+1;


