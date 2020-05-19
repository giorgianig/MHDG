function plotSolutionOnBoundary(u,X,T,exteriorFaces,refEl,B,nint,scale,proj)

np = size(refEl.NodesCoord,1);
np1d = size(refEl.NodesCoord1d,1);
elements = exteriorFaces(:,1);
faces = exteriorFaces(:,2);
faceNodes = refEl.faceNodes;
nDeg = refEl.degree;
Nf = size(exteriorFaces,1);

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

% Ordinate indeces for speed
indfix = zeros(size(ind));
indfixc = indfix;
[~,i] = min(ind(:)); 
i = (i-1)/size(ind,1)+1;
%check
if min(min(ind))~=ind(1,i),error('something wrong'),end
indfixc(:,1) = indc(:,i);
indfix(:,1) = ind(:,i);
for i = 2:Nf
    aux = indfixc(end,i-1)==indc(1,:);
    indfixc(:,i) = indc(:,aux);
    indfix(:,i) = ind(:,aux);
end

ind = indfix;
indc = indfixc;



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
point_x = zeros(nint,Nf);
point_y = point_x;

for iFace = 1:Nf

    % Integration weight
    xyDerNorm = sqrt(x_der(:,iFace).^2 + y_der(:,iFace).^2);

    % Tangent vector
    t = [ x_der(:,iFace)./xyDerNorm y_der(:,iFace)./xyDerNorm];
    n = [t(:,2) -t(:,1)]; 
    
    bn = b1e(:,iFace).*n(:,1)+b2e(:,iFace).*n(:,2);

    % p face
    pf = u_int(:,iFace);
    pf = pf*scale;
    
    if proj
        pf = pf.*bn;
    end
    % draw points
%     pos_pres = pf>=0;
%     vect_x =  n(:,1).*abs(pf);
%     vect_y =  n(:,2).*abs(pf);
    vect_x =  n(:,1).*pf;
    vect_y =  n(:,2).*pf;
    point_x(:,iFace) = x_int(:,iFace) + vect_x;
    point_y(:,iFace) = y_int(:,iFace) + vect_y;

    % plot bubble
%     plot(point_x(~pos_pres),point_y(~pos_pres),'b','LineWidth',2)
%     plot(point_x(pos_pres),point_y(pos_pres),'r','LineWidth',2)
%     plot(x_int(:,iFace),y_int(:,iFace),'k','LineWidth',3);
%     plot(point_x,point_y,'r','LineWidth',2)
%     plot(x_int(:,iFace),y_int(:,iFace),'k','LineWidth',3);

    % plot vectors
%     xv = x_int(ind_vect,iFace);
%     yv = y_int(ind_vect,iFace);
%     px = point_x(ind_vect);
%     py = point_y(ind_vect);
%     vx = vect_x(ind_vect);
%     vy = vect_y(ind_vect);
%     pp = pos_pres(ind_vect);
%     if any(~pp)
%         quiver(xv(~pp),yv(~pp), vx(~pp), vy(~pp),0,'b')
%     end
%     if any(pp)
%         quiver(px(pp),py(pp), -vx(pp), -vy(pp),0,'r')
%     end
end
plot(point_x(:),point_y(:),'r','LineWidth',2)
hold on
plot(x_int(:),y_int(:),'k','LineWidth',3);
axis equal