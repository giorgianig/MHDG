function plotCpOnAirfoil(press,X,T,infoFaces,refEl,nint)

exteriorFaces = [infoFaces.exteriorFaces_WALL_UP;infoFaces.exteriorFaces_WALL_DOWN];
np = size(refEl.NodesCoord,1);
np1d = size(refEl.NodesCoord1d,1);
elements = exteriorFaces(:,1);
faces = exteriorFaces(:,2);
faceNodes = refEl.faceNodes;
nDeg = refEl.degree;

% index for the pressure
ind = bsxfun(@plus,np*(elements'-1),faceNodes(faces,:)');

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

T_t = T';
x = X(:,1); y = X(:,2);
p_int = shapeFunctions(:,:,1)*press(ind);
x_int = shapeFunctions(:,:,1)*x(T_t(ind));
y_int = shapeFunctions(:,:,1)*y(T_t(ind));
x_der = shapeFunctions(:,:,2)*x(T_t(ind));
y_der = shapeFunctions(:,:,2)*y(T_t(ind));

% Loop in faces
Nf = size(exteriorFaces,1);
figure
hold on
ind_vect =1:nint/2:nint;
for iFace = 1:Nf

    % Integration weight
    xyDerNorm = sqrt(x_der(:,iFace).^2 + y_der(:,iFace).^2);

    % Tangent vector
    t = [ x_der(:,iFace)./xyDerNorm y_der(:,iFace)./xyDerNorm];
    n = [-t(:,2) t(:,1)]; % the normal is the opposite of the usual one

    % p face
    pf = p_int(:,iFace);

    % draw points
    pos_pres = pf>=0;
    vect_x =  n(:,1).*abs(pf);
    vect_y =  n(:,2).*abs(pf);
    point_x = x_int(:,iFace) + vect_x;
    point_y = y_int(:,iFace) + vect_y;

    % plot bubble
    plot(point_x(~pos_pres),point_y(~pos_pres),'b','LineWidth',2)
    plot(point_x(pos_pres),point_y(pos_pres),'r','LineWidth',2)
    plot(x_int(:,iFace),y_int(:,iFace),'k','LineWidth',3);

    % plot vectors
    xv = x_int(ind_vect,iFace);
    yv = y_int(ind_vect,iFace);
    px = point_x(ind_vect);
    py = point_y(ind_vect);
    vx = vect_x(ind_vect);
    vy = vect_y(ind_vect);
    pp = pos_pres(ind_vect);
    if any(~pp)
        quiver(xv(~pp),yv(~pp), vx(~pp), vy(~pp),0,'b')
    end
    if any(pp)
        quiver(px(pp),py(pp), -vx(pp), -vy(pp),0,'r')
    end
end
axis equal