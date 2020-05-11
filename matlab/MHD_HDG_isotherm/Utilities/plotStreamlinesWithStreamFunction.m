function plotStreamlinesWithStreamFunction(X,T,u,refEl,boundary,np,varargin)

if isempty(varargin)
   zoom = [];
else
    zoom = varargin{1};
end

Ne = size(T,1);     % number of elements
Nn = size(X,1);     % number of nodes
Nen = size(T,2);    % number of element nodes
allocation = Nen^2*Ne;
f = zeros(Nn,1);
I = zeros(allocation,1);
J = I;
K = I;
aux_ones = ones(1,Nen);

% velocity components
u = transpose(reshape(u,2,0.5*numel(u)));

% solve the Laplace problem for the stream function
for ielem = 1:Ne
    
    Te = T(ielem,:);
    Xe = X(Te,:);
    ue = u(Te,:);
    % elemental matrix
    [Ke fe] = computeElementalMatrix(Xe,refEl,ue);
    
    % Assembling
    f(Te) = f(Te) + fe;
    Te_transp = transpose(Te);
    aux_row = Te_transp(:,aux_ones);
    aux_col = Te(aux_ones,:);
    indK = (ielem-1)*Nen^2+1:ielem*Nen^2;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    K(indK) = Ke;
end

KK = sparse(I,J,K,Nn,Nn);

% boundary conditions
boundaryNames = fieldnames(boundary);
for ibound = 1:numel(boundaryNames)
    name = boundaryNames{ibound};
    Tb = boundary.(name);
    for iface = 1:size(Tb,1)
        Tf = Tb(iface,:);
        Xf = X(Tf,:);
        uf = u(Tf,:);
        ff = - computeElementalBoundaryVector(Xf,refEl,uf);
        f(Tf) = f(Tf) + ff; 
    end
end

% set one value of the stream function
KK(1,:) = 0;
KK(:,1) = 0;
KK(1,1) = 1;
f(1) = 0;

% stream function
SF = KK\f;

% plotting grid
if isempty(zoom)
    x = linspace(min(X(:,1)),max(X(:,1)),np);
    y = linspace(min(X(:,2)),max(X(:,2)),np);
else
    x = linspace(zoom(1),zoom(2),np);
    y = linspace(zoom(3),zoom(4),np);
end
[X_plot,Y_plot] = meshgrid(x,y);

% stream function in the plotting mesh
SF_plot = griddata(X(:,1),X(:,2),SF,X_plot,Y_plot,'cubic');

% plot streamlines
figure(18), clf
plotBoundary(X,boundary)
% plotMesh(X,T)
hold on
% plotSolution(X,T,sqrt(u(:,1).^2+u(:,2).^2),refEl)
contour(X_plot,Y_plot,SF_plot,np,'b')
axis equal

function [Ke fe] = computeElementalMatrix(Xe,refEl,ue)

%Information of the reference element
IPw = refEl.IPweights;
N = refEl.N;
Nxi = refEl.Nxi;
Neta = refEl.Neta;
ngauss = length(IPw);

% iniatilization
Nen = size(refEl.NodesCoord,1);
Ke = zeros(Nen,Nen);
fe = zeros(Nen,1);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

for g = 1:ngauss
    
    % Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    % Jacobian
    Jac = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    
    % Integration weight
    dvolu=IPw(g)*det(Jac);
    
    % x and y derivatives
    invJ = inv(Jac);  %Warning: the inverse is not usually computed!
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    
    % velocity derivatives
    uy = Ny_g * ue(:,1);
    vx = Nx_g * ue(:,2);
    
    %Contribution of the current integration point to the elemental matrix
    Ke = Ke + (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    fe = fe + N_g'*(uy-vx)*dvolu;
end

function fe = computeElementalBoundaryVector(Xf,refEl,uf)

Nf = size(Xf,1);
fe = zeros(Nf,1);

%Information of the reference element
IPw = refEl.IPweights1d; 
N = refEl.N1d; 
Nxi = refEl.N1dxi;

%Number of Gauss points
ngauss = length(IPw);

%Loop in integration points
for g = 1:ngauss
    
  %Shape functions and derivatives at the current integration point 
  N_g = N(g,:);
  Nxi_g = Nxi(g,:);
  
  %Integration weight
  xyDer_g = Nxi_g*Xf;
  xyDerNorm_g = norm(xyDer_g);
  dline=IPw(g)*xyDerNorm_g;
  
  % velocity derivatives
  u_g = N_g * uf(:,1);
  v_g = N_g * uf(:,2);
  
  %Unit normal to the boundary
  t_g = xyDer_g/xyDerNorm_g;
  n_g = [t_g(2) -t_g(1)];
  
  %Contribution of the current integration point to the elemental matrix
  fe = fe + N_g'*(-v_g*n_g(1) + u_g*n_g(2))*dline;
end
