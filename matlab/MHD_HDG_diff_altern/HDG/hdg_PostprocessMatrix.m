function [K Bt int_N M] = hdg_PostprocessMatrix(X,T,refElv)

% fine mesh informations
Ne = size(T,1);
Nv = size(T,2);

% allocating
K = zeros(2*Nv,2*Nv,Ne);
M = zeros(Nv,Nv,Ne);
Bt = zeros(2*Nv,4*Nv,Ne);
int_N = zeros(2*Nv,2,Ne);

% Loop in elements
for iElem = 1:Ne
    Te = T(iElem,:);
    Xe = X(Te,:);
    [Ke,Bte,int_Ne,Me] = ElementalMatrices(Xe,refElv);
    
    % store matrix
    K(:,:,iElem) = Ke;
    M(:,:,iElem) = Me;
    Bt(:,:,iElem) = Bte;
    int_N(:,:,iElem) = int_Ne;
        
end

function [K,Bt,int_N,Me] = ElementalMatrices(Xe,refElv)

Nv = size(refElv.NodesCoord,1);
K_loc = zeros(Nv);
Me = K_loc;
Bx = K_loc;
By = K_loc;
int_N_loc = zeros(Nv,1);

%Information of the reference element
IPw = refElv.IPweights;
N = refElv.N;
Nxi = refElv.Nxi;
Neta = refElv.Neta;

%Number of Gauss points in the interior
ngauss = length(IPw);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% VOLUME COMPUTATIONS: LOOP IN GAUSS POINTS
for g = 1:ngauss
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);

    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];

    %Integration weight
    dvolu=IPw(g)*det(J);

    %x and y derivatives
    invJ = inv(J);
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;

    %Contribution of the current integration point to the elemental matrix
    K_loc = K_loc + (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    Me = Me + N_g'*N_g*dvolu;
    Bx = Bx + Nx_g'*N_g*dvolu;
    By = By + Ny_g'*N_g*dvolu;
    int_N_loc = int_N_loc + N_g'*dvolu;
end

K = expandMatrixA(K_loc,2);
Bt = transpose(expandMatrixB(Bx',By'));
int_N = expandMatrixA(int_N_loc,2);

%% additional routines

function res = expandMatrixA(A,n)
% expand matrix A and M
%  [ A 0 0 0 
%    0 A 0 0
%    0 0 A 0
%    0 0 0 A ]
% dimension n
res = zeros([size(A) n n]);
res(:,:,1:n+1:n^2) = repmat(A, [1 1 n]);
res = permute(res, [3 1 4 2]);
res = reshape(res, n*size(A));


function res = expandMatrixB(Bx,By)
% expand matrix B
%   [ Bx  0
%     By  0
%     0  Bx
%     0  By ]
res = zeros([size(Bx) 4 2]);
res(:,:,[1 2 7 8]) = cat(3,Bx,By,Bx,By);
res = permute(res, [3 1 4 2]);
res = reshape(res, 4*size(Bx,1),2*size(Bx,2));