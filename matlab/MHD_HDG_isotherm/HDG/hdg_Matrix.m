function [KK L L0 Lro U U0 Uro P P0 Pro f] =...
    hdg_Matrix(X,T,Xp,Tp,Fcon,flipFace,refElv,refElp,h,velo,dt,BC)

% mesh data
Ne = size(T,1);                     % number of elements
Nv = size(T,2);                     % number of element nodes for the velocity
Np = size(Tp,2);                    % number of element nodes for the pressure
nOfFaces = max(max(Fcon));          % number of faces
nv = size(refElv.NodesCoord1d,1);   % number of face nodes for the velocity
nBC = size(BC.ind,1);

% allocation and initialization
v_unk = nOfFaces*nv*2;
p_unk = Ne;
dim = v_unk + p_unk + 1 + nBC;
allocation = 9*Ne*(2*nv)^2 + 3*2*2*nv*Ne + 2*Ne;
aux_ones_v = ones(1,6*nv);
index = 0;
I = zeros(allocation,1);
J = zeros(allocation,1);
K = zeros(allocation,1);
f = zeros(dim,1);
L = zeros(4*Nv,6*nv,Ne);
L0 = zeros(4*Nv,2*Nv,Ne);
Lro = zeros(4*Nv,1,Ne);
U = zeros(2*Nv,6*nv,Ne);
U0 = zeros(2*Nv,2*Nv,Ne);
Uro = zeros(2*Nv,1,Ne);
P = zeros(Np,6*nv,Ne);
P0 = zeros(Np,2*Nv,Ne);
Pro = zeros(Np,1,Ne);

% loop in elements
for iElem = 1:Ne
    
    ind = (iElem-1)*2*Nv + (1:2*Nv);
    Te = T(iElem,:);
    Xe = X(Te,:);
    Tep = Tp(iElem,:);
    Xep = Xp(Tep,:);
    Fcone = Fcon(iElem,:);
    flipFace_e = flipFace(iElem,:);
    velo_e = velo(ind);
    
    % elemental matrices
    [Le L0e Lroe Ue U0e Uroe Pe P0e Proe Ye ...
         Ffe Hfe Dfe Efe Lfe elemSize] = elementalMatrices(Xe,Xep,refElv,refElp,h,velo_e,dt);
    
    % global assembly indexes for the velocity
    ind_1_v_G = (Fcone(1)-1)*2*nv + (1:2*nv);
    ind_2_v_G = (Fcone(2)-1)*2*nv + (1:2*nv);
    ind_3_v_G = (Fcone(3)-1)*2*nv + (1:2*nv);
    
    % global assembly indexes for the pressure
    ind_p_G = v_unk + iElem;
    
    % assembly index for the pressure lagrangian
    ind_lambda_p = v_unk + Ne + 1;
    
    % assembly index for the pressure lagrangian
    ind_BC = v_unk + Ne + (2:2+nBC-1);
    
    % local assembly indexes
    ind_1_v_L = (1:2*nv);
    ind_2_v_L = 2*nv + (1:2*nv);
    ind_3_v_L = 4*nv + (1:2*nv);
    
    if flipFace_e(1)
        Le(:,ind_1_v_L) = fliplr2(Le(:,ind_1_v_L));
        Pe(:,ind_1_v_L) = fliplr2(Pe(:,ind_1_v_L));
        Ue(:,ind_1_v_L) = fliplr2(Ue(:,ind_1_v_L));
        Ye(ind_1_v_L) = fliplr2(Ye(ind_1_v_L));
        Dfe(ind_1_v_L,:) = flipud2(Dfe(ind_1_v_L,:));
        Ffe(ind_1_v_L,:) = flipud2(Ffe(ind_1_v_L,:));
        Lfe(ind_1_v_L,:) = flipud2(Lfe(ind_1_v_L,:));
    end
    if flipFace_e(2)
        Le(:,ind_2_v_L) = fliplr2(Le(:,ind_2_v_L));
        Pe(:,ind_2_v_L) = fliplr2(Pe(:,ind_2_v_L));
        Ue(:,ind_2_v_L) = fliplr2(Ue(:,ind_2_v_L));
        Ye(ind_2_v_L) = fliplr2(Ye(ind_2_v_L));
        Dfe(ind_2_v_L,:) = flipud2(Dfe(ind_2_v_L,:));
        Ffe(ind_2_v_L,:) = flipud2(Ffe(ind_2_v_L,:));
        Lfe(ind_2_v_L,:) = flipud2(Lfe(ind_2_v_L,:));
    end
    if  flipFace_e(3)
        Le(:,ind_3_v_L) = fliplr2(Le(:,ind_3_v_L));
        Pe(:,ind_3_v_L) = fliplr2(Pe(:,ind_3_v_L));
        Ue(:,ind_3_v_L) = fliplr2(Ue(:,ind_3_v_L));
        Ye(ind_3_v_L) = fliplr2(Ye(ind_3_v_L));
        Dfe(ind_3_v_L,:) = flipud2(Dfe(ind_3_v_L,:));
        Ffe(ind_3_v_L,:) = flipud2(Ffe(ind_3_v_L,:));
        Lfe(ind_3_v_L,:) = flipud2(Lfe(ind_3_v_L,:));
    end
    
    % faces matrices
    KKe_v = -Lfe*Le + Ffe*Pe + Dfe*Ue + Hfe-Efe;
    KKe_p = - Lfe*Lroe + Ffe*Proe + Dfe*Uroe;
    ffe = (Lfe*L0e - Ffe*P0e - Dfe*U0e)*velo_e;
    
    % assembly for velocity: momentum
    indv_m = [ind_1_v_G ind_2_v_G ind_3_v_G];
    indv_transp_m = transpose(indv_m);
    aux_row = indv_transp_m(:,aux_ones_v);
    aux_col = indv_m(aux_ones_v,:);
    index = index(end) + (1:9*4*nv^2);
    I(index) = aux_row(:);
    J(index) = aux_col(:);
    K(index) = KKe_v;
    f(indv_m) = f(indv_m) + ffe;
    
    % assembly for the velocity: continuity
    indv_c = ind_p_G;
    aux_col = indv_transp_m;    
    aux_row = indv_c(:,aux_ones_v);
    index = index(end) + (1:3*2*nv);
    I(index) = aux_row(:);
    J(index) = aux_col(:);
    K(index) = Ye;
    
    % assembly for pressure
    indp = ind_p_G;
    aux_col = indp(:,aux_ones_v);
    aux_row = indv_transp_m;
    index = index(end) + (1:3*2*nv);
    I(index) = aux_row(:);
    J(index) = aux_col(:);
    K(index) = KKe_p;
    
    % assembly the Lagrange multipliers for pressure
    aux_col = [ind_lambda_p ind_p_G];
    aux_row = [ind_p_G ind_lambda_p];
    index = index(end) + (1:2);
    I(index) = aux_row(:);
    J(index) = aux_col(:);
    K(index) = [elemSize elemSize];
    
    % store mapping
    L(:,:,iElem) = Le;
    L0(:,:,iElem) = L0e;
    Lro(:,:,iElem) = Lroe;
    U(:,:,iElem) = Ue;
    U0(:,:,iElem) = U0e;
    Uro(:,:,iElem) = Uroe;
    P(:,:,iElem) = Pe;
    P0(:,:,iElem) = P0e;
    Pro(:,:,iElem) = Proe;
end

% check allocation
if size(I,1)>allocation
    error('size overpassed')
end

% assembly BC
aux_row = [BC.ind; ind_BC'];
aux_col = [ind_BC'; BC.ind];
index = index(end) + (1:2*nBC);
I(index) = aux_row(:);
J(index) = aux_col(:);
K(index) = ones(2*nBC,1);
f(ind_BC) = f(ind_BC) + BC.val;

% create sparse matrix
KK = sparse(I(I~=0),J(I~=0),K(I~=0),dim,dim);

%% Elemental matrices
function [ LL LL0 LLro UU UU0 UUro PP PP0 PPro Y ...
           Ff Hf Df Ef Lf elemSize] = elementalMatrices(Xe,Xep,refElv,refElp,h,velo_e,dt)

% mesh data
Nv = size(refElv.NodesCoord,1);
nv = size(refElv.NodesCoord1d,1);
Np = size(refElp.NodesCoord,1);
np = size(refElp.NodesCoord1d,1);
faceNodesv = refElv.faceNodes;
faceNodesp = refElp.faceNodes;
nOfFaces = size(faceNodesv,1);

% initialize all the matrices
C = zeros(4*Nv,6*nv);
D = zeros(2*Nv,2*Nv);
Df = zeros(6*nv,2*Nv);
E = zeros(2*Nv,6*nv);
Ef = zeros(6*nv,6*nv);
F = zeros(2*Nv,Np);
Ff = zeros(6*nv,Np);
H = zeros(2*Nv,6*nv);
Hf = zeros(6*nv,6*nv);
L = zeros(2*Nv,4*Nv);
O = zeros(Np,6*nv);
Y = zeros(1,6*nv);
Mass = zeros(Nv);
Bx = Mass;
By = Mass;
Gx = zeros(Nv,Np);
Gy = Gx;
Rx = zeros(Np,Nv);
Ry = Rx;
Cvxx = Mass;
Cvxy = Mass;
Cvyx = Mass;
Cvyy = Mass;
intN = zeros(Np,1);
elemSize = 0;
elemPer = 0;
 
% Information of the reference element for the velocity
IPw = refElv.IPweights;                 % use the velocity gauss points to integrate
Niv = refElv.N;
Nxiv = refElv.Nxi;
Netav = refElv.Neta;
IPw_fv = refElv.IPweights1d;
N1dv = refElv.N1d;
Nx1dv = refElv.N1dxi;

% Information of the reference element for the pressure
Nip = refElp.N;
Nxip = refElp.Nxi;
Netap = refElp.Neta;
N1dp = refElp.N1d;

% Number of Gauss points in the interior 
ngauss = length(IPw);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);
xep = Xep(:,1); yep = Xep(:,2);

% reshape velocity 
velo_e = reshape(velo_e,2,0.5*numel(velo_e))';

%% VOLUME COMPUTATIONS
for g = 1:ngauss
    
    % Velocity shape functions and derivatives at the current integration point 
    Niv_g = Niv(g,:);
    Nxiv_g = Nxiv(g,:);
    Netav_g = Netav(g,:);
    
    % Pressure shape functions at the current integration point
    Nip_g = Nip(g,:);
    Nxip_g = Nxip(g,:);
    Netap_g = Netap(g,:);
    
    % Jacobian
    J = [Nxiv_g*xe	  Nxiv_g*ye
        Netav_g*xe  Netav_g*ye];
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end
    
    % Jacobian for the pressure element
    Jp = [Nxip_g*xep	  Nxip_g*yep
        Netap_g*xep  Netap_g*yep];
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end
    
    % Integration weight
    dvolu=IPw(g)*det(J);
    dvolup = IPw(g)*det(Jp);
    
    % x and y derivatives
    invJ = inv(J);  
    Nx_g = invJ(1,1)*Nxiv_g + invJ(1,2)*Netav_g;
    Ny_g = invJ(2,1)*Nxiv_g + invJ(2,2)*Netav_g;
    
    % x and y derivatives for the pressure
    invJp = inv(Jp);  
    Nxp_g = invJp(1,1)*Nxip_g + invJp(1,2)*Netap_g;
    Nyp_g = invJp(2,1)*Nxip_g + invJp(2,2)*Netap_g;
    
    % interpolate velocity in the current integration point
    velox_g = Niv_g*velo_e(:,1);
    veloy_g = Niv_g*velo_e(:,2);
    
    % Contribution of the current integration point to the elemental matrix
    Mass = Mass + Niv_g'*Niv_g*dvolu;
    Bx = Bx + Nx_g'*Niv_g*dvolu;
    By = By + Ny_g'*Niv_g*dvolu;
    Gx = Gx + Nx_g'*Nip_g*dvolu;
    Gy = Gy + Ny_g'*Nip_g*dvolu;
    Rx = Rx + Nxp_g'*Niv_g*dvolup;
    Ry = Ry + Nyp_g'*Niv_g*dvolup;
    Cvxx = Cvxx + velox_g*Nx_g'*Niv_g*dvolu;
    Cvxy = Cvxy + velox_g*Ny_g'*Niv_g*dvolu;
    Cvyx = Cvyx + veloy_g*Nx_g'*Niv_g*dvolu;
    Cvyy = Cvyy + veloy_g*Ny_g'*Niv_g*dvolu;
    intN = intN + Nip_g'*dvolu;
    elemSize = elemSize + dvolu;
end

% expand the matrices
A = expandMatrixA(Mass,4);
M = expandMatrixA(Mass,2);
B = expandMatrixB(Bx,By);
Cv = expandMatrixCv(Cvxx,Cvxy,Cvyx,Cvyy);
G = expandMatrixF(Gx,Gy);
R = transpose(expandMatrixF(Rx',Ry'));

% tilded pressure shape functions
b = 1/elemSize*intN;

%% FACES COMPUTATIONS:
ngauss_f = length(IPw_fv);
for iface = 1:nOfFaces
    
    % face nodes
    nodesv = faceNodesv(iface,:);
    nodesp = faceNodesp(iface,:);
    
    % indices for local assembly
    ind_face_2 = (iface-1)*2*nv + (1:2*nv);
    ind2 = reshape(bsxfun(@plus,(nodesv-1)*2,(1:2)'),2*nv,1); % assembly face to elem for velocity
    ind4 = reshape(bsxfun(@plus,(nodesv-1)*4,(1:4)'),4*nv,1); % assembly face to elem for velocity gradient
    
    xf = xe(nodesv);
    yf = ye(nodesv);
    velo_f = velo_e(nodesv,:);
    %     tau_f = 1/h^2;
    %     if iface~=1
    %         tau_f = 0;
    %     end
%     tau_f = 1/h;
    tau_f = 1e0;
    
    % initialize local matrices
    Cnx = zeros(nv);
    Cny = Cnx;
    Fx = zeros(nv,np);
    Fy = Fx;
    Hxx = Cnx;
    Hxy = Cnx;
    Hyx = Cnx;
    Hyy = Cnx;
    Massf = Cnx;
    Ox = Fx';
    Oy = Ox;
    Yx = zeros(1,nv);
    Yy = Yx;
    
    %  LOOP IN GAUSS POINTS
    for g = 1:ngauss_f 
        
        % Velocity shape functions and derivatives at the current integration point
        Nfv_g = N1dv(g,:);
        Nfxiv_g = Nx1dv(g,:);
        
        % interpolate velocity in the current integration point
        velofx_g = Nfv_g*velo_f(:,1);
        velofy_g = Nfv_g*velo_f(:,2);
        
        % Pressure shape functions at the current integration point
        Nfp_g = N1dp(g,:);
        Nfp_tilde_g = Nfp_g;% - bf';   %%% bf to be used just with P0 elements for the pressure?? 
        
        %Integration weight
        xyDer_g = Nfxiv_g*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline = IPw_fv(g)*xyDerNorm_g;
        
        %Unit normal to the boundary
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) -t_g(1)];
        
        %Contribution of the current integration point to the elemental matrix
        Cnx = Cnx + Nfv_g'*Nfv_g*n_g(1)*dline;
        Cny = Cny + Nfv_g'*Nfv_g*n_g(2)*dline;
        Fx = Fx + Nfv_g'*Nfp_g*n_g(1)*dline;
        Fy = Fy + Nfv_g'*Nfp_g*n_g(2)*dline;
        Hxx = Hxx + Nfv_g'*Nfv_g*n_g(1)*velofx_g*dline;
        Hxy = Hxy + Nfv_g'*Nfv_g*n_g(2)*velofx_g*dline;
        Hyx = Hyx + Nfv_g'*Nfv_g*n_g(1)*velofy_g*dline;
        Hyy = Hyy + Nfv_g'*Nfv_g*n_g(2)*velofy_g*dline;
        Massf = Massf + tau_f *(Nfv_g'*Nfv_g)*dline;        
        Ox = Ox + Nfp_tilde_g'*Nfv_g*n_g(1)*dline;
        Oy = Oy + Nfp_tilde_g'*Nfv_g*n_g(2)*dline;       
        Yx = Yx + Nfv_g*n_g(1)*dline;
        Yy = Yy + Nfv_g*n_g(2)*dline;
        elemPer = elemPer + dline;
    end
    
    % expand the matrices
    C_loc = expandMatrixB(Cnx,Cny);
    F_loc = expandMatrixF(Fx,Fy);
    H_loc = expandMatrixCv(Hxx,Hxy,Hyx,Hyy);
    D_loc = expandMatrixA(Massf,2);
    E_loc = expandMatrixA(Massf,2);
    O_loc = transpose(expandMatrixF(Ox',Oy'));
    Y_loc = transpose(expandMatrixF(Yx',Yy'));
     
    % elemental assembly
    C(ind4,ind_face_2) = C(ind4,ind_face_2) + C_loc;
    L(ind2,ind4) = L(ind2,ind4) + C_loc';
    F(ind2,nodesp) = F(ind2,nodesp) + F_loc;
    Ff(ind_face_2,nodesp) = Ff(ind_face_2,nodesp) + F_loc;
    H(ind2,ind_face_2) = H(ind2,ind_face_2) + H_loc;
    Hf(ind_face_2,ind_face_2) = Hf(ind_face_2,ind_face_2) + H_loc;
    D(ind2,ind2) = D(ind2,ind2) + D_loc;
    Df(ind_face_2,ind2) = Df(ind_face_2,ind2) + D_loc;
    E(ind2,ind_face_2) = E(ind2,ind_face_2) + E_loc;
    Ef(ind_face_2,ind_face_2) = Ef(ind_face_2,ind_face_2) + E_loc;
    O(nodesp,ind_face_2) = O(nodesp,ind_face_2) + O_loc;
    Y(ind_face_2) = Y(ind_face_2) + Y_loc;

end
Lf = C';
W = b;

%% matrix multiplication to get the mapping

% first set
invA = inv(A);
Mu = M/dt - Cv + D - (B'- L)*invA*B; 
Mp = F - G;
Mu_tilde = E - H - (B' - L)*invA*C;

% second set
invMu = inv(Mu); 
Gp = R*invMu*Mp;
Gu_tilde = R*invMu*Mu_tilde-O;
Gu = R*invMu*M/dt;

% third set
K = [Gp W; W' 0];
invK = inv(K);

% mapping for the pressure
PP = invK(1:end-1,1:end-1)*Gu_tilde;
PP0 = invK(1:end-1,1:end-1)*Gu;
PPro = invK(1:end-1,end);

% mapping for the velocity
UU = invMu*(Mu_tilde-Mp*PP);
UU0 = invMu*(M/dt-Mp*PP0);
UUro = -invMu*Mp*PPro;

% mapping for the velocity gradient
LL = invA*(C-B*UU);
LL0 = -invA*B*UU0;
LLro = -invA*B*UUro;

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

function res = expandMatrixCv(Cxx,Cxy,Cyx,Cyy)
% expand matrix Cv
%   [ Cxx Cxy
%     Cyx Cyy ]
res = zeros([size(Cxx) 2 2]);
res(:,:,[1 3 2 4]) = cat(3,Cxx,Cxy,Cyx,Cyy);
res = permute(res, [3 1 4 2]);
res = reshape(res, 2*size(Cxx,1),2*size(Cxx,2));

function res = expandMatrixF(Fx,Fy)
% expand matrix
%   [ Fx
%     Fy ]
res = zeros([size(Fx) 2 1]);
res(:,:,[1 2]) = cat(3,Fx,Fy);
res = permute(res, [3 1 4 2]);
res = reshape(res, 2*size(Fx,1),size(Fx,2));

function res = flipud2(A)
% [ 1 1 1      [ 5 5 5 
%   2 2 2        6 6 6
%   3 3 3  ===>  3 3 3
%   4 4 4        4 4 4
%   5 5 5        1 1 1
%   6 6 6 ]      2 2 2 ]
res = zeros(size(A));
res(end-1:-2:1,:) = A(1:2:end-1,:);
res(end:-2:2,:) = A(2:2:end,:);

function res = fliplr2(A)
% [ 1 1 1      [ 5 5 5 
%   2 2 2        6 6 6
%   3 3 3  ===>  3 3 3
%   4 4 4        4 4 4
%   5 5 5        1 1 1
%   6 6 6 ]      2 2 2 ]
res = zeros(size(A));
res(:,end-1:-2:1) = A(:,1:2:end-1);
res(:,end:-2:2) = A(:,2:2:end);




