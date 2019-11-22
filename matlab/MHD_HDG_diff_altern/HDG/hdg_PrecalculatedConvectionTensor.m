function [ Cv_tens H_tens ff_dir] =...
    hdg_PrecalculatedConvectionTensor(X,T,F_dir,refElv,Re)

% Mesh data
Ne = size(T,1);                     % number of elements
Nv = size(T,2);                     % number of element nodes for the velocity
Nf_dir = sum(sum(F_dir));          % number of Dirichlet faces
nv = size(refElv.NodesCoord1d,1);   % number of face nodes for the velocity
nd = 2; %xi,eta
nD = 2; %x,y

% Information of the reference element
IPw = refElv.IPweights;                % ng x 1
N = refElv.N;                          % ng x nOfElementNodes
Nxi = refElv.Nxi;                      % ng x nOfElementNodes
Neta = refElv.Neta;                    % ng x nOfElementNodes
Nf = refElv.N1d;                      % ng_f x nOfFaceNodes
Nxi_f = refElv.N1dxi;                  % ng_f x nOfFaceNodes
IPw_f = refElv.IPweights1d;            % ng_f x 1
faceNodes = refElv.faceNodes;          % nOfFaceNodes x nOfElementFaces
ngauss = length(IPw);
ngauss_f = length(IPw_f);
nf = size(faceNodes,1);

% Some reshapes and permutes
N = permute(N,[2 1]);
N = reshape(N,[Nv 1 ngauss 1]);        % nOfElementNodes x 1 x ngauss
Nxi = permute(Nxi,[2 1]);
Nxi = reshape(Nxi,[Nv 1 ngauss 1]);    % nOfElementNodes x 1 x ngauss
Neta = permute(Neta,[2 1]);
Neta = reshape(Neta,[Nv 1 ngauss 1]);  % nOfElementNodes x 1 x ngauss
IPw = reshape(IPw,[1 1 1 1 ngauss]);
Nf = permute(Nf,[2 1]);
IPw_f = reshape(IPw_f,[1 1 1 ngauss_f]);
Nxieta = [Nxi Neta];                   % np x nd x ngauss

% A_i_j_k_id_g = IPw_g * N_i_g * N_j_id_g * N_k_g
N_i_id_g = reshape(Nxieta,[Nv 1 1 nd ngauss]);
N_j_g = reshape(N,[1 Nv 1 1 ngauss]);
N_k_g = reshape(N,[1 1 Nv 1 ngauss]);
N_k_g = bsxfun(@times,N_k_g,IPw); % apply weights
N_i_j_k_id_g = bsxfun(@times,N_j_g,N_k_g);
N_i_j_k_id_g = bsxfun(@times,N_i_j_k_id_g,N_i_id_g);
A_tens = reshape(N_i_j_k_id_g,[Nv*Nv*Nv nd*ngauss]);

%%%% F_i_j_k_g = IPw_g * N_i_g * N_j_g * N_k_g
Nf_i = reshape(Nf,[nv 1 1 ngauss_f]);
Nf_j = reshape(Nf,[1 nv 1 ngauss_f]);
Nf_k = reshape(Nf,[1 1 nv ngauss_f]);
Nf_g = bsxfun(@times,Nf_i,IPw_f);
Nf_i_j_k_g = bsxfun(@times,Nf_g,Nf_j);
Nf_i_j_k_g = bsxfun(@times,Nf_i_j_k_g,Nf_k);
F_tens = reshape(Nf_i_j_k_g,[nv*nv*nv ngauss_f]);

% Tensor for elemental Jacobian
T_resh = transpose(T);                  % nOfElementNodes x nOfElements
T_resh = reshape(T_resh,1,Ne*Nv);       % 1 x nOfElementNodes * nOfElements
X_t = transpose(X);
X_t = reshape(X_t(:,T_resh),2,Nv,1,Ne); % nD x nOfElementNodes x 1 x e
Nxieta_r = reshape(Nxieta,[1 Nv nd*ngauss]);
J = bsxfun(@times,X_t,Nxieta_r);
J = reshape(sum(J,2),[nD nd ngauss Ne]);
invJ11 =  J(2,2,:,:);                   %Not 1/detJ (optimized)!
invJ12 = -J(1,2,:,:);
invJ21 = -J(2,1,:,:);
invJ22 =  J(1,1,:,:);
B_tens = zeros(nd,nD,ngauss,Ne);
B_tens(1,1,:,:) = invJ11;
B_tens(1,2,:,:) = invJ12;
B_tens(2,1,:,:) = invJ21;
B_tens(2,2,:,:) = invJ22;
B_tens = permute(B_tens,[1 3 2 4]);
B_tens = reshape(B_tens,[nd*ngauss nD*Ne]);

% Tensor for faces Jacobian
X_t = permute(X_t,[2 1 3 4]);                                 % nOfElementNodes x nD x 1 x e
X_f(:,:,1,:) = X_t(faceNodes(1,:),:,:,:);                     % nOfFaceNodes x nD x f x e
X_f(:,:,2,:) = X_t(faceNodes(2,:),:,:,:);
X_f(:,:,3,:) = X_t(faceNodes(3,:),:,:,:);
xyDer_g = bsxfun(@times,reshape(X_f,[1 nv nD nf Ne]),Nxi_f);
xyDer_g = reshape(sum(xyDer_g,2),[ngauss_f nD nf Ne]);             % ng x nD x f x e
normal_J = [xyDer_g(:,2,:,:) -xyDer_g(:,1,:,:)];                   % ng x nD x f x e
Bf_tens = reshape(normal_J,[ngauss_f,nD*nf*Ne]);

% Tensor for elemental convection matrix
Cv_tens = A_tens*B_tens;
Cv_tens = reshape(Cv_tens,[Nv Nv Nv 1 nD Ne]);

% Tensor for faces convection matrix
H_tens = F_tens*Bf_tens;
H_tens = reshape(H_tens,[nv nv nv 1 nD nf Ne]);

% Apply Dirichlet boundary conditions
X_f = permute(X_f,[1 2 4 3]);  % nOfFaceNodes x nD x e x f
normal_J = permute(normal_J,[1 2 4 3]); % nOfFaceNodes x nD x e x f
X_f_dir = reshape(X_f(:,:,F_dir),nv,2*Nf_dir); % nOfFaceNodes x nD*Nf_dir
Xg_f_dir =  refElv.N1d*X_f_dir;   % ng x nD*Nf_dir
Xg_f_dir = reshape(Xg_f_dir,ngauss_f,nD,Nf_dir); % ng x nD x Nf_dir
Xg_f_dir = reshape(permute(Xg_f_dir, [1 3 2]),ngauss_f*Nf_dir,nD);
velo_dir = reshape(setDirichletBoundaryConditions(Xg_f_dir,Re),ngauss_f,Nf_dir,nD);
velo_dir = permute(velo_dir,[1 3 2]);  % ng x nD x Nf_dir
Nf_g = reshape(Nf_g,nv,ngauss_f);
velo_dir_n = bsxfun(@times,velo_dir,normal_J(:,:,F_dir)); % ng x nD x Nf_dir
velo_dir_n = reshape(velo_dir_n,ngauss_f,nD,1,Nf_dir);  % ng x nD x 1 x Nf_dir
velo_dir = reshape(velo_dir,ngauss_f,1,nD,Nf_dir); % ng x 1 x nD x Nf_dir
velo_dir_2_n = bsxfun(@times,velo_dir_n,velo_dir); % ng x nD x nD x Nf_dir
velo_dir_2_n = reshape(sum(velo_dir_2_n, 2),ngauss_f,nD*Nf_dir);
fH_tens = reshape(Nf_g*velo_dir_2_n,nv,nD,Nf_dir);
ff_dir = reshape(permute(fH_tens,[2 1 3]),[2*nv Nf_dir]);   % 2*nv x Nf
H_tens(:,:,:,:,:,F_dir') = 0;





