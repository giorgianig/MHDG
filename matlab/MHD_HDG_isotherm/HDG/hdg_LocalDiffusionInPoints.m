function [Lf_locDiffPts,P_locDiffPts,Q_locDiffPts] =  ...
    hdg_LocalDiffusionInPoints(X,T,flipFace,refEl,Lf,P,Q)

% mesh data
global neq Mesh axisym          % number of equations (2 for isothermal)
nf  = size(refEl.faceNodes,1);    % number of faces in one element
Ne = size(T,1);                   % number of elements
nv = size(refEl.NodesCoord1d,1);    % number of face nodes

% local assembly indexes
ind_v_L = zeros(nf, neq*nv);
for iface = 1:nf
    ind_v_L(iface,:) = neq*nv*(iface-1)+ (1:neq*nv);
end

% compute permutation for flipping faces
perm = setperm(neq*nv,neq);

% plotMesh(X*Mesh.lscale,T)
% hold on
Lf_locDiffPts = zeros(size(Lf));
P_locDiffPts  = zeros(size(P));
Q_locDiffPts  = zeros(size(Q));

% loop in elements
for iElem = 1:Ne
    
    Te = T(iElem,:);
    Xe = X(Te,:);    
    
    [~,~,xy_corn] = computeLocalDiffusionInPoints(Xe*Mesh.lscale);
    
    bc = [(Xe(1:3,1)*Mesh.lscale)';(Xe(1:3,2)*Mesh.lscale)';ones(1,3)]...
        \[xy_corn(1,:);xy_corn(2,:);ones(1,size(xy_corn,2))];
    if ~any(all(bc>-1e-6,1)), continue, end
    
%     plot(X(T(iElem,:),1)*Mesh.lscale,X(T(iElem,:),2)*Mesh.lscale,'ro')
    
    flipFace_e = flipFace(iElem,:);
    
    % elemental matrices
    [Lfe,Pe,Qe] = elementalMatrices(Xe,refEl);
    
    for iface = 1:nf
        
        if flipFace_e(iface)
            Lfe(ind_v_L(iface,:),:) = Lfe(ind_v_L(iface,perm),:);
        end
    end
    
    % store matrices
    Lf_locDiffPts(:,:,iElem) =  Lfe;
    P_locDiffPts(:,:,iElem)  =  Pe;
    Q_locDiffPts(:,:,iElem)  =  Qe;
end


%% Elemental matrices
function [Lf,P,Q] = elementalMatrices(Xe,refEl)

% eps_nodal: nodal value of the artificial diffusion


global neq Mesh axisym

% mesh data
nd = size(Xe,2);
Nv = size(refEl.NodesCoord,1);
nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);

% initialize all the matrices
Lf_transp = zeros(neq*nd*Nv,nd*nf*nv);
Q = zeros(2*Nv,4*Nv);
P = Q;

% Information of the reference element for the velocity
IPw = refEl.IPweights;                 % use the velocity gauss points to integrate
Niv = refEl.N;
Nxiv = refEl.Nxi;
Netav = refEl.Neta;
IPw_fv = refEl.IPweights1d;
N1dv = refEl.N1d;
Nx1dv = refEl.N1dxi;

% Number of Gauss points in the interior
ngauss = length(IPw);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);
Xg = Niv*Xe;

[eps_point_n,eps_point_u] = computeLocalDiffusionInPoints(Xg*Mesh.lscale);


%% VOLUME COMPUTATIONS
for g = 1:ngauss
    
    % Velocity shape functions and derivatives at the current integration point
    Niv_g = Niv(g,:);
    Nxiv_g = Nxiv(g,:);
    Netav_g = Netav(g,:);
    
    % gauss point position
    %     xg = Niv_g*xe;
    %     yg = Niv_g*ye;
    
    % Jacobian
    J = [Nxiv_g*xe	  Nxiv_g*ye
        Netav_g*xe  Netav_g*ye];
    
    % x and y derivatives
    invJ = inv(J);
    Nx_g = invJ(1,1)*Nxiv_g + invJ(1,2)*Netav_g;
    Ny_g = invJ(2,1)*Nxiv_g + invJ(2,2)*Netav_g;
    
    % Integration weight
    dvolu=IPw(g)*det(J);
    if axisym
        xg = Niv_g*xe;
        dvolu=dvolu*xg;
    end
    % Contribution of the current integration point to the elemental matrix
    Px = Nx_g'*Niv_g*dvolu;
    Py = Ny_g'*Niv_g*dvolu;
    P = P + transpose(expandMatrixB(Px',Py',eps_point_n(g),eps_point_u(g)));
end

%% FACES COMPUTATIONS:
ngauss_f = length(IPw_fv);
for iface = 1:nf
    
    % face nodes
    nodesv = faceNodesv(iface,:);
    
    % indices for local assembly
    ind_face_2 = (iface-1)*neq*nv + (1:neq*nv);
    ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
    ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for velocity gradient
    
    xf = xe(nodesv);
    yf = ye(nodesv);
    xyfg = N1dv*[xf yf];
    [eps_face_n_gauss,eps_face_u_gauss] = computeLocalDiffusionInPoints(xyfg*Mesh.lscale);
    
    % initialize local matrices
    C_loc = zeros(4*nv,2*nv);
    
    %  LOOP IN GAUSS POINTS
    for g = 1:ngauss_f
        
        % Velocity shape functions and derivatives at the current integration point
        Nfv_g = N1dv(g,:);
        Nfxiv_g = Nx1dv(g,:);
        
        % Integration weight
        xyDer_g = Nfxiv_g*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline = IPw_fv(g)*xyDerNorm_g;
        if axisym
            x = Nfv_g*xf;
            dline = dline*x;
        end       
        % Unit normal to the boundary
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) -t_g(1)];
        
        % Contribution of the current integration point to the elemental matrix
        Cnx = Nfv_g'*Nfv_g*n_g(1)*dline;
        Cny = Nfv_g'*Nfv_g*n_g(2)*dline;
        C_loc = C_loc + expandMatrixB(Cnx,Cny,eps_face_n_gauss(g),eps_face_u_gauss(g));
    end
    
    % elemental assembly
    Lf_transp(ind4,ind_face_2) = Lf_transp(ind4,ind_face_2) + C_loc;
    Q(ind2,ind4) = Q(ind2,ind4) + C_loc';
end
Lf = Lf_transp';


%% additional routines
function res = expandMatrixB(Bx,By,eps_n,eps_u)
% expand matrix B
%   [ Bx  0
%     By  0
%     0  Bx
%     0  By ]
res = zeros([size(Bx) 4 2]);
res(:,:,[1 2 7 8]) = cat(3,eps_n*Bx,eps_n*By,eps_u*Bx,eps_u*By);
res = permute(res, [3 1 4 2]);
res = reshape(res, 4*size(Bx,1),2*size(Bx,2));


