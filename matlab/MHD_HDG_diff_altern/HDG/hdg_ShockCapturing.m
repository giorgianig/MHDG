function [Lf,P,Q] = hdg_ShockCapturing(X,T,u,flipFace,refEl,shock_st)

global Mesh shockcapt
% swt_sc: switch of the shock capturing parameter:
% 0 - elementwise
% 1 - linear interpolation in each element

% mesh data
global neq               % number of equations (2 for isothermal)
nf  = size(refEl.faceNodes,1);    % number of faces in one element
Ne = size(T,1);                   % number of elements
Nv = size(refEl.NodesCoord,1);  % number of nodes
nv = size(refEl.NodesCoord1d,1);    % number of face nodes
nd = size(X,2);

% shock capturing parameter in each element
eps_elem = findCoeffShockCaptur(u,X,T,refEl,-5.4);

% local assembly indexes
ind_v_L = zeros(nf, neq*nv);
for iface = 1:nf
   ind_v_L(iface,:) = neq*nv*(iface-1)+ (1:neq*nv);
end

% compute permutation for flipping faces
perm = setperm(neq*nv,neq);
eps_plot = zeros(size(X,1),1);

Lf = zeros(neq*nf*nv,neq*nd*Nv,Ne);
P = zeros(neq*Nv,nd*neq*Nv,Ne);
Q = zeros(neq*Nv,nd*neq*Nv,Ne);

% loop in elements
for iElem = 1:Ne
    
    Te = T(iElem,:);
    Xe = X(Te,:);
    % shock capturing parameter in each node
    switch shockcapt
        case 1
            if ~eps_elem(iElem), continue, end
            
            % constant value in each element
            eps_nodal_n = eps_elem(iElem)*ones(Nv,1);
            eps_nodal_u = eps_elem(iElem)*zeros(Nv,1);
            eps_plot(Te) = eps_nodal_n;
        case 2
            % linear interpolation
            eps_nodal_n = eps_elem(iElem)*zeros(Nv,1);
            eps_nodal_u = eps_elem(iElem)*zeros(Nv,1);            
            els = Mesh.N(Te(1:3),:); 
            aux_eps_nodal_n = zeros(3,1);
            for i = 1:size(els,1)
                ind = els(i,els(i,:)~=0);
                aux_eps_nodal_n(i) = max(eps_elem(ind));
            end
            if all(aux_eps_nodal_n==0), continue,end
            eps_nodal_n = shock_st.N*aux_eps_nodal_n;
            eps_plot(Te) = eps_nodal_n;
    end
%     disp('Artificial diffusion added')    
    flipFace_e = flipFace(iElem,:);
   
    % elemental matrices
    [Lfe,Pe,Qe] = ...
        elementalMatrices(Xe,refEl,eps_nodal_n,eps_nodal_u);
 
 for iface = 1:nf
     
     if flipFace_e(iface)
         Lfe(ind_v_L(iface,:),:) = Lfe(ind_v_L(iface,perm),:);
     end
 end
     
    % store matrices
    Lf(:,:,iElem) = Lfe;
    P(:,:,iElem)  = Pe;
    Q(:,:,iElem)  = Qe;
end
 

% figure, plotSolution(X,T,eps_plot,refEl), hold on, plotMesh(X,T);

%% Elemental matrices
function [Lf,P,Q] = elementalMatrices(Xe,refEl,eps_nodal_n,eps_nodal_u)

% eps_nodal: nodal value of the artificial diffusion


global neq axisym

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

eps_sc_n = Niv*eps_nodal_n; 
eps_sc_u = Niv*eps_nodal_u;
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
    P = P + transpose(expandMatrixB(Px',Py',eps_sc_n(g),eps_sc_u(g)));
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
    eps_face_n = eps_nodal_n(nodesv);
    eps_face_u = eps_nodal_u(nodesv);
%     xyfg = N1dv*[xf yf];
    eps_face_n_gauss = N1dv*eps_face_n;
    eps_face_u_gauss = N1dv*eps_face_u;
    
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


