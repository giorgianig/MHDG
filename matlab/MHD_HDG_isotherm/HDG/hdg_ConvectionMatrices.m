function [ C,H,Hdir,Hf] =...
    hdg_ConvectionMatrices(X,T,F,flipFace,refEl,sol,sol_hat,F_dir)

% mesh data
global neq
Ne = size(T,1);                     % number of elements
Nv = size(T,2);                     % number of element nodes for the velocity
nv = size(refEl.NodesCoord1d,1);   % number of face nodes for the velocity
nf  = size(refEl.faceNodes,1);    % number of faces in one element

% allocation and initialization
C = zeros(neq*Nv,neq*Nv,Ne);
H = zeros(neq*Nv,neq*nf*nv,Ne);
Hdir = zeros(neq*Nv,Ne);
Hf = zeros(neq*nf*nv,neq*nf*nv,Ne);

% local assembly indexes
ind_v_L = zeros(nf, neq*nv);
for iface = 1:nf
    ind_v_L(iface,:) = neq*nv*(iface-1)+ (1:neq*nv);
end
% compute permutation for flipping faces
perm = setperm(neq*nv,neq);

% loop in elements
for iElem = 1:Ne
    
    Fe = F(iElem,:);
    ind = (iElem-1)*neq*Nv + (1:neq*Nv);
    ind_uhat =  bsxfun(@plus,(Fe-1)*neq*nv,(1:neq*nv)');
    sol_hat_e = sol_hat(ind_uhat);
    Te = T(iElem,:);
    Xe = X(Te,:);
    sol_e = sol(ind);
    flipFace_e = flipFace(iElem,:);
    aux_dir = F_dir(iElem,:);
    
    % elemental matrices
    [Ce,He,Hdire,Hfe] = elementalMatrices(Xe,refEl,sol_e,sol_hat_e,aux_dir,flipFace_e,iElem,Fe);
        
    for iface = 1:nf
        
        if flipFace_e(iface)
            Hfe(ind_v_L(iface,:),:) = Hfe(ind_v_L(iface,perm),:);
            Hfe(:,ind_v_L(iface,:)) = Hfe(:,ind_v_L(iface,perm));
        end
    end
    
    % store matrices
    C(:,:,iElem) = Ce;
    H(:,:,iElem) = He;
    Hf(:,:,iElem) = Hfe;
    Hdir(:,iElem) = Hdire;
end

%% Elemental matrices
function [Ce,He,Hdire,Hf] = elementalMatrices(Xe,refEl,sol,sol_hat,aux_dir,flipFace_e,iElem,Fe)

% mesh data
global neq axisym
Nv = size(refEl.NodesCoord,1);
nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);

% initialize all the matrices
He = zeros(neq*Nv,neq*nf*nv);
Hf = zeros(neq*nf*nv,neq*nf*nv);
Hdire = zeros(neq*Nv,1);

% Information of the reference element
IPw = refEl.IPweights;                 % use the velocity gauss points to integrate
N = refEl.N;
Nxi = refEl.Nxi;
Neta = refEl.Neta;
IPw_fv = refEl.IPweights1d;
N1dv = refEl.N1d;
Nx1dv = refEl.N1dxi;

% Number of Gauss points in the interior
ngauss = length(IPw);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% Initialization
Cx11 = zeros(Nv); Cx12 = Cx11; 
Cx21 = Cx11; Cx22 = Cx11; 

Cy11 = zeros(Nv); Cy12 = Cy11; 
Cy21 = Cy11; Cy22 = Cy11; 

% solution at Gauss points
sol = reshape(sol,neq,numel(sol)/neq)';
sol_g = N*sol; % ng x nvar
prova = 0;
%% VOLUME COMPUTATIONS
for g = 1:ngauss
    
    % Velocity shape functions and derivatives at the current integration point
    Nig = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);

    % gauss point position
    xg = Nig*xe;
    yg = Nig*ye;
    
    % Jacobian of the element
    J = [Nxi_g*xe  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end
    
    % Integration weight
    dvolu=IPw(g)*det(J);
    
    if axisym        
        dvolu=dvolu*xg;
    end
    
    % x and y derivatives
    invJ = inv(J);
    Nxg = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Nyg = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
       
    % Jacobian matrices
    [Ax,Ay] = jacobianMatrices(xg,yg,sol_g(g,:),iElem,g);
    
    % Contribution of the current integration point to the elemental matrix
    Cx11 = Cx11 + Nxg'*Nig*Ax(1,1)*dvolu;
    Cx12 = Cx12 + Nxg'*Nig*Ax(1,2)*dvolu;
    Cx21 = Cx21 + Nxg'*Nig*Ax(2,1)*dvolu;
    Cx22 = Cx22 + Nxg'*Nig*Ax(2,2)*dvolu;
    
    Cy11 = Cy11 + Nyg'*Nig*Ay(1,1)*dvolu;
    Cy12 = Cy12 + Nyg'*Nig*Ay(1,2)*dvolu;
    Cy21 = Cy21 + Nyg'*Nig*Ay(2,1)*dvolu;
    Cy22 = Cy22 + Nyg'*Nig*Ay(2,2)*dvolu;
    
    prova = prova + dvolu;
end

% expand the matrices
Cx = expandMatrixCv(Cx11,Cx12,Cx21,Cx22);
Cy = expandMatrixCv(Cy11,Cy12,Cy21,Cy22);
Ce = Cx+Cy;


if any(isnan(Ce))
    stop
end

%% FACES COMPUTATIONS:

% compute permutation for flipping faces
perm = setperm(neq*nv,neq);

ngauss_f = length(IPw_fv);
for iface = 1:nf
    
    % face nodes
    nodesv = faceNodesv(iface,:);
    
    % indices for local assembly
    ind_face_2 = (iface-1)*neq*nv + (1:neq*nv);
    ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
    
    xf = xe(nodesv);
    yf = ye(nodesv);
    sol_f = sol_hat(:,iface);
    if flipFace_e(iface)
        sol_f = sol_f(perm);
    end
    sol_f = transpose(reshape(sol_f,neq,nv));
    
    % Gauss point position
    xyfg = N1dv*[xf yf];
    
    if aux_dir(iface)
        % exact velocity in the faces for Dirichlet boundary        
        sol_g = analyticalSolution(xyfg);
    else
        sol_g = N1dv*sol_f;
    end
    
    % Initialization    
    Hn11 = zeros(nv); Hn12 = Hn11;
    Hn21 = Hn11; Hn22 = Hn11;
    
    
    Hdirf = zeros(nv,neq);
    
    prova = 0;
    %  LOOP IN GAUSS POINTS
    for g = 1:ngauss_f
        
        % Shape functions and derivatives at the current integration point
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
        
        
        if flipFace_e(iface)
            igauss = ngauss_f-g+1;
        else
            igauss = g;
        end
        An = jacobianMatricesFace(xyfg(g,1),xyfg(g,2),n_g(1),n_g(2),sol_g(g,:),Fe(iface),igauss);
               
        Hn11 = Hn11 + Nfv_g'*Nfv_g*An(1,1)*dline;
        Hn12 = Hn12 + Nfv_g'*Nfv_g*An(1,2)*dline;
        Hn21 = Hn21 + Nfv_g'*Nfv_g*An(2,1)*dline;
        Hn22 = Hn22 + Nfv_g'*Nfv_g*An(2,2)*dline;

        Hdirf = Hdirf + Nfv_g'*transpose(An*sol_g(g,:)')*dline;
        
        prova = prova + dline;
        
    end
    
    % expand the matrices    
    Hloc = expandMatrixCv(Hn11,Hn12,Hn21,Hn22);
    
    % elemental assembly
    He(ind2,ind_face_2) = He(ind2,ind_face_2) + ~aux_dir(iface)*Hloc;
    Hf(ind_face_2,ind_face_2) = Hf(ind_face_2,ind_face_2) + ~aux_dir(iface)*Hloc;
    Hdire(ind2) = Hdire(ind2) + aux_dir(iface)*col(Hdirf');    
end

%% additional routines

function res = expandMatrixCv(Cxx,Cxy,Cyx,Cyy)
% expand matrix Cv
%   [ Cxx Cxy 
%     Cyx Cyy]
res = zeros([size(Cxx) 2 2]);
res(:,:,[1 3 2 4]) = cat(3,Cxx,Cxy,Cyx,Cyy);
res = permute(res, [3 1 4 2]);
res = reshape(res, 2*size(Cxx,1),2*size(Cxx,2));

function [Jx, Jy] = jacobianMatrices(x,y,U,iElem,ig)

global kT Magnetic testcase useThreshold
rho = U(1);
Gamma = U(2);
if useThreshold
    if U(1)<useThreshold,
        rho = useThreshold;
        Gamma = 0;
        
    end
end
if (testcase.n >= 50 && testcase.n<60)
    b = [Magnetic.bx(ig,iElem), Magnetic.by(ig,iElem)];
else
    b = defineMagneticField([x,y]);
end
Jx = [0, b(1); (-1*Gamma^2/rho^2+kT)*b(1), 2*Gamma/rho*b(1)];
Jy = [0, b(2); (-1*Gamma^2/rho^2+kT)*b(2), 2*Gamma/rho*b(2)];




function Jn = jacobianMatricesFace(x,y,nx,ny,U,iFace,ig)

global kT  Magnetic testcase useThreshold
rho = U(1);
Gamma = U(2);

if useThreshold
    if U(1)<useThreshold,
        rho = useThreshold;
        Gamma = 0;
        
    end
end
if (testcase.n >= 50 && testcase.n<60)
    b = [Magnetic.bxfaces(ig,iFace), Magnetic.byfaces(ig,iFace)];
else
    b = defineMagneticField([x,y]);
end
bn = b(1)*nx+b(2)*ny;
Jn = [0, bn; (-1*Gamma^2/rho^2+kT)*bn, 2*Gamma/rho*bn];

