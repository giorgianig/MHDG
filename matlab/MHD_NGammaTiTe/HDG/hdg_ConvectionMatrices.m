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
Cx11 = zeros(Nv); Cx12 = Cx11; Cx13 = Cx11; Cx14 = Cx11;
Cx21 = Cx11; Cx22 = Cx11; Cx23 = Cx11; Cx24 = Cx11;
Cx31 = Cx11; Cx32 = Cx11; Cx33 = Cx11; Cx34 = Cx11;
Cx41 = Cx11; Cx42 = Cx11; Cx43 = Cx11; Cx44 = Cx11;

Cy11 = zeros(Nv); Cy12 = Cy11; Cy13 = Cy11; Cy14 = Cy11;
Cy21 = Cy11; Cy22 = Cy11; Cy23 = Cy11; Cy24 = Cy11;
Cy31 = Cy11; Cy32 = Cy11; Cy33 = Cy11; Cy34 = Cy11;
Cy41 = Cy11; Cy42 = Cy11; Cy43 = Cy11; Cy44 = Cy11;

G11 = zeros(Nv); G12 = G11; G13 = G11; G14 = G11;
G21 = G11; G22 = G11; G23 = G11; G24 = G11;
G31 = G11; G32 = G11; G33 = G11; G34 = G11;
G41 = G11; G42 = G11; G43 = G11; G44 = G11;

% solution at Gauss points
sol = reshape(sol,neq,numel(sol)/neq)';
sol_g = N*sol; % ng x nvar

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

    G = GimplicitMatrix(xg,yg,sol_g(g,:),iElem,g);
    
    % Contribution of the current integration point to the elemental matrix
    Cx11 = Cx11 + Nxg'*Nig*Ax(1,1)*dvolu;
    Cx12 = Cx12 + Nxg'*Nig*Ax(1,2)*dvolu;
    Cx13 = Cx13 + Nxg'*Nig*Ax(1,3)*dvolu;
    Cx14 = Cx14 + Nxg'*Nig*Ax(1,4)*dvolu;
    Cx21 = Cx21 + Nxg'*Nig*Ax(2,1)*dvolu;
    Cx22 = Cx22 + Nxg'*Nig*Ax(2,2)*dvolu;
    Cx23 = Cx23 + Nxg'*Nig*Ax(2,3)*dvolu;
    Cx24 = Cx24 + Nxg'*Nig*Ax(2,4)*dvolu;
    Cx31 = Cx31 + Nxg'*Nig*Ax(3,1)*dvolu;
    Cx32 = Cx32 + Nxg'*Nig*Ax(3,2)*dvolu;
    Cx33 = Cx33 + Nxg'*Nig*Ax(3,3)*dvolu;
    Cx34 = Cx34 + Nxg'*Nig*Ax(3,4)*dvolu;
    Cx41 = Cx41 + Nxg'*Nig*Ax(4,1)*dvolu;
    Cx42 = Cx42 + Nxg'*Nig*Ax(4,2)*dvolu;
    Cx43 = Cx43 + Nxg'*Nig*Ax(4,3)*dvolu;
    Cx44 = Cx44 + Nxg'*Nig*Ax(4,4)*dvolu;

    
    
    Cy11 = Cy11 + Nyg'*Nig*Ay(1,1)*dvolu;
    Cy12 = Cy12 + Nyg'*Nig*Ay(1,2)*dvolu;
    Cy13 = Cy13 + Nyg'*Nig*Ay(1,3)*dvolu;
    Cy14 = Cy14 + Nyg'*Nig*Ay(1,4)*dvolu;    
    Cy21 = Cy21 + Nyg'*Nig*Ay(2,1)*dvolu;
    Cy22 = Cy22 + Nyg'*Nig*Ay(2,2)*dvolu;
    Cy23 = Cy23 + Nyg'*Nig*Ay(2,3)*dvolu;
    Cy24 = Cy24 + Nyg'*Nig*Ay(2,4)*dvolu;
    Cy31 = Cy31 + Nyg'*Nig*Ay(3,1)*dvolu;
    Cy32 = Cy32 + Nyg'*Nig*Ay(3,2)*dvolu;
    Cy33 = Cy33 + Nyg'*Nig*Ay(3,3)*dvolu;
    Cy34 = Cy34 + Nyg'*Nig*Ay(3,4)*dvolu;
    Cy41 = Cy41 + Nyg'*Nig*Ay(4,1)*dvolu;
    Cy42 = Cy42 + Nyg'*Nig*Ay(4,2)*dvolu;
    Cy43 = Cy43 + Nyg'*Nig*Ay(4,3)*dvolu;
    Cy44 = Cy44 + Nyg'*Nig*Ay(4,4)*dvolu;

    
    G11 = G11 + Nig'*Nig*G(1,1)*dvolu;
    G12 = G12 + Nig'*Nig*G(1,2)*dvolu;
    G13 = G13 + Nig'*Nig*G(1,3)*dvolu;
    G14 = G14 + Nig'*Nig*G(1,4)*dvolu;
    G21 = G21 + Nig'*Nig*G(2,1)*dvolu;
    G22 = G22 + Nig'*Nig*G(2,2)*dvolu;
    G23 = G23 + Nig'*Nig*G(2,3)*dvolu;
    G24 = G24 + Nig'*Nig*G(2,4)*dvolu;
    G31 = G31 + Nig'*Nig*G(3,1)*dvolu;
    G32 = G32 + Nig'*Nig*G(3,2)*dvolu;
    G33 = G33 + Nig'*Nig*G(3,3)*dvolu;
    G34 = G34 + Nig'*Nig*G(3,4)*dvolu;
    G41 = G41 + Nig'*Nig*G(4,1)*dvolu;
    G42 = G42 + Nig'*Nig*G(4,2)*dvolu;
    G43 = G43 + Nig'*Nig*G(4,3)*dvolu;
    G44 = G44 + Nig'*Nig*G(4,4)*dvolu;    
    
    
end

% expand the matrices
Cx = expandMatrixCv(Cx11,Cx12,Cx13,Cx14,Cx21,Cx22,Cx23,Cx24,Cx31,Cx32,Cx33,Cx34,Cx41,Cx42,Cx43,Cx44);
Cy = expandMatrixCv(Cy11,Cy12,Cy13,Cy14,Cy21,Cy22,Cy23,Cy24,Cy31,Cy32,Cy33,Cy34,Cy41,Cy42,Cy43,Cy44);
GG = expandMatrixCv(G11,G12,G13,G14,G21,G22,G23,G24,G31,G32,G33,G34,G41,G42,G43,G44);
Ce = Cx+Cy+GG;


if any(isnan(Ce))
    stop
end

if ~isreal(Ce)
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
    Hn11 = zeros(nv); Hn12 = Hn11; Hn13 = Hn11; Hn14 = Hn11;
    Hn21 = Hn11; Hn22 = Hn11; Hn23 = Hn11; Hn24 = Hn11;
    Hn31 = Hn11; Hn32 = Hn11; Hn33 = Hn11; Hn34 = Hn11;
    Hn41 = Hn11; Hn42 = Hn11; Hn43 = Hn11; Hn44 = Hn11;
    
    
    Hdirf = zeros(nv,neq);
    
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
        Hn13 = Hn13 + Nfv_g'*Nfv_g*An(1,3)*dline;
        Hn14 = Hn14 + Nfv_g'*Nfv_g*An(1,4)*dline;
        Hn21 = Hn21 + Nfv_g'*Nfv_g*An(2,1)*dline;
        Hn22 = Hn22 + Nfv_g'*Nfv_g*An(2,2)*dline;
        Hn23 = Hn23 + Nfv_g'*Nfv_g*An(2,3)*dline;
        Hn24 = Hn24 + Nfv_g'*Nfv_g*An(2,4)*dline;
        Hn31 = Hn31 + Nfv_g'*Nfv_g*An(3,1)*dline;
        Hn32 = Hn32 + Nfv_g'*Nfv_g*An(3,2)*dline;
        Hn33 = Hn33 + Nfv_g'*Nfv_g*An(3,3)*dline;
        Hn34 = Hn34 + Nfv_g'*Nfv_g*An(3,4)*dline;
        Hn41 = Hn41 + Nfv_g'*Nfv_g*An(4,1)*dline;
        Hn42 = Hn42 + Nfv_g'*Nfv_g*An(4,2)*dline;
        Hn43 = Hn43 + Nfv_g'*Nfv_g*An(4,3)*dline;
        Hn44 = Hn44 + Nfv_g'*Nfv_g*An(4,4)*dline;
        Hdirf = Hdirf + Nfv_g'*transpose(An*sol_g(g,:)')*dline;
        
    end
    
    % expand the matrices    
    Hloc = expandMatrixCv(Hn11,Hn12,Hn13,Hn14,Hn21,Hn22,Hn23,Hn24,Hn31,Hn32,Hn33,Hn34,Hn41,Hn42,Hn43,Hn44);

 
    % elemental assembly
    He(ind2,ind_face_2) = He(ind2,ind_face_2) + ~aux_dir(iface)*Hloc;
    Hf(ind_face_2,ind_face_2) = Hf(ind_face_2,ind_face_2) + ~aux_dir(iface)*Hloc;
    Hdire(ind2) = Hdire(ind2) + aux_dir(iface)*col(Hdirf');    
end

if ~isreal(He)
    stop
end
%% additional routines

function res = expandMatrixCv(C11,C12,C13,C14,C21,C22,C23,C24,C31,C32,C33,C34,C41,C42,C43,C44)
% expand matrix Cv
%   [ C11 C12 C13 C14
%     C21 C22 C23 C24
%     C31 C32 C33 C34
%     C41 C42 C43 C44]
res = zeros([size(C11) 4 4]);
res(:,:,[1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16]) = cat(3,C11,C12,C13,C14,C21,C22,C23,C24,C31,C32,C33,C34,C41,C42,C43,C44);
res = permute(res, [3 1 4 2]);
res = reshape(res, 4*size(C11,1),4*size(C11,2));

function [Jx, Jy] = jacobianMatrices(x,y,U,iElem,ig)

global Magnetic testcase useThreshold  decoupleEquations Mref
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
%
%
%
if ~decoupleEquations
    Jx = b(1)*[                         0,                                                                                             1,                                        0,            0; ...
                    -2/3*U(2)^2/U(1)^2                                                                         4/3*U(2)/U(1)                                      2/3,        2/3; ...
                    -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2                 5/3*U(2)/U(1),    0 ;   ...
                    -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                                     0,       5/3*U(2)/U(1)];

    
    Jy = b(2)*[                         0,                                                                                             1,                                        0,            0; ...
                    -2/3*U(2)^2/U(1)^2                                                                         4/3*U(2)/U(1)                                      2/3,        2/3; ...
                    -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2                  5/3*U(2)/U(1),    0 ;   ...
                    -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                                     0,       5/3*U(2)/U(1)];
else
    Jx = b(1)*[                         0,                                                                                             1,                                        0,            0; ...
                                (-1*Gamma^2/rho^2+Mref)                                                      2*Gamma/rho                                     0     ,        0; ...
                                -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                5/3*U(3)/U(1)-U(2)^2/U(1)^2                 5/3*U(2)/U(1),    0 ;   ...
                                -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                         0,       5/3*U(2)/U(1)];
    
    
    Jy = b(2)*[                         0,                                                                                             1,                                        0,            0; ...
                            (-1*Gamma^2/rho^2+Mref)                                                         2*Gamma/rho                                     0     ,        0; ...
                            -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2            5/3*U(2)/U(1),    0 ;   ...
                            -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                            0,       5/3*U(2)/U(1)];
    
end


function Jn = jacobianMatricesFace(x,y,nx,ny,U,iFace,ig)

global Magnetic testcase useThreshold  decoupleEquations Mref
rho = U(1);
Gamma = U(2);

if useThreshold
    if U(1)<useThreshold
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
if ~decoupleEquations
    Jn = bn*[                         0,                                                                                             1,                                        0,            0; ...
                    -2/3*U(2)^2/U(1)^2                                                                         4/3*U(2)/U(1)                                      2/3,        2/3; ...
                   -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2                 5/3*U(2)/U(1),    0 ;   ...
                   -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                                     0,       5/3*U(2)/U(1)];
else
    Jn = bn*[                         0,                                                                                             1,                                        0,            0; ...
                            (-1*Gamma^2/rho^2+Mref)                                                         2*Gamma/rho                                     0     ,        0; ...
                            -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2            5/3*U(2)/U(1),    0 ;   ...
                            -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                            0,       5/3*U(2)/U(1)];

end


         
         
         
function G = GimplicitMatrix(x,y,U,iElem,ig)

global  Magnetic testcase useThreshold  decoupleEquations Mref
rho = U(1);
Gamma = U(2);
if useThreshold
    if U(1)<useThreshold
        rho = useThreshold;
        Gamma = 0;
    end
end

% magnetic field
if (testcase.n >= 50 && testcase.n<60)
    divb = Magnetic.div(ig,iElem);
else
    [~,db] = defineMagneticField([x,y]);
    divb = sum(db);
end
if ~decoupleEquations
G = divb*[                0,                                 0,                                    0                               0; ...
                  1/3*U(2)^2/U(1)^2            -2/3*U(2)/U(1)                      2/3,                             2/3;...
                                0                                  0                                       0                              0; ...
                                0                                   0                                      0                              0];
else
    G = Mref*divb*[            0,                                 0                                      0                             0; ...
                                1                                  0                                      0                             0;...
                                0                                  0                                      0                             0; ...
                                0                                  0                                      0                             0];
end
