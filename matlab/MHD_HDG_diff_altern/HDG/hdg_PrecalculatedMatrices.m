function [A,B,C,C_dir,G,D,E,Df,Ef,Edir,invL,force,Lf,...
    L,P,Pb,Q,Qb,Qf] =...
    hdg_PrecalculatedMatrices(X,T,flipFace,refEl,tau,Fcon,F_dir)

% mesh data
global neq               % number of equations (2 for isothermal)
nf  = size(refEl.faceNodes,1);    % number of faces in one element
Ne = size(T,1);                    % number of elements
Nv = size(T,2);                    % number of element nodes for the velocity
nv = size(refEl.NodesCoord1d,1);   % number of face nodes for the velocity
nd = size(X,2);

% allocation and initialization
A  = zeros(neq*Nv,neq*Nv,Ne);
B = zeros(neq*nd*Nv,nd*Nv,Ne);
C = zeros(neq*nd*Nv,nf*nd*nv,Ne);
C_dir = zeros(neq*nd*Nv,Ne);
G = A;
D  = A;
E  = zeros(neq*Nv,neq*nf*nv,Ne);
Df = zeros(neq*nf*nv,neq*Nv,Ne);
Ef = zeros(neq*nf*nv,neq*nf*nv,Ne);
Edir = zeros(neq*Nv,Ne);
% PQinvL = zeros(neq*Nv,nd*neq*Nv,Ne);
invL  = zeros(neq*nd*Nv,neq*nd*Nv,Ne);
force = zeros(2*Nv,Ne);
Lf = zeros(neq*nf*nv,neq*nd*Nv,Ne);

L = zeros(neq*nd*Nv,neq*nd*Nv,Ne);
P = zeros(neq*Nv,nd*neq*Nv,Ne);
Pb = zeros(neq*Nv,nd*neq*Nv,Ne);
Q = zeros(neq*Nv,nd*neq*Nv,Ne);
Qb = zeros(neq*Nv,nd*neq*Nv,Ne);
Qf = Lf;

% local assembly indexes
ind_v_L = zeros(nf, neq*nv);
for iface = 1:nf
   ind_v_L(iface,:) = neq*nv*(iface-1)+ (1:neq*nv);
end

% compute permutation for flipping faces
perm = setperm(neq*nv,neq);

% 
% ind_1_v_L = (1:2*nv);
% ind_2_v_L = 2*nv + (1:2*nv);
% ind_3_v_L = 4*nv + (1:2*nv);

% loop in elements
for iElem = 1:Ne
    
    Fe = Fcon(iElem,:);
    Te = T(iElem,:);
    Xe = X(Te,:);
    flipFace_e = flipFace(iElem,:);
    aux_dir = F_dir(iElem,:);
    
    
    % elemental matrices
    [Ae,Be,Ce,C_dire,Ge,De,Ee,Dfe,Efe,Edire,invLe,Lfe,fe,...
        Le,Pe,Qe,Pbe,Qbe,Qfe] = ...
        elementalMatrices(Xe,refEl,tau,aux_dir,iElem,Fe,flipFace_e);
 
 for iface = 1:nf
     
     if flipFace_e(iface)
         Dfe(ind_v_L(iface,:),:) = Dfe(ind_v_L(iface,perm),:);
         Lfe(ind_v_L(iface,:),:) = Lfe(ind_v_L(iface,perm),:);
         Qfe(ind_v_L(iface,:),:) = Qfe(ind_v_L(iface,perm),:);
         Efe(ind_v_L(iface,:),:) = Efe(ind_v_L(iface,perm),:);
         Efe(:,ind_v_L(iface,:)) = Efe(:,ind_v_L(iface,perm));
     end
 end
     
    % store matrices
    A(:,:,iElem) = Ae;
    B(:,:,iElem) = Be;
    C(:,:,iElem) = Ce;
    C_dir(:,iElem) = C_dire;
    G(:,:,iElem) = Ge;
    D(:,:,iElem) = De;
    E(:,:,iElem) = Ee;
    Df(:,:,iElem) = Dfe;
    Ef(:,:,iElem) = Efe;
    Edir(:,iElem) = Edire;
%     PQinvL(:,:,iElem) = PQinvLe;
    invL(:,:,iElem) = invLe;
    force(:,iElem) = fe;
    Lf(:,:,iElem) = Lfe;
    
    L(:,:,iElem) = Le;
    P(:,:,iElem) = Pe;
    Q(:,:,iElem) = Qe;
    Pb(:,:,iElem) = Pbe;
    Qb(:,:,iElem) = Qbe;
    Qf(:,:,iElem) = Qfe;
end


%% Elemental matrices
function [A,B,C,C_dir,G,D,E,Df,Ef,Edir,invL,Lf,f,...
    L,P,Q,Pb,Qb,Qf] = ...
    elementalMatrices(Xe,refEl,tau,aux_dir,iElem,Fe,flipFace)

global kT neq diff_n diff_u Magnetic testcase axisym driftvel

% mesh data
nd = size(Xe,2);
Nv = size(refEl.NodesCoord,1);
nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);

% initialize all the matrices
C = zeros(neq*nd*Nv,nd*nf*nv);
C_dir = zeros(neq*nd*Nv,1);
Lf_transp = C;
Lf = Lf_transp';
D = zeros(neq*Nv,neq*Nv);
Df = zeros(neq*nf*nv,neq*Nv);
E = zeros(neq*Nv,neq*nf*nv);
Edir = zeros(neq*Nv,1);
Ef = zeros(neq*nf*nv,neq*nf*nv);
Q = zeros(2*Nv,4*Nv);
Qf = Lf_transp';
Qb = Q;
Qb_aux = Q;
Mass = zeros(Nv);
Ge = Mass;
fe = zeros(Nv,neq);
Bx = Mass;
By = Mass;
Px = Mass;
Py = Mass;
Pbx = Bx;
Pby = By;
if driftvel
    Dve = Mass;
end

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
prova = 0;
elarea = 0;
%% VOLUME COMPUTATIONS
for g = 1:ngauss
    
    % Velocity shape functions and derivatives at the current integration point
    Niv_g = Niv(g,:);
    Nxiv_g = Nxiv(g,:);
    Netav_g = Netav(g,:);
        
    % gauss point position
    xg = Niv_g*xe;
    yg = Niv_g*ye;
    
    % Jacobian
    J = [Nxiv_g*xe	  Nxiv_g*ye
        Netav_g*xe  Netav_g*ye];
    if det(J)<0
        error(['computeElementalMatrices: det(J)<0 in elem. ', num2str(iElem)])
    end
    
    % x and y derivatives
    invJ = inv(J);
    Nx_g = invJ(1,1)*Nxiv_g + invJ(1,2)*Netav_g;
    Ny_g = invJ(2,1)*Nxiv_g + invJ(2,2)*Netav_g;
    
    if axisym
        Nx_ax_g = Nx_g+1/xg*Niv_g;
    else
        Nx_ax_g = Nx_g;
    end
    
    
    % Integration weight
    dvolu=IPw(g)*det(J);
    
    if axisym
        dvolu=dvolu*xg;
    end
    
    % body force at integration point
    force = [0,0];
    switch testcase.n
        
        case 51 % Source between -0.88 and -0.90
            if (Magnetic.flux2D(g,iElem)<=-0.88 && Magnetic.flux2D(g,iElem)>=-0.90)
                force(1) = 3.20119388718018e-05;
            end
        case 52 % Source between -0.90 and -1
            if (Magnetic.flux2D(g,iElem)<=-0.90 && Magnetic.flux2D(g,iElem)>=-1)
                force(1) = 9.45155008295538e-06;
            end
        case 53 % Source from -0.90 
            if Magnetic.flux2D(g,iElem)<=-0.90
                force(1) = 7.24032211339971e-06;
            end
        case 54 % Source betwenn -0.98 and -1
            if (Magnetic.flux2D(g,iElem)<=-0.98  && Magnetic.flux2D(g,iElem)>=-1)
                force(1) = 5.78723047118297e-05;
            end        
        case 55 % Source from -1.03
            if Magnetic.flux2D(g,iElem)<=-1.03
                force(1) = 0.000115575293741846;
            end              
        otherwise
            force = bodyForce([xg,yg]); 
    end

    
    % magnetic field
    if (testcase.n >= 50 && testcase.n<60)
        b = [Magnetic.bx(g,iElem), Magnetic.by(g,iElem)];
        divb = Magnetic.div(g,iElem);
        if driftvel
            dv = [Magnetic.dvx(g,iElem),Magnetic.dvy(g,iElem)];
        end
    else
        if driftvel
            [b,db,dvx,dvy] = defineMagneticField([xg,yg]);
            dv = [dvx,dvy];
        else
            [b,db] = defineMagneticField([xg,yg]);
        end
        divb = sum(db);
    end
    
    
    % Contribution of the current integration point to the elemental matrix
    Mass = Mass + Niv_g'*Niv_g*dvolu;
    Bx = Bx + Nx_ax_g'*Niv_g*dvolu;
    By = By + Ny_g'*Niv_g*dvolu;   
    Px = Px + Nx_g'*Niv_g*dvolu;
    Py = Py + Ny_g'*Niv_g*dvolu;      
    Pbx = Pbx + b(1)*(b(1)*Nx_g'+b(2)*Ny_g')*Niv_g*dvolu;
    Pby = Pby + b(2)*(b(1)*Nx_g'+b(2)*Ny_g')*Niv_g*dvolu;
    Ge = Ge - kT*divb*(Niv_g'*Niv_g)*dvolu;
    fe = fe + Niv_g'*force*dvolu;
    prova = prova + Niv_g'*dvolu;
    elarea = elarea + dvolu;
    
    % Formulation 1 of the drift velocity term
    if driftvel
        Dve = Dve + Niv_g'*(dv(1)*Nx_g + dv(2)*Ny_g)*dvolu;
    end
end

% expand the matrices
L = expandMatrixA(Mass,neq*nd);
B = expandMatrixB(Bx,By);
P = transpose(expandMatrixB(Px',Py'));
Pb = transpose(expandMatrixB(Pbx',Pby'));
A = expandMatrixA(Mass,neq);
G = expandMatrixCv(zeros(size(Ge)),zeros(size(Ge)),Ge,zeros(size(Ge)));
f = col(fe');
% Formulation 1 of the drift velocity term
if driftvel
    Dv = expandMatrixA(Dve,neq);
    G = G+Dv;
end
%********************************************************
%% FACES COMPUTATIONS:
ngauss_f = length(IPw_fv);
for iface = 1:nf
    
    % face nodes
    nodesv = faceNodesv(iface,:);
    iFace = Fe(iface);
    
    % indices for local assembly
    ind_face_2 = (iface-1)*neq*nv + (1:neq*nv);
    ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
    ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for velocity gradient
    
    xf = xe(nodesv);
    yf = ye(nodesv);
    
    xyfg = N1dv*[xf yf];
    u_ex = analyticalSolution(xyfg);
    tau_f = tau;
    
    % initialize local matrices
    Massf = zeros(nv);
    E_dirf = zeros(nv,neq);
    Cnx = zeros(nv);
    Cny = Cnx;
    Qbx = Cnx;
    Qby = Cnx;
    Qbx_aux = Cnx; Qby_aux = Cnx;
    Cuxnx_dir = zeros(nv,1);
    Cuxny_dir = Cuxnx_dir;
    Cuynx_dir = Cuxnx_dir;
    Cuyny_dir = Cuxnx_dir;    
    
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
        
        % magnetic field
        if (testcase.n >= 50 && testcase.n<60)
            if flipFace(iface)
                igauss = ngauss_f-g+1;
            else
                igauss = g;
            end
            b = [Magnetic.bxfaces(igauss,iFace), Magnetic.byfaces(igauss,iFace)];
        else
            b = defineMagneticField([xyfg(g,1),xyfg(g,2)]);
        end
        
        % Contribution of the current integration point to the elemental matrix
        Cnx = Cnx + Nfv_g'*Nfv_g*n_g(1)*dline;
        Cny = Cny + Nfv_g'*Nfv_g*n_g(2)*dline;
        Cuxnx_dir = Cuxnx_dir + Nfv_g'*u_ex(g,1)*n_g(1)*dline;
        Cuxny_dir = Cuxny_dir + Nfv_g'*u_ex(g,1)*n_g(2)*dline;
        Cuynx_dir = Cuynx_dir + Nfv_g'*u_ex(g,2)*n_g(1)*dline;
        Cuyny_dir = Cuyny_dir + Nfv_g'*u_ex(g,2)*n_g(2)*dline;       
        Qbx = Qbx + b(1)*Nfv_g'*(b(1)*Nfv_g*n_g(1)+b(2)*Nfv_g*n_g(2))*dline;
        Qby = Qby + b(2)*Nfv_g'*(b(1)*Nfv_g*n_g(1)+b(2)*Nfv_g*n_g(2))*dline;            
        E_dirf = E_dirf + tau_f*(Nfv_g'*u_ex(g,:))*dline;
        Massf = Massf + tau_f *(Nfv_g'*Nfv_g)*dline;
     end
    
    % expand the matrices
    C_loc = expandMatrixB(Cnx,Cny);
    Qb_loc = transpose( expandMatrixB(Qbx',Qby') );
    Qb_loc_aux = transpose( expandMatrixB(Qbx_aux',Qby_aux') );
    C_loc_dir = expandMatrixF(expandMatrixF(Cuxnx_dir,Cuynx_dir),...
        expandMatrixF(Cuxny_dir,Cuyny_dir));
    D_loc = expandMatrixA(Massf,neq);
    E_loc = expandMatrixA(Massf,neq);

    % elemental assembly
    C(ind4,ind_face_2) = C(ind4,ind_face_2) + ~aux_dir(iface)*C_loc;
    C_dir(ind4) = C_dir(ind4) + aux_dir(iface)*C_loc_dir;
%     Lf_transp(ind4,ind_face_2) = Lf_transp(ind4,ind_face_2) + C_loc;
    Lf(ind_face_2,ind4) = Lf(ind_face_2,ind4) +  ~aux_dir(iface)*transpose( expandMatrixB(Cnx',Cny') );
    Q(ind2,ind4) = Q(ind2,ind4) + C_loc';
    D(ind2,ind2) = D(ind2,ind2) + D_loc;
    Df(ind_face_2,ind2) = Df(ind_face_2,ind2) + ~aux_dir(iface)*D_loc;
    E(ind2,ind_face_2) = E(ind2,ind_face_2) + ~aux_dir(iface)*E_loc;
    Edir(ind2) = Edir(ind2) + aux_dir(iface)*col(E_dirf');
    Ef(ind_face_2,ind_face_2) = Ef(ind_face_2,ind_face_2) + ~aux_dir(iface)*E_loc;
    Qb(ind2,ind4) = Qb(ind2,ind4) + Qb_loc;
    Qf(ind_face_2,ind4) = Qf(ind_face_2,ind4) + ~aux_dir(iface)*Qb_loc;
    Qb_aux(ind2,ind4) = Qb_aux(ind2,ind4) + Qb_loc_aux;
end

diffdiag = repmat(col(repmat([diff_n diff_u],nd,1)),Nv,1);
Diff = diag(diffdiag);
% Lf = Lf_transp'*Diff;
Lf = Lf*Diff;
Qf = Qf*Diff;

invL = L\eye(size(L));
% PQinvL = ( (P-Pb) - (Q-Qb))*Diff*invL;
P  = P*Diff;
Pb = Pb*Diff;
Q  = Q*Diff;
Qb = Qb*Diff;
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

% function res = expandMatrixB(Bx,By)
% % expand matrix B
% %   [ Bx  0
% %     By  0
% %     0  Bx
% %     0  By ]
% res = zeros([size(Bx) 4 2]);
% res(:,:,[1 2 7 8]) = cat(3,Bx,By,Bx,By);
% res = permute(res, [3 1 4 2]);
% res = reshape(res, 4*size(Bx,1),2*size(Bx,2));

function res = expandMatrixF(Fx,Fy)
% expand matrix
%   [ Fx
%     Fy ]
res = zeros([size(Fx) 2 1]);
res(:,:,[1 2]) = cat(3,Fx,Fy);
res = permute(res, [3 1 4 2]);
res = reshape(res, 2*size(Fx,1),size(Fx,2));


function res = expandMatrixCv(Cxx,Cxy,Cyx,Cyy)
% expand matrix Cv
%   [ Cxx Cxy 
%     Cyx Cyy]
res = zeros([size(Cxx) 2 2]);
res(:,:,[1 3 2 4]) = cat(3,Cxx,Cxy,Cyx,Cyy);
res = permute(res, [3 1 4 2]);
res = reshape(res, 2*size(Cxx,1),2*size(Cxx,2));

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