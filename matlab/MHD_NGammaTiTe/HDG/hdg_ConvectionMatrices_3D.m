function [ C,H,Hdir,Hf] =...
    hdg_ConvectionMatrices_3D(X,T,F,flipFace,refEl,sol,sol_hat,F_dir)

% mesh data
global ntor neq refElTor theta             % number of equations (2 for isothermal)
nf  = size(refEl.faceNodes,1);    % number of faces in one element
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Np2d      = size(refEl.NodesCoord,1);
Nv     = Np1dTor*Np2d;
N2d    = size(T,1);                    % number of elements
Nfl   = Np1dPol*Np1dTor;
Nfp  = Np2d*2+nf*Np1dPol*Np1dTor;
Ne = N2d*ntor;
Nf = max(max(F));
% nv = size(refEl.NodesCoord1d,1);   % number of face nodes for the velocity
% Ndim = 3;
tdiv = linspace(0,theta,ntor+1);
nDirFaces =  sum(sum(F_dir));

% allocation and initialization
C = zeros(neq*Nv,neq*Nv,Ne);
H = zeros(neq*Nv,neq*Nfp,Ne);
Hdir = zeros(neq*Nv,Ne);
Hf = zeros(neq*Nfp,neq*Nfp,Ne);

% compute permutation for flipping faces
perm =  setperm3d(Np1dPol,Np1dTor,neq);
ind_dim = neq*[Np2d,repmat(Np1dPol*Np1dTor,1,nf),Np2d];
ind_sta = [1, 1+cumsum(ind_dim)];
ind_uhat = zeros(Nfp*neq,1);

% loop in elements
for itor = 1:ntor
    tel = tdiv(itor)+0.5*(refElTor.NodesCoord1d+1)*(tdiv(itor+1)-tdiv(itor));
    for iel = 1:N2d
        Fe = F(iel,:);
        Te = T(iel,:);
        Xe = X(Te,:);
        flipFace_e = flipFace(iel,:);
        aux_dir = F_dir(iel,:);
        iElem = (itor-1)*N2d + iel;
        ind = (iElem-1)*neq*Nv + (1:neq*Nv);
        sol_e = sol(ind);
        delta = 1+(itor-1)*(N2d*Np2d+(Nf-nDirFaces)*Nfl)*neq+neq*[(iel-1)*Np2d,N2d*Np2d+...
            (Fe-1)*Nfl,N2d*Np2d+(Nf-nDirFaces)*Nfl+(iel-1)*Np2d];
        if itor==ntor
            delta(end) = 1+(iel-1)*Np2d*neq;
        end
        for iface = 1:nf+2
            ind_loc = ind_sta(iface)+(0:ind_dim(iface)-1);
            ind_uhat(ind_loc) = delta(iface)+(0:ind_dim(iface)-1);
        end
        sol_hat_e = sol_hat(ind_uhat);
        
        % elemental matrices
        [Ce,He,Hdire,Hfe] = elementalMatrices(Xe,tel,refEl,sol_e,sol_hat_e,aux_dir,flipFace_e,iElem,Fe);
        
        for iface = 1:nf
            if flipFace_e(iface)
                ind_v_L = Np2d*neq + (iface-1)*Nfl*neq+ (1:Nfl*neq);
                Hfe(ind_v_L,:) = Hfe(ind_v_L(perm),:);
                Hfe(:,ind_v_L) = Hfe(:,ind_v_L(perm));
            end
        end
        
        % store matrices
        C(:,:,iElem) = Ce;
        H(:,:,iElem) = He;
        Hf(:,:,iElem) = Hfe;
        Hdir(:,iElem) = Hdire;
    end
end

%% Elemental matrices
function [Ce,He,Hdire,Hf] = elementalMatrices(Xe,tel,refEl,sol,sol_hat,aux_dir,flipFace_e,iElem,Fe)

global neq axisym refElTor

% mesh data
htheta   = tel(end)-tel(1);

% Ndim = 3;
nf = size(refEl.faceNodes,1);
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);     % Number of 1d points for toroidal line
Np2d      = size(refEl.NodesCoord,1);
Nv = Np1dTor*Np2d;
Nfl   = Np1dPol*Np1dTor;
Nfp  = Np2d*2+nf*Np1dPol*Np1dTor;

% initialize all the matrices
He = zeros(neq*Nv,neq*Nfp);
Hf = zeros(neq*Nfp,neq*Nfp);
Hdire = zeros(neq*Nv,1);


% Information of the reference element
ipw2d  = refEl.IPweights;
ipw1dp = refEl.IPweights1d;
ipw1dt = refElTor.IPweights1d;
N2d       = refEl.N;
N2dxi    = refEl.Nxi;
N2deta  = refEl.Neta;
N1dPol  = refEl.N1d;
N1dxPol = refEl.N1dxi;
N1dTor  = refElTor.N1d;
N1dxTor= refElTor.N1dxi;
ngauss2d = length(ipw2d);
ngauss1dtor = length(ipw1dt);
ngauss1dpol = length(ipw1dp);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);
teg = N1dTor*tel;

ngausstor = ngauss1dtor;
ngausspol = ngauss2d;

% Initialization

Cx11 = zeros(Nv); Cx12 = Cx11; Cx13 = Cx11; Cx14 = Cx11;
Cx21 = Cx11; Cx22 = Cx11; Cx23 = Cx11; Cx24 = Cx11;
Cx31 = Cx11; Cx32 = Cx11; Cx33 = Cx11; Cx34 = Cx11;
Cx41 = Cx11; Cx42 = Cx11; Cx43 = Cx11; Cx44 = Cx11;

Cy11 = zeros(Nv); Cy12 = Cy11; Cy13 = Cy11; Cy14 = Cy11;
Cy21 = Cy11; Cy22 = Cy11; Cy23 = Cy11; Cy24 = Cy11;
Cy31 = Cy11; Cy32 = Cy11; Cy33 = Cy11; Cy34 = Cy11;
Cy41 = Cy11; Cy42 = Cy11; Cy43 = Cy11; Cy44 = Cy11;


Ct11 = zeros(Nv); Ct12 = Ct11; Ct13 = Ct11; Ct14 = Ct11;
Ct21 = Ct11; Ct22 = Ct11; Ct23 = Ct11; Ct24 = Ct11;
Ct31 = Ct11; Ct32 = Ct11; Ct33 = Ct11; Ct34 = Ct11;
Ct41 = Ct11; Ct42 = Ct11; Ct43 = Ct11; Ct44 = Ct11;


G11 = zeros(Nv); G12 = G11; G13 = G11; G14 = G11;
G21 = G11; G22 = G11; G23 = G11; G24 = G11;
G31 = G11; G32 = G11; G33 = G11; G34 = G11;
G41 = G11; G42 = G11; G43 = G11; G44 = G11;

% reshape solution
sol = reshape(sol,neq,numel(sol)/neq)';

prova = 0;

xyg = N2d*Xe;
b = defineMagneticField_3D([xyg(:,1),xyg(:,2)],teg);

%% VOLUME COMPUTATIONS
for igtor = 1:ngausstor
    
    N1_g   = N1dTor(igtor,:);                             % toroidal shape functions
    N1x_g = N1dxTor(igtor,:)*2/htheta;             % toroidal shape functions derivative
    tg       = teg(igtor);
    dvolu1d = 0.5*ipw1dt(igtor)*htheta;           % toroidal contribution to elemental volume
    
    % Loop in 2D Gauss points
    for igpol = 1:ngausspol
        
        g = (igtor-1)*ngausspol+igpol;
        
        % Velocity shape functions and derivatives at the current integration point
        N2_g = N2d(igpol,:);
        N2dxi_g = N2dxi(igpol,:);
        N2deta_g = N2deta(igpol,:);
        
        % gauss point position
        xg = N2_g*xe;
        yg = N2_g*ye;
        
        % Jacobian
        J = [N2dxi_g*xe	  N2dxi_g*ye
            N2deta_g*xe  N2deta_g*ye];
        if det(J)<0
            error(['computeElementalMatrices: det(J)<0 in elem. ', num2str(iel)])
        end
        
        % x and y derivatives
        invJ = inv(J);
        Nx_g = invJ(1,1)*N2dxi_g + invJ(1,2)*N2deta_g;
        Ny_g = invJ(2,1)*N2dxi_g + invJ(2,2)*N2deta_g;
        
        if axisym
            Nt_ax_g = 1/xg*N1x_g;
        else
            Nt_ax_g = N1x_g;
        end
        
        %% 3d shape functions
        Ni      = colt(N2_g'*N1_g);
        Nx      = colt(Nx_g'*N1_g);
        Ny      = colt(Ny_g'*N1_g);
        Nt      = colt(N2_g'*Nt_ax_g);
        
        % solution at Gauss points
        sol_g = Ni*sol; % ng x nvar
        
        % Integration weight
        dvolu=ipw2d(igpol)*det(J)*dvolu1d;
        
        if axisym
            dvolu=dvolu*xg;
        end
        % Jacobian matrices
        [Ax,Ay,At] = jacobianMatrices(xg,yg,tg,sol_g,iElem,g);
        G = GimplicitMatrix(xg,yg,tg,sol_g,iElem,g);
        
        
        
        
% A = jacobianMatricesNOMAGFIELD(sol_g);        
        % Contribution of the current integration point to the elemental matrix
        
        NiN = Ni'*Ni;
        NxN = Nx'*Ni;
        NyN = Ny'*Ni;
        NtN = Nt'*Ni;
        
        Cx11 = Cx11 + NxN*Ax(1,1)*dvolu;
        Cx12 = Cx12 + NxN*Ax(1,2)*dvolu;
        Cx13 = Cx13 + NxN*Ax(1,3)*dvolu;
        Cx14 = Cx14 + NxN*Ax(1,4)*dvolu;
        Cx21 = Cx21 + NxN*Ax(2,1)*dvolu;
        Cx22 = Cx22 + NxN*Ax(2,2)*dvolu;
        Cx23 = Cx23 + NxN*Ax(2,3)*dvolu;
        Cx24 = Cx24 + NxN*Ax(2,4)*dvolu;
        Cx31 = Cx31 + NxN*Ax(3,1)*dvolu;
        Cx32 = Cx32 + NxN*Ax(3,2)*dvolu;
        Cx33 = Cx33 + NxN*Ax(3,3)*dvolu;
        Cx34 = Cx34 + NxN*Ax(3,4)*dvolu;
        Cx41 = Cx41 + NxN*Ax(4,1)*dvolu;
        Cx42 = Cx42 + NxN*Ax(4,2)*dvolu;
        Cx43 = Cx43 + NxN*Ax(4,3)*dvolu;
        Cx44 = Cx44 + NxN*Ax(4,4)*dvolu;
        
        
        Cy11 = Cy11 + NyN*Ay(1,1)*dvolu;
        Cy12 = Cy12 + NyN*Ay(1,2)*dvolu;
        Cy13 = Cy13 + NyN*Ay(1,3)*dvolu;
        Cy14 = Cy14 + NyN*Ay(1,4)*dvolu;
        Cy21 = Cy21 + NyN*Ay(2,1)*dvolu;
        Cy22 = Cy22 + NyN*Ay(2,2)*dvolu;
        Cy23 = Cy23 + NyN*Ay(2,3)*dvolu;
        Cy24 = Cy24 + NyN*Ay(2,4)*dvolu;
        Cy31 = Cy31 + NyN*Ay(3,1)*dvolu;
        Cy32 = Cy32 + NyN*Ay(3,2)*dvolu;
        Cy33 = Cy33 + NyN*Ay(3,3)*dvolu;
        Cy34 = Cy34 + NyN*Ay(3,4)*dvolu;
        Cy41 = Cy41 + NyN*Ay(4,1)*dvolu;
        Cy42 = Cy42 + NyN*Ay(4,2)*dvolu;
        Cy43 = Cy43 + NyN*Ay(4,3)*dvolu;
        Cy44 = Cy44 + NyN*Ay(4,4)*dvolu;
        
        
        Ct11 = Ct11 + NtN*At(1,1)*dvolu;
        Ct12 = Ct12 + NtN*At(1,2)*dvolu;
        Ct13 = Ct13 + NtN*At(1,3)*dvolu;
        Ct14 = Ct14 + NtN*At(1,4)*dvolu;
        Ct21 = Ct21 + NtN*At(2,1)*dvolu;
        Ct22 = Ct22 + NtN*At(2,2)*dvolu;
        Ct23 = Ct23 + NtN*At(2,3)*dvolu;
        Ct24 = Ct24 + NtN*At(2,4)*dvolu;
        Ct31 = Ct31 + NtN*At(3,1)*dvolu;
        Ct32 = Ct32 + NtN*At(3,2)*dvolu;
        Ct33 = Ct33 + NtN*At(3,3)*dvolu;
        Ct34 = Ct34 + NtN*At(3,4)*dvolu;
        Ct41 = Ct41 + NtN*At(4,1)*dvolu;
        Ct42 = Ct42 + NtN*At(4,2)*dvolu;
        Ct43 = Ct43 + NtN*At(4,3)*dvolu;
        Ct44 = Ct44 + NtN*At(4,4)*dvolu;
        
        
        G11 = G11 + NiN*G(1,1)*dvolu;
        G12 = G12 + NiN*G(1,2)*dvolu;
        G13 = G13 + NiN*G(1,3)*dvolu;
        G14 = G14 + NiN*G(1,4)*dvolu;
        G21 = G21 + NiN*G(2,1)*dvolu;
        G22 = G22 + NiN*G(2,2)*dvolu;
        G23 = G23 + NiN*G(2,3)*dvolu;
        G24 = G24 + NiN*G(2,4)*dvolu;
        G31 = G31 + NiN*G(3,1)*dvolu;
        G32 = G32 + NiN*G(3,2)*dvolu;
        G33 = G33 + NiN*G(3,3)*dvolu;
        G34 = G34 + NiN*G(3,4)*dvolu;
        G41 = G41 + NiN*G(4,1)*dvolu;
        G42 = G42 + NiN*G(4,2)*dvolu;
        G43 = G43 + NiN*G(4,3)*dvolu;
        G44 = G44 + NiN*G(4,4)*dvolu;
        
        prova = prova + dvolu;
    end
end
% expand the matrices
Cx = expandMatrixCv(Cx11,Cx12,Cx13,Cx14,Cx21,Cx22,Cx23,Cx24,Cx31,Cx32,Cx33,Cx34,Cx41,Cx42,Cx43,Cx44);
Cy = expandMatrixCv(Cy11,Cy12,Cy13,Cy14,Cy21,Cy22,Cy23,Cy24,Cy31,Cy32,Cy33,Cy34,Cy41,Cy42,Cy43,Cy44);
Ct = expandMatrixCv(Ct11,Ct12,Ct13,Ct14,Ct21,Ct22,Ct23,Ct24,Ct31,Ct32,Ct33,Ct34,Ct41,Ct42,Ct43,Ct44);

GG = expandMatrixCv(G11,G12,G13,G14,G21,G22,G23,G24,G31,G32,G33,G34,G41,G42,G43,G44);
Ce = Cx+Cy+Ct+GG;

if any(isnan(Ce))
    stop
end

%% FACES COMPUTATIONS:

% compute permutation for flipping faces
perm =  setperm3d(Np1dPol,Np1dTor,neq);

%********************************************************
%% FACES COMPUTATIONS:
for iface = 1:(nf+2)
    % face nodes
    if iface==1
        nv = Np2d;
        ind_ff = 1:Np2d*neq;                              % assembly face-face
        ind_fe = 1:Np2d*neq;                              % assembly face-element (for var)
        xyfg     = N2d*Xe;                                % Gauss points coordinates of the face (x and y components)
        thetafg  = tel(1);                                % Gauss points coordinates of the face (theta component)
        Nfi      = N2d;                                   % Shape functions in the face
        J11 = N2dxi*Xe(:,1);                              % ng x 1
        J21 = N2dxi*Xe(:,2);                              % ng x 1
        J12 = N2deta*Xe(:,1);                             % ng x 1
        J22 = N2deta*Xe(:,2);                             % ng x 1
        detJ = J11.*J22 - J21.*J12;                       % determinant of the Jacobian
        ngausspol = ngauss2d;                             % number of Gauss points in the poloidal plane
        ngausstor = 1;                                    % number of Gauss points in the toroidal direction
        ngauss = ngausspol*ngausstor;
        ipw = ipw2d;
        dsurf = detJ.*ipw;                                % elemental surface
        n_g = repmat([0,0,-1],ngauss,1);                  % outer normal
        rg = xyfg(:,1);                                   % radius
        isdir = 0;                                        % to take into account Dirichlet boundary
    elseif iface==nf+2
        nv = Np2d;
        ind_ff   = Np2d*neq+nf*Nfl*neq + (1:Np2d*neq);
        ind_fe   = Np2d*(Np1dTor-1)*neq + (1:Np2d*neq);
        xyfg     = N2d*Xe;                                % Gauss points coordinates of the face
        thetafg  = tel(end);                              % Gauss points coordinates of the face (theta component)
        Nfi      = N2d;                                   % Shape functions in the face
        J11 = N2dxi*Xe(:,1);                              % ng x 1
        J21 = N2dxi*Xe(:,2);                              % ng x 1
        J12 = N2deta*Xe(:,1);                             % ng x 1
        J22 = N2deta*Xe(:,2);                             % ng x 1
        detJ = J11.*J22 - J21.*J12;                       % determinant of the Jacobian
        ipw = ipw2d;
        ngausspol = ngauss2d;                             % number of Gauss points in the poloidal plane
        ngausstor = 1;                                    % number of Gauss points in the toroidal direction
        ngauss = ngausspol*ngausstor;
        dsurf = detJ.*ipw;                                % elemental surface
        n_g = repmat([0,0,1],ngauss,1);                   % outer normal
        rg = xyfg(:,1);                                   % radius
        isdir = 0;                                        % to take into account Dirichlet boundary
    else
        nodesv      = refElTor.faceNodes3(iface-1,:);
        nodes2      = refEl.faceNodes(iface-1,:);
        nv          = Np1dPol*Np1dTor;                            % Number of nodes in this face
        ind_ff      = Np2d*neq + (iface-2)*nv*neq+ (1:nv*neq);
        ind_fe      = colt(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'));
        xyfg        = N1dPol*Xe(nodes2,:);                        % Gauss points coordinates of the face (xy component)
        thetafg     = N1dTor*tel;                                 % Gauss points coordinates of the face (theta component)
        Nfi         = refElTor.sfTor;                             % Shape functions in the face
        ngausstor   = ngauss1dtor;                                % number of Gauss points in the poloidal plane
        ngausspol   = ngauss1dpol;                                % number of Gauss points in the toroidal direction
        rg = col(repmat(xyfg(:,1),ngausstor,1));                  % radius
        isdir = aux_dir(iface-1);                                 % to take into account Dirichlet boundary
        
        % computation of dsurf for the toroidal faces
        xyd_g         = N1dxPol*Xe(nodes2,:);
        xydNorm_g = sqrt(xyd_g(:,1).^2+xyd_g(:,2).^2);
        %         dsurf = col((ipw1dt*0.5*htheta)*(ipw1dp.*xydNorm_g)');
        dsurf = col((ipw1dp.*xydNorm_g)*(ipw1dt*0.5*htheta)');
        
        % computation of the outer normal
        t_g = xyd_g./xydNorm_g;
        n_g  = repmat([t_g(:,2) -t_g(:,1)],ngausstor,1);
        n_g(:,3) = 0;
    end
    
    
    if axisym && iface>1 && iface<nf+2
        dsurf=dsurf.*rg;
    end
    
    if isdir
        % Exact solution at face Gauss points
        sol_g = analyticalSolution_3D(xyfg,thetafg);
    else
        sol_f = sol_hat(ind_ff);
        if (iface>1 && iface<nf+2)
            if flipFace_e(iface-1)
                sol_f = sol_f(perm);
            end
        end
        sol_g = Nfi*transpose(reshape(sol_f,neq,nv));
    end
    
    % Initialization
    Hn11 = zeros(nv); Hn12 = Hn11; Hn13 = Hn11; Hn14 = Hn11;
    Hn21 = Hn11; Hn22 = Hn11; Hn23 = Hn11; Hn24 = Hn11;
    Hn31 = Hn11; Hn32 = Hn11; Hn33 = Hn11; Hn34 = Hn11;
    Hn41 = Hn11; Hn42 = Hn11; Hn43 = Hn11; Hn44 = Hn11;
    Hdirf = zeros(nv,neq);
    prova = 0;
    
    %     for igpol = 1:ngausspol
    %         for igtor = 1:ngausstor
    for igtor = 1:ngausstor
        for igpol = 1:ngausspol
            
            %             g = (igpol-1)*ngausstor+igtor;
            g = (igtor-1)*ngausspol+igpol;
            
            % Velocity shape functions and derivatives at the current integration point
            Nf_g = Nfi(g,:);
            NN = Nf_g'*Nf_g;
            
            % Contribution of the current integration point to the elemental matrix
            %             igauss = g;
            %             if (iface>1 && iface<nf+2)
            %                 if flipFace_e(iface-1)
            %                     igauss = ngauss_f-g+1;
            %                 end
            %             end
            An = jacobianMatricesFace(xyfg(igpol,1),xyfg(igpol,2),...
                thetafg(igtor),n_g(g,1),n_g(g,2),n_g(g,3),sol_g(g,:));
            
            
            Hn11 = Hn11 + NN*An(1,1)*dsurf(g);
            Hn12 = Hn12 + NN*An(1,2)*dsurf(g);
            Hn13 = Hn13 + NN*An(1,3)*dsurf(g);
            Hn14 = Hn14 + NN*An(1,4)*dsurf(g);
            Hn21 = Hn21 + NN*An(2,1)*dsurf(g);
            Hn22 = Hn22 + NN*An(2,2)*dsurf(g);
            Hn23 = Hn23 + NN*An(2,3)*dsurf(g);
            Hn24 = Hn24 + NN*An(2,4)*dsurf(g);
            Hn31 = Hn31 + NN*An(3,1)*dsurf(g);
            Hn32 = Hn32 + NN*An(3,2)*dsurf(g);
            Hn33 = Hn33 + NN*An(3,3)*dsurf(g);
            Hn34 = Hn34 + NN*An(3,4)*dsurf(g);
            Hn41 = Hn41 + NN*An(4,1)*dsurf(g);
            Hn42 = Hn42 + NN*An(4,2)*dsurf(g);
            Hn43 = Hn43 + NN*An(4,3)*dsurf(g);
            Hn44 = Hn44 + NN*An(4,4)*dsurf(g);
            Hdirf = Hdirf + Nf_g'*transpose(An*sol_g(g,:)')*dsurf(g);
            
            prova = prova + dsurf(g);
        end
    end
    
    
    
    
    % expand the matrices
    Hloc = expandMatrixCv(Hn11,Hn12,Hn13,Hn14,Hn21,Hn22,Hn23,Hn24,Hn31,Hn32,Hn33,Hn34,Hn41,Hn42,Hn43,Hn44);
    
    % elemental assembly
    He(ind_fe,ind_ff) = He(ind_fe,ind_ff) + ~isdir*Hloc;
    Hf(ind_ff,ind_ff) = Hf(ind_ff,ind_ff) + ~isdir*Hloc;
    Hdire(ind_fe) = Hdire(ind_fe) + isdir*col(Hdirf');
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

function [Jx,Jy,Jt] = jacobianMatrices(x,y,t,U,iElem,ig)

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
    b = defineMagneticField_3D([x,y],t);
end
%
%
%
if ~decoupleEquations
    J = [                         0,                                                                                             1,                                        0,            0; ...
        -2/3*U(2)^2/U(1)^2                                                                         4/3*U(2)/U(1)                                      2/3,        2/3; ...
        -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2                 5/3*U(2)/U(1),    0 ;   ...
        -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                                     0,       5/3*U(2)/U(1)];
%     Jx = b(1)*[                         0,                                                                                             1,                                        0,            0; ...
%         -2/3*U(2)^2/U(1)^2                                                                         4/3*U(2)/U(1)                                      2/3,        2/3; ...
%         -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2                 5/3*U(2)/U(1),    0 ;   ...
%         -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                                     0,       5/3*U(2)/U(1)];
%     
%     
%     Jy = b(2)*[                         0,                                                                                             1,                                        0,            0; ...
%         -2/3*U(2)^2/U(1)^2                                                                         4/3*U(2)/U(1)                                      2/3,        2/3; ...
%         -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2                  5/3*U(2)/U(1),    0 ;   ...
%         -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                                     0,       5/3*U(2)/U(1)];
else
    J = [                         0,                                                                                             1,                                        0,            0; ...
        (-1*Gamma^2/rho^2+Mref)                                                      2*Gamma/rho                                     0     ,        0; ...
        -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                5/3*U(3)/U(1)-U(2)^2/U(1)^2                 5/3*U(2)/U(1),    0 ;   ...
        -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                         0,       5/3*U(2)/U(1)];
%     Jx = b(1)*[                         0,                                                                                             1,                                        0,            0; ...
%         (-1*Gamma^2/rho^2+Mref)                                                      2*Gamma/rho                                     0     ,        0; ...
%         -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                5/3*U(3)/U(1)-U(2)^2/U(1)^2                 5/3*U(2)/U(1),    0 ;   ...
%         -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                         0,       5/3*U(2)/U(1)];
%     
%     
%     Jy = b(2)*[                         0,                                                                                             1,                                        0,            0; ...
%         (-1*Gamma^2/rho^2+Mref)                                                         2*Gamma/rho                                     0     ,        0; ...
%         -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2            5/3*U(2)/U(1),    0 ;   ...
%         -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                            0,       5/3*U(2)/U(1)];
    
end

Jx = b(1)*J;
Jy = b(2)*J;
Jt = b(3)*J;

function Jn = jacobianMatricesFace(x,y,t,nx,ny,nt,U,iFace,ig)

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
    b = defineMagneticField_3D([x,y],t);
end
bn = b(1)*nx+b(2)*ny+b(3)*nt;
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









function  J= jacobianMatricesNOMAGFIELD(U)

global Magnetic testcase useThreshold  decoupleEquations Mref
rho = U(1);
Gamma = U(2);
if useThreshold
    if U(1)<useThreshold,
        rho = useThreshold;
        Gamma = 0;
    end
end

if ~decoupleEquations
    J = [                         0,                                                                                             1,                                        0,            0; ...
        -2/3*U(2)^2/U(1)^2                                                                         4/3*U(2)/U(1)                                      2/3,        2/3; ...
        -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                      5/3*U(3)/U(1)-U(2)^2/U(1)^2                 5/3*U(2)/U(1),    0 ;   ...
        -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                                     0,       5/3*U(2)/U(1)];
else
    J = [                         0,                                                                                             1,                                        0,            0; ...
        (-1*Gamma^2/rho^2+Mref)                                                      2*Gamma/rho                                     0     ,        0; ...
        -5/3*U(2)*U(3)/U(1)^2+2/3*U(2)^3/U(1)^3                5/3*U(3)/U(1)-U(2)^2/U(1)^2                 5/3*U(2)/U(1),    0 ;   ...
        -5/3*U(4)*U(2)/U(1)^2,                                                       5/3*U(4)/U(1),                                         0,       5/3*U(2)/U(1)];    
end








function G = GimplicitMatrix(x,y,t,U,iElem,ig)

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
    [~,db] = defineMagneticField_3D([x,y],t);
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