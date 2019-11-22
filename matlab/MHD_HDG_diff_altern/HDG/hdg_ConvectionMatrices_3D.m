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
        H(:,:,iElem) = He;if axisym
        Hf(:,:,iElem) = Hfe;
        Hdir(:,iElem) = Hdire;
    end
end

%% Elemental matrices
function [Ce,He,Hdire,Hf] = elementalMatrices(Xe,tel,refEl,sol,sol_hat,aux_dir,flipFace_e,iElem,Fe)

global neq axisym refElTor kT

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
Cx11 = zeros(Nv); Cx12 = Cx11;
Cx21 = Cx11; Cx22 = Cx11;

Cy11 = zeros(Nv); Cy12 = Cy11;
Cy21 = Cy11; Cy22 = Cy11;

Ct11 = zeros(Nv); Ct12 = Ct11;
Ct21 = Ct11; Ct22 = Ct11;

Ce = zeros(Nv*neq);
% reshape solution
sol = reshape(sol,neq,numel(sol)/neq)';
ind_ass = 0:neq:neq*(Nv-1);%(/ (i, i = 0, Neq*(Np-1),Neq) /)
prova = 0;
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
        sol_g = Ni*sol; % 1 x nvar
        
        % Integration weight
        dvolu=ipw2d(igpol)*det(J)*dvolu1d;
        
        if axisym
            dvolu=dvolu*xg;
        end
        
        
        
        
        
        
        
        % Jacobian matrices
%         [Ax,Ay,At] = jacobianMatrices(xg,yg,tg,sol_g,iElem,g);
        
        % Contribution of Nythe current integration point to the elemental matrix
%         Cx11 = Cx11 + Nx'*Ni*Ax(1,1)*dvolu;
%         Cx12 = Cx12 + Nx'*Ni*Ax(1,2)*dvolu;
%         Cx21 = Cx21 + Nx'*Ni*Ax(2,1)*dvolu;
%         Cx22 = Cx22 + Nx'*Ni*Ax(2,2)*dvolu;
%         
%         Cy11 = Cy11 + Ny'*Ni*Ay(1,1)*dvolu;
%         Cy12 = Cy12 + Ny'*Ni*Ay(1,2)*dvolu;
%         Cy21 = Cy21 + Ny'*Ni*Ay(2,1)*dvolu;
%         Cy22 = Cy22 + Ny'*Ni*Ay(2,2)*dvolu;
%         
%         Ct11 = Ct11 + Nt'*Ni*At(1,1)*dvolu;
%         Ct12 = Ct12 + Nt'*Ni*At(1,2)*dvolu;
%         Ct21 = Ct21 + Nt'*Ni*At(2,1)*dvolu;
%         Ct22 = Ct22 + Nt'*Ni*At(2,2)*dvolu;
        
        
        
        
            b = defineMagneticField3d([xg,yg],tg);
            A = [0, 1; (-1*sol_g(2)^2/sol_g(1)^2+kT), 2*sol_g(2)/sol_g(1)];
            NxNi = Nx'*Ni*dvolu;
            NyNi = Ny'*Ni*dvolu;
            NtNi = Nt'*Ni*dvolu;
            
            NNxy = b(1)*NxNi+b(2)*NyNi+b(3)*NtNi;
              for i=1:neq
                 for j=1:neq
                    Ce(i+ind_ass,j+ind_ass) = Ce(i+ind_ass,j+ind_ass) + A(i,j)*NNxy;
                 end
              end        
        
        
        
        prova = prova + dvolu;
    end
end
% expand the matrices
% Cx = expandMatrixCv(Cx11,Cx12,Cx21,Cx22);
% Cy = expandMatrixCv(Cy11,Cy12,Cy21,Cy22);
% Ct = expandMatrixCv(Ct11,Ct12,Ct21,Ct22);
% Ce = Cx+Cy+Ct;


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
        sol_g = analyticalSolution3d(xyfg,thetafg);
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
    Hn11 = zeros(nv); Hn12 = Hn11;
    Hn21 = Hn11; Hn22 = Hn11; 
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
            Hn21 = Hn21 + NN*An(2,1)*dsurf(g);
            Hn22 = Hn22 + NN*An(2,2)*dsurf(g);
            
            Hdirf = Hdirf + Nf_g'*transpose(An*sol_g(g,:)')*dsurf(g);
%             auxprova = (0.5-xyfg(igpol,1))/0.5;
%            if (xyfg(igpol,1)==max(xyfg(:,1)) && thetafg(igtor)==min(thetafg))
%                auxprova = 1;
%            else
%                auxprova = 0;
%            end
            prova = prova + dsurf(g);
        end
    end
    
    

    
    % expand the matrices
    Hloc = expandMatrixCv(Hn11,Hn12,Hn21,Hn22);
    
    % elemental assembly
    He(ind_fe,ind_ff) = He(ind_fe,ind_ff) + ~isdir*Hloc;
    Hf(ind_ff,ind_ff) = Hf(ind_ff,ind_ff) + ~isdir*Hloc;
    Hdire(ind_fe) = Hdire(ind_fe) + isdir*col(Hdirf');
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

function [Jx,Jy,Jt] = jacobianMatrices(x,y,t,U,iElem,ig)

global kT Magnetic testcase useThreshold
rho = U(1);
Gamma = U(2);
if useThreshold
    if U(1)<useThreshold
        rho = useThreshold;
        Gamma = 0;
        
    end
end
if (testcase.n >= 50 && testcase.n<60)
    b = [Magnetic.bx(ig,iElem), Magnetic.by(ig,iElem)];
else
    b = defineMagneticField3d([x,y],t);
end
Jx = [0, b(1); (-1*Gamma^2/rho^2+kT)*b(1), 2*Gamma/rho*b(1)];
Jy = [0, b(2); (-1*Gamma^2/rho^2+kT)*b(2), 2*Gamma/rho*b(2)];
Jt = [0, b(3); (-1*Gamma^2/rho^2+kT)*b(3), 2*Gamma/rho*b(3)];



% function Jn = jacobianMatricesFace(x,y,t,nx,ny,nt,U,iFace,ig)
% 
% global kT  Magnetic testcase useThreshold
% rho = U(1);
% Gamma = U(2);
% 
% if useThreshold
%     if U(1)<useThreshold
%         rho = useThreshold;
%         Gamma = 0;
%         
%     end
% end
% if (testcase.n >= 50 && testcase.n<60)
%     b = [Magnetic.bxfaces(ig,iFace), Magnetic.byfaces(ig,iFace)];
% else
%     b = defineMagneticField3d([x,y],t);
% end
% bn = b(1)*nx+b(2)*ny+b(3)*nt;
% Jn = [0, bn; (-1*Gamma^2/rho^2+kT)*bn, 2*Gamma/rho*bn];



function Jn = jacobianMatricesFace(x,y,t,nx,ny,nt,U)

global kT  Magnetic testcase useThreshold
rho = U(1);
Gamma = U(2);

if useThreshold
    if U(1)<useThreshold
        rho = useThreshold;
        Gamma = 0;
        
    end
end
% if (testcase.n >= 50 && testcase.n<60)
%     b = [Magnetic.bxfaces(ig,iFace), Magnetic.byfaces(ig,iFace)];
% else
    b = defineMagneticField3d([x,y],t);
% end
bn = b(1)*nx+b(2)*ny+b(3)*nt;
Jn = [0, bn; (-1*Gamma^2/rho^2+kT)*bn, 2*Gamma/rho*bn];
