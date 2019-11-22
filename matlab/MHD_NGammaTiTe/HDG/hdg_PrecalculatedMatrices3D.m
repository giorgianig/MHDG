function [A,B,C,C_dir,D,E,Df,Ef,Edir,invL,force,Lf,...
    L,P,Pb,Q,Qb,Qf] =...
    hdg_PrecalculatedMatrices3D(X,T,flipFace,refEl,tau,Fcon,F_dir)

% mesh data
global ntor neq refElTor theta             % number of equations (2 for isothermal)
nf  = size(refEl.faceNodes,1);    % number of faces in one element
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Np2d      = size(refEl.NodesCoord,1);
Nv       = Np1dTor*Np2d;
N2d      = size(T,1);                    % number of elements
Ne       = N2d*ntor;
Nfl   = Np1dPol*Np1dTor;
Nfp  = Np2d*2+nf*Np1dPol*Np1dTor;
% nv = size(refEl.NodesCoord1d,1);   % number of face nodes for the velocity
Ndim = 3;
tdiv = linspace(0,theta,ntor+1);

% allocation and initialization
A  = zeros(neq*Nv,neq*Nv,Ne);
B = zeros(neq*Ndim*Nv,neq*Nv,Ne);
C = zeros(neq*Ndim*Nv,neq*Nfp,Ne);
C_dir = zeros(neq*Ndim*Nv,Ne);
D  = A;
E  = zeros(neq*Nv,neq*Nfp,Ne);
Df = zeros(neq*Nfp,neq*Nv,Ne);
Ef = zeros(neq*Nfp,neq*Nfp,Ne);
Edir = zeros(neq*Nv,Ne);
% PQinvL = zeros(neq*Nv,Ndim*neq*Nv,Ne);
invL  = zeros(neq*Ndim*Nv,neq*Ndim*Nv,Ne);
force = zeros(neq*Nv,Ne);
Lf = zeros(neq*Nfp,neq*Ndim*Nv,Ne);

L = zeros(neq*Ndim*Nv,neq*Ndim*Nv,Ne);
P = zeros(neq*Nv,Ndim*neq*Nv,Ne);
Pb = zeros(neq*Nv,Ndim*neq*Nv,Ne);
Q = zeros(neq*Nv,Ndim*neq*Nv,Ne);
Qb = zeros(neq*Nv,Ndim*neq*Nv,Ne);
Qf = Lf;


% compute permutation for flipping faces
perm =  setperm3d(Np1dPol,Np1dTor,neq);

%
% ind_1_v_L = (1:2*nv);
% ind_2_v_L = 2*nv + (1:2*nv);
% ind_3_v_L = 4*nv + (1:2*nv);

% loop in elements
for itor = 1:ntor
    tel = tdiv(itor)+0.5*(refElTor.NodesCoord1d+1)*(tdiv(itor+1)-tdiv(itor));
    for iel = 1:N2d
        
        iElem = (itor-1)*N2d + iel;
        Fe = Fcon(iel,:);
        Te = T(iel,:);
        Xe = X(Te,:);
        flipFace_e = flipFace(iel,:);
        aux_dir = F_dir(iel,:);
        
        
        % elemental matrices
        [Ae,Be,Ce,C_dire,De,Ee,Dfe,Efe,Edire,invLe,Lfe,fe,...
            Le,Pe,Qe,Pbe,Qbe,Qfe] = ...
            elementalMatrices(Xe,tel,refEl,tau,aux_dir,iel,Fe,flipFace_e);
        
        for iface = 1:nf
            
            if flipFace_e(iface)
                ind_v_L = Np2d*neq + (iface-1)*Nfl*neq+ (1:Nfl*neq);
                Dfe(ind_v_L,:) = Dfe(ind_v_L(perm),:);
                Lfe(ind_v_L,:) = Lfe(ind_v_L(perm),:);
                Qfe(ind_v_L,:) = Qfe(ind_v_L(perm),:);
                Efe(ind_v_L,:) = Efe(ind_v_L(perm),:);
                Efe(:,ind_v_L) = Efe(:,ind_v_L(perm));
            end
        end
        
        % store matrices
        A(:,:,iElem) = Ae;
        B(:,:,iElem) = Be;
        C(:,:,iElem) = Ce;
        C_dir(:,iElem) = C_dire;
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
end

%% Elemental matrices
function [A,B,C,C_dir,D,E,Df,Ef,Edir,invL,Lf,f,...
    L,P,Q,Pb,Qb,Qf] = ...
    elementalMatrices(Xe,tel,refEl,tau,aux_dir,iel,Fe,flipFace)

global kT neq diff_n diff_u diff_ei diff_ee Magnetic testcase axisym driftvel refElTor

% mesh data
htheta   = tel(end)-tel(1);

Ndim = 3;
nf = size(refEl.faceNodes,1);
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);     % Number of 1d points for toroidal line
Np2d      = size(refEl.NodesCoord,1);
Nv = Np1dTor*Np2d;
Nfl   = Np1dPol*Np1dTor;
Nfp  = Np2d*2+nf*Np1dPol*Np1dTor;
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);

% initialize all the matrices
C = zeros(neq*Ndim*Nv,neq*Nfp);
C_dir = zeros(neq*Ndim*Nv,1);
Lf_transp = C;
Lf = Lf_transp';
D = zeros(neq*Nv,neq*Nv);
Df = zeros(neq*Nfp,neq*Nv);
E = zeros(neq*Nv,neq*Nfp);
Edir = zeros(neq*Nv,1);
Ef = zeros(neq*Nfp,neq*Nfp);
Q = zeros(neq*Nv,neq*Ndim*Nv);
Qf = Lf_transp';
Qb = Q;
Qb_aux = Q;
Mass = zeros(Nv);
fe = zeros(Nv,neq);
Bx = Mass;
By = Mass;
Bt = Mass;
Px = Mass;
Py = Mass;
Pt = Mass;
Pbx = Bx;
Pby = By;
Pbt = By;
if driftvel
    Dve = Mass;
end

% Information of the reference element for the velocity
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

prova = 0;
elarea = 0;
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
        
        %         Niv_g = Niv(g,:);
        %         Nxiv_g = Nxiv(g,:);
        %         Netav_g = Netav(g,:);
        
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
            Nx_ax_g = Nx_g+1/xg*N2_g;
        else
            Nx_ax_g = Nx_g;
        end
        if axisym
            Nt_ax_g = 1/xg*N1x_g;
        else
            Nt_ax_g = N1x_g;
        end
        
        %% 3d shape functions
        Ni      = colt(N2_g'*N1_g);
        Nx      = colt(Nx_g'*N1_g);
        Nxx     = colt(Nx_ax_g'*N1_g);
        Ny      = colt(Ny_g'*N1_g);
        Nt      = colt(N2_g'*Nt_ax_g);
        
        
        % Integration weight
        dvolu=ipw2d(igpol)*det(J)*dvolu1d;
        
        if axisym
            dvolu=dvolu*xg;
        end
        
        % body force at integration point
        force = [0,0];
        switch testcase.n
            
            case 51 % Source between -0.88 and -0.90
                if (Magnetic.flux2D(g,iel)<=-0.88 && Magnetic.flux2D(g,iel)>=-0.90)
                    force(1) = 3.20119388718018e-05;
                end
            case 52 % Source between -0.90 and -1
                if (Magnetic.flux2D(g,iel)<=-0.90 && Magnetic.flux2D(g,iel)>=-1)
                    force(1) = 9.45155008295538e-06;
                end
            case 53 % Source from -0.90
                if Magnetic.flux2D(g,iel)<=-0.90
                    force(1) = 7.24032211339971e-06;
                end
            case 54 % Source betwenn -0.98 and -1
                if (Magnetic.flux2D(g,iel)<=-0.98  && Magnetic.flux2D(g,iel)>=-1)
                    force(1) = 5.78723047118297e-05;
                end
            case 55 % Source from -1.03
                if Magnetic.flux2D(g,iel)<=-1.03
                    force(1) = 0.000115575293741846;
                end
            otherwise
                force = bodyForce_3D([xg,yg],tg);
        end
        
 
%  display(['xy: ', num2str(xg),'  ',num2str(yg)])       
%  display(['tg: ', num2str(tg)])
%  display(['force: ', num2str(force(1)),'  ',num2str(force(2))])       
        
        % magnetic field
        if (testcase.n >= 50 && testcase.n<60)
            b = [Magnetic.bx(g,iel), Magnetic.by(g,iel)];
            divb = Magnetic.div(g,iel);
            if driftvel
                dv = [Magnetic.dvx(g,iel),Magnetic.dvy(g,iel)];
            end
        else
            if driftvel
                [b,db,dvx,dvy] = defineMagneticField([xg,yg]);
                dv = [dvx,dvy];
            else
                [b,db] = defineMagneticField_3D([xg,yg],tg);
            end
            divb = sum(db);
        end
        
        
        % Contribution of the current integration point to the elemental matrix
        Mass = Mass + Ni'*Ni*dvolu;
        Bx = Bx + Nxx'*Ni*dvolu;
        By = By + Ny'*Ni*dvolu;
        Bt = Bt + Nt'*Ni*dvolu;
        Px = Px + Nx'*Ni*dvolu;
        Py = Py + Ny'*Ni*dvolu;
        Pt = Pt + Nt'*Ni*dvolu;
        Pbx = Pbx + b(1)*(b(1)*Nx'+b(2)*Ny'+b(3)*Nt')*Ni*dvolu;
        Pby = Pby + b(2)*(b(1)*Nx'+b(2)*Ny'+b(3)*Nt')*Ni*dvolu;
        Pbt = Pbt + b(3)*(b(1)*Nx'+b(2)*Ny'+b(3)*Nt')*Ni*dvolu;
        fe = fe + Ni'*force*dvolu;
        prova = prova + Ni'*dvolu;
        elarea = elarea + dvolu;
        
        
        
        
    end
end
% stop

% expand the matrices
L = expandMatrixA(Mass,neq*Ndim);
B = expandMatrixB(Bx,By,Bt);
P = transpose(expandMatrixB(Px',Py',Pt'));
Pb = transpose(expandMatrixB(Pbx',Pby',Pbt'));
A = expandMatrixA(Mass,neq);
f = col(fe');


%********************************************************
%% FACES COMPUTATIONS:
for iface = 1:(nf+2)
    % face nodes
    if iface==1
        nv = Np2d;
        ind_ff = 1:Np2d*neq;                                 % assembly face-face
        ind_fe = 1:Np2d*neq;                                 % assembly face-element (for var)
        ind_fq = 1:Ndim*Np2d*neq;                        % assembly face-element (for grad)
        xyfg     = N2d*Xe;                                 % Gauss points coordinates of the face (x and y components)
        thetafg  = tel(1);                                % Gauss points coordinates of the face (theta component)
        Nfi      = N2d;                                      % Shape functions in the face
        J11 = N2dxi*Xe(:,1);                                   % ng x 1
        J21 = N2dxi*Xe(:,2);                                   % ng x 1
        J12 = N2deta*Xe(:,1);                                 % ng x 1
        J22 = N2deta*Xe(:,2);                                 % ng x 1
        detJ = J11.*J22 - J21.*J12;                           % determinant of the Jacobian
        ngausspol = ngauss2d;                               % number of Gauss points in the poloidal plane
        ngausstor = 1;                                            % number of Gauss points in the toroidal direction
        ngauss = ngausspol*ngausstor;
        ipw = ipw2d;
        dsurf = detJ.*ipw;                                       % elemental surface
        n_g = repmat([0,0,-1],ngauss,1);                     % outer normal
        rg = xyfg(:,1);                                             % radius
        isdir = 0;                                       % to take into account Dirichlet boundary
    elseif iface==nf+2
        nv = Np2d;
        ind_ff   = Np2d*neq+nf*Nfl*neq + (1:Np2d*neq);
        ind_fe   = Np2d*(Np1dTor-1)*neq + (1:Np2d*neq);
        ind_fq   = Np2d*(Np1dTor-1)*Ndim*neq +(1:Ndim*Np2d*neq);
        xyfg     = N2d*Xe;                                % Gauss points coordinates of the face
        thetafg  = tel(end);                            % Gauss points coordinates of the face (theta component)
        Nfi      = N2d;                                      % Shape functions in the face
        J11 = N2dxi*Xe(:,1);                                   % ng x 1
        J21 = N2dxi*Xe(:,2);                                   % ng x 1
        J12 = N2deta*Xe(:,1);                                 % ng x 1
        J22 = N2deta*Xe(:,2);                                 % ng x 1
        detJ = J11.*J22 - J21.*J12;                           % determinant of the Jacobian
        ipw = ipw2d;
        ngausspol = ngauss2d;                               % number of Gauss points in the poloidal plane
        ngausstor = 1;                                            % number of Gauss points in the toroidal direction
        ngauss = ngausspol*ngausstor;
        dsurf = detJ.*ipw;                                        % elemental surface
        n_g = repmat([0,0,1],ngauss,1);                    % outer normal
        rg = xyfg(:,1);                                             % radius
        isdir = 0;                                       % to take into account Dirichlet boundary
    else
        nodesv      = refElTor.faceNodes3(iface-1,:);
        nodes2      = refEl.faceNodes(iface-1,:);
        nv          = Np1dPol*Np1dTor;                 % Number of nodes in this face
        ind_ff      = Np2d*neq + (iface-2)*nv*neq+ (1:nv*neq);
        ind_fe      = col(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'));
        ind_fq      = col(bsxfun(@plus,(nodesv-1)*Ndim*neq,(1:Ndim*neq)'));
        xyfg        = N1dPol*Xe(nodes2,:);            % Gauss points coordinates of the face (xy component)
        thetafg     = N1dTor*tel;                       % Gauss points coordinates of the face (theta component)
        Nfi         = refElTor.sfTor;                        % Shape functions in the face
        ngausstor   = ngauss1dtor;                              % number of Gauss points in the poloidal plane
        ngausspol   = ngauss1dpol;                              % number of Gauss points in the toroidal direction
        rg = col(repmat(xyfg(:,1),ngausstor,1));       % radius
        isdir = aux_dir(iface-1);                               % to take into account Dirichlet boundary
        
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
    
    
    % Exact solution at face Gauss points
    u_ex = analyticalSolution_3D(xyfg,thetafg);
    
    
     if (testcase.n >= 50 && testcase.n<60)
        if flipFace(iface)
            igauss = ngauss_f-g+1;
        else
            igauss = g;
        end
        b = [Magnetic.bxfaces(igauss,iFace), Magnetic.byfaces(igauss,iFace)];
    else
        b = defineMagneticField_3D([xyfg(:,1),xyfg(:,2)],thetafg(:));
    end   
    
    
    
    %     % face nodes
    %     nodesv = faceNodesv(iface,:);
    %     iFace = Fe(iface);
    %
    %     % indices for local assembly
    %     ind_face_2 = (iface-1)*neq*nv + (1:neq*nv);
    %     ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
    %     ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*Ndim,(1:neq*Ndim)'),neq*Ndim*nv,1); % assembly face to elem for velocity gradient
    %
    %     xf = xe(nodesv);
    %     yf = ye(nodesv);
    %
    %     xyfg = N1dv*[xf yf];
    %     u_ex = analyticalSolution(xyfg);
    %     tau_f = tau;
    
    
    
    
    
    
    
    
    
    
    % initialize local matrices
    Massf = zeros(nv);
    E_dirf = zeros(nv,neq);
    Cnx = zeros(nv);
    Cny = Cnx;
    Cnt = Cnx;
    Qbx = Cnx;
    Qby = Cnx;
    Qbt = Cnx;
    Qbx_aux = Cnx; Qby_aux = Cnx; Qbt_aux = Cnx;
     Cloc_dir = zeros(neq*Ndim*nv,1);   
%     Cu1nx_dir = zeros(nv,1);
%     Cu1ny_dir  = Cu1nx_dir ;
%     Cu1nt_dir  = Cu1nx_dir ;
%     Cu2nx_dir  = Cu1nx_dir ;
%     Cu2ny_dir  = Cu1nx_dir ;
%     Cu2nt_dir  = Cu1nt_dir ;
    
    %  LOOP IN GAUSS POINTS
%     for igpol = 1:ngausspol
%         for igtor = 1:ngausstor

%         g = (igpol-1)*ngausstor+igtor;

    for igtor = 1:ngausstor
        for igpol = 1:ngausspol 
            
%             
            g = (igtor-1)*ngausspol+igpol;
            
            % Velocity shape functions and derivatives at the current integration point
            Nf_g = Nfi(g,:);
             
           tau_f = tau;
            uexn = col(n_g(g,:)'*u_ex(g,:));
            
            % Contribution of the current integration point to the elemental matrix
            Cnx = Cnx + Nf_g'*Nf_g*n_g(g,1)*dsurf(g);
            Cny = Cny + Nf_g'*Nf_g*n_g(g,2)*dsurf(g);
            Cnt = Cnt + Nf_g'*Nf_g*n_g(g,3)*dsurf(g);
            
%             Cu1nx_dir = Cu1nx_dir + Nf_g'*u_ex(g,1)*n_g(1)*dsurf(g);
%             Cu1ny_dir  = Cu1ny_dir  + Nf_g'*u_ex(g,1)*n_g(2)*dsurf(g);
%             Cu1nt_dir  = Cu1nt_dir   + Nf_g'*u_ex(g,1)*n_g(3)*dsurf(g);
%             
%             Cu2nx_dir  = Cu2nx_dir  + Nf_g'*u_ex(g,2)*n_g(1)*dsurf(g);
%             Cu2ny_dir  = Cu2ny_dir  + Nf_g'*u_ex(g,2)*n_g(2)*dsurf(g);
%             Cu2nt_dir  = Cu2nt_dir   + Nf_g'*u_ex(g,2)*n_g(3)*dsurf(g);
%             check this!!
            Cloc_dir = Cloc_dir + col(uexn*Nf_g)*dsurf(g);
            
            Qbx = Qbx + b(g,1)*Nf_g'*(b(g,1)*Nf_g*n_g(g,1)+b(g,2)*Nf_g*n_g(g,2)+b(g,3)*Nf_g*n_g(g,3))*dsurf(g);
            Qby = Qby + b(g,2)*Nf_g'*(b(g,1)*Nf_g*n_g(g,1)+b(g,2)*Nf_g*n_g(g,2)+b(g,3)*Nf_g*n_g(g,3))*dsurf(g);
            Qbt = Qbt + b(g,3)*Nf_g'*(b(g,1)*Nf_g*n_g(g,1)+b(g,2)*Nf_g*n_g(g,2)+b(g,3)*Nf_g*n_g(g,3))*dsurf(g);
            E_dirf = E_dirf + tau_f*(Nf_g'*u_ex(g,:))*dsurf(g);
            Massf = Massf + tau_f *(Nf_g'*Nf_g)*dsurf(g);
        end
    end
    
    % expand the matrices
    C_loc = expandMatrixB(Cnx,Cny,Cnt);
    Qb_loc = transpose( expandMatrixB(Qbx',Qby',Qbt') );
    Qb_loc_aux = transpose( expandMatrixB(Qbx_aux',Qby_aux',Qbt_aux') );
%     C_loc_dir = expandMatrixF(expandMatrixF(Cu1nx_dir ,Cu2nx_dir ),...
%         expandMatrixF(Cu1ny_dir ,Cu2ny_dir )); % check this
    D_loc = expandMatrixA(Massf,neq);
    E_loc = expandMatrixA(Massf,neq);
    
%     Massf
    
    % elemental assembly
    C(ind_fq,ind_ff) = C(ind_fq,ind_ff) + ~isdir*C_loc;
    C_dir(ind_fq) = C_dir(ind_fq) + isdir*Cloc_dir;
    %     Lf_transp(ind_fq,ind_ff) = Lf_transp(ind_fq,ind_ff) + C_loc;
    Lf(ind_ff,ind_fq) = Lf(ind_ff,ind_fq) +  ~isdir*transpose( expandMatrixB(Cnx',Cny',Cnt') );
    Q(ind_fe,ind_fq) = Q(ind_fe,ind_fq) + C_loc';
    D(ind_fe,ind_fe) = D(ind_fe,ind_fe) + D_loc;
    Df(ind_ff,ind_fe) = Df(ind_ff,ind_fe) + ~isdir*D_loc;
    E(ind_fe,ind_ff) = E(ind_fe,ind_ff) + ~isdir*E_loc;
    Edir(ind_fe) = Edir(ind_fe) + isdir*col(E_dirf');
    Ef(ind_ff,ind_ff) = Ef(ind_ff,ind_ff) + ~isdir*E_loc;
    Qb(ind_fe,ind_fq) = Qb(ind_fe,ind_fq) + Qb_loc;
    Qf(ind_ff,ind_fq) = Qf(ind_ff,ind_fq) + ~isdir*Qb_loc;
    Qb_aux(ind_fe,ind_fq) = Qb_aux(ind_fe,ind_fq) + Qb_loc_aux;
end



% stop




diffdiag = repmat(col(repmat([diff_n diff_u diff_ei diff_ee],Ndim,1)),Nv,1);
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

% function res = expandMatrixB(Bx,By,Bt)
% % expand matrix B
% %   [ Bx  0
% %     By  0
% %     Bt   0
% %     0  Bx
% %     0  By
% %     0  Bt ]
% res = zeros([size(Bx) 6 2]);
% res(:,:,[1 2 3 10 11 12]) = cat(3,Bx,By,Bt,Bx,By,Bt);
% res = permute(res, [3 1 4 2]);
% res = reshape(res, 6*size(Bx,1),2*size(Bx,2));

function res = expandMatrixB(Bx,By,Bt)
% expand matrix B
%   [ Bx   0   0   0 
%     By   0   0   0
%     Bt   0   0   0
%     0   Bx   0   0
%     0   By   0   0  
%     0   Bt   0   0 
%     0    0  Bx   0 
%     0    0  By   0 
%     0    0  Bt   0
%     0    0   0  Bx 
%     0    0   0  By
%     0    0   0  Bt]

res = zeros([size(Bx) 12 4]);
res(:,:,[1 2 3 16 17 18 31 32 33 46 47 48]) = cat(3,Bx,By,Bt,Bx,By,Bt,Bx,By,Bt,Bx,By,Bt);
res = permute(res, [3 1 4 2]);
res = reshape(res, 12*size(Bx,1),4*size(Bx,2));