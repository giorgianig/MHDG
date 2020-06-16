function [Df,Ef,Hf,Lf,Qf,fH] = hdg_bc3d(Df,Ef,Hf,Lf,Qf,X,T,infoFaces,refEl,bn,F,F_dir,sol_hat)
                                     
global ntor neq testcase refElTor theta
%
Ne = size(Df,3);
nf = size(F,2);
N2d = size(F,1);
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Np2d    = size(refEl.NodesCoord,1);
Nfp     = Np2d*2+nf*Np1dPol*Np1dTor;
Nfl   = Np1dPol*Np1dTor;
tdiv = linspace(0,theta,ntor+1);
Nf        = max(max(F));                                              % Number of faces in the 2d plane
nDirFaces =  sum(sum(F_dir));

ind_dim = neq*[Np2d,repmat(Np1dPol*Np1dTor,1,nf),Np2d];
ind_sta = [1, 1+cumsum(ind_dim)];
ind_uhat = zeros(Nfp*neq,1);

% number of boundaries
nb = numel(bn);
fH = zeros(Nfp*neq,Ne);
for ibound = 1:nb
    
    name = bn{ibound}(4:end);
    
    % skip Dirichlet faces
    if testcase.n<50
        if any(strcmp(name,{'Diriclet','IN','OUT'})), continue ,end
    else
        if any(strcmp(name,{'Diriclet'})), continue ,end
    end
    faces = infoFaces.(['exteriorFaces_' name]);
    for itor = 1:ntor
        tel = tdiv(itor)+0.5*(refElTor.NodesCoord1d+1)*(tdiv(itor+1)-tdiv(itor));
        for iF = 1:size(faces,1)
            
            iel = faces(iF,1);
            ifac = faces(iF,2);
            Fe = F(iel,ifac);
            Te = T(iel,:);
            Xe = X(Te,:);
            
            iElem = (itor-1)*N2d+iel;
            Dfe = Df(:,:,iElem);
            Efe = Ef(:,:,iElem);
            Hfe = Hf(:,:,iElem);
            Lfe = Lf(:,:,iElem);
            Qfe = Qf(:,:,iElem);
            
            
            delta = 1+(itor-1)*(N2d*Np2d+(Nf-nDirFaces)*Nfl)*neq+neq*(N2d*Np2d+(Fe-1)*Nfl);
            ind_uhat = delta+(0: (Nfl*neq-1));
            solf = sol_hat(ind_uhat);
            fHe = fH(:,iElem);
            
            
            
            if any(strcmp(name,{'UP','DOWN','LEFT','RIGHT','WALL','ULIM'}))
                stop
                %                 [Dfe,Efe,Hfe,Lfe,Qfe] = elementalMatricesWall(Dfe,Efe,Hfe,Lfe,Qfe,Xe,tel,refEl,ifac);
            elseif any(strcmp(name,{'INLET','OUTLET'}))
                stop
                %                 [Dfe,Efe,Hfe,Qfe,fHe] = elementalMatricesInlOut(Dfe,Efe,Hfe,Qfe,Xe,tel,refEl,ifac,solf);
            elseif any(strcmp(name,{'Neumann'}))
                stop
                %                 [Dfe,Efe,Hfe,Qfe,fHe] = elementalMatricesNeumann(Dfe,Efe,Hfe,Qfe,Xe,tel,refEl,ifac);
            elseif any(strcmp(name,{'IN'}))
                [Dfe,Efe,Hfe,Lfe,Qfe,fHe] = elementalMatricesCore(Dfe,Efe,Hfe,Lfe,Qfe,Xe,tel,refEl,ifac,Fe);
            elseif any(strcmp(name,{'Dweak'}))
                [Dfe,Efe,Hfe,Lfe,Qfe,fHe] = elementalMatricesDirichletWeak(Dfe,Efe,Hfe,Lfe,Qfe,Xe,tel,refEl,ifac,Fe,iel);
                
            elseif any(strcmp(name,{'LIM','OUT'}))
                [Dfe,Efe,Hfe,Lfe,Qfe,fHe] = elementalMatricesBohm(Dfe,Efe,Hfe,Lfe,Qfe,Xe,tel,refEl,ifac,solf,Fe);
                %         elseif any(strcmp(name,{'Diriclet'}))
                %             [Dfe,Efe,Hfe,Lfe,Qfe,fHe] = elementalMatricesDirichlet(Dfe,Efe,Hfe,Lfe,Qfe,Xe,refEl,ifac,Fe);
            else
                error('Something wrong')
            end
            Df(:,:,iElem) = Dfe;
            Ef(:,:,iElem) = Efe;
            Hf(:,:,iElem) = Hfe;
            Lf(:,:,iElem) = Lfe;
            Qf(:,:,iElem) = Qfe;
            fH(:,iElem)   = fHe;
            
        end
    end
end





%% Elemental matrices
function [Df,Ef,Hf,Lf,Qf,fHf] = elementalMatricesBohm(Df,Ef,Hf,Lf,Qf,Xel,tel,refEl,ifac,solf,iFace)

% mesh data
global neq Magnetic testcase kT useThreshold limitRhoMin axisym refElTor

Ndim        = 3;
htheta      = tel(end)-tel(1);
nodesv      = refElTor.faceNodes3(ifac,:);
nodes2      = refEl.faceNodes(ifac,:);
Np2d        = size(refEl.NodesCoord,1);
Np1dPol     = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor     = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Nfl         = Np1dPol*Np1dTor;
ind_ff      = Np2d*neq + (ifac-1)*Nfl*neq+ (1:Nfl*neq);
ind_fe      = col(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'));
ind_fq      = col(bsxfun(@plus,(nodesv-1)*Ndim*neq,(1:Ndim*neq)'));
xyfg        = N1dPol*Xel(nodes2,:);      % Gauss points coordinates of the face (xy component)
thetafg     = N1dTor*tel;               % Gauss points coordinates of the face (theta component)
ipw1dp      = refEl.IPweights1d;
ipw1dt      = refElTor.IPweights1d;
ngausstor   = length(ipw1dt);
ngausspol   = length(ipw1dp);
N1dxPol     = refEl.N1dxi;

% computation of dsurf for the toroidal faces
xyd_g     = N1dxPol*Xel(nodes2,:);
xydNorm_g = sqrt(xyd_g(:,1).^2+xyd_g(:,2).^2);
dsurf     = col((ipw1dp.*xydNorm_g)*(ipw1dt*0.5*htheta)');

% computation of the outer normal
t_g  = xyd_g./xydNorm_g;
n_g  = repmat([t_g(:,2) -t_g(:,1)],ngausstor,1);
n_g(:,3) = 0;

% Shape functions on the face
Nfi = refElTor.sfTor;

% Radius at Gauss points
rg = col(repmat(xyfg(:,1),ngausstor,1));       % radius
if axisym
    dsurf = dsurf.*rg;
end

% solution at Gauss points
solf = transpose(reshape(solf,neq,Nfl));
sol_g = Nfi*solf;

%% FACES COMPUTATIONS:

% Initialization
B11 = zeros(Nfl); B12 = B11;
B21 = B11; B22 = B11;
Mf = B11;
fHf = zeros(neq*nf*nv,1);
fH = zeros(nv,neq);
if useThreshold
    sol_g(sol_g(:,1)<useThreshold,1) = useThreshold;
end
ug = sol_g(:,2)./sol_g(:,1);

%  LOOP IN GAUSS POINTS
for igtor = 1:ngausstor
    for igpol = 1:ngausspol
        g = (igtor-1)*ngausspol+igpol;
        
        if (testcase.n >= 50 && testcase.n<60)
            b = [Magnetic.bxfaces(g,iFace), Magnetic.byfaces(g,iFace)];
        else
            b = defineMagneticField3d([xyfg(igpol,1),xyfg(igpol,2)],thetafg(igtor));
        end
        
        if abs(dot(n_g,b/norm(b)))<1e-1
            sn = 1;
            delta = 0;
        else
            sn = sign(dot(n_g,b));
            delta = sn*ug(g)<sqrt(kT);
        end
        
        Wm = BohmMatrices(delta,delta_rho, sn);
        Nfg = Nfi(g,:);
        NN = Nfg'*Nfg*dsurf(g);
        Mf = Mf + (NN);
        B11 = B11 + Wm(1,1)*(NN);
        B12 = B12 + Wm(1,2)*(NN);
        B21 = B21 + Wm(2,1)*(NN);
        B22 = B22 + Wm(2,2)*(NN);
        fH = fH + Nfv_g'*[delta_rho*limitRhoMin,0]*dsurf(g);
    end
end

% expand the matrices
B = expandMatrixCv(B11,B12,B21,B22);
E_loc = expandMatrixA(Mf,neq);

% elemental assembly
Df(ind_ff,ind_fe) = B;
Ef(ind_ff,ind_ff) = E_loc;
Hf(ind_ff,ind_ff) = 0;
fHf(ind_ff) = -col(fH');







%% Elemental matrices
function [Df,Ef,Hf,Lf,Qf,fHf] = elementalMatricesCore(Df,Ef,Hf,Lf,Qf,Xel,tel,refEl,ifac,iFace)

% mesh data
global neq Magnetic testcase coreflux axisym refElTor

Ndim        = 3;
htheta      = tel(end)-tel(1);
nodesv      = refElTor.faceNodes3(ifac,:);
nodes2      = refEl.faceNodes(ifac,:);
Np2d        = size(refEl.NodesCoord,1);
Np1dPol     = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor     = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Nfl         = Np1dPol*Np1dTor;
ind_ff      = Np2d*neq + (ifac-1)*Nfl*neq+ (1:Nfl*neq);
ind_fe      = col(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'));
ind_fq      = col(bsxfun(@plus,(nodesv-1)*Ndim*neq,(1:Ndim*neq)'));
xyfg        = refEl.N1d*Xel(nodes2,:);      % Gauss points coordinates of the face (xy component)
thetafg     = refElTor.N1d*tel;               % Gauss points coordinates of the face (theta component)
ipw1dp      = refEl.IPweights1d;
ipw1dt      = refElTor.IPweights1d;
ngausstor   = length(ipw1dt);
ngausspol   = length(ipw1dp);
N1dxPol     = refEl.N1dxi;

% computation of dsurf for the toroidal faces
xyd_g     = N1dxPol*Xel(nodes2,:);
xydNorm_g = sqrt(xyd_g(:,1).^2+xyd_g(:,2).^2);
dsurf     = col((ipw1dp.*xydNorm_g)*(ipw1dt*0.5*htheta)');

% computation of the outer normal
t_g  = xyd_g./xydNorm_g;
n_g  = repmat([t_g(:,2) -t_g(:,1)],ngausstor,1);
n_g(:,3) = 0;

% Shape functions on the face
Nfi = refElTor.sfTor;

% Radius at Gauss points
rg = col(repmat(xyfg(:,1),ngausstor,1));       % radius
if axisym
    dsurf = dsurf.*rg;
end

% Initialization
B11 = zeros(nv); B12 = B11;
B21 = B11; B22 = B11;
fHf = zeros(neq*Nfp,1);
fH = zeros(nv,neq);

%  LOOP IN GAUSS POINTS
for igtor = 1:ngausstor
    for igpol = 1:ngausspol
        g = (igtor-1)*ngausspol+igpol;
        
        % Shape functions and derivatives at the current integration point
        Nfg = Nfi(g,:);
        NN = Nfg'*Nfg*dsurf(g);
        B11 = B11 + NN;
        B22 = B22 + NN;
        fH = fH + Nfg'*coreflux'*dsurf(g);
        
    end
end
% expand the matrices
E_loc = expandMatrixCv(B11, B12, B21,B22);

% elemental assembly
Df(ind_ff,ind_fe) = 0;
Hf(ind_ff,ind_ff) = 0;

if (testcase.n==51 || testcase.n==52 || testcase.n==53 || testcase.n==54 || testcase.n==55)
    Ef(ind_ff,ind_ff) = 0;
    Qf(ind_ff,ind4) = 0;
else
    
    Ef(ind_ff,ind_ff) = E_loc;
    Lf(ind_ff,ind4) = 0;
    Qf(ind_ff,ind4) = 0;
    fHf(ind_ff) = -col(fH');
end




%% Elemental matrices
function [Df,Ef,Hf,Lf,Qf,fHf] = elementalMatricesDirichletWeak(Df,Ef,Hf,Lf,Qf,Xel,tel,refEl,ifac,iFace,iel)

% mesh data
global neq   axisym refElTor
Ndim        = 3;
htheta      = tel(end)-tel(1);
nodesv      = refElTor.faceNodes3(ifac,:);
nodes2      = refEl.faceNodes(ifac,:);
Np2d        = size(refEl.NodesCoord,1);
Np1dPol     = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor     = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Nfl         = Np1dPol*Np1dTor;
ind_ff      = Np2d*neq + (ifac-1)*Nfl*neq+ (1:Nfl*neq);
ind_fe      = col(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'));
ind_fq      = col(bsxfun(@plus,(nodesv-1)*Ndim*neq,(1:Ndim*neq)'));
xyfg        = refEl.N1d*Xel(nodes2,:);      % Gauss points coordinates of the face (xy component)
thetafg     = refElTor.N1d*tel;               % Gauss points coordinates of the face (theta component)
ipw1dp      = refEl.IPweights1d;
ipw1dt      = refElTor.IPweights1d;
ngausstor   = length(ipw1dt);
ngausspol   = length(ipw1dp);
N1dxPol     = refEl.N1dxi;
nv          = Np1dPol*Np1dTor;
nf          = size(refEl.faceNodes,1);
Nfp         = Np2d*2+nf*Np1dPol*Np1dTor;


% computation of dsurf for the toroidal faces
xyd_g     = N1dxPol*Xel(nodes2,:);
xydNorm_g = sqrt(xyd_g(:,1).^2+xyd_g(:,2).^2);
dsurf     = col((ipw1dp.*xydNorm_g)*(ipw1dt*0.5*htheta)');

% Shape functions on the face
Nfi = refElTor.sfTor;

% Radius at Gauss points
rg = col(repmat(xyfg(:,1),ngausstor,1));       % radius
if axisym
    dsurf = dsurf.*rg;
end

% Initialization
Mf = zeros(nv); 
fHf = zeros(neq*Nfp,1);
fH = zeros(nv,neq);

sol_g = analyticalSolution3d(xyfg,thetafg);

%  LOOP IN GAUSS POINTS
for igtor = 1:ngausstor
    for igpol = 1:ngausspol
        g = (igtor-1)*ngausspol+igpol;
        
        % Shape functions and derivatives at the current integration point
        Nfg = Nfi(g,:);
        Mf = Mf+Nfg'*Nfg*dsurf(g);
        fH = fH + Nfg'*sol_g(g,:)*dsurf(g);
        
    end
end
% expand the matrices
E_loc = expandMatrixA(Mf,2);

% if (iel==15) 
%     stop
% end

% E_loc(E_loc==0) = 1e-8;

% elemental assembly
Df(ind_ff,ind_fe) = 0;
Hf(ind_ff,ind_ff) = 0;
Ef(ind_ff,ind_ff) = E_loc;
Lf(ind_ff,ind_fq) = 0;
Qf(ind_ff,ind_fq) = 0;
fHf(ind_ff) = -col(fH');

% Df(ind_ff,ind_fe) = 1e-8;




%% Elemental matrices
% function [Df,Ef,Hf,Lf,fHf] = elementalMatricesDirichlet(Df,Ef,Hf,Lf,X,refEl,ifac,iFace)
%
% % mesh data
% global neq Magnetic testcase coreflux axisym
%
% nv = size(refEl.NodesCoord1d,1);
% faceNodesv = refEl.faceNodes;
% nf = size(faceNodesv,1);
% nd = size(X,2);
%
% % Information of the reference element
% IPw_fv = refEl.IPweights1d;
% N1dv = refEl.N1d;
% Nx1dv = refEl.N1dxi;
%
% % % solution at Gauss points
% % solf = transpose(reshape(solf,neq,nv));
% % sol_g = N1dv*solf;
%
% %% FACES COMPUTATIONS:
% ngauss_f = length(IPw_fv);
%
% % face nodes
% nodesv = faceNodesv(ifac,:);
%
% % exact velocity in the face
% xe = X(:,1); ye = X(:,2);
% xf = xe(nodesv);
% yf = ye(nodesv);
% xyfg = N1dv*[xf yf];
% u_ex = analyticalSolution(xyfg);
%
% % indices for local assembly
% ind_ff = (ifac-1)*neq*nv + (1:neq*nv);
% ind_fe = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
% ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for velocity gradient
%
% % Initialization
% B11 = zeros(nv); B12 = B11;
% B21 = B11; B22 = B11;
% fHf = zeros(neq*nf*nv,1);
% fH = zeros(nv,neq);
%
% %  LOOP IN GAUSS POINTS
% for g = 1:ngauss_f
%
%     % Shape functions and derivatives at the current integration point
%     Nfv_g = N1dv(g,:);
%     Nfxiv_g = Nx1dv(g,:);
%
%     % Integration weight
%     xyDer_g = Nfxiv_g*[xf yf];
%     xyDerNorm_g = norm(xyDer_g);
%     dline = IPw_fv(g)*xyDerNorm_g;
%     if axisym
%         x = Nfv_g*xf;
%         dline = dline*x;
%     end
%     % Unit normal to the boundary
%     t_g = xyDer_g/xyDerNorm_g;
%     n_g = [t_g(2) -t_g(1)];
%
%     if (testcase.n >= 50 && testcase.n<60)
%         b = [Magnetic.bxfaces(g,iFace), Magnetic.byfaces(g,iFace)];
%     else
%         b = defineMagneticField([xyfg(g,1),xyfg(g,2)]);
%     end
%
%     B11 = B11 + (Nfv_g'*Nfv_g)*dline;
%     %     B12 = B12 + (Nfv_g'*Nfv_g)*dline;
%     %     B21 = B21 + (Nfv_g'*Nfv_g)*dline;
%     B22 = B22 + (Nfv_g'*Nfv_g)*dline;
%     fH = fH + Nfv_g'*u_ex(g,:)*dline;
%
% end
%
% % expand the matrices
% E_loc = expandMatrixCv(B11, zeros(size(B11)),zeros(size(B11)),B22);
%
% % elemental assembly
% Df(ind_ff,ind_fe) = 0;
% Ef(ind_ff,ind_ff) = E_loc;
% Lf(ind_ff,ind4) = 0;
% Hf(ind_ff,ind_ff) = 0;
% fHf(ind_ff) = -col(fH');







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

function res = expandMatrixCv(Cxx,Cxy,Cyx,Cyy)
% expand matrix Cv
%   [ Cxx Cxy
%     Cyx Cyy]
res = zeros([size(Cxx) 2 2]);
res(:,:,[1 3 2 4]) = cat(3,Cxx,Cxy,Cyx,Cyy);
res = permute(res, [3 1 4 2]);
res = reshape(res, 2*size(Cxx,1),2*size(Cxx,2));

% function res = wallMatrix(nx,ny)
%
% res = [1 0 0; 0 1-nx^2 -nx*ny; 0 -nx*ny 1-ny^2];
%
% function [An,Anp,Anm] = jacobianMatrices(x,y,U,nx,ny)
%
% global kT useThreshold
% rho = U(1);
% Gamma = U(2);
%
% if useThreshold
%     if U(1)<useThreshold,
%         rho = useThreshold;
%         Gamma = 0;
%
%     end
% end
%
% b = defineMagneticField([x,y]);
% bn = dot(b,[nx,ny]);
%
% An = [0, bn; (-1*Gamma^2/rho^2+kT)*bn, 2*Gamma/rho*bn];
% [V,D] = eig(An);
% invV = inv(V);
% Anp = 0.5*(An+V*abs(D)*invV);
% Anm = 0.5*(An-V*abs(D)*invV);

function B = BohmMatrices(delta,delta_rho,sn)

global kT

B = [1-delta_rho, 0; delta*sn*sqrt(kT), (1-delta)];


