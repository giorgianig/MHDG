function [Df,Ef,Hf,Lf,Qf,fH] = hdg_bc(Df,Ef,Hf,Lf,Qf,X,T,infoFaces,refEl,bn,F,sol_hat)

global neq testcase
%
Ne = size(Df,3);
nv = size(refEl.NodesCoord1d,1);   % number of face nodes for the velocity
nf = size(F,2);

% number of boundaries
nb = numel(bn);
fH = zeros(nf*nv*neq,Ne);
for ibound = 1:nb
    
    name = bn{ibound}(4:end);
    
    % skip Dirichlet faces
    if testcase.n<50
        if any(strcmp(name,{'Diriclet','IN','OUT'})), continue ,end
    else
        if any(strcmp(name,{'Diriclet'})), continue ,end
    end
    faces = infoFaces.(['exteriorFaces_' name]);
    
    for iF = 1:size(faces,1)
        
        iel = faces(iF,1);
        ifac = faces(iF,2);
        Fe = F(iel,ifac);
        Te = T(iel,:);
        Xe = X(Te,:);
        Dfe = Df(:,:,iel);
        Efe = Ef(:,:,iel);
        Hfe = Hf(:,:,iel);
        Lfe = Lf(:,:,iel);
        Qfe = Qf(:,:,iel);
        ind_uhat =  (Fe-1)*neq*nv+(1:neq*nv);
        solf = sol_hat(ind_uhat);
        fHe = fH(:,iel);
        if any(strcmp(name,{'UP','DOWN','LEFT','RIGHT','WALL','ULIM'}))
            [Dfe,Efe,Hfe,Lfe,Qfe] = elementalMatricesWall(Dfe,Efe,Hfe,Lfe,Qfe,Xe,refEl,ifac);
        elseif any(strcmp(name,{'INLET','OUTLET'}))
            [Dfe,Efe,Hfe,Qfe,fHe] = elementalMatricesInlOut(Dfe,Efe,Hfe,Qfe,Xe,refEl,ifac,solf);
        elseif any(strcmp(name,{'Neumann'}))
            [Dfe,Efe,Hfe,Qfe,fHe] = elementalMatricesNeumann(Dfe,Efe,Hfe,Qfe,Xe,refEl,ifac);
        elseif any(strcmp(name,{'IN'}))
            [Dfe,Efe,Hfe,Lfe,Qfe,fHe] = elementalMatricesCore(Dfe,Efe,Hfe,Lfe,Qfe,Xe,refEl,ifac,Fe);
        elseif any(strcmp(name,{'LIM','OUT'}))
            [Dfe,Efe,Hfe,Lfe,Qfe,fHe] = elementalMatricesBohm(Dfe,Efe,Hfe,Lfe,Qfe,Xe,refEl,ifac,solf,Fe);
            %         elseif any(strcmp(name,{'Diriclet'}))
            %             [Dfe,Efe,Hfe,Lfe,Qfe,fHe] = elementalMatricesDirichlet(Dfe,Efe,Hfe,Lfe,Qfe,Xe,refEl,ifac,Fe);
        else
            error('Something wrong')
        end
        Df(:,:,iel) = Dfe;
        Ef(:,:,iel) = Efe;
        Hf(:,:,iel) = Hfe;
        Lf(:,:,iel) = Lfe;
        Qf(:,:,iel) = Qfe;
        fH(:,iel)   = fHe;
        
    end
end



%% Elemental matrices
function [Df,Ef,Hf,Lf,Qfe] = elementalMatricesWall(Df,Ef,Hf,Lf,Qfe,X,refEl,ifac)

% mesh data
global neq axisym
nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nd = size(X,2);

% Information of the reference element
IPw_fv = refEl.IPweights1d;
N1dv = refEl.N1d;
Nx1dv = refEl.N1dxi;

% x and y coordinates of the element nodes
xe = X(:,1); ye = X(:,2);
% %% FACES COMPUTATIONS:
% ngauss_f = length(IPw_fv);
%
% % face nodes
% nodesv = faceNodesv(ifac,:);
%
% % indices for local assembly
% ind_face_2 = (ifac-1)*neq*nv + (1:neq*nv);
% ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
%
% xf = xe(nodesv);
% yf = ye(nodesv);
%
% % Initialization
% B11 = zeros(nv); B12 = B11;
% B21 = B11; B22 = B11;
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
%
%     % Unit normal to the boundary
%     t_g = xyDer_g/xyDerNorm_g;
%     n_g = [t_g(2) -t_g(1)];
%
%     Wm = wallMatrix(n_g(1),n_g(2));
%
%     B11 = B11 + Wm(1,1)*(Nfv_g'*Nfv_g)*dline;
%     B12 = B12 + Wm(1,2)*(Nfv_g'*Nfv_g)*dline;
%     B21 = B21 + Wm(2,1)*(Nfv_g'*Nfv_g)*dline;
%     B22 = B22 + Wm(2,2)*(Nfv_g'*Nfv_g)*dline;
% end
%
% % expand the matrices
% B = expandMatrixCv(B11,B12,B21,B22);
% E_loc = expandMatrixA(B11,neq);
%
% % elemental assembly
% Df(ind_face_2,ind2) = B;
% Ef(ind_face_2,ind_face_2) = E_loc;
% Hf(ind_face_2,ind_face_2) = 0;

%% FACES COMPUTATIONS:
ngauss_f = length(IPw_fv);

% face nodes
nodesv = faceNodesv(ifac,:);

% exact velocity in the face
xe = X(:,1); ye = X(:,2);
xf = xe(nodesv);
yf = ye(nodesv);
xyfg = N1dv*[xf yf];

% indices for local assembly
ind_face_2 = (ifac-1)*neq*nv + (1:neq*nv);
ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for unknown gradient

% Initialization
B11 = zeros(nv); B12 = B11;
B21 = B11; B22 = B11;
Mf = B11;


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
    %     t_g = xyDer_g/xyDerNorm_g;
    %     n_g = [t_g(2) -t_g(1)];
    %
    %
    %     delta = 0;
    %     sn = 1;
    %     delta = sn*ug(g)<sqrt(kT);
    
    %     Wm = BohmMatrices(delta, sn);
    
    Mf = Mf + (Nfv_g'*Nfv_g)*dline;
    B11 = B11 + (Nfv_g'*Nfv_g)*dline;
    %     B12 = B12 + Wm(1,2)*(Nfv_g'*Nfv_g)*dline;
    %     B21 = B21 + Wm(2,1)*(Nfv_g'*Nfv_g)*dline;
    B22 = B22 + (Nfv_g'*Nfv_g)*dline;
end

% expand the matrices
B = expandMatrixCv(B11,B12,B21,B22);
E_loc = expandMatrixA(Mf,neq);

% elemental assembly
Df(ind_face_2,ind2) = B;
Ef(ind_face_2,ind_face_2) = E_loc;
Hf(ind_face_2,ind_face_2) = 0;
% Lf(ind_face_2,ind4) = 0;

%% Elemental matrices
function [Df,Ef,Hf,Qf,fHf] = elementalMatricesInlOut(Df,Ef,Hf,Qf,X,refEl,ifac,solf)

% mesh data
global neq axisym
nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);

% Information of the reference element
IPw_fv = refEl.IPweights1d;
N1dv = refEl.N1d;
Nx1dv = refEl.N1dxi;

% solution at Gauss points
solf = transpose(reshape(solf,neq,nv));
sol_g = N1dv*solf;

%% FACES COMPUTATIONS:
ngauss_f = length(IPw_fv);

% face nodes
nodesv = faceNodesv(ifac,:);

% exact velocity in the face
xe = X(:,1); ye = X(:,2);
xf = xe(nodesv);
yf = ye(nodesv);
xyfg = N1dv*[xf yf];
sol_ex = setDirichletBoundaryConditions(xyfg);

% indices for local assembly
ind_face_2 = (ifac-1)*neq*nv + (1:neq*nv);
ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity

% Initialization
Hp11 = zeros(nv); Hp12 = Hp11;
Hp21 = Hp11; Hp22 = Hp11;

Hm11 = zeros(nv); Hm12 = Hm11;
Hm21 = Hm11; Hm22 = Hm11;
fHf = zeros(neq*nf*nv,1);
fH = zeros(nv,neq);

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
    
    [An,Anp,Anm] = jacobianMatrices(xyfg(g,1),xyfg(g,2),sol_g(g,:),n_g(1),n_g(2));
    
    Hp11 = Hp11 + Nfv_g'*Nfv_g*Anp(1,1)*dline;
    Hp12 = Hp12 + Nfv_g'*Nfv_g*Anp(1,2)*dline;
    Hp21 = Hp21 + Nfv_g'*Nfv_g*Anp(2,1)*dline;
    Hp22 = Hp22 + Nfv_g'*Nfv_g*Anp(2,2)*dline;
    
    Hm11 = Hm11 + Nfv_g'*Nfv_g*Anm(1,1)*dline;
    Hm12 = Hm12 + Nfv_g'*Nfv_g*Anm(1,2)*dline;
    Hm21 = Hm21 + Nfv_g'*Nfv_g*Anm(2,1)*dline;
    Hm22 = Hm22 + Nfv_g'*Nfv_g*Anm(2,2)*dline;
    
    fH = fH + Nfv_g'*transpose(Anm*sol_ex(g,:)' )*dline;
end

% expand the matrices
Hp = expandMatrixCv(Hp11,Hp12,Hp21,Hp22);
Hm = expandMatrixCv(Hm11,Hm12,Hm21,Hm22);

% elemental assembly
Df(ind_face_2,ind2) = Hp;
Ef(ind_face_2,ind_face_2) = Hp-Hm;
Hf(ind_face_2,ind_face_2) = 0;
fHf(ind_face_2) = col(fH');


%% Elemental matrices
function [Df,Ef,Hf,Qf,fHf] = elementalMatricesNeumann(Df,Ef,Hf,Qf,X,refEl,ifac)

% mesh data
global neq axisym
nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);

% Information of the reference element
IPw_fv = refEl.IPweights1d;
N1dv = refEl.N1d;
Nx1dv = refEl.N1dxi;

%% FACES COMPUTATIONS:
ngauss_f = length(IPw_fv);

% face nodes
nodesv = faceNodesv(ifac,:);

% exact velocity in the face
xe = X(:,1); ye = X(:,2);
xf = xe(nodesv);
yf = ye(nodesv);
xyfg = N1dv*[xf yf];
% exact velocity in the faces for Dirichlet boundary
sol_g = analyticalSolution(xyfg);

% indices for local assembly
ind_face_2 = (ifac-1)*neq*nv + (1:neq*nv);
ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity

% Initialization
fHf = zeros(neq*nf*nv,1);
fH = zeros(nv,neq);

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
    An = jacobianMatrices(xyfg(g,1),xyfg(g,2),sol_g(g,:),n_g(1),n_g(2));
    fH = fH + Nfv_g'*transpose(An*sol_g(g,:)'  )*dline;
end

% elemental assembly
% Df(ind_face_2,ind2) = 0;
% Ef(ind_face_2,ind_face_2) = 0;
% Hf(ind_face_2,ind_face_2) = 0;
fHf(ind_face_2) = col(fH');




%% Elemental matrices
function [Df,Ef,Hf,Lf,Qf,fHf] = elementalMatricesBohm(Df,Ef,Hf,Lf,Qf,X,refEl,ifac,solf,iFace)


% mesh data
global neq Magnetic testcase kT useThreshold limitRhoMin axisym

nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);
nd = size(X,2);

% Information of the reference element
IPw_fv = refEl.IPweights1d;
N1dv = refEl.N1d;
Nx1dv = refEl.N1dxi;

% solution at Gauss points
solf = transpose(reshape(solf,neq,nv));
sol_g = N1dv*solf;

%% FACES COMPUTATIONS:
ngauss_f = length(IPw_fv);

% face nodes
nodesv = faceNodesv(ifac,:);

% exact velocity in the face
xe = X(:,1); ye = X(:,2);
xf = xe(nodesv);
yf = ye(nodesv);
xyfg = N1dv*[xf yf];

% indices for local assembly
ind_face_2 = (ifac-1)*neq*nv + (1:neq*nv);
ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for unknown gradient
indLdens = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:2)'),neq*nv,1); % assembly face to elem for unknown gradient (only dens)
indLmome = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(3:4)'),neq*nv,1); % assembly face to elem for unknown gradient (only mom)

% Initialization
B11 = zeros(nv); B12 = B11;
B21 = B11; B22 = B11;
Mf = B11;
fHf = zeros(neq*nf*nv,1);
fH = zeros(nv,neq);
if useThreshold
    sol_g(sol_g(:,1)<useThreshold,1) = useThreshold;
end
ug = sol_g(:,2)./sol_g(:,1);

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
    
    if (testcase.n >= 50 && testcase.n<60)
        b = [Magnetic.bxfaces(g,iFace), Magnetic.byfaces(g,iFace)];
    else
        b = defineMagneticField([xyfg(g,1),xyfg(g,2)]);
    end
    
    if abs(dot(n_g,b/norm(b)))<1e-1;
        sn = 1;
        delta = 0;
    else
        sn = sign(dot(n_g,b));
        delta = sn*ug(g)<sqrt(kT);
    end
    
    if limitRhoMin
        delta_rho = sol_g(g,1)<limitRhoMin;
    else
        delta_rho = 0;
    end
    
    %% delta=1 impose the velocity to the speed of sound
    %     delta=1;
    
    Wm = BohmMatrices(delta,delta_rho, sn);
    
    Mf = Mf + (Nfv_g'*Nfv_g)*dline;
    B11 = B11 + Wm(1,1)*(Nfv_g'*Nfv_g)*dline;
    B12 = B12 + Wm(1,2)*(Nfv_g'*Nfv_g)*dline;
    B21 = B21 + Wm(2,1)*(Nfv_g'*Nfv_g)*dline;
    B22 = B22 + Wm(2,2)*(Nfv_g'*Nfv_g)*dline;
    fH = fH + Nfv_g'*[delta_rho*limitRhoMin,0]*dline;
end


% expand the matrices
B = expandMatrixCv(B11,B12,B21,B22);
E_loc = expandMatrixA(Mf,neq);

% elemental assembly
Df(ind_face_2,ind2) = B;
Ef(ind_face_2,ind_face_2) = E_loc;
Hf(ind_face_2,ind_face_2) = 0;
%  Lf(ind_face_2,indLmome) = 0;
% Lf(ind_face_2,ind4) = 0;
fHf(ind_face_2) = -col(fH');







%% Elemental matrices
function [Df,Ef,Hf,Lf,Qf,fHf] = elementalMatricesCore(Df,Ef,Hf,Lf,Qf,X,refEl,ifac,iFace)

% mesh data
global neq Magnetic testcase coreflux axisym

nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);
nd = size(X,2);

% Information of the reference element
IPw_fv = refEl.IPweights1d;
N1dv = refEl.N1d;
Nx1dv = refEl.N1dxi;

% % solution at Gauss points
% solf = transpose(reshape(solf,neq,nv));
% sol_g = N1dv*solf;

%% FACES COMPUTATIONS:
ngauss_f = length(IPw_fv);

% face nodes
nodesv = faceNodesv(ifac,:);

% exact velocity in the face
xe = X(:,1); ye = X(:,2);
xf = xe(nodesv);
yf = ye(nodesv);
xyfg = N1dv*[xf yf];

% indices for local assembly
ind_face_2 = (ifac-1)*neq*nv + (1:neq*nv);
ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for velocity gradient

% Initialization
B11 = zeros(nv); B12 = B11;
B21 = B11; B22 = B11;
fHf = zeros(neq*nf*nv,1);
fH = zeros(nv,neq);

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
    B11 = B11 + (Nfv_g'*Nfv_g)*dline;
    B22 = B22 + (Nfv_g'*Nfv_g)*dline;
    fH = fH + Nfv_g'*coreflux'*dline;
    
end

% expand the matrices
E_loc = expandMatrixCv(B11, B12, B21,B22);

% elemental assembly
Df(ind_face_2,ind2) = 0;
Hf(ind_face_2,ind_face_2) = 0;

if (testcase.n==51 || testcase.n==52 || testcase.n==53 || testcase.n==54 || testcase.n==55)
    Ef(ind_face_2,ind_face_2) = 0;
    Qf(ind_face_2,ind4) = 0;
else
    
    Ef(ind_face_2,ind_face_2) = E_loc;
    Lf(ind_face_2,ind4) = 0;
    Qf(ind_face_2,ind4) = 0;
    fHf(ind_face_2) = -col(fH');
end








%% Elemental matrices
function [Df,Ef,Hf,Lf,fHf] = elementalMatricesDirichlet(Df,Ef,Hf,Lf,X,refEl,ifac,iFace)

% mesh data
global neq Magnetic testcase coreflux axisym

nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);
nd = size(X,2);

% Information of the reference element
IPw_fv = refEl.IPweights1d;
N1dv = refEl.N1d;
Nx1dv = refEl.N1dxi;

% % solution at Gauss points
% solf = transpose(reshape(solf,neq,nv));
% sol_g = N1dv*solf;

%% FACES COMPUTATIONS:
ngauss_f = length(IPw_fv);

% face nodes
nodesv = faceNodesv(ifac,:);

% exact velocity in the face
xe = X(:,1); ye = X(:,2);
xf = xe(nodesv);
yf = ye(nodesv);
xyfg = N1dv*[xf yf];
u_ex = analyticalSolution(xyfg);

% indices for local assembly
ind_face_2 = (ifac-1)*neq*nv + (1:neq*nv);
ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for velocity gradient

% Initialization
B11 = zeros(nv); B12 = B11;
B21 = B11; B22 = B11;
fHf = zeros(neq*nf*nv,1);
fH = zeros(nv,neq);

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
    
    if (testcase.n >= 50 && testcase.n<60)
        b = [Magnetic.bxfaces(g,iFace), Magnetic.byfaces(g,iFace)];
    else
        b = defineMagneticField([xyfg(g,1),xyfg(g,2)]);
    end
    
    B11 = B11 + (Nfv_g'*Nfv_g)*dline;
    %     B12 = B12 + (Nfv_g'*Nfv_g)*dline;
    %     B21 = B21 + (Nfv_g'*Nfv_g)*dline;
    B22 = B22 + (Nfv_g'*Nfv_g)*dline;
    fH = fH + Nfv_g'*u_ex(g,:)*dline;
    
end

% expand the matrices
E_loc = expandMatrixCv(B11, zeros(size(B11)),zeros(size(B11)),B22);

% elemental assembly
Df(ind_face_2,ind2) = 0;
Ef(ind_face_2,ind_face_2) = E_loc;
Lf(ind_face_2,ind4) = 0;
Hf(ind_face_2,ind_face_2) = 0;
fHf(ind_face_2) = -col(fH');







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

function res = wallMatrix(nx,ny)

res = [1 0 0; 0 1-nx^2 -nx*ny; 0 -nx*ny 1-ny^2];

function [An,Anp,Anm] = jacobianMatrices(x,y,U,nx,ny)

global kT useThreshold
rho = U(1);
Gamma = U(2);

if useThreshold
    if U(1)<useThreshold,
        rho = useThreshold;
        Gamma = 0;
        
    end
end

b = defineMagneticField([x,y]);
bn = dot(b,[nx,ny]);

An = [0, bn; (-1*Gamma^2/rho^2+kT)*bn, 2*Gamma/rho*bn];
[V,D] = eig(An);
invV = inv(V);
Anp = 0.5*(An+V*abs(D)*invV);
Anm = 0.5*(An-V*abs(D)*invV);

function B = BohmMatrices(delta,delta_rho,sn)

global kT

B = [1-delta_rho, 0; delta*sn*sqrt(kT), (1-delta)];


