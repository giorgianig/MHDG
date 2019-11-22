function [Df,Ef,Hf,Lf,Qf,fH,TUhf,TQhf,Tfhf] = hdg_bc(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,X,T,infoFaces,refEl,bn,F,sol_hat,solQ)

global neq testcase useNeumannBC
%
Ne = size(Df,3);
nd = size(X,2);
Nv =size(refEl.NodesCoord,1);        % number of element nodes
nv = size(refEl.NodesCoord1d,1);   % number of face nodes
nf = size(F,2);

% number of boundaries
nb = numel(bn);
fH = zeros(nf*nv*neq,Ne);
for ibound = 1:nb
    
    name = bn{ibound}(4:end);
    
    if testcase.n ==69
        continue
    end
    % skip Dirichlet faces
    if testcase.n<50
        if testcase.n~=5
            if any(strcmp(name,{'Diriclet','IN','OUT','UP','DOWN','LEFT','RIGHT'})), continue ,end
        end
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
        TUhfe = TUhf(:,:,iel);
        TQhfe = TQhf(:,:,iel);
        Tfhfe   =  Tfhf(:,iel);
        ind_uhat =  (Fe-1)*neq*nv+(1:neq*nv);
        indQ = (iel-1)*neq*Nv*nd + (1:neq*Nv*nd);
        solf = sol_hat(ind_uhat);
        solQe = solQ(indQ);
        solQe    = reshape(solQe,neq*nd,numel(solQe)/neq/nd)';
        nodesv = refEl.faceNodes(ifac,:);
        solQf    = solQe(nodesv,:);
        
        fHe = fH(:,iel);
        if any(strcmp(name,{'UP','DOWN','WALL','ULIM'}))
            [Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe] = elementalMatricesWall(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,refEl,ifac);
        elseif any(strcmp(name,{'INLET','OUTLET'}))
            [Dfe,Efe,Hfe,Qfe,fHe] = elementalMatricesInlOut(Dfe,Efe,Hfe,Qfe,Xe,refEl,ifac,solf);
        elseif any(strcmp(name,{'Neumann'}))
            [Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe] = elementalMatricesNeumann(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,refEl,ifac);
        elseif any(strcmp(name,{'IN','LEFT'}))
            [Dfe,Efe,Hfe,Lfe,Qfe,fHe,TUhfe,TQhfe,Tfhfe] = elementalMatricesCore(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,refEl,ifac,Fe);
        elseif any(strcmp(name,{'LIM','RIGHT','OUT'}))
            if useNeumannBC~=2
                [Dfe,Efe,Hfe,Lfe,Qfe,fHe,TUhfe,TQhfe,Tfhfe] = elementalMatricesBohm(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,refEl,ifac,solf,solQf,Fe,iel);
            else
                [Dfe,Efe,Hfe,Lfe,Qfe,fHe,TUhfe,TQhfe,Tfhfe] = elementalMatricesBohmSimp(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,refEl,ifac,solf,solQf,Fe,iel);
%                 [Dfe,Efe,Hfe,Lfe,Qfe,fHe,TUhfe,TQhfe,Tfhfe] = elementalMatricesBohmTest(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,refEl,ifac,solf,solQf,Fe,iel);
            end
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
        TUhf(:,:,iel) = TUhfe;
        TQhf(:,:,iel) = TQhfe;
        Tfhf(:,iel)    = Tfhfe;
        
        
        
    end
end



%% Elemental matrices
function [Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf] = elementalMatricesWall(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,X,refEl,ifac)

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
    %     B11 = B11 + (Nfv_g'*Nfv_g)*dline;
    %     %     B12 = B12 + Wm(1,2)*(Nfv_g'*Nfv_g)*dline;
    %     %     B21 = B21 + Wm(2,1)*(Nfv_g'*Nfv_g)*dline;
    %     B22 = B22 + (Nfv_g'*Nfv_g)*dline;
end

% expand the matrices
B =  expandMatrixA(Mf,neq);
E_loc = expandMatrixA(Mf,neq);

% elemental assembly
Df(ind_face_2,ind2) = B;
Ef(ind_face_2,ind_face_2) = E_loc;
Hf(ind_face_2,ind_face_2) = 0;
% Lf(ind_face_2,ind4) = 0;
Qf(ind_face_2,ind4) = 0;
TUhf(ind_face_2,ind_face_2) = 0;
TQhf(ind_face_2,ind4)           = 0;
Tfhf(ind_face_2)                    = 0;

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
function [Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf] = elementalMatricesNeumann(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,X,refEl,ifac)

% mesh data
global neq axisym
nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);
nd = size(X,2);

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
ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for unknown gradient


% elemental assembly
Df(ind_face_2,ind2) = 0;
Ef(ind_face_2,ind_face_2) = 0;
Hf(ind_face_2,ind_face_2) = 0;
TUhf(ind_face_2,ind_face_2) = 0;
TQhf(ind_face_2,ind4)           = 0;
Tfhf(ind_face_2)                    = 0;
Qf(ind_face_2,ind4)               = 0;


% Lf(ind_face_2,ind4)               =Lf(ind_face_2,ind4)/10;



%% Elemental matrices Bohm
function [Df,Ef,Hf,Lf,Qf,fHf,TUhf,TQhf,Tfhf] = elementalMatricesBohm(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,X,refEl,ifac,solf,solQf,iFace,iel)


% mesh data
global neq Magnetic testcase useThreshold limitRhoMin axisym epn diff_pari diff_pare Mref useNeumannBC decoupleEquations tau
global  diff_n diff_u diff_ei diff_ee
mult = [1 1 1 1];

coefi = (2/(3*Mref) )^(1+epn)*diff_pari;
coefe = (2/(3*Mref) )^(1+epn)*diff_pare;



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
solQ_g = N1dv*solQf;

sol_g_phys = cons2phys(sol_g);

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
indLdens = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:2)'),nd*nv,1); % assembly face to elem for unknown gradient (only dens)
indLmome = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(3:4)'),nd*nv,1); % assembly face to elem for unknown gradient (only mom)
indLener = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(5:6)'),nd*nv,1); % assembly face to elem for unknown gradient (only energy)

% Initialization
B11 = zeros(nv); B12 = B11; B13 = B11; B14 = B11;
B21 = B11; B22 = B11; B23 = B11; B24 = B11;
B31 = B11; B32 = B11; B33 = B11; B34 = B11;
B41 = B11; B42 = B11; B43 = B11; B44 = B11;
Mf = B11;
fHf = zeros(neq*nf*nv,1);
fH = zeros(nv,neq);
if useThreshold
    sol_g(sol_g(:,1)<useThreshold,1) = useThreshold;
end
ug = sol_g(:,2)./sol_g(:,1);

Cnx = zeros(nv);Cny=Cnx;Qbx=Cnx;Qby=Cnx;
Hn11 = zeros(nv); Hn12 = Hn11; Hn13 = Hn11; Hn14 = Hn11;
Hn21 = Hn11; Hn22 = Hn11; Hn23 = Hn11; Hn24 = Hn11;
Hn31 = Hn11; Hn32 = Hn11; Hn33 = Hn11; Hn34 = Hn11;
Hn41 = Hn11; Hn42 = Hn11; Hn43 = Hn11; Hn44 = Hn11;


TUh11 = zeros(nv); TUh12 = TUh11; TUh13 = TUh11; TUh14 = TUh11;
TUh21 = TUh11; TUh22 = TUh11; TUh23 = TUh11; TUh24 = TUh11;
TUh31 = TUh11; TUh32 = TUh11; TUh33 = TUh11; TUh34 = TUh11;
TUh41 = TUh11; TUh42 = TUh11; TUh43 = TUh11; TUh44 = TUh11;

TQh11 = zeros(nv); TQh12 = TQh11; TQh13 = TQh11; TQh14 = TQh11; TQh15 = TQh11; TQh16 = TQh11; TQh17 = TQh11;  TQh18 = TQh11;
TQh21 = TQh11; TQh22 = TQh11; TQh23 = TQh11; TQh24 = TQh11; TQh25 = TQh11; TQh26 = TQh11; TQh27 = TQh11;  TQh28 = TQh11;
TQh31 = TQh11; TQh32 = TQh11; TQh33 = TQh11; TQh34 = TQh11; TQh35 = TQh11; TQh36 = TQh11; TQh37 = TQh11;  TQh38 = TQh11;
TQh41 = TQh11; TQh42 = TQh11; TQh43 = TQh11; TQh44 = TQh11; TQh45 = TQh11; TQh46 = TQh11; TQh47 = TQh11;  TQh48 = TQh11;

Tfh1          = zeros(nv,1); Tfh2 = Tfh1; Tfh3 = Tfh1; Tfh4 = Tfh1;
zer11 = zeros(nv); zer12 = zeros(nv); zer13 = zeros(nv); zer14 = zeros(nv);
zer21 = zeros(nv); zer22 = zeros(nv); zer23 = zeros(nv); zer24 = zeros(nv);
zer31 = zeros(nv); zer32 = zeros(nv); zer33 = zeros(nv); zer34 = zeros(nv);
zer41 = zeros(nv); zer42 = zeros(nv); zer43 = zeros(nv); zer44 = zeros(nv);

% analytical rhs for convergence test
fan = rhsbohm(xyfg);

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
    
    soundSpeed = sol_g_phys(g,9) ;
    if ~isreal(soundSpeed)
        soundSpeed = abs(soundSpeed);
    end
    
    if decoupleEquations
        soundSpeed = sqrt(Mref);
    end
    
    
    if abs(dot(n_g,b/norm(b)))<1e-1
        sn = 1;
        delta = 0;
        iota  =0;
    else
        sn = sign(dot(n_g,b));
        delta = sn*ug(g)<soundSpeed;
        iota = 1;
    end
    
    if limitRhoMin
        delta_rho = sol_g(g,1)<limitRhoMin;
    else
        delta_rho = 0;
    end
    
    %% delta=1 impose the velocity to the speed of sound
    %         delta=0;
    
    
    % Compute G^T^(k-1)
    G = reshape(solQ_g(g,:),nd,neq);
    
    % Compute V(U^(k-1))
    Vi = computeVi(sol_g(g,:));
    Ve = computeVe(sol_g(g,:));
    
    % Compute dV_dU (k-1)
    dVi_dU = compute_dVi_dU(sol_g(g,:));
    dVe_dU = compute_dVe_dU(sol_g(g,:));
    
    % Compute Alpha(U^(k-1))
    Alphai = computeAlphai(sol_g(g,:));
    Alphae = computeAlphae(sol_g(g,:));
    
    % Compute dAlpha/dU^(k-1)
    dAlphai_dU = compute_dAlphai_dU(sol_g(g,:));
    dAlphae_dU = compute_dAlphae_dU(sol_g(g,:));
    
%% NEW****************************************************************    
%     NNf =  Nfv_g'*dline*dot(n_g,b);
%     NN  =  Nfv_g'*Nfv_g*dline*dot(n_g,b);    
    NNf =  Nfv_g'*dline*dot(n_g,b);
    NN  =  Nfv_g'*Nfv_g*dline*dot(n_g,b);
%% NEW****************************************************************    




    gammai = dot(G*Vi,b);    % scalar
    gammae = dot(G*Ve,b);    % scalar
    Taui = G*dVi_dU;      % 2x3
    Taue = G*dVe_dU;      % 2x3
    
    Abohm = jacobianMatricesBohm(sol_g(g,:));
    
    
    %      Hn11 = Hn11 + Nfv_g'*Nfv_g*An(1,1)*dline;
    %      Hn12 = Hn12 + Nfv_g'*Nfv_g*An(1,2)*dline;
    %      Hn13 = Hn13 + Nfv_g'*Nfv_g*An(1,3)*dline;
    %      Hn21 = Hn21 + Nfv_g'*Nfv_g*An(2,1)*dline;
    %      Hn22 = Hn22 + Nfv_g'*Nfv_g*An(2,2)*dline;
    %      Hn23 = Hn23 + Nfv_g'*Nfv_g*An(2,3)*dline;
    
    if iota
%         Hn31 = Hn31 + Nfv_g'*Nfv_g*Abohm(3,1)*dline;
%         Hn32 = Hn32 + Nfv_g'*Nfv_g*Abohm(3,2)*dline;
%         Hn33 = Hn33 + Nfv_g'*Nfv_g*Abohm(3,3)*dline;
%         Hn34 = Hn34 + Nfv_g'*Nfv_g*Abohm(3,4)*dline;
%         
%         Hn41 = Hn41 + Nfv_g'*Nfv_g*Abohm(4,1)*dline;
%         Hn42 = Hn42 + Nfv_g'*Nfv_g*Abohm(4,2)*dline;
%         Hn43 = Hn43 + Nfv_g'*Nfv_g*Abohm(4,3)*dline;
%         Hn44 = Hn44 + Nfv_g'*Nfv_g*Abohm(4,4)*dline;
 
        
        
        Hn31 = Hn31 + Abohm(3,1)*NN;
        Hn32 = Hn32 + Abohm(3,2)*NN;
        Hn33 = Hn33 + Abohm(3,3)*NN;
        Hn34 = Hn34 + Abohm(3,4)*NN;
        
        Hn41 = Hn41 + Abohm(4,1)*NN;
        Hn42 = Hn42 + Abohm(4,2)*NN;
        Hn43 = Hn43 + Abohm(4,3)*NN;
        Hn44 = Hn44 + Abohm(4,4)*NN;        
        
        
% if iel==202
%     stop
% end
        % Contribution of the current integration point to the elemental matrix
        TUh31 = TUh31 + coefi*(gammai*dAlphai_dU(1)+Alphai*(Taui(1,1)*b(1)+Taui(2,1)*b(2)))*NN;
        TUh32 = TUh32 + coefi*(gammai*dAlphai_dU(2)+Alphai*(Taui(1,2)*b(1)+Taui(2,2)*b(2)))*NN;
        TUh33 = TUh33 + coefi*(gammai*dAlphai_dU(3)+Alphai*(Taui(1,3)*b(1)+Taui(2,3)*b(2)))*NN;
        TUh34 = TUh34 + coefi*(gammai*dAlphai_dU(4)+Alphai*(Taui(1,4)*b(1)+Taui(2,4)*b(2)))*NN;
        TUh41 = TUh41 + coefe*(gammae*dAlphae_dU(1)+Alphae*(Taue(1,1)*b(1)+Taue(2,1)*b(2)))*NN;
        TUh42 = TUh42 + coefe*(gammae*dAlphae_dU(2)+Alphae*(Taue(1,2)*b(1)+Taue(2,2)*b(2)))*NN;
        TUh43 = TUh43 + coefe*(gammae*dAlphae_dU(3)+Alphae*(Taue(1,3)*b(1)+Taue(2,3)*b(2)))*NN;
        TUh44 = TUh44 + coefe*(gammae*dAlphae_dU(4)+Alphae*(Taue(1,4)*b(1)+Taue(2,4)*b(2)))*NN;
        
        TQh31 = TQh31 + coefi*Alphai*Vi(1)*b(1)*NN;
        TQh32 = TQh32 + coefi*Alphai*Vi(1)*b(2)*NN;
        TQh33 = TQh33 + coefi*Alphai*Vi(2)*b(1)*NN;
        TQh34 = TQh34 + coefi*Alphai*Vi(2)*b(2)*NN;
        TQh35 = TQh35 + coefi*Alphai*Vi(3)*b(1)*NN;
        TQh36 = TQh36 + coefi*Alphai*Vi(3)*b(2)*NN;
        TQh37 = TQh37 + coefi*Alphai*Vi(4)*b(1)*NN;
        TQh38 = TQh38 + coefi*Alphai*Vi(4)*b(2)*NN;
        TQh41 = TQh41 + coefe*Alphae*Ve(1)*b(1)*NN;
        TQh42 = TQh42 + coefe*Alphae*Ve(1)*b(2)*NN;
        TQh43 = TQh43 + coefe*Alphae*Ve(2)*b(1)*NN;
        TQh44 = TQh44 + coefe*Alphae*Ve(2)*b(2)*NN;
        TQh45 = TQh45 + coefe*Alphae*Ve(3)*b(1)*NN;
        TQh46 = TQh46 + coefe*Alphae*Ve(3)*b(2)*NN;
        TQh47 = TQh47 + coefe*Alphae*Ve(4)*b(1)*NN;
        TQh48 = TQh48 + coefe*Alphae*Ve(4)*b(2)*NN;
        
        Tfh3    = Tfh3 + coefi*Alphai*sum( dot(Taui(1,:),sol_g(g,:))*b(1)+dot(Taui(2,:),sol_g(g,:))*b(2)   )*NNf;
        Tfh4    = Tfh4 + coefe*Alphae*sum( dot(Taue(1,:),sol_g(g,:))*b(1)+dot(Taue(2,:),sol_g(g,:))*b(2)   )*NNf;
    end
    
    Wm = BohmMatrices(delta,delta_rho, sn,soundSpeed);
    
    
    
    Mf = Mf + (Nfv_g'*Nfv_g)*dline;
    B11 = B11 + Wm(1,1)*(Nfv_g'*Nfv_g)*dline;
    B12 = B12 + Wm(1,2)*(Nfv_g'*Nfv_g)*dline;
    B21 = B21 + Wm(2,1)*(Nfv_g'*Nfv_g)*dline;
    B22 = B22 + Wm(2,2)*(Nfv_g'*Nfv_g)*dline;
    fH = fH + Nfv_g'*[delta_rho*limitRhoMin,0,-fan(g,3),-fan(g,4)]*dline;
    
    
        Cnx = Cnx + Nfv_g'*Nfv_g*n_g(1)*dline;
        Cny = Cny + Nfv_g'*Nfv_g*n_g(2)*dline;    
        Qbx = Qbx + b(1)*Nfv_g'*(b(1)*Nfv_g*n_g(1)+b(2)*Nfv_g*n_g(2))*dline;
        Qby = Qby + b(2)*Nfv_g'*(b(1)*Nfv_g*n_g(1)+b(2)*Nfv_g*n_g(2))*dline;            
    
    
    
    
end


% expand the matrices
% Df_loc = expandMatrixCv(mult(1)*B11,mult(1)*B12,B13,B14,mult(2)*B21,mult(2)*B22,B23,B24,B31,B32,mult(3)*Mf,B34,B41,B42,B43,mult(4)*Mf);
% Ef_loc = expandMatrixA_mult(Mf,neq,mult);


Df_loc = expandMatrixCv(B11,B12,B13,B14,B21,B22,B23,B24,B31,B32,Mf,B34,B41,...
                       B42,B43,Mf);
Ef_loc = expandMatrixCv(Mf,zer12,zer13,zer14,zer21,Mf,zer23,zer24,zer31,...
    zer32,Mf,zer34,zer41,zer42,zer43,Mf);


Hloc = expandMatrixCv(Hn11,Hn12,Hn13,Hn14,Hn21,Hn22,Hn23,Hn24,Hn31,Hn32,Hn33,Hn34,Hn41,Hn42,Hn43,Hn44);
TUhloc = expandMatrixCv(TUh11,TUh12,TUh13,TUh14,TUh21,TUh22,TUh23,TUh24,TUh31,TUh32,TUh33,TUh34,TUh41,TUh42,TUh43,TUh44);
TQhloc = expandMatrixTQ( TQh11,TQh12,TQh13,TQh14,TQh15,TQh16,TQh17,TQh18, ...
    TQh21,TQh22,TQh23,TQh24,TQh25,TQh26,TQh27,TQh28, ...
    TQh31,TQh32,TQh33,TQh34,TQh35,TQh36,TQh37,TQh38,...
    TQh41,TQh42,TQh43,TQh44,TQh45,TQh46,TQh47,TQh48);




% if iel==202
%     stop
% end


% % elemental assembly
% Df(ind_face_2,ind2) = Df_loc;
% Ef(ind_face_2,ind_face_2) = Ef_loc;
% Qf(ind_face_2,ind4) = 0;
% fHf(ind_face_2) = -col(fH')*iota;
% Hf(ind_face_2,ind_face_2)     = iota*Hloc;
% TUhf(ind_face_2,ind_face_2) = iota*TUhloc;
% TQhf(ind_face_2,ind4)           = iota*TQhloc;
% Tfhf(ind_face_2)                    = iota*col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));

% elemental assembly
Df(ind_face_2,ind2) = Df_loc;
Ef(ind_face_2,ind_face_2) = Ef_loc;
Qf(ind_face_2,ind4) = 0;
fHf(ind_face_2) = -col(fH');
Hf(ind_face_2,ind_face_2)     = Hloc;
TUhf(ind_face_2,ind_face_2) = TUhloc;
TQhf(ind_face_2,ind4)           = TQhloc;
Tfhf(ind_face_2)                    = col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));


%% NEW*********************************************
Lf(ind_face_2,ind4) =  transpose( expandMatrixB(Cnx',Cny') );
Qf(ind_face_2,ind4) =  transpose( expandMatrixB(Qbx',Qby') );
if iota
    diffdiag = repmat(col(repmat([diff_n diff_u diff_ei diff_ee],nd,1)),nv,1);
    Diff = diag(diffdiag);
    Lf(ind_face_2,ind4) = Lf(ind_face_2,ind4)*Diff;
    Qf(ind_face_2,ind4) = Qf(ind_face_2,ind4)*Diff;
end
%% NEW*********************************************


% if iel==143
%     stop
% end
%% Neumann bc
if useNeumannBC==1
    Hf(ind_face_2,ind_face_2)   = 0;
    TUhf(ind_face_2,ind_face_2) = 0;
    TQhf(ind_face_2,ind4)       = 0;
    Tfhf(ind_face_2)            = 0;
    fHf(ind_face_2)             = 0;
    Qf(ind_face_2,ind4)         = 0;
end













%% Elemental matrices Bohm
function [Df,Ef,Hf,Lf,Qf,fHf,TUhf,TQhf,Tfhf] = elementalMatricesBohmSimp(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,X,refEl,ifac,solf,solQf,iFace,iel)


% mesh data
global neq Magnetic testcase useThreshold limitRhoMin axisym epn diff_pari diff_pare Mref useNeumannBC decoupleEquations tau

mult = [1 1 1 1];

coefi = (2/(3*Mref) )^(1+epn)*diff_pari;
coefe = (2/(3*Mref) )^(1+epn)*diff_pare;



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
solQ_g = N1dv*solQf;

sol_g_phys = cons2phys(sol_g);

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
indLdens = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:2)'),nd*nv,1); % assembly face to elem for unknown gradient (only dens)
indLmome = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(3:4)'),nd*nv,1); % assembly face to elem for unknown gradient (only mom)
indLener = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(5:6)'),nd*nv,1); % assembly face to elem for unknown gradient (only energy)
indLTiTe = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(5:8)'),2*nd*nv,1); % assembly face to elem for unknown gradient (only TiTe)
indLTi = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(5:6)'),nd*nv,1); % assembly face to elem for unknown gradient (only Ti)
indLTe = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(7:8)'),nd*nv,1); % assembly face to elem for unknown gradient (only Te)



% Initialization
B11 = zeros(nv); B12 = B11; B13 = B11; B14 = B11;
B21 = B11; B22 = B11; B23 = B11; B24 = B11;
B31 = B11; B32 = B11; B33 = B11; B34 = B11;
B41 = B11; B42 = B11; B43 = B11; B44 = B11;
Mf = B11;
fHf = zeros(neq*nf*nv,1);
fH = zeros(nv,neq);
if useThreshold
    sol_g(sol_g(:,1)<useThreshold,1) = useThreshold;
end
ug = sol_g(:,2)./sol_g(:,1);




TUh11 = zeros(nv); TUh12 = TUh11; TUh13 = TUh11; TUh14 = TUh11;
TUh21 = TUh11; TUh22 = TUh11; TUh23 = TUh11; TUh24 = TUh11;
TUh31 = TUh11; TUh32 = TUh11; TUh33 = TUh11; TUh34 = TUh11;
TUh41 = TUh11; TUh42 = TUh11; TUh43 = TUh11; TUh44 = TUh11;

TQh11 = zeros(nv); TQh12 = TQh11; TQh13 = TQh11; TQh14 = TQh11; TQh15 = TQh11; TQh16 = TQh11; TQh17 = TQh11;  TQh18 = TQh11;
TQh21 = TQh11; TQh22 = TQh11; TQh23 = TQh11; TQh24 = TQh11; TQh25 = TQh11; TQh26 = TQh11; TQh27 = TQh11;  TQh28 = TQh11;
TQh31 = TQh11; TQh32 = TQh11; TQh33 = TQh11; TQh34 = TQh11; TQh35 = TQh11; TQh36 = TQh11; TQh37 = TQh11;  TQh38 = TQh11;
TQh41 = TQh11; TQh42 = TQh11; TQh43 = TQh11; TQh44 = TQh11; TQh45 = TQh11; TQh46 = TQh11; TQh47 = TQh11;  TQh48 = TQh11;

Tfh1          = zeros(nv,1); Tfh2 = Tfh1; Tfh3 = Tfh1; Tfh4 = Tfh1;

% analytical rhs for convergence test
fan = rhsbohm(xyfg);

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
    
    soundSpeed = sol_g_phys(g,9) ;
    if ~isreal(soundSpeed)
        soundSpeed = abs(soundSpeed);
    end
    
    if decoupleEquations
        soundSpeed = sqrt(Mref);
    end
    
    
    if abs(dot(n_g,b/norm(b)))<1e-1
        sn = 1;
        delta = 0;
        iota  =0;
    else
        sn = sign(dot(n_g,b));
        delta = sn*ug(g)<soundSpeed;
        iota = 1;
    end
    
    if limitRhoMin
        delta_rho = sol_g(g,1)<limitRhoMin;
    else
        delta_rho = 0;
    end
    
    %% delta=1 impose the velocity to the speed of sound
    %         delta=0;
    
    
    % Compute G^T^(k-1)
    G = reshape(solQ_g(g,:),nd,neq);
    
    % Compute V(U^(k-1))
    Vi = computeVi(sol_g(g,:));
    Ve = computeVe(sol_g(g,:));
    
    % Compute dV_dU (k-1)
    dVi_dU = compute_dVi_dU(sol_g(g,:));
    dVe_dU = compute_dVe_dU(sol_g(g,:));
    
    
    NNf =  Nfv_g'*dline;
    NN  =  Nfv_g'*Nfv_g*dline;
    
    Taui = G*dVi_dU;      % 2x3
    Taue = G*dVe_dU;      % 2x3
    
    
%     Vi = [0 0 1 0];
%     Ve = [0 0 0 1];
%     b = b/norm(b);
%     b = n_g;
    
    if iota
        % Contribution of the current integration point to the elemental matrix
        TUh31 = TUh31 + ((Taui(1,1)*b(1)+Taui(2,1)*b(2)))*NN;
        TUh32 = TUh32 + ((Taui(1,2)*b(1)+Taui(2,2)*b(2)))*NN;
        TUh33 = TUh33 + ((Taui(1,3)*b(1)+Taui(2,3)*b(2)))*NN;
        TUh34 = TUh34 + ((Taui(1,4)*b(1)+Taui(2,4)*b(2)))*NN;
        TUh41 = TUh41 + ((Taue(1,1)*b(1)+Taue(2,1)*b(2)))*NN;
        TUh42 = TUh42 + ((Taue(1,2)*b(1)+Taue(2,2)*b(2)))*NN;
        TUh43 = TUh43 + ((Taue(1,3)*b(1)+Taue(2,3)*b(2)))*NN;
        TUh44 = TUh44 + ((Taue(1,4)*b(1)+Taue(2,4)*b(2)))*NN;
                
        TQh31 = TQh31 + Vi(1)*b(1)*NN;
        TQh32 = TQh32 + Vi(1)*b(2)*NN;
        TQh33 = TQh33 + Vi(2)*b(1)*NN;
        TQh34 = TQh34 + Vi(2)*b(2)*NN;
        TQh35 = TQh35 + Vi(3)*b(1)*NN;
        TQh36 = TQh36 + Vi(3)*b(2)*NN;
        TQh37 = TQh37 + Vi(4)*b(1)*NN;
        TQh38 = TQh38 + Vi(4)*b(2)*NN;
        TQh41 = TQh41 + Ve(1)*b(1)*NN;
        TQh42 = TQh42 + Ve(1)*b(2)*NN;
        TQh43 = TQh43 + Ve(2)*b(1)*NN;
        TQh44 = TQh44 + Ve(2)*b(2)*NN;
        TQh45 = TQh45 + Ve(3)*b(1)*NN;
        TQh46 = TQh46 + Ve(3)*b(2)*NN;
        TQh47 = TQh47 + Ve(4)*b(1)*NN;
        TQh48 = TQh48 + Ve(4)*b(2)*NN;
        
        Tfh3    = Tfh3 + sum( dot(Taui(1,:),sol_g(g,:))*b(1)+dot(Taui(2,:),sol_g(g,:))*b(2)   )*NNf;
        Tfh4    = Tfh4 + sum( dot(Taue(1,:),sol_g(g,:))*b(1)+dot(Taue(2,:),sol_g(g,:))*b(2)   )*NNf;
        
    end
    
    
    Wm = BohmMatrices(delta,delta_rho, sn,soundSpeed);
    
    
    
    Mf = Mf + (Nfv_g'*Nfv_g)*dline;
    B11 = B11 + Wm(1,1)*(Nfv_g'*Nfv_g)*dline;
    B12 = B12 + Wm(1,2)*(Nfv_g'*Nfv_g)*dline;
    B21 = B21 + Wm(2,1)*(Nfv_g'*Nfv_g)*dline;
    B22 = B22 + Wm(2,2)*(Nfv_g'*Nfv_g)*dline;
    fH = fH + Nfv_g'*[delta_rho*limitRhoMin,0,-fan(g,3),-fan(g,4)]*dline;
end


% expand the matrices
Df_loc = expandMatrixCv(mult(1)*B11,mult(1)*B12,B13,B14,mult(2)*B21,mult(2)*B22,B23,B24,B31,B32,mult(3)*Mf,B34,B41,B42,B43,mult(4)*Mf);
Ef_loc = expandMatrixA_mult(Mf,neq,mult);


% Df_loc = expandMatrixCv(B11,B12,B13,B14,B21,B22,B23,B24,B31,B32,B33,B34,B41,B42,B43,B44);
% Ef_loc = expandMatrixCv(Mf,zer12,zer13,zer14,zer21,Mf,zer23,zer24,zer31,zer32,zer33,zer34,zer41,zer42,zer43,zer44);


TUhloc = expandMatrixCv(TUh11,TUh12,TUh13,TUh14,TUh21,TUh22,TUh23,TUh24,TUh31,TUh32,TUh33,TUh34,TUh41,TUh42,TUh43,TUh44);
TQhloc = expandMatrixTQ( TQh11,TQh12,TQh13,TQh14,TQh15,TQh16,TQh17,TQh18, ...
    TQh21,TQh22,TQh23,TQh24,TQh25,TQh26,TQh27,TQh28, ...
    TQh31,TQh32,TQh33,TQh34,TQh35,TQh36,TQh37,TQh38,...
    TQh41,TQh42,TQh43,TQh44,TQh45,TQh46,TQh47,TQh48);


% % elemental assembly
% Df(ind_face_2,ind2) = Df_loc;
% Ef(ind_face_2,ind_face_2) = Ef_loc;
% Qf(ind_face_2,ind4) = 0;
% fHf(ind_face_2) = -col(fH')*iota;
% Hf(ind_face_2,ind_face_2)     = iota*Hloc;
% TUhf(ind_face_2,ind_face_2) = iota*TUhloc;
% TQhf(ind_face_2,ind4)           = iota*TQhloc;
% Tfhf(ind_face_2)                    = iota*col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));

% elemental assembly
Df(ind_face_2,ind2)              = Df_loc;
Ef(ind_face_2,ind_face_2)     = Ef_loc;
fHf(ind_face_2)                     = 0;
Hf(ind_face_2,ind_face_2)     = 0;
TUhf(ind_face_2,ind_face_2) = TUhloc;
TQhf(ind_face_2,ind4)           = TQhloc;
Tfhf(ind_face_2)                    = col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));
Lf(ind_face_2,indLTiTe)         = 0;
Qf(ind_face_2,ind4)               = 0;


% prova1
% TQhf(ind_face_2,ind4)          = 0;
% TUhf(ind_face_2,ind_face_2) = 0;
% Tfhf(ind_face_2)                    = 0;

%prova2
% TQhf(ind_face_2,ind4) = 0;
% TUhf(ind_face_2,ind_face_2) = 0;
% Tfhf(ind_face_2)                    = 0;


% % TQhf(ind_face_2,ind4)           = 0;
% % 
% % 

% % 
% % TQhf(ind_face_2,ind4) = Lf(ind_face_2,ind4) ;
% Lf(ind_face_2,ind4) = 0;












function [Df,Ef,Hf,Lf,Qf,fHf,TUhf,TQhf,Tfhf] = elementalMatricesBohmTest(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,X,refEl,ifac,solf,solQf,iFace,iel)


% mesh data
global neq Magnetic testcase useThreshold limitRhoMin axisym epn diff_pari diff_pare Mref useNeumannBC decoupleEquations tau

mult = [1 1 1 1];




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
solQ_g = N1dv*solQf;

sol_g_phys = cons2phys(sol_g);

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
indLdens = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:2)'),nd*nv,1); % assembly face to elem for unknown gradient (only dens)
indLmome = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(3:4)'),nd*nv,1); % assembly face to elem for unknown gradient (only mom)
indLener = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(5:6)'),nd*nv,1); % assembly face to elem for unknown gradient (only energy)
indLTiTe = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(5:8)'),2*nd*nv,1); % assembly face to elem for unknown gradient (only TiTe)
indLTi = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(5:6)'),nd*nv,1); % assembly face to elem for unknown gradient (only Ti)
indLTe = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(7:8)'),nd*nv,1); % assembly face to elem for unknown gradient (only Te)



% Initialization
B11 = zeros(nv); B12 = B11; B13 = B11; B14 = B11;
B21 = B11; B22 = B11; B23 = B11; B24 = B11;
B31 = B11; B32 = B11; B33 = B11; B34 = B11;
B41 = B11; B42 = B11; B43 = B11; B44 = B11;
Mf = B11;

Cnx = zeros(nv);
Cny = Cnx;
    
fHf = zeros(neq*nf*nv,1);
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
    
    soundSpeed = sol_g_phys(g,9) ;
    if ~isreal(soundSpeed)
        soundSpeed = abs(soundSpeed);
    end
    
    if decoupleEquations
        soundSpeed = sqrt(Mref);
    end
    
    
    if abs(dot(n_g,b/norm(b)))<1e-1
        sn = 1;
        delta = 0;
        iota  =0;
    else
        sn = sign(dot(n_g,b));
        delta = sn*ug(g)<soundSpeed;
        iota = 1;
    end
    
      
    if limitRhoMin
        delta_rho = sol_g(g,1)<limitRhoMin;
    else
        delta_rho = 0;
    end
    
    Wm = BohmMatrices(delta,delta_rho, sn,soundSpeed);
    
%     
%     Cnx = Cnx + Nfv_g'*Nfv_g*n_g(1)*dline;
%     Cny = Cny + Nfv_g'*Nfv_g*n_g(2)*dline;

if ismember(iel, [204,202,208,206,104])
    
    
    Cnx = Cnx + Nfv_g'*Nfv_g*(-n_g(1))*dline;
    Cny = Cny + Nfv_g'*Nfv_g*n_g(2)*dline;
else

    Cnx = Cnx + Nfv_g'*Nfv_g*(n_g(1))*dline;
    Cny = Cny + Nfv_g'*Nfv_g*n_g(2)*dline;
end
%     
%     Cnx = Cnx + Nfv_g'*Nfv_g*(b(1)/norm(b))*dline;
%     Cny = Cny + Nfv_g'*Nfv_g*(b(2)/norm(b))*dline;

    
    
    Mf = Mf + (Nfv_g'*Nfv_g)*dline;
    B11 = B11 + Wm(1,1)*(Nfv_g'*Nfv_g)*dline;
    B12 = B12 + Wm(1,2)*(Nfv_g'*Nfv_g)*dline;
    B21 = B21 + Wm(2,1)*(Nfv_g'*Nfv_g)*dline;
    B22 = B22 + Wm(2,2)*(Nfv_g'*Nfv_g)*dline;
end


% expand the matrices
Df_loc = expandMatrixCv(mult(1)*B11,mult(1)*B12,B13,B14,mult(2)*B21,mult(2)*B22,B23,B24,B31,B32,mult(3)*Mf,B34,B41,B42,B43,mult(4)*Mf);
Ef_loc = expandMatrixA_mult(Mf,neq,mult);
aux =10*transpose( expandMatrixB(Cnx',Cny') );
aa = size(aux,2);
ind = col([5:8:aa-3;6:8:aa-2;7:8:aa-1;8:8:aa]);
Lf(ind_face_2,indLTiTe) = aux(:,ind);

% elemental assembly
Df(ind_face_2,ind2)              = Df_loc;
Ef(ind_face_2,ind_face_2)     = Ef_loc;
fHf(ind_face_2)                     =0;
Hf(ind_face_2,ind_face_2)     = 0;
TUhf(ind_face_2,ind_face_2) = 0;
TQhf(ind_face_2,ind4)           = 0;
Tfhf(ind_face_2)                    = 0;
Qf(ind_face_2,ind4)               = 0;






%% Elemental matrices
function [Df,Ef,Hf,Lf,Qf,fHf,TUhf,TQhf,Tfhf] = elementalMatricesCore(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,X,refEl,ifac,iFace)

% mesh data
global neq testcase coreflux axisym tau

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
B11 = zeros(nv); B12 = B11; B13 = B11; B14 = B11;
B21 = B11; B22 = B11; B23 = B11; B24 = B11;
B31 = B11; B32 = B11; B33 = B11; B34 = B11;
B41 = B11; B42 = B11; B43 = B11; B44 = B11;
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
    B11 = B11 + tau*(Nfv_g'*Nfv_g)*dline;
    fH = fH + tau*Nfv_g'*coreflux'*dline;
    
end

B22 = B11;
B33 = B11;
B44 = B11;
% expand the matrices
E_loc = expandMatrixCv(B11, B12, B13,B14,B21,B22,B23,B24,B31,B32,B33,B34,B41,B42,B43,B44);

% elemental assembly
Df(ind_face_2,ind2) = 0;
Hf(ind_face_2,ind_face_2) = 0;
TUhf(ind_face_2,ind_face_2) = 0;
TQhf(ind_face_2,ind4)           = 0;
Tfhf(ind_face_2)                    = 0;

if (testcase.n==51 || testcase.n==52 || testcase.n==53 || testcase.n==54 || testcase.n==55)
    Ef(ind_face_2,ind_face_2) = 0;
    Qf(ind_face_2,ind4) = 0;
else
    
    Ef(ind_face_2,ind_face_2) = E_loc;
    Lf(ind_face_2,ind4) = 0;
    Qf(ind_face_2,ind4) = 0;
    fHf(ind_face_2) = -col(fH');
    
    % Ef_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/Ef.txt');
    % fh_fort= readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/fh.txt');
    % max(abs(Ef_fort(:)-E_loc(:)))
    % max(abs(fh_fort(:)-  (-col(fH') ) ))
    % stop
    
    
end





%
%
%
% %% Elemental matrices
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
% ind_face_2 = (ifac-1)*neq*nv + (1:neq*nv);
% ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
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
% Df(ind_face_2,ind2) = 0;
% Ef(ind_face_2,ind_face_2) = E_loc;
% Lf(ind_face_2,ind4) = 0;
% Hf(ind_face_2,ind_face_2) = 0;
% fHf(ind_face_2) = -col(fH');
%






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


function res = expandMatrixA_mult(A,n,mult)
% expand matrix A and M
%  [ A 0 0 0
%    0 A 0 0
%    0 0 A 0
%    0 0 0 A ]
% dimension n
res = zeros([size(A) n n]);
aux = repmat(A, [1 1 n]);
aux(:,:,1) = aux(:,:,1)*mult(1);
aux(:,:,2) = aux(:,:,2)*mult(2);
aux(:,:,3) = aux(:,:,3)*mult(3);
aux(:,:,4) = aux(:,:,4)*mult(4);
res(:,:,1:n+1:n^2) = aux;
res = permute(res, [3 1 4 2]);
res = reshape(res, n*size(A));


% function res = expandMatrixCv(Cxx,Cxy,Cyx,Cyy)
% % expand matrix Cv
% %   [ Cxx Cxy
% %     Cyx Cyy]
% res = zeros([size(Cxx) 2 2]);
% res(:,:,[1 3 2 4]) = cat(3,Cxx,Cxy,Cyx,Cyy);
% res = permute(res, [3 1 4 2]);
% res = reshape(res, 2*size(Cxx,1),2*size(Cxx,2));




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
%     if U(1)<useThreshold
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

function B = BohmMatrices(delta,delta_rho,sn,soundSpeed)

B = [1-delta_rho, 0; delta*sn*soundSpeed, (1-delta)];



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


function res = expandMatrixTQ(C11,C12,C13,C14,C15,C16,C17,C18,...
    C21,C22,C23,C24,C25,C26,C27,C28,...
    C31,C32,C33,C34,C35,C36,C37,C38,...
    C41,C42,C43,C44,C45,C46,C47,C48)
% expand matrix TQ
%   [C11,C12,C13,C14,C15,C16,C17,C18,...
%    C21,C22,C23,C24,C25,C26,C27,C28,...
%    C31,C32,C33,C34,C35,C36,C37,C38,...
%    C41,C42,C43,C44,C45,C46,C47,C48 ]

res = zeros([size(C11) 4 8]);
res(:,:,[1 5  9  13 17 21 25 29 ...
    2 6 10 14 18 22 26 30 ...
    3 7 11 15 19 23 27 31 ...
    4 8 12 16 20 24 28  32 ]) = ...
    cat(3,C11,C12,C13,C14,C15,C16,C17,C18,...
    C21,C22,C23,C24,C25,C26,C27,C28,...
    C31,C32,C33,C34,C35,C36,C37,C38,...
    C41,C42,C43,C44,C45,C46,C47,C48);
res = permute(res, [3 1 4 2]);
res = reshape(res, 4*size(C11,1),8*size(C11,2));



function Jn = jacobianMatricesBohm(U)

global   useThreshold Gmbohmi Gmbohme
rho = U(1);
Gamma = U(2);

if useThreshold
    if U(1)<useThreshold
        rho = useThreshold;
        Gamma = 0;
        
    end
end
auxi = (5-2*Gmbohmi)/3;
auxe =  (5-2*Gmbohme)/3;
Jn31 = -auxi*(U(2)^3/U(1)^3-U(2)*U(3)/U(1)^2);
Jn32 =  auxi*(U(3)/U(1)-3/2*U(2)^2/U(1)^2) ;
Jn33 =  auxi*U(2)/U(1);
Jn34 = 0;
Jn41 = -auxe*U(2)*U(4)/U(1)^2;
Jn42 = auxe*U(4)/U(1);
Jn43 = 0;
Jn44 = auxe*U(2)/U(1);
Jn = [  0,          0,        0,         0; ...
    0,          0,        0,         0;...
    Jn31,    Jn32,    Jn33,    Jn34;...
    Jn41,    Jn42,    Jn43,    Jn44];







function res = computeVi(U)
% res:  [4x1]

res = [  U(2)^2/U(1)^3 - U(3)/U(1)^2; ...
    -U(2)/U(1)^2; ...
    1/U(1); ...
    0                             ];


function res = computeVe(U)
% res:  [4x1]

res = [      -U(4)/U(1)^2; ...
    0; ...
    0; ...
    1/U(1)                  ];

function res = compute_dVi_dU(U)
% res: [4 x 4]

res = [2*U(3)/U(1)^3-3*U(2)^2/U(1)^4,           2*U(2)/U(1)^3,               -1/U(1)^2,                  0;...
    2*U(2)/U(1)^3,                                    -1/U(1)^2,                               0,                0;...
    -1/U(1)^2,                                          0,                                          0,                 0;...
    0,                                                   0,                                          0,                 0];

function res = compute_dVe_dU(U)
% res: [4 x 4]

res = [2*U(4)/U(1)^3 ,                                     0,                      0,                  -1/U(1)^2;...
    0,                                                0,                      0,                        0;...
    0,                                                0,                      0,                        0;...
    -1/U(1)^2,                                         0,                      0,                       0];


% function res = computeAlphai(U)
% % res: scalar
% tol = 1e-5;
% global epn
%
% aux = U(3)/U(1) - 0.5*U(2)^2/U(1)^2;
% if aux<0, aux = tol; end
% res = ( aux)^epn;

function res = computeAlphai(U)
% res: scalar
tol = 1e-5;
global epn prolongateExponents

aux = U(3)/U(1) - 0.5*U(2)^2/U(1)^2;

if ~prolongateExponents
    if aux<0, aux = tol; end
    res = ( aux)^epn;
else
    res = prolExp(aux,epn);
end


% function res = computeAlphae(U)
% % res: scalar
% tol = 1e-5;
% global epn
%
% aux = U(4)/U(1);
% if aux<0, aux = tol; end
% res = ( aux)^epn;


function res = computeAlphae(U)
% res: scalar
tol = 1e-5;
global epn prolongateExponents

aux = U(4)/U(1);

if ~prolongateExponents
    if aux<0, aux = tol; end
    res = ( aux)^epn;
else
    res = prolExp(aux,epn);
end

% function res = compute_dAlphai_dU(U)
% % res = [4x1]
% global epn
% tol = 1e-5;
%
% aux = U(3)/U(1) - 0.5*U(2)^2/U(1)^2;
% if aux<0, aux = tol; end
%
% res = epn*(aux)^(epn-1) *...
%     [      -U(3)/U(1)^2+U(2)^2/U(1)^3; ...
%     -U(2)/U(1)^2; ...
%     1/U(1);...
%     0                        ];

function res = compute_dAlphai_dU(U)
% res = [4x1]
global epn prolongateExponents
tol = 1e-5;

aux = U(3)/U(1) - 0.5*U(2)^2/U(1)^2;
if ~prolongateExponents
    if aux<0, aux = tol; end
    res = (aux)^(epn-1);
else
    res = prolExp(aux,epn-1);
end
res = epn*res*...
    [      -U(3)/U(1)^2+U(2)^2/U(1)^3; ...
    -U(2)/U(1)^2; ...
    1/U(1);...
    0                        ];



% function res = compute_dAlphae_dU(U)
% % res = [4x1]
% global epn
% tol = 1e-5;
%
% aux = U(4)/U(1) ;
% if aux<0, aux = tol; end
%
% res = epn*(aux)^(epn-1) *...
%     [      -U(4)/U(1)^2; ...
%     0; ...
%     0;...
%     1/U(1)];

function res = compute_dAlphae_dU(U)
% res = [4x1]
global epn prolongateExponents
tol = 1e-5;
aux = U(4)/U(1) ;
if ~prolongateExponents
    if aux<0, aux = tol; end
    res = (aux)^(epn-1);
else
    res = prolExp(aux,epn-1);
end
res = epn*res *...
    [      -U(4)/U(1)^2; ...
    0; ...
    0;...
    1/U(1)];







function res = expandMatrixB(Bx,By)
% expand matrix B
%   [ Bx  0  0   0 
%     By  0  0   0
%     0  Bx  0   0
%     0  By  0   0
%     0   0  Bx  0
%     0   0  By  0
%     0   0   0  Bx 
%     0   0   0  By ]
res = zeros([size(Bx) 8 4]);
res(:,:,[1 2 11 12 21 22 31 32]) = cat(3,Bx,By,Bx,By,Bx,By,Bx,By);
res = permute(res, [3 1 4 2]);
res = reshape(res, 8*size(Bx,1),4*size(Bx,2));
