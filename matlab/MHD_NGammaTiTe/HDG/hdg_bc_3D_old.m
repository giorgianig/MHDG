function [Df,Ef,Hf,Lf,Qf,fH,TUhf,TQhf,Tfhf] = hdg_bc_3D(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,X,T,infoFaces,refEl,bn,           F,F_dir,sol_hat,solQ)

global neq testcase refElTor useNeumannBC theta ntor
%
Ne = size(Df,3);
nf = size(F,2);
N2d = size(F,1);
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);         % Number of 1d points for toroidal line
Np2d    = size(refEl.NodesCoord,1);
Np = Np2d*Np1dTor;
Nfp     = Np2d*2+nf*Np1dPol*Np1dTor;
Nfl     = Np1dPol*Np1dTor;
tdiv    = linspace(0,theta,ntor+1);
Nf        = max(max(F));                                              % Number of faces in the 2d plane
nDirFaces =  sum(sum(F_dir));
nd = 3;
ind_dim = neq*[Np2d,repmat(Np1dPol*Np1dTor,1,nf),Np2d];
ind_sta = [1, 1+cumsum(ind_dim)];
ind_uhat = zeros(Nfp*neq,1);

% number of boundaries
nb = numel(bn);
fH = zeros(Nfp*neq,Ne);
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
            TUhfe = TUhf(:,:,iElem);
            TQhfe = TQhf(:,:,iElem);
            Tfhfe   =  Tfhf(:,iElem);
            delta = 1+(itor-1)*(N2d*Np2d+(Nf-nDirFaces)*Nfl)*neq+neq*(N2d*Np2d+(Fe-1)*Nfl);
            ind_uhat = delta+(0: (Nfl*neq-1));
            indQ = (iElem-1)*neq*Np*nd + (1:neq*Np*nd);
            solf = sol_hat(ind_uhat);
            solQe = solQ(indQ);
            solQe    = reshape(solQe,neq*nd,numel(solQe)/neq/nd)';
            nodesv = refElTor.faceNodes3(ifac,:);
            solQf    = solQe(nodesv,:);
            
            fHe = fH(:,iElem);
            if any(strcmp(name,{'UP','DOWN','LEFT','RIGHT','WALL','ULIM'}))
                stop
                %             [Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe] = elementalMatricesWall(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,refEl,ifac);
            elseif any(strcmp(name,{'INLET','OUTLET'}))
                stop
                %             [Dfe,Efe,Hfe,Qfe,fHe] = elementalMatricesInlOut(Dfe,Efe,Hfe,Qfe,Xe,refEl,ifac,solf);
            elseif any(strcmp(name,{'Neumann'}))
                stop
                %             [Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe] = elementalMatricesNeumann(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,refEl,ifac);
            elseif any(strcmp(name,{'IN'}))
                [Dfe,Efe,Hfe,Lfe,Qfe,fHe,TUhfe,TQhfe,Tfhfe] = elementalMatricesCore(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,tel,refEl,ifac,Fe);
            elseif any(strcmp(name,{'Dweak'}))
                [Dfe,Efe,Hfe,Lfe,Qfe,fHe,TUhfe,TQhfe,Tfhfe] = elementalMatricesDirichletWeak(Dfe,Efe,Hfe,Lfe,Qfe,Xe,tel,refEl,ifac,Fe);
            elseif any(strcmp(name,{'LIM','OUT'}))
                [Dfe,Efe,Hfe,Lfe,Qfe,fHe,TUhfe,TQhfe,Tfhfe] = elementalMatricesBohm(Dfe,Efe,Hfe,Lfe,Qfe,TUhfe,TQhfe,Tfhfe,Xe,tel,refEl,ifac,solf,solQf,Fe,iElem);
            else
                error('Something wrong')
            end
            Df(:,:,iElem) = Dfe;
            Ef(:,:,iElem) = Efe;
            Hf(:,:,iElem) = Hfe;
            Lf(:,:,iElem) = Lfe;
            Qf(:,:,iElem) = Qfe;
            fH(:,iElem)   = fHe;
            TUhf(:,:,iElem) = TUhfe;
            TQhf(:,:,iElem) = TQhfe;
            Tfhf(:,iElem)    = Tfhfe;
            
            
            
        end
    end
end



%% Elemental matrices Bohm
function [Df,Ef,Hf,Lf,Qf,fHf,TUhf,TQhf,Tfhf] = elementalMatricesBohm(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,Xel,tel,refEl,ifac,solf,solQf,iFace,iel)


% mesh data
global neq Magnetic testcase useThreshold limitRhoMin axisym epn diff_pari diff_pare Mref useNeumannBC decoupleEquations tau
global refElTor

mult = [1 1 1 1];

coefi = (2/(3*Mref) )^(1+epn)*diff_pari;
coefe = (2/(3*Mref) )^(1+epn)*diff_pare;



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
solQ_g = Nfi*solQf;

sol_g_phys = cons2phys(sol_g);

%% FACES COMPUTATIONS:
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
%
% % indices for local assembly
% ind_face_2 = (ifac-1)*neq*nv + (1:neq*nv);
% ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for velocity
% ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for unknown gradient
% indLdens = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:2)'),nd*nv,1); % assembly face to elem for unknown gradient (only dens)
% indLmome = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(3:4)'),nd*nv,1); % assembly face to elem for unknown gradient (only mom)
% indLener = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(5:6)'),nd*nv,1); % assembly face to elem for unknown gradient (only energy)

% Initialization
zz = zeros(nv);

B11 = zz; B12 = zz; B13 = zz; B14 = zz;
B21 = zz; B22 = zz; B23 = zz; B24 = zz;
B31 = zz; B32 = zz; B33 = zz; B34 = zz;
B41 = zz; B42 = zz; B43 = zz; B44 = zz;
Mf = zz;
fHf = zeros(neq*Nfp,1);
fH = zeros(nv,neq);
if useThreshold
    sol_g(sol_g(:,1)<useThreshold,1) = useThreshold;
end
ug = sol_g(:,2)./sol_g(:,1);


Hn11 = zz; Hn12 = zz; Hn13 = zz; Hn14 = zz;
Hn21 = zz; Hn22 = zz; Hn23 = zz; Hn24 = zz;
Hn31 = zz; Hn32 = zz; Hn33 = zz; Hn34 = zz;
Hn41 = zz; Hn42 = zz; Hn43 = zz; Hn44 = zz;


TUh11 = zz; TUh12 = zz; TUh13 = zz; TUh14 = zz;
TUh21 = zz; TUh22 = zz; TUh23 = zz; TUh24 = zz;
TUh31 = zz; TUh32 = zz; TUh33 = zz; TUh34 = zz;
TUh41 = zz; TUh42 = zz; TUh43 = zz; TUh44 = zz;

TQh1_01 = zz; TQh1_02 = zz; TQh1_03 = zz; TQh1_04 = zz; TQh1_05 = zz; TQh1_06 = zz; TQh1_07 = zz; TQh1_08 = zz; TQh1_09 = zz; TQh1_10 = zz; TQh1_11 = zz; TQh1_12 = zz;
TQh2_01 = zz; TQh2_02 = zz; TQh2_03 = zz; TQh2_04 = zz; TQh2_05 = zz; TQh2_06 = zz; TQh2_07 = zz; TQh2_08 = zz; TQh2_09 = zz; TQh2_10 = zz; TQh2_11 = zz; TQh2_12 = zz;
TQh3_01 = zz; TQh3_02 = zz; TQh3_03 = zz; TQh3_04 = zz; TQh3_05 = zz; TQh3_06 = zz; TQh3_07 = zz; TQh3_08 = zz; TQh3_09 = zz; TQh3_10 = zz; TQh3_11 = zz; TQh3_12 = zz;
TQh4_01 = zz; TQh4_02 = zz; TQh4_03 = zz; TQh4_04 = zz; TQh4_05 = zz; TQh4_06 = zz; TQh4_07 = zz; TQh4_08 = zz; TQh4_09 = zz; TQh4_10 = zz; TQh4_11 = zz; TQh4_12 = zz;

Tfh1 = zeros(nv,1); Tfh2 = Tfh1; Tfh3 = Tfh1; Tfh4 = Tfh1;
zer11 = zz; zer12 = zz; zer13 = zz; zer14 = zz;
zer21 = zz; zer22 = zz; zer23 = zz; zer24 = zz;
zer31 = zz; zer32 = zz; zer33 = zz; zer34 = zz;
zer41 = zz; zer42 = zz; zer43 = zz; zer44 = zz;

% analytical rhs for convergence test
fan = rhsbohm(xyfg);

%  LOOP IN GAUSS POINTS
for igtor = 1:ngausstor
    for igpol = 1:ngausspol
        g = (igtor-1)*ngausspol+igpol;
        n = n_g(g,:);
        
							        
        if (testcase.n >= 50 && testcase.n<60)
            b = [Magnetic.bxfaces(g,iFace), Magnetic.byfaces(g,iFace)];
        else
            b = defineMagneticField_3D([xyfg(igpol,1),xyfg(igpol,2)],thetafg(igtor));
        end
        
        soundSpeed = sol_g_phys(g,9) ;
        if ~isreal(soundSpeed)
            soundSpeed = abs(soundSpeed);
        end
        
        if decoupleEquations
            soundSpeed = sqrt(Mref);
        end
        
        if abs(dot(n(1:2),b(1:2)/norm(b(1:2))))<1e-1
            sn = 1;
            delta = 0;
            iota  =0;
        else
            sn = sign(dot(n,b));
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
        G = reshape(solQ_g(g,:),Ndim,neq);
        
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
        
        %
        %      NNf =  Nfv_g'*dot(b,n_g)*dline;
        %      NN  =  Nfv_g'*dot(b,n_g)*Nfv_g*dline;
        
        Nfg = Nfi(g,:);
        NNf =  Nfg'*dsurf(g);
        NN  =  Nfg'*Nfg*dsurf(g);
        
        gammai = dot(G*Vi,b);    % scalar
        gammae = dot(G*Ve,b);    % scalar
        Taui = G*dVi_dU;      % 2x3
        Taue = G*dVe_dU;      % 2x3
        
        Abohm = jacobianMatricesBohm(sol_g(g,:));
        
        
        if iota
            Hn31 = Hn31 + NN*Abohm(3,1);
            Hn32 = Hn32 + NN*Abohm(3,2);
            Hn33 = Hn33 + NN*Abohm(3,3);
            Hn34 = Hn34 + NN*Abohm(3,4);
            
            Hn41 = Hn41 + NN*Abohm(4,1);
            Hn42 = Hn42 + NN*Abohm(4,2);
            Hn43 = Hn43 + NN*Abohm(4,3);
            Hn44 = Hn44 + NN*Abohm(4,4);
            
            % Contribution of the current integration point to the elemental matrix
            TUh31 = TUh31 + coefi*(gammai*dAlphai_dU(1)+Alphai*(Taui(1,1)*b(1)+Taui(2,1)*b(2)+Taui(3,1)*b(3)))*NN;
            TUh32 = TUh32 + coefi*(gammai*dAlphai_dU(2)+Alphai*(Taui(1,2)*b(1)+Taui(2,2)*b(2)+Taui(3,2)*b(3)))*NN;
            TUh33 = TUh33 + coefi*(gammai*dAlphai_dU(3)+Alphai*(Taui(1,3)*b(1)+Taui(2,3)*b(2)+Taui(3,3)*b(3)))*NN;
            TUh34 = TUh34 + coefi*(gammai*dAlphai_dU(4)+Alphai*(Taui(1,4)*b(1)+Taui(2,4)*b(2)+Taui(3,4)*b(3)))*NN;
            TUh41 = TUh41 + coefe*(gammae*dAlphae_dU(1)+Alphae*(Taue(1,1)*b(1)+Taue(2,1)*b(2)+Taue(3,1)*b(3)))*NN;
            TUh42 = TUh42 + coefe*(gammae*dAlphae_dU(2)+Alphae*(Taue(1,2)*b(1)+Taue(2,2)*b(2)+Taue(3,2)*b(3)))*NN;
            TUh43 = TUh43 + coefe*(gammae*dAlphae_dU(3)+Alphae*(Taue(1,3)*b(1)+Taue(2,3)*b(2)+Taue(3,3)*b(3)))*NN;
            TUh44 = TUh44 + coefe*(gammae*dAlphae_dU(4)+Alphae*(Taue(1,4)*b(1)+Taue(2,4)*b(2)+Taue(3,4)*b(3)))*NN;
            
            
            TQh3_01 = TQh3_01 + (coefi*Alphai*Vi(1)*NN )*b(1);
            TQh3_02 = TQh3_02 + (coefi*Alphai*Vi(1)*NN )*b(2);
            TQh3_03 = TQh3_03 + (coefi*Alphai*Vi(1)*NN )*b(3);
            TQh3_04 = TQh3_04 + (coefi*Alphai*Vi(2)*NN )*b(1);
            TQh3_05 = TQh3_05 + (coefi*Alphai*Vi(2)*NN )*b(2);
            TQh3_06 = TQh3_06 + (coefi*Alphai*Vi(2)*NN )*b(3);
            TQh3_07 = TQh3_07 + (coefi*Alphai*Vi(3)*NN )*b(1);
            TQh3_08 = TQh3_08 + (coefi*Alphai*Vi(3)*NN )*b(2);
            TQh3_09 = TQh3_09 + (coefi*Alphai*Vi(3)*NN )*b(3);
            TQh3_10 = TQh3_10 + (coefi*Alphai*Vi(4)*NN )*b(1);
            TQh3_11 = TQh3_11 + (coefi*Alphai*Vi(4)*NN )*b(2);
            TQh3_12 = TQh3_12 + (coefi*Alphai*Vi(4)*NN )*b(3);
            TQh4_01 = TQh4_01 + (coefe*Alphae*Ve(1)*NN )*b(1);
            TQh4_02 = TQh4_02 + (coefe*Alphae*Ve(1)*NN )*b(2);
            TQh4_03 = TQh4_03 + (coefe*Alphae*Ve(1)*NN )*b(3);
            TQh4_04 = TQh4_04 + (coefe*Alphae*Ve(2)*NN )*b(1);
            TQh4_05 = TQh4_05 + (coefe*Alphae*Ve(2)*NN )*b(2);
            TQh4_06 = TQh4_06 + (coefe*Alphae*Ve(2)*NN )*b(3);
            TQh4_07 = TQh4_07 + (coefe*Alphae*Ve(3)*NN )*b(1);
            TQh4_08 = TQh4_08 + (coefe*Alphae*Ve(3)*NN )*b(2);
            TQh4_09 = TQh4_09 + (coefe*Alphae*Ve(3)*NN )*b(3);
            TQh4_10 = TQh4_10 + (coefe*Alphae*Ve(4)*NN )*b(1);
            TQh4_11 = TQh4_11 + (coefe*Alphae*Ve(4)*NN )*b(2);
            TQh4_12 = TQh4_12 + (coefe*Alphae*Ve(4)*NN )*b(3);
            
            Tfh3    = Tfh3 + coefi*Alphai*sum( dot(Taui(1,:),sol_g(g,:))*b(1)+dot(Taui(2,:),sol_g(g,:))*b(2)+dot(Taui(3,:),sol_g(g,:))*b(3)    )*NNf;
            Tfh4    = Tfh4 + coefe*Alphae*sum( dot(Taue(1,:),sol_g(g,:))*b(1)+dot(Taue(2,:),sol_g(g,:))*b(2)+dot(Taue(3,:),sol_g(g,:))*b(3)   )*NNf;
        end
        
        Wm = BohmMatrices(delta,delta_rho, sn,soundSpeed);
        
        Mf = Mf + NN;
        B11 = B11 + Wm(1,1)*NN;
        B12 = B12 + Wm(1,2)*NN;
        B21 = B21 + Wm(2,1)*NN;
        B22 = B22 + Wm(2,2)*NN;
%         fH = fH + [delta_rho*limitRhoMin,0,-fan(g,3),-fan(g,4)]*NNf;
    end
end


% expand the matrices
% Df_loc = expandMatrixCv(mult(1)*B11,mult(1)*B12,B13,B14,mult(2)*B21,mult(2)*B22,B23,B24,B31,B32,mult(3)*Mf,B34,B41,B42,B43,mult(4)*Mf);
% Ef_loc = expandMatrixA_mult(Mf,neq,mult);


Df_loc = expandMatrixCv(B11,B12,B13,B14,B21,B22,B23,B24,B31,B32,B33,B34,B41,B42,B43,B44);
Ef_loc = expandMatrixCv(Mf,zer12,zer13,zer14,zer21,Mf,zer23,zer24,zer31,zer32,zer33,zer34,zer41,zer42,zer43,zer44);


Hloc = expandMatrixCv(Hn11,Hn12,Hn13,Hn14,Hn21,Hn22,Hn23,Hn24,Hn31,Hn32,Hn33,Hn34,Hn41,Hn42,Hn43,Hn44);
TUhloc = expandMatrixCv(TUh11,TUh12,TUh13,TUh14,TUh21,TUh22,TUh23,TUh24,TUh31,TUh32,TUh33,TUh34,TUh41,TUh42,TUh43,TUh44);
TQhloc = expandMatrixTQ(TQh1_01,TQh1_02,TQh1_03,TQh1_04,TQh1_05,TQh1_06,TQh1_07,TQh1_08,TQh1_09,TQh1_10,TQh1_11,TQh1_12, ...
    TQh2_01,TQh2_02,TQh2_03,TQh2_04,TQh2_05,TQh2_06,TQh2_07,TQh2_08,TQh2_09,TQh2_10,TQh2_11,TQh2_12, ...
    TQh3_01,TQh3_02,TQh3_03,TQh3_04,TQh3_05,TQh3_06,TQh3_07,TQh3_08,TQh3_09,TQh3_10,TQh3_11,TQh3_12, ...
    TQh4_01,TQh4_02,TQh4_03,TQh4_04,TQh4_05,TQh4_06,TQh4_07,TQh4_08,TQh4_09,TQh4_10,TQh4_11,TQh4_12);


% % elemental assembly
% Df(ind_face_2,ind2) = Df_loc;
% Ef(ind_face_2,ind_face_2) = Ef_loc;
% Qf(ind_face_2,ind4) = 0;
% fHf(ind_face_2) = -col(fH')*iota;xyfg        = refEl.N1d*Xel(nodes2,:);      % Gauss points coordinates of the face (xy component)
% thetafg     = refElTor.N1d*tel;               % Gauss points coordinates of the face (theta component)

% Hf(ind_face_2,ind_face_2)     = iota*Hloc;
% TUhf(ind_face_2,ind_face_2) = iota*TUhloc;
% TQhf(ind_face_2,ind4)           = iota*TQhloc;
% Tfhf(ind_face_2)                    = iota*col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));

% elemental assembly
Df(ind_ff,ind_fe)    = Df_loc;
Ef(ind_ff,ind_ff)    = Ef_loc;
Qf(ind_ff,ind_fq)    = 0;
fHf(ind_ff)          = -col(fH');
Hf(ind_ff,ind_ff)    = Hloc;
TUhf(ind_ff,ind_ff)  = TUhloc;
TQhf(ind_ff,ind_fq)  = TQhloc;
Tfhf(ind_ff)         = col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));

%% Neumann bc
if useNeumannBC==1
    Hf(ind_ff,ind_ff)   = 0;
    TUhf(ind_ff,ind_ff) = 0;
    TQhf(ind_ff,ind_fq) = 0;
    Tfhf(ind_ff)        = 0;
    fHf(ind_ff)         = 0;
    Qf(ind_ff,ind_fq)   = 0;
end






%% Elemental matrices
function [Df,Ef,Hf,Lf,Qf,fHf,TUhf,TQhf,Tfhf] = elementalMatricesCore(Df,Ef,Hf,Lf,Qf,TUhf,TQhf,Tfhf,Xel,tel,refEl,ifac,iFace)

% mesh data
global neq testcase coreflux axisym tau
global refElTor

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
B11 = zeros(nv); B12 = B11; B13 = B11; B14 = B11;
B21 = B11; B22 = B11; B23 = B11; B24 = B11;
B31 = B11; B32 = B11; B33 = B11; B34 = B11;
B41 = B11; B42 = B11; B43 = B11; B44 = B11;
fHf = zeros(neq*Nfp,1);
fH = zeros(nv,neq);

%  LOOP IN GAUSS POINTS
for igtor = 1:ngausstor
    for igpol = 1:ngausspol
        g = (igtor-1)*ngausspol+igpol;
        
        % Shape functions and derivatives at the current integration point
        Nfg = Nfi(g,:);
        NN = Nfg'*Nfg*dsurf(g);
        
        B11 = B11 + tau*NN;
        fH = fH + tau*Nfg'*coreflux'*dsurf(g);
        
    end
end

B22 = B11;
B33 = B11;
B44 = B11;
% expand the matrices
E_loc = expandMatrixCv(B11, B12, B13,B14,B21,B22,B23,B24,B31,B32,B33,B34,B41,B42,B43,B44);

% elemental assembly
Df(ind_ff,ind_fe)   = 0;
Hf(ind_ff,ind_ff)   = 0;
TUhf(ind_ff,ind_ff) = 0;
TQhf(ind_ff,ind_fq) = 0;
Tfhf(ind_ff)        = 0;

if (testcase.n==51 || testcase.n==52 || testcase.n==53 || testcase.n==54 || testcase.n==55)
    Ef(ind_ff,ind_ff) = 0;
    Qf(ind_ff,ind_fq) = 0;
else
    
    Ef(ind_ff,ind_ff) = E_loc;
    Lf(ind_ff,ind_fq) = 0;
    Qf(ind_ff,ind_fq) = 0;
    fHf(ind_ff) = -col(fH');
end






%% Elemental matrices
function [Df,Ef,Hf,Lf,Qf,fHf,TUhf,TQhf,Tfhf] = elementalMatricesDirichletWeak(Df,Ef,Hf,Lf,Qf,Xel,tel,refEl,ifac,iFace)

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


% E_loc(E_loc==0) = 1e-8;

% elemental assembly
Df(ind_ff,ind_fe) = 0;
Hf(ind_ff,ind_ff) = 0;
Ef(ind_ff,ind_ff) = E_loc;
Lf(ind_ff,ind_fq) = 0;
Qf(ind_ff,ind_fq) = 0;
fHf(ind_ff) = -col(fH');
TUhf(ind_ff,ind_ff) = 0;
TQhf(ind_ff,ind_fq) = 0;
Tfhf(ind_ff)        = 0;

% Df(ind_ff,ind_fe) = 1e-8;






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


function res = expandMatrixTQ(C11,C12,C13,C14,C15,C16,C17,C18,C19,C110,C111,C112,...
                              C21,C22,C23,C24,C25,C26,C27,C28,C29,C210,C211,C212,...
                              C31,C32,C33,C34,C35,C36,C37,C38,C39,C310,C311,C312,...
                              C41,C42,C43,C44,C45,C46,C47,C48,C49,C410,C411,C412)


res = zeros([size(C11) 4 12]);
res(:,:,[1 5 9  13 17 21 25 29 33 37 41 45 ...
         2 6 10 14 18 22 26 30 34 38 42 46 ...
         3 7 11 15 19 23 27 31 35 39 43 47 ...
         4 8 12 16 20 24 28 32 36 40 44 48 ]) = ...
    cat(3,C11,C12,C13,C14,C15,C16,C17,C18,C19,C110,C111,C112,...
          C21,C22,C23,C24,C25,C26,C27,C28,C29,C210,C211,C212,...
          C31,C32,C33,C34,C35,C36,C37,C38,C39,C310,C311,C312,...
          C41,C42,C43,C44,C45,C46,C47,C48,C49,C410,C411,C412);
res = permute(res, [3 1 4 2]);
res = reshape(res, 4*size(C11,1),12*size(C11,2));



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
