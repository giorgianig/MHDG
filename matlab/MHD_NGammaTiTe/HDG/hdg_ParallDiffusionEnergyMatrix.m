function [ TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf,TUhdir] =...
    hdg_ParallDiffusionEnergyMatrix(X,T,F,flipFace,refEl,solQ,solU,solUhat,F_dir)

% mesh data
global neq
Ne = size(T,1);                     % number of elements
Nv = size(T,2);                     % number of element nodes for the velocity
nv = size(refEl.NodesCoord1d,1);   % number of face nodes for the velocity
nf  = size(refEl.faceNodes,1);    % number of faces in one element
nd = size(X,2);

% allocation and initialization
TU       = zeros(neq*Nv,neq*Nv,Ne);
TUh     = zeros(neq*Nv,neq*nf*nv,Ne);
TQ       = zeros(neq*Nv,neq*Nv*nd,Ne);
TQh     = zeros(neq*Nv,nd*neq*Nv,Ne);
TUhf    = zeros(neq*nf*nv,neq*nf*nv,Ne);
TQhf    = zeros(neq*nf*nv,neq*nd*Nv,Ne);
Tf        = zeros(neq*Nv,Ne);
Tfh      = zeros(neq*Nv,Ne);
Tfhf     = zeros(neq*nf*nv,Ne);
TUhdir = zeros(neq*Nv,Ne);


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
    indQ = (iElem-1)*neq*Nv*nd + (1:neq*Nv*nd);
    indU = (iElem-1)*neq*Nv + (1:neq*Nv);
    ind_Uhat =  bsxfun(@plus,(Fe-1)*neq*nv,(1:neq*nv)');
    solQ_e      = solQ(indQ);
    solU_e      = solU(indU);
    solUhat_e = solUhat(ind_Uhat);
    Te = T(iElem,:);
    Xe = X(Te,:);
    flipFace_e = flipFace(iElem,:);
    aux_dir = F_dir(iElem,:);
    
    % elemental matrices
    [TUe,TQe,TUhe,TQhe,TUhfe,TQhfe,TUhdire,Tfe,Tfhe,Tfhfe] = elementalMatrices(Xe,refEl,solQ_e,solU_e,solUhat_e,aux_dir,flipFace_e,iElem,Fe);
    
    %  TQ_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/TQ.txt');
    %  TQh_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/TQh.txt');
    %
    % max(abs(col(TQe)-TQ_fort(:)))
    % max(abs(col(TQhe)-TQh_fort(:)))
    % stop
    
    
    for iface = 1:nf
        if flipFace_e(iface)
            TUhfe(ind_v_L(iface,:),:) = TUhfe(ind_v_L(iface,perm),:);
            TUhfe(:,ind_v_L(iface,:)) = TUhfe(:,ind_v_L(iface,perm));
            TQhfe(ind_v_L(iface,:),:) = TQhfe(ind_v_L(iface,perm),:);
            Tfhfe(ind_v_L(iface,:))     = Tfhfe(ind_v_L(iface,perm));
        end
    end
    
    
  
    
% if iElem==144
%     stop
% end

    
    % store matrices
    TU(:,:,iElem)   = TUe;           % Local problem, volume contribution, U part
    TUh(:,:,iElem) = TUhe;         % Local problem, face contribution, U part
    TQ(:,:,iElem)   = TQe;           % Local problem, volume contribution, Q part
    TQh(:,:,iElem) = TQhe;         % Local problem, face contribution, Q part
    TUhf(:,:,iElem)= TUhfe;        % Global problem, U part
    TQhf(:,:,iElem)= TQhfe;        % Global problem, Q part
    Tf(:,iElem)       = Tfe;            % Local problem, volume contribution, RHS
    Tfh(:,iElem)     = Tfhe;          % Local problem, face contribution, RHS
    Tfhf(:,iElem)    = Tfhfe;         % Global problem, RHS
    TUhdir(:,iElem)= TUhdire;     % Local problem, Dirichlet
end

%% Elemental matrices
function [TU,TQ,TUh,TQh,TUhf,TQhf,TUhdir,Tf,Tfh,Tfhf] = elementalMatrices(Xe,refEl,solQ,solU,solU_hat,aux_dir,flipFace_e,iElem,Fe)

% mesh data
global neq axisym testcase diff_pari diff_pare epn Magnetic Mref decoupleTiTe

coefi = (2/(3*Mref) )^(1+epn)*diff_pari;
coefe = (2/(3*Mref) )^(1+epn)*diff_pare;

Nv = size(refEl.NodesCoord,1);
nv = size(refEl.NodesCoord1d,1);
faceNodesv = refEl.faceNodes;
nf = size(faceNodesv,1);
nd = size(Xe,2);

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
TUh     = zeros(neq*Nv,neq*nf*nv);
TQh     = zeros(neq*Nv,nd*neq*Nv);
TUhf    = zeros(neq*nf*nv,neq*nf*nv);
TUhdir = zeros(neq*Nv,1);
Tfh      = zeros(neq*Nv,1);
TQhf    = zeros(neq*nf*nv,neq*nd*Nv);
Tfhf     = zeros(neq*nf*nv,1);

% solution at Gauss points
solQ    = reshape(solQ,neq*nd,numel(solQ)/neq/nd)';
solQ_g = N*solQ; % ng x nvar*nd
solU    = reshape(solU,neq,numel(solU)/neq)';
solU_g = N*solU; % ng x nvar

%% VOLUME COMPUTATIONS
% Initialization
TU11 = zeros(Nv); TU12 = TU11; TU13 = TU11; TU14 = TU11;
TU21 = TU11; TU22 = TU11; TU23 = TU11; TU24 = TU11;
TU31 = TU11; TU32 = TU11; TU33 = TU11; TU34 = TU11;
TU41 = TU11; TU42 = TU11; TU43 = TU11; TU44 = TU11;

TQ11 = zeros(Nv); TQ12 = TQ11; TQ13 = TQ11; TQ14 = TQ11; TQ15 = TQ11; TQ16 = TQ11; TQ17 = TQ11; TQ18 = TQ11;
TQ21 = TQ11; TQ22 = TQ11; TQ23 = TQ11; TQ24 = TQ11; TQ25 = TQ11; TQ26 = TQ11; TQ27 = TQ11; TQ28 = TQ11;
TQ31 = TQ11; TQ32 = TQ11; TQ33 = TQ11; TQ34 = TQ11; TQ35 = TQ11; TQ36 = TQ11; TQ37 = TQ11; TQ38 = TQ11;
TQ41 = TQ11; TQ42 = TQ11; TQ43 = TQ11; TQ44 = TQ11; TQ45 = TQ11; TQ46 = TQ11; TQ47 = TQ11; TQ48 = TQ11;

Tf1     = zeros(Nv,1); Tf2 = Tf1; Tf3 = Tf1; Tf4 = Tf1;

% gauss point position
xg = N*xe;
yg = N*ye;

% magnetic field
if (testcase.n >= 50 && testcase.n<60)
    bgauss = [Magnetic.bx(:,iElem), Magnetic.by(:,iElem)];
else
    bgauss = defineMagneticField([xg,yg]);
end

% aux1 = Tf1;  aux2 = Tf1; aux3 = Tf1; aux4 = Tf1; aux5 = Tf1; aux6  =Tf1;
for g = 1:ngauss
    
    % Velocity shape functions and derivatives at the current integration point
    Nig = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    
    
    b = bgauss(g,:);
    
    
    % Jacobian of the element
    J = [Nxi_g*xe  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end
    
    % Integration weight
    dvolu=IPw(g)*det(J);
    
    if axisym
        dvolu=dvolu*xg(g);
    end
    
    % x and y derivatives
    invJ = inv(J);
    Nxg = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Nyg = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    
    
    % Compute G^T^(k-1)
    G = reshape(solQ_g(g,:),nd,neq);
    
    % Compute V(U^(k-1))
    Vi = computeVi(solU_g(g,:));
    Ve = computeVe(solU_g(g,:));
    
%     **** parallel E terms *****************
    % Compute W(U^(k-1))   
    W = compute_W(solU_g(g,:));
    
    % Compute dW_dU (k-1)
    dW_dU = compute_dW_dU(solU_g(g,:));
%    **************************************

    % Compute dV_dU (k-1)
    dVi_dU = compute_dVi_dU(solU_g(g,:));
    dVe_dU = compute_dVe_dU(solU_g(g,:));
    
    % Compute Alpha(U^(k-1))
    Alphai = computeAlphai(solU_g(g,:));
    Alphae = computeAlphae(solU_g(g,:));
    
    % Compute dAlpha/dU^(k-1)
    dAlphai_dU = compute_dAlphai_dU(solU_g(g,:));
    dAlphae_dU = compute_dAlphae_dU(solU_g(g,:));
    
    % Compute ss(U^(k-1))
    ss = compute_s(solU_g(g,:));
    
    % Compute ds_dU (k-1)
    ds_dU = compute_ds_dU(solU_g(g,:));
  
    if decoupleTiTe
        W = zeros(size(W));
        dW_dU = zeros(size(dW_dU));
        ss = zeros(size(ss));
        ds_dU = zeros(size(ds_dU));
    end
    
    
    %
    NNf =  (Nxg'*b(1)+Nyg'*b(2))*dvolu;
    NN  =  (Nxg'*b(1)+Nyg'*b(2))*Nig*dvolu;
    NNi  =  Nig'*Nig*dvolu;
    
    gammai = dot(G*Vi,b);    % scalar
    gammae = dot(G*Ve,b);    % scalar
    Taui = G*dVi_dU;      % 2x4
    Taue= G*dVe_dU;      % 2x4
    Zet = G*dW_dU;        % 2x4
    
    
    % Contribution of the current integration point to the elemental matrix
    TU31 = TU31 + coefi*(gammai*dAlphai_dU(1)+Alphai*(Taui(1,1)*b(1)+Taui(2,1)*b(2)))*NN + (Zet(1,1)*b(1)+Zet(2,1)*b(2)+ds_dU(1))*NNi;
    TU32 = TU32 + coefi*(gammai*dAlphai_dU(2)+Alphai*(Taui(1,2)*b(1)+Taui(2,2)*b(2)))*NN + (Zet(1,2)*b(1)+Zet(2,2)*b(2)+ds_dU(2))*NNi;
    TU33 = TU33 + coefi*(gammai*dAlphai_dU(3)+Alphai*(Taui(1,3)*b(1)+Taui(2,3)*b(2)))*NN + (Zet(1,3)*b(1)+Zet(2,3)*b(2)+ds_dU(3))*NNi;
    TU34 = TU34 + coefi*(gammai*dAlphai_dU(4)+Alphai*(Taui(1,4)*b(1)+Taui(2,4)*b(2)))*NN + (Zet(1,4)*b(1)+Zet(2,4)*b(2)+ds_dU(4))*NNi;
    TU41 = TU41 + coefe*(gammae*dAlphae_dU(1)+Alphae*(Taue(1,1)*b(1)+Taue(2,1)*b(2)))*NN - (Zet(1,1)*b(1)+Zet(2,1)*b(2)+ds_dU(1))*NNi;
    TU42 = TU42 + coefe*(gammae*dAlphae_dU(2)+Alphae*(Taue(1,2)*b(1)+Taue(2,2)*b(2)))*NN - (Zet(1,2)*b(1)+Zet(2,2)*b(2)+ds_dU(2))*NNi;
    TU43 = TU43 + coefe*(gammae*dAlphae_dU(3)+Alphae*(Taue(1,3)*b(1)+Taue(2,3)*b(2)))*NN - (Zet(1,3)*b(1)+Zet(2,3)*b(2)+ds_dU(3))*NNi;
    TU44 = TU44 + coefe*(gammae*dAlphae_dU(4)+Alphae*(Taue(1,4)*b(1)+Taue(2,4)*b(2)))*NN - (Zet(1,4)*b(1)+Zet(2,4)*b(2)+ds_dU(4))*NNi;
    
    
    TQ31 = TQ31 + (coefi*Alphai*Vi(1)*NN + W(1)*NNi)*b(1);
    TQ32 = TQ32 + (coefi*Alphai*Vi(1)*NN + W(1)*NNi)*b(2);
    TQ33 = TQ33 + (coefi*Alphai*Vi(2)*NN + W(2)*NNi)*b(1);
    TQ34 = TQ34 + (coefi*Alphai*Vi(2)*NN + W(2)*NNi)*b(2);
    TQ35 = TQ35 + (coefi*Alphai*Vi(3)*NN + W(3)*NNi)*b(1);
    TQ36 = TQ36 + (coefi*Alphai*Vi(3)*NN + W(3)*NNi)*b(2);
    TQ37 = TQ37 + (coefi*Alphai*Vi(4)*NN + W(4)*NNi)*b(1);
    TQ38 = TQ38 + (coefi*Alphai*Vi(4)*NN + W(4)*NNi)*b(2);
    TQ41 = TQ41 + (coefe*Alphae*Ve(1)*NN - W(1)*NNi)*b(1);
    TQ42 = TQ42 + (coefe*Alphae*Ve(1)*NN - W(1)*NNi)*b(2);
    TQ43 = TQ43 + (coefe*Alphae*Ve(2)*NN - W(2)*NNi)*b(1);
    TQ44 = TQ44 + (coefe*Alphae*Ve(2)*NN - W(2)*NNi)*b(2);
    TQ45 = TQ45 + (coefe*Alphae*Ve(3)*NN - W(3)*NNi)*b(1);
    TQ46 = TQ46 + (coefe*Alphae*Ve(3)*NN - W(3)*NNi)*b(2);
    TQ47 = TQ47 + (coefe*Alphae*Ve(4)*NN - W(4)*NNi)*b(1);
    TQ48 = TQ48 + (coefe*Alphae*Ve(4)*NN - W(4)*NNi)*b(2);
    
    
    Tf3    = Tf3 + coefi*Alphai*( dot(Taui(1,:),solU_g(g,:))*b(1)+dot(Taui(2,:),solU_g(g,:))*b(2)   )*NNf + ss*Nig'*dvolu;
    Tf4    = Tf4 + coefe*Alphae*( dot(Taue(1,:),solU_g(g,:))*b(1)+dot(Taue(2,:),solU_g(g,:))*b(2)   )*NNf - ss*Nig'*dvolu;
end


% expand the matrices
TU = expandMatrixCv(TU11,TU12,TU13,TU14,TU21,TU22,TU23,TU24,TU31,TU32,TU33,TU34,TU41,TU42,TU43,TU44);
TQ = expandMatrixTQ( TQ11,TQ12,TQ13,TQ14,TQ15,TQ16,TQ17,TQ18, ...
                                     TQ21,TQ22,TQ23,TQ24,TQ25,TQ26,TQ27,TQ28,...
                                     TQ31,TQ32,TQ33,TQ34,TQ35,TQ36,TQ37,TQ38,...
                                     TQ41,TQ42,TQ43,TQ44,TQ45,TQ46,TQ47,TQ48);
Tf = col(transpose(cat(2,Tf1,Tf2,Tf3,Tf4)));


%% FACES COMPUTATIONS:

% compute permutation for flipping faces
perm = setperm(neq*nv,neq);

ngauss_f = length(IPw_fv);
for iface = 1:nf
    
    % face nodes
    nodesv = faceNodesv(iface,:);
    
    % indices for local assembly
    ind_face_2 = (iface-1)*neq*nv + (1:neq*nv);
    ind2 = reshape(bsxfun(@plus,(nodesv-1)*neq,(1:neq)'),neq*nv,1); % assembly face to elem for U
    ind4 = reshape(bsxfun(@plus,(nodesv-1)*neq*nd,(1:neq*nd)'),neq*nd*nv,1); % assembly face to elem for Q
    
    xf = xe(nodesv);
    yf = ye(nodesv);
    solUh_f = solU_hat(:,iface);
    solQf    = solQ(nodesv,:);
    if flipFace_e(iface)
        solUh_f = solUh_f(perm);
    end
    solUh_f = transpose(reshape(solUh_f,neq,nv));
    
    % Gauss point position
    xyfg = N1dv*[xf yf];
    
    if aux_dir(iface)
        % exact velocity in the faces for Dirichlet boundary
        solU_g = analyticalSolution(xyfg);
    else
        solU_g = N1dv*solUh_f;
    end
    solQ_g = N1dv*solQf;
    
    
    % Initialization
    TUh11 = zeros(nv); TUh12 = TUh11; TUh13 = TUh11; TUh14 = TUh11;
    TUh21 = TUh11; TUh22 = TUh11; TUh23 = TUh11; TUh24 = TUh11;
    TUh31 = TUh11; TUh32 = TUh11; TUh33 = TUh11; TUh34 = TUh11;
    TUh41 = TUh11; TUh42 = TUh11; TUh43 = TUh11; TUh44 = TUh11;
    
    TQh11 = zeros(nv); TQh12 = TQh11; TQh13 = TQh11; TQh14 = TQh11; TQh15 = TQh11; TQh16 = TQh11; TQh17 = TQh11;  TQh18 = TQh11;
    TQh21 = TQh11; TQh22 = TQh11; TQh23 = TQh11; TQh24 = TQh11; TQh25 = TQh11; TQh26 = TQh11; TQh27 = TQh11;  TQh28 = TQh11;
    TQh31 = TQh11; TQh32 = TQh11; TQh33 = TQh11; TQh34 = TQh11; TQh35 = TQh11; TQh36 = TQh11; TQh37 = TQh11;  TQh38 = TQh11;
    TQh41 = TQh11; TQh42 = TQh11; TQh43 = TQh11; TQh44 = TQh11; TQh45 = TQh11; TQh46 = TQh11; TQh47 = TQh11;  TQh48 = TQh11;
    
    Tfh1          = zeros(nv,1); Tfh2 = Tfh1; Tfh3 = Tfh1; Tfh4 = Tfh1;
    TUhdir1     = zeros(nv,1); TUhdir2 = Tfh1; TUhdir3 = Tfh1; TUhdir4 = Tfh1;
    
    
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
        
        % magnetic field
        if (testcase.n >= 50 && testcase.n<60)
            if flipFace_e(iface)
                igauss = ngauss_f-g+1;
            else
                igauss = g;
            end
            b = [Magnetic.bxfaces(igauss,Fe(iface)), Magnetic.byfaces(igauss,Fe(iface))];
        else
            b = defineMagneticField([xyfg(g,1),xyfg(g,2)]);
        end
        
        
        % Compute G^T^(k-1)
        G = reshape(solQ_g(g,:),nd,neq);
        
        % Compute V(U^(k-1))
        Vi = computeVi(solU_g(g,:));
        Ve = computeVe(solU_g(g,:));
        
        % Compute dV_dU (k-1)
        dVi_dU = compute_dVi_dU(solU_g(g,:));
        dVe_dU = compute_dVe_dU(solU_g(g,:));
        
        % Compute Alpha(U^(k-1))
        Alphai = computeAlphai(solU_g(g,:));
        Alphae = computeAlphae(solU_g(g,:));
        
        % Compute dAlpha/dU^(k-1)
        dAlphai_dU = compute_dAlphai_dU(solU_g(g,:));
        dAlphae_dU = compute_dAlphae_dU(solU_g(g,:));
        
        %
        NNf =  Nfv_g'*dot(b,n_g)*dline;
        NN  =  Nfv_g'*dot(b,n_g)*Nfv_g*dline;
        
        gammai = dot(G*Vi,b);    % scalar
        gammae = dot(G*Ve,b);    % scalar
        Taui = G*dVi_dU;      % 2x3
        Taue= G*dVe_dU;      % 2x3
        
        
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

        
        Tfh3    = Tfh3 + coefi*Alphai*sum( dot(Taui(1,:),solU_g(g,:))*b(1)+dot(Taui(2,:),solU_g(g,:))*b(2)   )*NNf;
        Tfh4    = Tfh4 + coefe*Alphae*sum( dot(Taue(1,:),solU_g(g,:))*b(1)+dot(Taue(2,:),solU_g(g,:))*b(2)   )*NNf;
        
        
        TUhdir3 = TUhdir3 + coefi*NNf*( (gammai*dAlphai_dU(1)+Alphai*(Taui(1,1)*b(1)+Taui(2,1)*b(2)))*solU_g(g,1)+...
                                                            (gammai*dAlphai_dU(2)+Alphai*(Taui(1,2)*b(1)+Taui(2,2)*b(2)))*solU_g(g,2)+...
                                                            (gammai*dAlphai_dU(3)+Alphai*(Taui(1,3)*b(1)+Taui(2,3)*b(2)))*solU_g(g,3)+...
                                                            (gammai*dAlphai_dU(4)+Alphai*(Taui(1,4)*b(1)+Taui(2,4)*b(2)))*solU_g(g,4) );
        TUhdir4 = TUhdir4 +coefe* NNf*( (gammae*dAlphae_dU(1)+Alphae*(Taue(1,1)*b(1)+Taue(2,1)*b(2)))*solU_g(g,1)+...
                                                             (gammae*dAlphae_dU(2)+Alphae*(Taue(1,2)*b(1)+Taue(2,2)*b(2)))*solU_g(g,2)+...
                                                             (gammae*dAlphae_dU(3)+Alphae*(Taue(1,3)*b(1)+Taue(2,3)*b(2)))*solU_g(g,3)+...
                                                             (gammae*dAlphae_dU(4)+Alphae*(Taue(1,4)*b(1)+Taue(2,4)*b(2)))*solU_g(g,4));
        
    end
    
    % expand the matrices
    TUhloc = expandMatrixCv(TUh11,TUh12,TUh13,TUh14,TUh21,TUh22,TUh23,TUh24,TUh31,TUh32,TUh33,TUh34,TUh41,TUh42,TUh43,TUh44);
    TQhloc = expandMatrixTQ( TQh11,TQh12,TQh13,TQh14,TQh15,TQh16,TQh17,TQh18, ...
                                              TQh21,TQh22,TQh23,TQh24,TQh25,TQh26,TQh27,TQh28, ...
                                              TQh31,TQh32,TQh33,TQh34,TQh35,TQh36,TQh37,TQh38,...
                                              TQh41,TQh42,TQh43,TQh44,TQh45,TQh46,TQh47,TQh48);
    
    
    
    % elemental assembly
    TUh(ind2,ind_face_2)        = TUh(ind2,ind_face_2) + ~aux_dir(iface)*TUhloc;
    TQh(ind2,ind4)              = TQh(ind2,ind4) + TQhloc;
    TUhf(ind_face_2,ind_face_2) = TUhf(ind_face_2,ind_face_2) + ~aux_dir(iface)*TUhloc;
    Tfh(ind2)                   = Tfh(ind2) + col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));
    TQhf(ind_face_2,ind4)       = TQhf(ind_face_2,ind4) + ~aux_dir(iface)*TQhloc;
    TUhdir(ind2)                = TUhdir(ind2) + aux_dir(iface)*col(transpose(cat(2,TUhdir1,TUhdir2,TUhdir3,TUhdir4)));
    Tfhf(ind_face_2)            = Tfhf(ind_face_2) + ~aux_dir(iface)*col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));
    
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




function res = compute_W(U)
% res = [4x1]
global pcourr

res =  pcourr*2/3*[ 0; ...
            0; ...
            0;...
    U(2)/U(1)];




function res = compute_dW_dU(U)
% res = [4x4]
global pcourr
res = pcourr*2/3* [ 0,               0,         0,        0; ...
                   0,                0,         0,        0; ...
                   0,                0,         0,        0; ...
            -U(2)/U(1)^2,  1/U(1),    0,        0];
        
        
        
function res = compute_s(U)
% res: scalar
global Mref tie  constantExchangeTemp 
tol = 1e-5;
if constantExchangeTemp
    coef = 2/(3*Mref*tie);
    res = coef* ((U(3)-U(4))/U(1)-0.5*U(2)^2/U(1)^2);
else
    coef = 1/tie*(2/3/Mref)^(-0.5);
    U1=U(1);
    U4=U(4);
    if U1<0, U1 = tol; end
    if U4<0, U4 = tol; end
    res = coef* (U1^(2.5)/U4^1.5)*(U(4)-U(3)+0.5*(U(2)^2/U(1)));
end

function res = compute_ds_dU(U)
% res = [4x1]
global Mref tie constantExchangeTemp
tol = 1e-5;
if constantExchangeTemp
    coef = 2/(3*Mref*tie);
    res =coef* [  (U(4)-U(3))/U(1)^2+U(2)^2/U(1)^3; ...
                                      -U(2)/U(1)^2;   ...
                                           1/U(1); ...
                                          -1/U(1)];
else
    coef = 1/tie*(2/3/Mref)^(-0.5);
    U1=U(1);
    U4=U(4);
    if U1<0, U1 = tol; end
    if U4<0, U4 = tol; end    
    res = coef* [ 2.5*(U1/U4)^1.5*(U(4)-U(3)+0.5*(U(2)^2/U(1))) - 0.5*U1^0.5*U(2)^2/U4^1.5;...
                                                          U(2)*(U1/U4)^1.5;...
                                                      -U1^2.5/U4^1.5;...
                                      -1.5*(U1/U4)^2.5*(U(4)-U(3)+0.5*U(2)^2/U(1))+U1^2.5/U4^1.5];
end