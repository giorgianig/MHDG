function [ TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf,TUhdir] =...
    hdg_ParallDiffusionEnergyMatrix_3D(X,T,F,flipFace,refEl,solQ,solU,solUhat,F_dir)

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
Ndim = 3;
% nv = size(refEl.NodesCoord1d,1);   % number of face nodes for the velocity
% Ndim = 3;
tdiv = linspace(0,theta,ntor+1);
nDirFaces =  sum(sum(F_dir));


% allocation and initialization
TU       = zeros(neq*Nv,neq*Nv,Ne);
TUh     = zeros(neq*Nv,neq*Nfp,Ne);
TQ       = zeros(neq*Nv,neq*Nv*Ndim,Ne);
TQh     = zeros(neq*Nv,Ndim*neq*Nv,Ne);
TUhf    = zeros(neq*Nfp,neq*Nfp,Ne);
TQhf    = zeros(neq*Nfp,neq*Ndim*Nv,Ne);
Tf        = zeros(neq*Nv,Ne);
Tfh      = zeros(neq*Nv,Ne);
Tfhf     = zeros(neq*Nfp,Ne);
TUhdir = zeros(neq*Nv,Ne);


% compute permutation for flipping faces
perm =  setperm3d(Np1dPol,Np1dTor,neq);
ind_dim = neq*[Np2d,repmat(Np1dPol*Np1dTor,1,nf),Np2d];
ind_sta = [1, 1+cumsum(ind_dim)];
ind_Uhat = zeros(Nfp*neq,1);


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
        indQ = (iElem-1)*neq*Nv*Ndim + (1:neq*Nv*Ndim);
        indU = (iElem-1)*neq*Nv + (1:neq*Nv);
        delta = 1+(itor-1)*(N2d*Np2d+(Nf-nDirFaces)*Nfl)*neq+neq*[(iel-1)*Np2d,N2d*Np2d+...
            (Fe-1)*Nfl,N2d*Np2d+(Nf-nDirFaces)*Nfl+(iel-1)*Np2d];
        if itor==ntor
            delta(end) = 1+(iel-1)*Np2d*neq;
        end
        for iface = 1:nf+2
            ind_loc = ind_sta(iface)+(0:ind_dim(iface)-1);
            ind_Uhat(ind_loc) = delta(iface)+(0:ind_dim(iface)-1);
        end
        solQ_e      = solQ(indQ);
        solU_e      = solU(indU);
        solUhat_e = solUhat(ind_Uhat);
        
        % elemental matrices
        [TUe,TQe,TUhe,TQhe,TUhfe,TQhfe,TUhdire,Tfe,Tfhe,Tfhfe] = elementalMatrices(Xe,tel,refEl,solQ_e,solU_e,solUhat_e,aux_dir,flipFace_e,iElem,Fe);
        
        %  TQ_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/TQ.txt');
        %  TQh_fort = readMatTxt('/home/giorgio/Dropbox/Fortran/MHDG/test/test_Energy/TQh.txt');
        %
        % max(abs(col(TQe)-TQ_fort(:)))
        % max(abs(col(TQhe)-TQh_fort(:)))
        % stop
        
        
        for iface = 1:nf
            if flipFace_e(iface)
                ind_v_L = Np2d*neq + (iface-1)*Nfl*neq+ (1:Nfl*neq);
                TUhfe(ind_v_L,:) = TUhfe(ind_v_L(perm),:);
                TUhfe(:,ind_v_L) = TUhfe(:,ind_v_L(perm));
                TQhfe(ind_v_L,:) = TQhfe(ind_v_L(perm),:);
                Tfhfe(ind_v_L)     = Tfhfe(ind_v_L(perm));
            end
        end
        
        
        
        
        
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
end

%% Elemental matrices
function [TU,TQ,TUh,TQh,TUhf,TQhf,TUhdir,Tf,Tfh,Tfhf] = elementalMatrices(Xe,tel,refEl,solQ,solU,solU_hat,aux_dir,flipFace_e,iElem,Fe)

% mesh data
global neq axisym testcase diff_pari diff_pare epn Magnetic Mref decoupleTiTe refElTor

htheta   = tel(end)-tel(1);

coefi = (2/(3*Mref) )^(1+epn)*diff_pari;
coefe = (2/(3*Mref) )^(1+epn)*diff_pare;

nf      = size(refEl.faceNodes,1);
Np1dPol = size(refEl.NodesCoord1d,1);          % Number of 1d points for poloidal line
Np1dTor = size(refElTor.NodesCoord1d,1);     % Number of 1d points for toroidal line
Np2d    = size(refEl.NodesCoord,1);
Nv      = Np1dTor*Np2d;
Nfl     = Np1dPol*Np1dTor;
Nfp     = Np2d*2+nf*Np1dPol*Np1dTor;
Ndim      = 3;


% Initialization
TUh     = zeros(neq*Nv,neq*Nfp);
TQh     = zeros(neq*Nv,Ndim*neq*Nv);
TUhf    = zeros(neq*Nfp,neq*Nfp);
TUhdir = zeros(neq*Nv,1);
Tfh      = zeros(neq*Nv,1);
TQhf    = zeros(neq*Nfp,neq*Ndim*Nv);
Tfhf     = zeros(neq*Nfp,1);


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

% solution at Gauss points
solQ    = reshape(solQ,neq*Ndim,numel(solQ)/neq/Ndim)';
solU    = reshape(solU,neq,numel(solU)/neq)';

%% VOLUME COMPUTATIONS
% Initialization
zz = zeros(Nv);
TU11 = zz; TU12 = zz; TU13 = zz; TU14 = zz;
TU21 = zz; TU22 = zz; TU23 = zz; TU24 = zz;
TU31 = zz; TU32 = zz; TU33 = zz; TU34 = zz;
TU41 = zz; TU42 = zz; TU43 = zz; TU44 = zz;

TQ1_01 = zz; TQ1_02 = zz; TQ1_03 = zz; TQ1_04 = zz; TQ1_05 = zz; TQ1_06 = zz; TQ1_07 = zz; TQ1_08 = zz; TQ1_09 = zz; TQ1_10 = zz; TQ1_11 = zz; TQ1_12 = zz;
TQ2_01 = zz; TQ2_02 = zz; TQ2_03 = zz; TQ2_04 = zz; TQ2_05 = zz; TQ2_06 = zz; TQ2_07 = zz; TQ2_08 = zz; TQ2_09 = zz; TQ2_10 = zz; TQ2_11 = zz; TQ2_12 = zz;
TQ3_01 = zz; TQ3_02 = zz; TQ3_03 = zz; TQ3_04 = zz; TQ3_05 = zz; TQ3_06 = zz; TQ3_07 = zz; TQ3_08 = zz; TQ3_09 = zz; TQ3_10 = zz; TQ3_11 = zz; TQ3_12 = zz;
TQ4_01 = zz; TQ4_02 = zz; TQ4_03 = zz; TQ4_04 = zz; TQ4_05 = zz; TQ4_06 = zz; TQ4_07 = zz; TQ4_08 = zz; TQ4_09 = zz; TQ4_10 = zz; TQ4_11 = zz; TQ4_12 = zz;
Tf1     = zeros(Nv,1); Tf2 = Tf1; Tf3 = Tf1; Tf4 = Tf1;

% aux1 = Tf1;  aux2 = Tf1; aux3 = Tf1; aux4 = Tf1; aux5 = Tf1; aux6  =Tf1;
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
        
        if (testcase.n >= 50 && testcase.n<60)
            b = [Magnetic.bx(g,iel), Magnetic.by(g,iel)];
        else
            b = defineMagneticField_3D([xg,yg],tg);
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
        solU_g = Ni*solU; % ng x nvar
        solQ_g = Ni*solQ; % ng x nvar*Ndim

        % Integration weight
        dvolu=ipw2d(igpol)*det(J)*dvolu1d;
        
        if axisym
            dvolu=dvolu*xg;
        end           
    
        % Compute G^T^(k-1)
        G = reshape(solQ_g,Ndim,neq);

        % Compute V(U^(k-1))
        Vi = computeVi(solU_g);
        Ve = computeVe(solU_g);

        %     **** parallel E terms *****************
        % Compute W(U^(k-1))
        W = compute_W(solU_g);

        % Compute dW_dU (k-1)
        dW_dU = compute_dW_dU(solU_g);
        %    **************************************

        % Compute dV_dU (k-1)
        dVi_dU = compute_dVi_dU(solU_g);
        dVe_dU = compute_dVe_dU(solU_g);

        % Compute Alpha(U^(k-1))
        Alphai = computeAlphai(solU_g);
        Alphae = computeAlphae(solU_g);

        % Compute dAlpha/dU^(k-1)
        dAlphai_dU = compute_dAlphai_dU(solU_g);
        dAlphae_dU = compute_dAlphae_dU(solU_g);

        % Compute ss(U^(k-1))
        ss = compute_s(solU_g);

        % Compute ds_dU (k-1)
        ds_dU = compute_ds_dU(solU_g);

        if decoupleTiTe
            W = zeros(size(W));
            dW_dU = zeros(size(dW_dU));
            ss = zeros(size(ss));
            ds_dU = zeros(size(ds_dU));
        end
    
    %
        NNf =  (Nx'*b(1)+Ny'*b(2)+Nt'*b(3))*dvolu;
        NN  =  (Nx'*b(1)+Ny'*b(2)+Nt'*b(3))*Ni*dvolu;
        NNi  =  Ni'*Ni*dvolu;
    
        gammai = dot(G*Vi,b);    % scalar
        gammae = dot(G*Ve,b);    % scalar
        Taui = G*dVi_dU;      % 2x4
        Taue= G*dVe_dU;      % 2x4
        Zet = G*dW_dU;        % 2x4
    
 

        
        % Contribution of the current integration point to the elemental matrix
        TU31 = TU31 + coefi*(gammai*dAlphai_dU(1)+Alphai*(Taui(1,1)*b(1)+Taui(2,1)*b(2)+Taui(3,1)*b(3)))*NN + (Zet(1,1)*b(1)+Zet(2,1)*b(2)+Zet(3,1)*b(3)+ds_dU(1))*NNi;
        TU32 = TU32 + coefi*(gammai*dAlphai_dU(2)+Alphai*(Taui(1,2)*b(1)+Taui(2,2)*b(2)+Taui(3,2)*b(3)))*NN + (Zet(1,2)*b(1)+Zet(2,2)*b(2)+Zet(3,2)*b(3)+ds_dU(2))*NNi;
        TU33 = TU33 + coefi*(gammai*dAlphai_dU(3)+Alphai*(Taui(1,3)*b(1)+Taui(2,3)*b(2)+Taui(3,3)*b(3)))*NN + (Zet(1,3)*b(1)+Zet(2,3)*b(2)+Zet(3,3)*b(3)+ds_dU(3))*NNi;
        TU34 = TU34 + coefi*(gammai*dAlphai_dU(4)+Alphai*(Taui(1,4)*b(1)+Taui(2,4)*b(2)+Taui(3,4)*b(3)))*NN + (Zet(1,4)*b(1)+Zet(2,4)*b(2)+Zet(3,4)*b(3)+ds_dU(4))*NNi;
        TU41 = TU41 + coefe*(gammae*dAlphae_dU(1)+Alphae*(Taue(1,1)*b(1)+Taue(2,1)*b(2)+Taue(3,1)*b(3)))*NN - (Zet(1,1)*b(1)+Zet(2,1)*b(2)+Zet(3,1)*b(3)+ds_dU(1))*NNi;
        TU42 = TU42 + coefe*(gammae*dAlphae_dU(2)+Alphae*(Taue(1,2)*b(1)+Taue(2,2)*b(2)+Taue(3,2)*b(3)))*NN - (Zet(1,2)*b(1)+Zet(2,2)*b(2)+Zet(3,2)*b(3)+ds_dU(2))*NNi;
        TU43 = TU43 + coefe*(gammae*dAlphae_dU(3)+Alphae*(Taue(1,3)*b(1)+Taue(2,3)*b(2)+Taue(3,3)*b(3)))*NN - (Zet(1,3)*b(1)+Zet(2,3)*b(2)+Zet(3,3)*b(3)+ds_dU(3))*NNi;
        TU44 = TU44 + coefe*(gammae*dAlphae_dU(4)+Alphae*(Taue(1,4)*b(1)+Taue(2,4)*b(2)+Taue(3,4)*b(3)))*NN - (Zet(1,4)*b(1)+Zet(2,4)*b(2)+Zet(3,4)*b(3)+ds_dU(4))*NNi;
    

%         if iElem==6
%    Alphai
%    coefi
%    Alphae
%    coefe
%     
%         end

        TQ3_01 = TQ3_01 + (coefi*Alphai*Vi(1)*NN + W(1)*NNi)*b(1);
        TQ3_02 = TQ3_02 + (coefi*Alphai*Vi(1)*NN + W(1)*NNi)*b(2);
        TQ3_03 = TQ3_03 + (coefi*Alphai*Vi(1)*NN + W(1)*NNi)*b(3);
        TQ3_04 = TQ3_04 + (coefi*Alphai*Vi(2)*NN + W(2)*NNi)*b(1);
        TQ3_05 = TQ3_05 + (coefi*Alphai*Vi(2)*NN + W(2)*NNi)*b(2);
        TQ3_06 = TQ3_06 + (coefi*Alphai*Vi(2)*NN + W(2)*NNi)*b(3);
        TQ3_07 = TQ3_07 + (coefi*Alphai*Vi(3)*NN + W(3)*NNi)*b(1);
        TQ3_08 = TQ3_08 + (coefi*Alphai*Vi(3)*NN + W(3)*NNi)*b(2);
        TQ3_09 = TQ3_09 + (coefi*Alphai*Vi(3)*NN + W(3)*NNi)*b(3);
        TQ3_10 = TQ3_10 + (coefi*Alphai*Vi(4)*NN + W(4)*NNi)*b(1);
        TQ3_11 = TQ3_11 + (coefi*Alphai*Vi(4)*NN + W(4)*NNi)*b(2);
        TQ3_12 = TQ3_12 + (coefi*Alphai*Vi(4)*NN + W(4)*NNi)*b(3);
        TQ4_01 = TQ4_01 + (coefe*Alphae*Ve(1)*NN - W(1)*NNi)*b(1);
        TQ4_02 = TQ4_02 + (coefe*Alphae*Ve(1)*NN - W(1)*NNi)*b(2);
        TQ4_03 = TQ4_03 + (coefe*Alphae*Ve(1)*NN - W(1)*NNi)*b(3);
        TQ4_04 = TQ4_04 + (coefe*Alphae*Ve(2)*NN - W(2)*NNi)*b(1);
        TQ4_05 = TQ4_05 + (coefe*Alphae*Ve(2)*NN - W(2)*NNi)*b(2);
        TQ4_06 = TQ4_06 + (coefe*Alphae*Ve(2)*NN - W(2)*NNi)*b(3);
        TQ4_07 = TQ4_07 + (coefe*Alphae*Ve(3)*NN - W(3)*NNi)*b(1);
        TQ4_08 = TQ4_08 + (coefe*Alphae*Ve(3)*NN - W(3)*NNi)*b(2);
        TQ4_09 = TQ4_09 + (coefe*Alphae*Ve(3)*NN - W(3)*NNi)*b(3);
        TQ4_10 = TQ4_10 + (coefe*Alphae*Ve(4)*NN - W(4)*NNi)*b(1);
        TQ4_11 = TQ4_11 + (coefe*Alphae*Ve(4)*NN - W(4)*NNi)*b(2);
        TQ4_12 = TQ4_12 + (coefe*Alphae*Ve(4)*NN - W(4)*NNi)*b(3);
            
    
        Tf3    = Tf3 + coefi*Alphai*( dot(Taui(1,:),solU_g)*b(1)+dot(Taui(2,:),solU_g)*b(2)+dot(Taui(3,:),solU_g)*b(3)   )*NNf + ss*Ni'*dvolu;
        Tf4    = Tf4 + coefe*Alphae*( dot(Taue(1,:),solU_g)*b(1)+dot(Taue(2,:),solU_g)*b(2)+dot(Taue(3,:),solU_g)*b(3)   )*NNf - ss*Ni'*dvolu;
    end
end


% if iElem==6
%     stop
% end

% expand the matrices
TU = expandMatrixCv(TU11,TU12,TU13,TU14,TU21,TU22,TU23,TU24,TU31,TU32,TU33,TU34,TU41,TU42,TU43,TU44);
TQ = expandMatrixTQ(TQ1_01,TQ1_02,TQ1_03,TQ1_04,TQ1_05,TQ1_06,TQ1_07,TQ1_08,TQ1_09,TQ1_10,TQ1_11,TQ1_12, ...
                    TQ2_01,TQ2_02,TQ2_03,TQ2_04,TQ2_05,TQ2_06,TQ2_07,TQ2_08,TQ2_09,TQ2_10,TQ2_11,TQ2_12, ...
                    TQ3_01,TQ3_02,TQ3_03,TQ3_04,TQ3_05,TQ3_06,TQ3_07,TQ3_08,TQ3_09,TQ3_10,TQ3_11,TQ3_12, ...
                    TQ4_01,TQ4_02,TQ4_03,TQ4_04,TQ4_05,TQ4_06,TQ4_07,TQ4_08,TQ4_09,TQ4_10,TQ4_11,TQ4_12);
Tf = col(transpose(cat(2,Tf1,Tf2,Tf3,Tf4)));


%% FACES COMPUTATIONS:

% compute permutation for flipping faces
perm =  setperm3d(Np1dPol,Np1dTor,neq);

for iface = 1:(nf+2)
    % face nodes
    if iface==1
        nv = Np2d;
        ind_ff = 1:Np2d*neq;                              % assembly face-face
        ind_fe = 1:Np2d*neq;                              % assembly face-element (for var)
        ind_ffe = 1:Np2d;
        ind_fq = 1:Ndim*Np2d*neq;                        % assembly face-element (for grad)
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
        ind_ffe   = Np2d*(Np1dTor-1) + (1:Np2d);
        ind_fq   = Np2d*(Np1dTor-1)*Ndim*neq +(1:Ndim*Np2d*neq);
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
        ind_ffe     = nodesv;
        ind_fq      = col(bsxfun(@plus,(nodesv-1)*Ndim*neq,(1:Ndim*neq)'));
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
        solU_g = analyticalSolution_3D(xyfg,thetafg);
    else
        solUh_f = solU_hat(ind_ff);
        if (iface>1 && iface<nf+2)
            if flipFace_e(iface-1)
                solUh_f = solUh_f(perm);
            end
        end
        solU_g = Nfi*transpose(reshape(solUh_f,neq,nv));
    end
    
    
     if (testcase.n >= 50 && testcase.n<60)
        if flipFace(iface)
            igauss = ngauss_f-g+1;
        else
            igauss = g;
        end
        b_g = [Magnetic.bxfaces(igauss,iFace), Magnetic.byfaces(igauss,iFace)];
    else
        b_g = defineMagneticField_3D([xyfg(:,1),xyfg(:,2)],thetafg(:));
     end   
    
    
    solQf   = solQ(ind_ffe,:);
    solQ_g = Nfi*solQf;
    
    
    % Initialization
    zz = zeros(nv);
    
    TUh11 = zz; TUh12 = zz; TUh13 = zz; TUh14 = zz;
    TUh21 = zz; TUh22 = zz; TUh23 = zz; TUh24 = zz;
    TUh31 = zz; TUh32 = zz; TUh33 = zz; TUh34 = zz;
    TUh41 = zz; TUh42 = zz; TUh43 = zz; TUh44 = zz;
    
   
    TQh1_01 = zz; TQh1_02 = zz; TQh1_03 = zz; TQh1_04 = zz; TQh1_05 = zz; TQh1_06 = zz; TQh1_07 = zz; TQh1_08 = zz; TQh1_09 = zz; TQh1_10 = zz; TQh1_11 = zz; TQh1_12 = zz;
    TQh2_01 = zz; TQh2_02 = zz; TQh2_03 = zz; TQh2_04 = zz; TQh2_05 = zz; TQh2_06 = zz; TQh2_07 = zz; TQh2_08 = zz; TQh2_09 = zz; TQh2_10 = zz; TQh2_11 = zz; TQh2_12 = zz;
    TQh3_01 = zz; TQh3_02 = zz; TQh3_03 = zz; TQh3_04 = zz; TQh3_05 = zz; TQh3_06 = zz; TQh3_07 = zz; TQh3_08 = zz; TQh3_09 = zz; TQh3_10 = zz; TQh3_11 = zz; TQh3_12 = zz;
    TQh4_01 = zz; TQh4_02 = zz; TQh4_03 = zz; TQh4_04 = zz; TQh4_05 = zz; TQh4_06 = zz; TQh4_07 = zz; TQh4_08 = zz; TQh4_09 = zz; TQh4_10 = zz; TQh4_11 = zz; TQh4_12 = zz;
    
    Tfh1          = zeros(nv,1); Tfh2 = Tfh1; Tfh3 = Tfh1; Tfh4 = Tfh1;
    TUhdir1     = zeros(nv,1); TUhdir2 = Tfh1; TUhdir3 = Tfh1; TUhdir4 = Tfh1;
    
    
    
    %  LOOP IN GAUSS POINTS
    for igtor = 1:ngausstor
        for igpol = 1:ngausspol 
            
%             g = (igpol-1)*ngausstor+igtor;
            g = (igtor-1)*ngausspol+igpol;
            b = b_g(g,:);
            n = n_g(g,:);
            
              % Velocity shape functions and derivatives at the current integration point
            Nf_g = Nfi(g,:);      
        
            % Compute G^T^(k-1)
            G = reshape(solQ_g(g,:),Ndim,neq);
            
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
            NNf =  Nf_g'*dot(b,n)*dsurf(g);
            NN  =  Nf_g'*dot(b,n)*Nf_g*dsurf(g);
            
            gammai = dot(G*Vi,b);    % scalar
            gammae = dot(G*Ve,b);    % scalar
            Taui = G*dVi_dU;      % 2x3
            Taue= G*dVe_dU;      % 2x3
            
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
        
 
            Tfh3    = Tfh3 + coefi*Alphai*sum( dot(Taui(1,:),solU_g(g,:))*b(1)+dot(Taui(2,:),solU_g(g,:))*b(2)+dot(Taui(3,:),solU_g(g,:))*b(3)    )*NNf;
            Tfh4    = Tfh4 + coefe*Alphae*sum( dot(Taue(1,:),solU_g(g,:))*b(1)+dot(Taue(2,:),solU_g(g,:))*b(2)+dot(Taue(3,:),solU_g(g,:))*b(3)   )*NNf;


            TUhdir3 = TUhdir3 + coefi*NNf*( (gammai*dAlphai_dU(1)+Alphai*(Taui(1,1)*b(1)+Taui(2,1)*b(2)+Taui(3,1)*b(3)))*solU_g(g,1)+...
                                            (gammai*dAlphai_dU(2)+Alphai*(Taui(1,2)*b(1)+Taui(2,2)*b(2)+Taui(3,2)*b(3)))*solU_g(g,2)+...
                                            (gammai*dAlphai_dU(3)+Alphai*(Taui(1,3)*b(1)+Taui(2,3)*b(2)+Taui(3,3)*b(3)))*solU_g(g,3)+...
                                            (gammai*dAlphai_dU(4)+Alphai*(Taui(1,4)*b(1)+Taui(2,4)*b(2)+Taui(3,4)*b(3)))*solU_g(g,4) );
            TUhdir4 = TUhdir4 +coefe* NNf*( (gammae*dAlphae_dU(1)+Alphae*(Taue(1,1)*b(1)+Taue(2,1)*b(2)+Taue(3,1)*b(3)))*solU_g(g,1)+...
                                            (gammae*dAlphae_dU(2)+Alphae*(Taue(1,2)*b(1)+Taue(2,2)*b(2)+Taue(3,2)*b(3)))*solU_g(g,2)+...
                                            (gammae*dAlphae_dU(3)+Alphae*(Taue(1,3)*b(1)+Taue(2,3)*b(2)+Taue(3,3)*b(3)))*solU_g(g,3)+...
                                            (gammae*dAlphae_dU(4)+Alphae*(Taue(1,4)*b(1)+Taue(2,4)*b(2)+Taue(3,4)*b(3)))*solU_g(g,4));
        
        end
    end
    
    % expand the matrices
  TUhloc = expandMatrixCv(TUh11,TUh12,TUh13,TUh14,TUh21,TUh22,TUh23,TUh24,TUh31,TUh32,TUh33,TUh34,TUh41,TUh42,TUh43,TUh44);  
  TQhloc = expandMatrixTQ(TQh1_01,TQh1_02,TQh1_03,TQh1_04,TQh1_05,TQh1_06,TQh1_07,TQh1_08,TQh1_09,TQh1_10,TQh1_11,TQh1_12, ...
                          TQh2_01,TQh2_02,TQh2_03,TQh2_04,TQh2_05,TQh2_06,TQh2_07,TQh2_08,TQh2_09,TQh2_10,TQh2_11,TQh2_12, ...
                          TQh3_01,TQh3_02,TQh3_03,TQh3_04,TQh3_05,TQh3_06,TQh3_07,TQh3_08,TQh3_09,TQh3_10,TQh3_11,TQh3_12, ...
                          TQh4_01,TQh4_02,TQh4_03,TQh4_04,TQh4_05,TQh4_06,TQh4_07,TQh4_08,TQh4_09,TQh4_10,TQh4_11,TQh4_12);
  
    
    % elemental assembly
    TUh(ind_fe,ind_ff)  = TUh(ind_fe,ind_ff) + ~isdir*TUhloc;
    TQh(ind_fe,ind_fq)  = TQh(ind_fe,ind_fq) + TQhloc;
    TUhf(ind_ff,ind_ff) = TUhf(ind_ff,ind_ff) + ~isdir*TUhloc;
    Tfh(ind_fe)         = Tfh(ind_fe) + col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));
    TQhf(ind_ff,ind_fq) = TQhf(ind_ff,ind_fq) + ~isdir*TQhloc;
    TUhdir(ind_fe)      = TUhdir(ind_fe) + isdir*col(transpose(cat(2,TUhdir1,TUhdir2,TUhdir3,TUhdir4)));
    Tfhf(ind_ff)        = Tfhf(ind_ff) + ~isdir*col(transpose(cat(2,Tfh1,Tfh2,Tfh3,Tfh4)));
    
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